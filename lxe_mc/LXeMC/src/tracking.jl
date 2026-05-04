"""
Stack-based Monte Carlo transport for photons and leptons in liquid xenon.

Design follows Geant4: a single LIFO particle stack, physics processes are
sampled at each interaction point, and secondaries are pushed onto the stack.
Charged-particle transport uses condensed history with explicit hard-brems
sampling.

# Public API

- [`simulate_event`](@ref): simulate one primary photon → list of energy deposits.
- [`cluster_deposits_in_z`](@ref): group deposits by z-proximity for SS/MS classification.
- [`is_single_site`](@ref): convenience wrapper.

# Geometry

Two geometry types are provided:
- [`InfiniteLXe`](@ref): no boundaries (all points inside).
- [`CylinderLXe`](@ref): finite cylinder along z.
"""

using Random


# =====================================================================
# Geometry
# =====================================================================

"""Abstract base for detector geometries."""
abstract type Geometry end

"""Infinite LXe volume (no boundaries)."""
struct InfiniteLXe <: Geometry end
is_inside(::InfiniteLXe, pos::Vector{Float64}) = true

"""
    CylinderLXe(radius_cm, half_height_cm)

Cylindrical LXe volume along z, centered at origin.
"""
struct CylinderLXe <: Geometry
    radius_cm::Float64
    half_height_cm::Float64
end
function is_inside(cyl::CylinderLXe, pos::Vector{Float64})::Bool
    pos[1]^2 + pos[2]^2 < cyl.radius_cm^2 && abs(pos[3]) < cyl.half_height_cm
end


# =====================================================================
# Track and stack
# =====================================================================

"""
    Track

A particle on the transport stack.

Fields:
- `kind`: `:gamma`, `:electron`, or `:positron`
- `energy`: kinetic energy [MeV] (for leptons) or total energy (for photons)
- `position`: 3-vector [cm]
- `direction`: unit 3-vector
- `parent_id`: unique id of parent track (-1 for primaries)
- `generation`: cascade depth (0 for primaries)
"""
struct Track
    kind::Symbol
    energy::Float64
    position::Vector{Float64}
    direction::Vector{Float64}
    parent_id::Int
    generation::Int
end


"""
    ParticleStack

LIFO stack of `Track` objects, backed by a `Vector`.
"""
mutable struct ParticleStack
    items::Vector{Track}
end
ParticleStack() = ParticleStack(Track[])
Base.push!(s::ParticleStack, t::Track) = push!(s.items, t)
Base.pop!(s::ParticleStack) = pop!(s.items)
Base.isempty(s::ParticleStack) = isempty(s.items)
Base.length(s::ParticleStack) = length(s.items)


"""
    Deposit

An energy deposit in the LXe volume.

Fields:
- `position`: 3-vector [cm]
- `energy`: deposited energy [MeV]
- `source`: origin label (`:electron`, `:positron`, `:photoelectric`,
  `:auger`, `:gamma_local`, `:brems_local`)
"""
struct Deposit
    position::Vector{Float64}
    energy::Float64
    source::Symbol
end


# =====================================================================
# Photon transport
# =====================================================================

"""
    transport_photon!(track::Track, geom::Geometry,
                      deposits::Vector{Deposit}, stack::ParticleStack,
                      nd::NISTData, cfg::SimConfig, rng::AbstractRNG)

Transport one gamma until escape, absorption, or energy falls below cutoff.

Pushes secondaries (electrons, positrons, photons) onto `stack` and
appends energy deposits to `deposits`. The photon is consumed by
pair production or photoelectric absorption, or continues after
Compton scattering with reduced energy.
"""
function transport_photon!(track::Track, geom::Geometry,
                           deposits::Vector{Deposit}, stack::ParticleStack,
                           nd::NISTData, cfg::SimConfig, rng::AbstractRNG)
    pos = copy(track.position)
    dir = copy(track.direction)
    E = track.energy
    tid = objectid(track)
    gen = track.generation

    while E >= cfg.Egamma_cut
        # Cross sections per atom
        sC  = sigma_compton_NIST(nd, E; per_atom=true)
        sP  = sigma_pair_NIST(nd, E; per_atom=true)
        sPh = sigma_phot_NIST(nd, E; per_atom=true)
        s_tot = sC + sP + sPh
        Σ_tot = cfg.n_atom * s_tot

        # Sample distance to next interaction
        s = sample_distance(Σ_tot, rng)
        pos .= pos .+ dir .* s
        is_inside(geom, pos) || return  # escapes

        # Sample which process
        proc = sample_process(sC/s_tot, sP/s_tot, sPh/s_tot, rng)

        if proc === :compton
            Egp, cos_t = sample_compton(E, cfg, rng)
            ϕ = 2π * rand(rng)

            # Push recoil electron
            n_e = compton_electron_direction(cos_t, ϕ, E, dir, cfg)
            T_e = E - Egp
            push!(stack, Track(:electron, T_e, copy(pos), n_e,
                               Int(tid % typemax(Int)), gen + 1))

            # Update scattered photon direction
            sin_t = sqrt(max(0.0, 1.0 - cos_t^2))
            local_vec = Float64[sin_t*cos(ϕ), sin_t*sin(ϕ), cos_t]
            dir = rotate_to_global(local_vec, dir)
            E = Egp

        elseif proc === :pair
            eps = sample_pair(E, cfg, rng)
            E_pos_total = eps * E
            E_ele_total = (1.0 - eps) * E
            T_pos = E_pos_total - cfg.me
            T_ele = E_ele_total - cfg.me

            θ_pos = pair_polar_angle(E_pos_total, cfg, rng)
            θ_ele = pair_polar_angle(E_ele_total, cfg, rng)
            ϕ = 2π * rand(rng)

            # Positron and electron, coplanar with incoming photon
            for (sign_val, θ, T, kind) in [(1.0, θ_pos, T_pos, :positron),
                                            (-1.0, θ_ele, T_ele, :electron)]
                ϕ_lep = ϕ + (sign_val < 0.0 ? π : 0.0)
                local_vec = Float64[sin(θ)*cos(ϕ_lep), sin(θ)*sin(ϕ_lep), cos(θ)]
                d_lep = rotate_to_global(local_vec, dir)
                if T > 0.0
                    push!(stack, Track(kind, T, copy(pos), d_lep,
                                       Int(tid % typemax(Int)), gen + 1))
                end
            end
            return  # photon consumed

        elseif proc === :photoelectric
            if E < cfg.EK
                # Below K-edge: deposit everything locally
                push!(deposits, Deposit(copy(pos), E, :photoelectric))
                return
            end

            # Above K-edge: deposit EK locally (relaxation cascade is
            # always local — fluorescence mfp ~0.05 mm, Auger range < 1 µm)
            push!(deposits, Deposit(copy(pos), cfg.EK, :photoelectric))

            # Photoelectron carries the rest
            T_e = E - cfg.EK
            if T_e > cfg.Te_cut
                θ_e = sample_photoelectron_angle(T_e, cfg, rng)
                ϕ_e = 2π * rand(rng)
                local_vec = Float64[sin(θ_e)*cos(ϕ_e), sin(θ_e)*sin(ϕ_e), cos(θ_e)]
                d_e = rotate_to_global(local_vec, dir)
                push!(stack, Track(:electron, T_e, copy(pos), d_e,
                                   Int(tid % typemax(Int)), gen + 1))
            else
                push!(deposits, Deposit(copy(pos), T_e, :photoelectric))
            end
            return  # photon consumed
        end
    end

    # Below cutoff: deposit residual energy locally
    push!(deposits, Deposit(pos, E, :gamma_local))
end


# =====================================================================
# Lepton transport (electron and positron)
# =====================================================================

"""
    transport_lepton!(track::Track, geom::Geometry,
                      deposits::Vector{Deposit}, stack::ParticleStack,
                      nd::NISTData, cfg::SimConfig, rng::AbstractRNG)

Transport an electron or positron through LXe using condensed-history
stepping with continuous collisional energy loss and explicit discrete
bremsstrahlung sampling.

**Stepping**: fixed step size `ds_step` (default 0.5 mm), clamped to
not overshoot the remaining kinetic energy. Collisional energy loss is
deposited at the midpoint of each step.

**Bremsstrahlung**: at each step, the probability of emitting a hard photon
(k > k_min) is P = n_atom × σ_brems(T) × ds, where σ_brems is looked up
from a pre-tabulated log-log table (no per-step numerical integration).
If triggered, the photon energy is sampled from the BH-Tsai spectrum and
pushed onto the stack.

**End of range**: residual energy is deposited locally. Positrons
annihilate at rest, producing two back-to-back 511 keV photons.
"""
function transport_lepton!(track::Track, geom::Geometry,
                           deposits::Vector{Deposit}, stack::ParticleStack,
                           nd::NISTData, cfg::SimConfig, rng::AbstractRNG)
    pos = copy(track.position)
    dir = copy(track.direction)
    T = track.energy
    tid = objectid(track)
    gen = track.generation

    while T >= cfg.Te_cut
        # Fixed step size, clamped to not overshoot
        dEdx_col = dEdx_collision_NIST(nd, T) * cfg.rho_LXe  # MeV/cm
        ds = cfg.ds_step
        dE_col = dEdx_col * ds
        if dE_col >= T
            ds = T / dEdx_col * 0.9
            dE_col = dEdx_col * ds
        end

        # Deposit collisional energy at step midpoint
        mid_pos = pos .+ dir .* (ds * 0.5)
        push!(deposits, Deposit(copy(mid_pos), dE_col, track.kind))
        T -= dE_col

        # Advance position
        pos .= pos .+ dir .* ds
        is_inside(geom, pos) || return

        T < cfg.Te_cut && break

        # Sample bremsstrahlung
        if T > cfg.k_min
            sig_b = sigma_brems_table(nd, T)
            P_brems = min(cfg.n_atom * sig_b * ds, 0.5)  # safety cap
            if rand(rng) < P_brems
                k = sample_brems(T, cfg.k_min, cfg, rng)
                if k !== nothing && k < T
                    θ_g = brems_photon_angle(T, cfg, rng)
                    ϕ_g = 2π * rand(rng)
                    local_vec = Float64[sin(θ_g)*cos(ϕ_g), sin(θ_g)*sin(ϕ_g), cos(θ_g)]
                    d_g = rotate_to_global(local_vec, dir)
                    push!(stack, Track(:gamma, k, copy(pos), d_g,
                                       Int(tid % typemax(Int)), gen + 1))
                    T -= k
                end
            end
        end
    end

    # End-of-range: deposit residual kinetic energy
    if T > 0.0 && is_inside(geom, pos)
        push!(deposits, Deposit(copy(pos), T, track.kind))
    end

    # Positron annihilation at rest: two back-to-back 511 keV photons
    if track.kind === :positron
        cos_t = -1.0 + 2.0 * rand(rng)
        ϕ = 2π * rand(rng)
        sin_t = sqrt(1.0 - cos_t^2)
        d_g = Float64[sin_t*cos(ϕ), sin_t*sin(ϕ), cos_t]
        for d in [d_g, -d_g]
            push!(stack, Track(:gamma, cfg.me, copy(pos), d,
                               Int(tid % typemax(Int)), gen + 1))
        end
    end
end


# =====================================================================
# Top-level event simulation
# =====================================================================

"""
    simulate_event(E_MeV::Float64, nd::NISTData, cfg::SimConfig;
                   position=(0.0, 0.0, 0.0),
                   direction=(0.0, 0.0, 1.0),
                   geom::Geometry=InfiniteLXe(),
                   rng::AbstractRNG=Random.default_rng()) -> Vector{Deposit}

Simulate one primary photon entering the calorimeter with energy `E_MeV`.
Returns a list of energy `Deposit` objects.

The simulation runs a LIFO stack loop: pop the next particle, transport it
(sampling interactions, depositing energy, pushing secondaries), repeat
until the stack is empty. A generation cap prevents runaway cascades.
"""
function simulate_event(E_MeV::Float64, nd::NISTData, cfg::SimConfig;
                        position::NTuple{3,Float64}=(0.0, 0.0, 0.0),
                        direction::NTuple{3,Float64}=(0.0, 0.0, 1.0),
                        geom::Geometry=InfiniteLXe(),
                        rng::AbstractRNG=Random.default_rng())::Vector{Deposit}
    deposits = Deposit[]
    stack = ParticleStack()
    push!(stack, Track(:gamma, E_MeV,
                       Float64[position...],
                       Float64[direction...],
                       -1, 0))

    while !isempty(stack)
        t = pop!(stack)
        t.generation > cfg.generation_cap && continue
        if t.kind === :gamma
            transport_photon!(t, geom, deposits, stack, nd, cfg, rng)
        else
            transport_lepton!(t, geom, deposits, stack, nd, cfg, rng)
        end
    end

    deposits
end


# =====================================================================
# Photon-only mode (Mode 2)
# =====================================================================

"""
    transport_photon_only!(track::Track, geom::Geometry,
                           deposits::Vector{Deposit}, stack::ParticleStack,
                           nd::NISTData, cfg::SimConfig, rng::AbstractRNG)

Photon-only transport: electrons are never tracked. At each interaction:

- **Compton**: recoil electron energy deposited locally at the interaction
  point; scattered photon continues with reduced energy and new direction.
- **Pair**: entire photon energy deposited locally (both leptons + rest mass
  are sub-mm range in LXe). Photon consumed.
- **Photoelectric**: same as full mode — deposit locally.

This mode is much faster but ignores the spatial extent of electron tracks
(~mm in LXe). Useful for studying the contribution of electron transport
to SS/MS topology.
"""
function transport_photon_only!(track::Track, geom::Geometry,
                                deposits::Vector{Deposit}, stack::ParticleStack,
                                nd::NISTData, cfg::SimConfig, rng::AbstractRNG)
    pos = copy(track.position)
    dir = copy(track.direction)
    E = track.energy
    gen = track.generation
    tid = objectid(track)

    while E >= cfg.Egamma_cut
        sC  = sigma_compton_NIST(nd, E; per_atom=true)
        sP  = sigma_pair_NIST(nd, E; per_atom=true)
        sPh = sigma_phot_NIST(nd, E; per_atom=true)
        s_tot = sC + sP + sPh
        Σ_tot = cfg.n_atom * s_tot

        s = sample_distance(Σ_tot, rng)
        pos .= pos .+ dir .* s
        is_inside(geom, pos) || return

        proc = sample_process(sC/s_tot, sP/s_tot, sPh/s_tot, rng)

        if proc === :compton
            Egp, cos_t = sample_compton(E, cfg, rng)
            T_e = E - Egp

            # Deposit electron energy locally
            push!(deposits, Deposit(copy(pos), T_e, :electron))

            # Update scattered photon direction
            ϕ = 2π * rand(rng)
            sin_t = sqrt(max(0.0, 1.0 - cos_t^2))
            local_vec = Float64[sin_t*cos(ϕ), sin_t*sin(ϕ), cos_t]
            dir = rotate_to_global(local_vec, dir)
            E = Egp

        elseif proc === :pair
            # Deposit entire photon energy locally
            push!(deposits, Deposit(copy(pos), E, :pair))
            return

        elseif proc === :photoelectric
            # Deposit entire photon energy locally
            push!(deposits, Deposit(copy(pos), E, :photoelectric))
            return
        end
    end

    # Below cutoff
    push!(deposits, Deposit(pos, E, :gamma_local))
end


"""
    simulate_event_photon_only(E_MeV::Float64, nd::NISTData, cfg::SimConfig;
                               position=(0.0, 0.0, 0.0),
                               direction=(0.0, 0.0, 1.0),
                               geom::Geometry=InfiniteLXe(),
                               rng::AbstractRNG=Random.default_rng()) -> Vector{Deposit}

Photon-only simulation (Mode 2): only photons are tracked on the stack.
All electron/positron energy is deposited locally at the interaction point.
No lepton transport, no bremsstrahlung, no positron annihilation photons.

Useful for comparing SS/MS topology with and without electron transport
to quantify the effect of finite electron range in LXe.
"""
function simulate_event_photon_only(E_MeV::Float64, nd::NISTData, cfg::SimConfig;
                                    position::NTuple{3,Float64}=(0.0, 0.0, 0.0),
                                    direction::NTuple{3,Float64}=(0.0, 0.0, 1.0),
                                    geom::Geometry=InfiniteLXe(),
                                    rng::AbstractRNG=Random.default_rng())::Vector{Deposit}
    deposits = Deposit[]
    stack = ParticleStack()
    push!(stack, Track(:gamma, E_MeV,
                       Float64[position...],
                       Float64[direction...],
                       -1, 0))

    while !isempty(stack)
        t = pop!(stack)
        t.generation > cfg.generation_cap && continue
        # Only photons on the stack in this mode
        transport_photon_only!(t, geom, deposits, stack, nd, cfg, rng)
    end

    deposits
end


# =====================================================================
# Deposit clustering and SS/MS classification
# =====================================================================

"""
    cluster_deposits_in_z(deposits::Vector{Deposit}, dz_cm::Float64;
                          E_min::Float64=0.0)
        -> Vector{Tuple{Float64, Float64}}

Group deposits into clusters whose z-extent is within `dz_cm`.
Clusters with total energy below `E_min` are discarded (they represent
sub-threshold deposits invisible to the detector).

Returns a list of `(z_centroid, total_energy)` tuples, sorted by z.
The centroid is energy-weighted.
"""
function cluster_deposits_in_z(deposits::Vector{Deposit},
                               dz_cm::Float64;
                               E_min::Float64=0.0)::Vector{Tuple{Float64,Float64}}
    isempty(deposits) && return Tuple{Float64,Float64}[]

    sorted = sort(deposits, by=d -> d.position[3])

    all_clusters = Tuple{Float64,Float64}[]
    cur_zs = Float64[sorted[1].position[3]]
    cur_es = Float64[sorted[1].energy]

    @inbounds for i in 2:length(sorted)
        z = sorted[i].position[3]
        e = sorted[i].energy
        if z - cur_zs[end] < dz_cm
            push!(cur_zs, z)
            push!(cur_es, e)
        else
            E_tot = sum(cur_es)
            z_cent = sum(z * e for (z, e) in zip(cur_zs, cur_es)) / E_tot
            push!(all_clusters, (z_cent, E_tot))
            cur_zs = Float64[z]
            cur_es = Float64[e]
        end
    end

    E_tot = sum(cur_es)
    z_cent = sum(z * e for (z, e) in zip(cur_zs, cur_es)) / E_tot
    push!(all_clusters, (z_cent, E_tot))

    # Filter by minimum cluster energy
    E_min > 0.0 ? filter(c -> c[2] >= E_min, all_clusters) : all_clusters
end


"""
    is_single_site(deposits::Vector{Deposit}, dz_cm::Float64;
                   E_min::Float64=0.0) -> Bool

Return `true` if all deposits cluster into a single site with
z-resolution `dz_cm`, after discarding clusters below `E_min`.
"""
function is_single_site(deposits::Vector{Deposit}, dz_cm::Float64;
                        E_min::Float64=0.0)::Bool
    length(cluster_deposits_in_z(deposits, dz_cm; E_min=E_min)) <= 1
end
