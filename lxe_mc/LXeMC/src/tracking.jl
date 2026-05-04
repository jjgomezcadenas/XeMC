"""
Stack-based Monte Carlo transport for photons and leptons.

Transport occurs inside a `PhysicalVolume` (typically a `PCyl` filled
with an active material). Material properties (cross sections, stopping
powers) come from the `Material` attached to the volume. Boundary
checking uses the volume's `is_inside` method.

# Public API

- [`simulate_event`](@ref): full simulation (photon + electron transport).
- [`simulate_event_photon_only`](@ref): photon-only mode (electrons deposit locally).
- [`cluster_deposits_in_z`](@ref): group deposits for SS/MS classification.
- [`is_single_site`](@ref): convenience wrapper.
"""

using Random


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

An energy deposit in the detector volume.

Fields:
- `position`: 3-vector [cm]
- `energy`: deposited energy [MeV]
- `source`: origin label (`:electron`, `:positron`, `:photoelectric`,
  `:pair`, `:gamma_local`)
"""
struct Deposit
    position::Vector{Float64}
    energy::Float64
    source::Symbol
end


# =====================================================================
# Photon transport (full mode)
# =====================================================================

"""
    transport_photon!(track::Track, vol::PhysicalVolume,
                      deposits::Vector{Deposit}, stack::ParticleStack,
                      cfg::SimConfig, rng::AbstractRNG)

Transport one gamma in the active volume until escape, absorption, or
energy falls below cutoff. Material properties come from `vol.material`.
"""
function transport_photon!(track::Track, vol::PhysicalVolume,
                           deposits::Vector{Deposit}, stack::ParticleStack,
                           cfg::SimConfig, rng::AbstractRNG)
    mat = vol.material
    pos = copy(track.position)
    dir = copy(track.direction)
    E = track.energy
    tid = objectid(track)
    gen = track.generation

    while E >= cfg.Egamma_cut
        sC  = sigma_compton(mat, E; per_atom=true)
        sP  = sigma_pair(mat, E; per_atom=true)
        sPh = sigma_phot(mat, E; per_atom=true)
        s_tot = sC + sP + sPh
        Σ_tot = mat.n_atom * s_tot

        s = sample_distance(Σ_tot, rng)
        pos .= pos .+ dir .* s
        is_inside(vol, pos) || return

        proc = sample_process(sC/s_tot, sP/s_tot, sPh/s_tot, rng)

        if proc === :compton
            Egp, cos_t = sample_compton(E, cfg, rng)
            ϕ = 2π * rand(rng)

            n_e = compton_electron_direction(cos_t, ϕ, E, dir, cfg)
            T_e = E - Egp
            push!(stack, Track(:electron, T_e, copy(pos), n_e,
                               Int(tid % typemax(Int)), gen + 1))

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
            return

        elseif proc === :photoelectric
            if E < mat.EK
                push!(deposits, Deposit(copy(pos), E, :photoelectric))
                return
            end

            push!(deposits, Deposit(copy(pos), mat.EK, :photoelectric))

            T_e = E - mat.EK
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
            return
        end
    end

    push!(deposits, Deposit(pos, E, :gamma_local))
end


# =====================================================================
# Lepton transport
# =====================================================================

"""
    transport_lepton!(track::Track, vol::PhysicalVolume,
                      deposits::Vector{Deposit}, stack::ParticleStack,
                      cfg::SimConfig, rng::AbstractRNG)

Transport an electron or positron using condensed-history stepping.
Material properties (stopping power, brems cross section) come from
`vol.material`.

**Stepping**: fixed step size `ds_step`, clamped to not overshoot.
**Bremsstrahlung**: probability from pre-tabulated σ_brems table.
**End of range**: deposit locally; positrons annihilate at rest.
"""
function transport_lepton!(track::Track, vol::PhysicalVolume,
                           deposits::Vector{Deposit}, stack::ParticleStack,
                           cfg::SimConfig, rng::AbstractRNG)
    mat = vol.material
    pos = copy(track.position)
    dir = copy(track.direction)
    T = track.energy
    tid = objectid(track)
    gen = track.generation

    while T >= cfg.Te_cut
        dEdx_col = dEdx_collision(mat, T) * mat.density  # MeV/cm
        ds = cfg.ds_step
        dE_col = dEdx_col * ds
        if dE_col >= T
            ds = T / dEdx_col * 0.9
            dE_col = dEdx_col * ds
        end

        mid_pos = pos .+ dir .* (ds * 0.5)
        push!(deposits, Deposit(copy(mid_pos), dE_col, track.kind))
        T -= dE_col

        pos .= pos .+ dir .* ds
        is_inside(vol, pos) || return

        T < cfg.Te_cut && break

        if T > cfg.k_min
            sig_b = sigma_brems(mat, T)
            P_brems = min(mat.n_atom * sig_b * ds, 0.5)
            if rand(rng) < P_brems
                k = sample_brems(T, cfg.k_min, mat.Z_eff, cfg, rng)
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

    if T > 0.0 && is_inside(vol, pos)
        push!(deposits, Deposit(copy(pos), T, track.kind))
    end

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
# Top-level event simulation (full mode)
# =====================================================================

"""
    simulate_event(E_MeV, vol::PhysicalVolume, cfg::SimConfig; ...) -> Vector{Deposit}

Simulate one primary photon in the active volume `vol`.
Full mode: photons and leptons are both tracked.
"""
function simulate_event(E_MeV::Float64, vol::PhysicalVolume, cfg::SimConfig;
                        position::NTuple{3,Float64}=(0.0, 0.0, 0.0),
                        direction::NTuple{3,Float64}=(0.0, 0.0, 1.0),
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
            transport_photon!(t, vol, deposits, stack, cfg, rng)
        else
            transport_lepton!(t, vol, deposits, stack, cfg, rng)
        end
    end

    deposits
end


# =====================================================================
# Photon-only mode (Mode 2)
# =====================================================================

"""
    transport_photon_only!(track::Track, vol::PhysicalVolume,
                           deposits::Vector{Deposit}, stack::ParticleStack,
                           cfg::SimConfig, rng::AbstractRNG)

Photon-only transport: electrons deposit locally. No lepton tracking.
"""
function transport_photon_only!(track::Track, vol::PhysicalVolume,
                                deposits::Vector{Deposit}, stack::ParticleStack,
                                cfg::SimConfig, rng::AbstractRNG)
    mat = vol.material
    pos = copy(track.position)
    dir = copy(track.direction)
    E = track.energy
    gen = track.generation
    tid = objectid(track)

    while E >= cfg.Egamma_cut
        sC  = sigma_compton(mat, E; per_atom=true)
        sP  = sigma_pair(mat, E; per_atom=true)
        sPh = sigma_phot(mat, E; per_atom=true)
        s_tot = sC + sP + sPh
        Σ_tot = mat.n_atom * s_tot

        s = sample_distance(Σ_tot, rng)
        pos .= pos .+ dir .* s
        is_inside(vol, pos) || return

        proc = sample_process(sC/s_tot, sP/s_tot, sPh/s_tot, rng)

        if proc === :compton
            Egp, cos_t = sample_compton(E, cfg, rng)
            T_e = E - Egp
            push!(deposits, Deposit(copy(pos), T_e, :electron))

            ϕ = 2π * rand(rng)
            sin_t = sqrt(max(0.0, 1.0 - cos_t^2))
            local_vec = Float64[sin_t*cos(ϕ), sin_t*sin(ϕ), cos_t]
            dir = rotate_to_global(local_vec, dir)
            E = Egp

        elseif proc === :pair
            push!(deposits, Deposit(copy(pos), E, :pair))
            return

        elseif proc === :photoelectric
            push!(deposits, Deposit(copy(pos), E, :photoelectric))
            return
        end
    end

    push!(deposits, Deposit(pos, E, :gamma_local))
end


"""
    simulate_event_photon_only(E_MeV, vol::PhysicalVolume, cfg::SimConfig; ...) -> Vector{Deposit}

Photon-only simulation: electrons deposit locally, no lepton transport.
"""
function simulate_event_photon_only(E_MeV::Float64, vol::PhysicalVolume, cfg::SimConfig;
                                    position::NTuple{3,Float64}=(0.0, 0.0, 0.0),
                                    direction::NTuple{3,Float64}=(0.0, 0.0, 1.0),
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
        transport_photon_only!(t, vol, deposits, stack, cfg, rng)
    end

    deposits
end


# =====================================================================
# Deposit clustering and SS/MS classification
# =====================================================================

"""
    cluster_deposits_in_z(deposits, dz_cm; E_min=0.0)

Group deposits into clusters whose z-extent is within `dz_cm`.
Clusters below `E_min` are discarded.
Returns `(z_centroid, total_energy)` tuples.
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

    E_min > 0.0 ? filter(c -> c[2] >= E_min, all_clusters) : all_clusters
end


"""
    is_single_site(deposits, dz_cm; E_min=0.0) -> Bool

True if deposits form at most one cluster above `E_min`.
"""
function is_single_site(deposits::Vector{Deposit}, dz_cm::Float64;
                        E_min::Float64=0.0)::Bool
    length(cluster_deposits_in_z(deposits, dz_cm; E_min=E_min)) <= 1
end
