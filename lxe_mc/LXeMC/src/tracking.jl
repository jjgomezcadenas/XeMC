"""
Stack-based Monte Carlo transport for photons and leptons.

Transport occurs inside a `PhysicalVolume` (typically a `PCyl` filled
with an active material). Material properties (cross sections, stopping
powers) come from the `Material` attached to the volume. Boundary
checking uses the volume's `is_inside` method.

# Public API

- [`simulate_event`](@ref): full simulation (photon + electron transport).
- [`simulate_event_photon_only`](@ref): photon-only mode (electrons deposit locally).
- [`is_single_site`](@ref): convenience wrapper.
"""

using Random

const TRANSPORT_BOUNDARY_TOL_CM = 1e-5
const TRANSPORT_BOUNDARY_PUSH_CM = 1e-4


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
- `track_id`: unique integer id (assigned by event loop counter)
- `parent_id`: track_id of parent (0 for primaries)
- `generation`: cascade depth (0 for primaries)
- `interaction`: what created this track (`:primary`, `:compton`, `:photoelectric`,
  `:pair`, `:bremsstrahlung`)
- `volume`: where this track was created (`:fv`, `:active`, `:passive`)
"""
struct Track
    kind::Symbol
    energy::Float64
    position::Vector{Float64}
    direction::Vector{Float64}
    track_id::Int
    parent_id::Int
    generation::Int
    interaction::Symbol
    volume::Symbol
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
- `interaction`: what interaction produced this deposit (`:compton`,
  `:photoelectric`, `:pair`, `:bremsstrahlung`, `:collisional`,
  `:escaped_gamma`, `:escaped_lepton`)
- `volume`: where the deposit occurred (`:fv`, `:active`, `:passive`)
"""
struct Deposit
    position::Vector{Float64}
    energy::Float64
    source::Symbol
    interaction::Symbol
    volume::Symbol
end


# =====================================================================
# Escape-point volume classification
# =====================================================================

"""
    _classify_escape_volume(fk, pos) -> Symbol

Classify the volume at position `pos` using the fast-kernel geometry.
Returns `:active` (Skin, TopActive, BarrelActive, BottomActive),
`:passive` (LXe_passive, FC_PTFE, FC_rings, or outside detector),
or `:fv` (should not happen for escapes, but included for completeness).
"""
@inline function _classify_escape_volume(fk::FastKernelGeometry, pos::Vector{Float64})::Symbol
    region = classify_fastkernel(fk, (pos[1], pos[2], pos[3]))
    region === nothing && return :passive
    rname = region.name
    rname == "FV" && return :fv
    rname == "Skin" && return :active
    rname in ("TopActive", "BarrelActive", "BottomActive") && return :active
    return :passive
end


# =====================================================================
# Photon transport (general, takes PhysicalVolume)
# =====================================================================

"""
    transport_photon!(track, vol, fk, deposits, stack, track_counter, cfg, rng)

Transport one gamma in `vol::PhysicalVolume` until escape, absorption,
or energy falls below `Egamma_cut`. Uses `is_inside(vol, pos)` for
boundary checking. On escape, records a deposit with the escaping
energy classified by `fk`.
"""
function transport_photon!(track::Track, vol::PhysicalVolume,
                           fk::FastKernelGeometry,
                           deposits::Vector{Deposit}, stack::ParticleStack,
                           track_counter::Ref{Int},
                           cfg::SimConfig, rng::AbstractRNG)
    mat = vol.material
    pos = copy(track.position)
    dir = copy(track.direction)
    E = track.energy
    tid = track.track_id
    gen = track.generation

    while E >= cfg.Egamma_cut
        sC, sP, sPh = sigma_three(mat, E)
        s_tot = sC + sP + sPh
        Σ_tot = mat.n_atom * s_tot

        s = sample_distance(Σ_tot, rng)
        pos .= pos .+ dir .* s
        if !is_inside(vol, pos)
            evol = _classify_escape_volume(fk, pos)
            push!(deposits, Deposit(copy(pos), E, :gamma, :escaped_gamma, evol))
            return
        end

        proc = sample_process(sC / s_tot, sP / s_tot, sPh / s_tot, rng)

        if proc === :compton
            Egp, cos_t = sample_compton(E, cfg, rng)
            ϕ = 2π * rand(rng)

            n_e = compton_electron_direction(cos_t, ϕ, E, dir, cfg)
            T_e = E - Egp
            track_counter[] += 1
            push!(stack, Track(:electron, T_e, copy(pos), n_e,
                               track_counter[], tid, gen + 1,
                               :compton, :fv))

            sin_t = sqrt(max(0.0, 1.0 - cos_t^2))
            local_vec = Float64[sin_t * cos(ϕ), sin_t * sin(ϕ), cos_t]
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
                local_vec = Float64[sin(θ) * cos(ϕ_lep), sin(θ) * sin(ϕ_lep), cos(θ)]
                d_lep = rotate_to_global(local_vec, dir)
                if T > 0.0
                    track_counter[] += 1
                    push!(stack, Track(kind, T, copy(pos), d_lep,
                                       track_counter[], tid, gen + 1,
                                       :pair, :fv))
                end
            end
            return

        elseif proc === :photoelectric
            if E < mat.EK
                push!(deposits, Deposit(copy(pos), E, :photoelectric, :photoelectric, :fv))
                return
            end

            push!(deposits, Deposit(copy(pos), mat.EK, :photoelectric, :photoelectric, :fv))

            T_e = E - mat.EK
            if T_e > cfg.Te_cut
                θ_e = sample_photoelectron_angle(T_e, cfg, rng)
                ϕ_e = 2π * rand(rng)
                local_vec = Float64[sin(θ_e) * cos(ϕ_e), sin(θ_e) * sin(ϕ_e), cos(θ_e)]
                d_e = rotate_to_global(local_vec, dir)
                track_counter[] += 1
                push!(stack, Track(:electron, T_e, copy(pos), d_e,
                                   track_counter[], tid, gen + 1,
                                   :photoelectric, :fv))
            else
                push!(deposits, Deposit(copy(pos), T_e, :photoelectric, :photoelectric, :fv))
            end
            return
        end
    end

    push!(deposits, Deposit(pos, E, :gamma_local, :gamma_local, :fv))
end


# =====================================================================
# Lepton transport (general, takes PhysicalVolume)
# =====================================================================

"""
    transport_lepton!(track, vol, fk, deposits, stack, track_counter, cfg, rng)

Transport an electron or positron through `vol::PhysicalVolume` using
condensed-history stepping. Uses `is_inside(vol, pos)` for boundary
checking. On escape, records a deposit with the remaining energy
classified by `fk`.
"""
function transport_lepton!(track::Track, vol::PhysicalVolume,
                           fk::FastKernelGeometry,
                           deposits::Vector{Deposit}, stack::ParticleStack,
                           track_counter::Ref{Int},
                           cfg::SimConfig, rng::AbstractRNG)
    mat = vol.material
    pos = copy(track.position)
    dir = copy(track.direction)
    T = track.energy
    tid = track.track_id
    gen = track.generation

    while T >= cfg.Te_cut
        dEdx_col = dEdx_collision(mat, T) * mat.density
        ds = cfg.ds_step
        dE_col = dEdx_col * ds
        if dE_col >= T
            ds = T / dEdx_col * 0.9
            dE_col = dEdx_col * ds
        end

        mid_pos = pos .+ dir .* (ds * 0.5)
        push!(deposits, Deposit(copy(mid_pos), dE_col, track.kind, :collisional, :fv))
        T -= dE_col

        pos .= pos .+ dir .* ds
        if !is_inside(vol, pos)
            evol = _classify_escape_volume(fk, pos)
            push!(deposits, Deposit(copy(pos), T, track.kind, :escaped_lepton, evol))
            return
        end

        T < cfg.Te_cut && break

        if T > cfg.k_min
            sig_b = sigma_brems(mat, T)
            P_brems = min(mat.n_atom * sig_b * ds, 0.5)
            if rand(rng) < P_brems
                M_rej = brems_rejection_M(mat, T)
                k = sample_brems(T, cfg.k_min, mat.Z_eff, M_rej, cfg, rng)
                if k !== nothing && k < T
                    θ_g = brems_photon_angle(T, cfg, rng)
                    ϕ_g = 2π * rand(rng)
                    local_vec = Float64[sin(θ_g) * cos(ϕ_g), sin(θ_g) * sin(ϕ_g), cos(θ_g)]
                    d_g = rotate_to_global(local_vec, dir)
                    track_counter[] += 1
                    push!(stack, Track(:gamma, k, copy(pos), d_g,
                                       track_counter[], tid, gen + 1,
                                       :bremsstrahlung, :fv))
                    T -= k
                end
            end
        end
    end

    if T > 0.0 && is_inside(vol, pos)
        push!(deposits, Deposit(copy(pos), T, track.kind, :collisional, :fv))
    elseif T > 0.0
        evol = _classify_escape_volume(fk, pos)
        push!(deposits, Deposit(copy(pos), T, track.kind, :escaped_lepton, evol))
    end

    if track.kind === :positron && is_inside(vol, pos)
        cos_t = -1.0 + 2.0 * rand(rng)
        ϕ = 2π * rand(rng)
        sin_t = sqrt(1.0 - cos_t^2)
        d_g = Float64[sin_t * cos(ϕ), sin_t * sin(ϕ), cos_t]
        for d in [d_g, -d_g]
            track_counter[] += 1
            push!(stack, Track(:gamma, cfg.me, copy(pos), d,
                               track_counter[], tid, gen + 1,
                               :pair, :fv))
        end
    end
end


"""
    propagate_gamma(E_MeV, vol::PhysicalVolume, fk::FastKernelGeometry, cfg;
                    position, direction, rng) -> Vector{Deposit}

Propagate one gamma inside `vol` with full stack physics (photon +
lepton cascade). Returns all deposits, including escape deposits
classified by volume (`:fv`, `:active`, `:passive`).
"""
function propagate_gamma(E_MeV::Float64, vol::PhysicalVolume,
                         fk::FastKernelGeometry, cfg::SimConfig;
                         position::NTuple{3,Float64}=(0.0, 0.0, 0.0),
                         direction::NTuple{3,Float64}=(0.0, 0.0, 1.0),
                         rng::AbstractRNG=Random.default_rng())::Vector{Deposit}
    deposits = Deposit[]
    stack = ParticleStack()
    track_counter = Ref(0)
    track_counter[] += 1
    push!(stack, Track(:gamma, E_MeV,
                       Float64[position...],
                       Float64[direction...],
                       track_counter[], 0, 0,
                       :primary, :fv))

    while !isempty(stack)
        t = pop!(stack)
        t.generation > cfg.generation_cap && continue
        if t.kind === :gamma
            transport_photon!(t, vol, fk, deposits, stack, track_counter, cfg, rng)
        else
            transport_lepton!(t, vol, fk, deposits, stack, track_counter, cfg, rng)
        end
    end

    deposits
end







