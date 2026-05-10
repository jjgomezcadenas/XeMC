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
"""
struct Track
    kind::Symbol
    energy::Float64
    position::Vector{Float64}
    direction::Vector{Float64}
    track_id::Int
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
# Photon transport (general, takes PhysicalVolume)
# =====================================================================

"""
    transport_photon!(track, vol, deposits, stack, track_counter, cfg, rng)

Transport one gamma in `vol::PhysicalVolume` until escape, absorption,
or energy falls below `Egamma_cut`. Uses `is_inside(vol, pos)` for
boundary checking.
"""
function transport_photon!(track::Track, vol::PhysicalVolume,
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
        is_inside(vol, pos) || return

        proc = sample_process(sC / s_tot, sP / s_tot, sPh / s_tot, rng)

        if proc === :compton
            Egp, cos_t = sample_compton(E, cfg, rng)
            ϕ = 2π * rand(rng)

            n_e = compton_electron_direction(cos_t, ϕ, E, dir, cfg)
            T_e = E - Egp
            track_counter[] += 1
            push!(stack, Track(:electron, T_e, copy(pos), n_e,
                               track_counter[], tid, gen + 1))

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
                                       track_counter[], tid, gen + 1))
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
                local_vec = Float64[sin(θ_e) * cos(ϕ_e), sin(θ_e) * sin(ϕ_e), cos(θ_e)]
                d_e = rotate_to_global(local_vec, dir)
                track_counter[] += 1
                push!(stack, Track(:electron, T_e, copy(pos), d_e,
                                   track_counter[], tid, gen + 1))
            else
                push!(deposits, Deposit(copy(pos), T_e, :photoelectric))
            end
            return
        end
    end

    push!(deposits, Deposit(pos, E, :gamma_local))
end


# =====================================================================
# Lepton transport (general, takes PhysicalVolume)
# =====================================================================

"""
    transport_lepton!(track, vol, deposits, stack, track_counter, cfg, rng)

Transport an electron or positron through `vol::PhysicalVolume` using
condensed-history stepping. Uses `is_inside(vol, pos)` for boundary
checking.
"""
function transport_lepton!(track::Track, vol::PhysicalVolume,
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
        push!(deposits, Deposit(copy(mid_pos), dE_col, track.kind))
        T -= dE_col

        pos .= pos .+ dir .* ds
        is_inside(vol, pos) || return

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
                                       track_counter[], tid, gen + 1))
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
        d_g = Float64[sin_t * cos(ϕ), sin_t * sin(ϕ), cos_t]
        for d in [d_g, -d_g]
            track_counter[] += 1
            push!(stack, Track(:gamma, cfg.me, copy(pos), d,
                               track_counter[], tid, gen + 1))
        end
    end
end


# =====================================================================
# Photon transport in FV (legacy, uses FVGeometry)
# =====================================================================

function transport_photon_in_fv!(track::Track, fv::FVGeometry,
                                 deposits::Vector{Deposit}, stack::ParticleStack,
                                 track_counter::Ref{Int},
                                 cfg::SimConfig, rng::AbstractRNG)
    mat = fv.material
    pos = copy(track.position)
    dir = copy(track.direction)
    E = track.energy
    tid = track.track_id
    gen = track.generation

    # Egamma_cut (10 keV).
    # Photons below 10 keV are absorbed (tracking stops).
    while E >= cfg.Egamma_cut
        sC, sP, sPh = sigma_three(mat, E)
        s_tot = sC + sP + sPh
        Σ_tot = mat.n_atom * s_tot

        s = sample_distance(Σ_tot, rng)
        pos .= pos .+ dir .* s
        is_inside_fv(fv, pos) || return

        proc = sample_process(sC / s_tot, sP / s_tot, sPh / s_tot, rng)

        if proc === :compton
            Egp, cos_t = sample_compton(E, cfg, rng)
            ϕ = 2π * rand(rng)

            n_e = compton_electron_direction(cos_t, ϕ, E, dir, cfg)
            T_e = E - Egp
            track_counter[] += 1
            push!(stack, Track(:electron, T_e, copy(pos), n_e,
                               track_counter[], tid, gen + 1))

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
                                       track_counter[], tid, gen + 1))
                end
            end
            return

        elseif proc === :photoelectric
            # EK = 0.034561 MeV (34.6 keV) for LXe.
            # the xenon K-shell binding energy
            if E < mat.EK
                push!(deposits, Deposit(copy(pos), E, :photoelectric))
                return
            end

            push!(deposits, Deposit(copy(pos), mat.EK, :photoelectric))

            T_e = E - mat.EK
            if T_e > cfg.Te_cut
                θ_e = sample_photoelectron_angle(T_e, cfg, rng)
                ϕ_e = 2π * rand(rng)
                local_vec = Float64[sin(θ_e) * cos(ϕ_e), sin(θ_e) * sin(ϕ_e), cos(θ_e)]
                d_e = rotate_to_global(local_vec, dir)
                track_counter[] += 1
                push!(stack, Track(:electron, T_e, copy(pos), d_e,
                                   track_counter[], tid, gen + 1))
            else
                push!(deposits, Deposit(copy(pos), T_e, :photoelectric))
            end
            return
        end
    end

    push!(deposits, Deposit(pos, E, :gamma_local))
end


# =====================================================================
# Lepton transport in FV
# =====================================================================

function transport_lepton_in_fv!(track::Track, fv::FVGeometry,
                                 deposits::Vector{Deposit}, stack::ParticleStack,
                                 track_counter::Ref{Int},
                                 cfg::SimConfig, rng::AbstractRNG)
    mat = fv.material
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
        push!(deposits, Deposit(copy(mid_pos), dE_col, track.kind))
        T -= dE_col

        pos .= pos .+ dir .* ds
        is_inside_fv(fv, pos) || return

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
                                       track_counter[], tid, gen + 1))
                    T -= k
                end
            end
        end
    end

    if T > 0.0 && is_inside_fv(fv, pos)
        push!(deposits, Deposit(copy(pos), T, track.kind))
    end

    if track.kind === :positron
        cos_t = -1.0 + 2.0 * rand(rng)
        ϕ = 2π * rand(rng)
        sin_t = sqrt(1.0 - cos_t^2)
        d_g = Float64[sin_t * cos(ϕ), sin_t * sin(ϕ), cos_t]
        for d in [d_g, -d_g]
            track_counter[] += 1
            push!(stack, Track(:gamma, cfg.me, copy(pos), d,
                               track_counter[], tid, gen + 1))
        end
    end
end


"""
    propagate_gamma(E_MeV, vol::PhysicalVolume, cfg; position, direction, rng) -> Vector{Deposit}

Propagate one gamma inside `vol` with full stack physics (photon +
lepton cascade). Returns all deposits inside the volume.
"""
function propagate_gamma(E_MeV::Float64, vol::PhysicalVolume, cfg::SimConfig;
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
                       track_counter[], 0, 0))

    while !isempty(stack)
        t = pop!(stack)
        t.generation > cfg.generation_cap && continue
        if t.kind === :gamma
            transport_photon!(t, vol, deposits, stack, track_counter, cfg, rng)
        else
            transport_lepton!(t, vol, deposits, stack, track_counter, cfg, rng)
        end
    end

    deposits
end


"""
    propagate_gamma_in_fv(E_MeV, fv::FVGeometry, cfg; position, direction, rng) -> Vector{Deposit}

Legacy wrapper. Propagate one gamma inside the FV cylinder described
by `FVGeometry`. Delegates to `propagate_gamma` via a PCyl.
"""
function propagate_gamma_in_fv(E_MeV::Float64, fv::FVGeometry, cfg::SimConfig;
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
                       track_counter[], 0, 0))

    while !isempty(stack)
        t = pop!(stack)
        t.generation > cfg.generation_cap && continue
        if t.kind === :gamma
            transport_photon_in_fv!(t, fv, deposits, stack, track_counter, cfg, rng)
        else
            transport_lepton_in_fv!(t, fv, deposits, stack, track_counter, cfg, rng)
        end
    end

    deposits
end


struct GammaPropagationResult
    status::Symbol
    interaction_type::Symbol
    deposit_E_MeV::Float64
    position::Vector{Float64}
    region::String
end


struct EventProcessingResult
    status::Symbol
    has_fv::Bool
    has_tpc_veto::Bool
    has_skin_veto::Bool
    n_processed::Int
    deposits::Vector{Deposit}  # FV deposits for :accepted events; empty otherwise
end


struct FastKernelCalibEventResult
    result::EventProcessingResult
    multiplicity::Int
    E1_MeV::Float64
    E2_MeV::Float64
end


"""
    FastGammaDeposit

One local gamma energy deposit used by the continued fast-kernel event
path. `coarse_region` is the transport region where the interaction
occurred, while `region` stores the analysis-refined region:

- `FV`
- `LXeTPC`
- `Skin`
- passive transport region name (`LXe_det`, `FC_PTFE`, `FC_rings`, ...)
"""
struct FastGammaDeposit
    Edep_MeV::Float64
    position::Vector{Float64}
    coarse_region::String
    region::String
end


"""
    FastGammaTrackResult

Result of continued fast-kernel photon transport for one gamma. The
transport follows repeated Compton scatters until a decisive analysis
outcome is reached or the gamma is terminated.

`status` values used by the current fast event path:

- `:escaped`
- `:below_cut`
- `:absorbed_passive`
- `:below_roi_fv`
- `:vetoed_tpc`
- `:vetoed_skin`
- `:ms_rejected`
"""
struct FastGammaTrackResult
    status::Symbol
    deposits::Vector{FastGammaDeposit}
    terminal_region::String
    energy_MeV::Float64
    position::Vector{Float64}
    direction::Vector{Float64}
end


function propagate_gamma_fastkernel(gamma,
                                    fk::FastKernelGeometry,
                                    cfg::SimConfig,
                                    rng::AbstractRNG;
                                    max_cm::Float64=400.0)::GammaPropagationResult

    pos = copy(gamma.position)
    dir = copy(gamma.direction)
    E = gamma.E_MeV
    traveled = 0.0

    while E >= cfg.Egamma_cut && traveled < max_cm
        region = classify_fastkernel(fk, (pos[1], pos[2], pos[3]))
        region === nothing && return GammaPropagationResult(:escaped, :none, 0.0, copy(pos), "MARS")
        region.name == "FV" && return GammaPropagationResult(:entered_fv, :none, 0.0, copy(pos), "FV")

        mat = region.material
        s_bnd = distance_to_boundary_fastkernel(region, (pos[1], pos[2], pos[3]), (dir[1], dir[2], dir[3]))
        if !isfinite(s_bnd)
            return GammaPropagationResult(:escaped, :none, 0.0, copy(pos), region.name)
        end

        if mat.density <= 0.0 || mat.xcom === nothing
            pos .= pos .+ dir .* (s_bnd + TRANSPORT_BOUNDARY_PUSH_CM)
            traveled += s_bnd + TRANSPORT_BOUNDARY_PUSH_CM
            continue
        end

        sC, sP, sPh = sigma_three(mat, E)
        s_tot = sC + sP + sPh

        if s_tot <= 0.0
            pos .= pos .+ dir .* (s_bnd + TRANSPORT_BOUNDARY_PUSH_CM)
            traveled += s_bnd + TRANSPORT_BOUNDARY_PUSH_CM
            continue
        end

        Σ_tot = mat.n_atom * s_tot
        s_int = sample_distance(Σ_tot, rng)

        if s_int + TRANSPORT_BOUNDARY_TOL_CM < s_bnd
            pos .= pos .+ dir .* s_int
            region_after = classify_fastkernel(fk, (pos[1], pos[2], pos[3]))
            region_name = region_after === nothing ? region.name : region_after.name

            proc = sample_process(sC/s_tot, sP/s_tot, sPh/s_tot, rng)
            if proc === :compton
                Egp, _ = sample_compton(E, cfg, rng)
                return GammaPropagationResult(:interacted, :compton, E - Egp, copy(pos), region_name)
            elseif proc === :pair
                return GammaPropagationResult(:interacted, :pair, E, copy(pos), region_name)
            else
                return GammaPropagationResult(:interacted, :photoelectric, E, copy(pos), region_name)
            end
        end

        pos .= pos .+ dir .* (s_bnd + TRANSPORT_BOUNDARY_PUSH_CM)
        traveled += s_bnd + TRANSPORT_BOUNDARY_PUSH_CM
    end

    GammaPropagationResult(:below_cut, :none, 0.0, copy(pos), "MARS")
end


@inline function _visible_threshold_MeV(cfg::SimConfig, region_name::String)::Float64
    region_name == "Skin" && return cfg.veto_skin
    region_name in ("TopActive", "BarrelActive", "BottomActive") && return cfg.veto_TPC
    region_name == "FV" && return cfg.veto_TPC
    return Inf
end


@inline function _is_passive_region(region_name::String)::Bool
    region_name == "LXe_passive" || region_name == "FC_PTFE" || region_name == "FC_rings"
end


function _replace_last_deposit!(deposits::Vector{FastGammaDeposit}, dep::FastGammaDeposit)
    deposits[end] = dep
    deposits
end


@inline function _fast_track_result(status::Symbol,
                                    deposits::Vector{FastGammaDeposit},
                                    terminal_region::String,
                                    E::Float64,
                                    pos::Vector{Float64},
                                    dir::Vector{Float64})
    FastGammaTrackResult(status, deposits, terminal_region, E, copy(pos), copy(dir))
end


"""
    transport_gamma_fastkernel(gamma, fk, cfg, rng; max_cm=400.0)
        -> FastGammaTrackResult

Continued fast-kernel photon transport used by the fast event path.

Unlike `propagate_gamma_fastkernel`, which stops at the first interaction
for geometry/transport validation, this routine continues after Compton
scatters and applies analysis-driven early exits:

- immediate veto in `Skin` and `LXeTPC \\ FV` when the local deposit is
  above threshold
- ROI-floor termination when the surviving photon energy falls below
  `cfg.E_roi_floor`
- local collapse of the remaining gamma energy at the current position
  when the ROI-floor condition is met

The returned deposits are later folded into event-level SS/MS logic.
"""
function transport_gamma_fastkernel(gamma,
                                    fk::FastKernelGeometry,
                                    cfg::SimConfig,
                                    rng::AbstractRNG;
                                    max_cm::Float64=400.0)::FastGammaTrackResult

    pos = copy(gamma.position)
    dir = copy(gamma.direction)
    E = gamma.E_MeV
    traveled = 0.0
    deposits = FastGammaDeposit[]

    while E >= cfg.Egamma_cut && traveled < max_cm
        coarse_region = classify_fastkernel(fk, (pos[1], pos[2], pos[3]))
        coarse_region === nothing && return _fast_track_result(:escaped, deposits, "MARS", E, pos, dir)

        if coarse_region.name == "FV"
            return _fast_track_result(:handoff_fv, deposits, "FV", E, pos, dir)
        end

        mat = coarse_region.material
        s_bnd = distance_to_boundary_fastkernel(coarse_region, (pos[1], pos[2], pos[3]), (dir[1], dir[2], dir[3]))
        if !isfinite(s_bnd)
            return _fast_track_result(:escaped, deposits, coarse_region.name, E, pos, dir)
        end

        if mat.density <= 0.0 || mat.xcom === nothing
            pos .= pos .+ dir .* (s_bnd + TRANSPORT_BOUNDARY_PUSH_CM)
            traveled += s_bnd + TRANSPORT_BOUNDARY_PUSH_CM
            continue
        end

        sC, sP, sPh = sigma_three(mat, E)
        s_tot = sC + sP + sPh

        if s_tot <= 0.0
            pos .= pos .+ dir .* (s_bnd + TRANSPORT_BOUNDARY_PUSH_CM)
            traveled += s_bnd + TRANSPORT_BOUNDARY_PUSH_CM
            continue
        end

        Σ_tot = mat.n_atom * s_tot
        s_int = sample_distance(Σ_tot, rng)

        if s_int + TRANSPORT_BOUNDARY_TOL_CM < s_bnd
            pos .= pos .+ dir .* s_int
            refined_region = coarse_region.name

            proc = sample_process(sC/s_tot, sP/s_tot, sPh/s_tot, rng)
            if proc === :compton
                Egp, cos_theta = sample_compton(E, cfg, rng)
                dep = E - Egp
                push!(deposits, FastGammaDeposit(dep, copy(pos), coarse_region.name, refined_region))

                if refined_region == "Skin" && dep >= cfg.veto_skin
                    return _fast_track_result(:vetoed_skin, deposits, refined_region, Egp, pos, dir)
                elseif refined_region in ("TopActive", "BarrelActive", "BottomActive") && dep >= cfg.veto_TPC
                    return _fast_track_result(:vetoed_tpc, deposits, refined_region, Egp, pos, dir)
                elseif refined_region == "FV" && dep >= cfg.veto_TPC
                    sin_theta = sqrt(max(0.0, 1.0 - cos_theta^2))
                    phi = 2π * rand(rng)
                    local_dir = Float64[sin_theta * cos(phi), sin_theta * sin(phi), cos_theta]
                    dir .= rotate_to_global(local_dir, dir)
                    E = Egp
                    if E < cfg.E_roi_floor
                        dep_total = deposits[end].Edep_MeV + E
                        _replace_last_deposit!(deposits,
                            FastGammaDeposit(dep_total, deposits[end].position,
                                             deposits[end].coarse_region,
                                             deposits[end].region))
                        return _fast_track_result(:below_roi_fv, deposits, refined_region, 0.0, pos, dir)
                    end
                    continue
                else
                    E = Egp
                    if E < cfg.E_roi_floor
                        dep_total = deposits[end].Edep_MeV + E
                        _replace_last_deposit!(deposits,
                            FastGammaDeposit(dep_total, deposits[end].position,
                                             deposits[end].coarse_region,
                                             deposits[end].region))

                        if refined_region == "Skin"
                            return _fast_track_result(:vetoed_skin, deposits, refined_region, 0.0, pos, dir)
                        elseif refined_region in ("TopActive", "BarrelActive", "BottomActive")
                            return _fast_track_result(:vetoed_tpc, deposits, refined_region, 0.0, pos, dir)
                        elseif refined_region == "FV"
                            return _fast_track_result(:below_roi_fv, deposits, refined_region, 0.0, pos, dir)
                        elseif _is_passive_region(refined_region)
                            return _fast_track_result(:absorbed_passive, deposits, refined_region, 0.0, pos, dir)
                        else
                            return _fast_track_result(:below_cut, deposits, refined_region, 0.0, pos, dir)
                        end
                    end

                    sin_theta = sqrt(max(0.0, 1.0 - cos_theta^2))
                    phi = 2π * rand(rng)
                    local_dir = Float64[sin_theta * cos(phi), sin_theta * sin(phi), cos_theta]
                    dir .= rotate_to_global(local_dir, dir)
                    continue
                end
            elseif proc === :pair
                push!(deposits, FastGammaDeposit(E, copy(pos), coarse_region.name, refined_region))
                if refined_region == "Skin"
                    return _fast_track_result(:vetoed_skin, deposits, refined_region, 0.0, pos, dir)
                elseif refined_region in ("TopActive", "BarrelActive", "BottomActive")
                    return _fast_track_result(:vetoed_tpc, deposits, refined_region, 0.0, pos, dir)
                elseif refined_region == "FV"
                    return _fast_track_result(:below_roi_fv, deposits, refined_region, 0.0, pos, dir)
                elseif _is_passive_region(refined_region)
                    return _fast_track_result(:absorbed_passive, deposits, refined_region, 0.0, pos, dir)
                else
                    return _fast_track_result(:below_cut, deposits, refined_region, 0.0, pos, dir)
                end
            else
                push!(deposits, FastGammaDeposit(E, copy(pos), coarse_region.name, refined_region))
                if refined_region == "Skin"
                    return _fast_track_result(:vetoed_skin, deposits, refined_region, 0.0, pos, dir)
                elseif refined_region in ("TopActive", "BarrelActive", "BottomActive")
                    return _fast_track_result(:vetoed_tpc, deposits, refined_region, 0.0, pos, dir)
                elseif refined_region == "FV"
                    return _fast_track_result(:below_roi_fv, deposits, refined_region, 0.0, pos, dir)
                elseif _is_passive_region(refined_region)
                    return _fast_track_result(:absorbed_passive, deposits, refined_region, 0.0, pos, dir)
                else
                    return _fast_track_result(:below_cut, deposits, refined_region, 0.0, pos, dir)
                end
            end
        end

        pos .= pos .+ dir .* (s_bnd + TRANSPORT_BOUNDARY_PUSH_CM)
        traveled += s_bnd + TRANSPORT_BOUNDARY_PUSH_CM
    end

    _fast_track_result(:below_cut, deposits, "MARS", E, pos, dir)
end


function _update_event_state(sel)
    has_fv = sel.class == :fv && sel.passes_threshold
    has_tpc_veto = sel.class == :tpc && sel.passes_threshold
    has_skin_veto = sel.class == :skin && sel.passes_threshold
    vetoed = has_tpc_veto || has_skin_veto
    (has_fv=has_fv, has_tpc_veto=has_tpc_veto, has_skin_veto=has_skin_veto, vetoed=vetoed)
end


function _fold_fast_deposit(dep::FastGammaDeposit,
                            has_fv::Bool,
                            has_tpc_veto::Bool,
                            has_skin_veto::Bool,
                            fv_z_ref::Float64,
                            cfg::SimConfig)
    if dep.region == "FV"
        if dep.Edep_MeV >= cfg.veto_TPC
            if !has_fv
                return (has_fv=true, has_tpc_veto=has_tpc_veto, has_skin_veto=has_skin_veto,
                        fv_z_ref=dep.position[3], vetoed=false, ms_rejected=false)
            elseif abs(dep.position[3] - fv_z_ref) > cfg.dz_resolution
                return (has_fv=true, has_tpc_veto=has_tpc_veto, has_skin_veto=has_skin_veto,
                        fv_z_ref=fv_z_ref, vetoed=true, ms_rejected=true)
            else
                return (has_fv=true, has_tpc_veto=has_tpc_veto, has_skin_veto=has_skin_veto,
                        fv_z_ref=fv_z_ref, vetoed=false, ms_rejected=false)
            end
        end
    elseif dep.region == "LXeTPC" && dep.Edep_MeV >= cfg.veto_TPC
        return (has_fv=has_fv, has_tpc_veto=true, has_skin_veto=has_skin_veto,
                fv_z_ref=fv_z_ref, vetoed=true, ms_rejected=false)
    elseif dep.region == "Skin" && dep.Edep_MeV >= cfg.veto_skin
        return (has_fv=has_fv, has_tpc_veto=has_tpc_veto, has_skin_veto=true,
                fv_z_ref=fv_z_ref, vetoed=true, ms_rejected=false)
    end

    (has_fv=has_fv, has_tpc_veto=has_tpc_veto, has_skin_veto=has_skin_veto,
     fv_z_ref=fv_z_ref, vetoed=false, ms_rejected=false)
end


function _classify_fv_stack_result(deposits::Vector{Deposit}, cfg::SimConfig)::Symbol
    clusters = cluster_deposits_in_z(deposits, cfg.dz_resolution; E_min=cfg.E_cluster_min)
    isempty(clusters) && return :no_fv
    length(clusters) == 1 && return :accepted
    :ms_rejected
end


function process_event(gammas, fk::FastKernelGeometry, fv::FVGeometry,
                       cfg::SimConfig, rng::AbstractRNG)
    has_fv = false
    has_tpc_veto = false
    has_skin_veto = false
    fv_z_ref = NaN
    n_processed = 0
    fv_deposits = Deposit[]

    empty_deps = Deposit[]

    for gamma in gammas
        result = transport_gamma_fastkernel(gamma, fk, cfg, rng)
        n_processed += 1

        for dep in result.deposits
            state = _fold_fast_deposit(dep, has_fv, has_tpc_veto, has_skin_veto, fv_z_ref, cfg)
            has_fv = state.has_fv
            has_tpc_veto = state.has_tpc_veto
            has_skin_veto = state.has_skin_veto
            fv_z_ref = state.fv_z_ref
            if state.vetoed
                return EventProcessingResult(:vetoed, has_fv, has_tpc_veto, has_skin_veto, n_processed, empty_deps)
            end
        end

        if result.status == :handoff_fv
            deps = propagate_gamma_in_fv(
                result.energy_MeV,
                fv,
                cfg;
                position=(result.position[1], result.position[2], result.position[3]),
                direction=(result.direction[1], result.direction[2], result.direction[3]),
                rng=rng,
            )
            fv_status = _classify_fv_stack_result(deps, cfg)
            if fv_status == :accepted
                has_fv = true
                append!(fv_deposits, deps)
            elseif fv_status == :ms_rejected
                return EventProcessingResult(:ms_rejected, true, has_tpc_veto, has_skin_veto, n_processed, empty_deps)
            end
        elseif result.status == :vetoed_tpc
            has_tpc_veto = true
            return EventProcessingResult(:vetoed, has_fv, has_tpc_veto, has_skin_veto, n_processed, empty_deps)
        elseif result.status == :vetoed_skin
            has_skin_veto = true
            return EventProcessingResult(:vetoed, has_fv, has_tpc_veto, has_skin_veto, n_processed, empty_deps)
        end
    end

    if has_fv
        return EventProcessingResult(:accepted, has_fv, has_tpc_veto, has_skin_veto, n_processed, fv_deposits)
    end
    EventProcessingResult(:no_fv, has_fv, has_tpc_veto, has_skin_veto, n_processed, empty_deps)
end


"""
    process_event(gammas, fk, vol::PhysicalVolume, cfg, rng) -> EventProcessingResult

Process a multi-gamma event through the fast kernel and full stack.
General version taking any PhysicalVolume for the FV region.
"""
function process_event(gammas, fk::FastKernelGeometry, vol::PhysicalVolume,
                       cfg::SimConfig, rng::AbstractRNG)
    has_fv = false
    has_tpc_veto = false
    has_skin_veto = false
    fv_z_ref = NaN
    n_processed = 0
    fv_deposits = Deposit[]

    empty_deps = Deposit[]

    for gamma in gammas
        result = transport_gamma_fastkernel(gamma, fk, cfg, rng)
        n_processed += 1

        for dep in result.deposits
            state = _fold_fast_deposit(dep, has_fv, has_tpc_veto, has_skin_veto, fv_z_ref, cfg)
            has_fv = state.has_fv
            has_tpc_veto = state.has_tpc_veto
            has_skin_veto = state.has_skin_veto
            fv_z_ref = state.fv_z_ref
            if state.vetoed
                return EventProcessingResult(:vetoed, has_fv, has_tpc_veto, has_skin_veto, n_processed, empty_deps)
            end
        end

        if result.status == :handoff_fv
            deps = propagate_gamma(
                result.energy_MeV,
                vol,
                cfg;
                position=(result.position[1], result.position[2], result.position[3]),
                direction=(result.direction[1], result.direction[2], result.direction[3]),
                rng=rng,
            )
            fv_status = _classify_fv_stack_result(deps, cfg)
            if fv_status == :accepted
                has_fv = true
                append!(fv_deposits, deps)
            elseif fv_status == :ms_rejected
                return EventProcessingResult(:ms_rejected, true, has_tpc_veto, has_skin_veto, n_processed, empty_deps)
            end
        elseif result.status == :vetoed_tpc
            has_tpc_veto = true
            return EventProcessingResult(:vetoed, has_fv, has_tpc_veto, has_skin_veto, n_processed, empty_deps)
        elseif result.status == :vetoed_skin
            has_skin_veto = true
            return EventProcessingResult(:vetoed, has_fv, has_tpc_veto, has_skin_veto, n_processed, empty_deps)
        end
    end

    if has_fv
        return EventProcessingResult(:accepted, has_fv, has_tpc_veto, has_skin_veto, n_processed, fv_deposits)
    end
    EventProcessingResult(:no_fv, has_fv, has_tpc_veto, has_skin_veto, n_processed, empty_deps)
end


function process_event_fastkernel_calib(fk::FastKernelGeometry, fv::FVGeometry, cfg::SimConfig,
                                        rng::AbstractRNG;
                                        E_MeV::Float64=2.615,
                                        x_cm::Float64=0.0,
                                        y_cm::Float64=0.0,
                                        z_cm::Float64=160.0,
                                        ux::Float64=1.0,
                                        uy::Float64=0.0,
                                        uz::Float64=0.0)
    gammas = sample_gammas("calib";
                           calib=true,
                           E_MeV=E_MeV,
                           x_cm=x_cm, y_cm=y_cm, z_cm=z_cm,
                           ux=ux, uy=uy, uz=uz,
                           rng=rng)
    multiplicity = length(gammas)

    E1 = multiplicity >= 1 ? gammas[1].E_MeV : 0.0
    E2 = multiplicity >= 2 ? gammas[2].E_MeV : 0.0
    FastKernelCalibEventResult(
        process_event(gammas, fk, fv, cfg, rng),
        multiplicity,
        E1,
        E2,
    )
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
