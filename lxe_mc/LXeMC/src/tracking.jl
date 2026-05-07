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
# Photon transport (full mode)
# =====================================================================

"""
    transport_photon!(track, vol, deposits, stack, track_counter, cfg, rng)

Transport one gamma in the active volume `vol` until escape, absorption,
or energy falls below `Egamma_cut`.

At each interaction point, cross sections are looked up via
[`sigma_three`](@ref) (single binary search for Compton + pair + photo).
The interaction channel is sampled, and secondaries are pushed onto
`stack` with IDs from `track_counter`. Material properties (cross
sections, n_atom, EK) come from `vol.material`.

Photoelectric absorption deposits E_K locally (relaxation cascade is
sub-mm) and only tracks the photoelectron if T_e > Te_cut.
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

        proc = sample_process(sC/s_tot, sP/s_tot, sPh/s_tot, rng)

        if proc === :compton
            Egp, cos_t = sample_compton(E, cfg, rng)
            ϕ = 2π * rand(rng)

            n_e = compton_electron_direction(cos_t, ϕ, E, dir, cfg)
            T_e = E - Egp
            track_counter[] += 1
            push!(stack, Track(:electron, T_e, copy(pos), n_e,
                               track_counter[], tid, gen + 1))

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
                local_vec = Float64[sin(θ_e)*cos(ϕ_e), sin(θ_e)*sin(ϕ_e), cos(θ_e)]
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
# Lepton transport
# =====================================================================

"""
    transport_lepton!(track, vol, deposits, stack, track_counter, cfg, rng)

Transport an electron or positron through `vol` using condensed-history
stepping with continuous collisional energy loss and discrete brems.

**Stepping**: fixed step `ds_step` [cm], clamped so dE < T. Collisional
energy loss dE = (dE/dx)_col × ρ × ds deposited at step midpoint.

**Bremsstrahlung**: per-step probability P = n_atom × σ_brems(T) × ds,
with σ_brems and rejection ceiling M from pre-tabulated log grids.
Hard photons (k > k_min) pushed onto stack.

**End of range**: residual T deposited locally. Positrons annihilate
at rest → two back-to-back 511 keV photons.
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
                M_rej = brems_rejection_M(mat, T)
                k = sample_brems(T, cfg.k_min, mat.Z_eff, M_rej, cfg, rng)
                if k !== nothing && k < T
                    θ_g = brems_photon_angle(T, cfg, rng)
                    ϕ_g = 2π * rand(rng)
                    local_vec = Float64[sin(θ_g)*cos(ϕ_g), sin(θ_g)*sin(ϕ_g), cos(θ_g)]
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
        d_g = Float64[sin_t*cos(ϕ), sin_t*sin(ϕ), cos_t]
        for d in [d_g, -d_g]
            track_counter[] += 1
            push!(stack, Track(:gamma, cfg.me, copy(pos), d,
                               track_counter[], tid, gen + 1))
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
    tid = track.track_id

    while E >= cfg.Egamma_cut
        sC, sP, sPh = sigma_three(mat, E)
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
                       1, 0, 0))

    while !isempty(stack)
        t = pop!(stack)
        t.generation > cfg.generation_cap && continue
        transport_photon_only!(t, vol, deposits, stack, cfg, rng)
    end

    deposits
end


# =====================================================================
# Fast veto pre-filter
# =====================================================================

"""
    PropagationResult

Result of [`propagate_to_fiducial`](@ref):
- `:vetoed` — gamma deposited visible energy outside FV (event killed by veto)
- `:lost` — gamma absorbed or fell below Egamma_cut without triggering veto
  (invisible to detector, no energy in FV)
- `:accepted` — gamma reached the FV. `energy` and `position` give its
  state at the FV entry point for full simulation.
"""
struct PropagationResult
    status::Symbol          # :vetoed, :lost, :accepted
    energy::Float64         # gamma energy at outcome point [MeV]
    position::Vector{Float64}  # position at outcome point [cm]
    direction::Vector{Float64} # direction at outcome point
    n_interactions::Int     # number of photon interactions before outcome
end


struct GammaPropagationV2Result
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
end


"""
    propagate_gamma_v2(gamma, det, cfg, rng; step_cm=0.05, max_cm=400.0)
        -> GammaPropagationV2Result

Propagate one sampled gamma through the Geometry V2 tracking tree until:

- it enters `FV` before interacting (`status = :entered_fv`)
- it has its first interaction in a tracking material (`status = :interacted`)
- it escapes the tracking geometry (`status = :escaped`)
- it falls below the photon cutoff (`status = :below_cut`)

This first V2 transport stage samples Compton, pair-production, and
photoelectric interactions, but still stops at the first interaction.
"""
function propagate_gamma_v2(gamma,
                            det::DetectorV2,
                            cfg::SimConfig,
                            rng::AbstractRNG;
                            step_cm::Float64=0.05,
                            max_cm::Float64=400.0)::GammaPropagationV2Result

    pos = copy(gamma.position)
    dir = copy(gamma.direction)
    E = gamma.E_MeV
    traveled = 0.0

    while E >= cfg.Egamma_cut && traveled < max_cm
        cls = classify_runtime_v2(det, (pos[1], pos[2], pos[3]))
        cls === nothing && return GammaPropagationV2Result(:escaped, :none, 0.0, copy(pos), "MARS")
        cls.fv_target && return GammaPropagationV2Result(:entered_fv, :none, 0.0, copy(pos), cls.name)

        mat = cls.node.lv.material
        s_bnd, _ = distance_to_node_change_v2(det, pos, dir; step_cm=step_cm, max_cm=max_cm - traveled)

        if !isfinite(s_bnd)
            return GammaPropagationV2Result(:escaped, :none, 0.0, copy(pos), cls.name)
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
            region = classify_runtime_v2(det, (pos[1], pos[2], pos[3]))
            region_name = region === nothing ? cls.name : region.name

            proc = sample_process(sC/s_tot, sP/s_tot, sPh/s_tot, rng)
            if proc === :compton
                Egp, _ = sample_compton(E, cfg, rng)
                return GammaPropagationV2Result(:interacted, :compton, E - Egp, copy(pos), region_name)
            elseif proc === :pair
                return GammaPropagationV2Result(:interacted, :pair, E, copy(pos), region_name)
            else
                return GammaPropagationV2Result(:interacted, :photoelectric, E, copy(pos), region_name)
            end
        end

        pos .= pos .+ dir .* (s_bnd + TRANSPORT_BOUNDARY_PUSH_CM)
        traveled += s_bnd + TRANSPORT_BOUNDARY_PUSH_CM
    end

    GammaPropagationV2Result(:below_cut, :none, 0.0, copy(pos), "MARS")
end


function propagate_gamma_fastkernel(gamma,
                                    fk::FastKernelGeometry,
                                    cfg::SimConfig,
                                    rng::AbstractRNG;
                                    max_cm::Float64=400.0)::GammaPropagationV2Result

    pos = copy(gamma.position)
    dir = copy(gamma.direction)
    E = gamma.E_MeV
    traveled = 0.0

    while E >= cfg.Egamma_cut && traveled < max_cm
        region = classify_fastkernel(fk, (pos[1], pos[2], pos[3]))
        region === nothing && return GammaPropagationV2Result(:escaped, :none, 0.0, copy(pos), "MARS")
        region.fv_target && return GammaPropagationV2Result(:entered_fv, :none, 0.0, copy(pos), region.name)

        mat = region.material
        s_bnd = distance_to_boundary_fastkernel(region, (pos[1], pos[2], pos[3]), (dir[1], dir[2], dir[3]))
        if !isfinite(s_bnd)
            return GammaPropagationV2Result(:escaped, :none, 0.0, copy(pos), region.name)
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
                return GammaPropagationV2Result(:interacted, :compton, E - Egp, copy(pos), region_name)
            elseif proc === :pair
                return GammaPropagationV2Result(:interacted, :pair, E, copy(pos), region_name)
            else
                return GammaPropagationV2Result(:interacted, :photoelectric, E, copy(pos), region_name)
            end
        end

        pos .= pos .+ dir .* (s_bnd + TRANSPORT_BOUNDARY_PUSH_CM)
        traveled += s_bnd + TRANSPORT_BOUNDARY_PUSH_CM
    end

    GammaPropagationV2Result(:below_cut, :none, 0.0, copy(pos), "MARS")
end


function _update_event_state(sel)
    has_fv = sel.class == :fv && sel.passes_threshold
    has_tpc_veto = sel.class == :tpc && sel.passes_threshold
    has_skin_veto = sel.class == :skin && sel.passes_threshold
    vetoed = has_tpc_veto || has_skin_veto
    (has_fv=has_fv, has_tpc_veto=has_tpc_veto, has_skin_veto=has_skin_veto, vetoed=vetoed)
end


function process_event_from_selections(selections)
    has_fv = false
    has_tpc_veto = false
    has_skin_veto = false

    for (i, sel) in enumerate(selections)
        state = _update_event_state(sel)
        has_fv |= state.has_fv
        has_tpc_veto |= state.has_tpc_veto
        has_skin_veto |= state.has_skin_veto

        if state.vetoed
            return EventProcessingResult(:vetoed, has_fv, has_tpc_veto, has_skin_veto, i)
        end
    end

    if has_fv
        return EventProcessingResult(:accepted, has_fv, has_tpc_veto, has_skin_veto, length(selections))
    end
    EventProcessingResult(:no_fv, has_fv, has_tpc_veto, has_skin_veto, length(selections))
end


function process_event(gammas, det::DetectorV2, cfg::SimConfig, rng::AbstractRNG)
    has_fv = false
    has_tpc_veto = false
    has_skin_veto = false
    n_processed = 0

    for gamma in gammas
        result = propagate_gamma_v2(gamma, det, cfg, rng)
        sel = select_interaction(result, det)
        n_processed += 1
        state = _update_event_state(sel)
        has_fv |= state.has_fv
        has_tpc_veto |= state.has_tpc_veto
        has_skin_veto |= state.has_skin_veto
        if state.vetoed
            return EventProcessingResult(:vetoed, has_fv, has_tpc_veto, has_skin_veto, n_processed)
        end
    end

    if has_fv
        return EventProcessingResult(:accepted, has_fv, has_tpc_veto, has_skin_veto, n_processed)
    end
    EventProcessingResult(:no_fv, has_fv, has_tpc_veto, has_skin_veto, n_processed)
end


function process_event_fastkernel(gammas, fk::FastKernelGeometry, det::DetectorV2,
                                  cfg::SimConfig, rng::AbstractRNG)
    has_fv = false
    has_tpc_veto = false
    has_skin_veto = false
    n_processed = 0

    for gamma in gammas
        result = propagate_gamma_fastkernel(gamma, fk, cfg, rng)
        sel = select_interaction(result, det)
        n_processed += 1
        state = _update_event_state(sel)
        has_fv |= state.has_fv
        has_tpc_veto |= state.has_tpc_veto
        has_skin_veto |= state.has_skin_veto
        if state.vetoed
            return EventProcessingResult(:vetoed, has_fv, has_tpc_veto, has_skin_veto, n_processed)
        end
    end

    if has_fv
        return EventProcessingResult(:accepted, has_fv, has_tpc_veto, has_skin_veto, n_processed)
    end
    EventProcessingResult(:no_fv, has_fv, has_tpc_veto, has_skin_veto, n_processed)
end


"""
    propagate_to_fiducial(E_MeV, position, direction, det, cfg, rng)
        -> PropagationResult

Lightweight photon-only propagation from a source surface toward the
fiducial volume. At each interaction, checks whether the energy deposit
would be visible to the veto:

- **Skin** (inside TPC volume but outside active, or in Skin volume):
  deposit > `veto_skin` → `:vetoed`
- **TPC active** (inside TPC but outside FV):
  deposit > `veto_TPC` → `:vetoed`
- **FV**: gamma arrived → `:accepted`

If the deposit is below threshold, the gamma is invisible to the
detector. For Compton, it continues with reduced energy. For
photoelectric/pair, it's absorbed → `:lost`.

No electron tracking, no secondaries — only photon steps.
Typical rejection rate: 95–98% of gammas from nearby surfaces.
"""
function propagate_to_fiducial(E_MeV::Float64,
                               position::NTuple{3,Float64},
                               direction::NTuple{3,Float64},
                               det::Detector,
                               cfg::SimConfig,
                               rng::AbstractRNG)::PropagationResult

    vol = active_volume(det)
    fv = fiducial_volume(det)
    mat = vol.material

    pos = Float64[position...]
    dir = Float64[direction...]
    E = E_MeV
    n_int = 0

    while E >= cfg.Egamma_cut
        sC, sP, sPh = sigma_three(mat, E)
        s_tot = sC + sP + sPh
        Σ_tot = mat.n_atom * s_tot

        s = sample_distance(Σ_tot, rng)
        pos .= pos .+ dir .* s
        n_int += 1

        # Escaped the active volume entirely
        if !is_inside(vol, pos)
            return PropagationResult(:lost, E, pos, dir, n_int)
        end

        # Reached the fiducial volume — accept for full simulation
        if is_inside(fv, pos)
            return PropagationResult(:accepted, E, pos, dir, n_int)
        end

        # Determine veto threshold based on region
        veto_threshold = cfg.veto_TPC

        # Sample interaction
        proc = sample_process(sC/s_tot, sP/s_tot, sPh/s_tot, rng)

        if proc === :compton
            Egp, cos_t = sample_compton(E, cfg, rng)
            T_e = E - Egp

            if T_e > veto_threshold
                return PropagationResult(:vetoed, E, pos, dir, n_int)
            end

            # Below veto: gamma continues, deposit is invisible
            ϕ = 2π * rand(rng)
            sin_t = sqrt(max(0.0, 1.0 - cos_t^2))
            local_vec = Float64[sin_t*cos(ϕ), sin_t*sin(ϕ), cos_t]
            dir = rotate_to_global(local_vec, dir)
            E = Egp

        elseif proc === :pair
            if E > veto_threshold
                return PropagationResult(:vetoed, E, pos, dir, n_int)
            end
            return PropagationResult(:lost, 0.0, pos, dir, n_int)

        elseif proc === :photoelectric
            if E > veto_threshold
                return PropagationResult(:vetoed, E, pos, dir, n_int)
            end
            return PropagationResult(:lost, 0.0, pos, dir, n_int)
        end
    end

    PropagationResult(:lost, E, pos, dir, n_int)
end


# =====================================================================
# Multi-volume propagation
# =====================================================================

"""
    classify_lxe_region(pos, det) -> Symbol

Classify the LXe region at `pos`:
- `:fv` — fiducial volume (analysis target)
- `:tpc` — TPC active (outside FV, veto threshold = veto_TPC)
- `:skin` — LXe skin (between field cage and ICV, veto_skin)
- `:passive` — dome or RFR (no veto, passive LXe)
- `:none` — not in any LXe region

Uses volume names from the detector JSON to identify regions.
"""
function classify_lxe_region(pos::Vector{Float64}, det::Detector)::Symbol
    fv = det.fiducial
    if fv !== nothing && is_inside(fv, pos)
        return :fv
    end

    vol = find_volume(det, pos)
    vol === nothing && return :none

    name = lowercase(vol.name)
    if name == "lxetpc"
        return :tpc
    elseif name == "skin"
        return :skin
    elseif name in ("rfr", "dome")
        return :passive
    elseif vol.material.name == "LXe"
        return :passive  # any other LXe volume is passive
    else
        return :none
    end
end


"""
    veto_threshold(region::Symbol, cfg::SimConfig) -> Float64

Return the veto energy threshold for a given LXe region.
`:tpc` → veto_TPC, `:skin` → veto_skin, `:passive` → Inf (no veto),
`:fv` → 0.0 (accept everything).
"""
function veto_threshold(region::Symbol, cfg::SimConfig)::Float64
    region === :tpc && return cfg.veto_TPC
    region === :skin && return cfg.veto_skin
    region === :fv && return 0.0
    region === :passive && return Inf  # no veto possible
    Inf
end


function _resolve_boundary_volume(det::Detector,
                                  pos::Vector{Float64},
                                  dir::Vector{Float64};
                                  tol_cm::Float64=TRANSPORT_BOUNDARY_TOL_CM,
                                  push_cm::Float64=TRANSPORT_BOUNDARY_PUSH_CM)::Union{PhysicalVolume,Nothing}
    vol = find_volume(det, pos)

    if vol === nothing
        if !is_inside(det.mars, pos)
            return nothing
        end
        next_vol, t_next = next_volume(pos, dir, det)
        if next_vol !== nothing && t_next <= tol_cm
            pos .= pos .+ dir .* push_cm
            return find_volume(det, pos)
        end
        return nothing
    end

    t_exit = distance_to_exit(pos, dir, vol.logical)
    if t_exit <= tol_cm
        pos .= pos .+ dir .* push_cm
        return find_volume(det, pos)
    end

    vol
end


"""
    propagate_to_lxe(E_MeV, position, direction, det, cfg, rng)
        -> PropagationResult

Multi-volume photon propagation from any point in MARS toward the
fiducial volume. Handles:

- **Material volumes** (Ti, PTFE): Compton scatter or absorb, with
  per-material cross sections.
- **Vacuum** (between OCV and ICV): straight-line propagation to the
  next volume boundary via `next_volume`.
- **LXe regions**: apply region-specific veto thresholds.
  TPC active: veto_TPC (10 keV). Skin: veto_skin (100 keV).
  Passive (dome/RFR): no veto (Inf threshold).
- **FV**: accept for full simulation.

Returns `PropagationResult` with status `:accepted`, `:vetoed`, or `:lost`.
"""
function propagate_to_lxe(E_MeV::Float64,
                          position::NTuple{3,Float64},
                          direction::NTuple{3,Float64},
                          det::Detector,
                          cfg::SimConfig,
                          rng::AbstractRNG)::PropagationResult

    pos = Float64[position...]
    dir = Float64[direction...]
    E = E_MeV
    n_int = 0

    while E >= cfg.Egamma_cut
        # Where are we?
        vol = _resolve_boundary_volume(det, pos, dir)

        if vol === nothing
            # In vacuum (MARS, between vessels): propagate to next volume
            if !is_inside(det.mars, pos)
                return PropagationResult(:lost, E, pos, dir, n_int)
            end
            next_vol, t_next = next_volume(pos, dir, det)
            if next_vol === nothing || t_next == Inf
                return PropagationResult(:lost, E, pos, dir, n_int)
            end
            # Advance to just inside the next volume
            pos .= pos .+ dir .* (t_next + TRANSPORT_BOUNDARY_PUSH_CM)
            continue
        end

        mat = vol.material

        # Vacuum material (shouldn't happen if volumes are correct, but handle it)
        if mat.density <= 0.0
            # Advance through this vacuum volume
            t_exit = distance_to_exit(pos, dir, vol.logical)
            pos .= pos .+ dir .* (t_exit + TRANSPORT_BOUNDARY_PUSH_CM)
            continue
        end

        # Check if this is an LXe region
        if mat.name == "LXe"
            region = classify_lxe_region(pos, det)

            # Reached FV: accept
            if region === :fv
                return PropagationResult(:accepted, E, pos, dir, n_int)
            end

            # LXe interaction with region-specific veto
            v_thresh = veto_threshold(region, cfg)

            sC, sP, sPh = sigma_three(mat, E)
            s_tot = sC + sP + sPh
            Σ_tot = mat.n_atom * s_tot
            s = sample_distance(Σ_tot, rng)
            pos .= pos .+ dir .* s
            n_int += 1

            # Left this LXe volume?
            if !is_inside(vol, pos)
                _resolve_boundary_volume(det, pos, dir)
                # Re-check where we are now
                continue
            end

            # Re-classify after moving (might have entered FV)
            region = classify_lxe_region(pos, det)
            if region === :fv
                return PropagationResult(:accepted, E, pos, dir, n_int)
            end

            v_thresh = veto_threshold(region, cfg)

            proc = sample_process(sC/s_tot, sP/s_tot, sPh/s_tot, rng)

            if proc === :compton
                Egp, cos_t = sample_compton(E, cfg, rng)
                T_e = E - Egp
                if T_e > v_thresh
                    return PropagationResult(:vetoed, E, pos, dir, n_int)
                end
                ϕ = 2π * rand(rng)
                sin_t = sqrt(max(0.0, 1.0 - cos_t^2))
                local_vec = Float64[sin_t*cos(ϕ), sin_t*sin(ϕ), cos_t]
                dir = rotate_to_global(local_vec, dir)
                E = Egp
            else
                if E > v_thresh
                    return PropagationResult(:vetoed, E, pos, dir, n_int)
                end
                return PropagationResult(:lost, 0.0, pos, dir, n_int)
            end

        else
            # Non-LXe material (Ti, PTFE): interact or traverse
            if mat.xcom === nothing
                # No cross section data: treat as transparent
                t_exit = distance_to_exit(pos, dir, vol.logical)
                pos .= pos .+ dir .* (t_exit + TRANSPORT_BOUNDARY_PUSH_CM)
                continue
            end

            sC, sP, sPh = sigma_three(mat, E)
            s_tot = sC + sP + sPh
            Σ_tot = mat.n_atom * s_tot
            s = sample_distance(Σ_tot, rng)
            pos .= pos .+ dir .* s
            n_int += 1

            # Left this volume? (passed through without interacting)
            if !is_inside(vol, pos)
                _resolve_boundary_volume(det, pos, dir)
                continue
            end

            # Interacted inside this material
            proc = sample_process(sC/s_tot, sP/s_tot, sPh/s_tot, rng)

            if proc === :compton
                Egp, cos_t = sample_compton(E, cfg, rng)
                ϕ = 2π * rand(rng)
                sin_t = sqrt(max(0.0, 1.0 - cos_t^2))
                local_vec = Float64[sin_t*cos(ϕ), sin_t*sin(ϕ), cos_t]
                dir = rotate_to_global(local_vec, dir)
                E = Egp
            else
                return PropagationResult(:lost, 0.0, pos, dir, n_int)
            end
        end
    end

    PropagationResult(:lost, E, pos, dir, n_int)
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
