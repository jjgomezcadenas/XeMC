"""
Fast-kernel gamma transport with analysis-driven shortcuts.

Transports gammas through the detector regions using the compiled
`FastKernelGeometry` partition. Unlike full-stack transport (tracking.jl),
this code applies veto thresholds and ROI-floor termination during
propagation for rapid event rejection.

Called by `process_event` (event.jl) for each gamma before FV handoff.
"""


"""
    FastGammaDeposit

One local gamma energy deposit produced during fast-kernel transport.
`region` is the detector region where the interaction occurred.
These deposits drive veto decisions internally; they are not returned
to the caller of `process_event`.
"""
struct FastGammaDeposit
    Edep_MeV::Float64
    position::Vector{Float64}
    region::String
end


"""
    FastGammaTrackResult

Result of fast-kernel photon transport for one gamma.

`status` values:
- `:handoff_fv` — gamma reached the FV, hand off to full stack
- `:escaped` — gamma left the detector
- `:below_cut` — gamma energy fell below Egamma_cut
- `:absorbed_passive` — absorbed in passive LXe or structural material
- `:below_roi_fv` — energy below ROI floor, collapsed locally
- `:vetoed_tpc` — deposit in TPC active region above veto threshold
- `:vetoed_skin` — deposit in Skin above veto threshold
"""
struct FastGammaTrackResult
    status::Symbol
    deposits::Vector{FastGammaDeposit}
    terminal_region::String
    energy_MeV::Float64
    position::Vector{Float64}
    direction::Vector{Float64}
end


@inline function _is_passive_region(tag::RegionTag)::Bool
    tag == TAG_PASSIVE_LXE || tag == TAG_STRUCTURAL
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


"""Return the appropriate terminal status for a gamma absorbed in the region with `tag`."""
@inline function _terminal_status(deposits::Vector{FastGammaDeposit},
                                   tag::RegionTag,
                                   name::String,
                                   pos::Vector{Float64},
                                   dir::Vector{Float64})
    status = if tag == TAG_SKIN
        :vetoed_skin
    elseif tag == TAG_TPC_ACTIVE
        :vetoed_tpc
    elseif _is_passive_region(tag)
        :absorbed_passive
    else
        :below_cut
    end
    _fast_track_result(status, deposits, name, 0.0, pos, dir)
end


"""
    transport_gamma_fastkernel(gamma, fk, cfg, rng; max_cm=400.0)
        -> FastGammaTrackResult

Fast-kernel photon transport with analysis-driven early exits.
Continues after Compton scatters (unlike the single-step diagnostic
`propagate_gamma_fastkernel` in tracking_util.jl).

Applies:
- immediate veto when a deposit in Skin or TPC active exceeds threshold
- ROI-floor termination when surviving energy < `cfg.E_roi_floor`
- FV handoff when the gamma enters the fiducial volume
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
        region = classify_fastkernel(fk, (pos[1], pos[2], pos[3]))
        region === nothing && return _fast_track_result(:escaped, deposits, "MARS", E, pos, dir)

        if region.tag == TAG_FV
            return _fast_track_result(:handoff_fv, deposits, region.name, E, pos, dir)
        end

        mat = region.material
        s_bnd = distance_to_boundary_fastkernel(region, (pos[1], pos[2], pos[3]), (dir[1], dir[2], dir[3]))
        if !isfinite(s_bnd)
            return _fast_track_result(:escaped, deposits, region.name, E, pos, dir)
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
            rname = region.name

            proc = sample_process(sC/s_tot, sP/s_tot, sPh/s_tot, rng)
            if proc === :compton
                Egp, cos_theta = sample_compton(E, cfg, rng)
                dep = E - Egp
                push!(deposits, FastGammaDeposit(dep, copy(pos), rname))

                # Veto check on Compton deposit
                if region.tag == TAG_SKIN && dep >= cfg.veto_skin
                    return _fast_track_result(:vetoed_skin, deposits, rname, Egp, pos, dir)
                elseif region.tag == TAG_TPC_ACTIVE && dep >= cfg.veto_TPC
                    return _fast_track_result(:vetoed_tpc, deposits, rname, Egp, pos, dir)
                end

                # Gamma survives with reduced energy
                E = Egp

                # ROI floor: gamma can no longer reach Q_ββ — collapse locally
                if E < cfg.E_roi_floor
                    dep_total = deposits[end].Edep_MeV + E
                    _replace_last_deposit!(deposits,
                        FastGammaDeposit(dep_total, deposits[end].position,
                                         deposits[end].region))
                    return _terminal_status(deposits, region.tag, rname, pos, dir)
                end

                # Update direction and continue
                sin_theta = sqrt(max(0.0, 1.0 - cos_theta^2))
                phi = 2π * rand(rng)
                local_dir = Float64[sin_theta * cos(phi), sin_theta * sin(phi), cos_theta]
                dir .= rotate_to_global(local_dir, dir)
                continue

            else
                # Pair or photoelectric: terminal, full energy deposited
                push!(deposits, FastGammaDeposit(E, copy(pos), rname))
                return _terminal_status(deposits, region.tag, rname, pos, dir)
            end
        end

        pos .= pos .+ dir .* (s_bnd + TRANSPORT_BOUNDARY_PUSH_CM)
        traveled += s_bnd + TRANSPORT_BOUNDARY_PUSH_CM
    end

    _fast_track_result(:below_cut, deposits, "MARS", E, pos, dir)
end
