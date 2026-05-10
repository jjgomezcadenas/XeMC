"""
Diagnostic and validation utilities for fast-kernel tracking.

These functions are NOT part of the production event pipeline. They
exist for geometry validation, transport debugging, and single-step
analysis. The production pipeline uses `transport_gamma_fastkernel`
(continued multi-scatter transport) called from `process_event`.

Placed in a separate file to keep `tracking.jl` focused on production
code.
"""


"""
    GammaPropagationResult

Result of single-step gamma propagation through the fast kernel.
Used by `propagate_gamma_fastkernel` for diagnostic purposes.
"""
struct GammaPropagationResult
    status::Symbol
    interaction_type::Symbol
    deposit_E_MeV::Float64
    position::Vector{Float64}
    region::String
end


"""
    propagate_gamma_fastkernel(gamma, fk, cfg, rng; max_cm=400.0)
        -> GammaPropagationResult

Single-step diagnostic transport: propagates one gamma through the
fast kernel until the FIRST interaction (Compton, pair, photoelectric),
FV entry, or escape. Unlike `transport_gamma_fastkernel` (production),
this function does NOT continue after Compton scatters.

Used for:
- Geometry validation (checking region classification)
- Transport debugging (verifying cross sections and distances)
- Single-interaction analysis scripts

Not used in the production pipeline (`process_event`).
"""
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


"""
    FastKernelCalibEventResult

Result of a calibration event processed through the fast kernel.
Wraps `EventProcessingResult` with calibration-specific metadata.
Used by `process_event_fastkernel_calib` for diagnostic/validation.
"""
struct FastKernelCalibEventResult
    result::EventProcessingResult
    multiplicity::Int
    E1_MeV::Float64
    E2_MeV::Float64
end


"""
    process_event_fastkernel_calib(fk, fv, cfg, rng; kwargs...) -> FastKernelCalibEventResult

Convenience wrapper for calibration events: samples gammas via
`sample_gammas("calib", ...)`, processes through `process_event`,
and packs the result with multiplicity and gamma energies.

Not used in the production pipeline.
"""
function process_event_fastkernel_calib(fk::FastKernelGeometry, fv::PhysicalVolume, cfg::SimConfig,
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
