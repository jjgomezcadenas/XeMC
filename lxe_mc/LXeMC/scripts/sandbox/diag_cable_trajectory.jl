"""
Diagnostic: trace cable gammas through the fast kernel + (if FV-handoff) the
full-stack transport. Dumps every deposit's (position, energy, region tag).

Run from LXeMC dir:
    julia --project=. scripts/sandbox/diag_cable_trajectory.jl [N] [seed]

Defaults: N=10 gammas, seed=42, position (R=72.75, z=160), direction 45deg
inward+down (dir=(-1,0,-1)/sqrt(2)).
"""

using LXeMC
using Printf
using Random
using LinearAlgebra


function trace_one(gamma::SampledGamma, fk::FastKernelGeometry,
                   fv_vol, cfg::SimConfig, rng::AbstractRNG)
    println("\n" * "="^80)
    @printf("INITIAL: r=%.2f z=%.2f phi=%.2f deg, E=%.3f MeV, dir=(%.3f, %.3f, %.3f)\n",
            hypot(gamma.position[1], gamma.position[2]),
            gamma.position[3],
            atand(gamma.position[2], gamma.position[1]),
            gamma.E_MeV,
            gamma.direction[1], gamma.direction[2], gamma.direction[3])

    region0 = classify_fastkernel(fk, (gamma.position[1], gamma.position[2], gamma.position[3]))
    if region0 === nothing
        println("  Starting position is OUTSIDE the fast-kernel envelope (MARS region)")
    else
        @printf("  Starting region: %s (tag=%s, material=%s)\n",
                region0.name, region0.tag, region0.material.name)
    end

    result = transport_gamma_fastkernel(gamma, fk, cfg, rng)
    @printf("FAST KERNEL STATUS: %s  terminal_region=%s  E_out=%.3f MeV\n",
            result.status, result.terminal_region, result.energy_MeV)
    @printf("  Final position: r=%.2f z=%.2f\n",
            hypot(result.position[1], result.position[2]),
            result.position[3])
    println("  Fast-kernel deposits (n=$(length(result.deposits))):")
    for (i, dep) in enumerate(result.deposits)
        r = hypot(dep.position[1], dep.position[2])
        @printf("    [%d] r=%.2f z=%.2f  E_dep=%.3f MeV  region=%s\n",
                i, r, dep.position[3], dep.Edep_MeV, dep.region)
    end

    fv_deposits = nothing
    if result.status == :handoff_fv
        println("  >>> Handoff to full-stack transport in FV <<<")
        fv_deposits = propagate_gamma(
            result.energy_MeV, fv_vol, fk, cfg;
            position=(result.position[1], result.position[2], result.position[3]),
            direction=(result.direction[1], result.direction[2], result.direction[3]),
            rng=rng,
        )
        println("  Full-stack FV deposits (n=$(length(fv_deposits))):")
        for (i, dep) in enumerate(fv_deposits)
            r = hypot(dep.position[1], dep.position[2])
            @printf("    [%d] r=%.2f z=%.2f  E_dep=%.3f MeV  src=%s  intr=%s  vol=%s\n",
                    i, r, dep.position[3], dep.energy,
                    dep.source, dep.interaction, dep.volume)
        end
    end

    (result, fv_deposits)
end


function main()
    N    = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 10
    seed = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 42

    cfg = default_config()
    mats = load_materials(cfg)
    det = load_tracking_detector(default_tracking_detector_path(), mats)
    fk = compile_fastkernel_geometry(det)
    fv_vol = compile_fv_volume(det)

    @printf("Fast-kernel geometry compiled: %d regions, FV at z=[%.1f, %.1f] r<=%.1f\n",
            length(fk.regions),
            fv_vol.logical.position[3] - fv_vol.logical.solid.half_height_cm,
            fv_vol.logical.position[3] + fv_vol.logical.solid.half_height_cm,
            fv_vol.logical.solid.radius_cm)
    @printf("Veto thresholds: TPC=%.3f MeV  Skin=%.3f MeV\n",
            cfg.veto_TPC, cfg.veto_skin)

    # Starting condition: cable shell at R=78, z=155, 45 deg inward+down
    R0     = 78.0
    z0     = 155.0
    E0     = 2.448  # Bi-214 line
    dir_in_down = [-1.0, 0.0, -1.0] / sqrt(2.0)
    @printf("\nDiagnostic: %d gammas at (R=%.2f, z=%.2f), E=%.3f MeV, dir=(%.3f, 0, %.3f) (45 deg inward+down)\n",
            N, R0, z0, E0, dir_in_down[1], dir_in_down[3])

    rng = MersenneTwister(seed)
    fv_event_count = 0
    handoff_count = 0
    fv_E_sum = 0.0

    for i in 1:N
        position = Float64[R0, 0.0, z0]
        gamma = SampledGamma(E0, copy(position), copy(dir_in_down))
        result, fv_deps = trace_one(gamma, fk, fv_vol, cfg, rng)

        if result.status == :handoff_fv
            handoff_count += 1
        end
        if fv_deps !== nothing
            for dep in fv_deps
                if dep.volume === :fv
                    fv_event_count += 1
                    fv_E_sum += dep.energy
                end
            end
        end
    end

    println("\n" * "="^80)
    @printf("SUMMARY: N=%d, handoff_fv=%d, FV deposits=%d, total FV energy=%.3f MeV\n",
            N, handoff_count, fv_event_count, fv_E_sum)
end


main()
