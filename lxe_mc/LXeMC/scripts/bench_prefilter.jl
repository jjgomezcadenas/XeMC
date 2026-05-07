#!/usr/bin/env julia
"""
bench_prefilter.jl — Legacy flat-detector benchmark for the old veto pre-filter.

Shoots gammas radially inward from the wall, downward from the top, and
upward from the bottom. Reports veto/acceptance rates and throughput.

This script exercises the historical `propagate_to_fiducial(...)` path,
not the canonical fast-kernel event workflow.

Usage:
    julia --project=.. scripts/bench_prefilter.jl
    julia --project=.. scripts/bench_prefilter.jl 5000000
"""

using LXeMC
using Random
using Printf

const N = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1_000_000
const E0 = 2.615  # Tl-208

function run_surface(gen_pos_dir::Function, label::String, N::Int,
                     det::Detector, cfg::SimConfig)
    rng = MersenneTwister(42)

    # Warmup
    for _ in 1:1000
        pos, dir = gen_pos_dir(rng)
        propagate_to_fiducial(E0, pos, dir, det, cfg, rng)
    end

    rng = MersenneTwister(42)
    n_vetoed = 0; n_accepted = 0; n_lost = 0
    GC.gc()
    t = @elapsed for _ in 1:N
        pos, dir = gen_pos_dir(rng)
        r = propagate_to_fiducial(E0, pos, dir, det, cfg, rng)
        if r.status === :vetoed
            n_vetoed += 1
        elseif r.status === :accepted
            n_accepted += 1
        else
            n_lost += 1
        end
    end

    @printf("  %-8s  %10d events  %6.2f s  %10.0f ev/s  vetoed=%5.1f%%  accepted=%5.2f%%  lost=%5.2f%%\n",
            label, N, t, N/t,
            100n_vetoed/N, 100n_accepted/N, 100n_lost/N)
end

function main()
    cfg = default_config()
    mats = load_materials(cfg)
    det = load_detector(default_detector_path(), mats)

    println("=" ^ 80)
    println("Veto pre-filter benchmark: $(N) gammas × 3 surfaces at $(E0) MeV")
    println("=" ^ 80)

    # Wall: radially inward from r=72, at TPC mid-height
    run_surface("WALL", N, det, cfg) do rng
        φ = 2π * rand(rng)
        ((72.0*cos(φ), 72.0*sin(φ), 72.8),
         (-cos(φ), -sin(φ), 0.0))
    end

    # Top: downward from z=145.0, uniform in disk r < 72.8
    run_surface("TOP", N, det, cfg) do rng
        r = 72.8 * sqrt(rand(rng))
        φ = 2π * rand(rng)
        ((r*cos(φ), r*sin(φ), 145.0),
         (0.0, 0.0, -1.0))
    end

    # Bottom: upward from z=0.5, uniform in disk r < 72.8
    run_surface("BOTTOM", N, det, cfg) do rng
        r = 72.8 * sqrt(rand(rng))
        φ = 2π * rand(rng)
        ((r*cos(φ), r*sin(φ), 0.5),
         (0.0, 0.0, 1.0))
    end

    println("=" ^ 80)
end

main()
