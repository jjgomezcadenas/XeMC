#!/usr/bin/env julia
"""
run_test.jl — Driver script for LXeMC simulation.

Simulates N gamma-ray events in a liquid xenon volume and reports
summary statistics (deposited energy, SS/MS fraction, timing).
Optionally saves per-event data in long-format CSV.

Usage:
    julia --project=.. scripts/run_test.jl -n 1000 -e 2.615
    julia --project=.. scripts/run_test.jl -n 500 -e 2.448 --save-events -o results/
    julia --project=.. scripts/run_test.jl -n 100 -e 2.615 --geometry cylinder,65,65
"""

using ArgParse
using LXeMC
using Random
using Statistics
using Printf

# =====================================================================
# CLI argument parsing
# =====================================================================

function parse_args()
    s = ArgParseSettings(
        description = "LXeMC: Monte Carlo simulation of gammas in liquid xenon",
        version = "0.1.0",
        add_version = true
    )

    @add_arg_table! s begin
        "--nevents", "-n"
            help = "Number of events to simulate"
            arg_type = Int
            required = true
        "--energy", "-e"
            help = "Primary gamma energy [MeV]"
            arg_type = Float64
            required = true
        "--seed", "-s"
            help = "RNG seed"
            arg_type = Int
            default = 42
        "--output-dir", "-o"
            help = "Output directory (default: current directory)"
            arg_type = String
            default = "."
        "--save-events"
            help = "Save per-event data to events.csv"
            action = :store_true
        "--geometry", "-g"
            help = "Geometry: 'infinite' or 'cylinder,R_cm,H_cm' (half-height)"
            arg_type = String
            default = "infinite"
        "--mode", "-m"
            help = "Simulation mode: 'full' (with electron transport) or 'photon-only'"
            arg_type = String
            default = "full"
    end

    return ArgParse.parse_args(s)
end

# =====================================================================
# Geometry parser
# =====================================================================

function build_geometry(spec::String)::Geometry
    spec = strip(spec)
    if lowercase(spec) == "infinite"
        return InfiniteLXe()
    end
    parts = split(spec, ",")
    if length(parts) == 3 && lowercase(strip(parts[1])) == "cylinder"
        R = parse(Float64, strip(parts[2]))
        H = parse(Float64, strip(parts[3]))
        return CylinderLXe(R, H)
    end
    error("Unknown geometry '$spec'. Use 'infinite' or 'cylinder,R,H'.")
end

# =====================================================================
# Main simulation loop
# =====================================================================

function run_simulation(args)
    N       = args["nevents"]
    E_MeV   = args["energy"]
    seed    = args["seed"]
    outdir  = args["output-dir"]
    save_ev = args["save-events"]
    geom    = build_geometry(args["geometry"])
    mode    = lowercase(strip(args["mode"]))
    mode in ("full", "photon-only") || error("Unknown mode '$mode'. Use 'full' or 'photon-only'.")
    photon_only = mode == "photon-only"

    mkpath(outdir)

    cfg = default_config()
    nd  = load_nist_data(cfg)
    rng = MersenneTwister(seed)

    # Pre-allocate per-event storage
    E_deps     = Vector{Float64}(undef, N)
    n_deposits = Vector{Int}(undef, N)
    n_clusters = Vector{Int}(undef, N)
    is_ss      = Vector{Bool}(undef, N)

    # Cluster-level storage (for CSV output)
    if save_ev
        ev_rows = Vector{NTuple{9, Any}}()  # will collect tuples
    end

    println("=" ^ 60)
    println("LXeMC Simulation")
    println("=" ^ 60)
    @printf("  Events:     %d\n", N)
    @printf("  Energy:     %.4f MeV\n", E_MeV)
    @printf("  Seed:       %d\n", seed)
    @printf("  Geometry:   %s\n", args["geometry"])
    @printf("  Mode:       %s\n", mode)
    @printf("  Output dir: %s\n", outdir)
    @printf("  Save CSV:   %s\n", save_ev ? "yes" : "no")
    println("-" ^ 60)

    t_start = time()

    for i in 1:N
        deposits = photon_only ?
            simulate_event_photon_only(E_MeV, nd, cfg; geom=geom, rng=rng) :
            simulate_event(E_MeV, nd, cfg; geom=geom, rng=rng)
        clusters = cluster_deposits_in_z(deposits, cfg.dz_resolution;
                                         E_min=cfg.E_cluster_min)

        E_dep = sum(d.energy for d in deposits; init=0.0)
        nc = length(clusters)
        ss = nc <= 1

        E_deps[i]     = E_dep
        n_deposits[i] = length(deposits)
        n_clusters[i] = nc
        is_ss[i]      = ss

        if save_ev
            for (ci, (z_cent, E_cl)) in enumerate(clusters)
                push!(ev_rows, (i, E_MeV, E_dep, length(deposits), nc, ss, ci, z_cent, E_cl))
            end
        end

        # Progress indicator
        if N >= 100 && i % (N ÷ 10) == 0
            @printf("  ... %d/%d events (%.0f%%)\n", i, N, 100.0 * i / N)
        end
    end

    t_elapsed = time() - t_start

    # ---- Summary statistics ----
    n_ss = count(is_ss)
    ss_frac = n_ss / N

    summary_lines = String[]
    push!(summary_lines, "=" ^ 60)
    push!(summary_lines, "Results")
    push!(summary_lines, "=" ^ 60)
    push!(summary_lines, @sprintf("  Events simulated:    %d", N))
    push!(summary_lines, @sprintf("  Primary energy:      %.4f MeV", E_MeV))
    push!(summary_lines, @sprintf("  Geometry:            %s", args["geometry"]))
    push!(summary_lines, @sprintf("  Mode:                %s", mode))
    push!(summary_lines, @sprintf("  Seed:                %d", seed))
    push!(summary_lines, "-" ^ 60)
    push!(summary_lines, @sprintf("  Mean E_deposited:    %.4f MeV", mean(E_deps)))
    push!(summary_lines, @sprintf("  Std  E_deposited:    %.4f MeV", std(E_deps)))
    push!(summary_lines, @sprintf("  Min  E_deposited:    %.4f MeV", minimum(E_deps)))
    push!(summary_lines, @sprintf("  Max  E_deposited:    %.4f MeV", maximum(E_deps)))
    push!(summary_lines, @sprintf("  Mean energy loss:    %.4f keV", (E_MeV - mean(E_deps)) * 1000))
    push!(summary_lines, "-" ^ 60)
    push!(summary_lines, @sprintf("  SS events:           %d / %d  (%.1f%%)", n_ss, N, ss_frac * 100))
    push!(summary_lines, @sprintf("  MS events:           %d / %d  (%.1f%%)", N - n_ss, N, (1 - ss_frac) * 100))
    push!(summary_lines, @sprintf("  Mean clusters/event: %.2f", mean(n_clusters)))
    push!(summary_lines, @sprintf("  Mean deposits/event: %.1f", mean(n_deposits)))
    push!(summary_lines, "-" ^ 60)
    push!(summary_lines, @sprintf("  Wall time:           %.2f s", t_elapsed))
    push!(summary_lines, @sprintf("  Events/sec:          %.0f", N / t_elapsed))
    push!(summary_lines, "=" ^ 60)

    # Print to stdout
    for line in summary_lines
        println(line)
    end

    # Write summary.txt
    summary_path = joinpath(outdir, "summary.txt")
    open(summary_path, "w") do io
        for line in summary_lines
            println(io, line)
        end
    end
    println("\nSummary written to: $summary_path")

    # ---- Per-event CSV ----
    if save_ev
        csv_path = joinpath(outdir, "events.csv")
        open(csv_path, "w") do io
            println(io, "event,E_MeV,E_dep,n_deposits,n_clusters,is_SS,cluster_idx,z_centroid,E_cluster")
            for row in ev_rows
                (ev, emev, edep, ndep, nc, ss, ci, zc, ecl) = row
                @printf(io, "%d,%.4f,%.6f,%d,%d,%s,%d,%.6f,%.6f\n",
                        ev, emev, edep, ndep, nc, ss ? "true" : "false", ci, zc, ecl)
            end
        end
        println("Events CSV written to: $csv_path  ($(length(ev_rows)) rows)")
    end
end

# =====================================================================
# Entry point
# =====================================================================

run_simulation(parse_args())
