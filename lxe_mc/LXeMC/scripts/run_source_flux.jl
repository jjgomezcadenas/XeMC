#!/usr/bin/env julia
"""
run_source_flux.jl — Legacy flat-detector source-flux generation.

Runs Stage 1 on the historical flat detector model:
decays in source material → photon-only propagation →
veto logic → binned flux table in (E, cos θ).

Outputs:
  - summary.txt: statistics
  - flux_table.csv: 2D histogram (E bins × u bins)
  - spectrum_E.csv: dN/dE marginal (25 bins)
  - spectrum_u.csv: dN/du marginal (10 bins)

Usage:
    julia --project=.. scripts/run_source_flux.jl --source ICV_barrel --decay Bi214 -n 1000000
    julia --project=.. scripts/run_source_flux.jl --source OCV_barrel --decay Tl208 -n 500000 -o results/
"""

using ArgParse
using LXeMC
using Random
using Printf

function parse_args()
    s = ArgParseSettings(
        description = "LXeMC: source flux generation (Stage 1)",
        version = "0.1.0",
        add_version = true
    )

    @add_arg_table! s begin
        "--nevents", "-n"
            help = "Number of decay events to generate"
            arg_type = Int
            required = true
        "--source", "-s"
            help = "Source volume name (as in detector JSON)"
            arg_type = String
            required = true
        "--decay", "-d"
            help = "Decay type: Bi214 or Tl208"
            arg_type = String
            required = true
        "--seed"
            help = "RNG seed"
            arg_type = Int
            default = 42
        "--output-dir", "-o"
            help = "Output directory"
            arg_type = String
            default = "."
        "--detector"
            help = "Path to detector JSON"
            arg_type = String
            default = ""
    end

    return ArgParse.parse_args(s)
end


function main()
    args = parse_args()

    N      = args["nevents"]
    src_name = args["source"]
    decay_name = args["decay"]
    seed   = args["seed"]
    outdir = args["output-dir"]
    det_path = args["detector"]

    if isempty(det_path)
        det_path = default_detector_path()
    end

    mkpath(outdir)

    cfg  = default_config()
    mats = load_materials(cfg)
    det  = load_detector(det_path, mats)
    decays = load_decays()

    # Find source volume
    by_name = Dict(v.name => v for v in det.volumes)
    haskey(by_name, src_name) || error("Source '$src_name' not found. Available: $(join(keys(by_name), ", "))")
    source = by_name[src_name]

    # Find decay scheme
    haskey(decays, decay_name) || error("Decay '$decay_name' not found. Available: $(join(keys(decays), ", "))")
    scheme = decays[decay_name]

    rng = MersenneTwister(seed)

    println("=" ^ 60)
    println("Source Flux Generation (Stage 1)")
    println("=" ^ 60)
    @printf("  Source:     %s (%s, %.1f kg)\n", src_name, source.material.name, mass(source)/1000)
    @printf("  Decay:      %s (%.3f MeV)\n", decay_name, scheme.gammas[1].E_MeV)
    @printf("  Events:     %d\n", N)
    @printf("  Seed:       %d\n", seed)
    @printf("  Output:     %s\n", outdir)
    println("-" ^ 60)

    show_progress(i, n) = @printf("  ... %d/%d (%.0f%%)\n", i, n, 100.0*i/n)
    t = @elapsed ft = generate_source_flux(N, source, scheme, det, cfg, rng;
                                            progress=show_progress)

    # Summary
    lines = String[]
    push!(lines, "=" ^ 60)
    push!(lines, "Results")
    push!(lines, "=" ^ 60)
    push!(lines, @sprintf("  Source:          %s", src_name))
    push!(lines, @sprintf("  Decay:           %s", decay_name))
    push!(lines, @sprintf("  Generated:       %d", ft.N_generated))
    push!(lines, @sprintf("  Surviving:       %d (%.4f%%)", ft.N_surviving, 100*survival_fraction(ft)))
    push!(lines, @sprintf("  Vetoed:          %d (%.1f%%)", ft.N_vetoed, 100*ft.N_vetoed/N))
    push!(lines, @sprintf("  Low energy:      %d (%.1f%%)", ft.N_low_energy, 100*ft.N_low_energy/N))
    push!(lines, @sprintf("  Absorbed:        %d (%.1f%%)", ft.N_absorbed, 100*ft.N_absorbed/N))
    push!(lines, @sprintf("  Invisible:       %d (%.1f%%)", ft.N_invisible, 100*ft.N_invisible/N))
    push!(lines, @sprintf("  Backward:        %d (%.1f%%)", ft.N_backward, 100*ft.N_backward/N))
    push!(lines, "-" ^ 60)
    E_line = scheme.gammas[1].E_MeV
    f_peak = peak_bin_fraction(ft, E_line)
    f_off  = off_peak_fraction(ft, E_line)
    f_surv = survival_fraction(ft)
    push!(lines, @sprintf("  Peak bin (%.3f MeV): %.4e per decay (%.1f%% of surviving)", E_line, f_peak, f_surv > 0 ? 100*f_peak/f_surv : 0))
    push!(lines, @sprintf("  Off-peak (in window): %.4e per decay (%.1f%% of surviving)", f_off, f_surv > 0 ? 100*f_off/f_surv : 0))
    push!(lines, @sprintf("  Total survival prob:  %.4e per decay", f_surv))
    push!(lines, "-" ^ 60)
    push!(lines, @sprintf("  E range:         [%.3f, %.3f] MeV (%d bins)", ft.E_min, ft.E_max, ft.n_E))
    push!(lines, @sprintf("  cos θ range:     [0, 1] (%d bins)", ft.n_u))
    push!(lines, @sprintf("  Wall time:       %.2f s", t))
    push!(lines, @sprintf("  Events/sec:      %.0f", N/t))
    push!(lines, "=" ^ 60)

    for l in lines; println(l); end

    # Write summary
    open(joinpath(outdir, "summary.txt"), "w") do io
        for l in lines; println(io, l); end
    end

    # Write flux table (2D)
    open(joinpath(outdir, "flux_table.csv"), "w") do io
        # Header: E_center, u1, u2, ..., u10
        u_c = u_centers(ft)
        print(io, "E_MeV")
        for u in u_c; @printf(io, ",u=%.2f", u); end
        println(io)
        E_c = E_centers(ft)
        for i in 1:ft.n_E
            @printf(io, "%.4f", E_c[i])
            for j in 1:ft.n_u
                @printf(io, ",%.6e", ft.pdf[i, j])
            end
            println(io)
        end
    end

    # Write dP/dE (marginal energy probability)
    open(joinpath(outdir, "spectrum_E.csv"), "w") do io
        println(io, "E_MeV,prob_per_decay")
        E_c = E_centers(ft)
        sp = spectrum_E(ft)
        for i in 1:ft.n_E
            @printf(io, "%.4f,%.6e\n", E_c[i], sp[i])
        end
    end

    # Write dP/du (marginal angular probability)
    open(joinpath(outdir, "spectrum_u.csv"), "w") do io
        println(io, "cos_theta,prob_per_decay")
        u_c = u_centers(ft)
        sp = spectrum_u(ft)
        for i in 1:ft.n_u
            @printf(io, "%.2f,%.6e\n", u_c[i], sp[i])
        end
    end

    println("\nOutput files:")
    println("  $(joinpath(outdir, "summary.txt"))")
    println("  $(joinpath(outdir, "flux_table.csv"))")
    println("  $(joinpath(outdir, "spectrum_E.csv"))")
    println("  $(joinpath(outdir, "spectrum_u.csv"))")
end

main()
