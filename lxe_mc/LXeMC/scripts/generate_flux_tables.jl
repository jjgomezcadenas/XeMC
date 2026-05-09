"""
Generate source flux tables for one (source, isotope) pair.

Usage:
    julia -t 8 --project=. scripts/generate_flux_tables.jl \\
        --source cryostat_barrel --isotope Bi214 \\
        --n 100000 --seed 42 \\
        --outdir results/fluxes/cryostat/barrel/Bi214

Parallelization: use `julia -t N` to set thread count.
The script splits N events across threads with independent RNGs.
"""

using LXeMC
using JSON
using Random
using Printf
using Base.Threads


# =====================================================================
# CLI
# =====================================================================

function print_help()
    println("""
Usage:
  julia -t N --project=. scripts/generate_flux_tables.jl [options]

Required:
  --source SOURCE    Source identifier ($(join(supported_sources(), ", ")))
  --isotope ISOTOPE  Isotope ($(join(supported_isotopes(), ", ")))
  --outdir DIR       Output directory

Optional:
  --n N              Number of decays to generate (default: 100000)
  --seed SEED        RNG seed (default: 42)
  -h, --help         Show this help
""")
end


function parse_cli(args)
    source = nothing
    isotope = nothing
    outdir = nothing
    N = 100_000
    seed = 42

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg in ("-h", "--help")
            print_help()
            exit(0)
        elseif arg == "--source"
            source = args[i+1]; i += 2; continue
        elseif arg == "--isotope"
            isotope = args[i+1]; i += 2; continue
        elseif arg == "--outdir"
            outdir = args[i+1]; i += 2; continue
        elseif arg == "--n"
            N = parse(Int, args[i+1]); i += 2; continue
        elseif arg == "--seed"
            seed = parse(Int, args[i+1]); i += 2; continue
        else
            error("Unknown argument '$arg'")
        end
        i += 1
    end

    source === nothing && error("--source is required")
    isotope === nothing && error("--isotope is required")
    outdir === nothing && error("--outdir is required")
    source in supported_sources() || error("Unknown source '$source'")
    isotope in supported_isotopes() || error("Unknown isotope '$isotope'")
    N > 0 || error("--n must be positive")

    (source=source, isotope=isotope, outdir=outdir, N=N, seed=seed)
end


# =====================================================================
# Main
# =====================================================================

function main()
    cli = parse_cli(ARGS)
    nt = nthreads()
    @printf("Source: %s, Isotope: %s, N: %d, Seed: %d, Threads: %d\n",
            cli.source, cli.isotope, cli.N, cli.seed, nt)

    t0 = time()

    cfg = default_config()
    mats = load_materials(cfg)
    src_path = normpath(joinpath(dirname(pathof(LXeMC)), "..", "..", "data",
                                 "source_geometry_lz_v1.json"))
    sg = load_source_geometry(src_path, mats)
    @printf("Setup time: %.2f s\n", time() - t0)

    # --- Parallel generation ---
    t1 = time()
    N_per_thread = cld(cli.N, nt)
    thread_results = Vector{Dict{Symbol,Any}}(undef, nt)

    @threads for tid in 1:nt
        rng = MersenneTwister(cli.seed + tid)
        thread_results[tid] = dispatch_source_flux(
            cli.source, cli.isotope, N_per_thread, sg, cfg, rng;
            verbose=(tid == 1))
    end

    merged = merge_dispatch_results(cli.isotope, thread_results)
    elapsed_gen = time() - t1
    N_actual = N_per_thread * nt
    @printf("Generation: %.2f s (%d events, %.0f events/s)\n",
            elapsed_gen, N_actual, N_actual / elapsed_gen)

    # --- Write output ---
    mkpath(cli.outdir)

    metadata = Dict{String,Any}(
        "source" => cli.source,
        "isotope" => cli.isotope,
        "N_requested" => cli.N,
        "N_actual" => N_actual,
        "seed" => cli.seed,
        "nthreads" => nt,
        "elapsed_generation_s" => elapsed_gen,
        "components" => Dict{String,Any}()
    )

    for (key, val) in merged
        prefix = string(key)
        if val isa SourceFluxBi214
            write_flux_bi214_csv(cli.outdir, prefix, val)
            metadata["components"][prefix] = flux_bi214_metadata(val)
            @printf("  %-20s  survival = %.6f  N_gen = %d\n",
                    prefix, sum(val.pdf), val.N_generated)
        elseif val isa SourceFluxTl208
            write_flux_tl208_csv(cli.outdir, prefix, val)
            metadata["components"][prefix] = flux_tl208_metadata(val)
            @printf("  %-20s  survival = %.6f  N_gen = %d\n",
                    prefix, sum(val.pdf_main), val.N_generated)
        elseif val isa SourceRateTable
            write_rate_table_csv(cli.outdir, prefix, val)
            metadata["components"][prefix] = rate_table_metadata(val)
            @printf("  %-20s  total_rate = %.6e gammas/s\n",
                    prefix, val.total_rate)
        end
    end

    # --- Activity summary ---
    activity_key = cli.isotope == "Bi214" ? "Bi214_mBq_per_kg" : "Tl208_mBq_per_kg"
    chain_name = cli.isotope == "Bi214" ? "U-238" : "Th-232"
    gamma_BR = cli.isotope == "Bi214" ? BR_BI214_2448 : BR_TL208_2615
    rate_key = cli.isotope == "Bi214" ? :bi214_rate : :tl208_rate

    if haskey(merged, rate_key)
        rt = merged[rate_key]
        activity_summary = Dict{String,Any}[]
        total_chain_Bq = 0.0
        total_gamma_emitted = 0.0

        for (i, cname) in enumerate(rt.component_names)
            sv = sg[cname]
            act_mBq = get(sv.activity, activity_key, 0.0)
            act_Bq = act_mBq * 1e-3 * sv.mass_kg
            gamma_emitted = act_Bq * gamma_BR
            total_chain_Bq += act_Bq
            total_gamma_emitted += gamma_emitted
            push!(activity_summary, Dict{String,Any}(
                "component" => cname,
                "mass_kg" => sv.mass_kg,
                "activity_mBq_per_kg" => act_mBq,
                "chain_activity_Bq" => act_Bq,
                "gamma_emitted_per_s" => gamma_emitted,
                "gamma_surviving_per_s" => rt.component_rates[i]
            ))
        end

        metadata["chain"] = chain_name
        metadata["gamma_BR"] = gamma_BR
        metadata["total_chain_activity_Bq"] = total_chain_Bq
        metadata["total_gamma_emitted_per_s"] = total_gamma_emitted
        metadata["total_gamma_surviving_per_s"] = rt.total_rate
        metadata["activity_breakdown"] = activity_summary

        @printf("\nActivity summary (%s chain, BR = %.4f):\n", chain_name, gamma_BR)
        @printf("  %-15s  %10s  %12s  %12s  %12s  %12s\n",
                "Component", "Mass [kg]", "Act [mBq/kg]", "Chain [Bq]",
                "Emitted [/s]", "Surviving [/s]")
        for d in activity_summary
            @printf("  %-15s  %10.2f  %12.3f  %12.6f  %12.6e  %12.6e\n",
                    d["component"], d["mass_kg"], d["activity_mBq_per_kg"],
                    d["chain_activity_Bq"], d["gamma_emitted_per_s"],
                    d["gamma_surviving_per_s"])
        end
        @printf("  %-15s  %10s  %12s  %12.6f  %12.6e  %12.6e\n",
                "TOTAL", "", "", total_chain_Bq, total_gamma_emitted, rt.total_rate)
    end

    write_flux_json(joinpath(cli.outdir, "metadata.json"), metadata)
    @printf("\nOutput written to %s\n", cli.outdir)
    @printf("Total time: %.2f s\n", time() - t0)
end


main()
