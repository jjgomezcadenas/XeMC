"""
Process source flux tables through the detector and collect FV deposits.

Reads pre-generated flux tables from --indir, samples gammas, runs them
through the fast kernel + FV stack, and writes:
  - fv_deposits.csv: all FV deposits with (event_id, x, y, z, energy, source)
  - statistics.csv: event status counts and fractions

SS/MS classification is done offline (in python) from the deposit CSV.

Usage:
    julia -t 8 --project=. scripts/run_source_backgrounds.jl \\
        --source cryostat_barrel --isotope Bi214 \\
        --indir results/fluxes/cryostat/barrel/Bi214 \\
        --outdir results/backgrounds/cryostat/barrel/Bi214 \\
        --n 1000000 --seed 42
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
  julia -t N --project=. scripts/run_source_backgrounds.jl [options]

Required:
  --source SOURCE    Source identifier ($(join(supported_sources(), ", ")))
  --isotope ISOTOPE  Isotope ($(join(supported_isotopes(), ", ")))
  --indir DIR        Input directory with flux tables
  --outdir DIR       Output directory for results

Optional:
  --n N              Number of events to sample (default: 100000)
  --seed SEED        RNG seed (default: 42)
  --fv_r R           FV radius [cm] (default: from detector JSON)
  --fv_zmin ZMIN     FV z_min [cm] (default: from detector JSON)
  --fv_zmax ZMAX     FV z_max [cm] (default: from detector JSON)
  -h, --help         Show this help
""")
end


function parse_cli(args)
    source = nothing
    isotope = nothing
    indir = nothing
    outdir = nothing
    N = 100_000
    seed = 42
    fv_r = nothing
    fv_zmin = nothing
    fv_zmax = nothing

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
        elseif arg == "--indir"
            indir = args[i+1]; i += 2; continue
        elseif arg == "--outdir"
            outdir = args[i+1]; i += 2; continue
        elseif arg == "--n"
            N = parse(Int, args[i+1]); i += 2; continue
        elseif arg == "--seed"
            seed = parse(Int, args[i+1]); i += 2; continue
        elseif arg == "--fv_r"
            fv_r = parse(Float64, args[i+1]); i += 2; continue
        elseif arg == "--fv_zmin"
            fv_zmin = parse(Float64, args[i+1]); i += 2; continue
        elseif arg == "--fv_zmax"
            fv_zmax = parse(Float64, args[i+1]); i += 2; continue
        else
            error("Unknown argument '$arg'")
        end
        i += 1
    end

    source === nothing && error("--source is required")
    isotope === nothing && error("--isotope is required")
    indir === nothing && error("--indir is required")
    outdir === nothing && error("--outdir is required")
    source in supported_sources() || error("Unknown source '$source'")
    isotope in supported_isotopes() || error("Unknown isotope '$isotope'")
    N > 0 || error("--n must be positive")

    (source=source, isotope=isotope, indir=indir, outdir=outdir, N=N, seed=seed,
     fv_r=fv_r, fv_zmin=fv_zmin, fv_zmax=fv_zmax)
end


# =====================================================================
# Per-thread processing
# =====================================================================

struct ThreadResult
    fv_events::Vector{Vector{Deposit}}
    n_fv::Int
    n_vetoed::Int
    n_vetoed_tpc::Int
    n_vetoed_skin::Int
    n_no_fv::Int
    n_total::Int
end


function process_batch(flux, surface_sampler::Function,
                        fk::FastKernelGeometry, fv::PhysicalVolume,
                        cfg::SimConfig, rng::AbstractRNG,
                        N::Int, isotope::String)::ThreadResult
    fv_events = FVEvent[]
    n_fv = 0
    n_vetoed = 0
    n_vetoed_tpc = 0
    n_vetoed_skin = 0
    n_no_fv = 0

    for _ in 1:N
        if isotope == "Bi214"
            g = sample_gamma_from_flux(flux, surface_sampler, rng)
            gammas = [g]
        else
            gammas = sample_gamma_from_flux(flux, surface_sampler, rng)
        end

        result = process_event(gammas, fk, fv, cfg, rng)

        if result.status == :fv
            n_fv += 1
            push!(fv_events, result.deposits)
        elseif result.status == :vetoed
            n_vetoed += 1
            result.has_tpc_veto && (n_vetoed_tpc += 1)
            result.has_skin_veto && (n_vetoed_skin += 1)
        else
            n_no_fv += 1
        end
    end

    ThreadResult(fv_events, n_fv, n_vetoed, n_vetoed_tpc, n_vetoed_skin, n_no_fv, N)
end


function merge_thread_results(results::Vector{ThreadResult})::ThreadResult
    ThreadResult(
        vcat([r.fv_events for r in results]...),
        sum(r.n_fv for r in results),
        sum(r.n_vetoed for r in results),
        sum(r.n_vetoed_tpc for r in results),
        sum(r.n_vetoed_skin for r in results),
        sum(r.n_no_fv for r in results),
        sum(r.n_total for r in results)
    )
end


# =====================================================================
# Output
# =====================================================================

function write_fv_deposits_csv(path::String, fv_events::Vector{Vector{Deposit}})
    open(path, "w") do io
        println(io, "event_id,x_cm,y_cm,z_cm,energy_MeV,source")
        for (eid, deps) in enumerate(fv_events)
            for d in deps
                @printf(io, "%d,%.6f,%.6f,%.6f,%.8e,%s\n",
                        eid, d.position[1], d.position[2], d.position[3],
                        d.energy, d.source)
            end
        end
    end
end


function write_statistics_csv(path::String, r::ThreadResult)
    open(path, "w") do io
        println(io, "quantity,value")
        @printf(io, "n_total,%d\n", r.n_total)
        @printf(io, "n_fv,%d\n", r.n_fv)
        @printf(io, "n_vetoed,%d\n", r.n_vetoed)
        @printf(io, "n_vetoed_tpc,%d\n", r.n_vetoed_tpc)
        @printf(io, "n_vetoed_skin,%d\n", r.n_vetoed_skin)
        @printf(io, "n_no_fv,%d\n", r.n_no_fv)
        @printf(io, "f_fv,%.8e\n", r.n_fv / r.n_total)
        @printf(io, "f_vetoed,%.8e\n", r.n_vetoed / r.n_total)
        @printf(io, "f_vetoed_tpc,%.8e\n", r.n_vetoed_tpc / r.n_total)
        @printf(io, "f_vetoed_skin,%.8e\n", r.n_vetoed_skin / r.n_total)
        @printf(io, "f_no_fv,%.8e\n", r.n_no_fv / r.n_total)
    end
end


# =====================================================================
# Main
# =====================================================================

function main()
    cli = parse_cli(ARGS)
    nt = nthreads()
    @printf("Source: %s, Isotope: %s, N: %d, Seed: %d, Threads: %d\n",
            cli.source, cli.isotope, cli.N, cli.seed, nt)
    @printf("Input:  %s\n", cli.indir)
    @printf("Output: %s\n", cli.outdir)

    t0 = time()

    # --- Setup ---
    cfg = default_config()
    mats = load_materials(cfg)
    det = load_tracking_detector(default_tracking_detector_path(), mats)
    fk = compile_fastkernel_geometry(det)

    # FV: use CLI overrides or detector JSON defaults
    fv_default = compile_fv_volume(det)
    if cli.fv_r !== nothing || cli.fv_zmin !== nothing || cli.fv_zmax !== nothing
        # Override: use CLI values, fall back to defaults for unspecified
        r = cli.fv_r !== nothing ? cli.fv_r : fv_default.logical.solid.radius_cm
        zmin = cli.fv_zmin !== nothing ? cli.fv_zmin :
               fv_default.logical.position[3] - fv_default.logical.solid.half_height_cm
        zmax = cli.fv_zmax !== nothing ? cli.fv_zmax :
               fv_default.logical.position[3] + fv_default.logical.solid.half_height_cm
        hh = (zmax - zmin) / 2.0
        zc = (zmax + zmin) / 2.0
        fv = PCyl("FV", LCyl(Cyl(r, hh), Float64[0.0, 0.0, zc]), fv_default.material)
        @printf("FV override: R=%.1f, z=[%.1f, %.1f]\n", r, zmin, zmax)
    else
        fv = fv_default
        r = fv.logical.solid.radius_cm
        zmin = fv.logical.position[3] - fv.logical.solid.half_height_cm
        zmax = fv.logical.position[3] + fv.logical.solid.half_height_cm
        @printf("FV default:  R=%.1f, z=[%.1f, %.1f]\n", r, zmin, zmax)
    end
    fv_r_used = fv.logical.solid.radius_cm
    fv_zmin_used = fv.logical.position[3] - fv.logical.solid.half_height_cm
    fv_zmax_used = fv.logical.position[3] + fv.logical.solid.half_height_cm

    src_path = normpath(joinpath(dirname(pathof(LXeMC)), "..", "..", "data",
                                 "source_geometry_lz_v1.json"))
    sg = load_source_geometry(src_path, mats)

    surface_sampler = make_surface_sampler(cli.source, sg)

    # --- Load flux tables ---
    isotope_prefix = cli.isotope == "Bi214" ? "bi214" : "tl208"
    rate_key = isotope_prefix * "_rate"
    rate_table = load_rate_table(cli.indir, rate_key)

    flux = if cli.isotope == "Bi214"
        total = sum(rate_table.pdf_rate)
        pdf_norm = total > 0 ? rate_table.pdf_rate ./ total : rate_table.pdf_rate
        SourceFluxBi214("combined", pdf_norm,
                        rate_table.E_min, rate_table.E_max,
                        rate_table.n_E, rate_table.n_u,
                        0, 0, 0, 0, 0)
    else
        # Load one Tl208 component for companion structure
        meta = JSON.parsefile(joinpath(cli.indir, "metadata.json"))
        tl_key = first(k for k in keys(meta["components"])
                       if startswith(k, "tl208") && !endswith(k, "rate"))
        load_flux_tl208(cli.indir, tl_key)
    end

    @printf("Setup time: %.2f s\n", time() - t0)
    @printf("Rate table: %.6e gammas/s\n", rate_table.total_rate)

    # --- Parallel processing ---
    t1 = time()
    N_per_thread = cld(cli.N, nt)
    thread_results = Vector{ThreadResult}(undef, nt)

    @threads for tid in 1:nt
        rng = MersenneTwister(cli.seed + tid)
        if tid == 1
            n_chunks = 10
            chunk_size = cld(N_per_thread, n_chunks)
            sub_results = ThreadResult[]
            for chunk in 1:n_chunks
                n_this = min(chunk_size, N_per_thread - (chunk - 1) * chunk_size)
                n_this <= 0 && break
                push!(sub_results, process_batch(flux, surface_sampler, fk, fv,
                                                  cfg, rng, n_this, cli.isotope))
                done = min(chunk * chunk_size, N_per_thread)
                @printf("  Progress: %d/%d (%.0f%%) [%.1fs]\n",
                        done, N_per_thread, 100.0 * done / N_per_thread, time() - t1)
            end
            thread_results[tid] = merge_thread_results(sub_results)
        else
            thread_results[tid] = process_batch(flux, surface_sampler, fk, fv,
                                                 cfg, rng, N_per_thread, cli.isotope)
        end
    end

    merged = merge_thread_results(thread_results)
    elapsed = time() - t1
    N_actual = N_per_thread * nt
    @printf("Processing: %.2f s (%d events, %.0f events/s)\n",
            elapsed, N_actual, N_actual / elapsed)

    # --- Summary ---
    @printf("\nResults:\n")
    @printf("  FV events:        %8d  (%.4e)\n", merged.n_fv, merged.n_fv / N_actual)
    @printf("  Vetoed:           %8d  (%.4e)\n", merged.n_vetoed, merged.n_vetoed / N_actual)
    @printf("    TPC veto:       %8d  (%.4e)\n", merged.n_vetoed_tpc, merged.n_vetoed_tpc / N_actual)
    @printf("    Skin veto:      %8d  (%.4e)\n", merged.n_vetoed_skin, merged.n_vetoed_skin / N_actual)
    @printf("  No FV:            %8d  (%.4e)\n", merged.n_no_fv, merged.n_no_fv / N_actual)
    n_deposits = sum(length(deps) for deps in merged.fv_events; init=0)
    @printf("  FV deposits:      %8d  (%.1f per FV event)\n",
            n_deposits, merged.n_fv > 0 ? n_deposits / merged.n_fv : 0.0)

    # --- Write output ---
    mkpath(cli.outdir)
    write_fv_deposits_csv(joinpath(cli.outdir, "fv_deposits.csv"), merged.fv_events)
    write_statistics_csv(joinpath(cli.outdir, "statistics.csv"), merged)

    metadata = Dict{String,Any}(
        "source" => cli.source,
        "isotope" => cli.isotope,
        "indir" => cli.indir,
        "N_requested" => cli.N,
        "N_actual" => N_actual,
        "seed" => cli.seed,
        "nthreads" => nt,
        "elapsed_processing_s" => elapsed,
        "rate_total_gammas_per_s" => rate_table.total_rate,
        "n_fv" => merged.n_fv,
        "n_vetoed" => merged.n_vetoed,
        "n_no_fv" => merged.n_no_fv,
        "n_fv_deposits" => n_deposits,
        "f_fv" => merged.n_fv / N_actual,
        "fv_r_cm" => fv_r_used,
        "fv_zmin_cm" => fv_zmin_used,
        "fv_zmax_cm" => fv_zmax_used
    )
    write_flux_json(joinpath(cli.outdir, "metadata.json"), metadata)

    @printf("\nOutput written to %s\n", cli.outdir)
    @printf("Total time: %.2f s\n", time() - t0)
end


main()
