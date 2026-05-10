"""
Process source flux tables through the detector and collect SS candidates.

Reads pre-generated flux tables from --indir, samples gammas, runs them
through the fast kernel + FV stack, and writes:
  - candidates.csv: SS candidates with (x, y, z, E_deposited)
  - statistics.csv: event status counts and fractions

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

    (source=source, isotope=isotope, indir=indir, outdir=outdir, N=N, seed=seed)
end


# =====================================================================
# Per-thread processing
# =====================================================================

struct CandidateEvent
    x::Float64
    y::Float64
    z::Float64
    E_deposited::Float64
end


struct ThreadResult
    candidates::Vector{CandidateEvent}
    n_accepted::Int
    n_vetoed::Int
    n_vetoed_tpc::Int
    n_vetoed_skin::Int
    n_ms_rejected::Int
    n_no_fv::Int
    n_total::Int
end


function process_batch(flux, surface_sampler::Function,
                        fk::FastKernelGeometry, fv::PhysicalVolume,
                        cfg::SimConfig, rng::AbstractRNG,
                        N::Int, isotope::String)::ThreadResult
    candidates = CandidateEvent[]
    n_accepted = 0
    n_vetoed = 0
    n_vetoed_tpc = 0
    n_vetoed_skin = 0
    n_ms_rejected = 0
    n_no_fv = 0

    for _ in 1:N
        # Sample gamma(s) from flux table
        if isotope == "Bi214"
            g = sample_gamma_from_flux(flux, surface_sampler, rng)
            gammas = [g]
        else
            gammas = sample_gamma_from_flux(flux, surface_sampler, rng)
        end

        # Process through fast kernel + FV stack
        result = process_event(gammas, fk, fv, cfg, rng)

        if result.status == :accepted
            n_accepted += 1
            # Compute energy-weighted centroid and total energy
            E_total = sum(d.energy for d in result.deposits)
            if E_total > 0.0
                x = sum(d.position[1] * d.energy for d in result.deposits) / E_total
                y = sum(d.position[2] * d.energy for d in result.deposits) / E_total
                z = sum(d.position[3] * d.energy for d in result.deposits) / E_total
            else
                x, y, z = 0.0, 0.0, 0.0
            end
            push!(candidates, CandidateEvent(x, y, z, E_total))
        elseif result.status == :vetoed
            n_vetoed += 1
            result.has_tpc_veto && (n_vetoed_tpc += 1)
            result.has_skin_veto && (n_vetoed_skin += 1)
        elseif result.status == :ms_rejected
            n_ms_rejected += 1
        else
            n_no_fv += 1
        end
    end

    ThreadResult(candidates, n_accepted, n_vetoed, n_vetoed_tpc, n_vetoed_skin,
                 n_ms_rejected, n_no_fv, N)
end


function merge_thread_results(results::Vector{ThreadResult})::ThreadResult
    ThreadResult(
        vcat([r.candidates for r in results]...),
        sum(r.n_accepted for r in results),
        sum(r.n_vetoed for r in results),
        sum(r.n_vetoed_tpc for r in results),
        sum(r.n_vetoed_skin for r in results),
        sum(r.n_ms_rejected for r in results),
        sum(r.n_no_fv for r in results),
        sum(r.n_total for r in results)
    )
end


# =====================================================================
# Output
# =====================================================================

function write_candidates_csv(path::String, candidates::Vector{CandidateEvent})
    open(path, "w") do io
        println(io, "x_cm,y_cm,z_cm,E_deposited_MeV")
        for c in candidates
            @printf(io, "%.6f,%.6f,%.6f,%.8e\n", c.x, c.y, c.z, c.E_deposited)
        end
    end
end


function write_statistics_csv(path::String, r::ThreadResult)
    open(path, "w") do io
        println(io, "quantity,value")
        @printf(io, "n_total,%d\n", r.n_total)
        @printf(io, "n_accepted,%d\n", r.n_accepted)
        @printf(io, "n_vetoed,%d\n", r.n_vetoed)
        @printf(io, "n_vetoed_tpc,%d\n", r.n_vetoed_tpc)
        @printf(io, "n_vetoed_skin,%d\n", r.n_vetoed_skin)
        @printf(io, "n_ms_rejected,%d\n", r.n_ms_rejected)
        @printf(io, "n_no_fv,%d\n", r.n_no_fv)
        @printf(io, "f_accepted,%.8e\n", r.n_accepted / r.n_total)
        @printf(io, "f_vetoed,%.8e\n", r.n_vetoed / r.n_total)
        @printf(io, "f_vetoed_tpc,%.8e\n", r.n_vetoed_tpc / r.n_total)
        @printf(io, "f_vetoed_skin,%.8e\n", r.n_vetoed_skin / r.n_total)
        @printf(io, "f_ms_rejected,%.8e\n", r.n_ms_rejected / r.n_total)
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
    fv = compile_fv_volume(det)
    src_path = normpath(joinpath(dirname(pathof(LXeMC)), "..", "..", "data",
                                 "source_geometry_lz_v1.json"))
    sg = load_source_geometry(src_path, mats)

    surface_sampler = make_surface_sampler(cli.source, sg)

    # --- Load flux tables ---
    isotope_prefix = cli.isotope == "Bi214" ? "bi214" : "tl208"
    flux_keys = [k for k in keys(
        JSON.parsefile(joinpath(cli.indir, "metadata.json"))["components"])
        if startswith(k, isotope_prefix) && !endswith(k, "rate")]

    # For sampling we use the summed rate table's PDF shape
    # (all components combined). Load it as the sampling source.
    rate_key = isotope_prefix * "_rate"
    rate_table = load_rate_table(cli.indir, rate_key)

    # For Bi-214: sample from rate table directly
    # For Tl-208: need the full Tl208 struct for companion sampling.
    #   Load the first non-rate component (they all have same shape after merging).
    flux = if cli.isotope == "Bi214"
        # Build a SourceFluxBi214 from the rate table PDF (treat as probability)
        total = sum(rate_table.pdf_rate)
        pdf_norm = total > 0 ? rate_table.pdf_rate ./ total : rate_table.pdf_rate
        SourceFluxBi214("combined", pdf_norm,
                        rate_table.E_min, rate_table.E_max,
                        rate_table.n_E, rate_table.n_u,
                        0, 0, 0, 0, 0)
    else
        # Load one of the Tl208 component tables for companion structure
        tl_key = first(k for k in flux_keys if !endswith(k, "rate"))
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
            # Thread 1: run in 10 sub-batches with progress
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
    @printf("  Accepted (SS):    %8d  (%.4e)\n", merged.n_accepted, merged.n_accepted / N_actual)
    @printf("  Vetoed:           %8d  (%.4e)\n", merged.n_vetoed, merged.n_vetoed / N_actual)
    @printf("    TPC veto:       %8d  (%.4e)\n", merged.n_vetoed_tpc, merged.n_vetoed_tpc / N_actual)
    @printf("    Skin veto:      %8d  (%.4e)\n", merged.n_vetoed_skin, merged.n_vetoed_skin / N_actual)
    @printf("  MS rejected:      %8d  (%.4e)\n", merged.n_ms_rejected, merged.n_ms_rejected / N_actual)
    @printf("  No FV:            %8d  (%.4e)\n", merged.n_no_fv, merged.n_no_fv / N_actual)
    @printf("  Candidates:       %8d\n", length(merged.candidates))

    # --- Write output ---
    mkpath(cli.outdir)
    write_candidates_csv(joinpath(cli.outdir, "candidates.csv"), merged.candidates)
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
        "n_accepted" => merged.n_accepted,
        "n_vetoed" => merged.n_vetoed,
        "n_ms_rejected" => merged.n_ms_rejected,
        "n_no_fv" => merged.n_no_fv,
        "f_accepted" => merged.n_accepted / N_actual
    )
    write_flux_json(joinpath(cli.outdir, "metadata.json"), metadata)

    @printf("\nOutput written to %s\n", cli.outdir)
    @printf("Total time: %.2f s\n", time() - t0)
end


main()
