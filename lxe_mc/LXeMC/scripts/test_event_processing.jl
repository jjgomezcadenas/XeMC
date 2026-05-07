using LXeMC
using Printf
using Random


function print_help()
    println("""
Usage:
  julia --project=. scripts/test_event_processing.jl [options]
  julia --project=. scripts/test_event_processing.jl [N] [SEED]

Options:
  -h, --help     Show this help message and exit
  --n N          Number of events to process (default: 1000)
  --seed SEED    RNG seed (default: 20260507)
  --path NAME    Event-processing path: v2 or fastkernel (default: v2)
  --outdir DIR   Write output tables into DIR

Examples:
  julia --project=. scripts/test_event_processing.jl
  julia --project=. scripts/test_event_processing.jl --n 10000 --seed 42 --path fastkernel --outdir /tmp/event_proc
  julia --project=. scripts/test_event_processing.jl 10000 42
""")
end


function parse_cli(args)
    N = 1000
    seed = 20260507
    path = "v2"
    outdir = nothing

    i = 1
    positional = String[]
    while i <= length(args)
        arg = args[i]
        if arg == "-h" || arg == "--help"
            print_help()
            exit(0)
        elseif arg == "--n"
            i == length(args) && error("--n requires an integer value")
            N = parse(Int, args[i + 1])
            i += 2
            continue
        elseif arg == "--seed"
            i == length(args) && error("--seed requires an integer value")
            seed = parse(Int, args[i + 1])
            i += 2
            continue
        elseif arg == "--path"
            i == length(args) && error("--path requires a value: v2 or fastkernel")
            path = lowercase(args[i + 1])
            i += 2
            continue
        elseif arg == "--outdir"
            i == length(args) && error("--outdir requires a directory path")
            outdir = args[i + 1]
            i += 2
            continue
        elseif startswith(arg, "--")
            error("Unknown option '$arg'")
        else
            push!(positional, arg)
        end
        i += 1
    end

    if !isempty(positional)
        N = parse(Int, positional[1])
    end
    if length(positional) >= 2
        seed = parse(Int, positional[2])
    end
    length(positional) <= 2 || error("Too many positional arguments")

    N > 0 || error("N must be positive")
    seed >= 0 || error("seed must be non-negative")
    path in ("v2", "fastkernel") || error("--path must be either 'v2' or 'fastkernel' (got '$path')")
    return N, seed, path, outdir
end


function print_counts(title::String, counts::Dict{Symbol,Int})
    println("\n$title")
    for key in sort(collect(keys(counts)); by=String)
        @printf("  %-14s %d\n", String(key), counts[key])
    end
end


function write_counts(path::String, counts::Dict{Symbol,Int})
    open(path, "w") do io
        println(io, "key,count")
        for key in sort(collect(keys(counts)); by=String)
            @printf(io, "%s,%d\n", String(key), counts[key])
        end
    end
end


function write_energy_status_table(path::String, counts::Dict{Tuple{Float64,Float64,Symbol},Int})
    open(path, "w") do io
        println(io, "E1_MeV,E2_MeV,status,count")
        for key in sort(collect(keys(counts)); by=k -> (k[1], k[2], String(k[3])))
            E1, E2, status = key
            @printf(io, "%.8f,%.8f,%s,%d\n", E1, E2, String(status), counts[key])
        end
    end
end


function main()
    N, seed, path, outdir = parse_cli(ARGS)

    cfg = default_config()
    mats = load_materials(cfg)
    det = load_detector_v2(default_detector_v2_path(), mats)
    fk = path == "fastkernel" ? compile_fastkernel_geometry(det) : nothing
    rng = MersenneTwister(seed)

    status_counts = Dict{Symbol,Int}()
    multiplicity_counts = Dict{Symbol,Int}()
    decisive_counts = Dict{Symbol,Int}()
    energy_status_counts = Dict{Tuple{Float64,Float64,Symbol},Int}()
    progress_step = max(1, cld(N, 10))
    t0 = time()

    for i in 1:N
        gammas = sample_event("Tl208", "calib"; calib=true, rng=rng)
        result = if path == "fastkernel"
            process_event_fastkernel(gammas, fk, det, cfg, rng)
        else
            process_event(gammas, det, cfg, rng)
        end

        status_counts[result.status] = get(status_counts, result.status, 0) + 1

        mult_key = Symbol(string(length(gammas)))
        multiplicity_counts[mult_key] = get(multiplicity_counts, mult_key, 0) + 1

        decisive_idx = result.status == :vetoed ? result.n_processed : 0
        dec_key = Symbol(string(decisive_idx))
        decisive_counts[dec_key] = get(decisive_counts, dec_key, 0) + 1

        energies = sort([g.E_MeV for g in gammas]; rev=true)
        E1 = length(energies) >= 1 ? energies[1] : 0.0
        E2 = length(energies) >= 2 ? energies[2] : 0.0
        tbl_key = (E1, E2, result.status)
        energy_status_counts[tbl_key] = get(energy_status_counts, tbl_key, 0) + 1

        if i % progress_step == 0 || i == N
            pct = 100.0 * i / N
            elapsed = time() - t0
            @printf("Processed %d / %d events (%.1f%%) in %.2f s\n", i, N, pct, elapsed)
        end
    end

    elapsed = time() - t0

    println("N = $N, seed = $seed, path = $path")
    @printf("Elapsed time = %.2f s\n", elapsed)
    @printf("Rate = %.2f events/s\n", N / elapsed)
    print_counts("Event status counts", status_counts)
    print_counts("Event multiplicity counts", multiplicity_counts)
    print_counts("First decisive gamma index counts", decisive_counts)

    if outdir !== nothing
        mkpath(outdir)
        println("\nOutput directory: $outdir")
        write_counts(joinpath(outdir, "event_status_counts.csv"), status_counts)
        write_counts(joinpath(outdir, "event_multiplicity_counts.csv"), multiplicity_counts)
        write_counts(joinpath(outdir, "first_decisive_gamma_counts.csv"), decisive_counts)
        write_energy_status_table(joinpath(outdir, "event_energy_status_table.csv"), energy_status_counts)
        open(joinpath(outdir, "summary.txt"), "w") do io
            @printf(io, "N = %d\n", N)
            @printf(io, "seed = %d\n", seed)
            @printf(io, "path = %s\n", path)
            @printf(io, "elapsed_s = %.8f\n", elapsed)
            @printf(io, "rate_events_per_s = %.8f\n", N / elapsed)
            @printf(io, "outdir = %s\n", outdir)
        end
    end
end


main()
