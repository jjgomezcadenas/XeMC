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
  --calib BOOL   Calibration sampling mode: true or false (default: true)
  --energy E     Calibration gamma energy in MeV (default: 2.615)
  --x X          Calibration source x in cm (default: 0.0)
  --y Y          Calibration source y in cm (default: 0.0)
  --z Z          Calibration source z in cm (default: 160.0)
  --ux UX        Calibration direction x-component (default: 1.0)
  --uy UY        Calibration direction y-component (default: 0.0)
  --uz UZ        Calibration direction z-component (default: 0.0)
  --outdir DIR   Write output tables into DIR

Examples:
  julia --project=. scripts/test_event_processing.jl
  julia --project=. scripts/test_event_processing.jl --n 10000 --seed 42 --calib true --outdir /tmp/event_proc
  julia --project=. scripts/test_event_processing.jl 10000 42
""")
end


function parse_bool_arg(s::AbstractString)::Bool
    v = lowercase(String(s))
    v == "true" && return true
    v == "false" && return false
    error("Expected true or false, got '$s'")
end


function parse_cli(args)
    N = 1000
    seed = 20260507
    calib = true
    energy = 2.615
    x_cm = 0.0
    y_cm = 0.0
    z_cm = 160.0
    ux = 1.0
    uy = 0.0
    uz = 0.0
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
        elseif arg == "--calib"
            i == length(args) && error("--calib requires true or false")
            calib = parse_bool_arg(args[i + 1])
            i += 2
            continue
        elseif arg == "--energy"
            i == length(args) && error("--energy requires a numeric value")
            energy = parse(Float64, args[i + 1])
            i += 2
            continue
        elseif arg == "--x"
            i == length(args) && error("--x requires a numeric value")
            x_cm = parse(Float64, args[i + 1])
            i += 2
            continue
        elseif arg == "--y"
            i == length(args) && error("--y requires a numeric value")
            y_cm = parse(Float64, args[i + 1])
            i += 2
            continue
        elseif arg == "--z"
            i == length(args) && error("--z requires a numeric value")
            z_cm = parse(Float64, args[i + 1])
            i += 2
            continue
        elseif arg == "--ux"
            i == length(args) && error("--ux requires a numeric value")
            ux = parse(Float64, args[i + 1])
            i += 2
            continue
        elseif arg == "--uy"
            i == length(args) && error("--uy requires a numeric value")
            uy = parse(Float64, args[i + 1])
            i += 2
            continue
        elseif arg == "--uz"
            i == length(args) && error("--uz requires a numeric value")
            uz = parse(Float64, args[i + 1])
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
    energy > 0.0 || error("--energy must be positive")
    calib || error("Only --calib true is supported in CleanupV3")
    if !calib
        cli_state_requested = energy != 2.615 || x_cm != 0.0 || y_cm != 0.0 || z_cm != 160.0 ||
                              ux != 1.0 || uy != 0.0 || uz != 0.0
        !cli_state_requested || error("With --calib false, do not pass --energy/--x/--y/--z/--ux/--uy/--uz")
    end
    return N, seed, calib, energy, x_cm, y_cm, z_cm, ux, uy, uz, outdir
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
    N, seed, calib, energy, x_cm, y_cm, z_cm, ux, uy, uz, outdir = parse_cli(ARGS)

    cfg = default_config()
    mats = load_materials(cfg)
    det_legacy = load_detector(default_detector_path(), mats)
    det = load_tracking_detector(default_tracking_detector_path(), mats)
    fk = compile_fastkernel_geometry(det)
    rng = MersenneTwister(seed)

    status_counts = Dict{Symbol,Int}()
    multiplicity_counts = Dict{Symbol,Int}()
    decisive_counts = Dict{Symbol,Int}()
    energy_status_counts = Dict{Tuple{Float64,Float64,Symbol},Int}()
    progress_step = max(1, cld(N, 10))
    t0 = time()

    for i in 1:N
        fused = process_event_fastkernel_calib(
            fk, det_legacy, cfg, rng;
            E_MeV=energy,
            x_cm=x_cm, y_cm=y_cm, z_cm=z_cm,
            ux=ux, uy=uy, uz=uz,
        )
        result = fused.result
        multiplicity = fused.multiplicity
        E1 = fused.E1_MeV
        E2 = fused.E2_MeV

        status_counts[result.status] = get(status_counts, result.status, 0) + 1

        mult_key = Symbol(string(multiplicity))
        multiplicity_counts[mult_key] = get(multiplicity_counts, mult_key, 0) + 1

        decisive_idx = result.status == :vetoed ? result.n_processed : 0
        dec_key = Symbol(string(decisive_idx))
        decisive_counts[dec_key] = get(decisive_counts, dec_key, 0) + 1

        tbl_key = (E1, E2, result.status)
        energy_status_counts[tbl_key] = get(energy_status_counts, tbl_key, 0) + 1

        if i % progress_step == 0 || i == N
            pct = 100.0 * i / N
            elapsed = time() - t0
            @printf("Processed %d / %d events (%.1f%%) in %.2f s\n", i, N, pct, elapsed)
        end
    end

    elapsed = time() - t0

    println("N = $N, seed = $seed, calib = $calib, path = fastkernel_fused")
    if calib
        @printf("Calibration source: E=%.3f MeV, pos=(%.3f, %.3f, %.3f) cm, dir=(%.3f, %.3f, %.3f)\n",
                energy, x_cm, y_cm, z_cm, ux, uy, uz)
    end
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
            @printf(io, "calib = %s\n", string(calib))
            @printf(io, "path = fastkernel_fused\n")
            if calib
                @printf(io, "energy_MeV = %.8f\n", energy)
                @printf(io, "x_cm = %.8f\n", x_cm)
                @printf(io, "y_cm = %.8f\n", y_cm)
                @printf(io, "z_cm = %.8f\n", z_cm)
                @printf(io, "ux = %.8f\n", ux)
                @printf(io, "uy = %.8f\n", uy)
                @printf(io, "uz = %.8f\n", uz)
            end
            @printf(io, "elapsed_s = %.8f\n", elapsed)
            @printf(io, "rate_events_per_s = %.8f\n", N / elapsed)
            @printf(io, "outdir = %s\n", outdir)
        end
    end
end


main()
