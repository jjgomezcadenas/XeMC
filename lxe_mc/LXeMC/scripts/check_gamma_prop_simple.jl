using LXeMC
using LinearAlgebra
using Printf
using Random


function print_help()
    println("""
Usage:
  julia --project=. scripts/check_gamma_prop_simple.jl [options]
  julia --project=. scripts/check_gamma_prop_simple.jl [N] [SEED]

Options:
  -h, --help     Show this help message and exit
  --n N          Number of gammas to shoot (default: 1000)
  --seed SEED    RNG seed (default: 20260507)
  --energy E     Gamma energy in MeV (default: 2.61)
  --x X          Initial x position in cm (default: 0.0)
  --y Y          Initial y position in cm (default: 0.0)
  --z Z          Initial z position in cm (default: 106.667)
  --ux UX        Direction cosine x (default: 0.0)
  --uy UY        Direction cosine y (default: 0.0)
  --uz UZ        Direction cosine z (default: -1.0)
  --outdir DIR   Write summary and histogram tables into DIR

Examples:
  julia --project=. scripts/check_gamma_prop_simple.jl
  julia --project=. scripts/check_gamma_prop_simple.jl --n 5000 --seed 42 --energy 2.615 --x 0 --y 0 --z 106.667 --ux 0 --uy 0 --uz -1
  julia --project=. scripts/check_gamma_prop_simple.jl --n 5000 --seed 42 --outdir /tmp/gamma_prop
  julia --project=. scripts/check_gamma_prop_simple.jl 5000 42
""")
end


function parse_cli(args)
    N = 1000
    seed = 20260507
    E0 = 2.61
    x0 = 0.0
    y0 = 0.0
    z0 = 106.67
    ux = 0.0
    uy = 0.0
    uz = -1.0
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
        elseif arg == "--energy"
            i == length(args) && error("--energy requires a floating-point value")
            E0 = parse(Float64, args[i + 1])
            i += 2
            continue
        elseif arg == "--x"
            i == length(args) && error("--x requires a floating-point value")
            x0 = parse(Float64, args[i + 1])
            i += 2
            continue
        elseif arg == "--y"
            i == length(args) && error("--y requires a floating-point value")
            y0 = parse(Float64, args[i + 1])
            i += 2
            continue
        elseif arg == "--z"
            i == length(args) && error("--z requires a floating-point value")
            z0 = parse(Float64, args[i + 1])
            i += 2
            continue
        elseif arg == "--ux"
            i == length(args) && error("--ux requires a floating-point value")
            ux = parse(Float64, args[i + 1])
            i += 2
            continue
        elseif arg == "--uy"
            i == length(args) && error("--uy requires a floating-point value")
            uy = parse(Float64, args[i + 1])
            i += 2
            continue
        elseif arg == "--uz"
            i == length(args) && error("--uz requires a floating-point value")
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
    E0 > 0.0 || error("energy must be positive")

    u = Float64[ux, uy, uz]
    norm_u = norm(u)
    norm_u > 0.0 || error("direction cosines must not all be zero")
    u ./= norm_u

    return N, seed, E0, Float64[x0, y0, z0], u, outdir
end


function hist1d(values::Vector{Float64}, lo::Float64, hi::Float64, nbins::Int)
    edges = collect(range(lo, hi; length=nbins + 1))
    counts = zeros(Int, nbins)
    width = (hi - lo) / nbins
    for v in values
        lo <= v <= hi || continue
        idx = v == hi ? nbins : clamp(fld(Int(floor((v - lo) / width)), 1) + 1, 1, nbins)
        counts[idx] += 1
    end
    edges, counts
end


function hist2d(xs::Vector{Float64},
                ys::Vector{Float64},
                xlo::Float64,
                xhi::Float64,
                ylo::Float64,
                yhi::Float64,
                nbx::Int,
                nby::Int)
    xedges = collect(range(xlo, xhi; length=nbx + 1))
    yedges = collect(range(ylo, yhi; length=nby + 1))
    counts = zeros(Int, nby, nbx)
    xw = (xhi - xlo) / nbx
    yw = (yhi - ylo) / nby

    for (x, y) in zip(xs, ys)
        xlo <= x <= xhi || continue
        ylo <= y <= yhi || continue
        ix = x == xhi ? nbx : clamp(fld(Int(floor((x - xlo) / xw)), 1) + 1, 1, nbx)
        iy = y == yhi ? nby : clamp(fld(Int(floor((y - ylo) / yw)), 1) + 1, 1, nby)
        counts[iy, ix] += 1
    end
    xedges, yedges, counts
end


function print_hist1d(title::String, edges, counts)
    println("\n$title")
    for i in eachindex(counts)
        @printf("  [%7.3f, %7.3f): %d\n", edges[i], edges[i + 1], counts[i])
    end
end


function print_counts(title::String, counts::Dict{Symbol,Int})
    println("\n$title")
    for key in sort(collect(keys(counts)); by=String)
        @printf("  %-14s %d\n", String(key), counts[key])
    end
end


function print_heatmap_rz(xedges, yedges, counts)
    println("\nR-Z heatmap counts (rows=z bins, cols=r bins)")
    header = join([@sprintf("%6.1f", 0.5 * (xedges[i] + xedges[i + 1])) for i in 1:length(xedges)-1], " ")
    println("        ", header)
    for iy in reverse(1:size(counts, 1))
        zmid = 0.5 * (yedges[iy] + yedges[iy + 1])
        row = join([@sprintf("%6d", counts[iy, ix]) for ix in 1:size(counts, 2)], " ")
        @printf("%6.1f  %s\n", zmid, row)
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


function write_hist1d(path::String, edges, counts)
    open(path, "w") do io
        println(io, "bin_lo,bin_hi,count")
        for i in eachindex(counts)
            @printf(io, "%.8f,%.8f,%d\n", edges[i], edges[i + 1], counts[i])
        end
    end
end


function write_heatmap(path::String, xedges, yedges, counts)
    open(path, "w") do io
        println(io, "z_bin_lo,z_bin_hi,r_bin_lo,r_bin_hi,count")
        for iy in 1:size(counts, 1), ix in 1:size(counts, 2)
            @printf(io, "%.8f,%.8f,%.8f,%.8f,%d\n",
                    yedges[iy], yedges[iy + 1], xedges[ix], xedges[ix + 1], counts[iy, ix])
        end
    end
end


fmt2(x::Real) = @sprintf("%.2f", Float64(x))


function cli_tag_dirname(E0, position0, direction0)
    join([
        "e_$(fmt2(E0))",
        "x_$(fmt2(position0[1]))",
        "y_$(fmt2(position0[2]))",
        "z_$(fmt2(position0[3]))",
        "ux_$(fmt2(direction0[1]))",
        "uy_$(fmt2(direction0[2]))",
        "uz_$(fmt2(direction0[3]))",
    ], "_")
end


function main()
    N, seed, E0, position0, direction0, outdir = parse_cli(ARGS)

    cfg = default_config()
    mats = load_materials(cfg)
    det = load_tracking_detector(default_tracking_detector_path(), mats)
    rng = MersenneTwister(seed)

    results = GammaPropagationV2Result[]
    for _ in 1:N
        gamma = SampledGamma(E0, copy(position0), copy(direction0))
        push!(results, propagate_gamma_v2(gamma, det, cfg, rng))
    end

    status_counts = Dict{Symbol,Int}()
    interaction_counts = Dict{Symbol,Int}()
    selection_counts = Dict{Symbol,Int}()
    deposits = Float64[]
    rs = Float64[]
    zs = Float64[]

    for r in results
        status_counts[r.status] = get(status_counts, r.status, 0) + 1
        interaction_counts[r.interaction_type] = get(interaction_counts, r.interaction_type, 0) + 1
        sel = select_interaction(r, det)
        sel_key = if sel.class == :fv
            :fv_pass
        elseif sel.class == :tpc
            sel.passes_threshold ? :tpc_pass : :tpc_fail
        elseif sel.class == :skin
            sel.passes_threshold ? :skin_pass : :skin_fail
        else
            :other
        end
        selection_counts[sel_key] = get(selection_counts, sel_key, 0) + 1
        if r.status == :interacted
            push!(deposits, r.deposit_E_MeV)
            push!(rs, hypot(r.position[1], r.position[2]))
            push!(zs, r.position[3])
        end
    end

    @printf("N = %d, seed = %d\n", N, seed)
    @printf("E0 = %.3f MeV\n", E0)
    @printf("Start position: (%.3f, %.3f, %.3f) cm\n", position0[1], position0[2], position0[3])
    @printf("Direction cosines: (%.6f, %.6f, %.6f)\n", direction0[1], direction0[2], direction0[3])

    if outdir !== nothing
        tagged_outdir = joinpath(outdir, cli_tag_dirname(E0, position0, direction0))
        mkpath(tagged_outdir)
        println("Base output directory: $outdir")
        println("Tagged output directory: $tagged_outdir")
        write_counts(joinpath(tagged_outdir, "status_counts.csv"), status_counts)
        write_counts(joinpath(tagged_outdir, "interaction_type_counts.csv"), interaction_counts)
        write_counts(joinpath(tagged_outdir, "selection_outcome_counts.csv"), selection_counts)
        open(joinpath(tagged_outdir, "summary.txt"), "w") do io
            @printf(io, "N = %d\n", N)
            @printf(io, "seed = %d\n", seed)
            @printf(io, "E0_MeV = %.8f\n", E0)
            @printf(io, "x0_cm = %.8f\n", position0[1])
            @printf(io, "y0_cm = %.8f\n", position0[2])
            @printf(io, "z0_cm = %.8f\n", position0[3])
            @printf(io, "ux = %.8f\n", direction0[1])
            @printf(io, "uy = %.8f\n", direction0[2])
            @printf(io, "uz = %.8f\n", direction0[3])
            @printf(io, "outdir = %s\n", tagged_outdir)
        end
    else
        print_counts("Status counts", status_counts)
        print_counts("Interaction-type counts", interaction_counts)
        print_counts("Selection outcome counts", selection_counts)
    end

    if !isempty(deposits)
        e_edges, e_counts = hist1d(deposits, 0.0, E0, 20)
        z_edges, z_counts = hist1d(zs, 0.0, 145.6, 20)
        r_edges, r_counts = hist1d(rs, 0.0, 82.1, 20)
        hz_r_edges, hz_z_edges, hz_counts = hist2d(rs, zs, 0.0, 82.1, 0.0, 145.6, 12, 12)

        if outdir !== nothing
            tagged_outdir = joinpath(outdir, cli_tag_dirname(E0, position0, direction0))
            write_hist1d(joinpath(tagged_outdir, "deposit_energy_hist.csv"), e_edges, e_counts)
            write_hist1d(joinpath(tagged_outdir, "interaction_z_hist.csv"), z_edges, z_counts)
            write_hist1d(joinpath(tagged_outdir, "interaction_r_hist.csv"), r_edges, r_counts)
            write_heatmap(joinpath(tagged_outdir, "interaction_rz_heatmap.csv"), hz_r_edges, hz_z_edges, hz_counts)
        else
            print_hist1d("Deposited energy [MeV]", e_edges, e_counts)
            print_hist1d("Interaction z [cm]", z_edges, z_counts)
            print_hist1d("Interaction r [cm]", r_edges, r_counts)
            print_heatmap_rz(hz_r_edges, hz_z_edges, hz_counts)
        end
    end
end


main()
