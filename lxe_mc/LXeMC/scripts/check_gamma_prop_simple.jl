using LXeMC
using Printf
using Random


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


function main()
    N = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1000
    seed = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 20260507

    cfg = default_config()
    mats = load_materials(cfg)
    det = load_detector_v2(default_detector_v2_path(), mats)
    rng = MersenneTwister(seed)

    E0 = 2.615
    sC, _, sPh = sigma_three(mats["LXe"], E0)
    λ = 1.0 / (mats["LXe"].n_atom * (sC + sPh))
    start_z = 96.0 + λ

    results = GammaPropagationV2Result[]
    for _ in 1:N
        gamma = SampledGamma(E0, Float64[0.0, 0.0, start_z], Float64[0.0, 0.0, -1.0])
        push!(results, propagate_gamma_v2(gamma, det, cfg, rng))
    end

    status_counts = Dict{Symbol,Int}()
    interaction_counts = Dict{Symbol,Int}()
    deposits = Float64[]
    rs = Float64[]
    zs = Float64[]

    for r in results
        status_counts[r.status] = get(status_counts, r.status, 0) + 1
        interaction_counts[r.interaction_type] = get(interaction_counts, r.interaction_type, 0) + 1
        if r.status == :interacted
            push!(deposits, r.deposit_E_MeV)
            push!(rs, hypot(r.position[1], r.position[2]))
            push!(zs, r.position[3])
        end
    end

    @printf("N = %d, seed = %d\n", N, seed)
    @printf("Start position: (0, 0, %.3f) cm, E0 = %.3f MeV\n", start_z, E0)
    @printf("LXe first-interaction mean free path (Compton+photo only): %.3f cm\n", λ)

    print_counts("Status counts", status_counts)
    print_counts("Interaction-type counts", interaction_counts)

    if !isempty(deposits)
        e_edges, e_counts = hist1d(deposits, 0.0, E0, 20)
        z_edges, z_counts = hist1d(zs, 0.0, 145.6, 20)
        r_edges, r_counts = hist1d(rs, 0.0, 82.1, 20)
        hz_r_edges, hz_z_edges, hz_counts = hist2d(rs, zs, 0.0, 82.1, 0.0, 145.6, 12, 12)

        print_hist1d("Deposited energy [MeV]", e_edges, e_counts)
        print_hist1d("Interaction z [cm]", z_edges, z_counts)
        print_hist1d("Interaction r [cm]", r_edges, r_counts)
        print_heatmap_rz(hz_r_edges, hz_z_edges, hz_counts)
    end
end


main()
