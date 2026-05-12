"""
Analyze FV deposit data: active veto, SS/MS classification, energy
smearing, ROI cut, and background rate estimate.

Reads fv_deposits.csv and metadata.json from --indir, writes CSV/JSON
results to --indir/analysis/ (or --outdir if specified).

Outputs:
  - event_summary.csv: per-event Emax, Etot, DZ, R, Z
  - classification.csv: per-event veto/SS/MS/E_cluster/E_smeared
  - esat.csv: energy of every non-Emax deposit
  - separated_deposits.csv: energy and dz of every separated deposit
  - etot_by_volume.csv: per-event energy split by volume
  - summary.json: all statistics, parameters, background rate

Usage:
    julia -t 8 --project=. scripts/analyze_backgrounds.jl \\
        --indir results/bfv/cryostat/barrel/Bi214
"""

using JSON
using Printf
using Random
using Base.Threads

# =====================================================================
# Constants
# =====================================================================

const Q_BB_KEV = 2458.07
const FWHM_FACTOR = 2.3548200
const SECONDS_PER_YEAR = 3.15576e7

# =====================================================================
# CLI
# =====================================================================

function print_help()
    println("""
Usage:
  julia -t N --project=. scripts/analyze_backgrounds.jl [options]

Required:
  --indir DIR        Input directory with fv_deposits.csv and metadata.json

Optional:
  --outdir DIR       Output directory (default: <indir>/analysis)
  --sigma SIGMA      Relative energy resolution sigma/E (default: 0.01)
  --dz DZ            Min z-separation for cluster splitting, mm (default: 3.0)
  --ks KS            Significance factor (default: 3.0)
  --eveto EVETO      Active veto threshold, keV (default: 10.0)
  --roi LOW HIGH     ROI bounds in keV (default: 1 FWHM centered on Qbb)
  --seed SEED        RNG seed for smearing (default: 42)
  -h, --help         Show this help
""")
end

function parse_cli(args)
    indir = nothing
    outdir = nothing
    sigma = 0.01
    dz_mm = 3.0
    ks = 3.0
    eveto = 10.0
    roi_low = NaN
    roi_high = NaN
    seed = 42

    i = 1
    while i <= length(args)
        a = args[i]
        if a in ("-h", "--help")
            print_help(); exit(0)
        elseif a == "--indir"
            indir = args[i+1]; i += 2; continue
        elseif a == "--outdir"
            outdir = args[i+1]; i += 2; continue
        elseif a == "--sigma"
            sigma = parse(Float64, args[i+1]); i += 2; continue
        elseif a == "--dz"
            dz_mm = parse(Float64, args[i+1]); i += 2; continue
        elseif a == "--ks"
            ks = parse(Float64, args[i+1]); i += 2; continue
        elseif a == "--eveto"
            eveto = parse(Float64, args[i+1]); i += 2; continue
        elseif a == "--roi"
            roi_low = parse(Float64, args[i+1])
            roi_high = parse(Float64, args[i+2])
            i += 3; continue
        elseif a == "--seed"
            seed = parse(Int, args[i+1]); i += 2; continue
        else
            error("Unknown argument '$a'")
        end
        i += 1
    end

    indir === nothing && error("--indir is required")
    if outdir === nothing
        outdir = joinpath(indir, "analysis")
    end

    sigma_abs = sigma * Q_BB_KEV
    fwhm = sigma_abs * FWHM_FACTOR
    if isnan(roi_low)
        roi_low = Q_BB_KEV - fwhm / 2.0
        roi_high = Q_BB_KEV + fwhm / 2.0
    end

    (indir=indir, outdir=outdir, sigma=sigma, dz_mm=dz_mm, ks=ks,
     eveto=eveto, roi_low=roi_low, roi_high=roi_high, seed=seed)
end

# =====================================================================
# CSV reader (fast, minimal)
# =====================================================================

struct Deposit
    event_id::Int
    x::Float64
    y::Float64
    z::Float64
    energy_keV::Float64
    source::String
    interaction::String
    volume::String
end

function read_deposits(path::String)::Vector{Deposit}
    lines = readlines(path)
    n = length(lines) - 1  # skip header
    deposits = Vector{Deposit}(undef, n)
    @inbounds for i in 1:n
        parts = split(lines[i+1], ',')
        deposits[i] = Deposit(
            parse(Int, parts[1]),
            parse(Float64, parts[2]),
            parse(Float64, parts[3]),
            parse(Float64, parts[4]),
            parse(Float64, parts[5]) * 1000.0,  # MeV -> keV
            String(parts[6]),
            String(parts[7]),
            String(parts[8]),
        )
    end
    deposits
end

# =====================================================================
# Group deposits by event_id (returns ranges into sorted array)
# =====================================================================

function group_by_event(deposits::Vector{Deposit})
    # deposits are already sorted by event_id from the Julia writer
    groups = UnitRange{Int}[]
    n = length(deposits)
    n == 0 && return groups
    start = 1
    current_id = deposits[1].event_id
    for i in 2:n
        if deposits[i].event_id != current_id
            push!(groups, start:i-1)
            start = i
            current_id = deposits[i].event_id
        end
    end
    push!(groups, start:n)
    groups
end

# =====================================================================
# Per-event analysis
# =====================================================================

struct EventSummary
    event_id::Int
    n_deposits::Int
    Emax_keV::Float64
    Etot_keV::Float64
    DZ_mm::Float64
    R_max_cm::Float64
    Z_max_cm::Float64
end

struct EventClassification
    event_id::Int
    active_vetoed::Bool
    n_separated::Int
    is_ss::Bool
    E_cluster_keV::Float64
    E_smeared_keV::Float64
end

struct EtotByVolume
    event_id::Int
    etot_all::Float64
    etot_fv::Float64
    etot_active::Float64
    etot_passive::Float64
end

struct SeparatedDeposit
    energy_keV::Float64
    dz_mm::Float64
end

function analyze_event(deps::AbstractVector{Deposit}, range::UnitRange{Int},
                       sigma_rel::Float64, dz_mm::Float64, ks::Float64,
                       eveto_keV::Float64, rng::AbstractRNG)
    eid = deps[range[1]].event_id
    n = length(range)

    # --- Etot by volume ---
    etot_all = 0.0; etot_fv = 0.0; etot_active = 0.0; etot_passive = 0.0
    for i in range
        e = deps[i].energy_keV
        etot_all += e
        v = deps[i].volume
        if v == "fv"
            etot_fv += e
        elseif v == "active"
            etot_active += e
        elseif v == "passive"
            etot_passive += e
        end
    end

    # --- Event summary (all deposits) ---
    idx_max = range[1]
    emax = deps[idx_max].energy_keV
    for i in range
        if deps[i].energy_keV > emax
            emax = deps[i].energy_keV
            idx_max = i
        end
    end

    z_max = deps[idx_max].z
    x_max = deps[idx_max].x
    y_max = deps[idx_max].y
    r_max = sqrt(x_max^2 + y_max^2)

    # DZ: Emax to energy-weighted centroid of rest
    sum_ez = 0.0; sum_e = 0.0
    esat_energies = Float64[]
    for i in range
        i == idx_max && continue
        e = deps[i].energy_keV
        push!(esat_energies, e)
        sum_ez += deps[i].z * e
        sum_e += e
    end
    dz_cm = sum_e > 0.0 ? abs(z_max - sum_ez / sum_e) : 0.0

    summary = EventSummary(eid, n, emax, etot_all, dz_cm * 10.0, r_max, z_max)
    evol = EtotByVolume(eid, etot_all, etot_fv, etot_active, etot_passive)

    # --- Classification ---
    # Step 1: active veto
    active_vetoed = false
    for i in range
        if deps[i].volume == "active" && deps[i].energy_keV >= eveto_keV
            active_vetoed = true
            break
        end
    end

    if active_vetoed
        clf = EventClassification(eid, true, 0, false, 0.0, 0.0)
        return summary, clf, evol, esat_energies, SeparatedDeposit[]
    end

    # Step 2: FV deposits only
    fv_idx = Int[]
    for i in range
        deps[i].volume == "fv" && push!(fv_idx, i)
    end

    if isempty(fv_idx)
        clf = EventClassification(eid, false, 0, true, 0.0, 0.0)
        return summary, clf, evol, esat_energies, SeparatedDeposit[]
    end

    # Step 3: SS/MS on FV deposits
    fv_idx_max = fv_idx[1]
    e1 = deps[fv_idx[1]].energy_keV
    for j in fv_idx
        if deps[j].energy_keV > e1
            e1 = deps[j].energy_keV
            fv_idx_max = j
        end
    end
    z1 = deps[fv_idx_max].z
    sigma_e1 = sigma_rel * e1
    e_threshold = ks * sigma_e1
    dz_threshold_cm = dz_mm / 10.0

    sep_deps = SeparatedDeposit[]
    attached_energy = 0.0
    for j in fv_idx
        j == fv_idx_max && continue
        ei = deps[j].energy_keV
        dz_i = abs(deps[j].z - z1)
        if dz_i > dz_threshold_cm && ei > e_threshold
            push!(sep_deps, SeparatedDeposit(ei, dz_i * 10.0))
        else
            attached_energy += ei
        end
    end

    n_separated = length(sep_deps)
    is_ss = n_separated == 0
    e_cluster = e1 + attached_energy
    sigma_cluster = sigma_rel * e_cluster
    e_smeared = e_cluster + sigma_cluster * randn(rng)

    clf = EventClassification(eid, false, n_separated, is_ss, e_cluster, e_smeared)
    return summary, clf, evol, esat_energies, sep_deps
end

# =====================================================================
# Parallel processing
# =====================================================================

function process_all(deposits::Vector{Deposit}, groups::Vector{UnitRange{Int}},
                     sigma_rel::Float64, dz_mm::Float64, ks::Float64,
                     eveto_keV::Float64, seed::Int)
    n_events = length(groups)
    summaries = Vector{EventSummary}(undef, n_events)
    classifications = Vector{EventClassification}(undef, n_events)
    etot_vols = Vector{EtotByVolume}(undef, n_events)
    esat_arrays = Vector{Vector{Float64}}(undef, n_events)
    sep_arrays = Vector{Vector{SeparatedDeposit}}(undef, n_events)

    @threads for i in 1:n_events
        rng = MersenneTwister(seed + i)
        s, c, ev, esat, sep = analyze_event(deposits, groups[i],
                                             sigma_rel, dz_mm, ks, eveto_keV, rng)
        summaries[i] = s
        classifications[i] = c
        etot_vols[i] = ev
        esat_arrays[i] = esat
        sep_arrays[i] = sep
    end

    summaries, classifications, etot_vols, esat_arrays, sep_arrays
end

# =====================================================================
# Writers
# =====================================================================

function write_event_summary(path::String, summaries::Vector{EventSummary})
    open(path, "w") do io
        println(io, "event_id,n_deposits,Emax_keV,Etot_keV,DZ_mm,R_max_cm,Z_max_cm")
        for s in summaries
            @printf(io, "%d,%d,%.4f,%.4f,%.4f,%.4f,%.4f\n",
                    s.event_id, s.n_deposits, s.Emax_keV, s.Etot_keV,
                    s.DZ_mm, s.R_max_cm, s.Z_max_cm)
        end
    end
end

function write_classification(path::String, clf::Vector{EventClassification})
    open(path, "w") do io
        println(io, "event_id,active_vetoed,n_separated,is_ss,E_cluster_keV,E_smeared_keV")
        for c in clf
            @printf(io, "%d,%s,%d,%s,%.4f,%.4f\n",
                    c.event_id, c.active_vetoed, c.n_separated, c.is_ss,
                    c.E_cluster_keV, c.E_smeared_keV)
        end
    end
end

function write_etot_by_volume(path::String, evols::Vector{EtotByVolume})
    open(path, "w") do io
        println(io, "event_id,etot_all,etot_fv,etot_active,etot_passive")
        for ev in evols
            @printf(io, "%d,%.4f,%.4f,%.4f,%.4f\n",
                    ev.event_id, ev.etot_all, ev.etot_fv, ev.etot_active, ev.etot_passive)
        end
    end
end

function write_esat(path::String, esat_arrays::Vector{Vector{Float64}})
    open(path, "w") do io
        println(io, "energy_keV")
        for arr in esat_arrays
            for e in arr
                @printf(io, "%.4f\n", e)
            end
        end
    end
end

function write_separated_deposits(path::String, sep_arrays::Vector{Vector{SeparatedDeposit}})
    open(path, "w") do io
        println(io, "energy_keV,dz_mm")
        for arr in sep_arrays
            for sd in arr
                @printf(io, "%.4f,%.4f\n", sd.energy_keV, sd.dz_mm)
            end
        end
    end
end

function write_summary_json(path::String, clf::Vector{EventClassification},
                            roi_low::Float64, roi_high::Float64,
                            sigma::Float64, dz_mm::Float64, ks::Float64,
                            eveto::Float64, metadata, elapsed::Float64)
    n_total = length(clf)
    n_vetoed = count(c -> c.active_vetoed, clf)
    n_fv_only = n_total - n_vetoed
    n_ss = count(c -> !c.active_vetoed && c.is_ss, clf)
    n_ms = n_fv_only - n_ss
    n_ss_roi = count(c -> !c.active_vetoed && c.is_ss &&
                          c.E_smeared_keV >= roi_low && c.E_smeared_keV <= roi_high, clf)
    n_ss_out = n_ss - n_ss_roi

    f(n) = n_total > 0 ? n / n_total : 0.0

    sigma_abs = sigma * Q_BB_KEV
    fwhm = sigma_abs * FWHM_FACTOR

    # Background rate from metadata
    rate_gammas_per_s = get(metadata, "rate_total_gammas_per_s", 0.0)
    f_fv_mc = get(metadata, "f_fv", 0.0)
    gamma_BR = get(metadata, "gamma_BR", 0.0)
    source = get(metadata, "source", "unknown")
    isotope = get(metadata, "isotope", "unknown")
    components = get(metadata, "source_components", [])

    rate_per_year = rate_gammas_per_s * SECONDS_PER_YEAR
    n_fv_per_year = rate_per_year * f_fv_mc
    f_fv_only = n_fv_only / n_total
    f_ss_of_fvonly = n_fv_only > 0 ? n_ss / n_fv_only : 0.0
    f_roi_of_ss = n_ss > 0 ? n_ss_roi / n_ss : 0.0
    n_ss_roi_per_year = n_fv_per_year * f_fv_only * f_ss_of_fvonly * f_roi_of_ss

    result = Dict{String,Any}(
        "parameters" => Dict{String,Any}(
            "sigma_rel" => sigma,
            "sigma_abs_keV" => sigma_abs,
            "fwhm_keV" => fwhm,
            "roi_low_keV" => roi_low,
            "roi_high_keV" => roi_high,
            "roi_width_keV" => roi_high - roi_low,
            "dz_mm" => dz_mm,
            "ks" => ks,
            "eveto_keV" => eveto,
            "Qbb_keV" => Q_BB_KEV,
        ),
        "counts" => Dict{String,Any}(
            "n_total" => n_total,
            "n_active_vetoed" => n_vetoed,
            "n_fv_only" => n_fv_only,
            "n_ms" => n_ms,
            "n_ss" => n_ss,
            "n_ss_outside_roi" => n_ss_out,
            "n_ss_in_roi" => n_ss_roi,
        ),
        "fractions" => Dict{String,Any}(
            "f_active_vetoed" => f(n_vetoed),
            "f_fv_only" => f(n_fv_only),
            "f_ms" => f(n_ms),
            "f_ss" => f(n_ss),
            "f_ss_outside_roi" => f(n_ss_out),
            "f_ss_in_roi" => f(n_ss_roi),
        ),
        "background_rate" => Dict{String,Any}(
            "source" => source,
            "isotope" => isotope,
            "gamma_BR" => gamma_BR,
            "rate_gammas_per_s" => rate_gammas_per_s,
            "rate_gammas_per_year" => rate_per_year,
            "f_fv_mc" => f_fv_mc,
            "n_fv_per_year" => n_fv_per_year,
            "f_no_active_veto" => f_fv_only,
            "f_ss_of_fvonly" => f_ss_of_fvonly,
            "f_roi_of_ss" => f_roi_of_ss,
            "n_ss_roi_per_year" => n_ss_roi_per_year,
            "source_components" => components,
        ),
        "elapsed_s" => elapsed,
    )

    open(path, "w") do io
        JSON.print(io, result, 2)
        println(io)
    end
end

# =====================================================================
# Main
# =====================================================================

function main()
    cli = parse_cli(ARGS)

    @printf("Input:  %s\n", cli.indir)
    @printf("Output: %s\n", cli.outdir)
    @printf("Threads: %d\n", nthreads())
    @printf("sigma=%.4f  dz=%.1f mm  ks=%.1f  eveto=%.1f keV\n",
            cli.sigma, cli.dz_mm, cli.ks, cli.eveto)
    @printf("ROI = [%.1f, %.1f] keV\n", cli.roi_low, cli.roi_high)

    # --- Read ---
    t0 = time()
    deposits_path = joinpath(cli.indir, "fv_deposits.csv")
    metadata_path = joinpath(cli.indir, "metadata.json")

    @printf("Reading %s ...\n", deposits_path)
    deposits = read_deposits(deposits_path)
    @printf("  %d deposits\n", length(deposits))

    metadata = if isfile(metadata_path)
        open(metadata_path) do io; JSON.parse(io) end
    else
        @printf("  Warning: %s not found\n", metadata_path)
        Dict{String,Any}()
    end

    groups = group_by_event(deposits)
    n_events = length(groups)
    @printf("  %d events\n", n_events)
    @printf("  Read time: %.2f s\n", time() - t0)

    # --- Process ---
    t1 = time()
    @printf("Processing ...\n")
    summaries, classifications, etot_vols, esat_arrays, sep_arrays =
        process_all(deposits, groups, cli.sigma, cli.dz_mm, cli.ks,
                    cli.eveto, cli.seed)
    elapsed = time() - t1
    @printf("  Processing time: %.2f s (%.0f events/s)\n", elapsed, n_events / elapsed)

    # --- Print summary ---
    n_vetoed = count(c -> c.active_vetoed, classifications)
    n_fv_only = n_events - n_vetoed
    n_ss = count(c -> !c.active_vetoed && c.is_ss, classifications)
    n_ms = n_fv_only - n_ss
    n_ss_roi = count(c -> !c.active_vetoed && c.is_ss &&
                          c.E_smeared_keV >= cli.roi_low &&
                          c.E_smeared_keV <= cli.roi_high, classifications)

    @printf("\nResults:\n")
    @printf("  Total events        : %8d\n", n_events)
    @printf("  Active-vetoed       : %8d  (%.6f)\n", n_vetoed, n_vetoed/n_events)
    @printf("  FV-only             : %8d  (%.6f)\n", n_fv_only, n_fv_only/n_events)
    @printf("  MS rejected         : %8d  (%.6f)\n", n_ms, n_ms/n_events)
    @printf("  SS events           : %8d  (%.6f)\n", n_ss, n_ss/n_events)
    @printf("  SS in ROI           : %8d  (%.6f)\n", n_ss_roi, n_ss_roi/n_events)

    # --- Write ---
    mkpath(cli.outdir)
    @printf("\nWriting output to %s ...\n", cli.outdir)

    write_event_summary(joinpath(cli.outdir, "event_summary.csv"), summaries)
    write_classification(joinpath(cli.outdir, "classification.csv"), classifications)
    write_etot_by_volume(joinpath(cli.outdir, "etot_by_volume.csv"), etot_vols)
    write_esat(joinpath(cli.outdir, "esat.csv"), esat_arrays)
    write_separated_deposits(joinpath(cli.outdir, "separated_deposits.csv"), sep_arrays)
    write_summary_json(joinpath(cli.outdir, "summary.json"), classifications,
                       cli.roi_low, cli.roi_high, cli.sigma, cli.dz_mm, cli.ks,
                       cli.eveto, metadata, elapsed)

    @printf("Done. Total time: %.2f s\n", time() - t0)
end

main()
