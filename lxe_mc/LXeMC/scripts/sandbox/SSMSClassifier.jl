module SSMSClassifier

"""
    SSMSClassifier

Single-site / multi-site event classification for ⁰νββ searches in
dual-phase LXe TPCs (XENON1T/nT, LZ, DARWIN).

This implements the simplified two-step algorithm:

  1. **FV-boundary veto**: any energy deposit outside the fiducial
     volume above `E_VETO` rejects the event immediately.

  2. **In-FV SS/MS classification**: a secondary cluster is tagged
     as MS if it satisfies BOTH:
        (a) |Δz| > Δz_min   (geometric peak-splitting limit)
        (b) E₂   > κ·σ(E₁)  (statistical significance above primary's
                              resolution)

Inside the FV the `√L` dependence of σ_L is averaged into a single
representative value, valid because the FV is centred in z and the
FV span in L is modest. The energy resolution σ_cal(E) is anchored
to the XENON1T-measured 0.80% at Q_ββ and scaled as 1/√E.

References
----------
- EXO-200, Phys. Rev. C 95, 025502 (2017): D_T = 55±4 cm²/s,
  v_d = 1.705 mm/μs at 380 V/cm.
- Sorensen, NIM A 635, 41 (2011): D_L ≈ 12 cm²/s at 730 V/cm.
- Capelli, PhD thesis, U. Zürich (2020): σ/E = 0.80±0.02% at
  2.46 MeV in XENON1T SR1; D_L = 29.35±0.05 cm²/s at 81 V/cm;
  v_d = 1.335±0.002 mm/μs.
- LZ Collaboration: Δz_min ≈ 3 mm operational separation.
"""
SSMSClassifier

using Printf

# ===========================================================================
# Physical constants (LXe at typical TPC operating conditions)
# ===========================================================================

"Q-value of ¹³⁶Xe 0νββ decay [keV]"
const Q_BB = 2458.07

"Drift velocity of electrons in LXe at ~380 V/cm [mm/μs] (EXO-200)"
const V_DRIFT = 1.705

"Longitudinal diffusion coefficient [cm²/s] (Sorensen, XENON10 reanalysis)"
const D_L = 12.0

"Transverse diffusion coefficient [cm²/s] (EXO-200)"
const D_T = 55.0

# ===========================================================================
# Detector-specific parameters
# ===========================================================================

"""
    DetectorParams

Container for detector-specific parameters needed by the classifier.

Fields
------
- `L_fv_mean::Float64`   : FV-averaged drift length [m]
- `L_fv_min::Float64`    : FV upper boundary (shortest drift) [m]
- `L_fv_max::Float64`    : FV lower boundary (longest drift) [m]
- `tau_e::Float64`       : electron lifetime [μs]
- `sigma_E_rel_at_Q::Float64` : relative energy resolution at Q_ββ
- `name::String`         : detector identifier
"""
struct DetectorParams
    L_fv_mean::Float64
    L_fv_min::Float64
    L_fv_max::Float64
    tau_e::Float64
    sigma_E_rel_at_Q::Float64
    name::String
end

# Detector presets. Numbers are FV-typical, not all from the same source;
# adjust for the specific analysis being projected.
const XENON1T = DetectorParams(0.50, 0.25, 0.75,  650.0, 0.0080, "XENON1T")
const XENONnT = DetectorParams(0.75, 0.40, 1.15, 2000.0, 0.0080, "XENONnT")
const LZ      = DetectorParams(0.70, 0.30, 1.15, 1000.0, 0.0080, "LZ")
const DARWIN  = DetectorParams(1.30, 0.60, 2.00, 5000.0, 0.0080, "DARWIN")

# ===========================================================================
# Algorithm-level constants (see manual for justification)
# ===========================================================================

"""
    Significance factor for declaring a secondary cluster real.

κ=3 corresponds to 3σ above the primary's energy-resolution fluctuation.
This empirically absorbs electron-attachment fluctuations, electronic
noise, AP/PI, and single-electron emission, since σ_cal is calibrated
on real data.
"""
const KAPPA = 3.0

"""
    Minimum-Δz multiplier for peak splitting.

Δz_min = N_SIGMA_DZ × σ_L. Below ~3σ_L, double-Gaussian fits become
degenerate between (E₁,E₂) and (E₁+δ, E₂−δ) and the secondary's
energy cannot be reliably assigned.
"""
const N_SIGMA_DZ = 3.0

"""
    Veto threshold for deposits OUTSIDE the FV [keV].

A deposit in the active volume but outside the FV with E > E_VETO
kills the event. ~10 keV is conservative; published LZ/XENONnT
analyses operate closer to a few keV. Scan this value as a systematic.
"""
const E_VETO = 10.0

"""
    Absolute floor on E₂ for in-FV MS tagging [keV].

Set by electronics noise + AP/PI floor; below this no secondary
cluster can be tagged regardless of statistical significance. From
operational experience in XENON1T S2SingleScatter ML classifier.
"""
const E2_FLOOR = 5.0

"""
    Energy ROI around Q_ββ for the SS sample [keV].

Default ±2σ at 0.8% resolution. Use a wider window (±4σ) for blinded
analyses; a narrower window for sensitivity-limit calculations.
"""
const ROI_LOW  = Q_BB - 2 * 0.008 * Q_BB   # ≈ 2418 keV
const ROI_HIGH = Q_BB + 2 * 0.008 * Q_BB   # ≈ 2497 keV

# ===========================================================================
# Derived quantities
# ===========================================================================

"""
    sigma_L(L_meters)

Longitudinal-diffusion RMS spread of the electron cloud after drifting
length `L` [m]. Returns σ_L in **mm**.

σ_L = √(2 D_L L / v_d), with D_L in cm²/s, v_d in cm/s, L in cm.
The √L scaling means the FV-averaged σ_L is well-defined for FVs
centred away from the cathode/gate.
"""
@inline function sigma_L(L_m::Float64)::Float64
    # Convert v_drift from mm/μs to cm/s:
    #   1.705 mm/μs = 0.1705 cm/μs = 0.1705 × 10⁶ cm/s = 1.705 × 10⁵ cm/s
    v_cm_per_s = V_DRIFT * 1e5     # mm/μs → cm/s
    L_cm = 100.0 * L_m
    sigma_cm = sqrt(2.0 * D_L * L_cm / v_cm_per_s)
    return 10.0 * sigma_cm   # cm → mm
end

"""
    sigma_E_cal(E, det)

Calibration energy resolution σ_cal(E) at energy `E` [keV] for
detector `det`. Anchored at Q_ββ and scaled as σ ∝ √E (photoelectron
counting limit). Already includes electron attachment, electronics
noise, AP/PI, and single-electron fluctuations as measured in
calibration data.
"""
@inline function sigma_E_cal(E_keV::Float64, det::DetectorParams)::Float64
    sigma_at_Q = det.sigma_E_rel_at_Q * Q_BB    # keV
    return sigma_at_Q * sqrt(E_keV / Q_BB)
end

"""
    dz_min(det)

Minimum z-separation [mm] at which two clusters can be split into
peaks with reliable energy assignment. Set to N_SIGMA_DZ × σ_L at
the FV-averaged drift length.
"""
@inline function dz_min(det::DetectorParams)::Float64
    return N_SIGMA_DZ * sigma_L(det.L_fv_mean)
end

"""
    e2_cut(E1, det)

Minimum E₂ [keV] for a secondary cluster to be tagged as a real
deposit. The dominant term is κ·σ_cal(E₁), with the floor E2_FLOOR
kicking in only at very small E₁.

Note: it's σ(E₁), not σ(E₂), because the relevant fluctuation for
detecting a small bump on top of a big primary is the resolution of
the primary, not of the secondary itself.
"""
@inline function e2_cut(E1_keV::Float64, det::DetectorParams)::Float64
    return max(KAPPA * sigma_E_cal(E1_keV, det), E2_FLOOR)
end

# ===========================================================================
# Event-level data structure
# ===========================================================================

"""
    Deposit

A single energy deposit from the underlying interaction MC.

Fields
------
- `x, y, z`    : position [m]; z is the drift coordinate
- `E`          : energy [keV]
- `in_fv`      : precomputed FV membership (set by the geometry stage)
"""
struct Deposit
    x::Float64
    y::Float64
    z::Float64
    E::Float64
    in_fv::Bool
end

"""
    EventResult

Outcome of classifying one event.

Fields
------
- `passed_veto::Bool`  : survived the FV-boundary veto
- `in_roi::Bool`       : total reconstructed energy in the [ROI_LOW, ROI_HIGH] window
- `is_ms::Bool`        : tagged as multi-site (rejected from 0νββ sample)
- `is_ss::Bool`        : tagged as single-site (passes to 0νββ sample)
- `E_total::Float64`   : total in-FV energy [keV]
- `E_primary::Float64` : energy of the dominant cluster [keV]
- `n_clusters::Int`    : number of distinct in-FV deposits
"""
struct EventResult
    passed_veto::Bool
    in_roi::Bool
    is_ms::Bool
    is_ss::Bool
    E_total::Float64
    E_primary::Float64
    n_clusters::Int
end

# ===========================================================================
# The classifier
# ===========================================================================

"""
    classify(deposits, det) -> EventResult

Apply the full SS/MS classification to a single MC event.

Steps
-----
1. **FV-boundary veto**: any out-of-FV deposit with E > E_VETO → reject.

2. **Energy ROI**: total in-FV energy must lie in [ROI_LOW, ROI_HIGH].

3. **MS tagging**: find the dominant in-FV cluster, then check whether
   any other in-FV cluster satisfies BOTH the geometric (Δz > dz_min)
   AND the energetic (E > κ·σ(E₁)) criteria.

The function returns the full `EventResult` so the caller can compute
SS/MS efficiencies, leakage probabilities, and depth-dependent
distributions in a single MC pass.
"""
function classify(deposits::AbstractVector{Deposit},
                  det::DetectorParams)::EventResult

    # --------------------------------------------------------------
    # Step 1: FV-boundary veto
    # --------------------------------------------------------------
    # A deposit outside the FV with E > E_VETO is fatal. This is the
    # cheap check that lets us reject ~most of 10⁸ events without
    # any per-event diffusion modelling.
    for d in deposits
        if !d.in_fv && d.E > E_VETO
            return EventResult(false, false, false, false, 0.0, 0.0, 0)
        end
    end

    # --------------------------------------------------------------
    # Collect in-FV deposits (the only ones that matter from here)
    # --------------------------------------------------------------
    fv_deposits = filter(d -> d.in_fv, deposits)
    n_clusters = length(fv_deposits)

    # No in-FV energy → nothing to classify
    if n_clusters == 0
        return EventResult(true, false, false, false, 0.0, 0.0, 0)
    end

    E_total = sum(d.E for d in fv_deposits)

    # --------------------------------------------------------------
    # Step 2: Energy ROI window
    # --------------------------------------------------------------
    in_roi = ROI_LOW ≤ E_total ≤ ROI_HIGH

    # Find the dominant cluster (anchors σ_cal for the threshold)
    E_primary = 0.0
    idx_primary = 0
    for (i, d) in enumerate(fv_deposits)
        if d.E > E_primary
            E_primary = d.E
            idx_primary = i
        end
    end

    # If outside ROI, no need to do MS tagging — the event is already
    # excluded from the 0νββ count. We still flag SS/MS for diagnostic
    # plots (e.g. background-component identification in sidebands).
    if !in_roi
        is_ms = _check_ms(fv_deposits, idx_primary, E_primary, det)
        return EventResult(true, false, is_ms, !is_ms,
                           E_total, E_primary, n_clusters)
    end

    # --------------------------------------------------------------
    # Step 3: MS tagging (in-ROI events)
    # --------------------------------------------------------------
    is_ms = _check_ms(fv_deposits, idx_primary, E_primary, det)

    return EventResult(true, true, is_ms, !is_ms,
                       E_total, E_primary, n_clusters)
end

"""
    _check_ms(fv_deposits, idx_primary, E_primary, det) -> Bool

Internal helper: returns `true` if any non-primary in-FV deposit
satisfies both the geometric (Δz) and energetic (κσ) criteria.
"""
@inline function _check_ms(fv_deposits::AbstractVector{Deposit},
                           idx_primary::Int,
                           E_primary::Float64,
                           det::DetectorParams)::Bool

    # Threshold values (constant within an event)
    Δz_threshold = dz_min(det)            # mm
    E_threshold  = e2_cut(E_primary, det) # keV

    z_primary = fv_deposits[idx_primary].z   # m

    @inbounds for (i, d) in enumerate(fv_deposits)
        i == idx_primary && continue
        # Δz in mm (positions stored in m → convert)
        Δz_mm = abs(d.z - z_primary) * 1000.0
        if Δz_mm > Δz_threshold && d.E > E_threshold
            return true
        end
    end
    return false
end

# ===========================================================================
# Convenience: batch processing for MC sweeps
# ===========================================================================

"""
    classify_batch(events, det) -> Vector{EventResult}

Apply `classify` to a vector of events. Each event is a vector of
`Deposit`s. Suitable for the inner loop of a 10⁸-event MC. For
parallelism, wrap this in `Threads.@threads` over event chunks.
"""
function classify_batch(events::AbstractVector{<:AbstractVector{Deposit}},
                        det::DetectorParams)::Vector{EventResult}
    results = Vector{EventResult}(undef, length(events))
    @inbounds for i in eachindex(events)
        results[i] = classify(events[i], det)
    end
    return results
end

# ===========================================================================
# Diagnostic / reporting
# ===========================================================================

"""
    summary(results) -> NamedTuple

Compute the standard summary statistics from a batch of event results:
veto rejection, ROI selection, SS/MS split. Useful for ROC scans.
"""
function summary(results::AbstractVector{EventResult})
    n_total       = length(results)
    n_passed_veto = count(r -> r.passed_veto, results)
    n_in_roi      = count(r -> r.in_roi, results)
    n_ss_in_roi   = count(r -> r.is_ss && r.in_roi, results)
    n_ms_in_roi   = count(r -> r.is_ms && r.in_roi, results)

    return (
        n_total       = n_total,
        n_passed_veto = n_passed_veto,
        n_in_roi      = n_in_roi,
        n_ss_in_roi   = n_ss_in_roi,
        n_ms_in_roi   = n_ms_in_roi,
        veto_eff      = 1 - n_passed_veto / n_total,
        ms_tag_eff    = n_ms_in_roi / max(1, n_in_roi),
        ss_leakage    = n_ss_in_roi / max(1, n_in_roi),
    )
end

"""
    print_config(det)

Print the algorithm-level numbers for a given detector. Useful for
sanity-checking before launching a long MC run.
"""
function print_config(det::DetectorParams)
    println("─"^60)
    println("Detector: $(det.name)")
    println("─"^60)
    @printf("  ⟨L⟩_FV         = %.2f m\n",  det.L_fv_mean)
    @printf("  L_FV range     = [%.2f, %.2f] m\n",
            det.L_fv_min, det.L_fv_max)
    @printf("  τ_e            = %.0f μs\n", det.tau_e)
    @printf("  σ_E/E at Q_ββ  = %.2f%%\n",  100 * det.sigma_E_rel_at_Q)
    println()
    @printf("  σ_L(⟨L⟩_FV)    = %.2f mm\n", sigma_L(det.L_fv_mean))
    @printf("  Δz_min         = %.2f mm  (= %.0f σ_L)\n",
            dz_min(det), N_SIGMA_DZ)
    @printf("  σ_cal(Q_ββ)    = %.1f keV\n", sigma_E_cal(Q_BB, det))
    @printf("  E_cut at Q_ββ  = %.1f keV (κ = %.0f)\n",
            e2_cut(Q_BB, det), KAPPA)
    @printf("  E_veto (out FV) = %.1f keV\n", E_VETO)
    @printf("  E2_floor        = %.1f keV\n", E2_FLOOR)
    @printf("  ROI            = [%.1f, %.1f] keV (±2σ)\n",
            ROI_LOW, ROI_HIGH)
    println("─"^60)
end

end  # module
