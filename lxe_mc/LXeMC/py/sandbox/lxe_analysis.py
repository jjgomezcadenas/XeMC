"""
ssms_classifier.py
==================

Single-site / multi-site event classification for 0νββ searches in
dual-phase LXe TPCs (XENON1T/nT, LZ, DARWIN).

This implements the simplified two-step algorithm:

  1. **FV-boundary veto**: any energy deposit outside the fiducial
     volume above ``E_VETO`` rejects the event immediately.

  2. **In-FV SS/MS classification**: a secondary cluster is tagged
     as MS if it satisfies BOTH:
        (a) |Δz| > Δz_min   (geometric peak-splitting limit)
        (b) E₂   > κ·σ(E₁)  (statistical significance above primary's
                              resolution)

Inside the FV the √L dependence of σ_L is averaged into a single
representative value, valid because the FV is centred in z and the
FV span in L is modest. The energy resolution σ_cal(E) is anchored
to the XENON1T-measured 0.80% at Q_ββ and scaled as σ ∝ √E.

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

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Iterable, List, Sequence


# ===========================================================================
# Physical constants (LXe at typical TPC operating conditions)
# ===========================================================================

Q_BB: float = 2458.07
"""Q-value of ¹³⁶Xe 0νββ decay [keV]."""

V_DRIFT: float = 1.705
"""Drift velocity of electrons in LXe at ~380 V/cm [mm/μs] (EXO-200)."""

D_L: float = 12.0
"""Longitudinal diffusion coefficient [cm²/s] (Sorensen, XENON10 reanalysis)."""

D_T: float = 55.0
"""Transverse diffusion coefficient [cm²/s] (EXO-200)."""


# ===========================================================================
# Detector-specific parameters
# ===========================================================================

@dataclass(frozen=True)
class DetectorParams:
    """
    Container for detector-specific parameters needed by the classifier.

    Parameters
    ----------
    L_fv_mean : float
        FV-averaged drift length [m].
    L_fv_min : float
        FV upper boundary (shortest drift) [m].
    L_fv_max : float
        FV lower boundary (longest drift) [m].
    tau_e : float
        Electron lifetime [μs].
    sigma_E_rel_at_Q : float
        Relative energy resolution at Q_ββ.
    name : str
        Detector identifier.
    """
    L_fv_mean: float
    L_fv_min: float
    L_fv_max: float
    tau_e: float
    sigma_E_rel_at_Q: float
    name: str


# Detector presets. Numbers are FV-typical, not all from the same source;
# adjust for the specific analysis being projected.
XENON1T = DetectorParams(0.50, 0.25, 0.75,  650.0, 0.0080, "XENON1T")
XENONnT = DetectorParams(0.75, 0.40, 1.15, 2000.0, 0.0080, "XENONnT")
LZ      = DetectorParams(0.70, 0.30, 1.15, 1000.0, 0.0080, "LZ")
DARWIN  = DetectorParams(1.30, 0.60, 2.00, 5000.0, 0.0080, "DARWIN")


# ===========================================================================
# Algorithm-level constants (see manual for justification)
# ===========================================================================

KAPPA: float = 3.0
"""
Significance factor for declaring a secondary cluster real.

κ=3 corresponds to 3σ above the primary's energy-resolution fluctuation.
This empirically absorbs electron-attachment fluctuations, electronic
noise, AP/PI, and single-electron emission, since σ_cal is calibrated
on real data.
"""

N_SIGMA_DZ: float = 3.0
"""
Minimum-Δz multiplier for peak splitting.

Δz_min = N_SIGMA_DZ × σ_L. Below ~3σ_L, double-Gaussian fits become
degenerate between (E₁, E₂) and (E₁ + δ, E₂ − δ) and the secondary's
energy cannot be reliably assigned.
"""

E_VETO: float = 10.0
"""
Veto threshold for deposits OUTSIDE the FV [keV].

A deposit in the active volume but outside the FV with E > E_VETO
kills the event. ~10 keV is conservative; published LZ/XENONnT
analyses operate closer to a few keV. Scan this value as a systematic.
"""

E2_FLOOR: float = 5.0
"""
Absolute floor on E₂ for in-FV MS tagging [keV].

Set by electronics noise + AP/PI floor; below this no secondary
cluster can be tagged regardless of statistical significance. From
operational experience in XENON1T S2SingleScatter ML classifier.
"""

ROI_LOW: float = Q_BB - 2 * 0.008 * Q_BB   # ≈ 2418 keV
ROI_HIGH: float = Q_BB + 2 * 0.008 * Q_BB  # ≈ 2497 keV
"""Energy ROI around Q_ββ for the SS sample [keV]. Default ±2σ at 0.8%."""


# ===========================================================================
# Derived quantities
# ===========================================================================

def sigma_L(L_m: float) -> float:
    """
    Longitudinal-diffusion RMS spread of the electron cloud.

    Parameters
    ----------
    L_m : float
        Drift length [m].

    Returns
    -------
    float
        σ_L in **mm**.

    Notes
    -----
    σ_L = √(2 D_L L / v_d). The √L scaling means the FV-averaged σ_L
    is well-defined for FVs centred away from the cathode/gate.
    """
    # Convert v_drift from mm/μs to cm/s:
    #   1.705 mm/μs = 0.1705 cm/μs = 1.705 × 10⁵ cm/s
    v_cm_per_s = V_DRIFT * 1e5      # mm/μs → cm/s
    L_cm = 100.0 * L_m
    sigma_cm = math.sqrt(2.0 * D_L * L_cm / v_cm_per_s)
    return 10.0 * sigma_cm           # cm → mm


def sigma_E_cal(E_keV: float, det: DetectorParams) -> float:
    """
    Calibration energy resolution at energy ``E_keV``.

    Anchored at Q_ββ and scaled as σ ∝ √E (photoelectron counting limit).
    Already includes electron attachment, electronics noise, AP/PI, and
    single-electron fluctuations as measured in calibration data.

    Returns σ_cal in keV.
    """
    sigma_at_Q = det.sigma_E_rel_at_Q * Q_BB
    return sigma_at_Q * math.sqrt(E_keV / Q_BB)


def dz_min(det: DetectorParams) -> float:
    """
    Minimum z-separation [mm] at which two clusters can be split into
    peaks with reliable energy assignment.

    Set to N_SIGMA_DZ × σ_L at the FV-averaged drift length.
    """
    return N_SIGMA_DZ * sigma_L(det.L_fv_mean)


def e2_cut(E1_keV: float, det: DetectorParams) -> float:
    """
    Minimum E₂ [keV] for a secondary cluster to be tagged as a real
    deposit. The dominant term is κ·σ_cal(E₁), with E2_FLOOR kicking in
    only at very small E₁.

    Note
    ----
    It's σ(E₁), not σ(E₂), because the relevant fluctuation for
    detecting a small bump on top of a big primary is the resolution
    of the primary, not of the secondary itself.
    """
    return max(KAPPA * sigma_E_cal(E1_keV, det), E2_FLOOR)


# ===========================================================================
# Event-level data structures
# ===========================================================================

@dataclass(frozen=True)
class Deposit:
    """
    A single energy deposit from the underlying interaction MC.

    Attributes
    ----------
    x, y, z : float
        Position [m]; z is the drift coordinate.
    E : float
        Energy [keV].
    in_fv : bool
        Precomputed FV membership (set by the geometry stage).
    """
    x: float
    y: float
    z: float
    E: float
    in_fv: bool


@dataclass(frozen=True)
class EventResult:
    """
    Outcome of classifying one event.

    Attributes
    ----------
    passed_veto : bool
        Survived the FV-boundary veto.
    in_roi : bool
        Total reconstructed energy in [ROI_LOW, ROI_HIGH].
    is_ms : bool
        Tagged as multi-site (rejected from 0νββ sample).
    is_ss : bool
        Tagged as single-site (passes to 0νββ sample).
    E_total : float
        Total in-FV energy [keV].
    E_primary : float
        Energy of the dominant cluster [keV].
    n_clusters : int
        Number of distinct in-FV deposits.
    """
    passed_veto: bool
    in_roi: bool
    is_ms: bool
    is_ss: bool
    E_total: float
    E_primary: float
    n_clusters: int


# ===========================================================================
# The classifier
# ===========================================================================

def classify(deposits: Sequence[Deposit],
             det: DetectorParams) -> EventResult:
    """
    Apply the full SS/MS classification to a single MC event.

    Steps
    -----
    1. **FV-boundary veto**: any out-of-FV deposit with E > E_VETO →
       reject.

    2. **Energy ROI**: total in-FV energy must lie in [ROI_LOW, ROI_HIGH].

    3. **MS tagging**: find the dominant in-FV cluster, then check
       whether any other in-FV cluster satisfies BOTH the geometric
       (Δz > dz_min) AND the energetic (E > κ·σ(E₁)) criteria.

    Returns the full :class:`EventResult` so the caller can compute
    SS/MS efficiencies, leakage probabilities, and depth-dependent
    distributions in a single MC pass.
    """
    # --------------------------------------------------------------
    # Step 1: FV-boundary veto
    # --------------------------------------------------------------
    # A deposit outside the FV with E > E_VETO is fatal. This is the
    # cheap check that lets us reject ~most of 10⁸ events without any
    # per-event diffusion modelling.
    for d in deposits:
        if (not d.in_fv) and d.E > E_VETO:
            return EventResult(
                passed_veto=False, in_roi=False,
                is_ms=False, is_ss=False,
                E_total=0.0, E_primary=0.0, n_clusters=0,
            )

    # --------------------------------------------------------------
    # Collect in-FV deposits (the only ones that matter from here on)
    # --------------------------------------------------------------
    fv_deposits = [d for d in deposits if d.in_fv]
    n_clusters = len(fv_deposits)

    # No in-FV energy → nothing to classify
    if n_clusters == 0:
        return EventResult(
            passed_veto=True, in_roi=False,
            is_ms=False, is_ss=False,
            E_total=0.0, E_primary=0.0, n_clusters=0,
        )

    E_total = sum(d.E for d in fv_deposits)

    # --------------------------------------------------------------
    # Step 2: Energy ROI window
    # --------------------------------------------------------------
    in_roi = ROI_LOW <= E_total <= ROI_HIGH

    # Find the dominant cluster (anchors σ_cal for the threshold)
    primary = max(fv_deposits, key=lambda d: d.E)
    E_primary = primary.E

    # --------------------------------------------------------------
    # Step 3: MS tagging
    # --------------------------------------------------------------
    # We compute is_ms regardless of in_roi: the SS/MS split outside the
    # ROI feeds background-component identification in the sidebands
    # (e.g., distinguishing Bi-214 shell from material backgrounds), as
    # done in the XENON1T 0νββ analysis.
    is_ms = _check_ms(fv_deposits, primary, det)

    return EventResult(
        passed_veto=True, in_roi=in_roi,
        is_ms=is_ms, is_ss=(not is_ms),
        E_total=E_total, E_primary=E_primary,
        n_clusters=n_clusters,
    )


def _check_ms(fv_deposits: Sequence[Deposit],
              primary: Deposit,
              det: DetectorParams) -> bool:
    """
    Return True if any non-primary in-FV deposit satisfies both the
    geometric (Δz) and energetic (κσ) criteria.
    """
    # Threshold values (constant within an event)
    dz_threshold_mm = dz_min(det)            # mm
    E_threshold_keV = e2_cut(primary.E, det) # keV

    z_primary_m = primary.z

    for d in fv_deposits:
        if d is primary:
            continue
        # Δz in mm (positions stored in m → convert)
        dz_mm = abs(d.z - z_primary_m) * 1000.0
        if dz_mm > dz_threshold_mm and d.E > E_threshold_keV:
            return True
    return False


# ===========================================================================
# Convenience: batch processing for MC sweeps
# ===========================================================================

def classify_batch(events: Iterable[Sequence[Deposit]],
                   det: DetectorParams) -> List[EventResult]:
    """
    Apply :func:`classify` to a sequence of events.

    Each event is itself a sequence of :class:`Deposit`s. Suitable
    for the inner loop of a 10⁸-event MC. For parallelism, split
    ``events`` into chunks and dispatch via ``multiprocessing.Pool``
    or use NumPy-vectorised variants for the hot path.
    """
    return [classify(ev, det) for ev in events]


# ===========================================================================
# Diagnostic / reporting
# ===========================================================================

@dataclass(frozen=True)
class Summary:
    """Standard summary statistics from a batch of events."""
    n_total: int
    n_passed_veto: int
    n_in_roi: int
    n_ss_in_roi: int
    n_ms_in_roi: int
    veto_eff: float
    ms_tag_eff: float
    ss_leakage: float


def summarize(results: Sequence[EventResult]) -> Summary:
    """
    Compute summary statistics from a batch of event results:
    veto rejection, ROI selection, SS/MS split. Useful for ROC scans.
    """
    n_total       = len(results)
    n_passed_veto = sum(1 for r in results if r.passed_veto)
    n_in_roi      = sum(1 for r in results if r.in_roi)
    n_ss_in_roi   = sum(1 for r in results if r.is_ss and r.in_roi)
    n_ms_in_roi   = sum(1 for r in results if r.is_ms and r.in_roi)

    return Summary(
        n_total=n_total,
        n_passed_veto=n_passed_veto,
        n_in_roi=n_in_roi,
        n_ss_in_roi=n_ss_in_roi,
        n_ms_in_roi=n_ms_in_roi,
        veto_eff=1 - n_passed_veto / max(1, n_total),
        ms_tag_eff=n_ms_in_roi / max(1, n_in_roi),
        ss_leakage=n_ss_in_roi / max(1, n_in_roi),
    )


def print_config(det: DetectorParams) -> None:
    """
    Print the algorithm-level numbers for a given detector. Useful
    for sanity-checking before launching a long MC run.
    """
    print("─" * 60)
    print(f"Detector: {det.name}")
    print("─" * 60)
    print(f"  ⟨L⟩_FV          = {det.L_fv_mean:.2f} m")
    print(f"  L_FV range      = [{det.L_fv_min:.2f}, {det.L_fv_max:.2f}] m")
    print(f"  τ_e             = {det.tau_e:.0f} μs")
    print(f"  σ_E/E at Q_ββ   = {100 * det.sigma_E_rel_at_Q:.2f}%")
    print()
    print(f"  σ_L(⟨L⟩_FV)     = {sigma_L(det.L_fv_mean):.2f} mm")
    print(f"  Δz_min          = {dz_min(det):.2f} mm  "
          f"(= {N_SIGMA_DZ:.0f} σ_L)")
    print(f"  σ_cal(Q_ββ)     = {sigma_E_cal(Q_BB, det):.1f} keV")
    print(f"  E_cut at Q_ββ   = {e2_cut(Q_BB, det):.1f} keV "
          f"(κ = {KAPPA:.0f})")
    print(f"  E_veto (out FV) = {E_VETO:.1f} keV")
    print(f"  E2_floor        = {E2_FLOOR:.1f} keV")
    print(f"  ROI             = [{ROI_LOW:.1f}, {ROI_HIGH:.1f}] keV (±2σ)")
    print("─" * 60)
