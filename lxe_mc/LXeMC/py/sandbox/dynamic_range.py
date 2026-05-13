"""
Dynamic range analysis for the LZ TPC at bb0nu energies.

Computes the S2 signal size at Q_bb and at lower energies, the ADC
dynamic range constraints on the minimum detectable energy deposit,
and the statistical floor for resolving a small S2 vertex next to
a Q_bb event.

All parameters from arXiv:1912.04248v2 (LZ bb0nu paper) and the
LZ TDR (arXiv:1703.09144).

Usage:
    python py/dynamic_range.py
    python py/dynamic_range.py --peak-factor 5
    python py/dynamic_range.py --drift-field 310
"""

import argparse
import math
import sys


# =====================================================================
# LZ detector parameters (from bb0nu paper + TDR)
# =====================================================================

W_EV = 13.7              # W-value of LXe [eV/quantum]
G2_TOTAL = 79.0           # S2 gain: detected photons per extracted e-
G2_BOTTOM = 27.0          # S2 gain in bottom PMT array only
EXTRACTION_EFF = 0.95     # electron extraction efficiency at liquid surface
E_QBB_KEV = 2458.0        # Xe-136 Q-value [keV]
N_PMT_TOP = 253           # top PMT array
N_PMT_BOT = 241           # bottom PMT array
S2_FWHM_US = 1.0          # S2 pulse FWHM for shallow events [us]
DRIFT_V_MM_US = 2.0       # maximum drift velocity [mm/us]
DZ_CUT_MM = 3.0           # single-scatter vertex separation cut [mm]


def ionization_fraction(E_keV: float) -> float:
    """Approximate ionization fraction for electron recoils in LXe at 310 V/cm.

    At high energy (~MeV), recombination is significant and f_ion ~ 0.65.
    At low energy (<100 keV), recombination is weaker and f_ion ~ 0.80.
    This is a simplified two-regime model; NEST gives a smooth curve.
    """
    if E_keV < 100:
        return 0.80
    return 0.65


def s2_signal(E_keV: float) -> dict:
    """Compute S2 signal parameters for an energy deposit E_keV."""
    f_ion = ionization_fraction(E_keV)
    n_quanta = E_keV * 1000.0 / W_EV
    n_electrons = n_quanta * f_ion
    s2_total = n_electrons * EXTRACTION_EFF * G2_TOTAL
    s2_bottom = n_electrons * EXTRACTION_EFF * G2_BOTTOM
    s2_top = s2_total - s2_bottom
    return dict(
        E_keV=E_keV,
        f_ion=f_ion,
        n_quanta=n_quanta,
        n_electrons=n_electrons,
        s2_total=s2_total,
        s2_bottom=s2_bottom,
        s2_top=s2_top,
        s2_per_pmt_bot=s2_bottom / N_PMT_BOT,
        s2_per_pmt_top=s2_top / N_PMT_TOP,
    )


def print_header(title: str) -> None:
    print()
    print("=" * 70)
    print(title)
    print("=" * 70)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="LZ dynamic range analysis at bb0nu energies")
    parser.add_argument("--peak-factor", type=float, default=8.0,
                        help="Peak PMT signal / average (geometric concentration)")
    parser.add_argument("--adc-bits", type=int, nargs="+", default=[12, 14, 16],
                        help="ADC bit depths to analyze")
    parser.add_argument("--adc-usage", type=float, default=0.8,
                        help="Fraction of ADC range used by peak PMT at Qbb")
    parser.add_argument("--min-adc-counts", type=int, default=5,
                        help="Minimum ADC counts for a resolvable S2 pulse")
    args = parser.parse_args()

    peak_factor = args.peak_factor

    # --- S2 at Q_bb ---
    print_header("S2 SIGNAL AT Q_bb (2458 keV)")

    qbb = s2_signal(E_QBB_KEV)
    print(f"  W-value          = {W_EV} eV")
    print(f"  g2 (total)       = {G2_TOTAL} phd/e-")
    print(f"  g2 (bottom)      = {G2_BOTTOM} phd/e-")
    print(f"  Extraction eff.  = {EXTRACTION_EFF}")
    print()
    print(f"  N_quanta         = {qbb['n_quanta']:.0f}")
    print(f"  N_electrons      = {qbb['n_electrons']:.0f}  (f_ion = {qbb['f_ion']})")
    print(f"  S2 total         = {qbb['s2_total']:.0f} phd")
    print(f"  S2 bottom        = {qbb['s2_bottom']:.0f} phd  ({N_PMT_BOT} PMTs)")
    print(f"  S2 top           = {qbb['s2_top']:.0f} phd  ({N_PMT_TOP} PMTs)")
    print(f"  S2/PMT (bot avg) = {qbb['s2_per_pmt_bot']:.0f} phd")
    print(f"  S2/PMT (top avg) = {qbb['s2_per_pmt_top']:.0f} phd")
    print(f"  Peak top PMT (~{peak_factor:.0f}x avg): "
          f"{qbb['s2_per_pmt_top'] * peak_factor:.0f} phd")

    # --- ADC dynamic range ---
    print_header("ADC DYNAMIC RANGE CONSTRAINTS")
    print(f"  Peak PMT concentration factor: {peak_factor:.0f}x average")
    print(f"  ADC usage at Qbb: {args.adc_usage:.0%} of full range")
    print(f"  Min resolvable: {args.min_adc_counts} ADC counts")

    peak_phd = qbb['s2_per_pmt_top'] * peak_factor

    for n_bits in args.adc_bits:
        adc_max = 2**n_bits - 1
        phd_per_count = peak_phd / (args.adc_usage * adc_max)
        min_phd = args.min_adc_counts * phd_per_count
        min_electrons = min_phd / (G2_TOTAL * EXTRACTION_EFF)
        f_ion_low = ionization_fraction(1.0)
        min_E_keV = min_electrons * W_EV / 1000.0 / f_ion_low

        print(f"\n  {n_bits}-bit ADC ({adc_max} counts):")
        print(f"    phd/count           = {phd_per_count:.2f}")
        print(f"    min S2 (~{args.min_adc_counts} counts)  = {min_phd:.0f} phd")
        print(f"    -> min electrons    = {min_electrons:.0f}")
        print(f"    -> min E deposit    = {min_E_keV:.1f} keV")

    # --- S2 at various deposits ---
    print_header("S2 SIGNAL vs ENERGY DEPOSIT")
    print(f"  {'E (keV)':>10s}  {'f_ion':>6s}  {'N_e':>8s}  "
          f"{'S2 total':>10s}  {'S2/PMT bot':>12s}  {'ratio to Qbb':>14s}")
    for E_keV in [0.5, 1, 2, 5, 10, 20, 50, 100, 500, 1000]:
        sig = s2_signal(E_keV)
        ratio = sig['s2_total'] / qbb['s2_total']
        print(f"  {E_keV:>10.1f}  {sig['f_ion']:>6.2f}  {sig['n_electrons']:>8.0f}  "
              f"{sig['s2_total']:>10.0f}  {sig['s2_per_pmt_bot']:>12.1f}  "
              f"{ratio:>14.1e}")

    # --- Vertex resolution ---
    print_header("VERTEX RESOLUTION AT 3 mm")

    dt_us = DZ_CUT_MM / DRIFT_V_MM_US
    resolved = dt_us > S2_FWHM_US
    print(f"  S2 FWHM:           {S2_FWHM_US} us")
    print(f"  Drift velocity:    <= {DRIFT_V_MM_US} mm/us")
    print(f"  {DZ_CUT_MM:.0f} mm separation -> {dt_us:.1f} us between S2 peaks")
    print(f"  Peaks are {'resolved' if resolved else 'overlapping'}"
          f"  (separation/FWHM = {dt_us / S2_FWHM_US:.1f})")

    print(f"\n  Statistical floor for resolving small S2 next to Qbb S2:")
    print(f"  The small peak sits on the tail of the large one.")
    print(f"  Need S2_small >> sqrt(S2_large) for statistical resolution.")
    sqrt_s2 = math.sqrt(qbb['s2_total'])
    n_e_floor = sqrt_s2 / (G2_TOTAL * EXTRACTION_EFF)
    f_ion_low = ionization_fraction(1.0)
    E_floor_keV = n_e_floor * W_EV / 1000.0 / f_ion_low
    print(f"  sqrt(S2 at Qbb)  = {sqrt_s2:.0f} phd")
    print(f"  -> {n_e_floor:.0f} electrons")
    print(f"  -> E_floor ~ {E_floor_keV:.1f} keV  (1-sigma statistical limit)")
    print(f"  -> 3-sigma floor ~ {3 * E_floor_keV:.1f} keV")
    print(f"  -> 5-sigma floor ~ {5 * E_floor_keV:.1f} keV")

    # --- Summary ---
    print_header("SUMMARY")
    print(f"  The ADC dynamic range is NOT the limiting factor:")
    print(f"  even a 12-bit ADC resolves sub-keV deposits.")
    print()
    print(f"  The statistical limit IS relevant: to resolve a small")
    print(f"  S2 vertex on the tail of a {E_QBB_KEV:.0f} keV S2 pulse,")
    print(f"  the deposit must exceed ~{3 * E_floor_keV:.0f} keV (3-sigma).")
    print()
    print(f"  Current LXeMC veto_TPC = 10 keV.")
    print(f"  Statistical floor (3-sigma) ~ {3 * E_floor_keV:.0f} keV.")
    print(f"  These are consistent: 10 keV > {E_floor_keV:.1f} keV (1-sigma)")
    print(f"  but below the 3-sigma level.")


if __name__ == "__main__":
    main()
