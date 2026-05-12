"""Plot analysis results from analyze_backgrounds.jl output.

Reads pre-computed CSV/JSON from <indir>/analysis/ and produces
the same panel plots as the original analyze_backgrounds.py.

Usage:
    python py/plot_analysis.py -i results/bfv/cryostat/barrel/Bi214/analysis
    python py/plot_analysis.py -i results/bfv/cryostat/barrel/Bi214/analysis --display
"""

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

Q_BB_KEV = 2458.07
E_GAMMA_BI214_KEV = 2447.86


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot analysis results from Julia analyze_backgrounds.jl output."
    )
    parser.add_argument(
        "-i", "--input", required=True, type=Path,
        help="Analysis directory (containing CSV/JSON from analyze_backgrounds.jl)"
    )
    parser.add_argument(
        "--display", action="store_true", help="Show plots interactively"
    )
    return parser.parse_args()


def plot_event_summary(summary: pd.DataFrame, esat: pd.DataFrame,
                       outdir: Path, display: bool) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))

    ax = axes[0, 0]
    ax.hist(summary["n_deposits"].values, bins=50, edgecolor="black", linewidth=0.3)
    ax.set_xlabel("Number of deposits")
    ax.set_ylabel("Counts")
    ax.set_title("Deposits per event")

    ax = axes[0, 1]
    ax.hist(summary["Emax_keV"].values, bins=100, edgecolor="black", linewidth=0.3)
    ax.set_xlabel("Emax (keV)")
    ax.set_ylabel("Counts")
    ax.set_title("Emax: largest single deposit")

    ax = axes[0, 2]
    ax.hist(esat["energy_keV"].values, bins=100, edgecolor="black", linewidth=0.3)
    ax.set_xlabel("Esat (keV)")
    ax.set_ylabel("Counts")
    ax.set_title("Esat: each non-Emax deposit")

    ax = axes[1, 0]
    ax.hist(summary["Etot_keV"].values, bins=100, edgecolor="black", linewidth=0.3)
    ax.set_xlabel("Etot (keV)")
    ax.set_ylabel("Counts")
    ax.set_title("Etot: total energy per event")

    ax = axes[1, 1]
    ax.hist(summary["DZ_mm"].values, bins=100, edgecolor="black", linewidth=0.3)
    ax.set_yscale("log")
    ax.set_xlabel("DZ (mm)")
    ax.set_ylabel("Counts")
    ax.set_title("DZ: Emax to satellite centroid")

    ax = axes[1, 2]
    ax.scatter(summary["Z_max_cm"].values, summary["R_max_cm"].values, s=1, alpha=0.3)
    ax.set_xlabel("Z (cm)")
    ax.set_ylabel("R (cm)")
    ax.set_title("Emax deposit: R vs Z")

    fig.tight_layout()
    outfile = outdir / "event_summary.png"
    fig.savefig(outfile, dpi=150)
    print(f"Saved {outfile}")
    if display:
        plt.show()
    plt.close(fig)


def plot_energy_by_volume(evol: pd.DataFrame, outdir: Path, display: bool) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    ax = axes[0, 0]
    ax.hist(evol["etot_all"].values, bins=100, edgecolor="black", linewidth=0.3)
    ax.axvline(E_GAMMA_BI214_KEV, color="green", ls="--", lw=0.8,
               label=f"Bi214={E_GAMMA_BI214_KEV:.0f}")
    ax.legend()
    ax.set_xlabel("Energy (keV)")
    ax.set_ylabel("Counts")
    ax.set_title("Etot: all deposits")

    ax = axes[0, 1]
    fv = evol["etot_fv"].values
    fv = fv[fv > 0]
    ax.hist(fv, bins=100, edgecolor="black", linewidth=0.3)
    ax.axvline(E_GAMMA_BI214_KEV, color="green", ls="--", lw=0.8,
               label=f"Bi214={E_GAMMA_BI214_KEV:.0f}")
    ax.legend()
    ax.set_xlabel("Energy (keV)")
    ax.set_ylabel("Counts")
    ax.set_title("Etot: FV deposits only")

    ax = axes[1, 0]
    passive = evol["etot_passive"].values
    passive = passive[passive > 0]
    if len(passive) > 0:
        ax.hist(passive, bins=100, edgecolor="black", linewidth=0.3)
    ax.set_xlabel("Energy (keV)")
    ax.set_ylabel("Counts")
    ax.set_title("Etot: passive deposits only")

    ax = axes[1, 1]
    active = evol["etot_active"].values
    active = active[active > 0]
    if len(active) > 0:
        ax.hist(active, bins=100, edgecolor="black", linewidth=0.3)
    ax.set_xlabel("Energy (keV)")
    ax.set_ylabel("Counts")
    ax.set_title("Etot: active deposits only")

    fig.tight_layout()
    outfile = outdir / "energy_by_volume.png"
    fig.savefig(outfile, dpi=150)
    print(f"Saved {outfile}")
    if display:
        plt.show()
    plt.close(fig)


def plot_ssms_analysis(clf: pd.DataFrame, sep: pd.DataFrame,
                       smry: dict, outdir: Path, display: bool) -> None:
    params = smry["parameters"]
    roi_low = params["roi_low_keV"]
    roi_high = params["roi_high_keV"]
    counts = smry["counts"]

    surv = clf[~clf["active_vetoed"]]
    ss_mask = surv["is_ss"].values
    ss_ecluster = surv.loc[ss_mask, "E_cluster_keV"].values
    ss_esmeared = surv.loc[ss_mask, "E_smeared_keV"].values

    n_total = counts["n_total"]
    n_vetoed = counts["n_active_vetoed"]
    n_surv = counts["n_fv_only"]

    fig, axes = plt.subplots(2, 3, figsize=(14, 8))

    # (0,0) Veto bar chart
    ax = axes[0, 0]
    f_fv = n_surv / n_total if n_total > 0 else 0
    f_vetoed = n_vetoed / n_total if n_total > 0 else 0
    ax.bar(["FV-only", "Active-vetoed"], [f_fv, f_vetoed],
           color=["steelblue", "tomato"], edgecolor="black", linewidth=0.5)
    ax.set_ylabel("Fraction")
    ax.set_title(f"FV-only={n_surv}  Vetoed={n_vetoed}")
    ax.set_ylim(0, 1)

    # (0,1) Cluster multiplicity
    ax = axes[0, 1]
    n_sep = surv["n_separated"].values
    max_sep = n_sep.max() if len(n_sep) > 0 else 1
    bins = list(range(0, max_sep + 1))
    counts_per_bin = [int((n_sep == k).sum()) for k in bins]
    bar_colors = ["steelblue" if k == 0 else "coral" for k in bins]
    ax.bar(bins, counts_per_bin, color=bar_colors, edgecolor="black", linewidth=0.3)
    ax.set_xlabel("N separated clusters (0 = SS)")
    ax.set_ylabel("Counts")
    ax.set_title("Cluster multiplicity (FV-only)")

    # (0,2) SS vs MS
    ax = axes[0, 2]
    n_ss = counts["n_ss"]
    n_ms = counts["n_ms"]
    ax.bar(["SS", "MS"], [n_ss, n_ms], color=["steelblue", "coral"],
           edgecolor="black", linewidth=0.5)
    ax.set_ylabel("Counts")
    ax.set_title(f"SS vs MS  (SS={n_ss}, MS={n_ms})")

    # (1,0) DZ of separated deposits
    ax = axes[1, 0]
    if len(sep) > 0:
        ax.hist(sep["dz_mm"].values, bins=100, edgecolor="black", linewidth=0.3)
    ax.set_xlabel("DZ (mm)")
    ax.set_ylabel("Counts")
    ax.set_title("DZ of separated deposits")

    # (1,1) E_cluster for SS
    ax = axes[1, 1]
    if len(ss_ecluster) > 0:
        ax.hist(ss_ecluster, bins=100, edgecolor="black", linewidth=0.3)
    ax.axvline(E_GAMMA_BI214_KEV, color="green", ls="--", lw=0.8,
               label=f"Bi214={E_GAMMA_BI214_KEV:.0f}")
    ax.legend()
    ax.set_xlabel("E_cluster (keV)")
    ax.set_ylabel("Counts")
    ax.set_title("E_cluster (SS events)")

    # (1,2) E_smeared zoomed
    ax = axes[1, 2]
    if len(ss_esmeared) > 0:
        ax.hist(ss_esmeared, bins=100, range=(2200, 2700),
                edgecolor="black", linewidth=0.3)
    ax.axvspan(roi_low, roi_high, alpha=0.2, color="red", label="ROI")
    ax.axvline(Q_BB_KEV, color="red", ls="--", lw=0.8, label=f"Qbb={Q_BB_KEV:.0f}")
    ax.axvline(E_GAMMA_BI214_KEV, color="green", ls="--", lw=0.8,
               label=f"Bi214={E_GAMMA_BI214_KEV:.0f}")
    ax.set_xlim(2200, 2700)
    ax.legend()
    ax.set_xlabel("E_smeared (keV)")
    ax.set_ylabel("Counts")
    ax.set_title("E_smeared (SS) - zoom around Qbb")

    fig.tight_layout()
    outfile = outdir / "ssms_analysis.png"
    fig.savefig(outfile, dpi=150)
    print(f"Saved {outfile}")
    if display:
        plt.show()
    plt.close(fig)


def plot_summary_table(smry: dict, outdir: Path, display: bool) -> None:
    params = smry["parameters"]
    counts = smry["counts"]
    fracs = smry["fractions"]

    table_data = [
        ["sigma/E", f"{params['sigma_rel']}", ""],
        ["sigma(Qbb)", f"{params['sigma_abs_keV']:.1f} keV", ""],
        ["FWHM(Qbb)", f"{params['fwhm_keV']:.1f} keV", ""],
        ["ROI", f"[{params['roi_low_keV']:.1f}, {params['roi_high_keV']:.1f}] keV",
         f"width = {params['roi_width_keV']:.1f} keV"],
        ["dz_min", f"{params['dz_mm']} mm", ""],
        ["ks", f"{params['ks']}", ""],
        ["E_veto", f"{params['eveto_keV']} keV", ""],
        ["", "", ""],
        ["Total events", f"{counts['n_total']}", f"{fracs['f_fv_only'] + fracs['f_active_vetoed']:.6f}"],
        ["Active-vetoed", f"{counts['n_active_vetoed']}", f"{fracs['f_active_vetoed']:.6f}"],
        ["FV-only events", f"{counts['n_fv_only']}", f"{fracs['f_fv_only']:.6f}"],
        ["Rejected (MS)", f"{counts['n_ms']}", f"{fracs['f_ms']:.6f}"],
        ["SS events", f"{counts['n_ss']}", f"{fracs['f_ss']:.6f}"],
        ["SS outside ROI", f"{counts['n_ss_outside_roi']}", f"{fracs['f_ss_outside_roi']:.6f}"],
        ["SS in ROI", f"{counts['n_ss_in_roi']}", f"{fracs['f_ss_in_roi']:.6f}"],
    ]
    col_labels = ["Parameter / Cut", "Value", "Detail / Fraction"]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.axis("off")
    tbl = ax.table(cellText=table_data, colLabels=col_labels,
                   loc="center", cellLoc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(10)
    tbl.scale(1.0, 1.5)
    for (row, _), cell in tbl.get_celld().items():
        if row == 0:
            cell.set_facecolor("#4472C4")
            cell.set_text_props(color="white", weight="bold")
        elif row == 8:
            cell.set_facecolor("#E0E0E0")
        elif row >= 9 and row % 2 == 1:
            cell.set_facecolor("#D9E2F3")

    ax.set_title("SS/MS Background Analysis Summary", fontsize=12, pad=12)
    fig.tight_layout()
    outfile = outdir / "summary_table.png"
    fig.savefig(outfile, dpi=150, bbox_inches="tight")
    print(f"Saved {outfile}")
    if display:
        plt.show()
    plt.close(fig)


def plot_background_rate(smry: dict, outdir: Path, display: bool) -> None:
    br = smry["background_rate"]
    components = br.get("source_components", [])

    table_data = [
        ["Source", br["source"], ""],
        ["Isotope", br["isotope"], ""],
        ["gamma BR", f"{br['gamma_BR']}", ""],
    ]
    for comp in components:
        table_data.append([
            comp["name"],
            f"{comp['mass_kg']:.3f} kg",
            f"{comp['activity_mBq_per_kg']:.4f} mBq/kg",
        ])
        table_data.append([
            "",
            f"geom. surv. = {comp['geometric_survival']:.6f}",
            f"rate = {comp['rate_gammas_per_s']:.4e} g/s",
        ])
    table_data += [
        ["", "", ""],
        ["Rate (gammas/s)", f"{br['rate_gammas_per_s']:.4e}", ""],
        ["Rate (gammas/year)", f"{br['rate_gammas_per_year']:.4e}", ""],
        ["f(FV) [MC]", f"{br['f_fv_mc']:.4e}", f"N(FV)/yr = {br['n_fv_per_year']:.2f}"],
        ["f(no active veto)", f"{br['f_no_active_veto']:.6f}", ""],
        ["f(SS | FV-only)", f"{br['f_ss_of_fvonly']:.6f}", ""],
        ["f(ROI | SS)", f"{br['f_roi_of_ss']:.6f}", ""],
        ["N(SS in ROI)/year", f"{br['n_ss_roi_per_year']:.4f}", ""],
    ]

    n_rows = len(table_data)
    fig_h = max(4, 0.4 * n_rows + 1)
    fig, ax = plt.subplots(figsize=(10, fig_h))
    ax.axis("off")
    col_labels = ["Quantity", "Value", "Detail"]
    tbl = ax.table(cellText=table_data, colLabels=col_labels,
                   loc="center", cellLoc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    tbl.scale(1.0, 1.4)
    for (row, _), cell in tbl.get_celld().items():
        if row == 0:
            cell.set_facecolor("#4472C4")
            cell.set_text_props(color="white", weight="bold")
        elif row == n_rows:
            cell.set_facecolor("#C6EFCE")
            cell.set_text_props(weight="bold")

    ax.set_title("Background Rate Estimate", fontsize=12, pad=12)
    fig.tight_layout()
    outfile = outdir / "background_rate.png"
    fig.savefig(outfile, dpi=150, bbox_inches="tight")
    print(f"Saved {outfile}")
    if display:
        plt.show()
    plt.close(fig)


def main() -> None:
    args = parse_args()
    indir = args.input

    print(f"Reading from {indir} ...")
    summary = pd.read_csv(indir / "event_summary.csv")
    clf = pd.read_csv(indir / "classification.csv")
    evol = pd.read_csv(indir / "etot_by_volume.csv")
    esat = pd.read_csv(indir / "esat.csv")
    sep = pd.read_csv(indir / "separated_deposits.csv")
    with open(indir / "summary.json") as f:
        smry = json.load(f)

    print(f"  {len(summary)} events")

    print("Plotting ...")
    plot_event_summary(summary, esat, indir, args.display)
    plot_energy_by_volume(evol, indir, args.display)
    plot_ssms_analysis(clf, sep, smry, indir, args.display)
    plot_summary_table(smry, indir, args.display)
    plot_background_rate(smry, indir, args.display)

    print("Done.")


if __name__ == "__main__":
    main()
