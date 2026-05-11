"""Analyze FV deposit data: per-event histograms of deposit multiplicity,
Emax, Esat, Etot, DZ, spatial distributions, and SS/MS classification.

All energies in keV.
"""

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

Q_BB_KEV = 2458.07
FWHM_FACTOR = 2.3548200  # 2 * sqrt(2 * ln(2))
E_GAMMA_BI214_KEV = 2447.86  # highest-energy Bi-214 gamma line


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Histogram analysis of FV deposits from background simulation."
    )
    parser.add_argument(
        "-i", "--input", required=True, type=Path,
        help="Input directory containing fv_deposits.csv and metadata.json"
    )
    parser.add_argument(
        "-o", "--output", type=Path, default=None,
        help="Output directory (default: <input>/analysis)"
    )
    parser.add_argument(
        "--display", action="store_true", help="Show plots interactively"
    )
    parser.add_argument(
        "--sigma", type=float, default=0.01,
        help="Relative energy resolution sigma/E (default: 0.01)"
    )
    parser.add_argument(
        "--roi", type=float, nargs=2, default=None, metavar=("LOW", "HIGH"),
        help="ROI bounds in keV (default: 1 FWHM centered on Qbb)"
    )
    parser.add_argument(
        "--dz", type=float, default=3.0,
        help="Minimum z-separation for cluster splitting (mm, default: 3.0)"
    )
    parser.add_argument(
        "--ks", type=float, default=3.0,
        help="Significance factor: E_i > ks * sigma(E1) (default: 3.0)"
    )
    parser.add_argument(
        "--eveto", type=float, default=10.0,
        help="Active veto threshold in keV (default: 10.0)"
    )
    return parser.parse_args()


# =========================================================================
# Stage 1: raw per-event summary (for event_summary.png)
# =========================================================================

def analyze(df: pd.DataFrame) -> tuple[pd.DataFrame, np.ndarray]:
    """Compute per-event summary and collect satellite deposit energies."""
    records = []
    esat_list: list[float] = []

    for event_id, grp in df.groupby("event_id"):
        n_deposits = len(grp)
        idx_max = grp["energy_keV"].idxmax()
        emax = grp.loc[idx_max, "energy_keV"]
        etot = grp["energy_keV"].sum()

        z_max = grp.loc[idx_max, "z_cm"]
        x_max = grp.loc[idx_max, "x_cm"]
        y_max = grp.loc[idx_max, "y_cm"]
        r_max = np.sqrt(x_max**2 + y_max**2)

        others = grp.drop(idx_max)
        if len(others) > 0:
            esat_list.extend(others["energy_keV"].values)
            z_centroid = np.average(others["z_cm"], weights=others["energy_keV"])
            dz = abs(z_max - z_centroid)
        else:
            dz = 0.0

        records.append(
            {
                "event_id": event_id,
                "n_deposits": n_deposits,
                "Emax_keV": emax,
                "Etot_keV": etot,
                "DZ_cm": dz,
                "R_max_cm": r_max,
                "Z_max_cm": z_max,
            }
        )

    summary = pd.DataFrame(records)
    esat_keV = np.array(esat_list)
    return summary, esat_keV


# =========================================================================
# Stage 2: active veto + SS/MS classification
# =========================================================================
#
# Selection steps (applied per event):
#
#   Step 1 — ACTIVE VETO: scan all deposits. If any deposit has
#            volume == "active" and energy >= eveto_keV, reject.
#
#   Step 2 — FV-ONLY: keep only deposits with volume == "fv".
#            (Passive deposits are lost energy, not vetoed.)
#
#   Step 3 — SS/MS CLASSIFICATION on FV deposits:
#            E1 = largest FV deposit, z1 = its z position.
#            For each other FV deposit E_i:
#              if |z_i - z1| > dz AND E_i > ks * sigma(E1) -> separated
#              else -> attached to E1
#            SS = 0 separated clusters, MS = >= 1.
#
#   Step 4 — E_CLUSTER (SS only): E1 + sum(attached E_i).
#
#   Step 5 — SMEARING: E_smeared = Gaussian(E_cluster, sigma * E_cluster).
#
#   Step 6 — ROI CUT: roi_low <= E_smeared <= roi_high.

def classify_events(
    df: pd.DataFrame,
    sigma_rel: float,
    dz_mm: float,
    ks: float,
    eveto_keV: float,
    rng: np.random.Generator,
) -> pd.DataFrame:
    """Apply active veto, then classify surviving events as SS/MS."""
    records = []

    for event_id, grp in df.groupby("event_id"):
        # --- Step 1: active veto ---
        active_deps = grp[grp["volume"] == "active"]
        active_vetoed = (active_deps["energy_keV"] >= eveto_keV).any()

        if active_vetoed:
            records.append(
                {
                    "event_id": event_id,
                    "active_vetoed": True,
                    "n_separated": 0,
                    "sep_energies": [],
                    "sep_dz_mm": [],
                    "is_ss": False,
                    "E_cluster_keV": 0.0,
                    "E_smeared_keV": 0.0,
                }
            )
            continue

        # --- Step 2: keep only FV deposits ---
        fv = grp[grp["volume"] == "fv"]
        if len(fv) == 0:
            records.append(
                {
                    "event_id": event_id,
                    "active_vetoed": False,
                    "n_separated": 0,
                    "sep_energies": [],
                    "sep_dz_mm": [],
                    "is_ss": True,
                    "E_cluster_keV": 0.0,
                    "E_smeared_keV": 0.0,
                }
            )
            continue

        # --- Step 3: SS/MS classification on FV deposits ---
        idx_max = fv["energy_keV"].idxmax()
        e1 = fv.loc[idx_max, "energy_keV"]
        z1_cm = fv.loc[idx_max, "z_cm"]
        sigma_e1 = sigma_rel * e1
        e_threshold = ks * sigma_e1
        dz_threshold_cm = dz_mm / 10.0

        others = fv.drop(idx_max)
        sep_energies: list[float] = []
        sep_dz: list[float] = []
        attached_energy = 0.0

        for _, row in others.iterrows():
            dz_i = abs(row["z_cm"] - z1_cm)
            ei = row["energy_keV"]
            if dz_i > dz_threshold_cm and ei > e_threshold:
                sep_energies.append(ei)
                sep_dz.append(dz_i * 10.0)  # cm -> mm
            else:
                attached_energy += ei

        n_separated = len(sep_energies)
        is_ss = n_separated == 0

        # --- Step 4: E_cluster ---
        e_cluster = e1 + attached_energy

        # --- Step 5: smearing ---
        sigma_cluster = sigma_rel * e_cluster
        e_smeared = rng.normal(e_cluster, sigma_cluster)

        records.append(
            {
                "event_id": event_id,
                "active_vetoed": False,
                "n_separated": n_separated,
                "sep_energies": sep_energies,
                "sep_dz_mm": sep_dz,
                "is_ss": is_ss,
                "E_cluster_keV": e_cluster,
                "E_smeared_keV": e_smeared,
            }
        )

    return pd.DataFrame(records)


# =========================================================================
# Plotting
# =========================================================================

def plot_event_summary(
    summary: pd.DataFrame, esat_keV: np.ndarray, outdir: Path, display: bool
) -> None:
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
    ax.hist(esat_keV, bins=100, edgecolor="black", linewidth=0.3)
    ax.set_xlabel("Esat (keV)")
    ax.set_ylabel("Counts")
    ax.set_title("Esat: each non-Emax deposit")

    ax = axes[1, 0]
    ax.hist(summary["Etot_keV"].values, bins=100, edgecolor="black", linewidth=0.3)
    ax.set_xlabel("Etot (keV)")
    ax.set_ylabel("Counts")
    ax.set_title("Etot: total energy per event")

    ax = axes[1, 1]
    dz_mm = summary["DZ_cm"].values * 10.0
    ax.hist(dz_mm, bins=100, edgecolor="black", linewidth=0.3)
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


def plot_energy_by_volume(df: pd.DataFrame, outdir: Path, display: bool) -> None:
    """2x2 panel: Etot (all), Etot(fv), Etot(passive), Etot(active) per event."""
    etot_all = df.groupby("event_id")["energy_keV"].sum()
    etot_fv = df[df["volume"] == "fv"].groupby("event_id")["energy_keV"].sum()
    etot_passive = df[df["volume"] == "passive"].groupby("event_id")["energy_keV"].sum()
    etot_active = df[df["volume"] == "active"].groupby("event_id")["energy_keV"].sum()

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    ax = axes[0, 0]
    ax.hist(etot_all.values, bins=100, edgecolor="black", linewidth=0.3)
    ax.axvline(E_GAMMA_BI214_KEV, color="green", ls="--", lw=0.8,
               label=f"Bi214={E_GAMMA_BI214_KEV:.0f}")
    ax.legend()
    ax.set_xlabel("Energy (keV)")
    ax.set_ylabel("Counts")
    ax.set_title("Etot: all deposits")

    ax = axes[0, 1]
    ax.hist(etot_fv.values, bins=100, edgecolor="black", linewidth=0.3)
    ax.axvline(E_GAMMA_BI214_KEV, color="green", ls="--", lw=0.8,
               label=f"Bi214={E_GAMMA_BI214_KEV:.0f}")
    ax.legend()
    ax.set_xlabel("Energy (keV)")
    ax.set_ylabel("Counts")
    ax.set_title("Etot: FV deposits only")

    ax = axes[1, 0]
    if len(etot_passive) > 0:
        ax.hist(etot_passive.values, bins=100, edgecolor="black", linewidth=0.3)
    ax.set_xlabel("Energy (keV)")
    ax.set_ylabel("Counts")
    ax.set_title("Etot: passive deposits only")

    ax = axes[1, 1]
    if len(etot_active) > 0:
        ax.hist(etot_active.values, bins=100, edgecolor="black", linewidth=0.3)
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


def plot_ssms_analysis(
    clf: pd.DataFrame, roi_low: float, roi_high: float, outdir: Path, display: bool
) -> None:
    # Work only with non-vetoed events for SS/MS plots
    surv = clf[~clf["active_vetoed"]]

    all_sep_dz: list[float] = []
    for _, row in surv.iterrows():
        all_sep_dz.extend(row["sep_dz_mm"])
    all_sep_dz_arr = np.array(all_sep_dz) if all_sep_dz else np.array([])

    ss_mask = surv["is_ss"].values
    ss_ecluster = surv.loc[ss_mask, "E_cluster_keV"].values
    ss_esmeared = surv.loc[ss_mask, "E_smeared_keV"].values

    n_total = len(clf)
    n_vetoed = int(clf["active_vetoed"].sum())
    n_surv = len(surv)

    fig, axes = plt.subplots(2, 3, figsize=(14, 8))

    # (0,0) Veto bar chart: FV-only vs active-vetoed (fractions)
    ax = axes[0, 0]
    f_fv = n_surv / n_total if n_total > 0 else 0
    f_vetoed = n_vetoed / n_total if n_total > 0 else 0
    ax.bar(["FV-only", "Active-vetoed"], [f_fv, f_vetoed],
           color=["steelblue", "tomato"], edgecolor="black", linewidth=0.5)
    ax.set_ylabel("Fraction")
    ax.set_title(f"FV-only={n_surv}  Vetoed={n_vetoed}")
    ax.set_ylim(0, 1)

    # (0,1) Cluster multiplicity (FV-only events), SS bar in different color
    ax = axes[0, 1]
    max_sep = surv["n_separated"].max() if len(surv) > 0 else 1
    n_sep = surv["n_separated"].values
    bins = list(range(0, max_sep + 1))
    counts_per_bin = [int((n_sep == k).sum()) for k in bins]
    bar_colors = ["steelblue" if k == 0 else "coral" for k in bins]
    ax.bar(bins, counts_per_bin, color=bar_colors,
           edgecolor="black", linewidth=0.3)
    ax.set_xlabel("N separated clusters (0 = SS)")
    ax.set_ylabel("Counts")
    ax.set_title("Cluster multiplicity (FV-only)")

    # (0,2) SS vs MS bar chart (FV-only events)
    ax = axes[0, 2]
    n_ss = int(ss_mask.sum())
    n_ms = int((~ss_mask).sum())
    ax.bar(["SS", "MS"], [n_ss, n_ms], color=["steelblue", "coral"],
           edgecolor="black", linewidth=0.5)
    ax.set_ylabel("Counts")
    ax.set_title(f"SS vs MS  (SS={n_ss}, MS={n_ms})")

    # (1,0) DZ of separated deposits
    ax = axes[1, 0]
    if len(all_sep_dz_arr) > 0:
        ax.hist(all_sep_dz_arr, bins=100, edgecolor="black", linewidth=0.3)
    ax.set_xlabel("DZ (mm)")
    ax.set_ylabel("Counts")
    ax.set_title("DZ of separated deposits")

    # (1,1) E_cluster for SS events
    ax = axes[1, 1]
    if len(ss_ecluster) > 0:
        ax.hist(ss_ecluster, bins=100, edgecolor="black", linewidth=0.3)
    ax.axvline(E_GAMMA_BI214_KEV, color="green", ls="--", lw=0.8,
               label=f"Bi214={E_GAMMA_BI214_KEV:.0f}")
    ax.legend()
    ax.set_xlabel("E_cluster (keV)")
    ax.set_ylabel("Counts")
    ax.set_title("E_cluster (SS events)")

    # (1,2) E_smeared for SS events - zoomed around Qbb
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


def write_summary(
    clf: pd.DataFrame,
    roi_low: float,
    roi_high: float,
    sigma_rel: float,
    dz_mm: float,
    ks: float,
    eveto_keV: float,
    outdir: Path,
    display: bool,
) -> None:
    n_total = len(clf)
    n_vetoed = int(clf["active_vetoed"].sum())
    surv = clf[~clf["active_vetoed"]]
    n_fv_only = len(surv)
    ss_mask = surv["is_ss"].values
    n_ss = int(ss_mask.sum())
    n_ms = n_fv_only - n_ss
    ss_esmeared = surv.loc[ss_mask, "E_smeared_keV"].values
    roi_mask = (ss_esmeared >= roi_low) & (ss_esmeared <= roi_high)
    n_ss_roi = int(roi_mask.sum())
    n_ss_out = n_ss - n_ss_roi

    sigma_abs = sigma_rel * Q_BB_KEV
    fwhm = sigma_abs * FWHM_FACTOR
    roi_width = roi_high - roi_low

    def frac(n: int) -> str:
        return f"{n / n_total:.6f}" if n_total > 0 else "0"

    lines = [
        "=" * 60,
        "SS/MS Background Analysis Summary",
        "=" * 60,
        f"  sigma/E       = {sigma_rel}",
        f"  sigma(Qbb)    = {sigma_abs:.1f} keV",
        f"  FWHM(Qbb)     = {fwhm:.1f} keV",
        f"  ROI           = [{roi_low:.1f}, {roi_high:.1f}] keV",
        f"  ROI width     = {roi_width:.1f} keV",
        f"  Qbb           = {Q_BB_KEV:.2f} keV",
        f"  dz_min        = {dz_mm} mm",
        f"  ks            = {ks}",
        f"  E_veto        = {eveto_keV} keV",
        "-" * 60,
        f"  Total events          : {n_total:>8d}   ({frac(n_total)})",
        f"  Active-vetoed         : {n_vetoed:>8d}   ({frac(n_vetoed)})",
        f"  FV-only events        : {n_fv_only:>8d}   ({frac(n_fv_only)})",
        f"  Rejected (MS)         : {n_ms:>8d}   ({frac(n_ms)})",
        f"  SS events             : {n_ss:>8d}   ({frac(n_ss)})",
        f"  SS outside ROI        : {n_ss_out:>8d}   ({frac(n_ss_out)})",
        f"  SS in ROI             : {n_ss_roi:>8d}   ({frac(n_ss_roi)})",
        "=" * 60,
    ]

    text = "\n".join(lines)
    print(text)

    outfile = outdir / "summary.txt"
    outfile.write_text(text + "\n")
    print(f"Saved {outfile}")

    # Summary table as PNG
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.axis("off")
    table_data = [
        ["sigma/E", f"{sigma_rel}", ""],
        ["sigma(Qbb)", f"{sigma_abs:.1f} keV", ""],
        ["FWHM(Qbb)", f"{fwhm:.1f} keV", ""],
        ["ROI", f"[{roi_low:.1f}, {roi_high:.1f}] keV", f"width = {roi_width:.1f} keV"],
        ["dz_min", f"{dz_mm} mm", ""],
        ["ks", f"{ks}", ""],
        ["E_veto", f"{eveto_keV} keV", ""],
        ["", "", ""],
        ["Total events", f"{n_total}", f"{frac(n_total)}"],
        ["Active-vetoed", f"{n_vetoed}", f"{frac(n_vetoed)}"],
        ["FV-only events", f"{n_fv_only}", f"{frac(n_fv_only)}"],
        ["Rejected (MS)", f"{n_ms}", f"{frac(n_ms)}"],
        ["SS events", f"{n_ss}", f"{frac(n_ss)}"],
        ["SS outside ROI", f"{n_ss_out}", f"{frac(n_ss_out)}"],
        ["SS in ROI", f"{n_ss_roi}", f"{frac(n_ss_roi)}"],
    ]
    col_labels = ["Parameter / Cut", "Value", "Detail / Fraction"]
    tbl = ax.table(
        cellText=table_data,
        colLabels=col_labels,
        loc="center",
        cellLoc="center",
    )
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


SECONDS_PER_YEAR = 3.15576e7


def write_background_rate(
    clf: pd.DataFrame,
    roi_low: float,
    roi_high: float,
    metadata: dict,
    outdir: Path,
    display: bool,
) -> None:
    """Compute and write absolute background rate per year."""
    n_total = len(clf)
    n_vetoed = int(clf["active_vetoed"].sum())
    n_fv_only = n_total - n_vetoed
    surv = clf[~clf["active_vetoed"]]
    ss_mask = surv["is_ss"].values
    n_ss = int(ss_mask.sum())
    ss_esmeared = surv.loc[ss_mask, "E_smeared_keV"].values
    n_ss_roi = int(((ss_esmeared >= roi_low) & (ss_esmeared <= roi_high)).sum())

    # Fractions from analysis (relative to FV events entering this script)
    f_fv_only = n_fv_only / n_total if n_total > 0 else 0.0
    f_ss_of_fvonly = n_ss / n_fv_only if n_fv_only > 0 else 0.0
    f_roi_of_ss = n_ss_roi / n_ss if n_ss > 0 else 0.0

    # From metadata
    source = metadata.get("source", "unknown")
    isotope = metadata.get("isotope", "unknown")
    rate_gammas_per_s = metadata.get("rate_total_gammas_per_s", 0.0)
    f_fv_mc = metadata.get("f_fv", 0.0)
    gamma_BR = metadata.get("gamma_BR", 0.0)
    components = metadata.get("source_components", [])

    rate_gammas_per_year = rate_gammas_per_s * SECONDS_PER_YEAR
    n_fv_per_year = rate_gammas_per_year * f_fv_mc
    n_ss_roi_per_year = n_fv_per_year * f_fv_only * f_ss_of_fvonly * f_roi_of_ss

    # --- Text file ---
    lines = [
        "=" * 70,
        "Background Rate Estimate",
        "=" * 70,
        f"  Source          : {source}",
        f"  Isotope         : {isotope}",
        f"  gamma BR        : {gamma_BR}",
        "-" * 70,
    ]
    for comp in components:
        lines.append(f"  Component: {comp['name']}")
        lines.append(f"    mass            = {comp['mass_kg']:.3f} kg")
        lines.append(f"    activity        = {comp['activity_mBq_per_kg']:.4f} mBq/kg")
        lines.append(f"    geom. survival  = {comp['geometric_survival']:.6f}")
        lines.append(f"    rate            = {comp['rate_gammas_per_s']:.6e} gammas/s")
    lines += [
        "-" * 70,
        f"  Total rate       : {rate_gammas_per_s:.6e} gammas/s",
        f"  Total rate       : {rate_gammas_per_year:.6e} gammas/year",
        f"  f(FV) [MC]       : {f_fv_mc:.6e}",
        f"  N(FV)/year       : {n_fv_per_year:.2f}",
        f"  f(no active veto): {f_fv_only:.6f}",
        f"  f(SS | FV-only)  : {f_ss_of_fvonly:.6f}",
        f"  f(ROI | SS)      : {f_roi_of_ss:.6f}",
        "-" * 70,
        f"  N(SS in ROI)/year: {n_ss_roi_per_year:.4f}",
        "=" * 70,
    ]

    text = "\n".join(lines)
    print(text)

    outfile = outdir / "background_rate.txt"
    outfile.write_text(text + "\n")
    print(f"Saved {outfile}")

    # --- Table PNG ---
    table_data = [
        ["Source", source, ""],
        ["Isotope", isotope, ""],
        ["gamma BR", f"{gamma_BR}", ""],
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
        ["Rate (gammas/s)", f"{rate_gammas_per_s:.4e}", ""],
        ["Rate (gammas/year)", f"{rate_gammas_per_year:.4e}", ""],
        ["f(FV) [MC]", f"{f_fv_mc:.4e}", f"N(FV)/yr = {n_fv_per_year:.2f}"],
        ["f(no active veto)", f"{f_fv_only:.6f}", ""],
        ["f(SS | FV-only)", f"{f_ss_of_fvonly:.6f}", ""],
        ["f(ROI | SS)", f"{f_roi_of_ss:.6f}", ""],
        ["N(SS in ROI)/year", f"{n_ss_roi_per_year:.4f}", ""],
    ]

    n_rows = len(table_data)
    fig_h = max(4, 0.4 * n_rows + 1)
    fig, ax = plt.subplots(figsize=(10, fig_h))
    ax.axis("off")
    col_labels = ["Quantity", "Value", "Detail"]
    tbl = ax.table(
        cellText=table_data,
        colLabels=col_labels,
        loc="center",
        cellLoc="center",
    )
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


# =========================================================================
# Main
# =========================================================================

def main() -> None:
    args = parse_args()
    indir = args.input
    if args.output is None:
        args.output = indir / "analysis"
    args.output.mkdir(parents=True, exist_ok=True)

    deposits_path = indir / "fv_deposits.csv"
    metadata_path = indir / "metadata.json"

    sigma_rel = args.sigma
    sigma_abs = sigma_rel * Q_BB_KEV
    fwhm = sigma_abs * FWHM_FACTOR
    if args.roi is not None:
        roi_low, roi_high = args.roi
    else:
        roi_low = Q_BB_KEV - fwhm / 2.0
        roi_high = Q_BB_KEV + fwhm / 2.0

    print(f"Reading {deposits_path} ...")
    df = pd.read_csv(deposits_path)
    df["energy_keV"] = df["energy_MeV"] * 1000.0
    n_events = df["event_id"].nunique()
    print(f"  {len(df)} deposits, {n_events} events")
    print(f"  sigma/E = {sigma_rel:.4f}  ->  sigma = {sigma_abs:.1f} keV, "
          f"FWHM = {fwhm:.1f} keV")
    print(f"  ROI = [{roi_low:.1f}, {roi_high:.1f}] keV")
    print(f"  dz = {args.dz} mm, ks = {args.ks}, eveto = {args.eveto} keV")

    # Stage 1: raw summary
    print("Computing per-event summary ...")
    summary, esat_keV = analyze(df)

    etot = summary["Etot_keV"].values
    print(f"\n  Etot distribution (all FV events):")
    print(f"    mean   = {etot.mean():.1f} keV")
    print(f"    median = {np.median(etot):.1f} keV")
    print(f"    min    = {etot.min():.1f} keV")
    print(f"    max    = {etot.max():.1f} keV")
    print(f"    Etot > 2400 keV : {(etot > 2400).sum():>6d}  ({(etot > 2400).sum()/len(etot):.4f})")
    print(f"    Etot in [2440,2460] : {((etot >= 2440) & (etot <= 2460)).sum():>6d}")
    print(f"    Etot < 2400 keV : {(etot < 2400).sum():>6d}  ({(etot < 2400).sum()/len(etot):.4f})")
    print(f"    Etot < 2000 keV : {(etot < 2000).sum():>6d}  ({(etot < 2000).sum()/len(etot):.4f})")
    print(f"    Etot < 1000 keV : {(etot < 1000).sum():>6d}  ({(etot < 1000).sum()/len(etot):.4f})")
    print()

    plot_event_summary(summary, esat_keV, args.output, args.display)

    # Energy by volume (2x2 panel)
    print("Plotting energy by volume ...")
    plot_energy_by_volume(df, args.output, args.display)

    # Stage 2: active veto + SS/MS classification
    print("Running active veto + SS/MS classification ...")
    rng = np.random.default_rng(42)
    clf = classify_events(df, sigma_rel, args.dz, args.ks, args.eveto, rng)
    plot_ssms_analysis(clf, roi_low, roi_high, args.output, args.display)

    # Summary statistics
    write_summary(clf, roi_low, roi_high, sigma_rel, args.dz, args.ks,
                  args.eveto, args.output, args.display)

    # Background rate estimate
    if metadata_path.exists():
        print(f"\nReading {metadata_path} ...")
        with open(metadata_path) as f:
            metadata = json.load(f)
        write_background_rate(clf, roi_low, roi_high, metadata,
                              args.output, args.display)
    else:
        print(f"\n  {metadata_path} not found, skipping rate estimate.")

    print("Done.")


if __name__ == "__main__":
    main()
