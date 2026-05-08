"""
Final background computation: smear, apply ROI cut, compute absolute rates.

Reads candidates.csv and metadata from --bgdir (background results) and
--fluxdir (flux tables). Applies Gaussian energy smearing and ROI window.
Computes absolute background rates.

Usage:
    python py/compute_backgrounds.py \\
        --source cryostat_barrel --isotope Bi214 \\
        --bgdir results/backgrounds/cryostat/barrel/Bi214 \\
        --fluxdir results/fluxes/cryostat/barrel/Bi214 \\
        --sigma 25 --roi 50 \\
        --display --save
"""

import sys
import os
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


SECONDS_PER_YEAR = 365.25 * 24 * 3600

# Source geometry: mass and activity per component for each cryostat source.
# These match source_geometry_lz_v1.json. Activity in mBq/kg.
SOURCE_INFO = {
    "cryostat_barrel": {
        "Bi214": {
            "components": [
                {"name": "OCV_barrel", "activity_mBq_per_kg": 0.08},
                {"name": "MLI", "activity_mBq_per_kg": 11.1},
                {"name": "ICV_barrel", "activity_mBq_per_kg": 0.08},
            ]
        },
        "Tl208": {
            "components": [
                {"name": "OCV_barrel", "activity_mBq_per_kg": 0.22},
                {"name": "MLI", "activity_mBq_per_kg": 7.79},
                {"name": "ICV_barrel", "activity_mBq_per_kg": 0.22},
            ]
        },
    },
    "cryostat_top": {
        "Bi214": {
            "components": [
                {"name": "OCV_top", "activity_mBq_per_kg": 0.08},
                {"name": "ICV_top", "activity_mBq_per_kg": 0.08},
            ]
        },
        "Tl208": {
            "components": [
                {"name": "OCV_top", "activity_mBq_per_kg": 0.22},
                {"name": "ICV_top", "activity_mBq_per_kg": 0.22},
            ]
        },
    },
    "cryostat_bottom": {
        "Bi214": {
            "components": [
                {"name": "OCV_bottom", "activity_mBq_per_kg": 0.08},
                {"name": "ICV_bottom", "activity_mBq_per_kg": 0.08},
            ]
        },
        "Tl208": {
            "components": [
                {"name": "OCV_bottom", "activity_mBq_per_kg": 0.22},
                {"name": "ICV_bottom", "activity_mBq_per_kg": 0.22},
            ]
        },
    },
}


def parse_cli(args):
    opts = {
        "source": None, "isotope": None,
        "bgdir": None, "fluxdir": None,
        "sigma": 25.0, "roi": 50.0,
        "display": False, "save": True,
    }
    i = 0
    while i < len(args):
        a = args[i]
        if a in ("-h", "--help"):
            print(__doc__)
            sys.exit(0)
        elif a == "--source":
            opts["source"] = args[i+1]; i += 2; continue
        elif a == "--isotope":
            opts["isotope"] = args[i+1]; i += 2; continue
        elif a == "--bgdir":
            opts["bgdir"] = args[i+1]; i += 2; continue
        elif a == "--fluxdir":
            opts["fluxdir"] = args[i+1]; i += 2; continue
        elif a == "--sigma":
            opts["sigma"] = float(args[i+1]); i += 2; continue
        elif a == "--roi":
            opts["roi"] = float(args[i+1]); i += 2; continue
        elif a == "--display":
            opts["display"] = True; i += 1; continue
        elif a == "--no-save":
            opts["save"] = False; i += 1; continue
        else:
            print(f"Unknown argument: {a}")
            sys.exit(1)
        i += 1

    for key in ("source", "isotope", "bgdir", "fluxdir"):
        if opts[key] is None:
            print(f"--{key} is required")
            sys.exit(1)
    return opts


def main():
    opts = parse_cli(sys.argv[1:])
    source = opts["source"]
    isotope = opts["isotope"]
    sigma_keV = opts["sigma"]
    roi_keV = opts["roi"]
    Q_bb_MeV = 2.458  # Xe-136

    sigma_MeV = sigma_keV / 1000.0
    roi_MeV = roi_keV / 1000.0
    roi_lo = Q_bb_MeV - roi_MeV
    roi_hi = Q_bb_MeV + roi_MeV

    # --- Load data ---
    candidates = pd.read_csv(os.path.join(opts["bgdir"], "candidates.csv"))
    bg_meta = json.load(open(os.path.join(opts["bgdir"], "metadata.json")))
    stats_df = pd.read_csv(os.path.join(opts["bgdir"], "statistics.csv"))
    stats = dict(zip(stats_df["quantity"], stats_df["value"]))

    flux_meta = json.load(open(os.path.join(opts["fluxdir"], "metadata.json")))

    E_true = candidates["E_deposited_MeV"].values
    N_sampled = int(float(bg_meta["N_actual"]))
    gamma_rate = float(bg_meta["rate_total_gammas_per_s"])

    n_accepted_ss = int(float(stats.get("n_accepted", 0)))
    n_ms_rejected = int(float(stats.get("n_ms_rejected", 0)))
    n_vetoed = int(float(stats.get("n_vetoed", 0)))
    n_no_fv = int(float(stats.get("n_no_fv", 0)))
    n_fv_total = n_accepted_ss + n_ms_rejected

    # --- Smear and ROI cut ---
    rng = np.random.default_rng(42)
    E_smeared = E_true + rng.normal(0, sigma_MeV, size=len(E_true))
    in_roi = (E_smeared >= roi_lo) & (E_smeared <= roi_hi)
    n_roi = int(np.sum(in_roi))

    # --- Absolute rates ---
    gammas_per_year = gamma_rate * SECONDS_PER_YEAR

    f_fv = n_fv_total / N_sampled if N_sampled > 0 else 0.0
    f_ss = n_accepted_ss / N_sampled if N_sampled > 0 else 0.0
    f_roi = n_roi / N_sampled if N_sampled > 0 else 0.0

    bg_fv_per_year = f_fv * gammas_per_year
    bg_ss_per_year = f_ss * gammas_per_year
    bg_ms_per_year = (n_ms_rejected / N_sampled) * gammas_per_year if N_sampled > 0 else 0.0
    bg_roi_per_year = f_roi * gammas_per_year

    ss_fraction_in_fv = n_accepted_ss / n_fv_total if n_fv_total > 0 else 0.0

    # --- Print summary ---
    print("=" * 60)
    print(f"  Background Summary: {source} / {isotope}")
    print("=" * 60)
    print(f"  Sigma:          {sigma_keV:.1f} keV")
    print(f"  ROI:            Q_bb +/- {roi_keV:.1f} keV = [{roi_lo*1000:.1f}, {roi_hi*1000:.1f}] keV")
    print(f"  Q_bb (Xe-136):  {Q_bb_MeV*1000:.1f} keV")
    print()
    print(f"  Gamma flux through ICV surface:  {gamma_rate:.6e} gammas/s")
    print(f"  Gammas per year:                 {gammas_per_year:.6e}")
    print()
    print(f"  N sampled:                       {N_sampled}")
    print(f"  N reaching FV (SS+MS):           {n_fv_total}")
    print(f"  N single-site (SS):              {n_accepted_ss}")
    print(f"  N multi-site (MS):               {n_ms_rejected}")
    print(f"  N vetoed:                        {n_vetoed}")
    print(f"  N no FV:                         {n_no_fv}")
    print(f"  SS fraction in FV:               {ss_fraction_in_fv:.4f}")
    print()
    print(f"  Background in FV:                {bg_fv_per_year:.4e} events/year")
    print(f"  Background SS:                   {bg_ss_per_year:.4e} events/year")
    print(f"  Background MS:                   {bg_ms_per_year:.4e} events/year")
    print(f"  N SS in ROI (after smearing):    {n_roi}")
    print(f"  Background SS in ROI:            {bg_roi_per_year:.4e} events/year")
    print("=" * 60)

    # --- Plots ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Left: true energy, no smearing
    bins = np.linspace(E_true.min() * 0.98, E_true.max() * 1.02, 60) if len(E_true) > 0 else np.linspace(2.0, 2.6, 60)
    ax1.hist(E_true, bins=bins, color="steelblue", edgecolor="black",
             linewidth=0.5, label="SS candidates (true E)")
    ax1.axvline(roi_lo, color="red", linestyle="--", linewidth=1.5, label=f"ROI [{roi_lo*1000:.0f}, {roi_hi*1000:.0f}] keV")
    ax1.axvline(roi_hi, color="red", linestyle="--", linewidth=1.5)
    ax1.axvspan(roi_lo, roi_hi, alpha=0.1, color="red")
    ax1.set_xlabel("E deposited [MeV]")
    ax1.set_ylabel("Counts")
    ax1.set_title("True energy (no smearing)")
    ax1.legend(fontsize=8)

    # Right: smeared energy with ROI
    ax2.hist(E_smeared, bins=bins, color="darkorange", edgecolor="black",
             linewidth=0.5, alpha=0.7, label="Smeared E")
    ax2.hist(E_smeared[in_roi], bins=bins, color="red", edgecolor="black",
             linewidth=0.5, alpha=0.5, label=f"In ROI ({n_roi} events)")
    ax2.axvline(roi_lo, color="red", linestyle="--", linewidth=1.5)
    ax2.axvline(roi_hi, color="red", linestyle="--", linewidth=1.5)
    ax2.axvspan(roi_lo, roi_hi, alpha=0.1, color="red")
    ax2.set_xlabel("E smeared [MeV]")
    ax2.set_ylabel("Counts")
    ax2.set_title(f"Smeared (sigma={sigma_keV:.0f} keV)")
    ax2.legend(fontsize=8)

    suptitle = f"{source} / {isotope} — BG in ROI: {bg_roi_per_year:.2e} events/year"
    fig.suptitle(suptitle, fontsize=12, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93])

    if opts["save"]:
        outpath = os.path.join(opts["bgdir"], "background_summary.png")
        fig.savefig(outpath, dpi=150)
        print(f"\nSaved: {outpath}")

        # Also save summary as JSON
        summary = {
            "source": source, "isotope": isotope,
            "sigma_keV": sigma_keV, "roi_keV": roi_keV,
            "Q_bb_MeV": Q_bb_MeV,
            "gamma_rate_per_s": gamma_rate,
            "gammas_per_year": gammas_per_year,
            "N_sampled": N_sampled,
            "n_fv_total": n_fv_total,
            "n_ss": n_accepted_ss,
            "n_ms": n_ms_rejected,
            "n_vetoed": n_vetoed,
            "n_no_fv": n_no_fv,
            "ss_fraction_in_fv": ss_fraction_in_fv,
            "n_ss_in_roi": n_roi,
            "bg_fv_per_year": bg_fv_per_year,
            "bg_ss_per_year": bg_ss_per_year,
            "bg_ms_per_year": bg_ms_per_year,
            "bg_roi_per_year": bg_roi_per_year,
        }
        json_path = os.path.join(opts["bgdir"], "background_summary.json")
        with open(json_path, "w") as f:
            json.dump(summary, f, indent=2)
        print(f"Saved: {json_path}")

    if opts["display"]:
        plt.show()
    else:
        plt.close(fig)


if __name__ == "__main__":
    main()
