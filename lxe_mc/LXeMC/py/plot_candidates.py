"""
Plot SS candidate distributions from background processing.

Reads candidates.csv and statistics.csv from a background results directory.
Produces a 2x2 figure: E distribution, x distribution, y distribution,
z distribution. Prints statistics summary to stdout.

Usage:
    python py/plot_candidates.py results/backgrounds/cryostat/barrel/Bi214
    python py/plot_candidates.py results/backgrounds/cryostat/barrel/Bi214 --display
    python py/plot_candidates.py results/backgrounds/cryostat/barrel/Bi214 --no-save
"""

import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def load_candidates(indir: str) -> pd.DataFrame:
    path = os.path.join(indir, "candidates.csv")
    if not os.path.exists(path):
        print(f"No candidates.csv found in {indir}")
        sys.exit(1)
    return pd.read_csv(path)


def load_statistics(indir: str) -> dict:
    path = os.path.join(indir, "statistics.csv")
    if not os.path.exists(path):
        return {}
    df = pd.read_csv(path)
    return dict(zip(df["quantity"], df["value"]))


def make_suptitle(indir: str) -> str:
    parts = os.path.normpath(indir).split(os.sep)
    if len(parts) >= 2:
        return f"{parts[-2]} / {parts[-1]}"
    return parts[-1]


def main():
    if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
        print("Usage: python py/plot_candidates.py <directory> [--display] [--no-save]")
        sys.exit(0)

    indir = sys.argv[1]
    display = "--display" in sys.argv
    save = "--no-save" not in sys.argv

    df = load_candidates(indir)
    stats = load_statistics(indir)

    # Print summary
    print(f"Directory: {indir}")
    print(f"Candidates: {len(df)}")
    if stats:
        for k, v in stats.items():
            print(f"  {k}: {v}")

    if len(df) == 0:
        print("No candidates to plot.")
        sys.exit(0)

    E = df["E_deposited_MeV"].values
    x = df["x_cm"].values
    y = df["y_cm"].values
    z = df["z_cm"].values
    r = np.sqrt(x**2 + y**2)

    fig, axes = plt.subplots(2, 2, figsize=(10, 8))

    # Energy distribution
    ax = axes[0, 0]
    ax.hist(E, bins=50, color="steelblue", edgecolor="black", linewidth=0.5)
    ax.set_xlabel("E deposited [MeV]")
    ax.set_ylabel("Counts")
    ax.set_title("Energy distribution")
    ax.axvline(2.458, color="red", linestyle="--", linewidth=1, label="Q_bb (Xe-136)")
    ax.legend(fontsize=8)

    # x distribution
    ax = axes[0, 1]
    ax.hist(x, bins=50, color="darkorange", edgecolor="black", linewidth=0.5)
    ax.set_xlabel("x [cm]")
    ax.set_ylabel("Counts")
    ax.set_title("x distribution")

    # y distribution
    ax = axes[1, 0]
    ax.hist(y, bins=50, color="seagreen", edgecolor="black", linewidth=0.5)
    ax.set_xlabel("y [cm]")
    ax.set_ylabel("Counts")
    ax.set_title("y distribution")

    # z distribution
    ax = axes[1, 1]
    ax.hist(z, bins=50, color="mediumpurple", edgecolor="black", linewidth=0.5)
    ax.set_xlabel("z [cm]")
    ax.set_ylabel("Counts")
    ax.set_title("z distribution")

    suptitle = make_suptitle(indir)
    n_cand = len(df)
    n_total = int(float(stats.get("n_total", 0)))
    if n_total > 0:
        suptitle += f" ({n_cand} candidates / {n_total} sampled)"
    else:
        suptitle += f" ({n_cand} candidates)"

    fig.suptitle(suptitle, fontsize=12, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.94])

    if save:
        outpath = os.path.join(indir, "candidates_plots.png")
        fig.savefig(outpath, dpi=150)
        print(f"Saved: {outpath}")
    if display:
        plt.show()
    else:
        plt.close(fig)

    # --- Statistics table ---
    if stats:
        plot_statistics_table(stats, indir, save, display)


def plot_statistics_table(stats: dict, indir: str, save: bool, display: bool):
    """Render a summary table of rejection fractions as a figure."""
    fraction_keys = [
        ("Accepted (SS)", "f_accepted"),
        ("Vetoed (total)", "f_vetoed"),
        ("  TPC veto", "f_vetoed_tpc"),
        ("  Skin veto", "f_vetoed_skin"),
        ("MS rejected", "f_ms_rejected"),
        ("No FV", "f_no_fv"),
    ]

    rows = []
    for label, key in fraction_keys:
        val = stats.get(key, None)
        if val is not None:
            rows.append([label, f"{float(val):.6e}"])
        else:
            rows.append([label, "—"])

    # SS and MS fractions among FV events
    f_acc = float(stats.get("f_accepted", 0))
    f_ms = float(stats.get("f_ms_rejected", 0))
    f_fv = f_acc + f_ms
    if f_fv > 0:
        rows.append(["", ""])
        rows.append(["SS fraction (in FV)", f"{f_acc / f_fv:.4f}"])
        rows.append(["MS fraction (in FV)", f"{f_ms / f_fv:.4f}"])

    fig, ax = plt.subplots(figsize=(5, 3))
    ax.axis("off")

    table = ax.table(
        cellText=rows,
        colLabels=["Category", "Fraction"],
        cellLoc="left",
        colLoc="left",
        loc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1.0, 1.6)

    # Style header
    for j in range(2):
        cell = table[0, j]
        cell.set_facecolor("#4472C4")
        cell.set_text_props(color="white", fontweight="bold")

    # Alternate row colors
    for i in range(len(rows)):
        for j in range(2):
            cell = table[i + 1, j]
            cell.set_facecolor("#D9E2F3" if i % 2 == 0 else "white")

    suptitle = make_suptitle(indir) + " — Rejection fractions"
    fig.suptitle(suptitle, fontsize=12, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.90])

    if save:
        outpath = os.path.join(indir, "statistics_table.png")
        fig.savefig(outpath, dpi=150, bbox_inches="tight")
        print(f"Saved: {outpath}")
    if display:
        plt.show()
    else:
        plt.close(fig)


if __name__ == "__main__":
    main()
