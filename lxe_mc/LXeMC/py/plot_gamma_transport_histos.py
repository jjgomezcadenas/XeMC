#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def read_counts_csv(path: Path):
    labels = []
    counts = []
    with path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            labels.append(row["key"])
            counts.append(int(row["count"]))
    return labels, np.asarray(counts, dtype=int)


def read_hist1d_csv(path: Path):
    lo = []
    hi = []
    counts = []
    with path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            lo.append(float(row["bin_lo"]))
            hi.append(float(row["bin_hi"]))
            counts.append(int(row["count"]))
    lo = np.asarray(lo, dtype=float)
    hi = np.asarray(hi, dtype=float)
    counts = np.asarray(counts, dtype=int)
    centers = 0.5 * (lo + hi)
    widths = hi - lo
    return lo, hi, centers, widths, counts


def read_heatmap_csv(path: Path):
    rows = []
    with path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(
                (
                    float(row["z_bin_lo"]),
                    float(row["z_bin_hi"]),
                    float(row["r_bin_lo"]),
                    float(row["r_bin_hi"]),
                    int(row["count"]),
                )
            )

    z_lo = sorted({r[0] for r in rows})
    z_hi = sorted({r[1] for r in rows})
    r_lo = sorted({r[2] for r in rows})
    r_hi = sorted({r[3] for r in rows})

    z_edges = np.asarray(z_lo + [max(z_hi)], dtype=float)
    r_edges = np.asarray(r_lo + [max(r_hi)], dtype=float)
    counts = np.zeros((len(z_lo), len(r_lo)), dtype=int)

    z_index = {v: i for i, v in enumerate(z_lo)}
    r_index = {v: i for i, v in enumerate(r_lo)}
    for zl, _, rl, _, c in rows:
        counts[z_index[zl], r_index[rl]] = c

    return r_edges, z_edges, counts


def read_summary(path: Path):
    data = {}
    if not path.exists():
        return data
    with path.open() as f:
        for line in f:
            if "=" not in line:
                continue
            key, value = line.split("=", 1)
            data[key.strip()] = value.strip()
    return data


def plot_counts(ax, labels, counts, title, color):
    x = np.arange(len(labels))
    ax.bar(x, counts, color=color, edgecolor="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=25, ha="right")
    ax.set_title(title)
    ax.set_ylabel("fraction")


def plot_hist(ax, centers, widths, counts, title, xlabel, color):
    ax.bar(centers, counts, width=widths, align="center", color=color, edgecolor="black", linewidth=0.6)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("count")


def make_figure(directory: Path):
    status_labels, status_counts = read_counts_csv(directory / "status_counts.csv")
    interaction_labels, interaction_counts = read_counts_csv(directory / "interaction_type_counts.csv")
    _, _, e_centers, e_widths, e_counts = read_hist1d_csv(directory / "deposit_energy_hist.csv")
    r_edges, z_edges, heatmap = read_heatmap_csv(directory / "interaction_rz_heatmap.csv")
    summary = read_summary(directory / "summary.txt")

    status_map = dict(zip(status_labels, status_counts))
    interaction_map = dict(zip(interaction_labels, interaction_counts))
    total = float(np.sum(status_counts))

    outcome_labels = ["compton", "pair", "photo", "inFV", "out"]
    outcome_counts = np.asarray(
        [
            interaction_map.get("compton", 0),
            interaction_map.get("pair", 0),
            interaction_map.get("photoelectric", 0),
            status_map.get("entered_fv", 0),
            status_map.get("escaped", 0),
        ],
        dtype=float,
    )
    if total > 0:
        outcome_counts /= total

    fig = plt.figure(figsize=(14, 10), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[0.9, 1.35])

    ax_outcome = fig.add_subplot(gs[0, 0])
    ax_e = fig.add_subplot(gs[0, 1])
    ax_heat = fig.add_subplot(gs[1, :])

    plot_counts(ax_outcome, outcome_labels, outcome_counts, "Interaction Outcome", "#F58518")

    title_bits = []
    for key in ["E0_MeV", "x0_cm", "y0_cm", "z0_cm", "ux", "uy", "uz", "N"]:
        if key in summary:
            title_bits.append(f"{key}={summary[key]}")
    if title_bits:
        ax_outcome.set_title("Interaction Outcome\n" + ", ".join(title_bits), fontsize=10)

    plot_hist(ax_e, e_centers, e_widths, e_counts, "Deposited Energy", "MeV", "#54A24B")

    mesh = ax_heat.pcolormesh(r_edges, z_edges, heatmap, shading="auto", cmap="viridis")
    ax_heat.set_title("Interaction Occupancy: r vs z")
    ax_heat.set_xlabel("r [cm]")
    ax_heat.set_ylabel("z [cm]")
    ax_heat.set_aspect("auto")
    fig.colorbar(mesh, ax=ax_heat, label="count")

    fig.suptitle("Gamma Transport Histograms", fontsize=16)
    return fig


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot gamma transport histogram CSVs from one output directory."
    )
    parser.add_argument(
        "directory",
        help="Directory containing summary.txt and the histogram CSV files.",
    )
    parser.add_argument(
        "--outfile",
        default="histos.png",
        help="PNG filename to write inside the directory (default: histos.png).",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    directory = Path(args.directory).expanduser().resolve()
    if not directory.is_dir():
        raise SystemExit(f"Not a directory: {directory}")

    fig = make_figure(directory)
    outpath = directory / args.outfile
    fig.savefig(outpath, dpi=180)
    print(f"Wrote {outpath}")
    plt.show(block=False)
    input("Press Return to close the figure...")
    plt.close(fig)


if __name__ == "__main__":
    main()
