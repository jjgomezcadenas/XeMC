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


def read_energy_status_table(path: Path):
    rows = []
    with path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(
                (
                    float(row["E1_MeV"]),
                    float(row["E2_MeV"]),
                    row["status"],
                    int(row["count"]),
                )
            )
    return rows


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


def plot_fraction_bars(ax, labels, counts, title, color):
    total = float(np.sum(counts))
    fractions = counts / total if total > 0 else np.zeros_like(counts, dtype=float)
    x = np.arange(len(labels))
    ax.bar(x, fractions, color=color, edgecolor="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=20, ha="right")
    ax.set_ylabel("fraction")
    ax.set_title(title)


def make_energy_status_table(rows):
    e1_vals = sorted({r[0] for r in rows})
    e2_vals = sorted({r[1] for r in rows})
    statuses = ["accepted", "vetoed", "no_fv"]

    e1_index = {v: i for i, v in enumerate(e1_vals)}
    e2_index = {v: i for i, v in enumerate(e2_vals)}
    counts = np.zeros((len(e1_vals), len(e2_vals), len(statuses)), dtype=int)

    for e1, e2, status, count in rows:
        if status not in statuses:
            continue
        counts[e1_index[e1], e2_index[e2], statuses.index(status)] += count

    return e1_vals, e2_vals, statuses, counts


def format_energy_label(e):
    return "0" if abs(e) < 1e-12 else f"{e:.2f}"


def make_figure(directory: Path):
    status_labels, status_counts = read_counts_csv(directory / "event_status_counts.csv")
    mult_labels, mult_counts = read_counts_csv(directory / "event_multiplicity_counts.csv")
    decisive_labels, decisive_counts = read_counts_csv(directory / "first_decisive_gamma_counts.csv")
    rows = read_energy_status_table(directory / "event_energy_status_table.csv")
    summary = read_summary(directory / "summary.txt")

    e1_vals, e2_vals, statuses, table_counts = make_energy_status_table(rows)
    total = np.sum(table_counts)
    frac_counts = table_counts / total if total > 0 else table_counts.astype(float)

    fig = plt.figure(figsize=(16, 10), constrained_layout=True)
    gs = fig.add_gridspec(2, 3, height_ratios=[0.9, 1.3])

    ax_status = fig.add_subplot(gs[0, 0])
    ax_mult = fig.add_subplot(gs[0, 1])
    ax_decisive = fig.add_subplot(gs[0, 2])
    ax_table = fig.add_subplot(gs[1, :])

    plot_fraction_bars(ax_status, status_labels, status_counts, "Event Status", "#4C78A8")
    plot_fraction_bars(ax_mult, mult_labels, mult_counts, "Multiplicity", "#F58518")
    plot_fraction_bars(ax_decisive, decisive_labels, decisive_counts, "First Decisive Gamma", "#54A24B")

    accepted = frac_counts[:, :, statuses.index("accepted")]
    vetoed = frac_counts[:, :, statuses.index("vetoed")]
    no_fv = frac_counts[:, :, statuses.index("no_fv")]

    x = np.arange(len(e1_vals))
    width = 0.8 / max(1, len(e2_vals))
    colors = {
        "accepted": "#4C78A8",
        "vetoed": "#E45756",
        "no_fv": "#72B7B2",
    }

    for j, e2 in enumerate(e2_vals):
        xj = x - 0.4 + width * (j + 0.5)
        base = np.zeros(len(e1_vals))
        ax_table.bar(
            xj, accepted[:, j], width=width, bottom=base,
            color=colors["accepted"], edgecolor="black", linewidth=0.4,
            label="accepted" if j == 0 else None,
        )
        base += accepted[:, j]
        ax_table.bar(
            xj, vetoed[:, j], width=width, bottom=base,
            color=colors["vetoed"], edgecolor="black", linewidth=0.4,
            label="vetoed" if j == 0 else None,
        )
        base += vetoed[:, j]
        ax_table.bar(
            xj, no_fv[:, j], width=width, bottom=base,
            color=colors["no_fv"], edgecolor="black", linewidth=0.4,
            label="no_fv" if j == 0 else None,
        )

    ax_table.set_xticks(x)
    ax_table.set_xticklabels([format_energy_label(e1) for e1 in e1_vals])
    ax_table.set_xlabel("E1 [MeV]")
    ax_table.set_ylabel("fraction")
    ax_table.set_title("Event Outcome by (E1, E2)")
    ax_table.legend(loc="upper right")

    e2_note = "  ".join([f"E2={format_energy_label(e)}" for e in e2_vals])
    ax_table.text(
        0.01, 0.98, e2_note,
        transform=ax_table.transAxes,
        ha="left", va="top",
        fontsize=10,
        bbox=dict(facecolor="white", alpha=0.85, edgecolor="none"),
    )

    title_parts = ["Event Processing Histograms"]
    if "path" in summary:
        title_parts.append(f"path={summary['path']}")
    if "N" in summary:
        title_parts.append(f"N={summary['N']}")
    if "energy_MeV" in summary:
        title_parts.append(f"E={float(summary['energy_MeV']):.2f} MeV")
    fig.suptitle(" | ".join(title_parts), fontsize=16)
    return fig


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot event-processing histogram CSVs from one output directory."
    )
    parser.add_argument(
        "directory",
        help="Directory containing summary.txt and the event-processing CSV files.",
    )
    parser.add_argument(
        "--outfile",
        default="event_histos.png",
        help="PNG filename to write inside the directory (default: event_histos.png).",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    directory = Path(args.directory).expanduser().resolve()
    if not directory.is_dir():
        raise SystemExit(f"Not a directory: {directory}")

    required = [
        "event_status_counts.csv",
        "event_multiplicity_counts.csv",
        "first_decisive_gamma_counts.csv",
        "event_energy_status_table.csv",
    ]
    missing = [name for name in required if not (directory / name).exists()]
    if missing:
        raise SystemExit(f"Missing required files in {directory}: {', '.join(missing)}")

    fig = make_figure(directory)
    outpath = directory / args.outfile
    fig.savefig(outpath, dpi=180)
    print(f"Wrote {outpath}")
    plt.show()
    plt.close(fig)


if __name__ == "__main__":
    main()
