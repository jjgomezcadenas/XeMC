#!/usr/bin/env python3
"""
plot_source_flux.py — Plot source flux results from run_source_flux.jl.

Reads spectrum_E.csv, spectrum_u.csv, flux_table.csv, summary.txt
from a results directory. Produces a 4-panel figure:
  1. dN/dE (log y)
  2. dN/d(cos θ)
  3. 2D heatmap E vs cos θ (log color)
  4. Event fractions (surviving, vetoed, absorbed, low-E, backward)

Usage:
    python py/plot_source_flux.py /path/to/results/
"""

import sys
import os
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


def parse_summary(path):
    """Extract counts from summary.txt."""
    info = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if "Generated:" in line:
                m = re.search(r"(\d+)", line.split("Generated:")[-1])
                if m:
                    info["generated"] = int(m.group(1))
            elif "Surviving:" in line:
                m = re.search(r"(\d+)", line.split("Surviving:")[-1])
                if m:
                    info["surviving"] = int(m.group(1))
            elif "Vetoed:" in line:
                m = re.search(r"(\d+)", line.split("Vetoed:")[-1])
                if m:
                    info["vetoed"] = int(m.group(1))
            elif "Absorbed:" in line:
                m = re.search(r"(\d+)", line.split("Absorbed:")[-1])
                if m:
                    info["absorbed"] = int(m.group(1))
            elif "Low energy:" in line:
                m = re.search(r"(\d+)", line.split("Low energy:")[-1])
                if m:
                    info["low_energy"] = int(m.group(1))
            elif "Invisible:" in line:
                m = re.search(r"(\d+)", line.split("Invisible:")[-1])
                if m:
                    info["invisible"] = int(m.group(1))
            elif "Backward:" in line:
                m = re.search(r"(\d+)", line.split("Backward:")[-1])
                if m:
                    info["backward"] = int(m.group(1))
            elif "Peak bin" in line:
                m = re.search(r"([\d.]+e[+-]?\d+)", line)
                if m:
                    info["peak_frac"] = float(m.group(1))
            elif "Off-peak" in line:
                m = re.search(r"([\d.]+e[+-]?\d+)", line)
                if m:
                    info["off_peak_frac"] = float(m.group(1))
            elif "Source:" in line and "generated" not in line.lower():
                info["source"] = line.split(":")[-1].strip()
            elif "Decay:" in line:
                info["decay"] = line.split(":")[-1].strip()
    # All counters are now per-event
    return info


def main():
    if len(sys.argv) < 2:
        print("Usage: python plot_source_flux.py <results_dir>")
        sys.exit(1)

    results_dir = sys.argv[1]
    if not os.path.isdir(results_dir):
        print(f"Error: directory '{results_dir}' not found")
        sys.exit(1)

    # Load data
    spec_E = np.genfromtxt(
        os.path.join(results_dir, "spectrum_E.csv"),
        delimiter=",", skip_header=1
    )
    spec_u = np.genfromtxt(
        os.path.join(results_dir, "spectrum_u.csv"),
        delimiter=",", skip_header=1
    )
    ft_raw = np.genfromtxt(
        os.path.join(results_dir, "flux_table.csv"),
        delimiter=",", skip_header=1
    )
    E_centers = ft_raw[:, 0]
    counts_2d = ft_raw[:, 1:]
    dE = E_centers[1] - E_centers[0] if len(E_centers) > 1 else 0.01

    # Parse summary
    summary_path = os.path.join(results_dir, "summary.txt")
    info = parse_summary(summary_path) if os.path.exists(summary_path) else {}

    title = f"{info.get('source', '?')} — {info.get('decay', '?')}"

    # Figure: 2x2
    fig, axs = plt.subplots(2, 2, figsize=(12, 9))

    # Panel 1: dN/dE (log y)
    ax = axs[0, 0]
    E_vals = spec_E[:, 0]
    N_vals = spec_E[:, 1]
    ax.bar(E_vals, N_vals, width=dE * 0.9, color="C0", edgecolor="C0")
    ax.set_xlabel("E [MeV]")
    ax.set_ylabel("prob / decay / bin")
    ax.set_yscale("log")
    ax.set_title("dP/dE")
    ax.set_xlim(E_vals[0] - dE, E_vals[-1] + dE)

    # Panel 2: dN/du
    ax = axs[0, 1]
    u_vals = spec_u[:, 0]
    N_u = spec_u[:, 1]
    du = u_vals[1] - u_vals[0] if len(u_vals) > 1 else 0.1
    ax.bar(u_vals, N_u, width=du * 0.9, color="C1", edgecolor="C1")
    ax.set_xlabel("cos θ")
    ax.set_ylabel("prob / decay / bin")
    ax.set_title("dP/d(cos θ)")
    ax.set_xlim(0, 1)

    # Panel 3: 2D heatmap (log color), peak bin excluded
    ax = axs[1, 0]
    n_u = counts_2d.shape[1]
    u_edges = np.linspace(0, 1, n_u + 1)
    E_edges = np.zeros(len(E_centers) + 1)
    E_edges[:-1] = E_centers - dE / 2
    E_edges[-1] = E_centers[-1] + dE / 2
    # Exclude peak bin (last bin for Tl-208, bin containing the line for Bi-214)
    # Find peak bin: the one with highest total probability
    row_sums = counts_2d.sum(axis=1)
    peak_idx = np.argmax(row_sums)
    plot_data = counts_2d.astype(float).copy()
    plot_data[peak_idx, :] = 0.0  # zero out peak bin
    plot_data[plot_data == 0] = np.nan
    if np.any(~np.isnan(plot_data)):
        vmin = np.nanmin(plot_data)
        vmax = np.nanmax(plot_data)
        im = ax.pcolormesh(u_edges, E_edges, plot_data, shading="flat",
                           cmap="hot", norm=LogNorm(vmin=vmin, vmax=vmax))
        fig.colorbar(im, ax=ax, label="prob/decay/bin")
    else:
        ax.text(0.5, 0.5, "no off-peak data", ha="center", va="center",
                transform=ax.transAxes)
    ax.set_xlabel("cos θ")
    ax.set_ylabel("E [MeV]")
    ax.set_title("off-peak flux (peak bin excluded) [log]")

    # Panel 4: Event fractions + peak/off-peak (horizontal bar)
    ax = axs[1, 1]
    N_gen = info.get("generated", 1)
    categories = [
        "peak (unscattered)",
        "off-peak (scattered)",
        "vetoed (companion)",
        "low energy",
        "backward",
        "absorbed/invisible",
    ]
    # Fractions: peak/off-peak are probabilities per decay (already fractional)
    # Others are integer counts → divide by N_gen
    fracs = [
        100.0 * info.get("peak_frac", 0),
        100.0 * info.get("off_peak_frac", 0),
        100.0 * info.get("vetoed", 0) / N_gen,
        100.0 * info.get("low_energy", 0) / N_gen,
        100.0 * info.get("backward", 0) / N_gen,
        100.0 * (info.get("absorbed", 0) + info.get("invisible", 0)) / N_gen,
    ]
    colors = ["C2", "C0", "C3", "C5", "C8", "C7"]
    bars = ax.barh(categories, fracs, color=colors, edgecolor="k", linewidth=0.5)
    ax.set_xlabel("fraction of generated [%]")
    ax.set_title("event fate")
    ax.set_xlim(0, max(fracs) * 1.2 if max(fracs) > 0 else 100)
    for bar, frac in zip(bars, fracs):
        if frac > 0.3:
            ax.text(bar.get_width() + 0.3, bar.get_y() + bar.get_height()/2,
                    f"{frac:.1f}%", va="center", fontsize=8)

    fig.suptitle(title, fontsize=13)
    fig.tight_layout()

    # Save
    out_path = os.path.join(results_dir, "flux_plot.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved: {out_path}")

    plt.show()


if __name__ == "__main__":
    main()
