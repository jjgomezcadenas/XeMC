#!/usr/bin/env python3
"""
plot_source_flux.py — Plot source flux results from run_source_flux.jl.

Reads spectrum_E.csv, spectrum_u.csv, flux_table.csv from a results
directory and produces a 3-panel figure (dN/dE, dN/du, 2D heatmap).

Usage:
    python py/plot_source_flux.py /path/to/results/
    python py/plot_source_flux.py results/flux_test
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt


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

    # Flux table: first column is E, rest are u bins
    ft_raw = np.genfromtxt(
        os.path.join(results_dir, "flux_table.csv"),
        delimiter=",", skip_header=1
    )
    E_centers = ft_raw[:, 0]
    counts_2d = ft_raw[:, 1:]

    # Read summary for title
    title = ""
    summary_path = os.path.join(results_dir, "summary.txt")
    if os.path.exists(summary_path):
        with open(summary_path) as f:
            for line in f:
                if "Source:" in line and "Generated" not in line:
                    src = line.strip().split(":")[-1].strip()
                elif "Decay:" in line:
                    decay = line.strip().split(":")[-1].strip()
                elif "Surviving:" in line:
                    surv = line.strip()
            title = f"{src} — {decay}"

    # Figure
    fig, axs = plt.subplots(1, 3, figsize=(14, 4))

    # Panel 1: dN/dE
    ax = axs[0]
    E_vals = spec_E[:, 0]
    N_vals = spec_E[:, 1]
    dE = E_vals[1] - E_vals[0] if len(E_vals) > 1 else 0.01
    ax.bar(E_vals, N_vals, width=dE * 0.9, color="C0", edgecolor="C0")
    ax.set_xlabel("E [MeV]")
    ax.set_ylabel("counts")
    ax.set_title("dN/dE")
    ax.set_xlim(E_vals[0] - dE, E_vals[-1] + dE)

    # Panel 2: dN/du
    ax = axs[1]
    u_vals = spec_u[:, 0]
    N_u = spec_u[:, 1]
    du = u_vals[1] - u_vals[0] if len(u_vals) > 1 else 0.1
    ax.bar(u_vals, N_u, width=du * 0.9, color="C1", edgecolor="C1")
    ax.set_xlabel("cos θ")
    ax.set_ylabel("counts")
    ax.set_title("dN/d(cos θ)")
    ax.set_xlim(0, 1)

    # Panel 3: 2D heatmap
    ax = axs[2]
    n_u = counts_2d.shape[1]
    u_edges = np.linspace(0, 1, n_u + 1)
    E_edges = np.zeros(len(E_centers) + 1)
    E_edges[:-1] = E_centers - dE / 2
    E_edges[-1] = E_centers[-1] + dE / 2
    im = ax.pcolormesh(u_edges, E_edges, counts_2d, shading="flat", cmap="hot")
    ax.set_xlabel("cos θ")
    ax.set_ylabel("E [MeV]")
    ax.set_title("flux (E, cos θ)")
    fig.colorbar(im, ax=ax, label="counts")

    if title:
        fig.suptitle(title, fontsize=12, y=1.02)

    fig.tight_layout()

    # Save
    out_path = os.path.join(results_dir, "flux_plot.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved: {out_path}")

    # Show
    plt.show()


if __name__ == "__main__":
    main()
