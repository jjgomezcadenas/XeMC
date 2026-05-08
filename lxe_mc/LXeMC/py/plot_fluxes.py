"""
Plot source flux tables as 2D colormaps.

Reads all .csv files in a directory (produced by generate_flux_tables.jl),
generates one colormap plot per file. Titles are derived from filenames.

Usage:
    python py/plot_fluxes.py results/fluxes/cryostat/barrel/Bi214
    python py/plot_fluxes.py results/fluxes/cryostat/barrel/Bi214 --save
"""

import sys
import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


def load_flux_csv(path: str):
    """Load a flux CSV file. Returns E centers, u centers, and the PDF matrix."""
    df = pd.read_csv(path, index_col=0)
    E_centers = df.index.values.astype(float)
    u_centers = np.array([float(c.split("=")[1]) for c in df.columns])
    pdf = df.values.astype(float)
    return E_centers, u_centers, pdf


def title_from_filename(fname: str) -> str:
    """Derive a human-readable title from the CSV filename."""
    name = os.path.splitext(fname)[0]
    # Replace underscores with spaces, capitalize words
    parts = name.replace("_", " ").split()
    return " ".join(p.capitalize() for p in parts)


def plot_flux(E, u, pdf, title: str, ax):
    """Plot a 2D flux table as a colormap."""
    # Mask zeros for log scale
    pdf_plot = np.where(pdf > 0, pdf, np.nan)
    vmin = np.nanmin(pdf_plot) if np.any(~np.isnan(pdf_plot)) else 1e-10
    vmax = np.nanmax(pdf_plot) if np.any(~np.isnan(pdf_plot)) else 1.0

    if vmin <= 0 or np.isnan(vmin):
        vmin = 1e-10
    if vmax <= vmin:
        vmax = vmin * 10

    im = ax.pcolormesh(u, E, pdf_plot, shading="nearest",
                       norm=LogNorm(vmin=vmin, vmax=vmax),
                       cmap="viridis")
    ax.set_xlabel(r"$u = \cos\theta$")
    ax.set_ylabel("E [MeV]")
    ax.set_title(title, fontsize=10)
    plt.colorbar(im, ax=ax, label="P / (decay · bin)")


def main():
    if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
        print("Usage: python py/plot_fluxes.py <directory> [--save]")
        sys.exit(0)

    indir = sys.argv[1]
    save = "--save" in sys.argv

    csv_files = sorted(glob.glob(os.path.join(indir, "*.csv")))
    if not csv_files:
        print(f"No .csv files found in {indir}")
        sys.exit(1)

    n = len(csv_files)
    ncols = min(n, 3)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows),
                             squeeze=False)

    for idx, path in enumerate(csv_files):
        row, col = divmod(idx, ncols)
        ax = axes[row][col]
        fname = os.path.basename(path)
        title = title_from_filename(fname)

        try:
            E, u, pdf = load_flux_csv(path)
            plot_flux(E, u, pdf, title, ax)
        except Exception as e:
            ax.text(0.5, 0.5, f"Error:\n{e}", transform=ax.transAxes,
                    ha="center", va="center", fontsize=8)
            ax.set_title(title, fontsize=10)

    # Hide unused axes
    for idx in range(n, nrows * ncols):
        row, col = divmod(idx, ncols)
        axes[row][col].set_visible(False)

    # Add directory as suptitle
    fig.suptitle(os.path.basename(os.path.normpath(indir)) + " — " +
                 os.path.basename(os.path.dirname(os.path.normpath(indir))),
                 fontsize=12, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    if save:
        outpath = os.path.join(indir, "flux_plots.pdf")
        fig.savefig(outpath, dpi=150)
        print(f"Saved to {outpath}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
