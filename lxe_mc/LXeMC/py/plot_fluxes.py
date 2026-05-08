"""
Plot source flux tables as 2D colormaps.

Reads all .csv files in a directory (produced by generate_flux_tables.jl).
For Bi-214: one figure with 2x2 panels (OCV, MLI, ICV, rate).
For Tl-208: one figure per source (main + 3 companions), plus rate.

Usage:
    python py/plot_fluxes.py results/fluxes/cryostat/barrel/Bi214
    python py/plot_fluxes.py results/fluxes/cryostat/barrel/Tl208 --display
    python py/plot_fluxes.py results/fluxes/cryostat/barrel/Bi214 --no-save --display
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


def plot_flux(E, u, pdf, title: str, ax):
    """Plot a 2D flux table as a colormap."""
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


def title_from_filename(fname: str) -> str:
    """Derive a human-readable title from the CSV filename."""
    name = os.path.splitext(fname)[0]
    parts = name.replace("_", " ").split()
    return " ".join(p.capitalize() for p in parts)


def find_csvs(indir: str):
    """Find and classify CSV files in directory."""
    files = sorted(glob.glob(os.path.join(indir, "*.csv")))
    by_name = {os.path.basename(f): f for f in files}
    return by_name


def is_tl208(by_name: dict) -> bool:
    return any("companion" in k for k in by_name)


def plot_bi214(by_name: dict, indir: str, save: bool, display: bool):
    """Plot Bi-214 flux tables: 2x2 grid."""
    order = ["bi214_ocv", "bi214_mli", "bi214_icv", "bi214_rate"]
    panels = []
    for prefix in order:
        fname = prefix + ".csv"
        if fname in by_name:
            panels.append((fname, by_name[fname]))

    if not panels:
        # Fallback: plot whatever CSVs we have
        panels = [(k, v) for k, v in by_name.items()]

    n = len(panels)
    ncols = min(n, 2)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows),
                             squeeze=False)

    for idx, (fname, path) in enumerate(panels):
        row, col = divmod(idx, ncols)
        E, u, pdf = load_flux_csv(path)
        plot_flux(E, u, pdf, title_from_filename(fname), axes[row][col])

    for idx in range(n, nrows * ncols):
        row, col = divmod(idx, ncols)
        axes[row][col].set_visible(False)

    suptitle = _make_suptitle(indir)
    fig.suptitle(suptitle, fontsize=12, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.94])

    if save:
        outpath = os.path.join(indir, "bi214_fluxes.png")
        fig.savefig(outpath, dpi=150)
        print(f"Saved: {outpath}")
    if display:
        plt.show()
    else:
        plt.close(fig)


def plot_tl208(by_name: dict, indir: str, save: bool, display: bool):
    """Plot Tl-208 flux tables: one figure per source + rate."""
    sources = ["ocv", "mli", "icv"]
    companion_suffixes = ["583keV", "511keV", "861keV"]

    for src in sources:
        prefix = f"tl208_{src}"
        main_file = f"{prefix}_main.csv"
        comp_files = [f"{prefix}_companion_{c}.csv" for c in companion_suffixes]

        panels = []
        if main_file in by_name:
            panels.append((main_file, by_name[main_file]))
        for cf in comp_files:
            if cf in by_name:
                panels.append((cf, by_name[cf]))

        if not panels:
            continue

        n = len(panels)
        ncols = min(n, 4)
        nrows = (n + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows),
                                 squeeze=False)

        for idx, (fname, path) in enumerate(panels):
            row, col = divmod(idx, ncols)
            E, u, pdf = load_flux_csv(path)
            plot_flux(E, u, pdf, title_from_filename(fname), axes[row][col])

        for idx in range(n, nrows * ncols):
            row, col = divmod(idx, ncols)
            axes[row][col].set_visible(False)

        suptitle = f"{_make_suptitle(indir)} — {src.upper()}"
        fig.suptitle(suptitle, fontsize=12, fontweight="bold")
        fig.tight_layout(rect=[0, 0, 1, 0.94])

        if save:
            outpath = os.path.join(indir, f"tl208_{src}_fluxes.png")
            fig.savefig(outpath, dpi=150)
            print(f"Saved: {outpath}")
        if display:
            plt.show()
        else:
            plt.close(fig)

    # Rate table
    rate_file = "tl208_rate.csv"
    if rate_file in by_name:
        fig, ax = plt.subplots(1, 1, figsize=(5, 4))
        E, u, pdf = load_flux_csv(by_name[rate_file])
        plot_flux(E, u, pdf, title_from_filename(rate_file), ax)
        fig.suptitle(f"{_make_suptitle(indir)} — Rate", fontsize=12, fontweight="bold")
        fig.tight_layout(rect=[0, 0, 1, 0.94])

        if save:
            outpath = os.path.join(indir, "tl208_rate_flux.png")
            fig.savefig(outpath, dpi=150)
            print(f"Saved: {outpath}")
        if display:
            plt.show()
        else:
            plt.close(fig)


def _make_suptitle(indir: str) -> str:
    parts = os.path.normpath(indir).split(os.sep)
    # Try to extract meaningful parts (e.g., "barrel / Bi214")
    if len(parts) >= 2:
        return f"{parts[-2]} / {parts[-1]}"
    return parts[-1]


def main():
    if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
        print("Usage: python py/plot_fluxes.py <directory> [--display] [--no-save]")
        sys.exit(0)

    indir = sys.argv[1]
    display = "--display" in sys.argv
    save = "--no-save" not in sys.argv

    by_name = find_csvs(indir)
    if not by_name:
        print(f"No .csv files found in {indir}")
        sys.exit(1)

    if is_tl208(by_name):
        plot_tl208(by_name, indir, save, display)
    else:
        plot_bi214(by_name, indir, save, display)


if __name__ == "__main__":
    main()
