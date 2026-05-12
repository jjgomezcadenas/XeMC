"""PMT background breakdown by source location (top/bottom) and component
(PMTs, bases, structure, R8778_dome).

Reads summary.json from the analysis directories of pmt_top and pmt_bottom,
decomposes the total background rate per component, and produces a cross-table
with row and column totals.

Usage:
    python py/pmt_background_summary.py \\
        --top results/bfv/pmt/top/Bi214/analysis \\
        --bottom results/bfv/pmt/bottom/Bi214/analysis \\
        -o results/bfv/pmt/Bi214_summary \\
        --display
"""

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="PMT background breakdown by location and component."
    )
    parser.add_argument(
        "--top", required=True, type=Path,
        help="Analysis directory for pmt_top (containing summary.json)"
    )
    parser.add_argument(
        "--bottom", required=True, type=Path,
        help="Analysis directory for pmt_bottom (containing summary.json)"
    )
    parser.add_argument(
        "-o", "--output", required=True, type=Path,
        help="Output directory for summary files"
    )
    parser.add_argument(
        "--display", action="store_true", help="Show plots interactively"
    )
    return parser.parse_args()


def load_summary(path: Path) -> dict:
    with open(path / "summary.json") as f:
        return json.load(f)


def component_breakdown(smry: dict) -> dict[str, float]:
    """Return {component_name: N_ss_roi_per_year} for each component."""
    br = smry["background_rate"]
    components = br.get("source_components", [])
    rate_total = br["rate_gammas_per_s"]
    n_ss_roi_total = br["n_ss_roi_per_year"]

    result = {}
    for comp in components:
        frac = comp["rate_gammas_per_s"] / rate_total if rate_total > 0 else 0.0
        result[comp["name"]] = frac * n_ss_roi_total
    return result


# Canonical component order and short names
TOP_COMPONENTS = [
    ("PMT_TOP_PMTs", "PMTs"),
    ("PMT_TOP_bases", "Bases"),
    ("PMT_TOP_structure", "Structure"),
]

BOT_COMPONENTS = [
    ("PMT_BOT_PMTs", "PMTs"),
    ("PMT_BOT_bases", "Bases"),
    ("PMT_BOT_structure", "Structure"),
    ("PMT_BOT_R8778_dome", "R8778 dome"),
]

# Unified column labels (superset)
ALL_SHORT = ["PMTs", "Bases", "Structure", "R8778 dome"]


def build_table(top_bd: dict[str, float], bot_bd: dict[str, float]):
    """Build the cross-table as a 2D array + row/col totals."""
    top_row = []
    for full, short in TOP_COMPONENTS:
        top_row.append(top_bd.get(full, 0.0))
    # Top has no R8778_dome
    top_row.append(0.0)

    bot_row = []
    for full, short in BOT_COMPONENTS:
        bot_row.append(bot_bd.get(full, 0.0))

    top_arr = np.array(top_row)
    bot_arr = np.array(bot_row)
    col_totals = top_arr + bot_arr
    top_total = top_arr.sum()
    bot_total = bot_arr.sum()
    grand_total = top_total + bot_total

    return top_arr, bot_arr, col_totals, top_total, bot_total, grand_total


def write_text(path: Path, top_smry: dict, bot_smry: dict,
               top_arr, bot_arr, col_totals,
               top_total, bot_total, grand_total) -> None:
    top_br = top_smry["background_rate"]
    bot_br = bot_smry["background_rate"]

    lines = [
        "=" * 80,
        "PMT Background Summary: N(SS in ROI) per year by component",
        "=" * 80,
        f"  Top source:    {top_br['source']}  isotope: {top_br['isotope']}",
        f"  Bottom source: {bot_br['source']}  isotope: {bot_br['isotope']}",
        "-" * 80,
        "",
        f"  {'':15s} {'PMTs':>12s} {'Bases':>12s} {'Structure':>12s} {'R8778 dome':>12s} {'Total':>12s}",
        f"  {'Top':15s} {top_arr[0]:12.4f} {top_arr[1]:12.4f} {top_arr[2]:12.4f} {top_arr[3]:12.4f} {top_total:12.4f}",
        f"  {'Bottom':15s} {bot_arr[0]:12.4f} {bot_arr[1]:12.4f} {bot_arr[2]:12.4f} {bot_arr[3]:12.4f} {bot_total:12.4f}",
        f"  {'Total':15s} {col_totals[0]:12.4f} {col_totals[1]:12.4f} {col_totals[2]:12.4f} {col_totals[3]:12.4f} {grand_total:12.4f}",
        "",
        "-" * 80,
        "  Per-component details (Top):",
    ]
    for comp in top_br.get("source_components", []):
        lines.append(f"    {comp['name']:25s}  mass={comp['mass_kg']:.3f} kg  "
                     f"act={comp['activity_mBq_per_kg']:.4f} mBq/kg  "
                     f"geom={comp['geometric_survival']:.6f}  "
                     f"rate={comp['rate_gammas_per_s']:.4e} g/s")

    lines.append("  Per-component details (Bottom):")
    for comp in bot_br.get("source_components", []):
        lines.append(f"    {comp['name']:25s}  mass={comp['mass_kg']:.3f} kg  "
                     f"act={comp['activity_mBq_per_kg']:.4f} mBq/kg  "
                     f"geom={comp['geometric_survival']:.6f}  "
                     f"rate={comp['rate_gammas_per_s']:.4e} g/s")

    lines += [
        "",
        f"  Top:    N(SS in ROI)/year = {top_br['n_ss_roi_per_year']:.4f}",
        f"  Bottom: N(SS in ROI)/year = {bot_br['n_ss_roi_per_year']:.4f}",
        f"  Total:  N(SS in ROI)/year = {grand_total:.4f}",
        "=" * 80,
    ]

    text = "\n".join(lines)
    print(text)
    path.write_text(text + "\n")
    print(f"Saved {path}")


def plot_table(top_arr, bot_arr, col_totals,
               top_total, bot_total, grand_total,
               outdir: Path, display: bool) -> None:
    def fmt(v: float) -> str:
        return f"{v:.4f}" if v > 0 else "-"

    table_data = [
        ["Top"] + [fmt(v) for v in top_arr] + [fmt(top_total)],
        ["Bottom"] + [fmt(v) for v in bot_arr] + [fmt(bot_total)],
        ["Total"] + [fmt(v) for v in col_totals] + [fmt(grand_total)],
    ]
    col_labels = ["", "PMTs", "Bases", "Structure", "R8778 dome", "Total"]

    fig, ax = plt.subplots(figsize=(10, 2.5))
    ax.axis("off")
    tbl = ax.table(cellText=table_data, colLabels=col_labels,
                   loc="center", cellLoc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(11)
    tbl.scale(1.0, 1.8)

    for (row, col), cell in tbl.get_celld().items():
        if row == 0:
            cell.set_facecolor("#4472C4")
            cell.set_text_props(color="white", weight="bold")
        elif col == 0:
            cell.set_text_props(weight="bold")
        elif col == 5:
            cell.set_facecolor("#D9E2F3")
        if row == 3:
            cell.set_facecolor("#E0E0E0")
            cell.set_text_props(weight="bold")
        if row == 3 and col == 5:
            cell.set_facecolor("#C6EFCE")
            cell.set_text_props(weight="bold")

    ax.set_title("PMT Background: N(SS in ROI)/year by component", fontsize=12, pad=12)
    fig.tight_layout()
    outfile = outdir / "pmt_background_table.png"
    fig.savefig(outfile, dpi=150, bbox_inches="tight")
    print(f"Saved {outfile}")
    if display:
        plt.show()
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    print(f"Reading top:    {args.top}")
    print(f"Reading bottom: {args.bottom}")
    top_smry = load_summary(args.top)
    bot_smry = load_summary(args.bottom)

    top_bd = component_breakdown(top_smry)
    bot_bd = component_breakdown(bot_smry)

    top_arr, bot_arr, col_totals, top_total, bot_total, grand_total = \
        build_table(top_bd, bot_bd)

    write_text(args.output / "pmt_background_summary.txt",
               top_smry, bot_smry,
               top_arr, bot_arr, col_totals,
               top_total, bot_total, grand_total)

    plot_table(top_arr, bot_arr, col_totals,
               top_total, bot_total, grand_total,
               args.output, args.display)

    print("Done.")


if __name__ == "__main__":
    main()
