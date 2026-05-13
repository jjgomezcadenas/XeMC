"""PMT background summary following LZ paper grouping.

Rows = component types (TPC PMTs, TPC PMT bases, TPC PMT structures,
TPC PMT cables, Skin PMTs + bases). Columns = mass, activities,
background per isotope, total.

Reads summary.json from six analysis directories: pmt_top, pmt_bottom,
pmt_top_cables, pmt_bottom_cables, pmt_skin_upper_ring, and
pmt_skin_lower_ring. Tolerant of missing isotopes (e.g., Tl208 not yet
computed).

Usage:
    python py/pmt_background_summary.py \\
        --top results/bfv/pmt/top/Bi214/analysis \\
        --bottom results/bfv/pmt/bottom/Bi214/analysis \\
        --top-cables results/bfv/pmt/top_cables/Bi214/analysis \\
        --bottom-cables results/bfv/pmt/bottom_cables/Bi214/analysis \\
        --skin-upper results/bfv/pmt/skin_upper_ring/Bi214/analysis \\
        --skin-lower results/bfv/pmt/skin_lower_ring/Bi214/analysis \\
        -o results/bfv/pmt/Bi214_summary \\
        --display

    # With Tl208 (when available): add the corresponding *-tl flags
    # pointing at the Tl208 analysis directories for each of the six
    # sources above.
"""

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="PMT background summary (LZ paper grouping)."
    )
    parser.add_argument("--top", required=True, type=Path, help="pmt_top Bi214 analysis dir")
    parser.add_argument("--bottom", required=True, type=Path, help="pmt_bottom Bi214 analysis dir")
    parser.add_argument("--top-cables", required=True, type=Path,
                        help="pmt_top_cables Bi214 analysis dir")
    parser.add_argument("--bottom-cables", required=True, type=Path,
                        help="pmt_bottom_cables Bi214 analysis dir")
    parser.add_argument("--skin-upper", required=True, type=Path,
                        help="pmt_skin_upper_ring Bi214 analysis dir")
    parser.add_argument("--skin-lower", required=True, type=Path,
                        help="pmt_skin_lower_ring Bi214 analysis dir")
    parser.add_argument("--top-tl", type=Path, default=None, help="pmt_top Tl208 analysis dir")
    parser.add_argument("--bottom-tl", type=Path, default=None, help="pmt_bottom Tl208 analysis dir")
    parser.add_argument("--top-cables-tl", type=Path, default=None,
                        help="pmt_top_cables Tl208 analysis dir")
    parser.add_argument("--bottom-cables-tl", type=Path, default=None,
                        help="pmt_bottom_cables Tl208 analysis dir")
    parser.add_argument("--skin-upper-tl", type=Path, default=None,
                        help="pmt_skin_upper_ring Tl208 analysis dir")
    parser.add_argument("--skin-lower-tl", type=Path, default=None,
                        help="pmt_skin_lower_ring Tl208 analysis dir")
    parser.add_argument("-o", "--output", required=True, type=Path, help="Output directory")
    parser.add_argument("--display", action="store_true", help="Show plots interactively")
    return parser.parse_args()


def load_summary(path: Path) -> dict | None:
    p = path / "summary.json"
    if not p.exists():
        return None
    with open(p) as f:
        return json.load(f)


def component_rates(smry: dict) -> dict[str, float]:
    """Return {component_name: N_ss_roi_per_year}."""
    br = smry["background_rate"]
    comps = br.get("source_components", [])
    rate_total = br["rate_gammas_per_s"]
    n_total = br["n_ss_roi_per_year"]
    result = {}
    for c in comps:
        frac = c["rate_gammas_per_s"] / rate_total if rate_total > 0 else 0.0
        result[c["name"]] = frac * n_total
    return result


def component_info(smry: dict, isotope: str = "bi") -> dict[str, dict]:
    """Return {component_name: {mass_kg, activity, rate_gammas_per_s}}."""
    br = smry["background_rate"]
    result = {}
    for c in br.get("source_components", []):
        result[c["name"]] = {
            "mass_kg": c["mass_kg"],
            f"activity_{isotope}_mBq_per_kg": c["activity_mBq_per_kg"],
            f"rate_{isotope}_gammas_per_s": c["rate_gammas_per_s"],
        }
    return result


# =====================================================================
# LZ-style grouping
# =====================================================================
# Each group: (label, list of full component names to sum)
GROUPS = [
    ("TPC PMTs",           ["PMT_TOP_PMTs", "PMT_BOT_PMTs"]),
    ("TPC PMT bases",      ["PMT_TOP_bases", "PMT_BOT_bases"]),
    ("TPC PMT structures", ["PMT_TOP_structure", "PMT_BOT_structure"]),
    ("TPC PMT cables",     ["PMT_TOP_cables", "PMT_BOT_cables"]),
    ("Skin PMTs + bases",  ["PMT_SKIN_UPPER_RING", "PMT_SKIN_LOWER_RING", "PMT_BOT_R8778_dome"]),
]


def gather_all(bi_summaries: dict[str, dict | None],
               tl_summaries: dict[str, dict | None]):
    """Build per-group totals for mass, activity, and background."""
    bi_rates: dict[str, float] = {}
    tl_rates: dict[str, float] = {}
    info: dict[str, dict] = {}

    for smry in bi_summaries.values():
        if smry is not None:
            bi_rates.update(component_rates(smry))
            for name, ci in component_info(smry, "bi").items():
                info.setdefault(name, {}).update(ci)

    for smry in tl_summaries.values():
        if smry is not None:
            tl_rates.update(component_rates(smry))
            for name, ci in component_info(smry, "tl").items():
                info.setdefault(name, {}).update(ci)

    seconds_per_year = 3.15576e7

    rows = []
    for label, members in GROUPS:
        mass = sum(info.get(m, {}).get("mass_kg", 0.0) for m in members)
        # Representative specific activities (members in a group share the same)
        act_bi = 0.0
        act_tl = 0.0
        for m in members:
            if m in info:
                act_bi = info[m].get("activity_bi_mBq_per_kg", 0.0)
                act_tl = info[m].get("activity_tl_mBq_per_kg", 0.0)
                break
        # Total activity = mass * specific_activity (mBq)
        total_act_bi = mass * act_bi
        total_act_tl = mass * act_tl
        # Gamma rate (gammas/year) from metadata, summed over group members
        rate_bi = sum(info.get(m, {}).get("rate_bi_gammas_per_s", 0.0) for m in members)
        rate_tl = sum(info.get(m, {}).get("rate_tl_gammas_per_s", 0.0) for m in members)
        rate_bi_yr = rate_bi * seconds_per_year
        rate_tl_yr = rate_tl * seconds_per_year
        bg_bi = sum(bi_rates.get(m, 0.0) for m in members)
        bg_tl = sum(tl_rates.get(m, 0.0) for m in members)
        bg_total = bg_bi + bg_tl
        rows.append({
            "label": label,
            "mass_kg": mass,
            "act_bi": act_bi,
            "act_tl": act_tl,
            "total_act_bi_mBq": total_act_bi,
            "total_act_tl_mBq": total_act_tl,
            "rate_bi_yr": rate_bi_yr,
            "rate_tl_yr": rate_tl_yr,
            "bg_bi": bg_bi,
            "bg_tl": bg_tl,
            "bg_total": bg_total,
        })
    return rows


def write_text(path: Path, rows: list[dict], has_tl: bool) -> None:
    w = 140 if has_tl else 130
    lines = [
        "=" * w,
        "PMT Background Summary: N(SS in ROI) per year",
        "=" * w,
        "",
    ]

    if has_tl:
        hdr = (f"  {'Component':22s} {'Mass':>8s} {'Bi214':>10s} {'Tl208':>10s} "
               f"{'Tot Act':>10s} {'Rate Bi':>12s} {'BG Bi214':>10s} {'BG Tl208':>10s} {'Total':>10s}")
        sub = (f"  {'':22s} {'(kg)':>8s} {'(mBq/kg)':>10s} {'(mBq/kg)':>10s} "
               f"{'(mBq)':>10s} {'(g/yr)':>12s} {'(evt/yr)':>10s} {'(evt/yr)':>10s} {'(evt/yr)':>10s}")
    else:
        hdr = (f"  {'Component':22s} {'Mass':>8s} {'Bi214':>10s} {'Tot Act':>10s} "
               f"{'Rate Bi':>12s} {'BG Bi214':>10s}")
        sub = (f"  {'':22s} {'(kg)':>8s} {'(mBq/kg)':>10s} {'(mBq)':>10s} "
               f"{'(g/yr)':>12s} {'(evt/yr)':>10s}")

    lines += [hdr, sub, "-" * w]

    total_bi = 0.0
    total_tl = 0.0
    total_all = 0.0
    for r in rows:
        total_bi += r["bg_bi"]
        total_tl += r["bg_tl"]
        total_all += r["bg_total"]
        if has_tl:
            lines.append(
                f"  {r['label']:22s} {r['mass_kg']:8.2f} {r['act_bi']:10.2f} "
                f"{r['act_tl']:10.2f} {r['total_act_bi_mBq']:10.2f} "
                f"{r['rate_bi_yr']:12.2e} {r['bg_bi']:10.4f} {r['bg_tl']:10.4f} "
                f"{r['bg_total']:10.4f}"
            )
        else:
            lines.append(
                f"  {r['label']:22s} {r['mass_kg']:8.2f} {r['act_bi']:10.2f} "
                f"{r['total_act_bi_mBq']:10.2f} {r['rate_bi_yr']:12.2e} "
                f"{r['bg_bi']:10.4f}"
            )

    lines.append("-" * w)
    if has_tl:
        lines.append(
            f"  {'TOTAL':22s} {'':>8s} {'':>10s} {'':>10s} {'':>10s} {'':>12s} "
            f"{total_bi:10.4f} {total_tl:10.4f} {total_all:10.4f}"
        )
    else:
        lines.append(
            f"  {'TOTAL':22s} {'':>8s} {'':>10s} {'':>10s} {'':>12s} "
            f"{total_bi:10.4f}"
        )
    lines.append("=" * w)

    text = "\n".join(lines)
    print(text)
    path.write_text(text + "\n")
    print(f"Saved {path}")


def plot_table(rows: list[dict], has_tl: bool,
               outdir: Path, display: bool) -> None:
    def fmt(v: float) -> str:
        return f"{v:.4f}" if v > 0 else "-"

    def fmtm(v: float) -> str:
        return f"{v:.2f}" if v > 0 else "-"

    def fmta(v: float) -> str:
        return f"{v:.2f}" if v > 0 else "-"

    def fmtr(v: float) -> str:
        return f"{v:.2e}" if v > 0 else "-"

    if has_tl:
        col_labels = ["Component", "Mass\n(kg)", "Bi214\n(mBq/kg)", "Tl208\n(mBq/kg)",
                       "Tot Act\n(mBq)", "Rate Bi\n(g/yr)",
                       "BG Bi214\n(evt/yr)", "BG Tl208\n(evt/yr)", "Total\n(evt/yr)"]
    else:
        col_labels = ["Component", "Mass\n(kg)", "Bi214\n(mBq/kg)", "Tot Act\n(mBq)",
                       "Rate Bi\n(g/yr)", "BG Bi214\n(evt/yr)"]

    table_data = []
    total_bi = 0.0
    total_tl = 0.0
    total_all = 0.0
    for r in rows:
        total_bi += r["bg_bi"]
        total_tl += r["bg_tl"]
        total_all += r["bg_total"]
        if has_tl:
            table_data.append([
                r["label"], fmtm(r["mass_kg"]), fmta(r["act_bi"]), fmta(r["act_tl"]),
                fmtm(r["total_act_bi_mBq"]), fmtr(r["rate_bi_yr"]),
                fmt(r["bg_bi"]), fmt(r["bg_tl"]), fmt(r["bg_total"]),
            ])
        else:
            table_data.append([
                r["label"], fmtm(r["mass_kg"]), fmta(r["act_bi"]),
                fmtm(r["total_act_bi_mBq"]), fmtr(r["rate_bi_yr"]),
                fmt(r["bg_bi"]),
            ])

    # Total row
    if has_tl:
        table_data.append(["TOTAL", "", "", "", "", "", fmt(total_bi), fmt(total_tl), fmt(total_all)])
    else:
        table_data.append(["TOTAL", "", "", "", "", fmt(total_bi)])

    n_cols = len(col_labels)
    n_rows = len(table_data)
    fig_w = 14 if has_tl else 10
    fig, ax = plt.subplots(figsize=(fig_w, 0.5 * n_rows + 1.5))
    ax.axis("off")
    tbl = ax.table(cellText=table_data, colLabels=col_labels,
                   loc="center", cellLoc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(10)
    tbl.scale(1.0, 1.8)

    for (row, col), cell in tbl.get_celld().items():
        if row == 0:
            cell.set_facecolor("#4472C4")
            cell.set_text_props(color="white", weight="bold")
        elif col == 0:
            cell.set_text_props(weight="bold")
        elif col == n_cols - 1 and row > 0:
            cell.set_facecolor("#D9E2F3")
        if row == n_rows:
            cell.set_facecolor("#E0E0E0")
            cell.set_text_props(weight="bold")
        if row == n_rows and col == n_cols - 1:
            cell.set_facecolor("#C6EFCE")
            cell.set_text_props(weight="bold")

    ax.set_title("PMT Background: N(SS in ROI)/year", fontsize=12, pad=12)
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

    bi_summaries = {
        "top": load_summary(args.top),
        "bottom": load_summary(args.bottom),
        "top_cables": load_summary(args.top_cables),
        "bottom_cables": load_summary(args.bottom_cables),
        "skin_upper": load_summary(args.skin_upper),
        "skin_lower": load_summary(args.skin_lower),
    }

    tl_summaries = {
        "top": load_summary(args.top_tl) if args.top_tl else None,
        "bottom": load_summary(args.bottom_tl) if args.bottom_tl else None,
        "top_cables": load_summary(args.top_cables_tl) if args.top_cables_tl else None,
        "bottom_cables": load_summary(args.bottom_cables_tl) if args.bottom_cables_tl else None,
        "skin_upper": load_summary(args.skin_upper_tl) if args.skin_upper_tl else None,
        "skin_lower": load_summary(args.skin_lower_tl) if args.skin_lower_tl else None,
    }

    has_tl = any(v is not None for v in tl_summaries.values())

    rows = gather_all(bi_summaries, tl_summaries)

    write_text(args.output / "pmt_background_summary.txt", rows, has_tl)
    plot_table(rows, has_tl, args.output, args.display)

    print("Done.")


if __name__ == "__main__":
    main()
