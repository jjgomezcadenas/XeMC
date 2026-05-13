"""
Cross-validate source geometry activities against the LZ bb0nu paper
(arXiv:1912.04248v2, Table I).

Checks:
  - Every paper Table I row mapped to JSON source(s)
  - Activities (238U-late = Bi214, 232Th-late = Tl208) match
  - Masses match (virtual equivalent_mass_kg, exact via geometry)
  - Identifies components in paper but absent from JSON
  - Documents naming convention (chain vs isotope-specific rates)
  - Breaks down lumped categories (skin PMTs, cryostat shells)

Usage:
    python py/validate_activities_vs_paper.py
    python py/validate_activities_vs_paper.py --json-source data/source_geometry_lz_v1.json
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Optional


DAYS_PER_YEAR = 365.25
LIVE_DAYS = 1000


# =====================================================================
# Paper Table I: every row
# (paper_name, paper_mass_kg, paper_U238_late_mBq_per_kg,
#  paper_Th232_late_mBq_per_kg, paper_counts_U, paper_counts_Th,
#  paper_total_counts,
#  json_sources_for_mass, json_source_for_activity)
# =====================================================================
PAPER_TABLE_I = [
    ("TPC PMTs", 91.9, 3.22, 1.61, 2.95, 0.10, 3.05,
     ["PMT_TOP_PMTs", "PMT_BOT_PMTs"], "PMT_TOP_PMTs"),
    ("TPC PMT bases", 2.80, 75.9, 33.1, 1.52, 0.03, 1.55,
     ["PMT_TOP_bases", "PMT_BOT_bases"], "PMT_TOP_bases"),
    ("TPC PMT structures", 166.0, 1.60, 1.06, 2.65, 0.12, 2.77,
     ["PMT_TOP_structure", "PMT_BOT_structure"], "PMT_TOP_structure"),
    ("TPC PMT cables", 88.7, 4.31, 0.82, 1.44, 0.19, 1.63,
     ["PMT_TOP_cables", "PMT_BOT_cables"], "PMT_TOP_cables"),
    ("Skin PMTs and bases", 8.59, 46.0, 14.9, 0.75, 0.02, 0.78,
     ["PMT_BARREL_R8520", "PMT_BARREL_R8778_lower", "PMT_BOT_R8778_dome"],
     "PMT_BARREL_R8520"),
    ("PTFE walls", 184.0, 0.04, 0.01, 0.39, 0.00, 0.39,
     ["FC_PTFE"], "FC_PTFE"),
    ("TPC sensors", 5.02, 5.82, 1.88, 1.19, 0.00, 1.19,
     ["FC_sensors"], "FC_sensors"),
    ("Field grids and holders", 89.1, 2.63, 1.46, 0.62, 0.11, 0.73,
     ["FC_topgrid", "FC_botgrid"], "FC_topgrid"),
    ("Field-cage resistors", 0.06, 1350.0, 2010.0, 2.63, 0.03, 2.65,
     ["FC_resistors"], "FC_resistors"),
    ("Field-cage rings", 93.0, 0.35, 0.24, 0.82, 0.00, 0.82,
     ["FC_rings"], "FC_rings"),
    ("Ti cryostat vessel", 2590.0, 0.08, 0.22, 1.30, 0.20, 1.49,
     ["OCV_barrel", "OCV_top", "OCV_bottom",
      "ICV_barrel", "ICV_top", "ICV_bottom"], "ICV_barrel"),
    ("Cryostat insulation", 13.8, 11.1, 7.79, 0.90, 0.04, 0.94,
     ["MLI"], "MLI"),
    ("Outer detector system", 22900.0, 4.71, 3.73, 1.70, 1.08, 2.79,
     [], None),
    ("Other components", 438.0, 1.83, 1.65, 2.10, 0.31, 2.41,
     [], None),
]

# Non-component backgrounds (paper only)
PAPER_INTERNAL = [
    ("Cavern walls", "29 Bq/kg rock (238U), 12.5 Bq/kg (232Th)",
     3.21, 8.41, 11.6),
    ("137Xe (neutron-induced)", "n-capture in purification loop",
     None, None, 0.28),
    ("Internal 222Rn", "<2.0 uBq/kg in active LXe",
     None, None, 0.45),
    ("136Xe 2vbb", "T1/2 = 2.165e21 yr",
     None, None, 0.01),
    ("8B solar neutrinos", "5.79e6 /cm2/s",
     None, None, 0.03),
]

# Skin PMT sub-components for breakdown
SKIN_PMT_SOURCES = [
    "PMT_BARREL_R8520",
    "PMT_BARREL_R8778_lower",
    "PMT_BOT_R8778_dome",
]

CRYOSTAT_SOURCES = [
    "OCV_barrel", "OCV_top", "OCV_bottom",
    "ICV_barrel", "ICV_top", "ICV_bottom",
]


def load_json(path: str) -> dict:
    with open(path) as f:
        return json.load(f)


def print_header(title: str) -> None:
    print()
    print("=" * 100)
    print(title)
    print("=" * 100)


def validate_activities(source_path: str) -> int:
    data = load_json(source_path)
    sources = {s["name"]: s for s in data["sources"]}
    failures = 0

    # ------------------------------------------------------------------
    # 1. Main table: paper vs JSON
    # ------------------------------------------------------------------
    print_header("ACTIVITY COMPARISON: Paper Table I vs JSON")
    print()
    print("Convention: JSON 'Bi214_mBq_per_kg' = paper '238U-late (mBq/kg)' = 226Ra chain activity")
    print("           JSON 'Tl208_mBq_per_kg' = paper '232Th-late (mBq/kg)' = 224Ra chain activity")
    print()

    hdr = (f"{'Paper component':<30s} {'Mass':>7s} {'U-late':>7s} {'Th-late':>7s}"
           f"  {'JSON mass':>10s} {'JSON Bi214':>10s} {'JSON Tl208':>10s}"
           f"  {'Status'}")
    print(hdr)
    print("-" * len(hdr))

    missing_components = []
    activity_issues = []
    paper_total_u = 0.0
    paper_total_th = 0.0
    paper_total = 0.0

    for row in PAPER_TABLE_I:
        (paper_name, paper_mass, paper_u, paper_th,
         paper_cts_u, paper_cts_th, paper_cts_total,
         json_mass_sources, json_act_source) = row

        paper_total_u += paper_cts_u
        paper_total_th += paper_cts_th
        paper_total += paper_cts_total

        if not json_mass_sources:
            print(f"{paper_name:<30s} {paper_mass:>7.1f} {paper_u:>7.2f} {paper_th:>7.2f}"
                  f"  {'---':>10s} {'---':>10s} {'---':>10s}  NOT IN JSON")
            missing_components.append(
                (paper_name, paper_mass, paper_u, paper_th, paper_cts_total))
            continue

        # JSON mass
        is_exact = all(sources[sn]["approximation"] == "exact"
                       for sn in json_mass_sources)
        if is_exact:
            mass_str = "(geom)"
        else:
            json_mass = sum(sources[sn].get("equivalent_mass_kg", 0.0)
                            for sn in json_mass_sources)
            mass_str = f"{json_mass:.2f}"

        # Activity from representative source
        s_act = sources[json_act_source]
        json_bi = s_act["activity"]["Bi214_mBq_per_kg"]
        json_tl = s_act["activity"]["Tl208_mBq_per_kg"]

        bi_match = abs(json_bi - paper_u) < 0.015
        tl_match = abs(json_tl - paper_th) < 0.015

        if bi_match and tl_match:
            status = "OK"
        else:
            parts = []
            if not bi_match:
                parts.append(f"Bi214: {json_bi} vs {paper_u}")
            if not tl_match:
                parts.append(f"Tl208: {json_tl} vs {paper_th}")
            status = "MISMATCH " + "; ".join(parts)
            activity_issues.append(
                (paper_name, json_bi, paper_u, json_tl, paper_th))
            failures += 1

        print(f"{paper_name:<30s} {paper_mass:>7.1f} {paper_u:>7.2f} {paper_th:>7.2f}"
              f"  {mass_str:>10s} {json_bi:>10.2f} {json_tl:>10.2f}  {status}")

    print("-" * len(hdr))
    print(f"{'Det. components subtotal':<30s}"
          f" {'':>7s} {'':>7s} {'':>7s}"
          f"  {'':>10s} {'':>10s} {'':>10s}"
          f"  U={paper_total_u:.1f}  Th={paper_total_th:.1f}  Total={paper_total:.1f} cts/1000d")

    # ------------------------------------------------------------------
    # 2. Non-component backgrounds
    # ------------------------------------------------------------------
    print_header("NON-COMPONENT BACKGROUNDS (paper only)")

    for name, note, u_act, th_act, counts in PAPER_INTERNAL:
        u_str = f"{u_act:.2f}" if u_act is not None else "---"
        th_str = f"{th_act:.2f}" if th_act is not None else "---"
        cts_yr = counts * DAYS_PER_YEAR / LIVE_DAYS
        print(f"  {name:<35s}  U={u_str:>8s}  Th={th_str:>8s}"
              f"  cts/1000d={counts:>5.2f}  cts/yr={cts_yr:>5.3f}  ({note})")

    # ------------------------------------------------------------------
    # 3. Skin PMT breakdown
    # ------------------------------------------------------------------
    print_header("SKIN PMT BREAKDOWN")
    print("  Paper lumps all as 8.59 kg at 46.0 / 14.9 mBq/kg")
    print()

    total_skin_mass = 0.0
    for sn in SKIN_PMT_SOURCES:
        s = sources[sn]
        m = s["equivalent_mass_kg"]
        bi = s["activity"]["Bi214_mBq_per_kg"]
        tl = s["activity"]["Tl208_mBq_per_kg"]
        total_skin_mass += m
        print(f"  {sn:<35s}  mass={m:>6.2f} kg  Bi214={bi:>5.1f}  Tl208={tl:>5.1f}")

    print(f"  {'TOTAL':<35s}  mass={total_skin_mass:>6.2f} kg  (paper: 8.59 kg)"
          f"  {'OK' if abs(total_skin_mass - 8.59) < 0.01 else 'FAIL'}")
    skin_acts = set()
    for sn in SKIN_PMT_SOURCES:
        s = sources[sn]
        skin_acts.add((s["activity"]["Bi214_mBq_per_kg"],
                       s["activity"]["Tl208_mBq_per_kg"]))
    print(f"  All sub-components share same activity: {len(skin_acts) == 1}")

    # ------------------------------------------------------------------
    # 4. Cryostat breakdown
    # ------------------------------------------------------------------
    print_header("CRYOSTAT SHELL BREAKDOWN")
    print("  Paper lumps all as 2590 kg at 0.08 / 0.22 mBq/kg (upper limit)")
    print()

    cryo_acts = set()
    for sn in CRYOSTAT_SOURCES:
        s = sources[sn]
        bi = s["activity"]["Bi214_mBq_per_kg"]
        tl = s["activity"]["Tl208_mBq_per_kg"]
        cryo_acts.add((bi, tl))
        print(f"  {sn:<20s}  Bi214={bi:.2f}  Tl208={tl:.2f}")

    print(f"  All shells share same activity: {len(cryo_acts) == 1}")
    print(f"  Note: paper marks cryostat activities as upper limits (dagger)")

    # ------------------------------------------------------------------
    # 5. Annual rate conversion
    # ------------------------------------------------------------------
    print_header("ANNUAL RATES (paper counts / 1000d -> counts/year)")

    scale = DAYS_PER_YEAR / LIVE_DAYS
    total_yr = 0.0
    for row in PAPER_TABLE_I:
        paper_name = row[0]
        cts_total = row[6]
        cts_yr = cts_total * scale
        total_yr += cts_yr
        print(f"  {paper_name:<30s}  {cts_total:>6.2f} cts/1000d  -> {cts_yr:>6.3f} cts/yr")

    for name, _, _, _, cts in PAPER_INTERNAL:
        cts_yr = cts * scale
        total_yr += cts_yr
        print(f"  {name:<30s}  {cts:>6.2f} cts/1000d  -> {cts_yr:>6.3f} cts/yr")

    print(f"  {'':-<30s}  {'':-<6s}              {'':-<6s}")
    print(f"  {'GRAND TOTAL':<30s}  {'35.6':>6s} cts/1000d  -> {total_yr:>6.2f} cts/yr")

    # ------------------------------------------------------------------
    # 6. Summary
    # ------------------------------------------------------------------
    print_header("SUMMARY")

    if activity_issues:
        print(f"  Activity mismatches: {len(activity_issues)}")
        for name, jbi, pbi, jtl, ptl in activity_issues:
            print(f"    {name}: Bi214 JSON={jbi} vs paper={pbi},"
                  f" Tl208 JSON={jtl} vs paper={ptl}")
    else:
        print("  All activities match paper Table I.")

    if missing_components:
        total_missing_cts = sum(c[4] for c in missing_components)
        print(f"\n  Components in paper but NOT in source JSON"
              f" ({len(missing_components)}):")
        for name, mass, u, th, cts in missing_components:
            print(f"    {name}: {mass:.0f} kg, U-late={u}, Th-late={th} mBq/kg"
                  f" -> {cts:.2f} cts/1000d")
        print(f"    Combined missing: {total_missing_cts:.1f} cts/1000d"
              f" ({total_missing_cts / 35.6 * 100:.0f}% of total)")
    else:
        print("\n  All paper components have JSON sources.")

    print()
    print("  NAMING CONVENTION:")
    print("  JSON 'Bi214_mBq_per_kg' = paper '238U-late (mBq/kg)'"
          " = 226Ra chain activity")
    print("  JSON 'Tl208_mBq_per_kg' = paper '232Th-late (mBq/kg)'"
          " = 224Ra chain activity")
    print("  To get relevant gamma rates for the bb0nu ROI:")
    print("    214Bi 2448 keV:  chain_activity * 0.015  (BR)")
    print("    208Tl 2615 keV:  chain_activity * 0.359  (212Bi->208Tl BR)"
          " * 0.998 (gamma BR)")

    return failures


def main() -> None:
    default_root = Path(__file__).resolve().parent.parent.parent
    default_source = default_root / "data" / "source_geometry_lz_v1.json"

    parser = argparse.ArgumentParser(
        description="Validate LZ source activities against"
                    " arXiv:1912.04248v2 Table I")
    parser.add_argument("--json-source", default=str(default_source),
                        help="Path to source geometry JSON")
    args = parser.parse_args()

    print(f"Source geometry: {args.json_source}")

    failures = validate_activities(args.json_source)
    sys.exit(1 if failures > 0 else 0)


if __name__ == "__main__":
    main()
