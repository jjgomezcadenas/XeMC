"""
Cross-validate source and tracking geometry JSON files against the
LZ bb0nu paper (arXiv:1912.04248v2, Table I).

Checks:
  - Cryostat shell masses (geometric mass from JSON dimensions vs paper)
  - Field cage masses (FC_PTFE, FC_rings)
  - Virtual source masses (PMTs, cables, MLI, resistors, sensors, grids)
  - Activities per component and isotope
  - Tracking geometry dimensions (FV, TPC active regions, skin, RFR)

Usage:
    python py/validate_geometry_vs_paper.py
    python py/validate_geometry_vs_paper.py --json-source data/source_geometry_lz_v1.json
    python py/validate_geometry_vs_paper.py --json-tracking data/detector_lz_v3.json
"""

import argparse
import json
import math
import os
import sys
from pathlib import Path


PI = math.pi

MATERIAL_DENSITY = {
    "Ti": 4.507,
    "PTFE": 2.20,
    "SS": 7.87,
    "LXe": 2.862,
    "Vacuum": 0.0,
    "HPGXe": 0.0,
}

DAYS_PER_YEAR = 365.25
LIVE_DAYS = 1000


def ellipsoid_half_area(radius_cm: float, aspect_ratio: float) -> float:
    a = radius_cm
    c = a / aspect_ratio
    if c < a:
        e = math.sqrt(1.0 - (c / a) ** 2)
        S_full = 2 * PI * a**2 + PI * (c**2 / e) * math.log((1 + e) / (1 - e))
    else:
        e = math.sqrt(1.0 - (a / c) ** 2)
        S_full = 2 * PI * a**2 + 2 * PI * a * c * math.asin(e) / e
    return S_full / 2.0


def geometric_mass_kg(source: dict, density: float) -> float:
    shape = source["shape"]
    if shape == "cylinder_shell":
        R_in = source["R_inner_cm"]
        t = source["wall_thickness_cm"]
        H = 2.0 * source["half_height_cm"]
        V = PI * ((R_in + t) ** 2 - R_in**2) * H
        return V * density / 1000.0
    elif shape == "cylinder":
        R = source["radius_cm"]
        H = 2.0 * source["half_height_cm"]
        return PI * R**2 * H * density / 1000.0
    elif shape == "disk":
        R = source["radius_cm"]
        t = source["wall_thickness_cm"]
        n = source["aspect_ratio"]
        S_half = ellipsoid_half_area(R, n)
        return S_half * t * density / 1000.0
    return float("nan")


def load_json(path: str) -> dict:
    with open(path) as f:
        return json.load(f)


def print_header(title: str) -> None:
    print()
    print("=" * 80)
    print(title)
    print("=" * 80)


def print_check(name: str, json_val: float, paper_val: float,
                unit: str = "kg", rtol: float = 0.05) -> bool:
    if paper_val == 0:
        ok = json_val == 0
    else:
        ok = abs(json_val - paper_val) / abs(paper_val) < rtol
    status = "OK" if ok else "FAIL"
    print(f"  {name:40s}  JSON={json_val:10.2f} {unit}  paper={paper_val:10.2f} {unit}  [{status}]")
    return ok


def validate_source_geometry(source_path: str) -> int:
    data = load_json(source_path)
    sources = {s["name"]: s for s in data["sources"]}
    failures = 0

    # ------------------------------------------------------------------
    # 1. Cryostat masses (exact geometry -> geometric mass)
    # ------------------------------------------------------------------
    print_header("CRYOSTAT MASS CHECKS (geometric mass vs paper)")

    paper_masses = {
        "OCV_barrel": 385.7,
        "OCV_top": 147.4,
        "OCV_bottom": 245.6,
        "ICV_barrel": 401.9,
        "ICV_top": 107.8,
        "ICV_bottom": 141.5,
    }

    cryo_total_json = 0.0
    for name, paper_m in paper_masses.items():
        s = sources[name]
        rho = MATERIAL_DENSITY[s["material"]]
        m = geometric_mass_kg(s, rho)
        cryo_total_json += m
        if not print_check(name, m, paper_m):
            failures += 1

    paper_cryo_total = 2590.0
    # Paper total includes more than just these 6 pieces (flanges, ports, etc.)
    print(f"\n  Cryostat sum (6 shells) = {cryo_total_json:.1f} kg")
    print(f"  Paper total (full vessel) = {paper_cryo_total:.0f} kg")
    print(f"  Ratio = {cryo_total_json / paper_cryo_total:.2%}")
    print(f"  (Difference is flanges, ports, and other structural elements)")

    # ------------------------------------------------------------------
    # 2. Field cage masses
    # ------------------------------------------------------------------
    print_header("FIELD CAGE MASS CHECKS")

    fc_checks = {"FC_PTFE": 184.0, "FC_rings": 93.0}
    for name, paper_m in fc_checks.items():
        s = sources[name]
        rho = MATERIAL_DENSITY[s["material"]]
        m = geometric_mass_kg(s, rho)
        if not print_check(name, m, paper_m):
            failures += 1

    # ------------------------------------------------------------------
    # 3. Virtual source masses (equivalent_mass_kg)
    # ------------------------------------------------------------------
    print_header("VIRTUAL / LUMPED SOURCE MASS CHECKS")

    paper_virtual_masses = {
        "MLI": ("Cryo insulation", 13.8),
        "FC_resistors": ("FC resistors", 0.06),
        "FC_sensors": ("TPC sensors", 5.02),
    }
    for json_name, (label, paper_m) in paper_virtual_masses.items():
        s = sources[json_name]
        json_m = s["equivalent_mass_kg"]
        if not print_check(label, json_m, paper_m):
            failures += 1

    # Aggregated PMT checks
    aggregates = [
        ("TPC PMTs (top+bot)",
         ["PMT_TOP_PMTs", "PMT_BOT_PMTs"], 91.9, "equivalent_mass_kg"),
        ("TPC PMT bases (top+bot)",
         ["PMT_TOP_bases", "PMT_BOT_bases"], 2.80, "equivalent_mass_kg"),
        ("TPC PMT structures (top+bot)",
         ["PMT_TOP_structure", "PMT_BOT_structure"], 166.0, "equivalent_mass_kg"),
        ("TPC PMT cables (top+bot)",
         ["PMT_TOP_cables", "PMT_BOT_cables"], 88.7, "equivalent_mass_kg"),
        ("Skin PMTs+bases (all)",
         ["PMT_BARREL_R8520", "PMT_BARREL_R8778_lower", "PMT_BOT_R8778_dome"], 8.59,
         "equivalent_mass_kg"),
        ("Field grids+holders (top+bot)",
         ["FC_topgrid", "FC_botgrid"], 89.1, "equivalent_mass_kg"),
    ]
    for label, names, paper_m, field in aggregates:
        total = sum(sources[n][field] for n in names)
        if not print_check(label, total, paper_m):
            failures += 1

    # ------------------------------------------------------------------
    # 4. Activities
    # ------------------------------------------------------------------
    print_header("ACTIVITY CHECKS (JSON vs paper Table I)")

    activity_checks = [
        ("Ti cryostat (ICV/OCV)", "ICV_barrel", 0.08, 0.22),
        ("PTFE walls", "FC_PTFE", 0.04, 0.01),
        ("FC rings", "FC_rings", 0.35, 0.24),
        ("FC resistors", "FC_resistors", 1350.0, 2010.0),
        ("TPC sensors", "FC_sensors", 5.82, 1.88),
        ("TPC PMTs", "PMT_TOP_PMTs", 3.22, 1.61),
        ("TPC PMT bases", "PMT_TOP_bases", 75.9, 33.1),
        ("TPC PMT structures", "PMT_TOP_structure", 1.60, 1.06),
        ("TPC PMT cables", "PMT_TOP_cables", 4.31, 0.82),
        ("Skin PMTs+bases", "PMT_BARREL_R8520", 46.0, 14.9),
        ("Cryo insulation (MLI)", "MLI", 11.1, 7.79),
        ("Field grids+holders", "FC_topgrid", 2.63, 1.46),
    ]

    for label, json_name, paper_bi214, paper_tl208 in activity_checks:
        s = sources[json_name]
        json_bi = s["activity"]["Bi214_mBq_per_kg"]
        json_tl = s["activity"]["Tl208_mBq_per_kg"]
        bi_ok = abs(json_bi - paper_bi214) < 0.01
        tl_ok = abs(json_tl - paper_tl208) < 0.01
        bi_status = "OK" if bi_ok else f"FAIL({json_bi})"
        tl_status = "OK" if tl_ok else f"FAIL({json_tl})"
        print(f"  {label:35s}  Bi214={bi_status:12s}  Tl208={tl_status:12s}")
        if not bi_ok:
            failures += 1
        if not tl_ok:
            failures += 1

    return failures


def validate_tracking_geometry(tracking_path: str) -> int:
    data = load_json(tracking_path)
    volumes = {v["name"]: v for v in data["volumes"]}
    world = data["world"]
    failures = 0

    print_header("TRACKING GEOMETRY CHECKS (JSON vs paper)")

    # Paper: TPC drift = 145.6 cm (cathode at z=0, gate at z=145.6)
    # Paper: TPC inner diameter = 145.6 cm -> radius = 72.8 cm
    # Paper: FV = 26 < z < 96, r < 39 cm, ~967 kg
    # Paper: RFR = 13.75 cm below cathode
    # Paper: skin = 4 cm at top, 8 cm at cathode level

    print("\n  --- FV ---")
    fv = volumes["FV"]
    fv_r = fv["radius_cm"]
    fv_hh = fv["half_height_cm"]
    fv_z0 = fv["position_cm"][2]
    fv_zmin = fv_z0 - fv_hh
    fv_zmax = fv_z0 + fv_hh
    print(f"  FV radius:  {fv_r} cm  (paper: 39.0 cm)  {'OK' if fv_r == 39.0 else 'FAIL'}")
    print(f"  FV z range: [{fv_zmin}, {fv_zmax}] cm  (paper: [26, 96] cm)  {'OK' if fv_zmin == 26.0 and fv_zmax == 96.0 else 'FAIL'}")
    if fv_r != 39.0 or fv_zmin != 26.0 or fv_zmax != 96.0:
        failures += 1

    print("\n  --- TPC active regions ---")
    for name in ["TopActive", "BarrelActive", "BottomActive"]:
        v = volumes[name]
        z0 = v["position_cm"][2]
        if v["shape"] == "cylinder":
            r = v["radius_cm"]
            hh = v["half_height_cm"]
            print(f"  {name:15s}  r={r} cm, z=[{z0-hh}, {z0+hh}] cm")
            if r != 72.8:
                print(f"    WARNING: radius {r} != 72.8 (TPC inner)")
                failures += 1
        elif v["shape"] == "cylinder_shell":
            ri = v["R_inner_cm"]
            t = v["wall_thickness_cm"]
            hh = v["half_height_cm"]
            print(f"  {name:15s}  R_in={ri}, R_out={ri+t} cm, z=[{z0-hh}, {z0+hh}] cm")

    # Check full drift coverage
    top = volumes["TopActive"]
    bot = volumes["BottomActive"]
    top_zmax = top["position_cm"][2] + top["half_height_cm"]
    bot_zmin = bot["position_cm"][2] - bot["half_height_cm"]
    print(f"\n  Drift coverage: z=[{bot_zmin}, {top_zmax}] cm  (paper: [0, 145.6] cm)")
    if abs(bot_zmin) > 0.1:
        print(f"    WARNING: bottom does not start at cathode (z=0)")
        failures += 1
    if abs(top_zmax - 145.6) > 0.1:
        print(f"    WARNING: top does not reach gate (z=145.6)")
        failures += 1

    print("\n  --- Field cage ---")
    ptfe = volumes["FC_PTFE"]
    rings = volumes["FC_rings"]
    ptfe_z0 = ptfe["position_cm"][2]
    ptfe_hh = ptfe["half_height_cm"]
    print(f"  FC_PTFE:  R=[{ptfe['R_inner_cm']}, {ptfe['R_inner_cm']+ptfe['wall_thickness_cm']:.3f}] cm, "
          f"z=[{ptfe_z0-ptfe_hh}, {ptfe_z0+ptfe_hh}] cm")
    print(f"  FC_rings: R=[{rings['R_inner_cm']}, {rings['R_inner_cm']+rings['wall_thickness_cm']:.3f}] cm, "
          f"z=[{rings['position_cm'][2]-rings['half_height_cm']}, {rings['position_cm'][2]+rings['half_height_cm']}] cm")

    rfr_bottom = ptfe_z0 - ptfe_hh
    print(f"\n  RFR bottom (FC z_min): {rfr_bottom} cm  (paper: -13.75 cm)  "
          f"{'OK' if abs(rfr_bottom - (-13.75)) < 0.01 else 'FAIL'}")
    if abs(rfr_bottom - (-13.75)) > 0.01:
        failures += 1

    print("\n  --- Skin ---")
    skin = volumes["Skin"]
    skin_rin = skin["R_inner_cm"]
    skin_rout = skin_rin + skin["wall_thickness_cm"]
    print(f"  Skin: R=[{skin_rin}, {skin_rout}] cm, thickness={skin['wall_thickness_cm']:.3f} cm")
    print(f"  Paper: 'skin 4 cm at top, 8 cm at cathode level'")
    print(f"  Note: JSON models Skin as uniform {skin['wall_thickness_cm']:.1f} cm shell")
    print(f"        (simplified vs paper's position-dependent thickness)")

    print("\n  --- Envelope ---")
    lz = volumes["LZ_detector"]
    print(f"  LZ_detector radius: {lz['radius_cm']} cm  (ICV inner radius)")
    print(f"  ICV_barrel R_inner (source JSON): 82.1 cm  -> {'OK' if lz['radius_cm'] == 82.1 else 'FAIL'}")

    return failures


def main() -> None:
    default_root = Path(__file__).resolve().parent.parent.parent
    default_source = default_root / "data" / "source_geometry_lz_v1.json"
    default_tracking = default_root / "data" / "detector_lz_v3.json"

    parser = argparse.ArgumentParser(
        description="Validate LZ geometry JSON files against arXiv:1912.04248v2")
    parser.add_argument("--json-source", default=str(default_source),
                        help="Path to source geometry JSON")
    parser.add_argument("--json-tracking", default=str(default_tracking),
                        help="Path to tracking geometry JSON")
    args = parser.parse_args()

    print(f"Source geometry:   {args.json_source}")
    print(f"Tracking geometry: {args.json_tracking}")

    total_failures = 0
    total_failures += validate_source_geometry(args.json_source)
    total_failures += validate_tracking_geometry(args.json_tracking)

    print_header("SUMMARY")
    if total_failures == 0:
        print("  All checks passed.")
    else:
        print(f"  {total_failures} check(s) failed.")
    sys.exit(1 if total_failures > 0 else 0)


if __name__ == "__main__":
    main()
