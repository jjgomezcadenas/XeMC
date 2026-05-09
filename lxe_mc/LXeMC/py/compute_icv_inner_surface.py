#!/usr/bin/env python3
"""Compute ICV inner surface dimensions from source geometry.

Compares with current tracking detector geometry to identify mismatches.
Used to guide geometry fixes for exact AirDome/LXe_passive placement.
"""

import json
import sys
from pathlib import Path


def main() -> None:
    data_dir = Path(__file__).resolve().parent.parent.parent / "data"

    src_path = data_dir / "source_geometry_lz_v1.json"
    det_path = data_dir / "detector_lz_v3.json"

    if not src_path.exists():
        print(f"Not found: {src_path}", file=sys.stderr)
        sys.exit(1)

    src = json.loads(src_path.read_text())
    det = json.loads(det_path.read_text())

    sources = {s["name"]: s for s in src["sources"]}
    volumes = {v["name"]: v for v in det["volumes"]}

    print("=" * 70)
    print("ICV source geometry (outer and inner surfaces)")
    print("=" * 70)

    # --- Barrel ---
    icv_b = sources["ICV_barrel"]
    Ri = icv_b["R_inner_cm"]
    t = icv_b["wall_thickness_cm"]
    hh = icv_b["half_height_cm"]
    z_b = icv_b["position_cm"][2]
    print(f"\nICV_barrel:")
    print(f"  Outer: R = {Ri + t:.1f}, z = [{z_b - hh:.3f}, {z_b + hh:.3f}]")
    print(f"  Inner: R = {Ri:.1f}, z = [{z_b - hh:.3f}, {z_b + hh:.3f}]")
    print(f"  Wall thickness: {t} cm")

    # --- Top head ---
    icv_t = sources["ICV_top"]
    R_t = icv_t["radius_cm"]
    t_t = icv_t["wall_thickness_cm"]
    ar_t = icv_t["aspect_ratio"]
    z_t = icv_t["position_cm"][2]
    depth_outer_t = R_t / ar_t
    R_inner_t = R_t - t_t
    # Inner dome: same aspect ratio, smaller radius
    depth_inner_t = R_inner_t / ar_t
    print(f"\nICV_top (orientation: {icv_t['orientation']}):")
    print(f"  Outer: R = {R_t:.1f}, z_equator = {z_t:.3f}, ar = {ar_t}")
    print(f"         dome depth = {depth_outer_t:.3f}, dome top = {z_t + depth_outer_t:.3f}")
    print(f"  Inner: R = {R_inner_t:.1f}, ar = {ar_t}")
    print(f"         dome depth = {depth_inner_t:.3f}, dome top = {z_t + depth_inner_t:.3f}")
    print(f"  Wall thickness: {t_t} cm")

    # --- Bottom head ---
    icv_bot = sources["ICV_bottom"]
    R_bot = icv_bot["radius_cm"]
    t_bot = icv_bot["wall_thickness_cm"]
    ar_bot = icv_bot["aspect_ratio"]
    z_bot = icv_bot["position_cm"][2]
    depth_outer_bot = R_bot / ar_bot
    R_inner_bot = R_bot - t_bot
    depth_inner_bot = R_inner_bot / ar_bot
    print(f"\nICV_bottom (orientation: {icv_bot['orientation']}):")
    print(f"  Outer: R = {R_bot:.1f}, z_equator = {z_bot:.3f}, ar = {ar_bot}")
    print(f"         dome depth = {depth_outer_bot:.3f}, dome bottom = {z_bot - depth_outer_bot:.3f}")
    print(f"  Inner: R = {R_inner_bot:.1f}, ar = {ar_bot}")
    print(f"         dome depth = {depth_inner_bot:.3f}, dome bottom = {z_bot - depth_inner_bot:.3f}")
    print(f"  Wall thickness: {t_bot} cm")

    print("\n" + "=" * 70)
    print("Current tracking detector geometry")
    print("=" * 70)

    lz = volumes["LZ_detector"]
    print(f"\nLZ_detector:")
    print(f"  R = {lz['radius_cm']}, barrel_hh = {lz['barrel_half_height_cm']}")
    print(f"  z_center = {lz['position_cm'][2]}")
    print(f"  top_ar = {lz['top_aspect_ratio']}, bot_ar = {lz['bottom_aspect_ratio']}")

    air = volumes["AirDome"]
    print(f"\nAirDome:")
    print(f"  R = {air['radius_cm']}, z_equator = {air['position_cm'][2]}, ar = {air['aspect_ratio']}")
    air_depth = air["radius_cm"] / air["aspect_ratio"]
    print(f"  dome top = {air['position_cm'][2] + air_depth:.3f}")

    lxe = volumes["LXe_passive"]
    print(f"\nLXe_passive:")
    print(f"  R = {lxe['radius_cm']}, barrel_hh = {lxe['barrel_half_height_cm']}")
    print(f"  z_center = {lxe['position_cm'][2]}")
    print(f"  bot_ar = {lxe['bottom_aspect_ratio']}")
    lxe_z = lxe["position_cm"][2]
    lxe_hh = lxe["barrel_half_height_cm"]
    lxe_bot_depth = lxe["radius_cm"] / lxe["bottom_aspect_ratio"]
    print(f"  barrel z = [{lxe_z - lxe_hh:.3f}, {lxe_z + lxe_hh:.3f}]")
    print(f"  bottom dome depth = {lxe_bot_depth:.3f}")
    print(f"  bottom dome bottom = {lxe_z - lxe_hh - lxe_bot_depth:.3f}")

    print("\n" + "=" * 70)
    print("Mismatches")
    print("=" * 70)

    # Barrel R
    if Ri != lz["radius_cm"]:
        print(f"\n  BARREL R: ICV inner = {Ri}, LZ_detector = {lz['radius_cm']}")
    else:
        print(f"\n  Barrel R: OK ({Ri})")

    # Top
    print(f"\n  TOP:")
    print(f"    AirDome R = {air['radius_cm']}, ICV_top inner R = {R_inner_t:.1f}")
    print(f"    AirDome z_eq = {air['position_cm'][2]}, ICV_top z_eq = {z_t}")
    print(f"    AirDome dome top = {air['position_cm'][2] + air_depth:.3f}, "
          f"ICV_top inner dome top = {z_t + depth_inner_t:.3f}")

    # Bottom
    print(f"\n  BOTTOM:")
    print(f"    LXe_passive R = {lxe['radius_cm']}, ICV_bottom inner R = {R_inner_bot:.1f}")
    print(f"    LXe_passive bot_ar = {lxe['bottom_aspect_ratio']}, ICV_bottom ar = {ar_bot}")
    lxe_bot_z = lxe_z - lxe_hh
    print(f"    LXe_passive barrel bottom = {lxe_bot_z:.3f}, ICV_bottom z_eq = {z_bot}")
    print(f"    LXe_passive dome bottom = {lxe_bot_z - lxe_bot_depth:.3f}, "
          f"ICV_bottom inner dome bottom = {z_bot - depth_inner_bot:.3f}")


if __name__ == "__main__":
    main()
