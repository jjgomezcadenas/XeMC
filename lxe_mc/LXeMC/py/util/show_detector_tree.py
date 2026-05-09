#!/usr/bin/env python3
"""Print the tracking detector geometry tree with positions and dimensions.

Reads the flat volumes list from detector_lz_v3.json and displays it as
an indented tree based on parent references.
"""

import argparse
import json
import sys
from pathlib import Path


def format_volume(vol: dict) -> str:
    """Format one volume entry for display."""
    name = vol.get("name", "?")
    shape = vol.get("shape", "?")
    tag = vol.get("tag", "")
    mat = vol.get("material", "")
    pos = vol.get("position_cm", [0, 0, 0])
    z = pos[2]

    base = f"{name:25s}  tag={tag:12s}  mat={mat:6s}  z={z:8.2f}"

    if shape == "cylinder":
        R = vol.get("radius_cm", 0)
        hh = vol.get("half_height_cm", 0)
        return f"{base}  R={R:6.1f}  hh={hh:7.2f}  z=[{z-hh:.2f}, {z+hh:.2f}]"

    elif shape == "cylinder_shell":
        Ri = vol.get("R_inner_cm", 0)
        t = vol.get("wall_thickness_cm", 0)
        hh = vol.get("half_height_cm", 0)
        return (
            f"{base}  Ri={Ri:6.1f}  t={t:5.3f}  hh={hh:7.2f}  "
            f"z=[{z-hh:.2f}, {z+hh:.2f}]"
        )

    elif shape == "domed_container":
        R = vol.get("radius_cm", 0)
        hh = vol.get("barrel_half_height_cm", 0)
        top_ar = vol.get("top_aspect_ratio", 0)
        bot_ar = vol.get("bottom_aspect_ratio", 0)
        return (
            f"{base}  R={R:6.1f}  hh={hh:7.2f}  "
            f"top_ar={top_ar}  bot_ar={bot_ar}  z=[{z-hh:.2f}, {z+hh:.2f}]"
        )

    elif shape == "capped_cylinder":
        R = vol.get("radius_cm", 0)
        hh = vol.get("barrel_half_height_cm", 0)
        bot_ar = vol.get("bottom_aspect_ratio", 0)
        return (
            f"{base}  R={R:6.1f}  hh={hh:7.2f}  "
            f"bot_ar={bot_ar}  z=[{z-hh:.2f}, {z+hh:.2f}]"
        )

    elif shape == "cap":
        R = vol.get("radius_cm", 0)
        ar = vol.get("aspect_ratio", 0)
        ori = vol.get("orientation", "?")
        depth = R / ar if ar > 0 else 0
        return f"{base}  R={R:6.1f}  ar={ar}  ori={ori}  depth={depth:.2f}"

    else:
        return f"{base}  shape={shape}"


def build_tree(world: dict, volumes: list) -> None:
    """Print the geometry as an indented tree."""
    # Build children map
    children: dict[str, list] = {}
    by_name: dict[str, dict] = {}

    by_name[world["name"]] = world
    for vol in volumes:
        by_name[vol["name"]] = vol
        parent = vol.get("parent", "")
        children.setdefault(parent, []).append(vol)

    def print_node(name: str, indent: int = 0) -> None:
        vol = by_name[name]
        prefix = "  " * indent
        if name == world["name"]:
            # World volume
            R = vol.get("radius_cm", 0)
            hh = vol.get("half_height_cm", 0)
            print(f"{prefix}{name:25s}  (world)  R={R}  hh={hh}")
        else:
            print(f"{prefix}{format_volume(vol)}")
        for child in children.get(name, []):
            print_node(child["name"], indent + 1)

    print_node(world["name"])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--detector",
        type=Path,
        default=Path(__file__).resolve().parent.parent.parent
        / "data"
        / "detector_lz_v3.json",
        help="Path to detector JSON (default: data/detector_lz_v3.json)",
    )
    args = parser.parse_args()

    if not args.detector.exists():
        print(f"File not found: {args.detector}", file=sys.stderr)
        sys.exit(1)

    with open(args.detector) as f:
        det = json.load(f)

    build_tree(det["world"], det["volumes"])


if __name__ == "__main__":
    main()
