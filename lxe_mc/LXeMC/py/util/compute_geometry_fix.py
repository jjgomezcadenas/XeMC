#!/usr/bin/env python3
"""Compute detector geometry with uniform 9 mm ICV wall thickness.

Derives the exact inner surface dimensions for barrel, top cap, and
bottom cap, and the dimensions of the AirCyl/LXeCyl gap fillers.
"""


def main() -> None:
    R_outer = 83.0  # ICV outer radius (all surfaces)
    t = 0.9         # uniform wall thickness [cm]
    R_inner = R_outer - t

    ar_top = 2.0
    ar_bot = 3.0
    z_eq_top = 148.5    # ICV top equator
    z_eq_bot = -41.33   # ICV bottom equator
    z_gate = 145.6      # liquid-gas boundary

    depth_inner_top = R_inner / ar_top
    depth_inner_bot = R_inner / ar_bot

    barrel_hh = (z_eq_top - z_eq_bot) / 2.0
    barrel_z_center = (z_eq_top + z_eq_bot) / 2.0

    print("=" * 60)
    print(f"Uniform ICV wall thickness: {t} cm")
    print(f"R_outer = {R_outer}, R_inner = {R_inner}")
    print("=" * 60)

    print(f"\n--- Barrel ---")
    print(f"  R_inner = {R_inner}")
    print(f"  z = [{z_eq_bot:.3f}, {z_eq_top:.3f}]")
    print(f"  half_height = {barrel_hh:.3f}")
    print(f"  z_center = {barrel_z_center:.3f}")

    print(f"\n--- Top cap (AirDome) ---")
    print(f"  shape: cap, orientation: up")
    print(f"  R = {R_inner}, ar = {ar_top}")
    print(f"  z_equator = {z_eq_top}")
    print(f"  dome depth = {depth_inner_top:.3f}")
    print(f"  dome top = {z_eq_top + depth_inner_top:.3f}")

    print(f"\n--- AirCyl (gas gap: gate to top equator) ---")
    print(f"  shape: cylinder")
    print(f"  R = {R_inner}")
    air_hh = (z_eq_top - z_gate) / 2.0
    air_z = (z_eq_top + z_gate) / 2.0
    print(f"  z = [{z_gate}, {z_eq_top}]")
    print(f"  half_height = {air_hh:.3f}")
    print(f"  z_center = {air_z:.3f}")

    print(f"\n--- Bottom cap (LXeDome) ---")
    print(f"  shape: cap, orientation: down")
    print(f"  R = {R_inner}, ar = {ar_bot}")
    print(f"  z_equator = {z_eq_bot}")
    print(f"  dome depth = {depth_inner_bot:.3f}")
    print(f"  dome bottom = {z_eq_bot - depth_inner_bot:.3f}")

    # LXeCyl: from ICV bottom equator up to ... where?
    # The cathode is at z=0. Below the cathode down to ICV_bottom equator
    # is passive LXe (no drift field). This is the LXeCyl.
    # But actually the barrel extends down to z_eq_bot = -41.33, and the
    # FC/active regions start at z=-13.75 (FC bottom). Below FC is passive LXe.
    # The LXeCyl fills the space between the bottom dome equator and the
    # bottom of the barrel region that is NOT covered by the dome.
    # Since the dome equator IS the barrel bottom, there's no gap — the dome
    # starts right where the barrel ends.
    #
    # Wait: for the top, AirCyl fills [gate, ICV_top_equator] because the
    # barrel LXe goes up to the gate, then there's a gas gap up to the dome.
    # For the bottom, the barrel LXe goes all the way down to z_eq_bot,
    # and the dome starts there. There's no gap equivalent to AirCyl.
    # The LXe just transitions from cylindrical barrel to domed bottom.

    print(f"\n--- LXeCyl (if needed) ---")
    print(f"  The bottom dome equator ({z_eq_bot}) coincides with the")
    print(f"  barrel bottom — no gap to fill. LXeCyl is NOT needed.")
    print(f"  The LXe_passive capped_cylinder handles this directly.")

    print(f"\n--- Summary: new LZ_detector ---")
    print(f"  shape: domed_container")
    print(f"  R = {R_inner}")
    print(f"  barrel_half_height = {barrel_hh:.3f}")
    print(f"  z_center = {barrel_z_center:.3f}")
    print(f"  top_ar = {ar_top}, bottom_ar = {ar_bot}")
    print(f"  z_range = [{z_eq_bot - depth_inner_bot:.3f}, {z_eq_top + depth_inner_top:.3f}]")

    print(f"\n--- Summary: new LXe_passive ---")
    lxe_hh = (z_gate - z_eq_bot) / 2.0
    lxe_z = (z_gate + z_eq_bot) / 2.0
    print(f"  shape: capped_cylinder")
    print(f"  R = {R_inner}")
    print(f"  barrel_half_height = {lxe_hh:.3f}")
    print(f"  z_center = {lxe_z:.3f}")
    print(f"  barrel z = [{z_eq_bot:.3f}, {z_gate:.3f}]")
    print(f"  bottom_ar = {ar_bot}")

    print(f"\n--- Source geometry fix ---")
    print(f"  ICV_top: change wall_thickness from 0.8 to {t}")
    print(f"  ICV_bottom: change wall_thickness from 1.2 to {t}")
    print(f"  All ICV inner R become {R_inner}")


if __name__ == "__main__":
    main()
