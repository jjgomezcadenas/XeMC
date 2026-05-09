#!/usr/bin/env python3
"""Compute ICV head dimensions that preserve mass while fixing R_inner = 82.1.

Strategy: keep the original mass of each ICV head, but set R_inner = 82.1
(matching the barrel). Solve for the wall thickness t such that
mass(R_outer = 82.1 + t, ar, t) = original mass.

This pushes R_outer slightly outward for each head.
"""

import math
from scipy.optimize import brentq


def ellipsoid_half_surface(R: float, aspect_ratio: float) -> float:
    """Half-surface area of an oblate ellipsoid (dome)."""
    a = R
    c = R / aspect_ratio
    if c < a:
        e = math.sqrt(1 - (c / a) ** 2)
        S_full = 2 * math.pi * a**2 + math.pi * (c**2 / e) * math.log((1 + e) / (1 - e))
    else:
        e = math.sqrt(1 - (a / c) ** 2)
        S_full = 2 * math.pi * a**2 + 2 * math.pi * a * c * math.asin(e) / e
    return S_full / 2.0


def dome_mass_kg(R_outer: float, aspect_ratio: float, t: float, rho: float) -> float:
    """Mass of a dome shell [kg] from half-ellipsoid surface x thickness x density."""
    return ellipsoid_half_surface(R_outer, aspect_ratio) * t * rho / 1000.0


def main() -> None:
    rho_Ti = 4.507  # g/cm^3
    R_inner = 82.1  # target inner radius (same as barrel)

    print("=" * 65)
    print("ICV head dimensions: preserve mass, fix R_inner = 82.1")
    print("=" * 65)

    heads = [
        ("ICV_top",    83.0, 0.8, 2.0),
        ("ICV_bottom", 83.0, 1.2, 3.0),
    ]

    for name, R_old, t_old, ar in heads:
        m_target = dome_mass_kg(R_old, ar, t_old, rho_Ti)

        print(f"\n--- {name} ---")
        print(f"  Original: R_outer = {R_old}, t = {t_old}, ar = {ar}")
        print(f"  Original mass: {m_target:.1f} kg")

        # Solve for t: dome_mass(82.1 + t, ar, t) = m_target
        def residual(t: float) -> float:
            R_out = R_inner + t
            return dome_mass_kg(R_out, ar, t, rho_Ti) - m_target

        t_new = brentq(residual, 0.01, 5.0)
        R_outer_new = R_inner + t_new
        m_check = dome_mass_kg(R_outer_new, ar, t_new, rho_Ti)

        print(f"  New:      R_outer = {R_outer_new:.3f}, t = {t_new:.3f}, ar = {ar}")
        print(f"  New R_inner = {R_inner}")
        print(f"  New mass: {m_check:.1f} kg (target: {m_target:.1f} kg)")
        print(f"  R_outer shift: {R_outer_new - R_old:.3f} cm")

    # Also show barrel for reference
    print(f"\n--- ICV_barrel (reference, unchanged) ---")
    print(f"  R_inner = {R_inner}, t = 0.9, R_outer = {R_inner + 0.9}")


if __name__ == "__main__":
    main()
