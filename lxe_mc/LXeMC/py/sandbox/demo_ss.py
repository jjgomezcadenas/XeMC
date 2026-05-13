#!/usr/bin/env python3
"""
demo.py
=======

Demonstrates the ssms_classifier module:
  1. Prints the derived constants for each detector preset.
  2. Runs a handful of synthetic events to verify the algorithm.
  3. Times a 10⁶-event batch as a back-of-the-envelope speed check.

Run with:
    python demo.py
"""
from __future__ import annotations

import random
import time

import ssms_classifier as sc


# ============================================================
# 1. Print configuration tables
# ============================================================
print("\n" + "=" * 60)
print("  SS/MS classifier — algorithm parameters by detector")
print("=" * 60)

for det in (sc.XENON1T, sc.XENONnT, sc.LZ, sc.DARWIN):
    sc.print_config(det)


# ============================================================
# 2. Sanity-check synthetic events
# ============================================================
print("\n" + "=" * 60)
print("  Sanity checks: synthetic events")
print("=" * 60)

det = sc.LZ

# ---- Event 1: clean SS at Q_ββ in FV ----
ev1 = [sc.Deposit(0.0, 0.0, 0.50, 2458.0, in_fv=True)]
r1 = sc.classify(ev1, det)
print("\nEvent 1: single 2458 keV deposit at z=0.5m in FV")
print(f"  → SS={r1.is_ss}, MS={r1.is_ms}, in_ROI={r1.in_roi}")
assert r1.is_ss and not r1.is_ms and r1.in_roi

# ---- Event 2: Bi-214-like, two clusters well-separated ----
ev2 = [
    sc.Deposit(0.0, 0.0, 0.50, 2000.0, in_fv=True),
    sc.Deposit(0.0, 0.0, 0.52,  448.0, in_fv=True),  # Δz=20mm, E2=448 keV
]
r2 = sc.classify(ev2, det)
print("\nEvent 2: two clusters Δz=20mm, E1=2000, E2=448 keV (Bi-214-like)")
print(f"  → SS={r2.is_ss}, MS={r2.is_ms}")
assert r2.is_ms and not r2.is_ss

# ---- Event 3: two clusters too close in z (geometric fail) ----
ev3 = [
    sc.Deposit(0.0, 0.0, 0.500, 2000.0, in_fv=True),
    sc.Deposit(0.0, 0.0, 0.501,  458.0, in_fv=True),  # Δz=1mm < 3mm
]
r3 = sc.classify(ev3, det)
print("\nEvent 3: two clusters Δz=1mm — too close to resolve")
print(f"  → SS={r3.is_ss}, MS={r3.is_ms}")
assert r3.is_ss and not r3.is_ms

# ---- Event 4: small E2, fails κσ test (energetic fail) ----
ev4 = [
    sc.Deposit(0.0, 0.0, 0.50, 2400.0, in_fv=True),
    sc.Deposit(0.0, 0.0, 0.55,   30.0, in_fv=True),  # Δz=50mm OK, E2 too small
]
r4 = sc.classify(ev4, det)
print("\nEvent 4: Δz=50mm OK but E2=30 keV (< κσ ≈ 60 keV)")
print(f"  → SS={r4.is_ss}, MS={r4.is_ms}")
assert r4.is_ss and not r4.is_ms

# ---- Event 5: deposit outside FV above veto ----
ev5 = [
    sc.Deposit(0.0, 0.0, 0.50, 2400.0, in_fv=True),
    sc.Deposit(0.0, 0.0, 1.45,   50.0, in_fv=False),  # outside FV, > veto
]
r5 = sc.classify(ev5, det)
print("\nEvent 5: 50 keV deposit outside FV (> 10 keV veto)")
print(f"  → passed_veto={r5.passed_veto}")
assert not r5.passed_veto

# ---- Event 6: small deposit outside FV (below veto) survives ----
ev6 = [
    sc.Deposit(0.0, 0.0, 0.50, 2458.0, in_fv=True),
    sc.Deposit(0.0, 0.0, 1.45,    5.0, in_fv=False),  # outside FV, < veto
]
r6 = sc.classify(ev6, det)
print("\nEvent 6: 5 keV deposit outside FV (< 10 keV veto, survives)")
print(f"  → passed_veto={r6.passed_veto}, SS={r6.is_ss}")
assert r6.passed_veto and r6.is_ss

print("\nAll sanity checks PASSED.")


# ============================================================
# 3. Speed test on 10⁶ random events
# ============================================================
print("\n" + "=" * 60)
print("  Speed test: 10⁶ random events")
print("=" * 60)

rng = random.Random(42)
N = 1_000_000

events = []
for _ in range(N):
    z1 = 0.3 + 0.7 * rng.random()                 # uniform within FV span
    E1 = 2458.0 - 500.0 * rng.random()
    deposits = [sc.Deposit(0.0, 0.0, z1, E1, in_fv=True)]
    if rng.random() < 0.3:                        # 30% have a satellite
        sign = 1 if rng.random() < 0.5 else -1
        dz = sign * (0.001 + 0.05 * rng.random()) # 1–51 mm offset
        E2 = 50.0 + 500.0 * rng.random()
        z2 = z1 + dz
        deposits.append(sc.Deposit(0.0, 0.0, z2, E2,
                                   in_fv=(0.30 < z2 < 1.00)))
    events.append(deposits)

# Warm up
sc.classify_batch(events[:1000], det)

t0 = time.perf_counter()
results = sc.classify_batch(events, det)
elapsed = time.perf_counter() - t0

s = sc.summarize(results)
print()
print(f"  Events processed     : {s.n_total:,}")
print(f"  Passed FV veto       : {s.n_passed_veto:,} "
      f"({100 * s.n_passed_veto / s.n_total:.1f}%)")
print(f"  In ROI               : {s.n_in_roi:,} "
      f"({100 * s.n_in_roi / s.n_total:.1f}%)")
print(f"  In ROI, tagged SS    : {s.n_ss_in_roi:,}")
print(f"  In ROI, tagged MS    : {s.n_ms_in_roi:,}")
print(f"  MS tagging eff (ROI) : {100 * s.ms_tag_eff:.1f}%")
print()
print(f"  Wall time            : {elapsed:.2f} s")
print(f"  Throughput           : {N / 1e6 / elapsed:.2f} Mevent/s")
print(f"  Projected for 10⁸    : {1e8 / (N / elapsed):.0f} s")
