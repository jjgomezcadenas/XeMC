#!/usr/bin/env julia
"""
demo.jl

Demonstrates the SSMSClassifier module:
  1. Prints the derived constants for each detector preset.
  2. Runs a handful of synthetic events to verify the algorithm.
  3. Times a 10⁶-event batch as a back-of-the-envelope speed check.
"""

include("SSMSClassifier.jl")
using .SSMSClassifier
using Printf
using Random

const SC = SSMSClassifier

println("\n" * "="^60)
println("  SS/MS classifier — algorithm parameters by detector")
println("="^60)

for det in (SC.XENON1T, SC.XENONnT, SC.LZ, SC.DARWIN)
    SC.print_config(det)
end

println("\n" * "="^60)
println("  Sanity checks: synthetic events")
println("="^60)

# Use LZ for the worked examples
det = SC.LZ

# ---- Event 1: clean SS at Q_ββ in FV ----
ev1 = [SC.Deposit(0.0, 0.0, 0.50, 2458.0, true)]
r1 = SC.classify(ev1, det)
println("\nEvent 1: single 2458 keV deposit at z=0.5m in FV")
println("  → SS=$(r1.is_ss), MS=$(r1.is_ms), in_ROI=$(r1.in_roi)")
@assert r1.is_ss && !r1.is_ms && r1.in_roi

# ---- Event 2: Bi-214-like, two clusters well-separated ----
ev2 = [
    SC.Deposit(0.0, 0.0, 0.50, 2000.0, true),
    SC.Deposit(0.0, 0.0, 0.52, 448.0,  true),  # Δz = 20 mm, E2 = 448 keV
]
r2 = SC.classify(ev2, det)
println("\nEvent 2: two clusters Δz=20mm, E1=2000, E2=448 keV (Bi-214-like)")
println("  → SS=$(r2.is_ss), MS=$(r2.is_ms)")
@assert r2.is_ms && !r2.is_ss

# ---- Event 3: two clusters too close in z (geometric fail) ----
ev3 = [
    SC.Deposit(0.0, 0.0, 0.500, 2000.0, true),
    SC.Deposit(0.0, 0.0, 0.501, 458.0,  true),  # Δz = 1 mm < 3 mm
]
r3 = SC.classify(ev3, det)
println("\nEvent 3: two clusters Δz=1mm — too close to resolve")
println("  → SS=$(r3.is_ss), MS=$(r3.is_ms)")
@assert r3.is_ss && !r3.is_ms

# ---- Event 4: small E2, fails κσ test (energetic fail) ----
ev4 = [
    SC.Deposit(0.0, 0.0, 0.50, 2400.0, true),
    SC.Deposit(0.0, 0.0, 0.55, 30.0,   true),  # Δz=50mm OK, E2=30 keV < 60 keV
]
r4 = SC.classify(ev4, det)
println("\nEvent 4: Δz=50mm OK but E2=30 keV (< κσ ≈ 60 keV)")
println("  → SS=$(r4.is_ss), MS=$(r4.is_ms)")
@assert r4.is_ss && !r4.is_ms

# ---- Event 5: deposit outside FV above veto ----
ev5 = [
    SC.Deposit(0.0, 0.0, 0.50,  2400.0, true),
    SC.Deposit(0.0, 0.0, 1.45,    50.0, false), # outside FV, E > E_VETO
]
r5 = SC.classify(ev5, det)
println("\nEvent 5: 50 keV deposit outside FV (> 10 keV veto)")
println("  → passed_veto=$(r5.passed_veto)")
@assert !r5.passed_veto

# ---- Event 6: small deposit outside FV (below veto) survives ----
ev6 = [
    SC.Deposit(0.0, 0.0, 0.50, 2458.0, true),
    SC.Deposit(0.0, 0.0, 1.45,    5.0, false), # outside FV, E < E_VETO
]
r6 = SC.classify(ev6, det)
println("\nEvent 6: 5 keV deposit outside FV (< 10 keV veto, survives)")
println("  → passed_veto=$(r6.passed_veto), SS=$(r6.is_ss)")
@assert r6.passed_veto && r6.is_ss

println("\nAll sanity checks PASSED.")

# ---- Speed test ----
println("\n" * "="^60)
println("  Speed test: 10⁶ random events")
println("="^60)

using Random
rng = MersenneTwister(42)

# Generate 10⁶ synthetic events: ~30% have a Compton satellite
N = 1_000_000
events = Vector{Vector{SC.Deposit}}(undef, N)
for i in 1:N
    z1 = 0.3 + 0.7 * rand(rng)             # uniform in FV
    E1 = 2458.0 - 500.0 * rand(rng)
    deposits = [SC.Deposit(0.0, 0.0, z1, E1, true)]
    if rand(rng) < 0.3                     # 30% MS
        Δz = (-1)^rand(rng, 0:1) * (0.001 + 0.05 * rand(rng))
        E2 = 50.0 + 500.0 * rand(rng)
        push!(deposits,
              SC.Deposit(0.0, 0.0, z1 + Δz, E2,
                         0.30 < z1 + Δz < 1.00))
    end
    events[i] = deposits
end

# Warm up
SC.classify_batch(events[1:1000], det)

# Timed run
t0 = time()
results = SC.classify_batch(events, det)
elapsed = time() - t0

s = SC.summary(results)
println()
@printf("  Events processed     : %d\n", s.n_total)
@printf("  Passed FV veto       : %d (%.1f%%)\n",
        s.n_passed_veto, 100 * s.n_passed_veto / s.n_total)
@printf("  In ROI               : %d (%.1f%%)\n",
        s.n_in_roi, 100 * s.n_in_roi / s.n_total)
@printf("  In ROI, tagged SS    : %d\n", s.n_ss_in_roi)
@printf("  In ROI, tagged MS    : %d\n", s.n_ms_in_roi)
@printf("  MS tagging eff (ROI) : %.1f%%\n", 100 * s.ms_tag_eff)
println()
@printf("  Wall time            : %.2f s\n", elapsed)
@printf("  Throughput           : %.2f Mevent/s\n", N / 1e6 / elapsed)
@printf("  Projected for 10⁸    : %.0f s\n", 1e8 / (N / elapsed))
