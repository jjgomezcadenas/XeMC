# Optimization History

## Measurement caveats

Early benchmarks (≤5k events from CLI) include Julia JIT compilation
overhead, which dominates for small N. Numbers below 50k events from
a cold session are unreliable. Warm-loop measurements (200k+ events,
post-warmup) are marked with *.

## Algorithm changes (affect precision/physics tradeoff)

| Change | Te_cut | k_min | Step | Events/sec | Notes |
|:-------|:-------|:------|:-----|:-----------|:------|
| Original (fraction-of-range, per-step brems integration) | 50 keV | 50 keV | 5% of range | ~655 (cold, 100 events) | JIT-dominated |
| Pre-tabulated σ_brems(T) + fixed 0.5 mm step | 50 keV | 50 keV | 0.5 mm | ~7,200 (cold, 5k) | Eliminated 100-point integration per step |
| Fixed 1 mm step | 50 keV | 50 keV | 1.0 mm | ~25,000 (cold, 5k) | Range at 400 keV = 0.74 mm, adequate for 3 mm resolution |
| Raised Te_cut to 100 keV | 100 keV | 100 keV | 1.0 mm | ~26,000 (cold, 5k) | |
| Te_cut scan: 400 keV optimal | 400 keV | 50 keV | 1.0 mm | ~46,000 (cold, 5k) | CSDA range at 400 keV = 0.74 mm ≪ 3 mm z-resolution |
| Simplified photoelectric (deposit EK locally) | 400 keV | 50 keV | 1.0 mm | ~26,000 (cold, 5k) | No speed gain (photo is 2% of events), but exact energy conservation |

**Key physics parameters (current):**
- Te_cut = 400 keV (CSDA range 0.74 mm)
- k_min = 50 keV (photon mfp 0.27 mm)
- Egamma_cut = 10 keV (mfp 0.027 mm)
- Step = 1 mm fixed
- All sub-cut: deposit locally. SS/MS classification at 3 mm z-resolution, 10 keV cluster threshold.

## Implementation optimizations (no physics change)

Measured with warm loop (200k–500k events, post-warmup).

| Change | Events/sec* | Speedup |
|:-------|:------------|:--------|
| Baseline (after algorithm changes + geometry refactor) | 162,000 | 1× |
| Pre-tabulated brems rejection ceiling M(T) | 220,000 | 1.36× |
| Combined sigma_three (1 binary search for 3 XCOM channels) | — | — |
| Pre-logged XCOM/ESTAR/brems grids (avoid log/exp on grid points) | — | — |
| Replace objectid with integer track counter | — | — |
| **All three combined** | **293,000** | **1.81×** |

## Current profile (500k events, 2.615 MeV in LZ TPC)

| Function | % CPU | Notes |
|:---------|:------|:------|
| transport_photon! | 59% | Cross-section lookups + Compton kinematics |
| ├ sigma_three (XCOM interp) | 22% | One search, three interpolations |
| ├ compton_electron_direction | 11% | Vector allocations |
| ├ sample_compton | 4% | KN rejection sampling |
| transport_lepton! | 38% | Stepping + brems |
| ├ dEdx_collision (ESTAR interp) | 7% | One lookup per step |
| ├ sigma_brems (table lookup) | 6% | One lookup per step |
| ├ brems sampling + stepping | 21% | |
| objectid (now removed) | 4% | Replaced by counter |

## Remaining opportunities

- **StaticArrays SVector{3,Float64}**: would eliminate ~9 KB/event of
  heap allocations from `Vector{Float64}` temporaries in position/direction
  vectors. Estimated 20-30% further gain. Requires changing Track, Deposit,
  and all transport functions.

## Performance summary

At 293k events/sec (warm), 10⁸ events takes ~5.7 minutes.
Photon-only mode runs at ~400k+ events/sec.
