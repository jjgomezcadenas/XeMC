# Tracking and Event Loop — Status (2026-05-10)

## Branch: `offline_fv_analysis`

Pushed and synced with `origin/offline_fv_analysis`. Tip: `fb93b32`
(Rewrite lxemc.tex with three-tier accuracy framing).


## What changed in this round

The transport code went from a single monolithic `tracking.jl` to a
five-file split that mirrors the three-tier accuracy model now documented
in `lxemc.tex`. Each file has a single responsibility and a clear
production-vs-diagnostic boundary.

| Commit | Change |
|:-------|:-------|
| `fbaa09b` | Move `cluster_deposits_in_z` and `is_single_site` to `analysis.jl` |
| `0c4be7c` | Move `propagate_gamma_fastkernel` (single-step diagnostic) to `tracking_util.jl` |
| `c334317` | Move `FastKernelCalibEventResult` and `process_event_fastkernel_calib` to `tracking_util.jl` |
| `7aaca36` | Split `tracking.jl` → `tracking.jl` + `tracking_fast.jl` + `event.jl` |
| `7f5ab71` | Add `tracking_and_transport.md`; rename design doc to `lxemc.tex` |
| `fb93b32` | Rewrite `lxemc.tex` around three-tier framing + new Physics section |


## Current file layout

| File | Lines | Responsibility |
|:-----|------:|:---------------|
| `tracking.jl`       | ~300 | Full-stack physics (Tier 1, FV-only): `Track`, `ParticleStack`, `Deposit`, `transport_photon!`, `transport_lepton!`, `propagate_gamma` |
| `tracking_fast.jl`  | ~200 | Fast-kernel transport (Tier 3): `FastGammaDeposit`, `FastGammaTrackResult`, `transport_gamma_fastkernel` |
| `event.jl`          |  ~80 | Event orchestration: `EventProcessingResult`, `process_event` |
| `tracking_util.jl`  | ~160 | Diagnostics only: single-step `propagate_gamma_fastkernel`, `FastKernelCalibEventResult`, `process_event_fastkernel_calib` |
| `analysis.jl`       |  ~60 | Offline helpers: `cluster_deposits_in_z`, `is_single_site` |

`include` order in `LXeMC.jl`: `tracking.jl` → `tracking_fast.jl` →
`event.jl` → `tracking_util.jl`. Diagnostics depend on production, not
the other way around.


## Production data flow

```
process_event(gammas, fk, vol, cfg, rng)            # event.jl
  for each gamma:
    result = transport_gamma_fastkernel(gamma, ...) # tracking_fast.jl
    |
    +-- :handoff_fv  --> propagate_gamma(...)       # tracking.jl
    |                      transport_photon! / transport_lepton!
    |                      append deposits to fv_deposits
    +-- :vetoed_tpc  --> short-circuit, return :vetoed
    +-- :vetoed_skin --> short-circuit, return :vetoed
    +-- other        --> continue (gamma terminated, no veto)
  return :fv | :no_fv | :vetoed
```

`transport_gamma_fastkernel` applies the two analysis-driven veto cuts
(`veto_TPC`, `veto_skin`) and the ROI floor (`E_roi_floor`) during
transport, so `process_event` only inspects the returned status.


## Status enums (stable, do not redefine elsewhere)

`FastGammaTrackResult.status` (in `tracking_fast.jl`):

| Value | Meaning |
|:------|:--------|
| `:handoff_fv`       | Reached FV; hand off to `propagate_gamma` |
| `:escaped`          | Left detector envelope |
| `:below_cut`        | Energy fell below `Egamma_cut` |
| `:absorbed_passive` | Absorbed in passive LXe / structural |
| `:below_roi_fv`     | ROI-floor termination after Compton |
| `:vetoed_tpc`       | Deposit in TPC active ≥ `veto_TPC` |
| `:vetoed_skin`      | Deposit in Skin ≥ `veto_skin` |

`EventProcessingResult.status` (in `event.jl`):

| Value     | Meaning |
|:----------|:--------|
| `:fv`     | At least one γ reached FV; `deposits` is the full per-deposit truth list |
| `:vetoed` | A `:vetoed_tpc` or `:vetoed_skin` short-circuited the event |
| `:no_fv`  | No γ reached FV (escaped, absorbed, or ROI-floor terminated) |


## Key design decisions

1. **Three accuracy tiers** drive the file split.
   Tier 1 (FV, full stack) lives in `tracking.jl`; Tier 2 (source-region
   photon-only) lives in `source_propagation.jl` (covered by
   `status_sources.md`); Tier 3 (fast vetoed transport) lives in
   `tracking_fast.jl`. `event.jl` is the orchestrator that ties Tiers 1
   and 3 together for detector-side events.

2. **Production vs. diagnostic, in separate files.**
   `transport_gamma_fastkernel` (multi-scatter, production) lives in
   `tracking_fast.jl`. `propagate_gamma_fastkernel` (single-step,
   diagnostic) lives in `tracking_util.jl`. Same for the calibration
   wrapper. This is enforced by file ownership, not just by docstring.

3. **Veto checks happen during fast-kernel transport, not after.**
   `transport_gamma_fastkernel` returns `:vetoed_tpc` / `:vetoed_skin`
   the moment a Compton deposit crosses threshold. `process_event` then
   short-circuits the rest of the event. There is no post-hoc "is this
   event vetoed?" pass over the deposit list.

4. **`veto_TPC` is a tunable parameter, not a fixed constant.**
   Default 10 keV, recommended range 5–10 keV for 0νββ analyses
   (XENONnT-style). The 1 keV thresholds quoted by DM searches would
   risk vetoing 0νββ signal whose resolution-broadened tail leaks into
   the active TPC outside the FV. `veto_skin` (default 100 keV) is
   LZ-specific and tracks the LZ outer-detector analysis.

5. **`process_event` returns raw deposits; SS/MS classification is
   offline.** No clustering, no resolution smearing, no ROI cuts during
   transport. The deposit list (`Vector{Deposit}` with `position`,
   `energy`, `source` symbol) is the handoff format to the offline FV
   analysis. `cluster_deposits_in_z` and `is_single_site` are convenience
   helpers in `analysis.jl`, not part of the event loop.

6. **FastGamma deposits are internal bookkeeping.**
   `FastGammaDeposit` records during fast-kernel transport drive the veto
   threshold check inside `transport_gamma_fastkernel`. They are *not*
   surfaced in `EventProcessingResult.deposits`; only Tier-1 FV deposits
   are returned to the caller, because only those are relevant to the
   FV-side topology classification.

7. **`PhysicalVolume`, not `FVGeometry`.**
   `compile_fv_volume(det)` returns a `PCyl` (a `PhysicalVolume`).
   `propagate_gamma` and `transport_photon!`/`transport_lepton!` accept
   any `PhysicalVolume`. The previous `FVGeometry` wrapper was removed:
   the FV is just a cylindrical active LXe volume, not a special type.

8. **Stack model is LIFO with a generation cap.**
   `ParticleStack` is a mutable `Vector{Track}`; secondaries are pushed
   during a step and resumed by `propagate_gamma` after the parent
   returns. `cfg.generation_cap = 100` discards anything deeper as a
   safety net (no runaway cascades observed, but the cap guards against
   bugs in future physics extensions).


## What's stable

- The five-file split. No further moves expected unless a new tier is
  added (e.g., a GXe full-MS transport).
- The two status enums (`FastGammaTrackResult.status`,
  `EventProcessingResult.status`). Downstream code reads these by symbol;
  changes here are breaking.
- The deposit output format (`Vector{Deposit}`). The offline analysis on
  this branch consumes it directly.
- Configuration parameter names in `sim_config.json` and `SimConfig`.


## What's pending

- **GXe physics extensions.** Multiple Coulomb scattering, δ-rays,
  in-flight positron annihilation kinematics, atomic relaxation cascade
  on photoelectric absorption, and possibly Rayleigh scattering. Needed
  if/when a high-pressure GXe configuration is added; not blocking for
  the current LXe FV studies. See `lxemc.tex` §9 (Limitations).
- **Richer SS/MS classification.** Current `is_single_site` is a
  one-dimensional z-clustering with a fixed `dz_resolution`. A real
  classifier will need 3D clustering plus a detector-response stage.
  The deposit-list output is designed to make this evolution possible
  without touching the transport code.
- **Energy-resolution smearing and ROI cuts** still live in the Python
  stage (`compute_backgrounds.py`). Could be folded into Julia
  post-processing once the resolution model stabilizes.
- **`E_roi_floor` configurability per analysis.** Default 2.3 MeV is
  tuned to the 0νββ ROI; analyses with wider windows or longer
  resolution tails should set their own.


## Cross-reference

- Prose architecture: `lxemc.tex` §1 (Introduction, three tiers),
  §4 (Full-stack), §5 (Fast kernel), §6 (Event processing).
- Compact code-level reference: `tracking_and_transport.md` (data
  types, transport flow ASCII diagram, parameter table).
- Source pipeline status: `status_sources.md`.
