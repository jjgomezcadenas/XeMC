# Tracking and Transport Architecture

Compact reference for the LXeMC transport code. Context-recovery document.


## File organization

| File | Lines | Content |
|:-----|:------|:--------|
| `tracking.jl` | ~360 | Full-stack physics: `transport_photon!`, `transport_lepton!`, `propagate_gamma`, `_classify_escape_volume` |
| `tracking_fast.jl` | ~200 | Fast-kernel transport: `transport_gamma_fastkernel` + helpers |
| `event.jl` | ~80 | Event processing: `process_event`, `EventProcessingResult` |
| `tracking_util.jl` | ~160 | Diagnostic tools: `propagate_gamma_fastkernel`, `process_event_fastkernel_calib` |
| `analysis.jl` | ~60 | Offline analysis: `cluster_deposits_in_z`, `is_single_site` |


## Data types

### tracking.jl

- `Track(kind, energy, position, direction, track_id, parent_id, generation, interaction, volume)` — one particle on the stack
- `ParticleStack` — struct wrapping `Vector{Track}` (LIFO stack)
- `Deposit(position, energy, source, interaction, volume)` — one energy deposit; `volume` is `:fv`, `:active`, or `:passive`
- `SampledGamma(E_MeV, position, direction)` — gamma ready for injection

### tracking_fast.jl

- `FastGammaDeposit(Edep_MeV, position, region)` — deposit during fast transport (internal bookkeeping, not returned to caller)
- `FastGammaTrackResult(status, deposits, terminal_region, energy_MeV, position, direction)` — result of fast-kernel transport

### event.jl

- `EventProcessingResult(status, has_tpc_veto, has_skin_veto, n_processed, deposits)`
  - `:fv` — gamma reached FV; `deposits` has full stack for offline SS/MS
  - `:vetoed` — veto threshold exceeded
  - `:no_fv` — no gamma reached FV

### tracking_util.jl (diagnostic only)

- `GammaPropagationResult(status, interaction_type, deposit_E_MeV, position, region)` — single-step result
- `FastKernelCalibEventResult(result, multiplicity, E1_MeV, E2_MeV)` — calibration wrapper


## Transport flow

```
process_event(gammas, fk, vol, cfg, rng)
  for each gamma:
    result = transport_gamma_fastkernel(gamma, fk, cfg, rng)
    |
    +-- :handoff_fv --> propagate_gamma(E, vol, fk, cfg; pos, dir, rng)
    |                     |
    |                     +-- transport_photon!(track, vol, deposits, stack, ...)
    |                     +-- transport_lepton!(track, vol, deposits, stack, ...)
    |                     +-- (loop until stack empty)
    |                     +-- returns Vector{Deposit}
    |                   append deposits to fv_deposits
    |
    +-- :vetoed_tpc --> return EventProcessingResult(:vetoed, ...)
    +-- :vetoed_skin --> return EventProcessingResult(:vetoed, ...)
    +-- other --> continue (gamma absorbed/escaped/below cut)

  if fv_deposits not empty:
    return EventProcessingResult(:fv, ..., fv_deposits)
  else:
    return EventProcessingResult(:no_fv, ...)
```


## Full-stack transport (tracking.jl)

### transport_photon!(track, vol, fk, deposits, stack, track_counter, cfg, rng)

Loop while `E >= Egamma_cut` (10 keV):
1. Sample interaction distance from total cross section
2. Advance; exit if outside `vol` (`is_inside(vol, pos)`)
3. Sample process (Compton / pair / photoelectric)
4. **Compton**: scatter gamma, push electron to stack
5. **Pair**: push positron + electron to stack, gamma dies
6. **Photoelectric**: deposit EK locally, push photoelectron if `T_e > Te_cut`

### transport_lepton!(track, vol, fk, deposits, stack, track_counter, cfg, rng)

Condensed-history stepping while `T >= Te_cut` (400 keV):
1. Fixed step `ds_step`, clamped so `dE < T`
2. Collisional dE/dx deposited at step midpoint
3. Exit if outside `vol`
4. Discrete bremsstrahlung: per-step probability, hard photon pushed to stack
5. End of range: residual T deposited locally
6. Positron annihilation: two 511 keV back-to-back gammas pushed to stack

### propagate_gamma(E, vol, fk, cfg; position, direction, rng)

`fk` is used by `_classify_escape_volume` to classify deposits that escape the FV boundary as `:active`, `:passive`, or `:fv`.

Full cascade loop:
1. Push initial gamma to stack
2. Pop particle, dispatch to `transport_photon!` or `transport_lepton!`
3. Repeat until stack empty or generation cap reached
4. Return all deposits


## Fast-kernel transport (tracking_fast.jl)

### transport_gamma_fastkernel(gamma, fk, cfg, rng)

Loop while `E >= Egamma_cut` and `traveled < max_cm`:
1. **Classify region**: `classify_fastkernel(fk, pos)` returns a `FastKernelRegion` or `nothing`
   - `nothing` --> `:escaped`
   - `region.name == "FV"` --> `:handoff_fv` (hand off to full stack)
2. **Vacuum/non-interacting**: advance to boundary, continue
3. **Sample interaction**: compare `s_int` vs `s_bnd`
   - Beyond boundary: advance to boundary, continue
   - Inside region: interact
4. **Compton**: deposit recoil, check veto thresholds
   - Skin above `veto_skin` --> `:vetoed_skin`
   - TopActive/BarrelActive/BottomActive above `veto_TPC` --> `:vetoed_tpc`
   - Below ROI floor (`E < 2.3 MeV`): collapse remaining energy, return terminal status
   - Otherwise: update direction, continue
5. **Pair/photoelectric**: terminal, deposit full energy, return status by region

### Passive region check

`_is_passive_region(name)` returns true for: `LXe_dome`, `LXe_below_FC`, `LXe_below_cathode`, `FC_PTFE`, `FC_rings`. These regions absorb energy without triggering a veto.

### Fast-kernel status values

| Status | Meaning |
|:-------|:--------|
| `:handoff_fv` | Reached FV; hand off to `propagate_gamma` |
| `:escaped` | Left the detector envelope |
| `:below_cut` | Energy fell below `Egamma_cut` |
| `:absorbed_passive` | Absorbed in passive LXe or structural material |
| `:below_roi_fv` | ROI-floor termination after Compton scatter |
| `:vetoed_tpc` | Deposit in active TPC outside FV exceeded `veto_TPC` |
| `:vetoed_skin` | Deposit in Skin exceeded `veto_skin` |

### Veto thresholds (from sim_config.json)

- `veto_skin`: Skin deposit threshold
- `veto_TPC`: TPC active deposit threshold

### ROI-floor optimization

If post-Compton energy < `E_roi_floor` (2.3 MeV), the gamma can no longer produce
a signal in the ROI (Q_bb = 2.458 MeV). Remaining energy is collapsed onto the
last deposit. This is an intentional approximation — it overestimates local deposits
but is safe for background estimation.


## Physics parameters

| Parameter | Value | Source |
|:----------|:------|:-------|
| `Egamma_cut` | 10 keV | Photon tracking cutoff |
| `Te_cut` | 400 keV | Lepton tracking cutoff |
| `E_roi_floor` | 2.3 MeV | Fast-kernel ROI floor |
| `EK` (Xe) | 34.561 keV | K-shell binding energy |
| `ds_step` | from config | Lepton condensed-history step |
| Cross sections | XCOM (photon), ESTAR (electron) | NIST tables, log-log interpolation |


## Key design decisions

1. **No SS/MS in event loop**: `process_event` returns raw FV deposits. SS/MS
   classification is done offline (python or Julia) from the deposit CSV.

2. **FVGeometry eliminated**: `compile_fv_volume` returns `PCyl` (a `PhysicalVolume`).
   All transport functions take `PhysicalVolume` + `FastKernelGeometry`. The `fk`
   parameter is used by `_classify_escape_volume` to classify deposits that escape
   the FV boundary as `:active`, `:passive`, or `:fv`.

3. **VirtualEnvelope from source geometry**: `make_virtual_envelope(source, sg)` for
   all sources. No `fk` dependency. Geometry fix ensures ICV inner surfaces match
   tracking detector boundary at R=82.1.

4. **Diagnostic vs production**: `propagate_gamma_fastkernel` (single-step, in
   `tracking_util.jl`) is for validation only. Production uses
   `transport_gamma_fastkernel` (multi-scatter, in `tracking_fast.jl`).

5. **Fast rejection**: veto decisions are made inside `transport_gamma_fastkernel`
   during transport. `process_event` just checks `result.status`.

6. **Lepton range approximation**: electrons/positrons are transported with
   condensed-history stepping (no multiple scattering). Valid because electron
   range in LXe is sub-mm at MeV energies and the FV analysis doesn't need
   sub-mm spatial precision.
