# Source Flux Implementation Plan

## Completed

All infrastructure is in place. The three-stage pipeline works end-to-end
for 5 source types x 2 isotopes = 10 combinations.

### Done: dense sources (cryostat)

- Data structures: `SourceFluxBi214`, `SourceFluxTl208`, `SourceRateTable`
- Source geometry loader: `load_source_geometry`
- Single-source flux: `generate_flux_bi214`, `generate_flux_tl208`
- Compound propagation: `propagate_through_layers`, `generate_flux_compound_*`
- Cryostat flux functions: `cryostat_barrel_flux`, `cryostat_top_flux`,
  `cryostat_bottom_flux` (in `cryostat_sources.jl`)
- MLI fix: compound propagation MLI -> ICV (not generated inside ICV)
- Sampling: `sample_from_flux`, `sample_from_rate_table`
- Surface samplers: `sample_barrel_point`, `sample_cap_point`,
  `sample_disk_point`, `reconstruct_direction`, `sample_gamma_from_flux`
- Flux I/O: CSV + JSON metadata, merge utilities
- Dispatch: `dispatch_source_flux`, `merge_dispatch_results`

### Done: transparent sources (PMTs)

- `pmt_sources.jl`: `pmt_top_flux`, `pmt_bottom_flux`
- `_merge_pmt_volume`: builds merged PDisk from sub-components
- Hemisphere filtering via PDisk orientation (:up/:down)
- Per-component rate bookkeeping (mass x activity)
- `sample_disk_point` with normal_z parameter

### Done: geometry fix

- ICV inner surfaces unified to R_inner = 82.1 (mass-preserving)
- AirDome moved to ICV_top equator (z=148.5), AirCyl added
- VirtualEnvelope unified on source geometry (`sg`), `fk` method removed
- `distance_to_exit(LDisk)` analytic for flat disks
- Vacuum handling in `propagate_in_source` (density <= 0)

### Done: scripts and analysis

- `generate_flux_tables.jl` (Stage 1)
- `run_source_backgrounds.jl` (Stage 2)
- `compute_backgrounds.py` (Stage 3)
- `plot_fluxes.py`, `plot_candidates.py`
- `running_recipe.md` with all 10 pipeline commands


## Remaining: PMT barrel sources

Three transparent shell sources along the FC outer wall:

| Source | R_inner (cm) | t (cm) | Mass (kg) | Bi214 (mBq/kg) | Tl208 (mBq/kg) |
|:-------|:-------------|:-------|:----------|:----------------|:----------------|
| PMT_BARREL_cables | 81.2 | 0.150 | 88.7 | 4.31 | 0.82 |
| PMT_BARREL_R8520 | 81.4 | 0.010 | 6.10 | 46.0 | 14.9 |
| PMT_BARREL_R8778_lower | 81.4 | 0.002 | 1.31 | 46.0 | 14.9 |

These are cylinder_shell sources (not flat disks like PMT_TOP/BOT).
All transparent (Vacuum, no self-shielding).

### Implementation plan

1. Add `pmt_barrel_flux` to `pmt_sources.jl`:
   - Merge the 3 components into a single transparent CylShell
   - Use existing `sample_barrel_point` for the VE (`:barrel` kind)
   - Per-component rate bookkeeping like PMT_TOP/BOT

2. Add `"pmt_barrel"` to `dispatch_source_flux` and `supported_sources`

3. Add `make_virtual_envelope` case for `pmt_barrel` in `source_dispatch.jl`

4. Run through full pipeline (flux tables + backgrounds + analysis)


## Remaining: field cage internal sources

Six sources inside or on the field cage, emitting directly into LXe:

| Source | Shape | Material | Transport | Mass (kg) | Bi214 | Tl208 |
|:-------|:------|:---------|:----------|:----------|:------|:------|
| FC_PTFE | cyl_shell | PTFE | KN | 184 (exact) | 0.04 | 0.01 |
| FC_rings | cyl_shell | Ti | KN | 93 (exact) | 0.35 | 0.24 |
| FC_resistors | cyl_shell | Vacuum | transparent | 0.06 | 1350 | 2010 |
| FC_sensors | cyl_shell | Vacuum | transparent | 5.02 | 5.82 | 1.88 |
| FC_topgrid | cyl_shell | SS | KN | 44.55 | 2.63 | 1.46 |
| FC_botgrid | cyl_shell | SS | KN | 22.28 | 2.63 | 1.46 |

### Key differences from cryostat/PMT sources

- FC_PTFE and FC_rings are **dense** sources with KN self-shielding,
  but emit directly into LXe (no compound propagation through layers)
- FC_topgrid and FC_botgrid are stainless steel, also KN
- FC_resistors and FC_sensors are transparent (vacuum)
- All are cylinder_shells at R ~ 72.6-80.3
- The VE is the source inner surface (barrel type)

### Implementation plan

1. Group by transport type:
   - KN sources (FC_PTFE, FC_rings, FC_topgrid, FC_botgrid): use
     `generate_flux_bi214/tl208` directly (single-source, no compound)
   - Transparent sources (FC_resistors, FC_sensors): same code path,
     vacuum propagation gives ~50% survival

2. Decision: treat individually or lump?
   - Each source has different R, z, material, activity
   - Probably treat individually (6 separate dispatch entries)
   - Or group by z-region: "fc_barrel" (PTFE, rings, resistors, sensors)
     vs "fc_topgrid" vs "fc_botgrid"

3. Add dispatch entries and VE support

4. Run through full pipeline


## Remaining: final background budget

After all sources are implemented:

1. Run all (source, isotope) combinations with high statistics
2. Collect `background_summary.json` from each
3. Sum all contributions into total counts/year per isotope
4. Compare with published LZ background budget
