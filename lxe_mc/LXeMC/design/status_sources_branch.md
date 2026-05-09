# Sources Branch — Status (2026-05-09)

## Branch: `sources`

Pushed and synced with `origin/sources`. Based on `main` (which includes
all FastGeometry work).


## What's implemented

### Library code (src/)

| File | Content |
|:-----|:--------|
| `source_propagation.jl` | `propagate_in_source`, `random_position_in_volume`, `cos_theta_to_lxe` |
| `source_flux.jl` | `SourceFlux` (abstract), `SourceFluxBi214`, `SourceFluxTl208`, `SourceRateTable`, `SourceVolumeInfo`, `load_source_geometry`, `generate_flux_bi214/tl208`, `propagate_through_layers`, `generate_flux_compound_bi214/tl208` |
| `source_sampling.jl` | `sample_from_flux`, `sample_from_rate_table`, `sample_barrel_point`, `sample_cap_point`, `reconstruct_direction`, `sample_gamma_from_flux` |
| `source_dispatch.jl` | `VirtualEnvelope`, `make_virtual_envelope`, `make_surface_sampler`, `dispatch_source_flux`, `merge_dispatch_results`, `supported_sources/isotopes` |
| `cryostat_sources.jl` | `cryostat_barrel_flux`, `cryostat_top_flux`, `cryostat_bottom_flux`, `_build_rate_table`, `_get_activity_Bq` (with verbose progress) |
| `flux_utils.jl` | `merge_flux_bi214/tl208`, `write_pdf_csv`, `write_flux_bi214/tl208_csv`, `write_rate_table_csv`, `load_pdf_csv`, `load_flux_bi214/tl208`, `load_rate_table`, `write_flux_json`, metadata helpers |
| `geometry_core.jl` | Analytic `distance_to_entry` for `LDisk` (ray-ellipsoid intersection, replaces step search) |
| `tracking.jl` | `EventProcessingResult` now includes `deposits::Vector{Deposit}` for `:accepted` events |

### Scripts

| Script | Purpose |
|:-------|:--------|
| `scripts/generate_flux_tables.jl` | Generate (E, u) flux tables for one (source, isotope) pair. Parallel via `julia -t N`. Outputs CSV + metadata.json. |
| `scripts/run_source_backgrounds.jl` | Sample from flux tables, process through fast kernel + FV stack, collect SS candidates. Outputs `candidates.csv` + `statistics.csv`. Parallel. |

### Python plots

| Script | Purpose |
|:-------|:--------|
| `py/plot_fluxes.py` | 2D colormaps of flux tables. Bi-214: 2x2 grid. Tl-208: per-source panels (main + 3 companions). |
| `py/plot_candidates.py` | E, x, y, z histograms of SS candidates + rejection fractions table. |
| `py/compute_backgrounds.py` | Energy smearing, ROI cut, absolute background rates. Summary table with FV params, sigma, FWHM, rates. |

### Design documents

| File | Content |
|:-----|:--------|
| `design/implement_sources.md` | 8-point implementation plan (all points complete) |
| `design/fast_geometry_design.tex` | Fast kernel + staged transport design note |


## Key design decisions

1. **Factored flux tables**: (E, u) PDFs per source layer, not per event.
   Bi-214: one table. Tl-208: main + 3 independent companion tables with
   survival fractions. Companions sampled via independent Bernoulli trials.

2. **No cross-contributions**: barrel sources → barrel flux only. Top/bottom
   likewise. No OCV_barrel → ICV_top coupling.

3. **Compound propagation**: OCV → vacuum → ICV as sequential layers via
   `propagate_through_layers` using `distance_to_entry` on logical volumes.

4. **Virtual envelope**: gammas placed on detector boundary (from
   `FastKernelGeometry`), not on source surface. Avoids geometric mismatch
   between source geometry and tracking detector. Formalized as
   `VirtualEnvelope` struct with explicit, inspectable parameters.

5. **Three-stage pipeline**:
   - Stage 1: `generate_flux_tables.jl` → CSV flux tables (expensive MC, run once)
   - Stage 2: `run_source_backgrounds.jl` → SS candidates (fast kernel + FV stack)
   - Stage 3: `compute_backgrounds.py` → smearing + ROI + absolute rates

6. **E_max = line + 10 keV**: tight binning to avoid empty bins above line energy.
   E_min = 2.2 MeV for main gammas (ROI floor).

7. **Parallelization**: both Julia scripts use `julia -t N` with per-thread
   RNGs. Results merged via `merge_flux_bi214/tl208` or `merge_thread_results`.


## Cryostat source geometry (from source_geometry_lz_v1.json)

### Barrel flux: 3 components
- OCV_barrel (CylShell, R=90.8, t=0.7, Ti) → vacuum → ICV_barrel (R=82.1, t=0.9, Ti)
- MLI (transparent, 13.8 kg) → ICV_barrel
- ICV_barrel self

### Top flux: 2 components
- OCV_top (Disk, R=91.5, t=0.9, 2:1, Ti) → vacuum → ICV_top (R=83, t=0.8, 2:1, Ti)
- ICV_top self

### Bottom flux: 2 components
- OCV_bottom (Disk, R=91.5, t=1.5, 2:1, Ti) → vacuum → ICV_bottom (R=83, t=1.2, 3:1, Ti)
- ICV_bottom self


## Virtual envelopes (from tracking detector)

| Surface | Kind | R [cm] | z range / z_eq [cm] | Aspect ratio |
|:--------|:-----|:-------|:---------------------|:-------------|
| Barrel | barrel | 82.09 | [-41.33, 148.5] | — |
| Top | cap_up | 82.09 | z_eq = 145.6 | 2.0 |
| Bottom | cap_down | 82.09 | z_eq = -41.33 | 3.0 |


## Tested end-to-end

- **Barrel / Bi-214**: full pipeline works. 10^6 events in ~1.3s (8 threads).
  Survival fractions stable to 3 significant figures at 10^6.
  Background script produces candidates + statistics.

- **Top / Bi-214**: flux generation works (after LDisk analytic fix and
  virtual envelope fix). Background processing confirmed to produce
  vetoed events (not all no_fv).

- **Bottom / Bi-214**: flux generation works. Very well shielded (~7 mfp
  of passive LXe). Gammas land inside detector.


## Performance

| Operation | Rate |
|:----------|:-----|
| Flux generation (barrel, 8 threads) | 1.67M events/s |
| Flux generation (top, 8 threads) | ~750k events/s |
| Background processing | ~3M events/s |


## What's next

1. Run remaining 5 combinations: barrel/Tl208, top/Tl208, bottom/Bi214, bottom/Tl208, top/Bi214 (top already has flux tables)
2. Internal sources: FC_PTFE, FC_rings, PMTs, etc. (different from cryostat — emit directly into LXe)
3. Final background budget: combine all sources into total counts/year
4. Possible: energy resolution parameterization, ROI optimization
