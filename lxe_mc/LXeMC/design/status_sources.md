# Sources Branch — Status (2026-05-13)

## Branch: `sources`

Pushed and synced with `origin/sources`.


## What's implemented

### Library code (src/)

| File | Content |
|:-----|:--------|
| `source_propagation.jl` | `propagate_in_source` (handles vacuum: straight-line exit for density <= 0), `random_position_in_volume` (PCyl, PCylShell, PDisk), `cos_theta_to_lxe` |
| `source_flux.jl` | `SourceFlux` (abstract), `SourceFluxBi214`, `SourceFluxTl208`, `SourceRateTable`, `SourceVolumeInfo`, `load_source_geometry`, `generate_flux_bi214/tl208`, `propagate_through_layers`, `generate_flux_compound_bi214/tl208` |
| `source_sampling.jl` | `sample_from_flux`, `sample_from_rate_table`, `sample_barrel_point`, `sample_disk_point(R, z, normal_z, rng)`, `sample_cap_point`, `reconstruct_direction`, `sample_gamma_from_flux` |
| `source_dispatch.jl` | `VirtualEnvelope` (kinds: `:barrel`, `:cap_up`, `:cap_down`, `:disk_flat`), `make_virtual_envelope(source, sg)` (unified, all sources from source geometry), `make_surface_sampler`, `dispatch_source_flux`, `merge_dispatch_results`, `supported_sources/isotopes` |
| `cryostat_sources.jl` | `cryostat_barrel_flux` (OCV compound + MLI compound + ICV self), `cryostat_top_flux`, `cryostat_bottom_flux`, `_build_rate_table`, `_get_activity_Bq` |
| `pmt_sources.jl` | `pmt_top_flux`, `pmt_bottom_flux`, `pmt_bottom_lxe_flux` (compound: PMT + passive LXe to cathode), `pmt_barrel_flux` (cables, R8520, R8778 lower-ring as merged PCylShell), `_merge_pmt_volume` |
| `flux_utils.jl` | `merge_flux_bi214/tl208`, CSV I/O, JSON metadata helpers |
| `geometry_core.jl` | Analytic `distance_to_entry(LDisk)`, `distance_to_exit(LDisk)` for flat disks |
| `tracking.jl` | `EventProcessingResult` includes `deposits::Vector{Deposit}` for `:accepted` events |

### Scripts

| Script | Purpose |
|:-------|:--------|
| `scripts/generate_flux_tables.jl` | Stage 1: generate (E, u) flux tables. Parallel via `julia -t N`. Works for all sources. |
| `scripts/run_source_backgrounds.jl` | Stage 2: sample from flux tables, fast kernel + FV stack, SS candidates. Uses `sg`-based sampler. |

### Python analysis (py/)

| Script | Purpose |
|:-------|:--------|
| `plot_fluxes.py` | 2D colormaps of flux tables |
| `plot_candidates.py` | SS candidate histograms + rejection fractions |
| `compute_backgrounds.py` | Energy smearing, ROI cut, absolute background rates |
| `plot_source_flux.py` | Source flux plotting |
| `plot_event_processing_histos.py` | Event processing histograms |
| `plot_gamma_transport_histos.py` | Gamma transport diagnostics |

### Python utilities (py/util/)

| Script | Purpose |
|:-------|:--------|
| `show_detector_tree.py` | Print tracking detector geometry tree |
| `compute_geometry_fix.py` | Compute detector dimensions with uniform ICV wall thickness |
| `compute_icv_inner_surface.py` | Compare ICV inner surfaces vs tracking detector |
| `compute_icv_masses.py` | Compute mass-preserving ICV head dimensions |
| `validate_activities_vs_paper.py` | Validate source activities against published values |
| `validate_geometry_vs_paper.py` | Validate source geometry against published values |


## Key design decisions

1. **Factored flux tables**: (E, u) PDFs per source layer, not per event.
   Bi-214: one table. Tl-208: main + 3 independent companion tables with
   survival fractions. Companions sampled via independent Bernoulli trials.

2. **No cross-contributions**: barrel sources -> barrel flux only. Top/bottom
   likewise. No OCV_barrel -> ICV_top coupling.

3. **Compound propagation**: OCV -> MLI (vacuum, straight-line) -> ICV as
   sequential layers via `propagate_through_layers`. MLI gammas born in MLI,
   propagate through MLI (vacuum) then ICV (KN) — not generated inside ICV.

4. **Vacuum handling**: `propagate_in_source` detects vacuum (density <= 0)
   and advances straight to exit. Unifies dense and transparent sources in
   the same code path.

5. **Source taxonomy**: dense sources (KN self-shielding, compound
   propagation) vs transparent/ideal sources (vacuum, ~50% survival from
   hemisphere cut). PMT sources are transparent — merged into a single
   volume per location, with per-component rate bookkeeping.

6. **Virtual envelope unified on source geometry**: `make_virtual_envelope`
   uses `sg` for all sources. The `fk`-based method was removed after the
   geometry fix ensured ICV inner surfaces match the tracking detector
   boundary exactly. Epsilon inset (`ENVELOPE_INSET_CM = 0.01`) avoids
   floating-point boundary issues.

7. **Geometry fix**: ICV head R_outer and wall_thickness adjusted so all
   inner surfaces share R_inner = 82.1 cm while preserving original dome
   masses (107.7 kg top, 141.4 kg bottom). AirDome equator moved from
   z=145.6 to z=148.5 (ICV_top equator). New AirCyl fills gas gap
   [145.6, 148.5]. The real ICV is not a perfect cylinder+dome; this is
   a deliberate simplification for consistent geometry.

8. **PMT merged volumes**: PMT_TOP (3 components) and PMT_BOT (4 components)
   lumped into single transparent PDisk volumes. Orientation-based hemisphere
   filtering (:up -> downward gammas, :down -> upward gammas). Z span
   computed from source geometry, no hardwired numbers. PMT_BARREL (3
   components: cables, R8520, R8778 lower-ring) merged into a single
   transparent PCylShell with inward radial filtering.

12. **PMT bottom with passive LXe**: `pmt_bottom_lxe` adds compound
    propagation through the passive LXe slab between the PMT array and the
    cathode. Flux recorded at cathode boundary.

9. **Three-stage pipeline**:
   - Stage 1: `generate_flux_tables.jl` -> CSV flux tables
   - Stage 2: `run_source_backgrounds.jl` -> SS candidates
   - Stage 3: `compute_backgrounds.py` -> smearing + ROI + absolute rates

10. **E_max = line + 10 keV**, E_min = 2.2 MeV (ROI floor).

11. **Parallelization**: `julia -t N`, per-thread RNGs, weight-averaged merging.


## Supported sources

| Source | Type | Components | Status |
|:-------|:-----|:-----------|:-------|
| `cryostat_barrel` | Dense | OCV, MLI, ICV | Complete |
| `cryostat_top` | Dense | OCV, ICV | Complete |
| `cryostat_bottom` | Dense | OCV, ICV | Complete |
| `pmt_top` | Transparent | PMTs, bases, structure | Complete |
| `pmt_bottom` | Transparent | PMTs, bases, structure, R8778_dome | Complete |
| `pmt_bottom_lxe` | Compound | PMT bottom + passive LXe slab to cathode | Complete |
| `pmt_barrel` | Transparent | cables, R8520 skin PMTs, R8778 lower-ring PMTs | Complete |


## Geometry (source_geometry_lz_v1.json)

### ICV wall thickness fix
ICV_top: R_outer = 82.902, t = 0.802 (was R=83, t=0.8) -> R_inner = 82.1, mass = 107.7 kg (preserved).
ICV_bottom: R_outer = 83.292, t = 1.192 (was R=83, t=1.2) -> R_inner = 82.1, mass = 141.4 kg (preserved).
ICV_barrel: unchanged (R_inner = 82.1, t = 0.9).

### Tracking detector (detector_lz_v3.json)
- AirDome: cap at z_eq = 148.5 (was 145.6), R = 82.1, ar = 2.0
- AirCyl: new cylinder R = 82.1, z = [145.6, 148.5] (gas gap)
- LXe_passive, LZ_detector: unchanged


## Pipeline test results (50k events each)

| Source | Isotope | SS / 50k | Vetoed | No FV |
|:-------|:--------|:---------|:-------|:------|
| pmt_top | Bi214 | 2 | 72% | 28% |
| pmt_top | Tl208 | 0 | 83% | 17% |
| pmt_bottom | Bi214 | 1 | 0.3% | 99.7% |
| pmt_bottom | Tl208 | 0 | 0.5% | 99.5% |


## Performance

| Operation | Rate |
|:----------|:-----|
| Flux generation barrel (8 threads) | ~1.7M events/s |
| Flux generation top (8 threads) | ~750k events/s |
| Flux generation PMT (trivial) | ~9k events/s |
| Background processing | ~30k events/s |
