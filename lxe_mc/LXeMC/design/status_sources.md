# Sources Branch — Status (2026-05-13)

## Branch: `sources`

Pushed and synced with `origin/sources`.

## Recent milestone — PMT source geometry refactor (complete)

The PMT source taxonomy was overhauled in May 2026 to match the LZ
geometry published in the LZ 0vbb paper and arXiv:1910.09124v2. The
old `pmt_barrel` endpoint (which lumped cables, R8520 Top-Skin PMTs
and R8778 Bottom-Skin lower-ring PMTs into a single full-barrel
shell) has been deleted entirely. In its place there are four new
endpoints: `pmt_top_cables`, `pmt_bottom_cables`,
`pmt_skin_upper_ring`, `pmt_skin_lower_ring`. See the "Supported
sources" table below for the current state. Background validation
against LZ paper Table I: cables 0.48 vs LZ 0.526 cts/yr;
Skin PMTs+bases 0.27 vs LZ 0.274. The full migration plan is archived
at `design/legacy/fix_pmt_cable_geometry.md`.


## What's implemented

### Library code (src/)

| File | Content |
|:-----|:--------|
| `source_propagation.jl` | `propagate_in_source` (handles vacuum: straight-line exit for density <= 0), `random_position_in_volume` (PCyl, PCylShell, PDisk), `cos_theta_to_lxe` |
| `source_flux.jl` | `SourceFlux` (abstract), `SourceFluxBi214`, `SourceFluxTl208`, `SourceRateTable`, `SourceVolumeInfo`, `load_source_geometry`, `generate_flux_bi214/tl208`, `propagate_through_layers`, `generate_flux_compound_bi214/tl208` |
| `source_sampling.jl` | `sample_from_flux`, `sample_from_rate_table`, `sample_barrel_point`, `sample_disk_point(R, z, normal_z, rng)`, `sample_cap_point`, `reconstruct_direction`, `sample_gamma_from_flux` |
| `source_dispatch.jl` | `VirtualEnvelope` (kinds: `:barrel`, `:cap_up`, `:cap_down`, `:disk_flat`), `make_virtual_envelope(source, sg)` (unified, all sources from source geometry), `make_surface_sampler`, `dispatch_source_flux`, `merge_dispatch_results`, `supported_sources/isotopes` |
| `cryostat_sources.jl` | `cryostat_barrel_flux` (OCV compound + MLI compound + ICV self), `cryostat_top_flux`, `cryostat_bottom_flux`, `_build_rate_table`, `_get_activity_Bq` |
| `pmt_sources.jl` | `pmt_top_flux`, `pmt_bottom_flux`, `pmt_bottom_lxe_flux` (compound: PMT + passive LXe to cathode), `pmt_top_cables_flux` / `pmt_bottom_cables_flux` (TPC PMT cables in upper/lower conduits), `pmt_skin_upper_ring_flux` (93 R8520 'Top Skin' PMTs just below LXe surface), `pmt_skin_lower_ring_flux` (20 R8778 'Bottom Skin' side PMTs near cathode), `_merge_pmt_volume`, `_merge_pmt_barrel_volume` |
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
   computed from source geometry, no hardwired numbers. TPC PMT cables and
   Skin PMT rings have their own dedicated endpoints (see "Supported
   sources" below), each implemented as a single-component transparent
   volume (PCyl for cables, PCylShell for rings) with the appropriate
   filter.

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
| `pmt_top_cables` | Transparent | TPC PMT cables in upper conduit (above top PMT array, in AirDome) | Complete |
| `pmt_bottom_cables` | Transparent | TPC PMT cables in lower conduit (below bottom PMT array, in passive LXe) | Complete |
| `pmt_skin_upper_ring` | Transparent | 93 R8520 'Top Skin' 1-inch PMTs just below LXe surface, attached to outer side of field cage | Complete |
| `pmt_skin_lower_ring` | Transparent | 20 R8778 'Bottom Skin' side PMTs at the bottom of the side Skin, attached to ICV | Complete |


## Remaining work: field-cage sources

Six field-cage sources are defined in `source_geometry_lz_v1.json` but not
yet implemented in the flux pipeline:

| Source | Type | Transport | Notes |
|:-------|:-----|:----------|:------|
| `FC_PTFE` | Dense | KN self-shielding | PTFE reflector walls, barrel VE |
| `FC_rings` | Dense | KN self-shielding | Ti rings, barrel VE |
| `FC_topgrid` | Dense | KN | SS grid holder, slab geometry |
| `FC_botgrid` | Dense | KN | SS grid holder, slab geometry |
| `FC_resistors` | Transparent | Vacuum | Mass-equivalent carrier |
| `FC_sensors` | Transparent | Vacuum | Mass-equivalent carrier |

All emit directly into the LXe (no compound propagation). Each needs:
1. A flux function in a new `fc_sources.jl` (or added to an existing file)
2. A dispatch entry in `source_dispatch.jl`
3. A `make_virtual_envelope` case (barrel-type VE for all six)

The infrastructure (`dispatch_source_flux`, `make_virtual_envelope`,
flux table I/O, `generate_flux_tables.jl`, `run_source_backgrounds.jl`)
is already in place.

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
