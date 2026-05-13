# LXeMC — Claude Onboarding

LXeMC is a staged Monte Carlo for gamma and lepton transport in a large
liquid-xenon TPC (LZ detector geometry). It simulates 0vbb backgrounds
at ~10^6 events/s by using three accuracy tiers: full physics inside the
fiducial volume, photon-only KN transport in dense source regions, and a
fast-kernel with veto cuts everywhere else.

## Branch and layout

Branch: `sources` (active development of source flux pipeline).

```
LXeMC/
  src/           21 Julia source files (the package)
  scripts/       Pipeline scripts + diagnostics
    sandbox/     Standalone exploration code (not part of pipeline)
  test/          Test suite (runtests.jl)
  data/          detector_lz_v3.json, source_geometry_lz_v1.json,
                 sim_config.json, materials/, nist_data/
  design/        Design docs (see below)
  py/            Python analysis scripts (Stage 3)
  py/util/       Python validation utilities
```

## Architecture (3-tier model)

1. **FV (full physics)**: photon + lepton cascade with bremsstrahlung,
   pair production, photoelectric — `tracking.jl`
2. **Source regions (photon-only)**: KN transport through dense source
   material, kerma approximation — `source_propagation.jl`
3. **Everywhere else (fast kernel)**: photon-only with veto cuts and
   ROI-floor termination — `tracking_fast.jl`

See `design/lxemc.tex` for the complete reference.

## Code map

| File | Responsibility |
|:-----|:---------------|
| `config.jl` | `SimConfig` struct, JSON loader, defaults |
| `nist_data.jl` | XCOM/ESTAR CSV loaders, log-log interpolation |
| `physics_utils.jl` | Bethe-Heitler-Tsai dsigma/dk, Coulomb correction |
| `materials.jl` | `Material` struct, cross-section methods, brems pre-tabulation |
| `geometry_core.jl` | Primitive solids/wrappers, containment, ray entry/exit |
| `geometry.jl` | `TrackingDetector`, `FastKernelGeometry` compiler, FV volume |
| `sampling.jl` | MC samplers (Compton, pair, photoelectric, brems, distance) |
| `tracking.jl` | Full-stack: `Track`, `Deposit`, `transport_photon!`, `transport_lepton!`, `propagate_gamma` |
| `tracking_fast.jl` | Fast kernel: `transport_gamma_fastkernel`, `FastGammaTrackResult` |
| `event.jl` | Event orchestration: `process_event`, `EventProcessingResult` |
| `tracking_util.jl` | Diagnostics: single-step `propagate_gamma_fastkernel`, calib wrapper |
| `decays.jl` | Decay schemes, calibration sampling |
| `source_propagation.jl` | Photon-only KN in source, position sampling, cos_theta |
| `source_flux.jl` | Flux tables, source geometry loader, single/compound generators |
| `source_sampling.jl` | 2D PDF sampler, surface point generators, `sample_gamma_from_flux` |
| `cryostat_sources.jl` | `cryostat_barrel_flux`, `cryostat_top_flux`, `cryostat_bottom_flux` |
| `pmt_sources.jl` | `pmt_top_flux`, `pmt_bottom_flux`, `pmt_bottom_lxe_flux`, `pmt_top_cables_flux`, `pmt_bottom_cables_flux`, `pmt_skin_upper_ring_flux`, `pmt_skin_lower_ring_flux` |
| `flux_utils.jl` | Flux table merging, CSV I/O, JSON metadata |
| `source_dispatch.jl` | `dispatch_source_flux`, `VirtualEnvelope`, `make_surface_sampler` |
| `analysis.jl` | Offline: `cluster_deposits_in_z`, `is_single_site` |

## Pipeline scripts

| Script | Role |
|:-------|:-----|
| `generate_flux_tables.jl` | Stage 1: flux table generation |
| `run_source_backgrounds.jl` | Stage 2: background preselection |
| `analyze_backgrounds.jl` | Stage 3: SS/MS + smearing + ROI (Julia) |
| `py/compute_backgrounds.py` | Stage 3: smearing + ROI (Python) |
| `check_gamma_prop_simple.jl` | Diagnostic: single-step gamma propagation |
| `test_event_processing.jl` | Diagnostic: calibration event processing |

## Design docs — reading order

Read in this order for fastest orientation:

1. **This file** — orientation and code map
2. `design/caveats.md` — operational gotchas (read before editing code)
3. `design/status_sources.md` — what sources are implemented, what's remaining
4. `design/tracking_and_transport.md` — struct definitions, function signatures, transport flow
5. `design/geometry_v3.md` — JSON geometry schema, tracking tree, mass anchoring
6. `design/running_recipe.md` — how to run the 3-stage pipeline
7. `design/lxemc.tex` — full authoritative reference (physics, geometry, transport, sources, limitations)

## Before ending a session

Update these files if you changed code or design:

- `design/status_sources.md` — source table, remaining work
- `design/caveats.md` — operational gotchas, struct fields, dispatch rules
- `design/tracking_and_transport.md` — struct definitions, function signatures
- `design/running_recipe.md` — pipeline commands (if new sources added)
- `design/lxemc.tex` — authoritative reference (physics, geometry, transport)
- `design/geometry_v3.md` — if JSON schema or tracking tree changed
- This file (`CLAUDE.md`) — if files were added/removed/renamed or code map changed

## Workflow rules

- Always check `git branch` before any edits
- Ask for explicit approval before editing any file
- Do not make changes proactively without user confirmation
- Julia 1.10+, formatted with JuliaFormatter (indent=4, margin=92)
- Run `Pkg.test()` before committing
