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
    legacy/      Archived plans + superseded docs (read-only history)
  docs/          LZ reference papers (LZ-TDR.txt, LZ-1910.09124v2.txt,
                 lz_bb0nu.txt + lz_bb0nu_summary.tex) and analysis docs
    analysis/    MC vs LZ paper comparison (background_status.md)
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
3. `design/status_sources.md` — what sources are implemented, recent milestone, what's remaining
4. `docs/analysis/background_status.md` — **current MC-vs-LZ-paper comparison table** (where each per-component rate stands; quick way to see today's state)
5. `design/tracking_and_transport.md` — struct definitions, function signatures, transport flow
6. `design/geometry_v3.md` — JSON geometry schema, tracking tree, mass anchoring
7. `design/running_recipe.md` — how to run the 3-stage pipeline
8. `design/lxemc.tex` — full authoritative reference (physics, geometry, transport, sources, limitations)

Archived design (do not edit, read for history only):

- `design/legacy/fix_transport.md` — applied: tag-based dispatch in transport helpers (commit 57895b0)
- `design/legacy/fix_pmt_cable_geometry.md` — applied: PMT cables + Skin PMT relocation, pmt_barrel deleted (commits 2dc8c7a..e3989c4)

## Active work / open issues

- **TPC PMTs / bases / structures all under-predict LZ paper by 12–36%**.
  Cables and Skin PMTs+bases now match LZ (~10%). The remaining
  under-prediction is shared across three components and likely has a
  single cause (PMT-array structural attenuation? analysis cuts?). Not
  yet investigated. See `docs/analysis/background_status.md` for the
  current comparison table.
- **Field-cage sources** (FC_PTFE, FC_rings, FC_topgrid, FC_botgrid,
  FC_resistors, FC_sensors): defined in JSON but no flux pipeline yet.
  See `design/status_sources.md` "Remaining work" section.

## Before ending a session

Update these files if you changed code or design:

- `design/status_sources.md` — source table, remaining work, recent milestones
- `docs/analysis/background_status.md` — refresh MC-vs-LZ comparison if rates changed
- `design/caveats.md` — operational gotchas, struct fields, dispatch rules
- `design/tracking_and_transport.md` — struct definitions, function signatures
- `design/running_recipe.md` — pipeline commands (if new sources added)
- `design/lxemc.tex` — authoritative reference (physics, geometry, transport)
- `design/geometry_v3.md` — if JSON schema or tracking tree changed
- This file (`CLAUDE.md`) — if files were added/removed/renamed, code map changed, or "Active work" list needs refresh

When a multi-commit refactor lands, **archive the plan doc**: move
`design/<plan>.md` → `design/legacy/<plan>.md` and add a one-line
pointer to it in this file's "Archived design" list. Keeps `design/`
to currently-active plans only.

## Workflow rules

- Always check `git branch` before any edits
- Ask for explicit approval before editing any file
- Do not make changes proactively without user confirmation
- Julia 1.10+, formatted with JuliaFormatter (indent=4, margin=92)
- Run `Pkg.test()` before committing
