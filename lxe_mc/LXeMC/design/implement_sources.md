# Source Flux Implementation Plan

## Overview

The source flux pipeline computes the (E, u) probability distributions of
gammas arriving at the inner surface of the ICV, where u = cos theta to the
inward surface normal. Three independent surfaces are considered: barrel,
top head, and bottom head. Each surface flux is the superposition of
contributions from different source layers (OCV, MLI, ICV), propagated
through all intervening material.

The tables store probability per generated decay per (E, u) bin. Absolute
rates are computed later by weighting with activity x mass x BR_from_chain.


## 1. Data structures

Location: `source.jl`

Already implemented:
- `SourceFlux` (abstract type)
- `SourceFluxBi214` -- single (E, u) PDF for one gamma per decay
- `SourceFluxTl208` -- main (E, u) PDF + 3 companion PDFs with BRs and
  survival fractions

To add:
- `SourceRateTable` -- activity-weighted sum of component flux tables for
  one surface, in units of gammas/sec per (E, u) bin. Contains:
  - `pdf_rate::Matrix{Float64}` (n_E x n_u, gammas/sec/bin)
  - binning parameters (E_min, E_max, n_E, n_u)
  - list of component names and their individual total rates
  - surface identifier (`:barrel`, `:top`, `:bottom`)

The summed rate table is computed as:

    pdf_rate[i,j] = sum_k ( A_k * m_k * pdf_k[i,j] )

where A_k is specific activity [Bq/kg], m_k is mass [kg], and pdf_k is
the probability-per-decay table for component k.


## 2. Source geometry loader

Location: `source.jl`

Function: `load_source_geometry(path, materials) -> Dict{String, ...}`

Reads `source_geometry_lz_v1.json` and builds a `PhysicalVolume` object
(`PCylShell` or `PDisk`) for each source entry. Each entry carries:
- the `PhysicalVolume` (for `random_position_in_volume`, `propagate_in_source`,
  `cos_theta_to_lxe`)
- material
- activity dict (Bi214_mBq_per_kg, Tl208_mBq_per_kg)
- mass (computed from geometry, or from `equivalent_mass_kg` for virtual sources)
- transport mode: `:KN` or `:transparent`
- source class: `shell_source` or `virtual_source`


## 3. Single-source flux generation

Location: `source.jl`

### `generate_flux_bi214(N, source_vol, cfg, rng; E_min, E_max, n_E, n_u)`

Generate N Bi-214 decays (single 2.448 MeV gamma each) inside `source_vol`.
For each decay:
1. Sample position via `random_position_in_volume`
2. Sample isotropic direction
3. Propagate through source material via `propagate_in_source` (KN)
4. If exited: check direction (u = cos_theta_to_lxe > 0)
5. If forward and in energy window: tally (E, u) bin

Returns `SourceFluxBi214` with pdf normalized to N_generated.

### `generate_flux_tl208(N, source_vol, cfg, rng; ...)`

Generate N Tl-208 decays. Each decay produces:
- 1 main gamma at 2.615 MeV (always)
- Up to 3 companions via independent Bernoulli trials:
  583 keV (85%), 511 keV (23%), 861 keV (12%)

Each gamma is propagated independently through the source. The main gamma
fills `pdf_main`. Each companion line fills its own `pdf_companion[i]`.
The survival fraction `companion_f[i]` is computed as:
N_surviving_companion_i / N_generated_companion_i (where N_generated is
the number of times that companion fired via its BR).

Returns `SourceFluxTl208`.


## 4. Compound propagation helper

Location: `source.jl`

### `propagate_through_layers(E, pos, dir, layers, cfg, rng)`

Propagate a gamma through a sequence of material layers. Each layer is a
`PhysicalVolume`. For each layer:
- If material is vacuum (density <= 0): straight-line to exit
  (`distance_to_exit` on the layer, advance position, no interaction)
- If material is real (Ti, etc.): call `propagate_in_source` (KN transport)
- If absorbed: return (:absorbed, ...)
- If exited: continue to next layer

Returns (:exited, E_final, pos_final, dir_final) or (:absorbed, ...).

This handles the OCV -> vacuum gap -> ICV compound propagation as a single
call with layers = [OCV_vol, vacuum_vol, ICV_vol].

### `generate_flux_compound_bi214(N, source_vol, layers, cfg, rng; ...)`

Generate N decays in `source_vol`, propagate each gamma through the layer
sequence, tally survivors at the final exit. Returns `SourceFluxBi214`.

### `generate_flux_compound_tl208(N, source_vol, layers, cfg, rng; ...)`

Same for Tl-208 (main + companions). Returns `SourceFluxTl208`.


## 5. Cryostat flux functions

Location: `source.jl` (or a new `cryostat_sources.jl` if source.jl grows too large)

### `cryostat_barrel_flux(N, source_geom, materials, cfg, rng)`

Compute barrel flux tables. Three contributing sources:

1. **OCV contribution**: generate decays in OCV_barrel, propagate through
   [OCV_barrel, vacuum_gap, ICV_barrel]. Tally at ICV inner surface.
2. **MLI contribution**: generate decays in MLI (transparent source, so
   gammas start at full decay energy on ICV outer surface), propagate
   through [ICV_barrel]. Tally at ICV inner surface.
3. **ICV contribution**: generate decays in ICV_barrel, propagate through
   [ICV_barrel]. Tally at ICV inner surface.

Returns per isotope (Bi214, Tl208):
- 3 component flux tables (OCV, MLI, ICV)
- 1 summed `SourceRateTable`

Total: 8 tables (4 x 2 isotopes).

### `cryostat_top_flux(N, source_geom, materials, cfg, rng)`

Two contributing sources:

1. **OCV contribution**: OCV_top -> vacuum -> ICV_top
2. **ICV contribution**: ICV_top self

Returns per isotope:
- 2 component flux tables
- 1 summed `SourceRateTable`

Total: 6 tables (3 x 2 isotopes).

### `cryostat_bottom_flux(N, source_geom, materials, cfg, rng)`

Two contributing sources:

1. **OCV contribution**: OCV_bottom -> vacuum -> ICV_bottom
2. **ICV contribution**: ICV_bottom self

Returns per isotope:
- 2 component flux tables
- 1 summed `SourceRateTable`

Total: 6 tables (3 x 2 isotopes).

Grand total: 20 flux tables.


## 6. Sampling from flux tables

Location: `source.jl`

### `sample_from_flux(flux::SourceFluxBi214, rng) -> (E, u)`

Sample one (E, u) pair from the 2D PDF using inverse CDF or rejection
sampling. Returns the energy and cos theta to the surface normal.

### `sample_from_flux(flux::SourceFluxTl208, rng) -> Vector{Tuple{Float64,Float64}}`

Sample main gamma (E, u) from `pdf_main` (always). For each companion i,
fire Bernoulli with probability `companion_BR[i] * companion_f[i]`; if it
fires, sample (E, u) from `pdf_companion[i]`. Returns 1--4 (E, u) tuples.

### `sample_from_rate_table(rate::SourceRateTable, rng) -> (E, u)`

Sample from the activity-weighted summed table. Used for the combined
cryostat flux.


## 7. Surface position and direction reconstruction

Location: `source.jl`

These functions convert (E, u) sampled from a flux table into a full
`SampledGamma` with 3D position and direction on the ICV inner surface.
The surface geometry comes from the source geometry JSON, not from the
flux table.

### `sample_barrel_point(R_inner, z_min, z_max, rng) -> (position, normal)`

Sample a uniform point on the ICV barrel inner surface:
- phi = 2pi * rand
- z = z_min + (z_max - z_min) * rand
- position = (R_inner * cos(phi), R_inner * sin(phi), z)
- normal = -r_hat (radially inward)

### `sample_cap_point(R, aspect_ratio, z_equator, orientation, rng) -> (position, normal)`

Sample a point on the ICV head inner surface (ellipsoidal cap):
- Sample r with correct area weighting on the ellipsoid
- phi = 2pi * rand
- z from the ellipsoid equation: z = z_equator + sgn * c * sqrt(1 - r^2/a^2)
- Compute local normal from the ellipsoid gradient

### `reconstruct_direction(u, normal, rng) -> Vector{Float64}`

Given u = cos theta to the local surface normal and a random azimuthal
angle, construct the 3D unit direction vector. Uses `rotate_to_global`
with the normal as the reference direction.

### `sample_gamma_from_flux(flux, surface_sampler, rng) -> SampledGamma`

Combines:
1. `sample_from_flux` -> (E, u)
2. `surface_sampler` -> (position, normal)
3. `reconstruct_direction(u, normal, rng)` -> direction

Returns a `SampledGamma(E, position, direction)` ready for `process_event`.

For Tl-208, returns `Vector{SampledGamma}` (1--4 gammas per event). All
gammas from one event share the same surface point (they entered from the
same decay, same cryostat location).


## 8. Driver script

Location: `scripts/run_cryostat_fluxes.jl`

Command-line script that:
1. Loads config, materials, source geometry
2. Calls `cryostat_barrel_flux`, `cryostat_top_flux`, `cryostat_bottom_flux`
3. Saves all 20 flux tables to disk (JSON)
4. Prints summary: survival fractions, peak bin fractions, total rates
   per source component

Usage:
```
julia --project=. scripts/run_cryostat_fluxes.jl --n 100000 --seed 42 --outdir results/fluxes/
```


## Implementation order (dense sources, points 1--8)

1. `SourceRateTable` struct (point 1)
2. `load_source_geometry` (point 2)
3. `generate_flux_bi214` / `generate_flux_tl208` for single source (point 3)
4. `propagate_through_layers` and compound variants (point 4)
5. `cryostat_barrel_flux` / `cryostat_top_flux` / `cryostat_bottom_flux` (point 5)
6. Sampling functions: `sample_from_flux`, `sample_from_rate_table` (point 6)
7. Surface reconstruction: `sample_barrel_point`, `sample_cap_point`,
   `reconstruct_direction`, `sample_gamma_from_flux` (point 7)
8. Driver script (point 8)

**All 8 points complete** for cryostat (dense) sources.

---

## 9. Source taxonomy: dense vs ideal

Sources fall into two categories based on their `transport_source` and
`source_class` fields in `source_geometry_lz_v1.json`:

### Dense sources

- `transport_source = "KN"`, `source_class = "shell_source"`
- Real material (Ti, PTFE, SS) with nonzero density
- Gammas undergo KN self-shielding inside the source volume
- May require compound propagation through intervening layers (OCV -> vacuum -> ICV)
- Flux tables are non-trivial: energy degradation + angular redistribution
- Virtual envelope sits on the **detector boundary** (ICV inner surface)
- Examples: OCV_barrel, ICV_top, FC_PTFE, FC_rings, FC_topgrid, FC_botgrid

### Ideal sources

- `transport_source = "transparent"`, `source_class = "virtual_source"`
- Material is Vacuum (no self-shielding, no energy loss)
- Mass given by `equivalent_mass_kg` (not geometric)
- Gammas emerge at full decay energy with isotropic angular distribution
- Flux tables are **trivial**: delta function in E, flat in u, survival = 1.0
- Virtual envelope is the **source geometry itself** (a thin disk or shell
  inside the detector volume)
- Examples: MLI, FC_resistors, FC_sensors, PMT_TOP_*, PMT_BOT_*, PMT_BARREL_*

The key difference for the pipeline: dense sources require expensive MC to
build flux tables; ideal sources build flux tables analytically in microseconds.
The output format (SourceFluxBi214 / SourceFluxTl208 structs, CSV files) is
identical for both, so the downstream pipeline (fast kernel + FV + smearing)
is unchanged.


## 10. Ideal source flux generation

### Flag: `--ideal_source`

The `generate_flux_tables.jl` script accepts `--ideal_source` (default `false`).
When `true`, the script skips MC propagation and constructs flux tables
analytically:

### Bi-214 (ideal)

- Single E bin centered at 2.448 MeV (n_E = 1 or minimal binning)
- Uniform in u (isotropic emission, all u bins equal)
- `n_survived = N`, `survival_fraction = 1.0`

### Tl-208 (ideal)

- Main gamma: single E bin at 2.615 MeV, flat in u, survival = 1.0
- Companion 1 (583 keV): single E bin, flat in u, survival = 1.0, BR = 0.85
- Companion 2 (511 keV): single E bin, flat in u, survival = 1.0, BR = 0.23
- Companion 3 (861 keV): single E bin, flat in u, survival = 1.0, BR = 0.12

All companion survival fractions are 1.0 (no material to attenuate).

### Library function

Location: `source_flux.jl`

```julia
generate_ideal_flux_bi214(N; E_gamma=2.448, n_u=10) -> SourceFluxBi214
generate_ideal_flux_tl208(N; E_main=2.615, n_u=10) -> SourceFluxTl208
```

These build the structs directly — no RNG, no propagation, no material.
The `N` parameter sets `n_generated` for rate normalization (activity × mass
× pdf_per_decay gives the correct absolute rate).


## 11. PMT_TOP sources

Three independent sources sharing a single virtual envelope:

| Source           | Mass (kg) | Bi214 (mBq/kg) | Tl208 (mBq/kg) |
|:-----------------|:----------|:----------------|:----------------|
| PMT_TOP_PMTs     | 47.07     | 3.22            | 1.61            |
| PMT_TOP_bases    | 1.434     | 75.9            | 33.1            |
| PMT_TOP_structure| 83.0      | 1.60            | 1.06            |

### Shared virtual envelope

- Kind: `:disk` (flat disk, not a dome cap)
- R = 72.8 cm
- z = 152.6 cm (bottom of the PMT stack, closest to LXe)
- Located inside the AirDome (gas phase), above the gate at z = 145.6

### Pipeline

1. **Flux tables** (trivial): `--ideal_source` flag, one run per source
   component (PMTs, bases, structure) × isotope (Bi214, Tl208). Six runs
   total, each instantaneous.
2. **Background processing**: sample from flux tables, place gammas on the
   shared VE disk, process through fast kernel + FV stack. Gammas travel
   downward through gas dome into LXe.
3. **Background rates**: smearing + ROI, per-component and summed.

### Source dispatch

`dispatch_source_flux` extended with:
- `"pmt_top_pmts"`, `"pmt_top_bases"`, `"pmt_top_structure"`

Each routes to `generate_ideal_flux_bi214` / `generate_ideal_flux_tl208`.


## 12. Virtual envelope generalization

The `VirtualEnvelope` struct and `make_virtual_envelope` are extended to
handle both dense and ideal sources:

### Dense sources (existing)

VE sits on the detector boundary (ICV inner surface):
- Barrel: cylinder at R = 82.09 cm (LZ_detector R - inset)
- Top: AirDome cap surface (R = 82.09, z_eq = 145.6, ar = 2.0)
- Bottom: LXe_passive bottom cap (R = 82.09, z_eq = -14.1, ar = 3.0)

### Ideal sources (new)

VE is the source geometry itself:
- PMT_TOP: flat disk at R = 72.8, z = 152.6
- PMT_BOT: flat disk at R = 72.8, z ≈ -15.75
- PMT_BARREL: cylindrical shell at R ≈ 81.2--81.4, z = [-13.75, 145.6]
- FC_resistors, FC_sensors: cylindrical shell at R ≈ 72.6--72.7

A new VE kind `:disk_flat` is added for flat disks (as opposed to `:cap_up`
/ `:cap_down` which are ellipsoidal). The surface sampler for `:disk_flat`
is simply uniform in (r^2, phi) at fixed z.

### `make_virtual_envelope` dispatch

```
make_virtual_envelope(source, fk)          # dense: from FastKernelGeometry
make_virtual_envelope(source, sg)          # ideal: from source geometry dict
```

The script or dispatch layer chooses which method to call based on the
source's `transport_source` field.


## Implementation order (ideal sources, points 9--12)

9.  Add `:disk_flat` VE kind + `sample_disk_point` surface sampler
10. `generate_ideal_flux_bi214` / `generate_ideal_flux_tl208` in `source_flux.jl`
11. `--ideal_source` flag in `generate_flux_tables.jl`
12. PMT_TOP dispatch entries + `make_virtual_envelope` from source geometry
13. Run PMT_TOP × {Bi214, Tl208} through full pipeline
14. Extend to PMT_BOT and PMT_BARREL (same pattern)
