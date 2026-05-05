# LXeMC — Current Status

## Architecture

| File | Role |
|:-----|:-----|
| `config.jl` | `SimConfig` — transport thresholds, stepping, veto thresholds, physical constants |
| `nist_data.jl` | XCOM/ESTAR CSV loaders, log-log interpolation (standard + pre-logged fast variant) |
| `physics_utils.jl` | Bethe-Heitler-Tsai bremsstrahlung dσ/dk, Coulomb correction |
| `materials.jl` | `Material` struct — density, composition, XCOM/ESTAR/brems data, pre-logged grids, physics methods |
| `geometry.jl` | `Cyl`/`CylShell`/`Box`/`Disk` → `LCyl`/`LCylShell`/`LBox`/`LDisk` → `PCyl`/`PCylShell`/`PBox`/`PDisk` → `RCyl`/`RCylShell`/`RDisk`, `Detector`, ray intersection |
| `sampling.jl` | MC samplers: Compton (KN), pair (BH), photoelectric (Sauter), brems (1/k rejection) |
| `tracking.jl` | Transport functions, event simulation, veto pre-filter, multi-volume propagation, clustering |
| `decays.jl` | Bi-214 and Tl-208 decay schemes, event generators |
| `source.jl` | Source flux generation (Stage 1), FluxTable, random position sampling |
| `LXeMC.jl` | Module wrapper, includes, exports |

## Physics implemented

### Photon interactions (from NIST XCOM tables)
- **Compton scattering**: Klein-Nishina sampled by composition+rejection (Butcher-Messel)
- **Pair production**: Bethe-Heitler energy split with flat+rejection
- **Photoelectric absorption**: above K-edge, E_K deposited locally; photoelectron tracked if T_e > Te_cut

### Electron/positron transport
- Condensed-history with fixed 1 mm step
- Continuous collisional energy loss from NIST ESTAR
- Discrete hard bremsstrahlung (k > k_min) via 1/k envelope rejection with pre-tabulated σ and M
- No multiple scattering (straight-line steps) — adequate for LXe
- Positron annihilation at rest → two 511 keV photons

### Current transport parameters
| Parameter | Value | Justification |
|:----------|:------|:--------------|
| Te_cut | 400 keV | CSDA range 0.74 mm ≪ 3 mm z-resolution |
| k_min | 50 keV | Photon mfp 0.27 mm ≪ 3 mm |
| Egamma_cut | 10 keV | mfp 0.027 mm |
| ds_step | 1 mm | Fixed, adequate for LXe |
| dz_resolution | 3 mm | SS/MS z-clustering |
| E_cluster_min | 10 keV | Sub-threshold cluster filter |
| veto_TPC | 10 keV | Detection threshold in TPC active region |
| veto_skin | 100 keV | Detection threshold in LXe skin |

## Simulation modes

| Mode | Function | Electrons | Use case |
|:-----|:---------|:----------|:---------|
| Full | `simulate_event` | Tracked (condensed-history) | Production: correct SS/MS topology |
| Photon-only | `simulate_event_photon_only` | Deposited locally | Fast studies, baseline comparison |

Key result at 2.615 MeV in LZ: SS fraction is 10% (full) vs 25% (photon-only).

## Geometry and materials

### Geometry hierarchy
```
Cyl / CylShell / Box / Disk   (geometric solid: dimensions)
  └─ LCyl / LCylShell / LBox / LDisk   (logical: solid + position in MARS)
      └─ PCyl / PCylShell / PBox / PDisk  (physical: logical + material)
          └─ RCyl / RCylShell / RDisk      (radioactive: physical + activities)

Detector = MARS (vacuum) + flat list of physical volumes + optional fiducial
```

### Volume tiling (no overlaps)
The LZ geometry is tiled so no two volumes in the list overlap:
- LXe regions (TPC, skin, RFR, dome) tile the interior by r and z ranges
- Ti volumes (ICV/OCV barrels, heads, flanges) surround the LXe
- FV is stored separately (not in volumes list), checked via `classify_lxe_region`
- `find_volume(det, pos)` linear scan works correctly because of this tiling

### Defined detector: LZ (detector_lz.json)
16 volumes + FV:
- **LXeTPC**: R=72.8, z=[0, 145.6], LXe, 7.16 t (active, cathode at z=0)
- **FieldCage**: PTFE shell R=[72.8, 74.3], z=[-13.75, 145.6]
- **Skin**: LXe shell R=[74.3, 82.1], z=[-13.75, 145.6], 1.80 t
- **RFR**: LXe R=72.8, z=[-13.75, 0], 0.68 t
- **Dome**: LXe R=82.1, z=[-69.0, -13.75], ~3.5 t (cylindrical approx)
- **ICV**: barrel + 2:1 top head + 3:1 bottom head, Ti, total 649 kg
- **OCV**: barrel + 2:1 top/bottom heads, Ti, total 777 kg
- **5 flanges**: Ti, total ~605 kg
- **FV**: R=39.0, z=[26, 96], 0.988 t (~1 tonne)

### Materials loaded (from materials.json)
| Material | Density | Type | XCOM | ESTAR | mfp at 2.615 MeV |
|:---------|:--------|:-----|:-----|:------|:-----------------|
| LXe | 2.953 g/cm³ | active+tracking | yes | yes | 8.9 cm |
| HPGXe | 0.089 g/cm³ | active+tracking | yes | yes | 297 cm |
| PTFE | 2.20 g/cm³ | passive | yes | no | 12.2 cm |
| Ti | 4.507 g/cm³ | passive | yes | no | 6.0 cm |
| Steel (Fe) | 7.99 g/cm³ | passive | yes | no | 3.3 cm |
| Vacuum | 0.0 | — | no | no | ∞ |

### Activities (from LZ assays, in detector_lz.json)
All Ti components: Bi-214 = 0.08 mBq/kg, Tl-208 = 0.22 mBq/kg.

## Ray intersection and multi-volume propagation

### Ray-volume intersection (geometry.jl)
- `distance_to_entry(pos, dir, logical)`: ray from outside → entry point. Analytic for LCyl, LCylShell. Step-search fallback for LBox, LDisk.
- `distance_to_exit(pos, dir, logical)`: ray from inside → exit point.
- `next_volume(pos, dir, det)`: scan all volumes, return nearest hit along ray.

### Multi-volume propagation (tracking.jl)
`propagate_to_lxe(E, pos, dir, det, cfg, rng)`:
- Steps through MARS, switching materials at volume boundaries
- **Vacuum**: straight-line to next volume via `next_volume`
- **Ti/PTFE**: KN scatter or absorb with per-material XCOM cross sections
- **LXe regions**: classify via `classify_lxe_region`, apply region-specific veto:
  - TPC active (outside FV): veto_TPC = 10 keV
  - Skin: veto_skin = 100 keV
  - Passive (dome, RFR): no veto (Inf threshold)
  - FV: accept for full simulation
- Returns `PropagationResult` (:accepted, :vetoed, :lost)

## Background simulation chain

### Stage 1 — Source flux (source.jl, run_source_flux.jl)
For each source volume (e.g., ICV barrel):
1. Generate N decays (Bi-214: single gamma; Tl-208: main + companion gammas)
2. All gammas from one decay originate at same random position in source
3. Propagate each gamma through source material (KN photon-only)
4. Direction filter: only gammas toward LXe (cos θ > 0)
5. Companion veto: if 2+ visible gammas (E > 100 keV) with |Δz| > 3 mm → discard
6. Energy window: [2.37, 2.62] MeV
7. Output: `FluxTable` — normalized pdf (probability per decay per (E, u) bin)

Results from ICV barrel (100k events):
- Bi-214: 50.7% backward, 37.9% surviving, 99.7% in peak bin
- Tl-208: 26.2% backward, 25.4% companion-vetoed, 20.2% surviving

### Stage 2a — Multi-volume LXe self-shielding veto (propagate_to_lxe)
Gamma from LXe surface → propagate through LXe (and possibly Ti/vacuum if from OCV).
LXe self-shielding distances to FV:
| Path | Distance | λ | Survival |
|:-----|:---------|:--|:---------|
| Bottom (cathode → FV) | 26.0 cm | 2.9 | ~5% |
| Radial (wall → FV) | 33.8 cm | 3.8 | ~2% |
| Top (gate → FV) | 49.6 cm | 5.5 | ~0.4% |

Veto pre-filter throughput: ~2.2M gammas/sec.

### Stage 2b — Full simulation in FV (NOT YET IMPLEMENTED)
Needed:
1. `sample_from_flux(ft, rng)` — sample (E, u) from flux table pdf
2. Convert to gamma at TPC surface with appropriate position/direction
3. `propagate_to_lxe` → accept/veto
4. For accepted: `simulate_event` in FV → deposits
5. Energy resolution smearing (Gaussian, σ/E configurable)
6. ROI cut: smeared energy in [Q_ββ - ΔE, Q_ββ + ΔE]
7. SS cut: `is_single_site`
8. Count: fraction surviving all cuts × decay rate = background/yr

### Absolute rates (source_rates.jl)
`rate [gammas/yr] = A_specific × mass × BR_from_chain × seconds_per_year × survival_fraction`

ICV barrel example:
- Bi-214: 15,125 decays/yr × 0.379 survival = 5,734 gammas/yr to LXe
- Tl-208: 8,319 decays/yr × 0.202 survival = 1,682 gammas/yr to LXe

## Performance

| Measurement | Events/sec |
|:------------|:-----------|
| Full simulation (warm, 500k events) | 293,000 |
| Photon-only mode | ~400,000+ |
| Veto pre-filter (propagate_to_fiducial) | 2,200,000 |
| Source flux generation | 210,000 |

10^8 source gammas: ~45s pre-filter + ~10s full sim = under 1 minute.

## Scripts

| Script | Purpose |
|:-------|:--------|
| `scripts/run_test.jl` | Run N events in detector, full/photon-only mode |
| `scripts/bench_prefilter.jl` | Benchmark veto pre-filter from TPC surfaces |
| `scripts/run_source_flux.jl` | Stage 1: generate flux table for a source volume |
| `scripts/source_rates.jl` | Compute annual rates from pre-computed flux tables |
| `py/plot_source_flux.py` | Plot flux results (4-panel: dP/dE, dP/du, heatmap, fractions) |

## Test suite

18 test sets covering: config, materials, detector geometry (masses validated against LZ reference), geometry solids (Cyl, CylShell, Disk), ray intersection, XCOM consistency, NIST interpolation, Coulomb correction, bremsstrahlung, Compton edge, pair split, brems sampling, energy conservation (full + photon-only), SS fraction, veto pre-filter (attenuation physics + interaction count), source flux (Bi-214 + Tl-208).

## What's next

1. **Tests for propagate_to_lxe**: region classification, multi-volume traversal (OCV → ICV → LXe), passive LXe (dome, no veto)
2. **Stage 2b implementation**: sample_from_flux → propagate_to_lxe → simulate_event → resolution → ROI + SS → final background count
3. **Run all source volumes**: generate flux tables for all 11 Ti components, compute total background budget
4. **Future**: insulation source (very radioactive), field cage (PTFE), source generators for surface/bulk contamination
