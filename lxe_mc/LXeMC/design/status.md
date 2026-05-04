# LXeMC — Current Status

## Architecture

| File | Role |
|:-----|:-----|
| `config.jl` | `SimConfig` — transport thresholds, stepping, physical constants |
| `nist_data.jl` | XCOM/ESTAR CSV loaders, log-log interpolation utilities |
| `physics_utils.jl` | Bethe-Heitler-Tsai bremsstrahlung dσ/dk, Coulomb correction |
| `materials.jl` | `Material` struct — density, composition, XCOM/ESTAR/brems data, physics methods |
| `geometry.jl` | `Cyl`/`Box` → `LCyl`/`LBox` → `PCyl`/`PBox` → `RCyl`, `Detector` |
| `sampling.jl` | MC samplers: Compton (KN), pair (BH), photoelectric (Sauter), brems (1/k rejection) |
| `tracking.jl` | `Track`, `Deposit`, `ParticleStack`, transport functions, event simulation, clustering |
| `LXeMC.jl` | Module wrapper, includes, exports |

## Physics implemented

### Photon interactions (from NIST XCOM tables)
- **Compton scattering**: Klein-Nishina sampled by composition+rejection (Butcher-Messel). Scattered photon continues; recoil electron pushed onto stack.
- **Pair production**: Bethe-Heitler energy split with flat+rejection sampling. Positron and electron pushed onto stack. Positrons annihilate at rest (two 511 keV back-to-back photons).
- **Photoelectric absorption**: above K-edge, E_K deposited locally (relaxation cascade is sub-mm); photoelectron tracked if T_e > Te_cut. Below K-edge, full energy deposited locally.

### Electron/positron transport
- Condensed-history with fixed step size (1 mm).
- Continuous collisional energy loss from NIST ESTAR tables.
- Discrete hard bremsstrahlung (k > k_min) via 1/k envelope rejection with pre-tabulated cross section and rejection ceiling.
- No multiple scattering (straight-line steps) — adequate for LXe where electron range ~ mm.
- No energy-loss straggling.

### What is NOT implemented
- Coherent (Rayleigh) scattering — negligible above 100 keV.
- Multiple scattering — needed for gas, not for LXe at current resolution.
- Delta-ray production.
- Energy-loss straggling (Landau/Vavilov fluctuations).
- Source generation (position, direction, spectrum sampling).
- Multi-volume transport (particle crossing volume boundaries).
- Passive shielding attenuation.

## Geometry and materials

### Geometry hierarchy
```
Cyl / Box          (geometric solid: dimensions)
  └─ LCyl / LBox   (logical volume: solid + position in MARS)
      └─ PCyl / PBox  (physical volume: logical + material)
          └─ RCyl      (radioactive volume: physical + specific activities)

Detector = MARS (vacuum) + list of physical volumes
```

### Defined detector: LZ
- MARS: cylinder R=110 cm, H_half=155 cm, Vacuum
- LXeTPC: cylinder R=72.8 cm, H_half=72.8 cm, LXe, centered at origin, mass ~7160 kg

### Materials loaded
| Material | Density [g/cm³] | Type | XCOM | ESTAR | mfp at 2.615 MeV |
|:---------|:----------------|:-----|:-----|:------|:-----------------|
| LXe | 2.953 | active+tracking | yes | yes | 8.9 cm |
| HPGXe | 0.089 | active+tracking | yes | yes | 297 cm |
| PTFE | 2.20 | passive | yes | no | 12.2 cm |
| Ti | 4.507 | passive | yes | no | 6.0 cm |
| Steel (Fe) | 7.99 | passive | yes | no | 3.3 cm |
| Vacuum | 0.0 | — | no | no | ∞ |

### Material capability flags
- `active`: XCOM + ESTAR + brems tables loaded. Can do photon + electron transport.
- `tracking`: condensed-history electron stepping (straight-line, no MS).
- `full_tracking` (future): + multiple scattering + straggling for gas.

## Source generation

**Current limitation**: primary photon is generated at a fixed position (default origin) with fixed direction (default +z) and fixed energy (caller-specified). No built-in source sampling.

**Needed**:
- Uniform generation within a volume (for bulk contamination).
- Surface generation (for external backgrounds).
- Isotropic direction sampling.
- Spectrum sampling (fixed gamma lines: 2.448 MeV Bi-214, 2.615 MeV Tl-208).
- Integration with `RCyl` radioactive volumes for flux-weighted generation.

## Simulation modes

| Mode | Function | Electrons | Use case |
|:-----|:---------|:----------|:---------|
| Full | `simulate_event` | Tracked (condensed-history) | Production: correct SS/MS topology |
| Photon-only | `simulate_event_photon_only` | Deposited locally | Fast studies, baseline comparison |

Key result at 2.615 MeV in LZ: SS fraction is 10% (full) vs 25% (photon-only), showing that mm-scale electron tracks drive MS topology in LXe.

## Output and analysis

- **Deposits**: list of `(position, energy, source_label)` per event.
- **Clustering**: `cluster_deposits_in_z(deposits, dz_cm; E_min)` groups deposits by z-proximity. Energy-weighted centroid. Sub-threshold clusters discarded.
- **SS/MS classification**: `is_single_site(deposits, dz, E_min)` — true if ≤1 cluster survives.
- **Driver script** (`scripts/run_test.jl`): CLI with `--nevents`, `--energy`, `--mode`, `--detector`, `--save-events`. Outputs summary to stdout + `summary.txt`; optionally writes long-format `events.csv` (one row per cluster, event fields repeated).

### Current parameters
| Parameter | Value | Justification |
|:----------|:------|:--------------|
| Te_cut | 400 keV | CSDA range 0.74 mm ≪ 3 mm z-resolution |
| k_min | 50 keV | Photon mfp 0.27 mm ≪ 3 mm |
| Egamma_cut | 10 keV | mfp 0.027 mm |
| ds_step | 1 mm | Fixed, adequate for LXe |
| dz_resolution | 3 mm | SS/MS z-clustering |
| E_cluster_min | 10 keV | Sub-threshold cluster filter |

## Performance

| Measurement | Events/sec |
|:------------|:-----------|
| Full mode (warm, 500k events) | 293,000 |
| Photon-only (warm) | ~400,000+ |
| Full mode (driver, 50k events) | ~58,000 |

At 293k events/sec, 10^8 events takes ~5.7 minutes.

Main CPU consumers (profiled): XCOM interpolation (22%), Compton electron direction (11%), ESTAR interpolation (7%), brems table lookup (6%).

## Test suite

14 test sets covering: config loading, material loading, detector geometry, XCOM consistency, NIST interpolation accuracy, Coulomb correction, bremsstrahlung dσ/dk, brems table vs direct integration, Compton edge, pair energy split, brems photon energy, energy conservation (full), SS fraction (Tl-208), energy conservation (photon-only).
