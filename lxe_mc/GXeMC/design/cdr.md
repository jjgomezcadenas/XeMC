# GXeMC — Conceptual Design Report

## Scope

GXeMC is a Monte Carlo simulation of photon and lepton transport in
high-pressure gaseous xenon (GXe) time projection chambers, targeting
NEXT-like detector configurations. It produces voxelized 3D energy
deposits suitable for track topology analysis.

Intended use cases:

- 0vbb signal simulation: two-electron track topology at Qbb = 2.458 MeV
- Gamma background characterisation: SS/MS discrimination from track
  structure
- Track imaging studies in pure Xe, pure He, and Xe/He mixtures


## Transport model

GXeMC uses a single-tier full-physics transport — every photon and
lepton is tracked with the same physics everywhere in the detector.

### Photon transport

Three interaction channels:

- **Compton scattering**: incoherent channel from XCOM. Klein-Nishina
  kinematics; recoil electron pushed to particle stack.
- **Pair production**: Bethe-Heitler (nuclear + electron). Electron and
  positron pushed to stack. Positron annihilation at end of range
  produces two back-to-back 511 keV gammas.
- **Photoelectric absorption**: total cross-section from XCOM. Shell
  treatment:
  - K-shell (E_K = 34.56 keV for Xe): the binding energy is emitted as
    a fluorescence X-ray (~30 keV) pushed to the stack as a secondary
    gamma. The photoelectron carries T_e = E_gamma - E_K.
  - L/M shells (< 5.5 keV): Auger electrons have sub-mm range in
    10 atm GXe. Deposited locally.

Rayleigh (coherent) scattering is not included in the initial version.
At MeV energies in Xe it is a ~1% effect on the total cross-section.
Can be added later if validation warrants it.

Photon tracking cutoff: Egamma_cut ~ 10 keV. Below this energy the
photon deposits its remaining energy locally.

### Lepton transport

Condensed-history stepping with three processes per step:

1. **Continuous collisional energy loss**: dE/dx from NIST ESTAR
   (collisional stopping power). Energy lost over step ds deposited at
   step midpoint.

2. **Multiple Coulomb scattering**: Highland formula applied per step.
   Scattering angle theta_0 = (13.6 MeV / beta*c*p) * sqrt(ds/X0) *
   [1 + 0.038 * ln(ds/X0)]. Random angular deflection applied to
   direction vector.

3. **Discrete bremsstrahlung**: per-step probability from total brems
   cross-section. On a hit, photon energy sampled from Bethe-Heitler-Tsai
   differential cross-section, pushed to stack.

4. **Delta-ray production**: hard knock-on electrons above a production
   threshold sampled from Moller cross-section and pushed to stack.

Step size: ds_step = 100 um (tunable). Shortened if nominal energy loss
would exceed the particle's kinetic energy.

Lepton tracking cutoff: Te_cut ~ 10 keV. Below this energy the residual
kinetic energy is deposited locally. Positrons produce two 511 keV gammas
at the termination point.

### Early termination

After each interaction, the maximum depositable energy is computed as:

    E_max = E_deposited + sum(E_i for all particles on the stack)

If E_max < Qbb - E_min (configurable, e.g. Qbb - 50 keV), the event
cannot contribute to the ROI and is killed immediately. This avoids
tracking electrons from background gammas that have already lost too much
energy to escaping photons.


## Cross-section data

All energy-dependent quantities are driven by tabulated NIST data:

- **XCOM**: photon mass attenuation coefficients (cm2/g) on a
  per-material energy grid, with seven channels (coherent, incoherent,
  photoelectric, pair-nuclear, pair-electron, total-with-coherent,
  total-without-coherent). GXeMC uses incoherent + photoelectric + pair
  for transport. Per-atom cross-sections obtained via
  sigma = (mu/rho) * A_eff / N_A.

- **ESTAR**: electron mass collision and radiative stopping powers
  (MeV cm2/g) and CSDA range (g/cm2). The collisional stopping power
  drives the continuous energy loss. The radiative stopping power
  provides the total bremsstrahlung cross-section for discrete sampling.

Both datasets are stored as CSV files and loaded at startup with log-log
interpolation for fast lookup at arbitrary energies.


## Gas model

Only Pure Xe for 0vbb studies


## Geometry

The detector is defined by a JSON file containing a tree of geometric
volumes. Each volume carries a shape, position, material, and a tag
that drives simulation policy (sensitive, passive, structural).

Supported analytic primitives:

- **Cyl**: filled cylinder (R, half_height)
- **CylShell**: cylindrical shell (R_inner, wall_thickness, half_height)
- **Box**: rectangular block (dx, dy, dz)
- **Cap**: ellipsoidal cap (R, aspect_ratio)

Each primitive supports:

- `is_inside(vol, pos)`: point containment test
- `distance_to_exit(vol, pos, dir)`: ray exit distance
- `distance_to_entry(vol, pos, dir)`: ray entry distance

The default geometry models a NEXT-like high-pressure GXe TPC. The
geometry is not hardwired — changing the detector is a JSON-only
operation.


## Output

Energy deposits are binned into 3D voxels of configurable size
(default: 2 mm). Each event produces a list of occupied voxels with
(ix, iy, iz, n_electrons) where n_electrons = E_deposited / W_gas and
the voxel center is reconstructed as:

    x = (ix + 0.5) * voxel_mm   [mm]

Output format is HDF5 with a flat layout:

```
/metadata/
    voxel_mm            float64
    pressure_bar        float64
    temperature_K       float64
    gas                 string
    Qbb_eV              float64
/events/
    event_id            (N,) int32
    n_voxels            (N,) int32
    n_electrons         (N,) int32
    status              (N,) string     # "accepted", "killed_energy", ...
    voxel_offset        (N,) int64      # index into /voxels
/voxels/
    ix                  (M,) int16
    iy                  (M,) int16
    iz                  (M,) int16
    n                   (M,) int16
```

Random access to event k: read voxels[voxel_offset[k]:voxel_offset[k]+n_voxels[k]].
