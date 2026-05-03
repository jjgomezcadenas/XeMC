# LXe Monte Carlo — Pedagogical implementation for SS/MS in 0νββ searches

A self-contained Monte Carlo simulation of MeV-scale γ-rays in liquid
xenon, with application to single-site / multi-site discrimination
near the ²³⁶Xe Q_ββ. Companion to the LaTeX manual `manual/lxe_manual.pdf`.

## Layout

```
lxe_mc/
├── README.md                      ← this file
├── manual/
│   ├── lxe_manual.tex             ← LaTeX source
│   ├── lxe_manual.pdf             ← compiled manual (18 pages)
│   └── figs/                      ← figures (PDF for LaTeX, PNG for direct view)
├── code/
│   ├── physics.py                 ← analytic cross sections (KN, BH-pair,
│   │                                Sandia, Berger-Seltzer, Tsai, BH-brems)
│   ├── nist_data.py               ← loads NIST XCOM and ESTAR tables, returns
│   │                                log-log interpolators (sigma_X_NIST,
│   │                                dEdx_*_NIST, csda_range_LXe_NIST)
│   ├── sampling.py                ← MC samplers: KN, BH-pair, BH-brems,
│   │                                Sauter, atomic relaxation
│   ├── tracking.py                ← Track, Stack, transport_photon,
│   │                                transport_lepton, simulate_event,
│   │                                cluster_deposits_in_z
│   ├── thresholds.py              ← all numerical parameters (E_cut, T_cut,
│   │                                k_min, f_range, etc.)
│   ├── make_figures.py            ← generates all manual figures and CSV summaries
│   ├── validation.py              ← 9 validation tests; all pass
│   └── tests/                     ← (placeholder for unit tests)
└── data/
    ├── xcom_xe.csv                ← NIST XCOM Xe table (verbatim, with
    │                                provenance header)
    ├── estar_xe.csv               ← NIST ESTAR Xe table (verbatim)
    ├── sandia_xe_coeffs.csv       ← Geant4 G4StaticSandiaData coefficients
    │                                for Z=54 (verbatim)
    ├── photon_xsec_summary.csv    ← summary table at key energies
    └── electron_stopping_summary.csv
```

## How to use

### Compile the manual

```bash
cd manual
pdflatex lxe_manual.tex
pdflatex lxe_manual.tex   # twice for cross-references
```

Requires standard LaTeX (TeX Live with `physics`, `siunitx`, `algorithm2e`,
`tcolorbox`, `cleveref`, `hyperref`, `booktabs`).

### Run the code

Each module is runnable and self-tests when run directly:

```bash
cd code
python physics.py        # prints cross-section / stopping-power values
python nist_data.py      # prints NIST-anchored values
python sampling.py       # not directly runnable; imported by tracking
python tracking.py       # runs 100 events of 2.615 MeV in infinite LXe
python make_figures.py   # regenerates all figures and CSVs
python validation.py     # runs the 9 validation tests
```

Requires Python 3.10+ and `numpy`, `matplotlib` only.

### Simulate one event

```python
import numpy as np
from tracking import simulate_event, cluster_deposits_in_z, is_single_site
rng = np.random.default_rng(42)
deposits = simulate_event(2.615, position=(0,0,0), direction=(0,0,1), rng=rng)
total_E = sum(d.energy for d in deposits)
ss = is_single_site(deposits, dz_cm=0.30)
print(f"Total E deposited: {total_E:.4f} MeV; SS: {ss}")
```

## Key results

From the validation suite (`python validation.py` → 9/9 PASS):

| quantity | value | source |
|----------|-------|--------|
| photoelectric branching at Q_ββ in Xe | 2.32% | NIST XCOM |
| Compton branching at Q_ββ in Xe | 85.0% | NIST XCOM |
| pair branching at Q_ββ in Xe | 12.7% | NIST XCOM |
| photon mfp at Q_ββ in LXe | 8.81 cm | NIST XCOM |
| electron CSDA range at 1 MeV in LXe | 2.42 mm | NIST ESTAR |
| electron CSDA range at 2 MeV in LXe | 5.12 mm | NIST ESTAR |
| radiative fraction at 1 MeV in LXe | 6.8% | NIST ESTAR |
| radiative fraction at 2 MeV in LXe | 11.8% | NIST ESTAR |
| MC SS fraction for ²⁰⁸Tl in inf. LXe | ~10% | this MC |

## Notable corrections from earlier draft

The previous version of the manual claimed photoelectric absorption
was "few-per-mille" at our energies. **This was wrong by a factor of
~30**: NIST XCOM gives 2.32% at Q_ββ. The current manual treats this
correctly and includes a dedicated section (§9, "Channel-by-channel
SS-mimic mechanism") explaining that photoelectric is the dominant
direct full-energy SS-mimic mode for ²¹⁴Bi.

## Validity caveats

- The Sandia photoelectric formula is used **only below 500 keV**; above
  it extrapolates badly (3.6× too high at 1 MeV in Xe). Above 500 keV
  we interpolate NIST XCOM directly. This dual approach is documented
  and tested.
- We omit multiple scattering (~20% effect on projected electron range),
  δ-rays, energy-loss straggling, Doppler broadening of Compton
  electrons, and Rayleigh scattering, all of which are sub-leading for
  SS/MS classification at 3 mm z-resolution. See manual §13.
- Density-effect correction in ESTAR is for elemental Xe; the
  difference for liquid Xe is < 2% at MeV energies.

## Provenance

- NIST XCOM: https://physics.nist.gov/xcom/ (SRD 8, 2010)
- NIST ESTAR: https://physics.nist.gov/Star/ (SRD 124, 2005)
- Geant4 PRM: https://geant4.web.cern.ch/documentation/dev/prm_html/ (v11.4, 2024)
- Biggs & Lighthill, Sandia SAND87-0070 (May 1990)

Data retrieved 2026-05-03.

## License

Public domain / CC0. Use freely, attribute if useful.
