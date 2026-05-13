# LZ Background Status — Annual Rates by Source Group

Annual background counts in the ±1σ ROI (2433–2482 keV) and inner 967 kg
fiducial volume, grouped by source category. Derived from Table I of the
LZ 0νββ paper (`tab:backgrounds` in `docs/lz_bb0nu_summary.tex`) by
scaling the 1000-day counts to annual rates (factor 365.25/1000 = 0.36525).
Activities marked † are upper limits. "Others" bundles the LZ outer
detector (liquid-scintillator veto, *not* the cryostat), miscellaneous
detector components, cavern walls, and intrinsic backgrounds
(neutron-induced ¹³⁷Xe, internal ²²²Rn, 2νββ, ⁸B solar ν).

Assumptions: 1.0% energy resolution at Q_ββ, 3 mm single-scatter
rejection in z. Mass and specific activities from `tab:backgrounds`;
per-year counts from `tab:rates`.

| Component | Mass (kg) | ²³⁸U-late (mBq/kg) | ²³²Th-late (mBq/kg) | U (cts/yr) | Th (cts/yr) | Total (cts/yr) |
|:----------|----------:|-------------------:|--------------------:|-----------:|------------:|---------------:|
| _PMTs_ | | | | | | |
| TPC PMTs              |   91.9 |    3.22 |    1.61 | 1.077 | 0.037 | 1.114 |
| TPC PMT bases         |    2.80 |   75.9  |   33.1  | 0.555 | 0.011 | 0.566 |
| TPC PMT structures    |  166   |    1.60 |    1.06 | 0.968 | 0.044 | 1.012 |
| TPC PMT cables        |   88.7 |    4.31 |    0.82 | 0.526 | 0.069 | 0.596 |
| Skin PMTs + bases     |    8.59 |   46.0  |   14.9  | 0.274 | 0.007 | 0.285 |
| **PMTs subtotal**     |        |         |         | **3.400** | **0.168** | **3.573** |
| _Field Cage_ | | | | | | |
| PTFE walls            |  184   |    0.04 |    0.01 | 0.142 | 0.000 | 0.142 |
| TPC sensors           |    5.02 |    5.82 |    1.88 | 0.435 | 0.000 | 0.435 |
| Field grids + holders |   89.1 |    2.63 |    1.46 | 0.226 | 0.040 | 0.267 |
| FC resistors          |    0.06 | 1350    | 2010    | 0.961 | 0.011 | 0.968 |
| FC rings†             |   93.0 |    0.35 |    0.24 | 0.300 | 0.000 | 0.300 |
| **Field Cage subtotal** |      |         |         | **2.064** | **0.051** | **2.112** |
| _Cryostat_ | | | | | | |
| Ti cryostat†          | 2590   |    0.08 |    0.22 | 0.475 | 0.073 | 0.544 |
| Cryo. insulation†     |   13.8 |   11.1  |    7.79 | 0.329 | 0.015 | 0.343 |
| **Cryostat subtotal** |        |         |         | **0.804** | **0.088** | **0.887** |
| _Others_ | | | | | | |
| Outer detector†       | 22900  |    4.71 |    3.73 | 0.621 | 0.395 | 1.019 |
| Other components      |  438   |    1.83 |    1.65 | 0.767 | 0.113 | 0.881 |
| Cavern walls          | 29000  |    3.21 |    8.41 | 1.170 | 3.060 | 4.240 |
| ¹³⁷Xe (neutron-induced) | —    | —       | —       | —     | —     | 0.102 |
| Internal ²²²Rn        | —      | —       | —       | —     | —     | 0.164 |
| ¹³⁶Xe 2νββ            | —      | —       | —       | —     | —     | 0.004 |
| ⁸B solar ν            | —      | —       | —       | —     | —     | 0.011 |
| **Others subtotal**   |        |         |         | **2.558** | **3.568** | **6.421** |
| **Grand total**       |        |         |         | **8.826** | **3.875** | **12.99** |

## Group ranking

By total annual rate in the ROI:

1. **Others**: 6.42 cts/yr — dominated by cavern walls (4.24 cts/yr) and
   the LZ outer-detector liquid scintillator (1.02 cts/yr).
2. **PMTs**: 3.57 cts/yr — spread roughly evenly across TPC PMTs,
   structures, and cables.
3. **Field Cage**: 2.11 cts/yr — dominated by FC resistors (0.97 cts/yr)
   despite their tiny mass (60 g), reflecting the very high specific
   activity (1350 mBq/kg ²³⁸U-late, 2010 mBq/kg ²³²Th-late).
4. **Cryostat**: 0.89 cts/yr — the cleanest detector subsystem by total
   contribution.

The grand total 12.99 cts/yr matches the LZ projection of ~13 background
events per year in the ±1σ ROI.

## LXeMC vs LZ paper — current PMT-group rates (2026-05-13)

After the PMT source-geometry refactor (branch `sources`, commits
`2dc8c7a`...`e3989c4`) the LXeMC predictions for the PMT-related groups
are:

| Component | LZ paper (cts/yr) | LXeMC (cts/yr) | MC / LZ |
|:----------|------------------:|---------------:|--------:|
| TPC PMTs | 1.077 | 0.684 | 0.64 |
| TPC PMT bases | 0.555 | 0.491 | 0.88 |
| TPC PMT structures | 0.968 | 0.622 | 0.64 |
| **TPC PMT cables** | **0.526** | **0.482** | **0.92** |
| **Skin PMTs + bases** | **0.274** | **0.270** | **0.99** |
| **PMT-group total** | **3.400** | **2.549** | **0.75** |

**Cables and Skin PMTs+bases match LZ paper to within 10%**, the
result of the May 2026 geometry refactor (cables relocated to upper/
lower conduits, Skin PMTs localised on the field-cage outer side and
near the cathode-level ring).

**TPC PMTs, bases, and structures all under-predict by 12–36%** with
similar fractions across the three components. This is a separate
open issue: all three have correct positions in the source geometry
and match LZ activities, so the discrepancy likely shares a single
cause (e.g. attenuation modelling of the PMT-array structural
material, or differing analysis cuts). Not yet investigated.
