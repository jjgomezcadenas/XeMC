# Fix: source-geometry for TPC PMT cables and lower-ring Skin PMTs

## Motivation

The MC over-predicts the LZ paper's per-component Bi-214 background rate
for two source components by a large factor:

| Component | LZ paper ²³⁸U-late (cts/yr) | MC Bi-214 (cts/yr) | MC/LZ |
|:----------|---------------------------:|-------------------:|------:|
| TPC PMT cables | 0.526 | 1.320 | 2.51 |
| Skin PMTs + bases | 0.274 | 1.368 | 4.99 |

Source activities (mBq/kg) and component masses match the LZ paper
exactly. The fast-kernel Skin veto (100 keV, `cfg.veto_skin`) and TPC
veto (10 keV, `cfg.veto_TPC`) match the LZ paper, are correctly wired
in `tracking_fast.jl`, and pass the tag-coverage test suite added in
the `refactor/tag-dispatch` branch.

The discrepancy is in the **source-side geometry** in
`data/source_geometry_lz_v1.json`, which places radioactive material
in positions that are *not* consistent with the LZ TDR. This document
proposes the changes to bring the source-geometry into agreement.

## Findings from the LZ TDR

### TPC PMT cables

LZ TDR Section 3.4.4 (line 4706 of `docs/LZ-TDR.txt`):

> "The cables from the PMTs associated with the upper and lower parts
> of the TPC are housed in separate conduits, so that **no cabling is
> routed through the side Skin region** (these conduits can be seen in
> Figure 1.2.1 in Chapter 1). This could interfere with the ability to
> hold a high voltage on the cathode and it would degrade the light
> collection efficiency in the Skin."

Cable inventory (LZ TDR Table 3.4.1, line 4716):

| Conduit | TPC PMT cables | Skin PMT cables | Other | Total |
|:--------|--------------:|----------------:|------:|------:|
| Upper | 506 | 186 | 109 | 801 |
| Lower | 482 | 76 |  67 | 625 |

Typical lengths: 12.8 m (upper), 11.6 m (lower). Total ≈ 17 km, mass
≈ 72 kg in the xenon space (LZ TDR line 4760). The LZ paper's 88.7 kg
"TPC PMT cables" entry includes additional ~17 kg of connector /
feedthrough / sheathing material.

### Skin PMTs

The LZ paper's "Skin PMTs + bases" entry (8.59 kg) is the *sum* of
three physically distinct populations (LZ TDR table at line 3754):

| Type | Count | Position | Mass (MC) |
|:-----|------:|:---------|----------:|
| R8520 "Top Skin" 1" | 93 | Side Skin, distributed cylindrically | 6.10 kg |
| R8778 "Bottom Skin" side 2" | 20 | **Cathode level only** (lower side ring) | 1.31 kg |
| R8778 "Dome" 2" | 18 | **Below cathode**, dome region (passive LXe) | 1.18 kg |

Side Skin thickness (LZ TDR table line 3789): **4 cm at top, 8 cm at
cathode**. Asymmetric.

## Current source-geometry (problematic)

From `data/source_geometry_lz_v1.json`:

### `PMT_BARREL_cables` (line 366)
```json
"shape": "cylinder_shell",
"R_inner_cm": 81.2,
"wall_thickness_cm": 0.150,
"half_height_cm": 79.675,
"position_cm": [0.0, 0.0, 65.925],
"equivalent_mass_kg": 88.7
```
Distributes 88.7 kg of cables as a thin shell at R=81.2 cm spanning
z ∈ [-13.75, 145.6] cm — i.e., **inside the side Skin region along the
entire barrel**. This is the opposite of the LZ TDR design.

### `PMT_BARREL_R8520` (line 384)
```json
"R_inner_cm": 81.4, "half_height_cm": 79.675, "position_cm": [0,0,65.925],
"equivalent_mass_kg": 6.10
```
Distributed cylindrically over the full side Skin z range. **Matches
LZ** (the 93 R8520 PMTs are distributed around the cylindrical side
Skin per TDR line 2350 / 4528).

### `PMT_BARREL_R8778_lower` (line 402)
```json
"R_inner_cm": 81.4, "half_height_cm": 79.675, "position_cm": [0,0,65.925],
"equivalent_mass_kg": 1.31
```
Same geometry as R8520, but should be **localized at the cathode level
only** (20 PMTs in a single lower ring near z ≈ -13.75 cm).

### `PMT_BOT_R8778_dome` (line 349)
```json
"shape": "cylinder", "radius_cm": 72.8, "half_height_cm": 0.0045,
"position_cm": [0.0, 0.0, -15.7545],
"equivalent_mass_kg": 1.18
```
Disk source at z = -15.75 cm. **Matches LZ** (18 dome PMTs below the
cathode). Note: this region is `TAG_PASSIVE_LXE` in the tracking
geometry, *not* `TAG_SKIN`, so the Skin veto cannot fire on its
gammas. This is the correct physical model — the dome PMT background
is reduced only by LXe path length and the TPC veto.

## Proposed changes

### 1. Split `PMT_BARREL_cables` into upper + lower conduits

Replace the single `PMT_BARREL_cables` shell with **two compact source
volumes** placed outside the side Skin, mimicking the LZ conduit
geometry.

**Mass split**: by cable-meters. Upper conduit carries 801 cables at
12.8 m → 10,253 cable-m. Lower carries 625 cables at 11.6 m → 7,250
cable-m. Ratio ≈ 58.6 / 41.4 (upper / lower). Total 88.7 kg:
- Upper conduit mass: 52.0 kg
- Lower conduit mass: 36.7 kg

**Position decisions needed** (placeholders below — confirm before
implementation):

Option A (simplest — compact disks near each PMT array):
```json
{
  "name": "PMT_TOP_cables",
  "shape": "cylinder",
  "radius_cm": 72.8,
  "half_height_cm": 1.0,
  "position_cm": [0.0, 0.0, <z_above_top_PMT_array>],
  "equivalent_mass_kg": 52.0,
  "activity": { "Bi214_mBq_per_kg": 4.31, "Tl208_mBq_per_kg": 0.82 }
},
{
  "name": "PMT_BOT_cables",
  "shape": "cylinder",
  "radius_cm": 72.8,
  "half_height_cm": 1.0,
  "position_cm": [0.0, 0.0, <z_below_bottom_PMT_array>],
  "equivalent_mass_kg": 36.7,
  "activity": { "Bi214_mBq_per_kg": 4.31, "Tl208_mBq_per_kg": 0.82 }
}
```

Option B (more faithful — narrow vertical conduit cylinders offset
radially outside the side Skin):
- Upper conduit: small cylinder at R ≈ 82.5 cm (outside Skin shell),
  z range above the top PMT array.
- Lower conduit: small cylinder at R ≈ 82.5 cm, z range below the
  bottom PMT array.

**Recommendation: Option A** for the first iteration. The cables exit
the cryostat through standoffs above the top array; their integrated
solid angle into FV is well approximated by a disk near the top PMT
array. Once Option A is calibrated against the LZ paper, Option B can
be tried if residual discrepancy demands it.

### 2. Localize `PMT_BARREL_R8778_lower` to the cathode-level ring

Reduce the z extent from full barrel to a narrow ring at the cathode
level. The 20 R8778 PMTs are mounted around the bottom of the side
Skin per the LZ TDR (Section 3.7, line 827).

```json
{
  "name": "PMT_BARREL_R8778_lower",
  "shape": "cylinder_shell",
  "R_inner_cm": 81.4,
  "wall_thickness_cm": 0.150,    // increase thickness to keep mass density realistic for a thin ring
  "half_height_cm": 5.0,         // ~10 cm tall ring (was 79.675)
  "position_cm": [0.0, 0.0, -8.75],  // center at z = -8.75 cm (range -13.75 to -3.75, just above cathode)
  "equivalent_mass_kg": 1.31,
  "activity": { "Bi214_mBq_per_kg": 46.0, "Tl208_mBq_per_kg": 14.9 }
}
```

Note the `wall_thickness_cm` increase from 0.00222 to 0.150 cm is
cosmetic — the source is "virtual_source"/"transparent" so the
material thickness doesn't affect transport. But the half-height
change from 79.675 to 5.0 cm is the load-bearing fix: it localizes
the radioactivity at z ≈ -8.75 ± 5 cm instead of smearing it over
z ∈ [-13.75, 145.6].

**Decision needed**: confirm cathode-ring z range. The cathode plane
is at z = -13.75 cm in our geometry (top of `LXe_below_cathode`). The
R8778 ring would sit just above the cathode, so a center at z ≈ -8.75
with half-height 5 cm (range -13.75 to -3.75 cm) is a reasonable
guess. Could also be tighter (half-height 2 cm) if the LZ TDR has a
more precise position.

### 3. `PMT_BOT_R8778_dome` — no change

The 18 dome PMTs at z = -15.75 cm are correctly placed. Their gammas
traverse `LXe_below_cathode` (`TAG_PASSIVE_LXE`), the cathode region,
and TPC active LXe before reaching FV. The Skin veto cannot fire on
them — this is correct physics, not a bug. The TPC active veto
(10 keV) applies in the active LXe region above the cathode.

If the dome PMTs turn out to dominate the residual "Skin PMTs"
discrepancy after fixes (1) and (2), investigate whether their
effective LXe path length matches LZ TDR — i.e., whether the
`LXe_below_cathode` region depth and the cathode-to-FV distance match
the TDR geometry.

### 4. `PMT_BARREL_R8520` — no change

The 93 R8520 PMTs are cylindrically distributed along the side Skin
in the LZ TDR. Our shell at R = 81.4, full barrel z, is consistent
with this.

## Validation plan

1. **Re-run the Bi-214 background pipeline** for pmt_top, pmt_bottom,
   pmt_barrel sources with the new geometry. Targets:
   - TPC PMT cables: 0.526 cts/yr (LZ paper)
   - Skin PMTs + bases: 0.274 cts/yr (LZ paper)
2. **Compare per-component breakdown** in
   `results/bfv/pmt/Bi214_summary/pmt_background_summary.txt`.
3. **Check the other PMT components are unchanged**:
   - TPC PMTs: 1.077 cts/yr (LZ) vs 0.684 cts/yr (current MC) — the
     ~36% low discrepancy here is a separate issue and should not be
     affected by these source-geometry changes.
   - TPC PMT bases: 0.555 vs 0.491 — ditto.
   - TPC PMT structures: 0.968 vs 0.622 — ditto.

## What this fix does NOT address

- The ~36% under-prediction for `TPC PMTs` and `TPC PMT structures`.
  These have correct (matching LZ) positions in the source geometry.
  The under-prediction is likely due to (a) different attenuation /
  ROI window assumptions, (b) different LXe path-length distribution
  between LZ TDR geometry and ours, or (c) a different "TPC PMTs"
  component mass attribution. A separate investigation.

- The LZ paper's veto model may differ from our threshold-only Skin
  veto in subtle ways (multi-deposit coincidence, time gating). We
  retracted this earlier as a likely explanation but it remains a
  small residual source of discrepancy.

## Test plan

After applying the JSON changes:
1. Re-run `Pkg.test()` — should be unaffected (source-geometry changes
   don't touch the Julia transport code).
2. Re-run flux table generation:
   `julia --project=. scripts/generate_flux_tables.jl --source pmt_top --isotope Bi214 ...`
   etc. for `pmt_bottom`, `pmt_barrel`.
3. Re-run background pipeline (Stage 1 + Stage 2 + Stage 3) for each
   PMT source × Bi-214.
4. Re-run `py/pmt_background_summary.py` and compare against the LZ
   paper numbers (target table above).
5. If totals agree to within ~30%, accept the fix. If not, iterate
   (try Option B for cables, tighten the R8778_lower ring z extent,
   etc.).

## Open decisions for user review

Before drafting the JSON diff:

1. **Cable conduit z positions** — what z values for the disk
   positions in Option A? Suggestion: top conduit at z just above
   `PMT_TOP_PMTs` position (currently at z ≈ 146 cm); bottom conduit
   just below `PMT_BOT_PMTs` (currently at z ≈ -17 cm). Need the user
   to confirm exact placement or provide TDR figure reference.

2. **Cables: Option A (compact disks) or Option B (offset cylinders
   outside Skin)** for the first iteration. Recommendation: A.

3. **R8778_lower ring z extent** — is the cathode-level half-height
   ~5 cm the right scale, or should it be tighter (~2 cm)?

4. **`PMT_BARREL_cables` removal**: delete the existing entry, or keep
   it commented out (`"_obsolete": true`) for traceability? Convention
   in this codebase appears to be deletion (per the
   `refactor/tag-dispatch` commit pattern).

## Scope

Source-geometry-only change. No Julia code touched. No tracking
geometry changes. Should be a single JSON edit + re-run of the
background pipeline.

---

## Step 2 update — `PMT_SKIN_LOWER_RING` (R8778 lower-ring localization)

**Applied**. The R8778 lower-ring component was previously named
`PMT_BARREL_R8778_lower` and was distributed as a thin cylindrical
shell over the entire side-Skin z range. The LZ TDR (table at line
3754 of `docs/LZ-TDR.txt`) labels these PMTs "Bottom Skin: 20 side +
18 dome", implying the 20 side PMTs are mounted in a single ring at
the bottom of the side Skin, not smeared along its height.

Changes made:

- `data/source_geometry_lz_v1.json`: renamed
  `PMT_BARREL_R8778_lower` → `PMT_SKIN_LOWER_RING`; relocated to a
  narrow cylindrical-shell ring (R=81.4, half_height=5 cm, center
  z=-8.75 cm).
- `src/pmt_sources.jl`: removed `PMT_BARREL_R8778_lower` from
  `pmt_barrel_flux` (now single-component R8520). Added
  `pmt_skin_lower_ring_flux`.
- `src/source_dispatch.jl`: new `pmt_skin_lower_ring` endpoint;
  removed the relocated component from the `pmt_barrel` VE
  composition.
- Python summary + validators updated.
- Tests: `pmt_barrel flux` testset now 1-component; new
  `pmt_skin_lower_ring flux` testset.

## **PENDING TDR REVIEW** — open geometry questions to resolve together

The following two questions should be resolved in a single review of
LZ TDR section 3.7 ("Mounting of dome and lower side Skin PMTs",
referenced at line 827 of `docs/LZ-TDR.txt`). They are deliberately
parked together because §3.7 is expected to cover both the lower ring
and the upper Skin PMTs:

### Q1. `PMT_SKIN_LOWER_RING` exact z range

Current values are an **initial estimate**:

| Parameter | Current | Source |
|:----------|--------:|:-------|
| `half_height_cm` | 5.0 | placeholder — 10 cm ring |
| `position_cm` (z) | -8.75 | placeholder — 5 cm above cathode (z=0), 5 cm below FC bottom (z=-13.75) |

What to verify in §3.7: actual mounting plane z, ring thickness, and
whether the 20 PMTs span the full circumference or a partial arc.

### Q2. `PMT_BARREL_R8520` ("Top Skin") z distribution

The MC currently models the 93 R8520 PMTs as cylindrically distributed
over the entire side-Skin z range (R=81.4 cm, `half_height_cm=79.675`).
The LZ TDR has *conflicting* statements:

- §1.3.5.1 (line 2350): "side skin is equipped with 90 R8520 looking
  downwards and another 90 looking upwards" — implies 180 PMTs in
  two rings (one looking down, one looking up).
- Spec table (line 3759): "Top Skin (R8520-406): 93" — implies a
  single population. The label "Top" suggests they sit at the *top*
  of the side Skin, not distributed along it.

If the 93 R8520 are actually localized at the top of the side Skin
(a ring near z ≈ 130–145 cm, where the Skin is thinnest — 4 cm per
the TDR), gammas would have a shorter LXe path to the FV through the
side, which would *increase* the predicted FV background — opposite
to the direction needed to close the 5× discrepancy. This may indicate
the discrepancy has a different cause (e.g., dome contribution
modelling, or the LZ paper applies additional cuts not represented in
our MC).

**Decision deferred** until §3.7 is read.

### Action items for the joint TDR review

1. Read LZ TDR section 3.7 in full.
2. Confirm or revise:
   - `PMT_SKIN_LOWER_RING` `half_height_cm` and `position_cm[3]`.
   - `PMT_BARREL_R8520` distribution (full barrel z vs. localized ring).
3. If R8520 needs localization, repeat the step-2 pattern:
   rename to e.g. `PMT_SKIN_UPPER_RING`, drop from `pmt_barrel_flux`,
   add `pmt_skin_upper_ring_flux` endpoint, update Python and tests.
4. Re-run the Bi-214 background pipeline and compare against LZ paper
   "Skin PMTs + bases" = 0.274 cts/yr.

---

## Step 3 update — `PMT_SKIN_UPPER_RING` (R8520 relocation, `pmt_barrel` deletion)

**Applied**. The TDR review was carried out and crucially supplemented
by the more recent LZ instrumentation paper
(arXiv:1910.09124v2, §2.2 "The Xe Skin Detector"). Key quotes:

> "The side region of the skin contains 4 cm of LXe at the top, widening
> to 8 cm at cathode level for increased standoff distance. This is
> viewed from above by 93 1-inch Hamamatsu R8520-406 PMTs. These are
> retained within PTFE structures **attached to the external side of
> the field cage**, located **below the liquid surface**."
>
> "At the bottom of the detector a **ring structure attached to the
> vessel** contains a further 20 2-inch Hamamatsu R8778 PMTs viewing
> upward into this lateral region."

The 93 R8520 PMTs are therefore **not** cylindrically distributed
along the side Skin. They are localised in a thin ring **on the
external side of the field cage** (R near 72.8 cm, not 81.4 cm), just
below the LXe surface (z ≈ 144 cm). The previous MC model placed all
6.10 kg of R8520 mass uniformly along R=81.4 cm over the full barrel
z range — wrong on both R and z.

Changes made:

- `data/source_geometry_lz_v1.json`: renamed `PMT_BARREL_R8520` →
  `PMT_SKIN_UPPER_RING`; relocated to a narrow cylindrical-shell ring
  at `R=74.224 cm` (just outside the FC outer wall at R=74.224),
  `half_height=1.5 cm`, `position=[0, 0, 144]` (1.6 cm below the LXe
  surface).
- `PMT_SKIN_LOWER_RING` also tightened to `half_height=1.5 cm` and
  `position=[0, 0, -12]` per the same paper ("ring structure attached
  to the vessel at the bottom of the detector").
- `src/pmt_sources.jl`: **deleted** `pmt_barrel_flux` (the barrel
  source now has zero components). Added `pmt_skin_upper_ring_flux`.
- `src/source_dispatch.jl`: **removed** `pmt_barrel` endpoint from
  `dispatch_source_flux`, `make_virtual_envelope`, and
  `supported_sources()`. Added `pmt_skin_upper_ring`.
- `src/LXeMC.jl`: removed `pmt_barrel_flux` export; added
  `pmt_skin_upper_ring_flux`.
- Tests, Python summary/validators, viewer, CLAUDE.md,
  status_sources.md, running_recipe.md: all references updated.

After this step, all the geometry-related "PENDING TDR REVIEW"
questions above are resolved. The MC's Skin-PMT geometry now matches
the published LZ description. Any residual MC-vs-paper discrepancy in
the "Skin PMTs + bases" rate after re-running the pipeline must come
from a different source (e.g., veto efficiency modelling, LZ paper
applying analysis cuts beyond the prompt veto, or modelling of the
PMT_BOT_R8778_dome attenuation path).

---

## Step 4 update — TPC PMT cables: use cold cable mass only

**Applied**. After Step 1, the MC predicted 1.32 cts/yr for TPC PMT
cables vs the LZ paper's 0.526 cts/yr — a 2.5x over-prediction. After
the Step 1 + 3 geometry refactors a fresh pipeline run gave 2.23
cts/yr (worse — 4.2x over). Root cause: the equivalent mass on the
cable cylinders was the *full* 88.7 kg from LZ paper Table I, which
includes both the cold cable segment inside the cryostat *and* the
warm cable segment outside.

LZ TDR §3.4.4 (line 4763) breaks the cable length down:

> "Rn emanation from cables is expected to be dominated by the warm
> region, which is ≈8 m long for the upper routing, and only ≈1 m
> for the lower routing."

With typical per-cable lengths 12.8 m upper, 11.6 m lower, the cold
fractions are:

| Routing | Per-cable cold length | Cable count | Cold cable-meters |
|:--------|----------------------:|------------:|------------------:|
| Upper   | 4.8 m                 | 506         | 2,429             |
| Lower   | 10.6 m                | 482         | 5,109             |
| **Total** | —                    | —           | **7,538**         |

Total cable-meters (warm + cold) = 12,068. Cold fraction = 7,538 /
12,068 = **62.5%**.

The warm cable mass sits *outside* the cryostat (and outside our MC
tracked volume entirely), shielded by the cryostat walls and OD
scintillator. Its contribution to FV background is negligible. So
the relevant mass for our MC is the cold portion only:

- Cold total: 88.7 x 0.625 = 55.4 kg
- Upper cold: 55.4 x (2429/7538) = **17.9 kg** (was 45.4 kg)
- Lower cold: 55.4 x (5109/7538) = **37.5 kg** (was 43.3 kg)

Note the redistribution: under cable count alone (Step 1) the upper
routing carried ~half the cable mass; with cold-meter accounting the
lower routing carries ~68% (because the lower routing has very little
warm cable).

Changes made:

- `data/source_geometry_lz_v1.json`: `PMT_TOP_cables.equivalent_mass_kg`
  45.4 -> 17.9; `PMT_BOT_cables.equivalent_mass_kg` 43.3 -> 37.5.
  `_doc` fields extended with the cold-fraction justification.
- `py/util/validate_geometry_vs_paper.py`: aggregate check for "TPC
  PMT cables" updated from 88.7 kg (LZ paper Table I, includes warm)
  to 55.4 kg (cold only) with explanatory comment.
- `py/util/validate_activities_vs_paper.py`: unchanged (validator
  only checks mBq/kg activities, not masses).
- Tests: unchanged (mass values are not asserted in the test suite).

After this step the cold cable mass is the only relevant component
in the MC and matches the physical setup: gammas from warm cable
outside the cryostat cannot reach the FV without losing most of their
energy to the cryostat + OD, so omitting them is the right
approximation.

---

## Step 5 update — `PMT_TOP_cables` radial position: TPC edge approximation

**Applied**. After Steps 1, 3, and 4 the MC predicted 2.23 cts/yr
for TPC PMT cables vs the LZ paper's 0.526 cts/yr (4.2x over). After
the cold-mass adjustment in Step 4 this scales down by 0.625 ->
predicted ~1.4 cts/yr (~2.6x over), still high. The residual cause is
the *radial* position of `PMT_TOP_cables`: the previous model placed
the upper cable bundle as a thin cylinder at R=10 cm on the central
axis, giving its gammas a clean straight-down line of sight to the
FV through ~75 cm of LXe.

LZ paper arXiv:1910.09124v2 §2.4 (line 558):

> "Two ports at the top head and the central port at the bottom are
> for the PMT and instrumentation cables."

The bottom cables therefore *do* exit through the central axis, but
the upper cables exit through two off-center ports. The cold cable
segment for the upper routing physically extends from the top of the
PMT array up to those off-center ports. The cables are routed along
the TPC periphery rather than down the central axis.

### Implementation choice: Option A (wide disk)

To represent this without adding new dispatch plumbing, the upper
cable cylinder is widened from `radius=10 cm` (central axis only) to
`radius=72.8 cm` (full TPC top region). The volume integration weight
`dV = r dr dphi dz` makes the average radial position
`<r> = (2/3) x 72.8 = 48.5 cm` — most of the cable radioactivity is
sampled at moderate-to-large r, which is closer to the physical
off-center port routing than the previous central-axis placement.

The half_height is reduced to 5.5 cm so the cylinder fits inside the
AirDome at r=72.8 (dome ceiling there is z=167.5). Position center
z=159.5 cm; z range [154, 165].

Limitations of Option A:

- The fill is *uniform* across r in [0, 72.8], so some cable mass is
  still placed near the central axis where cables don't physically
  run. The volume weighting biases against this but doesn't eliminate
  it.
- Azimuthal symmetry is enforced even though the real cables emerge
  from two specific ports, not a uniform ring. For the integrated FV
  rate this is the right average (FV itself is azimuthally
  symmetric), so the approximation is justified.

A more faithful Option B (thin cylindrical-shell ring at R=72.8 with
axial cos_theta_to_lxe filter) would require either a new
`cos_theta_to_lxe(::PCylShell, ...)` axial method or a custom flux
function bypassing the existing `_merge_pmt_volume` /
`_merge_pmt_barrel_volume` helpers. Deferred unless Option A still
over-predicts.

### Changes made

- `data/source_geometry_lz_v1.json`: `PMT_TOP_cables.radius_cm`
  10 -> 72.8; `half_height_cm` 17.5 -> 5.5;
  `position_cm[3]` 171.5 -> 159.5. The `_doc` field is extended with
  the LZ paper port reference and the Option A approximation
  justification.
- `PMT_BOT_cables`: unchanged. The LZ paper explicitly places the
  lower-conduit cable bundle on the central axis at the bottom flange
  of the ICV, matching our existing model.

Mass, activity, dispatch, and Julia/Python pipeline plumbing are all
unchanged from Step 4.
