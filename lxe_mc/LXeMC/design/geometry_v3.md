# Geometry V3

## Scope

V3 is the live geometry schema. It encodes:

- a partitioned, non-overlapping tracking tree in `data/detector_lz_v3.json`
- a separate, anatomically anchored source-side reference in
  `data/source_geometry_lz_v1.json`
- explicit mother/daughter topology with strict containment
- tag-based simulation policy + role-based descriptive metadata
- explicit modeling status via approximation labels

The previous V2 schema is archived in `legacy/`.

## Two files, one geometry

**Source and tracking live in independent JSON files**, but they
describe the **same anatomical geometry** wherever they overlap. The
two-file split exists for a different reason: each file is a *subset*
of the full set of detector volumes, optimized for what its consumer
needs. It is not a license to draw the same volume with different
dimensions on the two sides — and we don't.

The split rule is:

- a volume that **hosts radioactivity** must be in the **source file**
  (with mass + activities)
- a volume that **affects photon transport in the tracking stage**
  must be in the **tracking file** (with material + tag + role)
- volumes in **both** files share the same R, z, wall thickness, etc.
  so the geometric-mass invariant validates both at once

Why the split is allowed (vs. forcing one mixed file):

- The source file lists Ti cryostat shells, MLI, PMT structures and
  cables — none of which the photon transport needs to see, because
  source-side propagation has already absorbed their effect into the
  emission spectrum at the detector boundary. Putting them in the
  tracking file would just slow the per-step volume-finder for no
  physics gain.
- The tracking file lists analysis-policy regions (`FV`,
  `BarrelActive`, `Skin`, …) that are pure tags on LXe — they don't
  host radioactivity, they just route deposits to different cuts.
  Putting them in the source file would clutter it with volumes that
  have no `activity`.

### Source geometry — `data/source_geometry_lz_v1.json`

The reference for radioactivity. Lists every physical hardware piece
that hosts a measured ²¹⁴Bi / ²⁰⁸Tl activity, with **anatomical
dimensions** and the bb0nu paper masses.

Examples:

- OCV and ICV barrels and ellipsoidal heads (Ti)
- MLI insulation (mass-equivalent carrier on ICV outer wall)
- PTFE reflector walls, Ti field-cage rings, resistors, TPC sensors,
  top-grid and cathode-grid holders
- PMT arrays + bases + Ti structures + cables + skin PMTs

These objects host radioactivity and may self-shield source gammas
before those gammas reach the active LXe. The source-flux generator
samples (energy, position, direction) from each one with:

- **KN transport** (real Klein-Nishina, real μ_material) for sources
  flagged `approximation == "exact"` — used for the real-material Ti
  shells, PTFE, FC rings.
- **Transparent transport** (no source-side scattering or shielding)
  for `mass_equivalent_*` virtual carriers — used for sources where
  the physical hardware (PMT vacuum tubes, MLI blanket, sparse cables)
  cannot be drawn at full bulk density.

### Tracking geometry — `data/detector_lz_v3.json`

The reduced geometry used after the source-flux stage. Includes only
what the photon-transport needs:

- geometric acceptance to the fiducial volume
- hit formation in the sensitive regions
- veto logic in the skin

The Ti cryostat shells, MLI, PMT structures and PMT cables are
**absent** from the tracking tree because their effect has already
been absorbed into the source-flux stage. The tracker starts each
photon at the source's exit surface, not inside a Ti vessel.

### Algorithmic separation

The deeper reason for the split is that the two stages use different
propagation algorithms. A single mixed geometry would either
over-simplify the source stage (no self-shielding) or over-complicate
the tracking stage (irrelevant Ti cryostat in the photon transport).

**Source stage** propagation is specialized:

- photon-only KN-style transport inside the source bulk for `exact`
  shells, transparent transport for `mass_equivalent_*` virtuals
- no event stack — sample one (E, position, direction) per photon
  and apply real μ-attenuation through the source material
- termination at a virtual exit surface; the photon then enters the
  tracker with the post-shielding (E, position, direction)

**Tracking stage** is full event transport:

- multi-gamma event with optional Tl-208 cascade companion
- per-photon stack + Compton / pair-production cascade
- FV prefilter as the gate for event survival
- sensitive-region hit formation + veto deposit accumulation
- ROI smearing + SS/MS classification at the cluster level

### Volumes appearing in both files

`FC_PTFE`, `FC_rings`, and `Skin` appear in both source and tracking
because they both host radioactivity (PTFE walls, Ti rings, skin LXe)
*and* affect photon transport (the field-cage structure attenuates
photons traversing the FC annulus; the skin is the veto layer). When
a volume appears in both files it carries identical R, z, and wall
dimensions in both — anchored to the bb0nu / TDR measurement on the
source side and used as-is by the tracking partition.

| volume | R range (cm) | z range (cm) | t (cm) | mass anchor |
|---|---|---|---|---|
| `FC_PTFE` | [72.8, 73.947] | [−13.75, 145.6] | 1.147 | **184 kg** (bb0nu, geometric) |
| `FC_rings` | [73.947, 74.224] | [−13.75, 145.6] | 0.277 | **93 kg** (bb0nu, geometric) |
| `Skin` | [74.224, 82.1] | [−13.75, 145.6] | 7.876 | LXe (no activity — pure veto layer) |

The FC_PTFE outer wall touches FC_rings inner wall (no LXe sliver
between the two effective field-cage layers). Skin absorbs the
displacement: its inner radius shifted inward by 0.353 cm so its
outer wall stays at the envelope radius 82.1; its wall thickness
grew by the same 0.353 cm.

The strict source/tracking dimensional consistency is checked by
`@testset "Tracking/source consistency"` in `test/runtests.jl`.

## FastKernel design philosophy

Tracking is **not** a faithful Geant4-style geometry. It's a
**partition of space into named regions** where each region carries a
fast point-in-volume test (3–4 cheap arithmetic comparisons) and a
`tag` that drives policy (sensitive / not, ecut, veto). Photons
advance via "find next region boundary, jump" — no per-step recursion
through a hierarchy, no full event history maintained as a tree per
step. The whole design intent is to take fast classification decisions
cheaply and let the upper-level event logic (FV prefilter, SS/MS
clustering, veto) be where the physics complexity lives.

Two consequences flow from this:

1. **Volume tests are pure geometry.** `_is_inside_fastkernel_region`
   for `cylinder_shell` is literally `r²_min ≤ r² < r²_max ∧ z_min < z
   < z_max`. No coupling to other regions, no constants tied to
   `z_cathode = 0` or `z_RFR_bottom`. Whatever (R, z) bounds we put
   on a volume, the kernel tests them; nothing more.

2. **Tracking dimensions can match source dimensions.** Because the
   point-in-volume test doesn't care about anything outside its own
   bounds, anatomically correct dimensions (real wall thickness, real
   z extent) work just as well in tracking as simplified ones — and
   give better physics for free. We use the source-side anatomical
   dimensions for any volume that appears in both files.

The tracking file *can still be a subset* of the source file: Ti
cryostat shells, MLI, PMT structures and cables stay source-only
because the per-step photon transport doesn't need them once the
source-flux step has emitted gammas at the LXe boundary. That's a
different reason from "the tracking partition needs different
dimensions" — it's a "this volume isn't relevant to per-step
transport at all" reason.

## Mass anchoring (approximation values)

The `approximation` field on each volume controls how its absolute
activity is computed.

| approximation | activity scaling | typical use |
|---|---|---|
| `exact` | (geometric volume × ρ_material) × spec_act_mBq_per_kg | source-side real-material shells whose drawn dimensions reproduce the bb0nu mass within ~1% (OCV/ICV barrels, OCV/ICV heads, FC_PTFE, FC_rings) |
| `mass_equivalent_shell` | `equivalent_mass_kg` × spec_act; geometry sets only where photons emit | barrel-shaped carriers (MLI on ICV outer, PMT cables along the FC outer) |
| `mass_equivalent_slab` | same | thick annular slabs (top + cathode grid holders, R = [72.8, 80.3]) |
| `mass_equivalent_disk` | same | flat-disc carriers (PMT TOP/BOT PMTs, bases, structures, dome PMTs) |

For `mass_equivalent_*` volumes, `material` is `Vacuum` and
`transport_source` is `transparent`: the carrier is purely a
position-and-area handle, with no source-side self-shielding.

For `exact` volumes, `material` is real (Ti, PTFE, SS) and
`transport_source` is `KN`: emission samples real Klein-Nishina
kinematics inside the source bulk and applies real μ-attenuation to
the exit-angle distribution.

The geometric-mass invariant for `exact` source-side volumes is
checked by `@testset "Source geometric mass vs bb0nu"` in
`test/runtests.jl` (rtol = 5%).

## Live tracking tree (V3)

```
MARS                                  Vacuum    world (cylinder, R=120, H=440)
└── LZ_detector                       Vacuum    domed_container, R=82.1, ICV inner cavity
    ├── AirDome                       Vacuum    cap, gas dome (z_eq=148.5, ar=2.0)
    ├── AirCyl                        Vacuum    cylinder, gas gap (z=[145.6, 148.5])
    ├── LXe_dome                      LXe       cap, bottom dome (z∈[-68.7, -41.33])
    ├── LXe_below_FC                  LXe       cylinder, R≤82.1, z∈[-41.33, -13.75]
    ├── LXe_below_cathode             LXe       cylinder, R≤72.8, z∈[-13.75, 0]
    ├── TopActive                     LXe       cylinder, sensitive (above FV)
    ├── BarrelActive                  LXe       cylinder_shell, sensitive (around FV)
    ├── BottomActive                  LXe       cylinder, sensitive (below FV)
    ├── FV                            LXe       cylinder, sensitive + prefilter target
    ├── FC_PTFE                       PTFE      cylinder_shell, structural
    ├── FC_rings                      Ti        cylinder_shell, structural
    └── Skin                          LXe       cylinder_shell, sensitive (veto)
```

### ICV inner surface consistency

All tracking detector boundaries match the ICV inner surface at
R = 82.1 cm. The ICV head wall thicknesses were adjusted so that
R_inner = R_outer - t = 82.1 for all three surfaces (barrel, top,
bottom), while preserving the original dome masses. The real ICV
is not a perfect cylinder+dome; this is a deliberate simplification.

- ICV_barrel: R_outer = 83.0, t = 0.9 (unchanged)
- ICV_top: R_outer = 82.902, t = 0.802 (mass = 107.7 kg, preserved)
- ICV_bottom: R_outer = 83.292, t = 1.192 (mass = 141.4 kg, preserved)

The AirDome cap equator sits at z = 148.5 (ICV_top equator), not at
the gate plane z = 145.6. The cylindrical gas gap between the gate
and the dome equator is filled by AirCyl (R = 82.1, z = [145.6, 148.5]).

Containment invariants:

- `LZ_detector` ⊂ `MARS`.
- The twelve named daughters of `LZ_detector` form a strict partition
  of its interior: every point inside `LZ_detector` lies in exactly
  one daughter, no gaps, no overlaps. This is verified at runtime by
  the `FastKernel geometry partition` testset (sampling 10⁵ random
  points and asserting one match each).
- The three `passive_lxe` primitives (`LXe_dome`, `LXe_below_FC`,
  `LXe_below_cathode`) tile the LXe volume that the active and
  structural regions do not claim. They replaced the former
  `LXe_passive` capped-cylinder fallback, which had no well-defined
  boundary against its inner siblings and caused spurious interactions
  for gammas in transit through it.
- `TopActive`, `BarrelActive`, `BottomActive`, and `FV` together cover
  the active drift cylindrical region with no overlap (FV is the inner
  box at R≤39, z∈[26, 96]; the three "Active" volumes tile what's left
  of the active LXe outside the FV).
- `FC_PTFE`, `FC_rings`, and `Skin` form the radial partition for
  z ∈ [−13.75, 145.6] at R ∈ [72.8, 82.1]: PTFE inner at R=72.8 →
  rings inner at R=73.947 (touching) → Skin inner at R=74.224
  (touching) → envelope at R=82.1. No LXe gaps in the chain.

The Ti cryostat (OCV/ICV) is **not** in the tracking tree. It exists
only in `source_geometry_lz_v1.json`.

## V3 schema

Top-level fields: `name`, `version`, `units`, `world`, `volumes`.

Each entry in `volumes` carries:

| field | required | meaning |
|---|---|---|
| `name` | yes | unique within the file |
| `shape` | yes | one of `cylinder`, `cylinder_shell`, `disk`, `cap`, `domed_container`, `capped_cylinder`, `box` |
| shape dimensions | yes | per-shape — see *Solids* below |
| `position_cm` | yes | center of the shape in detector coords (cathode at z=0) |
| `orientation` | optional | `up`/`down` for `disk` and `cap` |
| `material` | yes | key into `materials.json` |
| `parent` | yes | name of the mother volume; `MARS` for world children |
| `tag` | yes | simulation policy (see *Tag set*) |
| `role` | yes | descriptive metadata, free-form |
| `approximation` | yes | `exact`, `mass_equivalent_shell`, `mass_equivalent_slab`, or `mass_equivalent_disk` |
| `ecut_keV`, `dz_mm` | optional | sensitive-region thresholds (only on `tpc_active`, `fv`, `skin`) |
| `_doc` | optional | one-line human comment |

The source-side schema additionally carries `source_class`,
`transport_source`, `equivalent_mass_kg`, and an `activity` block
(`Bi214_mBq_per_kg`, `Tl208_mBq_per_kg`) per volume.

## Solids

| solid | parameters | use |
|---|---|---|
| `Cyl` | `radius_cm`, `half_height_cm` | filled cylinder (FV, TopActive, PMT carriers) |
| `CylShell` | `R_inner_cm`, `wall_thickness_cm`, `half_height_cm` | shell (FC_PTFE, FC_rings, Skin, BarrelActive, PMT barrel carriers) |
| `Disk` | `R_cm`, `wall_thickness_cm`, `aspect_ratio`, `orientation` | flat or oblate-ellipsoidal head (cryostat heads) |
| `Cap` | `R_cm`, `aspect_ratio` | filled half-oblate cap (used inside `DomedContainer` / `CappedCylinder`) |
| `DomedContainer` | `radius_cm`, `barrel_half_height_cm`, top + bottom `Cap` | barrel + two coaxial caps (`LZ_detector`) |
| `CappedCylinder` | `radius_cm`, `barrel_half_height_cm`, optional top/bottom caps | barrel + 0–2 caps (no longer used by V3 — `LXe_passive` was replaced by three explicit primitives) |
| `Box` | `dx_cm`, `dy_cm`, `dz_cm` | rectangular block (not currently used in V3 but supported) |

`Cap` and `Disk` differ by purpose: `Disk` is a wall-shell with a
finite `wall_thickness_cm`; `Cap` is a filled solid used to extend a
cylinder into a dome.

## Tag set

| tag | semantics |
|---|---|
| `world` | the unique root region |
| `vacuum` | non-interacting medium (geometry navigation only) |
| `passive_lxe` | LXe but not sensitive — neither hits nor veto deposits |
| `tpc_active` | sensitive LXe in the active drift region; deposits become hits above `ecut_keV` |
| `fv` | sensitive AND the prefilter target — gates whether the event survives at all |
| `skin` | sensitive LXe in the veto layer — deposits feed skin-veto logic at a higher `ecut_keV` |
| `structural` | non-LXe material occupying space in the tracking tree (PTFE, Ti rings) — geometric obstruction only |

## Sensitivity vs material

Material alone does not define analysis behavior — the `tag` does.

- `tpc_active` LXe → sensitive, hits collected.
- `fv` LXe → sensitive AND target for the geometric prefilter.
- `skin` LXe → sensitive, deposits feed skin-veto logic.
- `passive_lxe` LXe → not sensitive (dome / RFR attenuate but don't hit).
- `structural` PTFE / Ti → never sensitive, purely a geometric obstacle.

This is why several LXe-filled volumes (`TopActive`, `BarrelActive`,
`BottomActive`, `FV`, `Skin`, and the three passive primitives
`LXe_dome` / `LXe_below_FC` / `LXe_below_cathode`) coexist as siblings
under `LZ_detector`: they encode different *analysis policies* on the
same material.

## Mother / daughter rule

`parent` defines registration into a mother volume.

Rule:

- The mother provides default material and default region semantics.
- Daughters displace the mother in the space they occupy.
- Navigation descends to the deepest containing daughter.
- If no daughter contains the point, the mother owns the point.

Required invariants:

- Every daughter is strictly contained in its declared mother
  (within float tolerance).
- Sibling daughters of the same mother do not overlap.

These are checked at load time. Strict mother-to-daughter descent is
the only navigation strategy used in production; there is no global
containing-node fallback.

## Validation

The V3 loader validates:

- unique names
- all `parent` references resolve
- exactly one root world (`MARS`)
- child containment inside parent
- sibling overlap
- required semantic fields present

The geometric-mass invariant for `exact` source-side volumes is
checked separately in the test suite.

## Runtime representation

V3 runtime model:

1. supported analytic solids (`Cyl`, `CylShell`, `Disk`, `Cap`,
   `DomedContainer`, `CappedCylinder`, `Box`)
2. logical volumes: solid + material + tag + role + approximation
3. placed nodes: logical volume + placement + parent + children
4. detector container: node set + root

The fast-kernel compiler (`compile_fastkernel_geometry`) flattens the
tree into a region array in JSON declaration order. Two distinguished
regions are identified by markers:
- **Envelope**: the unique region with `role == "tracking_envelope"`
- **Fallback**: the unique region with `tag == passive_lxe`

Classification loops over all regions; no volume name appears in the
dispatch logic.

## Adding a tracked region

Adding or removing a tracked region is a JSON-only change. No code edits
are required. The steps:

1. Add the volume entry to `data/detector_lz_v3.json` with:
   - `name`: any unique string (purely a label)
   - `shape` + dimensions: the geometric solid
   - `material`: key into `materials.json`
   - `parent`: the mother volume
   - `tag`: one of `world`, `vacuum`, `passive_lxe`, `tpc_active`,
     `fv`, `skin`, `structural`
   - `role`: a descriptive string (free-form)
   - Optional overrides: `ecut_keV`, `dz_mm`, `inFastKernel`, etc.

2. The `tag` drives all simulation behavior:
   - `_capabilities_for_tag(tag)` provides defaults for `inFastKernel`,
     `isXe`, `sensitive`, `ecut_keV`, `dz_mm`, `fv_target`
   - `select_interaction_fastkernel` classifies by tag, not by name
   - `veto_threshold` dispatches by tag
   - JSON fields override any default when present

3. The compiler validates structural invariants at load time:
   - Exactly one region with `role == "tracking_envelope"`
   - Exactly one region with `tag == passive_lxe`
   - Child containment and sibling non-overlap

To remove a region, delete its JSON entry. To rename one, change the
`name` field. Neither operation requires touching any `.jl` file.
