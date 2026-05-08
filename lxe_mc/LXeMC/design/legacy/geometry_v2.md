# Geometry V2

## Scope

Geometry V2 introduces:

- a new detector JSON schema in `data/detector_lz_v2.json`
- explicit mother/daughter topology
- explicit simulation semantics via tags
- explicit modeling status via approximation labels
- a separation between source geometry and tracking geometry

The current geometry path remains in place during the migration.

## Current model

Current code in `src/geometry.jl` uses:

- solids: `Cyl`, `CylShell`, `Box`, `Disk`
- logical volumes: `LCyl`, `LCylShell`, `LBox`, `LDisk`
- physical volumes: `PCyl`, `PCylShell`, `PBox`, `PDisk`
- detector: `mars`, flat `volumes`, separate `fiducial`

Current limitations:

- flat volume list
- special-case fiducial volume
- region classification from volume names
- containment not represented explicitly

## V2 schema

Top-level fields:

- `name`
- `version`
- `units`
- `world`
- `volumes`

`world` is the unique root region.

Each entry in `volumes` contains:

- `name`
- `shape`
- shape dimensions
- `position_cm`
- optional `orientation`
- `material`
- `parent`
- `tag`
- `role`
- `approximation`
- optional `activity`
- optional `_doc`

## Mother / daughter rule

`parent` defines registration into a mother volume.

The intended rule is:

- the mother provides default material and default region semantics
- daughters displace the mother in the space they occupy
- navigation descends to the deepest containing daughter
- if no daughter contains the point, the mother owns the point

This requires one geometric constraint:

- every daughter must be strictly contained in its declared mother, within tolerance

This is the key invariant of the V2 tree.

## Source geometry vs tracking geometry

This distinction is essential.

### Source geometry

Source geometry is used to model radioactive components and compute source fluxes.

Examples:

- Ti cryostat barrel and heads
- flanges
- MLI
- later PMT source objects

These objects are relevant because:

- they host radioactivity
- they may self-shield source gammas before those gammas enter the detector cavity

### Tracking geometry

Tracking geometry is the reduced geometry used after the source-flux stage.

It should include only what is needed for:

- geometric acceptance to `FV`
- hit formation in sensitive regions
- veto logic

Therefore, source objects such as Ti cryostat shells do not necessarily belong in the tracking tree if their effect has already been absorbed into the source-flux stage.

This is now the intended philosophy for V2.

### Consequence

The runtime tracking tree should be built from detector cavities and internal detector regions, not automatically from every radioactive shell object.

For current purposes, the desired high-level tracking structure is closer to:

- `LZ_det`
  - `OCV_void`
    - `ICV_void`
      - `LXe_detector`
        - `LXeTPC`
        - `Skin`
        - `RFR`
        - `Dome`
        - field-cage structures
        - later PMTs and internal objects as needed

The Ti shells remain source geometry unless a later transport use case requires them.

## Materials and containers

Container volumes are allowed.

Examples:

- `OCV`
- `ICV`
- `ICV_void`
- `LXe_detector`

These may be filled with vacuum or LXe and then have daughters registered inside them.

Examples:

- `ICV_void` may be vacuum
- `LXe_detector` registered inside `ICV_void` displaces vacuum with LXe
- PTFE, PMTs, grids, etc. registered inside `LXe_detector` displace LXe

So a mother does not need to remain homogeneous after placement. It only defines the default material of space not occupied by daughters.

## Tags and roles

### `tag`

`tag` defines simulation policy.

Current tag set:

- `world`
- `vacuum`
- `structural`
- `tpc_active`
- `fv`
- `skin`
- `passive_lxe`

### `role`

`role` is descriptive metadata.

Examples:

- `cryostat_cavity`
- `fiducial_region`
- `field_cage`
- `top_grid_holder`
- `support_flange`

### `approximation`

`approximation` records whether the object is exact or a proxy.

Current values:

- `exact`
- `cylindrical_proxy`
- `mass_equivalent_shell`
- `mass_equivalent_slab`
- `negligible_thickness_proxy`

## Sensitive regions

Material alone does not define analysis behavior.

Important distinction:

- sensitive LXe regions
- passive LXe regions

Current intended semantics:

- `LXeTPC`: sensitive
- `Skin`: sensitive
- passive LXe (`RFR`, `Dome`, etc.): not sensitive
- `FV`: sensitive and also the special fast-prefilter target region

`FV` is not just another sensitive region. It is the event-level target used by the geometric prefilter.

## Helper regions and fake volumes

Volumes such as:

- `LXeTPC`
- `FV`
- `Skin`
- `RFR`
- `DomeBarrel`
- `DomeBottomCap`

may be retained even when they are all LXe-filled, because they encode:

- sensitivity
- hit policy
- veto policy
- FV-target semantics

These are analysis/transport control regions, not merely material partitions.

## Region nodes

V2 allows explicit containment regions that are not part objects.

Examples:

- `OCV_void`
- `ICV_void`
- `LXe_detector`

These exist to encode space and topology, not necessarily to represent separate hardware pieces.

## Restricted container solids

The current preferred extension is not a generic union, but a restricted domed container solid.

Examples:

- `LZ_det = barrel + top cap + bottom cap`
- `OCV_void = barrel cavity + top cavity cap + bottom cavity cap`
- `ICV_void = barrel cavity + top cavity cap + bottom cavity cap`
- `LXe_detector = barrel + cap(s)` if needed

This is preferred over proxy cylindrical mothers when the daughters extend into head regions.

The intended scope is narrow:

- one cylinder plus one top cap plus one bottom cap
- coaxial only
- not general boolean geometry

## Runtime representation

V2 runtime model:

1. supported analytic solids
2. logical volumes: solid + material + tag
3. placed nodes: logical volume + placement + parent + children
4. detector container: node set + root

The current detector representation remains available during migration.

## Validation

The V2 loader should validate:

- unique names
- all `parent` references resolve
- exactly one root world
- child containment inside parent
- obvious sibling overlaps
- required semantic fields present

Validation should be local:

- child inside mother
- sibling overlap

Global search should only be a debug fallback, not the design.

## Current caveat

The present V2 file still contains some helper/proxy mothers that are not yet true geometric containers for all children.

Because of that, the current `find_node_v2` implementation includes a temporary global containing-node selection.

This is not the target design.

Target design:

- repair the hierarchy so strict mother-to-daughter descent is valid
- then remove the global fallback

## Migration

1. Add V2 structs.
2. Implement `load_detector_v2`.
3. Add validation and tree inspection utilities.
4. Add semantic helpers and point location.
5. Separate source geometry from tracking geometry.
6. Repair mother/daughter containment in the tracking tree so strict descent is valid.
7. Add V2 navigation.
8. Migrate `propagate_to_lxe` to V2 after comparison with the current path.
