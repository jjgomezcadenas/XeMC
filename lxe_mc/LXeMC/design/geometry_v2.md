# Geometry V2

## Scope

Geometry V2 introduces:

- a new detector JSON schema in `data/detector_lz_v2.json`
- explicit parent/child topology
- explicit simulation semantics via tags
- explicit modeling status via approximation labels

The current geometry path remains in place during the migration.


## Current Model

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


## V2 Schema

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


## Required Semantics

### `parent`

`parent` defines containment.

Examples:

- `FV.parent = LXeTPC`
- `LXeTPC.parent = ICV_LXe_interior`
- `ICV_barrel.parent = OCV_void`
- `OCV_void.parent = MARS`

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

- `fiducial_region`
- `field_cage`
- `cryostat_barrel`
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


## Region Nodes

V2 allows explicit containment regions that are not part objects.

Current examples:

- `OCV_void`
- `ICV_LXe_interior`

These exist to encode space and topology.


## Runtime Representation

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


## Migration

1. Add V2 structs.
2. Implement `load_detector_v2`.
3. Add validation and tree inspection utilities.
4. Use tags instead of name-based region classification in the V2 path.
5. Add V2 navigation.
6. Migrate `propagate_to_lxe` to V2 after comparison with the current path.
