# Geometry V2 Next Steps

Branch: `newGeometry`

## Current committed state

Recent commits:

- `2a7bc2f` `Add DomedContainer solid for Geometry V2`
- `6f789af` `Add Geometry V2 semantic helpers and point location`
- `b83504a` `Add Geometry V2 validation and figure tooling`
- `89c3336` `Add Geometry V2 loader and hierarchy types`
- `7c34a55` `Add geometry redesign baseline and LZ schema v2`

The working tree was clean when this note was updated.

## What is already in place

### 1. Geometry V2 data model

Implemented in [src/geometry2.jl](/Users/jjgomezcadenas/Projects/XeMC/lxe_mc/LXeMC/src/geometry2.jl):

- `RegionTag`
- `PlacementV2`
- `LogicalVolumeV2`
- `DetectorNode`
- `DetectorV2`
- `load_detector_v2`
- `node_by_name`, `root_node`, `child_nodes`, `detector_summary`
- semantic helpers
- `find_node_v2`
- validation helpers
- overlap reporting
- tree dump
- `DomedContainer`

The V2 path is loaded in [src/LXeMC.jl](/Users/jjgomezcadenas/Projects/XeMC/lxe_mc/LXeMC/src/LXeMC.jl).

### 2. Geometry V2 schema

Main file:

- [data/detector_lz_v2.json](/Users/jjgomezcadenas/Projects/XeMC/lxe_mc/data/detector_lz_v2.json)

This file already contains:

- `world`
- flat `volumes`
- `parent`
- `tag`
- `role`
- `approximation`

It also includes helper/container-style regions such as:

- `OCV_void`
- `ICV_LXe_interior`

### 3. New container solid

Implemented only in the V2 path:

- `DomedContainer`

This is the preferred restricted mother/container solid for:

- one cylinder
- one top cap
- one bottom cap
- coaxial layouts

### 4. New solid for filled lower head

Implemented only in the V2 path:

- `Cap`

The lower passive LXe region was split into:

- `DomeBarrel`
- `DomeBottomCap`

### 5. Validation and inspection

Implemented:

- child/parent consistency checks
- reachability from root
- containment sampling
- conservative sibling-overlap reporting
- exact-only overlap reporting

Official inspection script:

- [scripts/geometry_viewer.html](/Users/jjgomezcadenas/Projects/XeMC/lxe_mc/LXeMC/scripts/geometry_viewer.html)

This viewer can:

- show tracking and source geometry side by side
- inspect metadata by clicking on a volume
- switch the right-side legend between tracking and source color maps

### 6. Official geometry inspection

Canonical inspection path:

- [scripts/geometry_viewer.html](/Users/jjgomezcadenas/Projects/XeMC/lxe_mc/LXeMC/scripts/geometry_viewer.html)

### 7. Tests

Current V2 coverage in:

- [test/runtests.jl](/Users/jjgomezcadenas/Projects/XeMC/lxe_mc/LXeMC/test/runtests.jl)

The Julia test suite passed after the current V2 semantic / point-location changes with:

```bash
julia --project=. test/runtests.jl
```

Run this from:

- `/Users/jjgomezcadenas/Projects/XeMC/lxe_mc/LXeMC`

## Important current limitations

### 1. Source geometry vs tracking geometry is not yet reflected in the JSON

We now want V2 to distinguish:

- source geometry
- tracking geometry

This has not yet been encoded cleanly in `detector_lz_v2.json`.

In particular, Ti cryostat shells are source objects first, and may not belong in the reduced tracking tree at all.

### 2. V2 is not yet used by transport

The production transport path is still based on the old geometry model in:

- [src/geometry.jl](/Users/jjgomezcadenas/Projects/XeMC/lxe_mc/LXeMC/src/geometry.jl)
- [src/tracking.jl](/Users/jjgomezcadenas/Projects/XeMC/lxe_mc/LXeMC/src/tracking.jl)

In particular:

- `propagate_to_lxe`
- `find_volume`
- `classify_lxe_region`

are still old-path logic.

### 3. Current point location is temporary

`find_node_v2` currently uses a global containing-node selection:

- collect all containing nodes
- prefer deepest
- then prefer exact over proxy at equal depth
- then prefer the smaller region

This was introduced as a correctness fallback because the current helper/mother hierarchy is not yet a strict geometric containment tree.

This is not the intended final design.

### 4. Mother/daughter containment is still inconsistent

The present V2 file still mixes source-style objects and tracking-style mothers, and also contains helper/proxy mothers that are not true geometric containers for all children.

This must be repaired before implementing final navigation.

### 5. Semantics are only partially migrated

The V2 path now has semantic helpers, but the runtime transport/veto path has not yet been switched to use them.

### 6. Geometry figure is a support artifact, not physics validation

The figure is useful for human inspection only. It does not validate transport correctness.

## Recommended next implementation steps

The next work should correct the geometry model before adding more transport logic.

### Step 1. Separate tracking geometry from source geometry

Goal:

- define the reduced runtime tracking tree

Likely actions:

- introduce top-level tracking containers such as:
  - `LZ_det`
  - `OCV_void`
  - `ICV_void`
  - `LXe_detector`
- keep Ti shells and other radioactive shells conceptually in source geometry rather than forcing them into the runtime tracking tree

This is now the priority issue.

### Step 2. Repair mother/daughter containment in the tracking tree

Goal:

- make declared mothers true geometric containers

Likely actions:

- replace proxy cylindrical mothers where they fail containment
- use `DomedContainer` for barrel + head container volumes
- reparent LXe and internal detector objects under the reduced tracking-tree mothers

### Step 3. Restore strict descent-based point location

Goal:

- remove the temporary global `find_node_v2` fallback

Target behavior:

- descend the containment tree from the root
- only test candidate daughters of the current mother
- return the deepest containing node

### Step 4. Add classification comparison tests

Goal:

- compare old geometry semantics with the repaired V2 semantics

Create representative points in:

- world / outside detector
- OCV metal
- OCV void
- ICV metal
- LXeTPC
- FV
- Skin
- RFR
- DomeBarrel
- DomeBottomCap

For each point:

- record old-path classification where applicable
- record V2 node + tag
- confirm the new semantic interpretation is what is intended

### Step 5. Add V2 boundary / navigation skeleton

Goal:

- prepare V2 for actual stepping across materials/regions

Add basic helpers:

- distance to exit current node
- candidate child entry tests
- parent fallback on exit

Do not try to solve the full optimized navigator at once.

### Step 6. Add event-level FV prefilter logic

Goal:

- separate geometric rejection from later veto/hit decisions

Event-level rule:

- `0` gammas: reject
- `1` gamma: keep only if it reaches `FV`
- `2+` gammas: keep only if any gamma reaches `FV`

This should happen before full event stack development.

### Step 7. Add `propagate_to_lxe_v2`

Goal:

- run the same high-level logic as `propagate_to_lxe`, but backed by V2 geometry
- keep veto / hit logic after the geometric prefilter

Do this as a parallel function first.

### Step 8. Compare old and new behavior

Before retiring old logic, compare:

- accepted vs vetoed outcomes
- region classifications
- traversed materials
- representative ray tests

Only once this is stable should the old path be considered for replacement.

## What should not be done yet

- Do not delete or rewrite `src/geometry.jl`
- Do not replace `propagate_to_lxe` in place
- Do not merge V2 and old geometry into one file yet
- Do not optimize the navigator before there is a correct baseline
- Do not normalize the temporary global `find_node_v2` fallback into the final design

## Practical resume notes

When resuming work:

1. Start from branch `newGeometry`
2. Open:
   - [src/geometry2.jl](/Users/jjgomezcadenas/Projects/XeMC/lxe_mc/LXeMC/src/geometry2.jl)
   - [src/tracking.jl](/Users/jjgomezcadenas/Projects/XeMC/lxe_mc/LXeMC/src/tracking.jl)
   - [data/detector_lz_v2.json](/Users/jjgomezcadenas/Projects/XeMC/lxe_mc/data/detector_lz_v2.json)
   - [scripts/geometry_viewer.html](/Users/jjgomezcadenas/Projects/XeMC/lxe_mc/LXeMC/scripts/geometry_viewer.html)
   - [design/geometry_v2_prefilter.md](/Users/jjgomezcadenas/Projects/XeMC/lxe_mc/LXeMC/design/geometry_v2_prefilter.md)
   - [design/source_vs_tracking_geometry.md](/Users/jjgomezcadenas/Projects/XeMC/lxe_mc/LXeMC/design/source_vs_tracking_geometry.md)
3. First target:
   - separate tracking geometry from source geometry
   - repair mother/daughter containment in the tracking tree
   - then restore strict descent
4. Then add/adjust tests before transport integration

## Suggested next commit sequence

1. `Separate Geometry V2 tracking tree from source geometry`
2. `Repair Geometry V2 containment hierarchy`
3. `Restore strict Geometry V2 point location`
4. `Add Geometry V2 classification comparison tests`
5. `Add Geometry V2 FV prefilter`
6. `Add experimental propagate_to_lxe_v2`

This sequence is preferred over one large integration commit.
