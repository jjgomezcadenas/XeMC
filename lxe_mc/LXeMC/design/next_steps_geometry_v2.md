# Geometry V2 Next Steps

Branch: `newGeometry`

## Current committed state

Recent commits:

- `b83504a` `Add Geometry V2 validation and figure tooling`
- `89c3336` `Add Geometry V2 loader and hierarchy types`
- `7c34a55` `Add geometry redesign baseline and LZ schema v2`

The working tree was clean when this note was written.

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
- validation helpers
- overlap reporting
- tree dump

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

It also includes explicit containment/helper regions such as:

- `OCV_void`
- `ICV_LXe_interior`

### 3. New solid for filled lower head

Implemented only in the V2 path:

- `Cap`

This was introduced to represent filled LXe inside a head-like region rather than a shell.

The lower passive LXe region was split into:

- `DomeBarrel`
- `DomeBottomCap`

### 4. Validation and inspection

Implemented:

- child/parent consistency checks
- reachability from root
- containment sampling
- conservative sibling-overlap reporting
- exact-only overlap reporting

Official inspection script:

- [scripts/inspect_geometry_v2.jl](/Users/jjgomezcadenas/Projects/XeMC/lxe_mc/LXeMC/scripts/inspect_geometry_v2.jl)

This script can:

- print summary
- print tree
- print overlaps
- generate LaTeX/TikZ figures

### 5. Official figure source

Canonical figure source:

- [design/latex/lz_geometry_v2.tex](/Users/jjgomezcadenas/Projects/XeMC/lxe_mc/LXeMC/design/latex/lz_geometry_v2.tex)

Important design decision:

- this file is now the wrapper and the official source
- the old `lz_geometry_v2_png.tex` driver was removed
- LaTeX byproducts in `design/latex/` are ignored via repo `.gitignore`

### 6. Tests

Current V2 loader / validation coverage was added in:

- [test/runtests.jl](/Users/jjgomezcadenas/Projects/XeMC/lxe_mc/LXeMC/test/runtests.jl)

At the time of implementation, the Julia test suite passed with:

```bash
julia --project=. test/runtests.jl
```

Run this from:

- `/Users/jjgomezcadenas/Projects/XeMC/lxe_mc/LXeMC`

## Important current limitations

### 1. V2 is not yet used by transport

The production transport path is still based on the old geometry model in:

- [src/geometry.jl](/Users/jjgomezcadenas/Projects/XeMC/lxe_mc/LXeMC/src/geometry.jl)
- [src/tracking.jl](/Users/jjgomezcadenas/Projects/XeMC/lxe_mc/LXeMC/src/tracking.jl)

In particular:

- `propagate_to_lxe`
- `find_volume`
- `classify_lxe_region`

are still old-path logic.

### 2. Semantics are not yet fully migrated

The V2 JSON has `tag`/`role`, but the runtime transport/veto path has not yet been switched to semantic lookup from V2 nodes.

### 3. Navigation is not yet implemented for V2

What exists now:

- representation
- loading
- validation
- inspection

What does not yet exist:

- V2 point location for transport
- V2 boundary stepping
- V2 neighbor-aware traversal

### 4. Geometry figure is a support artifact, not physics validation

The figure is useful for human inspection only. It does not validate transport correctness.

## Recommended next implementation steps

The next work should move from representation to runtime use, but in a staged way.

### Step 1. Add semantic helpers for V2

Goal:

- stop relying on volume names in the new path

Add helpers in `geometry2.jl` or a small companion file if it becomes cleaner:

- `region_tag(node)`
- `is_fv(node)`
- `is_active_lxe(node)`
- `is_veto_lxe(node)`
- `is_passive_lxe(node)`
- `is_structural(node)`
- `is_vacuum(node)`

Also add a helper for policy mapping:

- `veto_threshold(tag, cfg)`

This should mirror the intent currently embedded in `classify_lxe_region` and the veto logic in `tracking.jl`, but without any string-name inspection.

Do not remove the old path yet.

### Step 2. Add V2 point-location primitives

Goal:

- determine which V2 node contains a point

Add a new API first, without touching old transport:

- `find_node_v2(det, pos)`
- maybe `find_deepest_node_v2(det, pos; start=root)`

Expected behavior:

- descend the V2 tree from the root
- at each level, check only candidate children
- return the deepest containing node

Initially, correctness matters more than speed.

Once working, add a second version later that accepts a current node hint for stateful navigation.

### Step 3. Add classification comparison tests

Goal:

- compare old geometry semantics with the new V2 semantics

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

This step is critical before trying to use V2 in transport.

### Step 4. Add V2 boundary / navigation skeleton

Goal:

- prepare V2 for actual stepping across materials/regions

Add basic helpers:

- distance to exit current node
- candidate child entry tests
- parent fallback on exit

Do not try to solve the full optimized navigator at once.

A reasonable first version is:

- correct
- tree-aware
- exact for the supported primitives

Speed optimization can come after functional equivalence.

### Step 5. Add `propagate_to_lxe_v2`

Goal:

- run the same high-level logic as `propagate_to_lxe`, but backed by V2 geometry

Do this as a parallel function first:

- do not replace `propagate_to_lxe` immediately

Expected behavior:

- identify current material/region from V2 nodes
- apply veto semantics from `tag`
- classify FV acceptance from `tag == TAG_FV` or equivalent
- preserve the external behavior and return style of the existing function as much as possible

### Step 6. Compare old and new behavior

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
- Do not use volume-name string parsing in new V2 logic except for temporary test comparisons

## Practical resume notes

When resuming work:

1. Start from branch `newGeometry`
2. Open:
   - [src/geometry2.jl](/Users/jjgomezcadenas/Projects/XeMC/lxe_mc/LXeMC/src/geometry2.jl)
   - [src/tracking.jl](/Users/jjgomezcadenas/Projects/XeMC/lxe_mc/LXeMC/src/tracking.jl)
   - [data/detector_lz_v2.json](/Users/jjgomezcadenas/Projects/XeMC/lxe_mc/data/detector_lz_v2.json)
   - [scripts/inspect_geometry_v2.jl](/Users/jjgomezcadenas/Projects/XeMC/lxe_mc/LXeMC/scripts/inspect_geometry_v2.jl)
3. First target:
   - semantic/tag helpers
   - V2 point location
4. Then add tests before transport integration

## Suggested next commit sequence

1. `Add Geometry V2 semantic helpers`
2. `Add Geometry V2 point location`
3. `Add Geometry V2 classification comparison tests`
4. `Add experimental propagate_to_lxe_v2`

This sequence is preferred over one large integration commit.
