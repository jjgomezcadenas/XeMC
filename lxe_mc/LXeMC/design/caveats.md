# Caveats and Operational Notes

Things a new Claude instance (or developer) should know before editing
the codebase. These are not bugs — they are deliberate design choices
that have non-obvious consequences.

## Name-based dispatch in transport helpers

The FastKernel geometry compiler (`compile_fastkernel_geometry`) and
classifier (`classify_fastkernel`) are fully data-driven — no volume
names in their logic. However, the **transport-level helpers** do
hardcode region names:

- `tracking_fast.jl::_visible_threshold_MeV` — checks `"TopActive"`,
  `"BarrelActive"`, `"BottomActive"` for TPC veto and `"Skin"` for
  skin veto.
- `tracking_fast.jl::_is_passive_region` — checks `"LXe_dome"`,
  `"LXe_below_FC"`, `"LXe_below_cathode"`, `"FC_PTFE"`, `"FC_rings"`.
- `tracking.jl::_classify_escape_volume` — classifies FV-escape
  positions as `:active` (Skin, TopActive, BarrelActive, BottomActive)
  or `:passive` (the three passive LXe + FC volumes).

**Consequence**: adding a new region to `detector_lz_v3.json` is NOT
a pure JSON change if the region should participate in veto or
escape-volume classification. You must also update these three helpers.

## ParticleStack is a struct, not a type alias

`ParticleStack` wraps `Vector{Track}` as a LIFO stack struct with
`items::Vector{Track}`. It is **not** `Vector{Track}` directly.

## Full-stack functions take `fk` parameter

`propagate_gamma`, `transport_photon!`, and `transport_lepton!` all
take a `FastKernelGeometry` parameter (`fk`) in addition to the
`PhysicalVolume` (`vol`). The `fk` is used by
`_classify_escape_volume` to tag deposits that escape the FV boundary.

Signatures:
- `propagate_gamma(E_MeV, vol, fk, cfg; position, direction, rng)`
- `transport_photon!(track, vol, fk, deposits, stack, track_counter, cfg, rng)`
- `transport_lepton!(track, vol, fk, deposits, stack, track_counter, cfg, rng)`

## Adding a new tracked region — checklist

1. Add volume entry to `data/detector_lz_v3.json` (tag, role, shape)
2. Update `_visible_threshold_MeV` in `tracking_fast.jl` if it needs a veto
3. Update `_is_passive_region` in `tracking_fast.jl` if it's passive
4. Update `_classify_escape_volume` in `tracking.jl` if it's active or passive
5. Run the partition test (`FastKernel geometry partition` testset)
6. If the region is also a source, add it to `source_geometry_lz_v1.json`
   and update the source/tracking consistency test
