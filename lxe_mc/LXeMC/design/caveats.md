# Caveats and Operational Notes

Things a new Claude instance (or developer) should know before editing
the codebase. These are not bugs — they are deliberate design choices
that have non-obvious consequences.

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

1. Add volume entry to `data/detector_lz_v3.json` with the correct
   `tag` (`world`/`vacuum`/`structural`/`tpc_active`/`fv`/`skin`/
   `passive_lxe`). Transport helpers dispatch on this tag — no
   Julia edits needed.
2. Run the partition test (`FastKernel geometry partition` testset)
   and the tag-coverage tests (`_is_passive_region`, `_terminal_status`,
   `_classify_escape_volume`).
3. If the region is also a source, add it to `source_geometry_lz_v1.json`
   and update the source/tracking consistency test.
