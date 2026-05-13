# Fix: replace name-based dispatch with tag-based dispatch in transport helpers

## Problem

The FastKernel geometry compiler and classifier are fully data-driven:
no volume names in `compile_fastkernel_geometry` or `classify_fastkernel`.
Adding a new region is supposed to be a JSON-only change.

But the transport-level helpers defeat this by hardcoding region names.
Adding a new `tpc_active` region (or a new `passive_lxe` region) requires
editing several Julia functions in addition to the JSON — and forgetting
to do so is a silent bug (the new region would be treated as neither
active nor passive, falling through to wrong behavior).

## Affected functions

### `tracking_fast.jl::_visible_threshold_MeV(cfg, region_name)`

Current:
```julia
region_name == "Skin" && return cfg.veto_skin
region_name in ("TopActive", "BarrelActive", "BottomActive") && return cfg.veto_TPC
region_name == "FV" && return cfg.veto_TPC
return Inf
```

### `tracking_fast.jl::_is_passive_region(region_name)`

Current:
```julia
region_name == "LXe_dome" || region_name == "LXe_below_FC" ||
    region_name == "LXe_below_cathode" ||
    region_name == "FC_PTFE" || region_name == "FC_rings"
```

### `tracking.jl::_classify_escape_volume(fk, pos)`

Current:
```julia
rname == "FV" && return :fv
rname == "Skin" && return :active
rname in ("TopActive", "BarrelActive", "BottomActive") && return :active
return :passive
```

### `tracking_fast.jl::transport_gamma_fastkernel` (main loop)

Current:
```julia
region.name == "FV"  # handoff check
rname in ("TopActive", "BarrelActive", "BottomActive")  # veto check
```

### `tracking_fast.jl::_terminal_status(deposits, rname, pos, dir)`

Current:
```julia
if rname == "Skin"
    _fast_track_result(:vetoed_skin, deposits, rname, 0.0, pos, dir)
elseif rname in ("TopActive", "BarrelActive", "BottomActive")
    _fast_track_result(:vetoed_tpc, deposits, rname, 0.0, pos, dir)
elseif _is_passive_region(rname)
    _fast_track_result(:absorbed_passive, deposits, rname, 0.0, pos, dir)
else
    _fast_track_result(:below_cut, deposits, rname, 0.0, pos, dir)
end
```

## Proposed fix

Replace all name-based checks with `RegionTag` enum comparisons. The
`FastKernelRegion` already carries `tag::RegionTag`, so no struct changes
are needed.

### `_visible_threshold_MeV(cfg, tag::RegionTag)`

```julia
@inline function _visible_threshold_MeV(cfg::SimConfig, tag::RegionTag)::Float64
    tag == TAG_SKIN && return cfg.veto_skin
    tag == TAG_TPC_ACTIVE && return cfg.veto_TPC
    tag == TAG_FV && return cfg.veto_TPC
    return Inf
end
```

### `_is_passive_region(tag::RegionTag)`

```julia
@inline function _is_passive_region(tag::RegionTag)::Bool
    tag == TAG_PASSIVE_LXE || tag == TAG_STRUCTURAL
end
```

### `_classify_escape_volume(fk, pos)` — use `region.tag`

```julia
@inline function _classify_escape_volume(fk::FastKernelGeometry, pos::Vector{Float64})::Symbol
    region = classify_fastkernel(fk, (pos[1], pos[2], pos[3]))
    region === nothing && return :passive
    tag = region.tag
    tag == TAG_FV && return :fv
    tag == TAG_SKIN && return :active
    tag == TAG_TPC_ACTIVE && return :active
    return :passive
end
```

### `transport_gamma_fastkernel` — use `region.tag`

Replace:
```julia
region.name == "FV"
```
with:
```julia
region.tag == TAG_FV
```

Replace:
```julia
rname in ("TopActive", "BarrelActive", "BottomActive") && dep >= cfg.veto_TPC
```
with:
```julia
region.tag == TAG_TPC_ACTIVE && dep >= cfg.veto_TPC
```

And pass `region.tag` instead of `region.name` to `_visible_threshold_MeV`
and `_is_passive_region`.

### `_terminal_status` — take `tag::RegionTag` (keep `rname` only for diagnostics)

```julia
@inline function _terminal_status(deposits::Vector{FastGammaDeposit},
                                   tag::RegionTag,
                                   rname::String,
                                   pos::Vector{Float64},
                                   dir::Vector{Float64})
    status = if tag == TAG_SKIN
        :vetoed_skin
    elseif tag == TAG_TPC_ACTIVE
        :vetoed_tpc
    elseif _is_passive_region(tag)
        :absorbed_passive
    else
        :below_cut
    end
    _fast_track_result(status, deposits, rname, 0.0, pos, dir)
end
```

Callers in `transport_gamma_fastkernel` pass `region.tag` alongside
`region.name` (the latter retained only as the `terminal_region`
diagnostic string in the returned `FastGammaTrackResult`).

## Relation to existing `veto_threshold`

`geometry.jl:156` already defines:

```julia
function veto_threshold(tag::RegionTag, cfg::SimConfig)::Float64
    tag == TAG_TPC_ACTIVE && return cfg.veto_TPC
    tag == TAG_SKIN && return cfg.veto_skin
    tag == TAG_FV && return 0.0
    tag == TAG_PASSIVE_LXE && return Inf
    ...
```

This is nearly identical to the proposed `_visible_threshold_MeV`. The
only semantic difference is `TAG_FV`:

| Tag | `veto_threshold` | proposed `_visible_threshold_MeV` |
|:----|:-----------------|:---------------------------------|
| `TAG_FV` | `0.0` | `cfg.veto_TPC` |

Before adding `_visible_threshold_MeV` as a new function, decide:

1. **Reuse**: replace `_visible_threshold_MeV` calls with `veto_threshold`
   and accept that `TAG_FV` returns `0.0` (any deposit triggers the
   threshold). This works if the fast kernel never actually calls the
   threshold function inside the FV (it shouldn't — FV is a handoff
   point, not a deposition region in the fast kernel).
2. **Diverge intentionally**: keep both functions, and document why the
   FV semantics differ. The current name-based code returns `cfg.veto_TPC`
   for `"FV"`, suggesting that the FV is treated as TPC-active for the
   purpose of any fallback threshold check — but this branch may be
   unreachable in practice.

Recommend option 1 (reuse `veto_threshold`) after verifying that no
fast-kernel code path calls the threshold function with `region.tag ==
TAG_FV`. If a code path does exist, that's likely a bug worth fixing
rather than preserving.

## What changes

| Before | After |
|:-------|:------|
| `_visible_threshold_MeV(cfg, region_name::String)` | `_visible_threshold_MeV(cfg, tag::RegionTag)` |
| `_is_passive_region(region_name::String)` | `_is_passive_region(tag::RegionTag)` |
| `_classify_escape_volume` checks `rname` strings | checks `region.tag` enum |
| `transport_gamma_fastkernel` checks `region.name` | checks `region.tag` |
| `_terminal_status` takes `rname::String` and dispatches on it | takes `tag::RegionTag` and dispatches on it (`rname` kept only for diagnostics) |

## What does NOT change

- `FastKernelRegion` struct (already has `tag::RegionTag`)
- `classify_fastkernel` (already returns `FastKernelRegion`)
- `FastGammaDeposit` (the `region` field can stay as the region name
  string for diagnostic purposes, or switch to tag — either works)
- JSON geometry files
- Test suite (tests check behavior, not internal dispatch)

## Benefits

1. Adding a new `tpc_active` region is truly JSON-only — no code edits
2. Adding a new `passive_lxe` region is truly JSON-only
3. No string comparisons in the hot path (integer enum instead)
4. The documented claim "no volume name in dispatch logic" becomes true
   for the entire codebase, not just the compiler/classifier

## Risk

Low. The tag enum already encodes exactly the semantic distinctions these
helpers need. The mapping is 1:1:

| Name-based check | Tag equivalent |
|:-----------------|:---------------|
| `"Skin"` | `TAG_SKIN` |
| `"TopActive"`, `"BarrelActive"`, `"BottomActive"` | `TAG_TPC_ACTIVE` |
| `"FV"` | `TAG_FV` |
| `"LXe_dome"`, `"LXe_below_FC"`, `"LXe_below_cathode"` | `TAG_PASSIVE_LXE` |
| `"FC_PTFE"`, `"FC_rings"` | `TAG_STRUCTURAL` |

## Test plan

1. Run existing test suite (`Pkg.test()`) — all tests should pass unchanged
2. Run `check_gamma_prop_simple.jl` and `test_event_processing.jl` —
   same output as before
3. Verify `FastKernel geometry partition` testset still passes
4. Optionally: add a test that confirms `_visible_threshold_MeV` and
   `_is_passive_region` cover all `RegionTag` values
