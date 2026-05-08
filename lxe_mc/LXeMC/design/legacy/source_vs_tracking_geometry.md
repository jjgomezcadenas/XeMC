# Source Geometry vs Tracking Geometry

## Purpose

This note records the current design decision that source modelling and tracking geometry are different problems and should be represented separately.

This is not a minor implementation detail. It determines:

- what geometry belongs in the tracking tree
- what geometry belongs in source modelling
- where navigation starts
- what kinds of propagation are performed in each stage

## Core decision

We do **not** want one single geometry model that mixes:

- radioactive source objects
- detector cavities
- sensitive LXe regions
- veto regions
- FV target regions

Instead, we want two separate geometry uses:

1. **source geometry**
2. **tracking geometry**

These may share dimensions and names where useful, but they serve different algorithms and different purposes.

## Source geometry

Source geometry is used to generate the gamma flux from radioactive components into the detector.

Examples of source families:

- cryostat titanium
  - outer cryostat barrel + heads
  - inner cryostat barrel + heads
  - flanges and related Ti structures
- PMTs
  - top PMTs
  - bottom PMTs
  - skin/barrel PMTs
- insulation / MLI
- any future surface or bulk contamination source classes

### What source geometry is used for

For a source such as the cryostat, a gamma may:

- start inside a Ti shell
- die locally
- go backward and be lost
- exit the shell
- cross vacuum
- possibly cross another shell
- eventually reach a virtual detector boundary

This is not the same problem as runtime event tracking in the detector.

Source geometry is used to compute a flux:

- `dPhi/dNdu` per event
- or more generally a flux table over energy and direction
- evaluated at the entry to the detector tracking stage

The relevant endpoint is not “full detector response”.
It is:

- the gamma flux entering the detector tracking mother volume

### Source propagation

Source propagation is specialized.

It is not full event transport.

Current intended/accepted simplifications:

- source-specific propagation
- no full stack development
- photon-only / KN-style propagation where appropriate
- termination at a virtual LXe or detector entry surface

Therefore, source objects are not merely a tagged subset of the runtime tracking geometry.

They belong to a separate source model.

## Tracking geometry

Tracking geometry begins only after the source-flux stage.

Its input is:

- one or more gammas already sampled at the detector entry
- with energy, direction, and normalization inherited from source modelling

Tracking geometry is then responsible for:

- FV reachability
- geometric prefilter logic
- sensitive-region hit formation
- veto logic
- later full-event tracking for surviving events

### Consequence

Tracking does **not** need to know anything about geometry outside the detector entry mother volume.

Once a source family has produced a flux into the detector, the external cryostat/source world is gone from the runtime problem.

This is the main simplification.

## The true tracking mother

The correct top-level mother for runtime propagation is not:

- OCV shell
- OCV cavity
- ICV shell
- ICV cavity

provided those are only needed to build source fluxes.

The true runtime mother is:

- `LZ_detector`

or whatever name is chosen for the detector tracking envelope.

This mother contains the entire region relevant to tracking after source propagation.

So the tracking tree should begin there.

## What belongs in the tracking tree

The tracking tree should contain only volumes relevant to runtime detector response.

Examples:

- `LZ_detector`
- sensitive and passive LXe regions
- `LXeTPC`
- `FV`
- `Skin`
- `RFR`
- `DomeBarrel`
- `DomeBottomCap`
- field cage and internal structures that matter for tracking or veto
- later PMTs only if they are needed as tracking obstacles or veto-relevant structures

What does **not** automatically belong in the tracking tree:

- Ti cryostat shells
- flanges
- MLI
- any source-only object whose effect is already absorbed into the source-flux stage

## Why a single mixed tree is the wrong abstraction

At first glance it may seem attractive to keep one geometry and attach labels such as:

- `source`
- `tracking`
- `shared`

This is not enough.

The problem is not merely classification. It is that source and tracking use different algorithms.

Source stage:

- special propagation
- source-local self-shielding
- virtual detector-entry boundary
- no full event stack

Tracking stage:

- event-level FV prefilter
- sensitive-region hit logic
- veto logic
- full event tracking only for surviving events

Because the propagation model differs, the right design is two representations or two sections, not one mixed navigation tree with labels.

## Recommended data organization

Two clean options are acceptable:

1. separate JSON files
   - one for source geometry
   - one for tracking geometry

2. one JSON with separate top-level sections
   - `"source_geometry": ...`
   - `"tracking_geometry": ...`

The important part is conceptual separation.

For the current codebase, a separate source JSON or source section is likely the cleaner mental model.

## Cryostat example

For the cryostat source:

- source geometry contains inner and outer cryostat titanium parts
- source propagation computes the flux entering the detector
- tracking starts only once the gamma has entered `LZ_detector`

So in runtime tracking:

- we never need to worry whether a gamma is outside `LZ_detector`
- we never need `OCV` or `ICV` in the tracking tree if they serve only source modelling

This is the central simplification.

## Relationship to sensitive volumes and FV logic

This geometry split works naturally with the previously defined hit and prefilter philosophy.

Tracking stage:

- `FV` is the special prefilter target
- `LXeTPC` and `Skin` are sensitive
- passive LXe is not sensitive

Event-level fast logic:

- `0` gammas: reject
- `1` gamma: keep only if it reaches `FV`
- `2+` gammas: keep only if any gamma reaches `FV`

This logic assumes that the input gammas are already inside the runtime detector mother volume.

That is another reason to start tracking at `LZ_detector`, not outside it.

## Practical implication for current V2 work

The next geometry work should not try to force `OCV` and `ICV` into the tracking tree just because they are physically present.

Instead:

1. define the reduced tracking tree beginning at `LZ_detector`
2. move sensitive/passive detector regions under that mother
3. keep cryostat shells and other radioactive source objects in source geometry
4. use source propagation to produce the flux entering `LZ_detector`

This should guide the next refactor of `detector_lz_v2.json`.

## Summary

- source geometry and tracking geometry are different representations
- source objects are not just tagged tracking objects
- runtime navigation should start at `LZ_detector`
- anything outside `LZ_detector` belongs to source modelling unless proven otherwise
- this separation is the intended design going forward
