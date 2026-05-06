# Geometry V2 Prefilter And Sensitivity

## Purpose

This note records the event-level logic that should be kept separate from the geometry/container discussion.

The important concepts are:

- source geometry vs tracking geometry
- mother/daughter material precedence
- sensitive volumes
- `FV` as a special fast-prefilter target

## Source geometry vs tracking geometry

The event-level fast logic assumes a reduced tracking geometry.

Radioactive source objects such as:

- Ti cryostat shells
- flanges
- MLI
- later PMT source objects

are primarily part of source geometry, not necessarily part of the runtime tracking tree.

Their effect should enter through the source-flux stage.

The tracking geometry should then start from the detector cavities and internal detector regions relevant for:

- FV reachability
- sensitive-region hits
- veto logic

## Material precedence

The intended rule is:

- the mother defines the default material
- daughters displace the mother in the space they occupy

Examples:

- `ICV_void` may be vacuum
- `LXe_detector` displaces that vacuum with LXe
- PTFE, Ti, PMTs, grids, etc. displace LXe where they are placed

This is the basis for navigation and transport.

## Sensitive volumes

Not all LXe subregions need the same analysis behavior.

Current intended classification:

- `LXeTPC`: sensitive
- `Skin`: sensitive
- passive LXe (`RFR`, `Dome`, etc.): not sensitive

Meaning:

- deposits/clusters in sensitive regions may become hits
- deposits/clusters in passive LXe do not become hits

## FV special role

`FV` is special in two ways:

1. it is a sensitive LXe subregion
2. it is the target of the event-level geometric prefilter

So `FV` should not be modeled only as “another sensitive region”.

It has a distinct role:

- `fv_target = true`

or an equivalent dedicated semantic concept.

## Hit logic

Intended separation:

1. transport produces deposits/clusters
2. sensitive regions convert qualifying clusters into hits
3. analysis decisions are made from hits

Examples:

- a cluster in `LXeTPC` above `veto_TPC` contributes to veto logic
- a cluster in `Skin` above `veto_skin` contributes to veto logic
- clusters in passive LXe do not become hits

SS/MS classification should operate on hits or hit-like clustered objects, not on raw transport history.

## Fast geometric prefilter

The code must reject most events before full stack development whenever possible.

Event-level rule:

### 0 gammas

- reject immediately

### 1 gamma

- does the gamma reach `FV`?
- if no: reject immediately
- if yes: keep the event and track normally

### 2 or more gammas

- does any gamma reach `FV`?
- if no: reject immediately
- if yes: keep the whole event and track all gammas normally

This is intentionally only a geometric prefilter.

It does **not** apply veto decisions yet.

## Why veto is later

The point of the prefilter is:

- decide the hard geometric rejection early

The point of the later hit logic is:

- preserve flexibility in veto thresholds and clustering cuts

So:

- prefilter asks: can any gamma reach `FV`?
- full event analysis later asks: do the sensitive-region hits reject the event?

## Consequence for implementation

Geometry V2 should eventually provide at least:

- material / region ownership by strict descent
- a way to identify `FV` as the prefilter target
- a way to identify sensitive regions and their hit policy
- a reduced tracking tree that excludes source-only geometry when appropriate

Transport should not need to know the full analysis policy beyond those semantics.
