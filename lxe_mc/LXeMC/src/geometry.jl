"""
Legacy stack/source geometry module for the LXe Monte Carlo.

This file still supports:
- the legacy flat `Detector` model used by source-transport utilities
- low-level geometry primitive tests
- historical stack-facing helpers that predate the canonical tracking
  geometry

The canonical detector workflow now lives in `geometry2.jl`, while this
module remains the support layer for the old flat detector representation.

It follows a Geant4-inspired hierarchy on top of the shared
geometry-core primitives:

1. **Logical volumes** (`LCyl`, `LCylShell`, `LBox`): solid + placement (position in MARS).
2. **Physical volumes** (`PCyl`, `PCylShell`, `PBox`): logical volume + material.
3. **Detector**: MARS (mother volume, Vacuum) + list of physical volumes.

All coordinates are in the MARS (Mother Absolute Reference System) frame.
Orientation is z-axis aligned for all volumes (rotation deferred).
"""

# =====================================================================
# Legacy flat detector: MARS + volumes
# =====================================================================

"""
    Detector

Complete detector geometry: MARS mother volume (Vacuum) containing
one or more physical volumes placed in MARS coordinates.

The optional `fiducial` field holds the fiducial volume (a virtual
volume used for analysis cuts, not transport). It is set automatically
by `load_detector` when a volume has `"fiducial": true` in the JSON.
"""
struct Detector
    name::String
    mars::PhysicalVolume
    volumes::Vector{PhysicalVolume}
    fiducial::Union{PhysicalVolume,Nothing}
end


"""
    active_volume(det) -> PhysicalVolume

Return the first volume whose material has `active=true`.
"""
function active_volume(det::Detector)::PhysicalVolume
    for v in det.volumes
        v.material.active && return v
    end
    error("No active volume found in detector '$(det.name)'")
end


"""
    fiducial_volume(det) -> PhysicalVolume

Return the fiducial volume. Errors if none defined.
"""
function fiducial_volume(det::Detector)::PhysicalVolume
    det.fiducial !== nothing && return det.fiducial
    error("No fiducial volume defined in detector '$(det.name)'")
end


"""
    find_volume(det, pos) -> Union{PhysicalVolume, Nothing}

Return the physical volume containing `pos`, or `nothing` if only in MARS.
Checks volumes in order; first match wins.
"""
function find_volume(det::Detector, pos::Vector{Float64})::Union{PhysicalVolume,Nothing}
    for v in det.volumes
        is_inside(v, pos) && return v
    end
    nothing
end


"""
    next_volume(pos, dir, det) -> (PhysicalVolume or nothing, Float64)

Find the nearest volume along ray (pos, dir) from a point in vacuum.
Returns (volume, distance) or (nothing, Inf) if ray escapes MARS.
"""
function next_volume(pos::Vector{Float64}, dir::Vector{Float64},
                     det::Detector)::Tuple{Union{PhysicalVolume,Nothing},Float64}
    best_vol = nothing
    best_t = Inf

    for v in det.volumes
        t = distance_to_entry(pos, dir, v.logical)
        if t < best_t
            best_t = t
            best_vol = v
        end
    end

    (best_vol, best_t)
end


# =====================================================================
# Detector loader
# =====================================================================

"""
    load_detector(path, materials) -> Detector

Load detector geometry from JSON. Supported shapes: `"cylinder"`,
`"cylinder_shell"`, `"box"`. Materials looked up by name.
"""
function load_detector(path::AbstractString,
                       materials::Dict{String,Material})::Detector
    raw = open(path, "r") do io
        JSON.parse(io)
    end

    det_name = raw["name"]

    md = raw["mars"]
    mars_mat = materials[md["material"]]
    mars = _build_volume("MARS", md, mars_mat)

    volumes = PhysicalVolume[]
    fid = nothing
    for vd in raw["volumes"]
        mat = materials[vd["material"]]
        vol = _build_volume(vd["name"], vd, mat)
        if Bool(get(vd, "fiducial", false))
            fid = vol
        else
            push!(volumes, vol)
        end
    end

    Detector(det_name, mars, volumes, fid)
end


function _build_volume(name::String, d, mat::Material)::PhysicalVolume
    shape = lowercase(d["shape"])
    pos = Float64.(get(d, "position_cm", [0.0, 0.0, 0.0]))

    if shape == "cylinder"
        solid = Cyl(Float64(d["radius_cm"]), Float64(d["half_height_cm"]))
        logical = LCyl(solid, pos)
        return PCyl(name, logical, mat)
    elseif shape == "cylinder_shell"
        solid = CylShell(Float64(d["R_inner_cm"]),
                         Float64(d["wall_thickness_cm"]),
                         Float64(d["half_height_cm"]))
        logical = LCylShell(solid, pos)
        return PCylShell(name, logical, mat)
    elseif shape == "box"
        solid = Box(Float64(d["half_x_cm"]), Float64(d["half_y_cm"]), Float64(d["half_z_cm"]))
        logical = LBox(solid, pos)
        return PBox(name, logical, mat)
    elseif shape == "disk"
        solid = Disk(Float64(d["radius_cm"]),
                     Float64(d["wall_thickness_cm"]),
                     Float64(get(d, "aspect_ratio", Inf)))
        orient = Symbol(get(d, "orientation", "up"))
        logical = LDisk(solid, pos, orient)
        return PDisk(name, logical, mat)
    else
        error("Unknown shape '$shape' for volume '$name'")
    end
end
