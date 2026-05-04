"""
Geometry module for the LXe Monte Carlo.

Follows a Geant4-inspired hierarchy:

1. **Geometric solids** (`Cyl`, `Box`): dimensions only.
2. **Logical volumes** (`LCyl`, `LBox`): solid + placement (position in MARS).
3. **Physical volumes** (`PCyl`, `PBox`): logical volume + material.
4. **Radioactive volumes** (`RCyl`): physical volume + specific activities.
5. **Detector**: MARS (mother volume, Vacuum) + list of physical volumes.

All coordinates are in the MARS (Mother Absolute Reference System) frame.
Orientation is z-axis aligned for all volumes (rotation deferred).
"""


# =====================================================================
# Geometric solids
# =====================================================================

"""
    Cyl

Cylindrical solid along z-axis, defined by radius and half-height.
"""
struct Cyl
    radius_cm::Float64
    half_height_cm::Float64
end

"""Volume of the cylinder [cm³]."""
volume(c::Cyl) = π * c.radius_cm^2 * 2.0 * c.half_height_cm

"""Lateral + end-cap surface area [cm²]."""
surface_area(c::Cyl) = 2π * c.radius_cm * (2.0 * c.half_height_cm + c.radius_cm)


"""
    Box

Rectangular box, defined by half-widths along x, y, z.
"""
struct Box
    half_x_cm::Float64
    half_y_cm::Float64
    half_z_cm::Float64
end

"""Volume of the box [cm³]."""
volume(b::Box) = 8.0 * b.half_x_cm * b.half_y_cm * b.half_z_cm

"""Surface area of the box [cm²]."""
surface_area(b::Box) = 8.0 * (b.half_x_cm * b.half_y_cm +
                               b.half_y_cm * b.half_z_cm +
                               b.half_z_cm * b.half_x_cm)


# =====================================================================
# Logical volumes (solid + placement)
# =====================================================================

"""
    LCyl

Logical cylinder: a `Cyl` placed at a position in MARS coordinates.
Axis is always along z (rotation deferred).
"""
struct LCyl
    solid::Cyl
    position::Vector{Float64}   # center [x, y, z] in MARS [cm]
end

"""Check if a point `pos` (MARS coordinates) is inside this logical cylinder."""
function is_inside(lc::LCyl, pos::Vector{Float64})::Bool
    dx = pos[1] - lc.position[1]
    dy = pos[2] - lc.position[2]
    dz = pos[3] - lc.position[3]
    dx^2 + dy^2 < lc.solid.radius_cm^2 && abs(dz) < lc.solid.half_height_cm
end


"""
    LBox

Logical box: a `Box` placed at a position in MARS coordinates.
"""
struct LBox
    solid::Box
    position::Vector{Float64}
end

"""Check if a point `pos` (MARS coordinates) is inside this logical box."""
function is_inside(lb::LBox, pos::Vector{Float64})::Bool
    abs(pos[1] - lb.position[1]) < lb.solid.half_x_cm &&
    abs(pos[2] - lb.position[2]) < lb.solid.half_y_cm &&
    abs(pos[3] - lb.position[3]) < lb.solid.half_z_cm
end


# =====================================================================
# Physical volumes (logical volume + material)
# =====================================================================

"""Abstract type for physical volumes."""
abstract type PhysicalVolume end

"""
    PCyl

Physical cylinder: a logical cylinder filled with a material.
"""
struct PCyl <: PhysicalVolume
    name::String
    logical::LCyl
    material::Material
end

"""Mass of the physical cylinder [g]."""
mass(pc::PCyl) = pc.material.density * volume(pc.logical.solid)

"""Check if a point is inside this physical volume."""
is_inside(pc::PCyl, pos::Vector{Float64}) = is_inside(pc.logical, pos)


"""
    PBox

Physical box: a logical box filled with a material.
"""
struct PBox <: PhysicalVolume
    name::String
    logical::LBox
    material::Material
end

"""Mass of the physical box [g]."""
mass(pb::PBox) = pb.material.density * volume(pb.logical.solid)

"""Check if a point is inside this physical volume."""
is_inside(pb::PBox, pos::Vector{Float64}) = is_inside(pb.logical, pos)


# =====================================================================
# Radioactive volumes
# =====================================================================

"""
    RCyl

Radioactive cylinder: a physical cylinder with specific activities
for U-238 chain (Bi-214, 2.448 MeV gamma) and Th-232 chain
(Tl-208, 2.615 MeV gamma).

Activities are in Bq/kg. Gamma fluxes are computed from activity × mass.
"""
struct RCyl
    pcyl::PCyl
    A_U238::Float64     # specific activity, U-238 chain [Bq/kg]
    A_Th232::Float64    # specific activity, Th-232 chain [Bq/kg]
end

"""Total U-238 activity [Bq]."""
activity_U238(rc::RCyl) = rc.A_U238 * mass(rc.pcyl) / 1000.0  # g → kg

"""Total Th-232 activity [Bq]."""
activity_Th232(rc::RCyl) = rc.A_Th232 * mass(rc.pcyl) / 1000.0

"""
    gamma_flux(rc::RCyl) -> (flux_Bi214, flux_Tl208)

Gamma flux [gammas/sec] for Bi-214 (2.448 MeV) and Tl-208 (2.615 MeV).
Assumes one gamma per decay (branching ratios can be applied externally).
"""
function gamma_flux(rc::RCyl)::Tuple{Float64,Float64}
    (activity_U238(rc), activity_Th232(rc))
end


# =====================================================================
# Detector: MARS + volumes
# =====================================================================

"""
    Detector

Complete detector geometry: a MARS mother volume (Vacuum) containing
one or more physical volumes. All volumes are placed in MARS coordinates.
"""
struct Detector
    name::String
    mars::PhysicalVolume
    volumes::Vector{PhysicalVolume}
end


"""
    active_volume(det::Detector) -> PhysicalVolume

Return the first active volume in the detector (the volume whose material
has `active=true`). Errors if none found.
"""
function active_volume(det::Detector)::PhysicalVolume
    for v in det.volumes
        v.material.active && return v
    end
    error("No active volume found in detector '$(det.name)'")
end


"""
    find_volume(det::Detector, pos::Vector{Float64}) -> Union{PhysicalVolume, Nothing}

Return the physical volume containing `pos`, or `nothing` if the point
is only in MARS (vacuum). Checks volumes in order; first match wins.
"""
function find_volume(det::Detector, pos::Vector{Float64})::Union{PhysicalVolume,Nothing}
    for v in det.volumes
        is_inside(v, pos) && return v
    end
    nothing
end


# =====================================================================
# Detector loader
# =====================================================================

"""
    load_detector(path::AbstractString, materials::Dict{String,Material}) -> Detector

Load a detector geometry from a JSON file. Volume materials are looked up
in the provided materials dictionary.
"""
function load_detector(path::AbstractString,
                       materials::Dict{String,Material})::Detector
    raw = open(path, "r") do io
        JSON.parse(io)
    end

    det_name = raw["name"]

    # Build MARS
    md = raw["mars"]
    mars_mat = materials[md["material"]]
    mars = _build_volume("MARS", md, mars_mat)

    # Build inner volumes
    volumes = PhysicalVolume[]
    for vd in raw["volumes"]
        mat = materials[vd["material"]]
        push!(volumes, _build_volume(vd["name"], vd, mat))
    end

    Detector(det_name, mars, volumes)
end


function _build_volume(name::String, d, mat::Material)::PhysicalVolume
    shape = lowercase(d["shape"])
    pos = Float64.(get(d, "position_cm", [0.0, 0.0, 0.0]))

    if shape == "cylinder"
        solid = Cyl(Float64(d["radius_cm"]), Float64(d["half_height_cm"]))
        logical = LCyl(solid, pos)
        return PCyl(name, logical, mat)
    elseif shape == "box"
        solid = Box(Float64(d["half_x_cm"]), Float64(d["half_y_cm"]), Float64(d["half_z_cm"]))
        logical = LBox(solid, pos)
        return PBox(name, logical, mat)
    else
        error("Unknown shape '$shape' for volume '$name'")
    end
end
