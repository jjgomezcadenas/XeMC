"""
Geometry module for the LXe Monte Carlo.

Follows a Geant4-inspired hierarchy:

1. **Geometric solids** (`Cyl`, `CylShell`, `Box`): dimensions only.
2. **Logical volumes** (`LCyl`, `LCylShell`, `LBox`): solid + placement (position in MARS).
3. **Physical volumes** (`PCyl`, `PCylShell`, `PBox`): logical volume + material.
4. **Radioactive volumes** (`RCyl`, `RCylShell`): physical volume + specific activities.
5. **Detector**: MARS (mother volume, Vacuum) + list of physical volumes.

All coordinates are in the MARS (Mother Absolute Reference System) frame.
Orientation is z-axis aligned for all volumes (rotation deferred).
"""


# =====================================================================
# Geometric solids
# =====================================================================

"""
    Cyl(radius_cm, half_height_cm)

Solid cylinder along z-axis. Used for active volumes (e.g., LXe TPC).
"""
struct Cyl
    radius_cm::Float64
    half_height_cm::Float64

    function Cyl(radius_cm::Real, half_height_cm::Real)
        radius_cm > 0 || error("Cyl: radius must be positive (got $radius_cm)")
        half_height_cm > 0 || error("Cyl: half_height must be positive (got $half_height_cm)")
        new(Float64(radius_cm), Float64(half_height_cm))
    end
end

"""Volume of the solid cylinder [cm³]: π R² × 2H."""
volume(c::Cyl) = π * c.radius_cm^2 * 2.0 * c.half_height_cm

"""Total surface area (lateral + two end caps) [cm²]."""
surface_area(c::Cyl) = 2π * c.radius_cm * (2.0 * c.half_height_cm + c.radius_cm)


"""
    CylShell(R_inner, wall_thickness, half_height)

Cylindrical shell (hollow cylinder) along z-axis. Inner surface at
radius `R_inner`, outer at `R_inner + wall_thickness`. Axial extent
±`half_height` from center. Used for passive components: cryostat
walls, field cage rings, PTFE reflectors.
"""
struct CylShell
    R_inner_cm::Float64
    wall_thickness_cm::Float64
    half_height_cm::Float64

    function CylShell(R_inner::Real, wall_thickness::Real, half_height::Real)
        R_inner > 0 || error("CylShell: R_inner must be positive (got $R_inner)")
        wall_thickness > 0 || error("CylShell: wall_thickness must be positive (got $wall_thickness)")
        half_height > 0 || error("CylShell: half_height must be positive (got $half_height)")
        new(Float64(R_inner), Float64(wall_thickness), Float64(half_height))
    end
end

"""Outer radius [cm]."""
R_outer(cs::CylShell) = cs.R_inner_cm + cs.wall_thickness_cm

"""Shell volume π(R_outer² − R_inner²) × 2H [cm³]."""
volume(cs::CylShell) = π * (R_outer(cs)^2 - cs.R_inner_cm^2) * 2.0 * cs.half_height_cm

"""Volume enclosed by the inner cylinder [cm³]."""
volume_inner(cs::CylShell) = π * cs.R_inner_cm^2 * 2.0 * cs.half_height_cm

"""Inner lateral surface area, 2π R_inner × 2H [cm²]."""
surface_area_inner(cs::CylShell) = 2π * cs.R_inner_cm * 2.0 * cs.half_height_cm

"""Outer lateral surface area, 2π R_outer × 2H [cm²]."""
surface_area_outer(cs::CylShell) = 2π * R_outer(cs) * 2.0 * cs.half_height_cm


"""
    Box(half_x_cm, half_y_cm, half_z_cm)

Rectangular box, defined by half-widths along x, y, z.
"""
struct Box
    half_x_cm::Float64
    half_y_cm::Float64
    half_z_cm::Float64

    function Box(half_x::Real, half_y::Real, half_z::Real)
        half_x > 0 || error("Box: half_x must be positive (got $half_x)")
        half_y > 0 || error("Box: half_y must be positive (got $half_y)")
        half_z > 0 || error("Box: half_z must be positive (got $half_z)")
        new(Float64(half_x), Float64(half_y), Float64(half_z))
    end
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
    LCyl(solid, position)

Logical cylinder: a `Cyl` placed at `position` [x,y,z] in MARS [cm].
Axis is along z (rotation deferred).
"""
struct LCyl
    solid::Cyl
    position::Vector{Float64}
end

"""True if `pos` (MARS coordinates) is inside this solid cylinder."""
function is_inside(lc::LCyl, pos::Vector{Float64})::Bool
    dx = pos[1] - lc.position[1]
    dy = pos[2] - lc.position[2]
    dz = pos[3] - lc.position[3]
    dx^2 + dy^2 < lc.solid.radius_cm^2 && abs(dz) < lc.solid.half_height_cm
end


"""
    LCylShell(solid, position)

Logical cylindrical shell: a `CylShell` placed at `position` in MARS [cm].
"""
struct LCylShell
    solid::CylShell
    position::Vector{Float64}
end

"""True if `pos` is inside the shell (between inner and outer radii, within height)."""
function is_inside(lcs::LCylShell, pos::Vector{Float64})::Bool
    dx = pos[1] - lcs.position[1]
    dy = pos[2] - lcs.position[2]
    dz = pos[3] - lcs.position[3]
    r2 = dx^2 + dy^2
    r2 >= lcs.solid.R_inner_cm^2 && r2 < R_outer(lcs.solid)^2 && abs(dz) < lcs.solid.half_height_cm
end


"""
    LBox(solid, position)

Logical box: a `Box` placed at `position` in MARS [cm].
"""
struct LBox
    solid::Box
    position::Vector{Float64}
end

"""True if `pos` is inside this box."""
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
    PCyl <: PhysicalVolume

Physical solid cylinder: logical cylinder + material.
"""
struct PCyl <: PhysicalVolume
    name::String
    logical::LCyl
    material::Material
end

"""Mass [g] = density × volume."""
mass(pc::PCyl) = pc.material.density * volume(pc.logical.solid)

"""True if `pos` is inside."""
is_inside(pc::PCyl, pos::Vector{Float64}) = is_inside(pc.logical, pos)


"""
    PCylShell <: PhysicalVolume

Physical cylindrical shell: logical shell + material.
Used for cryostat walls, field cage, reflectors.
"""
struct PCylShell <: PhysicalVolume
    name::String
    logical::LCylShell
    material::Material
end

"""Mass [g] = density × shell volume."""
mass(pcs::PCylShell) = pcs.material.density * volume(pcs.logical.solid)

"""True if `pos` is inside the shell."""
is_inside(pcs::PCylShell, pos::Vector{Float64}) = is_inside(pcs.logical, pos)


"""
    PBox <: PhysicalVolume

Physical box: logical box + material.
"""
struct PBox <: PhysicalVolume
    name::String
    logical::LBox
    material::Material
end

"""Mass [g] = density × volume."""
mass(pb::PBox) = pb.material.density * volume(pb.logical.solid)

"""True if `pos` is inside."""
is_inside(pb::PBox, pos::Vector{Float64}) = is_inside(pb.logical, pos)


# =====================================================================
# Radioactive volumes
# =====================================================================

"""
    RCyl

Radioactive solid cylinder. Wraps a `PCyl` with specific activities
for U-238 chain (Bi-214, 2.448 MeV) and Th-232 chain (Tl-208, 2.615 MeV).
Activities in Bq/kg.
"""
struct RCyl
    phys::PCyl
    A_U238::Float64     # Bq/kg
    A_Th232::Float64    # Bq/kg
end

"""Total U-238 activity [Bq]."""
activity_U238(rc::RCyl) = rc.A_U238 * mass(rc.phys) / 1000.0

"""Total Th-232 activity [Bq]."""
activity_Th232(rc::RCyl) = rc.A_Th232 * mass(rc.phys) / 1000.0

"""Gamma flux (Bi-214, Tl-208) [gammas/sec]. One gamma per decay assumed."""
gamma_flux(rc::RCyl) = (activity_U238(rc), activity_Th232(rc))


"""
    RCylShell

Radioactive cylindrical shell. Wraps a `PCylShell` with specific activities.
Typical use: contaminated cryostat wall, field cage ring, PTFE reflector.
"""
struct RCylShell
    phys::PCylShell
    A_U238::Float64     # Bq/kg
    A_Th232::Float64    # Bq/kg
end

"""Total U-238 activity [Bq]."""
activity_U238(rcs::RCylShell) = rcs.A_U238 * mass(rcs.phys) / 1000.0

"""Total Th-232 activity [Bq]."""
activity_Th232(rcs::RCylShell) = rcs.A_Th232 * mass(rcs.phys) / 1000.0

"""Gamma flux (Bi-214, Tl-208) [gammas/sec]."""
gamma_flux(rcs::RCylShell) = (activity_U238(rcs), activity_Th232(rcs))


# =====================================================================
# Detector: MARS + volumes
# =====================================================================

"""
    Detector

Complete detector geometry: MARS mother volume (Vacuum) containing
one or more physical volumes placed in MARS coordinates.
"""
struct Detector
    name::String
    mars::PhysicalVolume
    volumes::Vector{PhysicalVolume}
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
    else
        error("Unknown shape '$shape' for volume '$name'")
    end
end
