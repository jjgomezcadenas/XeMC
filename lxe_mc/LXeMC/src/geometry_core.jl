"""
Shared geometry-core primitives used by both the canonical tracking path
and the legacy flat-detector support layer.
"""


# =====================================================================
# Geometric solids
# =====================================================================

"""
    Cyl(radius_cm, half_height_cm)

Solid cylinder along z-axis.
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

volume(c::Cyl) = π * c.radius_cm^2 * 2.0 * c.half_height_cm
surface_area(c::Cyl) = 2π * c.radius_cm * (2.0 * c.half_height_cm + c.radius_cm)


"""
    CylShell(R_inner, wall_thickness, half_height)

Cylindrical shell (hollow cylinder) along z-axis.
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

R_outer(cs::CylShell) = cs.R_inner_cm + cs.wall_thickness_cm
volume(cs::CylShell) = π * (R_outer(cs)^2 - cs.R_inner_cm^2) * 2.0 * cs.half_height_cm
volume_inner(cs::CylShell) = π * cs.R_inner_cm^2 * 2.0 * cs.half_height_cm
surface_area_inner(cs::CylShell) = 2π * cs.R_inner_cm * 2.0 * cs.half_height_cm
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

volume(b::Box) = 8.0 * b.half_x_cm * b.half_y_cm * b.half_z_cm
surface_area(b::Box) = 8.0 * (b.half_x_cm * b.half_y_cm +
                               b.half_y_cm * b.half_z_cm +
                               b.half_z_cm * b.half_x_cm)


"""
    Disk(radius_cm, wall_thickness_cm, aspect_ratio)

Cryostat head or end cap. A thin shell whose shape ranges from flat to
ellipsoidal.
"""
struct Disk
    radius_cm::Float64
    wall_thickness_cm::Float64
    aspect_ratio::Float64

    function Disk(radius::Real, wall_thickness::Real, aspect_ratio::Real)
        radius > 0 || error("Disk: radius must be positive (got $radius)")
        wall_thickness >= 0 || error("Disk: wall_thickness must be ≥ 0 (got $wall_thickness)")
        aspect_ratio > 0 || error("Disk: aspect_ratio must be positive (got $aspect_ratio)")
        new(Float64(radius), Float64(wall_thickness), Float64(aspect_ratio))
    end
end

depth(d::Disk) = isinf(d.aspect_ratio) ? 0.0 : d.radius_cm / d.aspect_ratio
is_flat(d::Disk) = isinf(d.aspect_ratio)

function surface_area_inner(d::Disk)::Float64
    is_flat(d) && return π * d.radius_cm^2
    n = d.aspect_ratio
    n ≈ 1.0 && return 2π * d.radius_cm^2

    a = d.radius_cm
    c = d.radius_cm / n
    if c < a
        e = sqrt(1.0 - (c / a)^2)
        S_full = 2π * a^2 + π * (c^2 / e) * log((1 + e) / (1 - e))
    else
        e = sqrt(1.0 - (a / c)^2)
        S_full = 2π * a^2 + 2π * a * c * asin(e) / e
    end
    S_full / 2.0
end

volume(d::Disk) = surface_area_inner(d) * d.wall_thickness_cm


struct Cap
    radius_cm::Float64
    aspect_ratio::Float64

    function Cap(radius::Real, aspect_ratio::Real)
        radius > 0 || error("Cap: radius must be positive (got $radius)")
        aspect_ratio > 0 || error("Cap: aspect_ratio must be positive (got $aspect_ratio)")
        new(Float64(radius), Float64(aspect_ratio))
    end
end

depth(c::Cap) = c.radius_cm / c.aspect_ratio


struct DomedContainer
    radius_cm::Float64
    barrel_half_height_cm::Float64
    top_cap::Cap
    bottom_cap::Cap

    function DomedContainer(radius::Real,
                            barrel_half_height::Real,
                            top_cap::Cap,
                            bottom_cap::Cap)
        radius > 0 || error("DomedContainer: radius must be positive (got $radius)")
        barrel_half_height >= 0 || error("DomedContainer: barrel_half_height must be non-negative (got $barrel_half_height)")
        top_cap.radius_cm ≈ Float64(radius) || error("DomedContainer: top cap radius must match container radius")
        bottom_cap.radius_cm ≈ Float64(radius) || error("DomedContainer: bottom cap radius must match container radius")
        new(Float64(radius), Float64(barrel_half_height), top_cap, bottom_cap)
    end
end


struct CappedCylinder
    radius_cm::Float64
    barrel_half_height_cm::Float64
    top_cap::Union{Nothing,Cap}
    bottom_cap::Union{Nothing,Cap}

    function CappedCylinder(radius::Real,
                            barrel_half_height::Real;
                            top_cap::Union{Nothing,Cap}=nothing,
                            bottom_cap::Union{Nothing,Cap}=nothing)
        radius > 0 || error("CappedCylinder: radius must be positive (got $radius)")
        barrel_half_height >= 0 || error("CappedCylinder: barrel_half_height must be non-negative (got $barrel_half_height)")
        top_cap === nothing || top_cap.radius_cm ≈ Float64(radius) || error("CappedCylinder: top cap radius must match cylinder radius")
        bottom_cap === nothing || bottom_cap.radius_cm ≈ Float64(radius) || error("CappedCylinder: bottom cap radius must match cylinder radius")
        new(Float64(radius), Float64(barrel_half_height), top_cap, bottom_cap)
    end
end


# =====================================================================
# Logical volumes (solid + placement)
# =====================================================================

struct LCyl
    solid::Cyl
    position::Vector{Float64}
end

function is_inside(lc::LCyl, pos::Vector{Float64})::Bool
    dx = pos[1] - lc.position[1]
    dy = pos[2] - lc.position[2]
    dz = pos[3] - lc.position[3]
    dx^2 + dy^2 < lc.solid.radius_cm^2 && abs(dz) < lc.solid.half_height_cm
end


struct LCylShell
    solid::CylShell
    position::Vector{Float64}
end

function is_inside(lcs::LCylShell, pos::Vector{Float64})::Bool
    dx = pos[1] - lcs.position[1]
    dy = pos[2] - lcs.position[2]
    dz = pos[3] - lcs.position[3]
    r2 = dx^2 + dy^2
    r2 >= lcs.solid.R_inner_cm^2 && r2 < R_outer(lcs.solid)^2 && abs(dz) < lcs.solid.half_height_cm
end


struct LBox
    solid::Box
    position::Vector{Float64}
end

function is_inside(lb::LBox, pos::Vector{Float64})::Bool
    abs(pos[1] - lb.position[1]) < lb.solid.half_x_cm &&
    abs(pos[2] - lb.position[2]) < lb.solid.half_y_cm &&
    abs(pos[3] - lb.position[3]) < lb.solid.half_z_cm
end


struct LDisk
    solid::Disk
    position::Vector{Float64}
    orientation::Symbol

    function LDisk(solid::Disk, position::Vector{Float64}, orientation::Symbol)
        orientation in (:up, :down) || error("LDisk: orientation must be :up or :down (got $orientation)")
        new(solid, position, orientation)
    end
end

function is_inside(ld::LDisk, pos::Vector{Float64})::Bool
    dx = pos[1] - ld.position[1]
    dy = pos[2] - ld.position[2]
    dz = pos[3] - ld.position[3]
    r2 = dx^2 + dy^2
    R = ld.solid.radius_cm
    r2 > R^2 && return false

    if is_flat(ld.solid)
        sgn = ld.orientation === :up ? 1.0 : -1.0
        return 0.0 <= sgn * dz <= ld.solid.wall_thickness_cm
    end

    c_inner = depth(ld.solid)
    c_outer = c_inner + ld.solid.wall_thickness_cm
    sgn = ld.orientation === :up ? 1.0 : -1.0
    dz_rel = sgn * dz
    dz_rel < 0.0 && return false

    r_frac2 = r2 / R^2
    inside_outer = r_frac2 + (dz_rel / c_outer)^2 <= 1.0
    outside_inner = r_frac2 + (dz_rel / c_inner)^2 >= 1.0

    inside_outer && outside_inner
end


# =====================================================================
# Physical volumes (logical volume + material)
# =====================================================================

abstract type PhysicalVolume end

struct PCyl <: PhysicalVolume
    name::String
    logical::LCyl
    material::Material
end

mass(pc::PCyl) = pc.material.density * volume(pc.logical.solid)
is_inside(pc::PCyl, pos::Vector{Float64}) = is_inside(pc.logical, pos)


struct PCylShell <: PhysicalVolume
    name::String
    logical::LCylShell
    material::Material
end

mass(pcs::PCylShell) = pcs.material.density * volume(pcs.logical.solid)
is_inside(pcs::PCylShell, pos::Vector{Float64}) = is_inside(pcs.logical, pos)


struct PBox <: PhysicalVolume
    name::String
    logical::LBox
    material::Material
end

mass(pb::PBox) = pb.material.density * volume(pb.logical.solid)
is_inside(pb::PBox, pos::Vector{Float64}) = is_inside(pb.logical, pos)


struct PDisk <: PhysicalVolume
    name::String
    logical::LDisk
    material::Material
end

mass(pd::PDisk) = pd.material.density * volume(pd.logical.solid)
is_inside(pd::PDisk, pos::Vector{Float64}) = is_inside(pd.logical, pos)


# =====================================================================
# Ray intersections on shared logical / physical volume abstractions
# =====================================================================

function distance_to_exit(pos::Vector{Float64}, dir::Vector{Float64}, lc::LCyl)::Float64
    R = lc.solid.radius_cm
    H = lc.solid.half_height_cm
    cx, cy, cz = lc.position

    dx = pos[1] - cx
    dy = pos[2] - cy
    dz = pos[3] - cz

    t_min = Inf
    a = dir[1]^2 + dir[2]^2
    if a > 1e-20
        b = 2.0 * (dx * dir[1] + dy * dir[2])
        c = dx^2 + dy^2 - R^2
        disc = b^2 - 4 * a * c
        if disc >= 0.0
            sq = sqrt(disc)
            for t in [(-b + sq) / (2a), (-b - sq) / (2a)]
                if t > 1e-10 && abs(dz + t * dir[3]) < H
                    t_min = min(t_min, t)
                end
            end
        end
    end

    if abs(dir[3]) > 1e-20
        for z_face in [H, -H]
            t = (z_face - dz) / dir[3]
            if t > 1e-10
                rx = dx + t * dir[1]
                ry = dy + t * dir[2]
                if rx^2 + ry^2 < R^2
                    t_min = min(t_min, t)
                end
            end
        end
    end

    t_min
end


function distance_to_entry(pos::Vector{Float64}, dir::Vector{Float64}, lc::LCyl)::Float64
    R = lc.solid.radius_cm
    H = lc.solid.half_height_cm
    cx, cy, cz = lc.position

    dx = pos[1] - cx
    dy = pos[2] - cy
    dz = pos[3] - cz

    t_min = Inf
    a = dir[1]^2 + dir[2]^2
    if a > 1e-20
        b = 2.0 * (dx * dir[1] + dy * dir[2])
        c = dx^2 + dy^2 - R^2
        disc = b^2 - 4 * a * c
        if disc >= 0.0
            sq = sqrt(disc)
            for t in [(-b - sq) / (2a), (-b + sq) / (2a)]
                if t > 1e-10 && abs(dz + t * dir[3]) < H
                    t_min = min(t_min, t)
                end
            end
        end
    end

    if abs(dir[3]) > 1e-20
        for z_face in [H, -H]
            t = (z_face - dz) / dir[3]
            if t > 1e-10
                rx = dx + t * dir[1]
                ry = dy + t * dir[2]
                if rx^2 + ry^2 < R^2
                    t_min = min(t_min, t)
                end
            end
        end
    end

    t_min
end


function distance_to_entry(pos::Vector{Float64}, dir::Vector{Float64}, lcs::LCylShell)::Float64
    R_i = lcs.solid.R_inner_cm
    R_o = R_outer(lcs.solid)
    H = lcs.solid.half_height_cm
    cx, cy, cz = lcs.position

    dx = pos[1] - cx
    dy = pos[2] - cy
    dz = pos[3] - cz

    t_min = Inf
    a = dir[1]^2 + dir[2]^2

    for R in [R_o, R_i]
        if a > 1e-20
            b = 2.0 * (dx * dir[1] + dy * dir[2])
            c = dx^2 + dy^2 - R^2
            disc = b^2 - 4 * a * c
            if disc >= 0.0
                sq = sqrt(disc)
                for t in [(-b - sq) / (2a), (-b + sq) / (2a)]
                    if t > 1e-10 && abs(dz + t * dir[3]) < H
                        rx = dx + t * dir[1]
                        ry = dy + t * dir[2]
                        r2 = rx^2 + ry^2
                        if r2 >= R_i^2 - 1e-6 && r2 <= R_o^2 + 1e-6
                            t_min = min(t_min, t)
                        end
                    end
                end
            end
        end
    end

    if abs(dir[3]) > 1e-20
        for z_face in [H, -H]
            t = (z_face - dz) / dir[3]
            if t > 1e-10
                rx = dx + t * dir[1]
                ry = dy + t * dir[2]
                r2 = rx^2 + ry^2
                if r2 >= R_i^2 && r2 <= R_o^2
                    t_min = min(t_min, t)
                end
            end
        end
    end

    t_min
end


function distance_to_exit(pos::Vector{Float64}, dir::Vector{Float64}, lcs::LCylShell)::Float64
    R_i = lcs.solid.R_inner_cm
    R_o = R_outer(lcs.solid)
    H = lcs.solid.half_height_cm
    cx, cy, cz = lcs.position

    dx = pos[1] - cx
    dy = pos[2] - cy
    dz = pos[3] - cz

    t_min = Inf
    a = dir[1]^2 + dir[2]^2

    for R in [R_i, R_o]
        if a > 1e-20
            b = 2.0 * (dx * dir[1] + dy * dir[2])
            c = dx^2 + dy^2 - R^2
            disc = b^2 - 4 * a * c
            if disc >= 0.0
                sq = sqrt(disc)
                for t in [(-b - sq) / (2a), (-b + sq) / (2a)]
                    if t > 1e-10 && abs(dz + t * dir[3]) < H
                        t_min = min(t_min, t)
                    end
                end
            end
        end
    end

    if abs(dir[3]) > 1e-20
        for z_face in [H, -H]
            t = (z_face - dz) / dir[3]
            if t > 1e-10
                rx = dx + t * dir[1]
                ry = dy + t * dir[2]
                r2 = rx^2 + ry^2
                if r2 >= R_i^2 && r2 <= R_o^2
                    t_min = min(t_min, t)
                end
            end
        end
    end

    t_min
end


function distance_to_entry(pos::Vector{Float64}, dir::Vector{Float64}, lv::LBox)::Float64
    ds = 0.1
    for i in 1:10000
        t = i * ds
        test_pos = pos .+ dir .* t
        if is_inside(lv, test_pos)
            return t
        end
    end
    Inf
end


function distance_to_entry(pos::Vector{Float64}, dir::Vector{Float64}, ld::LDisk)::Float64
    ds = 0.1
    for i in 1:10000
        t = i * ds
        test_pos = pos .+ dir .* t
        if is_inside(ld, test_pos)
            return t
        end
    end
    Inf
end


function distance_to_entry(pos::Vector{Float64}, dir::Vector{Float64}, vol::PhysicalVolume)::Float64
    ds = 0.1
    for i in 1:10000
        t = i * ds
        test_pos = pos .+ dir .* t
        if is_inside(vol, test_pos)
            return t
        end
    end
    Inf
end


function distance_to_exit(pos::Vector{Float64}, dir::Vector{Float64}, vol::PhysicalVolume)::Float64
    ds = 0.1
    for i in 1:10000
        t = i * ds
        test_pos = pos .+ dir .* t
        if !is_inside(vol, test_pos)
            return t
        end
    end
    Inf
end
