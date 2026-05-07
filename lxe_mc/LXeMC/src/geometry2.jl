"""
Tracking geometry tree: explicit detector hierarchy and semantic metadata.

This file defines the canonical tracking-geometry representation used by
the current fast-kernel workflow. It parses the tracking JSON schema,
builds parent/child relationships, and exposes inspection and compiled
geometry helpers.
"""


@enum RegionTag begin
    TAG_WORLD
    TAG_VACUUM
    TAG_STRUCTURAL
    TAG_TPC_ACTIVE
    TAG_FV
    TAG_SKIN
    TAG_PASSIVE_LXE
end


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


struct Placement
    position_cm::Vector{Float64}
    orientation::Symbol
end


struct LogicalVolume
    name::String
    solid::Union{Cyl,CylShell,Box,Disk,Cap,DomedContainer,CappedCylinder}
    material::Material
    tag::RegionTag
    role::String
    inFastKernel::Bool
    isXe::Bool
    sensitive::Bool
    ecut_keV::Float64
    dz_mm::Float64
    fv_target::Bool
    approximation::String
    activity::Dict{String,Float64}
end


mutable struct DetectorNode
    id::Int
    parent_id::Int
    children_ids::Vector{Int}
    lv::LogicalVolume
    placement::Placement
end


struct TrackingDetector
    name::String
    nodes::Vector{DetectorNode}
    name_to_id::Dict{String,Int}
    root_id::Int
end


struct FastKernelRegion
    id::Int
    name::String
    role::String
    kind::Symbol
    material::Material
    isXe::Bool
    sensitive::Bool
    ecut_keV::Float64
    fv_target::Bool
    rmin_cm::Float64
    rmax_cm::Float64
    zmin_cm::Float64
    zmax_cm::Float64
    has_top_cap::Bool
    has_bottom_cap::Bool
    top_cap_radius_cm::Float64
    top_cap_aspect_ratio::Float64
    bottom_cap_radius_cm::Float64
    bottom_cap_aspect_ratio::Float64
end


struct FastKernelGeometry
    regions::Vector{FastKernelRegion}
    name_to_index::Dict{String,Int}
    detector_region::Int
    air_region::Int
    ptfe_region::Int
    rings_region::Int
    fv_region::Int
    top_active_region::Int
    barrel_active_region::Int
    bottom_active_region::Int
    skin_region::Int
    passive_region::Int
end

const GEOM_VALIDATE_TOL_CM = 0.1


region_tag(node::DetectorNode) = node.lv.tag

is_fv(node::DetectorNode)::Bool = region_tag(node) == TAG_FV

is_active_lxe(node::DetectorNode)::Bool = region_tag(node) in (TAG_TPC_ACTIVE, TAG_FV)

is_veto_lxe(node::DetectorNode)::Bool = region_tag(node) in (TAG_TPC_ACTIVE, TAG_SKIN)

is_passive_lxe(node::DetectorNode)::Bool = region_tag(node) == TAG_PASSIVE_LXE

is_structural(node::DetectorNode)::Bool = region_tag(node) == TAG_STRUCTURAL

is_vacuum(node::DetectorNode)::Bool = region_tag(node) in (TAG_WORLD, TAG_VACUUM)

is_sensitive(node::DetectorNode)::Bool = node.lv.sensitive

is_fv_target(node::DetectorNode)::Bool = node.lv.fv_target

ecut_keV(node::DetectorNode)::Float64 = node.lv.ecut_keV

dz_mm(node::DetectorNode)::Float64 = node.lv.dz_mm


function select_interaction_fastkernel(result, fk::FastKernelGeometry)
    idx = get(fk.name_to_index, result.region, 0)
    idx == 0 && return (
        class = :other,
        sensitive = false,
        passes_threshold = false,
        ecut_keV = 0.0,
        region = "MARS",
        interaction_type = result.interaction_type
    )

    region = fk.regions[idx]
    class = if region.name == "FV"
        :fv
    elseif region.name in ("TopActive", "BarrelActive", "BottomActive")
        :tpc
    elseif region.name == "Skin"
        :skin
    else
        :other
    end

    dep_keV = 1000.0 * result.deposit_E_MeV
    passes_threshold = region.sensitive && dep_keV >= region.ecut_keV

    (
        class = class,
        sensitive = region.sensitive,
        passes_threshold = passes_threshold,
        ecut_keV = region.ecut_keV,
        region = region.name,
        interaction_type = result.interaction_type
    )
end




function veto_threshold(tag::RegionTag, cfg::SimConfig)::Float64
    tag == TAG_TPC_ACTIVE && return cfg.veto_TPC
    tag == TAG_SKIN && return cfg.veto_skin
    tag == TAG_FV && return 0.0
    tag == TAG_PASSIVE_LXE && return Inf
    Inf
end


veto_threshold(node::DetectorNode, cfg::SimConfig)::Float64 = veto_threshold(region_tag(node), cfg)


function _parse_region_tag(raw)::RegionTag
    s = lowercase(String(raw))
    s == "world" && return TAG_WORLD
    s == "vacuum" && return TAG_VACUUM
    s == "structural" && return TAG_STRUCTURAL
    s == "tpc_active" && return TAG_TPC_ACTIVE
    s == "fv" && return TAG_FV
    s == "skin" && return TAG_SKIN
    s == "passive_lxe" && return TAG_PASSIVE_LXE
    error("Unknown RegionTag '$raw'")
end


function _parse_activity(d)::Dict{String,Float64}
    activity = Dict{String,Float64}()
    raw = get(d, "activity", nothing)
    raw === nothing && return activity
    for (k, v) in raw
        activity[String(k)] = Float64(v)
    end
    activity
end


function _build_tracking_solid(name::String, d)
    shape = lowercase(String(d["shape"]))
    if shape == "cylinder"
        return Cyl(Float64(d["radius_cm"]), Float64(d["half_height_cm"]))
    elseif shape == "cylinder_shell"
        return CylShell(Float64(d["R_inner_cm"]),
                        Float64(d["wall_thickness_cm"]),
                        Float64(d["half_height_cm"]))
    elseif shape == "box"
        return Box(Float64(d["half_x_cm"]),
                   Float64(d["half_y_cm"]),
                   Float64(d["half_z_cm"]))
    elseif shape == "disk"
        return Disk(Float64(d["radius_cm"]),
                    Float64(d["wall_thickness_cm"]),
                    Float64(get(d, "aspect_ratio", Inf)))
    elseif shape == "cap"
        return Cap(Float64(d["radius_cm"]),
                   Float64(d["aspect_ratio"]))
    elseif shape == "domed_container"
        radius = Float64(d["radius_cm"])
        return DomedContainer(
            radius,
            Float64(d["barrel_half_height_cm"]),
            Cap(radius, Float64(d["top_aspect_ratio"])),
            Cap(radius, Float64(d["bottom_aspect_ratio"]))
        )
    elseif shape == "capped_cylinder"
        radius = Float64(d["radius_cm"])
        top_cap = haskey(d, "top_aspect_ratio") ? Cap(radius, Float64(d["top_aspect_ratio"])) : nothing
        bottom_cap = haskey(d, "bottom_aspect_ratio") ? Cap(radius, Float64(d["bottom_aspect_ratio"])) : nothing
        return CappedCylinder(
            radius,
            Float64(d["barrel_half_height_cm"]);
            top_cap=top_cap,
            bottom_cap=bottom_cap
        )
    end
    error("Unknown shape '$shape' for volume '$name'")
end


function _default_tracking_semantics(name::String)
    if name == "MARS"
        return (tag=TAG_WORLD, role="world", inFastKernel=false, isXe=false, sensitive=false, ecut_keV=0.0, dz_mm=0.0, fv_target=false)
    elseif name in ("LZ_detector", "AirDome")
        return (tag=TAG_VACUUM, role=lowercase(name), inFastKernel=true, isXe=false, sensitive=false, ecut_keV=0.0, dz_mm=0.0, fv_target=false)
    elseif name == "LXe_passive"
        return (tag=TAG_PASSIVE_LXE, role="lxe_passive", inFastKernel=true, isXe=true, sensitive=false, ecut_keV=0.0, dz_mm=0.0, fv_target=false)
    elseif name in ("TopActive", "BarrelActive", "BottomActive")
        return (tag=TAG_TPC_ACTIVE, role=lowercase(name), inFastKernel=true, isXe=true, sensitive=true, ecut_keV=10.0, dz_mm=0.0, fv_target=false)
    elseif name == "FV"
        return (tag=TAG_FV, role="fiducial_region", inFastKernel=true, isXe=true, sensitive=true, ecut_keV=10.0, dz_mm=3.0, fv_target=true)
    elseif name == "Skin"
        return (tag=TAG_SKIN, role="veto_region", inFastKernel=true, isXe=true, sensitive=true, ecut_keV=100.0, dz_mm=0.0, fv_target=false)
    elseif name in ("FC_PTFE", "FC_rings")
        return (tag=TAG_STRUCTURAL, role=lowercase(name), inFastKernel=true, isXe=false, sensitive=false, ecut_keV=0.0, dz_mm=0.0, fv_target=false)
    end
    return (tag=TAG_STRUCTURAL, role=lowercase(name), inFastKernel=false, isXe=false, sensitive=false, ecut_keV=0.0, dz_mm=0.0, fv_target=false)
end


function _build_tracking_node(id::Int, d, material::Material; version::Int=2)::Tuple{DetectorNode,String}
    name = String(d["name"])
    solid = _build_tracking_solid(name, d)
    defaults = version >= 3 ? _default_tracking_semantics(name) : nothing
    tag = haskey(d, "tag") ? _parse_region_tag(d["tag"]) : defaults.tag
    role = String(get(d, "role", defaults === nothing ? "" : defaults.role))
    inFastKernel = Bool(get(d, "inFastKernel", defaults === nothing ? false : defaults.inFastKernel))
    isXe = Bool(get(d, "isXe", defaults === nothing ? false : defaults.isXe))
    sensitive = Bool(get(d, "sensitive", defaults === nothing ? false : defaults.sensitive))
    ecut_keV = Float64(get(d, "ecut_keV", defaults === nothing ? 0.0 : defaults.ecut_keV))
    dz_mm = Float64(get(d, "dz_mm", defaults === nothing ? 0.0 : defaults.dz_mm))
    fv_target = Bool(get(d, "fv_target", defaults === nothing ? false : defaults.fv_target))
    approximation = String(get(d, "approximation", "exact"))
    activity = _parse_activity(d)
    lv = LogicalVolume(name, solid, material, tag, role, inFastKernel, isXe, sensitive, ecut_keV, dz_mm, fv_target, approximation, activity)

    pos = Float64.(get(d, "position_cm", [0.0, 0.0, 0.0]))
    orientation = Symbol(get(d, "orientation", "none"))
    placement = Placement(pos, orientation)

    parent_name = String(get(d, "parent", ""))
    node = DetectorNode(id, 0, Int[], lv, placement)
    (node, parent_name)
end


function _logical_volume(node::DetectorNode)
    solid = node.lv.solid
    pos = node.placement.position_cm
    if solid isa Cyl
        return LCyl(solid, pos)
    elseif solid isa CylShell
        return LCylShell(solid, pos)
    elseif solid isa Box
        return LBox(solid, pos)
    elseif solid isa Disk
        return LDisk(solid, pos, node.placement.orientation)
    elseif solid isa Cap || solid isa DomedContainer || solid isa CappedCylinder
        return nothing
    end
    error("Unsupported solid type for node '$(node.lv.name)'")
end


function _is_inside_cap(cap::Cap, placement::Placement, pos::Vector{Float64})::Bool
    _is_inside_cap(cap,
                   placement.position_cm[1], placement.position_cm[2], placement.position_cm[3],
                   placement.orientation,
                   pos[1], pos[2], pos[3])
end


function _is_inside_cap(cap::Cap,
                        cx::Float64, cy::Float64, cz::Float64,
                        orientation::Symbol,
                        x::Float64, y::Float64, z::Float64)::Bool
    dx = x - cx
    dy = y - cy
    dz = z - cz
    r2 = dx^2 + dy^2
    R = cap.radius_cm
    r2 > R^2 && return false

    sgn = orientation === :down ? -1.0 : 1.0
    dz_rel = sgn * dz
    dz_rel < 0.0 && return false

    c = depth(cap)
    r_frac2 = r2 / R^2
    r_frac2 + (dz_rel / c)^2 <= 1.0
end


function _is_inside_domed_container(dc::DomedContainer, placement::Placement, pos::Vector{Float64})::Bool
    dx = pos[1] - placement.position_cm[1]
    dy = pos[2] - placement.position_cm[2]
    dz = pos[3] - placement.position_cm[3]
    r2 = dx^2 + dy^2
    R = dc.radius_cm
    r2 > R^2 && return false

    H = dc.barrel_half_height_cm
    if abs(dz) <= H
        return true
    elseif dz > H
        return _is_inside_cap(dc.top_cap,
                              placement.position_cm[1], placement.position_cm[2], placement.position_cm[3] + H,
                              :up,
                              pos[1], pos[2], pos[3])
    else
        return _is_inside_cap(dc.bottom_cap,
                              placement.position_cm[1], placement.position_cm[2], placement.position_cm[3] - H,
                              :down,
                              pos[1], pos[2], pos[3])
    end
end


function _is_inside_capped_cylinder(cc::CappedCylinder, placement::Placement, pos::Vector{Float64})::Bool
    dx = pos[1] - placement.position_cm[1]
    dy = pos[2] - placement.position_cm[2]
    dz = pos[3] - placement.position_cm[3]
    r2 = dx^2 + dy^2
    R = cc.radius_cm
    r2 > R^2 && return false

    H = cc.barrel_half_height_cm
    if abs(dz) <= H
        return true
    elseif dz > H
        cc.top_cap === nothing && return false
        return _is_inside_cap(cc.top_cap,
                              placement.position_cm[1], placement.position_cm[2], placement.position_cm[3] + H,
                              :up,
                              pos[1], pos[2], pos[3])
    else
        cc.bottom_cap === nothing && return false
        return _is_inside_cap(cc.bottom_cap,
                              placement.position_cm[1], placement.position_cm[2], placement.position_cm[3] - H,
                              :down,
                              pos[1], pos[2], pos[3])
    end
end


function is_inside(node::DetectorNode, pos::Vector{Float64})
    solid = node.lv.solid
    cx = node.placement.position_cm[1]
    cy = node.placement.position_cm[2]
    cz = node.placement.position_cm[3]
    x = pos[1]
    y = pos[2]
    z = pos[3]

    if solid isa Cap
        return _is_inside_cap(solid, cx, cy, cz, node.placement.orientation, x, y, z)
    elseif solid isa DomedContainer
        return _is_inside_domed_container(solid, node.placement, pos)
    elseif solid isa CappedCylinder
        return _is_inside_capped_cylinder(solid, node.placement, pos)
    elseif solid isa Cyl
        dx = x - cx
        dy = y - cy
        dz = z - cz
        return dx^2 + dy^2 < solid.radius_cm^2 && abs(dz) < solid.half_height_cm
    elseif solid isa CylShell
        dx = x - cx
        dy = y - cy
        dz = z - cz
        r2 = dx^2 + dy^2
        return r2 >= solid.R_inner_cm^2 && r2 < R_outer(solid)^2 && abs(dz) < solid.half_height_cm
    elseif solid isa Box
        return abs(x - cx) < solid.half_x_cm &&
               abs(y - cy) < solid.half_y_cm &&
               abs(z - cz) < solid.half_z_cm
    elseif solid isa Disk
        dx = x - cx
        dy = y - cy
        dz = z - cz
        r2 = dx^2 + dy^2
        R = solid.radius_cm
        r2 > R^2 && return false
        if is_flat(solid)
            sgn = node.placement.orientation === :up ? 1.0 : -1.0
            return 0.0 <= sgn * dz <= solid.wall_thickness_cm
        end

        c_inner = depth(solid)
        c_outer = c_inner + solid.wall_thickness_cm
        sgn = node.placement.orientation === :up ? 1.0 : -1.0
        dz_rel = sgn * dz
        dz_rel < 0.0 && return false

        r_frac2 = r2 / R^2
        inside_outer = r_frac2 + (dz_rel / c_outer)^2 <= 1.0
        outside_inner = r_frac2 + (dz_rel / c_inner)^2 >= 1.0
        return inside_outer && outside_inner
    end
    is_inside(_logical_volume(node), pos)
end


is_inside(node::DetectorNode, pos::NTuple{3,<:Real}) = is_inside(node, Float64[pos...])


function _sample_points(node::DetectorNode)::Vector{Vector{Float64}}
    solid = node.lv.solid
    cx, cy, cz = node.placement.position_cm

    if solid isa Cyl
        R = solid.radius_cm
        H = solid.half_height_cm
        return [
            Float64[cx, cy, cz],
            Float64[cx + 0.5R, cy, cz],
            Float64[cx - 0.5R, cy, cz],
            Float64[cx, cy + 0.5R, cz],
            Float64[cx, cy, cz + 0.999H],
            Float64[cx, cy, cz - 0.999H],
            Float64[cx + 0.999R, cy, cz]
        ]
    elseif solid isa CylShell
        Rmid = solid.R_inner_cm + 0.5 * solid.wall_thickness_cm
        Rin = solid.R_inner_cm + 1e-4
        Rout = R_outer(solid) - 1e-4
        H = solid.half_height_cm
        return [
            Float64[cx + Rmid, cy, cz],
            Float64[cx - Rmid, cy, cz],
            Float64[cx, cy + Rmid, cz],
            Float64[cx + Rin, cy, cz],
            Float64[cx + Rout, cy, cz],
            Float64[cx + Rmid, cy, cz + 0.999H],
            Float64[cx + Rmid, cy, cz - 0.999H]
        ]
    elseif solid isa Box
        x = 0.999 * solid.half_x_cm
        y = 0.999 * solid.half_y_cm
        z = 0.999 * solid.half_z_cm
        return [
            Float64[cx, cy, cz],
            Float64[cx + x, cy, cz],
            Float64[cx - x, cy, cz],
            Float64[cx, cy + y, cz],
            Float64[cx, cy, cz + z],
            Float64[cx, cy, cz - z]
        ]
    elseif solid isa Disk
        sgn = node.placement.orientation === :down ? -1.0 : 1.0
        if is_flat(solid)
            zmid = cz + sgn * 0.5 * solid.wall_thickness_cm
            return [
                Float64[cx, cy, zmid],
                Float64[cx + 0.5 * solid.radius_cm, cy, zmid]
            ]
        end
        z_apex = cz + sgn * (depth(solid) + 0.5 * solid.wall_thickness_cm)
        z_mid = cz + sgn * (0.5 * depth(solid) + 0.5 * solid.wall_thickness_cm)
        return [
            Float64[cx, cy, z_apex],
            Float64[cx, cy, z_mid],
            Float64[cx + 0.5 * solid.radius_cm, cy, z_mid]
        ]
    elseif solid isa Cap
        sgn = node.placement.orientation === :down ? -1.0 : 1.0
        c = depth(solid)
        z_apex = cz + sgn * (0.999 * c)
        z_mid = cz + sgn * (0.5 * c)
        return [
            Float64[cx, cy, cz],
            Float64[cx, cy, z_mid],
            Float64[cx, cy, z_apex],
            Float64[cx + 0.5 * solid.radius_cm, cy, z_mid]
        ]
    elseif solid isa DomedContainer
        H = solid.barrel_half_height_cm
        dt = depth(solid.top_cap)
        db = depth(solid.bottom_cap)
        return [
            Float64[cx, cy, cz],
            Float64[cx + 0.5 * solid.radius_cm, cy, cz],
            Float64[cx, cy, cz + 0.999 * H],
            Float64[cx, cy, cz - 0.999 * H],
            Float64[cx, cy, cz + H + 0.5 * dt],
            Float64[cx, cy, cz - H - 0.5 * db],
            Float64[cx + 0.5 * solid.radius_cm, cy, cz + H + 0.5 * dt],
            Float64[cx + 0.5 * solid.radius_cm, cy, cz - H - 0.5 * db]
        ]
    elseif solid isa CappedCylinder
        H = solid.barrel_half_height_cm
        pts = Vector{Float64}[
            Float64[cx, cy, cz],
            Float64[cx + 0.5 * solid.radius_cm, cy, cz],
            Float64[cx, cy, cz + 0.999 * H],
            Float64[cx, cy, cz - 0.999 * H]
        ]
        if solid.top_cap !== nothing
            dt = depth(solid.top_cap)
            push!(pts, Float64[cx, cy, cz + H + 0.5 * dt])
            push!(pts, Float64[cx + 0.5 * solid.radius_cm, cy, cz + H + 0.5 * dt])
        end
        if solid.bottom_cap !== nothing
            db = depth(solid.bottom_cap)
            push!(pts, Float64[cx, cy, cz - H - 0.5 * db])
            push!(pts, Float64[cx + 0.5 * solid.radius_cm, cy, cz - H - 0.5 * db])
        end
        return pts
    end
    error("Unsupported solid type for node '$(node.lv.name)'")
end


function _interior_sample_points(node::DetectorNode; tol_cm::Float64=GEOM_VALIDATE_TOL_CM)::Vector{Vector{Float64}}
    solid = node.lv.solid
    cx, cy, cz = node.placement.position_cm

    if solid isa Cyl
        R = max(solid.radius_cm - tol_cm, 0.5 * solid.radius_cm)
        H = max(solid.half_height_cm - tol_cm, 0.5 * solid.half_height_cm)
        return [
            Float64[cx, cy, cz],
            Float64[cx + 0.5R, cy, cz],
            Float64[cx - 0.5R, cy, cz],
            Float64[cx, cy + 0.5R, cz],
            Float64[cx, cy, cz + H],
            Float64[cx, cy, cz - H]
        ]
    elseif solid isa CylShell
        Rin = solid.R_inner_cm + tol_cm
        Rout = R_outer(solid) - tol_cm
        Rin < Rout || return Float64[Float64[cx + solid.R_inner_cm + 0.5 * solid.wall_thickness_cm, cy, cz]]
        Rmid = 0.5 * (Rin + Rout)
        H = max(solid.half_height_cm - tol_cm, 0.5 * solid.half_height_cm)
        return [
            Float64[cx + Rmid, cy, cz],
            Float64[cx - Rmid, cy, cz],
            Float64[cx, cy + Rmid, cz],
            Float64[cx + Rin, cy, cz],
            Float64[cx + Rout, cy, cz],
            Float64[cx + Rmid, cy, cz + H],
            Float64[cx + Rmid, cy, cz - H]
        ]
    elseif solid isa Box
        x = max(solid.half_x_cm - tol_cm, 0.5 * solid.half_x_cm)
        y = max(solid.half_y_cm - tol_cm, 0.5 * solid.half_y_cm)
        z = max(solid.half_z_cm - tol_cm, 0.5 * solid.half_z_cm)
        return [
            Float64[cx, cy, cz],
            Float64[cx + x, cy, cz],
            Float64[cx - x, cy, cz],
            Float64[cx, cy + y, cz],
            Float64[cx, cy, cz + z],
            Float64[cx, cy, cz - z]
        ]
    elseif solid isa Disk
        sgn = node.placement.orientation === :down ? -1.0 : 1.0
        if is_flat(solid)
            tmid = 0.5 * max(solid.wall_thickness_cm - tol_cm, 0.5 * solid.wall_thickness_cm)
            zmid = cz + sgn * tmid
            return [
                Float64[cx, cy, zmid],
                Float64[cx + 0.5 * max(solid.radius_cm - tol_cm, 0.5 * solid.radius_cm), cy, zmid]
            ]
        end
        c = depth(solid)
        z_apex = cz + sgn * max(c - tol_cm, 0.5 * c)
        z_mid = cz + sgn * 0.5 * c
        return [
            Float64[cx, cy, z_mid],
            Float64[cx, cy, z_apex],
            Float64[cx + 0.5 * max(solid.radius_cm - tol_cm, 0.5 * solid.radius_cm), cy, z_mid]
        ]
    elseif solid isa Cap
        sgn = node.placement.orientation === :down ? -1.0 : 1.0
        c = depth(solid)
        z_apex = cz + sgn * max(c - tol_cm, 0.5 * c)
        z_mid = cz + sgn * 0.5 * c
        return [
            Float64[cx, cy, z_mid],
            Float64[cx, cy, z_apex],
            Float64[cx + 0.5 * max(solid.radius_cm - tol_cm, 0.5 * solid.radius_cm), cy, z_mid]
        ]
    elseif solid isa DomedContainer
        H = max(solid.barrel_half_height_cm - tol_cm, 0.5 * solid.barrel_half_height_cm)
        dt = depth(solid.top_cap)
        db = depth(solid.bottom_cap)
        R = max(solid.radius_cm - tol_cm, 0.5 * solid.radius_cm)
        return [
            Float64[cx, cy, cz],
            Float64[cx + 0.5 * R, cy, cz],
            Float64[cx, cy, cz + H],
            Float64[cx, cy, cz - H],
            Float64[cx, cy, cz + solid.barrel_half_height_cm + 0.5 * max(dt - tol_cm, 0.5 * dt)],
            Float64[cx, cy, cz - solid.barrel_half_height_cm - 0.5 * max(db - tol_cm, 0.5 * db)]
        ]
    elseif solid isa CappedCylinder
        H = max(solid.barrel_half_height_cm - tol_cm, 0.5 * solid.barrel_half_height_cm)
        R = max(solid.radius_cm - tol_cm, 0.5 * solid.radius_cm)
        pts = Vector{Float64}[
            Float64[cx, cy, cz],
            Float64[cx + 0.5 * R, cy, cz],
            Float64[cx, cy, cz + H],
            Float64[cx, cy, cz - H]
        ]
        if solid.top_cap !== nothing
            dt = depth(solid.top_cap)
            push!(pts, Float64[cx, cy, cz + solid.barrel_half_height_cm + 0.5 * max(dt - tol_cm, 0.5 * dt)])
        end
        if solid.bottom_cap !== nothing
            db = depth(solid.bottom_cap)
            push!(pts, Float64[cx, cy, cz - solid.barrel_half_height_cm - 0.5 * max(db - tol_cm, 0.5 * db)])
        end
        return pts
    end
    error("Unsupported solid type for node '$(node.lv.name)'")
end


function _validate_reachability(det::TrackingDetector)
    seen = falses(length(det.nodes))
    stack = Int[det.root_id]
    while !isempty(stack)
        id = pop!(stack)
        seen[id] && continue
        seen[id] = true
        append!(stack, det.nodes[id].children_ids)
    end
    all(seen) || error("TrackingDetector contains unreachable nodes")
end


function _validate_child_containment(det::TrackingDetector)
    for node in det.nodes
        node.parent_id == 0 && continue
        parent = det.nodes[node.parent_id]
        parent.lv.approximation == "exact" || continue
        for p in _sample_points(node)
            is_inside(parent, p) || error("Node '$(node.lv.name)' is not contained in parent '$(parent.lv.name)'")
        end
    end
end


function _has_obvious_overlap(a::DetectorNode, b::DetectorNode)::Bool
    for p in _interior_sample_points(a)
        is_inside(b, p) && return true
    end
    for p in _interior_sample_points(b)
        is_inside(a, p) && return true
    end
    false
end


function _validate_sibling_overlaps(det::TrackingDetector)
    for parent in det.nodes
        ids = parent.children_ids
        for i in 1:length(ids)-1
            a = det.nodes[ids[i]]
            for j in i+1:length(ids)
                b = det.nodes[ids[j]]
                (a.lv.approximation == "exact" && b.lv.approximation == "exact") || continue
                _has_obvious_overlap(a, b) && error("Sibling nodes '$(a.lv.name)' and '$(b.lv.name)' overlap inside parent '$(parent.lv.name)'")
            end
        end
    end
end


function cap_volume(cap::Cap)::Float64
    a = cap.radius_cm
    c = depth(cap)
    2.0 / 3.0 * π * a^2 * c
end


function volume(dc::DomedContainer)::Float64
    2.0 * dc.barrel_half_height_cm * π * dc.radius_cm^2 +
    cap_volume(dc.top_cap) + cap_volume(dc.bottom_cap)
end


function volume(cc::CappedCylinder)::Float64
    V = 2.0 * cc.barrel_half_height_cm * π * cc.radius_cm^2
    cc.top_cap !== nothing && (V += cap_volume(cc.top_cap))
    cc.bottom_cap !== nothing && (V += cap_volume(cc.bottom_cap))
    V
end


function volume(node::DetectorNode)::Float64
    solid = node.lv.solid
    if solid isa Cap
        return cap_volume(solid)
    elseif solid isa DomedContainer
        return volume(solid)
    elseif solid isa CappedCylinder
        return volume(solid)
    end
    return volume(solid)
end


function sibling_overlaps(det::TrackingDetector; exact_only::Bool=true)
    overlaps = NamedTuple[]
    for parent in det.nodes
        ids = parent.children_ids
        for i in 1:length(ids)-1
            a = det.nodes[ids[i]]
            for j in i+1:length(ids)
                b = det.nodes[ids[j]]
                exact_only && !(a.lv.approximation == "exact" && b.lv.approximation == "exact") && continue
                _has_obvious_overlap(a, b) || continue
                push!(overlaps, (
                    parent = parent.lv.name,
                    a = a.lv.name,
                    b = b.lv.name,
                    a_approximation = a.lv.approximation,
                    b_approximation = b.lv.approximation
                ))
            end
        end
    end
    overlaps
end


function validate_tracking_detector(det::TrackingDetector)::Bool
    root = root_node(det)
    root.parent_id == 0 || error("Root node '$(root.lv.name)' must have parent_id = 0")
    root.lv.tag == TAG_WORLD || error("Root node '$(root.lv.name)' must have tag TAG_WORLD")
    _validate_reachability(det)
    _validate_child_containment(det)
    _validate_sibling_overlaps(det)
    true
end


function tree_dump(det::TrackingDetector)::String
    lines = String[]
    function visit(id::Int, depth::Int)
        node = det.nodes[id]
        indent = "  "^depth
        push!(lines, "$(indent)- $(node.lv.name) [$(node.lv.role), $(node.lv.tag)]")
        for child_id in node.children_ids
            visit(child_id, depth + 1)
        end
    end
    visit(det.root_id, 0)
    join(lines, "\n")
end


"""
    _load_tracking_tree(path, materials) -> TrackingDetector

Load the canonical tracking detector JSON. The file defines a unique
`world` entry and a flat list of volumes with explicit `parent` links.
"""
function _load_tracking_tree(path::AbstractString,
                             materials::Dict{String,Material};
                             validate::Bool=true)::TrackingDetector
    raw = open(path, "r") do io
        JSON.parse(io)
    end

    det_name = String(raw["name"])
    version = Int(get(raw, "version", 2))

    haskey(raw, "world") || error("Tracking geometry file '$path' is missing 'world'")
    haskey(raw, "volumes") || error("Tracking geometry file '$path' is missing 'volumes'")

    nodes = DetectorNode[]
    name_to_id = Dict{String,Int}()
    parent_names = String[]

    world_raw = raw["world"]
    world_name = String(world_raw["name"])
    world_mat_name = String(world_raw["material"])
    haskey(materials, world_mat_name) || error("Unknown material '$world_mat_name' for world '$world_name'")
    world_node, world_parent = _build_tracking_node(1, world_raw, materials[world_mat_name]; version=version)
    world_node.lv.tag == TAG_WORLD || error("World node '$world_name' must have tag 'world'")
    isempty(world_parent) || error("World node '$world_name' cannot define a parent")
    push!(nodes, world_node)
    name_to_id[world_name] = 1
    push!(parent_names, "")

    for d in raw["volumes"]
        name = String(d["name"])
        haskey(name_to_id, name) && error("Duplicate detector node name '$name'")

        mat_name = String(d["material"])
        haskey(materials, mat_name) || error("Unknown material '$mat_name' for node '$name'")

        id = length(nodes) + 1
        node, parent_name = _build_tracking_node(id, d, materials[mat_name]; version=version)
        isempty(parent_name) && error("Tracking geometry node '$name' is missing 'parent'")

        push!(nodes, node)
        name_to_id[name] = id
        push!(parent_names, parent_name)
    end

    for node in nodes[2:end]
        parent_name = parent_names[node.id]
        parent_id = get(name_to_id, parent_name, 0)
        parent_id == 0 && error("Tracking geometry node '$(node.lv.name)' references unknown parent '$parent_name'")
        node.parent_id = parent_id
        push!(nodes[parent_id].children_ids, node.id)
    end

    det = TrackingDetector(det_name, nodes, name_to_id, 1)
    validate && validate_tracking_detector(det)
    det
end


"""
    load_tracking_detector(path, materials; validate=true) -> TrackingDetector

Load the canonical tracking detector geometry. This expects the explicit
partitioned tracking schema used by the current tracking detector file.
"""
function load_tracking_detector(path::AbstractString,
                               materials::Dict{String,Material};
                               validate::Bool=true)::TrackingDetector
    _load_tracking_tree(path, materials; validate=validate)
end


root_node(det::TrackingDetector) = det.nodes[det.root_id]


function node_by_name(det::TrackingDetector, name::AbstractString)::DetectorNode
    id = get(det.name_to_id, String(name), 0)
    id == 0 && error("TrackingDetector has no node named '$name'")
    det.nodes[id]
end


function child_nodes(det::TrackingDetector, node::DetectorNode)::Vector{DetectorNode}
    [det.nodes[id] for id in node.children_ids]
end


function child_nodes(det::TrackingDetector, name::AbstractString)::Vector{DetectorNode}
    child_nodes(det, node_by_name(det, name))
end


function detector_summary(det::TrackingDetector)::String
    root = root_node(det)
    "TrackingDetector($(det.name)): $(length(det.nodes)) nodes, root=$(root.lv.name)"
end


function _compile_fastkernel_region(id::Int, node::DetectorNode)::FastKernelRegion
    solid = node.lv.solid
    x0, y0, z0 = node.placement.position_cm
    abs(x0) <= GEOM_VALIDATE_TOL_CM || error("FastKernel currently expects coaxial geometry; node '$(node.lv.name)' has x=$(x0)")
    abs(y0) <= GEOM_VALIDATE_TOL_CM || error("FastKernel currently expects coaxial geometry; node '$(node.lv.name)' has y=$(y0)")

    kind = :unknown
    rmin_cm = 0.0
    rmax_cm = 0.0
    zmin_cm = 0.0
    zmax_cm = 0.0
    has_top_cap = false
    has_bottom_cap = false
    top_cap_radius_cm = 0.0
    top_cap_aspect_ratio = 0.0
    bottom_cap_radius_cm = 0.0
    bottom_cap_aspect_ratio = 0.0

    if solid isa Cyl
        kind = :cylinder
        rmax_cm = solid.radius_cm
        zmin_cm = z0 - solid.half_height_cm
        zmax_cm = z0 + solid.half_height_cm
    elseif solid isa CylShell
        kind = :cylinder_shell
        rmin_cm = solid.R_inner_cm
        rmax_cm = solid.R_inner_cm + solid.wall_thickness_cm
        zmin_cm = z0 - solid.half_height_cm
        zmax_cm = z0 + solid.half_height_cm
    elseif solid isa Cap
        kind = :cap
        cap_depth = depth(solid)
        rmax_cm = solid.radius_cm
        if node.placement.orientation === :down
            has_bottom_cap = true
            bottom_cap_radius_cm = solid.radius_cm
            bottom_cap_aspect_ratio = solid.aspect_ratio
            zmin_cm = z0 - cap_depth
            zmax_cm = z0
        else
            has_top_cap = true
            top_cap_radius_cm = solid.radius_cm
            top_cap_aspect_ratio = solid.aspect_ratio
            zmin_cm = z0
            zmax_cm = z0 + cap_depth
        end
    elseif solid isa DomedContainer
        kind = :domed_container
        rmax_cm = solid.radius_cm
        zmin_cm = z0 - solid.barrel_half_height_cm - depth(solid.bottom_cap)
        zmax_cm = z0 + solid.barrel_half_height_cm + depth(solid.top_cap)
        has_top_cap = true
        has_bottom_cap = true
        top_cap_radius_cm = solid.top_cap.radius_cm
        top_cap_aspect_ratio = solid.top_cap.aspect_ratio
        bottom_cap_radius_cm = solid.bottom_cap.radius_cm
        bottom_cap_aspect_ratio = solid.bottom_cap.aspect_ratio
    elseif solid isa CappedCylinder
        kind = :capped_cylinder
        rmax_cm = solid.radius_cm
        zmin_cm = z0 - solid.barrel_half_height_cm
        zmax_cm = z0 + solid.barrel_half_height_cm
        if solid.top_cap !== nothing
            has_top_cap = true
            top_cap_radius_cm = solid.top_cap.radius_cm
            top_cap_aspect_ratio = solid.top_cap.aspect_ratio
            zmax_cm += depth(solid.top_cap)
        end
        if solid.bottom_cap !== nothing
            has_bottom_cap = true
            bottom_cap_radius_cm = solid.bottom_cap.radius_cm
            bottom_cap_aspect_ratio = solid.bottom_cap.aspect_ratio
            zmin_cm -= depth(solid.bottom_cap)
        end
    else
        error("Unsupported fast-kernel solid type '$(typeof(solid))' for node '$(node.lv.name)'")
    end

    FastKernelRegion(
        id,
        node.lv.name,
        node.lv.role,
        kind,
        node.lv.material,
        node.lv.isXe,
        node.lv.sensitive,
        node.lv.ecut_keV,
        node.lv.fv_target,
        rmin_cm,
        rmax_cm,
        zmin_cm,
        zmax_cm,
        has_top_cap,
        has_bottom_cap,
        top_cap_radius_cm,
        top_cap_aspect_ratio,
        bottom_cap_radius_cm,
        bottom_cap_aspect_ratio
    )
end


function _with_fastkernel_identity(id::Int, name::String, role::String,
                                   region::FastKernelRegion;
                                   kind::Symbol=region.kind,
                                   sensitive::Bool=region.sensitive,
                                   ecut_keV::Float64=region.ecut_keV,
                                   fv_target::Bool=region.fv_target,
                                   rmin_cm::Float64=region.rmin_cm,
                                   rmax_cm::Float64=region.rmax_cm,
                                   zmin_cm::Float64=region.zmin_cm,
                                   zmax_cm::Float64=region.zmax_cm,
                                   has_top_cap::Bool=region.has_top_cap,
                                   has_bottom_cap::Bool=region.has_bottom_cap,
                                   top_cap_radius_cm::Float64=region.top_cap_radius_cm,
                                   top_cap_aspect_ratio::Float64=region.top_cap_aspect_ratio,
                                   bottom_cap_radius_cm::Float64=region.bottom_cap_radius_cm,
                                   bottom_cap_aspect_ratio::Float64=region.bottom_cap_aspect_ratio)
    FastKernelRegion(
        id,
        name,
        role,
        kind,
        region.material,
        region.isXe,
        sensitive,
        ecut_keV,
        fv_target,
        rmin_cm,
        rmax_cm,
        zmin_cm,
        zmax_cm,
        has_top_cap,
        has_bottom_cap,
        top_cap_radius_cm,
        top_cap_aspect_ratio,
        bottom_cap_radius_cm,
        bottom_cap_aspect_ratio
    )
end


function _push_fastkernel_region!(regions::Vector{FastKernelRegion},
                                  name_to_index::Dict{String,Int},
                                  region::FastKernelRegion)
    push!(regions, region)
    name_to_index[region.name] = region.id
    region
end


function compile_fastkernel_geometry(det::TrackingDetector)::FastKernelGeometry
    raw_regions = Dict{String,FastKernelRegion}()

    for node in det.nodes
        node.lv.inFastKernel || continue
        region = _compile_fastkernel_region(length(raw_regions) + 1, node)
        raw_regions[region.name] = region
    end

    regions = FastKernelRegion[]
    name_to_index = Dict{String,Int}()
    for (idx, name) in enumerate(("LZ_detector", "AirDome", "FC_PTFE", "FC_rings", "FV", "TopActive", "BarrelActive", "BottomActive", "Skin", "LXe_passive"))
        haskey(raw_regions, name) || error("FastKernelGeometry requires canonical tracking region '$name'")
        region = raw_regions[name]
        _push_fastkernel_region!(regions, name_to_index,
            _with_fastkernel_identity(idx, name, region.role, region))
    end

    FastKernelGeometry(
        regions,
        name_to_index,
        name_to_index["LZ_detector"],
        name_to_index["AirDome"],
        name_to_index["FC_PTFE"],
        name_to_index["FC_rings"],
        name_to_index["FV"],
        name_to_index["TopActive"],
        name_to_index["BarrelActive"],
        name_to_index["BottomActive"],
        name_to_index["Skin"],
        name_to_index["LXe_passive"]
    )
end


function _is_inside_fastkernel_region(region::FastKernelRegion,
                                      x::Float64, y::Float64, z::Float64)::Bool
    r2 = x^2 + y^2
    if region.kind == :cylinder
        return r2 < region.rmax_cm^2 &&
               z > region.zmin_cm &&
               z < region.zmax_cm
    elseif region.kind == :cylinder_shell
        return r2 >= region.rmin_cm^2 &&
               r2 < region.rmax_cm^2 &&
               z > region.zmin_cm &&
               z < region.zmax_cm
    elseif region.kind == :cap
        r2 > region.rmax_cm^2 && return false
        if region.has_top_cap
            z < region.zmin_cm && return false
            dz_rel = z - region.zmin_cm
            c = region.top_cap_radius_cm / region.top_cap_aspect_ratio
            return r2 / region.top_cap_radius_cm^2 + (dz_rel / c)^2 <= 1.0
        else
            z > region.zmax_cm && return false
            dz_rel = region.zmax_cm - z
            c = region.bottom_cap_radius_cm / region.bottom_cap_aspect_ratio
            return r2 / region.bottom_cap_radius_cm^2 + (dz_rel / c)^2 <= 1.0
        end
    elseif region.kind == :domed_container
        r2 > region.rmax_cm^2 && return false
        top_depth = region.has_top_cap ? region.top_cap_radius_cm / region.top_cap_aspect_ratio : 0.0
        bottom_depth = region.has_bottom_cap ? region.bottom_cap_radius_cm / region.bottom_cap_aspect_ratio : 0.0
        barrel_zmin = region.zmin_cm + bottom_depth
        barrel_zmax = region.zmax_cm - top_depth
        if z >= barrel_zmin && z <= barrel_zmax
            return true
        elseif z > barrel_zmax
            dz_rel = z - barrel_zmax
            c = region.top_cap_radius_cm / region.top_cap_aspect_ratio
            return r2 / region.top_cap_radius_cm^2 + (dz_rel / c)^2 <= 1.0
        else
            dz_rel = barrel_zmin - z
            c = region.bottom_cap_radius_cm / region.bottom_cap_aspect_ratio
            return r2 / region.bottom_cap_radius_cm^2 + (dz_rel / c)^2 <= 1.0
        end
    elseif region.kind == :capped_cylinder
        r2 > region.rmax_cm^2 && return false
        top_depth = region.has_top_cap ? region.top_cap_radius_cm / region.top_cap_aspect_ratio : 0.0
        bottom_depth = region.has_bottom_cap ? region.bottom_cap_radius_cm / region.bottom_cap_aspect_ratio : 0.0
        barrel_zmin = region.zmin_cm + bottom_depth
        barrel_zmax = region.zmax_cm - top_depth
        if z >= barrel_zmin && z <= barrel_zmax
            return true
        elseif z > barrel_zmax
            region.has_top_cap || return false
            dz_rel = z - barrel_zmax
            c = region.top_cap_radius_cm / region.top_cap_aspect_ratio
            return r2 / region.top_cap_radius_cm^2 + (dz_rel / c)^2 <= 1.0
        else
            region.has_bottom_cap || return false
            dz_rel = barrel_zmin - z
            c = region.bottom_cap_radius_cm / region.bottom_cap_aspect_ratio
            return r2 / region.bottom_cap_radius_cm^2 + (dz_rel / c)^2 <= 1.0
        end
    end
    false
end


_is_inside_fastkernel_region(region::FastKernelRegion, pos::NTuple{3,<:Real}) =
    _is_inside_fastkernel_region(region, Float64(pos[1]), Float64(pos[2]), Float64(pos[3]))


function classify_fastkernel(fk::FastKernelGeometry, pos::NTuple{3,<:Real})
    x = Float64(pos[1])
    y = Float64(pos[2])
    z = Float64(pos[3])

    det_region = fk.regions[fk.detector_region]
    _is_inside_fastkernel_region(det_region, x, y, z) || return nothing

    air_region = fk.regions[fk.air_region]
    _is_inside_fastkernel_region(air_region, x, y, z) && return air_region

    ptfe_region = fk.regions[fk.ptfe_region]
    _is_inside_fastkernel_region(ptfe_region, x, y, z) && return ptfe_region

    rings_region = fk.regions[fk.rings_region]
    _is_inside_fastkernel_region(rings_region, x, y, z) && return rings_region

    skin_region = fk.regions[fk.skin_region]
    _is_inside_fastkernel_region(skin_region, x, y, z) && return skin_region

    fv_region = fk.regions[fk.fv_region]
    _is_inside_fastkernel_region(fv_region, x, y, z) && return fv_region

    top_active_region = fk.regions[fk.top_active_region]
    _is_inside_fastkernel_region(top_active_region, x, y, z) && return top_active_region

    barrel_active_region = fk.regions[fk.barrel_active_region]
    _is_inside_fastkernel_region(barrel_active_region, x, y, z) && return barrel_active_region

    bottom_active_region = fk.regions[fk.bottom_active_region]
    _is_inside_fastkernel_region(bottom_active_region, x, y, z) && return bottom_active_region

    passive_region = fk.regions[fk.passive_region]
    _is_inside_fastkernel_region(passive_region, x, y, z) && return passive_region

    det_region
end


function _min_positive(ts::Float64...)
    tmin = Inf
    for t in ts
        if t > 1e-10 && t < tmin
            tmin = t
        end
    end
    tmin
end


function _distance_to_cap_exit(region::FastKernelRegion,
                               x::Float64, y::Float64, z::Float64,
                               ux::Float64, uy::Float64, uz::Float64)::Float64
    R = region.has_top_cap ? region.top_cap_radius_cm : region.bottom_cap_radius_cm
    ar = region.has_top_cap ? region.top_cap_aspect_ratio : region.bottom_cap_aspect_ratio
    z_eq = region.has_top_cap ? region.zmin_cm : region.zmax_cm
    c = R / ar

    a = (ux^2 + uy^2) / R^2 + (uz^2) / c^2
    a <= 1e-20 && return Inf
    b = 2.0 * ((x * ux + y * uy) / R^2 + ((z - z_eq) * uz) / c^2)
    cc = (x^2 + y^2) / R^2 + ((z - z_eq)^2) / c^2 - 1.0
    disc = b^2 - 4.0 * a * cc
    disc < 0.0 && return Inf

    sq = sqrt(disc)
    tmin = Inf
    for t in ((-b - sq) / (2.0 * a), (-b + sq) / (2.0 * a))
        if t > 1e-10
            zh = z + uz * t
            if region.has_top_cap
                zh >= z_eq - 1e-10 || continue
            else
                zh <= z_eq + 1e-10 || continue
            end
            tmin = min(tmin, t)
        end
    end

    t_plane = Inf
    if region.has_top_cap
        if uz < -1e-20
            t = (z_eq - z) / uz
            if t > 1e-10
                xh = x + ux * t
                yh = y + uy * t
                if xh^2 + yh^2 <= R^2 + 1e-10
                    t_plane = t
                end
            end
        end
    else
        if uz > 1e-20
            t = (z_eq - z) / uz
            if t > 1e-10
                xh = x + ux * t
                yh = y + uy * t
                if xh^2 + yh^2 <= R^2 + 1e-10
                    t_plane = t
                end
            end
        end
    end

    min(tmin, t_plane)
end


function _distance_to_barrel_side_exit(rmax_cm::Float64,
                                       zbar_min::Float64, zbar_max::Float64,
                                       x::Float64, y::Float64, z::Float64,
                                       ux::Float64, uy::Float64, uz::Float64)::Float64
    a = ux^2 + uy^2
    a <= 1e-20 && return Inf
    b = 2.0 * (x * ux + y * uy)
    c = x^2 + y^2 - rmax_cm^2
    disc = b^2 - 4.0 * a * c
    disc < 0.0 && return Inf

    sq = sqrt(disc)
    tmin = Inf
    for t in ((-b - sq) / (2.0 * a), (-b + sq) / (2.0 * a))
        if t > 1e-10
            zh = z + uz * t
            if zh >= zbar_min - 1e-10 && zh <= zbar_max + 1e-10
                tmin = min(tmin, t)
            end
        end
    end
    tmin
end


function _distance_to_barrel_disk_exit(rmin_cm::Float64, rmax_cm::Float64,
                                       zface::Float64,
                                       x::Float64, y::Float64, z::Float64,
                                       ux::Float64, uy::Float64, uz::Float64)::Float64
    abs(uz) <= 1e-20 && return Inf
    t = (zface - z) / uz
    t <= 1e-10 && return Inf
    xh = x + ux * t
    yh = y + uy * t
    r2 = xh^2 + yh^2
    if r2 >= rmin_cm^2 - 1e-10 && r2 <= rmax_cm^2 + 1e-10
        return t
    end
    Inf
end


function distance_to_boundary_fastkernel(region::FastKernelRegion,
                                         pos::NTuple{3,<:Real},
                                         dir::NTuple{3,<:Real})::Float64
    x = Float64(pos[1]); y = Float64(pos[2]); z = Float64(pos[3])
    ux = Float64(dir[1]); uy = Float64(dir[2]); uz = Float64(dir[3])

    if region.kind == :cylinder
        return _min_positive(
            _distance_to_barrel_side_exit(region.rmax_cm, region.zmin_cm, region.zmax_cm, x, y, z, ux, uy, uz),
            _distance_to_barrel_disk_exit(0.0, region.rmax_cm, region.zmin_cm, x, y, z, ux, uy, uz),
            _distance_to_barrel_disk_exit(0.0, region.rmax_cm, region.zmax_cm, x, y, z, ux, uy, uz),
        )
    elseif region.kind == :cylinder_shell
        t_inner = begin
            a = ux^2 + uy^2
            if a <= 1e-20
                Inf
            else
                b = 2.0 * (x * ux + y * uy)
                c = x^2 + y^2 - region.rmin_cm^2
                disc = b^2 - 4.0 * a * c
                if disc < 0.0
                    Inf
                else
                    sq = sqrt(disc)
                    tmin = Inf
                    for t in ((-b - sq) / (2.0 * a), (-b + sq) / (2.0 * a))
                        if t > 1e-10
                            zh = z + uz * t
                            if zh >= region.zmin_cm - 1e-10 && zh <= region.zmax_cm + 1e-10
                                tmin = min(tmin, t)
                            end
                        end
                    end
                    tmin
                end
            end
        end
        return _min_positive(
            t_inner,
            _distance_to_barrel_side_exit(region.rmax_cm, region.zmin_cm, region.zmax_cm, x, y, z, ux, uy, uz),
            _distance_to_barrel_disk_exit(region.rmin_cm, region.rmax_cm, region.zmin_cm, x, y, z, ux, uy, uz),
            _distance_to_barrel_disk_exit(region.rmin_cm, region.rmax_cm, region.zmax_cm, x, y, z, ux, uy, uz),
        )
    elseif region.kind == :cap
        return _distance_to_cap_exit(region, x, y, z, ux, uy, uz)
    elseif region.kind == :domed_container || region.kind == :capped_cylinder
        top_depth = region.has_top_cap ? region.top_cap_radius_cm / region.top_cap_aspect_ratio : 0.0
        bottom_depth = region.has_bottom_cap ? region.bottom_cap_radius_cm / region.bottom_cap_aspect_ratio : 0.0
        barrel_zmin = region.zmin_cm + bottom_depth
        barrel_zmax = region.zmax_cm - top_depth

        t_side = _distance_to_barrel_side_exit(region.rmax_cm, barrel_zmin, barrel_zmax, x, y, z, ux, uy, uz)

        t_top = Inf
        if region.has_top_cap
            top_region = FastKernelRegion(
                0, "", "", :cap, region.material, region.isXe, region.sensitive, region.ecut_keV, region.fv_target,
                0.0, region.top_cap_radius_cm, barrel_zmax, barrel_zmax + top_depth,
                true, false, region.top_cap_radius_cm, region.top_cap_aspect_ratio, 0.0, 0.0
            )
            t_top = _distance_to_cap_exit(top_region, x, y, z, ux, uy, uz)
        else
            t_top = _distance_to_barrel_disk_exit(0.0, region.rmax_cm, barrel_zmax, x, y, z, ux, uy, uz)
        end

        t_bottom = Inf
        if region.has_bottom_cap
            bottom_region = FastKernelRegion(
                0, "", "", :cap, region.material, region.isXe, region.sensitive, region.ecut_keV, region.fv_target,
                0.0, region.bottom_cap_radius_cm, barrel_zmin - bottom_depth, barrel_zmin,
                false, true, 0.0, 0.0, region.bottom_cap_radius_cm, region.bottom_cap_aspect_ratio
            )
            t_bottom = _distance_to_cap_exit(bottom_region, x, y, z, ux, uy, uz)
        else
            t_bottom = _distance_to_barrel_disk_exit(0.0, region.rmax_cm, barrel_zmin, x, y, z, ux, uy, uz)
        end

        return _min_positive(t_side, t_top, t_bottom)
    end

    Inf
end


function _node_depth(det::TrackingDetector, node::DetectorNode)::Int
    depth = 0
    current = node
    while current.parent_id != 0
        depth += 1
        current = det.nodes[current.parent_id]
    end
    depth
end


function _child_preference(a::DetectorNode, b::DetectorNode)::DetectorNode
    a_exact = a.lv.approximation == "exact"
    b_exact = b.lv.approximation == "exact"
    if a_exact != b_exact
        return a_exact ? a : b
    end

    va = volume(a)
    vb = volume(b)
    va <= vb ? a : b
end


function find_tracking_node(det::TrackingDetector, pos::Vector{Float64})::Union{DetectorNode,Nothing}
    root = det.nodes[det.root_id]
    is_inside(root, pos) || return nothing

    current = root
    while true
        best_child = nothing
        for child_id in current.children_ids
            child = det.nodes[child_id]
            if is_inside(child, pos)
                best_child = best_child === nothing ? child : _child_preference(best_child, child)
            end
        end
        best_child === nothing && return current
        current = best_child
    end
end


find_tracking_node(det::TrackingDetector, pos::NTuple{3,<:Real}) = find_tracking_node(det, Float64[pos...])
