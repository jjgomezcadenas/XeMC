"""
Geometry V2: explicit detector hierarchy and semantic metadata.

This file adds a second geometry representation alongside the current
flat `Detector` model. It is loader-focused for now: it parses the
V2 JSON schema, builds parent/child relationships, and exposes basic
inspection helpers. Transport integration comes later.
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


struct PlacementV2
    position_cm::Vector{Float64}
    orientation::Symbol
end


struct LogicalVolumeV2
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
    lv::LogicalVolumeV2
    placement::PlacementV2
end


struct DetectorV2
    name::String
    nodes::Vector{DetectorNode}
    name_to_id::Dict{String,Int}
    root_id::Int
end


struct FastKernelVolume
    name::String
    material::Material
    role::String
    inFastKernel::Bool
    isXe::Bool
    sensitive::Bool
    ecut_keV::Float64
    fv_target::Bool
end


struct FastKernelGeometry
    volumes::Vector{FastKernelVolume}
    name_to_index::Dict{String,Int}
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


function classify_runtime_v2(det::DetectorV2, pos::NTuple{3,<:Real})
    node = find_node_v2(det, pos)
    node === nothing && return nothing
    (
        node = node,
        name = node.lv.name,
        sensitive = is_sensitive(node),
        ecut_keV = ecut_keV(node),
        dz_mm = dz_mm(node),
        fv_target = is_fv_target(node)
    )
end


function select_interaction(result, det::DetectorV2)
    cls = classify_runtime_v2(det, (result.position[1], result.position[2], result.position[3]))
    cls === nothing && return (
        class = :other,
        sensitive = false,
        passes_threshold = false,
        ecut_keV = 0.0,
        region = "MARS",
        interaction_type = result.interaction_type
    )

    class = if cls.fv_target
        :fv
    elseif cls.name == "LXeTPC"
        :tpc
    elseif cls.name == "Skin"
        :skin
    else
        :other
    end

    dep_keV = 1000.0 * result.deposit_E_MeV
    passes_threshold = cls.sensitive && dep_keV >= cls.ecut_keV

    (
        class = class,
        sensitive = cls.sensitive,
        passes_threshold = passes_threshold,
        ecut_keV = cls.ecut_keV,
        region = cls.name,
        interaction_type = result.interaction_type
    )
end


function _same_runtime_node(a, b)::Bool
    a === nothing && return b === nothing
    b === nothing && return false
    a.node.id == b.node.id
end


function _runtime_node_id(det::DetectorV2, x::Float64, y::Float64, z::Float64)
    node = find_node_v2(det, (x, y, z))
    node === nothing ? 0 : node.id
end


function distance_to_node_change_v2(det::DetectorV2,
                                    pos::Vector{Float64},
                                    dir::Vector{Float64};
                                    step_cm::Float64=0.05,
                                    max_cm::Float64=400.0,
                                    tol_cm::Float64=1e-4)
    x0, y0, z0 = pos
    dx, dy, dz = dir
    current_id = _runtime_node_id(det, x0, y0, z0)
    t_lo = 0.0
    t_hi = step_cm

    while t_lo < max_cm
        thi = min(t_hi, max_cm)
        cls_hi_id = _runtime_node_id(det, x0 + dx * thi, y0 + dy * thi, z0 + dz * thi)

        if cls_hi_id != current_id
            lo = t_lo
            hi = thi
            while hi - lo > tol_cm
                mid = 0.5 * (lo + hi)
                cls_mid_id = _runtime_node_id(det, x0 + dx * mid, y0 + dy * mid, z0 + dz * mid)
                if cls_mid_id == current_id
                    lo = mid
                else
                    hi = mid
                end
            end
            next_cls = classify_runtime_v2(det, (x0 + dx * hi, y0 + dy * hi, z0 + dz * hi))
            return hi, next_cls
        end

        t_lo = t_hi
        t_hi += step_cm
    end

    Inf, nothing
end


function geometric_prefilter_v2(det::DetectorV2, gammas;
                                cfg::SimConfig=default_config(),
                                rng::AbstractRNG=Random.default_rng())::Bool
    isempty(gammas) && return false
    any(g -> begin
        result = propagate_gamma_v2(g, det, cfg, rng)
        result.status == :entered_fv || result.region == "FV"
    end, gammas)
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


function _build_solid_v2(name::String, d)
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


function _build_node_v2(id::Int, d, material::Material)::Tuple{DetectorNode,String}
    name = String(d["name"])
    solid = _build_solid_v2(name, d)
    tag = _parse_region_tag(d["tag"])
    role = String(get(d, "role", ""))
    inFastKernel = Bool(get(d, "inFastKernel", false))
    isXe = Bool(get(d, "isXe", false))
    sensitive = Bool(get(d, "sensitive", false))
    ecut_keV = Float64(get(d, "ecut_keV", 0.0))
    dz_mm = Float64(get(d, "dz_mm", 0.0))
    fv_target = Bool(get(d, "fv_target", false))
    approximation = String(get(d, "approximation", "exact"))
    activity = _parse_activity(d)
    lv = LogicalVolumeV2(name, solid, material, tag, role, inFastKernel, isXe, sensitive, ecut_keV, dz_mm, fv_target, approximation, activity)

    pos = Float64.(get(d, "position_cm", [0.0, 0.0, 0.0]))
    orientation = Symbol(get(d, "orientation", "none"))
    placement = PlacementV2(pos, orientation)

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


function _is_inside_cap(cap::Cap, placement::PlacementV2, pos::Vector{Float64})::Bool
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


function _is_inside_domed_container(dc::DomedContainer, placement::PlacementV2, pos::Vector{Float64})::Bool
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


function _is_inside_capped_cylinder(cc::CappedCylinder, placement::PlacementV2, pos::Vector{Float64})::Bool
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


function _validate_reachability(det::DetectorV2)
    seen = falses(length(det.nodes))
    stack = Int[det.root_id]
    while !isempty(stack)
        id = pop!(stack)
        seen[id] && continue
        seen[id] = true
        append!(stack, det.nodes[id].children_ids)
    end
    all(seen) || error("DetectorV2 contains unreachable nodes")
end


function _validate_child_containment(det::DetectorV2)
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


function _validate_sibling_overlaps(det::DetectorV2)
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


function sibling_overlaps(det::DetectorV2; exact_only::Bool=true)
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


function validate_detector_v2(det::DetectorV2)::Bool
    root = root_node(det)
    root.parent_id == 0 || error("Root node '$(root.lv.name)' must have parent_id = 0")
    root.lv.tag == TAG_WORLD || error("Root node '$(root.lv.name)' must have tag TAG_WORLD")
    _validate_reachability(det)
    _validate_child_containment(det)
    _validate_sibling_overlaps(det)
    true
end


function tree_dump(det::DetectorV2)::String
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
    load_detector_v2(path, materials) -> DetectorV2

Load a Geometry V2 detector JSON. The file defines a unique `world`
entry and a flat list of volumes with explicit `parent` links.
"""
function load_detector_v2(path::AbstractString,
                          materials::Dict{String,Material};
                          validate::Bool=true)::DetectorV2
    raw = open(path, "r") do io
        JSON.parse(io)
    end

    det_name = String(raw["name"])

    haskey(raw, "world") || error("Geometry V2 file '$path' is missing 'world'")
    haskey(raw, "volumes") || error("Geometry V2 file '$path' is missing 'volumes'")

    nodes = DetectorNode[]
    name_to_id = Dict{String,Int}()
    parent_names = String[]

    world_raw = raw["world"]
    world_name = String(world_raw["name"])
    world_mat_name = String(world_raw["material"])
    haskey(materials, world_mat_name) || error("Unknown material '$world_mat_name' for world '$world_name'")
    world_node, world_parent = _build_node_v2(1, world_raw, materials[world_mat_name])
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
        node, parent_name = _build_node_v2(id, d, materials[mat_name])
        isempty(parent_name) && error("Geometry V2 node '$name' is missing 'parent'")

        push!(nodes, node)
        name_to_id[name] = id
        push!(parent_names, parent_name)
    end

    for node in nodes[2:end]
        parent_name = parent_names[node.id]
        parent_id = get(name_to_id, parent_name, 0)
        parent_id == 0 && error("Geometry V2 node '$(node.lv.name)' references unknown parent '$parent_name'")
        node.parent_id = parent_id
        push!(nodes[parent_id].children_ids, node.id)
    end

    det = DetectorV2(det_name, nodes, name_to_id, 1)
    validate && validate_detector_v2(det)
    det
end


root_node(det::DetectorV2) = det.nodes[det.root_id]


function node_by_name(det::DetectorV2, name::AbstractString)::DetectorNode
    id = get(det.name_to_id, String(name), 0)
    id == 0 && error("DetectorV2 has no node named '$name'")
    det.nodes[id]
end


function child_nodes(det::DetectorV2, node::DetectorNode)::Vector{DetectorNode}
    [det.nodes[id] for id in node.children_ids]
end


function child_nodes(det::DetectorV2, name::AbstractString)::Vector{DetectorNode}
    child_nodes(det, node_by_name(det, name))
end


function detector_summary(det::DetectorV2)::String
    root = root_node(det)
    "DetectorV2($(det.name)): $(length(det.nodes)) nodes, root=$(root.lv.name)"
end


function compile_fastkernel_geometry(det::DetectorV2)::FastKernelGeometry
    volumes = FastKernelVolume[]
    name_to_index = Dict{String,Int}()

    for node in det.nodes
        node.lv.inFastKernel || continue
        push!(volumes, FastKernelVolume(
            node.lv.name,
            node.lv.material,
            node.lv.role,
            node.lv.inFastKernel,
            node.lv.isXe,
            node.lv.sensitive,
            node.lv.ecut_keV,
            node.lv.fv_target
        ))
        name_to_index[node.lv.name] = length(volumes)
    end

    FastKernelGeometry(volumes, name_to_index)
end


function _node_depth(det::DetectorV2, node::DetectorNode)::Int
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


function find_node_v2(det::DetectorV2, pos::Vector{Float64})::Union{DetectorNode,Nothing}
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


find_node_v2(det::DetectorV2, pos::NTuple{3,<:Real}) = find_node_v2(det, Float64[pos...])
