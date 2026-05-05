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


struct PlacementV2
    position_cm::Vector{Float64}
    orientation::Symbol
end


struct LogicalVolumeV2
    name::String
    solid::Union{Cyl,CylShell,Box,Disk,Cap}
    material::Material
    tag::RegionTag
    role::String
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
    end
    error("Unknown shape '$shape' for volume '$name'")
end


function _build_node_v2(id::Int, d, material::Material)::Tuple{DetectorNode,String}
    name = String(d["name"])
    solid = _build_solid_v2(name, d)
    tag = _parse_region_tag(d["tag"])
    role = String(get(d, "role", ""))
    approximation = String(get(d, "approximation", "exact"))
    activity = _parse_activity(d)
    lv = LogicalVolumeV2(name, solid, material, tag, role, approximation, activity)

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
    elseif solid isa Cap
        return nothing
    end
    error("Unsupported solid type for node '$(node.lv.name)'")
end


function _is_inside_cap(cap::Cap, placement::PlacementV2, pos::Vector{Float64})::Bool
    dx = pos[1] - placement.position_cm[1]
    dy = pos[2] - placement.position_cm[2]
    dz = pos[3] - placement.position_cm[3]
    r2 = dx^2 + dy^2
    R = cap.radius_cm
    r2 > R^2 && return false

    sgn = placement.orientation === :down ? -1.0 : 1.0
    dz_rel = sgn * dz
    dz_rel < 0.0 && return false

    c = depth(cap)
    r_frac2 = r2 / R^2
    r_frac2 + (dz_rel / c)^2 <= 1.0
end


function is_inside(node::DetectorNode, pos::Vector{Float64})
    solid = node.lv.solid
    if solid isa Cap
        return _is_inside_cap(solid, node.placement, pos)
    end
    is_inside(_logical_volume(node), pos)
end


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
    for p in _sample_points(a)
        is_inside(b, p) && return true
    end
    for p in _sample_points(b)
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
