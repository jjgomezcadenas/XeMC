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


struct PlacementV2
    position_cm::Vector{Float64}
    orientation::Symbol
end


struct LogicalVolumeV2
    name::String
    solid::Union{Cyl,CylShell,Box,Disk}
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


"""
    load_detector_v2(path, materials) -> DetectorV2

Load a Geometry V2 detector JSON. The file defines a unique `world`
entry and a flat list of volumes with explicit `parent` links.
"""
function load_detector_v2(path::AbstractString,
                          materials::Dict{String,Material})::DetectorV2
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

    DetectorV2(det_name, nodes, name_to_id, 1)
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
