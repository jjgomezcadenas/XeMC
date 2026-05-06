#!/usr/bin/env julia
"""
draw_geometry.jl — Draw a Geometry V2 subtree as embedded TikZ.

Usage:
    julia --project=. scripts/draw_geometry.jl
    julia --project=. scripts/draw_geometry.jl --root LZ_detector --out design/latex/geometry_draw.tex
"""

using ArgParse
using LXeMC
using Printf


function parse_args()
    s = ArgParseSettings(
        description = "Draw a Geometry V2 subtree as a LaTeX/TikZ wrapper",
        version = "0.1.0",
        add_version = true
    )

    @add_arg_table! s begin
        "--detector", "-d"
            help = "Path to Geometry V2 detector JSON"
            arg_type = String
            default = ""
        "--root", "-r"
            help = "Root node name for the drawn subtree"
            arg_type = String
            default = "LZ_detector"
        "--include-children"
            help = "Draw the full subtree below --root instead of just the root node"
            action = :store_true
        "--out", "-o"
            help = "Output .tex path"
            arg_type = String
            default = "design/latex/geometry_draw.tex"
        "--scale"
            help = "Resize factor relative to linewidth in wrapper"
            arg_type = Float64
            default = 0.55
    end

    ArgParse.parse_args(s)
end


function z_bounds(node::DetectorNode)
    s = node.lv.solid
    z0 = node.placement.position_cm[3]
    if s isa Cyl
        return (z0 - s.half_height_cm, z0 + s.half_height_cm)
    elseif s isa CylShell
        return (z0 - s.half_height_cm, z0 + s.half_height_cm)
    elseif s isa Box
        return (z0 - s.half_z_cm, z0 + s.half_z_cm)
    elseif s isa Disk
        if is_flat(s)
            if node.placement.orientation === :down
                return (z0 - s.wall_thickness_cm, z0)
            end
            return (z0, z0 + s.wall_thickness_cm)
        end
        if node.placement.orientation === :down
            return (z0 - (depth(s) + s.wall_thickness_cm), z0)
        end
        return (z0, z0 + depth(s) + s.wall_thickness_cm)
    elseif s isa Cap
        if node.placement.orientation === :down
            return (z0 - depth(s), z0)
        end
        return (z0, z0 + depth(s))
    elseif s isa DomedContainer
        return (z0 - (s.barrel_half_height_cm + depth(s.bottom_cap)),
                z0 + s.barrel_half_height_cm + depth(s.top_cap))
    elseif s isa CappedCylinder
        zmin = z0 - s.barrel_half_height_cm
        zmax = z0 + s.barrel_half_height_cm
        s.bottom_cap !== nothing && (zmin -= depth(s.bottom_cap))
        s.top_cap !== nothing && (zmax += depth(s.top_cap))
        return (zmin, zmax)
    end
    error("Unsupported solid type for z_bounds")
end


function r_bounds(node::DetectorNode)
    s = node.lv.solid
    if s isa Cyl
        return (0.0, s.radius_cm)
    elseif s isa CylShell
        return (s.R_inner_cm, R_outer(s))
    elseif s isa Box
        return (0.0, max(s.half_x_cm, s.half_y_cm))
    elseif s isa Disk
        return (0.0, s.radius_cm)
    elseif s isa Cap
        return (0.0, s.radius_cm)
    elseif s isa DomedContainer
        return (0.0, s.radius_cm)
    elseif s isa CappedCylinder
        return (0.0, s.radius_cm)
    end
    error("Unsupported solid type for r_bounds")
end


function subtree_ids(det::DetectorV2, root_id::Int)
    ids = Int[]
    function visit(id::Int)
        push!(ids, id)
        for child_id in det.nodes[id].children_ids
            visit(child_id)
        end
    end
    visit(root_id)
    ids
end


function fill_style(node::DetectorNode)
    tag = node.lv.tag
    if tag == TAG_FV
        return "fill=blue!18, draw=none"
    elseif tag == TAG_TPC_ACTIVE
        return "fill=blue!8, draw=none"
    elseif tag == TAG_SKIN
        return "fill=blue!4, draw=none"
    elseif tag == TAG_PASSIVE_LXE
        return "fill=blue!5, draw=none"
    elseif tag == TAG_VACUUM || tag == TAG_WORLD
        return ""
    end
    return ""
end


function boundary_style(node::DetectorNode)
    if node.lv.tag == TAG_STRUCTURAL
        return "black, line width=0.9pt"
    elseif node.lv.tag == TAG_VACUUM
        return "blue!55!black, dashed, line width=0.8pt"
    elseif node.lv.tag == TAG_WORLD
        return "gray!60, dotted, line width=0.7pt"
    elseif node.lv.tag in (TAG_TPC_ACTIVE, TAG_FV, TAG_SKIN, TAG_PASSIVE_LXE)
        return "black, line width=0.8pt"
    end
    return "black, line width=0.8pt"
end


function label_pos(node::DetectorNode)
    s = node.lv.solid
    cx, _, cz = node.placement.position_cm
    if s isa Cyl
        return (cx, cz)
    elseif s isa CylShell
        return (s.R_inner_cm + 0.5 * s.wall_thickness_cm, cz)
    elseif s isa Box
        return (cx, cz)
    elseif s isa Disk
        zmin, zmax = z_bounds(node)
        return (0.0, 0.5 * (zmin + zmax))
    elseif s isa Cap
        zmin, zmax = z_bounds(node)
        return (0.0, 0.5 * (zmin + zmax))
    elseif s isa DomedContainer
        return (0.55 * s.radius_cm, cz + 0.7 * s.barrel_half_height_cm)
    elseif s isa CappedCylinder
        return (0.55 * s.radius_cm, cz)
    end
    return (0.0, cz)
end


function draw_head_shell(io, style::String, z0::Float64, R::Float64, d::Float64, t::Float64, orientation::Symbol)
    if orientation === :up
        @printf(io, "\\draw[%s] (-%.3f, %.3f) arc[start angle=180, end angle=0, x radius=%.3f, y radius=%.3f];\n",
                style, R, z0, R, d)
        @printf(io, "\\draw[%s] (-%.3f, %.3f) arc[start angle=180, end angle=0, x radius=%.3f, y radius=%.3f];\n",
                style, R, z0 + t, R, d + t)
    else
        @printf(io, "\\draw[%s] (-%.3f, %.3f) arc[start angle=180, end angle=360, x radius=%.3f, y radius=%.3f];\n",
                style, R, z0, R, d)
        @printf(io, "\\draw[%s] (-%.3f, %.3f) arc[start angle=180, end angle=360, x radius=%.3f, y radius=%.3f];\n",
                style, R, z0 - t, R, d + t)
    end
end


function draw_cap(io, style::String, fill::String, z0::Float64, R::Float64, d::Float64, orientation::Symbol)
    if orientation === :up
        !isempty(fill) && @printf(io, "\\path[%s] (0, %.3f) ellipse [x radius=%.3f, y radius=%.3f];\n", fill, z0, R, d)
        @printf(io, "\\draw[%s] (-%.3f, %.3f) arc[start angle=180, end angle=0, x radius=%.3f, y radius=%.3f];\n",
                style, R, z0, R, d)
    else
        !isempty(fill) && @printf(io, "\\path[%s] (0, %.3f) ellipse [x radius=%.3f, y radius=%.3f];\n", fill, z0, R, d)
        @printf(io, "\\draw[%s] (-%.3f, %.3f) arc[start angle=180, end angle=360, x radius=%.3f, y radius=%.3f];\n",
                style, R, z0, R, d)
    end
end


function draw_domed_container(io, style::String, fill::String, node::DetectorNode, dc::DomedContainer)
    cx, _, cz = node.placement.position_cm
    R = dc.radius_cm
    H = dc.barrel_half_height_cm
    z1 = cz - H
    z2 = cz + H
    !isempty(fill) && @printf(io, "\\path[%s] (-%.3f, %.3f) rectangle (%.3f, %.3f);\n", fill, R, z1, R, z2)
    @printf(io, "\\draw[%s] (-%.3f, %.3f) rectangle (%.3f, %.3f);\n", style, R, z1, R, z2)
    draw_cap(io, style, fill, z2, R, depth(dc.top_cap), :up)
    draw_cap(io, style, fill, z1, R, depth(dc.bottom_cap), :down)
end


function draw_capped_cylinder(io, style::String, fill::String, node::DetectorNode, cc::CappedCylinder)
    _, _, cz = node.placement.position_cm
    R = cc.radius_cm
    H = cc.barrel_half_height_cm
    z1 = cz - H
    z2 = cz + H
    !isempty(fill) && @printf(io, "\\path[%s] (-%.3f, %.3f) rectangle (%.3f, %.3f);\n", fill, R, z1, R, z2)
    @printf(io, "\\draw[%s] (-%.3f, %.3f) rectangle (%.3f, %.3f);\n", style, R, z1, R, z2)
    cc.top_cap !== nothing && draw_cap(io, style, fill, z2, R, depth(cc.top_cap), :up)
    cc.bottom_cap !== nothing && draw_cap(io, style, fill, z1, R, depth(cc.bottom_cap), :down)
end


function draw_node(io, node::DetectorNode)
    style = boundary_style(node)
    fill = fill_style(node)
    s = node.lv.solid
    zmin, zmax = z_bounds(node)
    rmin, rmax = r_bounds(node)

    if s isa Cyl
        !isempty(fill) && @printf(io, "\\path[%s] (-%.3f, %.3f) rectangle (%.3f, %.3f);\n", fill, rmax, zmin, rmax, zmax)
        @printf(io, "\\draw[%s] (-%.3f, %.3f) rectangle (%.3f, %.3f);\n", style, rmax, zmin, rmax, zmax)
    elseif s isa CylShell
        !isempty(fill) && begin
            @printf(io, "\\path[%s] (-%.3f, %.3f) rectangle (-%.3f, %.3f);\n", fill, rmax, zmin, rmin, zmax)
            @printf(io, "\\path[%s] (%.3f, %.3f) rectangle (%.3f, %.3f);\n", fill, rmin, zmin, rmax, zmax)
        end
        @printf(io, "\\draw[%s] (-%.3f, %.3f) rectangle (-%.3f, %.3f);\n", style, rmax, zmin, rmin, zmax)
        @printf(io, "\\draw[%s] (%.3f, %.3f) rectangle (%.3f, %.3f);\n", style, rmin, zmin, rmax, zmax)
    elseif s isa Disk
        z0 = node.placement.position_cm[3]
        if is_flat(s)
            @printf(io, "\\draw[%s] (-%.3f, %.3f) rectangle (%.3f, %.3f);\n", style, rmax, zmin, rmax, zmax)
        else
            draw_head_shell(io, style, z0, rmax, depth(s), s.wall_thickness_cm, node.placement.orientation)
        end
    elseif s isa Cap
        draw_cap(io, style, fill, node.placement.position_cm[3], rmax, depth(s), node.placement.orientation)
    elseif s isa DomedContainer
        draw_domed_container(io, style, fill, node, s)
    elseif s isa CappedCylinder
        draw_capped_cylinder(io, style, fill, node, s)
    elseif s isa Box
        !isempty(fill) && @printf(io, "\\path[%s] (-%.3f, %.3f) rectangle (%.3f, %.3f);\n", fill, rmax, zmin, rmax, zmax)
        @printf(io, "\\draw[%s] (-%.3f, %.3f) rectangle (%.3f, %.3f);\n", style, rmax, zmin, rmax, zmax)
    else
        error("Unsupported solid type for drawing")
    end

    x, z = label_pos(node)
    label = replace(node.lv.name, "_" => "\\_")
    @printf(io, "\\node[scale=0.42, fill=white, fill opacity=0.88, text opacity=1, inner sep=1.0pt] at (%.3f, %.3f) {%s};\n",
            x, z, label)
end


function write_tex(path::String, det::DetectorV2, root_name::String, scale::Float64; include_children::Bool=false)
    root = node_by_name(det, root_name)
    ids = include_children ? subtree_ids(det, root.id) : [root.id]
    nodes = [det.nodes[id] for id in ids]
    sort!(nodes, by=node -> volume(node), rev=true)

    open(path, "w") do io
        println(io, "% Auto-generated by scripts/draw_geometry.jl")
        println(io, "\\documentclass{article}")
        println(io, "\\usepackage[margin=1.5cm]{geometry}")
        println(io, "\\usepackage{tikz}")
        println(io, "\\usepackage{xcolor}")
        println(io, "\\usepackage{graphicx}")
        println(io, "\\pagestyle{empty}")
        println(io, "\\begin{document}")
        println(io, "\\begin{figure}[p]")
        println(io, "\\centering")
        @printf(io, "\\resizebox{%.2f\\linewidth}{!}{%%\n", scale)
        println(io, "\\begin{tikzpicture}[x=0.03cm,y=0.03cm]")
        for node in nodes
            draw_node(io, node)
        end
        println(io, "\\end{tikzpicture}")
        println(io, "}")
        println(io, "\\end{figure}")
        println(io, "\\end{document}")
    end
end


function main()
    args = parse_args()
    det_path = isempty(args["detector"]) ? default_detector_v2_path() : args["detector"]
    cfg = default_config()
    mats = load_materials(cfg)
    det = load_detector_v2(det_path, mats; validate=false)
    write_tex(args["out"], det, args["root"], args["scale"]; include_children=args["include-children"])
    println("Geometry drawing written to: $(args["out"])")
end


main()
