#!/usr/bin/env julia
"""
inspect_geometry_v2.jl — Inspect the Geometry V2 detector tree.

Reports:
- detector summary
- tree dump
- sibling overlaps (exact-only and all approximations)

Optional:
- write a TikZ R-z schematic with line styles keyed to approximation type

Usage:
    julia --project=.. scripts/inspect_geometry_v2.jl
    julia --project=.. scripts/inspect_geometry_v2.jl --tikz-out design/lz_geometry_v2.tex
    julia --project=.. scripts/inspect_geometry_v2.jl --detector ../data/detector_lz_v2.json
"""

using ArgParse
using LXeMC
using Printf


function parse_args()
    s = ArgParseSettings(
        description = "Inspect LXeMC Geometry V2 detector hierarchy and overlaps",
        version = "0.1.0",
        add_version = true
    )

    @add_arg_table! s begin
        "--detector", "-d"
            help = "Path to Geometry V2 detector JSON"
            arg_type = String
            default = ""
        "--tikz-out", "-t"
            help = "Write a TikZ R-z schematic to this file"
            arg_type = String
            default = ""
        "--tikz-profile"
            help = "TikZ drawing profile: 'core' or 'full'"
            arg_type = String
            default = "core"
        "--latex-mode"
            help = "Output mode for --tikz-out: 'bare', 'standalone', or 'wrapper'"
            arg_type = String
            default = "wrapper"
        "--standalone"
            help = "Deprecated alias for --latex-mode standalone"
            action = :store_true
        "--no-tree"
            help = "Do not print the detector tree"
            action = :store_true
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
    elseif s isa LXeMC.Cap
        if node.placement.orientation === :down
            return (z0 - LXeMC.depth(s), z0)
        end
        return (z0, z0 + LXeMC.depth(s))
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
    elseif s isa LXeMC.Cap
        return (0.0, s.radius_cm)
    end
    error("Unsupported solid type for r_bounds")
end


function approximation_style(approx::String)
    if approx == "exact"
        return "black, line width=0.9pt"
    elseif approx == "cylindrical_proxy"
        return "blue, dashed, line width=0.9pt"
    elseif approx == "mass_equivalent_shell"
        return "orange!80!black, dashed, line width=0.8pt"
    elseif approx == "mass_equivalent_slab"
        return "red!75!black, dashed, line width=0.8pt"
    elseif approx == "negligible_thickness_proxy"
        return "gray, dotted, line width=0.8pt"
    end
    return "black, dotted, line width=0.8pt"
end


function should_draw_node(node::DetectorNode, profile::String)::Bool
    profile == "full" && return true
    profile == "core" || error("Unknown tikz profile '$profile'")
    node.lv.name in (
        "LXeTPC", "FV", "FieldCage", "Skin", "RFR", "DomeBarrel", "DomeBottomCap",
        "OCV_barrel", "OCV_top", "OCV_bottom",
        "ICV_barrel", "ICV_top", "ICV_bottom"
    )
end


function fill_style(node::DetectorNode)
    name = node.lv.name
    if name == "FV"
        return "fill=blue!20, draw=none"
    elseif name == "LXeTPC"
        return "fill=blue!8, draw=none"
    elseif name == "Skin"
        return "fill=blue!4, draw=none"
    elseif name == "RFR"
        return "fill=blue!10, draw=none"
    elseif name in ("DomeBarrel", "DomeBottomCap")
        return "fill=blue!6, draw=none"
    end
    return ""
end


function boundary_style(node::DetectorNode)
    name = node.lv.name
    approx = node.lv.approximation
    if name in ("OCV_barrel", "OCV_top", "OCV_bottom")
        return "black, line width=1.2pt"
    elseif name in ("ICV_barrel", "ICV_top", "ICV_bottom")
        return "black, line width=1.1pt"
    elseif name in ("FieldCage", "Skin", "LXeTPC", "FV", "RFR", "DomeBarrel", "DomeBottomCap")
        return "black, line width=0.8pt"
    elseif name in ("OCV_void", "ICV_LXe_interior")
        return "blue!50, dashed, line width=0.5pt"
    end
    return approximation_style(approx)
end


function label_color(node::DetectorNode)
    name = node.lv.name
    if name in ("FV", "LXeTPC")
        return "blue!60!black"
    elseif name == "Skin"
        return "teal!60!black"
    elseif name in ("RFR", "DomeBarrel", "DomeBottomCap")
        return "blue!30!black"
    elseif name in ("FieldCage", "OCV_barrel", "OCV_top", "OCV_bottom", "ICV_barrel", "ICV_top", "ICV_bottom")
        return "black"
    end
    return "black"
end


function label_text(node::DetectorNode)
    name = node.lv.name
    if startswith(name, "OCV_")
        return "OCV"
    elseif startswith(name, "ICV_")
        return "ICV"
    elseif name == "DomeBarrel" || name == "DomeBottomCap"
        return "Dome"
    end
    name
end


function label_pos(node::DetectorNode)
    name = node.lv.name
    if name == "FV"
        return (0.0, 61.0)
    elseif name == "LXeTPC"
        return (0.0, 112.0)
    elseif name == "FieldCage"
        return (-73.55, 112.0)
    elseif name == "Skin"
        return (78.2, 112.0)
    elseif name == "RFR"
        return (0.0, -6.5)
    elseif name == "DomeBarrel"
        return (0.0, -28.0)
    elseif name == "DomeBottomCap"
        return nothing
    elseif name == "ICV_barrel"
        return (56.0, 158.0)
    elseif name == "OCV_barrel"
        return (42.0, -74.0)
    elseif name == "ICV_top"
        return nothing
    elseif name == "OCV_top"
        return nothing
    elseif name == "ICV_bottom"
        return nothing
    elseif name == "OCV_bottom"
        return nothing
    end
    return nothing
end


function label_style(node::DetectorNode)
    name = node.lv.name
    if name in ("FieldCage", "Skin")
        return "scale=0.42, rotate=90"
    elseif name in ("ICV_barrel", "OCV_barrel")
        return "scale=0.46"
    end
    return "scale=0.48"
end


function label_fill(node::DetectorNode)
    name = node.lv.name
    if name in ("ICV_barrel", "OCV_barrel")
        return "gray!18"
    end
    return "white"
end


function draw_head_shell(io, style::String, z0::Float64, R::Float64, d::Float64, t::Float64, orientation::Symbol)
    if orientation === :up
        @printf(io, "\\draw[%s] (%.3f, %.3f) arc[start angle=180, end angle=0, x radius=%.3f, y radius=%.3f];\n",
                style, -R, z0, R, d)
        @printf(io, "\\draw[%s] (%.3f, %.3f) arc[start angle=180, end angle=0, x radius=%.3f, y radius=%.3f];\n",
                style, -R, z0 + t, R, d + t)
    elseif orientation === :down
        @printf(io, "\\draw[%s] (%.3f, %.3f) arc[start angle=180, end angle=360, x radius=%.3f, y radius=%.3f];\n",
                style, -R, z0, R, d)
        @printf(io, "\\draw[%s] (%.3f, %.3f) arc[start angle=180, end angle=360, x radius=%.3f, y radius=%.3f];\n",
                style, -R, z0 - t, R, d + t)
    else
        error("Unsupported disk orientation '$orientation'")
    end
end


function draw_cap(io, style::String, fill::String, z0::Float64, R::Float64, d::Float64, orientation::Symbol)
    if orientation === :down
        !isempty(fill) && @printf(io, "\\path[%s] (0, %.3f) ellipse [x radius=%.3f, y radius=%.3f];\n", fill, z0, R, d)
        @printf(io, "\\draw[%s] (-%.3f, %.3f) arc[start angle=180, end angle=360, x radius=%.3f, y radius=%.3f];\n",
                style, R, z0, R, d)
    elseif orientation === :up
        !isempty(fill) && @printf(io, "\\path[%s] (0, %.3f) ellipse [x radius=%.3f, y radius=%.3f];\n", fill, z0, R, d)
        @printf(io, "\\draw[%s] (-%.3f, %.3f) arc[start angle=180, end angle=0, x radius=%.3f, y radius=%.3f];\n",
                style, R, z0, R, d)
    else
        error("Unsupported cap orientation '$orientation'")
    end
end


function draw_label(io, node::DetectorNode)
    pos = label_pos(node)
    pos === nothing && return
    x, y = pos
    color = label_color(node)
    fill = label_fill(node)
    text = replace(label_text(node), "_" => "\\_")
    style = label_style(node)
    @printf(io, "\\node[%s, text=%s, fill=%s, fill opacity=0.92, text opacity=1, inner sep=1.0pt] at (%.3f, %.3f) {%s};\n",
            style, color, fill, x, y, text)
end


function tikz_for_node(io, node::DetectorNode)
    style = boundary_style(node)
    fstyle = fill_style(node)
    zmin, zmax = z_bounds(node)
    rmin, rmax = r_bounds(node)
    s = node.lv.solid

    if s isa Cyl
        !isempty(fstyle) && @printf(io, "\\path[%s] (%.3f, %.3f) rectangle (%.3f, %.3f);\n", fstyle, -rmax, zmin, rmax, zmax)
        @printf(io, "\\draw[%s] (%.3f, %.3f) rectangle (%.3f, %.3f);\n", style, -rmax, zmin, rmax, zmax)
    elseif s isa CylShell
        !isempty(fstyle) && begin
            @printf(io, "\\path[%s] (%.3f, %.3f) rectangle (%.3f, %.3f);\n", fstyle, -rmax, zmin, -rmin, zmax)
            @printf(io, "\\path[%s] (%.3f, %.3f) rectangle (%.3f, %.3f);\n", fstyle, rmin, zmin, rmax, zmax)
        end
        @printf(io, "\\draw[%s] (%.3f, %.3f) rectangle (%.3f, %.3f);\n", style, -rmax, zmin, -rmin, zmax)
        @printf(io, "\\draw[%s] (%.3f, %.3f) rectangle (%.3f, %.3f);\n", style, rmin, zmin, rmax, zmax)
    elseif s isa Box
        !isempty(fstyle) && @printf(io, "\\path[%s] (%.3f, %.3f) rectangle (%.3f, %.3f);\n", fstyle, -rmax, zmin, rmax, zmax)
        @printf(io, "\\draw[%s] (%.3f, %.3f) rectangle (%.3f, %.3f);\n", style, -rmax, zmin, rmax, zmax)
    elseif s isa Disk
        z0 = node.placement.position_cm[3]
        d = depth(s)
        t = s.wall_thickness_cm
        if is_flat(s)
            @printf(io, "\\draw[%s] (%.3f, %.3f) rectangle (%.3f, %.3f);\n", style, -rmax, zmin, rmax, zmax)
        else
            draw_head_shell(io, style, z0, rmax, d, t, node.placement.orientation)
        end
    elseif s isa LXeMC.Cap
        z0 = node.placement.position_cm[3]
        d = LXeMC.depth(s)
        draw_cap(io, style, fstyle, z0, rmax, d, node.placement.orientation)
    end
end


function write_tikz(io, det::DetectorV2; profile::String="core")
    scale = profile == "core" ? 0.025 : 0.06
    @printf(io, "\\begin{tikzpicture}[x=%.3fcm,y=%.3fcm]\n", scale, scale)
    for node in det.nodes
        node.id == det.root_id && continue
        should_draw_node(node, profile) || continue
        tikz_for_node(io, node)
    end
    if profile == "core"
        for node in det.nodes
            node.id == det.root_id && continue
            should_draw_node(node, profile) || continue
            draw_label(io, node)
        end
    end
    println(io, "\\end{tikzpicture}")
end


function write_wrapper(io, det::DetectorV2; profile::String="core")
    println(io, "\\documentclass{article}")
    println(io, "\\usepackage[margin=1.5cm]{geometry}")
    println(io, "\\usepackage{tikz}")
    println(io, "\\usepackage{xcolor}")
    println(io, "\\usepackage{graphicx}")
    println(io, "\\pagestyle{empty}")
    println(io, "\\begin{document}")
    println(io, "\\begin{figure}[p]")
    println(io, "\\centering")
    println(io, "\\resizebox{0.45\\linewidth}{!}{%")
    write_tikz(io, det; profile=profile)
    println(io, "}")
    println(io, "\\end{figure}")
    println(io, "\\end{document}")
end


function write_tikz_file(path::String, det::DetectorV2; profile::String="core", mode::String="wrapper")
    open(path, "w") do io
        println(io, "% Auto-generated by scripts/inspect_geometry_v2.jl")
        if mode == "standalone"
            println(io, "\\documentclass[border=2pt,convert={density=300,outext=.png}]{standalone}")
            println(io, "\\usepackage{tikz}")
            println(io, "\\usepackage{xcolor}")
            println(io, "\\begin{document}")
            write_tikz(io, det; profile=profile)
            println(io, "\\end{document}")
        elseif mode == "wrapper"
            write_wrapper(io, det; profile=profile)
        elseif mode == "bare"
            write_tikz(io, det; profile=profile)
        else
            error("Unknown latex mode '$mode'")
        end
    end
end


function print_overlaps(det::DetectorV2)
    exact = sibling_overlaps(det; exact_only=true)
    allov = sibling_overlaps(det; exact_only=false)

    println("=" ^ 72)
    println("Exact-exact sibling overlaps")
    println("=" ^ 72)
    if isempty(exact)
        println("none")
    else
        for ov in exact
            println("$(ov.parent): $(ov.a) <-> $(ov.b)")
        end
    end

    println()
    println("=" ^ 72)
    println("All sibling overlaps including proxy geometry")
    println("=" ^ 72)
    if isempty(allov)
        println("none")
    else
        for ov in allov
            println("$(ov.parent): $(ov.a) <-> $(ov.b) [$(ov.a_approximation) / $(ov.b_approximation)]")
        end
    end
end


function main()
    args = parse_args()
    det_path = isempty(args["detector"]) ? default_detector_v2_path() : args["detector"]
    tikz_profile = lowercase(strip(args["tikz-profile"]))
    latex_mode = args["standalone"] ? "standalone" : lowercase(strip(args["latex-mode"]))

    cfg = default_config()
    mats = load_materials(cfg)
    det = load_detector_v2(det_path, mats; validate=true)

    println(detector_summary(det))
    println()

    if !args["no-tree"]
        println(tree_dump(det))
        println()
    end

    print_overlaps(det)

    tikz_out = args["tikz-out"]
    if !isempty(tikz_out)
        write_tikz_file(tikz_out, det; profile=tikz_profile, mode=latex_mode)
        println()
        println("TikZ schematic written to: $tikz_out")
    end
end


main()
