"""
Utility functions for source flux tables: merging per-thread results
and JSON serialization/deserialization.
"""


# =====================================================================
# Merging per-thread flux tables
# =====================================================================

"""
    merge_flux_bi214(tables::Vector{SourceFluxBi214}) -> SourceFluxBi214

Merge multiple Bi-214 flux tables (e.g., from parallel threads) into one.
PDFs are weight-averaged by N_generated; counters are summed.
"""
function merge_flux_bi214(tables::Vector{SourceFluxBi214})::SourceFluxBi214
    length(tables) == 1 && return tables[1]
    N_total = sum(t.N_generated for t in tables)
    pdf = sum(t.pdf .* t.N_generated for t in tables) ./ N_total
    SourceFluxBi214(
        tables[1].source_name, pdf,
        tables[1].E_min, tables[1].E_max, tables[1].n_E, tables[1].n_u,
        N_total,
        sum(t.N_surviving for t in tables),
        sum(t.N_absorbed for t in tables),
        sum(t.N_backward for t in tables),
        sum(t.N_low_energy for t in tables)
    )
end


"""
    merge_flux_tl208(tables::Vector{SourceFluxTl208}) -> SourceFluxTl208

Merge multiple Tl-208 flux tables into one.
"""
function merge_flux_tl208(tables::Vector{SourceFluxTl208})::SourceFluxTl208
    length(tables) == 1 && return tables[1]
    N_total = sum(t.N_generated for t in tables)
    pdf_main = sum(t.pdf_main .* t.N_generated for t in tables) ./ N_total

    n_comp = length(tables[1].companion_E_line)
    pdf_comp = [sum(t.pdf_companion[ic] .* t.N_generated for t in tables) ./ N_total
                for ic in 1:n_comp]

    comp_f = [sum(t.companion_f[ic] * t.N_generated for t in tables) / N_total
              for ic in 1:n_comp]

    SourceFluxTl208(
        tables[1].source_name,
        pdf_main, tables[1].E_min_main, tables[1].E_max_main, tables[1].n_E_main,
        pdf_comp,
        tables[1].companion_E_line,
        tables[1].companion_BR,
        comp_f,
        tables[1].E_min_companion,
        tables[1].E_max_companion,
        tables[1].n_E_companion,
        tables[1].n_u,
        N_total,
        sum(t.N_surviving_main for t in tables),
        sum(t.N_absorbed_main for t in tables),
        sum(t.N_backward_main for t in tables),
        sum(t.N_low_energy_main for t in tables)
    )
end


# =====================================================================
# CSV output for PDF matrices
# =====================================================================

"""
    write_pdf_csv(path, pdf, E_min, E_max, n_E, n_u)

Write a 2D PDF matrix to CSV. First row is a header with u bin centers.
First column of each data row is the E bin center. Remaining columns
are the PDF values.
"""
function write_pdf_csv(path::String, pdf::Matrix{Float64},
                        E_min::Float64, E_max::Float64,
                        n_E::Int, n_u::Int)
    dE = (E_max - E_min) / n_E
    du = 1.0 / n_u
    open(path, "w") do io
        # Header: E_center, u_center_1, u_center_2, ...
        print(io, "E_MeV")
        for j in 1:n_u
            @printf(io, ",u=%.4f", (j - 0.5) * du)
        end
        println(io)
        # Data rows
        for i in 1:n_E
            @printf(io, "%.6f", E_min + (i - 0.5) * dE)
            for j in 1:n_u
                @printf(io, ",%.8e", pdf[i, j])
            end
            println(io)
        end
    end
end


"""Write all CSV files for a Bi-214 flux table under `dir` with prefix `prefix`."""
function write_flux_bi214_csv(dir::String, prefix::String, ft::SourceFluxBi214)
    write_pdf_csv(joinpath(dir, prefix * ".csv"),
                  ft.pdf, ft.E_min, ft.E_max, ft.n_E, ft.n_u)
end


"""Write all CSV files for a Tl-208 flux table under `dir` with prefix `prefix`."""
function write_flux_tl208_csv(dir::String, prefix::String, ft::SourceFluxTl208)
    write_pdf_csv(joinpath(dir, prefix * "_main.csv"),
                  ft.pdf_main, ft.E_min_main, ft.E_max_main, ft.n_E_main, ft.n_u)
    for ic in 1:length(ft.companion_E_line)
        E_line_str = @sprintf("%.0fkeV", ft.companion_E_line[ic] * 1000)
        write_pdf_csv(joinpath(dir, prefix * "_companion_" * E_line_str * ".csv"),
                      ft.pdf_companion[ic],
                      ft.E_min_companion[ic], ft.E_max_companion[ic],
                      ft.n_E_companion[ic], ft.n_u)
    end
end


"""Write a rate table as CSV."""
function write_rate_table_csv(dir::String, prefix::String, rt::SourceRateTable)
    write_pdf_csv(joinpath(dir, prefix * ".csv"),
                  rt.pdf_rate, rt.E_min, rt.E_max, rt.n_E, rt.n_u)
end


# =====================================================================
# JSON for metadata only
# =====================================================================

"""Write a Dict as pretty-printed JSON."""
function write_flux_json(path::String, data)
    open(path, "w") do io
        JSON.print(io, data, 2)
    end
end


"""Build metadata dict for a Bi-214 flux table."""
function flux_bi214_metadata(ft::SourceFluxBi214)::Dict{String,Any}
    Dict{String,Any}(
        "type" => "SourceFluxBi214",
        "source_name" => ft.source_name,
        "E_min" => ft.E_min,
        "E_max" => ft.E_max,
        "n_E" => ft.n_E,
        "n_u" => ft.n_u,
        "N_generated" => ft.N_generated,
        "N_surviving" => ft.N_surviving,
        "N_absorbed" => ft.N_absorbed,
        "N_backward" => ft.N_backward,
        "N_low_energy" => ft.N_low_energy,
        "survival_fraction" => sum(ft.pdf)
    )
end


"""Build metadata dict for a Tl-208 flux table."""
function flux_tl208_metadata(ft::SourceFluxTl208)::Dict{String,Any}
    Dict{String,Any}(
        "type" => "SourceFluxTl208",
        "source_name" => ft.source_name,
        "E_min_main" => ft.E_min_main,
        "E_max_main" => ft.E_max_main,
        "n_E_main" => ft.n_E_main,
        "n_u" => ft.n_u,
        "companion_E_line" => ft.companion_E_line,
        "companion_BR" => ft.companion_BR,
        "companion_f" => ft.companion_f,
        "E_min_companion" => ft.E_min_companion,
        "E_max_companion" => ft.E_max_companion,
        "n_E_companion" => ft.n_E_companion,
        "N_generated" => ft.N_generated,
        "N_surviving_main" => ft.N_surviving_main,
        "N_absorbed_main" => ft.N_absorbed_main,
        "N_backward_main" => ft.N_backward_main,
        "N_low_energy_main" => ft.N_low_energy_main,
        "survival_fraction_main" => sum(ft.pdf_main)
    )
end


# =====================================================================
# Loading flux tables from CSV + metadata
# =====================================================================

"""
    load_pdf_csv(path) -> (E_centers, u_centers, pdf)

Read a PDF CSV file back into arrays.
"""
function load_pdf_csv(path::String)::Tuple{Vector{Float64},Vector{Float64},Matrix{Float64}}
    lines = readlines(path)
    header = split(lines[1], ",")
    n_u = length(header) - 1
    u_centers = [parse(Float64, split(h, "=")[2]) for h in header[2:end]]

    n_E = length(lines) - 1
    E_centers = Vector{Float64}(undef, n_E)
    pdf = Matrix{Float64}(undef, n_E, n_u)
    for i in 1:n_E
        cols = split(lines[i+1], ",")
        E_centers[i] = parse(Float64, cols[1])
        for j in 1:n_u
            pdf[i, j] = parse(Float64, cols[j+1])
        end
    end
    (E_centers, u_centers, pdf)
end


"""
    load_flux_bi214(dir, key) -> SourceFluxBi214

Load a Bi-214 flux table from `dir`. `key` is the component name
(e.g., `"bi214_ocv"`). Reads `<key>.csv` and metadata from `metadata.json`.
"""
function load_flux_bi214(dir::String, key::String)::SourceFluxBi214
    _, _, pdf = load_pdf_csv(joinpath(dir, key * ".csv"))
    meta = open(joinpath(dir, "metadata.json"), "r") do io
        JSON.parse(io)
    end
    m = meta["components"][key]
    SourceFluxBi214(
        String(m["source_name"]),
        pdf,
        Float64(m["E_min"]), Float64(m["E_max"]),
        Int(m["n_E"]), Int(m["n_u"]),
        Int(m["N_generated"]), Int(m["N_surviving"]),
        Int(m["N_absorbed"]), Int(m["N_backward"]),
        Int(m["N_low_energy"])
    )
end


"""
    load_flux_tl208(dir, key) -> SourceFluxTl208

Load a Tl-208 flux table from `dir`. `key` is the component name
(e.g., `"tl208_ocv"`). Reads `<key>_main.csv`, companion CSVs,
and metadata from `metadata.json`.
"""
function load_flux_tl208(dir::String, key::String)::SourceFluxTl208
    meta = open(joinpath(dir, "metadata.json"), "r") do io
        JSON.parse(io)
    end
    m = meta["components"][key]

    _, _, pdf_main = load_pdf_csv(joinpath(dir, key * "_main.csv"))

    companion_E_line = Float64.(m["companion_E_line"])
    companion_BR = Float64.(m["companion_BR"])
    companion_f = Float64.(m["companion_f"])
    E_min_comp = Float64.(m["E_min_companion"])
    E_max_comp = Float64.(m["E_max_companion"])
    n_E_comp = Int.(m["n_E_companion"])
    n_u = Int(m["n_u"])

    n_comp = length(companion_E_line)
    pdf_companion = Matrix{Float64}[]
    for ic in 1:n_comp
        E_line_str = @sprintf("%.0fkeV", companion_E_line[ic] * 1000)
        _, _, pdf_c = load_pdf_csv(joinpath(dir, key * "_companion_" * E_line_str * ".csv"))
        push!(pdf_companion, pdf_c)
    end

    SourceFluxTl208(
        String(m["source_name"]),
        pdf_main,
        Float64(m["E_min_main"]), Float64(m["E_max_main"]),
        Int(m["n_E_main"]),
        pdf_companion,
        companion_E_line,
        companion_BR,
        companion_f,
        E_min_comp,
        E_max_comp,
        n_E_comp,
        n_u,
        Int(m["N_generated"]),
        Int(m["N_surviving_main"]),
        Int(m["N_absorbed_main"]),
        Int(m["N_backward_main"]),
        Int(m["N_low_energy_main"])
    )
end


"""
    load_rate_table(dir, key) -> SourceRateTable

Load a rate table from `dir`. `key` is the component name
(e.g., `"bi214_rate"`). Reads `<key>.csv` and metadata from `metadata.json`.
"""
function load_rate_table(dir::String, key::String)::SourceRateTable
    _, _, pdf_rate = load_pdf_csv(joinpath(dir, key * ".csv"))
    meta = open(joinpath(dir, "metadata.json"), "r") do io
        JSON.parse(io)
    end
    m = meta["components"][key]
    SourceRateTable(
        Symbol(m["surface"]),
        pdf_rate,
        Float64(m["E_min"]), Float64(m["E_max"]),
        Int(m["n_E"]), Int(m["n_u"]),
        String.(m["component_names"]),
        Float64.(m["component_rates"]),
        Float64(m["total_rate"])
    )
end


"""Build metadata dict for a rate table."""
function rate_table_metadata(rt::SourceRateTable)::Dict{String,Any}
    Dict{String,Any}(
        "type" => "SourceRateTable",
        "surface" => String(rt.surface),
        "E_min" => rt.E_min,
        "E_max" => rt.E_max,
        "n_E" => rt.n_E,
        "n_u" => rt.n_u,
        "component_names" => rt.component_names,
        "component_rates" => rt.component_rates,
        "total_rate" => rt.total_rate
    )
end
