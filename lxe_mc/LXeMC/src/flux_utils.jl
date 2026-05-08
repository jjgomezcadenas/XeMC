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
