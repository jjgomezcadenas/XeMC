"""
Source dispatch: maps (source_name, isotope) pairs to the appropriate
flux generator function. Extensible for future source types.
"""


"""
    dispatch_source_flux(source, isotope, N, sg, cfg, rng) -> NamedTuple

Generate flux tables for one (source, isotope) pair.
Returns a named tuple whose keys are component symbols
(e.g., :bi214_ocv, :bi214_icv, :bi214_rate).

Supported sources:
- `"cryostat_barrel"`, `"cryostat_top"`, `"cryostat_bottom"`

Supported isotopes:
- `"Bi214"`, `"Tl208"`
"""
function dispatch_source_flux(source::String, isotope::String,
                               N::Int, sg::Dict{String,SourceVolumeInfo},
                               cfg::SimConfig, rng::AbstractRNG;
                               verbose::Bool=false)
    if source == "cryostat_barrel"
        result = cryostat_barrel_flux(N, sg, cfg, rng; verbose=verbose)
    elseif source == "cryostat_top"
        result = cryostat_top_flux(N, sg, cfg, rng; verbose=verbose)
    elseif source == "cryostat_bottom"
        result = cryostat_bottom_flux(N, sg, cfg, rng; verbose=verbose)
    else
        error("Unknown source '$source'. Supported: cryostat_barrel, cryostat_top, cryostat_bottom")
    end

    # Extract only the requested isotope's tables
    if isotope == "Bi214"
        _extract_bi214(source, result)
    elseif isotope == "Tl208"
        _extract_tl208(source, result)
    else
        error("Unknown isotope '$isotope'. Supported: Bi214, Tl208")
    end
end


function _extract_bi214(source::String, result)
    if source == "cryostat_barrel"
        Dict{Symbol,Any}(
            :bi214_ocv => result.bi214_ocv,
            :bi214_mli => result.bi214_mli,
            :bi214_icv => result.bi214_icv,
            :bi214_rate => result.bi214_rate
        )
    else
        Dict{Symbol,Any}(
            :bi214_ocv => result.bi214_ocv,
            :bi214_icv => result.bi214_icv,
            :bi214_rate => result.bi214_rate
        )
    end
end


function _extract_tl208(source::String, result)
    if source == "cryostat_barrel"
        Dict{Symbol,Any}(
            :tl208_ocv => result.tl208_ocv,
            :tl208_mli => result.tl208_mli,
            :tl208_icv => result.tl208_icv,
            :tl208_rate => result.tl208_rate
        )
    else
        Dict{Symbol,Any}(
            :tl208_ocv => result.tl208_ocv,
            :tl208_icv => result.tl208_icv,
            :tl208_rate => result.tl208_rate
        )
    end
end


"""
    merge_dispatch_results(isotope, results::Vector{Dict{Symbol,Any}}) -> Dict{Symbol,Any}

Merge per-thread dispatch results. Each element in `results` is a Dict
from `dispatch_source_flux`. Component flux tables are merged via
`merge_flux_bi214` or `merge_flux_tl208`. Rate tables are recomputed
from merged components (not merged directly).
"""
function merge_dispatch_results(isotope::String,
                                 results::Vector{Dict{Symbol,Any}})::Dict{Symbol,Any}
    length(results) == 1 && return results[1]

    merged = Dict{Symbol,Any}()
    keys_all = collect(keys(results[1]))

    for key in keys_all
        key == :bi214_rate && continue
        key == :tl208_rate && continue

        tables = [r[key] for r in results]
        if isotope == "Bi214"
            merged[key] = merge_flux_bi214(Vector{SourceFluxBi214}(tables))
        else
            merged[key] = merge_flux_tl208(Vector{SourceFluxTl208}(tables))
        end
    end

    # Rate table: take from first result (will be recomputed by caller if needed)
    rate_key = isotope == "Bi214" ? :bi214_rate : :tl208_rate
    if haskey(results[1], rate_key)
        merged[rate_key] = results[1][rate_key]
    end

    merged
end


"""
    supported_sources() -> Vector{String}

List of currently supported source identifiers.
"""
function supported_sources()::Vector{String}
    ["cryostat_barrel", "cryostat_top", "cryostat_bottom"]
end


"""
    supported_isotopes() -> Vector{String}

List of currently supported isotope identifiers.
"""
function supported_isotopes()::Vector{String}
    ["Bi214", "Tl208"]
end


"""
    make_surface_sampler(source, sg) -> Function

Return a surface sampler function `rng -> (position, normal)` for the
given source. The sampler places gammas on the ICV inner surface
appropriate for the source geometry.
"""
function make_surface_sampler(source::String,
                               sg::Dict{String,SourceVolumeInfo})::Function
    if source == "cryostat_barrel"
        icv = sg["ICV_barrel"]
        s = icv.volume.logical.solid
        c = icv.volume.logical.position
        R_inner = s.R_inner_cm
        z_min = c[3] - s.half_height_cm
        z_max = c[3] + s.half_height_cm
        return rng -> sample_barrel_point(R_inner, z_min, z_max, rng)
    elseif source == "cryostat_top"
        icv = sg["ICV_top"]
        s = icv.volume.logical.solid
        c = icv.volume.logical.position
        return rng -> sample_cap_point(s.radius_cm, s.aspect_ratio,
                                        c[3], :up, rng)
    elseif source == "cryostat_bottom"
        icv = sg["ICV_bottom"]
        s = icv.volume.logical.solid
        c = icv.volume.logical.position
        return rng -> sample_cap_point(s.radius_cm, s.aspect_ratio,
                                        c[3], :down, rng)
    else
        error("No surface sampler for source '$source'")
    end
end
