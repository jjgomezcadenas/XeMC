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
    make_surface_sampler(source, fk) -> Function

Return a surface sampler function `rng -> (position, normal)` for the
given source. The sampler places gammas on a virtual envelope that sits
just inside the tracking detector boundary, extracted from the compiled
`FastKernelGeometry`.

- **Barrel**: cylindrical shell at R = detector rmax, z from Skin region.
- **Top**: ellipsoidal cap matching the AirDome geometry.
- **Bottom**: ellipsoidal cap matching the LXe_passive bottom cap.
"""
function make_surface_sampler(source::String,
                               fk::FastKernelGeometry)::Function
    if source == "cryostat_barrel"
        # Barrel envelope: LZ_detector barrel region (between top and bottom caps)
        lz = fk.regions[fk.name_to_index["LZ_detector"]]
        top_depth = lz.has_top_cap ? lz.top_cap_radius_cm / lz.top_cap_aspect_ratio : 0.0
        bot_depth = lz.has_bottom_cap ? lz.bottom_cap_radius_cm / lz.bottom_cap_aspect_ratio : 0.0
        R = lz.rmax_cm - 0.01     # slightly inside boundary (82.09 cm)
        z_min = lz.zmin_cm + bot_depth  # top of bottom cap
        z_max = lz.zmax_cm - top_depth  # bottom of top cap
        return rng -> sample_barrel_point(R, z_min, z_max, rng)
    elseif source == "cryostat_top"
        # Top envelope: AirDome cap, slightly inside boundary
        air = fk.regions[fk.name_to_index["AirDome"]]
        R = air.top_cap_radius_cm - 0.01   # slightly inside (82.09 cm)
        ar = air.top_cap_aspect_ratio      # 2.0
        z_eq = air.zmin_cm                 # 145.6 cm (equator of the cap)
        return rng -> sample_cap_point(R, ar, z_eq, :up, rng)
    elseif source == "cryostat_bottom"
        # Bottom envelope: LXe_passive bottom cap, slightly inside boundary
        passive = fk.regions[fk.name_to_index["LXe_passive"]]
        R = passive.bottom_cap_radius_cm - 0.01   # slightly inside (82.09 cm)
        ar = passive.bottom_cap_aspect_ratio       # 3.0
        cap_depth = passive.bottom_cap_radius_cm / ar  # use full R for depth calc
        z_eq = passive.zmin_cm + cap_depth
        return rng -> sample_cap_point(R, ar, z_eq, :down, rng)
    else
        error("No surface sampler for source '$source'")
    end
end
