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
    VirtualEnvelope

Surface where sampled gammas are placed. For dense (cryostat) sources
the envelope sits just inside the detector boundary. For ideal
(transparent) sources the envelope is the source geometry itself.

Dense envelopes carry a 0.01 cm inset from the detector boundary to
avoid floating-point boundary issues with the fast-kernel classifier.

Fields:
- `kind`: `:barrel`, `:cap_up`, `:cap_down`, or `:disk_flat`
- `R_cm`: radius [cm]
- `z_min_cm`, `z_max_cm`: z range (barrel only; 0 for caps and disks)
- `z_equator_cm`: equator z (caps); disk z (`:disk_flat`); 0 for barrel
- `aspect_ratio`: cap aspect ratio (caps only; 0 for barrel and disk_flat)
"""
struct VirtualEnvelope
    kind::Symbol
    R_cm::Float64
    z_min_cm::Float64
    z_max_cm::Float64
    z_equator_cm::Float64
    aspect_ratio::Float64
end

const ENVELOPE_INSET_CM = 0.01  # inset from detector boundary


"""
    make_virtual_envelope(source, fk) -> VirtualEnvelope

Build the virtual envelope for a source from the tracking detector
geometry. The envelope shape matches the detector boundary for the
relevant surface.

- **Barrel**: cylindrical shell from the LZ_detector barrel region.
- **Top**: ellipsoidal cap from the AirDome.
- **Bottom**: ellipsoidal cap from the LXe_passive bottom.
"""
function make_virtual_envelope(source::String,
                                fk::FastKernelGeometry)::VirtualEnvelope
    if source == "cryostat_barrel"
        lz = fk.regions[fk.name_to_index["LZ_detector"]]
        top_depth = lz.has_top_cap ? lz.top_cap_radius_cm / lz.top_cap_aspect_ratio : 0.0
        bot_depth = lz.has_bottom_cap ? lz.bottom_cap_radius_cm / lz.bottom_cap_aspect_ratio : 0.0
        VirtualEnvelope(
            :barrel,
            lz.rmax_cm - ENVELOPE_INSET_CM,
            lz.zmin_cm + bot_depth,   # top of bottom cap
            lz.zmax_cm - top_depth,   # bottom of top cap
            0.0,
            0.0
        )
    elseif source == "cryostat_top"
        air = fk.regions[fk.name_to_index["AirDome"]]
        VirtualEnvelope(
            :cap_up,
            air.top_cap_radius_cm - ENVELOPE_INSET_CM,
            0.0,
            0.0,
            air.zmin_cm,              # equator of the AirDome cap
            air.top_cap_aspect_ratio
        )
    elseif source == "cryostat_bottom"
        passive = fk.regions[fk.name_to_index["LXe_passive"]]
        cap_depth = passive.bottom_cap_radius_cm / passive.bottom_cap_aspect_ratio
        VirtualEnvelope(
            :cap_down,
            passive.bottom_cap_radius_cm - ENVELOPE_INSET_CM,
            0.0,
            0.0,
            passive.zmin_cm + cap_depth,   # equator of the bottom cap
            passive.bottom_cap_aspect_ratio
        )
    else
        error("No virtual envelope for source '$source'")
    end
end


"""
    make_virtual_envelope(sv::SourceVolumeInfo) -> VirtualEnvelope

Build a virtual envelope directly from a source volume (ideal sources).
The envelope is the source geometry itself:
- `PCyl` (cylinder) → `:disk_flat` at the bottom face of the volume
- `PCylShell` → `:barrel` spanning the shell z range
"""
function make_virtual_envelope(sv::SourceVolumeInfo)::VirtualEnvelope
    vol = sv.volume
    if vol isa PCyl
        lv = vol.logical
        R = lv.solid.radius_cm
        z_bottom = lv.position[3] - lv.solid.half_height_cm
        VirtualEnvelope(:disk_flat, R, 0.0, 0.0, z_bottom, 0.0)
    elseif vol isa PCylShell
        lv = vol.logical
        R = lv.solid.R_inner_cm
        z_center = lv.position[3]
        hh = lv.solid.half_height_cm
        VirtualEnvelope(:barrel, R, z_center - hh, z_center + hh, 0.0, 0.0)
    else
        error("Cannot build VE from source volume type $(typeof(vol))")
    end
end


"""
    make_surface_sampler(env::VirtualEnvelope) -> Function

Return a surface sampler function `rng -> (position, normal)` from a
`VirtualEnvelope`.
"""
function make_surface_sampler(env::VirtualEnvelope)::Function
    if env.kind === :barrel
        return rng -> sample_barrel_point(env.R_cm, env.z_min_cm, env.z_max_cm, rng)
    elseif env.kind === :cap_up
        return rng -> sample_cap_point(env.R_cm, env.aspect_ratio,
                                        env.z_equator_cm, :up, rng)
    elseif env.kind === :cap_down
        return rng -> sample_cap_point(env.R_cm, env.aspect_ratio,
                                        env.z_equator_cm, :down, rng)
    elseif env.kind === :disk_flat
        return rng -> sample_disk_point(env.R_cm, env.z_equator_cm, rng)
    else
        error("Unknown envelope kind '$(env.kind)'")
    end
end


"""
    make_surface_sampler(source, fk) -> Function

Convenience: build virtual envelope and return sampler in one call.
"""
function make_surface_sampler(source::String,
                               fk::FastKernelGeometry)::Function
    env = make_virtual_envelope(source, fk)
    make_surface_sampler(env)
end
