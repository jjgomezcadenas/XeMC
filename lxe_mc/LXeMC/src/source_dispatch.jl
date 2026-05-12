"""
Source dispatch: maps (source_name, isotope) pairs to the appropriate
flux generator function. Extensible for future source types.
"""


"""
    dispatch_source_flux(source, isotope, N, sg, cfg, rng; mats, verbose) -> Dict

Generate flux tables for one (source, isotope) pair.
Returns a Dict whose keys are component symbols
(e.g., :bi214_ocv, :bi214_icv, :bi214_rate).

Supported sources:
- `"cryostat_barrel"`, `"cryostat_top"`, `"cryostat_bottom"`
- `"pmt_top"`, `"pmt_bottom"`, `"pmt_bottom_lxe"`

Supported isotopes:
- `"Bi214"`, `"Tl208"`

`mats` and `det` are required only for `"pmt_bottom_lxe"`.
"""
function dispatch_source_flux(source::String, isotope::String,
                               N::Int, sg::Dict{String,SourceVolumeInfo},
                               cfg::SimConfig, rng::AbstractRNG;
                               mats::Union{Dict{String,Material},Nothing}=nothing,
                               verbose::Bool=false)
    if source == "cryostat_barrel"
        result = cryostat_barrel_flux(N, sg, cfg, rng; verbose=verbose)
    elseif source == "cryostat_top"
        result = cryostat_top_flux(N, sg, cfg, rng; verbose=verbose)
    elseif source == "cryostat_bottom"
        result = cryostat_bottom_flux(N, sg, cfg, rng; verbose=verbose)
    elseif source == "pmt_top"
        result = pmt_top_flux(N, sg, cfg, rng; verbose=verbose)
    elseif source == "pmt_bottom"
        result = pmt_bottom_flux(N, sg, cfg, rng; verbose=verbose)
    elseif source == "pmt_bottom_lxe"
        mats === nothing && error("pmt_bottom_lxe requires mats keyword argument")
        result = pmt_bottom_lxe_flux(N, sg, mats, cfg, rng; verbose=verbose)
    else
        error("Unknown source '$source'. Supported: $(join(supported_sources(), ", "))")
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
    elseif startswith(source, "pmt_")
        Dict{Symbol,Any}(
            :bi214 => result.bi214,
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
    elseif startswith(source, "pmt_")
        Dict{Symbol,Any}(
            :tl208 => result.tl208,
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
    ["cryostat_barrel", "cryostat_top", "cryostat_bottom", "pmt_top", "pmt_bottom", "pmt_bottom_lxe"]
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

Surface where sampled gammas are placed, built from the source geometry.
For cryostat sources the envelope is the ICV inner surface; for PMT
sources it is the source volume itself. All envelopes carry a tiny
inset (`ENVELOPE_INSET_CM = 0.01 cm`) to avoid floating-point boundary
issues with the fast-kernel region classifier.

Fields:
- `kind`: `:barrel`, `:cap_up`, `:cap_down`, or `:disk_flat`
- `R_cm`: radius [cm]
- `z_min_cm`, `z_max_cm`: z range (barrel only; 0 for caps and disks)
- `z_equator_cm`: equator z (caps); disk z (`:disk_flat`); 0 for barrel
- `aspect_ratio`: cap aspect ratio (caps only; normal z-sign for `:disk_flat`; 0 for barrel)
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
    make_virtual_envelope(source, sg) -> VirtualEnvelope

Build the virtual envelope for any source from the source geometry.
The VE sits on the ICV inner surface (for cryostat sources) or on
the source itself (for PMT sources), with a small inset
(`ENVELOPE_INSET_CM`) to avoid floating-point boundary issues with
the fast-kernel region classifier.

- **cryostat_barrel**: ICV_barrel inner cylinder
- **cryostat_top**: ICV_top inner dome (:cap_up)
- **cryostat_bottom**: ICV_bottom inner dome (:cap_down)
- **pmt_top**: merged PMT disk (:disk_flat, downward)
- **pmt_bottom**: merged PMT disk (:disk_flat, upward)
"""
function make_virtual_envelope(source::String,
                                sg::Dict{String,SourceVolumeInfo})::VirtualEnvelope
    ε = ENVELOPE_INSET_CM

    if source == "cryostat_barrel"
        icv = sg["ICV_barrel"]
        R = icv.volume.logical.solid.R_inner_cm - ε
        z_eq_top = sg["ICV_top"].volume.logical.position[3]
        z_eq_bot = sg["ICV_bottom"].volume.logical.position[3]
        # Barrel z excludes dome regions
        top_depth = R / sg["ICV_top"].volume.logical.solid.aspect_ratio
        bot_depth = R / sg["ICV_bottom"].volume.logical.solid.aspect_ratio
        VirtualEnvelope(:barrel, R, z_eq_bot + bot_depth, z_eq_top - top_depth, 0.0, 0.0)

    elseif source == "cryostat_top"
        icv = sg["ICV_top"]
        R = icv.volume.logical.solid.radius_cm - icv.volume.logical.solid.wall_thickness_cm - ε
        z_eq = icv.volume.logical.position[3]
        ar = icv.volume.logical.solid.aspect_ratio
        VirtualEnvelope(:cap_up, R, 0.0, 0.0, z_eq, ar)

    elseif source == "cryostat_bottom"
        icv = sg["ICV_bottom"]
        R = icv.volume.logical.solid.radius_cm - icv.volume.logical.solid.wall_thickness_cm - ε
        z_eq = icv.volume.logical.position[3]
        ar = icv.volume.logical.solid.aspect_ratio
        VirtualEnvelope(:cap_down, R, 0.0, 0.0, z_eq, ar)

    elseif source == "pmt_top"
        merged, _ = _merge_pmt_volume(sg,
            ["PMT_TOP_PMTs", "PMT_TOP_bases", "PMT_TOP_structure"],
            "PMT_TOP_merged", :up)
        sv = SourceVolumeInfo("PMT_TOP_merged", merged, merged.material,
            Dict{String,Float64}(), 0.0, :transparent, "virtual_source", "merged")
        make_virtual_envelope(sv)

    elseif source == "pmt_bottom"
        merged, _ = _merge_pmt_volume(sg,
            ["PMT_BOT_PMTs", "PMT_BOT_bases", "PMT_BOT_structure", "PMT_BOT_R8778_dome"],
            "PMT_BOT_merged", :down)
        sv = SourceVolumeInfo("PMT_BOT_merged", merged, merged.material,
            Dict{String,Float64}(), 0.0, :transparent, "virtual_source", "merged")
        make_virtual_envelope(sv)

    elseif source == "pmt_bottom_lxe"
        # VE just above the cathode: flat disk at z_cathode + ε, gammas
        # going upward. The cathode plane (z=0) is the boundary between
        # LXe_below_cathode and BottomActive; with strict-inequality
        # containment, neither region claims z=0 exactly. Insetting by ε
        # places the launch point unambiguously inside BottomActive.
        z_cathode = _cathode_z(sg)
        merged, _ = _merge_pmt_volume(sg,
            ["PMT_BOT_PMTs", "PMT_BOT_bases", "PMT_BOT_structure", "PMT_BOT_R8778_dome"],
            "PMT_BOT_merged", :down)
        R = merged.logical.solid.radius_cm
        VirtualEnvelope(:disk_flat, R, 0.0, 0.0, z_cathode + ε, 1.0)

    else
        error("Unknown source '$source'. Supported: $(join(supported_sources(), ", "))")
    end
end


"""
    make_virtual_envelope(sv::SourceVolumeInfo) -> VirtualEnvelope

Build a virtual envelope directly from a source volume (ideal sources).
The envelope is the source geometry itself:
- `PCyl` (cylinder) → `:disk_flat` at the bottom face, normal_z = -1
- `PDisk` (flat disk) → `:disk_flat` with normal_z from orientation
  (:up → -1.0, toward LXe below; :down → +1.0, toward LXe above)
- `PCylShell` → `:barrel` spanning the shell z range

For `:disk_flat`, the `aspect_ratio` field stores the normal z-sign
(-1.0 or +1.0) used by `sample_disk_point`.
"""
function make_virtual_envelope(sv::SourceVolumeInfo)::VirtualEnvelope
    vol = sv.volume
    if vol isa PCyl
        lv = vol.logical
        R = lv.solid.radius_cm
        z_bottom = lv.position[3] - lv.solid.half_height_cm
        VirtualEnvelope(:disk_flat, R, 0.0, 0.0, z_bottom, -1.0)
    elseif vol isa PDisk
        lv = vol.logical
        R = lv.solid.radius_cm
        normal_z = lv.orientation === :up ? -1.0 : 1.0
        # VE at bottom face for :up, top face for :down
        if lv.orientation === :up
            z_ve = lv.position[3]  # bottom face
        else
            z_ve = lv.position[3]  # top face (position is at z_max for :down)
        end
        VirtualEnvelope(:disk_flat, R, 0.0, 0.0, z_ve, normal_z)
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
        nz = env.aspect_ratio  # normal z-sign stored in aspect_ratio field
        return rng -> sample_disk_point(env.R_cm, env.z_equator_cm, nz, rng)
    else
        error("Unknown envelope kind '$(env.kind)'")
    end
end


"""
    make_surface_sampler(source, sg) -> Function

Convenience: build virtual envelope from source geometry and return
sampler in one call. Works for all source types.
"""
function make_surface_sampler(source::String,
                               sg::Dict{String,SourceVolumeInfo})::Function
    env = make_virtual_envelope(source, sg)
    make_surface_sampler(env)
end
