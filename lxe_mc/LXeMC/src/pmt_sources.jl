"""
PMT flux functions: top and bottom.

Each function builds a merged transparent volume spanning the full z extent
of all PMT sub-components (PMTs, bases, structure, dome) from the source
geometry, generates flux tables via the standard vacuum propagation path,
and returns per-component rate breakdowns by bookkeeping (mass × activity).

The merged volume is a flat PDisk (aspect_ratio = Inf) with orientation
that determines the "toward LXe" hemisphere:
- TOP PMTs: orientation = :up  → cos_theta_to_lxe = -dir_z (downward = toward LXe)
- BOT PMTs: orientation = :down → cos_theta_to_lxe = +dir_z (upward = toward LXe)
"""


"""
    _merge_pmt_volume(sg, names, merged_name, orientation, materials)
        -> (PDisk, Vector{SourceVolumeInfo})

Build a merged flat disk volume spanning the full z extent and common
radius of the named PMT sub-components. The merged volume is vacuum
(transparent) with the given orientation for hemisphere filtering.

The z span is [min(z - hh), max(z + hh)] over all components — no
hardwired numbers, always computed from the source geometry.
"""
function _merge_pmt_volume(sg::Dict{String,SourceVolumeInfo},
                            names::Vector{String},
                            merged_name::String,
                            orientation::Symbol,
                            materials::Dict{String,Material})::Tuple{PDisk,Vector{SourceVolumeInfo}}
    components = [sg[n] for n in names]

    # Compute z span from all components (maximum extent)
    z_min = Inf
    z_max = -Inf
    R = 0.0
    for sv in components
        vol = sv.volume
        if vol isa PCyl
            lv = vol.logical
            z_lo = lv.position[3] - lv.solid.half_height_cm
            z_hi = lv.position[3] + lv.solid.half_height_cm
            R = lv.solid.radius_cm
        else
            error("PMT component '$(sv.name)' is not a PCyl")
        end
        z_min = min(z_min, z_lo)
        z_max = max(z_max, z_hi)
    end

    thickness = z_max - z_min
    # Position at z_min for :up (dome faces up, LXe below),
    # at z_max for :down (dome faces down, LXe above)
    z_pos = orientation === :up ? z_min : z_max
    solid = Disk(R, thickness, Inf)  # flat disk (aspect_ratio = Inf)
    vac = materials["Vacuum"]
    merged = PDisk(merged_name, LDisk(solid, Float64[0.0, 0.0, z_pos], orientation), vac)

    (merged, components)
end


"""
    pmt_top_flux(N, sg, materials, cfg, rng; kwargs...) -> NamedTuple

Compute top PMT flux tables. The three sub-components (PMTs, bases,
structure) are lumped into a single transparent volume spanning
z=[min, max] of all components. Gammas going downward pass the VE;
upward gammas are lost (filtered by cos_theta_to_lxe with :up orientation).

Returns flux tables and summed rate table for each isotope.
"""
function pmt_top_flux(N::Int, sg::Dict{String,SourceVolumeInfo},
                       materials::Dict{String,Material},
                       cfg::SimConfig, rng::AbstractRNG;
                       verbose::Bool=false, kwargs...)
    pmt_names = ["PMT_TOP_PMTs", "PMT_TOP_bases", "PMT_TOP_structure"]
    merged, components = _merge_pmt_volume(sg, pmt_names, "PMT_TOP_merged", :up, materials)

    t0 = time()
    bi = generate_flux_bi214(N, merged, cfg, rng; kwargs...)
    tl = generate_flux_tl208(N, merged, cfg, rng; kwargs...)
    verbose && @printf("  [pmt_top] flux generation done (%.1fs)\n", time() - t0)

    A_bi = "Bi214_mBq_per_kg"
    A_tl = "Tl208_mBq_per_kg"

    rate_bi = _build_rate_table(:top, BR_BI214_2448,
        [(sv.name, _get_activity_Bq(sv, A_bi), sv.mass_kg, bi) for sv in components],
        bi.E_min, bi.E_max, bi.n_E, bi.n_u)

    rate_tl = _build_rate_table(:top, BR_TL208_2615,
        [(sv.name, _get_activity_Bq(sv, A_tl), sv.mass_kg, tl) for sv in components],
        tl.E_min_main, tl.E_max_main, tl.n_E_main, tl.n_u)

    (bi214=bi, tl208=tl, bi214_rate=rate_bi, tl208_rate=rate_tl)
end


"""
    pmt_bottom_flux(N, sg, materials, cfg, rng; kwargs...) -> NamedTuple

Compute bottom PMT flux tables. Four sub-components (PMTs, bases,
structure, R8778_dome) lumped into one transparent volume. Gammas
going upward pass the VE; downward gammas are lost.
"""
function pmt_bottom_flux(N::Int, sg::Dict{String,SourceVolumeInfo},
                          materials::Dict{String,Material},
                          cfg::SimConfig, rng::AbstractRNG;
                          verbose::Bool=false, kwargs...)
    pmt_names = ["PMT_BOT_PMTs", "PMT_BOT_bases", "PMT_BOT_structure", "PMT_BOT_R8778_dome"]
    merged, components = _merge_pmt_volume(sg, pmt_names, "PMT_BOT_merged", :down, materials)

    t0 = time()
    bi = generate_flux_bi214(N, merged, cfg, rng; kwargs...)
    tl = generate_flux_tl208(N, merged, cfg, rng; kwargs...)
    verbose && @printf("  [pmt_bot] flux generation done (%.1fs)\n", time() - t0)

    A_bi = "Bi214_mBq_per_kg"
    A_tl = "Tl208_mBq_per_kg"

    rate_bi = _build_rate_table(:bottom, BR_BI214_2448,
        [(sv.name, _get_activity_Bq(sv, A_bi), sv.mass_kg, bi) for sv in components],
        bi.E_min, bi.E_max, bi.n_E, bi.n_u)

    rate_tl = _build_rate_table(:bottom, BR_TL208_2615,
        [(sv.name, _get_activity_Bq(sv, A_tl), sv.mass_kg, tl) for sv in components],
        tl.E_min_main, tl.E_max_main, tl.n_E_main, tl.n_u)

    (bi214=bi, tl208=tl, bi214_rate=rate_bi, tl208_rate=rate_tl)
end
