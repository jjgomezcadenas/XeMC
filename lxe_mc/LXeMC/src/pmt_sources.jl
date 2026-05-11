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
    _merge_pmt_volume(sg, names, merged_name, orientation)
        -> (PDisk, Vector{SourceVolumeInfo})

Build a merged flat disk volume spanning the full z extent and common
radius of the named PMT sub-components. The merged volume is vacuum
(transparent, material taken from first component) with the given
orientation for hemisphere filtering.

The z span is [min(z - hh), max(z + hh)] over all components — no
hardwired numbers, always computed from the source geometry.
"""
function _merge_pmt_volume(sg::Dict{String,SourceVolumeInfo},
                            names::Vector{String},
                            merged_name::String,
                            orientation::Symbol)::Tuple{PDisk,Vector{SourceVolumeInfo}}
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
    # All PMT components are vacuum — take material from first component
    vac = components[1].material
    merged = PDisk(merged_name, LDisk(solid, Float64[0.0, 0.0, z_pos], orientation), vac)

    (merged, components)
end


"""
    pmt_top_flux(N, sg, cfg, rng; kwargs...) -> NamedTuple

Compute top PMT flux tables. The three sub-components (PMTs, bases,
structure) are lumped into a single transparent volume spanning
z=[min, max] of all components. Gammas going downward pass the VE;
upward gammas are lost (filtered by cos_theta_to_lxe with :up orientation).

Returns flux tables and summed rate table for each isotope.
"""
function pmt_top_flux(N::Int, sg::Dict{String,SourceVolumeInfo},
                       cfg::SimConfig, rng::AbstractRNG;
                       verbose::Bool=false, kwargs...)
    pmt_names = ["PMT_TOP_PMTs", "PMT_TOP_bases", "PMT_TOP_structure"]
    merged, components = _merge_pmt_volume(sg, pmt_names, "PMT_TOP_merged", :up)

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
    pmt_bottom_flux(N, sg, cfg, rng; kwargs...) -> NamedTuple

Compute bottom PMT flux tables. Four sub-components (PMTs, bases,
structure, R8778_dome) lumped into one transparent volume. Gammas
going upward pass the VE; downward gammas are lost.
"""
function pmt_bottom_flux(N::Int, sg::Dict{String,SourceVolumeInfo},
                          cfg::SimConfig, rng::AbstractRNG;
                          verbose::Bool=false, kwargs...)
    pmt_names = ["PMT_BOT_PMTs", "PMT_BOT_bases", "PMT_BOT_structure", "PMT_BOT_R8778_dome"]
    merged, components = _merge_pmt_volume(sg, pmt_names, "PMT_BOT_merged", :down)

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


"""
    _cathode_z(sg) -> Float64

Return the cathode z position from the source geometry.
Derived as the top face of FC_botgrid (position_z + half_height).
"""
function _cathode_z(sg::Dict{String,SourceVolumeInfo})::Float64
    bg = sg["FC_botgrid"]
    vol = bg.volume
    if vol isa PCylShell
        vol.logical.position[3] + vol.logical.solid.half_height_cm
    elseif vol isa PCyl
        vol.logical.position[3] + vol.logical.solid.half_height_cm
    else
        error("FC_botgrid has unexpected volume type: $(typeof(vol))")
    end
end


"""
    pmt_bottom_lxe_flux(N, sg, mats, cfg, rng; kwargs...) -> NamedTuple

Compute bottom PMT flux tables at the cathode after propagation
through the passive LXe below the cathode.

Gammas are born at the PMT surface (transparent, hemisphere-filtered
upward), then propagated through a slab of passive LXe from
z_pmt_top to z_cathode. The flux is recorded at the cathode.

This allows comparing bottom PMT backgrounds with top PMTs on equal
footing: from the cathode onward, both see only active LXe before
the FV.
"""
function pmt_bottom_lxe_flux(N::Int, sg::Dict{String,SourceVolumeInfo},
                              mats::Dict{String,Material},
                              cfg::SimConfig, rng::AbstractRNG;
                              verbose::Bool=false, kwargs...)
    pmt_names = ["PMT_BOT_PMTs", "PMT_BOT_bases", "PMT_BOT_structure", "PMT_BOT_R8778_dome"]
    merged, components = _merge_pmt_volume(sg, pmt_names, "PMT_BOT_merged", :down)

    # Cathode z from source geometry (top face of FC_botgrid)
    z_cathode = _cathode_z(sg)
    z_pmt = merged.logical.position[3]  # top face of merged disk
    R_pmt = merged.logical.solid.radius_cm
    hh = (z_cathode - z_pmt) / 2.0
    z_center = (z_pmt + z_cathode) / 2.0

    verbose && @printf("  [pmt_bot_lxe] LXe slab: z=[%.2f, %.2f], R=%.1f\n",
                       z_pmt, z_cathode, R_pmt)

    lxe_mat = mats["LXe"]
    lxe_slab = PCyl("LXe_passive_slab",
                     LCyl(Cyl(R_pmt, hh), Float64[0.0, 0.0, z_center]),
                     lxe_mat)

    # Exit surface at cathode: orientation :down means
    # cos_theta_to_lxe = +dir_z (upward = toward active LXe)
    exit_disk = PDisk("cathode_surface",
                      LDisk(Disk(R_pmt, 0.01, Inf), Float64[0.0, 0.0, z_cathode], :down),
                      lxe_mat)

    layers = PhysicalVolume[lxe_slab]

    t0 = time()
    bi = generate_flux_compound_bi214(N, merged, layers, exit_disk, cfg, rng; kwargs...)
    tl = generate_flux_compound_tl208(N, merged, layers, exit_disk, cfg, rng; kwargs...)
    verbose && @printf("  [pmt_bot_lxe] flux generation done (%.1fs)\n", time() - t0)

    A_bi = "Bi214_mBq_per_kg"
    A_tl = "Tl208_mBq_per_kg"

    rate_bi_lxe = _build_rate_table(:bottom, BR_BI214_2448,
        [(sv.name, _get_activity_Bq(sv, A_bi), sv.mass_kg, bi) for sv in components],
        bi.E_min, bi.E_max, bi.n_E, bi.n_u)

    rate_tl_lxe = _build_rate_table(:bottom, BR_TL208_2615,
        [(sv.name, _get_activity_Bq(sv, A_tl), sv.mass_kg, tl) for sv in components],
        tl.E_min_main, tl.E_max_main, tl.n_E_main, tl.n_u)

    (bi214=bi, tl208=tl, bi214_rate=rate_bi_lxe, tl208_rate=rate_tl_lxe)
end
