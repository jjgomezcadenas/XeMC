"""
Cryostat flux functions: barrel, top, and bottom.

Each function computes the flux of gammas arriving at the ICV inner surface
from all contributing source layers (OCV, MLI, ICV), for both Bi-214 and
Tl-208. Returns component flux tables plus activity-weighted summed rate
tables.
"""


"""
    _build_rate_table(surface, gamma_BR, components, E_min, E_max, n_E, n_u)

Build a `SourceRateTable` by summing activity-weighted component flux tables.
Each component is `(name, A_Bq_per_kg, mass_kg, flux)`.

`gamma_BR` is the branching ratio from chain decay to the specific gamma
line simulated (e.g. `BR_BI214_2448` or `BR_TL208_2615`).  The chain
activity `A_Bq_per_kg` is multiplied by `gamma_BR` to obtain the
single-gamma rate.
"""
function _build_rate_table(surface::Symbol,
                            gamma_BR::Float64,
                            components::Vector{<:Tuple},
                            E_min::Float64, E_max::Float64,
                            n_E::Int, n_u::Int)::SourceRateTable
    pdf_rate = zeros(n_E, n_u)
    names = String[]
    rates = Float64[]

    for (name, A_Bq_per_kg, mass_kg, flux) in components
        pdf = flux isa SourceFluxBi214 ? flux.pdf : flux.pdf_main
        rate_scale = A_Bq_per_kg * mass_kg * gamma_BR
        if size(pdf) == (n_E, n_u)
            pdf_rate .+= rate_scale .* pdf
        end
        push!(names, name)
        push!(rates, rate_scale * sum(pdf))
    end

    SourceRateTable(surface, pdf_rate, E_min, E_max, n_E, n_u,
                    names, rates, sum(rates))
end


"""Convert mBq/kg activity to Bq/kg. Returns 0.0 if key is missing."""
function _get_activity_Bq(sv::SourceVolumeInfo, key::String)::Float64
    get(sv.activity, key, 0.0) * 1e-3
end


"""
    cryostat_barrel_flux(N, sg, cfg, rng; kwargs...)

Compute barrel flux tables. Three contributing sources:
1. OCV_barrel -> vacuum -> ICV_barrel (compound)
2. MLI -> ICV_barrel (transparent: propagate through ICV only)
3. ICV_barrel self

Returns named tuple with component tables and summed rate tables.
"""
function cryostat_barrel_flux(N::Int, sg::Dict{String,SourceVolumeInfo},
                               cfg::SimConfig, rng::AbstractRNG;
                               verbose::Bool=false, kwargs...)
    ocv = sg["OCV_barrel"]
    mli = sg["MLI"]
    icv = sg["ICV_barrel"]
    layers_icv = PhysicalVolume[icv.volume]

    t0 = time()
    bi_ocv = generate_flux_compound_bi214(N, ocv.volume, layers_icv, icv.volume, cfg, rng; kwargs...)
    tl_ocv = generate_flux_compound_tl208(N, ocv.volume, layers_icv, icv.volume, cfg, rng; kwargs...)
    verbose && @printf("  [barrel] OCV compound done (%.1fs)\n", time() - t0)

    t1 = time()
    bi_mli = generate_flux_bi214(N, icv.volume, cfg, rng; kwargs...)
    tl_mli = generate_flux_tl208(N, icv.volume, cfg, rng; kwargs...)
    verbose && @printf("  [barrel] MLI done (%.1fs)\n", time() - t1)

    t2 = time()
    bi_icv = generate_flux_bi214(N, icv.volume, cfg, rng; kwargs...)
    tl_icv = generate_flux_tl208(N, icv.volume, cfg, rng; kwargs...)
    verbose && @printf("  [barrel] ICV done (%.1fs)\n", time() - t2)

    A_bi = "Bi214_mBq_per_kg"
    A_tl = "Tl208_mBq_per_kg"

    rate_bi = _build_rate_table(:barrel, BR_BI214_2448,
        [("OCV_barrel", _get_activity_Bq(ocv, A_bi), ocv.mass_kg, bi_ocv),
         ("MLI",        _get_activity_Bq(mli, A_bi), mli.mass_kg, bi_mli),
         ("ICV_barrel", _get_activity_Bq(icv, A_bi), icv.mass_kg, bi_icv)],
        bi_ocv.E_min, bi_ocv.E_max, bi_ocv.n_E, bi_ocv.n_u)

    rate_tl = _build_rate_table(:barrel, BR_TL208_2615,
        [("OCV_barrel", _get_activity_Bq(ocv, A_tl), ocv.mass_kg, tl_ocv),
         ("MLI",        _get_activity_Bq(mli, A_tl), mli.mass_kg, tl_mli),
         ("ICV_barrel", _get_activity_Bq(icv, A_tl), icv.mass_kg, tl_icv)],
        tl_ocv.E_min_main, tl_ocv.E_max_main, tl_ocv.n_E_main, tl_ocv.n_u)

    (bi214_ocv=bi_ocv, bi214_mli=bi_mli, bi214_icv=bi_icv, bi214_rate=rate_bi,
     tl208_ocv=tl_ocv, tl208_mli=tl_mli, tl208_icv=tl_icv, tl208_rate=rate_tl)
end


"""
    cryostat_top_flux(N, sg, cfg, rng; kwargs...)

Compute top head flux tables. Two contributing sources:
1. OCV_top -> vacuum -> ICV_top (compound)
2. ICV_top self
"""
function cryostat_top_flux(N::Int, sg::Dict{String,SourceVolumeInfo},
                            cfg::SimConfig, rng::AbstractRNG;
                            verbose::Bool=false, kwargs...)
    ocv = sg["OCV_top"]
    icv = sg["ICV_top"]
    layers_icv = PhysicalVolume[icv.volume]

    t0 = time()
    bi_ocv = generate_flux_compound_bi214(N, ocv.volume, layers_icv, icv.volume, cfg, rng; kwargs...)
    tl_ocv = generate_flux_compound_tl208(N, ocv.volume, layers_icv, icv.volume, cfg, rng; kwargs...)
    verbose && @printf("  [top] OCV compound done (%.1fs)\n", time() - t0)

    t1 = time()
    bi_icv = generate_flux_bi214(N, icv.volume, cfg, rng; kwargs...)
    tl_icv = generate_flux_tl208(N, icv.volume, cfg, rng; kwargs...)
    verbose && @printf("  [top] ICV done (%.1fs)\n", time() - t1)

    A_bi = "Bi214_mBq_per_kg"
    A_tl = "Tl208_mBq_per_kg"

    rate_bi = _build_rate_table(:top, BR_BI214_2448,
        [("OCV_top", _get_activity_Bq(ocv, A_bi), ocv.mass_kg, bi_ocv),
         ("ICV_top", _get_activity_Bq(icv, A_bi), icv.mass_kg, bi_icv)],
        bi_ocv.E_min, bi_ocv.E_max, bi_ocv.n_E, bi_ocv.n_u)

    rate_tl = _build_rate_table(:top, BR_TL208_2615,
        [("OCV_top", _get_activity_Bq(ocv, A_tl), ocv.mass_kg, tl_ocv),
         ("ICV_top", _get_activity_Bq(icv, A_tl), icv.mass_kg, tl_icv)],
        tl_ocv.E_min_main, tl_ocv.E_max_main, tl_ocv.n_E_main, tl_ocv.n_u)

    (bi214_ocv=bi_ocv, bi214_icv=bi_icv, bi214_rate=rate_bi,
     tl208_ocv=tl_ocv, tl208_icv=tl_icv, tl208_rate=rate_tl)
end


"""
    cryostat_bottom_flux(N, sg, cfg, rng; kwargs...)

Compute bottom head flux tables. Two contributing sources:
1. OCV_bottom -> vacuum -> ICV_bottom (compound)
2. ICV_bottom self
"""
function cryostat_bottom_flux(N::Int, sg::Dict{String,SourceVolumeInfo},
                               cfg::SimConfig, rng::AbstractRNG;
                               verbose::Bool=false, kwargs...)
    ocv = sg["OCV_bottom"]
    icv = sg["ICV_bottom"]
    layers_icv = PhysicalVolume[icv.volume]

    t0 = time()
    bi_ocv = generate_flux_compound_bi214(N, ocv.volume, layers_icv, icv.volume, cfg, rng; kwargs...)
    tl_ocv = generate_flux_compound_tl208(N, ocv.volume, layers_icv, icv.volume, cfg, rng; kwargs...)
    verbose && @printf("  [bottom] OCV compound done (%.1fs)\n", time() - t0)

    t1 = time()
    bi_icv = generate_flux_bi214(N, icv.volume, cfg, rng; kwargs...)
    tl_icv = generate_flux_tl208(N, icv.volume, cfg, rng; kwargs...)
    verbose && @printf("  [bottom] ICV done (%.1fs)\n", time() - t1)

    A_bi = "Bi214_mBq_per_kg"
    A_tl = "Tl208_mBq_per_kg"

    rate_bi = _build_rate_table(:bottom, BR_BI214_2448,
        [("OCV_bottom", _get_activity_Bq(ocv, A_bi), ocv.mass_kg, bi_ocv),
         ("ICV_bottom", _get_activity_Bq(icv, A_bi), icv.mass_kg, bi_icv)],
        bi_ocv.E_min, bi_ocv.E_max, bi_ocv.n_E, bi_ocv.n_u)

    rate_tl = _build_rate_table(:bottom, BR_TL208_2615,
        [("OCV_bottom", _get_activity_Bq(ocv, A_tl), ocv.mass_kg, tl_ocv),
         ("ICV_bottom", _get_activity_Bq(icv, A_tl), icv.mass_kg, tl_icv)],
        tl_ocv.E_min_main, tl_ocv.E_max_main, tl_ocv.n_E_main, tl_ocv.n_u)

    (bi214_ocv=bi_ocv, bi214_icv=bi_icv, bi214_rate=rate_bi,
     tl208_ocv=tl_ocv, tl208_icv=tl_icv, tl208_rate=rate_tl)
end
