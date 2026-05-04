"""
Materials module for the LXe Monte Carlo.

Defines `Material` — a complete physics description of a detector material,
including density, composition, and (for active materials) photon cross
sections, electron stopping powers, CSDA range, and pre-tabulated
bremsstrahlung cross sections.

Materials are loaded from `data/materials.json` and `data/elements.json`.
XCOM and ESTAR tables are loaded from CSV files referenced in the JSON.

# Material capability levels

| Level          | Flags                          | Data loaded                    |
|:---------------|:-------------------------------|:-------------------------------|
| Passive        | `active=false`                 | XCOM only (photon attenuation) |
| Active         | `active=true, tracking=false`  | XCOM + ESTAR + brems table     |
| Tracking       | `active=true, tracking=true`   | same (straight-line e⁻ steps)  |
| Full tracking  | `full_tracking=true` (future)  | + MS, straggling for gas       |
"""


# =====================================================================
# Element lookup
# =====================================================================

"""
    ElementData

Basic element properties: atomic number and standard atomic weight.
"""
struct ElementData
    symbol::String
    Z::Int
    A::Float64
end


"""
    load_elements(path::AbstractString) -> Dict{String, ElementData}

Load element data from `elements.json`.
"""
function load_elements(path::AbstractString)::Dict{String,ElementData}
    raw = open(path, "r") do io
        JSON.parse(io)
    end
    elements = Dict{String,ElementData}()
    for (sym, d) in raw
        startswith(sym, "_") && continue
        elements[sym] = ElementData(sym, Int(d["Z"]), Float64(d["A"]))
    end
    elements
end


# =====================================================================
# Material struct
# =====================================================================

"""
    Material

Complete physics description of a detector material.

# Fields
- `name::String`: material identifier (e.g., "LXe", "PTFE")
- `density::Float64`: bulk density [g/cm³] (0 for vacuum)
- `components::Dict{String,Float64}`: element symbol → mass fraction
- `active::Bool`: has full transport data (XCOM + ESTAR + brems)
- `tracking::Bool`: electron transport with condensed-history stepping
- `full_tracking::Bool`: (future) multiple scattering + straggling
- `EK::Float64`: K-shell binding energy [MeV] (active materials, 0 otherwise)
- `Z_eff::Int`: effective atomic number (for single-element materials)
- `A_eff::Float64`: effective atomic weight [g/mol]
- `n_atom::Float64`: atom number density [cm⁻³] (derived)

# NIST data (loaded from CSV)
- `xcom`: photon cross-section table (all non-vacuum materials)
- `E_nudged`: K-edge-safe XCOM energy grid
- `estar`: electron stopping-power table (active materials only)
- `brems_T`, `brems_σ`: pre-tabulated brems cross section (active materials only)
"""
struct Material
    name::String
    density::Float64
    components::Dict{String,Float64}
    active::Bool
    tracking::Bool
    full_tracking::Bool
    EK::Float64
    Z_eff::Int
    A_eff::Float64
    n_atom::Float64

    # XCOM photon data
    xcom::Union{XCOMData,Nothing}
    E_nudged::Union{Vector{Float64},Nothing}

    # Pre-logged XCOM grids for fast interpolation
    log_E::Union{Vector{Float64},Nothing}          # log of E_nudged
    log_incoherent::Union{Vector{Float64},Nothing}
    log_photoelectric::Union{Vector{Float64},Nothing}
    log_pair_nuc::Union{Vector{Float64},Nothing}
    log_pair_ele::Union{Vector{Float64},Nothing}
    A_over_NA::Float64   # A_eff / N_A, pre-computed for per_atom conversion

    # ESTAR electron data (active only)
    estar::Union{ESTARData,Nothing}

    # Pre-logged ESTAR grids
    log_T::Union{Vector{Float64},Nothing}
    log_S_col::Union{Vector{Float64},Nothing}

    # Pre-tabulated brems cross section (active only)
    brems_T::Union{Vector{Float64},Nothing}
    brems_σ::Union{Vector{Float64},Nothing}
    brems_M::Union{Vector{Float64},Nothing}
    log_brems_T::Union{Vector{Float64},Nothing}
    log_brems_σ::Union{Vector{Float64},Nothing}
    log_brems_M::Union{Vector{Float64},Nothing}
end


# =====================================================================
# Material loader
# =====================================================================

"""
    load_materials(cfg::SimConfig; data_dir=nothing) -> Dict{String, Material}

Load all materials from `data/materials.json` and `data/elements.json`.
For each material, load XCOM (and ESTAR if active) from the referenced
CSV files. Pre-tabulates bremsstrahlung cross sections for active materials.

Missing CSV files produce a warning and a material with `nothing` data.
"""
function load_materials(cfg::SimConfig; data_dir::Union{String,Nothing}=nothing)::Dict{String,Material}
    if data_dir === nothing
        data_dir = normpath(joinpath(@__DIR__, "..", "..", "data"))
    end

    elements = load_elements(joinpath(data_dir, "elements.json"))

    raw = open(joinpath(data_dir, "materials.json"), "r") do io
        JSON.parse(io)
    end

    materials = Dict{String,Material}()
    for (name, d) in raw
        startswith(name, "_") && continue
        materials[name] = _build_material(name, d, elements, cfg, data_dir)
    end
    materials
end


function _build_material(name::String, d, elements::Dict{String,ElementData},
                         cfg::SimConfig, data_dir::String)::Material
    density = Float64(d["density"])
    comp_raw = d["components"]
    components = Dict{String,Float64}(String(k) => Float64(v) for (k, v) in comp_raw)

    is_active = Bool(get(d, "active", false))
    is_tracking = Bool(get(d, "tracking", false))
    is_full = Bool(get(d, "full_tracking", false))
    EK = Float64(get(d, "EK_MeV", 0.0))

    # Effective Z and A (mass-fraction weighted for compounds)
    Z_eff, A_eff = _effective_ZA(components, elements)

    # Atom number density
    n_atom = density > 0.0 ? cfg.N_A * density / A_eff : 0.0

    # Load XCOM
    xcom = nothing
    E_nudged = nothing
    xcom_file = get(d, "xcom", nothing)
    if xcom_file !== nothing
        xcom_path = joinpath(data_dir, xcom_file)
        if isfile(xcom_path)
            xcom = load_xcom(xcom_path)
            E_nudged = _prepare_xcom_energy(xcom)
        else
            @warn "XCOM file not found: $xcom_path (material $name)"
        end
    end

    # Load ESTAR (active materials only)
    estar = nothing
    if is_active
        estar_file = get(d, "estar", nothing)
        if estar_file !== nothing
            estar_path = joinpath(data_dir, estar_file)
            if isfile(estar_path)
                estar = load_estar(estar_path)
            else
                @warn "ESTAR file not found: $estar_path (material $name)"
            end
        end
    end

    # Pre-log XCOM grids
    log_E = E_nudged !== nothing ? log.(E_nudged) : nothing
    log_incoherent = xcom !== nothing ? prelog_data(xcom.incoherent) : nothing
    log_photoelectric = xcom !== nothing ? prelog_data(xcom.photoelectric) : nothing
    log_pair_nuc = xcom !== nothing ? prelog_data(xcom.pair_nuclear) : nothing
    log_pair_ele = xcom !== nothing ? prelog_data(xcom.pair_electron) : nothing
    A_over_NA = A_eff > 0.0 ? A_eff / cfg.N_A : 0.0

    # Pre-log ESTAR grids
    log_T_estar = (estar !== nothing) ? log.(estar.T_MeV) : nothing
    log_S_col = (estar !== nothing) ? prelog_data(estar.S_col) : nothing

    # Pre-tabulate brems (active materials only)
    brems_T = nothing
    brems_σ = nothing
    brems_M = nothing
    log_brems_T = nothing
    log_brems_σ = nothing
    log_brems_M = nothing
    if is_active && EK > 0.0
        brems_T, brems_σ, brems_M = _build_brems_table_for_material(cfg, Z_eff, A_eff)
        log_brems_T = log.(brems_T)
        log_brems_σ = prelog_data(brems_σ)
        log_brems_M = prelog_data(brems_M)
    end

    Material(name, density, components,
             is_active, is_tracking, is_full,
             EK, Z_eff, A_eff, n_atom,
             xcom, E_nudged,
             log_E, log_incoherent, log_photoelectric, log_pair_nuc, log_pair_ele,
             A_over_NA,
             estar, log_T_estar, log_S_col,
             brems_T, brems_σ, brems_M, log_brems_T, log_brems_σ, log_brems_M)
end


function _effective_ZA(components::Dict{String,Float64},
                       elements::Dict{String,ElementData})::Tuple{Int,Float64}
    isempty(components) && return (0, 0.0)

    # For single-element materials, return exact Z and A
    if length(components) == 1
        sym = first(keys(components))
        el = elements[sym]
        return (el.Z, el.A)
    end

    # For compounds: mass-fraction weighted average
    Z_sum = 0.0
    A_inv_sum = 0.0
    for (sym, frac) in components
        el = elements[sym]
        Z_sum += frac * el.Z
        A_inv_sum += frac / el.A
    end
    (round(Int, Z_sum), 1.0 / A_inv_sum)
end


"""
    _build_brems_table_for_material(cfg, Z, A) -> (T_grid, σ_grid, M_grid)

Pre-compute bremsstrahlung cross section σ(T) and rejection ceiling
M(T) = max{k × dσ/dk} × 1.05 on a log grid. M(T) is used by
`sample_brems` to avoid recomputing the 50-point grid at every call.
"""
function _build_brems_table_for_material(cfg::SimConfig, Z::Int, A::Float64)
    k_min = cfg.k_min
    T_max = 10.0
    n_grid = 200
    T_grid = exp.(range(log(k_min * 1.01), log(T_max), length=n_grid))
    σ_grid = Vector{Float64}(undef, n_grid)
    M_grid = Vector{Float64}(undef, n_grid)

    for i in 1:n_grid
        T = T_grid[i]
        σ_grid[i] = sigma_brems_above_kmin(T, k_min, cfg)

        # Rejection ceiling: max(k * dσ/dk) over k ∈ [k_min, T]
        k_pts = exp.(range(log(k_min), log(T * 0.9999), length=50))
        ds = dsigma_dk_brems(k_pts, T, Z, cfg)
        M_grid[i] = maximum(k_pts .* ds) * 1.05
    end

    (T_grid, σ_grid, M_grid)
end


# =====================================================================
# Physics methods on Material — photon cross sections
# =====================================================================

function _check_xcom(mat::Material)
    mat.xcom === nothing && error("Material '$(mat.name)' has no XCOM data")
end

function _check_active(mat::Material)
    mat.active || error("Material '$(mat.name)' is not active (no ESTAR/brems data)")
end


"""
    sigma_three(mat::Material, E_MeV) -> (σ_C, σ_P, σ_Ph)

Combined Compton, pair, and photoelectric cross sections [cm²/atom]
with a single binary search over the XCOM energy grid. This is the
hot path in photon transport — doing one search instead of three
saves ~15% of total CPU.
"""
function sigma_three(mat::Material, E_MeV::Float64)::Tuple{Float64,Float64,Float64}
    _check_xcom(mat)
    lx = log(E_MeV)
    n = length(mat.E_nudged)
    lo = searchsortedlast(mat.E_nudged, E_MeV)
    lo = clamp(lo, 1, n - 1)

    conv = mat.A_over_NA  # cm²/g → cm²/atom

    sC = interp_loglog_prelogged(lx, mat.log_E, mat.log_incoherent,
                                  mat.xcom.incoherent, lo) * conv
    sPh = interp_loglog_prelogged(lx, mat.log_E, mat.log_photoelectric,
                                   mat.xcom.photoelectric, lo) * conv

    # Pair = nuclear + electron
    sP_nuc = interp_loglog_prelogged(lx, mat.log_E, mat.log_pair_nuc,
                                      mat.xcom.pair_nuclear, lo) * conv
    sP_ele = interp_loglog_prelogged(lx, mat.log_E, mat.log_pair_ele,
                                      mat.xcom.pair_electron, lo) * conv

    (sC, sP_nuc + sP_ele, sPh)
end


"""
    sigma_compton(mat::Material, E_MeV; per_atom=true) -> Float64

Compton (incoherent) cross section from NIST XCOM.
"""
function sigma_compton(mat::Material, E_MeV::Float64; per_atom::Bool=true)::Float64
    _check_xcom(mat)
    _interp_xcom(mat.E_nudged, mat.xcom.incoherent, E_MeV;
                 per_atom=per_atom, A=mat.A_eff, N_A=6.02214076e23)
end


"""
    sigma_pair(mat::Material, E_MeV; per_atom=true) -> Float64

Pair production cross section from NIST XCOM.
"""
function sigma_pair(mat::Material, E_MeV::Float64; per_atom::Bool=true)::Float64
    _check_xcom(mat)
    nuc = _interp_xcom(mat.E_nudged, mat.xcom.pair_nuclear, E_MeV;
                       per_atom=per_atom, A=mat.A_eff, N_A=6.02214076e23)
    ele = _interp_xcom(mat.E_nudged, mat.xcom.pair_electron, E_MeV;
                       per_atom=per_atom, A=mat.A_eff, N_A=6.02214076e23)
    nuc + ele
end


"""
    sigma_phot(mat::Material, E_MeV; per_atom=true) -> Float64

Photoelectric cross section from NIST XCOM.
"""
function sigma_phot(mat::Material, E_MeV::Float64; per_atom::Bool=true)::Float64
    _check_xcom(mat)
    _interp_xcom(mat.E_nudged, mat.xcom.photoelectric, E_MeV;
                 per_atom=per_atom, A=mat.A_eff, N_A=6.02214076e23)
end


"""
    sigma_total(mat::Material, E_MeV; per_atom=true, include_coherent=false) -> Float64

Total photon cross section from NIST XCOM.
"""
function sigma_total(mat::Material, E_MeV::Float64;
                     per_atom::Bool=true, include_coherent::Bool=false)::Float64
    _check_xcom(mat)
    ch = include_coherent ? mat.xcom.total_w_coh : mat.xcom.total_no_coh
    _interp_xcom(mat.E_nudged, ch, E_MeV;
                 per_atom=per_atom, A=mat.A_eff, N_A=6.02214076e23)
end


"""
    mu_over_rho(mat::Material, E_MeV) -> Float64

Mass attenuation coefficient µ/ρ [cm²/g] from NIST XCOM.
"""
function mu_over_rho(mat::Material, E_MeV::Float64)::Float64
    sigma_total(mat, E_MeV; per_atom=false, include_coherent=false)
end


"""
    mfp(mat::Material, E_MeV) -> Float64

Photon mean free path [cm] in this material.
"""
function mfp(mat::Material, E_MeV::Float64)::Float64
    mat.density <= 0.0 && return Inf
    1.0 / (mu_over_rho(mat, E_MeV) * mat.density)
end


"""
    branching(mat::Material, E_MeV) -> (fC, fP, fPh)

Compton, pair, and photoelectric branching fractions.
"""
function branching(mat::Material, E_MeV::Float64)::Tuple{Float64,Float64,Float64}
    sC  = sigma_compton(mat, E_MeV; per_atom=false)
    sP  = sigma_pair(mat, E_MeV; per_atom=false)
    sPh = sigma_phot(mat, E_MeV; per_atom=false)
    st  = sC + sP + sPh
    (sC/st, sP/st, sPh/st)
end


# =====================================================================
# Physics methods on Material — electron stopping power and range
# =====================================================================

"""
    dEdx_collision(mat::Material, T_MeV) -> Float64

Mass collision stopping power [MeV cm²/g] from NIST ESTAR.
Only available for active materials.
"""
function dEdx_collision(mat::Material, T_MeV::Float64)::Float64
    _check_active(mat)
    n = length(mat.estar.T_MeV)
    lo = searchsortedlast(mat.estar.T_MeV, T_MeV)
    lo = clamp(lo, 1, n - 1)
    interp_loglog_prelogged(log(T_MeV), mat.log_T, mat.log_S_col,
                            mat.estar.S_col, lo)
end


"""
    dEdx_radiative(mat::Material, T_MeV) -> Float64

Mass radiative stopping power [MeV cm²/g] from NIST ESTAR.
"""
function dEdx_radiative(mat::Material, T_MeV::Float64)::Float64
    _check_active(mat)
    interp_loglog(T_MeV, mat.estar.T_MeV, mat.estar.S_rad)
end


"""
    dEdx_total(mat::Material, T_MeV) -> Float64

Total mass stopping power [MeV cm²/g] from NIST ESTAR.
"""
function dEdx_total(mat::Material, T_MeV::Float64)::Float64
    _check_active(mat)
    interp_loglog(T_MeV, mat.estar.T_MeV, mat.estar.S_tot)
end


"""
    csda_range_g_per_cm2(mat::Material, T_MeV) -> Float64

CSDA range [g/cm²] from NIST ESTAR.
"""
function csda_range_g_per_cm2(mat::Material, T_MeV::Float64)::Float64
    _check_active(mat)
    interp_loglog(T_MeV, mat.estar.T_MeV, mat.estar.R_csda)
end


"""
    csda_range_mm(mat::Material, T_MeV) -> Float64

CSDA range [mm] in this material.
"""
function csda_range_mm(mat::Material, T_MeV::Float64)::Float64
    mat.density <= 0.0 && return Inf
    csda_range_g_per_cm2(mat, T_MeV) / mat.density * 10.0
end


"""
    sigma_brems(mat::Material, T_MeV) -> Float64

Pre-tabulated total bremsstrahlung cross section [cm²/atom] for k > k_min.
"""
function sigma_brems(mat::Material, T_MeV::Float64)::Float64
    _check_active(mat)
    mat.brems_T === nothing && return 0.0
    n = length(mat.brems_T)
    lo = searchsortedlast(mat.brems_T, T_MeV)
    lo = clamp(lo, 1, n - 1)
    interp_loglog_prelogged(log(T_MeV), mat.log_brems_T, mat.log_brems_σ,
                            mat.brems_σ, lo)
end


"""
    brems_rejection_M(mat::Material, T_MeV) -> Float64

Pre-tabulated rejection ceiling M(T) for bremsstrahlung sampling.
"""
function brems_rejection_M(mat::Material, T_MeV::Float64)::Float64
    _check_active(mat)
    mat.brems_M === nothing && return 0.0
    n = length(mat.brems_T)
    lo = searchsortedlast(mat.brems_T, T_MeV)
    lo = clamp(lo, 1, n - 1)
    interp_loglog_prelogged(log(T_MeV), mat.log_brems_T, mat.log_brems_M,
                            mat.brems_M, lo)
end
