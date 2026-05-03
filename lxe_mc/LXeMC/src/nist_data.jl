"""
NIST XCOM and ESTAR data tables for xenon, with log-log interpolators.

These are the authoritative numerical values used in the Monte Carlo.
The XCOM table provides photon cross sections (cm²/g); the ESTAR table
provides electron stopping powers (MeV cm²/g).

# Data structures

- [`XCOMData`](@ref): holds XCOM energy grid and cross-section columns.
- [`ESTARData`](@ref): holds ESTAR energy grid, stopping powers, and
  CSDA range computed by numerical integration of 1/S_total.

# Interpolation

All interpolations use log-log linear interpolation via
[`interp_loglog`](@ref). This is exact for power-law regions and
accurate to < 1% between NIST grid points.
"""


# =====================================================================
# Data structures
# =====================================================================

"""
    XCOMData

NIST XCOM photon cross sections for Xe. All fields are `Vector{Float64}`
with units cm²/g (mass attenuation coefficients).
"""
struct XCOMData
    E_MeV::Vector{Float64}
    coherent::Vector{Float64}
    incoherent::Vector{Float64}
    photoelectric::Vector{Float64}
    pair_nuclear::Vector{Float64}
    pair_electron::Vector{Float64}
    total_w_coh::Vector{Float64}
    total_no_coh::Vector{Float64}
end


"""
    ESTARData

NIST ESTAR electron stopping powers for Xe and derived CSDA range.
Stopping powers in MeV cm²/g; CSDA range in g/cm².
"""
struct ESTARData
    T_MeV::Vector{Float64}
    S_col::Vector{Float64}
    S_rad::Vector{Float64}
    S_tot::Vector{Float64}
    delta::Vector{Float64}
    R_csda::Vector{Float64}   # g/cm², computed by integration of 1/S_tot
end


# =====================================================================
# CSV loaders
# =====================================================================

"""
    load_xcom(path::AbstractString) -> XCOMData

Read NIST XCOM CSV for xenon. Lines starting with `#` or the header
`E_MeV` are skipped. The K-edge duplicate at 0.034561 MeV is kept
as-is; the interpolator handles it by nudging.
"""
function load_xcom(path::AbstractString)::XCOMData
    rows = Vector{Vector{Float64}}()
    open(path, "r") do io
        for line in eachline(io)
            s = strip(line)
            (isempty(s) || startswith(s, "#") || startswith(s, "E_MeV")) && continue
            push!(rows, parse.(Float64, split(s, ",")))
        end
    end
    m = reduce(hcat, rows)'  # N×8 matrix
    XCOMData(m[:,1], m[:,2], m[:,3], m[:,4], m[:,5], m[:,6], m[:,7], m[:,8])
end


"""
    load_estar(path::AbstractString) -> ESTARData

Read NIST ESTAR CSV for xenon. Computes CSDA range R(T) [g/cm²] by
trapezoidal integration of 1/S_total over the tabulated grid.
"""
function load_estar(path::AbstractString)::ESTARData
    rows = Vector{Vector{Float64}}()
    open(path, "r") do io
        for line in eachline(io)
            s = strip(line)
            (isempty(s) || startswith(s, "#") || startswith(s, "T_MeV")) && continue
            push!(rows, parse.(Float64, split(s, ",")))
        end
    end
    m = reduce(hcat, rows)'
    T = m[:,1]; Sc = m[:,2]; Sr = m[:,3]; St = m[:,4]; delta = m[:,5]

    # CSDA range by trapezoidal integration of 1/S_tot
    R = zeros(length(T))
    @inbounds for i in 2:length(T)
        R[i] = R[i-1] + 0.5 * (1.0/St[i-1] + 1.0/St[i]) * (T[i] - T[i-1])
    end

    ESTARData(T, Sc, Sr, St, delta, R)
end


# =====================================================================
# Log-log interpolation
# =====================================================================

"""
    interp_loglog(x, xp, fp) -> Float64

Linear interpolation in log-log space. `xp` must be increasing and
positive. Zero or negative values in `fp` are handled by returning 0
when `x` falls in or below the zero region.

Used for cross sections and stopping powers, which are smooth power-law
functions between grid points.
"""
function interp_loglog(x::Float64, xp::Vector{Float64}, fp::Vector{Float64})::Float64
    x <= 0.0 && return 0.0

    # Find positive-value region
    n = length(xp)

    # Below table range
    if x < xp[1]
        return fp[1] > 0.0 ? fp[1] : 0.0
    end
    # Above table range
    if x > xp[end]
        return fp[end]
    end

    # Binary search for bracketing interval
    lo = searchsortedlast(xp, x)
    lo = clamp(lo, 1, n - 1)
    hi = lo + 1

    # Handle zero/negative values (e.g., pair below threshold)
    if fp[lo] <= 0.0 && fp[hi] <= 0.0
        return 0.0
    elseif fp[lo] <= 0.0
        return 0.0  # below first positive value
    elseif fp[hi] <= 0.0
        return fp[lo]
    end

    # Log-log linear interpolation
    lx = log(x)
    t = (lx - log(xp[lo])) / (log(xp[hi]) - log(xp[lo]))
    return exp(log(fp[lo]) + t * (log(fp[hi]) - log(fp[lo])))
end


# =====================================================================
# XCOM interpolation helpers
# =====================================================================

"""
    _prepare_xcom_energy(xc::XCOMData) -> Vector{Float64}

Return a copy of the XCOM energy grid with duplicate K-edge entries
nudged apart so that `searchsortedlast` works correctly.
"""
function _prepare_xcom_energy(xc::XCOMData)::Vector{Float64}
    E = copy(xc.E_MeV)
    @inbounds for i in 2:length(E)
        if E[i] <= E[i-1]
            E[i] = E[i-1] * (1.0 + 1e-9)
        end
    end
    E
end


"""
    _interp_xcom(xc::XCOMData, E_nudged::Vector{Float64},
                 channel::Vector{Float64}, E_MeV::Float64;
                 per_atom::Bool=true, cfg::SimConfig) -> Float64

Interpolate one XCOM channel at energy `E_MeV`.
Returns cm²/atom if `per_atom=true`, otherwise cm²/g.
"""
function _interp_xcom(E_nudged::Vector{Float64},
                      channel::Vector{Float64},
                      E_MeV::Float64;
                      per_atom::Bool=true,
                      A::Float64=131.293,
                      N_A::Float64=6.02214076e23)::Float64
    val = interp_loglog(E_MeV, E_nudged, channel)
    per_atom ? val * A / N_A : val
end


# =====================================================================
# Public photon cross-section API
# =====================================================================

"""
    NISTData

Bundle of loaded XCOM and ESTAR data with pre-computed energy grid
for efficient repeated interpolation. Includes a pre-tabulated
bremsstrahlung cross section σ_brems(T) for k > k_min, eliminating
the need for per-step numerical integration during transport.
"""
struct NISTData
    xcom::XCOMData
    estar::ESTARData
    E_nudged::Vector{Float64}   # K-edge-safe XCOM energy grid
    A::Float64
    N_A::Float64
    rho::Float64
    # Pre-tabulated brems cross section
    brems_T::Vector{Float64}    # log-spaced kinetic energy grid [MeV]
    brems_σ::Vector{Float64}    # σ_brems(T) for k > k_min [cm²/atom]
end


"""
    load_nist_data(cfg::SimConfig; data_dir=nothing) -> NISTData

Load all NIST tables from CSV files in `data_dir` (default: `../../data/`
relative to source). Pre-computes the K-edge-nudged energy grid.
"""
function load_nist_data(cfg::SimConfig; data_dir::Union{String,Nothing}=nothing)::NISTData
    if data_dir === nothing
        data_dir = normpath(joinpath(@__DIR__, "..", "..", "data"))
    end
    xcom = load_xcom(joinpath(data_dir, "xcom_xe.csv"))
    estar = load_estar(joinpath(data_dir, "estar_xe.csv"))
    E_nudged = _prepare_xcom_energy(xcom)

    # Pre-tabulate brems cross section σ(T) for k > k_min
    brems_T, brems_σ = _build_brems_table(cfg)

    NISTData(xcom, estar, E_nudged, cfg.A, cfg.N_A, cfg.rho_LXe, brems_T, brems_σ)
end


"""
    _build_brems_table(cfg::SimConfig) -> (T_grid, σ_grid)

Pre-compute the total bremsstrahlung cross section [cm²/atom] for photon
energies k > k_min on a 200-point log-spaced grid from k_min to 10 MeV.
Each point is computed by trapezoidal integration of the BH-Tsai
differential. This table is built once at load time and used via
log-log interpolation during transport.
"""
function _build_brems_table(cfg::SimConfig)
    T_max = 10.0  # MeV, well above our ROI
    k_min = cfg.k_min
    n_grid = 200
    T_grid = exp.(range(log(k_min * 1.01), log(T_max), length=n_grid))
    σ_grid = Vector{Float64}(undef, n_grid)

    for i in 1:n_grid
        σ_grid[i] = sigma_brems_above_kmin(T_grid[i], k_min, cfg)
    end

    (T_grid, σ_grid)
end


"""
    sigma_compton_NIST(nd::NISTData, E_MeV; per_atom=true) -> Float64

Compton (incoherent) scattering cross section from NIST XCOM.
Returns cm²/atom (default) or cm²/g.
"""
function sigma_compton_NIST(nd::NISTData, E_MeV::Float64; per_atom::Bool=true)::Float64
    _interp_xcom(nd.E_nudged, nd.xcom.incoherent, E_MeV;
                 per_atom=per_atom, A=nd.A, N_A=nd.N_A)
end


"""
    sigma_pair_NIST(nd::NISTData, E_MeV; per_atom=true) -> Float64

Pair production (nuclear + electron field) cross section from NIST XCOM.
"""
function sigma_pair_NIST(nd::NISTData, E_MeV::Float64; per_atom::Bool=true)::Float64
    nuc = _interp_xcom(nd.E_nudged, nd.xcom.pair_nuclear, E_MeV;
                       per_atom=per_atom, A=nd.A, N_A=nd.N_A)
    ele = _interp_xcom(nd.E_nudged, nd.xcom.pair_electron, E_MeV;
                       per_atom=per_atom, A=nd.A, N_A=nd.N_A)
    nuc + ele
end


"""
    sigma_phot_NIST(nd::NISTData, E_MeV; per_atom=true) -> Float64

Photoelectric cross section from NIST XCOM.
"""
function sigma_phot_NIST(nd::NISTData, E_MeV::Float64; per_atom::Bool=true)::Float64
    _interp_xcom(nd.E_nudged, nd.xcom.photoelectric, E_MeV;
                 per_atom=per_atom, A=nd.A, N_A=nd.N_A)
end


"""
    sigma_total_NIST(nd::NISTData, E_MeV; per_atom=true, include_coherent=false) -> Float64

Total photon cross section from NIST XCOM.
"""
function sigma_total_NIST(nd::NISTData, E_MeV::Float64;
                          per_atom::Bool=true,
                          include_coherent::Bool=false)::Float64
    channel = include_coherent ? nd.xcom.total_w_coh : nd.xcom.total_no_coh
    _interp_xcom(nd.E_nudged, channel, E_MeV;
                 per_atom=per_atom, A=nd.A, N_A=nd.N_A)
end


"""
    mfp_LXe(nd::NISTData, E_MeV) -> Float64

Photon mean free path [cm] in liquid xenon.
"""
function mfp_LXe(nd::NISTData, E_MeV::Float64)::Float64
    mu_rho = sigma_total_NIST(nd, E_MeV; per_atom=false)
    1.0 / (mu_rho * nd.rho)
end


"""
    branching_NIST(nd::NISTData, E_MeV) -> (fC, fP, fPh)

Compton, pair, and photoelectric branching fractions at energy `E_MeV`.
"""
function branching_NIST(nd::NISTData, E_MeV::Float64)::Tuple{Float64,Float64,Float64}
    sC  = sigma_compton_NIST(nd, E_MeV; per_atom=false)
    sP  = sigma_pair_NIST(nd, E_MeV; per_atom=false)
    sPh = sigma_phot_NIST(nd, E_MeV; per_atom=false)
    st  = sC + sP + sPh
    (sC/st, sP/st, sPh/st)
end


# =====================================================================
# Public electron stopping-power and range API
# =====================================================================

"""
    dEdx_collision_NIST(nd::NISTData, T_MeV) -> Float64

Mass collision stopping power [MeV cm²/g] from NIST ESTAR.
"""
function dEdx_collision_NIST(nd::NISTData, T_MeV::Float64)::Float64
    interp_loglog(T_MeV, nd.estar.T_MeV, nd.estar.S_col)
end


"""
    dEdx_radiative_NIST(nd::NISTData, T_MeV) -> Float64

Mass radiative stopping power [MeV cm²/g] from NIST ESTAR.
"""
function dEdx_radiative_NIST(nd::NISTData, T_MeV::Float64)::Float64
    interp_loglog(T_MeV, nd.estar.T_MeV, nd.estar.S_rad)
end


"""
    dEdx_total_NIST(nd::NISTData, T_MeV) -> Float64

Total mass stopping power [MeV cm²/g] from NIST ESTAR.
"""
function dEdx_total_NIST(nd::NISTData, T_MeV::Float64)::Float64
    interp_loglog(T_MeV, nd.estar.T_MeV, nd.estar.S_tot)
end


"""
    csda_range_g_per_cm2(nd::NISTData, T_MeV) -> Float64

CSDA range [g/cm²] for electrons of kinetic energy `T_MeV` in Xe,
computed by trapezoidal integration of 1/S_total over the ESTAR grid.
"""
function csda_range_g_per_cm2(nd::NISTData, T_MeV::Float64)::Float64
    interp_loglog(T_MeV, nd.estar.T_MeV, nd.estar.R_csda)
end


"""
    csda_range_LXe_mm(nd::NISTData, T_MeV) -> Float64

CSDA range in liquid xenon [mm].
"""
function csda_range_LXe_mm(nd::NISTData, T_MeV::Float64)::Float64
    csda_range_g_per_cm2(nd, T_MeV) / nd.rho * 10.0  # cm → mm
end


# =====================================================================
# Pre-tabulated bremsstrahlung cross section
# =====================================================================

"""
    sigma_brems_table(nd::NISTData, T_MeV) -> Float64

Total bremsstrahlung cross section [cm²/atom] for k > k_min at electron
kinetic energy `T_MeV`, from the pre-tabulated log-log interpolation.
Returns 0 for T below the table range.
"""
function sigma_brems_table(nd::NISTData, T_MeV::Float64)::Float64
    interp_loglog(T_MeV, nd.brems_T, nd.brems_σ)
end
