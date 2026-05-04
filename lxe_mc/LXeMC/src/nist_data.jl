"""
NIST data loaders and interpolation utilities.

Provides CSV loaders for NIST XCOM (photon cross sections) and ESTAR
(electron stopping powers) tables, plus log-log interpolation routines
used throughout the simulation.

# Data structures
- [`XCOMData`](@ref): XCOM energy grid and 7 cross-section columns [cm²/g].
- [`ESTARData`](@ref): ESTAR energy grid, stopping powers [MeV cm²/g],
  and CSDA range [g/cm²] computed by numerical integration of 1/S_total.

# Interpolation
- [`interp_loglog`](@ref): general-purpose log-log interpolation (with search).
- [`interp_loglog_prelogged`](@ref): fast variant using pre-computed log arrays
  and a known bracket index (avoids `searchsortedlast` and `log` on grid data).
- [`prelog_data`](@ref): utility to pre-compute log arrays at load time.
"""


# =====================================================================
# Data structures
# =====================================================================

"""
    XCOMData

NIST XCOM photon cross sections. All columns are `Vector{Float64}`
in units of cm²/g (mass attenuation coefficients µ/ρ).
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

NIST ESTAR electron stopping powers and derived CSDA range.
Stopping powers in MeV cm²/g; CSDA range in g/cm².
"""
struct ESTARData
    T_MeV::Vector{Float64}
    S_col::Vector{Float64}
    S_rad::Vector{Float64}
    S_tot::Vector{Float64}
    delta::Vector{Float64}
    R_csda::Vector{Float64}   # g/cm², by trapezoidal integration of 1/S_tot
end


# =====================================================================
# CSV loaders
# =====================================================================

"""
    load_xcom(path) -> XCOMData

Read a NIST XCOM table in the standard space-delimited format
(2 header lines + blank line + data rows with 8 numeric columns).
Columns: E [MeV], coherent, incoherent, photoelectric, pair_nuclear,
pair_electron, total_w_coh, total_no_coh — all in cm²/g.

K-edge duplicate energies are preserved; the interpolator handles
them via [`_prepare_xcom_energy`](@ref).
"""
function load_xcom(path::AbstractString)::XCOMData
    rows = Vector{Vector{Float64}}()
    open(path, "r") do io
        for line in eachline(io)
            s = strip(line)
            isempty(s) && continue
            c = s[1]
            (c >= '0' && c <= '9') || continue
            fields = split(s)
            length(fields) >= 8 || continue
            push!(rows, parse.(Float64, fields[1:8]))
        end
    end
    m = reduce(hcat, rows)'
    XCOMData(m[:,1], m[:,2], m[:,3], m[:,4], m[:,5], m[:,6], m[:,7], m[:,8])
end


"""
    load_estar(path) -> ESTARData

Read a NIST ESTAR table (comma- or space-delimited, 5 columns:
T [MeV], S_col, S_rad, S_tot, delta). Lines not starting with a
digit are skipped. CSDA range R(T) is computed by trapezoidal
integration of 1/S_total.
"""
function load_estar(path::AbstractString)::ESTARData
    rows = Vector{Vector{Float64}}()
    open(path, "r") do io
        for line in eachline(io)
            s = strip(line)
            isempty(s) && continue
            c = s[1]
            (c >= '0' && c <= '9') || continue
            fields = occursin(",", s) ? split(s, ",") : split(s)
            length(fields) >= 5 || continue
            push!(rows, parse.(Float64, fields[1:5]))
        end
    end
    m = reduce(hcat, rows)'
    T = m[:,1]; Sc = m[:,2]; Sr = m[:,3]; St = m[:,4]; delta = m[:,5]

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

Log-log linear interpolation of `fp` at `x`, given grid `xp` (must be
increasing, positive). Handles zero/negative `fp` values by returning 0
in the zero region.

Uses `searchsortedlast` for the bracket and computes `log`/`exp` on the
fly. For hot paths, prefer [`interp_loglog_prelogged`](@ref) which
avoids these costs.
"""
function interp_loglog(x::Float64, xp::Vector{Float64}, fp::Vector{Float64})::Float64
    x <= 0.0 && return 0.0
    n = length(xp)
    x < xp[1] && return (fp[1] > 0.0 ? fp[1] : 0.0)
    x > xp[end] && return fp[end]

    lo = clamp(searchsortedlast(xp, x), 1, n - 1)
    hi = lo + 1

    (fp[lo] <= 0.0 && fp[hi] <= 0.0) && return 0.0
    fp[lo] <= 0.0 && return 0.0
    fp[hi] <= 0.0 && return fp[lo]

    lx = log(x)
    t = (lx - log(xp[lo])) / (log(xp[hi]) - log(xp[lo]))
    exp(log(fp[lo]) + t * (log(fp[hi]) - log(fp[lo])))
end


"""
    interp_loglog_prelogged(lx, log_xp, log_fp, fp, lo) -> Float64

Fast log-log interpolation using pre-computed log arrays and a known
bracket index `lo` (from a prior `searchsortedlast`). Only one `exp`
call is needed; all `log` calls on grid data are eliminated.

This is the hot-path interpolator used by [`sigma_three`](@ref),
[`dEdx_collision`](@ref), [`sigma_brems`](@ref), and
[`brems_rejection_M`](@ref).

# Arguments
- `lx`: log of the query point (caller computes once, reuses for multiple channels)
- `log_xp`: pre-computed `log.(xp)`
- `log_fp`: pre-computed `log.(fp)` (with `-Inf` for non-positive entries)
- `fp`: raw data (for zero-checking)
- `lo`: bracket index such that `xp[lo] ≤ x < xp[lo+1]`
"""
function interp_loglog_prelogged(lx::Float64, log_xp::Vector{Float64},
                                  log_fp::Vector{Float64}, fp::Vector{Float64},
                                  lo::Int)::Float64
    hi = lo + 1
    @inbounds begin
        fp_lo = fp[lo]
        fp_hi = fp[hi]
    end

    (fp_lo <= 0.0 && fp_hi <= 0.0) && return 0.0
    fp_lo <= 0.0 && return 0.0
    fp_hi <= 0.0 && return fp_lo

    @inbounds begin
        t = (lx - log_xp[lo]) / (log_xp[hi] - log_xp[lo])
        exp(log_fp[lo] + t * (log_fp[hi] - log_fp[lo]))
    end
end


"""
    prelog_data(fp) -> Vector{Float64}

Compute `log.(fp)`, mapping non-positive values to `-Inf`.
Used at load time to build pre-logged arrays for
[`interp_loglog_prelogged`](@ref).
"""
function prelog_data(fp::Vector{Float64})::Vector{Float64}
    [f > 0.0 ? log(f) : -Inf for f in fp]
end


# =====================================================================
# XCOM energy grid preparation
# =====================================================================

"""
    _prepare_xcom_energy(xc::XCOMData) -> Vector{Float64}

Return a copy of the XCOM energy grid with duplicate K-edge entries
nudged apart (by a factor of 1+10⁻⁹) so that `searchsortedlast`
returns the correct bracket across the discontinuity.
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
