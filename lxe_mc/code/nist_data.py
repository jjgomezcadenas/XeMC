"""
nist_data.py
============

NIST XCOM (photon cross sections) and ESTAR (electron stopping power
and CSDA range) tables for xenon, plus log-log interpolators.

These are the *authoritative numerical values* used in the Monte Carlo
simulation and overlaid as markers in figures. The analytic formulas
in physics.py are cross-checked against these.

Usage:
    >>> from nist_data import (sigma_compton_NIST, sigma_pair_NIST,
    ...                        sigma_phot_NIST, dEdx_total_NIST,
    ...                        csda_range_LXe_NIST)
    >>> sigma_phot_NIST(1.5)         # cm^2/atom of Xe
    >>> dEdx_total_NIST(1.0)         # MeV cm^2/g
    >>> csda_range_LXe_NIST(1.0)     # mm in liquid Xe
"""

import os
import numpy as np
from physics import A_XE, NA, RHO_LXE, EK_XE_MEV

_HERE = os.path.dirname(os.path.abspath(__file__))
_DATA = os.path.normpath(os.path.join(_HERE, "..", "data"))


# =====================================================================
# XCOM loader
# =====================================================================

def _load_xcom():
    """Load NIST XCOM data for Xe; return dict of arrays (cm^2/g)."""
    rows = []
    with open(os.path.join(_DATA, "xcom_xe.csv")) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("E_MeV"):
                continue
            rows.append([float(x) for x in line.split(",")])
    a = np.array(rows)
    return {
        "E_MeV":         a[:, 0],
        "coherent":      a[:, 1],
        "incoherent":    a[:, 2],   # = Compton mass-attenuation
        "photoelectric": a[:, 3],
        "pair_nuclear":  a[:, 4],
        "pair_electron": a[:, 5],
        "total_w_coh":   a[:, 6],
        "total_no_coh":  a[:, 7],
    }

_XCOM = _load_xcom()


def _interp_loglog(x, xp, fp):
    """Linear interpolation of fp at x in log-log space.

    xp must be increasing positive. fp may have zeros (handled by
    falling back to linear interp where fp <= 0).
    Returns scalar if x is scalar, else array.
    """
    x_arr = np.atleast_1d(x).astype(float)
    out = np.zeros_like(x_arr)
    pos = fp > 0
    if pos.sum() >= 2:
        # Log-log interpolation in the positive region
        lx_pos = np.log(xp[pos])
        ly_pos = np.log(fp[pos])
        # only meaningful where x is within positive region
        for i, xi in enumerate(x_arr):
            if xi < xp[0]:
                out[i] = fp[0] if fp[0] > 0 else 0.0
            elif xi > xp[-1]:
                out[i] = fp[-1]
            elif xi < xp[pos][0]:
                # Below the first positive value: zero (e.g. pair below threshold)
                out[i] = 0.0
            else:
                out[i] = np.exp(np.interp(np.log(xi), lx_pos, ly_pos))
    else:
        out = np.interp(x_arr, xp, fp)
    return out if np.ndim(x) > 0 else float(out[0])


def _interp_xcom(channel, E_MeV, per_atom=True):
    """Interpolate one XCOM channel at energy E_MeV.

    Returns cross section in cm^2/atom (per_atom=True) or cm^2/g.
    Properly handles the K-edge discontinuity at 34.561 keV: the
    XCOM table has two entries at this energy (below-edge, above-edge);
    we treat them as a step.
    """
    Ex = _XCOM["E_MeV"].copy()
    Fx = _XCOM[channel].copy()
    # Handle the duplicate energy at K-edge: nudge the second one up
    for i in range(1, len(Ex)):
        if Ex[i] <= Ex[i-1]:
            Ex[i] = Ex[i-1] * (1.0 + 1e-9)
    out = _interp_loglog(E_MeV, Ex, Fx)   # cm^2/g, scalar or array
    if per_atom:
        out = out * (A_XE / NA)
    return out


# Public photon API ---------------------------------------------------

def sigma_compton_NIST(E_MeV, per_atom=True):
    """Compton (incoherent) cross section in Xe from NIST XCOM."""
    return _interp_xcom("incoherent", E_MeV, per_atom=per_atom)

def sigma_pair_NIST(E_MeV, per_atom=True):
    """Pair production (nuclear + electron) cross section in Xe from NIST XCOM."""
    nuc = _interp_xcom("pair_nuclear", E_MeV, per_atom=per_atom)
    ele = _interp_xcom("pair_electron", E_MeV, per_atom=per_atom)
    return nuc + ele

def sigma_phot_NIST(E_MeV, per_atom=True):
    """Photoelectric cross section in Xe from NIST XCOM."""
    return _interp_xcom("photoelectric", E_MeV, per_atom=per_atom)

def sigma_total_NIST(E_MeV, per_atom=True, include_coherent=False):
    """Total photon cross section in Xe from NIST XCOM."""
    chan = "total_w_coh" if include_coherent else "total_no_coh"
    return _interp_xcom(chan, E_MeV, per_atom=per_atom)

def mu_over_rho_NIST(E_MeV, include_coherent=False):
    """Mass attenuation coefficient mu/rho [cm^2/g] from NIST XCOM."""
    return sigma_total_NIST(E_MeV, per_atom=False, include_coherent=include_coherent)

def mfp_LXe_NIST(E_MeV, include_coherent=False):
    """Photon mean free path [cm] in liquid Xe from NIST XCOM."""
    mu = mu_over_rho_NIST(E_MeV, include_coherent=include_coherent) * RHO_LXE
    return 1.0 / mu

def branching_NIST(E_MeV):
    """Compton, pair, photoelectric branching fractions at energy E_MeV."""
    sC = sigma_compton_NIST(E_MeV, per_atom=False)
    sP = sigma_pair_NIST(E_MeV, per_atom=False)
    sPh = sigma_phot_NIST(E_MeV, per_atom=False)
    st = sC + sP + sPh
    return sC/st, sP/st, sPh/st


# =====================================================================
# ESTAR loader
# =====================================================================

def _load_estar():
    rows = []
    with open(os.path.join(_DATA, "estar_xe.csv")) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("T_MeV"):
                continue
            rows.append([float(x) for x in line.split(",")])
    a = np.array(rows)
    return {
        "T_MeV": a[:, 0],
        "S_col": a[:, 1],   # MeV cm^2/g
        "S_rad": a[:, 2],
        "S_tot": a[:, 3],
        "delta": a[:, 4],
    }

_ESTAR = _load_estar()


def dEdx_collision_NIST(T_MeV):
    """Mass collision stopping power [MeV cm^2/g] from NIST ESTAR for Xe."""
    return _interp_loglog(T_MeV, _ESTAR["T_MeV"], _ESTAR["S_col"])

def dEdx_radiative_NIST(T_MeV):
    """Mass radiative stopping power [MeV cm^2/g] from NIST ESTAR for Xe."""
    return _interp_loglog(T_MeV, _ESTAR["T_MeV"], _ESTAR["S_rad"])

def dEdx_total_NIST(T_MeV):
    """Total mass stopping power [MeV cm^2/g] from NIST ESTAR for Xe."""
    return _interp_loglog(T_MeV, _ESTAR["T_MeV"], _ESTAR["S_tot"])


# =====================================================================
# CSDA range (computed from ESTAR by integration of 1/S_tot)
# =====================================================================

def _build_csda_range_table():
    """Compute CSDA range R(T) [g/cm^2] by trapezoid integration of 1/S_tot
    over the ESTAR grid. R(T_min) defined as 0 at the lowest tabulated T.
    """
    T = _ESTAR["T_MeV"]
    S = _ESTAR["S_tot"]
    R = np.zeros_like(T)
    for i in range(1, len(T)):
        # Use log-log midpoint: dR = dT / S, with S evaluated at log-mean energy
        T_lo, T_hi = T[i-1], T[i]
        S_lo, S_hi = S[i-1], S[i]
        # Integrate 1/S using trapezoid in linear space; ESTAR grid is dense enough
        R[i] = R[i-1] + 0.5 * (1.0/S_lo + 1.0/S_hi) * (T_hi - T_lo)
    return T, R

_R_T_grid, _R_g_per_cm2 = _build_csda_range_table()


def csda_range_g_per_cm2_NIST(T_MeV):
    """CSDA range [g/cm^2] for electrons of kinetic energy T in Xe."""
    return _interp_loglog(T_MeV, _R_T_grid, _R_g_per_cm2)


def csda_range_LXe_NIST(T_MeV):
    """CSDA range in liquid Xe [mm] for electrons of kinetic energy T."""
    R_g = csda_range_g_per_cm2_NIST(T_MeV)
    return R_g / RHO_LXE * 10.0   # cm -> mm


if __name__ == "__main__":
    print("=== NIST XCOM photon cross sections in Xe ===")
    print(f"{'E[MeV]':>8} {'C [cm^2/g]':>12} {'P [cm^2/g]':>12} {'Ph [cm^2/g]':>13} {'mfp[cm]':>10}")
    for E in [0.05, 0.1, 0.5, 1.0, 1.5, 2.0, 2.4476, 2.4578, 2.6145, 3.0]:
        sC = sigma_compton_NIST(E, per_atom=False)
        sP = sigma_pair_NIST(E, per_atom=False)
        sPh = sigma_phot_NIST(E, per_atom=False)
        mfp = mfp_LXe_NIST(E)
        print(f"{E:8.4f} {sC:12.4e} {sP:12.4e} {sPh:13.4e} {mfp:10.3f}")

    print("\n=== Branching fractions (NIST) ===")
    for E, name in [(2.4476, "Bi-214"),
                    (2.4578, "Q_bb"),
                    (2.6145, "Tl-208")]:
        bC, bP, bPh = branching_NIST(E)
        print(f"{name} ({E} MeV): C={bC:.3f}  P={bP:.3f}  Ph={bPh:.4f}  Ph[%]={bPh*100:.2f}")

    print("\n=== Electron stopping power and CSDA range in LXe (NIST) ===")
    print(f"{'T[MeV]':>8} {'S_col':>10} {'S_rad':>10} {'rad%':>8} {'R[mm]':>10}")
    for T in [0.05, 0.1, 0.25, 0.5, 1.0, 1.5, 2.0, 2.4, 3.0]:
        Sc = dEdx_collision_NIST(T)
        Sr = dEdx_radiative_NIST(T)
        R = csda_range_LXe_NIST(T)
        print(f"{T:8.3f} {Sc:10.3f} {Sr:10.4f} {Sr/(Sc+Sr)*100:8.2f} {R:10.4f}")
