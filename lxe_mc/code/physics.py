"""
physics.py
==========

Analytic electromagnetic cross-section formulas for photons and electrons
in xenon, used by the pedagogical LXe Monte Carlo.

Each function returns a *numerical value* with documented units. Where two
versions exist (e.g. Sandia photoelectric vs NIST XCOM interpolation), use
nist_data.py for the interpolation and physics.py for the analytic formula;
the two are then compared in figures and validation tests.

Formulas implemented:
- sigma_compton       : Klein-Nishina, Z*KN per atom (above K-edge)
- sigma_pair          : Bethe-Heitler with screening + Coulomb correction
                        (Geant4 PRM v11.4 form)
- sigma_phot_sandia   : Sandia parameterization (Biggs-Lighthill 1990)
                        with Xe coefficients from data/sandia_xe_coeffs.csv
- dEdx_collision      : Berger-Seltzer mass collision stopping power (ICRU-37)
- dEdx_radiative      : Tsai radiative stopping power (complete-screening)
- dsig_dk_brems       : Bethe-Heitler bremsstrahlung differential
                        cross section (Tsai)

All cross sections returned as cm^2/atom unless stated. Energies in MeV.
"""

import numpy as np
import os

# =====================================================================
# Physical constants  (PDG 2024, CODATA-consistent)
# =====================================================================
ME_MEV   = 0.5109989461          # electron mass energy [MeV]
RE_CM    = 2.8179403262e-13      # classical electron radius [cm]
ALPHA_FS = 1.0 / 137.035999084   # fine-structure constant
NA       = 6.02214076e23         # Avogadro [/mol]
BARN_CM2 = 1.0e-24               # 1 barn = 1e-24 cm^2

# =====================================================================
# Xenon and LXe parameters
# =====================================================================
Z_XE          = 54
A_XE          = 131.293          # standard atomic weight
RHO_LXE       = 2.953            # g/cm^3 at ~165 K, 1 atm
NAT_LXE       = NA * RHO_LXE / A_XE   # atoms / cm^3 = 1.354e22

# Atomic data
EK_XE_MEV     = 0.034561         # K-shell binding energy
EL1_XE_MEV    = 0.005453
EL2_XE_MEV    = 0.005104
EL3_XE_MEV    = 0.004782
OMEGA_K_XE    = 0.89             # K-shell fluorescence yield
EK_ALPHA_MEV  = 0.02978          # dominant K-alpha line
I_XE_MEV      = 482e-6           # ICRU-37 mean excitation energy

# =====================================================================
# Compton scattering: Klein-Nishina
# =====================================================================
def sigma_compton(E_MeV):
    """Total Compton cross section per atom of Xe [cm^2].

    Above the K-edge: sigma_C = Z * sigma_KN, with sigma_KN the
    Klein-Nishina cross section per electron. The free-electron
    approximation is good to better than 1% for our energies
    (E_gamma >> E_K).
    """
    x = np.asarray(E_MeV) / ME_MEV
    pre = 2.0 * np.pi * RE_CM**2
    t1 = (1.0 + x) / x**2 * (2.0*(1.0 + x)/(1.0 + 2.0*x) - np.log(1.0 + 2.0*x)/x)
    t2 = np.log(1.0 + 2.0*x) / (2.0*x)
    t3 = -(1.0 + 3.0*x) / (1.0 + 2.0*x)**2
    sigma_KN = pre * (t1 + t2 + t3)
    return Z_XE * sigma_KN

def dsigma_dOmega_KN(theta_rad, E_MeV):
    """Klein-Nishina differential cross section per electron [cm^2/sr]."""
    a = E_MeV / ME_MEV
    ct = np.cos(theta_rad)
    Eg_p = E_MeV / (1.0 + a*(1.0 - ct))
    r = Eg_p / E_MeV
    return 0.5 * RE_CM**2 * r**2 * (r + 1.0/r - np.sin(theta_rad)**2)

def compton_kinematics(E_MeV, theta_rad):
    """Return (E_gamma_scattered, T_e, phi_e) for given incident E and angle.

    phi_e is the recoil-electron polar angle wrt the incoming photon
    direction.
    """
    a = E_MeV / ME_MEV
    Egp = E_MeV / (1.0 + a*(1.0 - np.cos(theta_rad)))
    Te = E_MeV - Egp
    cos_phi_e = (1.0 + a) / np.sqrt(1.0 + 2.0*a) * np.tan(theta_rad/2.0)
    phi_e = np.arctan(np.sqrt(1.0 - cos_phi_e**2) / np.maximum(cos_phi_e, 1e-12))
    return Egp, Te, phi_e

# =====================================================================
# Pair production: Bethe-Heitler with screening
# =====================================================================
def _bh_screening_phi(delta):
    """Geant4 PRM piecewise screening functions Phi_1, Phi_2."""
    delta = np.atleast_1d(delta).astype(float)
    p1 = np.where(
        delta <= 1.0,
        20.867 - 3.242*delta + 0.625*delta**2,
        21.12 - 4.184*np.log(delta + 0.952))
    p2 = np.where(
        delta <= 1.0,
        20.209 - 1.930*delta - 0.086*delta**2,
        21.12 - 4.184*np.log(delta + 0.952))
    return p1, p2

def _coulomb_correction_fc(Z):
    """Davies-Bethe-Maximon Coulomb correction f_c(Z)."""
    a = ALPHA_FS * Z
    return a**2 * (1.0/(1.0 + a**2) + 0.20206 - 0.0369*a**2
                   + 0.0083*a**4 - 0.0020*a**6)

def sigma_pair(E_MeV, calibration=1.0):
    """Pair production cross section per atom of Xe [cm^2].

    Numerical integration of the Geant4-PRM differential cross section
    (Bethe-Heitler with screening + Coulomb correction). Below 1.5 MeV
    the bare formula underestimates NIST XCOM by ~50% in Xe due to the
    near-threshold regime; this function takes an optional calibration
    factor (default 1.0 = no calibration).

    For accurate numbers below 5 MeV use nist_data.sigma_pair_NIST(E).
    """
    Eg = np.atleast_1d(E_MeV).astype(float)
    out = np.zeros_like(Eg)
    fc = _coulomb_correction_fc(Z_XE)
    F_Z = (8.0/3.0) * np.log(Z_XE)
    xi_Z = np.log(1440.0 / Z_XE**(2.0/3.0)) / (np.log(183.0/Z_XE**(1.0/3.0)) - fc)
    pre = ALPHA_FS * RE_CM**2 * Z_XE * (Z_XE + xi_Z)
    for i, Eg_i in enumerate(Eg):
        if Eg_i <= 2.0 * ME_MEV:
            out[i] = 0.0
            continue
        eps0 = ME_MEV / Eg_i
        if eps0 >= 0.5:
            out[i] = 0.0
            continue
        eps = np.linspace(eps0 + 1e-5, 0.5 - 1e-5, 300)
        delta = (136.0 / Z_XE**(1.0/3.0)) * eps0 / (eps * (1.0 - eps))
        phi1, phi2 = _bh_screening_phi(delta)
        f1 = np.maximum(phi1 - F_Z/2.0, 0.0)
        f2 = np.maximum(phi2 - F_Z/2.0, 0.0)
        integrand = ((eps**2 + (1.0 - eps)**2) * f1
                     + (2.0/3.0) * eps * (1.0 - eps) * f2)
        sigma = 2.0 * pre * np.trapezoid(integrand, eps)
        out[i] = max(0.0, calibration * sigma)
    return out if Eg.size > 1 else float(out[0])

# =====================================================================
# Photoelectric absorption: Sandia parameterization
# =====================================================================
def _load_sandia_table_xe():
    """Load Xe Sandia coefficients from data/sandia_xe_coeffs.csv.

    Returns array shape (N_intervals, 5) with columns
    [E_lower_keV, a1, a2, a3, a4]. Coefficients in cm^2 keV^i / g.
    """
    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(here, "..", "data", "sandia_xe_coeffs.csv")
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("E_lower"):
                continue
            parts = line.split(",")
            rows.append([float(x) for x in parts])
    return np.array(rows)

_SANDIA_XE = _load_sandia_table_xe()

def sigma_phot_sandia(E_MeV, per_atom=True):
    """Photoelectric cross section in Xe via Sandia parameterization.

    sigma_phot/rho [cm^2/g] = a1/E + a2/E^2 + a3/E^3 + a4/E^4
    with E in keV and a_i in cm^2 keV^i / g.

    By default returns per-atom cross section in cm^2 (multiply by
    rho * N_A / A internally, then convert by A/N_A: net effect is
    multiply mu/rho by A/N_A).

    *** VALIDITY DOMAIN ***
    Sandia agrees with NIST XCOM to ~1% below 500 keV in Xe. The
    high-energy interval [500 keV, infinity) is fit to a small number
    of points and extrapolates BADLY at MeV energies in Xe:
        E = 1 MeV : Sandia 1.56e-2,  XCOM 4.34e-3 cm^2/g  (off 3.6x)
        E = 3 MeV : Sandia 1.03e-3,  XCOM 6.42e-4 cm^2/g  (off 1.6x)
    Use sigma_phot_xcom() for E > 500 keV in production code.
    """
    Eg = np.atleast_1d(E_MeV).astype(float)
    out = np.zeros_like(Eg)
    E_keV = Eg * 1000.0
    edges = _SANDIA_XE[:, 0]
    coefs = _SANDIA_XE[:, 1:]
    for i, e_keV in enumerate(E_keV):
        # find the interval that contains e_keV: largest edge <= e_keV
        idx = np.searchsorted(edges, e_keV, side="right") - 1
        if idx < 0:
            out[i] = 0.0
            continue
        a1, a2, a3, a4 = coefs[idx]
        out[i] = a1/e_keV + a2/e_keV**2 + a3/e_keV**3 + a4/e_keV**4  # cm^2/g
        if out[i] < 0:
            out[i] = 0.0
    if per_atom:
        out = out * (A_XE / NA)   # cm^2/g -> cm^2/atom
    return out if Eg.size > 1 else float(out[0])

# =====================================================================
# Electron stopping power: Berger-Seltzer (collisional) + Tsai (radiative)
# =====================================================================
def _beta_gamma(T_MeV):
    g = 1.0 + np.asarray(T_MeV) / ME_MEV
    b = np.sqrt(1.0 - 1.0/g**2)
    return b, g

def dEdx_collision(T_MeV):
    """Mass collision stopping power for electrons in Xe [MeV cm^2/g].

    Berger-Seltzer (ICRU-37) with the Moller cross section for
    indistinguishable identical particles. Density-effect correction
    delta is set to zero (negligible below ~20 MeV in LXe).
    """
    T = np.atleast_1d(T_MeV).astype(float)
    b, g = _beta_gamma(T)
    tau = T / ME_MEV
    Fminus = (1.0 - b**2 +
              (tau**2/8.0 - (2.0*tau + 1.0)*np.log(2.0)) / (tau + 1.0)**2)
    K = 0.307075   # MeV cm^2/g  ( = 4 pi N_A r_e^2 m_e c^2 )
    delta = 0.0
    arg = (T / I_XE_MEV)**2 * (tau + 2.0) / 2.0
    S = (K/2.0) * (Z_XE / A_XE) / b**2 * (np.log(arg) + Fminus - delta)
    return S if T.size > 1 else float(S[0])

def dEdx_radiative(T_MeV, calibration=1.0):
    """Mass radiative stopping power for electrons in Xe [MeV cm^2/g].

    Tsai complete-screening form for bremsstrahlung losses:
       S_rad/rho = 4 alpha N_A/A Z(Z+1) r_e^2 (T+m_e) L_rad
    where L_rad = ln(183/Z^(1/3)) - f_c.

    The bare formula needs a low-energy turn-on factor; we use a
    smooth one tuned to match NIST ESTAR at T = 1 MeV in Xe to ~10%.
    For accurate numbers use nist_data.dEdx_radiative_NIST(T).
    """
    T = np.atleast_1d(T_MeV).astype(float)
    Etot = T + ME_MEV
    fc = _coulomb_correction_fc(Z_XE)
    Lrad = np.log(183.0 / Z_XE**(1.0/3.0)) - fc
    S_asymp = 4.0 * ALPHA_FS * (NA/A_XE) * Z_XE*(Z_XE + 1.0) * RE_CM**2 * Etot * Lrad
    factor_low = 1.0 - np.exp(-T / 0.4)
    S = S_asymp * factor_low * 0.43 * calibration
    return S if T.size > 1 else float(S[0])

# =====================================================================
# Bremsstrahlung differential cross section: Bethe-Heitler-Tsai
# =====================================================================
def dsigma_dk_brems(k_MeV, T_MeV):
    """Differential bremsstrahlung cross section in Xe [cm^2 / MeV / atom].

    Bethe-Heitler-Tsai with screening:
       dsigma/dk = 4 alpha r_e^2 / (3 k) {
          [y^2 + 2(1 + (1-y)^2)] [ Z^2 (F_el - f_c) + Z F_inel ]
          + (1-y) (Z^2 + Z) / 3
       }
    with y = k / (T + m_e). Form factors:
       F_el   = ln(184.15 / Z^(1/3))
       F_inel = ln(1194 / Z^(2/3))

    Accuracy ~10% in MeV range; for production use Seltzer-Berger
    tabulated cross sections (Geant4 G4SeltzerBergerModel).
    """
    Etot = T_MeV + ME_MEV
    k = np.atleast_1d(k_MeV).astype(float)
    out = np.zeros_like(k)
    mask = (k > 0) & (k < T_MeV)
    y = k[mask] / Etot
    fc = _coulomb_correction_fc(Z_XE)
    Fel = np.log(184.15 / Z_XE**(1.0/3.0))
    Finel = np.log(1194.0 / Z_XE**(2.0/3.0))
    out[mask] = (4.0 * ALPHA_FS * RE_CM**2 / (3.0 * k[mask])) * (
        (y**2 + 2.0*(1.0 + (1.0 - y)**2)) * (Z_XE**2*(Fel - fc) + Z_XE*Finel)
        + (1.0 - y) * (Z_XE**2 + Z_XE) / 3.0)
    return out if k.size > 1 else float(out[0])

# =====================================================================
# Total cross section per atom (sum of three channels)
# =====================================================================
def sigma_total(E_MeV, use_sandia=True, sandia_max_MeV=0.5):
    """Total photon cross section per atom of Xe [cm^2].

    Compton (always KN), pair (BH numerical), photoelectric: Sandia
    below sandia_max_MeV (default 500 keV), zero above (since Sandia
    extrapolates poorly there). For accurate photoelectric above 500
    keV use nist_data interpolation.
    """
    sC = sigma_compton(E_MeV)
    sP = sigma_pair(E_MeV)
    sPh = np.zeros_like(np.atleast_1d(E_MeV).astype(float))
    if use_sandia:
        Eg = np.atleast_1d(E_MeV).astype(float)
        below = Eg <= sandia_max_MeV
        if np.any(below):
            sPh_below = sigma_phot_sandia(Eg[below])
            sPh[below] = sPh_below
    return sC + sP + (sPh if sPh.size > 1 else float(sPh[0]))


if __name__ == "__main__":
    # sanity print
    print("Compton, pair, Sandia-phot per atom of Xe (in barns):")
    print(f"{'E[MeV]':>8} {'C [b]':>10} {'P [b]':>10} {'Ph_Sandia [b]':>15}")
    for E in [0.05, 0.1, 0.5, 1.0, 1.5, 2.0, 2.4476, 2.6145, 3.0]:
        sC = sigma_compton(E) / BARN_CM2
        sP = sigma_pair(E) / BARN_CM2
        sPh = sigma_phot_sandia(E) / BARN_CM2
        print(f"{E:8.4f} {sC:10.3f} {sP:10.3f} {sPh:15.4f}")

    print("\nElectron stopping power in Xe [MeV cm^2/g] (Berger-Seltzer / Tsai):")
    print(f"{'T[MeV]':>8} {'S_col':>10} {'S_rad':>10} {'frac_rad':>10}")
    for T in [0.1, 0.5, 1.0, 1.5, 2.0, 2.4]:
        Sc = dEdx_collision(T)
        Sr = dEdx_radiative(T)
        print(f"{T:8.3f} {Sc:10.3f} {Sr:10.4f} {Sr/(Sc+Sr):10.4f}")
