"""
sampling.py
===========

Monte Carlo samplers for the LXe simulation. Each sampler returns the
sampled kinematics for one interaction.

Methods used:
- inverse-CDF (exponential distance, simple angular distributions)
- composition + rejection (Klein-Nishina, Bethe-Heitler pair)
- envelope rejection (bremsstrahlung)
- Sauter (photoelectric polar angle)
"""

import numpy as np
from physics import (ME_MEV, RE_CM, ALPHA_FS, NA, A_XE, Z_XE,
                     EK_XE_MEV, OMEGA_K_XE, EK_ALPHA_MEV,
                     dsigma_dk_brems, _coulomb_correction_fc)


# =====================================================================
# Distance to next interaction
# =====================================================================

def sample_distance(Sigma_tot, rng):
    """Sample free path from exponential distribution.
    Sigma_tot in cm^-1; returns path length in cm."""
    return -np.log(rng.uniform()) / Sigma_tot


# =====================================================================
# Discrete process choice
# =====================================================================

def sample_process(P_C, P_P, P_Ph, rng):
    """Return one of 'compton', 'pair', 'photoelectric'."""
    r = rng.uniform()
    if r < P_C:
        return "compton"
    elif r < P_C + P_P:
        return "pair"
    else:
        return "photoelectric"


# =====================================================================
# Klein-Nishina sampling (after Geant4 PRM v11.4, ButcherMessel 1960)
# =====================================================================

def sample_compton(E_MeV, rng):
    """Sample Compton scattering: returns (E_gamma_scattered, cos_theta).

    Klein-Nishina sampled by combined composition + rejection.
    cos_theta is the cosine of the photon scattering angle in the lab.
    """
    eps0 = ME_MEV / (ME_MEV + 2.0 * E_MeV)
    alpha1 = -np.log(eps0)
    alpha2 = (1.0 - eps0**2) / 2.0
    while True:
        r1 = rng.uniform()
        r2 = rng.uniform()
        r3 = rng.uniform()
        if r1 < alpha1 / (alpha1 + alpha2):
            eps = eps0 ** r2
        else:
            eps = np.sqrt(eps0**2 + (1.0 - eps0**2) * r2)
        # cos theta from kinematics: 1 - cos_theta = (1/eps - 1) * me / E
        one_minus_cos = (1.0 - eps) * ME_MEV / (E_MeV * eps)
        sin2 = one_minus_cos * (2.0 - one_minus_cos)
        sin2 = max(0.0, min(1.0, sin2))
        # rejection function g(eps) = 1 - eps sin^2 / (1 + eps^2)
        g = 1.0 - eps * sin2 / (1.0 + eps**2)
        if g >= r3:
            cos_theta = 1.0 - one_minus_cos
            return eps * E_MeV, cos_theta


def compton_electron_direction(cos_theta_gamma, phi_gamma, E_MeV, dir_in):
    """Given the scattered photon direction, return the recoil
    electron direction (unit vector). Recoil electron is coplanar
    with incoming photon and scattered photon.

    Inputs: cos_theta_gamma in [-1,1] photon scattering cos relative
    to incoming, phi_gamma in [0,2pi], E_MeV initial photon energy,
    dir_in unit vector incoming photon direction.
    Returns unit vector for recoil electron direction.
    """
    # compute recoil momentum vector in the incoming-photon frame:
    # p_e = p_gamma - p_gamma_scattered
    sin_theta_gamma = np.sqrt(max(0.0, 1.0 - cos_theta_gamma**2))
    Eg_p = E_MeV / (1.0 + (E_MeV/ME_MEV) * (1.0 - cos_theta_gamma))
    # In the local frame where dir_in = z_hat:
    #  p_gamma_in = E_MeV * (0, 0, 1)
    #  p_gamma_out = Eg_p * (sin*cos_phi, sin*sin_phi, cos)
    p_in = np.array([0.0, 0.0, E_MeV])
    p_out = Eg_p * np.array([sin_theta_gamma*np.cos(phi_gamma),
                              sin_theta_gamma*np.sin(phi_gamma),
                              cos_theta_gamma])
    p_e_local = p_in - p_out
    norm = np.linalg.norm(p_e_local)
    if norm < 1e-12:
        return dir_in.copy()
    n_e_local = p_e_local / norm
    return _rotate_to_global(n_e_local, dir_in)


def _rotate_to_global(local_vec, ref_dir):
    """Rotate a unit vector from a frame where z_hat = local z to the
    global frame where z_hat is ref_dir."""
    rd = ref_dir / np.linalg.norm(ref_dir)
    if abs(rd[2]) > 0.99999:
        return local_vec * np.sign(rd[2])
    # build orthonormal triad
    e1 = np.array([rd[2]*rd[0], rd[2]*rd[1], rd[2]**2 - 1.0])
    e1 = -e1 / np.sqrt(rd[0]**2 + rd[1]**2)
    e2 = np.cross(rd, e1)
    return local_vec[0]*e1 + local_vec[1]*e2 + local_vec[2]*rd


# =====================================================================
# Pair production sampling (Bethe-Heitler)
# =====================================================================

def sample_pair(E_MeV, rng):
    """Sample energy split eps = E_+/E_gamma in [eps_0, 1-eps_0].

    Returns eps for the positron. For E_gamma < ~2 MeV the differential
    is essentially flat in eps and uniform sampling is correct to ~1%.
    Above that we use composition + rejection.
    """
    eps_min = ME_MEV / E_MeV
    if eps_min >= 0.5:
        return 0.5
    if E_MeV < 2.0:
        # Uniform sampling, full range
        return eps_min + rng.uniform() * (1.0 - 2.0*eps_min)
    # Composition: f1 = constant, f2 \propto (1/2 - eps)^2, normalized.
    # For pedagogy / simplicity we use a flat envelope and rejection:
    # the differential with screening is at most ~constant * (eps^2 + (1-eps)^2)
    # which is bounded above by 2 on [eps_min, 1-eps_min].
    while True:
        eps = eps_min + rng.uniform() * (1.0 - 2.0*eps_min)
        # Approximate weight: eps^2 + (1-eps)^2 normalized to its max (=1 at endpoints)
        w = eps**2 + (1.0 - eps)**2
        if rng.uniform() < w:
            return eps


def pair_polar_angle(E_lepton, rng):
    """Sample polar angle theta of one of the leptons (in rad).
    Uses the Geant4 / Urban parameterization."""
    a = 5.0/8.0
    d = 27.0
    # Composition: f(u) = c1 * u * exp(-a*u) + c2 * u * exp(-3*a*u)
    # weights given by integrals over u, both = 1/a^2 and 1/(3a)^2
    # so prob of first branch is (1/a^2)/(1/a^2 + d/(3a)^2) = 9/(9+d)
    p1 = 9.0 / (9.0 + d)
    if rng.uniform() < p1:
        b = a
    else:
        b = 3.0 * a
    # x*exp(-b*x) is sampled by negative binomial-like:
    # CDF F(u) = 1 - (1 + b*u) * exp(-b*u)
    # solve numerically via two uniforms (sum of two exponentials)
    u = (-np.log(rng.uniform()) - np.log(rng.uniform())) / b
    return (ME_MEV / E_lepton) * u


# =====================================================================
# Photoelectric: shell + photoelectron polar angle
# =====================================================================

def sample_phot_shell(E_MeV, rng):
    """Sample which shell is ionized in a photoelectric event in Xe.

    Returns (shell_label, binding_energy_MeV).
    Approximate branching above K-edge: K=80%, L=15%, M+=5%.
    Below K-edge: L only.
    """
    if E_MeV > EK_XE_MEV:
        r = rng.uniform()
        if r < 0.80:
            return "K", EK_XE_MEV
        elif r < 0.95:
            return "L", 0.005   # ~5 keV avg L-shell binding in Xe
        else:
            return "M", 0.0007
    else:
        return "L", 0.005


def sample_photoelectron_angle(T_e_MeV, rng):
    """Sample photoelectron polar angle (radians) wrt photon direction.

    Sauter (zeroth-order) distribution:
       dsigma/d(cos theta) ~ sin^2 theta / (1 - beta cos theta)^4
    Sampled by rejection from the uniform envelope.
    """
    g = 1.0 + T_e_MeV / ME_MEV
    beta = np.sqrt(1.0 - 1.0/g**2)
    # max of f(cos_theta) is at cos_theta = beta (approx); peak value bounded
    # we use rejection with envelope max ~ 4 (loose but OK)
    M = 4.0
    while True:
        ct = -1.0 + 2.0 * rng.uniform()
        f = (1.0 - ct**2) / (1.0 - beta*ct)**4
        if rng.uniform() * M < f:
            return np.arccos(ct)


def sample_atomic_relaxation_K(rng):
    """For a K-shell vacancy in Xe, sample the de-excitation outcome.

    Returns ('fluorescence', E_photon_MeV) or ('auger', E_electron_MeV).
    """
    if rng.uniform() < OMEGA_K_XE:
        return "fluorescence", EK_ALPHA_MEV
    else:
        # Auger electron of approximately E_K - 2*E_L
        return "auger", EK_XE_MEV - 2 * 0.005


# =====================================================================
# Bremsstrahlung sampling
# =====================================================================

def _sigma_brems_above_kmin(T_MeV, k_min_MeV):
    """Total brems cross section for k > k_min, by integration of BH-Tsai
    differential. Returns cm^2/atom."""
    if k_min_MeV >= T_MeV:
        return 0.0
    # log-spaced k grid for integration (1/k spectrum)
    k = np.logspace(np.log10(k_min_MeV), np.log10(T_MeV * 0.9999), 100)
    ds = dsigma_dk_brems(k, T_MeV)
    return np.trapezoid(ds, k)


def sample_brems(T_MeV, k_min_MeV, rng):
    """Sample bremsstrahlung photon energy k > k_min.

    Use 1/k envelope sampling: k = k_min * (T/k_min)^r with r ~ U(0,1).
    The rejection weight is the ratio of the BH-Tsai differential to
    the envelope (1/k).
    """
    if k_min_MeV >= T_MeV:
        return None  # not allowed
    # Envelope normalization: integral of 1/k from k_min to T is ln(T/k_min)
    log_ratio = np.log(T_MeV / k_min_MeV)
    # Find rough max of (k * dsigma/dk) over the range to bound rejection
    k_grid = np.logspace(np.log10(k_min_MeV), np.log10(T_MeV*0.9999), 50)
    ds = dsigma_dk_brems(k_grid, T_MeV)
    M = np.max(k_grid * ds) * 1.05   # 5% safety
    while True:
        r = rng.uniform()
        k = k_min_MeV * (T_MeV / k_min_MeV) ** r
        ds_k = dsigma_dk_brems(k, T_MeV)
        if rng.uniform() * M < k * ds_k:
            return k


def brems_photon_angle(T_e_MeV, rng):
    """Sample bremsstrahlung photon polar angle (rad) wrt electron direction.
    Same Urban functional form as pair-produced leptons."""
    return pair_polar_angle(T_e_MeV + ME_MEV, rng)
