"""
tracking.py
===========

Stack-based Monte Carlo tracking for photons and leptons in a large
LXe calorimeter.

Design follows Geant4: a single LIFO particle stack, processes are
sampled at each interaction, secondaries are pushed onto the stack.
Charged particle transport is condensed-history with explicit hard
bremsstrahlung sampling.

Public API:
    simulate_event(E_MeV, position, direction, geometry, rng)
        -> list of (x, y, z, dE) energy deposits

The geometry is currently a simple infinite-LXe (no boundary check)
or a cylindrical fiducial volume.
"""

from dataclasses import dataclass, field
from typing import List, Tuple
import numpy as np

from physics import (ME_MEV, NAT_LXE, RHO_LXE, EK_XE_MEV,
                     sigma_compton, dsigma_dk_brems)
from nist_data import (sigma_compton_NIST, sigma_pair_NIST, sigma_phot_NIST,
                       sigma_total_NIST, dEdx_collision_NIST,
                       csda_range_g_per_cm2_NIST)
import sampling
import thresholds as th


# =====================================================================
# Geometry: simple LXe-everywhere (infinite) or a cylinder
# =====================================================================

@dataclass
class InfiniteLXe:
    """Infinite LXe (no boundaries)."""
    def is_inside(self, pos): return True

@dataclass
class CylinderLXe:
    """Cylinder along z, radius R, half-height H, centered at origin."""
    radius_cm: float
    half_height_cm: float
    def is_inside(self, pos):
        return (pos[0]**2 + pos[1]**2 < self.radius_cm**2
                and abs(pos[2]) < self.half_height_cm)


# =====================================================================
# Track and stack
# =====================================================================

@dataclass
class Track:
    kind: str             # 'gamma', 'electron', 'positron'
    energy: float         # MeV (kinetic for leptons; total for gamma)
    position: np.ndarray  # cm
    direction: np.ndarray # unit vector
    parent_id: int = -1
    generation: int = 0


class Stack:
    """Plain LIFO stack of Track objects."""
    def __init__(self):
        self._items: List[Track] = []
    def push(self, t: Track):
        self._items.append(t)
    def pop(self) -> Track:
        return self._items.pop()
    def empty(self) -> bool:
        return len(self._items) == 0
    def __len__(self):
        return len(self._items)


@dataclass
class Deposit:
    """An energy deposit in the LXe."""
    position: np.ndarray   # cm
    energy: float          # MeV
    source: str            # 'electron', 'positron', 'photoelectric',
                           # 'auger', 'gamma_local', 'brems_local'


# =====================================================================
# Photon transport
# =====================================================================

def transport_photon(track: Track, geom, deposits: List[Deposit],
                     stack: Stack, rng: np.random.Generator) -> None:
    """Transport one gamma until escape, absorption, or below cutoff.

    May push additional secondaries (electrons, positrons, photons) onto
    the stack. Energy deposits are appended to deposits.
    """
    pos = track.position.copy()
    direction = track.direction.copy()
    E = track.energy

    while E >= th.EGAMMA_CUT_MEV:
        # NIST cross sections (per atom of Xe in cm^2)
        sC = sigma_compton_NIST(E, per_atom=True)
        sP = sigma_pair_NIST(E, per_atom=True)
        sPh = sigma_phot_NIST(E, per_atom=True)
        s_tot = sC + sP + sPh
        Sigma_tot = NAT_LXE * s_tot

        # Sample distance to next interaction
        s = sampling.sample_distance(Sigma_tot, rng)
        pos = pos + direction * s
        if not geom.is_inside(pos):
            return  # escapes

        # Sample which process
        proc = sampling.sample_process(sC/s_tot, sP/s_tot, sPh/s_tot, rng)

        if proc == "compton":
            # sample (Egp, cos_theta), update direction, push electron
            Egp, cos_t = sampling.sample_compton(E, rng)
            phi = 2*np.pi*rng.uniform()
            n_e = sampling.compton_electron_direction(cos_t, phi, E, direction)
            T_e = E - Egp
            # push electron
            stack.push(Track(kind="electron", energy=T_e, position=pos.copy(),
                             direction=n_e, parent_id=id(track),
                             generation=track.generation + 1))
            # update photon direction: rotate by cos_t, phi about old direction
            sin_t = np.sqrt(max(0.0, 1.0 - cos_t**2))
            local = np.array([sin_t*np.cos(phi), sin_t*np.sin(phi), cos_t])
            direction = sampling._rotate_to_global(local, direction)
            E = Egp
            # continue the while loop (this same gamma keeps going)

        elif proc == "pair":
            eps = sampling.sample_pair(E, rng)
            E_pos_total = eps * E
            E_ele_total = (1.0 - eps) * E
            T_pos = E_pos_total - ME_MEV
            T_ele = E_ele_total - ME_MEV
            # Sample directions (forward-peaked)
            theta_pos = sampling.pair_polar_angle(E_pos_total, rng)
            theta_ele = sampling.pair_polar_angle(E_ele_total, rng)
            phi = 2*np.pi*rng.uniform()
            # Place leptons coplanar with incoming photon
            for sign, theta, T, kind in [(1.0, theta_pos, T_pos, "positron"),
                                          (-1.0, theta_ele, T_ele, "electron")]:
                local = np.array([np.sin(theta)*np.cos(phi + (np.pi if sign < 0 else 0)),
                                   np.sin(theta)*np.sin(phi + (np.pi if sign < 0 else 0)),
                                   np.cos(theta)])
                d_lep = sampling._rotate_to_global(local, direction)
                if T > 0:
                    stack.push(Track(kind=kind, energy=T, position=pos.copy(),
                                     direction=d_lep, parent_id=id(track),
                                     generation=track.generation + 1))
            return  # photon consumed

        elif proc == "photoelectric":
            shell, E_bind = sampling.sample_phot_shell(E, rng)
            T_e = E - E_bind
            theta_e = sampling.sample_photoelectron_angle(T_e, rng)
            phi_e = 2*np.pi*rng.uniform()
            local = np.array([np.sin(theta_e)*np.cos(phi_e),
                               np.sin(theta_e)*np.sin(phi_e),
                               np.cos(theta_e)])
            d_e = sampling._rotate_to_global(local, direction)
            if T_e > 0:
                stack.push(Track(kind="electron", energy=T_e,
                                 position=pos.copy(), direction=d_e,
                                 parent_id=id(track),
                                 generation=track.generation + 1))
            # Atomic relaxation
            if shell == "K":
                outcome, E_relax = sampling.sample_atomic_relaxation_K(rng)
                if outcome == "fluorescence":
                    # Push a new photon, isotropic
                    cos_t = -1 + 2*rng.uniform()
                    phi = 2*np.pi*rng.uniform()
                    sin_t = np.sqrt(1 - cos_t**2)
                    d_g = np.array([sin_t*np.cos(phi), sin_t*np.sin(phi), cos_t])
                    stack.push(Track(kind="gamma", energy=E_relax,
                                     position=pos.copy(), direction=d_g,
                                     parent_id=id(track),
                                     generation=track.generation + 1))
                else:
                    # Auger electron deposits locally
                    deposits.append(Deposit(pos.copy(), E_relax, "auger"))
            else:
                # L/M binding: deposit locally
                deposits.append(Deposit(pos.copy(), E_bind, "auger"))
            return  # photon consumed

    # Below threshold: deposit residual energy locally
    deposits.append(Deposit(pos, E, "gamma_local"))


# =====================================================================
# Lepton transport (electron and positron)
# =====================================================================

def transport_lepton(track: Track, geom, deposits: List[Deposit],
                     stack: Stack, rng: np.random.Generator) -> None:
    """Step the electron/positron through LXe with continuous energy loss
    plus discrete bremsstrahlung. At end of range, positron annihilates."""
    pos = track.position.copy()
    direction = track.direction.copy()
    T = track.energy

    while T >= th.T_E_CUT_MEV:
        # determine step size from residual range
        R_cm = csda_range_g_per_cm2_NIST(T) / RHO_LXE  # cm
        ds = th.F_RANGE * R_cm
        ds = max(th.DS_FLOOR_CM, min(ds, th.DS_CEIL_CM))

        # collisional energy loss (continuous)
        dEdx_col = dEdx_collision_NIST(T) * RHO_LXE   # MeV/cm
        dE_col = dEdx_col * ds
        if dE_col >= T:
            # would over-shoot; reduce step
            ds = T / dEdx_col * 0.9
            dE_col = dEdx_col * ds

        # deposit collisional energy at midpoint of step
        mid_pos = pos + direction * (ds * 0.5)
        deposits.append(Deposit(mid_pos.copy(), dE_col, track.kind))
        T -= dE_col

        # advance position
        pos = pos + direction * ds
        if not geom.is_inside(pos):
            # exit detector; deposit nothing more
            return

        if T < th.T_E_CUT_MEV:
            break

        # Sample bremsstrahlung
        if T > th.K_MIN_MEV:
            sig_b = sampling._sigma_brems_above_kmin(T, th.K_MIN_MEV)
            P_brems = NAT_LXE * sig_b * ds
            if P_brems > 0.5:
                P_brems = 0.5  # safety cap; should not happen with our step sizes
            if rng.uniform() < P_brems:
                k = sampling.sample_brems(T, th.K_MIN_MEV, rng)
                if k is not None and k < T:
                    theta_g = sampling.brems_photon_angle(T, rng)
                    phi_g = 2*np.pi*rng.uniform()
                    local = np.array([np.sin(theta_g)*np.cos(phi_g),
                                       np.sin(theta_g)*np.sin(phi_g),
                                       np.cos(theta_g)])
                    d_g = sampling._rotate_to_global(local, direction)
                    stack.push(Track(kind="gamma", energy=k,
                                     position=pos.copy(), direction=d_g,
                                     parent_id=id(track),
                                     generation=track.generation + 1))
                    T -= k

    # End-of-range: deposit remaining T
    if T > 0 and geom.is_inside(pos):
        deposits.append(Deposit(pos.copy(), T, track.kind))

    # Positron: annihilation at rest
    if track.kind == "positron":
        # two 511 keV photons back-to-back, isotropic
        cos_t = -1 + 2*rng.uniform()
        phi = 2*np.pi*rng.uniform()
        sin_t = np.sqrt(1 - cos_t**2)
        d_g = np.array([sin_t*np.cos(phi), sin_t*np.sin(phi), cos_t])
        for d in [d_g, -d_g]:
            stack.push(Track(kind="gamma", energy=ME_MEV,
                             position=pos.copy(), direction=d,
                             parent_id=id(track),
                             generation=track.generation + 1))


# =====================================================================
# Top-level event simulation
# =====================================================================

def simulate_event(E_MeV: float, position=(0.0, 0.0, 0.0),
                   direction=(0.0, 0.0, 1.0),
                   geom=None, rng=None) -> List[Deposit]:
    """Simulate one primary photon entering the calorimeter.

    Returns a list of energy Deposit objects.
    """
    if geom is None:
        geom = InfiniteLXe()
    if rng is None:
        rng = np.random.default_rng()

    deposits: List[Deposit] = []
    stack = Stack()
    stack.push(Track(kind="gamma", energy=E_MeV,
                     position=np.asarray(position, dtype=float),
                     direction=np.asarray(direction, dtype=float),
                     parent_id=-1, generation=0))

    while not stack.empty():
        t = stack.pop()
        if t.generation > th.GEN_CAP:
            continue   # safety
        if t.kind == "gamma":
            transport_photon(t, geom, deposits, stack, rng)
        else:
            transport_lepton(t, geom, deposits, stack, rng)

    return deposits


def cluster_deposits_in_z(deposits: List[Deposit],
                           dz_cm: float = th.DZ_RESOLUTION_CM
                           ) -> List[Tuple[float, float]]:
    """Group deposits into clusters whose z-extent is within dz_cm.

    Returns list of (z_centroid, total_energy). Used for SS/MS classification.
    """
    if not deposits:
        return []
    items = sorted(deposits, key=lambda d: d.position[2])
    clusters = []
    cur_zs = [items[0].position[2]]
    cur_es = [items[0].energy]
    for d in items[1:]:
        if d.position[2] - cur_zs[-1] < dz_cm:
            cur_zs.append(d.position[2])
            cur_es.append(d.energy)
        else:
            E_tot = sum(cur_es)
            z_cent = sum(z*e for z, e in zip(cur_zs, cur_es)) / E_tot
            clusters.append((z_cent, E_tot))
            cur_zs = [d.position[2]]
            cur_es = [d.energy]
    E_tot = sum(cur_es)
    z_cent = sum(z*e for z, e in zip(cur_zs, cur_es)) / E_tot
    clusters.append((z_cent, E_tot))
    return clusters


def is_single_site(deposits: List[Deposit],
                    dz_cm: float = th.DZ_RESOLUTION_CM) -> bool:
    return len(cluster_deposits_in_z(deposits, dz_cm)) == 1


if __name__ == "__main__":
    rng = np.random.default_rng(42)
    print("Simulating 100 events of 2.615 MeV gammas in infinite LXe...")
    n_ss = 0
    total_E = []
    for i in range(100):
        deposits = simulate_event(2.615, rng=rng)
        E_dep = sum(d.energy for d in deposits)
        total_E.append(E_dep)
        if is_single_site(deposits):
            n_ss += 1
    print(f"  Mean total deposit: {np.mean(total_E):.4f} MeV")
    print(f"  std total deposit: {np.std(total_E):.4f} MeV")
    print(f"  Range: {np.min(total_E):.4f} to {np.max(total_E):.4f}")
    print(f"  SS fraction (3mm in z, infinite LXe): {n_ss}/100 = {n_ss/100*100:.1f}%")
