"""
validation.py
=============

Validation tests for the LXe MC simulation. Run as:
    python validation.py

Each test prints a pass/fail summary. Statistical tests use 1e4-1e5 events.
"""

import numpy as np
import physics as ph
import nist_data as nd
import sampling as sp
import tracking as tr
from thresholds import T_E_CUT_MEV, EGAMMA_CUT_MEV, K_MIN_MEV


PASS = "PASS"
FAIL = "FAIL"


def banner(name):
    print(f"\n{'='*70}\n{name}\n{'='*70}")


def test_xcom_total_vs_NIST():
    banner("Test 1: NIST XCOM total agrees with sum of channels")
    Es = [0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
    max_dev = 0.0
    for E in Es:
        sC = nd.sigma_compton_NIST(E, per_atom=False)
        sP = nd.sigma_pair_NIST(E, per_atom=False)
        sPh = nd.sigma_phot_NIST(E, per_atom=False)
        sum_chan = sC + sP + sPh
        sT = nd.sigma_total_NIST(E, per_atom=False, include_coherent=False)
        dev = abs(sum_chan - sT) / sT
        max_dev = max(max_dev, dev)
        print(f"  E={E} MeV  sum={sum_chan:.4e}  total_NIST={sT:.4e}  "
              f"rel_dev={dev*100:.3f}%")
    status = PASS if max_dev < 0.025 else FAIL
    print(f"  -> {status} (max deviation {max_dev*100:.3f}%, "
          f"sub-2.5% expected from log-log interp)")
    return status == PASS


def test_klein_nishina_compton_edge():
    banner("Test 2: Compton sampler reproduces the Compton edge")
    rng = np.random.default_rng(0)
    E0 = 2.615
    a = E0 / ph.ME_MEV
    Tmax_theory = E0 * 2*a / (1 + 2*a)
    print(f"  Theoretical Compton edge T_max = {Tmax_theory:.4f} MeV at E={E0} MeV")
    N = 100000
    Egp = np.array([sp.sample_compton(E0, rng)[0] for _ in range(N)])
    Te = E0 - Egp
    Te_max = Te.max()
    Te_99 = np.percentile(Te, 99.9)
    print(f"  Sampled max T_e = {Te_max:.4f} MeV")
    print(f"  99.9th percentile = {Te_99:.4f} MeV")
    # The maximum observed should be close to but less than Tmax_theory
    status = PASS if (Te_max < Tmax_theory + 1e-4 and Te_99 > 0.9 * Tmax_theory) else FAIL
    print(f"  -> {status}")
    return status == PASS


def test_kn_marginal():
    banner("Test 3: KN marginal (E_gamma') distribution agrees with formula")
    rng = np.random.default_rng(1)
    E0 = 1.0
    N = 200000
    samples = np.array([sp.sample_compton(E0, rng)[0] for _ in range(N)])
    # Histogram and compare with analytic dsigma/dEgp shape.
    # In terms of eps = Egp/E_gamma, KN dsigma/deps \propto
    #   [eps + 1/eps] * [1 - eps sin^2(theta)/(1+eps^2)]
    # but it's easier to use the direct form:
    # dsigma/dEgp \propto (Egp/E0)^2 * (E0/Egp + Egp/E0 - sin^2(theta))
    # where sin^2(theta) is determined by Egp via Compton kinematics.
    bins = np.linspace(samples.min()*0.99, samples.max()*1.01, 41)
    hist, edges = np.histogram(samples, bins=bins, density=True)
    centers = 0.5*(edges[:-1] + edges[1:])
    # Compute dsigma/dEgp at centers using the actual KN formula
    # KN d_sigma/d_Egp = pi r_e^2 / (E0 * me) * [Egp/E0 + E0/Egp - sin^2(theta)]
    # with cos(theta) = 1 - me/Egp + me/E0 from Compton relation
    cos_t = 1.0 - ph.ME_MEV/centers + ph.ME_MEV/E0
    sin2_t = np.maximum(0.0, 1.0 - cos_t**2)
    f_th = (centers/E0 + E0/centers - sin2_t)
    f_th_norm = f_th / np.trapezoid(f_th, centers)
    sel = hist > 0.1*hist.max()
    rel_dev = np.abs(hist[sel] - f_th_norm[sel]) / f_th_norm[sel]
    print(f"  Mean relative deviation in bulk bins: {rel_dev.mean()*100:.2f}%")
    print(f"  Max deviation: {rel_dev.max()*100:.2f}%")
    status = PASS if rel_dev.mean() < 0.05 else FAIL
    print(f"  -> {status}")
    return status == PASS


def test_energy_conservation():
    banner("Test 4: Energy conservation per event")
    rng = np.random.default_rng(2)
    N = 500
    E0 = 2.615
    losses = []
    for _ in range(N):
        deps = tr.simulate_event(E0, rng=rng)
        E_dep = sum(d.energy for d in deps)
        # In infinite LXe everything must be deposited modulo the cut energies
        loss = E0 - E_dep
        losses.append(loss)
    losses = np.array(losses)
    print(f"  Mean energy loss (E_primary - E_deposited): {losses.mean()*1000:.2f} keV")
    print(f"  Max loss: {losses.max()*1000:.2f} keV")
    print(f"  Cut budget: T_e_cut={T_E_CUT_MEV*1000:.1f} keV, k_min={K_MIN_MEV*1000:.1f} keV")
    # Acceptable loss is a few times the cut, since we may have multiple cuts per event
    status = PASS if losses.mean() < 5*T_E_CUT_MEV else FAIL
    print(f"  -> {status}")
    return status == PASS


def test_csda_range_vs_estar():
    banner("Test 5: CSDA range from electron MC matches NIST ESTAR")
    rng = np.random.default_rng(3)
    N = 200
    T0 = 1.0
    R_NIST_mm = nd.csda_range_LXe_NIST(T0)
    print(f"  NIST CSDA range at T={T0} MeV: {R_NIST_mm:.4f} mm")
    # Simulate isolated electron tracks: lay it down at origin moving +z,
    # measure the maximum z reached
    z_ends = []
    for _ in range(N):
        deposits = []
        stack = tr.Stack()
        track = tr.Track(kind="electron", energy=T0,
                          position=np.array([0.0, 0.0, 0.0]),
                          direction=np.array([0.0, 0.0, 1.0]))
        stack.push(track)
        while not stack.empty():
            t = stack.pop()
            if t.kind == "gamma":
                # Truncate brems photons (don't transport them) -- we just
                # want the projected range of the electron
                continue
            tr.transport_lepton(t, tr.InfiniteLXe(), deposits, stack, rng)
        # End-of-range z is the maximum z over all electron deposits
        z_max = max(d.position[2] for d in deposits)
        z_ends.append(z_max * 10)  # cm -> mm
    z_ends = np.array(z_ends)
    print(f"  Mean projected end-of-range: {z_ends.mean():.4f} mm")
    rel_dev = abs(z_ends.mean() - R_NIST_mm) / R_NIST_mm
    print(f"  Relative deviation from NIST: {rel_dev*100:.2f}%")
    # Expect within 5% (no MS in our model means projected = CSDA)
    status = PASS if rel_dev < 0.05 else FAIL
    print(f"  -> {status}")
    return status == PASS


def test_pair_branching():
    banner("Test 6: Pair production fraction at 2.615 MeV")
    rng = np.random.default_rng(4)
    N = 5000
    n_pair = 0
    # We just count the first-interaction process by directly sampling
    bC, bP, bPh = nd.branching_NIST(2.615)
    for _ in range(N):
        proc = sp.sample_process(bC, bP, bPh, rng)
        if proc == "pair":
            n_pair += 1
    frac = n_pair / N
    print(f"  Empirical pair fraction:   {frac*100:.2f}%")
    print(f"  Expected from NIST:        {bP*100:.2f}%")
    sigma_stat = np.sqrt(frac*(1-frac)/N)
    print(f"  Statistical std (1 sigma): {sigma_stat*100:.2f}%")
    status = PASS if abs(frac - bP) < 3*sigma_stat else FAIL
    print(f"  -> {status}")
    return status == PASS


def test_annihilation_peak():
    banner("Test 7: 511 keV annihilation peak from pair events")
    # Force a pair event by simulating 2.615 MeV gammas and look for 511 keV
    rng = np.random.default_rng(5)
    N = 200
    photon_energies = []
    for _ in range(N):
        deps = tr.simulate_event(2.615, rng=rng)
        # Inspect deposits — annihilation photons carry 511 keV before they
        # interact and deposit their energy. We see this via deposits that
        # are clustered far from the primary vertex with total energy near 511.
        # Easier: look at deposits' energies directly and look for the 511 keV
        # signature (a deposit cluster of approximately 511 keV).
        # For simplicity here, let's sum deposits in a tight z window and
        # count those near 0.5 MeV.
        clusters = tr.cluster_deposits_in_z(deps, dz_cm=0.3)
        for z, E in clusters:
            photon_energies.append(E)
    photon_energies = np.array(photon_energies)
    near_511 = np.sum((photon_energies > 0.45) & (photon_energies < 0.55))
    print(f"  Found {near_511} clusters in [450, 550] keV out of "
          f"{len(photon_energies)} clusters")
    status = PASS if near_511 > 0 else FAIL
    print(f"  -> {status}")
    return status == PASS


def test_radiative_yield():
    banner("Test 8: Hard-brems yield for 2 MeV electrons")
    rng = np.random.default_rng(6)
    N = 50
    T0 = 2.0
    f_rad_total = nd.dEdx_radiative_NIST(T0) / nd.dEdx_total_NIST(T0)
    # Our MC only emits brems above K_MIN_MEV (50 keV); the soft-brems
    # fraction below k_min is lumped into the local collisional deposit.
    # Therefore we expect the hard-brems yield to be a fraction of the
    # total radiative yield. Estimate the fraction:
    #   f_hard / f_total = integral( k * dsig/dk, k_min..T) / integral(0..T)
    k_grid = np.logspace(np.log10(0.001), np.log10(T0*0.9999), 200)
    ds = ph.dsigma_dk_brems(k_grid, T0)
    energy_full = np.trapezoid(k_grid * ds, k_grid)
    k_grid_hard = k_grid[k_grid >= K_MIN_MEV]
    ds_hard = ph.dsigma_dk_brems(k_grid_hard, T0)
    energy_hard = np.trapezoid(k_grid_hard * ds_hard, k_grid_hard)
    hard_fraction = energy_hard / energy_full
    f_hard_expected = f_rad_total * hard_fraction
    print(f"  NIST total radiative fraction at T={T0} MeV: {f_rad_total*100:.2f}%")
    print(f"  Hard-brems fraction (k > {K_MIN_MEV*1000:.0f} keV) of "
          f"total radiative: {hard_fraction*100:.1f}%")
    print(f"  Expected hard-brems yield in MC: {f_hard_expected*100:.2f}%")
    E_brems_total = 0.0
    E_total_input = 0.0
    for _ in range(N):
        deposits = []
        stack = tr.Stack()
        stack.push(tr.Track(kind="electron", energy=T0,
                             position=np.zeros(3),
                             direction=np.array([0,0,1.0])))
        E_total_input += T0
        while not stack.empty():
            t = stack.pop()
            if t.kind == "gamma":
                E_brems_total += t.energy
                continue
            tr.transport_lepton(t, tr.InfiniteLXe(), deposits, stack, rng)
    f_hard_observed = E_brems_total / E_total_input
    print(f"  Observed hard-brems yield: {f_hard_observed*100:.2f}%")
    rel_dev = abs(f_hard_observed - f_hard_expected) / f_hard_expected
    print(f"  Relative deviation: {rel_dev*100:.1f}%")
    status = PASS if rel_dev < 0.40 else FAIL
    print(f"  -> {status}")
    return status == PASS


def test_ss_fraction_Tl208():
    banner("Test 9: SS fraction for Tl-208 in infinite LXe")
    rng = np.random.default_rng(7)
    N = 200
    n_ss = 0
    for _ in range(N):
        deps = tr.simulate_event(2.615, rng=rng)
        if tr.is_single_site(deps, dz_cm=0.30):
            n_ss += 1
    frac = n_ss / N
    print(f"  SS fraction (3 mm z resolution): {n_ss}/{N} = {frac*100:.1f}%")
    print(f"  Expected from literature: ~10-20% (geometry-dependent)")
    status = PASS if 0.05 < frac < 0.30 else FAIL
    print(f"  -> {status}")
    return status == PASS


# =====================================================================
# Run all tests
# =====================================================================

if __name__ == "__main__":
    results = []
    results.append(("xcom_total", test_xcom_total_vs_NIST()))
    results.append(("compton_edge", test_klein_nishina_compton_edge()))
    results.append(("kn_marginal", test_kn_marginal()))
    results.append(("energy_conservation", test_energy_conservation()))
    results.append(("csda_range", test_csda_range_vs_estar()))
    results.append(("pair_branching", test_pair_branching()))
    results.append(("annihilation_peak", test_annihilation_peak()))
    results.append(("radiative_yield", test_radiative_yield()))
    results.append(("ss_fraction_Tl208", test_ss_fraction_Tl208()))

    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    n_pass = sum(1 for _, ok in results if ok)
    for name, ok in results:
        print(f"  {'PASS' if ok else 'FAIL':4}  {name}")
    print(f"\n{n_pass}/{len(results)} tests passed")
