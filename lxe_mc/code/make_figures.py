"""
make_figures.py
===============

Generate all figures and data tables for the LXe MC manual.

Lines = analytic formulas (from physics.py)
Markers = NIST data (XCOM, ESTAR; from nist_data.py)

Outputs go into ../manual/figs/ (PDF + PNG) and ../data/.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import physics as ph
import nist_data as nd

mpl.rcParams.update({
    "font.size": 11, "axes.grid": True, "grid.alpha": 0.3,
    "lines.linewidth": 1.6, "figure.dpi": 110, "savefig.dpi": 200,
    "savefig.bbox": "tight",
})

_HERE = os.path.dirname(os.path.abspath(__file__))
FIGS_DIR = os.path.normpath(os.path.join(_HERE, "..", "manual", "figs"))
DATA_DIR = os.path.normpath(os.path.join(_HERE, "..", "data"))
os.makedirs(FIGS_DIR, exist_ok=True)
os.makedirs(DATA_DIR, exist_ok=True)


def savefig(fig, name):
    pdf = os.path.join(FIGS_DIR, f"{name}.pdf")
    png = os.path.join(FIGS_DIR, f"{name}.png")
    fig.savefig(pdf)
    fig.savefig(png)
    plt.close(fig)
    print(f"  saved {name}.pdf and {name}.png")


# Energies of interest
Q_BB     = 2.4578
E_BI214  = 2.4476
E_TL208  = 2.6145


# =====================================================================
# Figure 1: photon cross sections in Xe (lines = formulas, points = NIST)
# =====================================================================
def fig_xsec():
    fig, ax = plt.subplots(figsize=(8, 5.5))
    E = np.logspace(np.log10(0.01), np.log10(20.0), 600)

    # Lines: analytic formulas (mass-attenuation)
    mu_C   = ph.sigma_compton(E)        * ph.NA / ph.A_XE
    # for pair, use NIST since BH numerical underestimates near threshold
    # but plot a thin line of analytic for completeness
    mu_P_NIST = nd.sigma_pair_NIST(E, per_atom=False)
    mu_Ph_S = ph.sigma_phot_sandia(E, per_atom=False)  # cm^2/g

    ax.loglog(E, mu_C, "-", color="C0", label="Compton (KN $\\times Z$)", lw=1.5)
    # Plot Sandia only over its valid domain solid, dashed elsewhere
    valid = E < 0.5
    invalid = E >= 0.5
    ax.loglog(E[valid], mu_Ph_S[valid], "-", color="C1",
              label="Photoelectric (Sandia, $E<500$ keV)", lw=1.5)
    ax.loglog(E[invalid], mu_Ph_S[invalid], ":", color="C1", lw=1.0,
              alpha=0.7, label="Sandia extrapolation (invalid)")
    ax.loglog(E[E > 2*ph.ME_MEV], mu_P_NIST[E > 2*ph.ME_MEV], "-", color="C2",
              label="Pair (NIST XCOM)", lw=1.5)

    # NIST points for each channel
    Ex = nd._XCOM["E_MeV"]
    ax.loglog(Ex, nd._XCOM["incoherent"], "o", color="C0", ms=4, mfc="white",
              label="NIST Compton")
    ax.loglog(Ex, nd._XCOM["photoelectric"], "s", color="C1", ms=4, mfc="white",
              label="NIST Photoelectric")
    pair_total_NIST = nd._XCOM["pair_nuclear"] + nd._XCOM["pair_electron"]
    pos = pair_total_NIST > 0
    ax.loglog(Ex[pos], pair_total_NIST[pos], "^", color="C2", ms=5, mfc="white",
              label="NIST Pair")
    ax.loglog(Ex, nd._XCOM["total_no_coh"], "kD", ms=3.5, mfc="black",
              label="NIST total")

    ax.axvspan(E_BI214, E_TL208, color="red", alpha=0.10)
    ax.axvline(Q_BB, color="red", ls=":", alpha=0.7)

    ax.set_xlabel("photon energy $E_\\gamma$ [MeV]")
    ax.set_ylabel("mass attenuation $\\mu/\\rho$ [cm$^2$/g]")
    ax.set_title("Photon cross sections in xenon ($Z=54$)")
    ax.legend(framealpha=0.95, fontsize=8.5, loc="best")
    ax.set_xlim(0.01, 20)
    ax.set_ylim(1e-4, 1e4)
    fig.tight_layout()
    savefig(fig, "xsec_xe")


# =====================================================================
# Figure 1b: Sandia validity — analytic vs NIST in the photoelectric channel
# =====================================================================
def fig_sandia_vs_nist():
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.6))
    E = np.logspace(np.log10(0.05), np.log10(10.0), 400)

    # Left: cross section
    ax = axs[0]
    sandia = ph.sigma_phot_sandia(E, per_atom=False)
    nist = nd.sigma_phot_NIST(E, per_atom=False)
    ax.loglog(E, sandia, "-", color="C1", label="Sandia (Biggs-Lighthill)", lw=1.5)
    ax.loglog(E, nist, "k--", label="NIST XCOM", lw=1.5)
    ax.axvline(0.5, color="gray", ls=":", alpha=0.7,
               label="$E=500$ keV (Sandia validity bound)")
    ax.axvspan(E_BI214, E_TL208, color="red", alpha=0.10, label="ROI")
    ax.set_xlabel("photon energy $E_\\gamma$ [MeV]")
    ax.set_ylabel("photoelectric $\\mu/\\rho$ [cm$^2$/g]")
    ax.set_title("Sandia vs. NIST photoelectric in Xe")
    ax.legend(framealpha=0.95, fontsize=9)

    # Right: ratio
    ax = axs[1]
    nist_arr = np.array([nd.sigma_phot_NIST(e, per_atom=False) for e in E])
    ratio = sandia / nist_arr
    ax.semilogx(E, ratio, "-", color="C1", lw=1.5)
    ax.axhline(1.0, color="k", ls="--", alpha=0.7)
    ax.axvline(0.5, color="gray", ls=":", alpha=0.7)
    ax.axvspan(E_BI214, E_TL208, color="red", alpha=0.10)
    ax.set_xlabel("$E_\\gamma$ [MeV]")
    ax.set_ylabel("Sandia / NIST")
    ax.set_title("Ratio Sandia / NIST (Xe)")
    ax.set_ylim(0, 5)
    fig.tight_layout()
    savefig(fig, "sandia_vs_nist")


# =====================================================================
# Figure 2: branching fractions (NIST)
# =====================================================================
def fig_branching():
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    E = np.linspace(0.5, 5.0, 300)
    sC = nd.sigma_compton_NIST(E, per_atom=False)
    sP = nd.sigma_pair_NIST(E, per_atom=False)
    sPh = nd.sigma_phot_NIST(E, per_atom=False)
    st = sC + sP + sPh
    ax.plot(E, sC/st, "-", color="C0", label="Compton", lw=1.7)
    ax.plot(E, sP/st, "-", color="C2", label="Pair", lw=1.7)
    ax.plot(E, sPh/st, "-", color="C1", label="Photoelectric", lw=1.7)
    ax.axvspan(E_BI214, E_TL208, color="red", alpha=0.10, label="ROI")
    ax.axvline(Q_BB, color="red", ls=":", alpha=0.7,
               label=f"$Q_{{\\beta\\beta}}={Q_BB:.3f}$ MeV")
    ax.set_xlabel("photon energy $E_\\gamma$ [MeV]")
    ax.set_ylabel("branching fraction")
    ax.set_title("Process branching in Xe (NIST XCOM)")
    ax.legend(framealpha=0.95, loc="center right")
    ax.set_ylim(0, 1)
    fig.tight_layout()
    savefig(fig, "branching_xe")


# =====================================================================
# Figure 3: Klein-Nishina differential cross section
# =====================================================================
def fig_kn_diff():
    fig, ax = plt.subplots(figsize=(7.5, 5))
    theta = np.linspace(0.001, np.pi - 0.001, 400)
    for E_MeV, c, lbl in [(0.1, "C0", "0.1 MeV"),
                           (0.5, "C1", "0.5 MeV"),
                           (1.0, "C2", "1.0 MeV"),
                           (2.4476, "C3", f"$E_{{Bi}}={E_BI214}$ MeV"),
                           (2.6145, "C4", f"$E_{{Tl}}={E_TL208}$ MeV")]:
        ds = ph.dsigma_dOmega_KN(theta, E_MeV) / ph.RE_CM**2
        ax.plot(np.degrees(theta), ds, color=c, label=lbl)
    ax.set_xlabel("photon scattering angle $\\theta$ [deg]")
    ax.set_ylabel("$d\\sigma/d\\Omega$ [$r_e^2$/sr]")
    ax.set_title("Klein-Nishina differential cross section")
    ax.set_xlim(0, 180)
    ax.set_yscale("log")
    ax.legend(framealpha=0.95)
    fig.tight_layout()
    savefig(fig, "kn_diff")


# =====================================================================
# Figure 4: Compton electron spectrum
# =====================================================================
def fig_compton_electron_spectrum():
    fig, ax = plt.subplots(figsize=(7.5, 5))
    me = ph.ME_MEV
    for E_MeV, c, lbl in [(0.5, "C1", "0.5 MeV"),
                           (1.0, "C2", "1.0 MeV"),
                           (2.4476, "C3", f"{E_BI214} MeV (Bi)"),
                           (2.6145, "C4", f"{E_TL208} MeV (Tl)")]:
        a = E_MeV / me
        Tmax = E_MeV * 2*a / (1 + 2*a)
        T = np.linspace(0.001, 0.9999*Tmax, 400)
        s = T / E_MeV
        # shape proportional to KN d sigma/dT
        y = 2.0 + s**2/(a**2*(1-s)**2) - 2.0*s/(a*(1-s))
        y = y / np.max(y)
        ax.plot(T, y, color=c, label=f"{lbl}, $T_e^{{\\max}}={Tmax:.3f}$")
        ax.axvline(Tmax, ls="--", color=c, alpha=0.4)
    ax.set_xlabel("Compton electron kinetic energy $T_e$ [MeV]")
    ax.set_ylabel("$dN/dT_e$ (peak normalized)")
    ax.set_title("Compton electron energy spectrum")
    ax.legend(loc="upper left", framealpha=0.95)
    fig.tight_layout()
    savefig(fig, "compton_electron_spectrum")


# =====================================================================
# Figure 5: stopping power and CSDA range (formulas vs NIST)
# =====================================================================
def fig_stopping_range():
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.6))

    # Stopping power
    ax = axs[0]
    T = np.logspace(np.log10(0.01), np.log10(10.0), 200)
    Sc_BS = ph.dEdx_collision(T)
    Sr_BH = ph.dEdx_radiative(T)
    Sc_NIST = nd.dEdx_collision_NIST(T)
    Sr_NIST = nd.dEdx_radiative_NIST(T)

    ax.loglog(T, Sc_BS, "-", color="C0", label="Berger-Seltzer (formula)", lw=1.4)
    ax.loglog(T, Sr_BH, "-", color="C2", label="Tsai (formula)", lw=1.4)
    ax.loglog(T, Sc_NIST, "kD", ms=3, mfc="white", label="NIST collisional")
    ax.loglog(T, Sr_NIST, "k^", ms=3, mfc="white", label="NIST radiative")
    ax.set_xlabel("electron kinetic energy $T$ [MeV]")
    ax.set_ylabel("$S/\\rho$ [MeV cm$^2$/g]")
    ax.set_title("Electron stopping power in Xe")
    ax.legend(framealpha=0.95, fontsize=9)

    # Range
    ax = axs[1]
    R_NIST = nd.csda_range_LXe_NIST(T)
    ax.loglog(T, R_NIST, "-", color="C3", lw=1.7, label="CSDA range (NIST)")
    ax.axhline(0.3, color="k", ls="--", alpha=0.5,
               label="$\\Delta z = 3$ mm")
    for Tval in [0.1, 0.5, 1.0, 2.0, 2.4]:
        Rv = nd.csda_range_LXe_NIST(Tval)
        ax.plot(Tval, Rv, "o", color="C3", ms=5)
        ax.annotate(f"{Rv:.2f} mm", (Tval, Rv), xytext=(8, -10),
                    textcoords="offset points", fontsize=9)
    ax.set_xlabel("electron kinetic energy $T$ [MeV]")
    ax.set_ylabel("CSDA range in LXe [mm]")
    ax.set_title("Electron CSDA range in LXe")
    ax.legend(framealpha=0.95)
    ax.set_xlim(0.01, 10)
    ax.set_ylim(1e-3, 100)
    fig.tight_layout()
    savefig(fig, "stopping_range_xe")


# =====================================================================
# Figure 6: bremsstrahlung spectrum
# =====================================================================
def fig_brems():
    fig, ax = plt.subplots(figsize=(7.5, 5))
    for T_e, c, lbl in [(0.5, "C1", "$T = 0.5$ MeV"),
                         (1.0, "C2", "$T = 1.0$ MeV"),
                         (2.0, "C3", "$T = 2.0$ MeV")]:
        k = np.linspace(0.005, 0.999*T_e, 200)
        ds = ph.dsigma_dk_brems(k, T_e)
        ax.plot(k/T_e, k*ds * 1e27, color=c, label=lbl, lw=1.5)
    ax.set_xlabel("photon energy fraction $k/T$")
    ax.set_ylabel("$k\\,d\\sigma/dk$ [mb] (per atom)")
    ax.set_title("Bremsstrahlung spectrum in Xe (Bethe-Heitler-Tsai)")
    ax.legend(framealpha=0.95)
    ax.set_xscale("log")
    ax.set_xlim(0.005, 1.0)
    fig.tight_layout()
    savefig(fig, "brems_spectrum")


# =====================================================================
# Figure 7: SS-faking channels (new pedagogical figure)
# =====================================================================
def fig_ss_faking():
    """Illustrate which channels can fake an SS event in the ROI from
    Bi-214 (2.448 MeV, essentially in ROI) and Tl-208 (2.615 MeV)."""
    fig, ax = plt.subplots(figsize=(8.5, 5))
    sources = ["$^{214}$Bi (2.448 MeV)", "$^{208}$Tl (2.615 MeV)"]
    Es = [E_BI214, E_TL208]

    # For each source, compute branching
    bC_arr = []
    bP_arr = []
    bPh_arr = []
    for E in Es:
        bC, bP, bPh = nd.branching_NIST(E)
        bC_arr.append(bC)
        bP_arr.append(bP)
        bPh_arr.append(bPh)

    x = np.arange(len(sources))
    width = 0.6
    p1 = ax.bar(x, bC_arr, width, color="C0", label="Compton")
    p2 = ax.bar(x, bP_arr, width, bottom=bC_arr, color="C2", label="Pair")
    p3 = ax.bar(x, bPh_arr, width,
                bottom=np.array(bC_arr) + np.array(bP_arr),
                color="C1", label="Photoelectric")

    ax.set_xticks(x)
    ax.set_xticklabels(sources, fontsize=11)
    ax.set_ylabel("first-interaction branching fraction (NIST)")
    ax.set_title("First-interaction branching for ROI gammas in Xe")
    ax.legend(framealpha=0.95)
    ax.set_ylim(0, 1.05)

    # Annotate the photoelectric fraction (the SS-mimic mode)
    for i, b in enumerate(bPh_arr):
        ax.text(i, 0.99, f"Ph: {b*100:.2f}%\n(direct full-E SS)",
                ha="center", va="top", color="black", fontsize=9,
                bbox=dict(facecolor="white", alpha=0.8, pad=2))

    fig.tight_layout()
    savefig(fig, "ss_faking_channels")


# =====================================================================
# Figure 8: topology schematic (no data)
# =====================================================================
def fig_topology():
    fig, axs = plt.subplots(1, 2, figsize=(11, 5))

    def draw(ax, with_e, title):
        ax.annotate("", xy=(2, 5.5), xytext=(2, 7.5),
                    arrowprops=dict(arrowstyle="->", color="gray", lw=1.6))
        verts = [(2, 5.5), (4, 3.0), (5.5, 1.0)]
        for v in verts:
            ax.plot(*v, "ks", ms=8)
        ax.annotate("", xy=verts[1], xytext=verts[0],
                    arrowprops=dict(arrowstyle="->", color="gray", lw=1.6))
        ax.annotate("", xy=verts[2], xytext=verts[1],
                    arrowprops=dict(arrowstyle="->", color="gray", lw=1.6))
        if not with_e:
            for x, y in verts:
                ax.plot(x, y, "o", ms=18, color="C3", alpha=0.5)
        else:
            ax.plot([2, 2.4, 2.7, 2.95, 3.0], [5.5, 5.3, 4.8, 4.3, 3.7],
                    "r-", lw=2.5)
            ax.plot([4, 4.1, 4.0], [3.0, 2.7, 2.5], "r-", lw=2.5)
            ax.plot([5.5, 5.7], [1.0, 0.9], "r-", lw=2.5)
            ax.annotate("", xy=(4.5, 5.1), xytext=(2.7, 5.2),
                        arrowprops=dict(arrowstyle="->", color="C2", lw=1.4))
            ax.plot(4.7, 5.1, "o", ms=10, color="C2", alpha=0.5)
            ax.text(4.85, 5.1, "brems\n$\\gamma$ deposit",
                    color="C2", fontsize=9)
        ax.set_xlim(0, 8)
        ax.set_ylim(0, 8)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(title, fontsize=12)
        ax.annotate("", xy=(0.4, 7.5), xytext=(0.4, 0.5),
                    arrowprops=dict(arrowstyle="->", color="black", alpha=0.6))
        ax.text(0.55, 7.3, "$z$", fontsize=12)

    draw(axs[0], with_e=False, title="naive Klein-Nishina (point deposits)")
    draw(axs[1], with_e=True, title="with electron transport + brems")
    fig.tight_layout()
    savefig(fig, "topology_schematic")


# =====================================================================
# Generate auxiliary CSV summaries
# =====================================================================
def write_summary_tables():
    # Photon summary
    rows = []
    for E in [0.05, 0.1, 0.5, 1.0, 1.5, 2.0, E_BI214, Q_BB, E_TL208, 5.0, 10.0]:
        sC = nd.sigma_compton_NIST(E, per_atom=False)
        sP = nd.sigma_pair_NIST(E, per_atom=False)
        sPh = nd.sigma_phot_NIST(E, per_atom=False)
        st = sC + sP + sPh
        mfp = 1.0 / (st * ph.RHO_LXE)
        rows.append([E, sC, sP, sPh, st, sC/st, sP/st, sPh/st, mfp])
    rows = np.array(rows)
    np.savetxt(os.path.join(DATA_DIR, "photon_xsec_summary.csv"), rows,
               header="E[MeV],sC[cm2/g],sP[cm2/g],sPh[cm2/g],sT[cm2/g],bC,bP,bPh,mfp[cm]",
               delimiter=",", fmt="%.6e", comments="")
    print(f"  saved photon_xsec_summary.csv")

    # Electron summary
    rows = []
    for T in [0.05, 0.1, 0.25, 0.5, 1.0, 1.5, 2.0, 2.382, 2.4, 2.448, 2.5, 2.6, 3.0]:
        Sc = nd.dEdx_collision_NIST(T)
        Sr = nd.dEdx_radiative_NIST(T)
        St = Sc + Sr
        R = nd.csda_range_LXe_NIST(T)
        rows.append([T, Sc, Sr, St, Sr/St, R])
    rows = np.array(rows)
    np.savetxt(os.path.join(DATA_DIR, "electron_stopping_summary.csv"), rows,
               header="T[MeV],S_col[MeVcm2/g],S_rad[MeVcm2/g],S_tot[MeVcm2/g],frac_rad,R_LXe[mm]",
               delimiter=",", fmt="%.6e", comments="")
    print(f"  saved electron_stopping_summary.csv")


# =====================================================================
# Main
# =====================================================================
if __name__ == "__main__":
    print("Generating figures...")
    fig_xsec()
    fig_sandia_vs_nist()
    fig_branching()
    fig_kn_diff()
    fig_compton_electron_spectrum()
    fig_stopping_range()
    fig_brems()
    fig_ss_faking()
    fig_topology()
    print("Generating data summaries...")
    write_summary_tables()
    print("Done.")
