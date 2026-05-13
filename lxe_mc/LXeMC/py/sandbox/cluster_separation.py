DETECTOR = "DARWIN"   # or "XENON1T", "XENONnT", "LZ"

# Detector-specific parameters
params = {
    "XENON1T":  dict(L_FV=0.50, tau_e=650e-6,  sigma_E_rel=0.008),
    "XENONnT":  dict(L_FV=0.75, tau_e=2000e-6, sigma_E_rel=0.008),
    "LZ":       dict(L_FV=0.70, tau_e=1000e-6, sigma_E_rel=0.008),
    "DARWIN":   dict(L_FV=1.30, tau_e=5000e-6, sigma_E_rel=0.008),
}[DETECTOR]

# Derived
sigma_L_FV = D_L_tilde * sqrt(params["L_FV"])         # for Δz_min
DZ_MIN     = 3 * sigma_L_FV                           # ≈ 3 mm for L_FV~0.7 m
E_CUT_at_Q = 3 * params["sigma_E_rel"] * 2458         # ≈ 60 keV
