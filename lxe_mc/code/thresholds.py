"""
thresholds.py
=============

Numerical parameters for the LXe Monte Carlo. Centralized here so that
sensitivity studies can be run by changing one file.

See manual section "Thresholds and numerical parameters" for justification.
"""

# Spatial resolution of the calorimeter [cm]
DZ_RESOLUTION_CM = 0.30   # 3 mm

# Photon transport thresholds
EGAMMA_CUT_MEV = 0.010    # below this, deposit residual energy locally
                          # (mfp at 10 keV in LXe is < 50 um -> always local)

# Lepton transport thresholds
T_E_CUT_MEV = 0.050       # below this, deposit residual energy locally
                          # CSDA range at 50 keV in LXe is ~30 um

# Bremsstrahlung threshold
K_MIN_MEV = 0.050         # photons below this energy are added to the local
                          # electron deposit (their mfp in LXe is ~3 mm, comparable
                          # to DZ resolution; not separable as a cluster)

# Stepping parameters
F_RANGE = 0.05            # max step = F_RANGE * residual range
DS_FLOOR_CM = 0.001       # 10 um, min step
DS_CEIL_CM = 0.10         # 1 mm, max step

# Generation cap (pathological-event protection)
GEN_CAP = 100

# Sandia validity domain for photoelectric
SANDIA_VALID_MAX_MEV = 0.50   # use Sandia below this; XCOM interpolation above
