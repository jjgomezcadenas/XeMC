"""
    LXeMC

Stack-based Monte Carlo simulation of photon and lepton transport in
liquid xenon (LXe). Designed for pedagogical study of single-site vs
multi-site event topology in xenon calorimeters.

# Quick start

```julia
using LXeMC, Random

cfg = default_config()
nd  = load_nist_data(cfg)
rng = MersenneTwister(42)

deposits = simulate_event(2.615, nd, cfg; rng=rng)
ss = is_single_site(deposits, cfg.dz_resolution)
```

# Module structure

- `config.jl`: `SimConfig` struct, JSON loader
- `nist_data.jl`: NIST XCOM/ESTAR tables, log-log interpolation
- `physics_utils.jl`: bremsstrahlung differential cross section
- `sampling.jl`: MC samplers for all interaction channels
- `tracking.jl`: particle transport, event simulation, clustering
"""
module LXeMC

using JSON
using Random

# Source files (order matters: each file may depend on earlier ones)
include("config.jl")
include("nist_data.jl")
include("physics_utils.jl")
include("sampling.jl")
include("tracking.jl")

# --- Config ---
export SimConfig, load_config, default_config

# --- NIST data ---
export NISTData, XCOMData, ESTARData
export load_nist_data
export sigma_compton_NIST, sigma_pair_NIST, sigma_phot_NIST, sigma_total_NIST
export mfp_LXe, branching_NIST
export dEdx_collision_NIST, dEdx_radiative_NIST, dEdx_total_NIST
export csda_range_g_per_cm2, csda_range_LXe_mm
export sigma_brems_table

# --- Physics utilities ---
export coulomb_correction_fc, dsigma_dk_brems

# --- Sampling ---
export sample_distance, sample_process
export sample_compton, compton_electron_direction, rotate_to_global
export sample_pair, pair_polar_angle
export sample_phot_shell, sample_photoelectron_angle, sample_atomic_relaxation_K
export sigma_brems_above_kmin, sample_brems, brems_photon_angle

# --- Tracking ---
export Geometry, InfiniteLXe, CylinderLXe, is_inside
export Track, ParticleStack, Deposit
export transport_photon!, transport_lepton!
export simulate_event, cluster_deposits_in_z, is_single_site

end # module LXeMC
