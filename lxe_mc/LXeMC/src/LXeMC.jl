"""
    LXeMC

Stack-based Monte Carlo simulation of photon and lepton transport in
xenon detectors. Designed for studying single-site vs multi-site event
topology.

# Quick start

```julia
using LXeMC, Random

cfg  = default_config()
mats = load_materials(cfg)
det  = load_detector(default_detector_path(), mats)
vol  = active_volume(det)
rng  = MersenneTwister(42)

deposits = simulate_event(2.615, vol, cfg; rng=rng)
ss = is_single_site(deposits, cfg.dz_resolution; E_min=cfg.E_cluster_min)
```

# Module structure

- `config.jl`: `SimConfig` — transport parameters
- `nist_data.jl`: XCOM/ESTAR CSV loaders, log-log interpolation
- `materials.jl`: `Material` — material properties + physics methods
- `geometry.jl`: `Cyl`/`LCyl`/`PCyl`/`RCyl`/`Detector` hierarchy
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
include("materials.jl")
include("geometry.jl")
include("sampling.jl")
include("tracking.jl")

# --- Config ---
export SimConfig, load_config, default_config

# --- NIST data (low-level, used internally by Material) ---
export XCOMData, ESTARData

# --- Materials ---
export Material, ElementData, load_elements, load_materials
export sigma_compton, sigma_pair, sigma_phot, sigma_total, sigma_three
export mu_over_rho, mfp, branching
export dEdx_collision, dEdx_radiative, dEdx_total
export csda_range_g_per_cm2, csda_range_mm, sigma_brems, brems_rejection_M

# --- Physics utilities ---
export coulomb_correction_fc, dsigma_dk_brems

# --- Geometry ---
export Cyl, Box, volume, surface_area
export LCyl, LBox, is_inside
export PhysicalVolume, PCyl, PBox, mass
export RCyl, activity_U238, activity_Th232, gamma_flux
export Detector, active_volume, find_volume, load_detector

# --- Sampling ---
export sample_distance, sample_process
export sample_compton, compton_electron_direction, rotate_to_global
export sample_pair, pair_polar_angle
export sample_photoelectron_angle
export sample_brems, brems_photon_angle

# --- Tracking ---
export Track, ParticleStack, Deposit
export transport_photon!, transport_lepton!
export simulate_event, simulate_event_photon_only
export cluster_deposits_in_z, is_single_site

# --- Convenience paths ---
"""Path to the default detector JSON (LZ)."""
function default_detector_path()
    normpath(joinpath(@__DIR__, "..", "..", "data", "detector_lz.json"))
end
export default_detector_path

"""Path to the default data directory."""
function default_data_dir()
    normpath(joinpath(@__DIR__, "..", "..", "data"))
end
export default_data_dir

end # module LXeMC
