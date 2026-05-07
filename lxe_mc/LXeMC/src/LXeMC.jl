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
- `geometry.jl`: `Cyl`/`LCyl`/`PCyl`/`Detector` hierarchy
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
include("geometry2.jl")
include("sampling.jl")
include("tracking.jl")
include("decays.jl")
include("source.jl")

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
export Cyl, CylShell, Box, Disk
export volume, volume_inner, surface_area, surface_area_inner, surface_area_outer
export R_outer, depth, is_flat
export LCyl, LCylShell, LBox, LDisk, is_inside
export PhysicalVolume, PCyl, PCylShell, PBox, PDisk, mass
export Detector, active_volume, fiducial_volume, find_volume, load_detector
export distance_to_entry, distance_to_exit, next_volume
export RegionTag, PlacementV2, LogicalVolumeV2, DetectorNode, DetectorV2
export FastKernelRegion, FastKernelGeometry
export Cap, DomedContainer, CappedCylinder
export TAG_WORLD, TAG_VACUUM, TAG_STRUCTURAL, TAG_TPC_ACTIVE
export TAG_FV, TAG_SKIN, TAG_PASSIVE_LXE
export load_detector_v2, root_node, node_by_name, child_nodes, detector_summary
export compile_fastkernel_geometry, classify_fastkernel, distance_to_boundary_fastkernel
export validate_detector_v2, tree_dump, sibling_overlaps
export region_tag, is_fv, is_active_lxe, is_veto_lxe, is_passive_lxe
export is_structural, is_vacuum, is_sensitive, is_fv_target, ecut_keV, dz_mm
export find_node_v2, classify_runtime_v2, select_interaction, geometric_prefilter_v2
export select_interaction_fastkernel
export propagate_gamma_fastkernel
export process_event_fastkernel

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
export PropagationResult, GammaPropagationV2Result, EventProcessingResult
export propagate_to_fiducial, propagate_to_lxe, propagate_gamma_v2
export process_event, process_event_from_selections
export classify_lxe_region, veto_threshold
export cluster_deposits_in_z, is_single_site

# --- Decays ---
export DecayGamma, DecayScheme, GammaEmission, SampledGamma
export load_decays, sample_decay, sample_event, sample_isotropic_direction
export FluxTable, spectrum_E, spectrum_u, E_centers, u_centers, survival_fraction
export peak_bin_fraction, off_peak_fraction
export propagate_in_source, random_position_in_volume, generate_source_flux

# --- Convenience paths ---
"""Path to the default detector JSON (LZ)."""
function default_detector_path()
    normpath(joinpath(@__DIR__, "..", "..", "data", "detector_lz.json"))
end
export default_detector_path

"""Path to the Geometry V2 detector JSON (LZ)."""
function default_detector_v2_path()
    normpath(joinpath(@__DIR__, "..", "..", "data", "detector_lz_v2.json"))
end
export default_detector_v2_path

"""Path to the default data directory."""
function default_data_dir()
    normpath(joinpath(@__DIR__, "..", "..", "data"))
end
export default_data_dir

end # module LXeMC
