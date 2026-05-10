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
det  = load_tracking_detector(default_tracking_detector_path(), mats)
fk   = compile_fastkernel_geometry(det)
fv   = compile_fv_geometry(det)
rng  = MersenneTwister(42)

deposits = propagate_gamma_in_fv(2.615, fv, cfg; rng=rng)
ss = is_single_site(deposits, cfg.dz_resolution; E_min=cfg.E_cluster_min)
```

# Module structure

- `config.jl`: `SimConfig` — transport parameters
- `nist_data.jl`: XCOM/ESTAR CSV loaders, log-log interpolation
- `materials.jl`: `Material` — material properties + physics methods
- `geometry.jl`: canonical tracking detector + fast-kernel geometry
- `physics_utils.jl`: bremsstrahlung differential cross section
- `sampling.jl`: MC samplers for all interaction channels
- `tracking.jl`: particle transport, event simulation, clustering
"""
module LXeMC

using JSON
using Printf
using Random

# Source files (order matters: each file may depend on earlier ones)
include("config.jl")
include("nist_data.jl")
include("physics_utils.jl")
include("materials.jl")
include("geometry_core.jl")
include("geometry.jl")
include("sampling.jl")
include("tracking.jl")
include("decays.jl")
include("source_propagation.jl")
include("source_flux.jl")
include("source_sampling.jl")
include("cryostat_sources.jl")
include("pmt_sources.jl")
include("flux_utils.jl")
include("source_dispatch.jl")

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
# Shared primitive solids and wrappers.
export Cyl, CylShell, Box, Disk
export volume, volume_inner, surface_area, surface_area_inner, surface_area_outer
export R_outer, depth, is_flat
export LCyl, LCylShell, LBox, LDisk, is_inside
export PhysicalVolume, PCyl, PCylShell, PBox, PDisk, mass
export distance_to_entry, distance_to_exit
# Canonical detector geometry and compiled fast-kernel support.
export RegionTag, Placement, LogicalVolume, DetectorNode, TrackingDetector
export FastKernelRegion, FastKernelGeometry, FVGeometry
export Cap, DomedContainer, CappedCylinder
export TAG_WORLD, TAG_VACUUM, TAG_STRUCTURAL, TAG_TPC_ACTIVE
export TAG_FV, TAG_SKIN, TAG_PASSIVE_LXE
export root_node, node_by_name, child_nodes, detector_summary
export load_tracking_detector
export compile_fastkernel_geometry, compile_fv_geometry, classify_fastkernel, distance_to_boundary_fastkernel
export validate_tracking_detector, tree_dump, sibling_overlaps
export region_tag, is_fv, is_active_lxe, is_veto_lxe, is_passive_lxe
export is_structural, is_vacuum, is_sensitive, is_fv_target, ecut_keV, dz_mm
export find_tracking_node
export is_inside_fv
export select_interaction_fastkernel
export propagate_gamma_fastkernel
export FastKernelCalibEventResult, process_event_fastkernel_calib
export FastGammaDeposit, FastGammaTrackResult, transport_gamma_fastkernel

# --- Sampling ---
export sample_distance, sample_process
export sample_compton, compton_electron_direction, rotate_to_global
export sample_pair, pair_polar_angle
export sample_photoelectron_angle
export sample_brems, brems_photon_angle

# --- Tracking ---
export Track, ParticleStack, Deposit
export transport_photon!, transport_lepton!
export propagate_gamma, propagate_gamma_in_fv
export EventProcessingResult
export GammaPropagationResult
export process_event
export cluster_deposits_in_z, is_single_site

# --- Decays ---
export DecayGamma, DecayScheme, GammaEmission, SampledGamma
export load_decays, sample_decay, sample_event, sample_gammas, sample_isotropic_direction
export propagate_in_source, random_position_in_volume
export SourceFlux, SourceFluxBi214, SourceFluxTl208, SourceRateTable
export BR_BI214_2448, BR_TL208_2615
export SourceVolumeInfo, load_source_geometry
export generate_flux_bi214, generate_flux_tl208
export propagate_through_layers
export generate_flux_compound_bi214, generate_flux_compound_tl208
export cryostat_barrel_flux, cryostat_top_flux, cryostat_bottom_flux
export pmt_top_flux, pmt_bottom_flux
export sample_from_flux, sample_from_rate_table
export sample_barrel_point, sample_disk_point, sample_cap_point, reconstruct_direction
export sample_gamma_from_flux
export merge_flux_bi214, merge_flux_tl208
export write_pdf_csv, write_flux_bi214_csv, write_flux_tl208_csv, write_rate_table_csv
export write_flux_json, flux_bi214_metadata, flux_tl208_metadata, rate_table_metadata
export load_pdf_csv, load_flux_bi214, load_flux_tl208, load_rate_table
export dispatch_source_flux, merge_dispatch_results
export supported_sources, supported_isotopes
export VirtualEnvelope, ENVELOPE_INSET_CM, make_virtual_envelope, make_surface_sampler

# --- Convenience paths ---
"""Path to the canonical detector-geometry JSON file."""
function default_tracking_detector_path()
    normpath(joinpath(@__DIR__, "..", "..", "data", "detector_lz_v3.json"))
end
export default_tracking_detector_path

"""Path to the default data directory."""
function default_data_dir()
    normpath(joinpath(@__DIR__, "..", "..", "data"))
end
export default_data_dir

end # module LXeMC
