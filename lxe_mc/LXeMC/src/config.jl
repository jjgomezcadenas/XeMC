"""
    SimConfig

Simulation transport parameters for the LXe Monte Carlo.

Material properties (Z, A, density, cross sections) are now in
`Material`; detector geometry lives in the canonical geometry layer
(`TrackingDetector`, `FastKernelGeometry`). This struct
holds only transport thresholds, stepping parameters, and physical constants.

Loaded from `data/sim_config.json`.
"""
struct SimConfig
    # --- Transport cuts ---
    Egamma_cut::Float64     # photon tracking cutoff [MeV]
    E_roi_floor::Float64    # gamma energy below which the photon cannot contribute to the ROI [MeV]
    Te_cut::Float64         # lepton tracking cutoff [MeV]
    k_min::Float64          # bremsstrahlung photon threshold [MeV]

    # --- Stepping ---
    ds_step::Float64        # fixed step size for lepton transport [cm]

    # --- Safety ---
    generation_cap::Int     # max cascade generations

    # --- Veto thresholds (detection, not physics) ---
    veto_TPC::Float64       # min visible deposit in TPC outside FV [MeV]
    veto_skin::Float64      # min visible deposit in LXe skin [MeV]

    # --- Clustering ---
    dz_resolution::Float64  # z-resolution for SS/MS classification [cm]
    E_cluster_min::Float64  # minimum cluster energy for SS/MS classification [MeV]

    # --- Physical constants ---
    me::Float64             # electron mass [MeV/c²]
    re::Float64             # classical electron radius [cm]
    alpha_fs::Float64       # fine-structure constant
    N_A::Float64            # Avogadro number [mol⁻¹]
end


"""
    load_config(path::AbstractString) -> SimConfig

Read a JSON configuration file and return a `SimConfig`.
"""
function load_config(path::AbstractString)::SimConfig
    raw = open(path, "r") do io
        JSON.parse(io)
    end

    pc  = raw["photon_cuts"]
    lc  = raw["lepton_cuts"]
    br  = raw["bremsstrahlung"]
    st  = raw["stepping"]
    sa  = raw["safety"]
    cl  = raw["clustering"]
    co  = raw["constants"]

    vt  = raw["veto"]

    SimConfig(
        Float64(pc["Egamma_cut_MeV"]),
        Float64(pc["E_roi_floor_MeV"]),
        Float64(lc["Te_cut_MeV"]),
        Float64(br["k_min_MeV"]),
        Float64(st["ds_step_cm"]),
        Int(sa["generation_cap"]),
        Float64(vt["veto_TPC_MeV"]),
        Float64(vt["veto_skin_MeV"]),
        Float64(cl["dz_resolution_cm"]),
        Float64(cl["E_cluster_min_MeV"]),
        Float64(co["me_MeV"]),
        Float64(co["re_cm"]),
        Float64(co["alpha_fs"]),
        Float64(co["N_A"])
    )
end


"""
    default_config() -> SimConfig

Load the default configuration from `data/sim_config.json`.
"""
function default_config()::SimConfig
    default_path = normpath(joinpath(@__DIR__, "..", "..", "data", "sim_config.json"))
    load_config(default_path)
end
