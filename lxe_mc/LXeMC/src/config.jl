"""
    SimConfig

Simulation configuration for the LXe Monte Carlo.

All energies are in MeV, lengths in cm. Loaded from a JSON file
(default: `data/sim_config.json`) so that sensitivity studies can be
run by editing the data file without touching source code.

# Sections
- **Xenon properties**: atomic number, mass, density, shell energies,
  fluorescence yield, mean excitation energy.
- **Transport cuts**: minimum energies for photon and lepton tracking;
  bremsstrahlung photon threshold.
- **Stepping parameters**: fixed step size for lepton transport.
- **Safety**: generation cap to prevent runaway cascades.
- **Clustering**: z-resolution for single-site / multi-site classification.
- **Physical constants**: electron mass, classical radius, fine-structure
  constant, Avogadro number.
"""
struct SimConfig
    # --- Xenon properties ---
    Z::Int                  # atomic number
    A::Float64              # standard atomic weight [g/mol]
    rho_LXe::Float64        # liquid xenon density [g/cm³]
    EK::Float64             # K-shell binding energy [MeV]
    EL1::Float64            # L1-shell binding energy [MeV]
    EL2::Float64            # L2-shell binding energy [MeV]
    EL3::Float64            # L3-shell binding energy [MeV]
    omega_K::Float64        # K-shell fluorescence yield
    EK_alpha::Float64       # dominant K-alpha fluorescence line [MeV]
    I_Xe::Float64           # mean excitation energy (ICRU-37) [MeV]

    # --- Transport cuts ---
    Egamma_cut::Float64     # photon tracking cutoff [MeV]
    Te_cut::Float64         # lepton tracking cutoff [MeV]
    k_min::Float64          # bremsstrahlung photon threshold [MeV]

    # --- Stepping ---
    ds_step::Float64        # fixed step size for lepton transport [cm]
    f_range::Float64        # (legacy) fraction-of-range step sizing
    ds_floor::Float64       # (legacy) minimum step size [cm]
    ds_ceil::Float64        # (legacy) maximum step size [cm]

    # --- Safety ---
    generation_cap::Int     # max cascade generations

    # --- Clustering ---
    dz_resolution::Float64  # z-resolution for SS/MS classification [cm]

    # --- Physical constants ---
    me::Float64             # electron mass [MeV/c²]
    re::Float64             # classical electron radius [cm]
    alpha_fs::Float64       # fine-structure constant
    N_A::Float64            # Avogadro number [mol⁻¹]

    # --- Derived ---
    n_atom::Float64         # atom number density in LXe [cm⁻³]
end


"""
    load_config(path::AbstractString) -> SimConfig

Read a JSON configuration file and return a `SimConfig`.
The JSON schema must match `data/sim_config.json`.
"""
function load_config(path::AbstractString)::SimConfig
    raw = open(path, "r") do io
        JSON.parse(io)
    end

    xe  = raw["xenon"]
    pc  = raw["photon_cuts"]
    lc  = raw["lepton_cuts"]
    br  = raw["bremsstrahlung"]
    st  = raw["stepping"]
    sa  = raw["safety"]
    cl  = raw["clustering"]
    co  = raw["constants"]

    Z   = Int(xe["Z"])
    A   = Float64(xe["A"])
    rho = Float64(xe["rho_LXe_g_cm3"])
    N_A = Float64(co["N_A"])

    SimConfig(
        Z, A, rho,
        Float64(xe["EK_MeV"]),
        Float64(xe["EL1_MeV"]),
        Float64(xe["EL2_MeV"]),
        Float64(xe["EL3_MeV"]),
        Float64(xe["omega_K"]),
        Float64(xe["EK_alpha_MeV"]),
        Float64(xe["I_MeV"]),
        Float64(pc["Egamma_cut_MeV"]),
        Float64(lc["Te_cut_MeV"]),
        Float64(br["k_min_MeV"]),
        Float64(st["ds_step_cm"]),
        Float64(st["f_range"]),
        Float64(st["ds_floor_cm"]),
        Float64(st["ds_ceil_cm"]),
        Int(sa["generation_cap"]),
        Float64(cl["dz_resolution_cm"]),
        Float64(co["me_MeV"]),
        Float64(co["re_cm"]),
        Float64(co["alpha_fs"]),
        N_A,
        N_A * rho / A   # n_atom: atom number density [cm⁻³]
    )
end


"""
    default_config() -> SimConfig

Load the default configuration from `data/sim_config.json`,
located relative to the package source directory.
"""
function default_config()::SimConfig
    default_path = normpath(joinpath(@__DIR__, "..", "..", "data", "sim_config.json"))
    load_config(default_path)
end
