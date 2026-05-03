"""
Minimal analytic physics formulas needed by the MC samplers.

Only the bremsstrahlung differential cross section and its prerequisite
(the Coulomb correction) are included here. All other cross sections
come from the NIST interpolation tables in `nist_data.jl`.

# Formulas

- [`coulomb_correction_fc`](@ref): Davies-Bethe-Maximon Coulomb correction f_c(Z).
- [`dsigma_dk_brems`](@ref): Bethe-Heitler-Tsai differential bremsstrahlung
  cross section dσ/dk [cm²/MeV/atom].
"""


"""
    coulomb_correction_fc(Z::Int, α::Float64) -> Float64

Davies-Bethe-Maximon Coulomb correction f_c(Z) for bremsstrahlung
and pair production in the field of a nucleus with atomic number `Z`.

The correction accounts for the distortion of the electron wave function
by the nuclear Coulomb field beyond Born approximation:

    f_c = (αZ)² [ 1/(1+(αZ)²) + 0.20206 - 0.0369(αZ)² + 0.0083(αZ)⁴ - 0.002(αZ)⁶ ]

For Xe (Z=54): f_c ≈ 0.0325.
"""
function coulomb_correction_fc(Z::Int, α::Float64)::Float64
    a = α * Z
    a2 = a * a
    a2 * (1.0/(1.0 + a2) + 0.20206 - 0.0369*a2 + 0.0083*a2^2 - 0.0020*a2^3)
end


"""
    dsigma_dk_brems(k::Float64, T::Float64, cfg::SimConfig) -> Float64

Differential bremsstrahlung cross section dσ/dk [cm²/MeV/atom] for an
electron of kinetic energy `T` [MeV] emitting a photon of energy `k` [MeV]
in xenon.

Uses the Bethe-Heitler-Tsai screened formula (Tsai, Rev. Mod. Phys. 46, 1974):

    dσ/dk = (4α r_e²) / (3k) × {
        [y² + 2(1 + (1−y)²)] × [Z²(F_el − f_c) + Z F_inel]
        + (1−y)(Z² + Z) / 3
    }

where:
- y = k / (T + m_e) is the photon energy fraction of total energy
- F_el  = ln(184.15 / Z^{1/3}) is the elastic atomic form factor
- F_inel = ln(1194 / Z^{2/3})  is the inelastic atomic form factor
- f_c is the Coulomb correction [`coulomb_correction_fc`](@ref)

Returns 0 for k ≤ 0 or k ≥ T (kinematically forbidden).

Accuracy: ~10% in the MeV range; for production use Seltzer-Berger tables.
"""
function dsigma_dk_brems(k::Float64, T::Float64, cfg::SimConfig)::Float64
    (k <= 0.0 || k >= T) && return 0.0

    Z = cfg.Z
    α = cfg.alpha_fs
    re = cfg.re
    me = cfg.me

    Etot = T + me
    y = k / Etot

    fc = coulomb_correction_fc(Z, α)
    Fel   = log(184.15 / Z^(1.0/3.0))
    Finel = log(1194.0 / Z^(2.0/3.0))

    Z2 = Float64(Z * Z)
    Zf = Float64(Z)

    (4.0 * α * re^2 / (3.0 * k)) * (
        (y^2 + 2.0*(1.0 + (1.0 - y)^2)) * (Z2*(Fel - fc) + Zf*Finel)
        + (1.0 - y) * (Z2 + Zf) / 3.0
    )
end


"""
    dsigma_dk_brems(k::Vector{Float64}, T::Float64, cfg::SimConfig) -> Vector{Float64}

Vectorized form: evaluate dσ/dk at each photon energy in `k`.
"""
function dsigma_dk_brems(k::Vector{Float64}, T::Float64, cfg::SimConfig)::Vector{Float64}
    [dsigma_dk_brems(ki, T, cfg) for ki in k]
end
