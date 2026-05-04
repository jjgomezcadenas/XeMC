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

    f_c = (αZ)² [ 1/(1+(αZ)²) + 0.20206 - 0.0369(αZ)² + 0.0083(αZ)⁴ - 0.002(αZ)⁶ ]
"""
function coulomb_correction_fc(Z::Int, α::Float64)::Float64
    a = α * Z
    a2 = a * a
    a2 * (1.0/(1.0 + a2) + 0.20206 - 0.0369*a2 + 0.0083*a2^2 - 0.0020*a2^3)
end


"""
    dsigma_dk_brems(k, T, Z, α, re, me) -> Float64

Differential bremsstrahlung cross section dσ/dk [cm²/MeV/atom].

Takes explicit physical parameters so it can be used with any material.
`Z` = atomic number, `α` = fine-structure constant, `re` = classical
electron radius [cm], `me` = electron mass [MeV].

Returns 0 for k ≤ 0 or k ≥ T.
"""
function dsigma_dk_brems(k::Float64, T::Float64,
                         Z::Int, α::Float64, re::Float64, me::Float64)::Float64
    (k <= 0.0 || k >= T) && return 0.0

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

"""Convenience: take Z from SimConfig (backwards compat for sampling)."""
function dsigma_dk_brems(k::Float64, T::Float64, Z::Int, cfg::SimConfig)::Float64
    dsigma_dk_brems(k, T, Z, cfg.alpha_fs, cfg.re, cfg.me)
end

"""Vectorized form."""
function dsigma_dk_brems(k::Vector{Float64}, T::Float64,
                         Z::Int, α::Float64, re::Float64, me::Float64)::Vector{Float64}
    [dsigma_dk_brems(ki, T, Z, α, re, me) for ki in k]
end

function dsigma_dk_brems(k::Vector{Float64}, T::Float64, Z::Int, cfg::SimConfig)::Vector{Float64}
    [dsigma_dk_brems(ki, T, Z, cfg) for ki in k]
end
