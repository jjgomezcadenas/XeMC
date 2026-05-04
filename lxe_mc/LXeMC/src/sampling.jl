"""
Monte Carlo samplers for photon and lepton interactions in xenon.

Each function samples the kinematics of one interaction, returning the
sampled quantities (scattered energy, angle, etc.). The sampling strategy
is documented in each function's docstring.

# Sampling strategies used

| Sampler                   | Method                                          |
|:--------------------------|:------------------------------------------------|
| `sample_distance`         | Inverse CDF of exponential distribution         |
| `sample_process`          | Direct discrete sampling                        |
| `sample_compton`          | Composition + rejection (Butcher-Messel / KN)   |
| `sample_pair`             | Flat envelope + rejection (BH shape)            |
| `pair_polar_angle`        | Composition of two gamma distributions (Urban)  |
| `sample_phot_shell`       | Direct discrete branching                       |
| `sample_photoelectron_angle` | Rejection from uniform (Sauter distribution) |
| `sample_atomic_relaxation_K` | Bernoulli (fluorescence yield)              |
| `sample_brems`            | 1/k envelope + rejection (BH-Tsai shape)        |
"""

using Random


# =====================================================================
# Distance to next interaction
# =====================================================================

"""
    sample_distance(Σ_tot::Float64, rng::AbstractRNG) -> Float64

Sample free-flight distance [cm] from an exponential distribution
with macroscopic cross section `Σ_tot` [cm⁻¹].

**Strategy**: inverse CDF.  If ξ ~ U(0,1), then s = -ln(ξ)/Σ.
"""
function sample_distance(Σ_tot::Float64, rng::AbstractRNG)::Float64
    -log(rand(rng)) / Σ_tot
end


# =====================================================================
# Discrete process selection
# =====================================================================

"""
    sample_process(P_C::Float64, P_P::Float64, P_Ph::Float64,
                   rng::AbstractRNG) -> Symbol

Select one of `:compton`, `:pair`, or `:photoelectric` with the given
probabilities. The three probabilities must sum to 1.

**Strategy**: direct partitioning of U(0,1).
"""
function sample_process(P_C::Float64, P_P::Float64, P_Ph::Float64,
                        rng::AbstractRNG)::Symbol
    r = rand(rng)
    r < P_C && return :compton
    r < P_C + P_P && return :pair
    return :photoelectric
end


# =====================================================================
# Compton scattering (Klein-Nishina)
# =====================================================================

"""
    sample_compton(E::Float64, cfg::SimConfig, rng::AbstractRNG)
        -> (E_scattered::Float64, cos_θ::Float64)

Sample Compton scattering of a photon with energy `E` [MeV].
Returns the scattered photon energy and cosine of the scattering angle.

**Strategy**: composition + rejection after Butcher & Messel (1960),
as implemented in Geant4 (PRM v11.4).

Two proposal distributions cover the Klein-Nishina shape:
- Branch 1 (weight α₁ = -ln ε₀): ε = ε₀^ξ  (samples the low-ε peak)
- Branch 2 (weight α₂ = (1-ε₀²)/2): ε = √(ε₀² + (1-ε₀²)ξ)  (samples the high-ε plateau)

where ε = E'/E ∈ [ε₀, 1], ε₀ = m_e/(m_e + 2E).

Rejection function: g(ε) = 1 - ε sin²θ/(1+ε²), with sin²θ from kinematics.
"""
function sample_compton(E::Float64, cfg::SimConfig, rng::AbstractRNG)::Tuple{Float64,Float64}
    me = cfg.me
    eps0 = me / (me + 2.0 * E)
    alpha1 = -log(eps0)
    alpha2 = (1.0 - eps0^2) / 2.0

    while true
        r1 = rand(rng)
        r2 = rand(rng)
        r3 = rand(rng)

        eps = if r1 < alpha1 / (alpha1 + alpha2)
            eps0^r2
        else
            sqrt(eps0^2 + (1.0 - eps0^2) * r2)
        end

        one_minus_cos = (1.0 - eps) * me / (E * eps)
        sin2 = clamp(one_minus_cos * (2.0 - one_minus_cos), 0.0, 1.0)
        g = 1.0 - eps * sin2 / (1.0 + eps^2)

        if g >= r3
            cos_theta = 1.0 - one_minus_cos
            return (eps * E, cos_theta)
        end
    end
end


"""
    compton_electron_direction(cos_θ_γ::Float64, ϕ_γ::Float64,
                               E::Float64, dir_in::Vector{Float64},
                               cfg::SimConfig) -> Vector{Float64}

Compute the recoil electron direction (unit vector) given the scattered
photon direction. The electron recoils in the plane defined by the
incoming and scattered photon momenta.

Uses momentum conservation: p_e = p_γ_in - p_γ_out in the local frame,
then rotates to the global frame.
"""
function compton_electron_direction(cos_θ_γ::Float64, ϕ_γ::Float64,
                                    E::Float64, dir_in::Vector{Float64},
                                    cfg::SimConfig)::Vector{Float64}
    me = cfg.me
    sin_θ_γ = sqrt(max(0.0, 1.0 - cos_θ_γ^2))
    Egp = E / (1.0 + (E / me) * (1.0 - cos_θ_γ))

    p_in  = Float64[0.0, 0.0, E]
    p_out = Egp .* Float64[sin_θ_γ*cos(ϕ_γ), sin_θ_γ*sin(ϕ_γ), cos_θ_γ]
    p_e = p_in .- p_out
    nrm = sqrt(p_e[1]^2 + p_e[2]^2 + p_e[3]^2)

    if nrm < 1e-12
        return copy(dir_in)
    end

    n_e_local = p_e ./ nrm
    rotate_to_global(n_e_local, dir_in)
end


# =====================================================================
# Frame rotation
# =====================================================================

"""
    rotate_to_global(local_vec::Vector{Float64},
                     ref_dir::Vector{Float64}) -> Vector{Float64}

Rotate a unit vector from a local frame (where ẑ = ref_dir) to the
global frame. Handles the degenerate case where ref_dir ≈ ±ẑ.
"""
function rotate_to_global(local_vec::Vector{Float64},
                          ref_dir::Vector{Float64})::Vector{Float64}
    rd = ref_dir ./ sqrt(ref_dir[1]^2 + ref_dir[2]^2 + ref_dir[3]^2)

    if abs(rd[3]) > 0.99999
        return local_vec .* sign(rd[3])
    end

    # Orthonormal triad: e1, e2, rd
    rho_perp = sqrt(rd[1]^2 + rd[2]^2)
    e1 = Float64[-rd[3]*rd[1], -rd[3]*rd[2], -(rd[3]^2 - 1.0)] ./ rho_perp
    e2 = Float64[-rd[2], rd[1], 0.0] ./ rho_perp

    local_vec[1] .* e1 .+ local_vec[2] .* e2 .+ local_vec[3] .* rd
end


# =====================================================================
# Pair production
# =====================================================================

"""
    sample_pair(E::Float64, cfg::SimConfig, rng::AbstractRNG) -> Float64

Sample the energy partition ε = E₊/E_γ ∈ [ε₀, 1-ε₀] for pair production,
where ε₀ = m_e/E_γ.

**Strategy**:
- Below 2 MeV: the Bethe-Heitler differential is nearly flat in ε,
  so uniform sampling is correct to ~1%.
- Above 2 MeV: flat envelope with rejection weight w(ε) = ε² + (1-ε)²,
  which captures the U-shaped BH differential. Maximum of w is 1 at ε=0 or 1,
  so acceptance is always > 50%.
"""
function sample_pair(E::Float64, cfg::SimConfig, rng::AbstractRNG)::Float64
    eps_min = cfg.me / E
    eps_min >= 0.5 && return 0.5

    if E < 2.0
        return eps_min + rand(rng) * (1.0 - 2.0 * eps_min)
    end

    while true
        eps = eps_min + rand(rng) * (1.0 - 2.0 * eps_min)
        w = eps^2 + (1.0 - eps)^2
        if rand(rng) < w
            return eps
        end
    end
end


"""
    pair_polar_angle(E_lepton::Float64, cfg::SimConfig,
                     rng::AbstractRNG) -> Float64

Sample the polar angle θ [rad] of a pair-produced lepton with total
energy `E_lepton` [MeV].

**Strategy**: Geant4/Urban parameterization. The distribution is a mixture
of two gamma-like components:

    f(u) ∝ c₁ u exp(-au) + c₂ u exp(-3au),  a = 5/8

sampled by choosing a branch (probability 9/(9+d), d=27) and then
drawing from u·exp(-bu) via the sum of two exponentials:
u = (-ln ξ₁ - ln ξ₂)/b.

The physical angle is θ = (m_e/E_lepton) × u, which is forward-peaked
as expected from relativistic kinematics.
"""
function pair_polar_angle(E_lepton::Float64, cfg::SimConfig,
                          rng::AbstractRNG)::Float64
    a = 5.0 / 8.0
    d = 27.0
    p1 = 9.0 / (9.0 + d)

    b = rand(rng) < p1 ? a : 3.0 * a
    u = (-log(rand(rng)) - log(rand(rng))) / b

    (cfg.me / E_lepton) * u
end


# =====================================================================
# Photoelectric effect
# =====================================================================

"""
    sample_photoelectron_angle(T_e::Float64, cfg::SimConfig,
                               rng::AbstractRNG) -> Float64

Sample the photoelectron polar angle θ [rad] relative to the photon
direction.

**Strategy**: rejection sampling from a uniform envelope over cos θ ∈ [-1, 1].

The target is the Sauter (zeroth-order) distribution:

    dσ/d(cos θ) ∝ sin²θ / (1 - β cos θ)⁴

where β = v/c of the photoelectron. The peak value is bounded by M = 4
(loose but safe). Acceptance rate is ~25-50% depending on energy.
"""
function sample_photoelectron_angle(T_e::Float64, cfg::SimConfig,
                                    rng::AbstractRNG)::Float64
    g = 1.0 + T_e / cfg.me
    beta = sqrt(1.0 - 1.0 / g^2)
    M = 4.0

    while true
        ct = -1.0 + 2.0 * rand(rng)
        f = (1.0 - ct^2) / (1.0 - beta * ct)^4
        if rand(rng) * M < f
            return acos(ct)
        end
    end
end


# =====================================================================
# Bremsstrahlung
# =====================================================================

"""
    sigma_brems_above_kmin(T, k_min, Z, cfg) -> Float64

Total bremsstrahlung cross section [cm²/atom] for photon energies
k > k_min, computed by trapezoidal integration of the BH-Tsai
differential on a log-spaced grid.
"""
function sigma_brems_above_kmin(T::Float64, k_min::Float64,
                                Z::Int, cfg::SimConfig)::Float64
    k_min >= T && return 0.0
    k = exp.(range(log(k_min), log(T * 0.9999), length=100))
    ds = dsigma_dk_brems(k, T, Z, cfg)
    s = 0.0
    @inbounds for i in 2:length(k)
        s += 0.5 * (ds[i-1] + ds[i]) * (k[i] - k[i-1])
    end
    s
end

# Legacy: no Z argument, uses Xe (Z=54) — for brems table building
function sigma_brems_above_kmin(T::Float64, k_min::Float64,
                                cfg::SimConfig)::Float64
    sigma_brems_above_kmin(T, k_min, 54, cfg)
end


"""
    sample_brems(T, k_min, Z, cfg, rng) -> Union{Float64, Nothing}

Sample bremsstrahlung photon energy k > k_min [MeV].

**Strategy**: 1/k envelope rejection. Acceptance rate typically > 60%.
"""
function sample_brems(T::Float64, k_min::Float64, Z::Int,
                      cfg::SimConfig, rng::AbstractRNG)::Union{Float64,Nothing}
    k_min >= T && return nothing

    k_grid = exp.(range(log(k_min), log(T * 0.9999), length=50))
    ds = dsigma_dk_brems(k_grid, T, Z, cfg)
    M = maximum(k_grid .* ds) * 1.05

    while true
        r = rand(rng)
        k = k_min * (T / k_min)^r
        ds_k = dsigma_dk_brems(k, T, Z, cfg)
        if rand(rng) * M < k * ds_k
            return k
        end
    end
end


"""
    brems_photon_angle(T_e::Float64, cfg::SimConfig,
                       rng::AbstractRNG) -> Float64

Sample bremsstrahlung photon polar angle [rad] relative to the electron
direction. Uses the same Urban functional form as [`pair_polar_angle`](@ref),
with E_lepton = T_e + m_e.
"""
function brems_photon_angle(T_e::Float64, cfg::SimConfig,
                            rng::AbstractRNG)::Float64
    pair_polar_angle(T_e + cfg.me, cfg, rng)
end
