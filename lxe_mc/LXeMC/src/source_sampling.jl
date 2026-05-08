"""
Sampling from source flux tables.

Provides functions to draw (E, u) pairs from the 2D probability
distributions stored in `SourceFluxBi214`, `SourceFluxTl208`, and
`SourceRateTable`.
"""


"""
    _sample_pdf_2d(pdf, E_min, E_max, n_E, n_u, rng) -> (E, u)

Sample one (E, u) pair from a 2D discrete PDF matrix.
Uses cumulative distribution + binary search on the flattened array.
Returns values with uniform jitter within the sampled bin.
"""
function _sample_pdf_2d(pdf::Matrix{Float64},
                         E_min::Float64, E_max::Float64,
                         n_E::Int, n_u::Int,
                         rng::AbstractRNG)::Tuple{Float64,Float64}
    total = sum(pdf)
    total <= 0.0 && error("Cannot sample from an empty PDF (sum = 0)")

    # Build flattened CDF
    flat = vec(pdf)
    cdf = cumsum(flat) ./ total

    # Sample
    r = rand(rng)
    idx = searchsortedfirst(cdf, r)
    idx = clamp(idx, 1, length(flat))

    # Convert flat index to (i_E, i_u) — Julia stores column-major
    i_E = mod1(idx, n_E)
    i_u = div(idx - 1, n_E) + 1

    dE = (E_max - E_min) / n_E
    du = 1.0 / n_u

    # Uniform jitter within bin
    E = E_min + (i_E - 0.5) * dE + (rand(rng) - 0.5) * dE
    u = (i_u - 0.5) * du + (rand(rng) - 0.5) * du

    (clamp(E, E_min, E_max), clamp(u, 0.0, 1.0))
end


"""
    sample_from_flux(flux::SourceFluxBi214, rng) -> (E, u)

Sample one (E, u) pair from a Bi-214 flux table.
"""
function sample_from_flux(flux::SourceFluxBi214, rng::AbstractRNG)::Tuple{Float64,Float64}
    _sample_pdf_2d(flux.pdf, flux.E_min, flux.E_max, flux.n_E, flux.n_u, rng)
end


"""
    sample_from_flux(flux::SourceFluxTl208, rng) -> Vector{Tuple{Float64,Float64}}

Sample one event from a Tl-208 flux table. Always returns the main gamma
(E, u). Each companion fires independently with probability
`companion_BR[i] * companion_f[i]`. Returns 1–4 (E, u) tuples.
"""
function sample_from_flux(flux::SourceFluxTl208, rng::AbstractRNG)::Vector{Tuple{Float64,Float64}}
    result = Tuple{Float64,Float64}[]

    # Main gamma (always)
    push!(result, _sample_pdf_2d(flux.pdf_main,
                                  flux.E_min_main, flux.E_max_main,
                                  flux.n_E_main, flux.n_u, rng))

    # Companions (independent Bernoulli)
    for i in 1:length(flux.companion_BR)
        p = flux.companion_BR[i] * flux.companion_f[i]
        rand(rng) >= p && continue
        sum(flux.pdf_companion[i]) <= 0.0 && continue
        push!(result, _sample_pdf_2d(flux.pdf_companion[i],
                                      flux.E_min_companion[i], flux.E_max_companion[i],
                                      flux.n_E_companion[i], flux.n_u, rng))
    end

    result
end


"""
    sample_from_rate_table(rate::SourceRateTable, rng) -> (E, u)

Sample one (E, u) pair from an activity-weighted summed rate table.
"""
function sample_from_rate_table(rate::SourceRateTable, rng::AbstractRNG)::Tuple{Float64,Float64}
    _sample_pdf_2d(rate.pdf_rate, rate.E_min, rate.E_max, rate.n_E, rate.n_u, rng)
end


# =====================================================================
# Surface point sampling and direction reconstruction
# =====================================================================

"""
    sample_barrel_point(R_inner, z_min, z_max, rng) -> (position, normal)

Sample a uniform point on a cylindrical barrel inner surface.
Returns position [cm] and inward normal (radially inward, unit vector).
"""
function sample_barrel_point(R_inner::Float64, z_min::Float64, z_max::Float64,
                              rng::AbstractRNG)::Tuple{Vector{Float64},Vector{Float64}}
    φ = 2π * rand(rng)
    z = z_min + (z_max - z_min) * rand(rng)
    cosφ = cos(φ)
    sinφ = sin(φ)
    pos = Float64[R_inner * cosφ, R_inner * sinφ, z]
    normal = Float64[-cosφ, -sinφ, 0.0]  # radially inward
    (pos, normal)
end


"""
    sample_cap_point(R, aspect_ratio, z_equator, orientation, rng) -> (position, normal)

Sample a point on an ellipsoidal cap inner surface with approximate
area weighting. The cap has semi-axes a = R (equatorial) and
c = R / aspect_ratio (polar).

`z_equator` is the z-coordinate of the cap equator.
`orientation` is `:up` (dome in +z) or `:down` (dome in -z).

Returns position [cm] and inward normal (unit vector pointing toward
the detector interior).
"""
function sample_cap_point(R::Float64, aspect_ratio::Float64,
                           z_equator::Float64, orientation::Symbol,
                           rng::AbstractRNG)::Tuple{Vector{Float64},Vector{Float64}}
    a = R
    c = R / aspect_ratio
    sgn = orientation === :up ? -1.0 : 1.0  # inward normal direction

    # Sample r with area weighting on ellipsoid.
    # Area element dA ∝ r × sqrt(1 + (dz/dr)^2) dr
    # where z(r) = c * sqrt(1 - r^2/a^2), dz/dr = -c*r / (a^2 * sqrt(1 - r^2/a^2))
    # dA ∝ r * sqrt(1 + c^2 * r^2 / (a^4 * (1 - r^2/a^2))) dr
    # For 2:1 aspect ratio the correction is modest; use rejection sampling.
    r_max_weight = sqrt(1.0 + (c / a)^2)  # upper bound on sqrt(1 + (dz/dr)^2) at r→a

    r = 0.0
    while true
        r_try = a * sqrt(rand(rng))  # uniform in r^2
        t = r_try / a
        if t >= 1.0
            continue
        end
        weight = sqrt(1.0 + c^2 * t^2 / (a^2 * (1.0 - t^2)))
        if rand(rng) * r_max_weight <= weight
            r = r_try
            break
        end
    end

    φ = 2π * rand(rng)
    cosφ = cos(φ)
    sinφ = sin(φ)

    # Position on ellipsoid
    z_local = c * sqrt(max(0.0, 1.0 - (r / a)^2))
    z_sign = orientation === :up ? 1.0 : -1.0
    pos = Float64[r * cosφ, r * sinφ, z_equator + z_sign * z_local]

    # Normal from ellipsoid gradient: ∇f = (2x/a², 2y/a², 2z_local/c²)
    # Inward normal points toward detector interior (sgn flips for up vs down)
    nx = r * cosφ / a^2
    ny = r * sinφ / a^2
    nz = z_local / c^2
    n_mag = sqrt(nx^2 + ny^2 + nz^2)
    if n_mag < 1e-12
        normal = Float64[0.0, 0.0, sgn]
    else
        normal = Float64[sgn * nx / n_mag, sgn * ny / n_mag, sgn * nz / n_mag]
    end

    (pos, normal)
end


"""
    reconstruct_direction(u, normal, rng) -> Vector{Float64}

Reconstruct a 3D unit direction vector from u = cos θ to the surface
normal and a random azimuthal angle. Uses `rotate_to_global` to
transform from the local frame (where ẑ = normal) to the global frame.
"""
function reconstruct_direction(u::Float64, normal::Vector{Float64},
                                rng::AbstractRNG)::Vector{Float64}
    sin_θ = sqrt(max(0.0, 1.0 - u^2))
    φ = 2π * rand(rng)
    local_dir = Float64[sin_θ * cos(φ), sin_θ * sin(φ), u]
    rotate_to_global(local_dir, normal)
end


"""
    sample_gamma_from_flux(flux::SourceFluxBi214, surface_sampler, rng) -> SampledGamma

Sample one gamma from a Bi-214 flux table and place it on the detector
surface. `surface_sampler` is a function `rng -> (position, normal)`.
"""
function sample_gamma_from_flux(flux::SourceFluxBi214,
                                 surface_sampler::Function,
                                 rng::AbstractRNG)::SampledGamma
    E, u = sample_from_flux(flux, rng)
    pos, normal = surface_sampler(rng)
    dir = reconstruct_direction(u, normal, rng)
    SampledGamma(E, pos, dir)
end


"""
    sample_gamma_from_flux(flux::SourceFluxTl208, surface_sampler, rng) -> Vector{SampledGamma}

Sample one Tl-208 event from the flux table and place all gammas on the
detector surface. All gammas from one event share the same surface point
(they originate from the same decay location on the cryostat).
Returns 1–4 `SampledGamma` objects.
"""
function sample_gamma_from_flux(flux::SourceFluxTl208,
                                 surface_sampler::Function,
                                 rng::AbstractRNG)::Vector{SampledGamma}
    eu_list = sample_from_flux(flux, rng)
    pos, normal = surface_sampler(rng)  # shared surface point

    gammas = SampledGamma[]
    for (E, u) in eu_list
        dir = reconstruct_direction(u, normal, rng)
        push!(gammas, SampledGamma(E, pos, dir))
    end
    gammas
end
