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
