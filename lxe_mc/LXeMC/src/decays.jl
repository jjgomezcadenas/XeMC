"""
Radioactive decay event generators for background gamma sources.

Loads decay schemes from `data/decays.json` and provides samplers
that return lists of gammas (energy + isotropic direction) for each
decay type.

# Decay types

- **Bi-214** (U-238 late chain): single 2.448 MeV gamma per event.
- **Tl-208** (Th-232 late chain): 2.615 MeV main gamma + companion
  gammas (583 keV at 85%, 511 keV at 23%, 861 keV at 12%).
  Companions are sampled independently; an event may have 0, 1, or
  multiple companions (each drawn by Bernoulli trial on its BR).

All gammas are emitted isotropically with independent directions.
"""

using Random


"""
    DecayGamma

A gamma line in a decay scheme: energy and branching ratio.
"""
struct DecayGamma
    E_MeV::Float64
    BR::Float64
end


"""
    DecayScheme

Complete decay scheme for one isotope.
- `name`: e.g., "Bi214", "Tl208"
- `parent_chain`: "U238" or "Th232"
- `BR_from_chain`: probability per parent chain decay
- `gammas`: main gamma lines (always emitted)
- `companions`: companion lines (each sampled independently by BR)
"""
struct DecayScheme
    name::String
    parent_chain::String
    BR_from_chain::Float64
    gammas::Vector{DecayGamma}
    companions::Vector{DecayGamma}
end


"""
    load_decays(; data_dir=nothing) -> Dict{String, DecayScheme}

Load decay schemes from `data/decays.json`.
"""
function load_decays(; data_dir::Union{String,Nothing}=nothing)::Dict{String,DecayScheme}
    if data_dir === nothing
        data_dir = normpath(joinpath(@__DIR__, "..", "..", "data"))
    end

    raw = open(joinpath(data_dir, "decays.json"), "r") do io
        JSON.parse(io)
    end

    schemes = Dict{String,DecayScheme}()
    for (name, d) in raw
        startswith(name, "_") && continue

        gammas = DecayGamma[]
        for g in get(d, "gammas", [])
            push!(gammas, DecayGamma(Float64(g["E_MeV"]), Float64(g["BR"])))
        end

        companions = DecayGamma[]
        for g in get(d, "companions", [])
            push!(companions, DecayGamma(Float64(g["E_MeV"]), Float64(g["BR"])))
        end

        schemes[name] = DecayScheme(
            name,
            String(d["parent_chain"]),
            Float64(d["BR_from_chain"]),
            gammas,
            companions
        )
    end
    schemes
end


"""
    GammaEmission

A single gamma emitted in a decay event: energy and direction (unit vector).
"""
struct GammaEmission
    E_MeV::Float64
    direction::Vector{Float64}
end


struct SampledGamma
    E_MeV::Float64
    position::Vector{Float64}
    direction::Vector{Float64}
end


"""
    sample_isotropic_direction(rng) -> Vector{Float64}

Sample a uniformly distributed unit vector on the sphere.
"""
function sample_isotropic_direction(rng::AbstractRNG)::Vector{Float64}
    cos_θ = -1.0 + 2.0 * rand(rng)
    sin_θ = sqrt(1.0 - cos_θ^2)
    φ = 2π * rand(rng)
    Float64[sin_θ*cos(φ), sin_θ*sin(φ), cos_θ]
end


"""
    sample_decay(scheme::DecayScheme, rng) -> Vector{GammaEmission}

Generate one decay event: sample which gammas are emitted and assign
each an independent isotropic direction.

- Main gammas (BR=1.0) are always emitted.
- Companions are sampled independently: each is emitted with
  probability = its BR (Bernoulli trial).
"""
function sample_decay(scheme::DecayScheme, rng::AbstractRNG)::Vector{GammaEmission}
    emissions = GammaEmission[]

    # Main gammas (always present)
    for g in scheme.gammas
        if rand(rng) < g.BR
            push!(emissions, GammaEmission(g.E_MeV, sample_isotropic_direction(rng)))
        end
    end

    # Companions (independent Bernoulli trials)
    for g in scheme.companions
        if rand(rng) < g.BR
            push!(emissions, GammaEmission(g.E_MeV, sample_isotropic_direction(rng)))
        end
    end

    emissions
end


function sample_event(isotope, source; calib::Bool=false, rng::AbstractRNG=Random.default_rng())::Vector{SampledGamma}
    calib || error("sample_event(...; calib=false) is not implemented yet")

    multiplicity = rand(rng, 0:2)
    gammas = SampledGamma[]
    for _ in 1:multiplicity
        push!(gammas, SampledGamma(
            2.615,
            Float64[0.0, 0.0, 160.0],
            sample_isotropic_direction(rng)
        ))
    end
    gammas
end
