"""
Source-side gamma propagation helpers.

This file covers source geometry sampling and photon transport inside a
radioactive source volume before handoff to the canonical detector
tracking path.
"""

using Random


# =====================================================================
# Source propagation
# =====================================================================

"""
    propagate_in_source(E, pos, dir, source_vol, cfg, rng)
        -> (status, E_exit, pos_exit, dir_exit)

Photon-only transport through the source material until the gamma
exits the volume or is absorbed. Uses Klein-Nishina Compton scattering.

Returns:
- `:exited` — gamma left the source volume with energy `E_exit`
- `:absorbed` — gamma was absorbed (photoelectric or pair) inside source
"""
function propagate_in_source(E::Float64, pos::Vector{Float64}, dir::Vector{Float64},
                             source_vol::PhysicalVolume, cfg::SimConfig,
                             rng::AbstractRNG)::Tuple{Symbol,Float64,Vector{Float64},Vector{Float64}}
    mat = source_vol.material

    while E >= cfg.Egamma_cut
        sC, sP, sPh = sigma_three(mat, E)
        s_tot = sC + sP + sPh
        Σ_tot = mat.n_atom * s_tot

        s = sample_distance(Σ_tot, rng)
        pos .= pos .+ dir .* s

        # Exited source volume
        if !is_inside(source_vol, pos)
            return (:exited, E, pos, dir)
        end

        # Interaction inside source
        proc = sample_process(sC/s_tot, sP/s_tot, sPh/s_tot, rng)

        if proc === :compton
            Egp, cos_t = sample_compton(E, cfg, rng)
            ϕ = 2π * rand(rng)
            sin_t = sqrt(max(0.0, 1.0 - cos_t^2))
            local_vec = Float64[sin_t*cos(ϕ), sin_t*sin(ϕ), cos_t]
            dir = rotate_to_global(local_vec, dir)
            E = Egp
        else
            # Photoelectric or pair: absorbed
            return (:absorbed, 0.0, pos, dir)
        end
    end

    # Below Egamma_cut
    (:absorbed, 0.0, pos, dir)
end


"""
    random_position_in_volume(vol::PhysicalVolume, rng) -> Vector{Float64}

Sample a uniform random position inside a physical volume.
Supports PCyl (solid cylinder) and PCylShell (hollow cylinder).
Uses rejection sampling for shells.
"""
function random_position_in_volume(vol::PCyl, rng::AbstractRNG)::Vector{Float64}
    s = vol.logical.solid
    c = vol.logical.position
    # Uniform in cylinder: r = R*sqrt(ξ), z uniform, φ uniform
    r = s.radius_cm * sqrt(rand(rng))
    φ = 2π * rand(rng)
    z = c[3] + s.half_height_cm * (-1.0 + 2.0 * rand(rng))
    Float64[c[1] + r*cos(φ), c[2] + r*sin(φ), z]
end

function random_position_in_volume(vol::PCylShell, rng::AbstractRNG)::Vector{Float64}
    s = vol.logical.solid
    c = vol.logical.position
    R_o = R_outer(s)
    R_i = s.R_inner_cm
    # Uniform in annular region: r² uniform in [R_i², R_o²]
    r = sqrt(R_i^2 + rand(rng) * (R_o^2 - R_i^2))
    φ = 2π * rand(rng)
    z = c[3] + s.half_height_cm * (-1.0 + 2.0 * rand(rng))
    Float64[c[1] + r*cos(φ), c[2] + r*sin(φ), z]
end

function random_position_in_volume(vol::PDisk, rng::AbstractRNG)::Vector{Float64}
    # Approximate: uniform in the thin shell region
    # For thin disks, sample uniformly on inner surface then offset by random t
    s = vol.logical.solid
    c = vol.logical.position
    orient = vol.logical.orientation
    # Sample r uniform in disk
    r = s.radius_cm * sqrt(rand(rng))
    φ = 2π * rand(rng)
    # For flat disk: z offset within thickness
    if is_flat(s)
        sgn = orient === :up ? 1.0 : -1.0
        dz = sgn * rand(rng) * s.wall_thickness_cm
    else
        # Ellipsoidal: sample along the dome surface (approximate)
        # Use parametric: at radius r, z = c * sqrt(1 - (r/a)^2)
        a = s.radius_cm
        c_depth = depth(s)
        z_inner = c_depth * sqrt(max(0.0, 1.0 - (r/a)^2))
        sgn = orient === :up ? 1.0 : -1.0
        dz = sgn * (z_inner + rand(rng) * s.wall_thickness_cm)
    end
    Float64[c[1] + r*cos(φ), c[2] + r*sin(φ), c[3] + dz]
end


# =====================================================================
# cos θ computation (angle to surface normal toward LXe)
# =====================================================================

"""
    cos_theta_to_lxe(pos_exit, dir, source_vol) -> Float64

Cosine of angle between gamma direction and the inward normal of the
source surface at the exit point. Positive = toward LXe (inward),
negative = away from LXe (backward).

- **CylShell** (barrel): inward normal is -r_hat at the exit position.
  cos θ = -(dx·dir_x + dy·dir_y) / r (radial inward projection).
- **Disk** (head): inward normal is -ẑ (if :up) or +ẑ (if :down).
  cos θ = ∓dir_z.
"""
function cos_theta_to_lxe(pos::Vector{Float64}, dir::Vector{Float64},
                          source_vol::PCylShell)::Float64
    cx = source_vol.logical.position[1]
    cy = source_vol.logical.position[2]
    dx = pos[1] - cx
    dy = pos[2] - cy
    r = sqrt(dx^2 + dy^2)
    r < 1e-10 && return 0.0
    # Inward radial component of direction (negative r_hat · dir)
    -(dx * dir[1] + dy * dir[2]) / r
end

function cos_theta_to_lxe(pos::Vector{Float64}, dir::Vector{Float64},
                          source_vol::PDisk)::Float64
    orient = source_vol.logical.orientation
    # Inward = toward the LXe: -z for :up heads, +z for :down heads
    orient === :up ? -dir[3] : dir[3]
end

# General fallback (solid cylinder sources)
function cos_theta_to_lxe(pos::Vector{Float64}, dir::Vector{Float64},
                          source_vol::PhysicalVolume)::Float64
    # Assume LXe is inward radially (same as shell logic)
    cx = source_vol.logical.position[1]
    cy = source_vol.logical.position[2]
    dx = pos[1] - cx
    dy = pos[2] - cy
    r = sqrt(dx^2 + dy^2)
    r < 1e-10 && return abs(dir[3])
    -(dx * dir[1] + dy * dir[2]) / r
end


# =====================================================================
# Source flux tables
# =====================================================================

"""
    SourceFlux

Abstract type for source flux tables. A flux table stores the (E, u)
probability distribution of gammas exiting a radioactive source volume,
where u = cos θ to the inward surface normal (toward the LXe).

Concrete subtypes:
- [`SourceFluxBi214`](@ref): single gamma per decay.
- [`SourceFluxTl208`](@ref): main gamma + up to 3 independent companions.
"""
abstract type SourceFlux end


"""
    SourceFluxBi214 <: SourceFlux

Flux table for Bi-214 (U-238 late chain): one 2.448 MeV gamma per decay.

The `pdf` matrix is normalized so that `sum(pdf) = N_surviving / N_generated`.
Each bin `pdf[i,j]` gives the probability per generated decay that a gamma
exits the source with energy in bin `i` and cos θ in bin `j`.
"""
struct SourceFluxBi214 <: SourceFlux
    source_name::String
    pdf::Matrix{Float64}        # n_E × n_u
    E_min::Float64
    E_max::Float64
    n_E::Int
    n_u::Int
    N_generated::Int
    N_surviving::Int
    N_absorbed::Int
    N_backward::Int
    N_low_energy::Int
end


"""
    SourceFluxTl208 <: SourceFlux

Flux table for Tl-208 (Th-232 late chain): one 2.615 MeV main gamma
plus up to three independent companion gammas per decay.

Companion lines (583 keV at 85%, 511 keV at 23%, 861 keV at 12%) are
sampled as independent Bernoulli trials. Each companion has its own
(E, u) table and survival fraction through the source material.

At sampling time, the main gamma always fires. Each companion i fires
independently with probability `companion_BR[i] × companion_f[i]`,
producing 1–4 gammas per event.
"""
struct SourceFluxTl208 <: SourceFlux
    source_name::String

    # Main gamma table
    pdf_main::Matrix{Float64}   # n_E_main × n_u
    E_min_main::Float64
    E_max_main::Float64
    n_E_main::Int

    # Companion tables (one per line: 583, 511, 861 keV)
    pdf_companion::Vector{Matrix{Float64}}  # 3 tables, each n_E_comp[i] × n_u
    companion_E_line::Vector{Float64}       # [0.583, 0.511, 0.861]
    companion_BR::Vector{Float64}           # [0.85, 0.23, 0.12]
    companion_f::Vector{Float64}            # survival fraction per line
    E_min_companion::Vector{Float64}
    E_max_companion::Vector{Float64}
    n_E_companion::Vector{Int}

    # Shared
    n_u::Int

    # Generation statistics
    N_generated::Int
    N_surviving_main::Int
    N_absorbed_main::Int
    N_backward_main::Int
    N_low_energy_main::Int
end


"""
    SourceRateTable

Activity-weighted sum of component flux tables for one ICV surface,
in units of gammas/sec per (E, u) bin.

    pdf_rate[i,j] = sum_k ( A_k [Bq/kg] × m_k [kg] × pdf_k[i,j] )

Used for sampling the combined cryostat flux from all contributing
source layers (OCV, MLI, ICV).
"""
struct SourceRateTable
    surface::Symbol             # :barrel, :top, or :bottom
    pdf_rate::Matrix{Float64}   # n_E × n_u, gammas/sec/bin
    E_min::Float64
    E_max::Float64
    n_E::Int
    n_u::Int
    component_names::Vector{String}
    component_rates::Vector{Float64}  # total rate per component [gammas/sec]
    total_rate::Float64               # sum of component_rates [gammas/sec]
end


# =====================================================================
# Source volume loader
# =====================================================================

"""
    SourceVolumeInfo

Loaded source volume with geometry, material, activity, and transport metadata.
Built by `load_source_geometry` from the source geometry JSON.
"""
struct SourceVolumeInfo
    name::String
    volume::PhysicalVolume
    material::Material
    activity::Dict{String,Float64}  # e.g. "Bi214_mBq_per_kg" => 0.08
    mass_kg::Float64                # from geometry or equivalent_mass_kg
    transport::Symbol               # :KN or :transparent
    source_class::String            # "shell_source" or "virtual_source"
    approximation::String
end


"""
    load_source_geometry(path, materials) -> Dict{String, SourceVolumeInfo}

Load source volumes from a source geometry JSON file. Builds a
`PhysicalVolume` for each entry and computes mass from geometry
(for shell sources) or from `equivalent_mass_kg` (for virtual sources).
"""
function load_source_geometry(path::AbstractString,
                              materials::Dict{String,Material})::Dict{String,SourceVolumeInfo}
    raw = open(path, "r") do io
        JSON.parse(io)
    end

    result = Dict{String,SourceVolumeInfo}()
    for d in raw["sources"]
        name = String(d["name"])
        mat_name = String(d["material"])
        haskey(materials, mat_name) || error("Unknown material '$mat_name' for source '$name'")
        mat = materials[mat_name]

        shape = lowercase(String(d["shape"]))
        pos = Float64.(d["position_cm"])

        vol = if shape == "cylinder_shell"
            solid = CylShell(Float64(d["R_inner_cm"]),
                             Float64(d["wall_thickness_cm"]),
                             Float64(d["half_height_cm"]))
            PCylShell(name, LCylShell(solid, pos), mat)
        elseif shape == "disk"
            solid = Disk(Float64(d["radius_cm"]),
                         Float64(d["wall_thickness_cm"]),
                         Float64(get(d, "aspect_ratio", Inf)))
            orient = Symbol(get(d, "orientation", "up"))
            PDisk(name, LDisk(solid, pos, orient), mat)
        elseif shape == "cylinder"
            solid = Cyl(Float64(d["radius_cm"]),
                        Float64(d["half_height_cm"]))
            PCyl(name, LCyl(solid, pos), mat)
        else
            error("Unsupported source shape '$shape' for '$name'")
        end

        transport = Symbol(lowercase(String(d["transport_source"])))
        source_class = String(d["source_class"])
        approximation = String(get(d, "approximation", "exact"))

        eq_mass = Float64(get(d, "equivalent_mass_kg", 0.0))
        mass_kg = if eq_mass > 0.0
            eq_mass
        else
            mass(vol) / 1000.0  # g -> kg
        end

        activity = Dict{String,Float64}()
        raw_act = get(d, "activity", nothing)
        if raw_act !== nothing
            for (k, v) in raw_act
                v !== nothing && (activity[String(k)] = Float64(v))
            end
        end

        result[name] = SourceVolumeInfo(name, vol, mat, activity, mass_kg,
                                        transport, source_class, approximation)
    end
    result
end


# =====================================================================
# Single-source flux generation
# =====================================================================

"""
    generate_flux_bi214(N, source_vol, cfg, rng; ...) -> SourceFluxBi214

Generate N Bi-214 decays (single 2.448 MeV gamma) inside `source_vol`,
propagate each through the source material (KN), and tally survivors
into an (E, u) probability table.

`E_min` is the minimum energy for a gamma to be relevant for the ROI.
Gammas that degrade below `E_min` are counted as `N_low_energy`.
`E_max` is the upper edge of the binning (no rejection).
"""
function generate_flux_bi214(N::Int, source_vol::PhysicalVolume,
                              cfg::SimConfig, rng::AbstractRNG;
                              E_gamma::Float64=2.448,
                              E_min::Float64=2.200, E_max::Float64=2.500,
                              n_E::Int=25, n_u::Int=10)::SourceFluxBi214
    dE = (E_max - E_min) / n_E
    counts = zeros(Int, n_E, n_u)
    n_surviving = 0
    n_absorbed = 0
    n_backward = 0
    n_low_energy = 0

    for _ in 1:N
        pos = random_position_in_volume(source_vol, rng)
        dir = sample_isotropic_direction(rng)

        status, E_out, pos_out, dir_out = propagate_in_source(
            E_gamma, pos, dir, source_vol, cfg, rng)

        if status === :absorbed
            n_absorbed += 1
            continue
        end

        u = cos_theta_to_lxe(pos_out, dir_out, source_vol)
        if u <= 0.0
            n_backward += 1
            continue
        end

        if E_out < E_min
            n_low_energy += 1
            continue
        end

        i_E = clamp(floor(Int, (E_out - E_min) / dE) + 1, 1, n_E)
        i_u = clamp(floor(Int, u * n_u) + 1, 1, n_u)
        counts[i_E, i_u] += 1
        n_surviving += 1
    end

    pdf = counts ./ N
    SourceFluxBi214(source_vol.name, pdf, E_min, E_max, n_E, n_u,
                    N, n_surviving, n_absorbed, n_backward, n_low_energy)
end


"""
    generate_flux_tl208(N, source_vol, cfg, rng; ...) -> SourceFluxTl208

Generate N Tl-208 decays inside `source_vol`. Each decay produces a
main 2.615 MeV gamma plus up to 3 companions via independent Bernoulli
trials (583 keV at 85%, 511 keV at 23%, 861 keV at 12%).

Each gamma is propagated independently through the source material.
The main gamma fills `pdf_main`; each companion line fills its own
`pdf_companion[i]`. Survival fractions are computed per line.

`E_min_main` is the minimum energy for the main gamma to be relevant
for the ROI. Companion energy windows are set independently.
"""
function generate_flux_tl208(N::Int, source_vol::PhysicalVolume,
                              cfg::SimConfig, rng::AbstractRNG;
                              E_main::Float64=2.615,
                              E_min_main::Float64=2.200, E_max_main::Float64=2.620,
                              n_E_main::Int=25,
                              companion_lines::Vector{Float64}=Float64[0.583, 0.511, 0.861],
                              companion_BRs::Vector{Float64}=Float64[0.85, 0.23, 0.12],
                              E_min_comp::Vector{Float64}=Float64[0.400, 0.350, 0.700],
                              E_max_comp::Vector{Float64}=Float64[0.650, 0.600, 0.950],
                              n_E_comp::Vector{Int}=Int[25, 25, 25],
                              n_u::Int=10)::SourceFluxTl208
    n_comp = length(companion_lines)
    dE_main = (E_max_main - E_min_main) / n_E_main
    dE_c = [(E_max_comp[i] - E_min_comp[i]) / n_E_comp[i] for i in 1:n_comp]

    counts_main = zeros(Int, n_E_main, n_u)
    counts_comp = [zeros(Int, n_E_comp[i], n_u) for i in 1:n_comp]

    n_surviving_main = 0
    n_absorbed_main = 0
    n_backward_main = 0
    n_low_energy_main = 0

    n_fired_comp = zeros(Int, n_comp)
    n_surviving_comp = zeros(Int, n_comp)

    for _ in 1:N
        origin = random_position_in_volume(source_vol, rng)

        # --- Main gamma (always) ---
        dir_main = sample_isotropic_direction(rng)
        status, E_out, pos_out, dir_out = propagate_in_source(
            E_main, copy(origin), dir_main, source_vol, cfg, rng)

        if status === :absorbed
            n_absorbed_main += 1
        else
            u = cos_theta_to_lxe(pos_out, dir_out, source_vol)
            if u <= 0.0
                n_backward_main += 1
            elseif E_out < E_min_main
                n_low_energy_main += 1
            else
                i_E = clamp(floor(Int, (E_out - E_min_main) / dE_main) + 1, 1, n_E_main)
                i_u = clamp(floor(Int, u * n_u) + 1, 1, n_u)
                counts_main[i_E, i_u] += 1
                n_surviving_main += 1
            end
        end

        # --- Companions (independent Bernoulli trials) ---
        for ic in 1:n_comp
            rand(rng) >= companion_BRs[ic] && continue
            n_fired_comp[ic] += 1

            dir_c = sample_isotropic_direction(rng)
            status_c, E_c, pos_c, dir_c_out = propagate_in_source(
                companion_lines[ic], copy(origin), dir_c, source_vol, cfg, rng)

            status_c === :absorbed && continue

            u_c = cos_theta_to_lxe(pos_c, dir_c_out, source_vol)
            u_c <= 0.0 && continue

            E_c < E_min_comp[ic] && continue

            i_E = clamp(floor(Int, (E_c - E_min_comp[ic]) / dE_c[ic]) + 1, 1, n_E_comp[ic])
            i_u = clamp(floor(Int, u_c * n_u) + 1, 1, n_u)
            counts_comp[ic][i_E, i_u] += 1
            n_surviving_comp[ic] += 1
        end
    end

    pdf_main = counts_main ./ N
    pdf_comp = [counts_comp[i] ./ N for i in 1:n_comp]
    comp_f = [n_fired_comp[i] > 0 ? n_surviving_comp[i] / n_fired_comp[i] : 0.0
              for i in 1:n_comp]

    SourceFluxTl208(
        source_vol.name,
        pdf_main, E_min_main, E_max_main, n_E_main,
        pdf_comp,
        companion_lines,
        companion_BRs,
        comp_f,
        E_min_comp,
        E_max_comp,
        n_E_comp,
        n_u,
        N, n_surviving_main, n_absorbed_main, n_backward_main, n_low_energy_main
    )
end


# =====================================================================
# Compound propagation through multiple layers
# =====================================================================

"""
    propagate_through_layers(E, pos, dir, layers, cfg, rng)
        -> (status, E_out, pos_out, dir_out)

Propagate a gamma through a sequence of material layers. Between layers,
the gamma travels in a straight line (vacuum gap) to the entry of the
next volume via `distance_to_entry` on the layer's logical volume.

Returns `(:exited, E, pos, dir)`, `(:absorbed, 0, pos, dir)`, or
`(:lost, E, pos, dir)` if the ray misses a layer.
"""
function propagate_through_layers(E::Float64, pos::Vector{Float64}, dir::Vector{Float64},
                                   layers::Vector{<:PhysicalVolume},
                                   cfg::SimConfig, rng::AbstractRNG
                                   )::Tuple{Symbol,Float64,Vector{Float64},Vector{Float64}}
    for layer in layers
        t_entry = distance_to_entry(pos, dir, layer.logical)
        if !isfinite(t_entry)
            return (:lost, E, pos, dir)
        end
        pos .= pos .+ dir .* (t_entry + 1e-4)

        status, E, pos, dir = propagate_in_source(E, pos, dir, layer, cfg, rng)
        if status === :absorbed
            return (:absorbed, 0.0, pos, dir)
        end
    end
    (:exited, E, pos, dir)
end


"""
    generate_flux_compound_bi214(N, source_vol, layers, exit_vol, cfg, rng; ...)
        -> SourceFluxBi214

Generate N Bi-214 decays in `source_vol`, propagate each gamma through
the source, then through the sequence of `layers`. `exit_vol` is used
for `cos_theta_to_lxe` at the final exit surface.
"""
function generate_flux_compound_bi214(N::Int, source_vol::PhysicalVolume,
                                       layers::Vector{<:PhysicalVolume},
                                       exit_vol::PhysicalVolume,
                                       cfg::SimConfig, rng::AbstractRNG;
                                       E_gamma::Float64=2.448,
                                       E_min::Float64=2.200, E_max::Float64=2.500,
                                       n_E::Int=25, n_u::Int=10)::SourceFluxBi214
    dE = (E_max - E_min) / n_E
    counts = zeros(Int, n_E, n_u)
    n_surviving = 0
    n_absorbed = 0
    n_backward = 0
    n_low_energy = 0

    for _ in 1:N
        pos = random_position_in_volume(source_vol, rng)
        dir = sample_isotropic_direction(rng)

        status, E_out, pos_out, dir_out = propagate_in_source(
            E_gamma, pos, dir, source_vol, cfg, rng)

        if status === :absorbed
            n_absorbed += 1
            continue
        end

        u_source = cos_theta_to_lxe(pos_out, dir_out, source_vol)
        if u_source <= 0.0
            n_backward += 1
            continue
        end

        status2, E_out2, pos_out2, dir_out2 = propagate_through_layers(
            E_out, pos_out, dir_out, layers, cfg, rng)

        if status2 === :absorbed || status2 === :lost
            n_absorbed += 1
            continue
        end

        u = cos_theta_to_lxe(pos_out2, dir_out2, exit_vol)
        if u <= 0.0
            n_backward += 1
            continue
        end

        if E_out2 < E_min
            n_low_energy += 1
            continue
        end

        i_E = clamp(floor(Int, (E_out2 - E_min) / dE) + 1, 1, n_E)
        i_u = clamp(floor(Int, u * n_u) + 1, 1, n_u)
        counts[i_E, i_u] += 1
        n_surviving += 1
    end

    pdf = counts ./ N
    SourceFluxBi214(source_vol.name, pdf, E_min, E_max, n_E, n_u,
                    N, n_surviving, n_absorbed, n_backward, n_low_energy)
end


"""
    generate_flux_compound_tl208(N, source_vol, layers, exit_vol, cfg, rng; ...)
        -> SourceFluxTl208

Generate N Tl-208 decays in `source_vol`, propagate main + companions
through the source and then through `layers`. Tally survivors at the
exit of the last layer.
"""
function generate_flux_compound_tl208(N::Int, source_vol::PhysicalVolume,
                                       layers::Vector{<:PhysicalVolume},
                                       exit_vol::PhysicalVolume,
                                       cfg::SimConfig, rng::AbstractRNG;
                                       E_main::Float64=2.615,
                                       E_min_main::Float64=2.200, E_max_main::Float64=2.620,
                                       n_E_main::Int=25,
                                       companion_lines::Vector{Float64}=Float64[0.583, 0.511, 0.861],
                                       companion_BRs::Vector{Float64}=Float64[0.85, 0.23, 0.12],
                                       E_min_comp::Vector{Float64}=Float64[0.400, 0.350, 0.700],
                                       E_max_comp::Vector{Float64}=Float64[0.650, 0.600, 0.950],
                                       n_E_comp::Vector{Int}=Int[25, 25, 25],
                                       n_u::Int=10)::SourceFluxTl208
    n_comp = length(companion_lines)
    dE_main = (E_max_main - E_min_main) / n_E_main
    dE_c = [(E_max_comp[i] - E_min_comp[i]) / n_E_comp[i] for i in 1:n_comp]

    counts_main = zeros(Int, n_E_main, n_u)
    counts_comp = [zeros(Int, n_E_comp[i], n_u) for i in 1:n_comp]

    n_surviving_main = 0
    n_absorbed_main = 0
    n_backward_main = 0
    n_low_energy_main = 0

    n_fired_comp = zeros(Int, n_comp)
    n_surviving_comp = zeros(Int, n_comp)

    for _ in 1:N
        origin = random_position_in_volume(source_vol, rng)

        # --- Main gamma ---
        dir_main = sample_isotropic_direction(rng)
        status, E_out, pos_out, dir_out = propagate_in_source(
            E_main, copy(origin), dir_main, source_vol, cfg, rng)

        if status === :absorbed
            n_absorbed_main += 1
        else
            u_source = cos_theta_to_lxe(pos_out, dir_out, source_vol)
            if u_source <= 0.0
                n_backward_main += 1
            else
                status2, E_out2, pos_out2, dir_out2 = propagate_through_layers(
                    E_out, pos_out, dir_out, layers, cfg, rng)

                if status2 === :absorbed || status2 === :lost
                    n_absorbed_main += 1
                else
                    u = cos_theta_to_lxe(pos_out2, dir_out2, exit_vol)
                    if u <= 0.0
                        n_backward_main += 1
                    elseif E_out2 < E_min_main
                        n_low_energy_main += 1
                    else
                        i_E = clamp(floor(Int, (E_out2 - E_min_main) / dE_main) + 1, 1, n_E_main)
                        i_u = clamp(floor(Int, u * n_u) + 1, 1, n_u)
                        counts_main[i_E, i_u] += 1
                        n_surviving_main += 1
                    end
                end
            end
        end

        # --- Companions ---
        for ic in 1:n_comp
            rand(rng) >= companion_BRs[ic] && continue
            n_fired_comp[ic] += 1

            dir_c = sample_isotropic_direction(rng)
            status_c, E_c, pos_c, dir_c_out = propagate_in_source(
                companion_lines[ic], copy(origin), dir_c, source_vol, cfg, rng)

            status_c === :absorbed && continue

            u_source_c = cos_theta_to_lxe(pos_c, dir_c_out, source_vol)
            u_source_c <= 0.0 && continue

            status2_c, E_c2, pos_c2, dir_c2 = propagate_through_layers(
                E_c, pos_c, dir_c_out, layers, cfg, rng)

            (status2_c === :absorbed || status2_c === :lost) && continue

            u_c = cos_theta_to_lxe(pos_c2, dir_c2, exit_vol)
            u_c <= 0.0 && continue

            E_c2 < E_min_comp[ic] && continue

            i_E = clamp(floor(Int, (E_c2 - E_min_comp[ic]) / dE_c[ic]) + 1, 1, n_E_comp[ic])
            i_u = clamp(floor(Int, u_c * n_u) + 1, 1, n_u)
            counts_comp[ic][i_E, i_u] += 1
            n_surviving_comp[ic] += 1
        end
    end

    pdf_main = counts_main ./ N
    pdf_comp = [counts_comp[i] ./ N for i in 1:n_comp]
    comp_f = [n_fired_comp[i] > 0 ? n_surviving_comp[i] / n_fired_comp[i] : 0.0
              for i in 1:n_comp]

    SourceFluxTl208(
        source_vol.name,
        pdf_main, E_min_main, E_max_main, n_E_main,
        pdf_comp,
        companion_lines,
        companion_BRs,
        comp_f,
        E_min_comp,
        E_max_comp,
        n_E_comp,
        n_u,
        N, n_surviving_main, n_absorbed_main, n_backward_main, n_low_energy_main
    )
end
