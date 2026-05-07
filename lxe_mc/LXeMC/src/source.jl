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

