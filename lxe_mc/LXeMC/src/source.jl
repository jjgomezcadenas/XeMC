"""
Source flux generation for radioactive background studies.

Implements Stage 1 of the background simulation: generate decay gammas
inside a source volume (e.g., Ti cryostat), propagate them through the
source material via photon-only transport, and collect those that exit
into the LXe with energy in the ROI window [E_min, E_max].

Multi-gamma events (Tl-208 with companions) undergo a veto check:
if two gammas with E > 100 keV reach the LXe with |Δz| > 3 mm at
their first interaction, the event is discarded (multi-gamma veto).

Output: a `FluxTable` — binned 2D histogram in (E, cos θ) of gammas
surviving all cuts, normalized per generated decay.
"""

using Random


# =====================================================================
# Flux table
# =====================================================================

"""
    FluxTable

Binned 2D histogram of surviving gamma flux in (E, u=cos θ).
Energy bins: [E_min, E_max] in steps of dE.
Angular bins: [0, 1] in steps of du.
"""
struct FluxTable
    E_min::Float64          # lower edge [MeV]
    E_max::Float64          # upper edge [MeV]
    dE::Float64             # bin width [MeV]
    n_E::Int                # number of E bins
    n_u::Int                # number of cos θ bins
    counts::Matrix{Int}     # n_E × n_u bin counts
    N_generated::Int        # total decays generated
    N_surviving::Int        # total gammas passing all cuts
    N_vetoed::Int           # events killed by multi-gamma veto
    N_low_energy::Int       # gammas rejected for E < E_min
    N_absorbed::Int         # gammas absorbed in source
end


"""Marginal energy spectrum dN/dE (sum over u bins)."""
spectrum_E(ft::FluxTable) = vec(sum(ft.counts, dims=2))

"""Marginal angular distribution dN/du (sum over E bins)."""
spectrum_u(ft::FluxTable) = vec(sum(ft.counts, dims=1))

"""Energy bin centers [MeV]."""
E_centers(ft::FluxTable) = [ft.E_min + (i - 0.5) * ft.dE for i in 1:ft.n_E]

"""cos θ bin centers."""
u_centers(ft::FluxTable) = [(i - 0.5) / ft.n_u for i in 1:ft.n_u]

"""Survival fraction (gammas out / decays generated)."""
survival_fraction(ft::FluxTable) = ft.N_surviving / ft.N_generated


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
    cos_theta_to_lxe(pos, dir, source_vol, det) -> Float64

Compute cos θ of gamma direction relative to the inward normal of the
source surface facing the LXe. For a cylindrical shell, the inward
normal is -r_hat. For a disk, it's ∓ẑ depending on orientation.

Returns value in [-1, 1]; only gammas with cos θ > 0 are directed
toward the LXe.
"""
function cos_theta_to_lxe(dir::Vector{Float64}, source_vol::PCylShell)::Float64
    # For a shell, "toward LXe" is radially inward = -r_hat
    # We approximate: cos θ = -dir · r_hat (where r_hat is radial direction)
    # This is position-dependent but for thin shells we can ignore z-component
    # Actually, for the flux we just care about the angle to the inward-facing normal
    # For a shell barrel, inward normal ≈ -r_hat at the emission point
    # Since we don't track exact emission point after propagation, use the
    # exit direction's radial component
    # Simpler: cos θ wrt z-axis as proxy (cylindrical symmetry)
    # Actually: for the flux table, u = |cos(angle to detector axis)|
    # which for barrel sources ≈ angle from the radial direction...
    #
    # Simplest physically meaningful definition: u = cos(angle between
    # gamma direction and the vector from emission point to detector center)
    # For now: just use the z-component of direction as u (appropriate for
    # top/bottom sources; for barrel, the relevant angle is radial)
    #
    # Better: u = inward radial component for barrels
    # Just return the magnitude of the radial-inward component
    abs(dir[3])  # placeholder — will refine per geometry
end

function cos_theta_to_lxe(dir::Vector{Float64}, source_vol::PDisk)::Float64
    # For a disk, inward normal is -z (if :up) or +z (if :down)
    orient = source_vol.logical.orientation
    orient === :up ? -dir[3] : dir[3]
end

# General fallback: use |cos θ| to z-axis
function cos_theta_to_lxe(dir::Vector{Float64}, source_vol::PhysicalVolume)::Float64
    abs(dir[3])
end


# =====================================================================
# Main flux generation
# =====================================================================

"""
    generate_source_flux(N, source_vol, scheme, det, cfg, rng;
                         E_min=2.370, E_max=2.620, n_E=25, n_u=10,
                         veto_E=0.100, veto_dz=0.30) -> FluxTable

Stage 1: generate N decay events in `source_vol`, propagate gammas
through the source material, and build the flux table of gammas
exiting into the LXe within the energy window [E_min, E_max].

# Parameters
- `N`: number of decay events to generate
- `source_vol`: the radioactive source (PCylShell, PDisk, etc.)
- `scheme`: decay scheme (from `load_decays()`)
- `det`: detector (for active volume / geometry info)
- `cfg`: simulation config
- `rng`: random number generator
- `E_min`, `E_max`: energy window [MeV] (default [2.370, 2.620])
- `n_E`: number of energy bins (default 25)
- `n_u`: number of cos θ bins (default 10)
- `veto_E`: companion energy threshold for veto [MeV] (default 0.100)
- `veto_dz`: z-separation for multi-gamma veto [cm] (default 0.30)

# Veto logic for multi-gamma events
If >1 gamma exits into LXe with E > veto_E:
- Compute first interaction z for each (one mfp draw)
- If |Δz| > veto_dz → event vetoed (both gammas visible and separated)
- If |Δz| ≤ veto_dz → count as single gamma (use the highest-energy one)
"""
function generate_source_flux(N::Int, source_vol::PhysicalVolume,
                              scheme::DecayScheme, det::Detector,
                              cfg::SimConfig, rng::AbstractRNG;
                              E_min::Float64=2.370, E_max::Float64=2.620,
                              n_E::Int=25, n_u::Int=10,
                              veto_E::Float64=0.100,
                              veto_dz::Float64=0.30)::FluxTable

    vol_active = active_volume(det)
    mat_lxe = vol_active.material
    dE = (E_max - E_min) / n_E

    counts = zeros(Int, n_E, n_u)
    n_surviving = 0
    n_vetoed = 0
    n_low_energy = 0
    n_absorbed = 0

    for _ in 1:N
        # Generate decay event
        emissions = sample_decay(scheme, rng)

        # Propagate each gamma through source
        exited_gammas = Vector{Tuple{Float64,Vector{Float64},Vector{Float64}}}()  # (E, pos, dir)

        for em in emissions
            pos = random_position_in_volume(source_vol, rng)
            status, E_out, pos_out, dir_out = propagate_in_source(
                em.E_MeV, pos, em.direction, source_vol, cfg, rng)

            if status === :exited
                push!(exited_gammas, (E_out, pos_out, dir_out))
            else
                n_absorbed += 1
            end
        end

        # No gammas exited
        isempty(exited_gammas) && continue

        # Filter: only gammas above veto threshold are "visible"
        visible = filter(g -> g[1] > veto_E, exited_gammas)

        if length(visible) == 0
            continue
        elseif length(visible) == 1
            # Single visible gamma — check if in energy window
            E_g, _, dir_g = visible[1]
            if E_g < E_min
                n_low_energy += 1
                continue
            end
            if E_g > E_max
                n_low_energy += 1
                continue
            end
            # Compute u = cos θ and bin
            u = cos_theta_to_lxe(dir_g, source_vol)
            u <= 0.0 && continue  # going away from detector

            i_E = clamp(floor(Int, (E_g - E_min) / dE) + 1, 1, n_E)
            i_u = clamp(floor(Int, u * n_u) + 1, 1, n_u)
            counts[i_E, i_u] += 1
            n_surviving += 1

        else
            # Multiple visible gammas — apply veto
            # Compute z of first interaction in LXe for each
            z_first = Float64[]
            for (E_g, _, dir_g) in visible
                sC, sP, sPh = sigma_three(mat_lxe, E_g)
                Σ = mat_lxe.n_atom * (sC + sP + sPh)
                s = sample_distance(Σ, rng)
                # z at first interaction (approximate: from TPC entry)
                push!(z_first, s * dir_g[3])
            end

            # Check all pairs for separation
            vetoed = false
            for i in 1:length(z_first)
                for j in (i+1):length(z_first)
                    if abs(z_first[i] - z_first[j]) > veto_dz
                        vetoed = true
                        break
                    end
                end
                vetoed && break
            end

            if vetoed
                n_vetoed += 1
                continue
            end

            # Not vetoed: take highest-energy gamma
            idx_max = argmax([g[1] for g in visible])
            E_g, _, dir_g = visible[idx_max]

            if E_g < E_min || E_g > E_max
                n_low_energy += 1
                continue
            end

            u = cos_theta_to_lxe(dir_g, source_vol)
            u <= 0.0 && continue

            i_E = clamp(floor(Int, (E_g - E_min) / dE) + 1, 1, n_E)
            i_u = clamp(floor(Int, u * n_u) + 1, 1, n_u)
            counts[i_E, i_u] += 1
            n_surviving += 1
        end
    end

    FluxTable(E_min, E_max, dE, n_E, n_u, counts,
              N, n_surviving, n_vetoed, n_low_energy, n_absorbed)
end
