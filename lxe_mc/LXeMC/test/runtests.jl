using Test
using LXeMC
using Random
using Statistics

const CFG  = default_config()
const MATS = load_materials(CFG)
const DET  = load_detector(default_detector_path(), MATS)
const VOL  = active_volume(DET)
const MAT  = VOL.material
# TPC center position (cathode at z=0, TPC spans [0, 145.6])
const TPC_CENTER = (0.0, 0.0, VOL.logical.position[3])


# =====================================================================
# Test 1: Config loading
# =====================================================================
@testset "Config" begin
    @test CFG.Egamma_cut > 0.0
    @test CFG.Te_cut > 0.0
    @test CFG.k_min > 0.0
    @test CFG.ds_step > 0.0
    @test CFG.generation_cap > 0
    @test CFG.dz_resolution > 0.0
    @test CFG.E_cluster_min >= 0.0
    @test CFG.me > 0.0
    @test CFG.re > 0.0
    @test CFG.alpha_fs > 0.0
    @test CFG.N_A > 0.0
end


# =====================================================================
# Test 2: Materials loading
# =====================================================================
@testset "Materials" begin
    @test haskey(MATS, "LXe")
    @test haskey(MATS, "Vacuum")
    @test MAT.active == true
    @test MAT.tracking == true
    @test MAT.density > 0.0
    @test MAT.Z_eff > 0
    @test MAT.n_atom > 0.0
    @test MAT.xcom !== nothing
    @test MAT.estar !== nothing
    @test MAT.brems_T !== nothing
    # Vacuum has zero density
    @test MATS["Vacuum"].density == 0.0
end


# =====================================================================
# Test 3: Detector geometry
# =====================================================================
@testset "Detector geometry" begin
    @test DET.name == "LZ"
    @test length(DET.volumes) == 11  # TPC + FC + Skin + RFR + Dome + 3 OCV + 3 ICV (FV separate)
    @test VOL.name == "LXeTPC"

    # TPC: centered at z=72.8 (cathode at z=0, drift upward)
    @test VOL.logical.solid.radius_cm ≈ 72.8
    @test VOL.logical.solid.half_height_cm ≈ 72.8
    @test is_inside(VOL, [0.0, 0.0, 72.8])    # center of TPC
    @test !is_inside(VOL, [0.0, 0.0, -1.0])   # below cathode
    @test !is_inside(VOL, [100.0, 0.0, 72.8])  # outside radius

    by_name = Dict(v.name => v for v in DET.volumes)

    # LXe masses vs reference table (tonnes)
    m_tpc  = mass(VOL) / 1e6                               # tonnes
    m_skin = mass(by_name["Skin"]) / 1e6
    m_rfr  = mass(by_name["RFR"]) / 1e6
    m_dome = mass(by_name["Dome"]) / 1e6
    m_lxe  = m_tpc + m_skin + m_rfr + m_dome
    @test m_tpc  ≈ 7.16  rtol=0.01    # ref 7.16 t
    @test m_skin ≈ 1.80  rtol=0.01    # ref 1.80 t
    @test m_rfr  ≈ 0.68  rtol=0.01    # ref 0.68 t
    @test m_dome ≈ 2.88  rtol=0.25    # cylindrical approx overestimates (ignores head curvature)
    @test m_lxe  ≈ 12.52 rtol=0.06   # total LXe, dome approx inflates by ~0.6 t

    # Ti cryostat masses (within 1%)
    m_ocv = sum(mass(by_name[n]) for n in ["OCV_barrel","OCV_top","OCV_bottom"]) / 1000.0
    m_icv = sum(mass(by_name[n]) for n in ["ICV_barrel","ICV_top","ICV_bottom"]) / 1000.0
    @test m_ocv ≈ 778.6  rtol=0.01
    @test m_icv ≈ 651.1  rtol=0.01

    # Field cage is PTFE
    @test by_name["FieldCage"].material.name == "PTFE"

    # Fiducial volume
    fv = fiducial_volume(DET)
    @test fv.name == "FV"
    @test fv.logical.solid.radius_cm ≈ 39.0
    @test fv.logical.solid.half_height_cm ≈ 35.0
    @test is_inside(fv, [0.0, 0.0, 61.0])      # FV center
    @test !is_inside(fv, [0.0, 0.0, 0.0])       # cathode, outside FV
    @test !is_inside(fv, [40.0, 0.0, 61.0])     # outside FV radius
    @test is_inside(VOL, fv.logical.position)    # FV center is inside TPC
    m_fv = mass(fv) / 1e6                        # tonnes
    @test 0.9 < m_fv < 1.1                       # ~1 tonne fiducial
end


# =====================================================================
# Test 3b: Geometric solids and CylShell
# =====================================================================
@testset "Geometry solids" begin
    # Solid cylinder
    c = Cyl(10.0, 20.0)
    @test volume(c) ≈ π * 100.0 * 40.0

    # CylShell: R_inner=50, wall=2, half_height=30
    cs = CylShell(50.0, 2.0, 30.0)
    @test R_outer(cs) ≈ 52.0
    @test volume(cs) ≈ π * (52.0^2 - 50.0^2) * 60.0
    @test volume_inner(cs) ≈ π * 50.0^2 * 60.0

    # LCylShell placement and is_inside
    lcs = LCylShell(cs, [0.0, 0.0, 0.0])
    @test is_inside(lcs, [51.0, 0.0, 0.0])   # in shell
    @test !is_inside(lcs, [49.0, 0.0, 0.0])  # inside bore
    @test !is_inside(lcs, [53.0, 0.0, 0.0])  # outside
    @test !is_inside(lcs, [51.0, 0.0, 31.0]) # above height

    # Offset placement
    lcs2 = LCylShell(cs, [0.0, 0.0, 100.0])
    @test is_inside(lcs2, [51.0, 0.0, 100.0])
    @test !is_inside(lcs2, [51.0, 0.0, 0.0])

    # Flat disk: R=50, t=1, aspect=Inf
    df = Disk(50.0, 1.0, Inf)
    @test is_flat(df)
    @test depth(df) ≈ 0.0
    @test surface_area_inner(df) ≈ π * 50.0^2
    @test volume(df) ≈ π * 50.0^2 * 1.0

    # Hemisphere: R=50, t=1, aspect=1
    dh = Disk(50.0, 1.0, 1.0)
    @test !is_flat(dh)
    @test depth(dh) ≈ 50.0
    @test surface_area_inner(dh) ≈ 2π * 50.0^2

    # 2:1 ellipsoidal: R=50, t=1, aspect=2
    de = Disk(50.0, 1.0, 2.0)
    @test depth(de) ≈ 25.0
    # Area should be between flat (πR²) and hemisphere (2πR²)
    @test π * 50.0^2 < surface_area_inner(de) < 2π * 50.0^2

    # LDisk is_inside: flat disk pointing up at origin
    ldf = LDisk(Disk(50.0, 2.0, Inf), [0.0, 0.0, 0.0], :up)
    @test is_inside(ldf, [10.0, 0.0, 1.0])    # inside disc
    @test !is_inside(ldf, [10.0, 0.0, -1.0])  # below equator
    @test !is_inside(ldf, [60.0, 0.0, 1.0])   # outside radius

    # LDisk is_inside: 2:1 ellipsoidal pointing up
    lde = LDisk(Disk(50.0, 2.0, 2.0), [0.0, 0.0, 0.0], :up)
    @test is_inside(lde, [0.0, 0.0, 25.5])     # near apex, in shell
    @test !is_inside(lde, [0.0, 0.0, 12.0])    # inside inner ellipsoid
    @test !is_inside(lde, [0.0, 0.0, 28.0])    # outside outer ellipsoid

    # Validation
    @test_throws ErrorException Cyl(-1.0, 5.0)
    @test_throws ErrorException CylShell(10.0, -1.0, 5.0)
    @test_throws ErrorException Box(0.0, 5.0, 5.0)
    @test_throws ErrorException Disk(-1.0, 1.0, 2.0)
end


# =====================================================================
# Test 4: XCOM total = sum of channels
# =====================================================================
@testset "XCOM total = sum of channels" begin
    for E in [0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
        sC  = sigma_compton(MAT, E; per_atom=false)
        sP  = sigma_pair(MAT, E; per_atom=false)
        sPh = sigma_phot(MAT, E; per_atom=false)
        sT  = sigma_total(MAT, E; per_atom=false)
        @test (sC + sP + sPh) ≈ sT  rtol=0.025
    end
end


# =====================================================================
# Test 5: NIST interpolation at grid points
# =====================================================================
@testset "NIST interpolation at grid points" begin
    @test sigma_compton(MAT, 1.0; per_atom=false) ≈ 5.211e-2  rtol=0.01
    @test dEdx_collision(MAT, 1.0) ≈ 1.122  rtol=0.01
end


# =====================================================================
# Test 6: Coulomb correction
# =====================================================================
@testset "Coulomb correction" begin
    fc = coulomb_correction_fc(54, CFG.alpha_fs)
    @test fc ≈ 0.16493  rtol=0.01
end


# =====================================================================
# Test 7: Bremsstrahlung dσ/dk
# =====================================================================
@testset "Bremsstrahlung dsigma/dk" begin
    @test dsigma_dk_brems(0.1, 1.0, 54, CFG) > 0.0
    @test dsigma_dk_brems(0.5, 1.0, 54, CFG) > 0.0
    @test dsigma_dk_brems(0.0, 1.0, 54, CFG) == 0.0
    @test dsigma_dk_brems(1.0, 1.0, 54, CFG) == 0.0
    @test dsigma_dk_brems(1.5, 1.0, 54, CFG) == 0.0
end


# =====================================================================
# Test 8: Brems table vs direct integration
# =====================================================================
@testset "Brems table vs direct integration" begin
    for T in [0.5, 1.0, 2.0, 5.0]
        T <= CFG.k_min && continue
        σ_table  = sigma_brems(MAT, T)
        σ_direct = LXeMC.sigma_brems_above_kmin(T, CFG.k_min, MAT.Z_eff, CFG)
        @test σ_table ≈ σ_direct  rtol=0.02
    end
end


# =====================================================================
# Test 9: Compton edge
# =====================================================================
@testset "Compton edge" begin
    rng = MersenneTwister(0)
    E0 = 2.615
    a = E0 / CFG.me
    Tmax_theory = E0 * 2a / (1 + 2a)

    N = 50_000
    Te_max = 0.0
    for _ in 1:N
        Egp, _ = sample_compton(E0, CFG, rng)
        Te = E0 - Egp
        Te_max = max(Te_max, Te)
    end
    @test Te_max < Tmax_theory + 1e-4
    @test Te_max > 0.9 * Tmax_theory
end


# =====================================================================
# Test 10: Pair energy split
# =====================================================================
@testset "Pair energy split" begin
    rng = MersenneTwister(1)
    E0 = 2.615
    eps_min = CFG.me / E0
    for _ in 1:1000
        eps = sample_pair(E0, CFG, rng)
        @test eps >= eps_min - 1e-10
        @test eps <= 1.0 - eps_min + 1e-10
    end
end


# =====================================================================
# Test 11: Brems photon energy
# =====================================================================
@testset "Brems photon energy" begin
    rng = MersenneTwister(2)
    T0 = 2.0
    for _ in 1:100
        M = brems_rejection_M(MAT, T0)
        k = sample_brems(T0, CFG.k_min, MAT.Z_eff, M, CFG, rng)
        @test k !== nothing
        @test k >= CFG.k_min
        @test k < T0
    end
end


# =====================================================================
# Test 12: Energy conservation (full mode)
# =====================================================================
@testset "Energy conservation" begin
    rng = MersenneTwister(3)
    E0 = 2.615
    N = 200
    losses = Float64[]
    for _ in 1:N
        deps = simulate_event(E0, VOL, CFG; position=TPC_CENTER, rng=rng)
        E_dep = sum(d.energy for d in deps)
        push!(losses, E0 - E_dep)
    end
    @test abs(mean(losses)) < 5 * CFG.Te_cut
end


# =====================================================================
# Test 13: SS fraction for Tl-208
# =====================================================================
@testset "SS fraction Tl-208" begin
    rng = MersenneTwister(7)
    N = 200
    n_ss = 0
    for _ in 1:N
        deps = simulate_event(2.615, VOL, CFG; position=TPC_CENTER, rng=rng)
        if is_single_site(deps, CFG.dz_resolution; E_min=CFG.E_cluster_min)
            n_ss += 1
        end
    end
    frac = n_ss / N
    @test 0.05 < frac < 0.30
end


# =====================================================================
# Test 14: Photon-only energy conservation
# =====================================================================
@testset "Photon-only energy conservation" begin
    rng = MersenneTwister(10)
    E0 = 2.615
    N = 200
    for _ in 1:N
        deps = simulate_event_photon_only(E0, VOL, CFG; position=TPC_CENTER, rng=rng)
        E_dep = sum(d.energy for d in deps)
        @test E_dep ≈ E0  rtol=1e-10
    end
end
