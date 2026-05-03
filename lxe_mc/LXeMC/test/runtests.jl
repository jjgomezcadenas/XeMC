using Test
using LXeMC
using Random
using Statistics

const CFG = default_config()
const ND  = load_nist_data(CFG)


# =====================================================================
# Test 1: Config loading
# =====================================================================
@testset "Config" begin
    @test CFG.Z == 54
    @test CFG.A ≈ 131.293
    @test CFG.rho_LXe ≈ 2.953
    @test CFG.EK ≈ 0.034561
    @test CFG.me ≈ 0.5109989461
    @test CFG.Egamma_cut ≈ 0.010
    @test CFG.Te_cut ≈ 0.050
    @test CFG.k_min ≈ 0.050
    @test CFG.generation_cap == 100
    @test CFG.dz_resolution ≈ 0.30
    # Derived quantity
    @test CFG.n_atom ≈ CFG.N_A * CFG.rho_LXe / CFG.A  rtol=1e-10
end


# =====================================================================
# Test 2: NIST XCOM total equals sum of channels
# =====================================================================
@testset "XCOM total = sum of channels" begin
    for E in [0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
        sC  = sigma_compton_NIST(ND, E; per_atom=false)
        sP  = sigma_pair_NIST(ND, E; per_atom=false)
        sPh = sigma_phot_NIST(ND, E; per_atom=false)
        sT  = sigma_total_NIST(ND, E; per_atom=false)
        @test (sC + sP + sPh) ≈ sT  rtol=0.025
    end
end


# =====================================================================
# Test 3: NIST interpolation at tabulated points
# =====================================================================
@testset "NIST interpolation at grid points" begin
    # Compton at 1.0 MeV: XCOM says 5.211e-2 cm²/g
    @test sigma_compton_NIST(ND, 1.0; per_atom=false) ≈ 5.211e-2  rtol=0.01

    # Collision stopping power at 1.0 MeV: ESTAR says 1.122 MeV cm²/g
    @test dEdx_collision_NIST(ND, 1.0) ≈ 1.122  rtol=0.01
end


# =====================================================================
# Test 4: Coulomb correction for Xe
# =====================================================================
@testset "Coulomb correction" begin
    fc = coulomb_correction_fc(54, CFG.alpha_fs)
    @test fc ≈ 0.16493  rtol=0.01  # verified against Python physics.py
end


# =====================================================================
# Test 5: Bremsstrahlung dσ/dk
# =====================================================================
@testset "Bremsstrahlung dsigma/dk" begin
    # Must be positive for valid k
    @test dsigma_dk_brems(0.1, 1.0, CFG) > 0.0
    @test dsigma_dk_brems(0.5, 1.0, CFG) > 0.0
    # Must be zero for invalid k
    @test dsigma_dk_brems(0.0, 1.0, CFG) == 0.0
    @test dsigma_dk_brems(1.0, 1.0, CFG) == 0.0  # k = T not allowed
    @test dsigma_dk_brems(1.5, 1.0, CFG) == 0.0  # k > T
end


# =====================================================================
# Test 5b: Brems table matches direct integration
# =====================================================================
@testset "Brems table vs direct integration" begin
    for T in [0.1, 0.5, 1.0, 2.0, 5.0]
        σ_table  = sigma_brems_table(ND, T)
        σ_direct = sigma_brems_above_kmin(T, CFG.k_min, CFG)
        @test σ_table ≈ σ_direct  rtol=0.02
    end
end


# =====================================================================
# Test 6: Compton sampler — Compton edge
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
# Test 7: Pair energy split in valid range
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
# Test 8: Brems photon energy in valid range
# =====================================================================
@testset "Brems photon energy" begin
    rng = MersenneTwister(2)
    T0 = 2.0
    for _ in 1:100
        k = sample_brems(T0, CFG.k_min, CFG, rng)
        @test k !== nothing
        @test k >= CFG.k_min
        @test k < T0
    end
end


# =====================================================================
# Test 9: Energy conservation in infinite LXe
# =====================================================================
@testset "Energy conservation" begin
    rng = MersenneTwister(3)
    E0 = 2.615
    N = 200
    losses = Float64[]
    for _ in 1:N
        deps = simulate_event(E0, ND, CFG; rng=rng)
        E_dep = sum(d.energy for d in deps)
        push!(losses, E0 - E_dep)
    end
    # Mean loss should be small (a few times the cut energies)
    @test mean(losses) < 5 * CFG.Te_cut
    @test mean(losses) >= 0.0  # should not gain energy
end


# =====================================================================
# Test 10: SS fraction for Tl-208 in plausible range
# =====================================================================
@testset "SS fraction Tl-208" begin
    rng = MersenneTwister(7)
    N = 200
    n_ss = 0
    for _ in 1:N
        deps = simulate_event(2.615, ND, CFG; rng=rng)
        if is_single_site(deps, CFG.dz_resolution)
            n_ss += 1
        end
    end
    frac = n_ss / N
    @test 0.05 < frac < 0.30
end
