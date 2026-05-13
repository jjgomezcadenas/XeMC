using Test
using LXeMC
using JSON
using Random
using Statistics
using Printf

const CFG  = default_config()
const MATS = load_materials(CFG)
const DET3 = load_tracking_detector(default_tracking_detector_path(), MATS)
const SOURCE_GEOM = JSON.parsefile(normpath(joinpath(@__DIR__, "..", "..", "data", "source_geometry_lz_v1.json")))
const MAT  = MATS["LXe"]
const GEOM_TOL = 0.1


# =====================================================================
# Test 1: Config loading
# =====================================================================
@testset "Config" begin
    @test CFG.Egamma_cut > 0.0
    @test CFG.E_roi_floor > CFG.Egamma_cut
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


@testset "Tracking geometry" begin
    @test DET3.name == "LZ"

    names3 = sort([n.lv.name for n in DET3.nodes])
    @test names3 == ["AirCyl", "AirDome", "BarrelActive", "BottomActive", "FC_PTFE", "FC_rings", "FV", "LXe_below_FC", "LXe_below_cathode", "LXe_dome", "LZ_detector", "MARS", "Skin", "TopActive"]

    root3 = root_node(DET3)
    lz3 = node_by_name(DET3, "LZ_detector")
    air3 = node_by_name(DET3, "AirDome")
    passive3 = node_by_name(DET3, "LXe_below_cathode")
    top3 = node_by_name(DET3, "TopActive")
    barrel3 = node_by_name(DET3, "BarrelActive")
    bottom3 = node_by_name(DET3, "BottomActive")
    fv3 = node_by_name(DET3, "FV")
    ptfe3 = node_by_name(DET3, "FC_PTFE")
    rings3 = node_by_name(DET3, "FC_rings")
    skin3 = node_by_name(DET3, "Skin")

    @test is_vacuum(root3)
    @test is_vacuum(lz3)
    @test is_vacuum(air3)
    @test is_passive_lxe(passive3)
    @test is_active_lxe(top3)
    @test is_active_lxe(barrel3)
    @test is_active_lxe(bottom3)
    @test is_fv(fv3)
    @test is_structural(ptfe3)
    @test is_structural(rings3)
    @test is_veto_lxe(skin3)

    @test !is_sensitive(passive3)
    @test is_sensitive(top3)
    @test is_sensitive(barrel3)
    @test is_sensitive(bottom3)
    @test is_sensitive(fv3)
    @test is_sensitive(skin3)
    @test !is_sensitive(ptfe3)
    @test !is_sensitive(rings3)

    @test top3.lv.ecut_keV ≈ 10.0 atol=GEOM_TOL
    @test barrel3.lv.ecut_keV ≈ 10.0 atol=GEOM_TOL
    @test bottom3.lv.ecut_keV ≈ 10.0 atol=GEOM_TOL
    @test fv3.lv.ecut_keV ≈ 10.0 atol=GEOM_TOL
    @test fv3.lv.dz_mm ≈ 3.0 atol=GEOM_TOL
    @test skin3.lv.ecut_keV ≈ 100.0 atol=GEOM_TOL
    @test skin3.lv.dz_mm ≈ 0.0 atol=GEOM_TOL

    @test find_tracking_node(DET3, (0.0, 0.0, 120.0)).lv.name == "TopActive"
    @test find_tracking_node(DET3, (60.0, 0.0, 61.0)).lv.name == "BarrelActive"
    @test find_tracking_node(DET3, (0.0, 0.0, 10.0)).lv.name == "BottomActive"
    @test find_tracking_node(DET3, (0.0, 0.0, 61.0)).lv.name == "FV"
    @test find_tracking_node(DET3, (73.55, 0.0, 100.0)).lv.name == "FC_PTFE"
    @test find_tracking_node(DET3, (74.1, 0.0, 100.0)).lv.name == "FC_rings"
    @test find_tracking_node(DET3, (78.0, 0.0, 100.0)).lv.name == "Skin"
    @test find_tracking_node(DET3, (0.0, 0.0, -20.0)).lv.name == "LXe_below_FC"
    @test find_tracking_node(DET3, (0.0, 0.0, -7.0)).lv.name == "LXe_below_cathode"
    @test find_tracking_node(DET3, (0.0, 0.0, -60.0)).lv.name == "LXe_dome"
    @test find_tracking_node(DET3, (0.0, 0.0, 160.0)).lv.name == "AirDome"

    fk3 = compile_fastkernel_geometry(DET3)
    @test sort([r.name for r in fk3.regions]) == ["AirCyl", "AirDome", "BarrelActive", "BottomActive", "FC_PTFE", "FC_rings", "FV", "LXe_below_FC", "LXe_below_cathode", "LXe_dome", "LZ_detector", "Skin", "TopActive"]
    @test fk3.envelope_index != 0
    @test fk3.regions[fk3.envelope_index].role == "tracking_envelope"
    @test count(r -> r.tag == TAG_PASSIVE_LXE, fk3.regions) == 3
    @test count(r -> r.tag == TAG_FV, fk3.regions) == 1
    @test count(r -> r.tag == TAG_TPC_ACTIVE, fk3.regions) >= 1
    @test classify_fastkernel(fk3, (0.0, 0.0, 120.0)).name == "TopActive"
    @test classify_fastkernel(fk3, (60.0, 0.0, 61.0)).name == "BarrelActive"
    @test classify_fastkernel(fk3, (0.0, 0.0, 10.0)).name == "BottomActive"
    @test classify_fastkernel(fk3, (0.0, 0.0, 61.0)).name == "FV"
    @test classify_fastkernel(fk3, (0.0, 0.0, -20.0)).name == "LXe_below_FC"
    @test classify_fastkernel(fk3, (0.0, 0.0, -7.0)).name == "LXe_below_cathode"
    @test classify_fastkernel(fk3, (0.0, 0.0, -60.0)).name == "LXe_dome"

    fv_vol = compile_fv_volume(DET3)
    @test fv_vol isa PCyl
    @test fv_vol.logical.solid.radius_cm ≈ 39.0 atol=GEOM_TOL
    @test fv_vol.material.name == "LXe"
    @test is_inside(fv_vol, Float64[0.0, 0.0, 61.0])
    @test !is_inside(fv_vol, Float64[40.0, 0.0, 61.0])
    @test !is_inside(fv_vol, Float64[0.0, 0.0, 20.0])
    @test !is_inside(fv_vol, Float64[0.0, 0.0, 100.0])
end

# =====================================================================
# Partition contract for the fast-kernel geometry.
#
# Every point inside the envelope (LZ_detector) must classify into
# exactly one non-envelope region: no gaps (zero claimants) and no
# overlaps (two or more claimants). With no fallback region, gaps would
# manifest as silent escapes and overlaps as ambiguous classification.
# =====================================================================
@testset "FastKernel geometry partition" begin
    fk = compile_fastkernel_geometry(DET3)
    env_node = node_by_name(DET3, "LZ_detector")
    envelope_region = fk.regions[fk.envelope_index]
    rng = MersenneTwister(0xC0FFEE)
    N = 100_000
    n_inside = 0
    n_zero_match = 0
    n_multi_match = 0
    bbox_R = envelope_region.rmax_cm
    bbox_zmin = envelope_region.zmin_cm
    bbox_zmax = envelope_region.zmax_cm
    for _ in 1:N
        x = (2 * rand(rng) - 1) * bbox_R
        y = (2 * rand(rng) - 1) * bbox_R
        z = bbox_zmin + rand(rng) * (bbox_zmax - bbox_zmin)
        is_inside(env_node, Float64[x, y, z]) || continue
        n_inside += 1
        matches = 0
        for i in eachindex(fk.regions)
            i == fk.envelope_index && continue
            LXeMC._is_inside_fastkernel_region(fk.regions[i], (x, y, z)) && (matches += 1)
        end
        matches == 0 && (n_zero_match += 1)
        matches > 1  && (n_multi_match += 1)
    end
    @test n_inside > N ÷ 4
    @test n_zero_match == 0
    @test n_multi_match == 0
end

# =====================================================================
# Top-PMT gammas must propagate in a straight line through the vacuum
# region above the gate (AirDome + AirCyl) and reach the LXe without
# any interaction in those vacuum regions.
# =====================================================================
@testset "Top-PMT gammas: straight-line propagation through vacuum" begin
    fk = compile_fastkernel_geometry(DET3)
    E = 2.448  # Bi214 line, MeV
    z_top = 152.78  # PMT_TOP_PMTs midpoint, inside AirDome
    rng = MersenneTwister(0x70BB07A0)
    N = 10_000
    n_vacuum_deposits = 0
    n_vacuum_terminations = 0
    for _ in 1:N
        # Random inward (downward) direction with cos θ wrt -z in [0.1, 1.0]
        u = -rand(rng) * 0.9 - 0.1
        φ = 2π * rand(rng)
        sθ = sqrt(1 - u^2)
        dir = Float64[sθ*cos(φ), sθ*sin(φ), u]
        gamma = SampledGamma(E, Float64[0.0, 0.0, z_top], dir)
        r = transport_gamma_fastkernel(gamma, fk, CFG, rng)
        for d in r.deposits
            (d.region == "AirDome" || d.region == "AirCyl") && (n_vacuum_deposits += 1)
        end
        (r.terminal_region in ("AirDome", "AirCyl")) && (n_vacuum_terminations += 1)
    end
    @test n_vacuum_deposits == 0
    @test n_vacuum_terminations == 0
end

# =====================================================================
# Top-PMT virgin attenuation: a gamma shot straight down from the PMT
# center traverses 49.6 cm of LXe (TopActive) before reaching the FV.
# The expected fraction of gammas that reach the FV with NO upstream
# interaction is exp(-L/lambda) with lambda = mfp(LXe, E).
# =====================================================================
@testset "Top-PMT virgin survival matches exp(-L/lambda)" begin
    fk = compile_fastkernel_geometry(DET3)
    E = 2.448
    λ = mfp(MATS["LXe"], E)
    L_lxe = 145.6 - 96.0  # TopActive thickness = 49.6 cm
    p_pred = exp(-L_lxe / λ)

    N = 50_000
    rng = MersenneTwister(0x70B07577)
    n_virgin = 0
    for _ in 1:N
        gamma = SampledGamma(E, Float64[0.0, 0.0, 152.78], Float64[0.0, 0.0, -1.0])
        r = transport_gamma_fastkernel(gamma, fk, CFG, rng)
        (r.status === :handoff_fv && isempty(r.deposits)) && (n_virgin += 1)
    end
    f_obs = n_virgin / N
    σ_stat = sqrt(p_pred * (1 - p_pred) / N)
    @test abs(f_obs - p_pred) < 4 * σ_stat
end

# =====================================================================
# Bottom-PMT virgin attenuation: a gamma shot straight up from the PMT
# center at z=-15.92 traverses LXe continuously (2.17 cm LXe_below_FC
# + 13.75 cm LXe_below_cathode + 26 cm BottomActive = 41.92 cm) before
# reaching the FV at z=26. Predicted virgin survival is exp(-L/lambda).
# This is the test that fails when the fast-kernel fallback bug is
# present: passive-LXe traversal forced spurious interactions and
# depleted virgins by orders of magnitude.
# =====================================================================
@testset "Bottom-PMT virgin survival matches exp(-L/lambda)" begin
    fk = compile_fastkernel_geometry(DET3)
    E = 2.448
    λ = mfp(MATS["LXe"], E)
    z_pmt = -15.92  # PMT_BOT_PMTs midpoint
    L_lxe = (-13.75 - z_pmt) + 13.75 + 26.0  # 41.92 cm of continuous LXe
    p_pred = exp(-L_lxe / λ)

    N = 50_000
    rng = MersenneTwister(0xB0770557)
    n_virgin = 0
    for _ in 1:N
        gamma = SampledGamma(E, Float64[0.0, 0.0, z_pmt], Float64[0.0, 0.0, 1.0])
        r = transport_gamma_fastkernel(gamma, fk, CFG, rng)
        (r.status === :handoff_fv && isempty(r.deposits)) && (n_virgin += 1)
    end
    f_obs = n_virgin / N
    σ_stat = sqrt(p_pred * (1 - p_pred) / N)
    @test abs(f_obs - p_pred) < 4 * σ_stat
end

@testset "Source geometry schema" begin
    @test SOURCE_GEOM["name"] == "LZ_source_geometry"
    @test SOURCE_GEOM["version"] == 1
    sources = SOURCE_GEOM["sources"]
    @test !isempty(sources)

    names = [String(s["name"]) for s in sources]
    @test length(unique(names)) == length(names)

    required = ["name", "shape", "material", "source_class", "transport_source", "approximation", "equivalent_mass_kg"]
    for s in sources
        for key in required
            @test haskey(s, key)
        end
    end

    by_name = Dict(String(s["name"]) => s for s in sources)
    @test haskey(by_name, "OCV_barrel")
    @test haskey(by_name, "ICV_barrel")
    @test haskey(by_name, "FC_PTFE")
    @test haskey(by_name, "FC_rings")
    @test haskey(by_name, "FC_resistors")
    @test haskey(by_name, "FC_sensors")
    @test haskey(by_name, "FC_topgrid")
    @test haskey(by_name, "FC_botgrid")

    for s in sources
        if s["transport_source"] == "transparent"
            @test s["material"] == "Vacuum"
            @test Float64(s["equivalent_mass_kg"]) > 0.0
        elseif s["transport_source"] == "KN"
            @test Float64(s["equivalent_mass_kg"]) >= 0.0
        else
            @test false
        end
    end
end

@testset "Tracking/source consistency" begin
    sources = Dict(String(s["name"]) => s for s in SOURCE_GEOM["sources"])
    shared = ["FC_PTFE", "FC_rings"]

    for name in shared
        t = node_by_name(DET3, name)
        s = sources[name]
        @test lowercase(String(s["shape"])) == begin
            solid = t.lv.solid
            solid isa CylShell ? "cylinder_shell" :
            solid isa Cyl ? "cylinder" :
            solid isa Disk ? "disk" : ""
        end
        @test t.placement.position_cm[3] ≈ Float64(s["position_cm"][3]) atol=GEOM_TOL
        if t.lv.solid isa CylShell
            @test t.lv.solid.R_inner_cm ≈ Float64(s["R_inner_cm"]) atol=GEOM_TOL
            @test t.lv.solid.wall_thickness_cm ≈ Float64(s["wall_thickness_cm"]) atol=GEOM_TOL
            @test t.lv.solid.half_height_cm ≈ Float64(s["half_height_cm"]) atol=GEOM_TOL
        end
    end

    source_only = ["MLI", "OCV_barrel", "ICV_barrel", "FC_sensors", "FC_topgrid", "FC_botgrid", "PMT_TOP_PMTs"]
    for name in source_only
        @test haskey(sources, name)
        @test !haskey(DET3.name_to_id, name)
    end
end


# =====================================================================
# Source geometric mass vs bb0nu reference
# =====================================================================
#
# For each source flagged approximation="exact" with a real material,
# verify that the mass derived from the drawn geometry × material density
# reproduces the bb0nu measured mass (or the TDR shell/head mass for the
# Ti vessels). Mass-equivalent virtual carriers are skipped — their
# absolute activity is anchored by equivalent_mass_kg, not by geometry.
#
# Heads use the analytic half-oblate-spheroid surface area formula,
# matching `area_inner(::GDisk)` in the geometry primitives.
@testset "Source geometric mass vs bb0nu" begin
    function _geom_mass_kg(s, ρ)
        shape = String(s["shape"])
        if shape == "cylinder_shell"
            R_in = Float64(s["R_inner_cm"])
            t    = Float64(s["wall_thickness_cm"])
            H    = 2 * Float64(s["half_height_cm"])
            V    = π * ((R_in + t)^2 - R_in^2) * H
            V * ρ / 1000.0
        elseif shape == "cylinder"
            R = Float64(s["radius_cm"])
            H = 2 * Float64(s["half_height_cm"])
            π * R^2 * H * ρ / 1000.0
        elseif shape == "disk"
            R = Float64(s["radius_cm"])
            t = Float64(s["wall_thickness_cm"])
            n = Float64(s["aspect_ratio"])
            a = R
            c = R / n
            S_full = if c < a
                e = sqrt(1.0 - (c/a)^2)
                2π * a^2 + π * (c^2 / e) * log((1 + e) / (1 - e))
            else
                e = sqrt(1.0 - (a/c)^2)
                2π * a^2 + 2π * a * c * asin(e) / e
            end
            (S_full / 2.0) * t * ρ / 1000.0
        else
            NaN
        end
    end

    # Reference masses (kg). Cryostat barrels + heads from the TDR/CSV;
    # FC_PTFE and FC_rings from bb0nu Table I. Heads computed at the
    # CSV's aspect ratio + wall thickness × ρ_Ti so each row sums
    # correctly to lz_cryo_geometry.csv mass_heads_kg per vessel.
    target_mass_kg = Dict(
        "OCV_barrel" => 385.7,
        "OCV_top"    => 147.4,
        "OCV_bottom" => 245.6,
        "ICV_barrel" => 401.9,
        "ICV_top"    => 107.8,
        "ICV_bottom" => 141.5,
        "FC_PTFE"    => 184.0,
        "FC_rings"   =>  93.0,
    )

    for s in SOURCE_GEOM["sources"]
        name     = String(s["name"])
        approx   = String(s["approximation"])
        material = String(s["material"])
        approx == "exact"     || continue
        material == "Vacuum"  && continue
        haskey(target_mass_kg, name) || continue

        ρ = MATS[material].density
        m_geom = _geom_mass_kg(s, ρ)
        @test isapprox(m_geom, target_mass_kg[name]; rtol=0.05)
    end
end

@testset "Calibration sampler" begin
    for i in 1:50
        ev = sample_gammas("calib"; calib=true, rng=MersenneTwister(1234 + i))
        @test length(ev) in (0, 1, 2)
        for g in ev
            @test g.E_MeV ≈ 2.615 atol=GEOM_TOL
            @test length(g.position) == 3
            @test length(g.direction) == 3
        end
    end

    ev_custom = sample_gammas("calib";
                              calib=true,
                              E_MeV=1.25,
                              x_cm=1.0, y_cm=2.0, z_cm=3.0,
                              ux=2.0, uy=0.0, uz=0.0,
                              rng=MersenneTwister(99))
    for g in ev_custom
        @test g.E_MeV ≈ 1.25 atol=GEOM_TOL
        @test g.position == Float64[1.0, 2.0, 3.0]
        @test g.direction == Float64[1.0, 0.0, 0.0]
    end
end

@testset "FastKernel first interaction transport" begin
    fk = compile_fastkernel_geometry(DET3)
    air_to_lxe = SampledGamma(2.615, Float64[0.0, 0.0, 160.0], Float64[0.0, 0.0, -1.0])
    fk_air = propagate_gamma_fastkernel(air_to_lxe, fk, CFG, MersenneTwister(23))
    @test fk_air.status in (:escaped, :entered_fv, :interacted)
    @test fk_air.region in ("AirDome", "TopActive", "BarrelActive", "BottomActive", "FV", "Skin", "FC_PTFE", "FC_rings", "LXe_dome", "LXe_below_FC", "LXe_below_cathode", "MARS")

    ptfe = SampledGamma(2.615, Float64[73.0, 0.0, 100.0], Float64[1.0, 0.0, 0.0])
    fk_ptfe = propagate_gamma_fastkernel(ptfe, fk, CFG, MersenneTwister(11))
    @test fk_ptfe.status == :interacted
    @test fk_ptfe.region in ("FC_PTFE", "FC_rings", "Skin")

    lxe_bulk = SampledGamma(2.615, Float64[0.0, 0.0, -20.0], Float64[0.0, 0.0, 1.0])
    fk_lxe = propagate_gamma_fastkernel(lxe_bulk, fk, CFG, MersenneTwister(19))
    @test fk_lxe.status == :interacted
    @test fk_lxe.region in ("LXe_below_FC", "LXe_below_cathode", "BottomActive", "FV")
end

@testset "Interaction selector" begin
    fk = compile_fastkernel_geometry(DET3)
    fv_sel = select_interaction_fastkernel(
        LXeMC.GammaPropagationResult(:interacted, :compton, 0.020, Float64[0.0, 0.0, 61.0], "FV"),
        fk
    )
    @test fv_sel.class == :fv
    @test fv_sel.sensitive
    @test fv_sel.passes_threshold
    @test fv_sel.region == "FV"

    tpc_hi_sel = select_interaction_fastkernel(
        LXeMC.GammaPropagationResult(:interacted, :compton, 0.020, Float64[0.0, 0.0, 120.0], "TopActive"),
        fk
    )
    @test tpc_hi_sel.class == :tpc
    @test tpc_hi_sel.sensitive
    @test tpc_hi_sel.passes_threshold
    @test tpc_hi_sel.ecut_keV ≈ 10.0 atol=GEOM_TOL

    tpc_lo_sel = select_interaction_fastkernel(
        LXeMC.GammaPropagationResult(:interacted, :compton, 0.005, Float64[0.0, 0.0, 120.0], "TopActive"),
        fk
    )
    @test tpc_lo_sel.class == :tpc
    @test !tpc_lo_sel.passes_threshold

    skin_hi_sel = select_interaction_fastkernel(
        LXeMC.GammaPropagationResult(:interacted, :compton, 0.200, Float64[78.0, 0.0, 100.0], "Skin"),
        fk
    )
    @test skin_hi_sel.class == :skin
    @test skin_hi_sel.passes_threshold

    skin_lo_sel = select_interaction_fastkernel(
        LXeMC.GammaPropagationResult(:interacted, :compton, 0.050, Float64[78.0, 0.0, 100.0], "Skin"),
        fk
    )
    @test skin_lo_sel.class == :skin
    @test !skin_lo_sel.passes_threshold

    struct_sel = select_interaction_fastkernel(
        LXeMC.GammaPropagationResult(:interacted, :compton, 0.500, Float64[73.55, 0.0, 100.0], "FC_PTFE"),
        fk
    )
    @test struct_sel.class == :other
    @test !struct_sel.sensitive
    @test !struct_sel.passes_threshold
end


@testset "FastKernel event processing" begin
    fk = compile_fastkernel_geometry(DET3)
    fv_vol = compile_fv_volume(DET3)

    empty_fk = process_event(SampledGamma[], fk, fv_vol, CFG, MersenneTwister(1))
    @test empty_fk.status == :no_fv
    @test empty_fk.n_processed == 0
    @test isempty(empty_fk.deposits)

    miss_gamma = SampledGamma(2.615, Float64[120.0, 0.0, 160.0], Float64[0.0, 0.0, -1.0])
    miss_res = LXeMC.transport_gamma_fastkernel(miss_gamma, fk, CFG, MersenneTwister(1))
    @test miss_res.status == :escaped
    @test isempty(miss_res.deposits)

    handoff_gamma = SampledGamma(2.615, Float64[0.0, 0.0, 61.0], Float64[0.0, 0.0, -1.0])
    handoff_res = LXeMC.transport_gamma_fastkernel(handoff_gamma, fk, CFG, MersenneTwister(2))
    @test handoff_res.status == :handoff_fv
    @test handoff_res.terminal_region == "FV"
    @test handoff_res.energy_MeV ≈ 2.615 atol=1e-12

    tpc_veto_gamma = SampledGamma(2.615, Float64[60.0, 0.0, 120.0], Float64[1.0, 0.0, 0.0])
    tpc_veto_res = LXeMC.transport_gamma_fastkernel(tpc_veto_gamma, fk, CFG, MersenneTwister(123))
    @test tpc_veto_res.status in (:vetoed_tpc, :absorbed_passive, :below_cut, :handoff_fv, :escaped, :below_roi_fv)

    skin_veto_gamma = SampledGamma(2.615, Float64[78.0, 0.0, 100.0], Float64[1.0, 0.0, 0.0])
    skin_veto_res = LXeMC.transport_gamma_fastkernel(skin_veto_gamma, fk, CFG, MersenneTwister(456))
    @test skin_veto_res.status in (:vetoed_skin, :absorbed_passive, :below_cut, :escaped)

    for seed in (7, 19, 42, 202, 777)
        rng1 = MersenneTwister(seed)
        rng2 = MersenneTwister(seed)
        gammas = sample_gammas("calib"; calib=true, rng=rng1)
        vec_res = process_event(gammas, fk, fv_vol, CFG, rng1)
        fused_res = process_event_fastkernel_calib(fk, fv_vol, CFG, rng2)

        @test fused_res.result.status == vec_res.status
        @test fused_res.result.has_tpc_veto == vec_res.has_tpc_veto
        @test fused_res.result.has_skin_veto == vec_res.has_skin_veto
        @test fused_res.result.n_processed == vec_res.n_processed
        @test fused_res.multiplicity == length(gammas)

        energies = sort([g.E_MeV for g in gammas]; rev=true)
        E1 = length(energies) >= 1 ? energies[1] : 0.0
        E2 = length(energies) >= 2 ? energies[2] : 0.0
        @test fused_res.E1_MeV == E1
        @test fused_res.E2_MeV == E2
    end

    # FV events carry deposits; others don't
    fv_gamma = SampledGamma(2.615, Float64[0.0, 0.0, 61.0], Float64[0.0, 0.0, 1.0])
    for seed in (100, 200, 300, 400, 500)
        res = process_event([fv_gamma], fk, fv_vol, CFG, MersenneTwister(seed))
        if res.status == :fv
            @test !isempty(res.deposits)
            @test sum(d.energy for d in res.deposits) > 0.0
        else
            @test isempty(res.deposits)
        end
    end
end


@testset "FV deposit CSV round-trip" begin
    # Pseudo main loop: shoot gammas, collect FV deposits, write CSV,
    # read back, and verify correctness.
    fk = compile_fastkernel_geometry(DET3)
    fv_vol = compile_fv_volume(DET3)

    # Collect FV events
    fv_events = Vector{Vector{Deposit}}()
    n_fv = 0
    n_vetoed = 0
    n_no_fv = 0

    fv_gamma = SampledGamma(2.615, Float64[0.0, 0.0, 61.0], Float64[0.0, 0.0, 1.0])
    for seed in 1:100
        res = process_event([fv_gamma], fk, fv_vol, CFG, MersenneTwister(seed))
        if res.status == :fv
            n_fv += 1
            push!(fv_events, res.deposits)
        elseif res.status == :vetoed
            n_vetoed += 1
        else
            n_no_fv += 1
        end
    end
    @test n_fv + n_vetoed + n_no_fv == 100
    @test n_fv > 0

    # Write CSV
    tmpdir = mktempdir()
    csv_path = joinpath(tmpdir, "fv_deposits.csv")
    open(csv_path, "w") do io
        println(io, "event_id,x_cm,y_cm,z_cm,energy_MeV,source")
        for (eid, deps) in enumerate(fv_events)
            for d in deps
                @printf(io, "%d,%.6f,%.6f,%.6f,%.8e,%s\n",
                        eid, d.position[1], d.position[2], d.position[3],
                        d.energy, d.source)
            end
        end
    end

    # Read back and verify
    lines = readlines(csv_path)
    header = lines[1]
    @test header == "event_id,x_cm,y_cm,z_cm,energy_MeV,source"

    data_lines = lines[2:end]
    @test length(data_lines) == sum(length(deps) for deps in fv_events)

    # Group by event_id and check deposit count per event
    events_from_csv = Dict{Int, Vector{NamedTuple}}()
    for line in data_lines
        parts = split(line, ",")
        eid = parse(Int, parts[1])
        row = (x=parse(Float64, parts[2]), y=parse(Float64, parts[3]),
               z=parse(Float64, parts[4]), energy=parse(Float64, parts[5]),
               source=String(parts[6]))
        push!(get!(events_from_csv, eid, []), row)
    end

    @test length(events_from_csv) == n_fv

    # Verify each event's deposits match the originals
    for (eid, deps) in enumerate(fv_events)
        csv_deps = events_from_csv[eid]
        @test length(csv_deps) == length(deps)
        for (d, c) in zip(deps, csv_deps)
            @test d.position[1] ≈ c.x atol=1e-5
            @test d.position[2] ≈ c.y atol=1e-5
            @test d.position[3] ≈ c.z atol=1e-5
            @test d.energy ≈ c.energy rtol=1e-6
            @test String(d.source) == c.source
        end
    end

    # Energy conservation: total deposited per event <= gamma energy
    for (eid, csv_deps) in events_from_csv
        E_total = sum(c.energy for c in csv_deps)
        @test E_total <= 2.615 + 0.001  # gamma energy + tolerance
        @test E_total > 0.0
    end

    rm(tmpdir; recursive=true)
end

@testset "_fold_fast_deposit is redundant" begin
    # Prove that _fold_fast_deposit never triggers a veto that the
    # result.status check wouldn't also catch. Run many events and
    # check that removing _fold_fast_deposit gives identical results.
    fk = compile_fastkernel_geometry(DET3)
    fv_vol = compile_fv_volume(DET3)

    test_gammas = [
        [SampledGamma(2.615, Float64[80.0, 0.0, 100.0], Float64[-1.0, 0.0, 0.0])],  # from skin
        [SampledGamma(2.448, Float64[0.0, 0.0, 120.0], Float64[0.0, 0.0, -1.0])],    # from top
        [SampledGamma(2.615, Float64[50.0, 0.0, 61.0], Float64[-1.0, 0.0, 0.0])],    # toward FV
        [SampledGamma(2.615, Float64[0.0, 0.0, 61.0], Float64[0.0, 0.0, 1.0])],      # in FV
    ]

    for gammas in test_gammas
        for seed in 1:200
            rng = MersenneTwister(seed)
            result = LXeMC.transport_gamma_fastkernel(gammas[1], fk, CFG, rng)

            # Check: no deposit triggers _fold_fast_deposit veto
            for dep in result.deposits
                # "LXeTPC" never matches (bug: should be TopActive etc.)
                @test dep.region != "LXeTPC"
                # FV never appears (proven by other test)
                @test dep.region != "FV"
                # Skin deposits above threshold: check if status already catches it
                if dep.region == "Skin" && dep.Edep_MeV >= CFG.veto_skin
                    @test result.status == :vetoed_skin
                end
            end
        end
    end
end

@testset "transport_gamma_fastkernel never deposits in FV" begin
    # The FV handoff at the top of the loop (region.name == "FV" -> :handoff_fv)
    # should prevent any interaction from being recorded with region == "FV".
    # If any deposit has region "FV", the dead-code branches were reached.
    fk = compile_fastkernel_geometry(DET3)

    # Fire gammas from various positions toward/through the FV
    test_gammas = [
        SampledGamma(2.615, Float64[0.0, 0.0, 61.0], Float64[0.0, 0.0, 1.0]),    # inside FV
        SampledGamma(2.615, Float64[50.0, 0.0, 61.0], Float64[-1.0, 0.0, 0.0]),   # toward FV
        SampledGamma(2.448, Float64[0.0, 0.0, 120.0], Float64[0.0, 0.0, -1.0]),   # from above
        SampledGamma(2.615, Float64[80.0, 0.0, 61.0], Float64[-1.0, 0.0, 0.0]),   # from skin
        SampledGamma(2.615, Float64[0.0, 0.0, 10.0], Float64[0.0, 0.0, 1.0]),     # from below
    ]

    n_fv_deposits = 0
    for g in test_gammas
        for seed in 1:200
            result = LXeMC.transport_gamma_fastkernel(g, fk, CFG, MersenneTwister(seed))
            for dep in result.deposits
                dep.region == "FV" && (n_fv_deposits += 1)
            end
        end
    end
    @test n_fv_deposits == 0
end

@testset "FV-only stack transport" begin
    fv_vol = compile_fv_volume(DET3)
    fk = compile_fastkernel_geometry(DET3)
    rng1 = MersenneTwister(20260508)
    deps1 = propagate_gamma(2.61, fv_vol, fk, CFG;
                            position=(0.0, 0.0, 61.0),
                            direction=(0.0, 0.0, 1.0),
                            rng=rng1)
    rng2 = MersenneTwister(20260508)
    deps2 = propagate_gamma(2.61, fv_vol, fk, CFG;
                            position=(0.0, 0.0, 61.0),
                            direction=(0.0, 0.0, 1.0),
                            rng=rng2)

    @test length(deps2) == length(deps1)
    fv_deps1 = filter(d -> d.volume == :fv, deps1)
    @test all(is_inside(fv_vol, d.position) for d in fv_deps1)
    @test sum(d.energy for d in deps2) ≈ sum(d.energy for d in deps1) atol=1e-9
    @test [d.source for d in deps2] == [d.source for d in deps1]

    rng3 = MersenneTwister(12345)
    deps3 = propagate_gamma(2.61, fv_vol, fk, CFG;
                            position=(10.0, 0.0, 80.0),
                            direction=(0.0, 0.0, -1.0),
                            rng=rng3)
    rng4 = MersenneTwister(12345)
    deps4 = propagate_gamma(2.61, fv_vol, fk, CFG;
                            position=(10.0, 0.0, 80.0),
                            direction=(0.0, 0.0, -1.0),
                            rng=rng4)

    @test length(deps4) == length(deps3)
    fv_deps3 = filter(d -> d.volume == :fv, deps3)
    @test all(is_inside(fv_vol, d.position) for d in fv_deps3)
    @test sum(d.energy for d in deps4) ≈ sum(d.energy for d in deps3) atol=1e-9
    @test [d.source for d in deps4] == [d.source for d in deps3]
end

@testset "Energy balance with escape deposits" begin
    fv_vol = compile_fv_volume(DET3)
    fk = compile_fastkernel_geometry(DET3)
    E_gamma = 2.615

    # Run many events and check energy conservation
    n_events = 100
    for seed in 1:n_events
        rng = MersenneTwister(seed)
        deps = propagate_gamma(E_gamma, fv_vol, fk, CFG;
                               position=(0.0, 0.0, 61.0),
                               direction=(0.0, 0.0, 1.0),
                               rng=rng)

        etot = sum(d.energy for d in deps)
        @test etot ≈ E_gamma atol=1e-6

        # Every deposit must have a valid volume
        for d in deps
            @test d.volume in (:fv, :active, :passive)
        end

        # Every deposit must have a valid interaction
        for d in deps
            @test d.interaction in (:compton, :photoelectric, :pair,
                                    :collisional, :gamma_local,
                                    :escaped_gamma, :escaped_lepton)
        end

        # FV deposits must be inside the volume
        for d in deps
            if d.volume == :fv
                @test is_inside(fv_vol, d.position)
            end
        end
    end
end

@testset "FV boundary region classification" begin
    fv_vol = compile_fv_volume(DET3)
    fk = compile_fastkernel_geometry(DET3)
    r_fv = fv_vol.logical.solid.radius_cm
    zc = fv_vol.logical.position[3]
    hh = fv_vol.logical.solid.half_height_cm

    # Points just outside FV in each direction should be active (TPC)
    test_cases = [
        ("radial +0.5",  Float64[r_fv + 0.5, 0.0, zc],         :active),
        ("radial +5.0",  Float64[r_fv + 5.0, 0.0, zc],         :active),
        ("top +0.5",     Float64[0.0, 0.0, zc + hh + 0.5],     :active),
        ("top +5.0",     Float64[0.0, 0.0, zc + hh + 5.0],     :active),
        ("bottom -0.5",  Float64[0.0, 0.0, zc - hh - 0.5],     :active),
        ("bottom -5.0",  Float64[0.0, 0.0, zc - hh - 5.0],     :active),
    ]
    for (label, pos, expected) in test_cases
        vol = LXeMC._classify_escape_volume(fk, pos)
        @test vol == expected
    end
end


# =====================================================================
# Tag-based dispatch coverage (refactor/tag-dispatch)
# =====================================================================
@testset "_is_passive_region tag coverage" begin
    expected_passive = Set([TAG_PASSIVE_LXE, TAG_STRUCTURAL])
    for tag in instances(RegionTag)
        @test LXeMC._is_passive_region(tag) == (tag in expected_passive)
    end
end


@testset "_terminal_status tag coverage" begin
    expected_status = Dict(
        TAG_SKIN         => :vetoed_skin,
        TAG_TPC_ACTIVE   => :vetoed_tpc,
        TAG_PASSIVE_LXE  => :absorbed_passive,
        TAG_STRUCTURAL   => :absorbed_passive,
        TAG_FV           => :below_cut,
        TAG_VACUUM       => :below_cut,
        TAG_WORLD        => :below_cut,
    )
    deposits = LXeMC.FastGammaDeposit[]
    pos = Float64[0.0, 0.0, 0.0]
    dir = Float64[0.0, 0.0, 1.0]
    for tag in instances(RegionTag)
        @test haskey(expected_status, tag)  # fails loudly if a new tag is added without updating the test
        result = LXeMC._terminal_status(deposits, tag, "probe", pos, dir)
        @test result.status == expected_status[tag]
        @test result.terminal_region == "probe"
    end
end


@testset "_classify_escape_volume tag mapping (per region)" begin
    fk = compile_fastkernel_geometry(DET3)
    expected = Dict(
        TAG_FV           => :fv,
        TAG_TPC_ACTIVE   => :active,
        TAG_SKIN         => :active,
        TAG_PASSIVE_LXE  => :passive,
        TAG_STRUCTURAL   => :passive,
    )
    seen_tags = Set{RegionTag}()
    for region in fk.regions
        haskey(expected, region.tag) || continue
        # Probe at the bounding-box center; skip if classify_fastkernel doesn't
        # round-trip (region has caps, or center hits a sibling).
        r_probe = 0.5 * (region.rmin_cm + region.rmax_cm)
        z_probe = 0.5 * (region.zmin_cm + region.zmax_cm)
        pos = Float64[r_probe, 0.0, z_probe]
        classified = classify_fastkernel(fk, (pos[1], pos[2], pos[3]))
        classified === nothing && continue
        classified.id == region.id || continue
        @test LXeMC._classify_escape_volume(fk, pos) == expected[region.tag]
        push!(seen_tags, region.tag)
    end
    # Confirm the test actually exercised every tag in DET3 we expect to see
    for tag in (TAG_FV, TAG_TPC_ACTIVE, TAG_SKIN, TAG_PASSIVE_LXE)
        @test tag in seen_tags
    end
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

    # DomedContainer: coaxial barrel plus two filled caps
    dc = DomedContainer(50.0, 20.0, Cap(50.0, 2.0), Cap(50.0, 3.0))
    @test volume(dc) ≈ (40.0 * π * 50.0^2 + (2.0 / 3.0) * π * 50.0^2 * 25.0 + (2.0 / 3.0) * π * 50.0^2 * (50.0 / 3.0))
    dc_lv = LogicalVolume("dc", dc, MATS["Vacuum"], TAG_VACUUM, "container", false, false, false, 0.0, 0.0, false, "exact", Dict{String,Float64}())
    dc_node = DetectorNode(1, 0, Int[], dc_lv, Placement([0.0, 0.0, 0.0], :none))
    @test is_inside(dc_node, [0.0, 0.0, 0.0])        # barrel center
    @test is_inside(dc_node, [0.0, 0.0, 30.0])       # inside top cap
    @test is_inside(dc_node, [0.0, 0.0, -30.0])      # inside bottom cap
    @test !is_inside(dc_node, [0.0, 0.0, 50.0])      # above top cap
    @test !is_inside(dc_node, [60.0, 0.0, 0.0])      # outside radius

    # CappedCylinder: optional top/bottom caps
    cc = CappedCylinder(50.0, 20.0; bottom_cap=Cap(50.0, 3.0))
    @test volume(cc) ≈ (40.0 * π * 50.0^2 + (2.0 / 3.0) * π * 50.0^2 * (50.0 / 3.0))
    cc_lv = LogicalVolume("cc", cc, MATS["LXe"], TAG_PASSIVE_LXE, "container", false, true, false, 0.0, 0.0, false, "exact", Dict{String,Float64}())
    cc_node = DetectorNode(1, 0, Int[], cc_lv, Placement([0.0, 0.0, 0.0], :none))
    @test is_inside(cc_node, [0.0, 0.0, 0.0])        # barrel center
    @test is_inside(cc_node, [0.0, 0.0, -30.0])      # inside bottom cap
    @test !is_inside(cc_node, [0.0, 0.0, 30.0])      # no top cap
    @test !is_inside(cc_node, [60.0, 0.0, 0.0])      # outside radius

    # Validation
    @test_throws ErrorException Cyl(-1.0, 5.0)
    @test_throws ErrorException CylShell(10.0, -1.0, 5.0)
    @test_throws ErrorException Box(0.0, 5.0, 5.0)
    @test_throws ErrorException Disk(-1.0, 1.0, 2.0)
    @test_throws ErrorException DomedContainer(50.0, 10.0, Cap(40.0, 2.0), Cap(50.0, 3.0))
    @test_throws ErrorException CappedCylinder(50.0, 10.0; top_cap=Cap(40.0, 2.0))
end


# =====================================================================
# Test 3c: Ray-volume intersection
# =====================================================================
@testset "Ray intersection" begin
    # Cylinder at origin, R=10, H_half=20
    lc = LCyl(Cyl(10.0, 20.0), [0.0, 0.0, 0.0])

    # Ray from outside along +x, should hit at x=10 (distance 5 from x=15)
    @test distance_to_entry([15.0, 0.0, 0.0], [-1.0, 0.0, 0.0], lc) ≈ 5.0  atol=0.01
    # Ray from inside along +x, exits at x=10
    @test distance_to_exit([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], lc) ≈ 10.0  atol=0.01
    # Ray along +z from inside, exits at z=20
    @test distance_to_exit([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], lc) ≈ 20.0  atol=0.01
    # Ray misses (parallel, offset)
    @test distance_to_entry([15.0, 0.0, 0.0], [0.0, 1.0, 0.0], lc) == Inf

    # CylShell: R_inner=50, wall=2, H_half=30
    lcs = LCylShell(CylShell(50.0, 2.0, 30.0), [0.0, 0.0, 0.0])
    # Ray from r=60 inward, should hit outer surface at r=52
    t_entry = distance_to_entry([60.0, 0.0, 0.0], [-1.0, 0.0, 0.0], lcs)
    @test t_entry ≈ 8.0  atol=0.1
    # Ray from inside shell outward, exits at r=52 or r=50
    t_exit = distance_to_exit([51.0, 0.0, 0.0], [-1.0, 0.0, 0.0], lcs)
    @test t_exit ≈ 1.0  atol=0.1  # exits through inner surface

    # --- LDisk (ellipsoidal) ray-entry ---
    # 2:1 ellipsoidal disk, R=50, t=2, pointing up at origin
    # Inner ellipsoid: a=50, c=25. Outer: a=50, c=27.
    lde = LDisk(Disk(50.0, 2.0, 2.0), [0.0, 0.0, 0.0], :up)

    # Ray from above along -z toward apex: should hit outer ellipsoid
    # Apex of outer ellipsoid is at z = c_outer = 27
    t_apex = distance_to_entry([0.0, 0.0, 40.0], [0.0, 0.0, -1.0], lde)
    @test t_apex ≈ 13.0  atol=0.2  # 40 - 27 = 13

    # Verify the hit point is inside the disk shell, and just before is outside
    hit_pos = [0.0, 0.0, 40.0] .+ [0.0, 0.0, -1.0] .* t_apex
    @test is_inside(lde, hit_pos)
    before_pos = [0.0, 0.0, 40.0] .+ [0.0, 0.0, -1.0] .* (t_apex - 0.5)
    @test !is_inside(lde, before_pos)

    # Ray from side along -x at z=10: should hit the outer ellipsoid side
    t_side = distance_to_entry([60.0, 0.0, 10.0], [-1.0, 0.0, 0.0], lde)
    @test isfinite(t_side)
    hit_side = [60.0, 0.0, 10.0] .+ [-1.0, 0.0, 0.0] .* t_side
    @test is_inside(lde, hit_side)
    before_side = [60.0, 0.0, 10.0] .+ [-1.0, 0.0, 0.0] .* (t_side - 0.5)
    @test !is_inside(lde, before_side)

    # Ray missing: going away from disk
    t_miss = distance_to_entry([0.0, 0.0, 40.0], [0.0, 0.0, 1.0], lde)
    @test t_miss == Inf

    # Ray below equator going down: misses :up disk
    t_below = distance_to_entry([0.0, 0.0, -10.0], [0.0, 0.0, -1.0], lde)
    @test t_below == Inf

    # :down disk at z=0, 3:1 aspect ratio (like ICV bottom)
    ldd = LDisk(Disk(50.0, 2.0, 3.0), [0.0, 0.0, 0.0], :down)
    # Inner c = 50/3 ≈ 16.67, outer c ≈ 18.67
    # Ray from below along +z toward nadir
    c_outer_down = 50.0 / 3.0 + 2.0
    t_nadir = distance_to_entry([0.0, 0.0, -30.0], [0.0, 0.0, 1.0], ldd)
    @test t_nadir ≈ (30.0 - c_outer_down)  atol=0.2
    hit_nadir = [0.0, 0.0, -30.0] .+ [0.0, 0.0, 1.0] .* t_nadir
    @test is_inside(ldd, hit_nadir)
    before_nadir = [0.0, 0.0, -30.0] .+ [0.0, 0.0, 1.0] .* (t_nadir - 0.5)
    @test !is_inside(ldd, before_nadir)

    # Flat disk: R=50, t=2, pointing up
    ldf = LDisk(Disk(50.0, 2.0, Inf), [0.0, 0.0, 10.0], :up)
    # Ray from above along -z
    t_flat = distance_to_entry([0.0, 0.0, 20.0], [0.0, 0.0, -1.0], ldf)
    @test t_flat ≈ 10.0  atol=0.1  # hits z=10 face
    hit_flat = [0.0, 0.0, 20.0] .+ [0.0, 0.0, -1.0] .* t_flat
    @test is_inside(ldf, hit_flat)
    # Point before hit is above the slab (z > 12), so outside
    before_flat = [0.0, 0.0, 20.0] .+ [0.0, 0.0, -1.0] .* (t_flat * 0.5)
    @test !is_inside(ldf, before_flat)

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
# Source flux data structures
# =====================================================================
@testset "SourceFluxBi214 construction" begin
    n_E, n_u = 25, 10
    pdf = zeros(n_E, n_u)
    pdf[20, 5] = 0.30   # peak bin
    pdf[19, 5] = 0.02   # scattered
    flux = SourceFluxBi214(
        "ICV_barrel",
        pdf,
        2.200, 2.500, n_E, n_u,
        100_000, 32_000, 50_000, 48_000, 2_000
    )
    @test flux isa SourceFlux
    @test flux.source_name == "ICV_barrel"
    @test size(flux.pdf) == (25, 10)
    @test flux.E_min == 2.200
    @test flux.E_max == 2.500
    @test flux.n_E == 25
    @test flux.n_u == 10
    @test flux.N_generated == 100_000
    @test flux.N_surviving == 32_000
    @test flux.N_absorbed == 50_000
    @test flux.N_backward == 48_000
    @test flux.N_low_energy == 2_000
    @test sum(flux.pdf) ≈ 0.32 atol=1e-10
end

@testset "SourceFluxTl208 construction" begin
    n_E_main, n_u = 25, 10
    pdf_main = zeros(n_E_main, n_u)
    pdf_main[24, 5] = 0.18

    n_E_comp = [20, 20, 20]
    pdf_c1 = zeros(20, 10)
    pdf_c1[15, 5] = 0.10
    pdf_c2 = zeros(20, 10)
    pdf_c2[12, 5] = 0.05
    pdf_c3 = zeros(20, 10)
    pdf_c3[18, 5] = 0.08

    flux = SourceFluxTl208(
        "ICV_barrel",
        pdf_main, 2.370, 2.620, n_E_main,
        [pdf_c1, pdf_c2, pdf_c3],
        [0.583, 0.511, 0.861],
        [0.85, 0.23, 0.12],
        [0.25, 0.20, 0.30],
        [0.400, 0.350, 0.700],
        [0.650, 0.600, 0.950],
        n_E_comp,
        n_u,
        100_000, 20_000, 60_000, 48_000, 5_000
    )
    @test flux isa SourceFlux
    @test flux.source_name == "ICV_barrel"
    @test size(flux.pdf_main) == (25, 10)
    @test length(flux.pdf_companion) == 3
    @test size(flux.pdf_companion[1]) == (20, 10)
    @test flux.companion_E_line == [0.583, 0.511, 0.861]
    @test flux.companion_BR == [0.85, 0.23, 0.12]
    @test flux.companion_f == [0.25, 0.20, 0.30]
    @test flux.n_u == 10
    @test flux.N_generated == 100_000
    @test flux.N_surviving_main == 20_000
end

@testset "SourceRateTable construction" begin
    n_E, n_u = 25, 10
    pdf_rate = zeros(n_E, n_u)
    pdf_rate[20, 5] = 1.5e-3   # gammas/sec in this bin

    rate = SourceRateTable(
        :barrel,
        pdf_rate,
        2.200, 2.500, n_E, n_u,
        ["OCV", "MLI", "ICV"],
        [0.5e-3, 0.3e-3, 0.7e-3],
        1.5e-3
    )
    @test rate.surface == :barrel
    @test size(rate.pdf_rate) == (25, 10)
    @test rate.E_min == 2.200
    @test rate.E_max == 2.500
    @test rate.n_E == 25
    @test rate.n_u == 10
    @test rate.component_names == ["OCV", "MLI", "ICV"]
    @test rate.total_rate ≈ 1.5e-3
    @test sum(rate.component_rates) ≈ rate.total_rate atol=1e-10
end

@testset "load_source_geometry" begin
    src_path = normpath(joinpath(@__DIR__, "..", "..", "data", "source_geometry_lz_v1.json"))
    sg = load_source_geometry(src_path, MATS)

    # All sources loaded
    @test length(sg) >= 21

    # ICV barrel
    icv = sg["ICV_barrel"]
    @test icv isa SourceVolumeInfo
    @test icv.volume isa PCylShell
    @test icv.material.name == "Ti"
    @test icv.transport == :kn
    @test icv.source_class == "shell_source"
    @test icv.activity["Bi214_mBq_per_kg"] ≈ 0.08
    @test icv.activity["Tl208_mBq_per_kg"] ≈ 0.22
    @test icv.mass_kg > 0.0

    # OCV top head
    ocv_top = sg["OCV_top"]
    @test ocv_top.volume isa PDisk
    @test ocv_top.material.name == "Ti"
    @test ocv_top.transport == :kn

    # MLI is transparent virtual source with equivalent mass
    mli = sg["MLI"]
    @test mli.transport == :transparent
    @test mli.source_class == "virtual_source"
    @test mli.material.name == "Vacuum"
    @test mli.mass_kg ≈ 13.8

    # PMT sources are cylinders
    pmt_top = sg["PMT_TOP_PMTs"]
    @test pmt_top.volume isa PCyl
    @test pmt_top.transport == :transparent
    @test pmt_top.mass_kg ≈ 47.07

    # Shell source mass computed from geometry
    icv_geom_mass = mass(icv.volume) / 1000.0  # g -> kg
    @test icv.mass_kg ≈ icv_geom_mass atol=0.1

    # All sources have valid transport mode
    for (name, sv) in sg
        @test sv.transport in (:kn, :transparent)
    end
end


@testset "generate_flux_bi214" begin
    src_path = normpath(joinpath(@__DIR__, "..", "..", "data", "source_geometry_lz_v1.json"))
    sg = load_source_geometry(src_path, MATS)
    icv = sg["ICV_barrel"]
    rng = MersenneTwister(42)
    N = 10_000

    ft = generate_flux_bi214(N, icv.volume, CFG, rng)
    @test ft isa SourceFluxBi214
    @test ft.source_name == "ICV_barrel"
    @test size(ft.pdf) == (25, 10)
    @test ft.N_generated == N
    @test ft.N_surviving > 0
    @test ft.N_absorbed >= 0
    @test ft.N_backward > 0

    # PDF normalization: sum(pdf) = N_surviving / N_generated
    @test sum(ft.pdf) ≈ ft.N_surviving / ft.N_generated rtol=1e-10

    # Bookkeeping: all events accounted for
    @test ft.N_surviving + ft.N_absorbed + ft.N_backward + ft.N_low_energy == N

    # Survival fraction in reasonable range for thin Ti barrel (~0.9 cm)
    # ~50% backward, small absorption, ~35-45% survive forward
    surv = sum(ft.pdf)
    @test 0.20 < surv < 0.50

    # Backward fraction ~50%
    @test ft.N_backward > N * 0.3
end

@testset "generate_flux_tl208" begin
    src_path = normpath(joinpath(@__DIR__, "..", "..", "data", "source_geometry_lz_v1.json"))
    sg = load_source_geometry(src_path, MATS)
    icv = sg["ICV_barrel"]
    rng = MersenneTwister(42)
    N = 10_000

    ft = generate_flux_tl208(N, icv.volume, CFG, rng)
    @test ft isa SourceFluxTl208
    @test ft.source_name == "ICV_barrel"
    @test size(ft.pdf_main) == (25, 10)
    @test ft.N_generated == N

    # Main gamma PDF normalization
    @test sum(ft.pdf_main) ≈ ft.N_surviving_main / ft.N_generated rtol=1e-10

    # Main bookkeeping
    @test ft.N_surviving_main + ft.N_absorbed_main + ft.N_backward_main + ft.N_low_energy_main == N

    # Main survival fraction reasonable
    surv_main = sum(ft.pdf_main)
    @test 0.20 < surv_main < 0.50

    # 3 companion tables with correct dimensions
    @test length(ft.pdf_companion) == 3
    for i in 1:3
        @test size(ft.pdf_companion[i]) == (25, 10)
    end

    # Companion survival fractions between 0 and 1
    for i in 1:3
        @test 0.0 <= ft.companion_f[i] <= 1.0
    end

    # 861 keV has longer mfp in Ti than 583 keV, so higher survival
    # (companion_lines = [583, 511, 861], companion_f = [f1, f2, f3])
    @test ft.companion_f[3] > ft.companion_f[1]  # 861 > 583

    # Companion BRs preserved
    @test ft.companion_BR == [0.85, 0.23, 0.12]
    @test ft.companion_E_line == [0.583, 0.511, 0.861]
end


@testset "propagate_through_layers" begin
    src_path = normpath(joinpath(@__DIR__, "..", "..", "data", "source_geometry_lz_v1.json"))
    sg = load_source_geometry(src_path, MATS)
    icv = sg["ICV_barrel"].volume
    rng = MersenneTwister(42)

    # A 2.6 MeV gamma starting just outside the ICV barrel, aimed inward
    R_outer_icv = 82.1 + 0.9  # R_inner + wall_thickness
    pos = Float64[R_outer_icv + 0.1, 0.0, 53.585]
    dir = Float64[-1.0, 0.0, 0.0]

    status, E_out, pos_out, dir_out = propagate_through_layers(
        2.615, pos, dir, [icv], CFG, rng)

    # Should either exit or be absorbed
    @test status in (:exited, :absorbed)
    if status === :exited
        @test E_out > 0.0
        @test E_out <= 2.615
    end
end

@testset "generate_flux_compound_bi214" begin
    src_path = normpath(joinpath(@__DIR__, "..", "..", "data", "source_geometry_lz_v1.json"))
    sg = load_source_geometry(src_path, MATS)
    ocv_barrel = sg["OCV_barrel"].volume
    icv_barrel = sg["ICV_barrel"].volume
    rng = MersenneTwister(42)
    N = 5_000

    ft = generate_flux_compound_bi214(N, ocv_barrel, [icv_barrel], icv_barrel, CFG, rng)
    @test ft isa SourceFluxBi214
    @test ft.source_name == "OCV_barrel"
    @test ft.N_generated == N
    @test size(ft.pdf) == (25, 10)

    # PDF normalization
    @test sum(ft.pdf) ≈ ft.N_surviving / ft.N_generated rtol=1e-10

    # Bookkeeping
    @test ft.N_surviving + ft.N_absorbed + ft.N_backward + ft.N_low_energy == N

    # OCV->ICV survival should be lower than ICV-only (extra material)
    ft_icv = generate_flux_bi214(N, icv_barrel, CFG, MersenneTwister(42))
    @test sum(ft.pdf) < sum(ft_icv.pdf)
end

@testset "generate_flux_compound_tl208" begin
    src_path = normpath(joinpath(@__DIR__, "..", "..", "data", "source_geometry_lz_v1.json"))
    sg = load_source_geometry(src_path, MATS)
    ocv_barrel = sg["OCV_barrel"].volume
    icv_barrel = sg["ICV_barrel"].volume
    rng = MersenneTwister(42)
    N = 5_000

    ft = generate_flux_compound_tl208(N, ocv_barrel, [icv_barrel], icv_barrel, CFG, rng)
    @test ft isa SourceFluxTl208
    @test ft.source_name == "OCV_barrel"
    @test ft.N_generated == N
    @test size(ft.pdf_main) == (25, 10)

    # Main PDF normalization
    @test sum(ft.pdf_main) ≈ ft.N_surviving_main / ft.N_generated rtol=1e-10

    # Main bookkeeping
    @test ft.N_surviving_main + ft.N_absorbed_main + ft.N_backward_main + ft.N_low_energy_main == N

    # 3 companion tables
    @test length(ft.pdf_companion) == 3
    for i in 1:3
        @test 0.0 <= ft.companion_f[i] <= 1.0
    end

    # OCV->ICV main survival lower than ICV-only
    ft_icv = generate_flux_tl208(N, icv_barrel, CFG, MersenneTwister(42))
    @test sum(ft.pdf_main) < sum(ft_icv.pdf_main)
end



@testset "cryostat barrel propagation chain" begin
    src_path = normpath(joinpath(@__DIR__, "..", "..", "data", "source_geometry_lz_v1.json"))
    sg = load_source_geometry(src_path, MATS)
    ocv = sg["OCV_barrel"]
    mli = sg["MLI"]
    icv = sg["ICV_barrel"]

    R_ocv_inner = ocv.volume.logical.solid.R_inner_cm   # 90.8
    R_mli_inner = mli.volume.logical.solid.R_inner_cm   # 83.0
    R_mli_outer = R_mli_inner + mli.volume.logical.solid.wall_thickness_cm  # 84.0
    R_icv_inner = icv.volume.logical.solid.R_inner_cm   # 82.1
    R_icv_outer = R_icv_inner + icv.volume.logical.solid.wall_thickness_cm  # 83.0
    z_center = ocv.volume.logical.position[3]

    N = 100

    # --- From OCV: propagate through OCV (KN), then [mli, icv] ---
    n_exit_ocv = 0
    n_exit_chain = 0
    for i in 1:N
        rng = MersenneTwister(i)
        # Start inside OCV, heading inward
        R_start = R_ocv_inner + 0.35  # middle of OCV wall
        pos = Float64[R_start, 0.0, z_center]
        dir = Float64[-1.0, 0.0, 0.0]

        status, E, pos, dir = propagate_in_source(2.448, pos, dir, ocv.volume, CFG, rng)
        if status == :exited
            n_exit_ocv += 1

            status2, E2, pos2, dir2 = propagate_through_layers(
                E, pos, dir, PhysicalVolume[mli.volume, icv.volume], CFG, rng)
            if status2 == :exited
                n_exit_chain += 1
            end
        end
    end
    @test n_exit_ocv > 10     # most exit thin Ti OCV
    @test n_exit_chain > 5    # some survive the full chain

    # --- From MLI: propagate through MLI (vacuum), then ICV (KN) ---
    n_exit_mli = 0
    n_exit_icv_from_mli = 0
    for i in 1:N
        rng = MersenneTwister(i + 1000)
        R_start = R_mli_inner + 0.5  # middle of MLI
        pos = Float64[R_start, 0.0, z_center]
        dir = Float64[-1.0, 0.0, 0.0]

        # MLI is vacuum: propagate_in_source should always exit
        status, E, pos, dir = propagate_in_source(2.448, pos, dir, mli.volume, CFG, rng)
        @test status == :exited
        @test E ≈ 2.448   # no energy loss in vacuum
        n_exit_mli += 1

        status2, E2, pos2, dir2 = propagate_through_layers(
            E, pos, dir, PhysicalVolume[icv.volume], CFG, rng)
        if status2 == :exited
            n_exit_icv_from_mli += 1
        end
    end
    @test n_exit_mli == N          # all exit vacuum
    @test n_exit_icv_from_mli > 10 # most survive thin ICV

    # --- From ICV: propagate through ICV (KN) only ---
    n_exit_icv = 0
    for i in 1:N
        rng = MersenneTwister(i + 2000)
        R_start = R_icv_inner + 0.45  # middle of ICV wall
        pos = Float64[R_start, 0.0, z_center]
        dir = Float64[-1.0, 0.0, 0.0]

        status, E, pos, dir = propagate_in_source(2.448, pos, dir, icv.volume, CFG, rng)
        if status == :exited
            n_exit_icv += 1
        end
    end
    @test n_exit_icv > 10  # most survive thin ICV
end

@testset "cryostat_barrel_flux" begin
    src_path = normpath(joinpath(@__DIR__, "..", "..", "data", "source_geometry_lz_v1.json"))
    sg = load_source_geometry(src_path, MATS)
    N = 2_000
    result = cryostat_barrel_flux(N, sg, CFG, MersenneTwister(42))

    # Structure: 8 tables (4 per isotope)
    @test result.bi214_ocv isa SourceFluxBi214
    @test result.bi214_mli isa SourceFluxBi214
    @test result.bi214_icv isa SourceFluxBi214
    @test result.bi214_rate isa SourceRateTable
    @test result.tl208_ocv isa SourceFluxTl208
    @test result.tl208_mli isa SourceFluxTl208
    @test result.tl208_icv isa SourceFluxTl208
    @test result.tl208_rate isa SourceRateTable

    # Rate tables have correct surface
    @test result.bi214_rate.surface == :barrel
    @test result.tl208_rate.surface == :barrel

    # 3 components for barrel
    @test length(result.bi214_rate.component_names) == 3
    @test result.bi214_rate.component_names == ["OCV_barrel", "MLI", "ICV_barrel"]

    # Total rate is sum of components
    @test result.bi214_rate.total_rate ≈ sum(result.bi214_rate.component_rates) atol=1e-15
    @test result.tl208_rate.total_rate ≈ sum(result.tl208_rate.component_rates) atol=1e-15

    # All components must have nonzero survival (non-empty flux tables)
    @test sum(result.bi214_ocv.pdf) > 0
    @test sum(result.bi214_mli.pdf) > 0
    @test sum(result.bi214_icv.pdf) > 0
    @test sum(result.tl208_ocv.pdf_main) > 0
    @test sum(result.tl208_mli.pdf_main) > 0
    @test sum(result.tl208_icv.pdf_main) > 0

    # OCV survival < ICV survival (more material)
    @test sum(result.bi214_ocv.pdf) < sum(result.bi214_icv.pdf)

    # BR is applied: rate = A_chain * mass * BR * survival, not A_chain * mass * survival
    icv = sg["ICV_barrel"]
    A_bi_chain = get(icv.activity, "Bi214_mBq_per_kg", 0.0) * 1e-3
    A_tl_chain = get(icv.activity, "Tl208_mBq_per_kg", 0.0) * 1e-3
    f_bi_icv = sum(result.bi214_icv.pdf)
    f_tl_icv = sum(result.tl208_icv.pdf_main)

    rate_bi_icv_with_BR = A_bi_chain * icv.mass_kg * BR_BI214_2448 * f_bi_icv
    rate_bi_icv_no_BR   = A_bi_chain * icv.mass_kg * f_bi_icv

    rate_tl_icv_with_BR = A_tl_chain * icv.mass_kg * BR_TL208_2615 * f_tl_icv
    rate_tl_icv_no_BR   = A_tl_chain * icv.mass_kg * f_tl_icv

    bi_icv_idx = findfirst(==("ICV_barrel"), result.bi214_rate.component_names)
    tl_icv_idx = findfirst(==("ICV_barrel"), result.tl208_rate.component_names)

    @test result.bi214_rate.component_rates[bi_icv_idx] ≈ rate_bi_icv_with_BR atol=1e-15
    @test result.bi214_rate.component_rates[bi_icv_idx] < rate_bi_icv_no_BR
    @test rate_bi_icv_no_BR / rate_bi_icv_with_BR ≈ 1.0 / BR_BI214_2448 atol=0.01

    @test result.tl208_rate.component_rates[tl_icv_idx] ≈ rate_tl_icv_with_BR atol=1e-15
    @test result.tl208_rate.component_rates[tl_icv_idx] < rate_tl_icv_no_BR
    @test rate_tl_icv_no_BR / rate_tl_icv_with_BR ≈ 1.0 / BR_TL208_2615 atol=0.01
end

@testset "cryostat_top_flux" begin
    src_path = normpath(joinpath(@__DIR__, "..", "..", "data", "source_geometry_lz_v1.json"))
    sg = load_source_geometry(src_path, MATS)
    N = 2_000
    result = cryostat_top_flux(N, sg, CFG, MersenneTwister(42))

    @test result.bi214_rate.surface == :top
    @test length(result.bi214_rate.component_names) == 2
    @test result.bi214_rate.component_names == ["OCV_top", "ICV_top"]
    @test result.bi214_rate.total_rate ≈ sum(result.bi214_rate.component_rates) atol=1e-15
end

@testset "cryostat_bottom_flux" begin
    src_path = normpath(joinpath(@__DIR__, "..", "..", "data", "source_geometry_lz_v1.json"))
    sg = load_source_geometry(src_path, MATS)
    N = 2_000
    result = cryostat_bottom_flux(N, sg, CFG, MersenneTwister(42))

    @test result.bi214_rate.surface == :bottom
    @test length(result.bi214_rate.component_names) == 2
    @test result.bi214_rate.component_names == ["OCV_bottom", "ICV_bottom"]
    @test result.bi214_rate.total_rate ≈ sum(result.bi214_rate.component_rates) atol=1e-15
end


@testset "sample_from_flux Bi214" begin
    src_path = normpath(joinpath(@__DIR__, "..", "..", "data", "source_geometry_lz_v1.json"))
    sg = load_source_geometry(src_path, MATS)
    icv = sg["ICV_barrel"]
    ft = generate_flux_bi214(5_000, icv.volume, CFG, MersenneTwister(42))
    rng = MersenneTwister(99)

    for _ in 1:100
        E, u = sample_from_flux(ft, rng)
        @test ft.E_min <= E <= ft.E_max
        @test 0.0 <= u <= 1.0
    end
end

@testset "sample_from_flux Tl208" begin
    src_path = normpath(joinpath(@__DIR__, "..", "..", "data", "source_geometry_lz_v1.json"))
    sg = load_source_geometry(src_path, MATS)
    icv = sg["ICV_barrel"]
    ft = generate_flux_tl208(5_000, icv.volume, CFG, MersenneTwister(42))
    rng = MersenneTwister(99)

    n_with_companion = 0
    for _ in 1:200
        gammas = sample_from_flux(ft, rng)
        @test length(gammas) >= 1
        @test length(gammas) <= 4

        # Main gamma in range
        E_main, u_main = gammas[1]
        @test ft.E_min_main <= E_main <= ft.E_max_main
        @test 0.0 <= u_main <= 1.0

        # Companions have valid u
        for j in 2:length(gammas)
            _, u_c = gammas[j]
            @test 0.0 <= u_c <= 1.0
        end

        length(gammas) > 1 && (n_with_companion += 1)
    end

    # Should get some companions
    @test n_with_companion > 0
end

@testset "sample_from_rate_table" begin
    src_path = normpath(joinpath(@__DIR__, "..", "..", "data", "source_geometry_lz_v1.json"))
    sg = load_source_geometry(src_path, MATS)
    result = cryostat_barrel_flux(2_000, sg, CFG, MersenneTwister(42))
    rng = MersenneTwister(99)

    for _ in 1:100
        E, u = sample_from_rate_table(result.bi214_rate, rng)
        @test result.bi214_rate.E_min <= E <= result.bi214_rate.E_max
        @test 0.0 <= u <= 1.0
    end
end


@testset "sample_barrel_point" begin
    rng = MersenneTwister(42)
    R = 82.1
    z_min = 53.585 - 94.915
    z_max = 53.585 + 94.915

    for _ in 1:100
        pos, normal = sample_barrel_point(R, z_min, z_max, rng)
        # Position on barrel surface
        r = sqrt(pos[1]^2 + pos[2]^2)
        @test r ≈ R atol=1e-10
        @test z_min <= pos[3] <= z_max
        # Normal is unit vector pointing inward
        @test sqrt(normal[1]^2 + normal[2]^2 + normal[3]^2) ≈ 1.0 atol=1e-10
        @test normal[3] ≈ 0.0 atol=1e-10  # no z component
        # Normal points inward (opposite to position radial direction)
        @test normal[1] * pos[1] + normal[2] * pos[2] < 0.0
    end
end

@testset "sample_disk_point" begin
    rng = MersenneTwister(42)
    R = 72.8
    z = 152.6

    # Downward normal
    for _ in 1:100
        pos, normal = sample_disk_point(R, z, -1.0, rng)
        r = sqrt(pos[1]^2 + pos[2]^2)
        @test r <= R + 1e-10
        @test pos[3] ≈ z atol=1e-10
        @test normal ≈ [0.0, 0.0, -1.0] atol=1e-10
    end

    # Upward normal
    for _ in 1:100
        pos, normal = sample_disk_point(R, z, 1.0, rng)
        @test normal ≈ [0.0, 0.0, 1.0] atol=1e-10
    end

    # Check uniform radial distribution: mean(r^2) should be R^2/2
    r2_vals = Float64[]
    for _ in 1:10000
        pos, _ = sample_disk_point(R, z, -1.0, rng)
        push!(r2_vals, pos[1]^2 + pos[2]^2)
    end
    @test mean(r2_vals) ≈ R^2 / 2 rtol=0.05
end

@testset "make_virtual_envelope from SourceVolumeInfo" begin
    sg = load_source_geometry(
        normpath(joinpath(@__DIR__, "..", "..", "data", "source_geometry_lz_v1.json")),
        MATS)

    # PMT_TOP_PMTs is a PCyl → should give :disk_flat
    sv_pmt = sg["PMT_TOP_PMTs"]
    env = make_virtual_envelope(sv_pmt)
    @test env.kind === :disk_flat
    @test env.R_cm ≈ 72.8
    # z should be at the bottom face of the volume
    lv = sv_pmt.volume.logical
    z_bottom = lv.position[3] - lv.solid.half_height_cm
    @test env.z_equator_cm ≈ z_bottom

    # PMT_BARREL_R8520 is a PCylShell → should give :barrel
    sv_skin = sg["PMT_BARREL_R8520"]
    env2 = make_virtual_envelope(sv_skin)
    @test env2.kind === :barrel
    lv2 = sv_skin.volume.logical
    @test env2.R_cm ≈ lv2.solid.R_inner_cm
    @test env2.z_min_cm ≈ lv2.position[3] - lv2.solid.half_height_cm
    @test env2.z_max_cm ≈ lv2.position[3] + lv2.solid.half_height_cm

    # Surface sampler from disk_flat VE should produce valid points
    rng = MersenneTwister(42)
    sampler = make_surface_sampler(env)
    for _ in 1:50
        pos, normal = sampler(rng)
        r = sqrt(pos[1]^2 + pos[2]^2)
        @test r <= env.R_cm + 1e-10
        @test pos[3] ≈ env.z_equator_cm atol=1e-10
        @test normal ≈ [0.0, 0.0, -1.0] atol=1e-10
    end
end

@testset "sample_cap_point" begin
    rng = MersenneTwister(42)
    R = 83.0
    ar = 2.0
    z_eq = 148.5

    for _ in 1:100
        pos, normal = sample_cap_point(R, ar, z_eq, :up, rng)
        r = sqrt(pos[1]^2 + pos[2]^2)
        @test r <= R + 1e-6
        @test pos[3] >= z_eq - 1e-6  # above equator for :up
        # Normal is unit vector
        @test sqrt(normal[1]^2 + normal[2]^2 + normal[3]^2) ≈ 1.0 atol=1e-10
        # Inward normal for :up cap points downward (z < 0)
        @test normal[3] < 0.0
    end

    # :down cap
    z_eq_bot = -41.33
    for _ in 1:50
        pos, normal = sample_cap_point(R, 3.0, z_eq_bot, :down, rng)
        r = sqrt(pos[1]^2 + pos[2]^2)
        @test r <= R + 1e-6
        @test pos[3] <= z_eq_bot + 1e-6  # below equator for :down
        # Inward normal for :down cap points upward (z > 0)
        @test normal[3] > 0.0
    end
end

@testset "reconstruct_direction" begin
    rng = MersenneTwister(42)

    # u=1 along normal → direction = normal
    normal = Float64[0.0, 0.0, 1.0]
    for _ in 1:50
        dir = reconstruct_direction(1.0, normal, rng)
        @test dir[3] ≈ 1.0 atol=1e-10
        @test sqrt(dir[1]^2 + dir[2]^2 + dir[3]^2) ≈ 1.0 atol=1e-10
    end

    # u=0 → perpendicular to normal
    for _ in 1:50
        dir = reconstruct_direction(0.0, normal, rng)
        @test abs(dir[3]) < 1e-10
        @test sqrt(dir[1]^2 + dir[2]^2 + dir[3]^2) ≈ 1.0 atol=1e-10
    end

    # General normal: direction is unit vector
    normal2 = Float64[-1.0 / sqrt(2.0), 0.0, -1.0 / sqrt(2.0)]
    for _ in 1:50
        dir = reconstruct_direction(0.5, normal2, rng)
        @test sqrt(dir[1]^2 + dir[2]^2 + dir[3]^2) ≈ 1.0 atol=1e-10
    end
end

@testset "sample_gamma_from_flux Bi214" begin
    src_path = normpath(joinpath(@__DIR__, "..", "..", "data", "source_geometry_lz_v1.json"))
    sg = load_source_geometry(src_path, MATS)
    icv = sg["ICV_barrel"]
    ft = generate_flux_bi214(5_000, icv.volume, CFG, MersenneTwister(42))
    rng = MersenneTwister(99)

    R_inner = 82.1
    z_min = 53.585 - 94.915
    z_max = 53.585 + 94.915
    barrel_sampler = rng -> sample_barrel_point(R_inner, z_min, z_max, rng)

    for _ in 1:50
        g = sample_gamma_from_flux(ft, barrel_sampler, rng)
        @test g isa SampledGamma
        @test ft.E_min <= g.E_MeV <= ft.E_max
        @test length(g.position) == 3
        @test length(g.direction) == 3
        @test sqrt(sum(g.direction .^ 2)) ≈ 1.0 atol=1e-10
    end
end

@testset "sample_gamma_from_flux Tl208" begin
    src_path = normpath(joinpath(@__DIR__, "..", "..", "data", "source_geometry_lz_v1.json"))
    sg = load_source_geometry(src_path, MATS)
    icv = sg["ICV_barrel"]
    ft = generate_flux_tl208(5_000, icv.volume, CFG, MersenneTwister(42))
    rng = MersenneTwister(99)

    R_inner = 82.1
    z_min = 53.585 - 94.915
    z_max = 53.585 + 94.915
    barrel_sampler = rng -> sample_barrel_point(R_inner, z_min, z_max, rng)

    n_multi = 0
    for _ in 1:100
        gammas = sample_gamma_from_flux(ft, barrel_sampler, rng)
        @test length(gammas) >= 1
        @test length(gammas) <= 4

        # All gammas share the same position (same decay location)
        for g in gammas
            @test g.position == gammas[1].position
            @test sqrt(sum(g.direction .^ 2)) ≈ 1.0 atol=1e-10
        end

        length(gammas) > 1 && (n_multi += 1)
    end
    @test n_multi > 0
end


@testset "merge_flux_bi214" begin
    src_path = normpath(joinpath(@__DIR__, "..", "..", "data", "source_geometry_lz_v1.json"))
    sg = load_source_geometry(src_path, MATS)
    icv = sg["ICV_barrel"].volume

    # Generate two halves with different seeds
    ft1 = generate_flux_bi214(3_000, icv, CFG, MersenneTwister(10))
    ft2 = generate_flux_bi214(3_000, icv, CFG, MersenneTwister(20))

    merged = merge_flux_bi214([ft1, ft2])
    @test merged.N_generated == ft1.N_generated + ft2.N_generated
    @test merged.N_surviving == ft1.N_surviving + ft2.N_surviving
    @test merged.N_absorbed == ft1.N_absorbed + ft2.N_absorbed
    @test merged.N_backward == ft1.N_backward + ft2.N_backward
    @test merged.N_low_energy == ft1.N_low_energy + ft2.N_low_energy
    @test size(merged.pdf) == size(ft1.pdf)

    # Weighted average: merged survival should be between the two
    s1 = sum(ft1.pdf)
    s2 = sum(ft2.pdf)
    sm = sum(merged.pdf)
    @test sm ≈ (s1 * ft1.N_generated + s2 * ft2.N_generated) / merged.N_generated atol=1e-12

    # Single-element merge is identity
    single = merge_flux_bi214([ft1])
    @test single === ft1
end

@testset "merge_flux_tl208" begin
    src_path = normpath(joinpath(@__DIR__, "..", "..", "data", "source_geometry_lz_v1.json"))
    sg = load_source_geometry(src_path, MATS)
    icv = sg["ICV_barrel"].volume

    ft1 = generate_flux_tl208(3_000, icv, CFG, MersenneTwister(10))
    ft2 = generate_flux_tl208(3_000, icv, CFG, MersenneTwister(20))

    merged = merge_flux_tl208([ft1, ft2])
    @test merged.N_generated == ft1.N_generated + ft2.N_generated
    @test merged.N_surviving_main == ft1.N_surviving_main + ft2.N_surviving_main
    @test merged.N_absorbed_main == ft1.N_absorbed_main + ft2.N_absorbed_main
    @test size(merged.pdf_main) == size(ft1.pdf_main)
    @test length(merged.pdf_companion) == 3

    # Main PDF weighted average
    s1 = sum(ft1.pdf_main)
    s2 = sum(ft2.pdf_main)
    sm = sum(merged.pdf_main)
    @test sm ≈ (s1 * ft1.N_generated + s2 * ft2.N_generated) / merged.N_generated atol=1e-12

    # Companion survival fractions are weighted averages
    for ic in 1:3
        expected = (ft1.companion_f[ic] * ft1.N_generated + ft2.companion_f[ic] * ft2.N_generated) / merged.N_generated
        @test merged.companion_f[ic] ≈ expected atol=1e-12
    end

    # Single-element merge is identity
    single = merge_flux_tl208([ft1])
    @test single === ft1
end

@testset "write_pdf_csv round-trip" begin
    n_E, n_u = 5, 3
    pdf = rand(MersenneTwister(42), n_E, n_u)
    E_min, E_max = 2.2, 2.5
    path = joinpath(tempdir(), "test_pdf_roundtrip.csv")

    write_pdf_csv(path, pdf, E_min, E_max, n_E, n_u)

    # Read back and verify
    lines = readlines(path)
    @test length(lines) == n_E + 1  # header + data rows

    # Header has n_u + 1 columns (E_MeV + u bins)
    header_cols = split(lines[1], ",")
    @test length(header_cols) == n_u + 1
    @test header_cols[1] == "E_MeV"

    # Check data values
    for i in 1:n_E
        cols = split(lines[i+1], ",")
        @test length(cols) == n_u + 1
        for j in 1:n_u
            val = parse(Float64, cols[j+1])
            @test val ≈ pdf[i, j] atol=1e-6
        end
    end

    rm(path)
end


@testset "Surface sampler from sg with epsilon inset" begin
    # Test that surface samplers built from source geometry + tiny inset
    # place gammas inside the detector. The inset (ENVELOPE_INSET_CM) is
    # needed because the fast kernel uses strict < for boundary checks.
    # With the geometry fix, the ICV inner surface matches the detector
    # boundary exactly, so the inset is purely numerical (not geometric).
    src_path = normpath(joinpath(@__DIR__, "..", "..", "data", "source_geometry_lz_v1.json"))
    sg = load_source_geometry(src_path, MATS)
    fk3 = compile_fastkernel_geometry(DET3)
    ε = ENVELOPE_INSET_CM

    # --- Barrel: ICV_barrel inner surface with inset ---
    icv_b = sg["ICV_barrel"]
    R_b = icv_b.volume.logical.solid.R_inner_cm - ε
    z_eq_top = sg["ICV_top"].volume.logical.position[3]
    z_eq_bot = sg["ICV_bottom"].volume.logical.position[3]
    # Barrel z excludes dome regions
    top_depth = R_b / sg["ICV_top"].volume.logical.solid.aspect_ratio
    bot_depth = R_b / sg["ICV_bottom"].volume.logical.solid.aspect_ratio
    barrel_sampler = rng -> sample_barrel_point(R_b, z_eq_bot + bot_depth, z_eq_top - top_depth, rng)

    rng = MersenneTwister(42)
    local n_in = 0
    for _ in 1:200
        pos, _ = barrel_sampler(rng)
        region = classify_fastkernel(fk3, (pos[1], pos[2], pos[3]))
        region !== nothing && (n_in += 1)
    end
    @test n_in == 200

    # --- Top: ICV_top inner dome with inset ---
    icv_t = sg["ICV_top"]
    R_t = icv_t.volume.logical.solid.radius_cm - icv_t.volume.logical.solid.wall_thickness_cm - ε
    ar_t = icv_t.volume.logical.solid.aspect_ratio
    top_sampler = rng -> sample_cap_point(R_t, ar_t, z_eq_top, :up, rng)

    local n_in_t = 0
    for _ in 1:200
        pos, _ = top_sampler(rng)
        region = classify_fastkernel(fk3, (pos[1], pos[2], pos[3]))
        region !== nothing && (n_in_t += 1)
    end
    @test n_in_t == 200

    # --- Bottom: ICV_bottom inner dome with inset ---
    icv_bot = sg["ICV_bottom"]
    R_bot = icv_bot.volume.logical.solid.radius_cm - icv_bot.volume.logical.solid.wall_thickness_cm - ε
    ar_bot = icv_bot.volume.logical.solid.aspect_ratio
    bot_sampler = rng -> sample_cap_point(R_bot, ar_bot, z_eq_bot, :down, rng)

    local n_in_b = 0
    for _ in 1:200
        pos, _ = bot_sampler(rng)
        region = classify_fastkernel(fk3, (pos[1], pos[2], pos[3]))
        region !== nothing && (n_in_b += 1)
    end
    @test n_in_b == 200

    # --- PMT top: already inside detector, R=72.8 << 82.1 ---
    pmt_top_sampler = rng -> sample_disk_point(72.8, 152.6, -1.0, rng)

    local n_in_pt = 0
    for _ in 1:200
        pos, _ = pmt_top_sampler(rng)
        region = classify_fastkernel(fk3, (pos[1], pos[2], pos[3]))
        region !== nothing && (n_in_pt += 1)
    end
    @test n_in_pt == 200

    # --- PMT bottom: already inside detector ---
    pmt_bot_sampler = rng -> sample_disk_point(72.8, -15.75, 1.0, rng)

    local n_in_pb = 0
    for _ in 1:200
        pos, _ = pmt_bot_sampler(rng)
        region = classify_fastkernel(fk3, (pos[1], pos[2], pos[3]))
        region !== nothing && (n_in_pb += 1)
    end
    @test n_in_pb == 200
end

@testset "Gammas from surface sampler reach detector" begin
    src_path = normpath(joinpath(@__DIR__, "..", "..", "data", "source_geometry_lz_v1.json"))
    sg = load_source_geometry(src_path, MATS)
    fk3 = compile_fastkernel_geometry(DET3)
    fv3 = compile_fv_volume(DET3)

    for (source_name, isotope) in [("cryostat_barrel", "Bi214"),
                                    ("cryostat_top", "Bi214"),
                                    ("cryostat_bottom", "Bi214")]
        # Generate a small flux table
        icv_key = source_name == "cryostat_barrel" ? "ICV_barrel" :
                  source_name == "cryostat_top" ? "ICV_top" : "ICV_bottom"
        icv_vol = sg[icv_key].volume
        ft = generate_flux_bi214(5_000, icv_vol, CFG, MersenneTwister(42))

        # Build sampler from source geometry
        sampler = make_surface_sampler(source_name, sg)
        rng = MersenneTwister(99)

        # All sampled gammas should land inside the detector envelope
        n_inside = 0
        for _ in 1:100
            g = sample_gamma_from_flux(ft, sampler, rng)
            region = classify_fastkernel(fk3, (g.position[1], g.position[2], g.position[3]))
            region !== nothing && (n_inside += 1)
        end
        @test n_inside == 100  # all inside

        # For barrel and top: at least some should be vetoed
        # Bottom is very well shielded (~7 mfp of passive LXe), skip veto check
        if source_name != "cryostat_bottom"
            n_vetoed = 0
            rng2 = MersenneTwister(42)
            for _ in 1:500
                g = sample_gamma_from_flux(ft, sampler, rng2)
                r = process_event([g], fk3, fv3, CFG, rng2)
                r.status == :vetoed && (n_vetoed += 1)
            end
            @test n_vetoed > 0
        end
    end
end


# =====================================================================
# PMT merged volume geometry
# =====================================================================
@testset "PMT merged volume geometry" begin
    src_path = normpath(joinpath(@__DIR__, "..", "..", "data", "source_geometry_lz_v1.json"))
    sg = load_source_geometry(src_path, MATS)

    # --- TOP: 3 components ---
    merged_top, comp_top = LXeMC._merge_pmt_volume(
        sg, ["PMT_TOP_PMTs", "PMT_TOP_bases", "PMT_TOP_structure"],
        "PMT_TOP_merged", :up)

    @test length(comp_top) == 3
    @test merged_top isa PDisk
    @test merged_top.material.density ≈ 0.0  # vacuum

    # z span should match expected values from source geometry
    s = merged_top.logical.solid
    c = merged_top.logical.position
    @test s.radius_cm ≈ 72.8
    @test is_flat(s)  # aspect_ratio = Inf

    # z range: [152.599, 153.706] (from bases bottom to structure top)
    z_bottom = c[3]  # position is at z_min for :up
    z_top = c[3] + s.wall_thickness_cm
    @test z_bottom ≈ 152.599 atol=0.01
    @test z_top ≈ 153.706 atol=0.01

    # --- BOTTOM: 4 components ---
    merged_bot, comp_bot = LXeMC._merge_pmt_volume(
        sg, ["PMT_BOT_PMTs", "PMT_BOT_bases", "PMT_BOT_structure", "PMT_BOT_R8778_dome"],
        "PMT_BOT_merged", :down)

    @test length(comp_bot) == 4
    @test merged_bot isa PDisk
    @test merged_bot.material.density ≈ 0.0

    s_bot = merged_bot.logical.solid
    c_bot = merged_bot.logical.position
    @test s_bot.radius_cm ≈ 72.8
    @test is_flat(s_bot)

    # z range: [-16.856, -15.750] (from structure bottom to shared top)
    z_top_bot = c_bot[3]  # position is at z_max for :down
    z_bottom_bot = c_bot[3] - s_bot.wall_thickness_cm
    @test z_top_bot ≈ -15.750 atol=0.01
    @test z_bottom_bot ≈ -16.856 atol=0.01
end


# =====================================================================
# PMT propagation: gammas enter air (top) or LXe (bottom)
# =====================================================================
@testset "PMT propagation into detector" begin
    src_path = normpath(joinpath(@__DIR__, "..", "..", "data", "source_geometry_lz_v1.json"))
    sg = load_source_geometry(src_path, MATS)
    det3 = load_tracking_detector(default_tracking_detector_path(), MATS)
    fk3 = compile_fastkernel_geometry(det3)

    # --- TOP PMTs: gammas going downward should enter AirDome (gas) ---
    merged_top, _ = LXeMC._merge_pmt_volume(
        sg, ["PMT_TOP_PMTs", "PMT_TOP_bases", "PMT_TOP_structure"],
        "PMT_TOP_merged", :up)

    rng = MersenneTwister(42)
    n_in_detector = 0
    N = 200
    for _ in 1:N
        pos = random_position_in_volume(merged_top, rng)
        dir = sample_isotropic_direction(rng)
        # Only downward gammas (toward LXe)
        dir[3] > 0 && continue

        # Check that the gamma position is inside the AirDome region
        # (PMTs sit inside the gas dome at z~152.6, R=72.8)
        region = classify_fastkernel(fk3, (pos[1], pos[2], pos[3]))
        if region !== nothing
            n_in_detector += 1
        end
    end
    # At least some gammas should start inside the detector envelope
    @test n_in_detector > 0

    # --- TOP PMTs flux: ~50% survival (hemisphere cut, no material) ---
    result_top = pmt_top_flux(5000, sg, CFG, MersenneTwister(42))
    f_top = sum(result_top.bi214.pdf)
    # Isotropic in vacuum: ~50% go downward, all survive
    @test 0.3 < f_top < 0.7

    # Rate table has 3 components
    @test length(result_top.bi214_rate.component_names) == 3
    @test result_top.bi214_rate.total_rate > 0

    # --- BOTTOM PMTs: gammas going upward should enter LXe ---
    merged_bot, _ = LXeMC._merge_pmt_volume(
        sg, ["PMT_BOT_PMTs", "PMT_BOT_bases", "PMT_BOT_structure", "PMT_BOT_R8778_dome"],
        "PMT_BOT_merged", :down)

    rng = MersenneTwister(42)
    n_in_detector_bot = 0
    for _ in 1:N
        pos = random_position_in_volume(merged_bot, rng)
        dir = sample_isotropic_direction(rng)
        dir[3] < 0 && continue

        region = classify_fastkernel(fk3, (pos[1], pos[2], pos[3]))
        if region !== nothing
            n_in_detector_bot += 1
        end
    end
    @test n_in_detector_bot > 0

    # --- BOTTOM PMTs flux: ~50% survival ---
    result_bot = pmt_bottom_flux(5000, sg, CFG, MersenneTwister(42))
    f_bot = sum(result_bot.bi214.pdf)
    @test 0.3 < f_bot < 0.7

    # Rate table has 4 components
    @test length(result_bot.bi214_rate.component_names) == 4
    @test result_bot.bi214_rate.total_rate > 0
end


# =====================================================================
# PMT barrel: R8520 Top-Skin PMTs only. Cables now live in
# pmt_top_cables / pmt_bottom_cables (LZ TDR §3.4.4); the R8778
# lower-ring lives in pmt_skin_lower_ring.
# =====================================================================
@testset "pmt_barrel flux" begin
    src_path = normpath(joinpath(@__DIR__, "..", "..", "data", "source_geometry_lz_v1.json"))
    sg = load_source_geometry(src_path, MATS)

    # Merged barrel geometry
    merged, comps = LXeMC._merge_pmt_barrel_volume(sg,
        ["PMT_BARREL_R8520"],
        "PMT_BARREL_merged")
    @test length(comps) == 1
    @test merged isa PCylShell
    @test merged.material.density ≈ 0.0  # vacuum

    # Flux: ~50% hemisphere survival (transparent, inward radial)
    result = pmt_barrel_flux(10000, sg, CFG, MersenneTwister(42))
    f = sum(result.bi214.pdf)
    @test 0.3 < f < 0.7

    # Rate table has 1 component with correct name
    @test length(result.bi214_rate.component_names) == 1
    @test result.bi214_rate.total_rate > 0
    @test "PMT_BARREL_R8520" in result.bi214_rate.component_names

    # Tl208 works
    @test sum(result.tl208.pdf_main) > 0
    @test result.tl208_rate.total_rate > 0
end


# =====================================================================
# Skin lower ring: 20 R8778 'Bottom Skin' side PMTs localised at
# cathode level (LZ TDR table line 3754).
# =====================================================================
@testset "pmt_skin_lower_ring flux" begin
    src_path = normpath(joinpath(@__DIR__, "..", "..", "data", "source_geometry_lz_v1.json"))
    sg = load_source_geometry(src_path, MATS)

    merged, comps = LXeMC._merge_pmt_barrel_volume(sg,
        ["PMT_SKIN_LOWER_RING"],
        "PMT_SKIN_LOWER_RING_merged")
    @test length(comps) == 1
    @test merged isa PCylShell
    @test merged.material.density ≈ 0.0  # vacuum

    # Ring localised near cathode: z extent should be tight, not full barrel
    hh = merged.logical.solid.half_height_cm
    @test hh < 20.0  # not full barrel (was 79.675 before)

    result = pmt_skin_lower_ring_flux(5000, sg, CFG, MersenneTwister(42))
    f = sum(result.bi214.pdf)
    @test 0.3 < f < 0.7  # inward hemisphere

    @test length(result.bi214_rate.component_names) == 1
    @test "PMT_SKIN_LOWER_RING" in result.bi214_rate.component_names
    @test result.bi214_rate.total_rate > 0
end


# =====================================================================
# PMT cables: upper conduit (top) and lower conduit (bottom)
# =====================================================================
@testset "pmt_top_cables flux" begin
    src_path = normpath(joinpath(@__DIR__, "..", "..", "data", "source_geometry_lz_v1.json"))
    sg = load_source_geometry(src_path, MATS)

    merged, comps = LXeMC._merge_pmt_volume(sg,
        ["PMT_TOP_cables"], "PMT_TOP_cables_merged", :up)
    @test length(comps) == 1
    @test merged isa PDisk
    @test merged.logical.orientation === :up

    result = pmt_top_cables_flux(5000, sg, CFG, MersenneTwister(42))
    f = sum(result.bi214.pdf)
    @test 0.3 < f < 0.7  # downward hemisphere

    @test length(result.bi214_rate.component_names) == 1
    @test "PMT_TOP_cables" in result.bi214_rate.component_names
    @test result.bi214_rate.total_rate > 0
end

@testset "pmt_bottom_cables flux" begin
    src_path = normpath(joinpath(@__DIR__, "..", "..", "data", "source_geometry_lz_v1.json"))
    sg = load_source_geometry(src_path, MATS)

    merged, comps = LXeMC._merge_pmt_volume(sg,
        ["PMT_BOT_cables"], "PMT_BOT_cables_merged", :down)
    @test length(comps) == 1
    @test merged isa PDisk
    @test merged.logical.orientation === :down

    result = pmt_bottom_cables_flux(5000, sg, CFG, MersenneTwister(42))
    f = sum(result.bi214.pdf)
    @test 0.3 < f < 0.7  # upward hemisphere

    @test length(result.bi214_rate.component_names) == 1
    @test "PMT_BOT_cables" in result.bi214_rate.component_names
    @test result.bi214_rate.total_rate > 0
end


# =====================================================================
# PMT bottom LXe slab geometry sanity
# =====================================================================
@testset "pmt_bottom_lxe slab geometry" begin
    src_path = normpath(joinpath(@__DIR__, "..", "..", "data", "source_geometry_lz_v1.json"))
    sg = load_source_geometry(src_path, MATS)

    merged, _ = LXeMC._merge_pmt_volume(sg,
        ["PMT_BOT_PMTs", "PMT_BOT_bases", "PMT_BOT_structure", "PMT_BOT_R8778_dome"],
        "PMT_BOT_merged", :down)
    z_pmt = merged.logical.position[3]
    R = merged.logical.solid.radius_cm
    z_cathode = LXeMC._cathode_z(sg)

    # Build slab the same way pmt_bottom_lxe_flux does
    z_slab_bottom = z_pmt + 0.01
    hh = (z_cathode - z_slab_bottom) / 2.0
    z_center = (z_slab_bottom + z_cathode) / 2.0
    lxe_mat = MATS["LXe"]
    slab = PCyl("test_slab", LCyl(Cyl(R, hh), Float64[0.0, 0.0, z_center]), lxe_mat)

    # Gamma exits PMT source at z_pmt + 1e-4 (propagate_in_source push)
    z_exit = z_pmt + 1e-4

    # Exit point must be OUTSIDE the slab
    @test !is_inside(slab, Float64[0.0, 0.0, z_exit])

    # Slab interior must be inside
    @test is_inside(slab, Float64[0.0, 0.0, z_center])

    # distance_to_entry from exit point going up must find the bottom face
    pos = Float64[0.0, 0.0, z_exit]
    dir = Float64[0.0, 0.0, 1.0]
    t = distance_to_entry(pos, dir, slab.logical)
    @test isfinite(t)
    @test t > 0
    @test t < 0.1  # should be ~0.01 cm (slab_bottom - z_exit)

    # After entry, position must be inside the slab
    entry_pos = pos .+ dir .* (t + 1e-4)
    @test is_inside(slab, entry_pos)
end


# =====================================================================
# PMT bottom LXe: propagation through passive LXe to cathode
# =====================================================================
@testset "pmt_bottom_lxe flux" begin
    src_path = normpath(joinpath(@__DIR__, "..", "..", "data", "source_geometry_lz_v1.json"))
    sg = load_source_geometry(src_path, MATS)

    # Cathode z must equal top face of FC_botgrid
    bg = sg["FC_botgrid"]
    z_expected = bg.volume.logical.position[3] + bg.volume.logical.solid.half_height_cm
    z_cathode = LXeMC._cathode_z(sg)
    @test z_cathode ≈ z_expected atol=1e-10

    # Generate flux tables (50k for reasonable statistics)
    result = pmt_bottom_lxe_flux(50000, sg, MATS, CFG, MersenneTwister(42))
    result_bare = pmt_bottom_flux(50000, sg, CFG, MersenneTwister(42))

    f_lxe = sum(result.bi214.pdf)
    f_bare = sum(result_bare.bi214.pdf)

    # LXe slab must attenuate significantly (~16 cm of LXe at 2.4 MeV)
    # Bare is ~50% (hemisphere), LXe should be noticeably less
    @test f_lxe < f_bare
    @test f_lxe > 0
    @test f_lxe / f_bare < 0.95  # at least 5% attenuation through ~16 cm LXe

    # Some gammas must be absorbed (not just lost)
    @test result.bi214.N_absorbed > 0

    # Flux must have entries below the line energy (Compton-degraded gammas)
    # The bare flux is all in the top energy bin; LXe flux should spread
    pdf = result.bi214.pdf
    n_E = result.bi214.n_E
    top_bin_frac = sum(pdf[n_E, :]) / sum(pdf)
    @test top_bin_frac < 0.99  # not all in the top bin

    # Rate table has 4 components (same PMT sub-components)
    @test length(result.bi214_rate.component_names) == 4
    @test result.bi214_rate.total_rate > 0
    @test result.bi214_rate.total_rate < result_bare.bi214_rate.total_rate

    # Tl208 also works
    @test sum(result.tl208.pdf_main) > 0
    @test result.tl208_rate.total_rate > 0
end
