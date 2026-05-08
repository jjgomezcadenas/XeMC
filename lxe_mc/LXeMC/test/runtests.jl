using Test
using LXeMC
using JSON
using Random
using Statistics

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
    @test names3 == ["AirDome", "BarrelActive", "BottomActive", "FC_PTFE", "FC_rings", "FV", "LXe_passive", "LZ_detector", "MARS", "Skin", "TopActive"]

    root3 = root_node(DET3)
    lz3 = node_by_name(DET3, "LZ_detector")
    air3 = node_by_name(DET3, "AirDome")
    passive3 = node_by_name(DET3, "LXe_passive")
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
    @test find_tracking_node(DET3, (74.45, 0.0, 100.0)).lv.name == "FC_rings"
    @test find_tracking_node(DET3, (78.0, 0.0, 100.0)).lv.name == "Skin"
    @test find_tracking_node(DET3, (0.0, 0.0, -20.0)).lv.name == "LXe_passive"
    @test find_tracking_node(DET3, (0.0, 0.0, 160.0)).lv.name == "AirDome"

    fk3 = compile_fastkernel_geometry(DET3)
    @test sort([r.name for r in fk3.regions]) == ["AirDome", "BarrelActive", "BottomActive", "FC_PTFE", "FC_rings", "FV", "LXe_passive", "LZ_detector", "Skin", "TopActive"]
    @test fk3.top_active_region == fk3.name_to_index["TopActive"]
    @test fk3.barrel_active_region == fk3.name_to_index["BarrelActive"]
    @test fk3.bottom_active_region == fk3.name_to_index["BottomActive"]
    @test fk3.passive_region == fk3.name_to_index["LXe_passive"]
    @test classify_fastkernel(fk3, (0.0, 0.0, 120.0)).name == "TopActive"
    @test classify_fastkernel(fk3, (60.0, 0.0, 61.0)).name == "BarrelActive"
    @test classify_fastkernel(fk3, (0.0, 0.0, 10.0)).name == "BottomActive"
    @test classify_fastkernel(fk3, (0.0, 0.0, 61.0)).name == "FV"
    @test classify_fastkernel(fk3, (0.0, 0.0, -20.0)).name == "LXe_passive"

    fvgeom = compile_fv_geometry(DET3)
    @test fvgeom.radius_cm ≈ 39.0 atol=GEOM_TOL
    @test fvgeom.zmin_cm ≈ 26.0 atol=GEOM_TOL
    @test fvgeom.zmax_cm ≈ 96.0 atol=GEOM_TOL
    @test fvgeom.material.name == "LXe"
    @test is_inside_fv(fvgeom, (0.0, 0.0, 61.0))
    @test !is_inside_fv(fvgeom, (40.0, 0.0, 61.0))
    @test !is_inside_fv(fvgeom, (0.0, 0.0, 20.0))
    @test !is_inside_fv(fvgeom, (0.0, 0.0, 100.0))
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
    @test fk_air.region in ("AirDome", "TopActive", "BarrelActive", "BottomActive", "FV", "Skin", "FC_PTFE", "FC_rings", "LXe_passive", "MARS")

    ptfe = SampledGamma(2.615, Float64[73.0, 0.0, 100.0], Float64[1.0, 0.0, 0.0])
    fk_ptfe = propagate_gamma_fastkernel(ptfe, fk, CFG, MersenneTwister(11))
    @test fk_ptfe.status == :interacted
    @test fk_ptfe.region in ("FC_PTFE", "FC_rings", "Skin")

    lxe_bulk = SampledGamma(2.615, Float64[0.0, 0.0, -20.0], Float64[0.0, 0.0, 1.0])
    fk_lxe = propagate_gamma_fastkernel(lxe_bulk, fk, CFG, MersenneTwister(19))
    @test fk_lxe.status == :interacted
    @test fk_lxe.region in ("LXe_passive", "BottomActive", "FV")
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
    fvgeom = compile_fv_geometry(DET3)

    empty_fk = process_event(SampledGamma[], fk, fvgeom, CFG, MersenneTwister(1))
    @test empty_fk.status == :no_fv
    @test empty_fk.n_processed == 0

    miss_gamma = SampledGamma(2.615, Float64[120.0, 0.0, 160.0], Float64[0.0, 0.0, -1.0])
    miss_res = LXeMC.transport_gamma_fastkernel(miss_gamma, fk, CFG, MersenneTwister(1))
    @test miss_res.status == :escaped
    @test isempty(miss_res.deposits)

    handoff_gamma = SampledGamma(2.615, Float64[0.0, 0.0, 61.0], Float64[0.0, 0.0, -1.0])
    handoff_res = LXeMC.transport_gamma_fastkernel(handoff_gamma, fk, CFG, MersenneTwister(2))
    @test handoff_res.status == :handoff_fv
    @test handoff_res.terminal_region == "FV"
    @test handoff_res.energy_MeV ≈ 2.615 atol=1e-12

    ss_deps = Deposit[
        Deposit(Float64[0.0, 0.0, 61.0], 1.0, :electron),
        Deposit(Float64[0.0, 0.0, 61.2], 1.0, :electron),
    ]
    ms_deps = Deposit[
        Deposit(Float64[0.0, 0.0, 61.0], 1.0, :electron),
        Deposit(Float64[0.0, 0.0, 62.0], 1.0, :electron),
    ]
    @test LXeMC._classify_fv_stack_result(Deposit[], CFG) == :no_fv
    @test LXeMC._classify_fv_stack_result(ss_deps, CFG) == :accepted
    @test LXeMC._classify_fv_stack_result(ms_deps, CFG) == :ms_rejected

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
        vec_res = process_event(gammas, fk, fvgeom, CFG, rng1)
        fused_res = process_event_fastkernel_calib(fk, fvgeom, CFG, rng2)

        @test fused_res.result.status == vec_res.status
        @test fused_res.result.has_fv == vec_res.has_fv
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
end


@testset "FV-only stack transport" begin
    fvgeom = compile_fv_geometry(DET3)
    rng1 = MersenneTwister(20260508)
    deps1 = propagate_gamma_in_fv(2.61, fvgeom, CFG;
                                  position=(0.0, 0.0, 61.0),
                                  direction=(0.0, 0.0, 1.0),
                                  rng=rng1)
    rng2 = MersenneTwister(20260508)
    deps2 = propagate_gamma_in_fv(2.61, fvgeom, CFG;
                                  position=(0.0, 0.0, 61.0),
                                  direction=(0.0, 0.0, 1.0),
                                  rng=rng2)

    @test length(deps2) == length(deps1)
    @test all(is_inside_fv(fvgeom, d.position) for d in deps1)
    @test sum(d.energy for d in deps2) ≈ sum(d.energy for d in deps1) atol=1e-9
    @test [d.source for d in deps2] == [d.source for d in deps1]

    rng3 = MersenneTwister(12345)
    deps3 = propagate_gamma_in_fv(2.61, fvgeom, CFG;
                                  position=(10.0, 0.0, 80.0),
                                  direction=(0.0, 0.0, -1.0),
                                  rng=rng3)
    rng4 = MersenneTwister(12345)
    deps4 = propagate_gamma_in_fv(2.61, fvgeom, CFG;
                                  position=(10.0, 0.0, 80.0),
                                  direction=(0.0, 0.0, -1.0),
                                  rng=rng4)

    @test length(deps4) == length(deps3)
    @test all(is_inside_fv(fvgeom, d.position) for d in deps3)
    @test sum(d.energy for d in deps4) ≈ sum(d.energy for d in deps3) atol=1e-9
    @test [d.source for d in deps4] == [d.source for d in deps3]
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

