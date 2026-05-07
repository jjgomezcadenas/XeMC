using Test
using LXeMC
using JSON
using Random
using Statistics

const CFG  = default_config()
const MATS = load_materials(CFG)
const DET  = load_detector(default_detector_path(), MATS)
const DET3 = load_tracking_detector(default_tracking_detector_path(), MATS)
const SOURCE_GEOM = JSON.parsefile(normpath(joinpath(@__DIR__, "..", "..", "data", "source_geometry_lz_v1.json")))
const VOL  = active_volume(DET)
const MAT  = VOL.material
const GEOM_TOL = 0.1
# TPC center position (cathode at z=0, TPC spans [0, 145.6])
const TPC_CENTER = (0.0, 0.0, VOL.logical.position[3])


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


# =====================================================================
# Test 3: Detector geometry
# =====================================================================
@testset "Detector geometry" begin
    @test DET.name == "LZ"
    @test length(DET.volumes) == 21  # legacy geometry plus 5 added field-cage/grid-holder proxy volumes (FV separate)
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
    fvvol = fiducial_volume(DET)

    rng1 = MersenneTwister(20260508)
    deps_old = simulate_event(2.61, fvvol, CFG;
                              position=(0.0, 0.0, 61.0),
                              direction=(0.0, 0.0, 1.0),
                              rng=rng1)
    rng2 = MersenneTwister(20260508)
    deps_new = propagate_gamma_in_fv(2.61, fvgeom, CFG;
                                     position=(0.0, 0.0, 61.0),
                                     direction=(0.0, 0.0, 1.0),
                                     rng=rng2)

    @test length(deps_new) == length(deps_old)
    @test all(is_inside_fv(fvgeom, d.position) for d in deps_new)
    @test sum(d.energy for d in deps_new) ≈ sum(d.energy for d in deps_old) atol=1e-9
    @test [d.source for d in deps_new] == [d.source for d in deps_old]

    rng3 = MersenneTwister(12345)
    deps_old_2 = simulate_event(2.61, fvvol, CFG;
                                position=(10.0, 0.0, 80.0),
                                direction=(0.0, 0.0, -1.0),
                                rng=rng3)
    rng4 = MersenneTwister(12345)
    deps_new_2 = propagate_gamma_in_fv(2.61, fvgeom, CFG;
                                       position=(10.0, 0.0, 80.0),
                                       direction=(0.0, 0.0, -1.0),
                                       rng=rng4)

    @test length(deps_new_2) == length(deps_old_2)
    @test all(is_inside_fv(fvgeom, d.position) for d in deps_new_2)
    @test sum(d.energy for d in deps_new_2) ≈ sum(d.energy for d in deps_old_2) atol=1e-9
    @test [d.source for d in deps_new_2] == [d.source for d in deps_old_2]
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

    # next_volume: ray from outside MARS toward ICV barrel
    by_name = Dict(v.name => v for v in DET.volumes)
    icv = by_name["ICV_barrel"]
    # From r=100, aimed at center, should hit ICV barrel
    vol, t = next_volume([100.0, 0.0, 53.585], [-1.0, 0.0, 0.0], DET)
    @test vol !== nothing
    @test t < 100.0  # should hit something before traveling 100 cm
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


# =====================================================================
# Test 15: Veto pre-filter (propagate_to_fiducial)
# =====================================================================
@testset "Veto pre-filter" begin
    rng = MersenneTwister(20)
    E0 = 2.615
    N = 1000

    # Gammas from TPC wall (r=72.8, aimed inward at z=72.8 = TPC center height)
    n_vetoed = 0; n_accepted = 0; n_lost = 0
    n_int_vetoed = Int[]
    for _ in 1:N
        φ = 2π * rand(rng)
        pos = (72.0 * cos(φ), 72.0 * sin(φ), 72.8)
        dir_vec = (-cos(φ), -sin(φ), 0.0)

        result = propagate_to_fiducial(E0, pos, dir_vec, DET, CFG, rng)
        if result.status === :vetoed
            n_vetoed += 1
            push!(n_int_vetoed, result.n_interactions)
        elseif result.status === :accepted
            n_accepted += 1
        else
            n_lost += 1
        end
    end

    frac_vetoed  = n_vetoed / N
    frac_accepted = n_accepted / N

    # Radial shielding: 33.8 cm / 8.9 cm mfp ≈ 3.8 λ
    # Unscattered survival ≈ exp(-3.8) ≈ 2.2%
    # Accepted is slightly higher (~3%) due to forward Compton below veto threshold
    @test 0.90 < frac_vetoed < 0.99     # ~97% vetoed
    @test 0.01 < frac_accepted < 0.07   # ~2-5% reach FV
    @test n_lost == 0                    # radial gammas don't escape

    # Nearly all vetoed events are killed at the first interaction
    # (99.8% at interaction 1, 0.2% at interaction 2, from 10k-event study)
    frac_first = count(==(1), n_int_vetoed) / length(n_int_vetoed)
    @test frac_first > 0.95              # >95% vetoed at first interaction
    @test maximum(n_int_vetoed) <= 3     # never needs more than ~2-3 steps
end


# =====================================================================
# Test 17: Source flux generation
# =====================================================================
@testset "Source flux" begin
    decays = load_decays()
    by_name = Dict(v.name => v for v in DET.volumes)
    source = by_name["ICV_barrel"]
    rng = MersenneTwister(42)
    N = 10000

    # Bi-214: single gamma, no companion veto
    ft_bi = generate_source_flux(N, source, decays["Bi214"], DET, CFG, rng)
    @test ft_bi.N_generated == N
    @test ft_bi.N_vetoed == 0                      # no companions → no veto
    @test ft_bi.N_surviving > 0
    # ~38% survive (forward + in energy window); ~51% backward
    @test 0.25 < survival_fraction(ft_bi) < 0.50
    @test ft_bi.N_backward > N * 0.3               # ~50% go backward
    # Peak bin: most surviving Bi-214 are unscattered (2.448 in bin 8)
    f_pk = peak_bin_fraction(ft_bi, 2.448)
    f_off = off_peak_fraction(ft_bi, 2.448)
    @test f_pk > f_off                             # Bi-214: peak dominates (thin source)
    # Flux table dimensions
    @test size(ft_bi.pdf) == (25, 10)
    @test sum(ft_bi.pdf) ≈ ft_bi.N_surviving / ft_bi.N_generated  rtol=1e-10

    # Tl-208: companion gamma causes veto
    rng2 = MersenneTwister(42)
    ft_tl = generate_source_flux(N, source, decays["Tl208"], DET, CFG, rng2)
    @test ft_tl.N_generated == N
    @test ft_tl.N_vetoed > N * 0.15                 # ~25% vetoed (both gammas inward + visible)
    @test ft_tl.N_surviving > 0
    @test ft_tl.N_backward > N * 0.15              # ~26% backward (no gamma toward LXe)
    @test sum(ft_tl.pdf) ≈ ft_tl.N_surviving / ft_tl.N_generated  rtol=1e-10
    # Tl-208: peak (2.615) should be large (most exit unscattered from thin Ti)
    # but off-peak (scattered) is the dangerous background
    f_pk_tl = peak_bin_fraction(ft_tl, 2.615)
    f_off_tl = off_peak_fraction(ft_tl, 2.615)
    @test f_pk_tl > 0.0
    @test f_off_tl > 0.0
end
