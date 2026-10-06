@testset "Threshold treatments" begin
    hard, smooth = treatments
    @test max(hard, 0.2, 0.3, moisture_scale(hard)) == 0.3
    @test min(hard, 0.2, 0.3, moisture_scale(hard)) == 0.2
    @test clamp(hard, 1.5, 0.0, 1.0, factor_scale(hard)) == 1.0
    # The smooth clamp never drops below the lower bound, and exceeds the upper by < s/2
    s = factor_scale(smooth)
    for x in (-1.0, 0.0, 0.5, 1.0, 2.0)
        @test 0 <= clamp(smooth, x, 0.0, 1.0, s) <= 1 + s / 2
    end
    # Smoothed f_wet has a finite slope at w_r = 0 (the hard 2/3 power does not)
    w_rmax = 0.6
    @test isfinite(
        derivative(w_r -> fraction_wet_vegetation(w_r, w_rmax, smooth), AutoForwardDiff(), 0.0)
    )
    # HardThresholds is the default treatment
    forcings = constant_forcings(P)
    u0 = [0.2, 0.32, 0.3]
    @test compute_diagnostics(u0, p_model, t_unix, forcings) ==
        compute_diagnostics(u0, p_model, t_unix, forcings, hard)
end

@testset "Surface layer bounded at saturation" begin
    @test surface_infiltration_factor(w_sat, w_sat) == 0
    @test surface_infiltration_factor(0.0, w_sat) ≈ 1
    @test surface_infiltration_factor(w_sat - 0.05, w_sat) > 0.99
    # Heavy rain (20 mm/h) on a saturated surface layer during the day: no further wetting
    P_heavy = 20 / 3600 # kg / (m2 * s)
    forcings = constant_forcings(P_heavy)
    for thresholds in treatments
        du = zeros(3)
        compute_tendencies!(du, [w_sat, 0.4, 0.7], p_model, t_unix, forcings, thresholds)
        @test du[1] <= 0
        # 12 h of heavy rain on a wet surface over a drier root zone (so most rain
        # infiltrates) drives w_1 to the bound, but not above it
        t_span = (t_unix, t_unix + 12 * 3600.0)
        model = ProcessBasedModel{Float64}(;
            forcings=forcings,
            parameters=p_model,
            t_span=t_span,
            u0=[0.35, 0.25, 0.5],
            saveat=collect(t_span[1]:600.0:t_span[2]),
            thresholds=thresholds,
        )
        EvaporationModel.initialize!(model)
        # Explicit solver: the ForwardDiff Jacobian is NaN for w_1 > w_fc, where r_ss = 0
        EvaporationModel.solve!(model; alg=Tsit5(), abstol=1e-8, reltol=1e-8)
        @test OrdinaryDiffEq.SciMLBase.successful_retcode(model.sol)
        @test maximum(u -> u[1], model.sol.u) <= w_sat + 1e-6 # w_1 bounded
        @test maximum(u -> u[2], model.sol.u) <= w_sat + 1e-6 # w_2 bounded
        @test maximum(u -> u[1], model.sol.u) > w_sat - 0.01 # w_1 approaches bound
    end
end

@testset "Check evaporation sum" begin
    R_nc, R_ns = net_radiation_partitioning(Rn, f_veg)
    A, A_c, A_s = available_energy_partitioning(R_nc, R_ns, G)
    λE_tot, λE_tot_p = total_evaporation(
        T_a, p_a, VPD_a, A, A_c, A_s, r_aa, r_ac, r_as, r_sc, r_ss, f_wet
    )
    VPD_m = vpd_veg_source_height(VPD_a, T_a, p_a, A, λE_tot, r_aa)
    ET_t, λE_t = transpiration(T_a, p_a, VPD_m, A_c, r_ac, r_sc, f_wet)
    ET_i, λE_i = interception_loss(T_a, p_a, VPD_m, A_c, r_ac, f_wet)
    ET_s, λE_s = soil_evaporation(T_a, p_a, VPD_m, A_s, r_as, r_ss)
    @test λE_tot ≈ λE_t + λE_i + λE_s
end
