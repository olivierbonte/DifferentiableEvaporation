@testset "Mass conservation of tendencies" begin
    forcings = constant_forcings(P)
    for thresholds in treatments,
        u0 in map(u0_from_fractions, u0_fractions)

        du = zeros(3)
        compute_tendencies!(du, u0, p_model, t_unix, forcings, thresholds)
        d = compute_diagnostics(u0, p_model, t_unix, forcings, thresholds)
        @test all(isfinite, du)
        @test d.λE_tot ≈ d.λE_t + d.λE_i + d.λE_s
        # Water balance of the column, d/dt(ρ_w d_2 w_2 + w_r) = inputs - outputs, with both
        # sides in kg m⁻² s⁻¹ (mm s⁻¹).
        lhs = ρ_w * p_model.d_2 * du[2] + du[3]
        rhs = P - d.Q_s - d.E_s - d.E_t - d.E_i - ρ_w * p_model.d_2 * d.K_2
        @test lhs ≈ rhs atol = 1e-12
    end
end

@testset "Energy conservation of diagnostics" begin
    # Daytime (10:00, R_n > 0) and nighttime (02:00, R_n < 0) conditions
    t_night = datetime2unix(DateTime(2010, 7, 1, 2))
    for (R_n, t) in ((Rn, t_unix), (-50.0, t_night)),
        thresholds in treatments,
        u0 in map(u0_from_fractions, u0_fractions)

        forcings = merge(constant_forcings(P), (; R_n=Returns(R_n)))
        d = compute_diagnostics(u0, p_model, t, forcings, thresholds)
        λE_c = d.λE_t + d.λE_i
        # 1) Total energy balance
        @test all(isfinite, (d.R_nc, d.R_ns, d.G, d.λE_tot, d.H_tot))
        @test d.R_nc + d.R_ns ≈ R_n
        @test abs(d.G) < abs(d.R_ns) # G is a fraction of R_ns
        @test d.λE_tot + d.H_tot ≈ R_n - d.G
        # 2) Component balances sum to the total
        # Soil
        @test d.λE_s + d.H_s ≈ d.R_ns - d.G
        # Canopy
        @test λE_c + d.H_c ≈ d.R_nc
        # Total
        @test d.λE_tot ≈ λE_c + d.λE_s
        @test d.H_tot ≈ d.H_c + d.H_s
    end
end

@testset "ProcessBasedModel initialize! and solve!" begin
    # Test exported solvers
    algs = (
        Tsit5(),
        Heun(),
        ImplicitEuler(; autodiff=AutoForwardDiff()),
        ImplicitEuler(; autodiff=AutoEnzyme(; function_annotation=Enzyme.Duplicated)),
    )
    for thresholds in treatments, alg in algs
        model = toy_model(; thresholds)
        @test model.solver_kwargs == (; abstol=1e-6, reltol=1e-6)
        EvaporationModel.solve!(model; alg=alg)
        @test OrdinaryDiffEq.SciMLBase.successful_retcode(model.sol)
        @test model.sol.t == saveat_toy
        @test all(u -> all(isfinite, u), model.sol.u)
    end
    # Keyword arguments of solve! take precedence over solver_kwargs
    model_default, model_tight = toy_model(), toy_model()
    EvaporationModel.solve!(model_default; alg=Tsit5())
    EvaporationModel.solve!(model_tight; alg=Tsit5(), abstol=1e-10, reltol=1e-10)
    # tighter tolerances require more accepted steps
    @test model_tight.sol.stats.naccept > model_default.sol.stats.naccept
end

@testset "Accumulated diagnostics" begin
    # Any diagnostic can be accumulated: the window means close the water balances of the
    # canopy and of the column, also with the (backward Euler consistent) ImplicitEuler
    water = (:P_c, :E_i, :D_c, :I_s, :E_s, :E_t, :K_2)
    tight = (; abstol=1e-10, reltol=1e-10)
    solvers = ((Tsit5(), tight), (ImplicitEuler(; autodiff=AutoForwardDiff()), (;)))
    for thresholds in treatments, (alg, tol) in solvers
        model = toy_model(; thresholds, accumulated=water)
        EvaporationModel.solve!(model; alg, tol...)
        m = model.diagnostics
        Δt = diff(saveat_toy)
        dw_r_dt = diff(model.sol[3, :]) ./ Δt
        dw_2_dt = diff(model.sol[2, :]) ./ Δt
        @test dw_r_dt ≈ m.P_c .- m.E_i .- m.D_c rtol = 1e-6 atol = 1e-12
        @test dw_2_dt ≈ (m.I_s .- m.E_s .- m.E_t) ./ (ρ_w * p_model.d_2) .- m.K_2 rtol = 1e-6 atol =
            1e-14
    end

    # Accumulators do not change the time steps or the physical states (implicit solvers: up
    # to the Newton tolerance, as the Newton solve includes the accumulators)
    for (alg, rtol) in ((Tsit5(), 1e-12), (ImplicitEuler(; autodiff=AutoForwardDiff()), 1e-5))
        with, without = toy_model(), toy_model(; accumulated=())
        EvaporationModel.solve!(with; alg)
        EvaporationModel.solve!(without; alg)
        @test with.sol.stats.naccept == without.sol.stats.naccept
        @test with.sol[1:3, :] ≈ Array(without.sol) rtol = rtol
        @test without.diagnostics == NamedTuple()
    end

    # Invalid input
    @test_throws ArgumentError toy_model(; accumulated=(:not_a_diagnostic,))
    @test_throws ArgumentError toy_model(; u0=[u0_toy; 0.0])
    @test_throws ArgumentError EvaporationModel.solve!(toy_model(); abstol=fill(1e-6, 6))
end
