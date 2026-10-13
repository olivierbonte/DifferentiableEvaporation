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
        @test length(model.diagnostics.saveval) == length(saveat_toy)
        @test all(u -> all(isfinite, u), model.sol.u)
    end
    # Keyword arguments of solve! take precedence over solver_kwargs
    model_default, model_tight = toy_model(), toy_model()
    EvaporationModel.solve!(model_default; alg=Tsit5())
    EvaporationModel.solve!(model_tight; alg=Tsit5(), abstol=1e-10, reltol=1e-10)
    # tighter tolerances require more accepted steps
    @test model_tight.sol.stats.naccept > model_default.sol.stats.naccept
end
