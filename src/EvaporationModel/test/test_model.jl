@testset "Mass conservation of tendencies" begin
    forcings = constant_forcings(P)
    for thresholds in treatments,
        u0 in ([0.2, 0.32, 0.3], [0.05, 0.1, 0.0], [0.4, 0.44, 0.7], [w_sat, w_sat, 0.7])

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
