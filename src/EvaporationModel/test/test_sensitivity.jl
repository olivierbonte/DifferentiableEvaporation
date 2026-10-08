# Gradients of a toy loss through the ODE solve. The solver kwargs (tolerances) are passed to
# `solve`, so that the adjoint solves use them as well.

function toy_loss(p, prob, saveat, y_obs, sensealg, solver_kwargs)
    sol = solve(remake(prob; p=p), Tsit5(); saveat, sensealg, solver_kwargs...)
    return sum(abs2, Array(sol) .- y_obs) / length(y_obs)
end

# For differentiating directly through the solver with Enzyme: the problem is built inside the
# loss, with `FullSpecialize` (no FunctionWrappers around `f`), instead of a `remake` of a
# `Constant` problem, and `SensitivityADPassThrough` skips the SciMLSensitivity adjoint rules
function toy_loss_direct(p, f, u0, tspan, saveat, y_obs, solver_kwargs)
    prob = ODEProblem{true,SciMLBase.FullSpecialize}(f, copy(u0), tspan, p)
    sensealg = DiffEqBase.SensitivityADPassThrough()
    sol = solve(prob, Tsit5(); saveat, sensealg, solver_kwargs...)
    return sum(abs2, Array(sol) .- y_obs) / length(y_obs)
end

# Error norm of the adaptive step size control, ignored by Enzyme, so that the step sizes are
# treated as constants (as ForwardDiff does)
inactive_norm(u, t) = Enzyme.ignore_derivatives(DiffEqBase.ODE_DEFAULT_NORM(u, t))

@testset "Sensitivities of a toy loss" verbose = true begin
    # Reference = slightly perturbed parameters
    model = toy_model()
    p_obs = copy(p_model)
    p_obs.r_smin *= 1.2
    y_obs = Array(
        solve(remake(model.prob; p=p_obs), Tsit5(); saveat=saveat_toy, model.solver_kwargs...)
    )
    # Every gradient call gets its own copy of the problem, so that a gradient cannot change
    # the problem used by the next one (see the `prob.u0` test below)
    fresh_prob() = remake(model.prob; u0=copy(model.u0))
    context(sensealg, solver_kwargs=model.solver_kwargs) = (
        Constant(fresh_prob()), Constant(saveat_toy), Constant(y_obs), Constant(sensealg),
        Constant(solver_kwargs),
    )

    # Forward mode: the reference gradient, once sanity checked against finite differences.
    # Thight tolerances to minimalise the effect of the solver tolerances on the gradient.
    tight = (; abstol=1e-12, reltol=1e-12)
    g_ref = gradient(
        toy_loss, AutoForwardDiff(), p_model, context(ForwardDiffSensitivity(), tight)...
    )
    @testset "Forward mode (ForwardDiff)" begin
        @test all(isfinite, g_ref)
        g_fd = gradient(
            toy_loss, AutoFiniteDiff(), p_model, context(ForwardDiffSensitivity(), tight)...
        )
        @test isapprox(g_ref, g_fd; rtol=1e-3)
    end

    @testset "Reverse mode (Enzyme + GaussAdjoint with EnzymeVJP)" begin
        enzyme_reverse = AutoEnzyme(; mode=Enzyme.set_runtime_activity(Enzyme.Reverse))
        sensealg = GaussAdjoint(; autojacvec=EnzymeVJP())
        g = gradient(toy_loss, enzyme_reverse, p_model, context(sensealg, tight)...)
        @test isapprox(g, g_ref; rtol=1e-3)
        prob = fresh_prob()
        gradient(toy_loss, enzyme_reverse, p_model, Constant(prob), Constant(saveat_toy),
            Constant(y_obs), Constant(sensealg), Constant(model.solver_kwargs))
        # Enzyme reverse mode used to write into `prob.u0` although the problem is a
        # `Constant`, so a later forward solve started from changed initial conditions
        @test prob.u0 == model.u0
    end

    # Enzyme reverse mode through the solver internals (discretize-then-differentiate),
    # following the Enzyme.jl SciML integration tests:
    # https://github.com/EnzymeAD/Enzyme.jl/blob/d11e92e353b2469d5fe513b3c356e52baf2b47b5/test/integration/SciML/runtests.jl#L45
    @testset "Reverse mode (Enzyme, direct through the solver)" begin
        enzyme_reverse = AutoEnzyme(; mode=Enzyme.set_runtime_activity(Enzyme.Reverse))
        # Without `inactive_norm`, Enzyme also differentiates the adaptive step size control,
        # which gives a wrong gradient (relative error ~1)
        solver_kwargs = (; tight..., internalnorm=inactive_norm)
        g = gradient(
            toy_loss_direct, enzyme_reverse, p_model, Constant(model.f), Constant(model.u0),
            Constant(model.t_span), Constant(saveat_toy), Constant(y_obs),
            Constant(solver_kwargs),
        )
        @test isapprox(g, g_ref; rtol=1e-3)
    end
end
