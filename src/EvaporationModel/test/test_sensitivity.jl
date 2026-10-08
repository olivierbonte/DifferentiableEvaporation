# Gradients of a toy loss through the ODE solve. The solver kwargs (tolerances) are passed to
# `solve`, so that the adjoint solves use them as well.
const PassThrough = SciMLSensitivity.SensitivityADPassThrough

function toy_loss(p, prob, saveat, y_obs, sensealg, solver_kwargs)
    sol = solve(remake(prob; p=p), Tsit5(); saveat, sensealg, solver_kwargs...)
    return sum(abs2, Array(sol) .- y_obs) / length(y_obs)
end

@testset "Sensitivities of a toy loss" verbose = true begin
    # Reference = slightly perturbed parameters
    model = toy_model()
    p_obs = copy(p_model)
    p_obs.r_smin *= 1.2
    y_obs = Array(
        solve(remake(model.prob; p=p_obs), Tsit5(); saveat=saveat_toy, model.solver_kwargs...)
    )
    # Enzyme reverse mode writes into `prob.u0` although the problem is a `Constant`, so
    # every gradient call gets its own copy of the problem
    fresh_prob() = remake(model.prob; u0=copy(model.u0))
    context(sensealg) = (
        Constant(fresh_prob()), Constant(saveat_toy), Constant(y_obs), Constant(sensealg),
        Constant(model.solver_kwargs),
    )

    # Forward mode: the reference gradient once sanity checked against finite differences
    g_ref = gradient(toy_loss, AutoForwardDiff(), p_model, context(ForwardDiffSensitivity())...)
    @testset "Forward mode (ForwardDiff)" begin
        @test all(isfinite, g_ref)
        # Finite differences through an adaptive solve are noisy, so only a loose check
        g_fd = gradient(toy_loss, AutoFiniteDiff(), p_model, context(ForwardDiffSensitivity())...)
        @test isapprox(g_ref, g_fd; rtol=5e-2) #TODO: very loose tolerance...
    end

    @testset "Reverse mode (Enzyme + GaussAdjoint with EnzymeVJP)" begin
        enzyme_reverse = AutoEnzyme(; mode=Enzyme.set_runtime_activity(Enzyme.Reverse))
        sensealg = GaussAdjoint(; autojacvec=EnzymeVJP())
        g = gradient(toy_loss, enzyme_reverse, p_model, context(sensealg)...)
        @test isapprox(g, g_ref; rtol=1e-3)
        prob = fresh_prob()
        gradient(toy_loss, enzyme_reverse, p_model, Constant(prob), Constant(saveat_toy),
            Constant(y_obs), Constant(sensealg), Constant(model.solver_kwargs))
        @test_broken prob.u0 == model.u0
        # Enzyme reverse mode writes into `prob.u0` although the problem is a `Constant`
        # So if you solve forward again, initial conditions have changed...
    end

    # Enzyme through the solver internals. Forward mode runs but its gradients differ from
    # ForwardDiff by about 1 % (Enzyme v0.13.209). Reverse mode is not tested: it aborts
    # Julia ("LLVM ERROR: augmented function failed verification").
    @testset "SensitivityADPassThrough" begin
        enzyme_forward = AutoEnzyme(;
            mode=Enzyme.set_runtime_activity(Enzyme.set_strong_zero(Enzyme.Forward))
        )
        @test_broken isapprox(
            gradient(toy_loss, enzyme_forward, p_model, context(PassThrough())...), g_ref;
            rtol=1e-6,
        )
    end
end
