# Implicit solver on the full model, with the Jacobian from ForwardDiff or Enzyme
@testset "ImplicitEuler with AD Jacobian" begin
    for thresholds in treatments
        reference = toy_model(; thresholds)
        EvaporationModel.solve!(reference; alg=Tsit5(), abstol=1e-10, reltol=1e-10)
        u_ref = reduce(hcat, reference.sol.u)
        for autodiff in (AutoForwardDiff(), AutoEnzyme(; function_annotation=Enzyme.Duplicated))
            model = toy_model(; thresholds)
            EvaporationModel.solve!(model; alg=ImplicitEuler(; autodiff))
            @test OrdinaryDiffEq.SciMLBase.successful_retcode(model.sol)
            # First-order method at the default tolerances (1e-6)
            @test maximum(abs, reduce(hcat, model.sol.u) .- u_ref) < 1e-3
        end
    end
end
