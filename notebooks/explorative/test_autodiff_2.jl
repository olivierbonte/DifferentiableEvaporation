# %% Imports
using DrWatson
@quickactivate "DifferentiableEvaporation"

using Revise
using Bigleaf, EvaporationModel
using DiffEqBase, DifferentiationInterface, OrdinaryDiffEq, StaticArrays
using Enzyme: Enzyme

# Test dx/dt for gradients
include("example_input.jl")
test_model = ProcessBasedModel{FT}(;
    forcings=forcings,
    parameters=param_test,
    t_span=t_span_test,
    u0=u0_test,
    saveat=saveat_test,
)
EvaporationModel.initialize!(test_model)

# Try to construct scimlproblem from scratch
# easiest example = dy[1,end]/du0 to be calculated
f_tendencies = test_model.f
const saveat_const = SA[collect(saveat_test)...]
function f_simple(u0)
    prob = ODEProblem{true,SciMLBase.FullSpecialize}(
        f_tendencies, u0, t_span_test, param_test
    )
    sol = DiffEqBase.solve(
        prob,
        Tsit5();
        abstol=1e-6,
        reltol=1e-6,
        saveat=saveat_const,
        sensealg=DiffEqBase.SensitivityADPassThrough(),
    )
    return sol[1, end]
end
f_simple(u0_test) # test run
# Using ForwardDiff
val_fd, grad_fd = value_and_gradient(f_simple, AutoForwardDiff(), u0_test)
# Works in nansafe mode of ForwardDiff

# Using Enzyme
val_test, jac_test = value_and_gradient(
    f_simple,
    AutoEnzyme(; mode=Enzyme.set_runtime_activity(Enzyme.set_strong_zero(Enzyme.Forward))),
    u0_test,
)
