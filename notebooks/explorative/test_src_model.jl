# %% Imports
using DrWatson
@quickactivate "DifferentiableEvaporation"

using Revise
using BenchmarkTools
using Bigleaf
using DataFrames
using DifferentiationInterface
using EvaporationModel
using OrdinaryDiffEq
using Plots
using YAXArrays
include("example_input.jl")

@benchmark compute_diagnostics($u0_test, $param_test, $t_span_test[1], $forcings)
@code_warntype compute_diagnostics(u0_test, param_test, t_span_test[1], forcings)

du_test = similar(u0_test)
compute_tendencies!(du_test, u0_test, param_test, t_span_test[1], forcings)
@benchmark compute_tendencies!($du_test, $u0_test, $param_test, $t_span_test[1], $forcings)
# function gradient_tendencies_wrapper(p)
#     return compute_tendencies!(du_test, u0_test, p, t_span_test[1], forcings)
# end
# @enter gradient(gradient_tendencies_wrapper, AutoForwardDiff(), param_test)
using FastClosures
const test_func = @closure (du, u, p, t) -> compute_tendencies!(du, u, p, t, forcings)
@benchmark test_func($du_test, $u0_test, $param_test, $t_span_test[1]) #zero allocs
# Not using FastClosures
const test_func_2 = (du, u, p, t) -> compute_tendencies!(du, u, p, t, forcings)
@benchmark test_func_2($du_test, $u0_test, $param_test, $t_span_test[1])
# Trying to make a constant forcings
const forcings_const = forcings
test_func_3(du, u, p, t) = compute_tendencies!(du, u, p, t, forcings_const)
@benchmark test_func_3($du_test, $u0_test, $param_test, $t_span_test[1])
# Using let block
test_func_4 = let forcings = forcings
    @closure (du, u, p, t) -> compute_tendencies!(du, u, p, t, forcings)
end
@benchmark test_func_4($du_test, $u0_test, $param_test, $t_span_test[1])

test_model = ProcessBasedModel{FT}(;
    forcings=forcings,
    parameters=param_test,
    t_span=t_span_test,
    u0=u0_test,
    saveat=saveat_test,
)
EvaporationModel.initialize!(test_model)
@code_warntype test_model.f(du_test, u0_test, param_test, t_span_test[1])
@benchmark test_model.f($du_test, $u0_test, $param_test, $t_span_test[1])
@code_warntype test_model.prob.f(du_test, u0_test, param_test, t_span_test[1])
@benchmark test_model.prob.f($du_test, $u0_test, $param_test, $t_span_test[1])
EvaporationModel.solve!(test_model; alg=ImplicitEuler(; autodiff=AutoForwardDiff()))
plot(test_model.sol)

# %% Real forcing data test
include("real_forcings.jl")
test_model_BE_Bra = ProcessBasedModel{FT}(;
    forcings=forcings_real,
    parameters=param_test,
    t_span=t_span_real,
    u0=u0_test,
    saveat=t_unix,
)
EvaporationModel.initialize!(test_model_BE_Bra)
@benchmark test_model_BE_Bra.f($du_test, $u0_test, $param_test, $t_unix[1])
# outside model source code definition
sol_outside = solve(test_model_BE_Bra.prob; tstops=t_unix)
# inside model source code definition (testing Rosenbrock method)
EvaporationModel.solve!(test_model_BE_Bra; alg=Rosenbrock23())
plot(test_model_BE_Bra.sol)

# Benchmarking a solution
@benchmark EvaporationModel.solve!(test_model_BE_Bra; alg=Heun())
@benchmark EvaporationModel.solve!(test_model_BE_Bra; alg=Tsit5())

# Test Ensemble simulation
function prob_func(prob, i, repeat)
    @. prob.u0 = max(0.2, rand()) * prob.u0
    return prob
end
ens_prob = EnsembleProblem(test_model_BE_Bra.prob; prob_func=prob_func)
@time ensemble_sim = solve(
    ens_prob;
    alg=Heun(),
    ensemblealg=EnsembleThreads(),
    tstops=t_unix,
    saveat=t_unix,
    trajectories=1000,
); #around 5s

# %% Experiment with saving diagnostics
# Manual
test_model_BE_Bra.diagnostics.saveval
# DataFrame
df_diagnostics = DataFrame(test_model_BE_Bra.diagnostics.saveval)
df_prognostics = DataFrame(test_model_BE_Bra.sol)
df_all = hcat(df_diagnostics, df_prognostics)
df_all_no_time = df_all[:, names(df_all) .!= "timestamp"]
time_saved = unix2datetime.(test_model_BE_Bra.diagnostics.t)
axlist = (YAXArrays.time(time_saved), Variables(names(df_all_no_time)))
