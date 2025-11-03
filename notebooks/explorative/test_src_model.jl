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
using SciMLSensitivity
using Plots
using YAXArrays
using Zygote: Zygote
include("example_input.jl")

@benchmark compute_diagnostics($u0_test, $param_test, $t_span_test[1], $forcings)
@code_warntype compute_diagnostics(u0_test, param_test, t_span_test[1], forcings)
# Zero allocs, type stable

du_test = similar(u0_test)
compute_tendencies!(du_test, u0_test, param_test, t_span_test[1], forcings)
@benchmark compute_tendencies!($du_test, $u0_test, $param_test, $t_span_test[1], $forcings)
# Zero allocs, type stable

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
# Two allocs (<-> zero allocs above...)
# Compare to SciML version of f inside problem
@code_warntype test_model.prob.f(du_test, u0_test, param_test, t_span_test[1])
@benchmark test_model.prob.f($du_test, $u0_test, $param_test, $t_span_test[1])
# Three allocs

# Key Q: does this matter in solving the ODE?
const forcings_const = forcings
test_func_no_alloc(du, u, p, t) = compute_tendencies!(du, u, p, t, forcings_const)
@benchmark test_func_no_alloc($du_test, $u0_test, $param_test, $t_span_test[1])
# zero allocs for this closure with forcings as const (impossible inside of struct)
test_prob_no_alloc = ODEProblem(test_func_no_alloc, u0_test, t_span_test, param_test)
@benchmark solve($test_prob_no_alloc) # 925 allocs
@benchmark solve($test_model.prob) # 926 allocs
# Conclusion: the 2 vs 0 allocs in f does not seem to matter for the overall solve allocs

# Test solve
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
cols_prognostics = filter(x -> x != "timestamp", names(df_prognostics))
rename!(df_prognostics, map(=>, cols_prognostics, ["w_1", "w_2", "w_r"]))
df_all = hcat(df_diagnostics, df_prognostics)
df_all_no_time = df_all[:, names(df_all) .!= "timestamp"]
time_saved = unix2datetime.(test_model_BE_Bra.diagnostics.t)
axlist = (YAXArrays.time(time_saved), Variables(names(df_all_no_time)))
test_yax_array = YAXArray(axlist, Array(df_all_no_time))

# Implemented in source code
plot(test_model_BE_Bra.output)

## Experiment with Differentiating the code
# NOTE: TIMESTAMP ARE NOT YET CORRECTLY ALLIGNED!
y_obs = ds_ec_sel.Qle[x=1, y=1]
function loss_function(param, y_obs, forcings, t_span, u0, saveat; kwargs...)
    model = ProcessBasedModel{FT}(;
        forcings=forcings, parameters=param, t_span=t_span, u0=u0, saveat=saveat
    )
    EvaporationModel.initialize!(model)
    EvaporationModel.solve!(model; AD=true, kwargs...)
    # y_pred = model.output[Variables=At("λE_tot")] Use of YAXArray and dataframes problematic with AD
    y_pred = [
        model.diagnostics.saveval[i].λE_tot for i in 1:length(model.diagnostics.saveval)
    ]
    return sum((y_pred .- y_obs) .^ 2)
end
loss_function(param_test, y_obs, forcings_real, t_span_real, u0_test, t_unix)
function loss_function_wrapper(param; kwargs...)
    return loss_function(
        param, y_obs, forcings_real, t_span_real, u0_test, t_unix; kwargs...
    )
end
# Direct AD
@benchmark gradient(loss_function_wrapper, AutoForwardDiff(), param_test) #35 ms
@benchmark gradient(loss_function_wrapper, AutoFiniteDiff(), param_test) # 350 ms
# Continuous AD
loss_function_gauss(param) = loss_function_wrapper(param; sensealg=GaussAdjoint())
gradient(loss_function_gauss, AutoReverseDiff(), param_test) # wrong gradient
#gradient(loss_function_gauss, AutoZygote(), param_test) # Fails

## Easier AD experiments: just sum the states!
function sum_of_solutions(model, p; kwargs...)
    _prob = remake(model.prob; p=p)
    return sum(solve(_prob, Heun(); kwargs...))
end
test_sum = sum_of_solutions(test_model_BE_Bra, param_test; tstops=t_unix)

sum_wrapper(p) = sum_of_solutions(test_model_BE_Bra, p; tstops=t_unix)
@benchmark gradient(sum_wrapper, AutoForwardDiff(), $param_test)
@benchmark gradient(sum_wrapper, AutoFiniteDiff(), $param_test)
# Note: these are fast, but this is because forwarddiff is choisen as a default sensealg
# I think, see
# https://github.com/SciML/SciMLSensitivity.jl/blob/5d955d387beb90c10bd7a3d85f15ed68baed50f5/src/concrete_solve.jl#L322
@benchmark gradient(sum_wrapper, AutoReverseDiff(), $param_test)
@benchmark gradient(sum_wrapper, AutoZygote(), $param_test)

function sum_wrapper_adjoint(p)
    return sum_of_solutions(test_model_BE_Bra, p; sensealg=GaussAdjoint(), tstops=t_unix)
end
@benchmark gradient(sum_wrapper_adjoint, AutoReverseDiff(), $param_test)
@benchmark gradient(sum_wrapper_adjoint, AutoZygote(), $param_test)

function sum_wrapper_reverse(p)
    return sum_of_solutions(
        test_model_BE_Bra, p; sensealg=ReverseDiffAdjoint(), tstops=t_unix
    )
end
gradient(sum_wrapper_reverse, AutoReverseDiff(), param_test) # All Nan
# function sum_wrapper_zygote(p)
#     return sum_of_solutions(test_model_BE_Bra, p; sensealg=ZygoteAdjoint(), tstops=t_unix)
# end
# gradient(sum_wrapper_zygote, AutoZygote(), param_test)
