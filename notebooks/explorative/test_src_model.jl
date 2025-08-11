# %% Imports
using DrWatson
@quickactivate "DifferentiableEvaporation"

using Revise
using BenchmarkTools
using Bigleaf, EvaporationModel
using OrdinaryDiffEq
include("example_input.jl")

# @kwdef mutable struct Diagnostics{FT}
#     w_rmax::FT = zero(FT)
#     λE_tot::FT = zero(FT)
#     VPD_m::FT = zero(FT)
#     E_t::FT = zero(FT)
#     λE_t::FT = zero(FT)
#     E_i::FT = zero(FT)
#     λE_i::FT = zero(FT)
#     E_s::FT = zero(FT)
#     λE_s::FT = zero(FT)
#     D_c::FT = zero(FT)
#     P_s::FT = zero(FT)
#     Q_s::FT = zero(FT)
#     D_1::FT = zero(FT)
#     K_2::FT = zero(FT)
#     I_s::FT = zero(FT)
# end
# test_diagnostics = Diagnostics{FT}()

function create_diagnostics(FT)
    fields = (
        :w_r,
        :w_rmax,
        :C_1,
        :f_veg,
        :λE_tot,
        :VPD_m,
        :E_t,
        :λE_t,
        :E_i,
        :λE_i,
        :E_s,
        :λE_s,
        :D_c,
        :P_s,
        :Q_s,
        :D_1,
        :K_2,
        :I_s,
    )
    values = ntuple(i -> zero(FT), length(fields))
    return ComponentArray(NamedTuple{fields}(values))
end
test_diagnostics = create_diagnostics(FT)
@benchmark compute_diagnostics!(
    $test_diagnostics, $u0_test, $param_test, $t_span_test[1], $forcings
)
@code_warntype compute_diagnostics!(
    test_diagnostics, u0_test, param_test, t_span_test[1], forcings
)

du_test = similar(u0_test)
compute_tendencies!(
    du_test, u0_test, param_test, t_span_test[1], forcings, test_diagnostics
)
@benchmark compute_tendencies!(
    $du_test, $u0_test, $param_test, $t_span_test[1], $forcings, $test_diagnostics
)
test_func(du, u, p, t) = compute_tendencies!(du, u, p, t, forcings, test_diagnostics)
@benchmark test_func($du_test, $u0_test, $param_test, $t_span_test[1])

test_model = ProcessBasedModel{FT}(;
    forcings=forcings, parameters=param_test, t_span=t_span_test, u0=u0_test
)
EvaporationModel.initialize!(test_model)
@enter test_model.f(du_test, u0_test, param_test, t_span_test[1])
@code_warntype test_model.f(du_test, u0_test, param_test, t_span_test[1])
@benchmark test_model.f($du_test, $u0_test, $param_test, $t_span_test[1])
@benchmark test_model.prob.f($du_test, $u0_test, $param_test, $t_span_test[1])

# abstract type AbstractModel end
# @kwdef struct ProcessBasedModel{FT} <: AbstractModel
#     forcings::NamedTuple
#     parameters::AbstractArray
#     t_span
#     u0
#     diagnostics::ComponentArray = create_diagnostics(FT)
# end
