using DrWatson
@quickactivate "DifferentiableEvaporation"
using Bigleaf, EvaporationModel
using ComponentArrays
using OrdinaryDiffEq
using DifferentiationInterface
using Enzyme: Enzyme
using ForwardDiff: ForwardDiff
using Zygote: Zygote
using SciMLSensitivity
using BenchmarkTools
using Revise
# Note that DifferentiatioIterface exports ADTypes

# %% Test potential_et from Bigleaf with ForwardDiff
Tair = 30.0 # [°C]
pressure = 101.325 # [kPa]
Rn = 400.0 # [W m-2]
VPD = 0.5 # [kPa]
Ga = 1 / 30.0 # [m s-1]
approach = PenmanMonteith()

test_et = potential_ET(approach, Tair, pressure, Rn, VPD, Ga)

# wrapper function
f(x) = potential_ET(PenmanMonteith(), x[1], x[2], x[3], x[4], x[5])[2]

x_trial = [Tair, pressure, Rn, VPD, Ga]

f(x_trial)

# try to differentiate
gradient(f, AutoForwardDiff(), x_trial)
gradient(f, AutoEnzyme(), x_trial)

@benchmark gradient($f, $AutoForwardDiff(), $x_trial)
@benchmark gradient($f, $AutoEnzyme(), $x_trial)

# pre-allocate gradient

grad = similar(x_trial)
@benchmark gradient!($f, $grad, $AutoForwardDiff(), $x_trial)
@benchmark gradient!($f, $grad, $AutoEnzyme(), $x_trial)

# prepare gradient (for Enzyme)
prep = prepare_gradient(f, AutoEnzyme(), zero(x_trial))
grad = similar(x_trial)
@benchmark gradient!($f, $grad, $prep, $AutoEnzyme(), $x_trial)

# %% Test functions from EvaporationModel
x_trial_bis = [Tair + 273.15, pressure * 1000, VPD * 1000, Rn, 1 / Ga, 1 / Ga * 4]
penman_monteith(Tair + 273.15, pressure * 1000, VPD * 1000, Rn, 1 / Ga, 1 / Ga * 4)

f_bis = x -> penman_monteith(x[1], x[2], x[3], x[4], x[5], x[6])[2]
gradient(f_bis, AutoForwardDiff(), x_trial_bis)

function f_ode(u, p, t)
    return -1 / (1000 * 1.5) * f_bis(p)
end
u0 = 0.5
tspan = (0.0, 20.0)
prob = ODEProblem(f_ode, u0, tspan, x_trial_bis)
@benchmark solve(
    prob, Rosenbrock23(; autodiff=AutoEnzyme(; function_annotation=Enzyme.Duplicated))
)
@benchmark solve(prob, Rosenbrock23(; autodiff=AutoForwardDiff()))

## Try to fix AD for full set of conservation equations
# Fixed/synthetic forcings as an easy start
P_test = (t) -> 5e-6 #kg/(m² s)
T_a_test = (t) -> 275.0 + 5.0 * sin(2.0 * π * t / 86400 - π / 2);
u_a_test = (t) -> 3.0 # m /s
p_a_test = (t) -> 101325.0 # Pa
VPD_a_test = (t) -> 200.0 # Pa
SW_in_test = (t) -> max(1361.0 * sin(2π * t / 86400 - π / 2), 0);
R_n_test = (t) -> SW_in_test(t) * 0.3 # W/m²
LAI_test = (t) -> 3.0

param_test = ComponentArray(;
    h=21.0, # [m] height of the canopy
    LAI=4.0, # [-] leaf area index
    z_0ms=0.01, # [m] roughness length for soil
    w_sat=0.5, # [-]  saturation water content
    a=0.15,
    p_soil=6.0,
    b=6.1,
    C_1sat=1.9,
    C_2ref=0.83,
    d_1=0.01,
    z_obs=39.0, # [m]
)
struct ModelParamsTest{T<:AbstractFloat}
    h::T # [m] height of the canopy
    LAI::T # [-] leaf area index
    z_0ms::T # [m] roughness length for soil
    w_sat::T # [-] saturation water content
    a::T
    p_soil::T
    b::T
    C_1sat::T
    C_2ref::T
    d_1::T
    z_obs::T # [m]
end
model_params_test = ModelParamsTest(
    param_test.h,
    param_test.LAI,
    param_test.z_0ms,
    param_test.w_sat,
    param_test.a,
    param_test.p_soil,
    param_test.b,
    param_test.C_1sat,
    param_test.C_2ref,
    param_test.d_1,
    param_test.z_obs,
)
# Test to reduce allocations
mutable struct TestDynParamModel
    f_veg
end
test_dyn_param_model = TestDynParamModel(0.4)
function fractional_vegetation_cover!(dyn_ar::TestDynParamModel, LAI)
    # This is a dummy function to test if AD works with dynamic parameters
    dyn_ar.f_veg = fractional_vegetation_cover(LAI)
end
fractional_vegetation_cover!(test_dyn_param_model, LAI_test(0.0))
# In-place form f!(du, u, p, t) should be fastest
u0_test = [0.2, 0.2, 0.001]
t_end = 6 * 86400.0 # 1 day in seconds
t_span_test = (0.0, t_end) # 1 day in seconds
# GOAL: get code below optimised!
function calculate_fluxes_test!(du, u, p, t)
    w_1, w_2, w_r = u
    # This adds 11 allocations, not wanted...
    @unpack h, LAI, z_0ms, w_sat, a, p_soil, b, C_1sat, C_2ref, d_1, z_obs = p
    d_c, z_0mc = Bigleaf.roughness_parameters(RoughnessCanopyHeightLAI(), h, LAI; hs=z_0ms)
    f_veg = fractional_vegetation_cover(LAI)
    f_wet = fraction_wet_vegetation(w_r, LAI)
    w_1eq = w_geq(w_2, w_sat, a, p_soil) #no allocs
    C_1 = c_1(w_1, w_sat, b, C_1sat) # no allocs
    C_2 = c_2(w_2, w_sat, C_2ref) # no allocs

    G = ground_heat_flux(Allen07(), R_n_test(t), T_a_test(t))
    A, A_c, A_s = available_energy_partioning(R_n_test(t), G, f_veg)

    ustar = ustar_from_u(u_a_test(t), z_obs, d_c, z_0mc)
    # To avoid issues with types in autodiff...
    T = eltype(u)
    du[1] = C_1 / (ρ_w * d_1) * P_test(t)
    du[2] = zero(T)
    du[3] = zero(T)
    return nothing
end
# Alternative: make a mutating function
# Some benchmarking on this function. Goal to get this close to zero allocations!
du0_test = similar(u0_test)
@benchmark calculate_fluxes_test!($du0_test, $u0_test, $param_test, $t_span_test[1])
@code_warntype calculate_fluxes_test!(du0_test, u0_test, param_test, t_span_test[1])
# Test if AD works with this function
prob_test = ODEProblem(calculate_fluxes_test!, u0_test, t_span_test, param_test)
# Test if AD works within solver
sol = solve(prob_test, Rosenbrock23(; autodiff=AutoForwardDiff()))
sol = solve(
    prob_test, Rosenbrock23(; autodiff=AutoEnzyme(; function_annotation=Enzyme.Duplicated))
)

# Test if AD works for parameter optimisation purpose (i.e. take derivaties)
t_obs_test = collect(t_span_test[1]:1:t_span_test[2])
observed_data = rand(3, length(t_obs_test)) # Simulated observed data for testing
function loss_function_test(p; sensealg=AutoForwardDiff())
    simul = Array(
        solve(
            prob_test,
            Rosenbrock23(; autodiff=AutoForwardDiff());#AutoForwardDiff());
            saveat=t_obs_test,
            sensealg=sensealg,
            p=p,
        ),
    )
    return mean(abs2, simul .- observed_data)
end
@benchmark loss_function_test(param_test)
loss_function_test(param_test)
gradient(loss_function_test, AutoForwardDiff(), param_test)
@benchmark gradient(loss_function_test, AutoForwardDiff(), param_test)
@benchmark gradient(loss_function_test, AutoFiniteDiff(), param_test)

# Test adjoint equations with Zygote for reverse mode AD
loss_function_adjoint(p) = loss_function_test(p; sensealg=BacksolveAdjoint(;autojacvec=EnzymeVJP()))
@benchmark gradient(loss_function_adjoint, AutoZygote(), param_test)
# Much slower...
# BenchmarkTools.Trial: 1 sample with 1 evaluation per sample.
#   Single result which took 128.607 s (5.03% GC) to evaluate,
#   with a memory estimate of 31.54 GiB, over 1116636507 allocations.

#gradient(loss_function_test, AutoEnzyme(), param_test)
ForwardDiff.gradient(loss_function_test, param_test)
