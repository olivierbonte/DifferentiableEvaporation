# %% Imports
using DrWatson
@quickactivate "DifferentiableEvaporation"

using Revise
using Bigleaf, EvaporationModel
using ComponentArrays
using OrdinaryDiffEq
using DifferentiationInterface
using Enzyme: Enzyme
using ForwardDiff: ForwardDiff
using ReverseDiff: ReverseDiff
using Mooncake: Mooncake
using Zygote: Zygote
using SciMLSensitivity
using BenchmarkTools
using Preferences
set_preferences!(ForwardDiff, "nansafe_mode" => true)
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

# %% Gradient benchmarking
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
function P_test(t::T) where {T}
    return T(5e-6)
end
function T_a_test(t::T) where {T}
    return T(275.0) #+ T(5.0) * sin(T(2) * π * t / T(86400) - T(π / 2))
end
function u_a_test(t::T) where {T}
    return T(3.0)
end
function p_a_test(t::T) where {T}
    return T(101325.0)
end
function VPD_a_test(t::T) where {T}
    return T(200.0)
end
function SW_in_test(t::T) where {T}
    return T(1000) #max(T(1361.0) * sin(T(2) * π * t / T(86400) - π / T(2)), zero(t))
end
function R_n_test(t::T) where {T}
    return T(SW_in_test(t)) * T(0.3)
end
function LAI_test(t::T) where {T}
    return T(3.0)
end

param_test = ComponentArray(;
    h=21.0, # [m] height of the canopy
    z_0ms=0.01, # [m] roughness length for soil
    w_sat=0.5, # [-]  saturation water content
    a=0.15,
    p_soil=6.0,
    b=6.1,
    C_1sat=1.9,
    C_2ref=0.83,
    C_3=0.25,
    d_1=0.01,
    d_2=1.3,
    w_res=0.04,
    w_wp=0.08,
    w_fc=0.3,
    z_obs=39.0, # [m]
    kB⁻¹=log(10),
    g_d=0.0003,
    r_smin=395.0,
)
# In-place form f!(du, u, p, t) should be fastest
u0_test = [0.2, 0.2, 0.001]
t_end = 1 * 86400.0 # 1 day in seconds
t_span_test = (0.0, t_end) # 1 day in seconds
# GOAL: get code below optimised!
function calculate_fluxes_test!(du, u, p, t)
    w_1, w_2, w_r = u
    # This adds 11 allocations, not wanted...
    @unpack h,
    z_0ms,
    w_sat,
    a,
    p_soil,
    b,
    w_res,
    w_wp,
    w_fc,
    C_1sat,
    C_2ref,
    C_3,
    d_1,
    d_2,
    z_obs,
    kB⁻¹,
    g_d,
    r_smin = p
    d_c, z_0mc = Bigleaf.roughness_parameters(
        RoughnessCanopyHeightLAI(), h, LAI_test(t); hs=z_0ms
    )
    f_veg = fractional_vegetation_cover(LAI_test(t))
    w_rmax = max_canopy_capacity(LAI_test(t))
    f_wet = fraction_wet_vegetation(w_r, w_rmax)

    w_1eq = w_geq(w_2, w_sat, a, p_soil) #no allocs
    C_1 = c_1(w_1, w_sat, b, C_1sat) # no allocs
    C_2 = c_2(w_2, w_sat, C_2ref) # no allocs

    G = ground_heat_flux(Allen07(), R_n_test(t), LAI_test(t))
    A, A_c, A_s = available_energy_partioning(R_n_test(t), G, f_veg)

    # Resistances
    ustar = ustar_from_u(u_a_test(t), z_obs, d_c, z_0mc)
    r_aa = Bigleaf.compute_Ram(ResistanceWindZr(), ustar, u_a_test(t))
    r_ac = (Bigleaf.Gb_constant_kB1(ustar, kB⁻¹))^-1
    r_as = soil_aerodynamic_resistance(Choudhury1988soil(), ustar, h, d_c, z_0mc, z_0ms)
    r_sc = surface_resistance(
        JarvisStewart(),
        SW_in_test(t),
        VPD_a_test(t),
        T_a_test(t),
        w_2,
        w_fc,
        w_wp,
        LAI_test(t),
        g_d,
        r_smin,
    )
    β = soil_evaporation_efficiency(Pielke92(), w_1, w_fc)
    r_ss = beta_to_r_ss(β, r_as)

    # Turbulent fluxes calculations
    λE_tot, λE_tot_p = total_evaporation(
        T_a_test(t),
        p_a_test(t),
        VPD_a_test(t),
        A,
        A_c,
        A_s,
        r_aa,
        r_ac,
        r_as,
        r_sc,
        r_ss,
        f_wet,
    )
    VPD_m = vpd_veg_source_height(VPD_a_test(t), T_a_test(t), p_a_test(t), A, λE_tot, r_aa)
    E_t, λE_t = transpiration(T_a_test(t), p_a_test(t), VPD_m, A_c, r_ac, r_sc, f_wet)
    E_i, λE_i = interception(T_a_test(t), p_a_test(t), VPD_m, A_c, r_ac, f_wet)
    E_s, λE_s = soil_evaporation(T_a_test(t), p_a_test(t), VPD_m, A_s, r_as, r_ss)

    D_c = canopy_drainage(P_test(t), w_r, f_veg)
    P_s = precip_below_canopy(P_test(t), f_veg, D_c)
    Q_s = surface_runoff(StaticInfiltration(), P_test(t), w_2, w_fc)
    D_1 = diffusion_layer_1(w_1, w_1eq, C_2)
    K_2 = vertical_drainage_layer_2(w_2, w_fc, C_3, d_2)
    I_s = P_s - Q_s
    # Define smoothing parameter for canopy water content
    m_can = w_rmax / 100
    # ALLOC FREE UP UNITL HERE
    du[1] = C_1 / (ρ_w * d_1) * (I_s - E_s) - D_1
    du[2] = 1 / (ρ_w * d_2) * (I_s - E_s - E_t) - K_2
    du[3] =
        f_veg * P_test(t) * smoothing_kernel(UpperBound(), w_r, w_rmax, m_can) -
        E_i * smoothing_kernel(LowerBound(), w_r, zero(w_r), m_can) - D_c
    return nothing
end
# Alternative: make a mutating function
# Some benchmarking on this function. Goal to get this close to zero allocations!
du0_test = similar(u0_test)
@benchmark calculate_fluxes_test!($du0_test, $u0_test, $param_test, $t_span_test[1])
@code_warntype calculate_fluxes_test!(du0_test, u0_test, param_test, t_span_test[1])

# Test AD on the RHS of the ODE
function rhs_diff_test(u)
    du = similar(u)
    calculate_fluxes_test!(du, u, param_test, t_span_test[1])
    return du
end
rhs_diff_test(u0_test)
jac_adf = jacobian(rhs_diff_test, AutoForwardDiff(), u0_test)
jac_fd = jacobian(rhs_diff_test, AutoFiniteDiff(), u0_test)
jac_enz = jacobian(rhs_diff_test, AutoEnzyme(), u0_test)
isapprox.(jac_adf, jac_fd; rtol=1e-6, atol=1e-6) # Check if the jacobians are approximately equal

# Test if AD works with this function
prob_test = ODEProblem{true,SciMLBase.FullSpecialize}(
    calculate_fluxes_test!, u0_test, t_span_test, param_test
)
# Test if AD works within solver
sol = solve(prob_test, Rosenbrock23(; autodiff=AutoFiniteDiff())) #Works
sol = solve(prob_test, Rosenbrock23(; autodiff=AutoForwardDiff()))
sol = solve(
    prob_test, Rosenbrock23(; autodiff=AutoEnzyme(; function_annotation=Enzyme.Duplicated))
)
t_obs_test = collect(t_span_test[1]:1800:t_span_test[2]) # Save every 30 minutes

# %% Benchmarking ODE solving code
@benchmark solve($prob_test, Rosenbrock23(; autodiff=AutoForwardDiff()))
@benchmark solve(prob_test, Heun())
@benchmark solve(prob_test, Tsit5())
@benchmark solve(prob_test, Rodas5P(; autodiff=AutoForwardDiff()))
# See the effect of saving at every time step
@benchmark solve($prob_test, $Heun(), save_everystep=false) # faster of the 3
@benchmark solve($prob_test, $Heun(), tstops=$t_obs_test) #slowest of the 3
@benchmark solve($prob_test, $Heun(), saveat=$t_obs_test) # middle

# %% Test if AD works for parameter optimisation purpose (i.e. take derivaties)
observed_data = rand(3, length(t_obs_test)) # Simulated observed data for testing
function loss_function_test(p; t_obs_test=t_obs_test, sensealg=ForwardDiffSensitivity())
    simul = Array(solve(prob_test, Tsit5(); saveat=t_obs_test, sensealg=sensealg, p=p))
    return mean(abs2, simul .- observed_data)
end
loss_function_test(param_test)
@benchmark loss_function_test(param_test)
# function loss_function_interpol_test(p; sensealg=AutoForwardDiff())
#     sol = solve(
#         prob_test, Rosenbrock23(; autodiff=AutoForwardDiff()); sensealg=sensealg, p=p
#     )
#     simul = Array(sol(t_obs_test))
#     return mean(abs2, simul .- observed_data)
# end
# @benchmark loss_function_interpol_test(param_test)
gradient(loss_function_test, AutoForwardDiff(), param_test)
gradient(loss_function_test, AutoFiniteDiff(), param_test)
@benchmark gradient($loss_function_test, $AutoForwardDiff(), $param_test)
@benchmark gradient($loss_function_test, $AutoFiniteDiff(), $param_test)
# # %% Test direct reverse AD trhough the solver
# function loss_function_zygote_AD(p)
#     return loss_function_test(p; sensealg=ZygoteAdjoint())
# end
# gradient(loss_function_zygote_AD, AutoZygote(), param_test)

# %% Test adjoint equations with Zygote for reverse mode AD
function loss_function_adjoint(p)
    return loss_function_test(p; sensealg=BacksolveAdjoint(; autojacvec=EnzymeVJP()))
end
gradient(loss_function_adjoint, AutoZygote(), param_test)
@benchmark gradient(loss_function_adjoint, AutoZygote(), param_test)
# Much slower...
# BenchmarkTools.Trial: 3471 samples with 1 evaluation per sample.
#  Range (min … max):  1.023 ms … 12.045 ms  ┊ GC (min … max): 0.00% … 86.09%
#  Time  (median):     1.169 ms              ┊ GC (median):    0.00%
#  Time  (mean ± σ):   1.434 ms ±  1.125 ms  ┊ GC (mean ± σ):  9.66% ± 10.92%
# This seems to be related to Zygote!! Becuase when you use AutoForwardDiff
@benchmark gradient($loss_function_adjoint, $AutoForwardDiff(), $param_test)
# NOTE: AutoForwardDif IGNOGRES the adjoint and just does forwarddiff I think...?
# The benchmark is too similar the one with ForwardDiffSensitivity
# BenchmarkTools.Trial: 9178 samples with 1 evaluation per sample.
# Range (min … max):  500.522 μs …  16.718 ms  ┊ GC (min … max): 0.00% … 94.20%
# Time  (median):     522.883 μs               ┊ GC (median):    0.00%
# Time  (mean ± σ):   542.809 μs ± 286.628 μs  ┊ GC (mean ± σ):  0.90% ±  1.70%
# Dual numbers on just solves, no adjoint is called!
# @enter gradient(loss_function_adjoint, AutoForwardDiff(), param_test)
# Very different from using AutoZygote, which does call the adjoint (debugger crashes)
# @enter gradient(loss_function_adjoint, AutoZygote(), param_test)

@benchmark gradient($loss_function_adjoint, $AutoReverseDiff(), $param_test)
# ReverseDiff also quite slow...
# BenchmarkTools.Trial: 3673 samples with 1 evaluation per sample.
# BenchmarkTools.Trial: 3767 samples with 1 evaluation per sample.
#  Range (min … max):  992.173 μs … 32.616 ms  ┊ GC (min … max):  0.00% … 94.00%
#  Time  (median):       1.170 ms              ┊ GC (median):     0.00%
#  Time  (mean ± σ):     1.324 ms ±  2.041 ms  ┊ GC (mean ± σ):  10.50% ±  6.54%
# @benchmark gradient($loss_function_adjoint, $AutoMooncake(), $param_test)
# Fails as of now, not trivial to get reverse mode AD working fast...

# Add sensitivity equations with reverse mode AD as a working option
function loss_function_sensitivity_equations(p)
    return loss_function_test(p; sensealg=ForwardSensitivity(; autodiff=true))
end
gradient(loss_function_sensitivity_equations, AutoZygote(), param_test)
@benchmark gradient($loss_function_sensitivity_equations, $AutoZygote(), $param_test)
# Direct sensitivities
# https://docs.sciml.ai/SciMLSensitivity/stable/tutorials/direct_sensitivity/
prob_sens = ODEForwardSensitivityProblem(
    calculate_fluxes_test!, u0_test, t_span_test, param_test
)
sol_sens = solve(prob_sens, Tsit5())
@benchmark solve($prob_sens, $Tsit5()) # quite slow
x, jac_sol = extract_local_sensitivities(sol_sens)
sens_rsmin = jac_sol[end]
plot(sol_sens.t, sens_rsmin'; ylabel="∂y/∂rₛₘᵢₙ", xlabel="t")

# Direct reverse mode AD
loss_function_discrete(p) = loss_function_test(p; sensealg=ReverseDiffAdjoint())
@benchmark gradient($loss_function_discrete, $AutoReverseDiff(), $param_test)

# %% Easier loss function + solver
function loss_easy_choice(p, sensealg)
    return sum(solve(prob_test, Heun(); p=p, sensealg=sensealg))
end
loss_easy_forwardiff(p) = loss_easy_choice(p, ForwardDiffSensitivity())
loss_easy_gauss_adjoint(p) = loss_easy_choice(p, GaussAdjoint())
loss_easy_backsolve_adjoint(p) = loss_easy_choice(p, BacksolveAdjoint())
@benchmark gradient($loss_easy_forwardiff, $AutoForwardDiff(), $param_test)
# BenchmarkTools.Trial: 10000 samples with 1 evaluation per sample.
#  Range (min … max):  19.300 μs …   8.577 ms  ┊ GC (min … max): 0.00% … 98.72%
#  Time  (median):     21.000 μs               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   26.399 μs ± 143.777 μs  ┊ GC (mean ± σ):  9.33% ±  1.71%
@benchmark gradient($loss_easy_gauss_adjoint, $AutoZygote(), $param_test)
# BenchmarkTools.Trial: 5310 samples with 1 evaluation per sample.
#  Range (min … max):  804.304 μs …  12.673 ms  ┊ GC (min … max): 0.00% … 87.86%
#  Time  (median):     874.905 μs               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   937.636 μs ± 697.328 μs  ┊ GC (mean ± σ):  4.79% ±  5.94%
@benchmark gradient($loss_easy_backsolve_adjoint, $AutoZygote(), $param_test)
# BenchmarkTools.Trial: 5246 samples with 1 evaluation per sample.
#  Range (min … max):  803.305 μs …  12.753 ms  ┊ GC (min … max): 0.00% … 86.09%
#  Time  (median):     875.256 μs               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   949.442 μs ± 738.863 μs  ┊ GC (mean ± σ):  5.41% ±  6.42%
# Enzyme.gradient(Reverse, loss_easy_gauss_adjoint, param_test)

# %% Experimenting with Enzyme
function loss_function_enzyme(p, t_obs_test, observed_data, sensealg=BacksolveAdjoint())
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
loss_function_enzyme(param_test, t_obs_test, observed_data)
# gradient(
#     loss_function_enzyme,
#     AutoEnzyme(; mode=Enzyme.set_runtime_activity(Enzyme.Forward)),
#     param_test,
#     Constant(t_obs_test),
#     Constant(observed_data),
# )
# Craches as of now

# %% A simpler Enzyme experiment
# Only works with explicit solvers!
# https://sciml.ai/news/2024/08/25/rootfinding_enzyme/#direct_ode_support_with_enzyme
using StaticArrays
function lorenz!(du, u, p, t)
    du[1] = 10.0(u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    return du[3] = u[1] * u[2] - (8 / 3) * u[3]
end

const _saveat = SA[0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0]

function f(u0::Array{Float64})
    tspan = (0.0, 3.0)
    prob = ODEProblem{true,SciMLBase.FullSpecialize}(lorenz!, u0, tspan)
    sol = DiffEqBase.solve(
        prob, Tsit5(); saveat=_saveat, sensealg=DiffEqBase.SensitivityADPassThrough()
    )
    return sum(sol)
end;
u0 = [1.0; 0.0; 0.0]
d_u0 = zeros(3)
y = zeros(13)
dy = zeros(13)
f(u0)

grad_enzyme = gradient(f, AutoEnzyme(), u0)

# Make a "normal" version for use regular functioanality
function f_normal(u0)
    tspan = (0.0, 3.0)
    prob = ODEProblem{true,SciMLBase.FullSpecialize}(lorenz!, u0, tspan)
    sol = DiffEqBase.solve(prob, Tsit5(); saveat=_saveat, sensealg=ForwardDiffSensitivity())
    return sum(sol)
end
grad_forward = gradient(f_normal, AutoForwardDiff(), u0)
grad_finite_diff = gradient(f_normal, AutoFiniteDiff(), u0)
grad_forward ≈ grad_enzyme #NOT THE SAME! Enzyme wrong I think...
# Benchmark the difference
@benchmark gradient(f, AutoEnzyme(), u0) # 282 µs
@benchmark gradient(f_normal, AutoForwardDiff(), u0) # 8.5 µs
