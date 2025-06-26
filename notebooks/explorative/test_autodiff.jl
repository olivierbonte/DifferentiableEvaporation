using DrWatson
@quickactivate "DifferentiableEvaporation"
using Bigleaf, EvaporationModel
using DifferentialEquations
using DifferentiationInterface
import Enzyme, ForwardDiff
using BenchmarkTools

# %% Test potential_et from Bigleaf with ForwardDiff
Tair = 30.0 # [°C]
pressure = 101.325 # [kPa]
Rn = 400.0 # [W m-2]
VPD = 0.5 # [kPa]
Ga = 1/30.0 # [m s-1]
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
penman_monteith(Tair + 273.15, pressure * 1000 , VPD * 1000, Rn, 1 / Ga, 1 / Ga * 4)

f_bis = x -> penman_monteith(x[1], x[2], x[3], x[4], x[5], x[6])[2]
gradient(f_bis, AutoForwardDiff(), x_trial_bis)

function f_ode(u, p, t)
    return -1 / (1000 * 1.5) * f_bis(p)
end
u0 = 0.5
tspan = (0.0, 20.0)
prob = ODEProblem(f_ode, u0, tspan, x_trial_bis)
@benchmark solve(prob, Rosenbrock23(;autodiff = AutoEnzyme(;function_annotation=Enzyme.Duplicated)))
@benchmark solve(prob, Rosenbrock23(;autodiff = AutoForwardDiff()))
