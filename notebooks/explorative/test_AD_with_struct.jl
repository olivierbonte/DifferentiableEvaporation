using DrWatson
@quickactivate "DifferentiableEvaporation"
using OrdinaryDiffEq, DataInterpolations, SciMLSensitivity, DifferentiationInterface

# Example forcing
tdata = 0:10
ydata = rand(length(tdata))
forcing = LinearInterpolation(ydata, tdata)
extra = (forcing=forcing,)
p_numeric = [1.0, 2.0]
t_span = (0.0, 10.0)
u0 = [0.0]

# NamedTuple input
p_nt = (parameters=p_numeric, forcing=forcing)
function f_normal!(du, u, p, t)
    α, β = p.parameters
    F = p.forcing(t)
    return du[1] = α * u[1] + β * F
end
prob_nt = ODEProblem(f_normal!, u0, t_span, p_nt)
sol_nt = solve(prob_nt)

function loss_nt(p)
    prob_tmp = remake(prob_nt; p=p)
    sol = solve(prob_tmp, Tsit5(); sensealg=InterpolatingAdjoint(; autojacvec=EnzymeVJP()))
    return sum(sol[end])
end
gradient(loss_nt, AutoZygote(), p_nt)
# Key error
# ERROR: Adjoint sensitivity analysis functionality requires being able to solve
# a differential equation defined by the parameter struct `p`. Thus while
# DifferentialEquations.jl can support any parameter struct type, usage
# with adjoint sensitivity analysis requires that `p` could be a valid
# type for being the initial condition `u0` of an array. This means that
# many simple types, such as `Tuple`s and `NamedTuple`s, will work as
# parameters in normal contexts but will fail during adjoint differentiation.

# To work around this issue for complicated cases like nested structs, look
# into defining `p` using `AbstractArray` libraries such as RecursiveArrayTools.jl
# or ComponentArrays.jl so that `p` is an `AbstractArray` with a concrete element type.

# If you have a non-standard type you wish to work with adjoint differentiation, you need
# to define the SciMLStructures.jl interface on that type. For more information, check out
# https://docs.sciml.ai/SciMLStructures/stable/example/ for an example.
#
function f_extra!(du, u, p, t, extra)
    α, β = p
    F = extra.forcing(t)
    return du[1] = α * u[1] + β * F
end
extra = (forcing=forcing,)

prob_extra = ODEProblem(
    (du, u, p, t) -> f_extra!(du, u, p, t, extra), u0, t_span, p_numeric
)
sol_extra = solve(prob_extra)

function loss_extra(p)
    prob_tmp = remake(prob_extra; p=p)
    sol = solve(prob_tmp, Tsit5(); sensealg=InterpolatingAdjoint(; autojacvec=EnzymeVJP()))
    return sum(sol[end])
end
gradient(loss_extra, AutoZygote(), p_numeric)

# Doing it this way, gradient calculation works just fine :)
