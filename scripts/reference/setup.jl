# Shared setup for the Δt study (PR 0 of the Terrarium port): reference cases and runs.
using DrWatson
@quickactivate "DifferentiableEvaporation"
using CSV
using DataFrames
using Dates
using EvaporationModel
using OrdinaryDiffEq

include(srcdir("site_data.jl")) # site_forcing, site_parameters
include(joinpath(@__DIR__, "synthetic_forcing.jl"))

const REFERENCE_DIR = datadir("reference")
const SAVE_DT = 1800.0 # s

site_u0(p) = [0.6 * p.w_fc, 0.6 * p.w_fc, 1e-4]

"""
    reference_cases()

The reference runs: name => (forcings, parameters, u0, tspan, tstops). The `tstops` lie just
before each jump of the piecewise-constant forcing (see `tstops_before`).
"""
function reference_cases()
    cases = Dict{String,Any}("synthetic_1day" => synthetic_case())
    p = site_parameters("BE-Bra")
    for start_date in (DateTime(2010, 3, 15), DateTime(2010, 7, 7))
        # piecewise constant over each 30-min forcing window
        forcings, t = site_forcing(site_ec_file("BE-Bra"), start_date, start_date + Day(10))
        name = "BE-Bra_$(Date(start_date))_10day"
        cases[name] = (
            forcings=forcings, p=p, u0=site_u0(p), tspan=(t[1], t[end]),
            tstops=tstops_before(t[2:end]),
        )
    end
    return cases
end

save_times(case) = collect(case.tspan[1]:SAVE_DT:case.tspan[2])

function build_model(case; thresholds=HardThresholds())
    model = ProcessBasedModel{Float64}(;
        forcings=case.forcings,
        parameters=case.p,
        t_span=case.tspan,
        u0=copy(case.u0),
        saveat=save_times(case),
        tstops=case.tstops,
        thresholds=thresholds,
    )
    EvaporationModel.initialize!(model)
    return model
end

"""
    run_case(case, alg; thresholds=HardThresholds(), kwargs...)

Solve one reference case with `alg` and the threshold treatment `thresholds`, saving states
and diagnostics every `SAVE_DT`. Returns the solved model; `kwargs` are passed to `solve`.
"""
function run_case(case, alg; thresholds=HardThresholds(), kwargs...)
    model = build_model(case; thresholds)
    EvaporationModel.solve!(model; alg=alg, maxiters=10^8, kwargs...)
    return model
end

"""
    to_dataframe(model)

States, diagnostics and forcings at the save times.
"""
function to_dataframe(model)
    sol = model.sol
    t = sol.t
    f = model.forcings
    df = DataFrame(; time=unix2datetime.(t), t_unix=t)
    df.w_1 = sol[1, :]
    df.w_2 = sol[2, :]
    df.w_r = sol[3, :]
    diag = DataFrame(model.diagnostics.saveval)
    df = hcat(df, diag)
    for name in (:P, :T_a, :u_a, :p_a, :VPD_a, :SW_in, :R_n, :LAI)
        df[!, Symbol("forcing_", name)] = getproperty(f, name).(t)
    end
    return df
end

function git_info()
    root = projectdir()
    hash = readchomp(`git -C $root rev-parse HEAD`)
    dirty = !isempty(readchomp(`git -C $root status --porcelain -- src`))
    return hash, dirty
end
