# Shared setup for the Δt study (PR 0 of the Terrarium port): the reference case and its runs.
using DrWatson
@quickactivate "DifferentiableEvaporation"
using CSV
using DataFrames
using Dates
using EvaporationModel
using OrdinaryDiffEq

include(srcdir("site_data.jl")) # site_forcing, site_parameters

const SITE = "DE-Tha"
const EVAL_YEAR = 2018 # drought year: exercises the soil moisture stress thresholds
const SPINUP = Day(14) # simulated before EVAL_YEAR and left out of the error metrics
# Diagnostic integrated by the model, for its window means
const ACCUMULATED = (:λE_tot,)

"""
    reference_case()

One year of `SITE` forcing (`EVAL_YEAR`) after a `SPINUP`, as `(forcings, p, u0, tspan, saveat)`.
The states are saved at the 30-min window bounds `saveat`, which are also the bounds of the
forcing windows: every evaluation within a window, also at its end (where backward Euler
evaluates the forcing), uses that window's forcing (`EvaporationModel.window_callback`).
"""
function reference_case()
    p = site_parameters(SITE)
    eval_start = DateTime(EVAL_YEAR, 1, 1)
    # piecewise constant over each 30-min forcing window
    forcings, times, _ = site_forcing(
        site_ec_file(SITE), eval_start - SPINUP, DateTime(EVAL_YEAR, 12, 31, 23, 30)
    )
    w_init = 0.6 * p.w_fc # both soil layers start at 60% of field capacity
    u0 = [w_init, w_init, 1e-4] # w_1, w_2 [m³ m⁻³], w_r [kg m⁻²]: a nearly dry canopy
    return (forcings=forcings, p=p, u0=u0, tspan=times.span, saveat=times.bounds)
end

"""
    run_case(case, alg; thresholds=HardThresholds(), kwargs...)

Solve `case` with `alg` and the threshold treatment `thresholds`, saving the states at
`case.saveat` and the window means of the `ACCUMULATED` diagnostics between them, both also in
the YAXArrays `Dataset` `model.output`. Returns the solved model; `kwargs` are passed to `solve`.
"""
function run_case(case, alg; thresholds=HardThresholds(), kwargs...)
    model = ProcessBasedModel{Float64}(;
        forcings=case.forcings,
        parameters=case.p,
        t_span=case.tspan,
        u0=copy(case.u0),
        saveat=case.saveat,
        accumulated=ACCUMULATED,
        thresholds=thresholds,
    )
    EvaporationModel.initialize!(model)
    EvaporationModel.solve!(model; alg=alg, maxiters=10^8, yaxarray_output=true, kwargs...)
    return model
end

"""
    evaluation_period(model)

The output `Dataset` of the solved `model` (the states at the window bounds and the window means
of λE_tot, see `EvaporationModel.output_dataset`) over `EVAL_YEAR`, without the `SPINUP`.
"""
evaluation_period(model) = model.output[time=DateTime(EVAL_YEAR, 1, 1) .. DateTime(EVAL_YEAR + 1, 1, 1)]
