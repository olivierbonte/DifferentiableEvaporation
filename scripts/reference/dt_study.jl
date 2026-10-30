# Δt study and time-stepper comparison (PR 0 of the Terrarium port), see
# docs/dev/2026-10-01-PLAN_dt_study_implicit_euler_smoothing.md.
# Every run is compared with a Kvaerno5 "truth" (abstol = reltol = 1e-10) of the same threshold
# treatment: the numerical error, as against the exact solution of the same (smoothed) model
# equations in Clark & Kavetski (2010, §6.1). The smoothing bias, the smoothed truth against the
# hard truth, is a separate row. Saves the results (one row per run) in
# data/sims/reference/dt_study_*.jld2 (`produce_or_load`, see `rerun`) and as
# data/sims/reference/dt_study.csv, and plots the relative RMSE against the step size
# (plots/reference/dt_study_rrmse.png) and against the cost (plots/reference/dt_study_rrmse_cost.png).
using DrWatson
@quickactivate "DifferentiableEvaporation"
include(scriptsdir("reference", "setup.jl"))
using CairoMakie
using OrdinaryDiffEqSDIRK: ImplicitEuler, KenCarp4, Kvaerno5
using Statistics

## Settings

const WINDOW = 1800.0 # s, forcing window
const DTS = WINDOW ./ [1, 2, 4, 8, 16, 32] # fixed steps: 1 to 32 sub-steps per window
# SUMMA-SUNDIALS production tolerance; also the Newton tolerance of fixed-step ImplicitEuler
const TOL = 1e-6
const TRUTH_TOL = 1e-10

# Accuracy criterion: RMSE of the states at the window bounds and of the window means of λE_tot,
# relative error of the cumulative latent heat
const MAX_RMSE_W = 0.004 # m³ m⁻³
const MAX_RMSE_LE = 5.0 # W m⁻²
const MAX_REL_CUM = 0.01 # -
# Reference line of the plot, not a criterion: SUMMA-SUNDIALS at 1e-6 (Spiteri et al., 2024)
const REFERENCE_RRMSE = 0.1

const TREATMENTS = ["hard" => HardThresholds(), "smooth" => KavetskiSmoothing()]

# Fixed-step backward Euler (N sub-steps per window), with the Newton settings of
# `fixed_step_newton`
backward_euler() = ImplicitEuler(; nlsolve=fixed_step_newton(), autodiff=AutoForwardDiff())

# The runs: fixed-step solvers have a step `dt` [s], adaptive solvers a tolerance `tol`
# (abstol = reltol); the other is missing. `kwargs` are passed to `solve`.
const RUNS = []
for dt in DTS
    push!(RUNS, (solver="Heun", dt=dt, tol=missing, alg=Heun(), kwargs=(adaptive=false, dt=dt)))
end
for dt in DTS
    backward_euler_kwargs = (adaptive=false, dt=dt, abstol=TOL, reltol=TOL)
    push!(RUNS, (solver="ImplicitEuler", dt=dt, tol=missing, alg=backward_euler(), kwargs=backward_euler_kwargs))
end
const ADAPTIVE_RUNS = [
    # (solver, tolerance, algorithm, extra kwargs)
    ("ImplicitEuler-adaptive", TOL, ImplicitEuler(; autodiff=AutoForwardDiff()), (dtmax=WINDOW,)),
    ("KenCarp4", TOL, KenCarp4(; autodiff=AutoForwardDiff()), (dtmax=WINDOW,)),
    ("Heun-adaptive", TOL, Heun(), (;)),
    ("Tsit5", TOL, Tsit5(), (;)),
]
for (solver, tol, alg, extra_kwargs) in ADAPTIVE_RUNS
    kwargs = merge((abstol=tol, reltol=tol), extra_kwargs)
    push!(RUNS, (solver=solver, dt=missing, tol=tol, alg=alg, kwargs=kwargs))
end

# Columns of the results (and dt_study.csv). Missing where they do not apply: the solver statistics and errors of
# runs that fail, and the solver statistics of the "smoothing bias" row (smoothed truth against
# the hard truth)
const CSV_COLUMNS = [
    :treatment, :solver, :dt, :tol, :retcode,
    # Cost: accepted and rejected steps, RHS evaluations, Jacobians, Newton iterations, the median
    # accepted step [s] of the adaptive solvers, and the RHS evaluations per forcing window
    :naccept, :nreject, :nf, :njacs, :nnonliniter, :median_dt, :nf_per_window,
    # Accuracy criterion and the errors it is based on
    :accurate, :rmse_w_1, :rmse_w_2, :rmse_λE_tot, :rel_cum_λE,
    # Relative RMSE (plotted)
    :rrmse_w_1, :rrmse_λE_tot,
]

## Errors

# The window means of λE_tot are NaN at the last time stamp, where no window starts
without_nan(x) = filter(!isnan, collect(x))
rmse(output, truth) = metric_skipnan(RMSE(), output, truth)
rrmse(output, truth) = metric_skipnan(RRMSE(), output, truth)
# Cumulative latent heat [J m⁻²] over the evaluation period, from the window means of λE_tot. As
# cumulative ET, but in energy: ET itself would need λ, which depends on the air temperature.
cumulative_latent_heat(output) = sum(without_nan(output.λE_tot)) * WINDOW
rel_error(c, c_truth) = abs(c - c_truth) / abs(c_truth)

"""
    errors(output, truth)

Errors of the run `output` against `truth`, both output `Dataset`s over the evaluation period
([`evaluation_period`](@ref)): the RMSE of w₁ and w₂ at the window bounds and of the window means
of λE_tot, the relative error of the cumulative latent heat, whether these meet the accuracy
criterion, and the relative RMSE of w₁ and λE_tot.
"""
function errors(output, truth)
    e = (
        rmse_w_1=rmse(output.w_1, truth.w_1),
        rmse_w_2=rmse(output.w_2, truth.w_2),
        rmse_λE_tot=rmse(output.λE_tot, truth.λE_tot),
        rel_cum_λE=rel_error(cumulative_latent_heat(output), cumulative_latent_heat(truth)),
        rrmse_w_1=rrmse(output.w_1, truth.w_1),
        rrmse_λE_tot=rrmse(output.λE_tot, truth.λE_tot),
    )
    accurate =
        all(isfinite, e) &&
            e.rmse_w_1 < MAX_RMSE_W &&
            e.rmse_w_2 < MAX_RMSE_W &&
            e.rmse_λE_tot < MAX_RMSE_LE &&
            e.rel_cum_λE < MAX_REL_CUM
    return merge((accurate=accurate,), e)
end

## Cost

"""
    nf_per_window(stats, case)

RHS evaluations per 30-min forcing window, the cost measure of Clark & Kavetski (2010, §6.1.3).
`stats.nf` counts every call of the RHS, those of the ForwardDiff Jacobians and Newton
iterations included (OrdinaryDiffEq.jl issue #2537); the linear solves are left out.
"""
function nf_per_window(stats, case)
    n_windows = length(case.saveat) - 1
    return stats.nf / n_windows
end

## Runs

"""
    truth(case, thresholds)

Kvaerno5 at abstol = reltol = `TRUTH_TOL` over the evaluation period. Throws if it fails, as no
run can be evaluated without it.
"""
function truth(case, thresholds)
    model = run_case(
        case, Kvaerno5(; autodiff=AutoForwardDiff()); thresholds, abstol=TRUTH_TOL, reltol=TRUTH_TOL
    )
    sol = model.sol
    if sol.retcode != ReturnCode.Success
        error("No truth for $thresholds: $(sol.retcode) at $(unix2datetime(sol.t[end]))")
    end
    @info "Truth for $thresholds: $(sol.stats.naccept) steps"
    return evaluation_period(model)
end

"""
    step_recorder(step_sizes)

Callback that pushes the size of every accepted step to `step_sizes`, without saving the
solution at every step. The recording happens in the condition, which never triggers the
(empty) affect.
"""
function step_recorder(step_sizes)
    function record_step(u, t, integrator)
        push!(step_sizes, t - integrator.tprev)
        return false
    end
    do_nothing(integrator) = nothing
    return DiscreteCallback(record_step, do_nothing; save_positions=(false, false))
end

"""
    evaluate(case, truth, run; thresholds)

One row of the results: `run` with its solver statistics and, if it succeeds, its cost per window
and errors against `truth`. A run that throws keeps the start of the error message as `retcode`.
"""
function evaluate(case, truth, run; thresholds)
    run_info = (solver=run.solver, dt=run.dt, tol=run.tol)
    is_adaptive = !ismissing(run.tol)
    step_sizes = Float64[]
    # Merged by `solve` with the `window_callback` of the problem
    callback_kwargs = is_adaptive ? (; callback=step_recorder(step_sizes)) : (;)
    model = try
        run_case(case, run.alg; thresholds, callback_kwargs..., run.kwargs...)
    catch err
        message = first(sprint(showerror, err), 200)
        return merge(run_info, (retcode="Error: " * message,))
    end

    sol = model.sol
    stats = sol.stats
    solver_statistics = (
        retcode=string(sol.retcode),
        naccept=stats.naccept,
        nreject=stats.nreject,
        nf=stats.nf,
        njacs=stats.njacs,
        nnonliniter=stats.nnonliniter,
        median_dt=isempty(step_sizes) ? missing : median(step_sizes),
    )
    row = merge(run_info, solver_statistics)
    succeeded = sol.retcode == ReturnCode.Success
    if !succeeded
        return row # no cost per window (it stopped early) and no errors
    end
    run_output = evaluation_period(model)
    cost = (nf_per_window=nf_per_window(stats, case),)
    return merge(row, cost, errors(run_output, truth))
end

## Study

# Only runs the study if `rerun` or if no saved results exist for this site, year and tolerances
# (set `rerun = true` after changing e.g. the runs or the reference case)
rerun = false
study_data, _ = produce_or_load(
    (; site=SITE, year=EVAL_YEAR, tol=TOL, truth_tol=TRUTH_TOL), datadir("sims", "reference");
    prefix="dt_study", force=rerun,
) do _
    case = reference_case()
    truths = Dict()
    for (treatment, thresholds) in TREATMENTS
        truths[treatment] = truth(case, thresholds)
    end

    # Rows have different columns (failed runs have no errors): collect with cols=:union, then
    # order the columns
    results = DataFrame()
    for (treatment, thresholds) in TREATMENTS
        for run in RUNS
            @info "$treatment: $(run.solver), dt = $(run.dt), tol = $(run.tol)"
            row = evaluate(case, truths[treatment], run; thresholds)
            push!(results, merge((treatment=treatment,), row); cols=:union)
        end
    end
    # The smoothing bias: the smoothed truth against the hard truth
    bias = errors(truths["smooth"], truths["hard"])
    push!(results, merge((treatment="smooth", solver="smoothing bias"), bias); cols=:union)
    select!(results, CSV_COLUMNS)
    # Also as CSV, to read without Julia
    CSV.write(datadir("sims", "reference", "dt_study.csv"), results)
    return @strdict(results)
end
@unpack results = study_data

## Plot

# Physical size, font size and resolution of the saved figures, see
# https://docs.makie.org/stable/how-to/match-figure-size-font-sizes-and-dpi
inch = 96 # figure units (CSS px) per inch
pt = 4 / 3
cm = inch / 2.54
dpi = 300
# AGU full page width (two-column figure) with ~8 pt lettering
figure_size = (17cm, 10.5cm)
fontsize = 8pt

# A color per time stepper, a marker per threshold treatment (Wong palette)
const SOLVERS = unique([run.solver for run in RUNS])
const FIXED_STEP_SOLVERS = ["Heun", "ImplicitEuler"]
const SOLVER_COLORS = Dict(zip(SOLVERS, Makie.wong_colors()))
const TREATMENT_MARKERS = Dict("hard" => :rect, "smooth" => :circle)

"""
    plot_rrmse(results, xaxis, file)

Relative RMSE of λE_tot and w₁ against the step size (`xaxis = :step`: Δt for fixed-step
solvers, the median Δt for adaptive ones) or the cost (`xaxis = :cost`: RHS evaluations per
forcing window), saved as `plots/reference/<file>`. Fixed-step solvers are lines over their
steps, labelled with N (Δt = WINDOW / N) against the cost; adaptive runs are markers. The
dashed line is `REFERENCE_RRMSE`.
"""
function plot_rrmse(results, xaxis, file)
    runs = results[results.solver .!= "smoothing bias", :]
    metrics = [
        :rrmse_λE_tot => "RRMSE λE_tot [-]",
        :rrmse_w_1 => "RRMSE w₁ [-]",
    ]
    if xaxis == :step
        xlabel = "Δt [s] (adaptive: median)"
    else
        xlabel = "RHS evaluations per 30-min window"
    end

    fig = Figure(; size=figure_size, fontsize)
    # One panel per metric, with a shared x-axis labelled under the last panel
    panels = []
    for (row, (_, ylabel)) in enumerate(metrics)
        ax = CairoMakie.Axis(fig[row, 1]; xscale=log10, yscale=log10, ylabel) # Axis is also exported by ComponentArrays
        push!(panels, ax)
    end
    linkxaxes!(panels...)
    for ax in panels[1:(end-1)]
        hidexdecorations!(ax; grid=false)
    end
    panels[end].xlabel = xlabel

    for (ax, (column, _)) in zip(panels, metrics)
        for series in groupby(dropmissing(runs, column), [:solver, :treatment, :tol])
            solver = first(series.solver)
            treatment = first(series.treatment)
            tol = first(series.tol)
            is_fixed_step = ismissing(tol)
            if xaxis == :cost
                x = Float64.(series.nf_per_window)
            elseif is_fixed_step
                x = Float64.(series.dt)
            else
                x = Float64.(series.median_dt)
            end
            y = series[!, column]
            color = SOLVER_COLORS[solver]
            marker = TREATMENT_MARKERS[treatment]

            if is_fixed_step
                scatterlines!(ax, x, y; color, marker, markersize=9, linewidth=1)
                # The treatments overlap, so label the N of the hard thresholds only
                if xaxis == :cost && treatment == "hard"
                    labels = ["N=$(round(Int, WINDOW / dt))" for dt in series.dt]
                    text!(ax, x, y; text=labels, color, fontsize=6pt, offset=(4, 4))
                end
            else
                scatter!(ax, x, y; color, marker, markersize=11)
            end
        end
        hlines!(ax, [REFERENCE_RRMSE]; color=:gray, linestyle=:dash)
    end

    # Legend in two groups, as the styles combine per series
    solver_elements = []
    solver_labels = String[]
    for solver in SOLVERS
        color = SOLVER_COLORS[solver]
        if solver in FIXED_STEP_SOLVERS
            push!(solver_elements, [LineElement(; color), MarkerElement(; marker=:circle, color)])
            push!(solver_labels, "$solver (fixed Δt)")
        else
            push!(solver_elements, MarkerElement(; marker=:circle, color, markersize=11))
            # The fixed-step runs are labelled as such, so no suffix is needed for the adaptive ones
            push!(solver_labels, replace(solver, "-adaptive" => ""))
        end
    end
    treatment_elements = []
    treatment_labels = String[]
    for (treatment, _) in TREATMENTS
        push!(treatment_elements, MarkerElement(; marker=TREATMENT_MARKERS[treatment], color=:gray40, markersize=11))
        push!(treatment_labels, treatment)
    end
    Legend(
        fig[length(metrics)+1, 1],
        [solver_elements, treatment_elements],
        [solver_labels, treatment_labels],
        ["Time stepper", "Threshold treatment"];
        orientation=:horizontal, tellheight=true, labelsize=fontsize, titlesize=fontsize, nbanks=2,
    )
    safesave(plotsdir("reference", file), fig; px_per_unit=dpi / inch)
    return nothing
end

plot_rrmse(results, :step, "dt_study_rrmse.png")
plot_rrmse(results, :cost, "dt_study_rrmse_cost.png")
