# Δt study and time-stepper comparison (PR 0 of the Terrarium port), see
# docs/dev/2026-10-01-PLAN_dt_study_implicit_euler_smoothing.md.
# Fixed-step Heun and backward Euler (BE-N, BackTracking Newton), adaptive ImplicitEuler,
# KenCarp4 and Tsit5, for hard and smoothed thresholds, against a Kvaerno5 truth at
# abstol = reltol = 1e-10 per treatment. Writes data/reference/dt_study.{csv,md,png}.
include(joinpath(@__DIR__, "setup.jl"))
using CairoMakie
using ForwardDiff
using LinearAlgebra: I, cond
using LineSearches: BackTracking, LineSearchException
using OrdinaryDiffEqNonlinearSolve: NLNewton
using OrdinaryDiffEqSDIRK: ImplicitEuler, KenCarp4, Kvaerno5
using Printf
using Statistics

const WINDOW = 1800.0 # s, forcing window
const HEUN_DTS = [1800.0, 900.0, 600.0, 300.0, 120.0, 60.0, 30.0, 10.0, 5.0, 2.0, 1.0]
const BE_NS = [1, 2, 4, 8, 16, 64] # sub-steps per forcing window
const TOL = 1e-6 # SUMMA-SUNDIALS production tolerance; also sets the Newton tolerance of BE-N
const TRUTH_TOL = 1e-10
const MAX_NEWTON_ITER = 30

# Accuracy criterion (Decision 1): RMSE at the save times, relative error of cumulative fluxes
const MAX_RMSE_W = 0.004 # m³ m⁻³
const MAX_RMSE_LE = 5.0 # W m⁻²
const MAX_REL_CUM = 0.01 # -
const CUM_FLOOR = 0.01 # kg m⁻², floor of the denominator of the relative cumulative error

const TREATMENTS = [
    "hard" => HardThresholds(),
    "smooth_small" => KavetskiSmoothing(; s_w=1e-3, s_r=0.01, s_E=2.0),
    "smooth_large" => KavetskiSmoothing(; s_w=5e-3, s_r=0.05, s_E=10.0),
]

be_newton() = NLNewton(; max_iter=MAX_NEWTON_ITER, relax=BackTracking())
backward_euler() = ImplicitEuler(; nlsolve=be_newton(), autodiff=AutoForwardDiff())

# Adaptive runs: (solver, abstol = reltol, alg, extra kwargs)
const ADAPTIVE = [
    ("ImplicitEuler", TOL, ImplicitEuler(; autodiff=AutoForwardDiff()), (dtmax=WINDOW,)),
    ("KenCarp4", TOL, KenCarp4(; autodiff=AutoForwardDiff()), (dtmax=WINDOW,)),
    ("Tsit5", 1e-4, Tsit5(), (;)),
    ("Tsit5", 1e-6, Tsit5(), (;)),
]

"""
    solver_configs()

`(solver, setting, alg, kwargs)` per run; `setting` is Δt [s] for fixed-step solvers and
abstol = reltol for adaptive ones.
"""
solver_configs() = vcat(
    [("Heun", dt, Heun(), (adaptive=false, dt=dt)) for dt in HEUN_DTS],
    [
        ("BE-N", WINDOW / N, backward_euler(), (adaptive=false, dt=WINDOW / N, abstol=TOL, reltol=TOL))
        for N in BE_NS
    ],
    [(solver, tol, alg, (; abstol=tol, reltol=tol, kw...)) for (solver, tol, alg, kw) in ADAPTIVE],
)

const RESULT_COLUMNS = (
    :retcode, :accurate, :wall_time, :naccept, :nreject, :nf, :nnonliniter, :nnonlinconvfail,
    :newton_per_step, :min_dt, :median_dt,
    :rmse_w_1, :rmse_w_2, :rmse_λE_tot, :rmse_H_tot, :rel_cum_ET, :rel_cum_runoff_drainage,
    :max_w_1, :max_w_2, :max_w_r, :max_λE_tot,
)

# Records the accepted step sizes without saving the solution at every step
step_recorder(dts) = DiscreteCallback(
    (u, t, integrator) -> (push!(dts, t - integrator.tprev); false),
    integrator -> nothing;
    save_positions=(false, false),
)

empty_result() = Dict{Symbol,Any}(c => missing for c in RESULT_COLUMNS)

trapz(t, y) = sum((t[i + 1] - t[i]) * (y[i] + y[i + 1]) / 2 for i in 1:(length(t) - 1))

"""
    cumulative_fluxes(df, p)

Cumulative evapotranspiration and runoff + drainage [kg m⁻²] over the run.
"""
function cumulative_fluxes(df, p)
    ET = df.E_t .+ df.E_i .+ df.E_s
    runoff_drainage = df.Q_s .+ ρ_w * p.d_2 .* df.K_2
    return trapz(df.t_unix, ET), trapz(df.t_unix, runoff_drainage)
end

rmse(a, b) = sqrt(mean(abs2, a .- b))
rel_cum(c, c_truth) = abs(c - c_truth) / max(abs(c_truth), CUM_FLOOR)

"""
    compare(df, truth_df, p)

Error metrics of `df` against `truth_df` at the save times, and whether they meet the criterion.
"""
function compare(df, truth_df, p)
    r = Dict{Symbol,Any}()
    for var in (:w_1, :w_2, :λE_tot, :H_tot)
        r[Symbol(:rmse_, var)] = rmse(df[!, var], truth_df[!, var])
    end
    for var in (:w_1, :w_2, :w_r, :λE_tot)
        r[Symbol(:max_, var)] = maximum(abs, df[!, var] .- truth_df[!, var])
    end
    ET, RD = cumulative_fluxes(df, p)
    ET_truth, RD_truth = cumulative_fluxes(truth_df, p)
    r[:rel_cum_ET] = rel_cum(ET, ET_truth)
    r[:rel_cum_runoff_drainage] = rel_cum(RD, RD_truth)
    r[:accurate] =
        all(isfinite, values(r)) &&
        r[:rmse_w_1] < MAX_RMSE_W &&
        r[:rmse_w_2] < MAX_RMSE_W &&
        r[:rmse_λE_tot] < MAX_RMSE_LE &&
        r[:rel_cum_ET] < MAX_REL_CUM &&
        r[:rel_cum_runoff_drainage] < MAX_REL_CUM
    return r
end

function evaluate(case, truth_df, alg, kwargs; thresholds)
    result = empty_result()
    dts = Float64[]
    local model
    try
        result[:wall_time] = @elapsed model = run_case(
            case, alg; thresholds, callback=step_recorder(dts), kwargs...
        )
    catch err
        result[:retcode] = "Error: " * first(sprint(showerror, err), 80)
        return result
    end
    sol = model.sol
    s = sol.stats
    result[:retcode] = string(sol.retcode)
    result[:naccept], result[:nreject], result[:nf] = s.naccept, s.nreject, s.nf
    result[:nnonliniter], result[:nnonlinconvfail] = s.nnonliniter, s.nnonlinconvfail
    nsteps = s.naccept + s.nreject
    result[:newton_per_step] = nsteps > 0 && s.nnonliniter > 0 ? s.nnonliniter / nsteps : missing
    dts = filter(>(1e-2), dts) # drop the ε-step from the last tstop to the end of the run
    if !isempty(dts)
        result[:min_dt], result[:median_dt] = minimum(dts), median(dts)
    end
    finished = sol.retcode == ReturnCode.Success && length(sol.t) == nrow(truth_df)
    finished || return result
    merge!(result, compare(to_dataframe(model), truth_df, case.p))
    return result
end

"""
    truth(case, thresholds)

Kvaerno5 at abstol = reltol = 1e-10, with a ForwardDiff Jacobian, or a finite-difference
Jacobian if that fails (the hard f_wet has an infinite slope at w_r = 0).
"""
function truth(case, thresholds)
    for autodiff in (AutoForwardDiff(), AutoFiniteDiff())
        alg = Kvaerno5(; autodiff)
        model = try
            run_case(case, alg; thresholds, abstol=TRUTH_TOL, reltol=TRUTH_TOL)
        catch err
            @warn "Truth with $(nameof(typeof(autodiff))) failed" exception = err
            continue
        end
        sol = model.sol
        if sol.retcode == ReturnCode.Success
            info = "Kvaerno5, $(nameof(typeof(autodiff))) Jacobian, $(sol.stats.naccept) steps"
            return to_dataframe(model), info
        end
        @warn "Truth with $(nameof(typeof(autodiff))) returned $(sol.retcode)"
    end
    return error("No converged truth")
end

"""
    active_thresholds(model, u, p, t)

Thresholds of `model` near which the state `u` lies at time `t`.
"""
function active_thresholds(model, u, p, t)
    d = model.f_diagnostics(u, p, t)
    w_1, w_2, w_r = u
    near(x, x_0, s) = abs(x - x_0) < s
    flags = (
        "w_r ≈ 0" => w_r < 0.01 * d.w_rmax,
        "w_1 ≈ w_wp" => near(w_1, p.w_wp, 5e-3),
        "w_2 ≈ w_fc" => near(w_2, p.w_fc, 5e-3),
        "w ≥ w_sat" => max(w_1, w_2) > p.w_sat - 5e-3,
        "λE_pc ≈ 0 (condensation)" => near(d.λE_pc, 0, 10),
        "R_ns ≈ 0 (G ratio)" => near(d.R_ns, 0, 10),
    )
    return join([name for (name, on) in flags if on], ", ")
end

"""
    diagnose_be(case, thresholds; N=8)

Steps BE-N with the BackTracking Newton through `case` and stops at the first failing step.
Returns a summary with the Newton iteration counts per step and, on failure, the state,
the active thresholds and the condition numbers of the ODE Jacobian J and of I - Δt J.
"""
function diagnose_be(case, thresholds; N=8)
    dt = WINDOW / N
    model = build_model(case; thresholds)
    integrator = init(
        model.prob, backward_euler();
        adaptive=false, dt, abstol=TOL, reltol=TOL, tstops=model.tstops,
        save_everystep=false, maxiters=10^8,
    )
    iters = Int[]
    failure = nothing
    while integrator.t < case.tspan[2]
        t_prev, u_prev = integrator.t, copy(integrator.u)
        n_prev = integrator.stats.nnonliniter
        retcode = try
            step!(integrator)
            OrdinaryDiffEq.SciMLBase.check_error(integrator)
        catch err
            err isa LineSearchException || rethrow()
            "LineSearchException"
        end
        push!(iters, integrator.stats.nnonliniter - n_prev)
        if retcode != ReturnCode.Success
            J = ForwardDiff.jacobian((du, u) -> model.f(du, u, case.p, t_prev), zeros(3), u_prev)
            failure = (
                retcode=string(retcode),
                time=unix2datetime(t_prev),
                u=u_prev,
                newton_iter=iters[end],
                thresholds=active_thresholds(model, u_prev, case.p, t_prev),
                cond_J=cond(J),
                cond_W=cond(I - dt * J),
            )
            break
        end
    end
    return (
        N=N,
        dt=dt,
        nsteps=length(iters),
        mean_iter=mean(iters),
        max_iter=maximum(iters),
        n_near_max=count(>=(0.8 * MAX_NEWTON_ITER), iters),
        failure=failure,
    )
end

function main()
    mkpath(REFERENCE_DIR)
    cases = reference_cases()
    rows = DataFrame([c => Any[] for c in (:case, :treatment, :solver, :setting, RESULT_COLUMNS...)])
    truth_info = Dict{Tuple{String,String},String}()
    diagnostics = Dict{Tuple{String,String},Any}()
    for name in sort(collect(keys(cases)))
        case = cases[name]
        truths = Dict{String,DataFrame}()
        for (treatment, thresholds) in TREATMENTS
            @info "Truth for $name, $treatment"
            truths[treatment], truth_info[(name, treatment)] = truth(case, thresholds)
            if treatment != "hard" # smoothing bias: smoothed truth against hard truth
                r = empty_result()
                merge!(r, compare(truths[treatment], truths["hard"], case.p))
                push!(rows, merge(r, Dict(
                    :case => name, :treatment => treatment, :solver => "smoothing bias",
                    :setting => missing, :retcode => "Success",
                )))
            end
            @info "BE-8 diagnostic for $name, $treatment"
            diagnostics[(name, treatment)] = diagnose_be(case, thresholds)
            for (solver, setting, alg, kwargs) in solver_configs()
                @info "$name, $treatment: $solver $setting"
                r = evaluate(case, truths[treatment], alg, kwargs; thresholds)
                push!(rows, merge(r, Dict(
                    :case => name, :treatment => treatment, :solver => solver, :setting => setting,
                )))
            end
        end
    end
    CSV.write(joinpath(REFERENCE_DIR, "dt_study.csv"), rows)
    write_markdown(rows, truth_info, diagnostics)
    plot_errors(rows)
    return rows, diagnostics
end

fmt(x::Real) = @sprintf("%.2e", x)
fmt(x::Integer) = string(x)
fmt(x::Bool) = x ? "yes" : "no"
fmt(::Missing) = "–"
fmt(x) = string(x)

const FIXED_STEP = ("Heun", "BE-N")

"""
Largest accurate Δt per fixed-step solver, and whether each adaptive run is accurate.
"""
function summary_table(io, rows)
    header = vcat(["$s Δt [s]" for s in FIXED_STEP], ["$s $tol" for (s, tol) in ADAPTIVE])
    println(io, "| case | treatment | ", join(header, " | "), " |")
    println(io, "|", repeat("---|", 2 + length(header)))
    for g in groupby(rows, [:case, :treatment]; sort=false)
        best(solver) = begin
            ok = g[(g.solver .== solver) .& coalesce.(g.accurate .=== true, false), :setting]
            isempty(ok) ? "none" : fmt(maximum(ok))
        end
        adaptive(solver, tol) = begin
            r = g[(g.solver .== solver) .& (g.setting .=== tol), :]
            nrow(r) == 0 ? "–" : fmt(r.accurate[1] === true)
        end
        cells = vcat(
            [best(s) for s in FIXED_STEP],
            [adaptive(s, tol) for (s, tol) in ADAPTIVE],
        )
        println(io, "| $(g.case[1]) | $(g.treatment[1]) | ", join(cells, " | "), " |")
    end
    return println(io)
end

function write_markdown(rows, truth_info, diagnostics)
    cols = [:treatment, :solver, :setting, (c for c in RESULT_COLUMNS if c != :nf)...]
    io = IOBuffer()
    println(io, "# Δt study and time-stepper comparison\n")
    println(io, "Generated by `scripts/reference/dt_study.jl` on $(today()).")
    hash, dirty = git_info()
    println(io, "Model source at `$hash`$(dirty ? " (uncommitted changes in src)" : "").\n")
    println(
        io,
        "Errors are against Kvaerno5 at abstol = reltol = $TRUTH_TOL for the same threshold ",
        "treatment, at the 30-min save times. `accurate`: RMSE(w₁), RMSE(w₂) < $MAX_RMSE_W m³ m⁻³, ",
        "RMSE(λE_tot) < $MAX_RMSE_LE W m⁻², and relative error of cumulative ET and of cumulative ",
        "runoff + drainage < $(100MAX_REL_CUM)% (denominator floored at $CUM_FLOOR kg m⁻²). ",
        "`max_*` are maximum absolute errors (reported, not gating). `setting` is Δt [s] for ",
        "fixed-step solvers and abstol = reltol for adaptive ones. BE-N uses ",
        "`NLNewton(max_iter=$MAX_NEWTON_ITER, relax=BackTracking())` with abstol = reltol = $TOL. ",
        "`smoothing bias` rows compare the smoothed truth with the hard truth.\n",
    )
    println(io, "## Summary: largest accurate Δt\n")
    summary_table(io, rows)
    println(io, "## Truth solutions\n")
    for ((case, treatment), info) in sort(collect(truth_info))
        println(io, "- $case, $treatment: $info")
    end
    println(io, "\n## BE-8 diagnostic (BackTracking Newton)\n")
    println(io, "| case | treatment | steps | mean Newton iter. | max | steps ≥ $(round(Int, 0.8MAX_NEWTON_ITER)) | first failure |")
    println(io, "|---|---|---|---|---|---|---|")
    for ((case, treatment), d) in sort(collect(diagnostics))
        f = d.failure
        failure = isnothing(f) ? "none" :
            "$(f.retcode) at $(f.time), u = $(round.(f.u; sigdigits=4)), $(f.newton_iter) iter., " *
            "near: $(isempty(f.thresholds) ? "–" : f.thresholds), " *
            "cond(J) = $(fmt(f.cond_J)), cond(I − ΔtJ) = $(fmt(f.cond_W))"
        println(io, "| $case | $treatment | $(d.nsteps) | $(fmt(d.mean_iter)) | $(d.max_iter) | $(d.n_near_max) | $failure |")
    end
    println(io)
    for case in unique(rows.case)
        sub = rows[rows.case .== case, :]
        println(io, "## $case\n")
        println(io, "| ", join(string.(cols), " | "), " |")
        println(io, "|", repeat("---|", length(cols)))
        for r in eachrow(sub)
            println(io, "| ", join([fmt(r[c]) for c in cols], " | "), " |")
        end
        println(io)
    end
    write(joinpath(REFERENCE_DIR, "dt_study.md"), String(take!(io)))
    return nothing
end

function plot_errors(rows)
    cases = unique(rows.case)
    metrics = [(:rmse_w_1, "RMSE w₁ [m³ m⁻³]", MAX_RMSE_W), (:rmse_λE_tot, "RMSE λE [W m⁻²]", MAX_RMSE_LE)]
    linestyles = Dict("hard" => :solid, "smooth_small" => :dash, "smooth_large" => :dot)
    colors = Dict("Heun" => Makie.wong_colors()[1], "BE-N" => Makie.wong_colors()[2])
    fig = Figure(; size=(450 * length(cases), 350 * length(metrics)))
    for (i, (metric, label, criterion)) in enumerate(metrics), (j, case) in enumerate(cases)
        ax = CairoMakie.Axis( # Axis is also exported by ComponentArrays
            fig[i, j]; xscale=log10, yscale=log10, xlabel="Δt [s]", ylabel=label,
            title=i == 1 ? case : "",
        )
        for solver in FIXED_STEP, (treatment, _) in TREATMENTS
            sub = rows[(rows.case .== case) .& (rows.solver .== solver) .& (rows.treatment .== treatment), :]
            ok = .!ismissing.(sub[!, metric])
            any(ok) || continue
            scatterlines!(
                ax, Float64.(sub.setting[ok]), Float64.(sub[ok, metric]) .+ 1e-12;
                color=colors[solver], linestyle=linestyles[treatment], label="$solver, $treatment",
            )
        end
        hlines!(ax, [criterion]; color=:gray, linestyle=:dash)
        i == 1 && j == length(cases) && axislegend(ax; position=:rb, labelsize=9, merge=true)
    end
    save(joinpath(REFERENCE_DIR, "dt_study.png"), fig)
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
