# Reproduction of Figure 4 of Kronenberg et al. (2026), https://doi.org/10.5194/gmd-19-5261-2026
# (subdaily BROOK90), with the model from src/EvaporationModel: energy balance components
# measured at ICOS DE-Tha and simulated, for 29 June - 2 July and 4 - 7 December 2024.
# As in the paper, 2023 is used as spin-up and 2024 for the evaluation.
#
# Forcing and evaluation data are made by
# scripts/data_processing/experimental/icos_hh_process.py (run that first).
using DrWatson
@quickactivate "DifferentiableEvaporation"
using ADTypes
using CairoMakie
using ComponentArrays
using Dates
using ErrorMetrics
using EvaporationModel
using LineSearches: BackTracking
using Optim
using OrdinaryDiffEq
using OrdinaryDiffEq: SciMLBase
using OrdinaryDiffEqNonlinearSolve: NLNewton
using Statistics
using YAXArrays

include(srcdir("site_data.jl")) # site_forcing, site_parameters

FT = Float64
site = "DE-Tha"
file_ec = datadir("exp_pro", "eddy_covariance", "experimental", site * "_ICOS_HH_2023-2024.nc")
# Forcings: one value per half-hour window
forcings, times, ds_ec = site_forcing(file_ec; FT=FT) # local standard time
time_ec = ds_ec.time.val # start of the half-hours windows

## Simulation of 2023 (spin-up) and 2024 with adaptive implicit Euler
param = site_parameters(site)
u0 = [param.w_fc, param.w_fc, 0.0] # w_1, w_2, w_r: wet soil on 1 January, spin-up of 1 year
model = ProcessBasedModel{FT}(;
    forcings=forcings,
    parameters=param,
    t_span=times.span,
    u0=u0,
    saveat=times.bounds,
    thresholds=KavetskiSmoothing(),
)
initialize!(model)
# Two years of half-hourly `tstops` exceed the default `maxiters` of 1e5
@time EvaporationModel.solve!(
    model; alg=ImplicitEuler(; autodiff=AutoForwardDiff()), maxiters=10^7
)
println(model.sol.retcode, ", number of steps: ", length(model.sol.t), " saved, ",
    model.sol.stats.naccept, " accepted")
@assert SciMLBase.successful_retcode(model.sol) "Solver failed: $(model.sol.retcode)"

## Energy balance components of the model: means over the half-hours
λE_mod = model.diagnostics.λE_tot
H_mod = model.diagnostics.H_tot
# Available energy Rₙ - G, with Rₙ constant over each half-hour
A_mod = forcings.R_n .- model.diagnostics.G

## Observed energy balance components
# Uncorrected (F_MDS) and energy balance corrected (CORR) fluxes
obs(x) = FT.(collect(x)) # plain Vector, without the time lookup
λE_obs, H_obs = obs(ds_ec.Qle), obs(ds_ec.Qh)
λE_obs_cor, H_obs_cor = obs(ds_ec.Qle_cor), obs(ds_ec.Qh_cor)
λE_qc, H_qc = obs(ds_ec.Qle_qc), obs(ds_ec.Qh_qc)
A_obs = collect(ds_ec.Rnet .- ds_ec.Qg)
λE_obs_25, λE_obs_75 = collect(ds_ec.Qle_cor_25), collect(ds_ec.Qle_cor_75)
H_obs_25, H_obs_75 = collect(ds_ec.Qh_cor_25), collect(ds_ec.Qh_cor_75)

# Mean fluxes for June and December 2024 (compare with Table 1 of the paper)
for m in (6, 12)
    sel = (year.(time_ec) .== 2024) .& (month.(time_ec) .== m)
    println(
        monthname(m),
        ": Rn-G model = ", round(mean(A_mod[sel]); digits=1),
        ", Rn-G ICOS = ", round(mean(A_obs[sel]); digits=1),
        ", λE model = ", round(mean(λE_mod[sel]); digits=1),
        ", λE ICOS raw = ", round(mean(λE_obs[sel]); digits=1),
        ", λE ICOS corr = ", round(mean(λE_obs_cor[sel]); digits=1),
        ", H model = ", round(mean(H_mod[sel]); digits=1),
        ", H ICOS raw = ", round(mean(H_obs[sel]); digits=1),
        ", H ICOS corr = ", round(mean(H_obs_cor[sel]); digits=1)
    )
end
# Based on figures from this Table: paper uses energy balance corrected values.

## Figure 4 of Kronenberg et al. (2026)
figdir(args...) = plotsdir("2026kronenberg", args...)
mkpath(figdir())
cm = 37.8 # 1cm = 37.8 px
colors = (λE=:blue, H=:chocolate2, A=:seagreen) # as in the paper
periods = (
    (DateTime(2024, 6, 29), DateTime(2024, 7, 3)),
    (DateTime(2024, 12, 4), DateTime(2024, 12, 8)),
)

# Observed fluxes are energy balance corrected, consistent with their 25-75th quantile band.
# Time in days since the start of the period, since `band!` does not support `DateTime`.
# Issue: https://github.com/MakieOrg/Makie.jl/issues/3932
# `Makie.Axis` since ComponentArrays also exports an `Axis`
function period_axis(pos, period, ylim)
    n_days = Day(period[2] - period[1]).value
    ticks = (0:n_days, Dates.format.(period[1] .+ Day.(0:n_days), "u dd"))
    return Makie.Axis(
        pos; ylabel="Energy flux [W m⁻²]", xticks=ticks, limits=((0, n_days), ylim)
    )
end

function period_selection(period)
    sel = period[1] .<= time_ec .<= period[2]
    ms_per_day = Dates.value(Millisecond(Day(1)))
    ms_since_start = time_ec[sel] .- period[1] #
    ms_since_start = Dates.value.(ms_since_start)
    t = ms_since_start ./ ms_per_day
    return sel, t
end

"""
    row_ylim(period; step=50)

Y limits of a row: the minimum and maximum of all plotted data in `period`, rounded
outward to a multiple of `step` [W m⁻²].
"""
function row_ylim(period; step=50)
    sel, _ = period_selection(period)
    data = (
        λE_obs_25, λE_obs_75, H_obs_25, H_obs_75, λE_obs_cor, H_obs_cor, A_obs,
        λE_mod, H_mod, A_mod,
    )
    y_min = minimum(minimum(x[sel]) for x in data)
    y_max = maximum(maximum(x[sel]) for x in data)
    return (floor(y_min / step) * step, ceil(y_max / step) * step)
end
ylims = row_ylim.(periods)

function plot_observed!(pos, period, ylim)
    sel, t = period_selection(period)
    ax = period_axis(pos, period, ylim)
    band!(ax, t, λE_obs_25[sel], λE_obs_75[sel]; color=colors.λE, alpha=0.35,
        label="λE ICOS 25-75th quantile")
    band!(ax, t, H_obs_25[sel], H_obs_75[sel]; color=colors.H, alpha=0.35,
        label="H ICOS 25-75th quantile")
    lines!(ax, t, λE_obs_cor[sel]; color=colors.λE, label="λE ICOS corrected")
    lines!(ax, t, H_obs_cor[sel]; color=colors.H, label="H ICOS corrected")
    lines!(ax, t, A_obs[sel]; color=colors.A, label="Rn-G ICOS")
    return ax
end

function plot_modelled!(pos, period, ylim)
    sel, t = period_selection(period)
    ax = period_axis(pos, period, ylim)
    lines!(ax, t, λE_mod[sel]; color=colors.λE, label="λE model")
    lines!(ax, t, H_mod[sel]; color=colors.H, label="H model")
    lines!(ax, t, A_mod[sel]; color=colors.A, label="Rn-G model")
    return ax
end

# Layout: rows are the periods, columns observed (ICOS) and modelled, legends below.
# The right padding leaves room for the last x tick label, which is not a protrusion.
fig = Figure(; size=(30cm, 30cm), figure_padding=(10, 30, 10, 10))
axes = Matrix{Any}(undef, 2, 2)
for (i, (period, ylim)) in enumerate(zip(periods, ylims))
    axes[i, 1] = plot_observed!(fig[i, 1], period, ylim)
    axes[i, 2] = plot_modelled!(fig[i, 2], period, ylim)
    linkaxes!(axes[i, 1], axes[i, 2])
    hideydecorations!(axes[i, 2]; grid=false) # same y axis as the observed panel
end
axes[2, 1].xlabel = "date"
axes[2, 2].xlabel = "date"
for (j, title) in enumerate(["ICOS DE-Tha", "Model"])
    Label(fig[1, j, Top()], title; font=:bold, padding=(0, 0, 5, 0), valign=:bottom)
end
# Panel labels in the top left protrusion, so they do not add a row or column
for (label, (i, j)) in zip(["a)", "b)", "c)", "d)"], [(1, 1), (1, 2), (2, 1), (2, 2)])
    Label(fig[i, j, TopLeft()], label; font=:bold, padding=(0, 5, 5, 0), halign=:right)
end
# tellwidth=false: the legends do not widen their column
for j in 1:2
    Legend(fig[3, j], axes[2, j]; orientation=:horizontal, nbanks=2, tellwidth=false)
end
colgap!(fig.layout, 45) # room for the overhanging first and last x tick labels
save(figdir("2026_Kronenberg_fig4_reproduction.png"), fig)
fig

## Performance in 2024: NSE, KGE and modified KGE
# Only measured half-hours (QC = 0), as in Sect. 4.2 of the paper.
in_2024 = year.(time_ec) .== 2024
comparisons = (
    ("λE", λE_mod, λE_obs, λE_qc),
    ("λE corrected", λE_mod, λE_obs_cor, λE_qc),
    ("H", H_mod, H_obs, H_qc),
    ("H corrected", H_mod, H_obs_cor, H_qc),
)
function print_metrics(comparisons)
    for (name, ŷ, y, qc) in comparisons
        sel = in_2024 .& (qc .== 0)
        nse, kge, kge_mod = (metric(m, ŷ[sel], y[sel]) for m in (NSE(), KGE(), KGEmod()))
        println(
            rpad(name, 13), ": NSE = ", round(nse; digits=2), ", KGE = ",
            round(kge; digits=2), ", KGE' = ", round(kge_mod; digits=2),
            " (n = ", count(sel), ")",
        )
    end
end
print_metrics(comparisons)

##

## Calibration of `calibrated`: BFGS on the RMSE of λE in 2024 (no validation data)
# Loss on the energy balance corrected λE of the measured half-hours of 2024, with 2023 as
# spin-up, as for the metrics above. Gradients with ForwardDiff through the solver.
sel_cal = in_2024 .& (λE_qc .== 0)
# Fixed-step implicit Euler (Δt = 30/8 min): a smooth loss it seems, because unlike adaptive
# time stepping with ImplicitEuler (tried before) accepted/rejected steps
# change with the parameters. Newton settings as in scripts/calibration/synthetic.jl.
solver = ImplicitEuler(;
    nlsolve=NLNewton(; max_iter=50, always_new=true, relax=BackTracking(), check_div=false),
    autodiff=AutoForwardDiff(),
)
options = solve_kwargs(model; dt=1800 / 8, adaptive=false, maxiters=10^7)
calibrated = (:r_smin, :kB⁻¹) # names in `param` of the calibrated parameters
θ_default = [param[name] for name in calibrated]

# x = log(θ / θ_default): positive parameters of order 1 for BFGS, default values at x = 0
function calibration_parameters(x)
    p = similar(param, eltype(x)) .= param
    for (name, θ) in zip(calibrated, θ_default .* exp.(x))
        setproperty!(p, name, θ)
    end
    return p
end

# Window means of the accumulated diagnostics, nothing if the solve failed
function simulate(x)
    prob = with_parameters(model.prob, calibration_parameters(x))
    sol = solve(prob, solver; options...)
    return SciMLBase.successful_retcode(sol) ? window_means(sol, model.accumulated) : nothing
end

rmse_λE(means) = metric(RMSE(), means.λE_tot[sel_cal], λE_obs_cor[sel_cal])

function loss(x)
    means = simulate(x)
    return isnothing(means) ? convert(eltype(x), Inf) : rmse_λE(means)
end

x_default = zeros(length(calibrated))
@time result = optimize(
    loss,
    x_default,
    BFGS(),
    Optim.Options(; g_tol=1e-3, iterations=30, show_trace=true);
    autodiff=AutoForwardDiff(),
)
x_cal = Optim.minimizer(result)
println(result)
for (name, θ_0, θ) in zip(calibrated, θ_default, θ_default .* exp.(x_cal))
    println(name, ": ", round(θ_0; sigdigits=3), " → ", round(θ; sigdigits=3))
end

# Performance in 2024 of the default and calibrated parameters, with the calibration solver
for (label, x) in (("default", x_default), ("calibrated", x_cal))
    means = simulate(x)
    println(label, ": RMSE λE corrected = ", round(rmse_λE(means); digits=2), " W m⁻²")
    print_metrics((
        ("λE corrected", means.λE_tot, λE_obs_cor, λE_qc),
        ("H corrected", means.H_tot, H_obs_cor, H_qc),
    ))
end
