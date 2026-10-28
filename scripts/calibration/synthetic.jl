using DrWatson
@quickactivate "DifferentiableEvaporation"
using CairoMakie
using OrdinaryDiffEq
using OrdinaryDiffEqSDIRK
using ErrorMetrics: metric
using LineSearches: BackTracking
using OrdinaryDiffEq: SciMLBase
using OrdinaryDiffEqNonlinearSolve: NLNewton
using Random: Xoshiro


include(srcdir("site_data.jl"))

FT = Float64
site = "DE-Tha"
file_ec = datadir("exp_pro", "eddy_covariance", site * ".nc")

start_date = DateTime(2009, 01, 01, 00, 00)
end_date = DateTime(2010, 03, 31, 23, 30)
forcings, times, ds_ec = site_forcing(
    file_ec, start_date, end_date; FT=FT
)
param = site_parameters(site)
u0 = [param.w_fc, param.w_fc, 0]
model = ProcessBasedModel{FT}(;
    forcings=forcings,
    parameters=param,
    t_span=times.span,
    u0=u0,
    saveat=times.bounds,
    thresholds=KavetskiSmoothing(),
)
initialize!(model)
## Spin-up testing: How long until initial effect wears out
# Figure setup
fig = Figure()
ax1 = Makie.Axis(fig[1, 1]; ylabel="w₁ [-]")
ax2 = Makie.Axis(fig[2, 1]; ylabel="w₂ [-]")
linkaxes!(ax1, ax2)
hidexdecorations!(ax1)
time_plot = unix2datetime.(times.bounds)

fracs_fc = [0.05, 1/2, 1]
for frac_ in fracs_fc
    model.u0 = frac_ * [param.w_fc, param.w_fc, 0]
    initialize!(model)
    EvaporationModel.solve!(model; alg=Tsit5())

    lines!(ax1, time_plot, model.sol[1, :])
    lines!(ax2, time_plot, model.sol[2, :])
end
fig

# Conclusions that a year is sufficient as a spin up period
w1_end, w2_end = model.sol[1, end], model.sol[2, end]

## Synthetic truth generation
# Use end of 2010-03 value as initial value!
# 2010-04 as the period of interest
start_date = DateTime(2010, 04, 01, 00, 00)
end_date = DateTime(2010, 04, 30, 12, 30)
forcings, times, ds_ec = site_forcing(file_ec, start_date, end_date; FT=FT)
model.forcings = forcings
model.t_span = times.span
model.u0 = [w1_end, w2_end, 0]
model.saveat = times.bounds
initialize!(model)
# Note: use a high tolerance and a different solver than what we would
# typically use here.
EvaporationModel.solve!(model; alg=Kvaerno5(), abstol=1e-8, reltol=1e-8, yaxarray_output=true)
λE_syn_truth = model.output.λE_tot

## Synthetic observations: truth + constant Gaussian noise
p_true = copy(model.parameters)
λE_truth = copy(model.diagnostics.λE_tot)
σ_noise = 10.0 # W m⁻²
λE_noisy = λE_truth .+ σ_noise .* randn(Xoshiro(42), length(λE_truth))

## Loss landscapes (RMSE of λE_tot) in the (kB⁻¹, r_smin) plane
# Conceptually as Fig. 1 of Kavetski & Clark (2010)

# Solvers to compare
# Newton fix for fixed-step implicit Euler, see https://github.com/SciML/OrdinaryDiffEq.jl/issues/4905
OrdinaryDiffEq.OrdinaryDiffEqCore.isadaptive(::ImplicitEuler) = false
backward_euler = ImplicitEuler(;
    nlsolve=NLNewton(; max_iter=50, always_new=true, relax=BackTracking(), check_div=false),
    autodiff=AutoForwardDiff(),
)
solvers = [
    (name="Tsit5 (adaptive, tol = 10⁻⁶)", alg=Tsit5(), kwargs=(;)),
    (name="Kvaerno5 (adaptive, tol = 10⁻⁶)", alg=Kvaerno5(), kwargs=(;)),
    (
        name="Implicit Euler (Δt = 30/8 min)",
        alg=backward_euler,
        kwargs=(; dt=225.0, adaptive=false),
    ),
    (
        name="Explicit Euler (Δt = 30/8 min), clipping",
        alg=Euler(),
        kwargs=(; dt=225.0, adaptive=false, clamp_states=true),
    ),
]

# Parameter grid: 0.5 to 1.5 times the true values
n_grid = 41
kB⁻¹_values = range(0.5, 1.5; length=n_grid) .* p_true.kB⁻¹
r_smin_values = range(0.5, 1.5; length=n_grid) .* p_true.r_smin
n_runs = n_grid * n_grid
# Run k of the ensemble is grid point (i, j): kB⁻¹_values[i], r_smin_values[j]
grid_points = CartesianIndices((n_grid, n_grid))

# Ensemble of runs, one per grid point
function set_parameters(prob, ctx)
    i, j = Tuple(grid_points[ctx.sim_id])
    p = copy(p_true)
    p.kB⁻¹ = kB⁻¹_values[i]
    p.r_smin = r_smin_values[j]
    return with_parameters(prob, p)
end

# Keep only the RMSEs of a run (NaN if the solve failed), not the full solution
function compute_rmses(sol, ctx)
    rerun_run = false
    if !SciMLBase.successful_retcode(sol)
        return (NaN, NaN), rerun_run
    end
    λE_model = window_means(sol, model.accumulated).λE_tot
    rmse_truth = metric(RMSE(), λE_model, λE_truth)
    rmse_noisy = metric(RMSE(), λE_model, λE_noisy)
    return (rmse_truth, rmse_noisy), rerun_run
end

# safetycopy=false: `set_parameters` does not change `prob`, so no need to copy it (and the
# forcings) for every run; `with_parameters` gives every run its own current forcing window
ensemble = EnsembleProblem(
    model.prob; prob_func=set_parameters, output_func=compute_rmses, safetycopy=false
)

# Only runs the ensembles if `rerun` or if no saved results exist for this grid size, noise
# level and threshold treatment (set `rerun = true` after changing e.g. the solvers)
rerun = false
thresholds = nameof(typeof(model.thresholds))
landscape_data, _ = produce_or_load(
    (; n_grid, σ_noise, thresholds), datadir("sims", "calibration");
    prefix="loss_landscape_kB⁻¹_r_smin", force=rerun,
) do _
    solver_names = String[]
    rmse_truth = Matrix{Float64}[] # one n_grid × n_grid matrix per solver
    rmse_noisy = Matrix{Float64}[]
    for solver in solvers
        solve_options = solve_kwargs(model; solver.kwargs...)
        sim = solve(ensemble, solver.alg, EnsembleThreads(); trajectories=n_runs,
            solve_options...)
        rmses = reshape(sim.u, n_grid, n_grid) # (rmse_truth, rmse_noisy) per grid point
        push!(solver_names, solver.name)
        push!(rmse_truth, first.(rmses))
        push!(rmse_noisy, last.(rmses))
        @info "$(solver.name) done" n_failed = count(isnan, rmse_truth[end])
    end
    return @strdict(
        kB⁻¹_values, r_smin_values, solver_names, rmse_truth, rmse_noisy, λE_noisy, σ_noise
    )
end
@unpack solver_names, rmse_truth, rmse_noisy = landscape_data

"""
    plot_landscapes(rmse_per_solver, title; colorscale=identity)

2 × 2 heatmaps with contours of the RMSE landscapes of the 4 solvers on a shared colour
scale, with the true parameters and the grid minimum marked.
"""
function plot_landscapes(rmse_per_solver, title; colorscale=identity)
    fig = Figure(; size=(950, 850))
    Label(fig[0, 1:2], title; fontsize=18, font=:bold)

    all_values = reduce(vcat, vec.(rmse_per_solver))
    colorrange = extrema(filter(isfinite, all_values))

    panel_positions = [(1, 1), (1, 2), (2, 1), (2, 2)] # (row, column) of the 4 solvers
    for k in 1:4
        rmse = rmse_per_solver[k]
        row, col = panel_positions[k]
        ax = Makie.Axis(fig[row, col]; title=solver_names[k], xlabel="kB⁻¹ [-]",
            ylabel="r_smin [s m⁻¹]")
        heatmap!(ax, kB⁻¹_values, r_smin_values, rmse; colorrange, colorscale,
            colormap=:viridis)
        contour!(ax, kB⁻¹_values, r_smin_values, colorscale.(rmse); levels=15,
            color=:white, linewidth=0.6)

        scatter!(ax, [p_true.kB⁻¹], [p_true.r_smin]; marker=:star5, color=:red,
            markersize=15, label="true")
        i_min, j_min = Tuple(argmin(replace(rmse, NaN => Inf)))
        scatter!(ax, [kB⁻¹_values[i_min]], [r_smin_values[j_min]]; marker=:xcross,
            color=:orange, markersize=13, label="grid minimum")
    end
    axislegend(content(fig[1, 1]); position=:rt)
    Colorbar(fig[1:2, 3]; colorrange, colormap=:viridis, scale=colorscale,
        label="RMSE λE_tot [W m⁻²]")
    return fig
end

fig_truth = plot_landscapes(
    rmse_truth, "Loss-landscapes of the noise-free synthetic truth"; colorscale=log10
)
safesave(plotsdir("calibration", "loss_landscape_kB⁻¹_r_smin.png"), fig_truth)

fig_noisy = plot_landscapes(
    rmse_noisy, "Loss-landscapes of the synthetic truth + Gaussian noise (σ = $σ_noise W m⁻²)"
)
safesave(plotsdir("calibration", "loss_landscape_kB⁻¹_r_smin_noisy.png"), fig_noisy)

# Most important conclusion here for now is that Tsit5 is most noisy because changing
# a parameter changes a low how stepsize is chosen (which steps rejected or not) and this give
# the most noisy landscape. This is not the case for adaptive Kvaerno5. Suprisingy,
# fixed-step explicit Euler is smooth here.
