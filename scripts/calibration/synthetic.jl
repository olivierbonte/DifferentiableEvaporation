using DrWatson
@quickactivate "DifferentiableEvaporation"
using CairoMakie
using DifferentiationInterface: value_and_derivative
using OrdinaryDiffEq
using OrdinaryDiffEqSDIRK
using ErrorMetrics: metric
using OrdinaryDiffEq: SciMLBase
using Random: Xoshiro
seed = 42 # of the noise of the synthetic observations
truth_tol = 1e-10


include(srcdir("site_data.jl"))

FT = Float64
site = "DE-Tha"
file_ec = datadir("exp_pro", "eddy_covariance", site * ".nc")

start_date = DateTime(2009, 01, 01, 00, 00)
end_date = DateTime(2010, 06, 30, 23, 30)
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
# Use end of 2010-06 value as initial value!
# 2010-07 as the period of interest: dry (w₁ near w_wp ~40% of the time), with ~100 mm of rain
start_date = DateTime(2010, 07, 01, 00, 00)
end_date = DateTime(2010, 07, 31, 23, 30)
forcings, times, _ = site_forcing(file_ec, start_date, end_date; FT=FT)
model.forcings = forcings
model.t_span = times.span
model.u0 = [w1_end, w2_end, 0]
model.saveat = times.bounds
initialize!(model)
# Note: use a high tolerance and a different solver than what we would
# typically use here.
# Recommendation from https://docs.sciml.ai/DiffEqDocs/stable/api/ordinarydiffeq/implicit/SDIRK/#Low-tolerance-(1e-8)
EvaporationModel.solve!(model; alg=Kvaerno5(), abstol=truth_tol, reltol=truth_tol)

p_true = copy(model.parameters)
period = Dates.format(start_date, "yyyy-mm") # 1 month



## Synthetic observations: truth + constant Gaussian noise
"""
    synthetic_target(truth_sol, name, label, unit, observe_function, σ_noise)

Calibration target `name`: `observe_function(sol)` gives the observed quantity of a solution
, here of the synthetic truth `truth_sol`, with Gaussian noise of standard deviation
`σ_noise` [`unit`] added for the noisy observations.
"""
function synthetic_target(truth_sol, name, label, unit, observe_function, σ_noise)
    truth = observe_function(truth_sol)
    noisy = truth .+ σ_noise .* randn(Xoshiro(seed), length(truth))
    return (; name, label, unit, observe_function, truth, noisy, σ_noise)
end
observe_λE_tot(sol) = window_means(sol, model.accumulated).λE_tot
observe_w_1(sol) = sol[1, :] # states at `saveat`
targets = (
    λE_tot=synthetic_target(model.sol, :λE_tot, "λE_tot", "W m⁻²", observe_λE_tot, 10.0),
    w_1=synthetic_target(model.sol, :w_1, "w₁", "m³ m⁻³", observe_w_1, 0.02),
)

"""
    target_rmses(target, sol)

RMSEs of the solution `sol` against the noise-free truth and the noisy observations of
`target`.
"""
function target_rmses(target, sol)
    simulated = target.observe_function(sol)
    rmse_vs_truth = metric(RMSE(), simulated, target.truth)
    rmse_vs_noisy = metric(RMSE(), simulated, target.noisy)
    return [rmse_vs_truth, rmse_vs_noisy]
end

## Loss landscapes: RMSE of λE_tot in the (kB⁻¹, r_smin) plane, RMSE of w₁ in the
## (b, C_2ref) plane
# Conceptually as Fig. 1 of Kavetski & Clark (2010)

# Solvers to compare
# Newton fix for fixed-step implicit Euler, see https://github.com/SciML/OrdinaryDiffEq.jl/issues/4905
OrdinaryDiffEq.OrdinaryDiffEqCore.isadaptive(::ImplicitEuler) = false
backward_euler = ImplicitEuler(;
    nlsolve=fixed_step_newton(; max_iter=50),
    autodiff=AutoForwardDiff(),
)
# Adaptive solvers at a tight and a loose tolerance (microscale noise from changes
# in the step-size sequence, Fig. 3 of Kavetski & Clark (2010)), fixed-step solvers with one
# step per forcing window (a bias of the optimum that grows with Δt). Explicit Euler without
# clamping is unstable at Δt = 30 min
solvers = [
    (name="Tsit5 (adaptive, tol = 10⁻⁶)", alg=Tsit5(), kwargs=(;)),
    (
        name="Tsit5 (adaptive, tol = 10⁻³)",
        alg=Tsit5(),
        kwargs=(; abstol=1e-3, reltol=1e-3),
    ),
    (
        name="KenCarp4 (adaptive, tol = 10⁻³)",
        alg=KenCarp4(),
        kwargs=(; abstol=1e-3, reltol=1e-3),
    ),
    (
        name="Implicit Euler (Δt = 30 min)",
        alg=backward_euler,
        kwargs=(; dt=1800.0, adaptive=false),
    ),
    # Stable at Δt = 30 min by clamping the states to their physical range after every step,
    # as in classic hydrological models.
    (
        name="Explicit Euler (Δt = 30 min), clamping",
        alg=Euler(),
        kwargs=(; dt=1800.0, adaptive=false, clamp_states=true),
    ),
]

# Parameter planes: the original (kB⁻¹, r_smin) and the soil parameters of the fast surface
# layer (b, C_1sat), whose stiffness grows when w₁ is dry
planes = [
    (x=:kB⁻¹, y=:r_smin, xlabel="kB⁻¹ [-]", ylabel="r_smin [s m⁻¹]"),
    (x=:b, y=:C_1sat, xlabel="b [-]", ylabel="C₁sat [-]"),
]

# Parameter grid: 0.5 to 1.5 times the true values
n_grid = 41
n_runs = n_grid * n_grid
# Run k of the ensemble is grid point (i, j): x_values[i], y_values[j]
grid_points = CartesianIndices((n_grid, n_grid))
grid_values(name) = range(0.5, 1.5; length=n_grid) .* getproperty(p_true, name)

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

# Only runs the ensembles if `rerun` or if no saved results exist for this plane, period, grid
# size, noise level and threshold treatment (set `rerun = true` after changing e.g. the solvers)
rerun = true
thresholds = nameof(typeof(model.thresholds))
period = Dates.format(start_date, "yyyy-mm")

"""
    compute_landscapes(plane)

RMSE landscapes of all `solvers` on the `n_grid × n_grid` grid of the parameters `plane.x`
and `plane.y`, loaded from `datadir("sims", "calibration")` if already computed.
"""
function compute_landscapes(plane)
    x_values, y_values = grid_values(plane.x), grid_values(plane.y)
    # Ensemble of runs, one per grid point
    function set_parameters(prob, ctx)
        i, j = Tuple(grid_points[ctx.sim_id])
        p = copy(p_true)
        setproperty!(p, plane.x, x_values[i])
        setproperty!(p, plane.y, y_values[j])
        return with_parameters(prob, p)
    end
    # safetycopy=false: `set_parameters` does not change `prob`, so no need to copy it (and the
    # forcings) for every run; `with_parameters` gives every run its own current forcing window
    ensemble = EnsembleProblem(
        model.prob; prob_func=set_parameters, output_func=compute_rmses, safetycopy=false
    )
    landscape_data, _ = produce_or_load(
        (; n_grid, σ_noise, thresholds, period), datadir("sims", "calibration");
        prefix="loss_landscape_$(plane.x)_$(plane.y)", force=rerun,
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
            @info "$(plane.x), $(plane.y): $(solver.name) done" n_failed = count(
                isnan, rmse_truth[end]
            )
        end
        return @strdict(
            x_values, y_values, solver_names, rmse_truth, rmse_noisy, λE_noisy, σ_noise
        )
    end
    return landscape_data
end

"""
    plot_landscapes(plane, landscape_data, rmse_per_solver, title; colorscale=identity)

3 × 2 heatmaps with contours of the RMSE landscapes of the 6 solvers in the parameter `plane`
on a shared colour scale, with the true parameters and the grid minimum marked.
"""
function plot_landscapes(plane, landscape_data, rmse_per_solver, title; colorscale=identity)
    @unpack x_values, y_values, solver_names = landscape_data
    fig = Figure(; size=(950, 1250))
    Label(fig[0, 1:2], title; fontsize=18, font=:bold)

    all_values = reduce(vcat, vec.(rmse_per_solver))
    colorrange = extrema(filter(isfinite, all_values))

    # (row, column) of the 6 solvers: one row per solver pair of `solvers`
    panel_positions = [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1), (3, 2)]
    for k in 1:6
        rmse = rmse_per_solver[k]
        row, col = panel_positions[k]
        ax = Makie.Axis(fig[row, col]; title=solver_names[k], xlabel=plane.xlabel,
            ylabel=plane.ylabel)
        heatmap!(ax, x_values, y_values, rmse; colorrange, colorscale, colormap=:viridis)
        contour!(ax, x_values, y_values, colorscale.(rmse); levels=15, color=:white,
            linewidth=0.6)

        scatter!(ax, [getproperty(p_true, plane.x)], [getproperty(p_true, plane.y)];
            marker=:star5, color=:red, markersize=15, label="true")
        i_min, j_min = Tuple(argmin(replace(rmse, NaN => Inf)))
        scatter!(ax, [x_values[i_min]], [y_values[j_min]]; marker=:xcross,
            color=:orange, markersize=13, label="grid minimum")
    end
    axislegend(content(fig[1, 1]); position=:rt)
    Colorbar(fig[1:3, 3]; colorrange, colormap=:viridis, scale=colorscale,
        label="RMSE λE_tot [W m⁻²]")
    return fig
end

for plane in planes
    landscape_data = compute_landscapes(plane)
    @unpack rmse_truth, rmse_noisy = landscape_data
    name = "loss_landscape_$(plane.x)_$(plane.y)_$(period)"

    fig_truth = plot_landscapes(plane, landscape_data, rmse_truth,
        "Loss-landscapes of the noise-free synthetic truth"; colorscale=log10)
    safesave(plotsdir("calibration", "$(name).png"), fig_truth)

    fig_noisy = plot_landscapes(plane, landscape_data, rmse_noisy,
        "Loss-landscapes of the synthetic truth + Gaussian noise (σ = $σ_noise W m⁻²)")
    safesave(plotsdir("calibration", "$(name)_noisy.png"), fig_noisy)
end
