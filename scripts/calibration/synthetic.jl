# Run multithreaded, e.g. `julia -t 8 scripts/calibration/synthetic.jl`:
# ` the landscape ensembles (EnsembleThreads) and the slice points
# (Threads.@threads) run in parallel over the threads, sequentially with one thread
using DrWatson
@quickactivate "DifferentiableEvaporation"
using CairoMakie
using DifferentiationInterface: AutoFiniteDiff, derivative, value_and_derivative
using OrdinaryDiffEq
using OrdinaryDiffEqSDIRK
using ErrorMetrics: metric
using OrdinaryDiffEq: SciMLBase
using Random: Xoshiro
FT = Float64
seed = 42 # of the noise of the synthetic observations
truth_tol = 1e-10
forced_rerun = false


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

# Runs that leave the physical range of the states are considered diverged
solve_options(solver) = solve_kwargs(model; unstable_check=unphysical_states, solver.kwargs...)

"""
    parameters_with(changes...)

Copy of `p_true` with the parameters in `changes` (`name => value` pairs) set. The
`ComponentArray` constructor promotes the number type to that of the new values (e.g.
ForwardDiff dual numbers).
"""
parameters_with(changes...) = ComponentArray(merge(NamedTuple(p_true), (; changes...)))

# Parameter planes with their calibration target. The soil parameters of the fast surface
# layer are calibrated on w₁ itself: they only act on λE_tot through the small soil
# evaporation. b sets the force term of w₁, C_2ref the restore term
planes = [(x=:kB⁻¹, y=:r_smin, target=targets.λE_tot), (x=:b, y=:C_2ref, target=targets.w_1)]
param_labels = (kB⁻¹="kB⁻¹ [-]", r_smin="r_smin [s m⁻¹]", b="b [-]", C_2ref="C₂ref [-]")

# Parameter grids: (1 ± rel_range) times the true values
rel_range = 0.5
n_grid = 41
# Run k of the ensemble is grid point (i, j): x_values[i], y_values[j]
grid_points = CartesianIndices((n_grid, n_grid))
function grid_values(name, n=n_grid)
    return range(1 - rel_range, 1 + rel_range; length=n) .* getproperty(p_true, name)
end

threshold_treatment = nameof(typeof(model.thresholds))

"""
    compute_landscapes(plane)

RMSE landscapes of `plane.target` of all `solvers` on the `n_grid × n_grid` grid of the
parameters `plane.x` and `plane.y`, loaded from `datadir("sims", "calibration", "synthetic")`
if already computed.
"""
function compute_landscapes(plane)
    (; target) = plane
    (; σ_noise) = target
    x_values, y_values = grid_values(plane.x), grid_values(plane.y)
    # Keep only the RMSEs of a run (NaN if the solve failed or diverged), not the full solution
    function compute_rmses(sol, ctx)
        if !SciMLBase.successful_retcode(sol)
            return [NaN, NaN], false
        end
        return target_rmses(target, sol), false
    end
    # Ensemble of runs, one per grid point
    function set_parameters(prob, ctx)
        i, j = Tuple(grid_points[ctx.sim_id])
        p = parameters_with(plane.x => x_values[i], plane.y => y_values[j])
        return with_parameters(prob, p)
    end
    # safetycopy=false: `set_parameters` does not change `prob`, so no need to copy it (and the
    # forcings) for every run; `with_parameters` gives every run its own current forcing window
    ensemble_problem = EnsembleProblem(
        model.prob; prob_func=set_parameters, output_func=compute_rmses, safetycopy=false
    )
    landscape_data, _ = produce_or_load(
        (; target=target.name, n_grid, rel_range, σ_noise, thresholds=threshold_treatment,
            period),
        datadir("sims", "calibration", "synthetic");
        prefix="loss_landscape_$(plane.x)_$(plane.y)", force=forced_rerun,
    ) do _
        solver_names = String[]
        rmse_truth = Matrix{FT}[] # one n_grid × n_grid matrix per solver
        rmse_noisy = Matrix{FT}[]
        for solver in solvers
            sim = solve(ensemble_problem, solver.alg, EnsembleThreads();
                trajectories=length(grid_points), solve_options(solver)...)
            rmses = reshape(sim.u, n_grid, n_grid) # [vs truth, vs noisy] per grid point
            solver_rmse_truth = first.(rmses)
            push!(solver_names, solver.name)
            push!(rmse_truth, solver_rmse_truth)
            push!(rmse_noisy, last.(rmses))
            n_failed = count(isnan, solver_rmse_truth)
            @info "$(plane.x), $(plane.y): $(solver.name) done" n_failed
        end
        return @strdict(x_values, y_values, solver_names, rmse_truth, rmse_noisy, σ_noise)
    end
    return landscape_data
end

"""
    plot_landscapes(plane, landscape_data, rmse_per_solver, title; colorscale=identity)

2 × 3 heatmaps with contours of the RMSE landscapes of the 5 solvers in the parameter `plane`
on a shared colour scale, with the true parameters and the grid minimum marked. The colour
bar and legend go in the empty sixth panel.
"""
function plot_landscapes(plane, landscape_data, rmse_per_solver, title; colorscale=identity)
    @unpack x_values, y_values, solver_names = landscape_data
    (; label, unit) = plane.target
    fig = Figure(; size=(1100, 680))
    Label(fig[0, 1:3], title; fontsize=18, font=:bold)

    all_values = reduce(vcat, vec.(rmse_per_solver))
    colorrange = extrema(filter(isfinite, all_values))

    # (row, column) of the 5 solvers, (2, 3) is left for the colour bar and legend
    panel_positions = [(1, 1), (1, 2), (1, 3), (2, 1), (2, 2)]
    axes = Makie.Axis[]
    for k in 1:5
        rmse = rmse_per_solver[k]
        row, col = panel_positions[k]
        ax = Makie.Axis(fig[row, col]; title=solver_names[k], xlabel=param_labels[plane.x],
            ylabel=param_labels[plane.y])
        push!(axes, ax)
        heatmap!(ax, x_values, y_values, rmse; colorrange, colorscale, colormap=:viridis)
        contour!(ax, x_values, y_values, colorscale.(rmse); levels=15, color=:white,
            linewidth=0.6)

        scatter!(ax, [getproperty(p_true, plane.x)], [getproperty(p_true, plane.y)];
            marker=:star5, color=:red, markersize=15, label="true")
        i_min, j_min = Tuple(argmin(replace(rmse, NaN => Inf))) # failed runs are NaN
        scatter!(ax, [x_values[i_min]], [y_values[j_min]]; marker=:xcross,
            color=:orange, markersize=13, label="grid minimum")
    end
    # Shared axes: tick labels only in the left column and below the bottom panels
    linkaxes!(axes...)
    hidexdecorations!.(axes[1:2])
    hideydecorations!.(axes[[2, 3, 5]])

    # tellwidth/tellheight=false so the colour bar and legend don't set the size of their panel
    Colorbar(fig[2, 3][1, 1]; colorrange, colormap=:viridis, scale=colorscale,
        label="RMSE $label [$unit]", vertical=false, flipaxis=false, tellwidth=false,
        tellheight=false, valign=:bottom)
    Legend(fig[2, 3][2, 1], axes[1]; framevisible=false, orientation=:horizontal,
        tellwidth=false, tellheight=false, valign=:top)
    return fig
end

for plane in planes
    landscape_data = compute_landscapes(plane)
    @unpack rmse_truth, rmse_noisy = landscape_data
    (; name, label, unit, σ_noise) = plane.target
    figure_name = "loss_landscape_$(plane.x)_$(plane.y)_$(name)_$(period)"
    # The noise-free RMSE goes to zero at the truth: log colour scale
    variants = [
        (suffix="", rmse=rmse_truth, colorscale=log10,
            description="the noise-free synthetic truth"),
        (suffix="_noisy", rmse=rmse_noisy, colorscale=identity,
            description="the synthetic truth + Gaussian noise (σ = $σ_noise $unit)"),
    ]
    for (; suffix, rmse, colorscale, description) in variants
        fig = plot_landscapes(plane, landscape_data, rmse,
            "Loss-landscapes ($label) of $description"; colorscale)
        safesave(plotsdir("calibration", "$(figure_name)$(suffix).png"), fig)
    end
end

## Loss and its forward mode AD derivative along slices: b and C_2ref on w₁, r_smin on λE_tot
# Compared with central finite differences of the loss, which the microscale noise of the
# step-size sequence makes noisy (Fig. 3 of Kavetski & Clark, 2010).
slices = [
    (param=:b, target=targets.w_1),
    (param=:C_2ref, target=targets.w_1),
    # On λE_tot: r_smin controls the transpiration, the largest part of λE_tot at this forest site
    (param=:r_smin, target=targets.λE_tot),
]
n_slice = 201

"""
    slice_rmses(slice, solver, param_value)

RMSEs of the model run with `solver` and `slice.param` set to `param_value` against the
noise-free truth and the noisy observations of `slice.target` (NaN, also for the derivative,
if the solve failed or diverged).
"""
function slice_rmses(slice, solver, param_value)
    p = parameters_with(slice.param => param_value)
    prob = with_parameters(model.prob, p)
    sol = solve(prob, solver.alg; solve_options(solver)...)
    if !SciMLBase.successful_retcode(sol)
        # param_value * NaN: NaN partials too (with NaN of type param_value they would be 0)
        return fill(param_value * NaN, 2)
    end
    return target_rmses(slice.target, sol)
end

"""
    compute_slice(slice)

RMSEs and their derivatives (forward mode AD and central finite differences, with the default
step of FiniteDiff) along `n_slice` values of `slice.param` for all `solvers`, against the
noise-free truth and the noisy observations, loaded from
`datadir("sims", "calibration", "synthetic")` if already computed. The derivatives of failed
runs are NaN.
"""
function compute_slice(slice)
    (; target) = slice
    (; σ_noise) = target
    slice_data, _ = produce_or_load(
        (; target=target.name, n_slice, rel_range, σ_noise, thresholds=threshold_treatment,
            period),
        datadir("sims", "calibration", "synthetic");
        prefix="loss_slice_$(slice.param)", force=forced_rerun,
    ) do _
        param_values = collect(grid_values(slice.param, n_slice))
        solver_names = [solver.name for solver in solvers]
        # One vector along the slice per solver
        rmse_truth, rmse_noisy = Vector{FT}[], Vector{FT}[]
        derivative_truth, derivative_noisy = Vector{FT}[], Vector{FT}[]
        derivative_fd_truth, derivative_fd_noisy = Vector{FT}[], Vector{FT}[]
        for solver in solvers
            rmses = Matrix{FT}(undef, 2, n_slice) # [vs truth, vs noisy] per slice point
            derivatives = Matrix{FT}(undef, 2, n_slice)
            derivatives_fd = Matrix{FT}(undef, 2, n_slice)
            Threads.@threads for i in 1:n_slice
                rmses_at(param_value) = slice_rmses(slice, solver, param_value)
                rmses[:, i], derivatives[:, i] = value_and_derivative(
                    rmses_at, AutoForwardDiff(), param_values[i]
                )
                derivatives_fd[:, i] = derivative(
                    rmses_at, AutoFiniteDiff(; fdtype=Val(:central)), param_values[i]
                )
            end
            solver_rmse_truth = rmses[1, :]
            push!(rmse_truth, solver_rmse_truth)
            push!(rmse_noisy, rmses[2, :])
            push!(derivative_truth, derivatives[1, :])
            push!(derivative_noisy, derivatives[2, :])
            push!(derivative_fd_truth, derivatives_fd[1, :])
            push!(derivative_fd_noisy, derivatives_fd[2, :])
            n_failed = count(isnan, solver_rmse_truth)
            @info "Slice $(slice.param): $(solver.name) done" n_failed
        end
        return @strdict(
            param_values, solver_names, rmse_truth, rmse_noisy, derivative_truth,
            derivative_noisy, derivative_fd_truth, derivative_fd_noisy, σ_noise
        )
    end
    return slice_data
end

"""
    derivative_ylims(derivatives, reference; max_factor=2)

y-limits covering all finite `derivatives` (vectors) up to `max_factor` times the largest
absolute value of the `reference` derivative: leaves out the exploding derivatives of a
failing solve. Pass only the AD derivatives, so the noisy finite differences (e.g. of the
loose tolerance adaptive solvers) don't hide them: finite differences outside are clipped.
"""
function derivative_ylims(derivatives, reference; max_factor=2)
    bound = max_factor * maximum(abs, filter(isfinite, reference))
    in_bounds = filter(d -> isfinite(d) && abs(d) <= bound, reduce(vcat, derivatives))
    low, high = extrema(in_bounds)
    margin = 0.1 * (high - low)
    return low - margin, high + margin
end

"""
    plot_slice(slice, slice_data, rmse_per_solver, derivative_per_solver,
        derivative_fd_per_solver, title)

3 × 2 panels with the RMSE (left y-axis) and its derivative (right y-axis) along the slice of
the 5 solvers: forward mode AD (line) and central finite differences (dots). All panels share
the x-axis, the RMSE and the derivative axes each share their range across the solvers. The
legend goes in the empty sixth panel.
"""
function plot_slice(
    slice, slice_data, rmse_per_solver, derivative_per_solver, derivative_fd_per_solver, title
)
    @unpack param_values, solver_names = slice_data
    (; unit) = slice.target
    fig = Figure(; size=(1000, 830))
    Label(fig[0, 1:2], title; fontsize=18, font=:bold)
    xlabel = param_labels[slice.param]
    color_rmse, color_ad, color_fd = :black, Makie.wong_colors()[1], :orange
    # Limits from the AD derivatives only, the first solver is the reference
    derivative_axis_limits = derivative_ylims(derivative_per_solver, derivative_per_solver[1])
    param_true = getproperty(p_true, slice.param)

    # (row, column) of the 5 solvers, (3, 2) is left for the legend
    panel_positions = [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1)]
    axes_rmse, axes_derivative = Makie.Axis[], Makie.Axis[]
    for k in 1:5
        row, col = panel_positions[k]
        # Axis labels only on the outer sides: RMSE left in the left column, derivative
        # right in the right column (tick labels everywhere, coloured per axis)
        ax_rmse = Makie.Axis(fig[row, col]; title=solver_names[k], xlabel,
            ylabel=col == 1 ? "RMSE [$unit]" : "")
        ax_derivative = Makie.Axis(fig[row, col]; yaxisposition=:right,
            ylabel=col == 2 ? "∂RMSE/∂$(slice.param)" : "", ylabelcolor=color_ad,
            yticklabelcolor=color_ad, ytickcolor=color_ad)
        hidespines!(ax_derivative)
        hidexdecorations!(ax_derivative)
        push!(axes_rmse, ax_rmse)
        push!(axes_derivative, ax_derivative)

        hlines!(ax_derivative, 0; color=color_ad, linestyle=:dot, linewidth=0.8)
        vlines!(ax_rmse, param_true; color=:red, linestyle=:dash, linewidth=0.8)
        scatter!(ax_derivative, param_values, derivative_fd_per_solver[k]; color=color_fd,
            markersize=4)
        lines!(ax_derivative, param_values, derivative_per_solver[k]; color=color_ad)
        lines!(ax_rmse, param_values, rmse_per_solver[k]; color=color_rmse, linewidth=2)
    end
    # x tick labels only at the bottom of each column
    hidexdecorations!.(axes_rmse[1:3]; grid=false)
    linkxaxes!(axes_rmse..., axes_derivative...)
    linkyaxes!(axes_rmse...)
    linkyaxes!(axes_derivative...)
    ylims!(axes_derivative[1], derivative_axis_limits)

    legend_entries = [
        LineElement(; color=color_rmse, linewidth=2),
        LineElement(; color=color_ad),
        MarkerElement(; color=color_fd, marker=:circle, markersize=6),
        LineElement(; color=:red, linestyle=:dash),
    ]
    legend_labels = [
        "RMSE (left axis)", "∂RMSE/∂$(slice.param), forward mode AD (right axis)",
        "∂RMSE/∂$(slice.param), finite differences (right axis)", "true value",
    ]
    # tellwidth/tellheight=false so the legend doesn't set the size of its panel
    Legend(fig[3, 2], legend_entries, legend_labels; framevisible=false, tellwidth=false,
        tellheight=false)
    return fig
end

for slice in slices
    slice_data = compute_slice(slice)
    @unpack rmse_truth, rmse_noisy, derivative_truth, derivative_noisy = slice_data
    @unpack derivative_fd_truth, derivative_fd_noisy = slice_data
    (; name, label, unit, σ_noise) = slice.target
    figure_name = "loss_slice_$(slice.param)_$(name)_$(period)"
    variants = [
        (suffix="", rmse=rmse_truth, derivative_ad=derivative_truth,
            derivative_fd=derivative_fd_truth, description="noise-free synthetic truth"),
        (suffix="_noisy", rmse=rmse_noisy, derivative_ad=derivative_noisy,
            derivative_fd=derivative_fd_noisy,
            description="truth + noise (σ = $σ_noise $unit)"),
    ]
    for (; suffix, rmse, derivative_ad, derivative_fd, description) in variants
        fig = plot_slice(slice, slice_data, rmse, derivative_ad, derivative_fd,
            "RMSE ($label) and its derivative along $(slice.param), $description")
        safesave(plotsdir("calibration", "$(figure_name)$(suffix).png"), fig)
    end
end
