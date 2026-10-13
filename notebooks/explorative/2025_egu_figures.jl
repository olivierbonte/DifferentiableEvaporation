# Figures of the EGU poster: adaptive implicit Euler (IEₐ) versus explicit Euler (EE) for
# BE-Bra, 2010-03-15 to 2010-03-25, with the model from src/EvaporationModel.
using DrWatson
@quickactivate "DifferentiableEvaporation"
using ADTypes
using CairoMakie
using ComponentArrays
using Dates
using EvaporationModel
using OrdinaryDiffEq
using Statistics
using YAXArrays

include(projectdir("notebooks", "explorative", "real_forcings.jl")) # forcings_real, t_unix, ds_ec_sel

"""
    site_parameters(site; vegtype)

Model parameters of a site, with ISBA coefficients from the clay pedotransfer functions.
"""
function site_parameters(site; vegtype=EvergreenNeedleleafTrees(), FT=Float64)
    ds_soil = open_dataset(datadir("exp_pro", "soil", site * "_total_agg.nc"))
    soil = Dict(v => FT(ds_soil[v][1]) for v in propertynames(ds_soil))
    ds_ec = open_dataset(datadir("exp_pro", "eddy_covariance", site * ".nc"))
    veg = VegetationParameters(; vegtype=vegtype)
    clay = soil[:mean_clay_percentage]
    d_1 = FT(0.01) # normalization depth of the surface layer [m]
    return ComponentArray(;
        h=FT(ds_ec.canopy_height[1]),
        z_0ms=FT(0.01),
        w_sat=soil[:w_sat],
        a=compute_a(clay),
        p_soil=compute_p(clay),
        b=compute_b(Clay(), clay),
        C_1sat=c_1sat(NoilhanMahfouf96(), clay, d_1),
        C_2ref=c_2ref(clay),
        C_3=c_3(clay),
        d_1=d_1,
        d_2=soil[:root_depth],
        w_res=soil[:w_res],
        w_wp=soil[:w_wp],
        w_fc=soil[:w_fc],
        z_obs=FT(ds_ec.reference_height[1]),
        kB⁻¹=FT(log(10)),
        g_d=FT(veg.g_d),
        r_smin=FT(veg.r_smin),
        k_ext=FT(0.5),
    )
end

## Simulations
param = site_parameters(site)
u0 = [param.w_fc / 3, param.w_fc / 3, 1e-4] # w_1, w_2, w_r
function site_model()
    model = ProcessBasedModel{FT}(;
        forcings=forcings_real,
        parameters=param,
        t_span=(t_unix[1], t_unix[end]),
        u0=u0,
        saveat=t_unix,
        thresholds=KavetskiSmoothing()
    )
    initialize!(model)
    return model
end

# Adaptive implicit Euler, at the default tolerances of the model (solver_kwargs)
model_ie = site_model()
EvaporationModel.solve!(model_ie; alg=ImplicitEuler(; autodiff=AutoForwardDiff()))

# Explicit Euler at the time step of the forcing data, with clipping of the states
"""
    clipping_callback(model)

Clip the states to their physical range after every step, ``0 ≤ w_1, w_2 ≤ w_{sat}`` and
``0 ≤ w_r ≤ w_{r,max}``, as in classic fixed-step hydrological models (not mass conserving).
"""
function clipping_callback(model)
    affect! = let forcings = model.forcings
        function (integrator)
            u, p, t = integrator.u, integrator.p, integrator.t
            u[1] = clamp(u[1], zero(u[1]), p.w_sat)
            u[2] = clamp(u[2], zero(u[2]), p.w_sat)
            u[3] = clamp(u[3], zero(u[3]), max_canopy_capacity(forcings.LAI(t)))
            # `saveat` values are stored before the callbacks are applied
            sol = integrator.sol
            if !isempty(sol.t) && sol.t[end] == t
                sol.u[end] .= u
            end
            return nothing
        end
    end
    return DiscreteCallback(Returns(true), affect!; save_positions=(false, false))
end
dt = t_unix[2] - t_unix[1]
model_ee = site_model()
EvaporationModel.solve!(
    model_ee; alg=Euler(), dt=dt, adaptive=false, callback=clipping_callback(model_ee)
)
for model in (model_ie, model_ee)
    println(model.sol.alg, ": ", model.sol.retcode)
end

λE_ie = [d.λE_tot for d in model_ie.diagnostics.saveval]
λE_ee = [d.λE_tot for d in model_ee.diagnostics.saveval]
λE_obs = collect(ds_ec_sel.Qle_cor[:])
println("cor(λE IEₐ, λE observed) = ", cor(λE_ie, λE_obs))

## PLOTS FOR EGU 2025
figdir(args...) = projectdir("figures", args...)
mkpath(figdir())
cm = 37.8 #1cm = 37.8 px
time_plot = collect(unix2datetime.(t_unix)) # plain Vector: t_unix is a DimVector
precip_plot_mm_h = collect(forcings_real.P.(t_unix)) * 3600 #kg/(m²s) -> mm/h
state_labels = ["w₁ [-]", "w₂ [-]", "wᵣ [kg/m²]"]
state_ylims = (0, maximum(model_ie.sol[1, :]) * 1.5)
colors = Makie.wong_colors()

# States of a solution, with the precipitation as an inverted bar on a second y-axis.
# `Makie.Axis` since ComponentArrays also exports an `Axis`.
function plot_states!(pos, sol, title; legend=true)
    ax = Makie.Axis(pos; title=title, xlabel="Time", limits=(nothing, state_ylims))
    for (i, label) in enumerate(state_labels)
        lines!(ax, time_plot, sol[i, :]; label=label)
    end
    legend && axislegend(ax)
    ax_precip = Makie.Axis(
        pos;
        yaxisposition=:right,
        yreversed=true,
        ylabel="Precipitation [mm/h]",
        yticks=([0, 10, 20], ["0", "10", "20"]),
        limits=(nothing, (0, maximum(precip_plot_mm_h) * 2.5)),
        backgroundcolor=:transparent,
    )
    hidexdecorations!(ax_precip)
    hidespines!(ax_precip, :l, :t, :b)
    linkxaxes!(ax, ax_precip)
    barplot!(ax_precip, time_plot, precip_plot_mm_h; color=(:gray, 0.6), gap=0)
    return ax
end

function plot_fluxes!(pos)
    ax = Makie.Axis(pos; xlabel="Time", ylabel="λE [W/m²]")
    lines!(ax, time_plot, λE_obs; label="Observation", color=:black)
    lines!(ax, time_plot, λE_ee; label="EE", color=colors[1])
    lines!(ax, time_plot, λE_ie; label="IEₐ", color=colors[2])
    axislegend(ax)
    return ax
end

function plot_flux_difference!(pos)
    ax = Makie.Axis(pos; xlabel="Time", ylabel="λE(IEₐ) - λE(EE) [W/m²]")
    lines!(ax, time_plot, λE_ie - λE_ee; color=colors[3])
    return ax
end

# Figure 1: states for stable adapative implicit euler
fig_implicit = Figure()
plot_states!(fig_implicit[1, 1], model_ie.sol, "Adaptive Implicit Euler (IEₐ)")
save(figdir("IE_a_states.png"), fig_implicit)

# Figure 2: states for explicit euler
fig_explicit = Figure()
plot_states!(fig_explicit[1, 1], model_ee.sol, "Explicit Euler (EE)")
save(figdir("EE_states.png"), fig_explicit)

# Figure 3: fluxes for both methods
fig_implicit_fluxes = Figure()
plot_fluxes!(fig_implicit_fluxes[1, 1])
save(figdir("IE_a_fluxes.png"), fig_implicit_fluxes)

# Figure 4: difference in fluxes between implicit and explicit euler
fig_fluxes_diff = Figure()
plot_flux_difference!(fig_fluxes_diff[1, 1])
save(figdir("IEa_diff_EE_fluxes.png"), fig_fluxes_diff)

# Combine the plots in one
title_fontsize = 18
tick_fontsize = title_fontsize - 4
function combined_figure(; background=:white)
    theme = Theme(;
        backgroundcolor=background,
        Axis=(
            backgroundcolor=background,
            titlesize=title_fontsize,
            xlabelsize=title_fontsize - 4,
            ylabelsize=title_fontsize - 4,
            xticklabelsize=tick_fontsize,
            yticklabelsize=tick_fontsize,
        ),
        Legend=(labelsize=12,),
        Lines=(linewidth=2.5,),
    )
    return with_theme(theme) do
        fig = Figure(; size=(30cm, 29cm))
        ax_ie = plot_states!(fig[1, 1], model_ie.sol, "Adaptive Implicit Euler (IEₐ)")
        ax_ee = plot_states!(fig[1, 2], model_ee.sol, "Explicit Euler (EE)"; legend=false)
        ax_fluxes = plot_fluxes!(fig[2, 1])
        ax_diff = plot_flux_difference!(fig[2, 2])
        hidexdecorations!(ax_ie; grid=false)
        hidexdecorations!(ax_ee; grid=false)
        linkxaxes!(ax_ie, ax_ee, ax_fluxes, ax_diff)
        return fig
    end
end
fig_combined = combined_figure()
save(figdir("combined.png"), fig_combined)
save(figdir("combined.svg"), fig_combined)
save(figdir("combined_transparent.svg"), combined_figure(; background=:transparent))
