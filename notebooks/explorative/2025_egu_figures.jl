# Figures of the EGU poster: adaptive implicit Euler (IEₐ) versus explicit Euler (EE) for
# BE-Bra, 2010-03-15 to 2010-03-25, with the model from src/EvaporationModel.
using DrWatson
@quickactivate "DifferentiableEvaporation"
using ADTypes
using ComponentArrays
using Dates
using EvaporationModel
using OrdinaryDiffEq
using Plots
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
gr()
figdir(args...) = projectdir("figures", args...)
mkpath(figdir())
cm = 37.8 #1cm = 37.8 px
time_plot = unix2datetime.(t_unix)
date_ticks = [unix2datetime(t_unix[50]), unix2datetime(t_unix[end-50])]
xticks = (date_ticks, Dates.format.(date_ticks, "yyyy-mm-dd"))
precip_plot_mm_h = forcings_real.P.(t_unix) * 3600 #kg/(m²s) -> mm/h
y_ticks_precip = ([0, 10, 20], ["0", "10", "20"])
state_labels = ["w₁ [-]" "w₂ [-]" "wᵣ [kg/m²]"]

function add_precipitation!(fig)
    plot!(
        twinx(fig),
        time_plot,
        precip_plot_mm_h;
        fill=(0, :gray),
        color=:gray,
        yflip=true,
        ylabel="Precipitation [mm/h]",
        legend=:none,
        xticks=xticks,
        yticks=y_ticks_precip,
        ylims=(0, maximum(precip_plot_mm_h) * 2.5),
    )
    return fig
end

# Figure 1: states for stable adapative implicit euler
fig_implicit = plot(
    time_plot,
    Array(model_ie.sol)';
    label=state_labels,
    xlabel="Time",
    ylims=(0, maximum(model_ie.sol[1, :]) * 1.5),
    title="Adaptive Implicit Euler (IEₐ)",
    xticks=xticks,
)
add_precipitation!(fig_implicit)
savefig(fig_implicit, figdir("IE_a_states.png"))

# Figure 2: states for explicit euler
fig_explicit = plot(
    time_plot,
    Array(model_ee.sol)';
    label=state_labels,
    xlabel="Time",
    ylims=ylims(fig_implicit),
    title="Explicit Euler (EE)",
)
add_precipitation!(fig_explicit)
savefig(fig_explicit, figdir("EE_states.png"))

# Figure 3: fluxes for both methods
fig_implicit_fluxes = plot(
    time_plot,
    λE_obs;
    label="Observation",
    color=:black,
    ylabel="λE [W/m²]",
    xlabel="Time",
    xticks=xticks,
    framestyle=:box,
)
plot!(time_plot, λE_ee; label="EE", color=palette(:default)[1])
plot!(time_plot, λE_ie; label="IEₐ", color=palette(:default)[2])
savefig(fig_implicit_fluxes, figdir("IE_a_fluxes.png"))

# Figure 4: difference in fluxes between implicit and explicit euler
fig_fluxes_diff = plot(
    time_plot,
    λE_ie - λE_ee;
    label=:none,
    ylabel="λE(IEₐ) - λE(EE) [W/m²]",
    xlabel="Time",
    xticks=xticks,
    framestyle=:box,
    color=palette(:default)[3],
)
savefig(fig_fluxes_diff, figdir("IEa_diff_EE_fluxes.png"))

# Combine the plots in one
l = @layout [_ a{0.485w} b{0.485w} _; _ c{0.485w} d{0.485w} _]
title_fontsize = 18
tick_fontsize = title_fontsize - 4
fig_implicit_combine = plot(
    fig_implicit;
    xlabel="",
    xtickfontcolor=:white,
    xtickfontsize=1,
    ytickfontsize=tick_fontsize,
)
fig_explicit_combine = plot(
    fig_explicit;
    xlabel="",
    xtickfontcolor=:white,
    xtickfontsize=1,
    ytickfontsize=tick_fontsize,
    legend=false,
)
fig_implicit_fluxes_combine = plot(fig_implicit_fluxes; tickfontsize=tick_fontsize)
fig_fluxes_diff_combine = plot(fig_fluxes_diff; tickfontsize=tick_fontsize)
fig_combined = plot(
    fig_implicit_combine,
    fig_explicit_combine,
    fig_implicit_fluxes_combine,
    fig_fluxes_diff_combine;
    layout=l,
    size=(30cm, 29cm),
    link=:x,
    legendfontsize=12,
    titlefontsize=title_fontsize,
    guidefontsize=title_fontsize - 4,
    line=2.5,
)
savefig(fig_combined, figdir("combined.png"))
savefig(fig_combined, figdir("combined.svg"))
fig_combined_transparent = plot(fig_combined; background_color=:transparent)
savefig(fig_combined_transparent, figdir("combined_transparent.svg"))
