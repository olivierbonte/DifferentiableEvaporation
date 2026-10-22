using DrWatson
@quickactivate "DifferentiableEvaporation"
using CairoMakie
using OrdinaryDiffEq
using OrdinaryDiffEqSDIRK


include(srcdir("site_data.jl"))

FT = Float64
site = "DE-Tha"
file_ec = datadir("exp_pro", "eddy_covariance", site * ".nc")

start_date = DateTime(2009, 01, 01, 00, 00)
end_date = DateTime(2009, 12, 31, 23, 30)
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
# Use end of 2009 value as initial value!
# 2010 as the period of interest
start_date = DateTime(2010, 01, 01, 00, 00)
end_date = DateTime(2010, 12, 31, 12, 30)
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
