## Idea of defining the ODE model here
using DrWatson
@quickactivate "DifferentiableEvaporation"
using Revise
using Plots, Dates, Statistics, Parameters
using EvaporationModel, Bigleaf, DifferentialEquations
using YAXArrays, NetCDF, ComponentArrays, DimensionalData
gr()



# Test data for this example
# PLUMBER2 dataset: https://dx.doi.org/10.25914/5fdb0902607e1,
# server: https://thredds.nci.org.au/thredds/catalog/ks32/CLEX_Data/PLUMBER2/v1-0/catalog.html
# Be-Bra = example site! mixed forest, for simplicity now take Needle leaf trees

# Download via OpenDAP
# Soil is loamy-sand here -> take sandy-loam as most similar... 
#flux_data_url = "https://thredds.nci.org.au/thredds/dodsC/ks32/CLEX_Data/PLUMBER2/v1-0/Flux/BE-Bra_2004-2014_FLUXNET2015_Flux.nc"
#meteo_data_url = "https://thredds.nci.org.au/thredds/dodsC/ks32/CLEX_Data/PLUMBER2/v1-0/Met/BE-Bra_2004-2014_FLUXNET2015_Met.nc"
#flux_data = NCDataset(flux_data_url)
#meteo_data = NCDataset(meteo_data_url)
#Swith from NCDataset to NetCDF and YAXArray for convenience -> donwload of the data 
#Download the data to a folder
ds_meteo = open_dataset(datadir("exp_raw", "eddy_covariance", "BE-Bra_2004-2014_FLUXNET2015_Met.nc"));
ds_flux = open_dataset(datadir("exp_raw", "eddy_covariance", "BE-Bra_2004-2014_FLUXNET2015_Flux.nc"));

ds_meteo_LM = open_dataset(datadir("exp_raw", "eddy_covariance", "ES-LM1_2014-2020_FLUXDATAKIT_Met.nc"));
#now make a shorter subset of the data
start_date = DateTime(2006, 04, 01, 00)
end_date = DateTime(2006, 07, 01, 00)
ds_meteo_sub = ds_meteo[time=DimensionalData.Between(start_date, end_date)]
ds_flux_sub = ds_flux[time=DimensionalData.Between(start_date, end_date)]

plot(collect(ds_meteo_sub.Ti), ds_meteo_sub["Precip"][:])

# Define parameters for initial experiment
# All defined in Float64 for type stability 
# Inspired by CLASS: short grass, sandy loam (p.129-131)
# Table 10.6: Needle leaf trees
# z0m = 0.02 # m
# z0h = Bigleaf.roughness_z0h(0.02, log(10.0)) # for now: use /10
# d = 10.0 * 0.02 * 2.0 / 3.0 # assuming displacemengt height is 2/3h and z0m = 1/10h
# LAI = 2.0 # [-] -> dynamic LAI from data is take instead!
gD = 0.03
rsmin = 250.0 # s/m (see ECMWF)
fveg = 0.9

# Table 10.7: Sandy Loam
w_res = 0.041 # [-] hihydrosoil
w_crit = 0.2 # [-] temp value
w_sat = 0.472 # [-]
w_fc = 0.323 # [-]
w_wp = 0.171 # [-]
C1sat = 0.132 # [-]
C2ref = 1.8 # [-]
a_ch = 0.219
b_ch = 4.9
p_ch = 4.0
#Get X_Clay from C_2ref
X_clay = 10 # [%]
c_3 = EvaporationModel.compute_c_3(X_clay)

d_1 = 0.1 # [m]
d_2 = 2 # [m], based on description in ds_meteo_sub.properties["soil_type"]

#test soil parameters assuming fixed SM
c_1_test = compute_c_1(0.35, w_sat, b_ch, C1sat)
c_2_test = compute_c_2(0.35, w_sat, C2ref)
w_geq_test = compute_w_geq(0.35, w_sat, a_ch, p_ch)

#site info
z_measur = Float64(ds_flux_sub.reference_height.data[:, :][1])
h_veg = Float64(ds_flux_sub.canopy_height.data[:, :][1])
z0m = 0.1 * h_veg
z0h = Bigleaf.roughness_z0h(z0m, 2) #8 according to Ridgen & Salvucci, 2015
#this is quite high though... no difference for 
d = 2 / 3 * h_veg


#Test Penman-Monteith 
# aerodynamic_resistance
param_ra = ComponentArray(
    d=d,
    z0m=z0m,
    z0h=z0h,
    z_measur=z_measur
)
r_a = EvaporationModel.aerodynamic_resistance(
    ds_meteo_sub["Wind"][:],
    param_ra
)

time_sub = collect(ds_meteo_sub.Ti)
plot(time_sub, r_a, xformatter=x -> Dates.format(x, "Y-m-d H:M:S"))
plot(time_sub, ds_meteo_sub.Wind[:])

#surface resistance
param_rs = ComponentArray(
    w_fc=w_fc,
    w_wilt=w_wp,
    gd=gD,
    r_smin=rsmin
)
# forcing_selection = NCDataset.@select(meteo_data, "SWdown", "VPD", "Tair", "LAI")
#make a time selection
forcing = ComponentArray(
    s_in=ds_meteo_sub["SWdown"][:],
    w_2=0.25, #constant SM for namd
    vpd=ds_meteo_sub["VPD"][:], #in hPa
    t_air=ds_meteo_sub["Tair"][:],
    lai=ds_meteo_sub["LAI"][:]
)
r_s = EvaporationModel.jarvis_stewart(forcing, param_rs)
plot(time_sub, r_s)
g_s_mol = ms_to_mol.(1.0 ./ r_s, ds_meteo_sub.Tair[:], ds_meteo_sub.Psurf[:])
plot(time_sub, g_s_mol)
# using ComponentArrays, Parameters
# test = ComponentArray(a = [3,4], b = [5,6])
# @unpack a,b = test
# et_test = Bigleaf.potential_ET.(
#     Val(:PenmanMonteith),
#     ds_meteo_sub.Tair[:] .- 273.15, #K to °C
#     ds_meteo_sub.Psurf[:]/1000, #kPa
#     ds_flux_sub.Rnet[:], 
#     ds_meteo_sub.VPD[:], #hPa -> kPa
#     1.0 ./r_a;
#     Gs_pot = 1.0 ./10, #eenheden aanpassen! mol m-2 s-1 moet het zijn 
# )
#Custom function to wrap Bigleaf function
function calculate_evaporation(Tair, Psurf, Rnet, VPD, r_a, r_s; kwargs...)
    et = zeros(length(Tair))
    le = zeros(length(Tair))
    for i in 1:length(Tair)
        et[i], le[i] = Bigleaf.potential_ET(
            Val(:PenmanMonteith),
            Tair[i], Psurf[i],
            Rnet[i], VPD[i],
            1.0 ./ r_a[i],
            G=0.05 * Rnet[i];
            Gs_pot=ms_to_mol(1.0 ./ r_s[i], Tair[i] + 273.15, Psurf[i] * 1000), #eenheden aanpassen! mol m-2 s-1 moet het zijn 
            kwargs...
        )
    end
    return et, le
end
et_pred, le_pred = calculate_evaporation(
    ds_meteo_sub.Tair[:] .- 273.15, #K to °C
    ds_meteo_sub.Psurf[:] / 1000, #Pa -> kPa
    ds_flux_sub.Rnet[:], #W/m²
    ds_meteo_sub.VPD[:] / 10, #hPa -> kPa
    r_a, r_s
)
corr_test = Statistics.cor(le_pred, ds_flux_sub.Qle_cor[:])
plot(time_sub, le_pred, label="predicted")
plot!(time_sub, ds_flux_sub.Qle_cor[:], label="observed")
#plot!(time_sub, ds_flux_sub.Rnet[:], label = "Rₙ")
yaxis!("Latent heat [W/m²]")
using Printf
title!(Printf.@sprintf("Correlation: %.4f", corr_test))
#title!("Correlation: $corr_test%.4f")
# NCDataset.@select(meteo_data["Tair"], start_date .≤ time .≤ end_date)

scatter(le_pred, ds_flux_sub.Qle_cor[:])
xaxis!("Predicted LE")
yaxis!("Observed LE")
plot!(ds_flux_sub.Qle_cor[:], ds_flux_sub.Qle_cor[:])

## Experiment with ODE for w_g
using DataInterpolations, Dates
#convert to numeric ms for interpolation 
interp_Rnet = ConstantInterpolation(
    ds_flux_sub["Rnet"][:],
    Dates.datetime2epochms.(collect(ds_flux_sub.Ti)),
    dir=:left, extrapolate=true
)
t_plot = collect(range(start_date, end_date, step=Minute(5)))
scatter(collect(ds_flux_sub.Ti), ds_flux_sub["Rnet"][:], label="Data")
plot!(t_plot, interp_Rnet(Dates.datetime2epochms.(t_plot)); label="Constant interpolation")
#keep w_2 as constant for now
#same interpolation everywhere
#t_interp = Float64.(Dates.datetime2epochms.(collect(ds_flux_sub.Ti)))
t_interp = Dates.value.(Dates.convert.(Second, ds_flux_sub.Ti .- ds_flux_sub.Ti[1])) #numer of seconds since start in Int64
@with_kw struct struct_forcing_wg
    r_net::DataInterpolations.AbstractInterpolation = ConstantInterpolation(ds_flux_sub["Rnet"][:], t_interp, dir=:left)
    vpd::DataInterpolations.AbstractInterpolation = ConstantInterpolation(ds_meteo_sub["VPD"][:], t_interp, dir=:left)
    t_air::DataInterpolations.AbstractInterpolation = ConstantInterpolation(ds_meteo_sub["Tair"][:], t_interp, dir=:left)
    p_surf::DataInterpolations.AbstractInterpolation = ConstantInterpolation(ds_meteo_sub["Psurf"][:], t_interp, dir=:left)
    precip::DataInterpolations.AbstractInterpolation = ConstantInterpolation(ds_meteo_sub["Precip"][:], t_interp, dir=:left)
    wind::DataInterpolations.AbstractInterpolation = ConstantInterpolation(ds_meteo_sub["Wind"][:], t_interp, dir=:left)
    s_in::DataInterpolations.AbstractInterpolation = ConstantInterpolation(ds_meteo_sub["SWdown"][:], t_interp, dir=:left)
    lai::DataInterpolations.AbstractInterpolation = ConstantInterpolation(ds_meteo_sub["LAI"][:], t_interp, dir=:left)
end
#ALTERNATIVE: JUST DEFINE THESE FUNCTIONS OUTSIDE OF THE 
# @with_kw struct struct_test
#     test::DataInterpolations.AbstractInterpolation
# end
param_all = ComponentArray(
    d_1=d_1,
    d_2=d_2,
    τ=24.0*3600, #[s]
    w_sat=w_sat,
    c1_sat=C1sat,
    c2_ref=C2ref,
    c_3=c_3,
    a=a_ch,
    b=b_ch,
    p=p_ch,
    w_res=w_res,
    w_crit = w_crit,
    w_fc=w_fc,
    w_wilt=w_wp,
    gd=gD,
    r_smin=rsmin,
    d=d,
    z0m=z0m,
    z0h=z0h,
    z_measur=z_measur
)
#p_total_wg = ComponentArray(forcing = forcing_wg, parameters = param_wg)
@with_kw struct struct_total_pwg
    forcings::struct_forcing_wg = struct_forcing_wg()
    params::ComponentArray = param_all
end

function compute_penman(Tair, Psurf, Rnet, VPD, r_a, r_s; kwargs...)
    et, le = Bigleaf.potential_ET(
        Val(:PenmanMonteith), Tair, Psurf, Rnet, VPD, 1.0 / r_a,
        G=0.05 * Rnet; Gs_pot=ms_to_mol(1.0 ./ r_s, Tair + 273.15, Psurf * 1000), #eenheden aanpassen! mol m-2 s-1 moet het zijn 
        kwargs...
    )
    return et, le
end

using StaticArrays
using SimulationLogs
function wg_conservation(u, p, t)
    #t in ms because this easiest to work with for now
    #@unpack forcings, params = p
    @unpack d_1, d_2, τ, w_sat, c1_sat, c2_ref, c_3, a, b, p, d, z0m, z0h, z_measur,
    w_res, w_crit, w_fc, w_wilt, gd, r_smin = p
    #@unpack r_net, vpd, t_air, p_surf, precip, wind, s_in, lai = forcings
    # DO NOT UNPACK FUNCTION HERE, THIS IS SUPER SLOW! 
    #w_2 = 0.25 #temporarily
    r_a = EvaporationModel.aerodynamic_resistance(wind(t), d, z0m, z0h, z_measur)
    r_s = EvaporationModel.jarvis_stewart(s_in(t), u[2], vpd(t), t_air(t), lai(t), w_fc, w_wilt, gd, r_smin)
    c_1 = compute_c_1(u[1], w_sat, b, c1_sat)
    c_2 = compute_c_2(u[2], w_sat, c2_ref)
    w_geq = compute_w_geq(u[2], w_sat, a, p) #0.25 is temporarily
    f_veg = 0.95 #temp value
    f_wet = 0.0 #temp value  
    e_tr = EvaporationModel.compute_transpiration(
        compute_penman(
            t_air(t) - 273.15, p_surf(t) / 1000, r_net(t), vpd(t) / 10, r_a, r_s
        )[1], #kg/(m^2 s)
        f_veg, f_wet
    )
    e_g = EvaporationModel.compute_bare_soil_evaporation(
        compute_penman(
            t_air(t) - 273.15, p_surf(t) / 1000, r_net(t), vpd(t) / 10, r_a, r_s
        )[1], #kg/(m^2 s)
        f_veg,
        EvaporationModel.compute_soil_evaporation_stress(
            u[1], w_crit, w_res
        ),
    )
    #correct fluxes:
    # - Bare soil evaporation to zero when below residual SM
    # - Transpirtation to zero when w_2 below wilting point
    w_min = min(u[1], u[2])
    e_g = max(w_min - w_res, 0.0) * e_g
    e_tr = max(u[2] - w_wilt, 0.0) * e_tr
    @log e_g_wm2 = Bigleaf.ET_to_LE(e_g, t_air(t))
    @log e_tr_wm2 = Bigleaf.ET_to_LE(e_tr, t_air(t))
    dwg = c_1 / (d_1 * ρ_w) * (precip(t) - e_g) - c_2 / τ * (u[1] - w_geq)
    dw2 = 1 / (ρ_w * d_2) * (precip(t) - e_g - e_tr) - c_3 / (d_2 * τ) * max(u[2] - w_fc, 0.0)
    SA[dwg, dw2]
end
#input_struct = struct_total_pwg()
@unpack r_net, vpd, t_air, p_surf, precip, wind, s_in, lai = struct_forcing_wg()
using BenchmarkTools
@benchmark wg_conservation(SA[0.4, 0.25], param_all, t_interp[1])

w_init = w_fc
prob = ODEProblem(wg_conservation, SA[w_init, w_init], (t_interp[1], t_interp[end]), param_all)
sol = solve(prob, ImplicitEuler(), saveat=t_interp, dt=(Second(Minute(30))).value, adaptive=false)
out = get_log(sol)
scope(sol, [:e_g_wm2, :e_tr_wm2])

#try implementing saving callback
using DiffEqCallbacks
saved_values = SavedValues(Float64, Tuple{Float64, Float64})
# saved_values = SavedValues(Float64, NamedTuple{(:value1, :value2), Tuple{Float64, Float64}})
function save_func(u, t , integrator)
    @unpack d, z0m, z0h, z_measur,
    w_res, w_crit, w_fc, w_wilt, gd, r_smin = integrator.p
    r_a = EvaporationModel.aerodynamic_resistance(wind(t), d, z0m, z0h, z_measur)
    r_s = EvaporationModel.jarvis_stewart(s_in(t), u[2], vpd(t), t_air(t), lai(t), w_fc, w_wilt, gd, r_smin) 
    f_veg = 0.95 #temp value
    f_wet = 0.0 #temp value  
    e_tr = EvaporationModel.compute_transpiration(
        compute_penman(
            t_air(t) - 273.15, p_surf(t) / 1000, r_net(t), vpd(t) / 10, r_a, r_s
        )[2], #W/m^2
        f_veg, f_wet
    )
    e_g = EvaporationModel.compute_bare_soil_evaporation(
        compute_penman(
            t_air(t) - 273.15, p_surf(t) / 1000, r_net(t), vpd(t) / 10, r_a, r_s
        )[2], #W/m^2
        f_veg,
        EvaporationModel.compute_soil_evaporation_stress(
            u[1], w_crit, w_res
        ),
    )
    #correct fluxes:
    # - Bare soil evaporation to zero when below residual SM
    # - Transpirtation to zero when w_2 below wilting point
    w_min = min(u[1], u[2])
    e_g = max(w_min - w_res, 0.0) * e_g
    e_tr = max(u[2] - w_wilt, 0.0) * e_tr
    return (e_g, e_tr)   
end
cb = SavingCallback(save_func, saved_values, saveat = t_interp)

w_init = w_fc
prob = ODEProblem(wg_conservation, SA[w_init, w_init], (t_interp[1], t_interp[end]), param_all)
sol = solve(prob, ImplicitEuler(), saveat=t_interp, dt=(Second(Minute(30))).value, adaptive=false, callback = cb)
display(sol.destats)

#plot
color1, color2 = :blue, :black
p = plot(collect(ds_flux_sub.Ti), sol[1, :], label="w_g", color=color1, dpi = 400)
plot!(
    collect(ds_flux_sub.Ti), sol[2, :], label="w_2", legend=:topleft, linestyle=:dash,
    color=color1, yguidefontcolor=color1, ylabel="soil moisture [-]", yforeground_color_axis=color1,
    ytickfontcolor=color1
)
plot!(
    twinx(), collect(ds_flux_sub.Ti), ds_meteo_sub["Precip"][:], color=:black,
    linestyle=:dot, ylabel="precipitation [kg/(m² s)]", label=""
)
savefig(projectdir("plots","test_picture.png"))
# plot 

#plot saved values
E_total = map(x -> x[1], saved_values.saveval) + map(x -> x[2], saved_values.saveval)
plot(collect(ds_flux_sub.Ti), ds_flux_sub.Qle_cor[:], label = "observed", ylabel = "LE [W/m²]", dpi = 400)
plot!(collect(ds_flux_sub.Ti), E_total, label = "predicted", linestyle = :dash)
corr_model = Statistics.cor(E_total, ds_flux_sub.Qle_cor[:])
title!(Printf.@sprintf("Correlation: %.4f", corr_model))
savefig(projectdir("plots","test_picture_LE.png"))

## Try different solvers for the same problem
sol_implicit = solve(prob, ImplicitEuler(), saveat=t_interp, dt=(Second(Minute(30))).value, adaptive=false)
sol_explicit = solve(prob, Euler(), saveat=t_interp, dt=(Second(Minute(30))).value)
sol_stiff = solve(prob, KenCarp4(), saveat=t_interp)

#benchmark
# using LineSearches
# display(@benchmark solve(prob, ImplicitEuler(nlsolve = NLNewton(relax=BackTracking())), saveat=t_interp, 
#     dt=(Second(Minute(30))).value, adaptive=false))
display(@benchmark solve(prob, ImplicitEuler(), saveat=t_interp, 
    dt=(Second(Minute(30))).value, adaptive=false))
display(@benchmark solve(prob, Euler(), saveat=t_interp, dt=(Second(Minute(30))).value))
# display(@benchmark solve(prob, QNDF(), saveat=t_interp))
display(@benchmark solve(prob, Rodas5P(), saveat=t_interp))
#Rodas4P() also good (more robust)

# test out of domain callback!
# More info on how to constrain the solution, see https://doi.org/10.1016/j.amc.2004.12.011
sol_Tsit5 = solve(prob, Tsit5(), cb = PositiveDomain(), abstol = 1e-8, reltol = 1e-8)
plot(sol_Tsit5)
# Works somewhat, but this solver can not handle the stiffness...
sol_Heun = solve(prob, Heun(), cb = PositiveDomain(), abstol = 1e-6, reltol = 1e-5)
plot(sol_Heun)
# https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/
# The stiff solver Rodas5P with reltol set to 1e-4
# This mean accurate op to 0.0001 for soil moisture
# so max 0.0001 difference between the two orders used 

sol = solve(prob, alg_hints = [:stiff])
plot(sol)