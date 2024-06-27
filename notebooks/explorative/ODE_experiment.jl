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
ds_meteo = open_dataset(datadir("exp_raw","BE-Bra_2004-2014_FLUXNET2015_Met.nc"));
ds_flux = open_dataset(datadir("exp_raw","BE-Bra_2004-2014_FLUXNET2015_Flux.nc"));
#now make a shorter subset of the data
start_date = DateTime(2006,04,01,00)
end_date = DateTime(2006,07,01,00)
ds_meteo_sub = ds_meteo[time = DimensionalData.Between(start_date, end_date)]
ds_flux_sub = ds_flux[time = DimensionalData.Between(start_date, end_date)]

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
C3 = 5.327 * X_clay^(-1.043)

d_1 = 0.1 # [m]
d_2 = 2 # [m], based on description in ds_meteo_sub.properties["soil_type"]

#test soil parameters assuming fixed SM
c_1_test = compute_c_1(0.35, w_sat, b_ch, C1sat)
c_2_test = compute_c_2(0.35, w_sat, C2ref)
w_geq_test = compute_w_geq(0.35, w_sat, a_ch, p_ch)

#site info
z_measur = Float64(ds_flux_sub.reference_height.data[:,:][1])
h_veg = Float64(ds_flux_sub.canopy_height.data[:,:][1])
z0m = 0.1 * h_veg
z0h = Bigleaf.roughness_z0h(z0m, 8) #8 according to Ridgen & Salvucci, 2015
d = 2/3 * h_veg


#Test Penman-Monteith 
# aerodynamic_resistance
param_ra = ComponentArray(
    d = d,
    z0m = z0m,
    z0h = z0h,
    z_measur = z_measur
)
r_a = EvaporationModel.aerodynamic_resistance(
    ds_meteo_sub["Wind"][:],
    param_ra
)

time_sub = collect(ds_meteo_sub.Ti)
plot(time_sub, r_a, xformatter = x -> Dates.format(x, "Y-m-d H:M:S"))
plot(time_sub, ds_meteo_sub.Wind[:])

#surface resistance
param_rs = ComponentArray(
    w_fc = w_fc,
    w_wilt = w_wp,
    gd = gD,
    r_smin = rsmin
)
# forcing_selection = NCDataset.@select(meteo_data, "SWdown", "VPD", "Tair", "LAI")
#make a time selection
forcing = ComponentArray(
    s_in = ds_meteo_sub["SWdown"][:],
    w_2 = 0.25, #constant SM for namd
    vpd = ds_meteo_sub["VPD"][:], #in hPa
    t_air = ds_meteo_sub["Tair"][:],
    lai = ds_meteo_sub["LAI"][:]
)
r_s = EvaporationModel.jarvis_stewart(forcing, param_rs)
plot(time_sub, r_s)
g_s_mol = ms_to_mol.(1.0 ./r_s, ds_meteo_sub.Tair[:], ds_meteo_sub.Psurf[:])
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
function calculate_evaporation(Tair, Psurf, Rnet, VPD, r_a, r_s; kwargs... )
    et = zeros(length(Tair))
    le = zeros(length(Tair))
    for i in 1:length(Tair)
        et[i], le[i] = Bigleaf.potential_ET(
            Val(:PenmanMonteith),
            Tair[i] ,Psurf[i],
            Rnet[i], VPD[i],
            1.0 ./r_a[i],
            G = 0.05*Rnet[i];
            Gs_pot = ms_to_mol(1.0 ./r_s[i], Tair[i] + 273.15, Psurf[i]*1000), #eenheden aanpassen! mol m-2 s-1 moet het zijn 
            kwargs...
        )
    end
    return et, le
end
et_pred, le_pred = calculate_evaporation(
    ds_meteo_sub.Tair[:] .- 273.15, #K to °C
    ds_meteo_sub.Psurf[:]/1000, #Pa -> kPa
    ds_flux_sub.Rnet[:], #W/m²
    ds_meteo_sub.VPD[:]/10, #hPa -> kPa
    r_a, r_s
)
corr_test = Statistics.cor(le_pred, ds_flux_sub.Qle_cor[:])
plot(time_sub, le_pred, label = "predicted")
plot!(time_sub, ds_flux_sub.Qle_cor[:], label = "observed")
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
    dir = :left, extrapolate = true
)
t_plot = collect(range(start_date, end_date, step = Minute(5)))
scatter(collect(ds_flux_sub.Ti), ds_flux_sub["Rnet"][:], label = "Data")
plot!(t_plot, interp_Rnet(Dates.datetime2epochms.(t_plot)); label = "Constant interpolation")

#keep w_2 as constant for now
#same interpolation everywhere
t_interp = Float64.(Dates.datetime2epochms.(collect(ds_flux_sub.Ti)))

@with_kw struct struct_forcing_wg 
    r_net = ConstantInterpolation(ds_flux_sub["Rnet"][:], t_interp, dir = :left)
    vpd = ConstantInterpolation(ds_meteo_sub["VPD"][:], t_interp, dir = :left)
    t_air = ConstantInterpolation(ds_meteo_sub["Tair"][:], t_interp, dir = :left)
    p_surf = ConstantInterpolation(ds_meteo_sub["Psurf"][:], t_interp, dir = :left)
    precip = ConstantInterpolation(ds_meteo_sub["Precip"][:], t_interp, dir = :left)
    wind = ConstantInterpolation(ds_meteo_sub["Wind"][:], t_interp, dir = :left)
    s_in = ConstantInterpolation(ds_meteo_sub["SWdown"][:], t_interp, dir = :left)
    lai = ConstantInterpolation(ds_meteo_sub["LAI"][:], t_interp, dir = :left)
end
param_all = ComponentArray(
    d_1 = d_1,
    d_2 = d_2,
    τ = 24, #[h]
    w_sat = w_sat,
    c1_sat = C1sat,
    c2_ref = C2ref,
    c_3 = C3, 
    a = a_ch,
    b = b_ch,
    p = p_ch, 
    w_fc = w_fc,
    w_wilt = w_wp,
    gd = gD,
    r_smin = rsmin,
    d = d,
    z0m = z0m,
    z0h = z0h,
    z_measur = z_measur
)
#p_total_wg = ComponentArray(forcing = forcing_wg, parameters = param_wg)
 @with_kw struct struct_total_pwg
    forcings = struct_forcing_wg() 
    params = param_all
 end

 function compute_penman(Tair, Psurf, Rnet, VPD, r_a, r_s; kwargs...)
    et, le = Bigleaf.potential_ET(
            Val(:PenmanMonteith), Tair ,Psurf, Rnet, VPD, 1.0 /r_a,
            G = 0.05*Rnet; Gs_pot = ms_to_mol(1.0 ./r_s, Tair + 273.15, Psurf*1000), #eenheden aanpassen! mol m-2 s-1 moet het zijn 
            kwargs...
        )
    return et, le 
 end

 function wg_conservation(du,u,p,t)
    #t in ms because this easiest to work with for now
    @unpack forcings, params = p
    @unpack d_1, d_2, τ, w_sat, c1_sat, c2_ref, c_3, a, b, p, d, z0m, z0h, z_measur, 
            w_fc, w_wilt, gd, r_smin = params
    @unpack r_net, vpd, t_air, p_surf, precip, wind, s_in, lai = forcings
    #w_2 = 0.25 #temporarily
    r_a = EvaporationModel.aerodynamic_resistance(wind(t), d, z0m, z0h, z_measur)
    r_s = EvaporationModel.jarvis_stewart(s_in(t), u[2], vpd(t), t_air(t), lai(t), w_fc, w_wilt, gd, r_smin)
    c_1 = compute_c_1(u[1], w_sat, b, c1_sat)
    c_2 = compute_c_2(u[2], w_sat, c2_ref)
    w_geq = compute_w_geq(u[2], w_sat, a, p) #0.25 is temporarily
    e_g = compute_penman(t_air(t) - 273.15, p_surf(t)/1000, r_net(t), vpd(t)/10, r_a, r_s)[1] #kg/(m^2 s)
    du[1] = 1/1000*(c_1/(d_1 * ρ_w) * (precip(t) - e_g) - c_2 / 86400 * (u[1] - w_geq))
    du[2] = 1/1000*(1 / (ρ_w * d_2) * (precip(t) - e_g) - c_3 / (d_2 * 86400) * max(u[2] - w_fc, 0.0)) 
    return nothing 
 end 

wg_conservation(zeros(2), [0.4, 0.25], struct_total_pwg(), t_interp[1])

w_init = w_sat
prob = ODEProblem(wg_conservation, [w_sat, w_sat], (t_interp[1],t_interp[end]), struct_total_pwg())
sol = solve(prob, ImplicitEuler(), saveat = t_interp, dt = (Millisecond(Minute(30))).value, adaptive = false)
plot(collect(ds_flux_sub.Ti), sol[1,:], label = "w_g")
plot!(collect(ds_flux_sub.Ti), sol[2,:], label = "w_2")
plot!(twinx(), collect(ds_flux_sub.Ti), ds_meteo_sub["Precip"][:], color = :red, linestyle = :dot)
display(sol.destats)