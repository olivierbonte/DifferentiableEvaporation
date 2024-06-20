## Idea of defining the ODE model here
using DrWatson
@quickactivate "DifferentiableEvaporation"
using Revise
using Plots, Dates, Statistics
using EvaporationModel, Bigleaf
using YAXArrays, NetCDF, ComponentArrays, DimensionalData
plotlyjs()



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
start_date = DateTime(2006,07,01,00)
end_date = DateTime(2006,07,06,00)
ds_meteo_sub = ds_meteo[time = DimensionalData.Between(start_date, end_date)]
ds_flux_sub = ds_flux[time = DimensionalData.Between(start_date, end_date)]

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

#test soil parameters assuming fixed SM
c_1_test = C_1(0.35, w_sat, b_ch, C1sat)
c_2_test = C_2(0.35, w_sat, C2ref)
w_geq_test = w_geq(0.35, w_sat, a_ch, p_ch)

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