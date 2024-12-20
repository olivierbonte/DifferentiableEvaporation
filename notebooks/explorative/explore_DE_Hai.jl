using DrWatson
@quickactivate "DifferentiableEvaporation"
using YAXArrays, DimensionalData, NetCDF
using Glob
using Plots; plotlyjs()

site = "DE-Tha"#"BE-Bra"#"ES-LM1"
ds_meteo = open_dataset(glob(site*"*FLUXDATAKIT_Met.nc", datadir("exp_raw", "eddy_covariance"))[1]);
ds_flux = open_dataset(glob(site*"*FLUXDATAKIT_Flux.nc", datadir("exp_raw", "eddy_covariance"))[1]);

plot(collect(ds_meteo.time), ds_meteo.Tair[:] .- 273.15)
ylabel!("T [°C]")

# Inspect shortwave upgoing radiation
plot(collect(ds_meteo.time), ds_flux.SWup[:], ylims = (-10, maximum(ds_flux.SWup[:]) + 10))
ylabel!("SWᵤₚ[W/m²]")

#Filter out snow from rain
bool_snow = ds_meteo.Tair[:] .< 273.15
snow = ds_meteo.Precip[:]
snow[.!bool_snow] .= 0.0
plot(collect(ds_meteo.time), snow * 3600)
ylabel!(" Snowfall [mm/h]")

