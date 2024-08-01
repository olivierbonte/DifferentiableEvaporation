using DrWatson
@quickactivate "DifferentiableEvaporation"
using Revise
using YAXArrays, DimensionalData

# Site info: https://www.bgc-jena.mpg.de/en/servicegroups/fieldexperiements/locations/majadas
# %% Load the data
ds_meteo = open_dataset(datadir("exp_raw", "eddy_covariance", "ES-LM1_2014-2020_FLUXDATAKIT_Met.nc"));
ds_flux = open_dataset(datadir("exp_raw", "eddy_covariance", "ES-LM1_2014-2020_FLUXDATAKIT_Flux.nc"));
ds_soil = open_dataset(datadir("exp_pro", "soil", "ES-LM1_soil_horizontal_agg.nc"))
ds_soil = Cube(ds_soil) #easier for computations!

# %% Parameter definitions
# Initial simplification: just take average over the depths
import Statistics: mean 
ds_soil_avg = mapslices(mean, ds_soil, dims = "depth")
# Table 8.2 of 
#TEST: dynamic f_veg
ds_meteo.cubes[:fveg] = 1 .- exp.(-0.7 .* ds_meteo["LAI"])

# %% make subset
start_date = DateTime(2014, 04, 01, 00)
end_date = DateTime(2014, 07, 01, 00)
ds_meteo_sub = ds_meteo[time=DimensionalData.Between(start_date, end_date)]
ds_flux_sub = ds_flux[time=DimensionalData.Between(start_date, end_date)]
