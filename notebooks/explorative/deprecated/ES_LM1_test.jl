using DrWatson
@quickactivate "DifferentiableEvaporation"
using Revise
using YAXArrays, DimensionalData, NetCDF
using EvaporationModel
using DataFrames
using Glob
using Bigleaf

# Site info ES-LM1: https://www.bgc-jena.mpg.de/en/servicegroups/fieldexperiements/locations/majadas
# Site info BE-Bra: https://meta.icos-cp.eu/resources/stations/ES_BE-Bra
# %% Load the data
site = "DE-Hai"#"BE-Bra"#"ES-LM1"
ds_meteo = open_dataset(
    glob(site*"*FLUXDATAKIT_Met.nc", datadir("exp_raw", "eddy_covariance"))[1]
);
ds_flux = open_dataset(
    glob(site*"*FLUXDATAKIT_Flux.nc", datadir("exp_raw", "eddy_covariance"))[1]
);
ds_soil = open_dataset(glob(site*"_soil_horizontal_agg.nc", datadir("exp_pro", "soil"))[1])
ds_veg = open_dataset(glob(site*".nc", datadir("exp_pro", "vegetation"))[1])
ds_soil = Cube(ds_soil) #easier for computations!

# %% Test

# %% Parameter definitions
# Initial simplification: just take average over the depths
import Statistics: mean
ds_soil_avg = mapslices(mean, ds_soil; dims="depth")
# # Table 8.2 of 
# #TEST: dynamic f_veg
# k_ext = 0.7 #see https://github.com/gee-hydro/gee_PML/blob/stable/src/pkg_PMLV2_v0.1.5.js#L422C32-L422C134
# ds_meteo.cubes[:fveg] = 1 .- exp.(-k_ext .* ds_meteo["LAI"])

# Vegetation parameters 
veg_par = EvaporationModel.df_veg
IGBP_land_cover = strip(join(collect(ds_meteo.IGBP_veg_long.data)))
print("Vegetation type: " * IGBP_land_cover)
# Assume Interruped forest here as closest from the ECMWF table for Savanna
# Temporary dict for linking these two land cover classifications
link_dict = Dict(
    "Savannas" => "Interrupted forest", "Mixed Forests" => "Mixed forest/woodland"
)
land_cover_type = link_dict[IGBP_land_cover]
veg_par_sub = veg_par[veg_par.vegetation_type .== land_cover_type, :]
g_d = veg_par_sub.g_d[1]
r_smin = veg_par_sub.r_smin[1]

# Soil parameters
d_1 = 0.1 # [m] #imposed
d_2 = 2 # [m], based on description in ds_meteo_sub.properties["soil_type"]
#this info only present in file from PLUMBER2 (not FLUXDATAKIT)

w_res = ds_soil_avg[Variable = At("w_res")].data[1]
w_wp = ds_soil_avg[Variable = At("w_wp")].data[1]
w_crit = ds_soil_avg[Variable = At("w_crit")].data[1]
w_fc = ds_soil_avg[Variable = At("w_fc")].data[1]
w_sat = ds_soil_avg[Variable = At("w_sat")].data[1]
perc_clay = ds_soil_avg[Variable = At("mean_clay_percentage")].data[1]
# Compute the other parameters
a_ch = compute_a(perc_clay)
b_ch = compute_b(Val(:clay), perc_clay)
p_ch = compute_p(perc_clay)

# %% make subset
start_date = DateTime(2014, 04, 01, 00)
end_date = DateTime(2014, 07, 01, 00)
ds_meteo_sub = ds_meteo[time = DimensionalData.Between(start_date, end_date)]
ds_flux_sub = ds_flux[time = DimensionalData.Between(start_date, end_date)]

#site info
z_measur = Float64(ds_flux.reference_height.data[:, :][1])
h_veg = Float64(ds_flux.canopy_height.data[:, :][1])
z0m = 0.1 * h_veg
kB⁻¹ = 8 #calibrate later 
z0h = Bigleaf.roughness_z0h(z0m, kB⁻¹) #8 according to Ridgen & Salvucci, 2015
#this is quite high though... no difference for 
d = 2 / 3 * h_veg
