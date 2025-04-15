using DrWatson
@quickactivate "DifferentiableEvaporation"
using Revise
using Plots, Dates, Statistics, Parameters
using EvaporationModel, Bigleaf, DifferentialEquations
using YAXArrays, NetCDF, ComponentArrays, DimensionalData

FT = Float64
## Site of interest: BE-Bra
site = "BE-Bra"
ds_ec = open_dataset(datadir("exp_pro", "eddy_covariance", site * ".nc"))
ds_soil = open_dataset(datadir("exp_pro", "soil", site * "_total_agg.nc"))

## Define static parameters
# In a later stage, we would want to move this to some form of 
# intialisation function 
# Vegetation characteristics, using info from https://meta.icos-cp.eu/resources/stations/ES_BE-Bra
veg_param = VegetationParameters(vegtype = EvergreenNeedleleafTrees())
kB⁻¹ = log(10) 
z_obs = FT(ds_ec.reference_height[1]) # m, height of the observations
h = FT(ds_ec.canopy_height[1]) # m, height of the canopy
# Static soil parameters based on PFT 
c1sat= c_1sat(ds_soil.mean_clay_percentage[1])
c2ref= c_2ref(ds_soil.mean_clay_percentage[1])
c3= c_3(ds_soil.mean_clay_percentage[1])
d1 = 0.01 # m, depth of the first layer 
d2 = ds_soil.root_depth[1] # m, depth of the second layer
z_0ms = 0.01 # m, roughness length for momentum transfer of soil

## Test Peman-Monteith equation at one time step
# Variable canopy parameters
rough_dict = Bigleaf.roughness_parameters(RoughnessCanopyHeightLAI(), h, ds_ec.LAI[1]; hs = z_0ms)
d_c = rough_dict.d
z_0mc = rough_dict.z0m
w_1 = ds_soil.w_fc[1] # Both soil layers at field capacity as test
w_2 = ds_soil.w_fc[1]
ustar = ustar_from_u(FT(ds_ec.Wind[1]), z_obs, d_c, z_0mc)
r_aa = Bigleaf.compute_Ram(ResistanceWindZr(), ustar, FT(ds_ec.Wind[1]))
r_ac = (Bigleaf.Gb_constant_kB1(ustar, kB⁻¹))^-1
# r_as = soil_aerodynamic_resistance(Choudhury1988soil(), ustar