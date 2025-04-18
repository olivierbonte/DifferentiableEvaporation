using DrWatson
@quickactivate "DifferentiableEvaporation"
using Revise
using Plots, Dates, Statistics, Parameters
using EvaporationModel, Bigleaf, DifferentialEquations
using YAXArrays, NetCDF, ComponentArrays, DimensionalData

FT = Float64
## Site of interest: BE-Bra
site = "BE-Bra"
# Use readcubedata to load into memory!
# See https://juliadatacubes.github.io/YAXArrays.jl/stable/UserGuide/read.html#readcubedata
ds_ec = readcubedata(open_dataset(datadir("exp_pro", "eddy_covariance", site * ".nc")))
# Readcubedata fails here (becuase no dims...) -> convert to simple dict
ds_soil = open_dataset(datadir("exp_pro", "soil", site * "_total_agg.nc"))
dict_soil = Dict()
soil_variables = propertynames(ds_soil)
for variable in soil_variables
    dict_soil[variable] = FT(ds_soil[variable][1])
end

## Define static parameters
# In a later stage, we would want to move this to some form of
# intialisation function
# Vegetation characteristics, using info from https://meta.icos-cp.eu/resources/stations/ES_BE-Bra
veg_param = VegetationParameters(; vegtype=EvergreenNeedleleafTrees())
kB⁻¹ = log(10)
z_obs = FT(ds_ec.reference_height[1]) # m, height of the observations
h = FT(ds_ec.canopy_height[1]) # m, height of the canopy
# Static soil parameters based on PFT
c1sat = c_1sat(dict_soil[:mean_clay_percentage])
c2ref = c_2ref(dict_soil[:mean_clay_percentage])
c3 = c_3(dict_soil[:mean_clay_percentage])
d1 = 0.01 # m, depth of the first layer
d2 = dict_soil[:root_depth] # m, depth of the second layer
z_0ms = 0.01 # m, roughness length for momentum transfer of soil

## Test Peman-Monteith equation at one time step at at time
## in a for loop over several days
## Also test the Shuttleworth&Wallace like model
start_date = DateTime(2010, 3, 5)
end_date = DateTime(2010, 3, 15)
ds_ec_sel = ds_ec[time=start_date .. end_date]
time_indices = 1:length(ds_ec_sel.time)[1]
le_array = similar(collect(time_indices), FT)
le_total_array = similar(le_array)
le_total_p_array = similar(le_array)
for i in time_indices
    # Variable canopy parameters
    rough_dict = Bigleaf.roughness_parameters(
        RoughnessCanopyHeightLAI(), h, ds_ec_sel.LAI[i]; hs=z_0ms
    )
    d_c = rough_dict.d
    z_0mc = rough_dict.z0m
    w_1 = dict_soil[:w_fc] * 2/3# Both soil layers at 2/3 field capacity as test
    w_2 = dict_soil[:w_fc] * 2/3
    ustar = ustar_from_u(FT(ds_ec_sel.Wind[i]), z_obs, d_c, z_0mc)
    r_aa = Bigleaf.compute_Ram(ResistanceWindZr(), ustar, FT(ds_ec_sel.Wind[i]))
    r_ac = (Bigleaf.Gb_constant_kB1(ustar, kB⁻¹))^-1
    r_as = soil_aerodynamic_resistance(Choudhury1988soil(), ustar, h, d_c, z_0mc, z_0ms)
    r_sc = surface_resistance(
        JarvisStewart(),
        FT(ds_ec_sel.SWdown[i]),
        FT(ds_ec_sel.VPD[i]) * 100,
        FT(ds_ec_sel.Tair[i]),
        w_2,
        dict_soil[:w_fc],
        dict_soil[:w_wp],
        FT(ds_ec_sel.LAI[i]),
        veg_param.g_d,
        veg_param.r_smin,
    )
    beta = soil_evaporation_efficiency(Pielke92(), w_1, dict_soil[:w_fc])
    beta_test = soil_evaporation_efficiency(
        Martens17(), w_1, dict_soil[:w_res], dict_soil[:w_crit]
    )
    r_ss = beta_to_r_ss(beta, r_as)
    # Normal Penman Monteith
    et_test, le_test = penman_monteith(
        FT(ds_ec_sel.Tair[i]),
        FT(ds_ec_sel.Psurf[i]),
        FT(ds_ec_sel.VPD[i]) * 100,
        FT(ds_ec_sel.Rnet[i]),
        r_aa + r_ac,
        r_sc,
    )
    le_array[i] = le_test
    # Shuttleworth&Wallace like model
    f_veg = 1 - exp(-0.5 * FT(ds_ec_sel.LAI[i]))
    G = FT(0.05) * FT(ds_ec_sel.Rnet[i]) #temp
    R_nc = f_veg * FT(ds_ec_sel.Rnet[i])
    R_ns = (1 - f_veg) * FT(ds_ec_sel.Rnet[i])
    A_c = R_nc
    A_s = R_ns - G
    f_wet = FT(0.01) # temp
    le_total, le_total_p = total_evaporation(
        FT(ds_ec_sel.Tair[i]),
        FT(ds_ec_sel.Psurf[i]),
        FT(ds_ec_sel.VPD[i]) * 100,
        A_c,
        A_c,
        A_s,
        r_aa,
        r_ac,
        r_as,
        r_sc,
        r_ss,
        f_wet,
    )
    le_total_array[i] = le_total
    le_total_p_array[i] = le_total_p
end
plot(collect(ds_ec_sel.time), le_array; label="λE Bigleaf")
plot!(collect(ds_ec_sel.time), le_total_array; label="λE Shuttleworth&Wallace")
plot!(collect(ds_ec_sel.time), le_total_p_array; label="λE Shuttleworth&Wallace Potential")
plot!(collect(ds_ec_sel.time), ds_ec_sel.Qle_cor[:]; label="λE observed")
