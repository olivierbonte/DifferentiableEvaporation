# Forcing and parameters of an eddy-covariance site for EvaporationModel, from the outputs of
# scripts/data_processing in data/exp_pro. Load with `include(srcdir("site_data.jl"))` after
# `@quickactivate "DifferentiableEvaporation"`.
using ComponentArrays
using DataInterpolations
using Dates
using EvaporationModel
using NetCDF
using YAXArrays

site_ec_file(site) = datadir("exp_pro", "eddy_covariance", site * ".nc")
site_soil_file(site) = datadir("exp_pro", "soil", site * "_total_agg.nc")

"""
    site_forcing(file, start_date=typemin(DateTime), end_date=typemax(DateTime); FT=Float64)

Forcing of an eddy-covariance `file` between `start_date` and `end_date` (inclusive; the
whole file by default), with constant piecewise interpolation (`ConstantInterpolation`,
`dir=:left`: the value stamped `t_k` holds on ``[t_k, t_{k+1})``). Returns `(forcings, t, ds)`:
the forcing NamedTuple, the forcing timestamps (local standard time, unix seconds) and the
selected dataset (for observations).
"""
function site_forcing(
    file, start_date=typemin(DateTime), end_date=typemax(DateTime); FT=Float64
)
    ds = readcubedata(open_dataset(file)[time=start_date..end_date])
    t = datetime2unix.(collect(ds.time))
    forcing(x) = ConstantInterpolation(FT.(vec(x)), t; dir=:left)
    forcings = (
        R_n=forcing(ds.Rnet[:]),
        T_a=forcing(ds.Tair[:]),
        p_a=forcing(ds.Psurf[:]),
        u_a=forcing(ds.Wind[:]),
        VPD_a=forcing(FT.(ds.VPD[:]) .* 100), # hPa -> Pa
        P=forcing(ds.Precip[:]),
        SW_in=forcing(ds.SWdown[:]),
        LAI=forcing(ds.LAI[:]),
        lon=FT(ds.longitude[1]),
        utc_offset=FT(NetCDF.ncgetatt(file, "time", "UTC_offset")),
    )
    return forcings, t, ds
end

"""
    model_parameters(; x_clay, w_sat, w_wp, w_fc, w_res, d_2, h, z_obs,
        d_1=of_value_type(x_clay, 0.01), kB⁻¹=log(10), z_0ms=of_value_type(x_clay, 0.01),
        k_ext=of_value_type(x_clay, 0.5),
        vegtype=EvergreenNeedleleafTrees(), FT=Float64)

Model parameters from site properties: ISBA coefficients from the clay percentage `x_clay`
[%] (pedotransfer functions), soil moisture at saturation, wilting point, field capacity and
residual [m³ m⁻³], root-zone depth `d_2` [m], canopy height `h` [m], observation height
`z_obs` [m] and the Jarvis-Stewart parameters of `vegtype`.
"""
function model_parameters(;
    x_clay, w_sat, w_wp, w_fc, w_res, d_2, h, z_obs, d_1=of_value_type(x_clay, 0.01),
    kB⁻¹=log(10), z_0ms=of_value_type(x_clay, 0.01), k_ext=of_value_type(x_clay, 0.5),
    vegtype=EvergreenNeedleleafTrees(), FT=Float64,
)
    x_clay = FT(x_clay)
    veg = VegetationParameters(; vegtype=vegtype)
    d_1 = FT(d_1) # normalization depth of the surface layer [m]
    return ComponentArray(;
        h=FT(h),
        z_0ms=FT(z_0ms),
        w_sat=FT(w_sat),
        a=compute_a(x_clay),
        p_soil=compute_p(x_clay),
        b=compute_b(Clay(), x_clay),
        C_1sat=c_1sat(NoilhanMahfouf96(), x_clay, d_1),
        C_2ref=c_2ref(x_clay),
        C_3=c_3(x_clay),
        d_1=d_1,
        d_2=FT(d_2),
        w_res=FT(w_res),
        w_wp=FT(w_wp),
        w_fc=FT(w_fc),
        z_obs=FT(z_obs),
        kB⁻¹=FT(kB⁻¹),
        # Jarvis-Stewart canopy resistance r_sc = r_smin / (LAI f_1 f_2 f_3 f_4), IFS Table 8.1
        g_d=FT(veg.g_d), # [Pa⁻¹] coefficient of the VPD stress factor f_3 = exp(-g_d VPD)
        r_smin=FT(veg.r_smin), # [s m⁻¹] minimum stomatal resistance
        k_ext=FT(k_ext),
    )
end

"""
    site_parameters(site; vegtype=EvergreenNeedleleafTrees(), FT=Float64,
        soil_file=site_soil_file(site), ec_file=site_ec_file(site))

[`model_parameters`](@ref) of a site: soil properties and root depth from the processed soil
data, canopy and observation heights from the eddy-covariance file.
"""
function site_parameters(
    site; vegtype=EvergreenNeedleleafTrees(), FT=Float64,
    soil_file=site_soil_file(site), ec_file=site_ec_file(site),
)
    ds_soil = open_dataset(soil_file)
    soil = Dict(v => ds_soil[v][1] for v in propertynames(ds_soil))
    ds_ec = open_dataset(ec_file)
    return model_parameters(;
        x_clay=soil[:mean_clay_percentage],
        w_sat=soil[:w_sat],
        w_wp=soil[:w_wp],
        w_fc=soil[:w_fc],
        w_res=soil[:w_res],
        d_2=soil[:root_depth],
        h=ds_ec.canopy_height[1],
        z_obs=ds_ec.reference_height[1],
        vegtype=vegtype,
        FT=FT,
    )
end
