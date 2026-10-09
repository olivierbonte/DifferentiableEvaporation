using DrWatson
@quickactivate "DifferentiableEvaporation"
using DataInterpolations
using Dates
using NetCDF
using CairoMakie
using YAXArrays

FT = Float64
# Read in
site = "BE-Bra"
file_ec = datadir("exp_pro", "eddy_covariance", site * ".nc")
ds_ec = readcubedata(open_dataset(file_ec))
start_date = DateTime(2010, 3, 15)
end_date = DateTime(2010, 3, 25)
ds_ec_sel = ds_ec[time=start_date..end_date]
t_unix = datetime2unix.(ds_ec_sel.time)
t_span_real = (t_unix[1], t_unix[end])
lines(ds_ec_sel.Precip[x=1, y=1])

#Constant Piecewise interpolation
R_n = ConstantInterpolation(FT.(ds_ec_sel.Rnet[:]), t_unix; dir=:left)
T_a = ConstantInterpolation(FT.(ds_ec_sel.Tair[:]), t_unix; dir=:left)
p_a = ConstantInterpolation(FT.(ds_ec_sel.Psurf[:]), t_unix; dir=:left)
u_a = ConstantInterpolation(FT.(ds_ec_sel.Wind[:]), t_unix; dir=:left)
VPD_a = ConstantInterpolation(FT.(ds_ec_sel.VPD[:]) .* 100, t_unix; dir=:left) #hPa -> Pa
P = ConstantInterpolation(FT.(ds_ec_sel.Precip[:]), t_unix; dir=:left)
SW_in = ConstantInterpolation(FT.(ds_ec_sel.SWdown[:]), t_unix; dir=:left)
LAI = ConstantInterpolation(FT.(ds_ec_sel.LAI[:]), t_unix; dir=:left)

# NamedTuple
forcings_real = (
    R_n=R_n,
    T_a=T_a,
    p_a=p_a,
    u_a=u_a,
    VPD_a=VPD_a,
    P=P,
    SW_in=SW_in,
    LAI=LAI,
    lon=FT(ds_ec.longitude[1]),
    utc_offset=FT(NetCDF.ncgetatt(file_ec, "time", "UTC_offset")),
)
