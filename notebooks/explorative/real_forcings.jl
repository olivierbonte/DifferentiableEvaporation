using DrWatson
@quickactivate "DifferentiableEvaporation"
using YAXArrays, NetCDF

# Read in
site = "BE-Bra"
ds_ec = readcubedata(open_dataset(datadir("exp_pro", "eddy_covariance", site * ".nc")))
start_date = DateTime(2010, 3, 15)
end_date = DateTime(2010, 3, 25)
ds_ec_sel = ds_ec[time=start_date .. end_date]

# Interpolate
