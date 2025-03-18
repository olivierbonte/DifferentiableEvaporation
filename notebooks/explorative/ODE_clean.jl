using DrWatson
@quickactivate "DifferentiableEvaporation"
using Revise
using Plots, Dates, Statistics, Parameters
using EvaporationModel, Bigleaf, DifferentialEquations
using YAXArrays, NetCDF, ComponentArrays, DimensionalData

# Site of interest: BE-Bra
site = "BE-Bra"
ds_ec = open_dataset(datadir("exp_pro", "eddy_covariance", site * ".nc"))
