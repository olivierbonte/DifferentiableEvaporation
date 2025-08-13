using DrWatson
@quickactivate "DifferentiableEvaporation"
using Revise
using YAXArrays, DimensionalData, NetCDF
using EvaporationModel
using Bigleaf
using Glob
using CSV
using Dates
using DataFrames
using Plots
using LinearAlgebra
using LinearSolve
using Statistics
using LaTeXStrings
using ComponentArrays
using Parameters

# %% Read in data
# FluxDataKit data
site = "DE-Hai"#"BE-Bra"#"ES-LM1"
ds_ec = open_dataset(glob(site * ".nc", datadir("exp_pro", "eddy_covariance"))[1])
ds_soil = open_dataset(
    glob(site * "_soil_horizontal_agg.nc", datadir("exp_pro", "soil"))[1]
)

# Seperate LW_out data 
folder_fluxnet2015 = glob(
    "*" * site * "*FLUXNET2015*", datadir("exp_raw", "fluxnet_for_soil_moisture")
)[1]
df = DataFrame(
    CSV.File(
        glob("*.csv", datadir("exp_raw", "fluxnet_for_soil_moisture", folder_fluxnet2015))[1],
    ),
)
df.TIMESTAMP_START = string.(df.TIMESTAMP_START)
df.TIMESTAMP_START = DateTime.(df.TIMESTAMP_START, dateformat"yyyymmddHHMM") # Set to Date
df.LW_OUT = Float32.(df.LW_OUT)
df_lw_out = df[:, [:TIMESTAMP_START, :LW_OUT, :TS_F_MDS_1]]
replace!(df_lw_out.LW_OUT, -9999 => NaN)

# Add LW_out data to ds_flux dataframe and make common selection
common_timestamps = intersect(collect(ds_ec.Ti), df_lw_out.TIMESTAMP_START)
start_time = minimum(common_timestamps)
stop_time = maximum(common_timestamps)
ds_ec_sel = ds_ec[Ti = start_time .. stop_time]
df_lw_out_sel = filter(row -> start_time <= row.TIMESTAMP_START <= stop_time, df_lw_out)
axlist = (
    Ti(df_lw_out_sel.TIMESTAMP_START),
    # Dim{:variable}("LWup")
)
da_lw_out_sel = YAXArray(axlist, df_lw_out_sel.LW_OUT)
da_tsoil1 = YAXArray(axlist, df_lw_out_sel.TS_F_MDS_1)
# cubes = [ds_meteo_sel.LWdown[x =1,y=1],da_lw_out_sel]
# var_axis = Dim{:variable}(["LWdown","LWup"])
# ds_LW = concatenatecubes(cubes, var_axis)
arrays = Dict(
    :LWdown => ds_ec_sel.LWdown[x = 1, y = 1],
    :LWup => da_lw_out_sel,
    :G => ds_ec_sel.Qg[x = 1, y = 1],
    :LAI => ds_ec_sel.LAI[x = 1, y = 1],
    :f_veg => ds_ec_sel.f_veg,
    :SWC => ds_ec_sel.SWC,
    :Rn => ds_ec_sel.Rnet,
    :Tsoil1 => da_tsoil1,
)
ds_LW = Dataset(; arrays...)
# Start with Fixed emissivity
const ϵ = 0.98f0 # Maes et al., 2019
function calc_Ts(LW_out::T, LW_in::T, ϵ::T) where {T}
    return ((LW_out - (T(1) - ϵ)LW_in) / (T(bigleaf_constants()[:sigma])ϵ))^(T(0.25))
end

# Ts = calc_Ts.(ds_LW[variable = At("LWup")], ds_LW[variable = At("LWdown")], ϵ)
da_Ts = calc_Ts.(ds_LW.LWup, ds_LW.LWdown, ϵ)
arrays[:Ts] = da_Ts
ds_comb = Dataset(; arrays...)

# %% Select a random day to start
# day_sel = Date("2012-08-10")
# ds_comb_day = ds_comb[Ti = Where(x -> Date(x) == day_sel)]
ds_comb_day = ds_comb[Ti = DateTime(2015, 06, 16, 00, 00) .. DateTime(2015, 06, 16, 23, 30)]
#2013-04-15 for max for DE-Hai

# Number of terms in the Fourier series
M_terms = 12 # as in Verhoef (2004)
ω = 2π / (24 * 60 * 60) # radial velocity in rad/s, period of 1 day
t = Dates.value.(Second.(ds_comb_day.Ti .- ds_comb_day.Ti[1]))
Ts_true = collect(ds_comb_day.Ts)

# Fit fourier coefficients
coeffs = fit_fourier_coefficients(t, Ts_true, M_terms, ω)
# Make the Fourier predictions
Ts_fourier = fourier_series.(t, Ref(coeffs), ω)
# Plot the data
plot(ds_comb_day.Ts; label="Original data", xrotation=10)
plot!(collect(ds_comb_day.Ti), Ts_fourier; label="Fourier series: M = $M_terms")
ylabel!("Surface temperature [K]")

# %% Adapted Fourier Coefficents
@unpack a0, an, bn = coeffs
a_bn, ϕ = compute_amplitude_and_phase(an, bn)
Ts_test = a0 .+ sum(a_bn[n] .* sin.(n * ω * t .+ ϕ[n]) for n in 1:M_terms)
plot!(collect(ds_comb_day.Ti), Ts_test; label="Test Fourier")

# %% Calculate the J_B term (harmonics terms surface temperature) of Verhoef (2007) (eq. 8)
# Note that the equation is adapted to account for a delay as explained in Verhoef (2012) 
# and Bhattacharya (2022)
f_s = 1 .- collect(ds_comb_day.f_veg) ./ 100 # soil fraction
Δt = (1 .- f_s) * 1.5 # hours , see Bhattacharya (2022) eq. 4
J_B_no_shift = sum(
    @. A_bn[n] * √(n * ω) * sin(n * ω * t + ϕ_bn[n] + π / 4) for n in 1:M_terms
)
J_B = sum(
    @. A_bn[n] * √(n * ω) * sin(n * ω * t + ϕ_bn[n] + π / 4 - π * Δt / 12) for
    n in 1:M_terms
)
# Illustrate the effect of delaying the shift 
plot(collect(ds_comb_day.Ti), J_B_no_shift; label="no shift")
plot!(collect(ds_comb_day.Ti), J_B; label="shift")

# Calculate estimated summation of harmonic terms J_s
J_s = (1 / 2 * f_s .+ 1 / 2) .* J_B

# %% Calculate thermal inertia using the near surface values
# choice of near surface because of use in Verhoef (2012)
min_measur_depth = ds_ec.SWC.depth[1] * 100 #m -> cm
ssm = ds_comb_day.SWC[depth = At(min_measur_depth / 100)]
info_depths = parse.(
    Float32, hcat(split.([x[1:(end - 2)] for x in collect(ds_soil.depth)], "-")...)
)
start_layer, end_layer = info_depths[1, :], info_depths[2, :]
layer_select = (min_measur_depth .> start_layer) .& (min_measur_depth .≤ end_layer)
categorial_layer = ds_soil.depth[layer_select][1]
ds_soil_top = ds_soil[depth = At(categorial_layer)]
θsat = ds_soil_top.w_sat[1]
f_sand = ds_soil_top.mean_sand_percentage[1] / 100 #% -> fraction
Γ_0 = -1062.4θsat + 1010.8 # Eq. 20a Murray & Verhoef 2007 (I)
Γ_star = 788.2θsat^-1.29 # Eq. 20b Murray & Verhoef 2007 (I)

"""
Values for δ and γ come from Verhoest (2012) (https://doi.org/10.1016/j.agrformet.2011.08.003)
"""
function Kersten(θ::T, θsat::T, f_sand::T) where {T}
    if f_sand > T(0.8)
        # δacc = T(1.78)
        # γacc = T(2.0) / δacc
        δ = T(2.0)
        γ = T(1.78)
    elseif f_sand > T(0.4)
        # δacc = T(0.93)
        # γacc = T(1.5) / δacc
        δ = T(1.5)
        γ = T(0.93)
    else
        # δacc = T(3.84)
        # γacc = T(4.0) / δacc
        δ = T(4.0)
        γ = T(3.8)
    end
    S_r = min(θ / θsat, T(1))
    K_e = exp(γ * (1 - S_r^(γ - δ)))
    return K_e
end

function soil_thermal_inertia(K_e, Γ_star, Γ_0)
    return K_e * (Γ_star - Γ_0) + Γ_0
end
# Try to recreate Fig 2 of Murray and Verhoef 2007 (I)
ssm_range = range(0, θsat, 100)
K_e_range = Kersten.(ssm_range, θsat, f_sand)
plot(ssm_range ./ θsat, K_e_range; xlabel="θ/θsat", ylabel=L"$K_e$ [-]")
Γ_range = soil_thermal_inertia.(K_e, Γ_star, Γ_0)
plot(
    ssm_range ./ θsat,
    Γ_range;
    label=missing,
    xlabel="θ/θsat",
    ylabel=L"$\Gamma$ (J m$^{-2}$ K$^{-1}$ s$^{-0.5}$)",
)

# Apply to chosen day 
K_e = Kersten.(collect(ssm) ./ 100, θsat, f_sand)
Γ = @. K_e * (Γ_star - Γ_0) + Γ_0 # Eq. 19 Murray & Verhoef 2007 (I)
plot(collect(ds_comb_day.Ti), Γ; xrotation=10, label="Γ")
ylabel!("Soil thermal intertia Γ (J m^-2 K-1 s^-0.5)")

# %% Ground heat flux calculation
gr()
G_remote = Γ .* J_s
rmse = √(mean((G_remote - collect(ds_comb_day.G)) .^ 2))
plot(ds_comb_day.G; xrotation=10, label="Observed", ylabel="G [W/m²]")
plot!(collect(ds_comb_day.Ti), G_remote; label="Modelled")
formatted_rmse = round(rmse; digits=3)
title!("RMSE : $formatted_rmse W/m²")
plot!(0.05 .* ds_comb_day.Rn[x = At(1), y = At(1)]; label="Modelled GLEAM", legend=:topleft)
plot!(twinx(), ds_comb_day.Ts; label="Surface temperature", colour=:black)
