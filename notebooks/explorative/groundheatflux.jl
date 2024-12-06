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

# %% Read in data
# FluxDataKit data
site = "DE-Hai"#"BE-Bra"#"ES-LM1"
ds_meteo = open_dataset(glob(site*"*FLUXDATAKIT_Met.nc", datadir("exp_raw", "eddy_covariance"))[1]);
ds_flux = open_dataset(glob(site*"*FLUXDATAKIT_Flux.nc", datadir("exp_raw", "eddy_covariance"))[1]);
ds_soil = open_dataset(glob(site*"_soil_horizontal_agg.nc", datadir("exp_pro", "soil"))[1])
ds_veg = open_dataset(glob(site*".nc", datadir("exp_pro", "vegetation"))[1])

# Seperate LW_out data 
folder_fluxnet2015 = glob("*" * site * "*FLUXNET2015*", datadir("exp_raw","fluxnet_for_soil_moisture"))[1]
df = DataFrame(CSV.File(glob("*.csv", datadir("exp_raw","fluxnet_for_soil_moisture",folder_fluxnet2015))[1]))
df.TIMESTAMP_START = string.(df.TIMESTAMP_START)
df.TIMESTAMP_START = DateTime.(df.TIMESTAMP_START, dateformat"yyyymmddHHMM") # Set to Date
df.LW_OUT = Float32.(df.LW_OUT)
df_lw_out = df[:,[:TIMESTAMP_START,:LW_OUT]]
replace!(df_lw_out.LW_OUT, -9999 => NaN)

# Add LW_out data to ds_flux dataframe and make common selection
common_timestamps = intersect(collect(ds_meteo.Ti), df_lw_out.TIMESTAMP_START)
start_time = minimum(common_timestamps) 
stop_time = maximum(common_timestamps)
ds_meteo_sel = ds_meteo[Ti = start_time .. stop_time]
ds_flux_sel = ds_flux[Ti = start_time .. stop_time]
df_lw_out_sel = filter(row -> start_time <= row.TIMESTAMP_START <= stop_time, df_lw_out)
axlist = (
    Ti(df_lw_out_sel.TIMESTAMP_START),
    # Dim{:variable}("LWup")
)
da_lw_out_sel = YAXArray(axlist, df_lw_out_sel.LW_OUT)
# cubes = [ds_meteo_sel.LWdown[x =1,y=1],da_lw_out_sel]
# var_axis = Dim{:variable}(["LWdown","LWup"])
# ds_LW = concatenatecubes(cubes, var_axis)
arrays = Dict(
    :LWdown => ds_meteo_sel.LWdown[x =1,y=1],
    :LWup => da_lw_out_sel,
    :G => ds_flux_sel.Qg[x=1,y=1],
    :LAI => ds_meteo_sel.LAI[x=1,y=1]
)
ds_LW = Dataset(;arrays...)
# Start with Fixed emissivity
const ϵ = 0.98f0 # Maes et al., 2019
function calc_Ts(LW_out, LW_in, ϵ)
    return ((LW_out - (1.0f0-ϵ)LW_in)/(Float32(bigleaf_constants()[:sigma])ϵ))^(0.25f0)
end

# Ts = calc_Ts.(ds_LW[variable = At("LWup")], ds_LW[variable = At("LWdown")], ϵ)
da_Ts = calc_Ts.(ds_LW.LWup, ds_LW.LWdown, ϵ)
arrays[:Ts] = da_Ts
ds_comb = Dataset(;arrays...)

# %% Select a random day to start
# day_sel = Date("2012-08-10")
# ds_comb_day = ds_comb[Ti = Where(x -> Date(x) == day_sel)]
ds_comb_day = ds_comb[Ti = DateTime(2012,08,15,00,00)..DateTime(2012,08,16,00,00)]


# %% Define the Fourier series function
function fourier_series(t, coeffs, ω)
    k = length(coeffs) - 1
    M = Int(k/2)
    a0, a, b = coeffs[1], coeffs[2:M+1], coeffs[M+2:end]
    series = a0 .+ sum(a[n] .* cos.(n * ω * t) for n in 1:M) + sum(b[n] .* sin.(n * ω * t) for n in 1:M)
    return series, (a0, a, b)
end

# Number of terms in the Fourier series
M_terms = 12 # as in Verhoef (2004)
ω = 2π/(24*60*60) # radial velocity in rad/s, period of 1 day
t = Dates.value.(Second.(ds_comb_day.Ti .- ds_comb_day.Ti[1]))
Ts_true = collect(ds_comb_day.Ts)

# Create the design matrix
X = hcat(ones(length(t)), [cos.(n * ω * t) for n in 1:M_terms]..., [sin.(n * ω * t) for n in 1:M_terms]...)

# Solve the linear system
coeffs = X \ Ts_true

# Print the coefficients
println("Coefficients: ", coeffs)
# Make the Fourier predictions
Ts_fourier, coeff_tuple = fourier_series(t, coeffs, ω)
T_mean = coeff_tuple[1]
A_n = coeff_tuple[2]
B_n = coeff_tuple[3]
# Plot the data
plot(ds_comb_day.Ts, label = "Original data", xrotation = 10)
plot!(collect(ds_comb_day.Ti), Ts_fourier, label = "Fourier series: M = $M_terms" )
ylabel!("Surface temperature [K]")

# %% Adapted Fourier Coefficents
A_bn = @. √(A_n^2 + B_n^2)
ϕ_bn = @. atan(A_n, B_n) #https://en.wikipedia.org/wiki/Atan2
# CRUCIAL to use Atan2 to get angle between positive x axis and
# point defined by cos and sin combination on unit circle! 
#Test to get same fourier series
Ts_test = T_mean .+ sum(A_bn[n] .* sin.(n * ω * t .+ ϕ_bn[n]) for n in 1:M_terms)
plot!(collect(ds_comb_day.Ti), Ts_test, label = "Test Fourier")

# Calculate fractional surface (soil) coverage assuming 0° viewing angle
τ = 0  
β = 0.5 # Radiation extinction coefficient
f_s = exp.(-β*collect(ds_comb_day.LAI)/cos(τ))

# Calculate the J_B term (harmonics terms surface temperature) of Verhoef (2007) (eq. 8)
# Note that the equation is adapted to account for a delay as explained in Verhoef (2012) 
# and Bhattacharya (2022)
Δt = (1 .- f_s) * 1.5 # hours , see Bhattacharya (2022) eq. 4
J_B_no_shift = sum(@. A_bn[n] * √(n * ω) * sin(n * ω * t + ϕ_bn[n] + π/4) for n in 1:M_terms)
J_B = sum(@. A_bn[n] * √(n * ω) * sin(n * ω * t + ϕ_bn[n] + π/4 - π*Δt/12) for n in 1:M_terms)
# Illustrate the effect of delaying the shift 
plot(collect(ds_comb_day.Ti), J_B_no_shift, label = "no shift")
plot!(collect(ds_comb_day.Ti), J_B, label = "shift")


# Calculate estimated summation of harmonic terms J_s
J_s = (0.5 * f_s .+ 0.5) .* J_B

# Calculate thermal inertia using the near surface values
# choice of near surface because of use in Verhoef (2012)
ds_soil_top = ds_soil[depth = At("0-5cm")]
θsat = collect(ds_soil_top.w_sat)[1]
f_sand = collect(ds_soil.mean_sand_percentage[depth = At("0-5cm")])
Γ_0 = -1062.4θsat + 1010.8 # Eq. 20a Murray & Verhoef 2007 (I)
Γ_star = 788.2θsat^1.29 # Eq. 20b Murray & Verhoef 2007 (I)
function Kersten(θ, θsat, f_sand)
    if f_sand > 0.8
        δacc = 1.78
        γacc = 2.0/δ'
    elseif f_sand > 0.4
        δacc = 0.93
        γacc = 1.5/δ'
    else
        δacc = 3.84
        γacc = 4.0/δ'
    end
    S_r = θ/θsat
    K_e = exp(γacc(1-S_r^(γacc-δacc)))
    return K_e
end
K_e = Kersten()
