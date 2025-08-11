using DrWatson
@quickactivate "DifferentiableEvaporation"
using Revise
using Plots, Dates, Statistics, Parameters
using BenchmarkTools
using YAXArrays, NetCDF, DimensionalData
using ADTypes,
    ComponentArrays, DataInterpolations, OrdinaryDiffEq, DiffEqCallbacks, StaticArrays
using Bigleaf
using EvaporationModel

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
C_1sat = c_1sat(dict_soil[:mean_clay_percentage])
C_2ref = c_2ref(dict_soil[:mean_clay_percentage])
C_3 = c_3(dict_soil[:mean_clay_percentage])
a = compute_a(dict_soil[:mean_clay_percentage])
p = compute_p(dict_soil[:mean_clay_percentage])
b = compute_b(Clay(), dict_soil[:mean_clay_percentage])
d_1 = 0.01 # m, depth of the first layer
d_2 = dict_soil[:root_depth] # m, depth of the second layer
z_0ms = 0.01 # m, roughness length for momentum transfer of soil

## Test Peman-Monteith equation at one time step at at time
## in a for loop over several days
## Also test the Shuttleworth&Wallace like model
# Goede data voor geen errors: 2009-03-15 - 2009-03-20,
# bodemvocht = 1/2 aan veldcapaciteit

start_date = DateTime(2010, 3, 15)
end_date = DateTime(2010, 3, 25)
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
    w_1 = dict_soil[:w_fc] * 2 / 3# Both soil layers at 2/3 field capacity as test
    w_2 = dict_soil[:w_fc] * 2 / 3
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

## Test one timestep of the complete model
# REPEAT FROM ABOVE
i = 1
# Variable canopy parameters
rough_dict = Bigleaf.roughness_parameters(
    RoughnessCanopyHeightLAI(), h, ds_ec_sel.LAI[i]; hs=z_0ms
)
d_c = rough_dict.d
z_0mc = rough_dict.z0m
w_1 = dict_soil[:w_fc] * 3 / 4# Both soil layers at 2/3 field capacity as test
w_2 = dict_soil[:w_fc] * 3 / 4
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
# END REPREAT
w_r = FT(0.2) * FT(ds_ec_sel.LAI[i]) * 1 / 6
f_veg = fractional_vegetation_cover(FT(ds_ec_sel.LAI[i]))
f_wet = fraction_wet_vegetation(w_r, FT(ds_ec_sel.LAI[i]))
G = ground_heat_flux(Allen07(), FT(ds_ec_sel.Rnet[i]), FT(ds_ec_sel.Tair[i]))
A, A_c, A_s = available_energy_partioning(FT(ds_ec_sel.Rnet[i]), G, f_veg)
λE_tot, λE_tot_p = total_evaporation(
    FT(ds_ec_sel.Tair[i]),
    FT(ds_ec_sel.Psurf[i]),
    FT(ds_ec_sel.VPD[i]) * 100,
    A,
    A_c,
    A_s,
    r_aa,
    r_ac,
    r_as,
    r_sc,
    r_ss,
    f_wet,
)
VPD_m = vpd_veg_source_height(
    FT(ds_ec_sel.VPD[i]) * 100,
    FT(ds_ec_sel.Tair[i]),
    FT(ds_ec_sel.Psurf[i]),
    A,
    λE_tot,
    r_aa,
)
E_t, λE_t = transpiration(
    FT(ds_ec_sel.Tair[i]), FT(ds_ec_sel.Psurf[i]), VPD_m, A_c, r_ac, r_sc, f_wet
)
E_i, λE_i = interception(
    FT(ds_ec_sel.Tair[i]), FT(ds_ec_sel.Psurf[i]), VPD_m, A_c, r_ac, f_wet
)
E_s, λE_s = soil_evaporation(
    FT(ds_ec_sel.Tair[i]), FT(ds_ec_sel.Psurf[i]), VPD_m, A_s, r_as, r_ss
)
D_c = canopy_drainage(FT(ds_ec_sel.Precip[i]), w_r, f_veg)
P_s = precip_below_canopy(FT(ds_ec_sel.Precip[i]), f_veg, D_c)
# Soil fluxes
Q_s = surface_runoff(StaticInfiltration(), P_s, w_2, FT(dict_soil[:w_fc]))
w_1eq = w_geq(w_2, dict_soil[:w_sat], a, p)
C_1 = c_1(w_1, dict_soil[:w_sat], b, C_1sat)
C_2 = c_2(w_2, dict_soil[:w_sat], C_2ref)
D_1 = diffusion_layer_1(w_1, w_1eq, C_2)
K_2 = vertical_drainage_layer_2(w_2, dict_soil[:w_fc], C_3, d_2)
I_s = P_s - Q_s
dw1dt = C_1 / (ρ_w * d_1) * (I_s - E_s) - D_1
dw2dt = 1 / (ρ_w * d_2) * (I_s - E_s - E_t) - K_2
dwrdt = f_veg * FT(ds_ec_sel.Precip[i]) - E_i - D_c

## Final experiment: constructing a proper ODE system
# Interpolate the forcings (using constant interpolation for now)
t_unix = datetime2unix.(ds_ec_sel.time[:]) # in seconds!
# R_n = ConstantInterpolation(FT.(ds_ec_sel.Rnet[:]), t_unix; dir=:left)
# T_a = ConstantInterpolation(FT.(ds_ec_sel.Tair[:]), t_unix; dir=:left)
# p_a = ConstantInterpolation(FT.(ds_ec_sel.Psurf[:]), t_unix; dir=:left)
# u_a = ConstantInterpolation(FT.(ds_ec_sel.Wind[:]), t_unix; dir=:left)
# VPD_a = ConstantInterpolation(FT.(ds_ec_sel.VPD[:]) .* 100, t_unix; dir=:left) #hPa -> Pa
# P = ConstantInterpolation(FT.(ds_ec_sel.Precip[:]), t_unix; dir=:left)
# SW_in = ConstantInterpolation(FT.(ds_ec_sel.SWdown[:]), t_unix; dir=:left)
# LAI = ConstantInterpolation(FT.(ds_ec_sel.LAI[:]), t_unix; dir=:left)

# R_n = PCHIPInterpolation(FT.(ds_ec_sel.Rnet[:]), t_unix)
# T_a = PCHIPInterpolation(FT.(ds_ec_sel.Tair[:]), t_unix)
# p_a = PCHIPInterpolation(FT.(ds_ec_sel.Psurf[:]), t_unix)
# u_a = PCHIPInterpolation(FT.(ds_ec_sel.Wind[:]), t_unix)
# VPD_a = PCHIPInterpolation(FT.(ds_ec_sel.VPD[:]) .* 100, t_unix) #hPa -> Pa
# P = PCHIPInterpolation(FT.(ds_ec_sel.Precip[:]), t_unix)
# SW_in = PCHIPInterpolation(FT.(ds_ec_sel.SWdown[:]), t_unix)
# LAI = PCHIPInterpolation(FT.(ds_ec_sel.LAI[:]), t_unix)

# Alternative: preserver integral by using smoothed constant interpolation
d_max = diff(t_unix)[1] * 0.2 # The maximum smoothing distance in the t direction from the data
# points in seconds, here as a fraction of the time interval
R_n = SmoothedConstantInterpolation(FT.(ds_ec_sel.Rnet[:]), t_unix; d_max=d_max)
T_a = SmoothedConstantInterpolation(FT.(ds_ec_sel.Tair[:]), t_unix; d_max=d_max)
p_a = SmoothedConstantInterpolation(FT.(ds_ec_sel.Psurf[:]), t_unix; d_max=d_max)
u_a = SmoothedConstantInterpolation(FT.(ds_ec_sel.Wind[:]), t_unix; d_max=d_max)
VPD_a = SmoothedConstantInterpolation(FT.(ds_ec_sel.VPD[:]) .* 100, t_unix; d_max=d_max) #hPa -> Pa
P = SmoothedConstantInterpolation(FT.(ds_ec_sel.Precip[:]), t_unix; d_max=d_max)
SW_in = SmoothedConstantInterpolation(FT.(ds_ec_sel.SWdown[:]), t_unix; d_max=d_max)
LAI = SmoothedConstantInterpolation(FT.(ds_ec_sel.LAI[:]), t_unix; d_max=d_max)

# Pack alle the parameters into component array
param = ComponentArray(;
    d_1=d_1,
    d_2=d_2,
    C_1sat=C_1sat,
    C_2ref=C_2ref,
    C_3=C_3,
    w_sat=dict_soil[:w_sat],
    w_fc=dict_soil[:w_fc],
    w_wp=dict_soil[:w_wp],
    w_res=dict_soil[:w_res],
    x_clay=dict_soil[:mean_clay_percentage],
    a=a,
    p=p,
    b=b,
    h=h,
    z_a=z_obs,
    z_0ms=z_0ms,
    kB⁻¹=kB⁻¹,
    r_smin=veg_param.r_smin,
    g_d=veg_param.g_d,
)

# function conservation_equations!(du, u, p, t)

# end

# Construct the ODE system function
function calculate_fluxes(u, p, t)
    ## Unpack the state variables
    w_1, w_2, w_r = u
    ## Unpack the static parameters
    @unpack d_1,
    d_2,
    C_1sat,
    C_2ref,
    C_3,
    w_sat,
    w_fc,
    w_wp,
    w_res,
    x_clay,
    p,
    b,
    h,
    z_a,
    z_0ms,
    kB⁻¹,
    r_smin,
    g_d = p
    ## Compute the dynamic parameters
    # Vegtation characteristics
    d_c, z_0mc = Bigleaf.roughness_parameters(
        RoughnessCanopyHeightLAI(), h, LAI(t); hs=z_0ms
    )
    f_veg = fractional_vegetation_cover(LAI(t))
    w_rmax = max_canopy_capacity(LAI(t))
    f_wet = fraction_wet_vegetation(w_r, w_rmax)

    # Soil characteristics
    w_1eq = w_geq(w_2, w_sat, a, p)
    C_1 = c_1(w_1, w_sat, b, C_1sat)
    C_2 = c_2(w_2, w_sat, C_2ref)

    # Energy partitioning
    G = ground_heat_flux(Allen07(), R_n(t), LAI(t))
    A, A_c, A_s = available_energy_partioning(R_n(t), G, f_veg)

    # Resistances
    ustar = ustar_from_u(u_a(t), z_obs, d_c, z_0mc)
    r_aa = Bigleaf.compute_Ram(ResistanceWindZr(), ustar, u_a(t))
    r_ac = (Bigleaf.Gb_constant_kB1(ustar, kB⁻¹))^-1
    r_as = soil_aerodynamic_resistance(Choudhury1988soil(), ustar, h, d_c, z_0mc, z_0ms)
    r_sc = surface_resistance(
        JarvisStewart(), SW_in(t), VPD_a(t), T_a(t), w_2, w_fc, w_wp, LAI(t), g_d, r_smin
    )
    β = soil_evaporation_efficiency(Pielke92(), w_1, w_fc)
    r_ss = beta_to_r_ss(β, r_as)

    # Turbulent fluxes calculations
    λE_tot, λE_tot_p = total_evaporation(
        T_a(t), p_a(t), VPD_a(t), A, A_c, A_s, r_aa, r_ac, r_as, r_sc, r_ss, f_wet
    )
    VPD_m = vpd_veg_source_height(VPD_a(t), T_a(t), p_a(t), A, λE_tot, r_aa)
    E_t, λE_t = transpiration(T_a(t), p_a(t), VPD_m, A_c, r_ac, r_sc, f_wet)
    E_i, λE_i = interception(T_a(t), p_a(t), VPD_m, A_c, r_ac, f_wet)
    E_s, λE_s = soil_evaporation(T_a(t), p_a(t), VPD_m, A_s, r_as, r_ss)

    # Soil and canopy fluxes
    D_c = canopy_drainage(P(t), w_r, f_veg)
    P_s = precip_below_canopy(P(t), f_veg, D_c)
    Q_s = surface_runoff(StaticInfiltration(), P(t), w_2, w_fc)
    D_1 = diffusion_layer_1(w_1, w_1eq, C_2)
    K_2 = vertical_drainage_layer_2(w_2, w_fc, C_3, d_2)
    I_s = P_s - Q_s
    return (
        d_c=d_c,
        z_0mc=z_0mc,
        w_rmax=w_rmax,
        f_veg=f_veg,
        f_wet=f_wet,
        w_1eq=w_1eq,
        C_1=C_1,
        C_2=C_2,
        G=G,
        A=A,
        A_c=A_c,
        A_s=A_s,
        ustar=ustar,
        r_aa=r_aa,
        r_ac=r_ac,
        r_as=r_as,
        r_sc=r_sc,
        β=β,
        r_ss=r_ss,
        λE_tot=λE_tot,
        λE_tot_p=λE_tot_p,
        VPD_m=VPD_m,
        E_t=E_t,
        λE_t=λE_t,
        E_i=E_i,
        λE_i=λE_i,
        E_s=E_s,
        λE_s=λE_s,
        D_c=D_c,
        P_s=P_s,
        Q_s=Q_s,
        D_1=D_1,
        K_2=K_2,
        I_s=I_s,
    )
end

function conservation_equations(u, p, t)
    fluxes_out = calculate_fluxes(u, p, t)
    # Unpack the static parameters needed
    @unpack d_1, d_2 = p
    # Unpack the fluxes
    @unpack D_c, I_s, D_1, K_2, E_s, E_t, E_i, w_rmax = fluxes_out
    dw1dt = C_1 / (ρ_w * d_1) * (I_s - E_s) - D_1
    dw2dt = 1 / (ρ_w * d_2) * (I_s - E_s - E_t) - K_2
    dwrdt = f_veg * P(t) - E_i - D_c
    return SA[dw1dt, dw2dt, dwrdt]
end

function conservation_equations!(du, u, p, t)
    fluxes_out = calculate_fluxes(u, p, t)
    @unpack d_1 = p
    # Unpack the fluxes
    # Define smoothing parameter for canopy water content
    m_can = w_rmax / 100
    @unpack D_c, I_s, D_1, K_2, E_s, E_t, E_i, w_rmax = fluxes_out
    du[1] = C_1 / (ρ_w * d_1) * (I_s - E_s) - D_1
    du[2] = 1 / (ρ_w * d_2) * (I_s - E_s - E_t) - K_2
    du[3] =
        f_veg * P(t) * smoothing_kernel(UpperBound(), w_r, w_rmax, m_can) -
        E_i * smoothing_kernel(LowerBound(), w_r, zero(w_r), m_can) - D_c
    return nothing
end

# Test the RHS of the ODE system function
# Define the initial conditions
u0 = SA[
    dict_soil[:w_fc] * 1 / 3, # w_1
    dict_soil[:w_fc] * 1 / 3, # w_2
    FT(0.0001),#FT(0.2) * FT(ds_ec_sel.LAI[i]) * 1 / 6, # w_r
]
# Definfe test time stamp
test_timestamp = ds_ec_sel.time[1] + Hour(12)
test_timestamp_unix = datetime2unix(test_timestamp)
# Compute RHS of the ODE system function
test_fluxes_output = calculate_fluxes(u0, param, test_timestamp_unix)
conservation_equations(u0, param, test_timestamp_unix)
#@benchmark conservation_equations($u0, $param, $test_timestamp_unix)

# Solve the ODE system
t_span = (t_unix[1], t_unix[end])
prob = ODEProblem(conservation_equations, u0, t_span, param)
dt = Dates.value(Second(Minute(30)))

# Hand made explicit euler
u_euler = zeros(length((t_unix)), length(u0))
fluxes_euler = Vector(undef, length(t_unix))
for j in 1:length(t_unix)
    if j == 1
        u_euler[j, :] = u0
        fluxes_euler[j] = calculate_fluxes(u0, param, t_unix[j])
    else
        # if t_unix[j] ≥ 1.268106722334931e9
        #     println("Problematic time step: ", t_unix[j])
        # end
        # Correct states below zero to be ≥ 0
        u_euler[j - 1, :] = max.(u_euler[j - 1, :], 0)
        w_rmax = FT(0.2) * FT(ds_ec_sel.LAI[i])
        u_euler[j - 1, 3] = min(u_euler[j - 1, 3], w_rmax)
        # Get the fluxes at the previous time step
        fluxes = calculate_fluxes(u_euler[j - 1, :], param, t_unix[j - 1])
        fluxes_euler[j] = fluxes
        # Update the state variables using explicit Euler method
        state_derivaties = conservation_equations(u_euler[j - 1, :], param, t_unix[j - 1])
        u_euler[j, :] = u_euler[j - 1, :] + dt * state_derivaties
    end
end
λE_tot_explicit = [fluxes_euler[i].λE_tot for i in 1:length(t_unix)]

# The default "Hydrology solver"
sol_explicit_euler = solve(prob, Euler(); dt=dt / 10)
plot(sol_explicit_euler)

# integrator interface
integrator = init(prob, Euler(); dt=dt / 10)
step!(integrator)

# Classic solver: ImplicitEuler fixed timestep at resolution of data
sol_implicit_euler = solve(
    prob, ImplicitEuler(; autodiff=AutoForwardDiff()); adaptive=false, dt=dt / 100
)

# Implicit Euler with adapative timestepping
sol_implicit_euler_adaptive = solve(
    prob,
    ImplicitEuler(; autodiff=AutoForwardDiff());
    adaptive=true,
    saveat=t_unix,
    #callback=PositiveDomain(),
    abstol=1e-5,
    reltol=1e-5,
    #isoutofdomain=(u, p, t) -> any(x -> x < 0, u),
)
plot(sol_implicit_euler_adaptive)

# Compare static array versions with in-place versions
@benchmark solve(
    prob,
    ImplicitEuler(; autodiff=AutoForwardDiff());
    adaptive=true,
    saveat=t_unix,
    #callback=PositiveDomain(),
    abstol=1e-5,
    reltol=1e-5,
) # gives around 150 ms
u0_inplace = [
    dict_soil[:w_fc] * 1 / 3, # w_1
    dict_soil[:w_fc] * 1 / 3, # w_2
    FT(0.0001),
]
prob_inplace = ODEProblem(conservation_equations!, u0_inplace, t_span, param)
solve(
    prob_inplace,
    ImplicitEuler(; autodiff=AutoFiniteDiff());
    adaptive=true,
    saveat=t_unix,
    #callback=PositiveDomain(),
    abstol=1e-5,
    reltol=1e-5,
) # gives around 60 ms -> faster option!

# A recommende solver for stiff problems: Rosenbrock23()
# sol_2 = solve(prob, Heun(); callback=PositiveDomain(), abstol=1e-6, reltol=1e-5
# reltol and abstol needed if ConstantInterpolation, not for PCHIPInterpolation,
# Also for ConstantInterpolation better to work at higerh abstol values
sol_rosenbrock = solve(
    prob, Rosenbrock23(; autodiff=AutoFiniteDiff()); saveat=t_unix, abstol=1e-5, reltol=1e-5
)
plot(sol_rosenbrock)

# Hinting that stiff
# Again, give tolerances if using ConstantInterpolation
# sol_alg_hints = solve(
#     prob; alg_hints=:stif, abstol=1e-5, reltol=1e-5
# )
# plot(sol_alg_hints)
# possible_algs = sol_alg_hints.alg.algs
# choice_stiff = unique(sol_alg_hints.alg_choice)
# for i in choice_stiff
#     println("Stiff algorithm choice: ", possible_algs[i])
# end

# No hints, full auto
# sol_auto = solve(prob, saveat=t_unix, reltol=1e-5)
# plot(sol_auto)
# possible_algs = sol_auto.alg.algs
# choice_auto = unique(sol_auto.alg_choice)
# for i in choice_auto
#     println("Auto algorithm choice: ", possible_algs[i])
# end
# Note: interpolation DOES have an effect!
# With PCHIP -> auto solving works
# With ConstantInterpolation -> auto solving fails

## do constant interpolation BUT force the solver to STOP at the points where there are discontinuities
R_n = ConstantInterpolation(FT.(ds_ec_sel.Rnet[:]), t_unix; dir=:left)
T_a = ConstantInterpolation(FT.(ds_ec_sel.Tair[:]), t_unix; dir=:left)
p_a = ConstantInterpolation(FT.(ds_ec_sel.Psurf[:]), t_unix; dir=:left)
u_a = ConstantInterpolation(FT.(ds_ec_sel.Wind[:]), t_unix; dir=:left)
VPD_a = ConstantInterpolation(FT.(ds_ec_sel.VPD[:]) .* 100, t_unix; dir=:left) #hPa -> Pa
P = ConstantInterpolation(FT.(ds_ec_sel.Precip[:]), t_unix; dir=:left)
SW_in = ConstantInterpolation(FT.(ds_ec_sel.SWdown[:]), t_unix; dir=:left)
LAI = ConstantInterpolation(FT.(ds_ec_sel.LAI[:]), t_unix; dir=:left)

# e.g. using recommend explicit Heun method from La Follette et al. 2021
sol_tstops = solve(prob, Heun(); tstops=t_unix)
plot(sol_tstops)

# The default non-stiff stolver
sol_tstops_tsit5 = solve(prob, Tsit5(); tstops=t_unix)

# even give explicit euler a shot!
# condition(u, t, integrator) = any(u .≤ 0) # Positive domain condition
# affect!(integrator) = integrator.u[integrator.u .≤ 0] .= 0
u0_bis = [
    dict_soil[:w_fc] * 1 / 3, # w_1
    dict_soil[:w_fc] * 1 / 3, # w_2
    FT(0.0001),#FT(0.2) * FT(ds_ec_sel.LAI[i]) * 1 / 6, # w_r
]
prob_bis = ODEProblem(conservation_equations, u0_bis, t_span, param)
condition(u, t, integrator) = true#(integrator.u[3] ≤ 0.0) | (u[3] ≤ 0.0)
affect!(integrator) = integrator.u = max.(0, integrator.u) # Set w_r to 0 if it goes below 0
cb = DiscreteCallback(condition, affect!)
sol_ee_pos = solve(prob_bis, Euler(); dt=dt / 10, callback=cb)#, tstops = t_unix)
plot(sol_ee_pos)

## Experiment
## Key to do: Implement a SavingCallback to save the fluxes!
# https://docs.sciml.ai/DiffEqCallbacks/stable/output_saving/#DiffEqCallbacks.SavingCallback
# Tutorial: https://nextjournal.com/sosiris-de/ode-diffeq
# IDEA: try using saveat, does this make a difference?
function save_flux_data(u, t, integrator)
    # Reuse the same calculations function
    fluxes = calculate_fluxes(u, integrator.p, t)

    # Return the specific values you want to save
    return (
        λE_t=fluxes.λE_t,
        λE_i=fluxes.λE_i,
        λE_s=fluxes.λE_s,
        λE_tot=fluxes.λE_tot,
        E_t=fluxes.E_t,
        E_i=fluxes.E_i,
        E_s=fluxes.E_s,
        f_veg=fluxes.f_veg,
        f_wet=fluxes.f_wet,
        # Add any other values you'd like to record
    )
end
save_cb = SavingCallback(save_flux_data, saved_flux_data; saveat=t_unix)
sol_cb_save = solve(
    prob,
    ImplicitEuler(; autodiff=AutoFiniteDiff());
    adaptive=true,
    saveat=t_unix,
    #callback=PositiveDomain(),
    abstol=1e-5,
    reltol=1e-5,
    callback=save_cb,
)
plosaved_flux_data = SavedValues(FT, NamedTuple)
t(sol_cb_save)
λE_tot_saved = [data.λE_tot for data in saved_flux_data.saveval]
λE_t_saved = [data.λE_t for data in saved_flux_data.saveval]

## PLOTS FOR EGU
# Plot the states
gr()
figdir(args...) = projectdir("figures", args...)
if ~isdir(figdir())
    mkdir(figdir())
end
cm = 37.8 #1cm = 37.8 px
mm = cm / 10
start_date = unix2datetime(t_unix[50])
end_date = unix2datetime(t_unix[end - 50])
date_ticks = [start_date, end_date]

# Figure 1: states for stable adapative implicit euler
plot(
    unix2datetime.(sol_cb_save.t),
    Array(sol_cb_save)';
    label=["w₁ [-]" "w₂ [-]" "wᵣ [kg/m²]"],
    #ylabel="Soil moisture [-]/ vegetation water content[kg/m²]",
    xlabel="Time",
    ylims=(0, maximum(sol_cb_save[1, :]) * 1.5),
    title="Adaptive Implicit Euler (IEₐ)",
    xticks=(date_ticks, Dates.format.(date_ticks, "yyyy-mm-dd")),
)
#p2 = twinx(p1)
precip_plot_mm_h = collect(P(t_unix)[:]) * 3600 #kg/(m²s) -> mm/h
y_ticks_precip = ([0, 10, 20], ["0", "10", "20"])
fig_implicit = plot!(
    twinx(),
    unix2datetime.(t_unix),
    precip_plot_mm_h;
    # linetype=:bar,
    fill=(0, :gray),
    color=:gray,
    yflip=true,
    ylabel="Precipitation [mm/h]",
    legend=:none,
    xticks=(date_ticks, Dates.format.(date_ticks, "yyyy-mm-dd")),
    yticks=y_ticks_precip,
    ylims=(0, maximum(precip_plot_mm_h) * 2.5),
)
savefig(fig_implicit, figdir("IE_a_states.png"))

# Figure 2: states for explicit euler without corrections
fig_explicit = plot(
    unix2datetime.(t_unix),
    u_euler;
    label=["w₁ [-]" "w₂ [-]" "wᵣ [kg/m²]"],
    #ylabel="Soil moisture [-]/ vegetation water content[kg/m²]",
    xlabel="Time",
    ylims=ylims(fig_implicit),
    title="Explicit Euler (EE)",
)
plot!(
    twinx(),
    unix2datetime.(t_unix),
    precip_plot_mm_h;
    # linetype=:bar,
    fill=(0, :gray),
    color=:gray,
    yflip=true,
    ylabel="Precipitation [mm/h]",
    legend=:none,
    xticks=(date_ticks, Dates.format.(date_ticks, "yyyy-mm-dd")),
    yticks=y_ticks_precip,
    ylims=(0, maximum(precip_plot_mm_h) * 2.5),
)

savefig(fig_explicit, figdir("EE_states.png"))

# Figure 3: fluxes for both methods
λE_obs = collect(ds_ec_sel.Qle_cor[:])
corr = cor(λE_tot_saved, λE_obs)
fig_implicit_fluxes = plot(
    collect(ds_ec_sel.time),
    ds_ec_sel.Qle_cor[:];
    label="Observation",
    color=:black,
    #ylabel=L"\lambda E \quad [\mathrm{W}/ \mathrm{m}^2]",
    ylabel="λE [W/m²]",
    xlabel="Time",
    # title="Nothing",
    # titlefontcolor=:white,
    xticks=(date_ticks, Dates.format.(date_ticks, "yyyy-mm-dd")),
    framestyle=:box,
)
plot!(
    unix2datetime.(t_unix),
    λE_tot_explicit;
    label="EE",
    color=collect(palette(:default))[1],#[6]
)
plot!(
    unix2datetime.(saved_flux_data.t),
    λE_tot_saved;
    label="IEₐ",
    color=collect(palette(:default))[2],#[5],
)
#plot!(unix2datetime.(saved_flux_data.t), λE_t_saved; label="λE_t")

savefig(fig_implicit_fluxes, figdir("IE_a_fluxes.png"))

# Figure 4: Compare fluxes for explicit and implicit euler
λE_tot_diff = λE_tot_saved - λE_tot_explicit
fig_fluxes_diff = plot(
    unix2datetime.(saved_flux_data.t),
    λE_tot_diff;
    label=:none,
    ylabel="λE(IEₐ) - λE(EE) [W/m²]",
    xlabel="Time",
    xticks=(date_ticks, Dates.format.(date_ticks, "yyyy-mm-dd")),
    framestyle=:box,
    color=palette(:default)[3],
)
savefig(fig_fluxes_diff, figdir("IEa_diff_EE_fluxes.png"))

# Combine the plots in one
l = @layout [_ a{0.485w} b{0.485w} _; _ c{0.485w} d{0.485w} _]
title_fontsize = 18
tick_fontsize = title_fontsize - 4
fig_implicit_combine = plot(
    fig_implicit;
    xlabel="",
    xtickfontcolor=:white,
    xtickfontsize=1,
    ytickfontsize=tick_fontsize,
)
fig_explicit_combine = plot(
    fig_explicit;
    xlabel="",
    xtickfontcolor=:white,
    xtickfontsize=1,
    ytickfontsize=tick_fontsize,
    legend=false,
)
fig_implicit_fluxes_combine = plot(fig_implicit_fluxes; tickfontsize=tick_fontsize)
fig_fluxes_diff_combine = plot(fig_fluxes_diff; tickfontsize=tick_fontsize)
fig_combined = plot(
    fig_implicit_combine,
    fig_explicit_combine,
    fig_implicit_fluxes_combine,
    fig_fluxes_diff_combine;
    layout=l,
    size=(30cm, 29cm),
    link=:x,
    legendfontsize=12,
    titlefontsize=title_fontsize,
    guidefontsize=title_fontsize - 4,
    line=2.5,
)
savefig(fig_combined, figdir("combined.png"))
savefig(fig_combined, figdir("combined.svg"))
fig_combined_transparent = plot(fig_combined; background_color=:transparent)
savefig(fig_combined_transparent, figdir("combined_transparent.svg"))

# Experiment Two: modellingtoolkit
