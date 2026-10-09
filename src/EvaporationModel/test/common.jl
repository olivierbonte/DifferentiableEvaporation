# Shared test inputs, included once by runtests.jl
println("Defining test inputs")
r_aa = 50.0 # s/m
r_ac = 30.0 # s/m
r_as = 100.0 # s/m
r_sc = 1000.0 # s/m
r_ss = 100.0 # s/m
f_wet = 0.1
f_veg = 0.5
SW_in = 1000.0 # W/m2
w_2 = 0.2
w_fc = 0.3
w_wp = 0.1
g_d = 0.003
r_smin = 200.0 # s/m
Rn = 400.0 # W/m2
G = 50.0 # W/m2
T_a = 300.0 # K
p_a = 101325.0 # Pa
VPD_a = 2000.0 # Pa
LAI = 3.0 # m2/m2
h = 21.0 # m
d_c = 2 / 3 * h # m
z_0mc = 0.1 * h # m
z_obs = 39.0 # m
z_0ms = 0.01 # m
u_star = 3.0 # m/s
u = 5.0 # m/s
P = 2.0e-5 # gross precipitation, kg / (m2 * s)
P_s = 1.5e-5 # precipitation below the canopy, kg / (m2 * s)
w_sat = 0.45

p_model = ComponentArray(;
    h=h, z_0ms=z_0ms, w_sat=w_sat, a=0.15, p_soil=6.0, b=6.1, w_res=0.04,
    w_wp=w_wp, w_fc=w_fc, C_1sat=0.019, C_2ref=0.83, C_3=0.25, d_1=0.01, d_2=1.3,
    z_obs=z_obs, kB⁻¹=log(10), g_d=3e-4, r_smin=395.0, k_ext=0.5,
)
# `Returns` captures the values: closures over (non-const) globals give silently wrong
# Enzyme derivatives unless runtime activity is enabled
constant_forcings(P) = (
    P=Returns(P), T_a=Returns(T_a), u_a=Returns(u), p_a=Returns(p_a), VPD_a=Returns(VPD_a),
    SW_in=Returns(SW_in), R_n=Returns(Rn), LAI=Returns(LAI), lon=4.52, utc_offset=1,
)
t_unix = datetime2unix(DateTime(2010, 7, 1, 10)) # model time [s since Unix epoch], 10:00 local
treatments = (HardThresholds(), KavetskiSmoothing())

# Toy case: one day of constant forcing with light rain, output every 30 minutes.
# w_1 < w_fc, where the Jacobian of the hard thresholds is finite (r_ss = 0 above w_fc).
u0_toy = [0.2, 0.2, 0.001]
t_span_toy = (t_unix, t_unix + 86400.0)
saveat_toy = collect(t_span_toy[1]:1800.0:t_span_toy[2])
forcings_toy = constant_forcings(5e-6) # kg/(m² * s)

function toy_model(; thresholds=HardThresholds(), parameters=p_model, u0=u0_toy)
    model = ProcessBasedModel{Float64}(;
        forcings=forcings_toy,
        parameters=parameters,
        t_span=t_span_toy,
        u0=copy(u0),
        saveat=saveat_toy,
        thresholds=thresholds,
    )
    initialize!(model)
    return model
end
