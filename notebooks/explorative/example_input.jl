using DrWatson
@quickactivate "DifferentiableEvaporation"
using ComponentArrays

# Fixed/synthetic forcings as an easy start
function P_test(t::T) where {T}
    return T(5e-6)
end
function T_a_test(t::T) where {T}
    return T(275.0) #+ T(5.0) * sin(T(2) * π * t / T(86400) - T(π / 2))
end
function u_a_test(t::T) where {T}
    return T(3.0)
end
function p_a_test(t::T) where {T}
    return T(101325.0)
end
function VPD_a_test(t::T) where {T}
    return T(200.0)
end
function SW_in_test(t::T) where {T}
    return T(1000) #max(T(1361.0) * sin(T(2) * π * t / T(86400) - π / T(2)), zero(t))
end
function R_n_test(t::T) where {T}
    return T(SW_in_test(t)) * T(0.3)
end
function LAI_test(t::T) where {T}
    return T(3.0)
end

forcings = (
    P=P_test,
    T_a=T_a_test,
    u_a=u_a_test,
    p_a=p_a_test,
    VPD_a=VPD_a_test,
    SW_in=SW_in_test,
    R_n=R_n_test,
    LAI=LAI_test,
)

# Model parameters
param_test = ComponentArray(;
    h=21.0, # [m] height of the canopy
    z_0ms=0.01, # [m] roughness length for soil
    w_sat=0.5, # [-]  saturation water content
    a=0.15,
    p_soil=6.0,
    b=6.1,
    C_1sat=1.9,
    C_2ref=0.83,
    C_3=0.25,
    d_1=0.01,
    d_2=1.3,
    w_res=0.04,
    w_wp=0.08,
    w_fc=0.3,
    z_obs=39.0, # [m]
    kB⁻¹=log(10),
    g_d=0.0003,
    r_smin=395.0,
)
FT = Float64

# Simulation settings
u0_test = [0.2, 0.2, 0.001]
t_end = 1 * 86400.0 # 1 day in seconds
t_span_test = (0.0, t_end) # 1 day in seconds
dt = 1800 #s, 30"
saveat_test = 0:dt:t_end
