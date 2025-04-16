abstract type SurfaceResistanceMethod end
struct JarvisStewart <: SurfaceResistanceMethod end

abstract type SoilAerodynamicResistanceMethod end
struct Choudhury1988soil <: SoilAerodynamicResistanceMethod end

"""
    surface_resistance(::JarvisStewart, s_in, vpd, t_air, w_2, w_fc, w_wilt, lai, g_d, r_smin)

Calculate surface resistance

# Arguments
- `approach`: The approach to use for calculating surface resistance
- `s_in`: Incoming solar radiation [W/m²]
- `vpd`: Vapor pressure deficit [Pa]
- `t_air`: Air temperature [K]
- `w_2`: Root zone soil moisture [m³/m³]
- `w_fc`: Soil moisture at field capacity [m³/m³]
- `w_wilt`: Soil moisture at wilting point [m³/m³]
- `lai`: Leaf area index [m²/m²]
- `g_d`: Parameter relating `vpd` to stomatal conductance [Pa⁻¹]
- `r_smin`: Minimum stomatal resistance [s/m]

# Returns
- `r_s`: Surface resistance [s/m]

# Details

With `approach = JarvisStewart()`, the Jarvis-Stewart model is used as given in equations 
(8.9), (8.10), (8.11) and (8.14) of the 
[IFS Cy47r1 documentation Part IV: Physical processes](https://doi.org/10.21957/eyrpir4vj). 
The temperature constraint (`f_4`) comes from 
[Noilhan and Planton (1989)](https://doi.org/10.1175/1520-0493(1989)117%3C0536:ASPOLS%3E2.0.CO;2)
equation (37), with the default value of ``T_{opt}`` set to 298.0 K
"""
function surface_resistance(
    approach::JarvisStewart,
    s_in::T,
    vpd::T,
    t_air::T,
    w_2::T,
    w_fc::T,
    w_wilt::T,
    lai::T,
    g_d::T,
    r_smin::T;
    T_opt::T=T(298.0),
) where {T}
    f_1 = clamp((T(0.004) * s_in + T(0.05)) / (T(0.81) * (T(0.004) * s_in + T(1.0))), 0, 1)
    f_2 = clamp((w_2 - w_wilt) / (w_fc - w_wilt), 0, 1)
    f_3 = clamp(exp(-g_d * vpd), 0, 1)
    f_4 = clamp(1 - T(0.0016) * (T_opt - t_air)^2, 0, 1)
    r_s = r_smin / lai * (f_1 * f_2 * f_3 * f_4)^(-1)
    return r_s
end
"""
    soil_aerodynamic_resistance(::Choudhury1988soil, ustar, h, d_c, z_0mc, z_0ms, η=3)

Calculate soil aerodynamic resistance, which is between the soil surface and the 
canopy source height (``z_m = z_{0mc} + d_{c}``)

# Arguments
- approach: The approach to use for calculating soil aerodynamic resistance
- ustar: Friction velocity [m/s]
- h: Canopy height [m]
- d_c: Canopy displacement height [m]
- z_0mc: Roughness length for momentum transfer for canopy [m]
- z_0ms: Roughness length for momentum transfer for soil [m]
- η: Extinction coefficient of ``K`` in canopy, default = 3 [-]

# Returns
- r_as: Soil aerodynamic resistance [s/m]

# Details
With `approach = Choudhury1988soil()`, equation (25) of 
[Choudhury and Monteith (1988)](https://doi.org/10.1002/qj.49711448006) is used.  
"""
function soil_aerodynamic_resistance(
    approach::Choudhury1988soil, ustar::T, h::T, d_c::T, z_0mc::T, z_0ms::T, η::T=T(3)
) where {T}
    Kh = T(BigleafConstants().k) * ustar * (h - d_c)
    r_as = (h * exp(η) / (η * Kh)) * (exp(-η * z_0ms / h) - exp(-η * (z_0mc + d_c) / h))
    return r_as
end

function jarvis_stewart(forcing::ComponentArray, parameters::ComponentArray)
    @unpack s_in, w_2, vpd, t_air, lai = forcing
    @unpack w_fc, w_wilt, gd, r_smin = parameters
    return jarvis_stewart.(s_in, w_2, vpd, t_air, lai, w_fc, w_wilt, gd, r_smin)
end

function jarvis_stewart(s_in, w_2, vpd, t_air, lai, w_fc, w_wilt, gd, r_smin)
    f_1 = min(1.0, (0.004 * s_in + 0.05) / 0.81(0.004 * s_in + 1.0))
    f_2 = max(0.0, min((w_2 - w_wilt) / (w_fc - w_wilt), 1.0))
    f_3 = exp(-gd * vpd)
    #f_4 = @. max(1.0 - 0.0016*(298.0 - t_air)^2, 0.0)
    r_s = r_smin ./ lai * (f_1 * f_2 * f_3)^(-1)
    #max r_s of 50e3 s/max
    r_s = max(50e3, r_s)
    return r_s
end

function aerodynamic_resistance(u::AbstractArray, parameters::ComponentArray)
    @unpack d, z0m, z0h, z_measur = parameters
    return aerodynamic_resistance.(u, d, z0m, z0h, z_measur)
end

function aerodynamic_resistance(u, d, z0m, z0h, z_measur)
    u = max(u, 0.1)
    #I use 0.1 instead of 1 here to allow higher values! https://doi.org/10.1016/S0168-1923(00)00164-7
    r_a = (log((z_measur - d) / z0m) * log((z_measur - d) / z0h)) ./ (0.41^2 * u)
    return r_a
end

function ustar_from_u(u::T, z_obs::T, d::T, z_0m::T, ψ_m::T=0.0) where {T}
    return T(BigleafConstants().k) * (u + ψ_m) / log((z_obs - d) / z_0m)
end

# r_a,a: https://earthyscience.github.io/Bigleaf.jl/dev/aerodynamic_conductance/#Bigleaf.compute_Ram

# r_a,c: https://earthyscience.github.io/Bigleaf.jl/dev/boundary_layer_conductance/#Bigleaf.Gb_constant_kB1
