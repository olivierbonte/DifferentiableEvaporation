abstract type SurfaceResistanceMethod end
struct JarvisStewart <: SurfaceResistanceMethod end

abstract type SoilAerodynamicResistanceMethod end
struct Choudhury1988soil <: SoilAerodynamicResistanceMethod end

abstract type SoilEvaporationEfficiencyMethod end
struct Martens17 <: SoilEvaporationEfficiencyMethod end
struct Pielke92 <: SoilEvaporationEfficiencyMethod end

"""
    surface_resistance(::JarvisStewart, SW_in, VPD, T_a, w_2, w_fc, w_wilt, LAI, g_d, r_smin)

Calculate surface resistance

# Arguments
- `approach`: The approach to use for calculating surface resistance
- `SW_in`: Incoming solar radiation [W/m²]
- `VPD`: Vapor pressure deficit [Pa]
- `T_a`: Air temperature [K]
- `w_2`: Root zone soil moisture [m³/m³]
- `w_fc`: Soil moisture at field capacity [m³/m³]
- `w_wilt`: Soil moisture at wilting point [m³/m³]
- `LAI`: Leaf area index [m²/m²]
- `g_d`: Parameter relating `VPD` to surface conductance [Pa⁻¹]
- `r_smin`: Minimum surface resistance [s/m]
- `T_opt`: Optimum temperature for stomatal conductance [K], default = 298.0 K
- `r_smax`: Maximum surface resistance [s/m], default = 500_000 s/m

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
    SW_in::T,
    VPD::T,
    T_a::T,
    w_2::T,
    w_fc::T,
    w_wilt::T,
    LAI::T,
    g_d::T,
    r_smin::T;
    T_opt::T=T(298.0),
    r_smax::T=T(500_000),
) where {T}
    f_1 = clamp((T(0.004) * SW_in + T(0.05)) / (T(0.81) * (T(0.004) * SW_in + T(1.0))), 0, 1)
    f_2 = clamp((w_2 - w_wilt) / (w_fc - w_wilt), 0, 1)
    f_3 = clamp(exp(-g_d * VPD), 0, 1)
    f_4 = clamp(1 - T(0.0016) * (T_opt - T_a)^2, 0, 1)
    r_s = min(r_smin / LAI * (f_1 * f_2 * f_3 * f_4)^(-1), r_smax)
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
- r_as: Aerodynamic resistance between soil and canopy source height [s/m]

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
    @unpack SW_in, w_2, VPD, T_a, LAI = forcing
    @unpack w_fc, w_wilt, gd, r_smin = parameters
    return jarvis_stewart.(SW_in, w_2, VPD, T_a, LAI, w_fc, w_wilt, gd, r_smin)
end

function jarvis_stewart(SW_in, w_2, VPD, T_a, LAI, w_fc, w_wilt, gd, r_smin)
    f_1 = min(1.0, (0.004 * SW_in + 0.05) / 0.81(0.004 *  SW_in + 1.0))
    f_2 = max(0.0, min((w_2 - w_wilt) / (w_fc - w_wilt), 1.0))
    f_3 = exp(-gd * VPD)
    #f_4 = @. max(1.0 - 0.0016*(298.0 - T_a)^2, 0.0)
    r_s = r_smin ./ LAI * (f_1 * f_2 * f_3)^(-1)
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

@doc raw"""
    ustar_from_u(u, z_obs, d, z_0m, ψ_m=0)

Calculate the friction velocity from the wind speed at a given height assuming a logarithmic wind profile.

# Arguments
- `u`: Wind speed at the measurement height [m/s]
- `z_obs`: Height of the wind speed measurement [m]
- `d`: Displacement height [m]
- `z_0m`: Roughness length for momentum transfer [m]
- `ψ_m`: Stability correction for momentum transfer [m], default = 0

# Returns
- `u_star`: Friction velocity [m/s]

# Details
The friction velocity is calculated using the logarithmic wind profile equation, given by:

``u(z) = \frac{u^*}{k} \ln \left( \frac{z - d}{z_{0m}} \right) - \psi_m``

For more info on how to calculate ``\psi_m``, see the
[Bigleaf package documentation](https://earthyscience.github.io/Bigleaf.jl/dev/stability_correction/#Bigleaf.stability_parameter)
"""
function ustar_from_u(u::T, z_obs::T, d::T, z_0m::T, ψ_m::T=T(0)) where {T}
    return T(BigleafConstants().k) * (u + ψ_m) / log((z_obs - d) / z_0m)
end

"""
    soil_evaporation_efficiency(approach::Pielke92, w_1, w_fc)
    soil_evaporation_efficiency(approach::Martens17, w_1, w_res, w_c)

Calculate soil evaporation efficiency β, a factor between 0 and 1 which scales the
potential soil evaporation

# Arguments
- `approach`: The approach to use for calculating soil evaporation efficiency, a subtype of
    `SoilEvaporationEfficiencyMethod`
- `w_1`: Soil moisture from first layer [m³/m³]

With `approach = Pielke92()`, the following inputs are

- `w_fc`: Soil moisture at field capacity [m³/m³]

With `approach = Martens17()`, the following inputs are

- `w_res`: Residual soil moisture [m³/m³]
- `w_c`: Critical soil moisture [m³/m³]

# Details

With `approach = Pielke92()`, β is calculated using euqation (7) of
[Lee & Pielke (1992)](https://doi.org/10.1175/1520-0450(1992)031%3C0480:ETSSSH%3E2.0.CO;2)

With `approach = Martens17()`, β is calculated using equation (6) of
[Martens et al. (2017)](https://doi.org/10.5194/gmd-10-1903-2017).

# Examples

```jldoctest
using EvaporationModel
w_fc = 0.35
w_1 = w_fc / 2
β_p = soil_evaporation_efficiency(Pielke92(), w_1, w_fc)
β_p ≈ 0.25

# output

true
```

"""
function soil_evaporation_efficiency(approach::Pielke92, w_1::T, w_fc::T) where {T}
    return 1 / 4 * (1 - cos(π * w_1 / w_fc))^2
end

function soil_evaporation_efficiency(
    approach::Martens17, w_1::T, w_res::T, w_c::T
) where {T}
    return clamp((w_1 - w_res) / (w_c - w_res), 0, 1)
end

"""
    beta_to_r_ss(beta, r_as)

Calculate soil surface resistance ``r_{ss}`` from soil evaporation efficiency β

# Arguments
- `beta`: Soil evaporation efficiency [-]
- `r_as`: Aerodynamic resistance between soil and canopy source height [s/m]

# Details

Equation used derived from equivalency between equation (6) and (7) of
[Merlin et al. (2016)](https://doi.org/10.1002/2015WR018233)

"""
function beta_to_r_ss(beta::T, r_as::T) where {T}
    return r_as / beta - r_as
end

"""
    r_ss_to_beta(r_ss, r_as)

Calculate soil evaporation efficiency β soil surface resistance ``r_{ss}``

# Arguments
- `r_ss`: Soil surface resistance [s/m]
- `r_as`: Aerodynamic resistance between soil and canopy source height [s/m]

# Details

Equation used derived from equivalency between equation (6) and (7) of
[Merlin et al. (2016)](https://doi.org/10.1002/2015WR018233)

"""
function r_ss_to_beta(r_ss::T, r_as::T) where {T}
    return r_as / (r_as + r_ss)
end

# r_a,a: https://earthyscience.github.io/Bigleaf.jl/dev/aerodynamic_conductance/#Bigleaf.compute_Ram

# r_a,c: https://earthyscience.github.io/Bigleaf.jl/dev/boundary_layer_conductance/#Bigleaf.Gb_constant_kB1
