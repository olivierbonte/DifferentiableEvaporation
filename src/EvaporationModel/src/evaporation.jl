import Bigleaf

"""
    penman_monteith(t_air, p_surf, r_net, vpd, r_a, r_s; g = 0.0, kwargs...)

Compute evaporation (ET) and latent heat flux (LE) 
using the Penman-Monteith equation.

# Arguments
- `t_air`: Air temperature [K].
- `p_surf`: Surface pressure [Pa].
- `r_net`: Net radiation [W/m²].
- `vpd`: Vapor pressure deficit [Pa].
- `r_a`: Aerodynamic resistance [s/m].
- `r_s`: Surface resistance [s/m].
- `g`: Soil heat flux [W/m²]. Default is 0.0.
- `kwargs`: Additional keyword arguments.

# Returns
- `et`: Potential evapotranspiration [kg/(m^2 * s)].
- `le`: Latent heat flux [W/m²].
"""
function penman_monteith(t_air, p_surf, r_net, vpd, r_a, r_s; g = 0.0, kwargs...)
    cp = Bigleaf.bigleaf_constants()
    et, le = Bigleaf.potential_ET(
        Val(:PenmanMonteith), t_air - cp[:Kelvin], p_surf * cp[:Pa2kPa], r_net, 
        vpd * cp[:Pa2kPa], 1.0 / r_a, G = g, 
        Gs_pot = Bigleaf.ms_to_mol(1.0 / r_s, t_air - cp[:Kelvin], p_surf * cp[:Pa2kPa])
    )
    return et, le 
end

"""
    compute_g_from_r_net(r_net, lai)

Compute the ground heat flux [W/m²] from net radiation (and LAI)

# Arguments
- `r_net`: The net radiation [W/m²].
- `lai`: The leaf area index [m²/m²].

# Details
This implementation follows the approach of the METRIC model.
See Equation 27 of 
[Allen et al., 2007](https://doi.org/10.1061/(ASCE)0733-9437(2007)133:4(380))

"""
function compute_g_from_r_net(r_net, lai)
    if lai < 0.5
        @error "Ground heat flux from net radiation for LAI < 0.5
        not yet implemented"
    else
        g = r_net * (0.05 + 0.18 * exp(-0.521 * lai))
    end
    return g
end