"""
    penman_monteith(t_air, p_surf, r_net, vpd, r_a, r_s; g = 0.0, kwargs...)

Compute evaporation (ET) and latent heat flux (LE)  using the Penman-Monteith equation.

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
- `et`: Potential evapotranspiration [kg/(m² * s)].
- `le`: Latent heat flux [W/m²].

# See also
[`Bigleaf.potential_ET`](https://earthyscience.github.io/Bigleaf.jl/dev/evapotranspiration/#Bigleaf.potential_ET)
 for more details on the Penman-Monteith equation.
"""
function penman_monteith(t_air::T, p_surf::T, r_net::T, vpd::T, r_a::T, r_s::T; g=0, kwargs...) where {T}
    con = Bigleaf.BigleafConstants()
    et, le = Bigleaf.potential_ET(
        PenmanMonteith(),
        t_air - T(con.Kelvin),
        p_surf * T(con.Pa2kPa),
        r_net,
        vpd * T(con.Pa2kPa),
        1 / r_a;
        G=g,
        Gs_pot=Bigleaf.ms_to_mol(1 / r_s, t_air - con.Kelvin, p_surf * con.Pa2kPa),
    )
    return et, le
end

function compute_soil_evaporation_stress(w_g, w_crit, w_res)
    return max(1.0, (w_g - w_res) / (w_crit - w_res))
end

function compute_bare_soil_evaporation(et_pot, f_veg, stress)
    return (1.0 - f_veg) * stress * et_pot
end

function compute_transpiration(et, f_veg, f_wet)
    return f_veg * (1 - f_wet) * et
end
