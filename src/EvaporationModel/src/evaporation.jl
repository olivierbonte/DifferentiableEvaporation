"""
    penman_monteith(t_a, p_a, r_net, vpd, r_a, r_s; g = 0.0, kwargs...)

Compute evaporation (ET) and latent heat flux (LE)  using the Penman-Monteith equation.

# Arguments
- `t_a`: Air temperature [K].
- `p_a`: Surface pressure [Pa].
- `r_net`: Net radiation [W/m²].
- `vpd`: Vapor pressure deficit [Pa].
- `r_a`: Aerodynamic resistance [s/m].
- `r_s`: Surface resistance [s/m].
- `g`: Soil heat flux [W/m²]. Default is 0.0.
- `kwargs`: Additional keyword arguments.

# Returns
- `et`: Potential evapotranspiration [kg/(m² * s)].
- `λe`: Latent heat flux [W/m²].

# See also
[`Bigleaf.potential_ET`](https://earthyscience.github.io/Bigleaf.jl/dev/evapotranspiration/#Bigleaf.potential_ET)
 for more details on the Penman-Monteith equation.
"""
function penman_monteith(
    t_a::T, p_a::T, r_net::T, vpd::T, r_a::T, r_s::T; g=0, kwargs...
) where {T}
    con = Bigleaf.BigleafConstants()
    et, λe = Bigleaf.potential_ET(
        PenmanMonteith(),
        t_a - T(con.Kelvin),
        p_a * T(con.Pa2kPa),
        r_net,
        vpd * T(con.Pa2kPa),
        1 / r_a;
        G=g,
        Gs_pot=Bigleaf.ms_to_mol(1 / r_s, t_a - con.Kelvin, p_a * con.Pa2kPa),
    )
    return et, λe
end

function total_evaporation(
    t_a::T,
    p_a::T,
    vpd::T,
    A::T,
    A_c::T,
    A_s::T,
    r_aa::T,
    r_ac::T,
    r_as::T,
    r_sc::T,
    r_ss::T,
    f_wet::T,
) where {T}
    con = Bigleaf.BigleafConstants()
    Δ = Bigleaf.Esat_from_Tair_deriv(t_a - T(con.Kelvin))
    γ = Bigleaf.psychrometric_constant(t_a - T(con.Kelvin), p_a * T(con.Pa2kPa))
    R_c = r_sc + (1 + Δ / γ) * r_ac
    R_s = r_ss + (1 + Δ / γ) * r_as
    R_a = (1 + Δ / γ) * r_aa
    R_i = (1 + Δ / γ) * r_ac
    DE = R_c * R_i * R_s + R_a * ((1 - f_wet) * R_i * R_s + f_wet * R_c * R_s + R_c * R_i)
    P_c = r_aa * (1 - f_wet) * R_i * R_s / DE
    P_i = r_aa * f_wet * R_c * R_s / DE
    P_s = r_aa * R_c * R_i / DE
    et_p, λe_p = penman_monteith(t_a, p_a, A, vpd, r_aa, T(0)) # r_s = 0 -> Penman
    λe =
        (Δ + γ) / γ * (P_c + P_i + P_s) * λe_p +
        Δ / (γ * r_aa) * (P_c * A_c * r_ac + P_i * A_c * r_ac + P_s * A_s * r_as)
    return λe, λe_p
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
