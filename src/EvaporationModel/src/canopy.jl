function fractional_vegetation_cover(LAI, k_ext=oftype(LAI,0.5))
    return 1 - exp(-k_ext * LAI)
end

function available_energy_partioning(R_n, G, f_veg)
    A = R_n - G
    R_nc = R_n * f_veg
    R_ns = R_n * (1 - f_veg)
    A_c = R_nc
    A_s = R_ns - G
    return A, A_c, A_s
end

function max_canopy_capacity(LAI, c=oftype(LAI, 0.2))
    return c * LAI
end

function fraction_wet_vegetation(w_r, w_rmax)
    w_r = smooth_max(w_r, zero(w_r), w_rmax / 1000)
    return (w_r / w_rmax)^oftype(w_r, 2 / 3)
end

function canopy_drainage(P, w_r, f_veg, c=oftype(f_veg, 0.2))
    k = (1 - f_veg) * P
    b = 1 / (2 * c)
    D_c = k * (exp(b * w_r) - 1)
    return D_c
end

function precip_below_canopy(P, f_veg, D_c)
    P_s = P * (1 - f_veg) + D_c
    return P_s
end

function vpd_veg_source_height(VPD_a, T_a, p_a, A, λE, r_aa)
    T = typeof(r_aa)
    con = Bigleaf.BigleafConstants()
    Δ = Bigleaf.Esat_from_Tair_deriv(T_a - T(con.Kelvin)) * T(con.kPa2Pa)
    γ =
        Bigleaf.psychrometric_constant(T_a - T(con.Kelvin), p_a * T(con.Pa2kPa)) *
        T(con.kPa2Pa)
    ρ_a = Bigleaf.air_density(T_a - T(con.Kelvin), p_a * T(con.Pa2kPa))
    c_p = T(Bigleaf.BigleafConstants().cp)
    VPD_m = VPD_a + (Δ * A - (Δ + γ) * λE) * r_aa / (ρ_a * c_p)
    return VPD_m
end
