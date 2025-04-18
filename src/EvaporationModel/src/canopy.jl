function fractional_vegetation_cover(LAI::T, k_ext::T=0.5) where {T}
    f_veg = 1 - exp(-k_ext * LAI)
    return f_veg
end

function available_energy_partioning(R_n::T, G::T, f_veg::T) where {T}
    A = R_n - G
    R_nc = R_n * f_veg
    R_ns = R_n * (1 - f_veg)
    A_c = R_nc
    A_s = R_ns - G
    return A, A_c, A_s
end

function fraction_wet_vegetation(w_r::T, LAI::T, c::T=T(0.2)) where {T}
    w_rmax = c * LAI
    f_wet = min(T(1), (w_r / w_rmax)^(2 / 3))
    return f_wet
end

function canopy_drainage(P::T, w_r::T, f_veg::T, c::T=T(0.2)) where {T}
    k = (1 - f_veg) * P
    b = 1 / (2 * c)
    D_c = k * (exp(b * w_r) - 1)
    return D_c
end

function vpd_veg_source_height(VPD_a::T, T_a::T, p_a::T, A::T, λE::T, r_aa::T) where {T}
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
