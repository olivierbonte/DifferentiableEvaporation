"""
    fractional_vegetation_cover(LAI, k_ext=0.5)

Fractional vegetation cover ``f_{veg} = 1 - e^{-k_{ext} \\mathrm{LAI}}`` [-]
written as Beer's law in LAI, as in
[De Ridder (2001)](https://doi.org/10.1029/2001JD900128) (``k_{ext} = 1/2``).
"""
function fractional_vegetation_cover(LAI, k_ext=of_value_type(LAI, 0.5))
    return 1 - exp(-k_ext * LAI)
end

"""
    net_radiation_partitioning(R_n, f_veg)

Split the net radiation `R_n` [W m⁻²] into the part absorbed by the canopy,
``R_{nc} = f_{veg} R_n``, and the part reaching the soil, ``R_{ns} = (1 - f_{veg}) R_n``.
"""
function net_radiation_partitioning(R_n, f_veg)
    R_nc = R_n * f_veg
    R_ns = R_n * (1 - f_veg)
    return R_nc, R_ns
end

"""
    available_energy_partitioning(R_nc, R_ns, G)

Available energy of the canopy, ``A_c = R_{nc}``, of the soil, ``A_s = R_{ns} - G``, and in
total, ``A = A_c + A_s`` [W m⁻²].
"""
function available_energy_partitioning(R_nc, R_ns, G)
    A_c = R_nc
    A_s = R_ns - G
    A = A_c + A_s
    return A, A_c, A_s
end

"""
    max_canopy_capacity(LAI, c=0.2)

Maximum canopy water storage ``w_{rmax} = c \\, \\mathrm{LAI}`` [kg m⁻², i.e. mm], with
``c = 0.2`` kg m⁻² (0.2 mm per unit LAI) from
[Noilhan & Planton (1989)](https://doi.org/10.1175/1520-0493(1989)117%3C0536:ASPOLS%3E2.0.CO;2).
"""
function max_canopy_capacity(LAI, c=of_value_type(LAI, 0.2))
    return c * LAI
end

"""
    fraction_wet_vegetation(w_r, w_rmax, t=HardThresholds())

Fraction of the foliage covered with water ``f_{wet} = (w_r / w_{rmax})^{2/3}`` [-], from
[Deardorff (1978)](https://doi.org/10.1029/JC083iC04p01889), with `w_r` clipped at zero and
``f_{wet}`` at one.

With [`KavetskiSmoothing`](@ref), the clip at zero is a smooth maximum (so the slope of the
2/3 power stays finite) times a logistic lower-bound kernel that brings ``f_{wet}`` to zero
below ``w_r = 0``, and the cap at one is a smooth minimum.
"""
function fraction_wet_vegetation(w_r, w_rmax, t=HardThresholds())
    s = storage_scale(t, w_rmax)
    w_r⁺ = threshold_max(t, w_r, zero(w_r), s)
    f_wet = (w_r⁺ / w_rmax)^of_value_type(w_r, 2 / 3) * lower_bound_kernel(t, w_r, zero(w_r), s)
    return threshold_min(t, f_wet, one(f_wet), factor_scale(t))
end

"""
    canopy_input(P, f_veg)

Canopy input ``P_c = f_{veg} P``, the rainfall captured by the canopy (the "canopy input" of the
Rutter model, see [Valente et al. (1997)](https://doi.org/10.1016/S0022-1694(96)03066-1), Fig. 1),
here with canopy and trunk lumped as one.
"""
function canopy_input(P, f_veg)
    return f_veg * P
end

"""
    canopy_drainage(P, w_r, f_veg, k_ext, c=0.2)

Canopy drainage ``D_c = (1 - f_{veg}) P (e^{b w_r} - 1)`` with ``b = k_{ext}/c``, after
De Ridder (2001) Eqs. 39-41, which assume
``k_{ext} = 1/2`` (uniform leaf angle distribution); `b = k_ext/c` extends this to other `k_ext`.
"""
function canopy_drainage(P, w_r, f_veg, k_ext, c=of_value_type(f_veg, 0.2))
    k = (1 - f_veg) * P
    b = k_ext / c
    D_c = k * (exp(b * w_r) - 1)
    return D_c
end

"""
    precip_below_canopy(P, P_c, D_c)

Net precipitation ``P_s = P - P_c + D_c``, i.e. direct throughfall plus canopy drainage
(indirect throughfall), see
[Van Dijk & Bruijnzeel (2001)](https://doi.org/10.1016/S0022-1694(01)00392-4). No stemflow.
"""
function precip_below_canopy(P, P_c, D_c)
    P_s = P - P_c + D_c
    return P_s
end

function vpd_veg_source_height(VPD_a, T_a, p_a, A, λE, r_aa)
    T = value_type(r_aa)
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
