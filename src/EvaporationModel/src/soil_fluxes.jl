abstract type InfiltrationMethod end
struct StaticInfiltration <: InfiltrationMethod end
struct VegetationInfiltration <: InfiltrationMethod end

"""
    surface_runoff(::StaticInfiltration, P_s, w_2, w_sat, p_inf=2)
    surface_runoff(::VegetationInfiltration, P_s, w_2, w_sat, f_veg, s_inf=3)

Bergström β-function surface runoff ``Q_s = (w_2 / w_{sat})^{p} P_s`` on the precipitation
reaching the soil `P_s`, see Eq. 1 of [Trautmann et al. (2022)](https://doi.org/10.5194/hess-26-1089-2022).
With `VegetationInfiltration`, ``p = s_{inf} f_{veg}``.
"""
function surface_runoff(
    approach::StaticInfiltration, P_s, w_2, w_sat, p_inf=of_value_type(w_sat, 2)
)
    Q_s = (w_2 / w_sat)^p_inf * P_s
    return Q_s
end

function surface_runoff(
    approach::VegetationInfiltration, P_s, w_2, w_sat, f_veg, s_inf=of_value_type(f_veg, 3)
)
    p_inf = f_veg * s_inf
    Q_s = surface_runoff(StaticInfiltration(), P_s, w_2, w_sat, p_inf)
    return Q_s
end

"""
    surface_infiltration_factor(w_1, w_sat, m_1=0.01)

Share of the infiltration `I_s` that wets the surface layer,
``f_1 = 1 - e^{(w_1 - w_{sat})/m_1}`` [-],`
It keeps ``w_1 \\le w_{sat}`` without clipping the state:

This is the continuous counterpart of the clip of `w_g` at `w_sat` in SURFEX
([`hydro_soil.F90`](https://github.com/joewkr/open-SURFEX/blob/70d23957e90ac9dfe4f076669f78e967a5e234e3/src/SURFEX/hydro_soil.F90#L513-L526)),
which produces no runoff: `w_1` is not a water store, it only represent as thin surface layer.
`w_2` still receives all of `I_s`and the water balance is unchanged.
The kernel is the exponential smoothing kernel of Eq. 20 in
[Kavetski & Kuczera (2007)](https://doi.org/10.1029/2006WR005195) (see [`smoothing_kernel`](@ref))
`m_1` [m³ m⁻³] sets the width of the transition.
"""
function surface_infiltration_factor(w_1, w_sat, m_1=of_value_type(w_1, 0.01))
    return smoothing_kernel(UpperBound(), w_1, w_sat, m_1)
end

function diffusion_layer_1(w_1, w_1eq, C_2)
    D_1 = C_2 / τ * (w_1 - w_1eq)
    return D_1
end

function vertical_drainage_layer_2(w_2, w_fc, C_3, d_2, t=HardThresholds())
    K_2 = C_3 / (d_2 * τ) * max(t, w_2 - w_fc, zero(w_2), moisture_scale(t))
    return K_2
end
