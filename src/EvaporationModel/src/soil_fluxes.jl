abstract type InfiltrationMethod end
struct StaticInfiltration <: InfiltrationMethod end
struct VegetationInfiltration <: InfiltrationMethod end

function surface_runoff(
    approach::StaticInfiltration, P_s::T, w_2::T, w_fc::T, p_inf::T=T(2)
) where {T}
    Q_s = (w_2 / w_fc)^p_inf * P_s
    return Q_s
end

function surface_runoff(
    approach::VegetationInfiltration, P_s::T, w_2::T, w_fc::T, f_veg::T, s_inf::T=T(3)
) where {T}
    p_inf = f_veg * s_inf
    Q_s = surface_runoff(StaticInfiltration(), P_s, w_2, w_fc, p_inf)
    return Q_s
end

function diffusion_layer_1(w_1::T, w_1eq::T, C_2::T) where {T}
    D_1 = C_2 / τ * (w_1 - w_1eq)
    return D_1
end

function vertical_drainage_layer_2(w_2, w_fc, C_3, d_2)
    K_2 = C_3 / (d_2 * τ) * max(zero(w_2), w_2 - w_fc)
    return K_2
end
