function C1(w_g::AbstractFloat, w_sat::AbstractFloat, b::AbstractFloat, C1_sat::AbstractFloat)
    return C1_sat * (w_g / w_sat) ^ (b / 2 + 1)
end

function C_2(w_2::AbstractFloat, w_sat::AbstractFloat, C2_ref::AbstractFloat)
    return C2_ref * (w_2 / (w_sat - w_2))
end

function w_geq(w_2::AbstractFloat, w_sat::AbstractFloat, a::AbstractFloat, p::AbstractFloat)
    return w_2 - a * w_sat * (w_2 / w_sat) ^ p * (1 - (w_2 / w_sat)^(8 * p))
end

