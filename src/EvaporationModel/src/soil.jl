function compute_c_1(w_g, w_sat, b, C1_sat)
    return C1_sat * (w_g / w_sat)^(b / 2 + 1)
end

function compute_c_2(w_2, w_sat, C2_ref)
    return C2_ref * (w_2 / (w_sat - w_2 + 0.01))
end

function compute_w_geq(w_2, w_sat, a, p)
    return w_2 - a * w_sat * (w_2 / w_sat)^p * (1 - (w_2 / w_sat)^(8 * p))
end
