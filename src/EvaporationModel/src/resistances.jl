function jarvis_stewart(forcing::ComponentArray, parameters::ComponentArray)
    @unpack s_in, w_2, vpd, t_air, lai = forcing
    @unpack w_fc, w_wilt, gd, r_smin = parameters
    return jarvis_stewart.(s_in, w_2, vpd, t_air, lai, w_fc, w_wilt, gd, r_smin)
end

function jarvis_stewart(s_in, w_2, vpd, t_air, lai, w_fc, w_wilt, gd, r_smin)
    f_1 = min(1.0, (0.004 * s_in + 0.05) / 0.81(0.004 * s_in + 1.0))
    f_2 = max(0.0, min((w_2 - w_wilt) / (w_fc - w_wilt), 1.0))
    f_3 = exp(-gd * vpd)
    #f_4 = @. max(1.0 - 0.0016*(298.0 - t_air)^2, 0.0)
    r_s = r_smin ./ lai * (f_1 * f_2 * f_3)^(-1)
    #max r_s of 50e3 s/max
    r_s = max(50e3, r_s)
    return r_s
end

function aerodynamic_resistance(u::AbstractArray, parameters::ComponentArray)
    @unpack d, z0m, z0h, z_measur = parameters
    return aerodynamic_resistance.(u, d, z0m, z0h, z_measur)
end

function aerodynamic_resistance(u, d, z0m, z0h, z_measur)
    u = max(u, 0.1)
    #I use 0.1 instead of 1 here to allow higher values! https://doi.org/10.1016/S0168-1923(00)00164-7
    r_a = (log((z_measur - d) / z0m) * log((z_measur - d) / z0h)) ./ (0.41^2 * u)
    return r_a
end
