function jarvis_stewart(forcing::ComponentArray, parameters::ComponentArray)
    @unpack s_in, w_2, vpd, t_air, lai = forcing
    @unpack w_fc, w_wilt, gd, r_smin = parameters
    f_1 = @. min(1.0, (0.004*s_in + 0.05)/0.81(0.004*s_in + 1))
    f_2 = @. max(0.0, min((w_2 - w_wilt)/(w_fc - w_wilt), 1.0))
    f_3 = @. exp(-gd * vpd)
    f_4 = @. 1.0 - 0.0016*(298.0 - t_air)^2
    r_s = @. r_smin./lai * (f_1 * f_2 * f_3 * f_4)^(-1)
    return r_s
end

function aerodynamic_resistances(u::AbstractArray, parameters::ComponentArray)
    @unpack d, z0m, z0h, z_measur = parameters
    r_a = (log((z_measur - d)/z0m)*log((z_measur - d)/z0h))/(0.41^2 * u)
    return r_a 
end