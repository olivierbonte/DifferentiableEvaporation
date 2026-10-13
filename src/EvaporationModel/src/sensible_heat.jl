"""Calculate sensible heat flux from total energy balance."""
function sensible_heat(A, λE)
    return A - λE
end

"""Calculate sensible heat flux from soil energy balance."""
function sensible_heat_soil(A_s, λE_s)
    return A_s - λE_s
end

"""Calculate sensible heat flux from canopy energy balance."""
function sensible_heat_canopy(A_c, λE_c)
    return A_c - λE_c
end
