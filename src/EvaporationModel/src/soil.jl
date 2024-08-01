"""
    compute_c_1(w_g, w_sat, b, c_1_sat)

Compute force coefficient `c_1` of force restore framework for soil mositure.
See equation 20 of [Noilhan & Mahfouf, 1996](https://doi.org/10.1016/0921-8181(95)00043-7).

# Arguments
- `w_g`: Surface soil moisutre [m³ m⁻³]
- `w_sat`: Saturated soil moisture [m³ m⁻³]
- `b`: the Brooks-Corey/Clapp-Hornberger parameter, see [compute_b](@ref compute_b)
- `c_1_sat`: See [compute_c_1_sat](@ref compute_c_1_sat)

"""
function compute_c_1(w_g, w_sat, b, c_1_sat)
    return c_1_sat * (w_g / w_sat)^(b / 2. + 1.)
end

"""
    compute_c_2(w_2, w_sat, c2_ref)

Compute restore coefficient `c_2` of force restore framework for soil moisture
See equation 21 of [Noilhan & Mahfouf, 1996](https://doi.org/10.1016/0921-8181(95)00043-7).

# Arguments
- `w_2`: The second layer soil mositure [m³ m⁻³]
- `w_sat`: The saturated soil moisture [m³ m⁻³]
- `c2_ref`: See [compute_c_2_ref](@ref compute_c_2_ref)

"""
function compute_c_2(w_2, w_sat, c2_ref)
    return c2_ref * (w_2 / (w_sat - w_2 + 0.01))
end

"""
    compute_w_geq(w_2, w_sat, a, p)

Compute `w_geq`, the equilibrium surface soil moisture (i.e. when capillary and 
gravitational forces are in equilibrium).
See equation 19 of [Noilhan & Mahfouf, 1996](https://doi.org/10.1016/0921-8181(95)00043-7).

# Arguments
- `w_2`: The second layer soil moisture [m³ m⁻³].
- `w_sat`: The saturated soil moisture [m³ m⁻³].
- `a`: Clapp-Hornberger parameter `a` (see [compute_a](@ref compute_a))
- `p`: Clapp-Hornberger parameter `a` (see [compute_p](@ref compute_p))
"""
function compute_w_geq(w_2, w_sat, a, p)
    return w_2 - a * w_sat * (w_2 / w_sat)^p * (1. - (w_2 / w_sat)^(8. * p))
end

"""
    compute_b(::Val{:clay}, perc_clay)
    compute_b(::Val{:van_genuchten}, n)

Compute `b`, the Brooks-Corey/Clapp-Hornberger parameter  (see equation 1
 of [Clapp & Hornberger](https://doi.org/10.1029/WR014i004p00601) 
for its definition), based on percentage clay or the van Genuchten paramter `n`.	

# Arguments
- `approach`: Either `Val(:clay)` or `Val(:van_genuchten)`.
- `perc_clay`: The percentage of clay in the soil [%], only for `approach = Val(:clay)`.
- `n`: The Van Genuchten parameter `n`, only for `approach = Val(:van_genuchten)`.	

# Details 
For the `Val(:clay)` approach, Equation (30) of 
[Noilhan & Mahfouf, 1996](https://doi.org/10.1016/0921-8181(95)00043-7) is used.

For the `Val(:van_genuchten)` approach, the parameter equivalence between the Brooks-Corey
and van Genuchten, is based on [Morel-Seytoux et al., 1996](https://doi.org/10.1029/96WR00069). 
Note that in this paper, `M` is equivalent to `b`. 
"""
function compute_b(approach, perc_clay)
    return compute_b(approach, perc_clay)
end

function compute_b(::Val{:clay}, perc_clay)
    return 0.137 * perc_clay + 3.501
end

function compute_b(::Val{:van_genuchten}, n)
    return -1. + n / (n - 1.)
end

"""
    compute_c_1_sat(perc_clay)

Compute the value for force-restore coefficient `c_1` when `w_g =  w_sat`, 
`c_2_ref`  based on the percentage of clay in the soil.
See equation 32 of [Noilhan & Mahfouf, 1996](https://doi.org/10.1016/0921-8181(95)00043-7).

# Arguments
- `perc_clay`: The percentage of clay in the soil.

# Returns
- `c_1_sat`: The computed value of `c_1_sat`.

"""
function compute_c_1_sat(perc_clay)
    return (5.58 * perc_clay + 84.88) * 1e-2
end

"""
    compute_c_2_ref(perc_clay)

Compute the value for force-restore coefficient `c_2` when `w_2 = 0.5 w_sat`, 
`c_2_ref`  based on the percentage of clay in the soil.
See equation 33 of [Noilhan & Mahfouf, 1996](https://doi.org/10.1016/0921-8181(95)00043-7).

# Arguments
- `perc_clay`: The percentage of clay in the soil [%].

"""
function compute_c_2_ref(perc_clay)
    return 13.815 * perc_clay^(-0.954)
end

"""
    compute_c_3(perc_clay)

Compute the coefficient for graviational drainage `c_3` based on the percentage of clay in the soil.
See equation 34 of [Noilhan & Mahfouf, 1996](https://doi.org/10.1016/0921-8181(95)00043-7).

# Arguments
- `perc_clay`: The percentage of clay in the soil [%].

"""
function compute_c_3(perc_clay)
    return 5.327 * perc_clay^(-1.043)
end

"""
    compute_a(perc_clay)

Compute a, a parameter for for w_geq calculation,  of clay in the soil.
See equation 35 of [Noilhan & Mahfouf, 1996](https://doi.org/10.1016/0921-8181(95)00043-7).

# Arguments
- `perc_clay`: The percentage of clay in the soil [%].

"""
function compute_a(perc_clay)
    return 732.43e-3 * perc_clay^(-0.539)
end

"""
    compute_p(perc_clay)

Compute `p`, a parameter for for w_geq calculation, based on the percentage of clay in the soil.
See equation 36 of [Noilhan & Mahfouf, 1996](https://doi.org/10.1016/0921-8181(95)00043-7).

# Arguments
- `perc_clay`: The percentage of clay in the soil [%].

"""
function compute_p(perc_clay)
    return 0.134 * perc_clay + 3.4
end
