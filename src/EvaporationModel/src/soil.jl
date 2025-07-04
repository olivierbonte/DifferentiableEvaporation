abstract type bMethod end
struct Clay <: bMethod end
struct VanGenuchten <: bMethod end

"""
    c_1(w_1, w_sat, b, c_1sat)

Compute force coefficient `c_1` [-] of force restore framework for soil mositure.
See equation 20 of [Noilhan & Mahfouf, 1996](https://doi.org/10.1016/0921-8181(95)00043-7).

# Arguments
- `w_1`: Surface soil moisutre [m³ m⁻³]
- `w_sat`: Saturated soil moisture [m³ m⁻³]
- `b`: the Brooks-Corey/Clapp-Hornberger parameter, see [`compute_b`](@ref compute_b)
- `c_1sat`: See [c_1sat](@ref c_1sat)

"""
function c_1(w_1, w_sat, b, c_1sat)
    return c_1sat * (w_1 / w_sat)^(b / 2 + 1)
end

"""
    c_2(w_2, w_sat, c2_ref)

Compute restore coefficient `c_2` [-] of force restore framework for soil moisture
See equation 21 of [Noilhan & Mahfouf, 1996](https://doi.org/10.1016/0921-8181(95)00043-7).

# Arguments
- `w_2`: The second layer soil mositure [m³ m⁻³]
- `w_sat`: The saturated soil moisture [m³ m⁻³]
- `c_2ref`: See [`c_2ref`](@ref c_2ref)

"""
function c_2(w_2, w_sat, c_2ref)
    return c_2ref * (w_2 / (w_sat - w_2 + oftype(c2_ref, 0.01)))
end

"""
    w_geq(w_2, w_sat, a, p)

Compute `w_geq` [m³ m⁻³], the equilibrium surface soil moisture (i.e. when capillary and
gravitational forces are in equilibrium).
See equation 19 of [Noilhan & Mahfouf, 1996](https://doi.org/10.1016/0921-8181(95)00043-7).

# Arguments
- `w_2`: The second layer soil moisture [m³ m⁻³].
- `w_sat`: The saturated soil moisture [m³ m⁻³].
- `a`: Clapp-Hornberger parameter `a` (see [`compute_a`](@ref compute_a))
- `p`: Clapp-Hornberger parameter `a` (see [`compute_p`](@ref compute_p))
"""
function w_geq(w_2, w_sat, a, p)
    return w_2 - a * w_sat * (w_2 / w_sat)^p * (1 - (w_2 / w_sat)^(8 * p))
end

"""
    compute_b(approach::Clay, x_clay)
    compute_b(approach::VanGenuchten, n)

Compute `b` [-], the Brooks-Corey/Clapp-Hornberger parameter  (see equation 1
of [Clapp & Hornberger](https://doi.org/10.1029/WR014i004p00601)
for its definition), based on percentage clay or the van Genuchten paramter `n`.

# Arguments
- `approach`: calculation approach, subtype of `bMethod`.

With `approach = Clay()`:
- `x_clay`: The percentage of clay in the soil [%]

With `approach = VanGenuchten()`:
- `n`: The Van Genuchten parameter `n` [-]

# Details
For `approach = Clay()`, Equation (30) of
[Noilhan & Mahfouf, 1996](https://doi.org/10.1016/0921-8181(95)00043-7) is used.

For`approach = VanGenuchten()`, the parameter equivalence between the Brooks-Corey
and van Genuchten, is based on
[Morel-Seytoux et al., 1996](https://doi.org/10.1029/96WR00069).
Note that in this paper, `M` is equivalent to `b`.
"""
function compute_b(approach::Clay, x_clay::T) where {T}
    return T(0.137) * x_clay + T(3.501)
end

function compute_b(approach::VanGenuchten, n::T) where {T}
    return -1 + n / (n - 1)
end

"""
    c_1sat(x_clay)

Compute the value for force-restore coefficient `c_1` when `w_g =  w_sat`.
See equation 32 of [Noilhan & Mahfouf, 1996](https://doi.org/10.1016/0921-8181(95)00043-7).

# Arguments
- `x_clay`: The percentage of clay in the soil [%]

# Returns
- `c_1sat`: The computed value of `c_1sat` [-]

"""
function c_1sat(x_clay::T) where {T}
    return (T(5.58) * x_clay + T(84.88)) * T(1e-2)
end

"""
    c_2ref(x_clay)

Compute the value for force-restore coefficient `c_2` [-] when `w_2 = 0.5 w_sat`,
`c_2ref` based on the percentage of clay in the soil.
See equation 33 of [Noilhan & Mahfouf, 1996](https://doi.org/10.1016/0921-8181(95)00043-7).

# Arguments
- `x_clay`: The percentage of clay in the soil [%].

"""
function c_2ref(x_clay::T) where {T}
    return T(13.815) * x_clay^(T(-0.954))
end

"""
    c_3(x_clay)

Compute the coefficient for graviational drainage `c_3` [m] based on the percentage of clay
in the soil.
See equation 34 of [Noilhan & Mahfouf, 1996](https://doi.org/10.1016/0921-8181(95)00043-7).

# Arguments
- `x_clay`: The percentage of clay in the soil [%].

"""
function c_3(x_clay::T) where {T}
    return T(5.327) * x_clay^(T(-1.043))
end

"""
    compute_a(x_clay)

Compute `a` [-], a parameter for for `w_geq` calculation, based on percentage of
clay in the soil.
See equation 35 of [Noilhan & Mahfouf, 1996](https://doi.org/10.1016/0921-8181(95)00043-7).

# Arguments
- `x_clay`: The percentage of clay in the soil [%].

"""
function compute_a(x_clay::T) where {T}
    return T(732.43e-3) * x_clay^(T(-0.539))
end

"""
    compute_p(x_clay)

Compute `p` [-], a parameter for for `w_geq` calculation, based on the percentage of clay
in the soil.
See equation 36 of [Noilhan & Mahfouf, 1996](https://doi.org/10.1016/0921-8181(95)00043-7).

# Arguments
- `x_clay`: The percentage of clay in the soil [%].

"""
function compute_p(x_clay::T) where {T}
    return T(0.134) * x_clay + T(3.4)
end
