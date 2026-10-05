abstract type bMethod end
struct Clay <: bMethod end
struct VanGenuchten <: bMethod end

"""
    c_1(w_1, w_sat, b, c_1sat, w_wp, t=HardThresholds())

Compute force coefficient `c_1` [-] of force restore framework for soil mositure,
``C_1 = C_{1sat} (w_{sat} / \\max(w_1, w_{wp}))^{b/2 + 1}``.
See equation 20 of [Noilhan & Mahfouf, 1996](https://doi.org/10.1016/0921-8181(95)00043-7).
`w_1` is floored at the wilting point, the lower validity bound of the formula, with
a `max` under the threshold treatment `t`.

# Arguments
- `w_1`: Surface soil moisutre [m³ m⁻³]
- `w_sat`: Saturated soil moisture [m³ m⁻³]
- `b`: the Brooks-Corey/Clapp-Hornberger parameter, see [`compute_b`](@ref compute_b)
- `c_1sat`: ``C_1`` at saturation for the chosen `d_1`, see [`c_1sat`](@ref c_1sat)
- `w_wp`: Soil moisture at wilting point [m³ m⁻³]

"""
function c_1(w_1, w_sat, b, c_1sat, w_wp, t=HardThresholds())
    return c_1sat * (w_sat / max(t, w_1, w_wp, moisture_scale(t)))^(b / 2 + 1)
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
    return c_2ref * (w_2 / (w_sat - w_2 + of_value_type(c_2ref, 0.01)))
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

abstract type C1satMethod end
struct NoilhanMahfouf96 <: C1satMethod end
struct NoilhanPlanton89 <: C1satMethod end

"""
    c_1sat(::NoilhanMahfouf96, x_clay, d_1)
    c_1sat(::NoilhanPlanton89, w_sat, b, ψ_sat, K_sat, d_1)

Force coefficient ``C_{1sat}`` [-], the value of [`c_1`](@ref c_1) at ``w_1 = w_{sat}``, for a
surface layer of normalization depth `d_1` [m].

``C_1`` is dimensionless and proportional to `d_1`: Eq. A4 of
[Noilhan & Planton (1989)](https://doi.org/10.1175/1520-0493(1989)117%3C0536:ASPOLS%3E2.0.CO;2)
gives ``C_1 = 2 d_1 / d``, with ``d`` the penetration depth of the diurnal soil moisture wave.
Only ``C_1 / d_1`` enters the force-restore equation, so any `d_1` gives the same dynamics as
long as ``C_{1sat}`` is computed for that same `d_1`.

- `NoilhanMahfouf96`: regression on the clay percentage `x_clay` [%], Eq. 32 of
  [Noilhan & Mahfouf (1996)](https://doi.org/10.1016/0921-8181(95)00043-7). The regression
  is normalized to ``d_1 = 1`` m (as used in the SURFEX and WRF Pleim–Xiu code), so it is
  scaled by `d_1 / 1 m` here.
- `NoilhanPlanton89`: from the soil hydraulic properties, Eq. A5 of Noilhan & Planton (1989),
  ``C_{1sat} = 2 \\sqrt{\\pi} d_1 \\sqrt{w_{sat} / (b |ψ_{sat}| K_{sat} τ)}``, with `w_sat`
  [m³ m⁻³], the Clapp–Hornberger `b` [-], the saturated matric potential `ψ_sat` [m], the
  saturated hydraulic conductivity `K_sat` [m s⁻¹] and ``τ`` = 1 day. Representative values
  per soil texture are given in Table 2 of
  [Clapp & Hornberger (1978)](https://doi.org/10.1029/WR014i004p00601) (there ``ψ_{sat}`` in
  cm and ``K_{sat}`` in 10⁻⁴ cm s⁻¹; convert to m and m s⁻¹).

# Examples

With the loam row of Table 2 of
[Clapp & Hornberger (1978)](https://doi.org/10.1029/WR014i004p00601) (``w_{sat}`` = 0.451,
``b`` = 5.39, ``ψ_{sat}`` = 47.8 cm, ``K_{sat}`` = 6.95 × 10⁻⁴ cm s⁻¹) and ``d_1 = 0.1`` m,
Eq. A5 gives the Noilhan & Planton (1989) table value of 0.191, close to the clay regression
at 19% clay:

```jldoctest
using EvaporationModel
c_np = c_1sat(NoilhanPlanton89(), 0.451, 5.39, 0.478, 6.95e-6, 0.1)
c_nm = c_1sat(NoilhanMahfouf96(), 19.0, 0.1)
round(c_np; digits=3), round(c_nm; digits=3)

# output

(0.191, 0.191)
```
"""
function c_1sat(::NoilhanMahfouf96, x_clay, d_1)
    T = value_type(x_clay)
    return (T(5.58) * x_clay + T(84.88)) * T(1e-2) * d_1 # regression is for d_1 = 1 m
end

function c_1sat(::NoilhanPlanton89, w_sat, b, ψ_sat, K_sat, d_1)
    return 2 * √(of_value_type(w_sat, π)) * d_1 * √(w_sat / (b * abs(ψ_sat) * K_sat * τ))
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
