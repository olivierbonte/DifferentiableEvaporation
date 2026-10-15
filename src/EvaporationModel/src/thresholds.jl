# Thresholds (floors, caps, steps, storage bounds) and their smoothing, after
# [Kavetski & Kuczera (2007)](https://doi.org/10.1029/2006WR005195)

## 1. Smooth approximants

"""
    BoundSide

Side of a bound for [`exponential_kernel`](@ref) and [`bound_kernel`](@ref):
[`LowerBound`](@ref) (``x \\ge`` bound) or [`UpperBound`](@ref) (``x \\le`` bound).
"""
abstract type BoundSide end
"Lower bound, see [`BoundSide`](@ref)."
struct LowerBound <: BoundSide end
"Upper bound, see [`BoundSide`](@ref)."
struct UpperBound <: BoundSide end

"""
    smooth_max(a, b, m)

Smooth maximum, Eq. 11 of [Kavetski & Kuczera (2007)](https://doi.org/10.1029/2006WR005195),
with `m` in the units of `a` squared. It exceeds `max(a, b)` by at most ``\\sqrt{m}/2``.
"""
function smooth_max(a::FT, b, m) where {FT}
    return convert(FT, 1 / 2) * (a + b + √((a - b)^2 + m))
end

"""
    smooth_min(a, b, m)

Smooth minimum, Eq. 18 of [Kavetski & Kuczera (2007)](https://doi.org/10.1029/2006WR005195),
with `m` in the units of `a` squared.
"""
function smooth_min(a::FT, b, m) where {FT}
    return convert(FT, 1 / 2) * (a + b - √((a - b)^2 + m))
end

"""
    smooth_clamp(x, lower, upper, m)

Smooth clamp, [`smooth_min`](@ref) first and [`smooth_max`](@ref) last, so the result never
drops below `lower`.
"""
function smooth_clamp(x::FT, lower, upper, m) where {FT}
    return smooth_max(smooth_min(x, upper, m), lower, m)
end

"""
    logistic_step(x, x_0, m)

Logistic step from 0 (``x \\ll x_0``) to 1 (``x \\gg x_0``), ``1 / (1 + e^{-(x - x_0)/m})``,
Eq. 8 of [Kavetski & Kuczera (2007)](https://doi.org/10.1029/2006WR005195), for step
discontinuities in constitutive functions, with `m` in the units of `x`. It stays in
``[0, 1]`` and equals 1/2 at ``x_0``.
"""
logistic_step(x, x_0, m) = 1 / (1 + exp(-(x - x_0) / m))

"""
    exponential_kernel(::LowerBound, x, bound, m)
    exponential_kernel(::UpperBound, x, bound, m)

Exponential smoothing kernel, Eq. 20 of [Kavetski & Kuczera (2007)](https://doi.org/10.1029/2006WR005195).
Zero at the bound, close to one more than a few `m` inside it and
negative beyond it. Multiplying the outflux of a store by it bounds the store (Eq. 19).
"""
exponential_kernel(::LowerBound, x, bound, m) = 1 - exp(-(x - bound) / m)
exponential_kernel(::UpperBound, x, bound, m) = 1 - exp((x - bound) / m)

## 2. Threshold treatments and their scales

"""
    ThresholdTreatment

How the thresholds of the model (floors, caps, steps, storage bounds) are evaluated, for the
whole model through the `thresholds` field of `ProcessBasedModel`.
"""
abstract type ThresholdTreatment end

"""
    HardThresholds()

Exact thresholds.
"""
struct HardThresholds <: ThresholdTreatment end

"""
    KavetskiSmoothing(; s_w=1e-3, s_r=0.01, s_f=0.01, s_E=2.0)

Smoothed thresholds, following the rationale of [Kavetski & Kuczera (2007)](https://doi.org/10.1029/2006WR005195),
with scales in the units of the thresholded quantity:

- `s_w`: soil moisture [m³ m⁻³]
- `s_r`: canopy storage, as a fraction of `w_rmax` [-]
- `s_f`: dimensionless factors in ``[0, 1]`` (Jarvis stress factors, `f_wet`) [-]
- `s_E`: energy fluxes (soil net radiation, potential canopy latent heat flux) [W m⁻²]
"""
Base.@kwdef struct KavetskiSmoothing{T} <: ThresholdTreatment
    s_w::T = 1e-3
    s_r::T = 0.01
    s_f::T = 0.01
    s_E::T = 2.0
end

moisture_scale(t::KavetskiSmoothing) = t.s_w
moisture_scale(::HardThresholds) = false
storage_scale(t::KavetskiSmoothing, w_rmax) = t.s_r * w_rmax
storage_scale(::HardThresholds, w_rmax) = false
factor_scale(t::KavetskiSmoothing) = t.s_f
factor_scale(::HardThresholds) = false
flux_scale(t::KavetskiSmoothing) = t.s_E
flux_scale(::HardThresholds) = false

## 3. Treatment-aware thresholds

"""
    max(t::ThresholdTreatment, a, b, s)
    min(t::ThresholdTreatment, a, b, s)
    clamp(t::ThresholdTreatment, x, lo, hi, s)

`max`, `min` and `clamp` for [`HardThresholds`](@ref); [`smooth_max`](@ref),
[`smooth_min`](@ref) and [`smooth_clamp`](@ref) with ``m = s^2`` for
[`KavetskiSmoothing`](@ref), with the scale `s` in the units of the arguments.
"""
Base.max(::HardThresholds, a, b, s) = max(a, b)
Base.max(::KavetskiSmoothing, a, b, s) = smooth_max(a, b, s^2)
Base.min(::HardThresholds, a, b, s) = min(a, b)
Base.min(::KavetskiSmoothing, a, b, s) = smooth_min(a, b, s^2)
Base.clamp(::HardThresholds, x, lo, hi, s) = clamp(x, lo, hi)
Base.clamp(::KavetskiSmoothing, x, lo, hi, s) = smooth_clamp(x, lo, hi, s^2)

"""
    step_function(t::ThresholdTreatment, x, x_0, s)

Step from 0 (``x < x_0``) to 1 (``x > x_0``): the Heaviside step for
[`HardThresholds`](@ref), the [`logistic_step`](@ref) with ``m = s`` for
[`KavetskiSmoothing`](@ref).
"""
step_function(::HardThresholds, x, x_0, s) = ifelse(x > x_0, one(x), zero(x))
step_function(::KavetskiSmoothing, x, x_0, s) = logistic_step(x, x_0, s)

"""
    nonneg_pow(t::ThresholdTreatment, x, exponent, s)

``\\max(x, 0)^{exponent}`` for `exponent > 0`, so an unphysical negative `x` (e.g. a
negative soil moisture ratio) gives zero instead of a `DomainError`. The clip at zero is a
[`smooth_max`](@ref) with scale `s` for [`KavetskiSmoothing`](@ref), which keeps the power
defined with a finite slope, and is ``(s/2)^{exponent}`` at ``x = 0``.
"""
nonneg_pow(t::ThresholdTreatment, x, exponent, s) = max(t, x, zero(x), s)^exponent

"""
    bound_kernel(t::ThresholdTreatment, side::BoundSide, x, bound, s)

Factor on the outflux of a store `x` that keeps it on the `side` of `bound`
([Kavetski & Kuczera, 2007](https://doi.org/10.1029/2006WR005195), Eq. 19): one for
[`HardThresholds`](@ref) (the hard model bounds the store itself), the
[`exponential_kernel`](@ref) with ``m = s`` for [`KavetskiSmoothing`](@ref), which is zero at
the bound, so the outflux stops there. Unlike the [`logistic_step`](@ref), which is 1/2 at
its threshold.
"""
bound_kernel(::HardThresholds, side::BoundSide, x, bound, s) = one(x)
bound_kernel(::KavetskiSmoothing, side::BoundSide, x, bound, s) =
    exponential_kernel(side, x, bound, s)
