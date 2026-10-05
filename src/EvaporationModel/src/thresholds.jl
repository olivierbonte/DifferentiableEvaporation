abstract type KernelMethod end
struct LowerBound <: KernelMethod end
struct UpperBound <: KernelMethod end

"""
    smooth_max(a, b, m)

Smooth maximum Eq. 11 of [Kavetski & Kuczera (2007)](https://doi.org/10.1029/2006WR005195),
with `m` in the units of `a` squared. It exceeds `max(a, b)` by at most ``\\sqrt{m}/2``.
"""
function smooth_max(a::FT, b, m) where {FT}
    return convert(FT, 1 / 2) * (a + b + √((a - b)^2 + m))
end

"""
    smooth_min(a, b, m)

Smooth minimum Eq. 18 of [Kavetski & Kuczera (2007)](https://doi.org/10.1029/2006WR005195).
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
    smoothing_kernel(::LowerBound, x, threshold, m)
    smoothing_kernel(::UpperBound, x, threshold, m)

Exponential smoothing kernel Eq. 20 of [Kavetski & Kuczera (2007)](https://doi.org/10.1029/2006WR005195):
zero at the threshold, close to one more than a few `m` away from it, with `m` in the units of `x`.
"""
function smoothing_kernel(approach::LowerBound, x, threshold, m)
    return 1 - exp(-(x - threshold) / m)
end

function smoothing_kernel(approach::UpperBound, x, threshold, m)
    return 1 - exp((x - threshold) / m)
end

"""
    ThresholdTreatment

How the thresholds of the model (floors, caps, storage bounds) are evaluated, for the whole
model through the `thresholds` field of `ProcessBasedModel`.
"""
abstract type ThresholdTreatment end

"""
    HardThresholds()

Exact thresholds.
"""
struct HardThresholds <: ThresholdTreatment end

"""
    KavetskiSmoothing(; s_w=1e-3, s_r=0.01, s_f=0.01)

Smoothed thresholds, following the rational of [Kavetski & Kuczera (2007)](https://doi.org/10.1029/2006WR005195),
with scales in the units of the thresholded quantity:

- `s_w`: soil moisture [m³ m⁻³]
- `s_r`: canopy storage, as a fraction of `w_rmax` [-]
- `s_f`: dimensionless factors in ``[0, 1]`` (Jarvis stress factors, `f_wet`) [-]
"""
Base.@kwdef struct KavetskiSmoothing{T} <: ThresholdTreatment
    s_w::T = 1e-3
    s_r::T = 0.01
    s_f::T = 0.01
end

moisture_scale(t::KavetskiSmoothing) = t.s_w
moisture_scale(::HardThresholds) = false
storage_scale(t::KavetskiSmoothing, w_rmax) = t.s_r * w_rmax
storage_scale(::HardThresholds, w_rmax) = false
factor_scale(t::KavetskiSmoothing) = t.s_f
factor_scale(::HardThresholds) = false

"""
    max(t::ThresholdTreatment, a, b, s)

`max(a, b)` for [`HardThresholds`](@ref), [`smooth_max`](@ref) with ``m = s^2`` for
[`KavetskiSmoothing`](@ref), with the scale `s` in the units of `a` and `b`. `min` and
`clamp` take the treatment and scale in the same way.
"""
Base.max(::HardThresholds, a, b, s) = max(a, b)
Base.max(::KavetskiSmoothing, a, b, s) = smooth_max(a, b, s^2)
Base.min(::HardThresholds, a, b, s) = min(a, b)
Base.min(::KavetskiSmoothing, a, b, s) = smooth_min(a, b, s^2)
Base.clamp(::HardThresholds, x, lo, hi, s) = clamp(x, lo, hi)
Base.clamp(::KavetskiSmoothing, x, lo, hi, s) = smooth_clamp(x, lo, hi, s^2)

"""
    lower_bound_kernel(t::ThresholdTreatment, x, x_min, s)

One for [`HardThresholds`](@ref). For [`KavetskiSmoothing`](@ref) the logistic step
smoother of Eq. 8 of [Kavetski & Kuczera (2007)](https://doi.org/10.1029/2006WR005195),
``1 / (1 + e^{-(x - x_{min})/s})``, which stays in ``[0, 1]`` (unlike
[`smoothing_kernel`](@ref)) and is therefore useful for implicit timestepping.
"""
lower_bound_kernel(::HardThresholds, x, x_min, s) = one(x)
lower_bound_kernel(::KavetskiSmoothing, x, x_min, s) = 1 / (1 + exp(-(x - x_min) / s))
