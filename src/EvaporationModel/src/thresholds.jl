"""
    ThresholdTreatment

How the thresholds in the model (floors, caps, storage bounds) are evaluated:
[`HardThresholds`](@ref) or [`KavetskiSmoothing`](@ref). Set it for a whole model through
the `thresholds` field of `ProcessBasedModel`.
"""
abstract type ThresholdTreatment end

"""
    HardThresholds()

Evaluate thresholds exactly with `max`, `min` and `clamp` (the kernels are one).
"""
struct HardThresholds <: ThresholdTreatment end

"""
    KavetskiSmoothing(; s_w=1e-3, s_r=0.01, s_f=0.01)

Smooth all thresholds following
[Kavetski & Kuczera (2007)](https://doi.org/10.1029/2006WR005195): `max`/`min` become
`smooth_max`/`smooth_min` (their Eqs. 11 and 18) and storage bounds the logistic step
smoother of their Eq. 8. Each smoothing scale has the units of the thresholded quantity:

- `s_w`: soil moisture scale [m³ m⁻³]
- `s_r`: canopy storage scale, as a fraction of `w_rmax` [-]
- `s_f`: scale for dimensionless factors in [0, 1] (Jarvis stress factors, `f_wet`) [-]
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
    threshold_max(t::ThresholdTreatment, a, b, s)
    threshold_min(t::ThresholdTreatment, a, b, s)

`max(a, b)` / `min(a, b)` for [`HardThresholds`](@ref); `smooth_max` / `smooth_min` with
``m = s^2`` for [`KavetskiSmoothing`](@ref), where `s` has the units of `a` and `b` (the
smooth result is offset by `s/2` at `a = b`).
"""
threshold_max(::HardThresholds, a, b, s) = max(a, b)
threshold_max(::KavetskiSmoothing, a, b, s) = smooth_max(a, b, s^2)
threshold_min(::HardThresholds, a, b, s) = min(a, b)
threshold_min(::KavetskiSmoothing, a, b, s) = smooth_min(a, b, s^2)

"""
    threshold_clamp(t::ThresholdTreatment, x, lo, hi, s)

Clamp `x` to `[lo, hi]`. The minimum is applied first and the maximum last, so the smooth
version never drops below `lo` (it may exceed `hi` by at most `s/2`).
"""
threshold_clamp(t::ThresholdTreatment, x, lo, hi, s) =
    threshold_max(t, threshold_min(t, x, hi, s), lo, s)

"""
    lower_bound_kernel(t::ThresholdTreatment, x, x_min, s)

Kernel that switches off a flux as a storage `x` approaches its lower bound `x_min`. It is
one for [`HardThresholds`](@ref). For [`KavetskiSmoothing`](@ref) it is the logistic step
smoother of Eq. 8 in Kavetski & Kuczera (2007), ``1 / (1 + e^{-(x - x_{min})/s})``: one far
above the bound, 1/2 at the bound, and zero below it. Unlike the exponential kernel of their
Eq. 20, it stays in [0, 1], so a trial step past the bound cannot produce an unbounded flux.
"""
lower_bound_kernel(::HardThresholds, x, x_min, s) = one(x)
lower_bound_kernel(::KavetskiSmoothing, x, x_min, s) = 1 / (1 + exp(-(x - x_min) / s))
