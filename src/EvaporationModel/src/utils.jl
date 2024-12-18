"""
    fourier_series(t, coeffs, ω)

Function that gives the Fourier series, as presented in equation
5 of [Murray & Verhoef](https://doi.org/10.1016/j.agrformet.2003.11.005) 

``f(t) = a_0 + \\sum_{n=1}^M \\left( a_n \\cos(n \\omega t) + 
b_n \\sin(n \\omega t) \\right)``

# Arguments

- `t`: Timestamp in seconds 
- `coeffs::ComponentArray`: array holding the coefficients
    - `a0`: Mean component
    - `an`: Cosine coefficients
    - `bn`: Sine coefficients
- `ω`: Radial frequency [rad/s]
"""
function fourier_series(t, coeffs::ComponentArray, ω)
    @unpack a0, an, bn = coeffs
    M = length(an)
    series =
        a0 +
        sum(an[n] * cos(n * ω * t) for n in 1:M) +
        sum(bn[n] * sin(n * ω * t) for n in 1:M)
    return series
end

function fit_fourier_coefficients(t::AbstractVector, x::AbstractVector, M::Int, ω::Real)
    X = hcat(
        ones(length(t)),
        [cos.(n * ω * t) for n in 1:M]...,
        [sin.(n * ω * t) for n in 1:M]...,
    )
    prob = LinearProblem(X, x)
    linsolve = init(prob)
    sol = solve(linsolve)
    coeffs = sol.u
    return ComponentArray(; a0=coeffs[1], an=coeffs[2:(M + 1)], bn=coeffs[(M + 2):end])
end

"""
    compute_amplitude_and_phase(an::AbstractVector, bn::AbstractVector)

Translate the fourier coefficients from sin-cos form (see 
[fourier_series](@ref fourier_series)) to amplitude phase form:

``f(t) = a_0 + \\sum_{n=1}^M \\left( a_{bn} \\sin(n \\omega t + \\phi) \\right)``

Specifically, it returns both `a_{bn}` and `\\phi` vectors.
FYI: proof of conversion found [here](http://wrean.ca/cazelais/math252/lc-trig.pdf).
Crucial to use atan2 function
"""
function compute_amplitude_and_phase(an::AbstractVector, bn::AbstractVector)
    a_bn = @. √(an^2 + bn^2)
    ϕ = @. atan(an, bn)
    return a_bn, ϕ
end
