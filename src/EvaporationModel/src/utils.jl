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
