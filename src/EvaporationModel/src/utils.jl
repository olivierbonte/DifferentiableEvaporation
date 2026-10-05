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
    return ComponentArray(; a0=coeffs[1], an=coeffs[2:(M+1)], bn=coeffs[(M+2):end])
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

# NOAA equation of time [min] for the fractional year γ [rad]
function equation_of_time(γ::T) where {T}
    return T(229.18) * (
        T(0.000075) + T(0.001868) * cos(γ) - T(0.032077) * sin(γ) -
            T(0.014615) * cos(2γ) - T(0.040849) * sin(2γ)
    )
end

"""
    local_to_solar_time(time, timezone, lon)

# Arguments
- `time`: DateTime object given the current time
- `timezone`: Timezone in hours ahead of UTC (e.g. +1 for Brussels)
- `lon`: Longitude in degrees

# Returns
- `t_sol`: DateTime object of the local solar time

# Details
Caluculations based on equations from
[NOAA](https://gml.noaa.gov/grad/solcalc/solareqns.PDF) and
[pveducation](https://www.pveducation.org/pvcdrom/properties-of-sunlight/solar-time).
For the equation ot time, the equation from NOAA is used.

# Examples
Below an example to check if calculation matches calculator provided
[here](https://www.pveducation.org/pvcdrom/properties-of-sunlight/solar-time)

```jldoctest
using Dates
using EvaporationModel
lon = 150 # °
timezone = 10 # UTC+10
time = DateTime(2003, 1, 5, 12, 30)
t_sol = local_to_solar_time(time, timezone, lon)
true_t_sol_hour = 12
true_t_sol_minute = 24
# Check if correct within minute
[hour(t_sol) - true_t_sol_hour, minute(t_sol) - true_t_sol_minute] ≤ [0, 1]

# output

true
```
"""
function local_to_solar_time(time::DateTime, timezone::Integer, lon)
    doy = dayofyear(time)
    # B = 350/365 * (doy - 81)
    # EoT = 9.87 * sin(2B) - 7.53 * cos(B) - 1.5 * sin(B) # Equation of Time, in minutes
    days_in_year = daysinyear(year(time))
    γ = 2 * π / days_in_year * (doy - 1 + (hour(time) - 12) / 24) #fractional year
    EoT = equation_of_time(γ)
    LSTM = 4 * (lon - 15 * timezone) + EoT # Local Standard Time Meridian , in minutes
    LSTM_ms = round(LSTM * 60 * 1000) # Local Standard Time Meridian, in milliseconds
    t_sol = time + Millisecond(LSTM_ms) # Local solar time
    return t_sol
end

"""
    seconds_since_solar_noon(t_sol)

Given the local solar time `t_sol`, this function returns `t_diff`, which is the number
of seconds since solar noon, which takes places at 12:00:00 of that day (in local solar time).
"""
function seconds_since_solar_noon(t_sol::DateTime)
    date = Date(t_sol)
    time = Time(12, 0, 0)
    solar_noon = DateTime(date, time)
    # Calculate the difference in seconds between t_sol and solar noon
    t_diff = round(t_sol - solar_noon, Dates.Second)
    return t_diff
end

"""
    seconds_since_solar_noon(t::Real, lon, utc_offset)

Seconds since solar noon [s], in [-43200, 43200), for `t` the local standard time in
seconds since 1970-01-01T00:00, `lon` the longitude [°E] and `utc_offset` the offset of
local standard time to UTC [h]. Negative before solar noon. The NOAA equation of time is
evaluated with a mean year length, so no `DateTime` is needed (AD-compatible).

# Examples

```jldoctest
using EvaporationModel, Dates
lon, timezone = 4.52, 1
time = DateTime(2010, 7, 1, 12, 30)
t_ref = seconds_since_solar_noon(local_to_solar_time(time, timezone, lon)).value
t = seconds_since_solar_noon(datetime2unix(time), lon, timezone)
abs(t - t_ref) < 60

# output

true
```
"""
function seconds_since_solar_noon(t::Real, lon, utc_offset)
    T = value_type(t)
    day = T(86_400)
    t_utc = t - 3600 * utc_offset
    γ = 2 * T(π) * mod(t_utc / day - T(0.5), T(365.2425)) / T(365.2425) # fractional year
    EoT = equation_of_time(γ)
    t_solar = t_utc + 60 * (4 * lon + EoT) # time offset in minutes
    return mod(t_solar, day) - day / 2
end

@inline function value_type(x::T) where {T}
    if T <: ForwardDiff.Dual
        return typeof(ForwardDiff.value(x))
    else
        return T
    end
end

@inline function of_value_type(x, y)
    return convert(value_type(x), y)
end
