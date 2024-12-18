"""
    compute_g_from_r_net(r_net, lai)

Compute the ground heat flux [W/m²] from net radiation (and LAI)

# Arguments
- `r_net`: The net radiation [W/m²].
- `lai`: The leaf area index [m²/m²].

# Details
This implementation follows the approach of the METRIC model.
See Equation 27 of 
[Allen et al., 2007](https://doi.org/10.1061/(ASCE)0733-9437(2007)133:4(380))

"""
function compute_g_from_r_net(r_net, lai)
    if lai < 0.5
        @error "Ground heat flux from net radiation for LAI < 0.5
        not yet implemented"
    else
        g = r_net * (0.05 + 0.18 * exp(-0.521 * lai))
    end
    return g
end

"""
    compute_harmonic_sum(t::Real, a_bn::AbstractVector, ϕ::AbstractVector, 
    ω::AbstractVector, Δt::Int)

Compute the sum of harmonic terms `\\Gamma_s` as defined in equation 1 of
[Murray and Verhoef (2007)](https://doi.org/10.1016/j.agrformet.2007.06.009)
"""
function compute_harmonic_sum(
    t::Real, a_bn::AbstractVector, ϕ::AbstractVector, ω::Real, Δt::Real
)
    return sum(
        a_bn[n] * √(n * ω) * sin(n * ω * t + ϕ[n] + π / 4 - π * Δt / 12) for n in 1:M_terms
    )
end

"""
    compute_harmonic_sum(t::AbstractVector, a_bn::AbstractVector, ϕ::AbstractVector, 
    ω::Real, Δt::Real)

Applies broadcasting of the function for when `t::AbstractVector` 
"""
function compute_harmonic_sum(
    t::AbstractVector, a_bn::AbstractVector, ϕ::AbstractVector, ω::Real, Δt::Real
)
    return compute_harmonic_sum.(t, Ref(a_bn), Ref(ϕ), ω, Δt)
end
