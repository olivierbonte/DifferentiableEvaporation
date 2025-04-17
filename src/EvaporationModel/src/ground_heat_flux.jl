abstract type GroundHeatFluxMethod end
struct Allen07 <: GroundHeatFluxMethod end
struct SantanelloFriedl03 <: GroundHeatFluxMethod end

function ground_heat_flux(method::Allen07, R_n::T, LAI::T) where {T}
    if LAI < 0.5
        @error "Ground heat flux from net radiation for LAI < 0.5
        not yet implemented"
    else
        G = R_n * (0.05 + 0.18 * exp(-0.521 * LAI))
    end
    return G
end

function ground_heat_flux(
    method::SantanelloFriedl03,
    R_ns::T,
    w_1::T,
    w_1sat::T,
    t_soln::T,
    c_gmin::T=T(0.31),
    c_gmax::T=T(0.35),
    t_gmin::T=T(74_000),
    t_gmax::T=T(100_000),
) where {T}
    Θ = w_1 / w_1sat
    c_g = c_gmax * (1 - Θ) + c_gmin * Θ
    t_g = t_gmax * (1 - Θ) + t_gmin * Θ
    G = c_g * cos(2π * (t_soln  + 10_800) / t_g) * R_ns
    return G
end

"""
    compute_g_from_r_n(R_n, lai)

Compute the ground heat flux [W/m²] from net radiation (and LAI)

# Arguments
- `R_n`: The net radiation [W/m²].
- `lai`: The leaf area index [m²/m²].

# Details
This implementation follows the approach of the METRIC model.
See Equation 27 of
[Allen et al., 2007](https://doi.org/10.1061/(ASCE)0733-9437(2007)133:4(380))

"""
function compute_g_from_r_n(R_n, lai)
    if lai < 0.5
        @error "Ground heat flux from net radiation for LAI < 0.5
        not yet implemented"
    else
        g = R_n * (0.05 + 0.18 * exp(-0.521 * lai))
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
