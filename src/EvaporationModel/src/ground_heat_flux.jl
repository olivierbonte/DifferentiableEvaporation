abstract type GroundHeatFluxMethod end
struct Allen07 <: GroundHeatFluxMethod end
struct SantanelloFriedl03 <: GroundHeatFluxMethod end

function ground_heat_flux(method::Allen07, R_n, LAI)
    if LAI < 0.5
        @error "Ground heat flux from net radiation for LAI < 0.5
        not yet implemented"
    else
        G = R_n * (0.05 + 0.18 * exp(-0.521 * LAI))
    end
    return G
end

"""
    ground_heat_flux(::SantanelloFriedl03, R_ns, w_1, w_sat, t_sol; c_gmin=0.31, c_gmax=0.35,
    t_gmin=74_000, t_gmax=100_000, phase_shift=10_800)

Ground heat flux ``G = c_g \\cos(2π (t_{sol} + t_{shift}) / t_g) R_{ns}`` [W m⁻²], with
`t_sol` the time since solar noon [s] - negative before noon - and `R_ns` the net radiation
at the soil surface [W m⁻²].

Equation 4 of [Santanello & Friedl (2003)](https://doi.org/10.1175/1520-0450(2003)042%3C0851:DCISHF%3E2.0.CO;2),
applied to `R_ns` as in ALEXI ([Anderson et al. 2018](https://lpdaac.usgs.gov/documents/332/ECO3ETALEXIU_ATBD_V1.pdf), Eq. 15),
with ``c_g`` and ``t_g`` interpolated linearly in ``Θ = w_1 / w_{sat}`` as in
[Mallick et al. (2022)](https://doi.org/10.1029/2021GL097568), Eqs. S1.19-S1.20.
"""
function ground_heat_flux(
    method::SantanelloFriedl03,
    R_ns,
    w_1,
    w_sat,
    t_sol;
    c_gmin=of_value_type(w_1, 0.31),
    c_gmax=of_value_type(w_1, 0.35),
    t_gmin=of_value_type(w_1, 74_000),
    t_gmax=of_value_type(w_1, 100_000),
    phase_shift=of_value_type(w_1, 10_800),
)
    Θ = w_1 / w_sat
    c_g = c_gmax * (1 - Θ) + c_gmin * Θ
    t_g = t_gmax * (1 - Θ) + t_gmin * Θ
    G = c_g * cos(2π * (t_sol + phase_shift) / t_g) * R_ns
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
