abstract type AbstractModel end
@kwdef mutable struct ProcessBasedModel{FT} <: AbstractModel
    forcings::NamedTuple
    parameters::AbstractArray
    t_span::Tuple
    u0
    saveat::AbstractArray
    tstops::AbstractArray = saveat
    f = nothing
    f_diagnostics = nothing
    prob = nothing
    sol = nothing
    diagnostics = SavedValues(FT, NamedTuple)
    output = nothing
end

function initialize!(model::ProcessBasedModel)
    model.f = create_rhs(model)
    model.f_diagnostics = create_f_diagnostics(model)
    model.prob = ODEProblem(model.f, model.u0, model.t_span, model.parameters)
    return nothing
end

function create_rhs(model::ProcessBasedModel)
    return (du, u, p, t) -> compute_tendencies!(du, u, p, t, model.forcings)
end

function create_f_diagnostics(model::ProcessBasedModel)
    return (u, p, t) -> compute_diagnostics(u, p, t, model.forcings)
end

function solve!(model::ProcessBasedModel; kwargs...)
    cb = SavingCallback(
        (u, t, integrator) -> model.f_diagnostics(u, integrator.p, t),
        model.diagnostics;
        saveat=model.saveat,
    )
    model.sol = solve(
        model.prob; callback=cb, saveat=model.saveat, tstops=model.tstops, kwargs...
    )
    # df_diagnostics = DataFrame(model.diagnostics.saveval)
    # df_prognostics = DataFrame(model.sol)
    # axlist = ()
    return nothing
end

@inline function compute_diagnostics(u, p::AbstractArray, t, forcings::NamedTuple)
    w_1, w_2, w_r = u
    @unpack h,
    z_0ms,
    w_sat,
    a,
    p_soil,
    b,
    w_res,
    w_wp,
    w_fc,
    C_1sat,
    C_2ref,
    C_3,
    d_1,
    d_2,
    z_obs,
    kB⁻¹,
    g_d,
    r_smin = p
    d_c, z_0mc = Bigleaf.roughness_parameters(
        RoughnessCanopyHeightLAI(), h, forcings.LAI(t); hs=z_0ms
    )
    f_veg = fractional_vegetation_cover(forcings.LAI(t))
    w_rmax = max_canopy_capacity(forcings.LAI(t))
    f_wet = fraction_wet_vegetation(w_r, w_rmax)

    w_1eq = w_geq(w_2, w_sat, a, p_soil) #no allocs
    C_1 = c_1(w_1, w_sat, b, C_1sat) # no allocs
    C_2 = c_2(w_2, w_sat, C_2ref) # no allocs

    G = ground_heat_flux(Allen07(), forcings.R_n(t), forcings.LAI(t))
    A, A_c, A_s = available_energy_partioning(forcings.R_n(t), G, f_veg)

    # Resistances
    ustar = ustar_from_u(forcings.u_a(t), z_obs, d_c, z_0mc)
    r_aa = Bigleaf.compute_Ram(ResistanceWindZr(), ustar, forcings.u_a(t))
    r_ac = (Bigleaf.Gb_constant_kB1(ustar, kB⁻¹))^-1
    r_as = soil_aerodynamic_resistance(Choudhury1988soil(), ustar, h, d_c, z_0mc, z_0ms)
    r_sc = surface_resistance(
        JarvisStewart(),
        forcings.SW_in(t),
        forcings.VPD_a(t),
        forcings.T_a(t),
        w_2,
        w_fc,
        w_wp,
        forcings.LAI(t),
        g_d,
        r_smin,
    )
    β = soil_evaporation_efficiency(Pielke92(), w_1, w_fc)
    r_ss = beta_to_r_ss(β, r_as)

    # Turbulent fluxes calculations
    λE_tot, λE_tot_p = total_evaporation(
        forcings.T_a(t),
        forcings.p_a(t),
        forcings.VPD_a(t),
        A,
        A_c,
        A_s,
        r_aa,
        r_ac,
        r_as,
        r_sc,
        r_ss,
        f_wet,
    )
    VPD_m = vpd_veg_source_height(
        forcings.VPD_a(t), forcings.T_a(t), forcings.p_a(t), A, λE_tot, r_aa
    )
    E_t, λE_t = transpiration(
        forcings.T_a(t), forcings.p_a(t), VPD_m, A_c, r_ac, r_sc, f_wet
    )
    E_i, λE_i = interception(forcings.T_a(t), forcings.p_a(t), VPD_m, A_c, r_ac, f_wet)
    E_s, λE_s = soil_evaporation(forcings.T_a(t), forcings.p_a(t), VPD_m, A_s, r_as, r_ss)

    D_c = canopy_drainage(forcings.P(t), w_r, f_veg)
    P_s = precip_below_canopy(forcings.P(t), f_veg, D_c)
    Q_s = surface_runoff(StaticInfiltration(), forcings.P(t), w_2, w_fc)
    D_1 = diffusion_layer_1(w_1, w_1eq, C_2)
    K_2 = vertical_drainage_layer_2(w_2, w_fc, C_3, d_2)
    I_s = P_s - Q_s
    return (
        w_rmax=w_rmax,
        C_1=C_1,
        f_veg=f_veg,
        λE_tot=λE_tot,
        VPD_m=VPD_m,
        E_t=E_t,
        λE_t=λE_t,
        E_i=E_i,
        λE_i=λE_i,
        E_s=E_s,
        λE_s=λE_s,
        D_c=D_c,
        P_s=P_s,
        Q_s=Q_s,
        D_1=D_1,
        K_2=K_2,
        I_s=I_s,
    )
end

function compute_tendencies!(du, u, p::AbstractArray, t, forcings::NamedTuple)
    diagnostics = compute_diagnostics(u, p, t, forcings)
    _, _, w_r = u
    @unpack d_1, d_2 = p
    @unpack D_c, I_s, D_1, K_2, E_s, E_t, E_i, w_rmax, f_veg, C_1 = diagnostics
    # Define smoothing parameter for canopy water content
    m_can = w_rmax / 100
    # Define ODE
    du[1] = C_1 / (ρ_w * d_1) * (I_s - E_s) - D_1
    du[2] = 1 / (ρ_w * d_2) * (I_s - E_s - E_t) - K_2
    du[3] =
        f_veg * forcings.P(t) * smoothing_kernel(UpperBound(), w_r, w_rmax, m_can) -
        E_i * smoothing_kernel(LowerBound(), w_r, zero(w_r), m_can) - D_c
    return nothing
end
