abstract type AbstractModel end

"""
    ProcessBasedModel{FT}(; forcings, parameters, t_span, u0, saveat, kwargs...)

Process-based model of the water balance of soil and canopy, solved with [`solve!`](@ref)
after [`initialize!`](@ref).

`solver_kwargs` are passed to `OrdinaryDiffEq.solve` by [`solve!`](@ref). The default
tolerances, `abstol = reltol = 1e-6`, following SUMMA with SUNDIALS
[Spiteri et al., 2024](https://doi.org/10.1029/2024MS004256).
Keyword arguments given to [`solve!`](@ref) take precedence over `solver_kwargs`.
"""
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
    thresholds::ThresholdTreatment = HardThresholds()
    solver_kwargs::NamedTuple = (; abstol=1e-6, reltol=1e-6)
end

function initialize!(model::ProcessBasedModel)
    model.f = create_rhs(model)
    model.f_diagnostics = create_f_diagnostics(model)
    model.prob = ODEProblem(model.f, model.u0, model.t_span, model.parameters)
    return nothing
end

function create_rhs(model::ProcessBasedModel)
    return let forcings = model.forcings, thresholds = model.thresholds
        (du, u, p, t) -> compute_tendencies!(du, u, p, t, forcings, thresholds)
    end
end

function create_f_diagnostics(model::ProcessBasedModel)
    return let forcings = model.forcings, thresholds = model.thresholds
        (u, p, t) -> compute_diagnostics(u, p, t, forcings, thresholds)
    end
end

"""
    solve!(model::ProcessBasedModel; AD=false, kwargs...)

Solve the model, saving the diagnostics at `model.saveat`. `kwargs` are passed to
`OrdinaryDiffEq.solve` and take precedence over `model.solver_kwargs`. A `callback` is
applied before the diagnostics are saved. yaxarray_output=true saves the output in a
`YAXArray datacube, but is false by default because it interferes with automatic
differentiation application.

The solver is set with `alg`, e.g. `Tsit5()`, `Heun()` or `ImplicitEuler()`. The Jacobian of
an implicit solver is computed with ForwardDiff, `ImplicitEuler(; autodiff=AutoForwardDiff())`,
or with Enzyme,
`ImplicitEuler(; autodiff=AutoEnzyme(; function_annotation=EvaporationModel.Enzyme.Duplicated))`.
`Duplicated` is needed because the right-hand side is a closure over the forcings.

Any OrdinaryDiffEq solver can be used, but only the options above are exported by this package.
For other solvers, the user has to manage installation of required packages.
"""
function solve!(model::ProcessBasedModel; yaxarray_output=false, callback=nothing, kwargs...)
    cb = SavingCallback(
        (u, t, integrator) -> model.f_diagnostics(u, integrator.p, t),
        model.diagnostics;
        saveat=model.saveat,
    )
    model.sol = solve(
        model.prob;
        callback=isnothing(callback) ? cb : CallbackSet(callback, cb),
        saveat=model.saveat,
        tstops=model.tstops,
        model.solver_kwargs...,
        kwargs...,
    )

    # Save data in datacube
    if yaxarray_output
        df_diagnostics = DataFrame(model.diagnostics.saveval)
        df_prognostics = DataFrame(model.sol)
        cols_prognostics = filter(x -> x != "timestamp", names(df_prognostics))
        rename!(df_prognostics, map(=>, cols_prognostics, ["w_1", "w_2", "w_r"]))
        df_all = hcat(df_diagnostics, df_prognostics)
        df_all_no_time = df_all[:, names(df_all) .!= "timestamp"]
        axlist = (
            YAXArrays.time(unix2datetime.(model.saveat)), Variables(names(df_all_no_time))
        )
        model.output = YAXArray(axlist, Array(df_all_no_time))
    end
    return nothing
end

@inline function compute_diagnostics(
    u, p::AbstractArray, t, forcings::NamedTuple, thresholds::ThresholdTreatment=HardThresholds()
)
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
    r_smin,
    k_ext = p
    P, T_a, u_a, p_a, VPD_a, SW_in, R_n, LAI = (
        forcings.P(t),
        forcings.T_a(t),
        forcings.u_a(t),
        forcings.p_a(t),
        forcings.VPD_a(t),
        forcings.SW_in(t),
        forcings.R_n(t),
        forcings.LAI(t),
    )
    d_c, z_0mc = Bigleaf.roughness_parameters(
        RoughnessCanopyHeightLAI(), h, LAI; hs=z_0ms
    )
    f_veg = fractional_vegetation_cover(LAI, k_ext)
    w_rmax = max_canopy_capacity(LAI)
    f_wet = fraction_wet_vegetation(w_r, w_rmax, thresholds)

    w_1eq = w_geq(w_2, w_sat, a, p_soil) #no allocs
    C_1 = c_1(w_1, w_sat, b, C_1sat, w_wp, thresholds) # no allocs
    C_2 = c_2(w_2, w_sat, C_2ref) # no allocs

    t_sol = seconds_since_solar_noon(t, forcings.lon, forcings.utc_offset)
    R_nc, R_ns = net_radiation_partitioning(R_n, f_veg)
    G = ground_heat_flux(SantanelloFriedl03(), R_ns, w_1, w_sat, t_sol)
    A, A_c, A_s = available_energy_partitioning(R_nc, R_ns, G)

    # Resistances
    ustar = ustar_from_u(u_a, z_obs, d_c, z_0mc)
    r_aa = Bigleaf.compute_Ram(ResistanceWindZr(), ustar, u_a)
    r_ac = (Bigleaf.Gb_constant_kB1(ustar, kB⁻¹))^-1
    r_as = soil_aerodynamic_resistance(Choudhury1988soil(), ustar, h, d_c, z_0mc, z_0ms)
    r_sc = surface_resistance(
        JarvisStewart(),
        SW_in,
        VPD_a,
        T_a,
        w_2,
        w_fc,
        w_wp,
        LAI,
        g_d,
        r_smin;
        thresholds,
    )
    β = soil_evaporation_efficiency(Pielke92(), w_1, w_fc)
    r_ss = beta_to_r_ss(β, r_as)

    # Turbulent fluxes calculations
    λE_tot, λE_tot_p = total_evaporation(
        T_a,
        p_a,
        VPD_a,
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
        VPD_a, T_a, p_a, A, λE_tot, r_aa
    )
    E_t, λE_t = transpiration(
        T_a, p_a, VPD_m, A_c, r_ac, r_sc, f_wet
    )
    E_i, λE_i = interception_loss(T_a, p_a, VPD_m, A_c, r_ac, f_wet)
    E_s, λE_s = soil_evaporation(T_a, p_a, VPD_m, A_s, r_as, r_ss)

    P_c = canopy_input(P, f_veg)
    D_c = canopy_drainage(P, w_r, f_veg, k_ext)
    P_s = precip_below_canopy(P, P_c, D_c)
    Q_s = surface_runoff(StaticInfiltration(), P_s, w_2, w_sat)
    D_1 = diffusion_layer_1(w_1, w_1eq, C_2)
    K_2 = vertical_drainage_layer_2(w_2, w_fc, C_3, d_2, thresholds)
    I_s = P_s - Q_s
    f_1 = surface_infiltration_factor(w_1, w_sat) # bounds w_1 at w_sat
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
        P_c=P_c,
        D_c=D_c,
        P_s=P_s,
        Q_s=Q_s,
        D_1=D_1,
        K_2=K_2,
        I_s=I_s,
        f_1=f_1,
    )
end

function compute_tendencies!(
    du, u, p::AbstractArray, t, forcings::NamedTuple, thresholds::ThresholdTreatment=HardThresholds()
)
    diagnostics = compute_diagnostics(u, p, t, forcings, thresholds)
    @unpack d_1, d_2 = p
    @unpack P_c, D_c, I_s, f_1, D_1, K_2, E_s, E_t, E_i, C_1 = diagnostics
    du[1] = C_1 / (ρ_w * d_1) * (f_1 * I_s - E_s) - D_1
    du[2] = 1 / (ρ_w * d_2) * (I_s - E_s - E_t) - K_2
    du[3] = P_c - E_i - D_c
    return nothing
end
