abstract type AbstractModel end

"Number of physical states: `w_1`, `w_2`, `w_r`"
const N_STATES = 3

"""
    ProcessBasedModel{FT}(; forcings, parameters, t_span, u0, saveat, kwargs...)

Process-based model of the soil and canopy water balance, solved with [`solve!`](@ref) after
[`initialize!`](@ref).

The ODE state is `u0 = [w_1, w_2, w_r]` plus one accumulator per diagnostic in `accumulated`
(default `(:λE_tot, :H_tot, :G)`): its time integral, a pure quadrature variable as in
[CVODES](https://sundials.readthedocs.io/en/latest/cvodes/Mathematics_link.html#pure-quadrature-integration),
left out of the error control (see [`state_tolerances`](@ref)).

The states are saved at ``t_k = `` `saveat[k]`, ``k = 1, …, n + 1``, and the accumulated
diagnostics are averaged over the windows ``[t_k, t_{k+1}]``, ``k = 1, …, n``: no window
starts at the last time ``t_{n+1}`` (e.g. one forcing time step after the last forcing time
stamp).

`solver_kwargs` (default `abstol = reltol = 1e-6`, as SUMMA-SUNDIALS,
[Spiteri et al., 2024](https://doi.org/10.1029/2024MS004256)) are stored in `model.prob` and
used by [`solve!`](@ref); keyword arguments of `solve!` take precedence.
"""
@kwdef mutable struct ProcessBasedModel{FT} <: AbstractModel
    forcings::NamedTuple
    parameters::AbstractArray
    t_span::Tuple
    u0
    saveat::AbstractArray
    tstops::AbstractArray = saveat
    accumulated::Tuple{Vararg{Symbol}} = (:λE_tot, :H_tot, :G)
    f = nothing
    f_diagnostics = nothing
    prob = nothing
    sol = nothing
    diagnostics = nothing
    output = nothing
    thresholds::ThresholdTreatment = HardThresholds()
    solver_kwargs::NamedTuple = (; abstol=1e-6, reltol=1e-6)
end

"""
    initialize!(model::ProcessBasedModel)

Build the right-hand side, the diagnostics function and the `ODEProblem` (with the
accumulators appended to `u0`) of `model`.
"""
function initialize!(model::ProcessBasedModel)
    length(model.u0) == N_STATES ||
        throw(ArgumentError("`u0` must hold the $N_STATES physical states `[w_1, w_2, w_r]`"))
    model.f = create_rhs(model)
    model.f_diagnostics = create_f_diagnostics(model)
    check_accumulated(model)
    u0 = vcat(model.u0, zeros(eltype(model.u0), length(model.accumulated)))
    model.prob = ODEProblem(
        model.f, u0, model.t_span, model.parameters;
        state_tolerances(model, model.solver_kwargs)...,
    )
    return nothing
end

function check_accumulated(model::ProcessBasedModel)
    d = model.f_diagnostics(model.u0, model.parameters, model.t_span[1])
    for name in model.accumulated
        hasproperty(d, name) || throw(
            ArgumentError("`$name` in `accumulated` is not a diagnostic of `compute_diagnostics`")
        )
    end
    return nothing
end

function create_rhs(model::ProcessBasedModel)
    return let forcings = model.forcings, thresholds = model.thresholds, accumulated = Val(model.accumulated)

        (du, u, p, t) -> compute_tendencies!(du, u, p, t, forcings, thresholds, accumulated)
    end
end

function create_f_diagnostics(model::ProcessBasedModel)
    return let forcings = model.forcings, thresholds = model.thresholds
        (u, p, t) -> compute_diagnostics(u, p, t, forcings, thresholds)
    end
end

"""
    state_tolerances(model::ProcessBasedModel, solver_kwargs::NamedTuple)

`solver_kwargs` with the scalar `abstol` and `reltol` for the physical states replaced by one
tolerance per state. Only scalar tolerances are allowed.

OrdinaryDiffEq scales the error of each state ``i`` as
``e_i / (\\mathrm{abstol}_i + \\max(|u_{\\mathrm{prev},i}|, |u_i|)\\,\\mathrm{reltol}_i)`` and takes the
root mean square over all ``n`` states
([Timestepping Method Descriptions](https://docs.sciml.ai/DiffEqDocs/stable/extras/timestepping/)).
The same norm is used in the Newton convergence test of implicit solvers. The accumulators get
`abstol = Inf`, so their scaled error is 0 and they are left out of the error control, as
[CVODES](https://sundials.readthedocs.io/en/latest/cvodes/Mathematics_link.html#pure-quadrature-integration)
does by default for pure quadratures. The physical-state tolerances are scaled by
``\\sqrt{3/n}``, so that the norm equals the root mean square over the 3 physical states only:
the time steps and physical states are then the same as without accumulators (for implicit
solvers up to the Newton tolerance, as the Newton solve includes the accumulators). Per-state
tolerances are passed to the solver as vectors, see
[Common Solver Options](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/).
Give both `abstol` and `reltol`: a missing one keeps its unscaled solver default.
"""
function state_tolerances(model::ProcessBasedModel, solver_kwargs::NamedTuple)
    # Tolerances rather than a custom `internalnorm`, which would clash with a user's (e.g.
    # Enzyme-inactive) norm and could reach SciMLSensitivity's adjoint solves
    n_accumulators = length(model.accumulated)
    n_total = N_STATES + n_accumulators
    # Compensates for the n_accumulators zero errors in the root mean square over n_total
    norm_scale = sqrt(N_STATES / n_total)

    tolerances = (;)
    if haskey(solver_kwargs, :abstol)
        abstol = solver_kwargs.abstol
        abstol isa Real || throw(ArgumentError("`abstol` must be a scalar (physical states)"))
        abstol_states = fill(abstol * norm_scale, N_STATES)
        abstol_accumulators = fill(Inf, n_accumulators) # left out of the error control
        tolerances = merge(tolerances, (; abstol=vcat(abstol_states, abstol_accumulators)))
    end
    if haskey(solver_kwargs, :reltol)
        reltol = solver_kwargs.reltol
        reltol isa Real || throw(ArgumentError("`reltol` must be a scalar (physical states)"))
        reltol_states = fill(reltol * norm_scale, N_STATES)
        # Irrelevant with abstol = Inf, but must have the same length
        reltol_accumulators = fill(reltol, n_accumulators)
        tolerances = merge(tolerances, (; reltol=vcat(reltol_states, reltol_accumulators)))
    end
    return merge(solver_kwargs, tolerances)
end

"""
    solve!(model::ProcessBasedModel; yaxarray_output=false, kwargs...)

Solve the model, saving the states at `model.saveat`, and store the [`window_means`](@ref)
of the accumulated diagnostics in `model.diagnostics`. `kwargs` (e.g. `alg`, `callback`,
`d_discontinuities`) are passed to `OrdinaryDiffEq.solve`. `yaxarray_output=true` stores the
states and window means in a YAXArrays `Dataset` in `model.output`; off by default as it
interferes with AD.

Implicit solvers take their Jacobian from ForwardDiff, `ImplicitEuler(; autodiff=AutoForwardDiff())`,
or Enzyme, `ImplicitEuler(; autodiff=AutoEnzyme(; function_annotation=EvaporationModel.Enzyme.Duplicated))`
(`Duplicated` because the right-hand side is a closure over the forcings).
"""
function solve!(model::ProcessBasedModel; yaxarray_output=false, kwargs...)
    solver_kwargs = state_tolerances(model, merge(model.solver_kwargs, NamedTuple(kwargs)))
    model.sol = solve(model.prob; saveat=model.saveat, tstops=model.tstops, solver_kwargs...)
    model.diagnostics = window_means(model)

    if yaxarray_output
        model.output = output_dataset(model)
    end
    return nothing
end

"""
    output_dataset(model::ProcessBasedModel)

Datacube (YAXArrays `Dataset`) of the solution: one variable per physical state and
accumulated diagnostic, on the time axis of the saved times ``t_k`` (`model.saveat`). Each
variable has a `description` property of how it relates to ``t_k``.
"""
function output_dataset(model::ProcessBasedModel)
    time_axis = YAXArrays.time(unix2datetime.(model.sol.t))
    variables = Pair{Symbol,YAXArray}[]

    # Physical states: values at the time stamps
    for (i, name) in enumerate((:w_1, :w_2, :w_r))
        values = model.sol[i, :]
        properties = Dict{String,Any}("description" => "value at the time stamp t_k")
        push!(variables, name => YAXArray((time_axis,), values, properties))
    end

    # Accumulated diagnostics: means over the window that starts at the time stamp. No
    # window starts at the last time stamp, hence a NaN there.
    for (name, means) in pairs(model.diagnostics)
        values = vcat(means, NaN)
        properties = Dict{String,Any}(
            "description" =>
                "mean over the window [t_k, t_k+1] starting at the time stamp t_k " *
                    "(NaN at the last time stamp)",
        )
        push!(variables, name => YAXArray((time_axis,), values, properties))
    end

    return Dataset(; variables...)
end

"""
    window_means(sol, accumulated)
    window_means(model::ProcessBasedModel)

Means of the `accumulated` diagnostics (the last states of `sol`) over the windows between
consecutive saved times: entry `k`, for `k = 1, …, length(sol.t) - 1`, is the mean over
`[sol.t[k], sol.t[k+1]]`.
"""
function window_means(sol, accumulated::Tuple{Vararg{Symbol}})
    first_state = sol.u[1] # all saved states have the same length
    n_physical_states = length(first_state) - length(accumulated)
    n_windows = length(sol.t) - 1
    T = eltype(first_state) # e.g. Dual numbers under ForwardDiff

    all_means = Vector{T}[]
    for i in eachindex(accumulated)
        accumulator_index = n_physical_states + i # position in the state vector
        means = Vector{T}(undef, n_windows)
        for k in 1:n_windows
            window_start, window_end = sol.t[k], sol.t[k+1]
            # The accumulator is the time integral of the diagnostic since t_span[1]
            integral_over_window = sol.u[k+1][accumulator_index] - sol.u[k][accumulator_index]
            means[k] = integral_over_window / (window_end - window_start)
        end
        push!(all_means, means)
    end
    return NamedTuple{accumulated}(Tuple(all_means))
end
window_means(model::ProcessBasedModel) = window_means(model.sol, model.accumulated)

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
    G = ground_heat_flux(SantanelloFriedl03(), R_ns, w_1, w_sat, t_sol, thresholds)
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
    β = soil_evaporation_efficiency(Pielke92(), w_1, w_fc, thresholds)
    r_ss = beta_to_r_ss(β, r_as)

    # Turbulent fluxes calculations
    # First solve: detects condensation on the canopy, with the sign of the canopy flux at
    # zero surface resistance (r_sc only enters the denominator of Penman-Monteith)
    λE_tot1, _ = total_evaporation(
        T_a, p_a, VPD_a, A, A_c, A_s, r_aa, r_ac, r_as, r_sc, r_ss, f_wet
    )
    VPD_m1 = vpd_veg_source_height(VPD_a, T_a, p_a, A, λE_tot1, r_aa)
    _, λE_pc1 = penman_monteith(T_a, p_a, VPD_m1, A_c, r_ac, of_value_type(f_wet, 0))
    g_cond = step_function(thresholds, λE_pc1, zero(λE_pc1), flux_scale(thresholds))
    f_wet_eff = fraction_wet_vegetation_with_condensation(f_wet, λE_pc1, thresholds)
    # Second solve, always evaluated (avoid branching): equal to the first when f_wet_eff = f_wet
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
        f_wet_eff,
    )
    VPD_m = vpd_veg_source_height(
        VPD_a, T_a, p_a, A, λE_tot, r_aa
    )
    E_t, λE_t = transpiration(
        T_a, p_a, VPD_m, A_c, r_ac, r_sc, f_wet_eff
    )
    E_i, λE_i = interception_loss(T_a, p_a, VPD_m, A_c, r_ac, f_wet_eff)
    E_s, λE_s = soil_evaporation(T_a, p_a, VPD_m, A_s, r_as, r_ss)

    # Sensible heat via energy balance
    H_tot = sensible_heat(A, λE_tot)
    H_s = sensible_heat_soil(A_s, λE_s)
    H_c = sensible_heat_canopy(A_c, λE_t + λE_i)

    # Remaining water fluxes
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
        H_tot=H_tot,
        H_s=H_s,
        H_c=H_c,
        R_nc=R_nc,
        R_ns=R_ns,
        G=G,
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

"""
    compute_tendencies!(du, u, p, t, forcings, thresholds=HardThresholds(), accumulated=Val(()))

Tendencies of `[w_1, w_2, w_r]`, followed by the `accumulated` diagnostics (a `Val` of names).
"""
function compute_tendencies!(
    du,
    u,
    p::AbstractArray,
    t,
    forcings::NamedTuple,
    thresholds::ThresholdTreatment=HardThresholds(),
    accumulated::Val=Val(()),
)
    diagnostics = compute_diagnostics(u, p, t, forcings, thresholds)
    @unpack d_1, d_2 = p
    @unpack P_c, D_c, I_s, f_1, D_1, K_2, E_s, E_t, E_i, C_1 = diagnostics
    du[1] = C_1 / (ρ_w * d_1) * (f_1 * I_s - E_s) - D_1
    du[2] = 1 / (ρ_w * d_2) * (I_s - E_s - E_t) - K_2
    du[3] = P_c - E_i - D_c
    accumulate_diagnostics!(du, diagnostics, accumulated)
    return nothing
end

"""
    accumulate_diagnostics!(du, diagnostics, accumulated::Val{names})

Tendencies of the accumulator states: `du[3 + i]` is the diagnostic `names[i]`.
"""
@inline function accumulate_diagnostics!(
    du, diagnostics, accumulated::Val{names}
) where {names}
    accumulated_diagnostics = NamedTuple{names}(diagnostics)
    accumulated_values = Tuple(accumulated_diagnostics)
    n_accumulators = length(names)

    ntuple(Val(n_accumulators)) do i
        accumulator_index = N_STATES + i
        du[accumulator_index] = accumulated_values[i]
        return nothing
    end
    return nothing
end
