"""
    ModelParameters(θ, forcing_table, window_starts, window=1)

The `p` of the model's `ODEProblem`: the model parameters `θ` (field `tunable`), the forcing
value of every window ([`forcing_table`](@ref)) and the index of the current window, kept up
to date by [`window_callback`](@ref).

A [SciMLStructures.jl](https://docs.sciml.ai/SciMLStructures/stable/) structure, like
`InversionBinder` in
[ODINN.jl](https://github.com/ODINN-SciML/ODINN.jl/blob/4a037e1e7004e0ba7a90fd6e1109fdad36949d4a/src/simulations/inversions/sciml_utils.jl).
`θ` is the `Tunable` portion; set it with [`with_parameters`](@ref). The field has to be
called `tunable` (as in ModelingToolkit's `MTKParameters`): DiffEqBase's Enzyme rule for
`solve` only adds the gradient to a field with that name, otherwise reverse-mode gradients
are zero. The window index is the `Discrete` portion, so the adjoint solve restores it at
every window change.
"""
mutable struct ModelParameters{P<:AbstractVector,F<:NamedTuple,B<:AbstractVector}
    tunable::P
    forcing_table::F
    window_starts::B
    window::Int
end

function ModelParameters(θ, forcing_table, window_starts)
    return ModelParameters(θ, forcing_table, window_starts, 1)
end

SciMLStructures.isscimlstructure(::ModelParameters) = true
SciMLStructures.ismutablescimlstructure(::ModelParameters) = true
SciMLStructures.hasportion(::SciMLStructures.Tunable, ::ModelParameters) = true
SciMLStructures.hasportion(::SciMLStructures.Discrete, ::ModelParameters) = true

function SciMLStructures.canonicalize(::SciMLStructures.Tunable, p::ModelParameters)
    buffer = Vector(p.tunable)
    repack(new_θ) = SciMLStructures.replace(SciMLStructures.Tunable(), p, new_θ)
    aliases = false # `buffer` is a copy of θ
    return buffer, repack, aliases
end

# Every copy gets its own window index, so that copies can be solved at the same time
function SciMLStructures.replace(::SciMLStructures.Tunable, p::ModelParameters, new_θ)
    θ = similar(p.tunable, eltype(new_θ)) # keeps the names of a ComponentArray
    copyto!(θ, new_θ)
    return ModelParameters(θ, p.forcing_table, p.window_starts, p.window)
end

function SciMLStructures.replace!(::SciMLStructures.Tunable, p::ModelParameters, new_θ)
    p.tunable .= new_θ
    return p
end

# The window index, as a vector of one element
function SciMLStructures.canonicalize(::SciMLStructures.Discrete, p::ModelParameters)
    buffer = [p.window]
    repack(new_window) = SciMLStructures.replace(SciMLStructures.Discrete(), p, new_window)
    aliases = false
    return buffer, repack, aliases
end

function SciMLStructures.replace(::SciMLStructures.Discrete, p::ModelParameters, new_window)
    window = only(new_window)
    return ModelParameters(p.tunable, p.forcing_table, p.window_starts, window)
end

function SciMLStructures.replace!(::SciMLStructures.Discrete, p::ModelParameters, new_window)
    p.window = only(new_window)
    return p
end

# Used by SciMLSensitivity to differentiate through callbacks (as in ODINN.jl)
function Base.zero(p::ModelParameters)
    return ModelParameters(zero(p.tunable), p.forcing_table, p.window_starts, p.window)
end

function Base.copyto!(dest::ModelParameters, src::ModelParameters)
    dest.tunable .= src.tunable
    return dest
end

"""
    with_parameters(prob, θ)

`prob` with the model parameters in its [`ModelParameters`](@ref) replaced by `θ`.
"""
function with_parameters(prob, θ)
    p = SciMLStructures.replace(SciMLStructures.Tunable(), prob.p, θ)
    return remake(prob; p)
end

"""
    forcing_window_starts(t_span, tstops)

Start times of the forcing windows: `t_span[1]` and the window bounds `tstops` before
`t_span[2]`.
"""
function forcing_window_starts(t_span, tstops)
    t_start, t_end = t_span
    bounds = sort(unique([t_start; tstops]))
    is_window_start(t) = t_start <= t < t_end
    return filter(is_window_start, bounds)
end

"""
    forcing_table(forcings, window_starts)

Value of each forcing in each window: a vector with one value per window, or a number for
a constant forcing (and `lon`, `utc_offset`).
"""
function forcing_table(forcings::NamedTuple, window_starts)
    return map(forcing -> window_values(forcing, window_starts), forcings)
end

# The values of one forcing in all windows
window_values(constant::Number, window_starts) = constant
function window_values(forcing, window_starts)
    throw(ArgumentError("a forcing must be a vector or a number, not $(typeof(forcing))"))
end
function window_values(forcing::AbstractVector, window_starts)
    n_values = length(forcing)
    n_windows = length(window_starts)
    if n_values != n_windows
        throw(ArgumentError("a forcing vector has $n_values values for $n_windows windows"))
    end
    return forcing
end

"""
    window_index(window_starts, t)

Index of the window containing `t`; at a bound, the window starting there. Clamped to the
first and last window.
"""
function window_index(window_starts, t)
    window = searchsortedlast(window_starts, t) # 0 before the first start
    n_windows = length(window_starts)
    return clamp(window, 1, n_windows)
end

"""
    window_forcing(forcing_table, window)

Forcing values of window number `window`, as passed to [`compute_diagnostics`](@ref).
"""
function window_forcing(forcing_table::NamedTuple, window)
    return map(values -> value_in_window(values, window), forcing_table)
end

# One forcing in window number `window`
value_in_window(constant::Number, window) = constant
value_in_window(values::AbstractVector, window) = values[window]

# The forcing of the current window of `p`, used by the right-hand side
current_forcing(p::ModelParameters) = window_forcing(p.forcing_table, p.window)

"""
    window_callback()

Keeps the window index of the [`ModelParameters`](@ref) up to date. It moves to the next
window only after the step ending at a window bound is accepted, so every stage of that step
still uses the old window's forcing. At the start of the solve it adds all window starts as
`tstops`, so no step crosses a bound.

A window change counts as a derivative discontinuity (the default for a callback), so FSAL
methods like `Tsit5` recompute their first stage with the new forcing and SciMLSensitivity
records the event. Ribasim does the same in
[`update_basin!`](https://github.com/Deltares/Ribasim/blob/642eccef4ed4d1873584eb705e43cc6e93afb252/core/src/callback.jl#L849-L865).
Added to `model.prob` by [`initialize!`](@ref).
"""
function window_callback()
    function at_next_window_start(u, t, integrator)
        p = integrator.p
        n_windows = length(p.window_starts)
        if p.window == n_windows # no next window
            return false
        end
        next_window_start = p.window_starts[p.window + 1]
        return t == next_window_start
    end
    # Also replayed by SciMLSensitivity on a stand-in integrator, which only has `u`, `p`, `t`
    function next_window!(integrator)
        integrator.p.window += 1
        return nothing
    end
    function first_window(callback, u, t, integrator)
        p = integrator.p
        p.window = window_index(p.window_starts, t)
        later_window_starts = p.window_starts[(p.window + 1):end]
        for window_start in later_window_starts
            add_tstop!(integrator, window_start)
        end
        derivative_discontinuity!(integrator, false) # nothing changed: no re-evaluation
        return nothing
    end
    return DiscreteCallback(
        at_next_window_start, next_window!;
        initialize=first_window, save_positions=(false, false),
    )
end
