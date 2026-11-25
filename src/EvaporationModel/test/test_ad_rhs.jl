# AD of the right-hand side of the ODE system (no solver involved)
function rhs(u, thresholds, accumulated=Val(()))
    du = similar(u)
    compute_tendencies!(du, u, p_model, t_unix, constant_forcing_toy, thresholds, accumulated)
    return du
end

@testset "Jacobian of the tendencies" begin
    backends = (
        AutoForwardDiff(),
        AutoEnzyme(; mode=Enzyme.Forward),
        # Without runtime activity, reverse mode silently gives wrong rows for du[1:2]
        AutoEnzyme(; mode=Enzyme.set_runtime_activity(Enzyme.Reverse)),
    )
    # Without and with accumulators of the main energy fluxes
    accumulated_cases = ((), (:λE_tot, :H_tot, :G))
    for thresholds in treatments, accumulated in accumulated_cases
        n_accumulators = length(accumulated)
        u = [u0_toy; zeros(n_accumulators)]
        f = u -> rhs(u, thresholds, Val(accumulated))
        # Check ForwardDiff (AD) against a finite difference Jacobian, as e.g. in
        # https://mc-stan.org/math/md_doxygen_2contributor__help__pages_2autodiff__test__guide.html
        J_fd = jacobian(f, AutoFiniteDiff(), u)
        J_fwd = jacobian(f, AutoForwardDiff(), u)
        @test isapprox(J_fwd, J_fd; rtol=1e-5)
        # Accumulators do not feed back into the tendencies: zero columns
        @test all(iszero, J_fwd[:, 4:end])
        # Given that ForwardDiff AD is correct, check that other AD backends give the same Jacobian
        for backend in backends
            J = jacobian(f, backend, u)
            @test all(isfinite, J)
            @test isapprox(J, J_fwd; rtol=1e-10)
        end
    end
end
