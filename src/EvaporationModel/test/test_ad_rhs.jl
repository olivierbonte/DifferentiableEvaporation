# AD of the right-hand side of the ODE system (no solver involved)
function rhs(u, thresholds)
    du = similar(u)
    compute_tendencies!(du, u, p_model, t_unix, forcings_toy, thresholds)
    return du
end

@testset "Jacobian of the tendencies" begin
    backends = (
        AutoForwardDiff(),
        AutoEnzyme(; mode=Enzyme.Forward),
        # Without runtime activity, reverse mode silently gives wrong rows for du[1:2]
        AutoEnzyme(; mode=Enzyme.set_runtime_activity(Enzyme.Reverse)),
    )
    for thresholds in treatments
        f = u -> rhs(u, thresholds)
        # Check ForwardDiff (AD) against a finite difference Jacobian, as e.g. in
        # https://mc-stan.org/math/md_doxygen_2contributor__help__pages_2autodiff__test__guide.html
        J_fd = jacobian(f, AutoFiniteDiff(), u0_toy)
        J_fwd = jacobian(f, AutoForwardDiff(), u0_toy)
        @test isapprox(J_fwd, J_fd; rtol=1e-5)
        # Given that ForwardDiff AD is correct, check that other AD backends give the same Jacobian
        for backend in backends
            J = jacobian(f, backend, u0_toy)
            @test all(isfinite, J)
            @test isapprox(J, J_fwd; rtol=1e-10)
        end
    end
end
