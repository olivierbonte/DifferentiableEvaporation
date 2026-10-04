# Import Packages to test
using Bigleaf
using ComponentArrays
using Dates
using EvaporationModel
# Import packages used for testing
using AllocCheck
using BenchmarkTools
using DifferentiationInterface
using Documenter
import ForwardDiff # backend for DifferentiationInterface
using OrdinaryDiffEq
using Test
# Here you include files using `srcdir`
# include(srcdir("file.jl"))

DocMeta.setdocmeta!(EvaporationModel, :DocTestSetup, :(using Dates))

# Run test suite
println("Starting tests")
ti = time()

println("Defining test inputs")
r_aa = 50.0 # s/m
r_ac = 30.0 # s/m
r_as = 100.0 # s/m
r_sc = 1000.0 # s/m
r_ss = 100.0 # s/m
f_wet = 0.1
f_veg = 0.5
SW_in = 1000.0 # W/m2
w_2 = 0.2
w_fc = 0.3
w_wp = 0.1
g_d = 0.003
r_smin = 200.0 # s/m
Rn = 400.0 # W/m2
G = 50.0 # W/m2
T_a = 300.0 # K
p_a = 101325.0 # Pa
VPD_a = 2000.0 # Pa
LAI = 3.0 # m2/m2
h = 21.0 # m
d_c = 2 / 3 * h # m
z_0mc = 0.1 * h # m
z_obs = 39.0 # m
z_0ms = 0.01 # m
u_star = 3.0 # m/s
u = 5.0 # m/s
P = 2.0e-5 # gross precipitation, kg / (m2 * s)
P_s = 1.5e-5 # precipitation below the canopy, kg / (m2 * s)
w_sat = 0.45

p_model = ComponentArray(;
    h=h, z_0ms=z_0ms, w_sat=w_sat, a=0.15, p_soil=6.0, b=6.1, w_res=0.04,
    w_wp=w_wp, w_fc=w_fc, C_1sat=0.019, C_2ref=0.83, C_3=0.25, d_1=0.01, d_2=1.3,
    z_obs=z_obs, kB⁻¹=log(10), g_d=3e-4, r_smin=395.0, k_ext=0.5,
)
constant_forcings(P) = (
    P=t -> P, T_a=t -> T_a, u_a=t -> u, p_a=t -> p_a, VPD_a=t -> VPD_a,
    SW_in=t -> SW_in, R_n=t -> Rn, LAI=t -> LAI, lon=4.52, utc_offset=1,
)
t_unix = datetime2unix(DateTime(2010, 7, 1, 10)) # model time [s since Unix epoch], 10:00 local
treatments = (HardThresholds(), KavetskiSmoothing())

@testset "Mass conservation of tendencies" begin
    forcings = constant_forcings(P)
    for thresholds in treatments,
        u0 in ([0.2, 0.32, 0.3], [0.05, 0.1, 0.0], [0.4, 0.44, 0.7], [w_sat, w_sat, 0.7])

        du = zeros(3)
        compute_tendencies!(du, u0, p_model, t_unix, forcings, thresholds)
        d = compute_diagnostics(u0, p_model, t_unix, forcings, thresholds)
        @test all(isfinite, du)
        @test d.λE_tot ≈ d.λE_t + d.λE_i + d.λE_s
        # Water balance of the column, d/dt(ρ_w d_2 w_2 + w_r) = inputs - outputs, with both
        # sides in kg m⁻² s⁻¹ (mm s⁻¹).
        lhs = ρ_w * p_model.d_2 * du[2] + du[3]
        rhs = P - d.Q_s - d.E_s - d.E_t - d.E_i - ρ_w * p_model.d_2 * d.K_2
        @test lhs ≈ rhs atol = 1e-12
    end
end

@testset "Threshold treatments" begin
    hard, smooth = treatments
    @test max(hard, 0.2, 0.3, moisture_scale(hard)) == 0.3
    @test min(hard, 0.2, 0.3, moisture_scale(hard)) == 0.2
    @test clamp(hard, 1.5, 0.0, 1.0, factor_scale(hard)) == 1.0
    # The smooth clamp never drops below the lower bound, and exceeds the upper by < s/2
    s = factor_scale(smooth)
    for x in (-1.0, 0.0, 0.5, 1.0, 2.0)
        @test 0 <= clamp(smooth, x, 0.0, 1.0, s) <= 1 + s / 2
    end
    # Smoothed f_wet has a finite slope at w_r = 0 (the hard 2/3 power does not)
    @test isfinite(
        derivative(w_r -> fraction_wet_vegetation(w_r, 0.6, smooth), AutoForwardDiff(), 0.0)
    )
    # HardThresholds is the default treatment
    forcings = constant_forcings(P)
    u0 = [0.2, 0.32, 0.3]
    @test compute_diagnostics(u0, p_model, t_unix, forcings) ==
        compute_diagnostics(u0, p_model, t_unix, forcings, hard)
end

@testset "Surface layer bounded at saturation" begin
    @test surface_infiltration_factor(w_sat, w_sat) == 0
    @test surface_infiltration_factor(0.0, w_sat) ≈ 1
    @test surface_infiltration_factor(w_sat - 0.05, w_sat) > 0.99
    @test surface_infiltration_factor(w_sat + 0.01, w_sat) < 0
    # Heavy rain (20 mm/h) on a saturated surface layer during the day: no further wetting
    P_heavy = 20 / 3600 # kg / (m2 * s)
    forcings = constant_forcings(P_heavy)
    for thresholds in treatments
        du = zeros(3)
        compute_tendencies!(du, [w_sat, 0.4, 0.7], p_model, t_unix, forcings, thresholds)
        @test du[1] <= 0
        # 12 h of heavy rain on a wet surface over a drier root zone (so most rain
        # infiltrates) drives w_1 to the bound, but not above it
        t_span = (t_unix, t_unix + 12 * 3600.0)
        model = ProcessBasedModel{Float64}(;
            forcings=forcings,
            parameters=p_model,
            t_span=t_span,
            u0=[0.35, 0.25, 0.5],
            saveat=collect(t_span[1]:600.0:t_span[2]),
            thresholds=thresholds,
        )
        EvaporationModel.initialize!(model)
        # Explicit solver: the ForwardDiff Jacobian is NaN for w_1 > w_fc, where r_ss = 0
        EvaporationModel.solve!(model; AD=true, alg=Tsit5(), abstol=1e-8, reltol=1e-8)
        @test OrdinaryDiffEq.SciMLBase.successful_retcode(model.sol)
        @test maximum(u -> u[1], model.sol.u) <= w_sat + 1e-6
        @test maximum(u -> u[2], model.sol.u) <= w_sat + 1e-6
        @test maximum(u -> u[1], model.sol.u) > w_sat - 0.01 # the bound is active
    end
end

@testset "Check evaporation sum" begin
    R_nc, R_ns = net_radiation_partitioning(Rn, f_veg)
    A, A_c, A_s = available_energy_partitioning(R_nc, R_ns, G)
    λE_tot, λE_tot_p = total_evaporation(
        T_a, p_a, VPD_a, A, A_c, A_s, r_aa, r_ac, r_as, r_sc, r_ss, f_wet
    )
    VPD_m = vpd_veg_source_height(VPD_a, T_a, p_a, A, λE_tot, r_aa)
    ET_t, λE_t = transpiration(T_a, p_a, VPD_m, A_c, r_ac, r_sc, f_wet)
    ET_i, λE_i = interception_loss(T_a, p_a, VPD_m, A_c, r_ac, f_wet)
    ET_s, λE_s = soil_evaporation(T_a, p_a, VPD_m, A_s, r_as, r_ss)
    @test λE_tot ≈ λE_t + λE_i + λE_s
end

# Test idea from https://modernjuliaworkflows.org/optimizing/#memory_management
@testset "Check functions on having no allocations" begin
    FT = Float64
    # @test (@ballocations fractional_vegetation_cover($LAI)) == 0
    println("Testing functions from canopy.jl")
    @test isempty(check_allocs(fractional_vegetation_cover, (FT,)))
    @test isempty(check_allocs(net_radiation_partitioning, (FT, FT)))
    @test isempty(check_allocs(available_energy_partitioning, (FT, FT, FT)))
    @test isempty(check_allocs(fraction_wet_vegetation, (FT, FT)))
    @test isempty(check_allocs(canopy_input, (FT, FT)))
    @test isempty(check_allocs(canopy_drainage, (FT, FT, FT, FT)))
    @test isempty(check_allocs(precip_below_canopy, (FT, FT, FT)))
    @test isempty(check_allocs(vpd_veg_source_height, (FT, FT, FT, FT, FT, FT)))

    println("Testing function from evaporation.jl")
    @test isempty(check_allocs(penman_monteith, (FT, FT, FT, FT, FT, FT)))
    @test isempty(
        check_allocs(total_evaporation, (FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT))
    )
    @test isempty(check_allocs(transpiration, (FT, FT, FT, FT, FT, FT, FT)))
    @test isempty(check_allocs(interception_loss, (FT, FT, FT, FT, FT, FT)))
    @test isempty(check_allocs(soil_evaporation, (FT, FT, FT, FT, FT, FT)))

    println("Testing functions from resistances.jl")
    @test isempty(check_allocs(ustar_from_u, (FT, FT, FT, FT)))
    @test (@ballocations soil_aerodynamic_resistance(
        Choudhury1988soil(), $u_star, $h, $d_c, $z_0mc, $z_0ms
    )) == 0
    @test (@ballocations surface_resistance(
        JarvisStewart(), $SW_in, $VPD_a, $T_a, $w_2, $w_fc, $w_wp, $LAI, $g_d, $r_smin
    )) == 0

    println("Testing functions from soil_fluxes.jl")
    @test (@ballocations surface_runoff(StaticInfiltration(), $P_s, $w_2, $w_fc)) == 0
    @test (@ballocations surface_runoff(
        VegetationInfiltration(), $P_s, $w_2, $w_fc, $f_veg
    )) == 0
    @test isempty(check_allocs(c_1, (FT, FT, FT, FT, FT)))
    @test isempty(check_allocs(surface_infiltration_factor, (FT, FT)))
    @test isempty(check_allocs(diffusion_layer_1, (FT, FT, FT)))
    @test isempty(check_allocs(vertical_drainage_layer_2, (FT, FT, FT, FT)))

    println("Testing Bigleaf functions")
    @test (@ballocations Bigleaf.roughness_parameters(
        RoughnessCanopyHeightLAI(), $h, $LAI; hs=$z_0ms
    )) == 0
    rough_dict = Bigleaf.roughness_parameters(RoughnessCanopyHeightLAI(), h, LAI; hs=z_0ms)
    @test (@ballocations Bigleaf.compute_Ram(ResistanceWindZr(), $u_star, $u)) == 0
    @test isempty(check_allocs(Bigleaf.Gb_constant_kB1, (FT, FT)))
end

@testset "Check test from doctest" begin
    Documenter.doctest(EvaporationModel)
end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti / 60; digits=3), " minutes")
