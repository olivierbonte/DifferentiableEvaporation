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

    println("Testing the model equations")
    forcings = constant_forcings(P)
    u0 = [0.2, 0.32, 0.3]
    du = zeros(3)
    for thresholds in treatments
        @test (@ballocations compute_diagnostics(
            $u0, $p_model, $t_unix, $forcings, $thresholds
        )) == 0
        @test (@ballocations compute_tendencies!(
            $du, $u0, $p_model, $t_unix, $forcings, $thresholds
        )) == 0
    end

    println("Testing Bigleaf functions")
    @test (@ballocations Bigleaf.roughness_parameters(
        RoughnessCanopyHeightLAI(), $h, $LAI; hs=($z_0ms)
    )) == 0
    rough_dict = Bigleaf.roughness_parameters(RoughnessCanopyHeightLAI(), h, LAI; hs=z_0ms)
    @test (@ballocations Bigleaf.compute_Ram(ResistanceWindZr(), $u_star, $u)) == 0
    @test isempty(check_allocs(Bigleaf.Gb_constant_kB1, (FT, FT)))
end
