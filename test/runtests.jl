using DrWatson, Test
@quickactivate "DifferentiableEvaporation"
using EvaporationModel
using AllocCheck
# Here you include files using `srcdir`
# include(srcdir("file.jl"))

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
Rn = 400.0 # W/m2
G = 50.0 # W/m2
T_a = 300.0 # K
p_a = 101325.0 # Pa
VPD_a = 2000.0 # Pa
LAI = 3.0 # m2/m2

@testset "Check evaporation sum" begin
    A, A_c, A_s = available_energy_partioning(Rn, G, f_veg)
    λE_tot, λE_tot_p = total_evaporation(
        T_a, p_a, VPD_a, A, A_c, A_s, r_aa, r_ac, r_as, r_sc, r_ss, f_wet
    )
    VPD_m = vpd_veg_source_height(VPD_a, T_a, p_a, A, λE_tot, r_aa)
    ET_t, λE_t = transpiration(T_a, p_a, VPD_m, A_c, r_ac, r_sc, f_wet)
    ET_i, λE_i = interception(T_a, p_a, VPD_m, A_c, r_ac, f_wet)
    ET_s, λE_s = soil_evaporation(T_a, p_a, VPD_m, A_s, r_as, r_ss)
    @test λE_tot ≈ λE_t + λE_i + λE_s
end

# Test idea from https://modernjuliaworkflows.org/optimizing/#memory_management
@testset "Check functions on having no allocations" begin
    # @test (@ballocations fractional_vegetation_cover($LAI)) == 0
    println("Testing functions from canopy.jl")
    @test isempty(check_allocs(fractional_vegetation_cover, (Float64,)))
    @test isempty(check_allocs(available_energy_partioning, (Float64, Float64, Float64)))
    @test isempty(check_allocs(fraction_wet_vegetation, (Float64, Float64)))
    @test isempty(check_allocs(canopy_drainage, (Float64, Float64, Float64)))
    @test isempty(check_allocs(precip_below_canopy, (Float64, Float64, Float64)))
    @test isempty(
        check_allocs(
            vpd_veg_source_height, (Float64, Float64, Float64, Float64, Float64, Float64)
        ),
    )
end
ti = time() - ti
println("\nTest took total time of:")
println(round(ti / 60; digits=3), " minutes")
