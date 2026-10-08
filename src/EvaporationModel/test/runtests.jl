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
using Enzyme: Enzyme # backends for DifferentiationInterface
using FiniteDiff: FiniteDiff
using ForwardDiff: ForwardDiff
using OrdinaryDiffEq
using SciMLSensitivity
using Test

DocMeta.setdocmeta!(EvaporationModel, :DocTestSetup, :(using Dates))

# Run test suite
println("Starting tests")
ti = time()

include("common.jl")

@testset "EvaporationModel" verbose = true begin
    include("test_processes.jl")
    include("test_model.jl")
    include("test_allocations.jl")
    include("test_ad_rhs.jl")
    include("test_ad_solver.jl")
    include("test_sensitivity.jl")
    include("test_doctests.jl")
end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti / 60; digits=3), " minutes")
