# Import Packages to test
using Bigleaf
using EvaporationModel

# Import AD backend that are in EvaporationModel
using EvaporationModel: Enzyme
using EvaporationModel: ForwardDiff

# Import packages used for testing
using AllocCheck
using BenchmarkTools
using ComponentArrays
using Dates
using DifferentiationInterface
using Documenter
using FiniteDiff: FiniteDiff # backend for DifferentiationInterface, reference
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
