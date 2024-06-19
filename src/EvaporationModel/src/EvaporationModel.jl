module EvaporationModel

using ComponentArrays
using Parameters

include("resistances.jl")
export jarvis_stewart, aerodynamic_resistance

greet() = print("Hello World!")

end # module EvaporationModel
