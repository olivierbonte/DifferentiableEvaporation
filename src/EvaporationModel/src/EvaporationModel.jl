module EvaporationModel

using ComponentArrays
using Parameters

include("resistances.jl")
export jarvis_stewart, aerodynamic_resistance

include("soil.jl")
export C_1, C_2, w_geq

greet() = print("Hello World!")

end # module EvaporationModel
