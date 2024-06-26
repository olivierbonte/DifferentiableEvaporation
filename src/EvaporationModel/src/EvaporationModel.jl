module EvaporationModel

using ComponentArrays
using Parameters

include("resistances.jl")
export jarvis_stewart, aerodynamic_resistance

include("soil.jl")
export compute_c_1, compute_c_2, compute_w_geq

include("constants.jl")
export ρ_w

greet() = print("Hello World!")

end # module EvaporationModel
