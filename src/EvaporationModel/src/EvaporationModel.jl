module EvaporationModel

using ComponentArrays
using Parameters

include("constants.jl")
export ρ_w

include("evaporation.jl")
export penman_monteith,
    compute_g_from_r_net,
    compute_soil_evaporation_stress,
    compute_bare_soil_evaporation,
    compute_transpiration

include("resistances.jl")
export jarvis_stewart, aerodynamic_resistance

include("soil.jl")
export compute_c_1, compute_c_2, compute_c_3, compute_w_geq, compute_b

end # module EvaporationModel
