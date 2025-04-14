module EvaporationModel

using Bigleaf
using ComponentArrays
using LinearSolve
using Parameters

include("config.jl")
export VegetationParameters, VegetationType, Crops, ShortGrass, EvergreenNeedleleafTrees,
    DeciduousNeedleleafTrees, DeciduousBroadleafTrees, EvergreenBroadleafTrees,
    TallGrass, Desert, Tundra, IrrigatedCrops, Semidesert, BogsAndMarshes,
    EvergreenShrubs, DeciduousShrubs, MixedForestWoodland, InterruptedForest,
    WaterAndLandMixtures

include("constants.jl")
export ρ_w

include("evaporation.jl")
export penman_monteith,
    compute_soil_evaporation_stress,
    compute_bare_soil_evaporation,
    compute_transpiration

include("groundheatlux.jl")
export compute_g_from_r_net, compute_harmonic_sum

include("resistances.jl")
export jarvis_stewart, aerodynamic_resistance

include("soil.jl")
export compute_c_1, compute_c_2, compute_c_3, compute_w_geq, compute_a, compute_b, compute_p

include("utils.jl")
export compute_amplitude_and_phase, fourier_series, fit_fourier_coefficients

end # module EvaporationModel
