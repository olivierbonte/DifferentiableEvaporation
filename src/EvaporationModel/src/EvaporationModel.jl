module EvaporationModel

using Bigleaf
using ComponentArrays
using LinearSolve
using Parameters

include("config.jl")
export VegetationParameters,
    VegetationType,
    Crops,
    ShortGrass,
    EvergreenNeedleleafTrees,
    DeciduousNeedleleafTrees,
    DeciduousBroadleafTrees,
    EvergreenBroadleafTrees,
    TallGrass,
    Desert,
    Tundra,
    IrrigatedCrops,
    SemiDesert,
    BogsAndMarshes,
    EvergreenShrubs,
    DeciduousShrubs,
    MixedForestWoodland,
    InterruptedForest,
    WaterAndLandMixtures

include("constants.jl")
export ρ_w, τ

include("evaporation.jl")
export penman_monteith,
    total_evaporation,
    compute_soil_evaporation_stress,
    compute_bare_soil_evaporation,
    compute_transpiration

include("groundheatlux.jl")
export compute_g_from_r_net, compute_harmonic_sum

include("resistances.jl")
export jarvis_stewart,
    aerodynamic_resistance,
    soil_aerodynamic_resistance,
    ustar_from_u,
    surface_resistance,
    soil_evaporation_efficiency,
    beta_to_r_ss,
    r_ss_to_beta,
    SurfaceResistanceMethod,
    JarvisStewart,
    SoilAerodynamicResistanceMethod,
    Choudhury1988soil,
    SoilEvaporationEfficiencyMethod,
    Martens17,
    Pielke92

include("soil.jl")
export c_1, c_2, c_3, c_1sat, c_2ref, w_geq, compute_a, compute_b, compute_p

include("utils.jl")
export compute_amplitude_and_phase, fourier_series, fit_fourier_coefficients

end # module EvaporationModel
