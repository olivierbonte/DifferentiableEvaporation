module EvaporationModel

using Bigleaf
using ComponentArrays
using Dates
using LinearSolve
using Parameters

include("canopy.jl")
export fractional_vegetation_cover,
    available_energy_partioning,
    fraction_wet_vegetation,
    canopy_drainage,
    vpd_veg_source_height

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
    transpiration,
    interception,
    bare_soil_evaporation

include("ground_heat_flux.jl")
export ground_heat_flux,
    compute_harmonic_sum, GroundHeatFluxMethod, Allen07, SantanelloFriedl03

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

include("soil_fluxes.jl")
export surface_runoff,
    diffusion_layer_1,
    vertical_drainage_layer_2,
    precip_below_canopy,
    InfiltrationMethod,
    StaticInfiltration,
    VegetationInfiltration

include("utils.jl")
export compute_amplitude_and_phase,
    fourier_series, fit_fourier_coefficients, local_to_solar_time, seconds_since_solar_noon

end # module EvaporationModel
