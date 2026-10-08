module EvaporationModel

using ADTypes: AutoEnzyme, AutoForwardDiff
using Bigleaf
using ComponentArrays
using Dates
using DataFrames
using DiffEqCallbacks
using Enzyme: Enzyme
using ForwardDiff
using OrdinaryDiffEq
using OrdinaryDiffEqLowOrderRK: Euler, Heun
using OrdinaryDiffEqSDIRK: ImplicitEuler
using Parameters
using YAXArrays

# Solvers and AD backends for the Jacobian of implicit solvers
export Tsit5, Euler, Heun, ImplicitEuler, AutoForwardDiff, AutoEnzyme

include("thresholds.jl")
export smooth_min,
    smooth_max,
    smooth_clamp,
    smoothing_kernel,
    KernelMethod,
    LowerBound,
    UpperBound,
    ThresholdTreatment,
    HardThresholds,
    KavetskiSmoothing,
    lower_bound_kernel,
    moisture_scale,
    storage_scale,
    factor_scale

include("canopy.jl")
export fractional_vegetation_cover,
    net_radiation_partitioning,
    available_energy_partitioning,
    max_canopy_capacity,
    fraction_wet_vegetation,
    canopy_drainage,
    precip_below_canopy,
    canopy_input,
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
export penman_monteith, total_evaporation, transpiration, interception_loss, soil_evaporation

include("ground_heat_flux.jl")
export ground_heat_flux,
    compute_harmonic_sum, GroundHeatFluxMethod, Allen07, SantanelloFriedl03

include("model.jl")
export AbstractModel,
    ProcessBasedModel, compute_diagnostics, compute_tendencies!, initialize!, solve!

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
export c_1,
    c_2,
    c_3,
    c_1sat,
    C1satMethod,
    NoilhanMahfouf96,
    NoilhanPlanton89,
    c_2ref,
    w_geq,
    compute_a,
    compute_b,
    compute_p,
    bMethod,
    Clay,
    VanGenuchten

include("soil_fluxes.jl")
export surface_runoff,
    surface_infiltration_factor,
    diffusion_layer_1,
    vertical_drainage_layer_2,
    InfiltrationMethod,
    StaticInfiltration,
    VegetationInfiltration

include("utils.jl")
export compute_amplitude_and_phase,
    fourier_series,
    fit_fourier_coefficients,
    local_to_solar_time,
    seconds_since_solar_noon,
    value_type,
    of_value_type

end # module EvaporationModel
