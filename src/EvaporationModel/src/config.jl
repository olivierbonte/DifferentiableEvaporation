## Vegetation parameters
# Define vegeteation types as Julia structs
"""
    VegetationType 

Abstract type defined to represent different types of vegetation.

Below example provided on how to check all the subtypes (which are structs):
```julia
subtypes(VegetationType)
```
"""
abstract type VegetationType end
struct Crops <: VegetationType end
struct ShortGrass <: VegetationType end
struct EvergreenNeedleleafTrees <: VegetationType end
struct DeciduousNeedleleafTrees <: VegetationType end
struct DeciduousBroadleafTrees <: VegetationType end
struct EvergreenBroadleafTrees <: VegetationType end
struct TallGrass <: VegetationType end
struct Desert <: VegetationType end
struct Tundra <: VegetationType end
struct IrrigatedCrops <: VegetationType end
struct SemiDesert <: VegetationType end
struct BogsAndMarshes <: VegetationType end
struct EvergreenShrubs <: VegetationType end
struct DeciduousShrubs <: VegetationType end
struct MixedForestWoodland <: VegetationType end
struct InterruptedForest <: VegetationType end
struct WaterAndLandMixtures <: VegetationType end

"""
    get_default_value(vegtype::VegetationType)

Assign default values for minimal stomatal resistance (`r_smin`) and coefficient relating 
vapour pressure deficit to stomatal resistance (`g_d`) based on vegetation type. 

The values come from Table 8.1 of the 
[IFS Cy49r1 documentation Part IV: Phyiscal processes](https://doi.org/10.21957/c731ee1102).
Not part of public API, only used in [`VegetationParameters`](@ref VegetationParameters).
"""
function get_default_value(vegtype::VegetationType)
    params = Dict(
        Crops() => (125.0, 0.0),
        ShortGrass() => (80.0, 0.0),
        EvergreenNeedleleafTrees() => (395.0, 3e-4),
        DeciduousNeedleleafTrees() => (320.0, 3e-4),
        DeciduousBroadleafTrees() => (215.0, 3e-4),
        EvergreenBroadleafTrees() => (320.0, 3e-4),
        TallGrass() => (100.0, 0.0),
        Desert() => (250.0, 0.0),
        Tundra() => (45.0, 0.0),
        IrrigatedCrops() => (110.0, 0.0),
        SemiDesert() => (45.0, 0.0),
        BogsAndMarshes() => (130.0, 0.0),
        EvergreenShrubs() => (230.0, 0.0),
        DeciduousShrubs() => (110.0, 0.0),
        MixedForestWoodland() => (180.0, 3e-4),
        InterruptedForest() => (175.0, 3e-4),
        WaterAndLandMixtures() => (150.0, 3e-4),
    )
    if haskey(params, vegtype)
        r_smin, g_d = params[vegtype]
        return (r_smin, g_d)
    else
        error("Unknown vegetation type: $vegtype")
    end
end

"""
    VegetationParameters(;...)

Default vegetation parameters for the Jarvis-Stewart model.
The fields are:

- `vegtype`: The type of vegetation (e.g., `Crops`, `ShortGrass`, etc.)
- `r_smin`: Minimum stomatal resistance [s/m]
- `g_d`: Coefficient relating vapour pressure deficit to stomatal resistance [Pa⁻¹]

Based on [`get_default_value`](@ref get_default_value) function, values of `r_smin` and `g_d` 
are set based on `vegtype`. Default values can be overridden by passing them as keyword arguments.

## Examples
```jldoctest
using EvaporationModel
# Get default parameters
params_default = VegetationParameters(vegtype=Crops())
# Adapt a default value
params = VegetationParameters(vegtype=Crops(), r_smin=300.0) 
params.r_smin == 300.0

# output

true
```
"""
@with_kw struct VegetationParameters{FT}
    vegtype::VegetationType
    r_smin::FT = get_default_value(vegtype)[1]
    g_d::FT = get_default_value(vegtype)[2]
end
