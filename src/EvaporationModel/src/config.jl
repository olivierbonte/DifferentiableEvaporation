## Vegetation parameters
# Define vegeteation types as Julia abstract types
"""
    VegetationType 

Abstract type defined to represent different types of vegetation.

Below example provided on how to check all the subtypes:
```julia
subtypes(VegetationType)
```
"""
abstract type VegetationType end
abstract type Crops <: VegetationType end
abstract type ShortGrass <: VegetationType end
abstract type EvergreenNeedleleafTrees <: VegetationType end
abstract type DeciduousNeedleleafTrees <: VegetationType end
abstract type DeciduousBroadleafTrees <: VegetationType end
abstract type EvergreenBroadleafTrees <: VegetationType end
abstract type TallGrass <: VegetationType end
abstract type Desert <: VegetationType end
abstract type Tundra <: VegetationType end
abstract type IrrigatedCrops <: VegetationType end
abstract type Semidesert <: VegetationType end
abstract type BogsAndMarshes <: VegetationType end
abstract type EvergreenShrubs <: VegetationType end
abstract type DeciduousShrubs <: VegetationType end
abstract type MixedForestWoodland <: VegetationType end
abstract type InterruptedForest <: VegetationType end
abstract type WaterAndLandMixtures <: VegetationType end


"""
    get_default_value(vegtype::Type{<:VegetationType})

Assign default values for minimal stomatal resistance (rsmin) and coefficient relating 
vapour pressure deficit to stomatal resistance (gd) based on vegetation type. 

The values come from Table 8.1 of the 
[IFS Cy49r1 documentation Part IV: Phyiscal processes](https://doi.org/10.21957/c731ee1102).
Not part of public API, only used in [`VegetationParameters`](
@ref VegetationParameters).
"""
function get_default_value(vegtype::Type{<:VegetationType})
    params = Dict(
        Crops => (125.0, 0.0),
        ShortGrass => (80.0, 0.0),
        EvergreenNeedleleafTrees => (395.0, 3.0),
        DeciduousNeedleleafTrees => (320.0, 3.0),
        DeciduousBroadleafTrees => (215.0, 3.0),
        EvergreenBroadleafTrees => (320.0, 3.0),
        TallGrass => (100.0, 0.0),
        Desert => (250.0, 0.0),
        Tundra => (45.0, 0.0),
        IrrigatedCrops => (110.0, 0.0),
        Semidesert => (45.0, 0.0),
        BogsAndMarshes => (130.0, 0.0),
        EvergreenShrubs => (230.0, 0.0),
        DeciduousShrubs => (110.0, 0.0),
        MixedForestWoodland => (180.0, 3.0),
        InterruptedForest => (175.0, 3.0),
        WaterAndLandMixtures => (150.0, 3.0),
    )
    if haskey(params, vegtype)
        rsmin, gd = params[vegtype]
        return (rsmin, gd)
    else
        error("Unknown vegetation type: $vegtype")
    end
end

"""
    VegetationParameters(;...)

Default vegetation parameters for the Jarvis-Stewart model.
The fields are:

- `vegtype`: The type of vegetation (e.g., `Crops`, `ShortGrass`, etc.)
- `rsmin`: Minimum stomatal resistance [s/m]
- `gd`: Coefficient relating vapour pressure deficit to stomatal resistance [Pa⁻¹]

Based on [`get_default_value`](@ref get_default_value) function, values of `rsmin` and `gd` 
are set based on `vegtype`. Default values can be overridden by passing them as keyword arguments.

## Examples
```jldoctest
# Get default parameters
params_Default = VegetationParameters(vegtype=Crops)
# Adapt a default value
params = VegetationParameters(vegtype=Crops, rsmin=300.0) 
params.rsmin == 300.0
# output
true
```
"""
@with_kw struct VegetationParameters{FT}
    vegtype::Type{<:VegetationType}
    rsmin::FT = get_default_value(vegtype)[1]
    gd::FT = get_default_value(vegtype)[2]
end
