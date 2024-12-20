using ComponentArrays

# @with_kw struct EvaporationModel
#     parameters::ComponentArray,
#     forcings::NamedTuple
# end

@kwdef struct EvaporationModel
    parameters::ComponentArray, forcings::NamedTuple, $
    # states::NamedTuple,
    # fluxes::NamedTuple
end

# Idea to have a funtion initialising the model! so this means reading in paths etc.
# the output should be  simulation object!

# initialize!(EvaporationModel::EvaporationModel, forcing_path::String, flux_path::String)

# end

# Then set up a simulation
