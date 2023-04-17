module ACEapi

using ACE1
using AtomsBase
using CellListMap
using Folds
using FLoops
#using SparseArrays
using StaticArrays
using Unitful

export ace_energy
export ace_forces
export ace_virial
export load_ace_model

# Combination calls prefer these over individuals
export ace_energy_forces
export ace_energy_forces_virial
export ace_forces_virial

include("api.jl")
include("backend.jl")
include("utils.jl")
include("experimental.jl")

end
