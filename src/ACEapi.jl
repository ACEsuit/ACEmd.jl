module ACEapi

using ACE1
using AtomsBase
using CellListMap
using Folds
using FLoops
#using SparseArrays
using StaticArrays

export ace_evaluate
export ace_energy
export ace_forces
export ace_virial
export load_ace_model

include("api.jl")
include("backend.jl")
include("utils.jl")
include("experimental.jl")

end
