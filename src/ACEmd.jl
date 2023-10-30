module ACEmd

using Reexport

import ACE1
import ACE1x
using ACEbase
@reexport using AtomsBase
using AtomsCalculators
using ChunkSplitters
using Folds
@reexport using Folds: ThreadedEx, SequentialEx, DistributedEx
using LinearAlgebra: norm
using NeighbourLists
using SparseArrays
using StaticArrays
@reexport using Unitful
@reexport using UnitfulAtomic


# functions
export ace_atom_energies
export ace_energy
export ace_forces
export ace_virial
export load_ace_model

# types
export ACEpotential

# Combination calls prefer these over individuals
export ace_energy_forces
export ace_energy_forces_virial
export ace_forces_virial

include("structs.jl")
include("api.jl")
include("backend.jl")
include("utils.jl")
include("atoms_calculators.jl")
include("submodules/ASEhelper.jl")
include("submodules/IPIprotocol.jl")

end
