module ACE_Molly_ext

using ACE1
using ACEapi
using CellListMap
using Folds
using Molly



function Molly.forces(
        acp::ACEpotential,
        sys,
        neighbors=nothing;
        n_threads=Threads.nthreads())
    # sys has atoms for atom identy data
    # coords for coordinates
    # and boundary for boundary conditions
    return ace_forces(acp, sys)
end


function Molly.potential_energy(
        acp::ACEpotential,
        sys,
        neighbors=nothing;
        n_threads=Threads.nthreads())
    # sys has atoms for atom identy data
    # coords for coordinates
    # and boundary for boundary conditions
    return ace_energy(acp, sys)
end
    

function CellListMap.neighborlist(sys::MollySystem, cutoff; kwargs...)
    cell = sys.boundary.side_lengths
    list = CellListMap.neighborlist(sys.coords, cutoff; unitcell=cell)
    return list
end

function _atomic_number(sys::System, i) 
    return ACE1.AtomicNumber( sys.atoms_data[i] )
end

end # Module