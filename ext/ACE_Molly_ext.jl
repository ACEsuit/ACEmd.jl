module ACE_Molly_ext

using ACE1
using ACEapi
using CellListMap
using Folds
using Molly
using Unitful
using UnitfulAtomic

const Nₐ = 6.022140857e23u"mol^-1"

function Molly.forces(
        acp::ACEpotential,
        sys,
        neighbors=nothing;
        n_threads=Threads.nthreads())
    # sys has atoms for atom identy data
    # coords for coordinates
    # and boundary for boundary conditions
    return ace_forces(acp, sys) *u"hartree/Å" * Nₐ
end


function Molly.potential_energy(
        acp::ACEpotential,
        sys,
        neighbors=nothing;
        n_threads=Threads.nthreads())
    # sys has atoms for atom identy data
    # coords for coordinates
    # and boundary for boundary conditions

    return ace_energy(acp, sys) * u"hartree" * Nₐ
end
    

function CellListMap.neighborlist(sys::Molly.System, cutoff; kwargs...)
    cell = sys.boundary.side_lengths
    list = CellListMap.neighborlist(sys.coords, cutoff; unitcell=cell)
    return list
end

function ACEapi._atomic_number(sys::Molly.System, i) 
    return ACE1.AtomicNumber( sys.atoms_data[i] )
end

end # Module