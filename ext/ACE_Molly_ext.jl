module ACE_Molly_ext

using ACE1
using ACEapi
using CellListMap
using Folds
using Molly
using Unitful
using UnitfulAtomic


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

    # ACE has hartree units
    return ace_energy(acp, sys) * 2625.5u"kJ/mol"
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