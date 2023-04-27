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
        n_threads=nothing,
        executor=ThreadedEx())
    # sys has atoms for atom identy data
    # coords for coordinates
    # and boundary for boundary conditions
    return ace_forces(acp, sys; executor=executor) *u"hartree/Å" * Nₐ
end


function Molly.potential_energy(
        acp::ACEpotential,
        sys,
        neighbors=nothing;
        n_threads=nothing,
        executor=ThreadedEx()
        )
    # sys has atoms for atom identy data
    # coords for coordinates
    # and boundary for boundary conditions

    return ace_energy(acp, sys; executor=executor) * u"hartree" * Nₐ
end
 

function ACEapi._atomic_number(sys::Molly.System, i) 
    return ACE1.AtomicNumber( sys.atoms_data[i].Z )
end

end # Module