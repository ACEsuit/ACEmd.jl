module ACE_Molly_ext

using ACE1
using ACEmd
using CellListMap
using Folds
using Molly
using Unitful
using UnitfulAtomic


function Molly.forces(
        acp::ACEpotential,
        sys,
        neighbors=nothing;
        n_threads=nothing,
        executor=ThreadedEx()
    )
    return ace_forces(acp, sys; executor=executor) *u"hartree/Å"
end


function Molly.potential_energy(
        acp::ACEpotential,
        sys,
        neighbors=nothing;
        n_threads=nothing,
        executor=ThreadedEx()
    )
    return ace_energy(acp, sys; executor=executor) * u"hartree"
end
 

function ACEmd._atomic_number(sys::Molly.System, i) 
    return ACE1.AtomicNumber( sys.atoms_data[i].Z )
end

end # Module
