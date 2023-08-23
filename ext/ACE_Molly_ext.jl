module ACE_Molly_ext

import ACE1: AtomicNumber
using ACEmd
using AtomsBase
using LinearAlgebra: isdiag
using Molly


function Molly.forces(
        acp::ACEpotential,
        sys,
        neighbors=nothing;
        n_threads=nothing,
        executor=ThreadedEx()
    )
    return ace_forces(acp, sys; executor=executor) 
end


function Molly.potential_energy(
        acp::ACEpotential,
        sys,
        neighbors=nothing;
        n_threads=nothing,
        executor=ThreadedEx()
    )
    return ace_energy(acp, sys; executor=executor) 
end
 

function ACEmd._atomic_number(sys::Molly.System, i) 
    return AtomicNumber( sys.atoms_data[i].Z )
end

function Molly.System(
    sys::AbstractSystem,
    pot::ACEpotential;
    kwargs...
    )
    atoms = [Molly.Atom( index=i, mass=atomic_mass(sys, i) ) for i in 1:length(sys) ]

    boundary = begin
        box = bounding_box(sys)
        if isdiag( hcat(box...) )
           tmp = CubicBoundary(box[1][1], box[2][2], box[3][3])
        else
           tmp = TriclinicBoundary(box...)
        end
        tmp
    end

    atom_data = [ (; :Z=>z,:element=>s)  for (z,s) in zip(atomic_number(sys), atomic_symbol(sys))  ]
    
    return Molly.System(
        atoms=atoms,
        atoms_data = atom_data,
        coords= map(sys) do r
            SVector(position(r)...)
        end,
        general_inters = (pot,),
        boundary=boundary,
        energy_units=pot.energy_unit,
        force_units=pot.energy_unit/pot.length_unit,
        velocities = map(sys) do a
            SVector(velocity(a)...)
        end,
        kwargs...
    )
end

end # Module
