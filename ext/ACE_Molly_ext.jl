module ACE_Molly_ext

import ACE1: AtomicNumber
using ACEmd
using AtomsBase
using LinearAlgebra: isdiag
using Molly

 

function Molly.System(
        sys::AbstractSystem,
        pot::ACEpotential;
        energy_units = pot.energy_unit,
        force_units = pot.energy_unit/pot.length_unit,
        velocity_units = pot.length_unit/u"fs",
        kwargs...
    )
    @assert dimension(energy_units) == dimension(u"J") "energy unit has wrong dimenstions"
    @assert dimension(force_units) == dimension(u"kg*m/s^2")  "force unit has wrong dimenstions"
    @assert dimension(velocity_units) == dimension(u"m/s")  "velocity unit has wrong dimenstions"

    atoms = [Molly.Atom(; index=i, mass=atomic_mass(sys, i) ) for i in 1:length(sys) ]
    atoms_data = [ Molly.AtomData(; element=String(atomic_symbol(sys,i))) for i in 1:length(sys)]

    boundary = begin
        box = bounding_box(sys)
        if isdiag( hcat(box...) )
           tmp = CubicBoundary(box[1][1], box[2][2], box[3][3])
        else
           tmp = TriclinicBoundary(box...)
        end
        tmp
    end
    
    return Molly.System(;
        atoms=atoms,
        atoms_data=atoms_data,
        coords= map(sys) do r
            SVector(position(r)...)
        end,
        general_inters = (pot,),
        boundary=boundary,
        energy_units=energy_units,
        force_units=force_units,
        velocities = map(sys) do a
            tmp = SVector(velocity(a)...)
            uconvert.(velocity_units, tmp)
        end,
        kwargs...
    )
end

end # Module
