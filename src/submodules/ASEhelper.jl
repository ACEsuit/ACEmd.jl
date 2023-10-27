module ASEhelper

using ..ACEmd
using AtomsBase
using StaticArrays
using Unitful


export ase_energy
export ase_forces
export ase_virial
export load_ace_model


function make_system(
    atom_symbols,
    positions,
    cell,
    pbc
)
    atoms = [ Atom(Symbol(s), SVector(r...) *u"Å" ) for (s, r) in zip(atom_symbols, eachrow(positions))  ]
    box = [ [cell[1], 0., 0.], [0., cell[2], 0.], [0., 0., cell[3]]  ]u"Å"
    periodic = [ x ? Periodic() : DirichletZero() for x in pbc ]
    sys = FlexibleSystem(atoms, box, periodic)
    return sys
end

function ase_energy(
    potential,
    atom_symbols,
    positions,
    cell,
    pbc,  
)
    sys = make_system(atom_symbols, positions, cell, pbc)
    e = ace_energy(potential, sys)
    return ustrip(u"eV", e)
end

function ase_forces(
    potential,
    atom_symbols,
    positions,
    cell,
    pbc,  
)
    sys = make_system(atom_symbols, positions, cell, pbc)
    f = ace_forces(potential, sys)
    unit(f[1][1])
    fm = reinterpret(reshape, typeof(1.0*unit(f[1][1])), f)
    return ustrip.(u"eV/Å", fm)'  # we need transpose for python
end

function ase_virial(
    potential,
    atom_symbols,
    positions,
    cell,
    pbc,  
)
    sys = make_system(atom_symbols, positions, cell, pbc)
    v = ace_virial(potential, sys)
    return ustrip.(u"eV", v)'  # we need transpose for python
end


end