module ASEhelper

using ..ACEmd
using AtomsBase
using Distributed
using StaticArrays
using Unitful


export ase_energy
export ase_forces
export ase_virial
export ase_energy_forces_virial
export load_ace_model


function make_system(
    atom_symbols,
    positions,
    cell,
    pbc
)
    atoms = [ Atom(Symbol(s), SVector(r...) *u"Å" ) for (s, r) in zip(atom_symbols, eachrow(positions))  ]
    box = collect( eachrow(Matrix(cell)*u"Å") )
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
    e = fetch( @spawn ace_energy(potential, sys) )
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
    f = fetch( @spawn ace_forces(potential, sys) )
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
    v = fetch( @spawn ace_virial(potential, sys) )
    return ustrip.(u"eV", v)'  # we need transpose for python
end

function ase_energy_forces_virial(
    potential,
    atom_symbols,
    positions,
    cell,
    pbc,
)
    sys = make_system(atom_symbols, positions, cell, pbc)
    tmp = fetch( @spawn ace_energy_forces_virial(potential, sys) )
    f = tmp["forces"]
    fm = reinterpret(reshape, typeof(1.0*unit(f[1][1])), f)
    e = tmp["energy"]
    v = tmp["virial"]
    return ustrip(u"eV", e), ustrip.(u"eV/Å", fm), ustrip.(u"eV", v)
end


end
