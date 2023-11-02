```@meta
CurrentModule = ACEmd
```

# ACEmd

[ACEmd](https://github.com/acesuit/ACEmd.jl) is molecular dynamics interface for Julia ACE. All MD use cases should use this package.

ACEmd is fully [AtomsBase](https://github.com/JuliaMolSim/AtomsBase.jl) compatable and you should use AtomsBase over older JuLIP. You can still use JuLIP `Atoms` type as an input, while in a long run JuLIP will be decomissioned.

## Example use case

```@example 1
using ACEmd # Reexports AtomsBase and Unitful
using ExtXYZ

# Load potential and atomic structure
fname_ace = joinpath(pkgdir(ACEmd), "data", "TiAl.json")
fname_xyz = joinpath(pkgdir(ACEmd), "data", "TiAl-big.xyz")

data = ExtXYZ.load(fname_xyz)
pot = load_ace_model(fname_ace)

# use ACE to calculate...
ace_energy(pot, data)
ace_forces(pot, data)
ace_virial(pot, data)
```

## Units

ACEmd supports [Unitful](https://github.com/PainterQubits/Unitful.jl) and it is reexported by default.

Defaults units are

```@repl 1
ACEmd.default_energy
ACEmd.default_length
```

Units can be changes when loading potential

```@example 1
pot_ev = load_ace_model(fname_ace;
    energy_unit=u"eV",
    length_unit=u"Å",
    cutoff_unit=u"Å"
)

ace_energy(pot_ev, data)
```

Alternatively you can overload the default unit by giving it explicitely

```@repl 1
ace_energy(pot_ev, data; energy_unit=u"hartree")
```

Following units are taken into account

- `energy_unit` defines the energy unit of the potential
- `length_unit` defines the length unit for the potential
- `cutoff_unit` defines unit for the cutoff that the potential is using


## Energy of each atom

You can get ACE energy of each atom by

```@repl 1
ace_atom_energies(pot, data)
```


## Index

```@index
```

```@autodocs
Modules = [ACEmd, ACEmd.IPIprotocol]
```
