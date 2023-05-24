```@meta
CurrentModule = ACEMD
```

# ACEMD

Documentation for [ACEMD](https://github.com/tjjarvinen/ACEMD.jl).

The main structure is now functional, but testing is still needed.

### Example use case

```julia
using ACEMD
using AtomsBase
using ExtXYZ

# Load potential and atomic structure
fname_ace = joinpath(pkgdir(ACEMD), "data", "TiAl.json")
fname_xyz = joinpath(pkgdir(ACEMD), "data", "TiAl-big.xyz")

data = FastSystem(ExtXYZ.Atoms(read_frame(fname_xyz)))
pot = load_ace_model(fname_ace)

# use ACE to calculate...
ace_energy(pot, data)
ace_forces(pot, data)
ace_virial(pot, data)
```


```@index
```

```@autodocs
Modules = [ACEMD]
```
