# ACE in Molly

ACE support for [Molly](https://github.com/JuliaMolSim/Molly.jl) is loaded with Julia package extensions.
To use you need to use at least v1.9 of Julia. You also need to have both Molly and ACEapi added to the environment you are using, for extension to be loaded.

There are couple of notes that need to be understood.

- Default unit for energy in Molly is kJ/mol. You need to change this to anything that is not per mole (same for force units)
- As of writing this Molly does not have fully compatible AtomsBase interface. So building `System` structure is a bit clunky. Most of all, `atoms_data` needs to have defined `element` for atomic symbol and `Z` for nuclear charge. As these are needed by ACE.

## Example

```julia
using Molly
using ACEapi

using AtomsBase
using ExtXYZ
using Unitful


# Load ACE data and initial structure
fname_ace = joinpath(pkgdir(ACEapi), "data", "TiAl.json")
fname_xyz = joinpath(pkgdir(ACEapi), "data", "TiAl-big.xyz")

data = FastSystem(ExtXYZ.Atoms(read_frame(fname_xyz)))
pot = load_ace_model(fname_ace)

# Prepare data to Molly compatible format
atoms = [Molly.Atom( index=i, mass=atomic_mass(data, i) ) for i in 1:length(data) ]

boundary = begin
    box = bounding_box(data)
    CubicBoundary(box[1][1], box[2][2], box[3][3])
end

atom_data = [ (; :Z=>z,:element=>s)  for (z,s) in zip(atomic_number(data), atomic_symbol(data))  ]

# Set up temperature and velocities
temp = 298.0u"K"
velocities = [random_velocity(m, temp) for m in atomic_mass(data)]

# Set up simulator
simulator = VelocityVerlet(
    dt=1.0u"fs",
    coupling=AndersenThermostat(temp, 1.0u"ps"),
)

# Set up Molly system
sys = System(
           atoms=atoms,
           atoms_data = atom_data,
           coords=position(data),
           velocities=velocities,
           general_inters = (pot,),
           boundary=boundary,
           loggers=(temp=TemperatureLogger(100),),
           energy_units=u"eV",
           force_units=u"eV/Å",
       )


# Perform MD
simulate!(sys, simulator, 1000)
```
