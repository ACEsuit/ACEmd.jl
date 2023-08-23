# ACE in Molly

ACE support for [Molly](https://github.com/JuliaMolSim/Molly.jl) is loaded with Julia package extensions.
To use you need to use at least v1.9 of Julia. You also need to have both Molly and ACEapi added to the environment you are using, for extension to be loaded.

There are couple of notes that need to be understood.

- Default unit for energy in Molly is kJ/mol. You need to change this to anything that is not per mole (same for force units)
- As of writing this Molly does not have fully compatible AtomsBase interface. So building `System` structure is a bit clunky. Most of all, `atoms_data` needs to have defined `element` for atomic symbol and `Z` for nuclear charge. As these are needed by ACE.

## Example

```julia
using Molly
using ACEmd

using ExtXYZ
using Unitful


# Load ACE data and initial structure
fname_ace = joinpath(pkgdir(ACEmd), "data", "TiAl.json")
fname_xyz = joinpath(pkgdir(ACEmd), "data", "TiAl-big.xyz")

data = ExtXYZ.Atoms(read_frame(fname_xyz))
pot = load_ace_model(fname_ace)

# Pack data to Molly compatible format
# data is AtomsBase system type
sys = Molly.System(data, pot)

# Set up temperature and velocities
temp = 298.0u"K"
vel = random_velocities(sys, tmp)

# Add velocities and loggers
sys = Molly.System(
    sys;
    velocities = vel,
    loggers=(temp=TemperatureLogger(100),)
)

# Set up simulator
simulator = VelocityVerlet(
    dt=1.0u"fs",
    coupling=AndersenThermostat(temp, 1.0u"ps"),
)


# Perform MD
simulate!(sys, simulator, 1000)
```
