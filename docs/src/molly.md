# ACE in Molly

ACE support for [Molly](https://github.com/JuliaMolSim/Molly.jl) is loaded with Julia package extensions.
To use you need to use at least v1.9 of Julia. You also need to have both Molly and ACEapi added to the environment you are using, for extension to be loaded.

You need to understand that by default Molly uses kJ/mol for energy. This is not compatible with ACEmd and AtomsCalculators, so you need to change the energy unit. ACEmd adds convenience functions to do this for you and at the same time initialize the system.

## Example

```julia
using Molly
using ACEmd

using ExtXYZ
using Unitful


# Load ACE model and initial structure
fname_ace = joinpath(pkgdir(ACEmd), "data", "TiAl.json")
fname_xyz = joinpath(pkgdir(ACEmd), "data", "TiAl-big.xyz")

data = ExtXYZ.Atoms(read_frame(fname_xyz))
pot = load_ace_model(fname_ace)

# Insert the potential to Molly and change the units
# to reflect this. Also initializes system structure.
sys = Molly.System(data, pot)

# Set up temperature and velocities
temp = 298.0u"K"
vel = random_velocities(sys, temp)

# Add initial velocities and loggers
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
