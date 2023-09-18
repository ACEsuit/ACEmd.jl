# Fitting Potentials With ACEfit

ACEmd has an extension fot ACEfit to perform fitting. It is only basic features. For more fitting methods look for [ACEpotentials](https://github.com/ACEsuit/ACEpotentials.jl).

Here is a basic fitting example
 
```@example fit
using ACE1x
using ACEfit
using ACEmd
using ExtXYZ

# Load training data
fname_train = joinpath(pkgdir(ACEmd), "data", "TiAl-train.xyz")
data = ExtXYZ.Atoms.( ExtXYZ.read_frames(fname_train) )

# Generate ACE basis to be trained
basis = ACE1x.ace_basis(
    elements = [:Ti, :Al],
    order = 3,
    totaldegree = 6,
    rcut = 5.5
);

# Generate one body potential
Vref = OneBody(:Ti => -1586.0195, :Al => -105.5954);

# Assemble training data
A, Y, W = ACEfit.assemble(data, basis; energy_default_weight=5, energy_ref=Vref)

# Smoothness prior
P = smoothness_prior(basis; p = 4)

# Solver
solver = ACEfit.LSQR(damp = 1e-2, atol = 1e-6, P = P)

# Fit the potential
results = ACEfit.solve(solver, W .* A, W .* Y)

# Form the final potential
ACEpotential(basis, results["C"], Vref)
```

## Customize Assembly

### Weights

There are two types of weights default weights for energy, force and virial that apply for all structures, and an individual weight for a specific structure.

General weights are given as a keywords for `assemble` command

```julia
ACEfit.assemble(data, basis;
    energy_default_weight=5,
    force_default_weight=2,
    virial_default_weight=0.5
)
```

Individual weight is given for a specific structure in training data.

```julia
data[1] = FlexibleSystem(data[1]; weight=4)
```

Also energy, force and virial can be given individual weight,

```julia
data[1] = FlexibleSystem(data[1];
    energy_weight=4,
    force_weight=2,
    virial_weight=1.5
)
```

### Control What is Used in Fitting

Training data needs to have comparison data in order for it to be used.
The keys for training data are controlled by keyword arguments:

```julia
ACEfit.assemble(data, basis;
    energy_key=:energy,
    force_key=:force,
    virial_key=:virial
)
```

Giving different values allows you to control what names are used for training.

You can test, if individual training point has the property by

```julia
julia> haskey(data[1], :virial)
true
```

You can disable energy, force or virial from being used in training, even though it is present in data, by giving `assemble` corresponding keyword

```julia
ACEfit.assemble(data, basis;
    energy=true,
    force=true,
    virial=false  # disable virial from training
)
```

### Passing Commands to Calculators

You can pass commands to `ace_energy`, `ace_force` and `ace_virial` calculators as keywords for `assebly`. E.g.

```julia
ACEfit.assemble(data, basis;
    executor=SequentialEx()   # Execute sequentially instead of multithreading
)
```

### One Body Potential

One Body Potential is added as a keyword `energy_ref` for `assemble` command. Not giving it means it is not used in training.

```julia
ACEfit.assemble(data, basis; 
    energy_ref=OneBody(:Ti => -1586.0195, :Al => -105.5954)
)
```

### Multithreading and Parallel Processing

By default Multithreading is in use and assembly uses all threads. Multiprocessing is also supported and used, if more than one process is present.