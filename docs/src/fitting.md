# Fitting Potentials With ACEfit

ACEmd has an extension fot ACEfit to perform fitting. It is only basic features. For more fitting methods look for [ACEpotenstials](https://github.com/ACEsuit/ACEpotentials.jl)

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

# Assemble training dataa
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