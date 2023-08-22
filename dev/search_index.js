var documenterSearchIndex = {"docs":
[{"location":"benchmark/#Benchmarking","page":"Benchmarks","title":"Benchmarking","text":"","category":"section"},{"location":"benchmark/","page":"Benchmarks","title":"Benchmarks","text":"To use the benchmarks you need to install additional packages","category":"page"},{"location":"benchmark/","page":"Benchmarks","title":"Benchmarks","text":"using Pkg\nPkg.add(\"ExtXYZ\")\nPkg.add(\"AtomsBase\")\nPkg.add(\"ACE1\")\nPkg.add(\"BenchmarkTools\")\nPkg.add(\"PkgBenchmark\")","category":"page"},{"location":"benchmark/","page":"Benchmarks","title":"Benchmarks","text":"After that you calculate the benchmarks by calling","category":"page"},{"location":"benchmark/","page":"Benchmarks","title":"Benchmarks","text":"using ACEmd\nusing BenchmarkTools\nusing PkgBenchmark\n\nbench = benchmarkpkg(ACEmd)\nresults = bench.benchmarkgroup","category":"page"},{"location":"benchmark/","page":"Benchmarks","title":"Benchmarks","text":"You can inspect the results using @tagged macro from BenchmarkTools","category":"page"},{"location":"benchmark/","page":"Benchmarks","title":"Benchmarks","text":"julia> results[@tagged \"energy\" && \"AtomsBase\"]\n1-element BenchmarkTools.BenchmarkGroup:\n  tags: []\n  \"AtomsBase\" => 1-element BenchmarkTools.BenchmarkGroup:\n          tags: []\n          \"energy\" => 2-element BenchmarkTools.BenchmarkGroup:\n                  tags: []\n                  \"big\" => Trial(26.091 ms)\n                  \"huge\" => Trial(544.891 ms)","category":"page"},{"location":"benchmark/#System-size","page":"Benchmarks","title":"System size","text":"","category":"section"},{"location":"benchmark/","page":"Benchmarks","title":"Benchmarks","text":"In the benchmarks there are two system sizes \"big\" and \"huge\", which correspond to 1024 and 16000 atoms.","category":"page"},{"location":"benchmark/#Controlling-number-of-threads-and-other-parameters","page":"Benchmarks","title":"Controlling number of threads and other parameters","text":"","category":"section"},{"location":"benchmark/","page":"Benchmarks","title":"Benchmarks","text":"The default config uses current Julia parameters to run the config. You can create your own configs with BenchmarkConfig","category":"page"},{"location":"benchmark/","page":"Benchmarks","title":"Benchmarks","text":"t2 = BenchmarkConfig(env = Dict(\"JULIA_NUM_THREADS\" => 2))\nbenchmarkpkg(ACEmd, t2)","category":"page"},{"location":"benchmark/#Comparing-results","page":"Benchmarks","title":"Comparing results","text":"","category":"section"},{"location":"benchmark/","page":"Benchmarks","title":"Benchmarks","text":"You can use judge function to compare benchmarks with different setups or to other branches.","category":"page"},{"location":"benchmark/","page":"Benchmarks","title":"Benchmarks","text":"t4 = BenchmarkConfig(env = Dict(\"JULIA_NUM_THREADS\" => 4))\nt8 = BenchmarkConfig(env = Dict(\"JULIA_NUM_THREADS\" => 8))\n\n# Compare how much changing from 4-threads to 8 improves the performance\nj = judge(ACEmd, t8, t4)\n\nshow(j.benchmarkgroup)","category":"page"},{"location":"benchmark/","page":"Benchmarks","title":"Benchmarks","text":"Compare current branch to \"origin/main\"","category":"page"},{"location":"benchmark/","page":"Benchmarks","title":"Benchmarks","text":"Note, you need to commit all the changes you made!","category":"page"},{"location":"benchmark/","page":"Benchmarks","title":"Benchmarks","text":"j = judge(ACEmd, \"origin/main\")\nshow(j.benchmarkgroup)","category":"page"},{"location":"benchmark/","page":"Benchmarks","title":"Benchmarks","text":"Note that BenchmarkConfig takes parameter id, which controls wich branch is used for the benchmark.","category":"page"},{"location":"benchmark/#CI-Benchmarks-on-Pull-Requests","page":"Benchmarks","title":"CI Benchmarks on Pull Requests","text":"","category":"section"},{"location":"benchmark/","page":"Benchmarks","title":"Benchmarks","text":"When making PRs you have an option to add \"Run Benchmarks\" label to make CI run the benchmark. This should be done every time, when change is made on how calculations are run in ACEmd.","category":"page"},{"location":"benchmark/#Benchmark-Scripts","page":"Benchmarks","title":"Benchmark Scripts","text":"","category":"section"},{"location":"benchmark/","page":"Benchmarks","title":"Benchmarks","text":"In scripts/ folder has a script file that can be used to benchmark scaling. To use it just","category":"page"},{"location":"benchmark/","page":"Benchmarks","title":"Benchmarks","text":"using ACEmd\ninclude( joinpath(pkgdir(ACEmd), \"scripts\", \"benchmark.jl\") )\n\n# For 1, 2, 4 and 8 threads\nresults = bench_scaling([1,2,4,8]) \n\n# Plot results using UnicodePlots\np = plot_scaling(results)","category":"page"},{"location":"benchmark/","page":"Benchmarks","title":"Benchmarks","text":"You can also use judge for the results","category":"page"},{"location":"benchmark/","page":"Benchmarks","title":"Benchmarks","text":"# Compare 2 and 4 threads \nj = judge(results[3], results[2])\nshow(j.benchmarkgroup)","category":"page"},{"location":"benchmark/","page":"Benchmarks","title":"Benchmarks","text":"It is also possible to do the benchmark on a specific branch. This is done by givin bench_scaling an extra parameter for branch","category":"page"},{"location":"benchmark/","page":"Benchmarks","title":"Benchmarks","text":"# same as above but for \"my-branch\" branch\nresults = bench_scaling([1,2,4,8]; branch=\"my-branch\") ","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = ACEmd","category":"page"},{"location":"#ACEmd","page":"Home","title":"ACEmd","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"ACEmd is molecular dynamics interface for Julia ACE. All MD use cases should use this package.","category":"page"},{"location":"","page":"Home","title":"Home","text":"ACEmd is fully AtomsBase compatable and you should use AtomsBase over older JuLIP. You can still use JuLIP Atoms type as an input, while in a long run JuLIP will be decomissioned.","category":"page"},{"location":"#Example-use-case","page":"Home","title":"Example use case","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"using ACEmd # Reexports AtomsBase and Unitful\nusing ExtXYZ\n\n# Load potential and atomic structure\nfname_ace = joinpath(pkgdir(ACEmd), \"data\", \"TiAl.json\")\nfname_xyz = joinpath(pkgdir(ACEmd), \"data\", \"TiAl-big.xyz\")\n\ndata = FastSystem(ExtXYZ.Atoms(read_frame(fname_xyz)))\npot = load_ace_model(fname_ace)\n\n# use ACE to calculate...\nace_energy(pot, data)\nace_forces(pot, data)\nace_virial(pot, data)","category":"page"},{"location":"#Units","page":"Home","title":"Units","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"ACEmd supports Unitful and it is reexported by default.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Defaults units are","category":"page"},{"location":"","page":"Home","title":"Home","text":"ACEmd.default_energy\nACEmd.default_length","category":"page"},{"location":"","page":"Home","title":"Home","text":"Units can be changes when loading potential","category":"page"},{"location":"","page":"Home","title":"Home","text":"pot_ev = load_ace_model(fname_ace;\n    energy_unit=u\"eV\",\n    length_unit=u\"Å\",\n    cutoff_unit=u\"Å\"\n)\n\nace_energy(pot_ev, data)","category":"page"},{"location":"","page":"Home","title":"Home","text":"Alternatively you can overload the default unit by giving it explicitely","category":"page"},{"location":"","page":"Home","title":"Home","text":"ace_energy(pot_ev, data; energy_unit=u\"hartree\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"Following units are taken into account","category":"page"},{"location":"","page":"Home","title":"Home","text":"energy_unit defines the energy unit of the potential\nlength_unit defines the length unit for the potential\ncutoff_unit defines unit for the cutoff that the potential is using","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [ACEmd]","category":"page"},{"location":"#ACEmd.ace_energy-Tuple{Any, Any}","page":"Home","title":"ACEmd.ace_energy","text":"ace_energy(potential, ACE1.Atoms, Kwargs)\n\nCalculates ACE potential energy for atomic system. The ACE1.Atoms object needs to be in ACE1.AtomsBase compatable format. The returned energy has a unit as defined by Unitful.\n\nParallel execution is done with Transducers.jl and there is an option to use different executors. Look for ThreadedEx for more details on how to control it. \n\nKwargs\n\ndomain=1:length(ACE1.Atoms)  :  choose subset of ACE1.Atoms to which energy is calculated\nexecutor=ThreadedEx()   :  used to control multithreading using Transducers.jl\nenergy_unit  :   used to override energy unit for the calculation\nlength_unit  :   used to override lenght unit for the calculation\ncutoff_unit  :   used to override unit that cutoff radius is defined\n\n\n\n\n\n","category":"method"},{"location":"#ACEmd.ace_forces-Tuple{Any, Any}","page":"Home","title":"ACEmd.ace_forces","text":"ace_forces(potential, ACE1.Atoms, Kwargs)\n\nCalculates forces for ACE potential for given atomic system. The ACE1.Atoms object needs to be in ACE1.AtomsBase compatable format. The returned energy has a unit as defined by Unitful.\n\nParallel execution is done with Transducers.jl and there is an option to use different executors. Look for ThreadedEx for more details on how to control it. \n\nKwargs\n\ndomain=1:length(ACE1.Atoms)  :  choose subset of ACE1.Atoms to which energy is calculated\nexecutor=ThreadedEx()   :  used to control multithreading using Transducers.jl\nenergy_unit  :   used to override energy unit for the calculation\nlength_unit  :   used to override lenght unit for the calculation\ncutoff_unit  :   used to override unit that cutoff radius is defined\n\n\n\n\n\n","category":"method"},{"location":"#ACEmd.ace_virial-Tuple{Any, Any}","page":"Home","title":"ACEmd.ace_virial","text":"ace_virial(potential, ACE1.Atoms, Kwargs)\n\nCalculates virial for ACE potential for given atomic system. The ACE1.Atoms object needs to be in ACE1.AtomsBase compatable format. The returned energy has a unit as defined by Unitful.\n\nParallel execution is done with Transducers.jl and there is an option to use different executors. Look for ThreadedEx for more details on how to control it. \n\nKwargs\n\ndomain=1:length(ACE1.Atoms)  :  choose subset of ACE1.Atoms to which energy is calculated\nexecutor=ThreadedEx()   :  used to control multithreading using Transducers.jl\nenergy_unit  :   used to override energy unit for the calculation\nlength_unit  :   used to override lenght unit for the calculation\ncutoff_unit  :   used to override unit that cutoff radius is defined\n\n\n\n\n\n","category":"method"},{"location":"#ACEmd.load_ace_model-Tuple{AbstractString}","page":"Home","title":"ACEmd.load_ace_model","text":"load_ace_model(fname::AbstractString; Kwargs)\nload_ace_model(potential: Kwargs)\n\nUsed to load potential from json or yml files. By default this adds units for the potential. You can use Kwargs to control what units are used. The default units are defined in ACEmd.default_energy and ACEmd.default_length.\n\nYou can also use this to convert non unit holding potentials to potentials that have units. Like e.g. when you have fitted a new potential it does not have unit. You can then use this function to wrap units for it.\n\nKwargs\n\nenergy_unit=default_energy  :  energy unit for the potential\nlength_unit=default_length  :  lenght unit used for the force and atomic system\ncutoff_unit=default_length  :  unit in which cuttoff radius has been defined\n\n\n\n\n\n","category":"method"},{"location":"molly/#ACE-in-Molly","page":"Using ACE in Molly","title":"ACE in Molly","text":"","category":"section"},{"location":"molly/","page":"Using ACE in Molly","title":"Using ACE in Molly","text":"ACE support for Molly is loaded with Julia package extensions. To use you need to use at least v1.9 of Julia. You also need to have both Molly and ACEapi added to the environment you are using, for extension to be loaded.","category":"page"},{"location":"molly/","page":"Using ACE in Molly","title":"Using ACE in Molly","text":"There are couple of notes that need to be understood.","category":"page"},{"location":"molly/","page":"Using ACE in Molly","title":"Using ACE in Molly","text":"Default unit for energy in Molly is kJ/mol. You need to change this to anything that is not per mole (same for force units)\nAs of writing this Molly does not have fully compatible AtomsBase interface. So building System structure is a bit clunky. Most of all, atoms_data needs to have defined element for atomic symbol and Z for nuclear charge. As these are needed by ACE.","category":"page"},{"location":"molly/#Example","page":"Using ACE in Molly","title":"Example","text":"","category":"section"},{"location":"molly/","page":"Using ACE in Molly","title":"Using ACE in Molly","text":"using Molly\nusing ACEmd\n\nusing AtomsBase\nusing ExtXYZ\nusing Unitful\n\n\n# Load ACE data and initial structure\nfname_ace = joinpath(pkgdir(ACEmd), \"data\", \"TiAl.json\")\nfname_xyz = joinpath(pkgdir(ACEmd), \"data\", \"TiAl-big.xyz\")\n\ndata = FastSystem(ExtXYZ.Atoms(read_frame(fname_xyz)))\npot = load_ace_model(fname_ace)\n\n# Prepare data to Molly compatible format\natoms = [Molly.Atom( index=i, mass=atomic_mass(data, i) ) for i in 1:length(data) ]\n\nboundary = begin\n    box = bounding_box(data)\n    CubicBoundary(box[1][1], box[2][2], box[3][3])\nend\n\natom_data = [ (; :Z=>z,:element=>s)  for (z,s) in zip(atomic_number(data), atomic_symbol(data))  ]\n\n# Set up temperature and velocities\ntemp = 298.0u\"K\"\nvelocities = [random_velocity(m, temp) for m in atomic_mass(data)]\n\n# Set up simulator\nsimulator = VelocityVerlet(\n    dt=1.0u\"fs\",\n    coupling=AndersenThermostat(temp, 1.0u\"ps\"),\n)\n\n# Set up Molly system\nsys = System(\n           atoms=atoms,\n           atoms_data = atom_data,\n           coords=position(data),\n           velocities=velocities,\n           general_inters = (pot,),\n           boundary=boundary,\n           loggers=(temp=TemperatureLogger(100),),\n           energy_units=u\"eV\",\n           force_units=u\"eV/Å\",\n       )\n\n\n# Perform MD\nsimulate!(sys, simulator, 1000)","category":"page"}]
}
