# Benchmarking

To use the benchmarks you need to install additional packages

```julia
using Pkg
Pkg.add("ExtXYZ")
Pkg.add("AtomsBase")
Pkg.add("ACE1")
Pkg.add("BenchmarkTools")
Pkg.add("PkgBenchmark")
```

After that you calculate the benchmarks by calling

```julia
using ACEmd
using BenchmarkTools
using PkgBenchmark

bench = benchmarkpkg(ACEmd)
results = bench.benchmarkgroup
```

You can inspect the results using [@tagged](https://juliaci.github.io/BenchmarkTools.jl/stable/manual/#Indexing-into-a-BenchmarkGroup-using-@tagged) macro from [BenchmarkTools](https://github.com/JuliaCI/BenchmarkTools.jl)

```julia
julia> results[@tagged "energy" && "AtomsBase"]
1-element BenchmarkTools.BenchmarkGroup:
  tags: []
  "AtomsBase" => 1-element BenchmarkTools.BenchmarkGroup:
          tags: []
          "energy" => 2-element BenchmarkTools.BenchmarkGroup:
                  tags: []
                  "big" => Trial(26.091 ms)
                  "huge" => Trial(544.891 ms)
```

### System size

In the benchmarks there are two system sizes "big" and "huge", which correspond to 1_024 and 16_000 atoms.


## Controlling number of threads and other parameters

The default config uses current Julia parameters to run the config.
You can create your own configs with [BenchmarkConfig](https://juliaci.github.io/PkgBenchmark.jl/stable/run_benchmarks/#PkgBenchmark.BenchmarkConfig)

```julia
t2 = BenchmarkConfig(env = Dict("JULIA_NUM_THREADS" => 2))
benchmarkpkg(ACEmd, t2)
```



## Comparing results

You can use `judge` function to compare benchmarks with different setups or to other branches.

```julia
t4 = BenchmarkConfig(env = Dict("JULIA_NUM_THREADS" => 4))
t8 = BenchmarkConfig(env = Dict("JULIA_NUM_THREADS" => 8))

# Compare how much changing from 4-threads to 8 improves the performance
j = judge(ACEmd, t8, t4)

show(j.benchmarkgroup)
```


Compare current branch to "origin/main"

Note, you need to commit all the changes you made!


```julia
j = judge(ACEmd, "origin/main")
show(j.benchmarkgroup)
```

Note that `BenchmarkConfig` takes parameter `id`, which controls wich branch is used for the benchmark.


## CI Benchmarks on Pull Requests

When making PRs you have an option to add "Run Benchmarks" label to make CI run the benchmark. This should be done every time, when change is made on how calculations are run in `ACEmd`.


## Benchmark Scripts

In `scripts/` folder has a script file that can be used to benchmark scaling. To use it just

```julia
using ACEmd
include( joinpath(pkgdir(ACEmd), "scripts", "benchmark.jl") )

# For 1, 2, 4 and 8 threads
results = bench_scaling([1,2,4,8]) 

# Plot results using UnicodePlots
p = plot_scaling(results)
```

You can also use `judge` for the results

```julia
# Compare 2 and 4 threads 
j = judge(results[3], results[2])
show(j.benchmarkgroup)
```

It is also possible to do the benchmark on a specific branch. This is done by givin `bench_scaling` an extra parameter for branch

```julia
# same as above but for "my-branch" branch
results = bench_scaling([1,2,4,8]; branch="my-branch") 
```