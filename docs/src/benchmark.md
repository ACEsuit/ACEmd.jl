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

bench = benchmark(ACEmd)
results = bench.benchmarkgroup
```

You can inspect the results using [@tagged](https://juliaci.github.io/BenchmarkTools.jl/stable/manual/#Indexing-into-a-BenchmarkGroup-using-@tagged) macro from [BenchmarkTools](https://github.com/JuliaCI/BenchmarkTools.jl)

```julia
julia> results[@tagged "energy"][@tagged "AtomsBase"]
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



## Comparing results to other branches

You can use [BenchmarkCI](https://github.com/JuliaCI/PkgBenchmark.jl) to compare benchmarks to other branches. This is useful only for `ACEmd` development.

Note, you need to commit all the changes you made!

```julia
BenchmarkCI
judge(ACEmd, "origin/main")
displayjudgement()
```

When making PRs you have an option to add "Run Benchmarks" label to make CI run the benchmark. This should be done every time, when change is made on how calculations are run in `ACEmd`.