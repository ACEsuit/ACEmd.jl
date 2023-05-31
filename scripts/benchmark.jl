using ACE1
using ACEmd
using AtomsBase
using ExtXYZ
using BenchmarkTools
using PkgBenchmark
using UnicodePlots
using Unitful

##

function bench_scaling(n_threads; branch="")
    scaling = map(n_threads) do n
        @info "Benchmark $(n) threads"
        conf = BenchmarkConfig(id=branch, env = Dict("JULIA_NUM_THREADS" => n))
        bench = benchmarkpkg(ACEmd, conf)
        bench
    end
    return scaling
end

##

function chop_results(results)
    energy = Dict()
    energy["AB big"] = map( results  ) do res
        tmp = res.benchmarkgroup["AtomsBase"]["energy"]["big"]
        mean(tmp.times) * u"ns" |> u"ms"
    end
    energy["AB huge"] = map( results  ) do res
        tmp = res.benchmarkgroup["AtomsBase"]["energy"]["huge"]
        mean(tmp.times) * u"ns" |> u"ms"
    end
    energy["JL huge"] = map( results  ) do res
        tmp = res.benchmarkgroup["JuLIP"]["energy"]["huge"]
        mean(tmp.times) * u"ns" |> u"ms"
    end
    energy["JL big"] = map( results  ) do res
        tmp = res.benchmarkgroup["JuLIP"]["energy"]["big"]
        mean(tmp.times) * u"ns" |> u"ms"
    end

    forces = Dict()
    forces["AB big"] = map( results  ) do res
        tmp = res.benchmarkgroup["AtomsBase"]["forces"]["big"]
        mean(tmp.times) * u"ns" |> u"ms"
    end
    forces["AB huge"] = map( results  ) do res
        tmp = res.benchmarkgroup["AtomsBase"]["forces"]["huge"]
        mean(tmp.times) * u"ns" |> u"ms"
    end
    forces["JL huge"] = map( results  ) do res
        tmp = res.benchmarkgroup["JuLIP"]["forces"]["huge"]
        mean(tmp.times) * u"ns" |> u"ms"
    end
    forces["JL big"] = map( results  ) do res
        tmp = res.benchmarkgroup["JuLIP"]["forces"]["big"]
        mean(tmp.times) * u"ns" |> u"ms"
    end

    virial = Dict()
    virial["AB big"] = map( results  ) do res
        tmp = res.benchmarkgroup["AtomsBase"]["virial"]["big"]
        mean(tmp.times) * u"ns" |> u"ms"
    end
    virial["AB huge"] = map( results  ) do res
        tmp = res.benchmarkgroup["AtomsBase"]["virial"]["huge"]
        mean(tmp.times) * u"ns" |> u"ms"
    end
    virial["JL huge"] = map( results  ) do res
        tmp = res.benchmarkgroup["JuLIP"]["virial"]["huge"]
        mean(tmp.times) * u"ns" |> u"ms"
    end
    virial["JL big"] = map( results  ) do res
        tmp = res.benchmarkgroup["JuLIP"]["virial"]["big"]
        mean(tmp.times) * u"ns" |> u"ms"
    end

    t = map(results) do res
        res.benchmarkconfig.env["JULIA_NUM_THREADS"]
    end
    return (;
        :nthreads=>t,
        :energy=>energy,
        :forces=>forces,
        :virial=>virial
    )
end


##

function plot_scaling(results)
    e = chop_results( collect(results) )

    x = e.nthreads.^-1


    base_line = e.energy["JL huge"][1]
    mi = min( e.energy["JL huge"]..., e.energy["AB huge"]...)
    ma = max( e.energy["JL huge"]..., e.energy["AB huge"]...)
    ylims = (mi, ma) ./ base_line
    p = lineplot(
        x,
        e.energy["AB huge"]./base_line;
        name="AtomsBase",
        title="AtomsBase energy 16_000 atoms",
        xlabel="1/nthreads",
        ylabel="relative time",
        ylim=ylims,
        xflip=true
    )
    lineplot!(p,
        x,
        e.energy["JL huge"]./base_line;
        name="JuLIP"
    )

    base_line = e.energy["JL big"][1]
    mi = min( e.energy["JL big"]..., e.energy["AB big"]...)
    ma = max( e.energy["JL big"]..., e.energy["AB big"]...)
    ylims = (mi, ma) ./ base_line
    pb = lineplot(
        x,
        e.energy["AB big"]./base_line;
        name="AtomsBase",
        title="AtomsBase energy 1024 atoms",
        xlabel="1/nthreads",
        ylabel="relative time",
        ylim=ylims,
        xflip=true
    )
    lineplot!(pb,
        x,
        e.energy["JL big"]./base_line;
        name="JuLIP"
    )

    base_line = e.forces["JL big"][1]
    mi = min( e.forces["JL big"]..., e.forces["AB big"]...)
    ma = max( e.forces["JL big"]..., e.forces["AB big"]...)
    ylims = (mi, ma) ./ base_line
    pfb = lineplot(
        x,
        e.forces["AB big"]./base_line;
        name="AtomsBase",
        title="AtomsBase forces 1024 atoms",
        xlabel="1/nthreads",
        ylabel="relative time",
        ylim=ylims,
        xflip=true
    )
    lineplot!(pfb,
        x,
        e.forces["JL big"]./base_line;
        name="JuLIP"
    )

    base_line = e.forces["JL huge"][1]
    mi = min( e.forces["JL huge"]..., e.forces["AB huge"]...)
    ma = max( e.forces["JL huge"]..., e.forces["AB huge"]...)
    ylims = (mi, ma) ./ base_line
    pfh = lineplot(
        x,
        e.forces["AB huge"]./base_line;
        name="AtomsBase",
        title="AtomsBase forces 16_000 atoms",
        xlabel="1/nthreads",
        ylabel="relative time",
        ylim=ylims,
        xflip=true
    )
    lineplot!(pfh,
        x,
        e.forces["JL huge"]./base_line;
        name="JuLIP",
    )

    base_line = e.virial["JL big"][1]
    mi = min( e.virial["JL big"]..., e.virial["AB big"]...)
    ma = max( e.virial["JL big"]..., e.virial["AB big"]...)
    ylims = (mi, ma) ./ base_line
    pvb = lineplot(
        x,
        e.virial["AB big"]./base_line;
        name="AtomsBase",
        title="AtomsBase virial 1024 atoms",
        xlabel="1/nthreads",
        ylabel="relative time",
        ylim=ylims,
        xflip=true
    )
    lineplot!(pvb,
        x,
        e.virial["JL big"]./base_line;
        name="JuLIP"
    )

    base_line = e.virial["JL huge"][1]
    mi = min( e.virial["JL huge"]..., e.virial["AB huge"]...)
    ma = max( e.virial["JL huge"]..., e.virial["AB huge"]...)
    ylims = (mi, ma) ./ base_line
    pvh = lineplot(
        x,
        e.virial["AB huge"]./base_line;
        name="AtomsBase",
        title="AtomsBase virial 16_000 atoms",
        xlabel="1/nthreads",
        ylabel="relative time",
        ylim=ylims,
        xflip=true
    )
    lineplot!(pvh,
        x,
        e.virial["JL huge"]./base_line;
        name="JuLIP",
    )

    return [pb, p, pfb, pfh, pvb, pvh]
end