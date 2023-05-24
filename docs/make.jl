using ACEMD
using Documenter

DocMeta.setdocmeta!(ACEMD, :DocTestSetup, :(using ACEMD); recursive=true)

makedocs(;
    modules=[ACEMD],
    authors="Teemu Järvinen <teemu.j.jarvinen@gmail.com> and contributors",
    repo="https://github.com/tjjarvinen/ACEMD.jl/blob/{commit}{path}#{line}",
    sitename="ACEMD.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://tjjarvinen.github.io/ACEMD.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Using ACE in Molly" => "molly.md",
    ],
)

deploydocs(;
    repo="github.com/tjjarvinen/ACEMD.jl",
    devbranch="main",
)
