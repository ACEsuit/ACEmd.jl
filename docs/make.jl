using ACEapi
using Documenter

DocMeta.setdocmeta!(ACEapi, :DocTestSetup, :(using ACEapi); recursive=true)

makedocs(;
    modules=[ACEapi],
    authors="Teemu Järvinen <teemu.j.jarvinen@gmail.com> and contributors",
    repo="https://github.com/tjjarvinen/ACEapi.jl/blob/{commit}{path}#{line}",
    sitename="ACEapi.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://tjjarvinen.github.io/ACEapi.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/tjjarvinen/ACEapi.jl",
    devbranch="main",
)
