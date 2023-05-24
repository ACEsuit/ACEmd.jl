using ACEmd
using Documenter

DocMeta.setdocmeta!(ACEmd, :DocTestSetup, :(using ACEmd); recursive=true)

makedocs(;
    modules=[ACEmd],
    authors="Teemu Järvinen <teemu.j.jarvinen@gmail.com> and contributors",
    repo="https://github.com/acesuit/ACEmd.jl/blob/{commit}{path}#{line}",
    sitename="ACEmd.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://acesuit.github.io/ACEmd.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Using ACE in Molly" => "molly.md",
    ],
)

deploydocs(;
    repo="github.com/acesuit/ACEmd.jl",
    devbranch="main",
)
