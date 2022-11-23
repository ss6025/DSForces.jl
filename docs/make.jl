using DSForces
using Documenter

DocMeta.setdocmeta!(DSForces, :DocTestSetup, :(using DSForces); recursive=true)

makedocs(;
    modules=[DSForces],
    authors="Andrea Di Gioacchino",
    repo="https://github.com/adigioacchino/DSForces.jl/blob/{commit}{path}#{line}",
    sitename="DSForces.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://adigioacchino.github.io/DSForces.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/adigioacchino/DSForces.jl",
    devbranch="main",
)
