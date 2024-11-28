using GeneString
using Documenter

DocMeta.setdocmeta!(GeneString, :DocTestSetup, :(using GeneString); recursive=true)

makedocs(;
    modules=[GeneString],
    authors="Jonathan Miller jonathan.miller@fieldofnodes.com",
    sitename="GeneString.jl",
    format=Documenter.HTML(;
        canonical="https://fieldofnodes.github.io/GeneString.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/fieldofnodes/GeneString.jl",
    devbranch="main",
)
