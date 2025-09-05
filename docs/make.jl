using ComodoFerrite
using Documenter

DocMeta.setdocmeta!(ComodoFerrite, :DocTestSetup, :(using ComodoFerrite); recursive=true)

makedocs(;
    modules=[ComodoFerrite],
    authors="Aminofa70 <amin.alibakhshi@upm.es> and contributors",
    sitename="ComodoFerrite.jl",
    format=Documenter.HTML(;
        canonical="https://COMODO-research.github.io/ComodoFerrite.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/COMODO-research/ComodoFerrite.jl",
    devbranch="main",
)
