using TMatrix
using Documenter

DocMeta.setdocmeta!(TMatrix, :DocTestSetup, :(using TMatrix); recursive=true)

makedocs(;
    modules=[TMatrix],
    authors="Gabriel Wu <wuzihua@pku.edu.cn> and contributors",
    repo="https://github.com/lucifer1004/TMatrix.jl/blob/{commit}{path}#{line}",
    sitename="TMatrix.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://lucifer1004.github.io/TMatrix.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/lucifer1004/TMatrix.jl",
    devbranch="main",
)
