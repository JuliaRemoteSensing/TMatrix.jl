using TMatrix
using Documenter

DocMeta.setdocmeta!(TMatrix, :DocTestSetup, :(using TMatrix); recursive = true)

makedocs(;
         modules = [TMatrix],
         authors = "Gabriel Wu <wuzihua@pku.edu.cn> and contributors",
         repo = "https://github.com/JuliaRemoteSensing/TMatrix.jl/blob/{commit}{path}#{line}",
         sitename = "TMatrix.jl",
         format = Documenter.HTML(;
                                  prettyurls = get(ENV, "CI", "false") == "true",
                                  canonical = "https://juliaremotesensing.github.io/TMatrix.jl",
                                  assets = String[]),
         pages = ["Home" => "index.md"
                  "API" => "api.md"])

deploydocs(; repo = "github.com/JuliaRemoteSensing/TMatrix.jl", devbranch = "main")
