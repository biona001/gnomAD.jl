using EasyLD
using Documenter

makedocs(
    sitename = "EasyLD.jl",
    format = Documenter.HTML(),
    modules = [EasyLD],
    pages = [
        "API" => "man/api.md",
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo   = "github.com/biona001/EasyLD.jl.git",
    target = "build"
)
