using AstroEOMs
using Documenter

DocMeta.setdocmeta!(AstroEOMs, :DocTestSetup, :(using AstroEOMs); recursive=true)

makedocs(;
    modules=[AstroEOMs],
    authors="Grant Hecht",
    repo="https://github.com/GrantHecht/AstroEOMs.jl/blob/{commit}{path}#{line}",
    sitename="AstroEOMs.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://GrantHecht.github.io/AstroEOMs.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/GrantHecht/AstroEOMs.jl",
    devbranch="main",
)
