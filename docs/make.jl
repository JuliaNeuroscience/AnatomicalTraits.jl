using AnatomicalTraits
using Documenter

DocMeta.setdocmeta!(AnatomicalTraits, :DocTestSetup, :(using AnatomicalTraits); recursive=true)

makedocs(;
    modules=[AnatomicalTraits],
    authors="Zachary P. Christensen <zchristensen7@gmail.com> and contributors",
    repo="https://github.com/Tokazama/AnatomicalTraits.jl/blob/{commit}{path}#{line}",
    sitename="AnatomicalTraits.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Tokazama.github.io/AnatomicalTraits.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Tokazama/AnatomicalTraits.jl",
)
