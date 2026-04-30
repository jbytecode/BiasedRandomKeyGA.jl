using Documenter, BiasedRandomKeyGA

makedocs(
    format=Documenter.HTML(
        prettyurls=get(ENV, "CI", nothing) == "true",
        collapselevel=2,
        # assets = ["assets/favicon.ico", "assets/extra_styles.css"],
    ),
    sitename="BiasedRandomKeyGA.jl",
    authors="Mehmet Hakan Satman",
    pages=[
        "API Reference" => "apireference.md",
    ],
)


deploydocs(repo="github.com/jbytecode/BiasedRandomKeyGA.jl")