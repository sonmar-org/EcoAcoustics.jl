using EcoAcoustics
using Documenter

DocMeta.setdocmeta!(EcoAcoustics, :DocTestSetup, :(using EcoAcoustics); recursive=true)

makedocs(;
    modules=[EcoAcoustics],
    authors="Sonmar <info@sonmar.org>",
    sitename="EcoAcoustics.jl",
    format=Documenter.HTML(;
        edit_link="dev",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
