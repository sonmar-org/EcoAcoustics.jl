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
        "Filename parsing and metadata" => "metadata.md",
        "Audio sources and time-based access" => "sources.md",
        "Calibration" => [
            "Rockhopper" => "calibration/rockhopper.md",
        ],
        "Explanations" => [
            "Calibration"  => "explanations/calibration.md",
            "Spectrogram"  => "explanations/spectrogram.md",
            "PSD"          => "explanations/psd.md",
            "Glossary"     => "explanations/glossary.md",
        ],
    ],
)

deploydocs(;
    repo = "github.com/sonmar-org/EcoAcoustics.jl.git",
    branch = "gh-pages",
    devbranch = "dev",
)


