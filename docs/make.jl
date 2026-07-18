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
        # The single @autodocs API page (api.md) grows with every exported
        # symbol and exceeds Documenter's 200 KiB default. Raised to keep the
        # build green as the package grows. FOLLOW-UP: split api.md into
        # per-area autodocs pages (Pages = [...]) for better reader UX.
        size_threshold = 400 * 1024,
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
            "SPL"          => "explanations/spl.md",
            "SPL bands"    => "explanations/spl_bands.md",
            "LTSA"         => "explanations/ltsa.md",
            "Band metrics" => "explanations/band_metrics.md",
            "Glossary"     => "explanations/glossary.md",
        ],
        "API Reference" => "api.md"
    ],
)

deploydocs(;
    repo = "github.com/sonmar-org/EcoAcoustics.jl.git",
    branch = "gh-pages",
    devbranch = "dev",
)


