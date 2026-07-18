# Process a whole SM3M deployment into a per-cycle candidate-band metrics CSV.
#
# One row per 5-minute duty-cycle file; columns are the Chapter 1 candidate bands
# (CANDIDATE_ACOUSTIC_BANDS.org) × metrics. Threaded across files, one shared FFT
# plan. Writes a plain CSV (no extra dependencies) ready for R.
#
# Run threaded (use all cores):
#   julia -t auto --project=. examples/process_deployment.jl <dir> [out.csv]
# e.g.
#   julia -t 24 --project=. examples/process_deployment.jl \
#         /home/robert/dissertation/data/T1C  T1C_bands.csv
#
# Calibration: SM3M profile default (-153 dB, the pre-Apr-2018 / 12 dB-gain era).
# For post-Apr-2018 files pass cal = ScalarCalibration(-165f0) to band_table.

using EcoAcoustics
using DataFrames
using Dates
using Base.Threads

# ── Candidate band set (CANDIDATE_ACOUSTIC_BANDS.org) ─────────────────────────
# Standard shipping TOLs (exact ISO/base-10 centres via narrow generator ranges)
# plus four broad exploratory bands. :b1_4k overlaps the whistle band — treat as
# an exploratory vessel indicator only.
const BANDS = merge(
    tol_bands(62.0, 64.0), tol_bands(124.0, 126.0), tol_bands(249.0, 251.0),
    Dict(:b100_200 => (100.0, 200.0), :b200_400 => (200.0, 400.0),
         :b400_800 => (400.0, 800.0), :b1_4k    => (1000.0, 4000.0)),
)

# ── Minimal CSV writer (avoids a CSV.jl dependency) ───────────────────────────
_fmt(x::DateTime)      = string(x)
_fmt(x::AbstractFloat) = string(round(x; digits = 4))
_fmt(x)                = string(x)

function write_csv(path::AbstractString, df::DataFrame)
    cols = names(df)
    open(path, "w") do io
        println(io, join(cols, ","))
        for row in eachrow(df)
            println(io, join((_fmt(row[c]) for c in cols), ","))
        end
    end
end

# ── Main ──────────────────────────────────────────────────────────────────────
function main()
    dir     = length(ARGS) >= 1 ? ARGS[1] : error(
        "usage: julia -t auto --project=. examples/process_deployment.jl <dir> [out.csv]")
    out_csv = length(ARGS) >= 2 ? ARGS[2] : "deployment_bands.csv"

    # Sorted so the SM3M yyyymmdd_HHMMSS filenames yield a time-ordered table.
    files = sort(filter(f -> endswith(lowercase(f), ".wav"),
                        readdir(dir; join = true)))
    @info "process_deployment: $(length(files)) files in $dir; threads = $(nthreads())"
    isempty(files) && error("no .wav files found in $dir")

    # One shared FFT plan for the whole run (all files are 48 kHz SM3M).
    fs   = Float64(read_audio(files[1]; recorder = "sm3m").fs)
    plan = make_spectrogram_plan(fs, 1.0)

    t = @elapsed df = band_table(files;
                                 reader          = f -> read_audio(f; recorder = "sm3m"),
                                 bands           = BANDS,
                                 window_seconds  = 1.0,
                                 fft_plan        = plan,
                                 parallel        = :threads,
                                 on_error        = :skip,
                                 progress        = true)

    write_csv(out_csv, df)
    @info "process_deployment: wrote $(nrow(df)) rows × $(ncol(df)) cols to $out_csv " *
          "in $(round(t; digits = 1)) s ($(round(t / max(nrow(df), 1); digits = 3)) s/file)"
    @info "process_deployment: time span $(minimum(df.start_time)) … $(maximum(df.start_time))"
end

main()
