# Candidate-band metrics table over a duty-cycled deployment.
#
# Builds a one-row-per-cycle table of soundscape band metrics for the Chapter 1
# exploratory analysis, then writes it to CSV for R. The band set is the
# candidate list from CANDIDATE_ACOUSTIC_BANDS.org: three standard shipping-noise
# one-third-octave levels (63/125/250 Hz) plus four broader exploratory bands.
#
# This is a general demonstration of `band_table`; the specific bands and
# parameters are the Chapter 1 test-run choices, not package defaults.
#
# Run:  julia --project=. examples/band_table_candidate_bands.jl <deployment_dir> [out.csv]

using EcoAcoustics
using DataFrames

# ── Candidate band set ────────────────────────────────────────────────────────
# Standard shipping TOLs. Narrow ranges pick exactly one ANSI-preferred centre
# each (63, 125, 250 Hz), so the edges are the canonical one-third-octave edges
# rather than hand-typed numbers.
tols = merge(tol_bands(62.0, 64.0),
             tol_bands(124.0, 126.0),
             tol_bands(249.0, 251.0))          # => :tol_63, :tol_125, :tol_250

# Broader exploratory bands (hand-picked edges). :b1_4k overlaps the dolphin
# whistle band (~2–20 kHz) — treat it as an exploratory vessel indicator only
# (self-contamination caveat, see CANDIDATE_ACOUSTIC_BANDS.org).
broad = Dict(:b100_200 => (100.0, 200.0),
             :b200_400 => (200.0, 400.0),
             :b400_800 => (400.0, 800.0),
             :b1_4k    => (1000.0, 4000.0))

bands = merge(tols, broad)                     # 7 candidate bands

# ── Inputs ────────────────────────────────────────────────────────────────────
deployment_dir = length(ARGS) >= 1 ? ARGS[1] : error(
    "usage: julia --project=. examples/band_table_candidate_bands.jl <deployment_dir> [out.csv]")
out_csv = length(ARGS) >= 2 ? ARGS[2] : "candidate_bands.csv"

# Every 5-min duty-cycle file in the deployment tree. `readdir(...; join=true)`
# returns full paths; filter to the audio extension you have.
files = filter(f -> endswith(lowercase(f), ".wav"),
               readdir(deployment_dir; join = true))
@info "band_table: found $(length(files)) file(s) in $deployment_dir"

# ── Per-era SM3M calibration ──────────────────────────────────────────────────
# The SM3M recorder profile currently encodes ONE era (−153 dB, pre-Apr-2018).
# Until per-era resolution lands (see ROADMAP.org), pass the correct scalar
# calibration explicitly. For a single-era batch, set it once here; for a mixed
# batch, branch inside `reader` on the file's timestamp.
#   pre-Apr-2018 : ScalarCalibration(-153f0)  (or omit `cal` to use the profile)
#   from Apr-2018 : ScalarCalibration(-165f0)
era_cal = ScalarCalibration(-153.0f0)

reader = f -> read_audio(f; recorder = "sm3m")

# ── Build the table ───────────────────────────────────────────────────────────
# One row per cycle. Percentiles are taken over 1 s PSD frames within each file.
# on_error=:skip drops any unreadable/corrupt file with a warning rather than
# aborting the whole 3-month run.
df = band_table(files;
                reader         = reader,
                bands          = bands,
                window_seconds = 1.0,
                cal            = era_cal,
                on_error       = :skip)

@info "band_table: produced $(nrow(df)) row(s), $(ncol(df)) column(s)"
show(first(df, 5); allcols = true)
println()

# ── Write for R ───────────────────────────────────────────────────────────────
# Requires CSV.jl in your environment: `using CSV; CSV.write(out_csv, df)`.
# Left commented so this example has no extra dependency.
#
# using CSV
# CSV.write(out_csv, df)
# @info "band_table: wrote $out_csv"
