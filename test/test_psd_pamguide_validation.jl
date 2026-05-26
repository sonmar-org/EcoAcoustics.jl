using Test
using EcoAcoustics
using Logging
using Statistics
using DelimitedFiles

# ─── PAMGuide PSD cross-validation ───────────────────────────────────────────
#
# Compares EcoAcoustics.jl PSD output against PAMGuide reference output for
# two DMON2 recordings (2 kHz, full spectrum 0–1000 Hz) and one synthetic
# white-noise file (48 kHz; PAMGuide CSV covers 10–1000 Hz subset).
#
# Normalization reference: ADEON DPS §2.2.1 (Ainslie et al. 2018). The DPS
# formula uses a pre-normalised Hann window (sqrt(8/3) factor); this is
# algebraically identical to the Merchant 2015 Eq.(1) convention used here.
# PAMGuide implements the same convention and serves as the independent check.
# See DD-24 for the 30 Hz lower-frequency cut and docs/design_decisions.md for
# the full normalization equivalence argument.
#
# Required files in test/validation/pamguide/ (not version-controlled):
#   230306_152540_VesselPassage_dmon2.wav
#   230306_201311_CallingPeriod_dmon2.wav
#   WhiteNoise_10s_48kHz_+-0.5.wav
#   VesselPassage_dmon2_PSD.csv
#   CallingPeriod_dmon2_PSD.csv
#   WhiteNoise_10s_48kHz_+-0.5_PSD.csv
#
# See test/validation/pamguide/README.md for PAMGuide settings and tolerances.

const _VAL_DIR = joinpath(@__DIR__, "validation", "pamguide")

const _VP_WAV  = joinpath(_VAL_DIR, "230306_152540_VesselPassage_dmon2.wav")
const _CP_WAV  = joinpath(_VAL_DIR, "230306_201311_CallingPeriod_dmon2.wav")
const _WN_WAV  = joinpath(_VAL_DIR, "WhiteNoise_10s_48kHz_+-0.5.wav")
const _VP_CSV  = joinpath(_VAL_DIR, "VesselPassage_dmon2_PSD.csv")
const _CP_CSV  = joinpath(_VAL_DIR, "CallingPeriod_dmon2_PSD.csv")
const _WN_CSV  = joinpath(_VAL_DIR, "WhiteNoise_10s_48kHz_+-0.5_PSD.csv")

let missing_files = filter(!isfile, [_VP_WAV, _CP_WAV, _WN_WAV,
                                     _VP_CSV, _CP_CSV, _WN_CSV])
    if !isempty(missing_files)
        @info "PAMGuide PSD validation: skipping — missing files: " *
              join(basename.(missing_files), ", ")
    end
    global _PAMGUIDE_PSD_SKIP = !isempty(missing_files)
end

if !_PAMGUIDE_PSD_SKIP

# ─── Tolerances ───────────────────────────────────────────────────────────────

const _MEAN_THRESHOLD = 0.05    # dB — fail if masked Δmean exceeds this
const _MAX_THRESHOLD  = 0.5     # dB — threshold for per-bin anomaly classification
const _MAX_EXCEED_PCT = 1.0     # %  — fail if more than this fraction exceeds _MAX_THRESHOLD
const _NOISE_FLOOR_DB = -120.0  # dB — bins below this are excluded from statistics

# ─── Calibration ──────────────────────────────────────────────────────────────

# DMON2: derive system sensitivity from the registered profile components so
# the arithmetic is auditable here. Source: src/recorders/dmon2.jl.
const _DMON2_SENS_DB = let cp = EcoAcoustics.CALIBRATION_PROFILES["dmon2"]
    cp.sensitivity + cp.preamp_gain + cp.board_gain + 20.0 * log10(1.0 / cp.Vadc_0pk)
end
const _DMON2_CAL = ScalarCalibration(Float32(_DMON2_SENS_DB))

# "test" recorder for the synthetic white-noise file — not a registered package
# profile. Components: sensitivity=-200, preamp=0, board=0, Vadc_0pk=2.0 V.
const _TEST_CAL = ScalarCalibration(Float32(-200.0 + 20.0 * log10(1.0 / 2.0)))

# ─── Helpers ──────────────────────────────────────────────────────────────────

# Load a PAMGuide PSD CSV. Header row format:
#   <metadata_field>, 10, 11, ..., 1000
# Data rows: first column is frame-centre time (s); remaining columns are
# PSD in dB re 1 µPa²/Hz at the labelled integer-Hz bin centres.
#
# Returns (pam_dB, freq_hz) where:
#   pam_dB  :: Matrix{Float64}  shape freq × time (transposed from CSV time × freq)
#   freq_hz :: Vector{Int}      Hz label for each row of pam_dB
function _load_pamguide_psd(path::String)
    data, hdr = readdlm(path, ',', header = true)
    # hdr is a 1×N matrix. First element is a metadata field ("2211",
    # significance unknown); remaining elements are integer Hz labels.
    freq_hz = parse.(Int, vec(hdr)[2:end])
    # Drop the time column (column 1); remaining columns are PSD dB values.
    pam_dB = Float64.(data[:, 2:end])'   # transpose: time×freq → freq×time
    return pam_dB, freq_hz
end

# Return the UnitRange of row indices into psd_freqs that spans hz_labels.
# Requires hz_labels to be consecutive integers and psd_freqs to have 1 Hz
# resolution — both guaranteed for a 1 s analysis window at any integer fs.
function _freq_slice(psd_freqs::AbstractVector, hz_labels::Vector{Int})
    i_lo = searchsortedfirst(psd_freqs, Float64(first(hz_labels)))
    i_hi = searchsortedfirst(psd_freqs, Float64(last(hz_labels)))
    return i_lo:i_hi
end

# Compare ea_dB (freq × time) against pam_dB element-wise.
# freq_hz must be the integer-Hz labels for each row of ea_dB / pam_dB.
# Only bins with freq_hz[i] >= freq_lo are included in all statistics.
# freq_lo defaults to 0.0 (include all bins); pass 30.0 to exclude the
# lowest-frequency bins where symmetric vs periodic Hann window differences
# produce large per-frame spectral-leakage residuals (see DD-24).
# Applies a noise-floor mask (_NOISE_FLOOR_DB) before computing statistics.
# Logs a summary line via @info; samples 100 random interior bins from the
# retained rows and classifies any anomaly > _MAX_THRESHOLD as isolated
# (0 neighbors also elevated) or systematic (≥ 1 neighbor above
# _MEAN_THRESHOLD) via @warn.
# Returns (Δmean, Δstd, Δmax, percent_exceed) on the masked residuals.
function _compare_psd(label::String,
                      ea_dB::AbstractMatrix{Float64},
                      pam_dB::AbstractMatrix{Float64},
                      freq_hz::Vector{Int};
                      freq_lo::Float64 = 0.0)
    @assert size(ea_dB) == size(pam_dB) "Size mismatch for $label: " *
        "EcoAcoustics $(size(ea_dB)) vs PAMGuide $(size(pam_dB))"
    @assert length(freq_hz) == size(ea_dB, 1) "freq_hz length $(length(freq_hz)) " *
        "must match row count $(size(ea_dB, 1))"

    row_mask  = freq_hz .>= freq_lo        # Bool vector over frequency rows
    ea_sub    = ea_dB[row_mask, :]
    pam_sub   = pam_dB[row_mask, :]
    freq_sub  = freq_hz[row_mask]

    full_Δ = abs.(ea_sub .- pam_sub)

    mask           = (ea_sub  .>= _NOISE_FLOOR_DB) .& (pam_sub .>= _NOISE_FLOOR_DB)
    masked_Δ       = full_Δ[mask]
    Δmean          = mean(masked_Δ)
    Δstd           = std(masked_Δ)
    Δmax           = maximum(masked_Δ)
    n_exceed       = count(>(_MAX_THRESHOLD), masked_Δ)
    percent_exceed = 100.0 * n_exceed / length(masked_Δ)

    @info "PSD vs PAMGuide [$label]" freq_lo=freq_lo Δmean=round(Δmean; digits=4) Δstd=round(Δstd; digits=4) Δmax=round(Δmax; digits=4) percent_exceeding=round(percent_exceed; digits=2)

    # Neighborhood analysis: sample 100 random interior bins from the retained
    # rows, staying 10 bins away from each edge so every sampled bin has a
    # full 3×3 neighborhood.
    n_freq, n_time = size(full_Δ)
    if n_freq >= 21 && n_time >= 21
        sample_f = rand(11:(n_freq - 10), 100)
        sample_t = rand(11:(n_time - 10), 100)
        for (fi, ti) in zip(sample_f, sample_t)
            δ = full_Δ[fi, ti]
            δ > _MAX_THRESHOLD || continue
            n_neighbors = sum(
                full_Δ[fi + df, ti + dt] > _MEAN_THRESHOLD
                for df in -1:1, dt in -1:1 if !(df == 0 && dt == 0)
            )
            if n_neighbors == 0
                @warn "[$label] isolated anomaly at ($fi, $ti): " *
                      "Δ=$(round(δ; digits=3)) dB  " *
                      "EA=$(round(ea_sub[fi,ti]; digits=2)) dB  " *
                      "PG=$(round(pam_sub[fi,ti]; digits=2)) dB"
            else
                @warn "[$label] systematic anomaly at ($fi, $ti): " *
                      "Δ=$(round(δ; digits=3)) dB  " *
                      "$n_neighbors/8 neighbors also exceed $(_MEAN_THRESHOLD) dB"
            end
        end
    end

    return Δmean, Δstd, Δmax, percent_exceed
end

# ─── VesselPassage ───────────────────────────────────────────────────────────

@testset "PAMGuide PSD validation: VesselPassage_dmon2.wav" begin
    audio = @test_logs min_level=Logging.Error read_audio(_VP_WAV; recorder = "dmon2")
    psd   = compute_psd(audio; window_seconds = 1.0, overlap_fraction = 0.5,
                        cal = _DMON2_CAL)
    ea_dB = to_dB(psd)

    pam_dB, freq_hz = _load_pamguide_psd(_VP_CSV)
    ea_slice = ea_dB[_freq_slice(psd.freqs, freq_hz), :]

    Δmean, Δstd, Δmax, pct = _compare_psd("VesselPassage_dmon2", ea_slice, pam_dB,
                                           freq_hz; freq_lo = 30.0)

    @test Δmean <= _MEAN_THRESHOLD
    @test pct   <= _MAX_EXCEED_PCT
end

# ─── CallingPeriod ────────────────────────────────────────────────────────────

@testset "PAMGuide PSD validation: CallingPeriod_dmon2.wav" begin
    audio = @test_logs min_level=Logging.Error read_audio(_CP_WAV; recorder = "dmon2")
    psd   = compute_psd(audio; window_seconds = 1.0, overlap_fraction = 0.5,
                        cal = _DMON2_CAL)
    ea_dB = to_dB(psd)

    pam_dB, freq_hz = _load_pamguide_psd(_CP_CSV)
    ea_slice = ea_dB[_freq_slice(psd.freqs, freq_hz), :]

    Δmean, Δstd, Δmax, pct = _compare_psd("CallingPeriod_dmon2", ea_slice, pam_dB,
                                           freq_hz; freq_lo = 30.0)

    @test Δmean <= _MEAN_THRESHOLD
    @test pct   <= _MAX_EXCEED_PCT
end

# ─── WhiteNoise ───────────────────────────────────────────────────────────────
# 48 kHz synthetic signal; EcoAcoustics PSD spans 0–24 000 Hz while the
# PAMGuide CSV covers 10–1 000 Hz. _freq_slice selects the matching rows.

@testset "PAMGuide PSD validation: WhiteNoise_10s_48kHz_+-0.5.wav" begin
    audio = @test_logs min_level=Logging.Error read_audio(_WN_WAV)
    psd   = compute_psd(audio; window_seconds = 1.0, overlap_fraction = 0.5,
                        cal = _TEST_CAL)
    ea_dB = to_dB(psd)

    pam_dB, freq_hz = _load_pamguide_psd(_WN_CSV)
    ea_slice = ea_dB[_freq_slice(psd.freqs, freq_hz), :]

    Δmean, Δstd, Δmax, pct = _compare_psd("WhiteNoise", ea_slice, pam_dB,
                                           freq_hz; freq_lo = 30.0)

    @test Δmean <= _MEAN_THRESHOLD
    @test pct   <= _MAX_EXCEED_PCT
end

end  # !_PAMGUIDE_PSD_SKIP
