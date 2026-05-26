using Test
using EcoAcoustics
using Logging
using Statistics
using DelimitedFiles

# ─── PAMGuide SPL cross-validation ───────────────────────────────────────────
#
# Compares EcoAcoustics.jl SPL output against two independent references:
#
#   DMON2 recordings — reference is derived by integrating the PAMGuide PSD
#   CSV files to broadband SPL via _pam_spl, bypassing PAMGuide's buggy
#   broadband SPL path (PG_DFT.m line 127, which omits the 1/B ENBW factor).
#   Each CSV contains per-frame PSD in dB re 1 µPa²/Hz; integrating with
#   df = 1 Hz and taking the energetic mean over frames gives a reference
#   independent of EcoAcoustics's PSD pipeline.
#
#   WhiteNoise — reference is the time-domain RMS of the WAV signal, converted
#   to SPL using the calibration formula directly. Bypasses both the PSD
#   pipeline and PAMGuide CSV entirely. For a Hann-windowed PSD, Parseval's
#   theorem guarantees that sum(psd) × df = mean(signal²), so PSD integration
#   and time-domain RMS should agree to within numerical noise.
#
#   VesselPassage FLAC — same PAMGuide CSV reference as the WAV test, verifying
#   that FLAC and WAV decoding produce identical SPL.
#
# Required files in test/validation/pamguide/ (not version-controlled):
#   230306_152540_VesselPassage_dmon2.wav
#   230306_152540_VesselPassage_dmon2.flac
#   230306_201311_CallingPeriod_dmon2.wav
#   WhiteNoise_10s_48kHz_+-0.5.wav
#   VesselPassage_dmon2_PSD.csv
#   CallingPeriod_dmon2_PSD.csv
#
# See test/validation/pamguide/README.md for calibration components,
# tolerance rationale, and the PAMGuide ENBW bug documentation.

const _SPL_VAL_DIR = joinpath(@__DIR__, "validation", "pamguide")

const _SPL_VP_WAV  = joinpath(_SPL_VAL_DIR, "230306_152540_VesselPassage_dmon2.wav")
const _SPL_VP_FLAC = joinpath(_SPL_VAL_DIR, "230306_152540_VesselPassage_dmon2.flac")
const _SPL_CP_WAV  = joinpath(_SPL_VAL_DIR, "230306_201311_CallingPeriod_dmon2.wav")
const _SPL_WN_WAV  = joinpath(_SPL_VAL_DIR, "WhiteNoise_10s_48kHz_+-0.5.wav")
const _SPL_VP_CSV  = joinpath(_SPL_VAL_DIR, "VesselPassage_dmon2_PSD.csv")
const _SPL_CP_CSV  = joinpath(_SPL_VAL_DIR, "CallingPeriod_dmon2_PSD.csv")

let missing_files = filter(!isfile, [_SPL_VP_WAV, _SPL_VP_FLAC,
                                     _SPL_CP_WAV, _SPL_WN_WAV,
                                     _SPL_VP_CSV, _SPL_CP_CSV])
    if !isempty(missing_files)
        @info "PAMGuide SPL validation: skipping — missing files: " *
              join(basename.(missing_files), ", ")
    end
    global _PAMGUIDE_SPL_SKIP = !isempty(missing_files)
end

# Integrate a PAMGuide PSD CSV (freq × time, dB re µPa²/Hz) to broadband SPL.
# df = 1 Hz is guaranteed by the 1-second analysis window at any integer fs.
# Restricts to bins f_lo ≤ freq_hz ≤ f_hi, then returns the energetic mean
# (linear average of per-frame power, then converted to dB) over all frames.
# Directly comparable to BandSPL.mean_dB from compute_spl (DD-22).
function _pam_spl(pam_dB::AbstractMatrix{Float64}, freq_hz::Vector{Int};
                  f_lo::Int = first(freq_hz), f_hi::Int = last(freq_hz))
    rows = (freq_hz .>= f_lo) .& (freq_hz .<= f_hi)
    sub  = pam_dB[rows, :]
    per_frame_power = [sum(10 .^ (sub[:, j] ./ 10)) for j in axes(sub, 2)]
    return 10.0 * log10(mean(per_frame_power))
end

if !_PAMGUIDE_SPL_SKIP

# ─── Calibration ──────────────────────────────────────────────────────────────

# DMON2: derive from registered profile components so the arithmetic is auditable.
const _SPL_DMON2_CAL = let cp = EcoAcoustics.CALIBRATION_PROFILES["dmon2"]
    ScalarCalibration(Float32(
        cp.sensitivity + cp.preamp_gain + cp.board_gain + 20.0 * log10(1.0 / cp.Vadc_0pk)
    ))
end

# "test" recorder for white-noise synthetic signal — not a registered profile.
# Components: sensitivity=-200, preamp=0, board=0, Vadc_0pk=2.0 V.
const _SPL_TEST_CAL = ScalarCalibration(Float32(-200.0 + 20.0 * log10(1.0 / 2.0)))

# ─── Shared PSD computation ───────────────────────────────────────────────────
#
# Each audio file is loaded once and its PSD computed once. Multiple testsets
# reuse the same PSDResult (cheap band integration) rather than re-reading
# and re-transforming the full audio signal per testset.

audio_vp_wav  = @test_logs min_level=Logging.Error read_audio(_SPL_VP_WAV;  recorder = "dmon2")
audio_vp_flac = @test_logs min_level=Logging.Error read_audio(_SPL_VP_FLAC; recorder = "dmon2")
audio_cp      = @test_logs min_level=Logging.Error read_audio(_SPL_CP_WAV;  recorder = "dmon2")
audio_wn      = @test_logs min_level=Logging.Error read_audio(_SPL_WN_WAV)

psd_vp_wav  = compute_psd(audio_vp_wav;  window_seconds = 1.0, overlap_fraction = 0.5,
                           cal = _SPL_DMON2_CAL)
psd_vp_flac = compute_psd(audio_vp_flac; window_seconds = 1.0, overlap_fraction = 0.5,
                           cal = _SPL_DMON2_CAL)
psd_cp      = compute_psd(audio_cp;      window_seconds = 1.0, overlap_fraction = 0.5,
                           cal = _SPL_DMON2_CAL)
psd_wn      = compute_psd(audio_wn;      window_seconds = 1.0, overlap_fraction = 0.5,
                           cal = _SPL_TEST_CAL)

# ─── PAMGuide CSV references (DMON2 only) ─────────────────────────────────────
#
# Load once here; each testset calls _pam_spl with an appropriate band.
# The WhiteNoise CSV covers only 10–1000 Hz of the 48 kHz spectrum, so it
# cannot serve as a full-broadband reference — that test uses RMS instead.

_pam_vp_dB, _pam_vp_hz = let
    data, hdr = readdlm(_SPL_VP_CSV, ',', header = true)
    hz = parse.(Int, vec(hdr)[2:end])
    dB = Float64.(data[:, 2:end])'
    dB, hz
end

_pam_cp_dB, _pam_cp_hz = let
    data, hdr = readdlm(_SPL_CP_CSV, ',', header = true)
    hz = parse.(Int, vec(hdr)[2:end])
    dB = Float64.(data[:, 2:end])'
    dB, hz
end

# ─── VesselPassage WAV — broadband ───────────────────────────────────────────

@testset "PAMGuide SPL: VesselPassage_dmon2.wav broadband" begin
    spl    = compute_spl(psd_vp_wav)
    actual = spl.bands[:broadband].mean_dB
    ref    = _pam_spl(_pam_vp_dB, _pam_vp_hz)
    @info "SPL [VesselPassage_dmon2.wav broadband] EA vs PAMGuide CSV" ea_dB=round(actual; digits=3) pam_csv_dB=round(ref; digits=3) Δ=round(abs(actual - ref); digits=3)
    @test abs(actual - ref) <= 0.05
end

# ─── VesselPassage FLAC — broadband ──────────────────────────────────────────
# Uses the same CSV reference as the WAV test. Agreement between WAV and FLAC
# confirms format-independent decoding; agreement with the CSV confirms
# correct normalization.

@testset "PAMGuide SPL: VesselPassage_dmon2.flac broadband" begin
    spl    = compute_spl(psd_vp_flac)
    actual = spl.bands[:broadband].mean_dB
    ref    = _pam_spl(_pam_vp_dB, _pam_vp_hz)
    @info "SPL [VesselPassage_dmon2.flac broadband] EA vs PAMGuide CSV" ea_dB=round(actual; digits=3) pam_csv_dB=round(ref; digits=3) Δ=round(abs(actual - ref); digits=3)
    @test abs(actual - ref) <= 0.05
end

# ─── VesselPassage WAV — 500–1000 Hz band ────────────────────────────────────

@testset "PAMGuide SPL: VesselPassage_dmon2.wav 500–1000 Hz" begin
    spl    = compute_spl(psd_vp_wav; bands = Dict(:band => (500.0, 1000.0)))
    actual = spl.bands[:band].mean_dB
    ref    = _pam_spl(_pam_vp_dB, _pam_vp_hz; f_lo = 500, f_hi = 1000)
    @info "SPL [VesselPassage_dmon2.wav 500–1000 Hz] EA vs PAMGuide CSV" ea_dB=round(actual; digits=3) pam_csv_dB=round(ref; digits=3) Δ=round(abs(actual - ref); digits=3)
    @test abs(actual - ref) <= 0.1
end

# ─── CallingPeriod WAV — broadband ───────────────────────────────────────────

@testset "PAMGuide SPL: CallingPeriod_dmon2.wav broadband" begin
    spl    = compute_spl(psd_cp)
    actual = spl.bands[:broadband].mean_dB
    ref    = _pam_spl(_pam_cp_dB, _pam_cp_hz)
    @info "SPL [CallingPeriod_dmon2.wav broadband] EA vs PAMGuide CSV" ea_dB=round(actual; digits=3) pam_csv_dB=round(ref; digits=3) Δ=round(abs(actual - ref); digits=3)
    @test abs(actual - ref) <= 0.05
end

# ─── WhiteNoise — broadband ───────────────────────────────────────────────────
# 48 kHz file; broadband = 10 Hz to Nyquist (24 000 Hz).
# The PAMGuide CSV covers only 10–1000 Hz and cannot serve as a broadband
# reference. Instead, compare against the time-domain RMS of the decoded
# signal: SPL_td = 20·log10(RMS_normalized) − S. For a Hann-windowed PSD,
# Parseval's theorem guarantees sum(psd)·df = mean(signal²), so PSD
# integration and RMS should agree to within numerical noise (~0.01 dB).

@testset "PAMGuide SPL: WhiteNoise_10s_48kHz_+-0.5.wav broadband" begin
    spl    = compute_spl(psd_wn)
    actual = spl.bands[:broadband].mean_dB
    S      = Float64(_SPL_TEST_CAL.system_sensitivity_dB)
    rms_ref = 20.0 * log10(sqrt(mean(audio_wn.sig .^ 2))) - S
    @info "SPL [WhiteNoise broadband] EA (PSD) vs time-domain RMS" psd_dB=round(actual; digits=3) rms_dB=round(rms_ref; digits=3) Δ=round(abs(actual - rms_ref); digits=3)
    @test abs(actual - rms_ref) <= 0.1
end

end  # !_PAMGUIDE_SPL_SKIP
