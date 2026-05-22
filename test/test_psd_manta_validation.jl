using Test
using Dates
using EcoAcoustics
using Statistics

# ─── MANTA / PAMGuide cross-validation ───────────────────────────────────────
#
# Compares EcoAcoustics.jl PSD output against MANTA reference output for three
# one-minute windows extracted from two 120-second Rockhopper clips (197368 Hz).
#
# Shipped audio clips (in test/validation/rockhopper/):
#   rh_clip_004600Z.flac   — 120 s starting 17-Sep-2023 00:46:00 UTC
#   rh_clip_184620Z.flac   — 120 s starting 17-Sep-2023 18:46:20 UTC
#
# Required but NOT shipped (must be placed in test/validation/rockhopper/):
#   139635MD01_RH416_A7M.1.197368_20230917_DAILY_MILLIDEC_MinRes.csv
#   Rockhopper_TF_calibration.csv
#   Rockhopper_TF.xlsx
#
# See test/validation/rockhopper/README.md for MANTA parameters and tolerances.

const _VAL_DIR    = joinpath(@__DIR__, "validation", "rockhopper")
const _CLIP1_PATH = joinpath(_VAL_DIR, "rh_clip_004600Z.flac")
const _CLIP2_PATH = joinpath(_VAL_DIR, "rh_clip_184620Z.flac")
const _CSV_PATH   = joinpath(_VAL_DIR, "139635MD01_RH416_A7M.1.197368_20230917_DAILY_MILLIDEC_MinRes.csv")
const _TF_CSV_PATH = joinpath(_VAL_DIR, "Rockhopper_TF_calibration.csv")
const _TF_XLS_PATH = joinpath(_VAL_DIR, "Rockhopper_TF.xlsx")

let missing_files = filter(!isfile, [_CLIP1_PATH, _CLIP2_PATH, _CSV_PATH,
                                     _TF_CSV_PATH, _TF_XLS_PATH])
    if !isempty(missing_files)
        @info "MANTA validation: skipping — missing files: " *
              join(basename.(missing_files), ", ")
    end
    global _MANTA_SKIP = !isempty(missing_files)
end

if !_MANTA_SKIP

# ─── helpers ─────────────────────────────────────────────────────────────────

# Construct a new Audiodata from a contiguous sample range of an existing clip,
# preserving all metadata. dt_offset shifts starttime by the window's sample
# offset converted to wall time.
function _audio_window(clip::Audiodata, range::UnitRange{Int}, dt_offset::Period)
    Audiodata(clip.sig[range], clip.fs, clip.starttime + dt_offset;
              calibration = clip.calibration,
              recorder    = clip.metadata.recorder,
              recorder_id = clip.metadata.recorder_id,
              site_id     = clip.metadata.site_id,
              lat         = clip.metadata.lat,
              lon         = clip.metadata.lon,
              timezone    = clip.metadata.timezone)
end

# Load three rows from the MANTA daily CSV by exact timestamp string match.
# Returns a Dict{String, Vector{Float64}} mapping each timestamp to 391 dB
# values covering integer Hz 10–400 (CSV columns 13–403, 0-indexed Hz = col-3).
# Aborts if any timestamp is missing or its N-seconds value is not 60.
function _load_manta_rows(path, timestamps)
    lines = filter(!isempty ∘ strip, readlines(path))
    # Build timestamp → row-parts lookup from data rows (skip header, line 1).
    row_map = Dict{String, Vector{String}}()
    for line in lines[2:end]
        parts = String.(strip.(split(line, ',')))
        isempty(parts[1]) && continue
        row_map[parts[1]] = parts
    end
    result = Dict{String, Vector{Float64}}()
    for ts in timestamps
        haskey(row_map, ts) || error(
            "MANTA validation: timestamp \"$ts\" not found in $(basename(path))")
        parts = row_map[ts]
        n_sec = parse(Float64, parts[2])
        n_sec == 60.0 || error(
            "MANTA validation: N_seconds = $n_sec for \"$ts\"; " *
            "expected 60 — partial-minute aggregation invalidates comparison.")
        # Columns 13–403 (1-indexed) = integer Hz 10–400.
        result[ts] = parse.(Float64, parts[13:403])
    end
    return result
end

# ─── Load clips and extract three one-minute windows ─────────────────────────

clip1 = read_audio(_CLIP1_PATH; recorder = "rockhopper")
clip2 = read_audio(_CLIP2_PATH; recorder = "rockhopper")

fs   = clip1.fs                       # 197368.0f0
n60  = round(Int, Float64(fs) * 60)   # 11_842_080 — exactly representable
n120 = round(Int, Float64(fs) * 120)  # 23_684_160 — even, exactly representable
n40  = round(Int, Float64(fs) * 40)   #  7_894_720 — exactly representable
n100 = round(Int, Float64(fs) * 100)  # 19_736_800 — even, exactly representable

# Clip 1 starts 00:46:00 UTC.  First minute = samples 1:n60.
audio_0046 = _audio_window(clip1, 1:n60,       Second(0))
# Second minute = samples n60+1:n120.
audio_0047 = _audio_window(clip1, n60+1:n120,  Second(60))
# Clip 2 starts 18:46:20 UTC.  Minute 18:47 begins 40 s in.
audio_1847 = _audio_window(clip2, n40+1:n100,  Second(40))

# ─── Compute per-minute PSDs (no explicit cal — auto-resolves Rockhopper TF) ─

psd_0046 = compute_psd(audio_0046; window_seconds = 1.0, overlap_fraction = 0.5)
psd_0047 = compute_psd(audio_0047; window_seconds = 1.0, overlap_fraction = 0.5)
psd_1847 = compute_psd(audio_1847; window_seconds = 1.0, overlap_fraction = 0.5)

# Timestamp strings must match the CSV's first-column format exactly.
# MANTA daily CSV uses "DD-Mon-YYYY HH:MM:SS".
const _TS_0046 = "17-Sep-2023 00:46:00"
const _TS_0047 = "17-Sep-2023 00:47:00"
const _TS_1847 = "17-Sep-2023 18:47:00"

# Load all three MANTA reference rows once; each testset reads from this dict.
manta_rows = _load_manta_rows(_CSV_PATH, [_TS_0046, _TS_0047, _TS_1847])

# ─── Testsets ────────────────────────────────────────────────────────────────

@testset "MANTA validation: shipped TF matches validation copy" begin
    shipped_tf = RockhopperProfile().tf
    val_tf = load_tf_calcurves(
        joinpath(_VAL_DIR, "Rockhopper_TF_calibration.csv");
        vmax_peak_V = 5.0
    )
    @test length(shipped_tf.tf_dB) == length(val_tf.tf_dB)
    @test shipped_tf.frequency == val_tf.frequency
    @test maximum(abs, shipped_tf.tf_dB .- val_tf.tf_dB) < 1e-10
end

@testset "MANTA validation: minute 00:46:00" begin
    ours  = average_psd(psd_0046)
    diffs = (10 .* log10.(ours[11:401])) .- manta_rows[_TS_0046]
    @info "MANTA 00:46:00" maximum_abs=maximum(abs, diffs) mean_abs=mean(abs, diffs) std=std(diffs)
    @test maximum(abs, diffs) < 0.1
    @test mean(abs, diffs)    < 0.05
    @test std(diffs)          < 0.05
end

@testset "MANTA validation: minute 00:47:00" begin
    ours  = average_psd(psd_0047)
    diffs = (10 .* log10.(ours[11:401])) .- manta_rows[_TS_0047]
    @info "MANTA 00:47:00" maximum_abs=maximum(abs, diffs) mean_abs=mean(abs, diffs) std=std(diffs)
    @test maximum(abs, diffs) < 0.1
    @test mean(abs, diffs)    < 0.05
    @test std(diffs)          < 0.05
end

@testset "MANTA validation: minute 18:47:00" begin
    ours  = average_psd(psd_1847)
    diffs = (10 .* log10.(ours[11:401])) .- manta_rows[_TS_1847]
    @info "MANTA 18:47:00" maximum_abs=maximum(abs, diffs) mean_abs=mean(abs, diffs) std=std(diffs)
    @test maximum(abs, diffs) < 0.1
    @test mean(abs, diffs)    < 0.05
    @test std(diffs)          < 0.05
end

end  # !_MANTA_SKIP
