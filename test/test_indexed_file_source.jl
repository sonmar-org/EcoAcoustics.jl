using Test
using EcoAcoustics
using Dates
import Arrow
import WAV

# ─── Shared test fixture ──────────────────────────────────────────────────────
#
# Two synthetic WAV files, each 1 second long at 1000 Hz (1000 samples):
#
#   file1: starts at T1 = 2023-01-01 12:00:00, signal value = 0.5
#   file2: starts at T2 = 2023-01-01 12:00:02, signal value = -0.5
#
# The gap between file1's end (12:00:01) and file2's start (12:00:02) is
# exactly 1 second. The total span is 3 seconds.

const _T1  = DateTime(2023, 1, 1, 12, 0, 0)
const _T2  = DateTime(2023, 1, 1, 12, 0, 2)
const _FS  = 1000  # Hz (Int; converted to Float32 by build_index)
const _N   = 1000  # samples per file (= 1 second at _FS Hz)

const _IFS_DIR    = mktempdir()
const _IFS_FILE1  = joinpath(_IFS_DIR, "T1-C__0__20230101_120000.wav")
const _IFS_FILE2  = joinpath(_IFS_DIR, "T1-C__0__20230101_120002.wav")
const _IFS_ARROW  = joinpath(_IFS_DIR, "index.arrow")

# Write files. Float64 samples round-trip through WAV.wavwrite / wavread exactly.
WAV.wavwrite(fill( 0.5, _N), _IFS_FILE1; Fs = _FS)
WAV.wavwrite(fill(-0.5, _N), _IFS_FILE2; Fs = _FS)

# Build the index once; reuse across all testsets in this file.
const _IFS_TBL = EcoAcoustics.build_index(_IFS_DIR; recorder = "sm3m",
                                           output_path = _IFS_ARROW)

# Convenience: file1 end time and gap end (= file2 start)
const _T1_END = _T1 + Nanosecond(round(Int, 1e9 * _N / _FS))  # 12:00:01
const _T2_END = _T2 + Nanosecond(round(Int, 1e9 * _N / _FS))  # 12:00:03

# ─── IndexedFileSource construction ──────────────────────────────────────────

@testset "IndexedFileSource: from Arrow path" begin
    src = EcoAcoustics.IndexedFileSource(_IFS_ARROW)
    @test src isa EcoAcoustics.IndexedFileSource
    @test src.root == dirname(abspath(_IFS_ARROW))
end

@testset "IndexedFileSource: from Arrow.Table" begin
    src = EcoAcoustics.IndexedFileSource(_IFS_TBL; root = _IFS_DIR)
    @test src isa EcoAcoustics.IndexedFileSource
    @test src.root == _IFS_DIR
end

@testset "IndexedFileSource: missing file throws" begin
    @test_throws ArgumentError EcoAcoustics.IndexedFileSource("/no/such/file.arrow")
end

@testset "IndexedFileSource: empty table throws" begin
    empty_tbl = let
        buf = IOBuffer()
        Arrow.write(buf, (
            file_path           = String[],
            start_time          = DateTime[],
            end_time            = DateTime[],
            fs                  = Float32[],
            nsamples            = Int64[],
            recorder            = String[],
            recorder_id         = Union{String,Missing}[],
            site_id             = Union{String,Missing}[],
            hydrophone_id       = Union{String,Missing}[],
            time_uncertainty_ms = Union{Float64,Missing}[],
            notes               = Union{String,Missing}[],
        ))
        seekstart(buf)
        Arrow.Table(buf)
    end
    @test_throws ArgumentError EcoAcoustics.IndexedFileSource(empty_tbl; root = _IFS_DIR)
end

# ─── time_range ───────────────────────────────────────────────────────────────

@testset "time_range" begin
    src = EcoAcoustics.IndexedFileSource(_IFS_ARROW)
    t_start, t_stop = EcoAcoustics.time_range(src)
    @test t_start == _T1      # earliest start_time across all rows
    @test t_stop  == _T2_END  # latest   end_time   across all rows
end

# ─── coverage_fraction ────────────────────────────────────────────────────────

@testset "coverage_fraction" begin
    src = EcoAcoustics.IndexedFileSource(_IFS_ARROW)

    # Full 3-second window: two 1-second files separated by a 1-second gap → 2/3.
    frac_full = EcoAcoustics.coverage_fraction(src, _T1, _T2_END)
    @test frac_full ≈ 2.0 / 3.0 atol = 1e-10

    # Window that exactly covers one file → 1.0.
    @test EcoAcoustics.coverage_fraction(src, _T1, _T1_END) ≈ 1.0 atol = 1e-10
    @test EcoAcoustics.coverage_fraction(src, _T2, _T2_END) ≈ 1.0 atol = 1e-10

    # Window entirely in the gap between files → 0.0.
    @test EcoAcoustics.coverage_fraction(src, _T1_END, _T2) ≈ 0.0 atol = 1e-10

    # Degenerate window (stop <= start) → 0.0.
    @test EcoAcoustics.coverage_fraction(src, _T1, _T1) == 0.0
end

# ─── read_audio_range ─────────────────────────────────────────────────────────

@testset "read_audio_range: start >= stop throws" begin
    src = EcoAcoustics.IndexedFileSource(_IFS_ARROW)
    @test_throws ArgumentError EcoAcoustics.read_audio_range(src, _T1, _T1)
    @test_throws ArgumentError EcoAcoustics.read_audio_range(src, _T2, _T1)
end

@testset "read_audio_range: invalid gap_handling throws" begin
    src = EcoAcoustics.IndexedFileSource(_IFS_ARROW)
    @test_throws ArgumentError EcoAcoustics.read_audio_range(
        src, _T1, _T2_END; gap_handling = :skip)
end

@testset "read_audio_range: single file — full extent" begin
    src   = EcoAcoustics.IndexedFileSource(_IFS_ARROW)
    chunk = EcoAcoustics.read_audio_range(src, _T1, _T1_END)

    @test chunk isa EcoAcoustics.Audiodata
    @test EcoAcoustics.nsamples(chunk) == _N
    @test chunk.starttime == _T1
    @test all(isapprox.(chunk.sig,  0.5; atol = 1e-10))
    @test chunk.is_calibrated == false
end

@testset "read_audio_range: second file — full extent" begin
    src   = EcoAcoustics.IndexedFileSource(_IFS_ARROW)
    chunk = EcoAcoustics.read_audio_range(src, _T2, _T2_END)

    @test EcoAcoustics.nsamples(chunk) == _N
    @test chunk.starttime == _T2
    @test all(isapprox.(chunk.sig, -0.5; atol = 1e-10))
end

@testset "read_audio_range: span two files with gap" begin
    # Window covers both files and the 1-second gap between them.
    # Expected: [1000 × 0.5, 1000 × 0.0, 1000 × -0.5]
    src   = EcoAcoustics.IndexedFileSource(_IFS_ARROW)
    chunk = EcoAcoustics.read_audio_range(src, _T1, _T2_END)

    n_out = round(Int, (_T2_END - _T1).value / 1000.0 * _FS)  # = 3000
    @test EcoAcoustics.nsamples(chunk) == n_out
    @test chunk.starttime == _T1

    @test all(isapprox.(chunk.sig[1        : _N],      0.5; atol = 1e-10))  # file 1
    @test all(chunk.sig[_N + 1   : 2 * _N]  .== 0.0)                        # gap (zeros)
    @test all(isapprox.(chunk.sig[2*_N + 1 : 3 * _N], -0.5; atol = 1e-10)) # file 2
end

@testset "read_audio_range: window entirely in gap" begin
    # Window falls in the 1-second gap between the two files.
    src   = EcoAcoustics.IndexedFileSource(_IFS_ARROW)
    chunk = EcoAcoustics.read_audio_range(src, _T1_END, _T2)

    @test EcoAcoustics.nsamples(chunk) == 0
end

@testset "read_audio_range: gap_handling=:error raises" begin
    src = EcoAcoustics.IndexedFileSource(_IFS_ARROW)

    # Full 3-second window has only 2/3 coverage → raises.
    @test_throws ArgumentError EcoAcoustics.read_audio_range(
        src, _T1, _T2_END; gap_handling = :error)

    # Window with full coverage → no throw.
    @test_nowarn EcoAcoustics.read_audio_range(
        src, _T1, _T1_END; gap_handling = :error)
end

@testset "read_audio_range: metadata from first contributing file" begin
    src = EcoAcoustics.IndexedFileSource(_IFS_ARROW)

    chunk = EcoAcoustics.read_audio_range(src, _T1, _T2_END)
    @test chunk.metadata.recorder == "sm3m"
    @test chunk.metadata.site_id  == "T1-C"
    @test chunk.is_calibrated     == false
end

@testset "read_audio_range: partial overlap at file head" begin
    # Window starts 0.5 s into file 1 — only the second half of file 1 overlaps.
    src    = EcoAcoustics.IndexedFileSource(_IFS_ARROW)
    offset = Millisecond(500)
    chunk  = EcoAcoustics.read_audio_range(src, _T1 + offset, _T1_END)

    expected_n = round(Int, 0.5 * _FS)  # 500 samples
    @test EcoAcoustics.nsamples(chunk) == expected_n
    @test all(isapprox.(chunk.sig, 0.5; atol = 1e-10))
end

# ─── coverage_fraction: binary-search correctness ────────────────────────────
#
# Verifies that the binary-search implementation of coverage_fraction produces
# identical results to the original O(N) linear scan on a 120-file index.
# Covers all edge cases from the finding: window before all files, after all
# files, in a gap, spanning a file boundary, and mid-deployment.
#
# The index is synthetic (no files on disk) — coverage_fraction reads only
# idx.start_time and idx.end_time, never file paths.

@testset "coverage_fraction: binary search matches linear scan (120-file index)" begin
    T0          = DateTime(2023, 6, 1)
    n_files     = 120
    file_dur_ms = 10_000   # 10 s per file
    gap_dur_ms  =  5_000   #  5 s gap between files
    stride_ms   = file_dur_ms + gap_dur_ms  # 15 s per slot

    starts = [T0 + Millisecond((i - 1) * stride_ms)               for i in 1:n_files]
    ends   = [T0 + Millisecond((i - 1) * stride_ms + file_dur_ms) for i in 1:n_files]

    buf = IOBuffer()
    Arrow.write(buf, (
        file_path           = ["file$i.wav"                       for i in 1:n_files],
        start_time          = starts,
        end_time            = ends,
        fs                  = fill(Float32(1000),  n_files),
        nsamples            = fill(Int64(10_000),  n_files),
        recorder            = fill("sm3m",          n_files),
        recorder_id         = Union{String,Missing}[missing for _ in 1:n_files],
        site_id             = Union{String,Missing}[missing for _ in 1:n_files],
        hydrophone_id       = Union{String,Missing}[missing for _ in 1:n_files],
        time_uncertainty_ms = Union{Float64,Missing}[missing for _ in 1:n_files],
        notes               = Union{String,Missing}[missing for _ in 1:n_files],
    ))
    seekstart(buf)
    src = EcoAcoustics.IndexedFileSource(Arrow.Table(buf); root = tempdir())

    # Reference: original O(N) linear scan — used only inside this testset.
    function _linear_coverage(src, t_start, t_stop)
        t_stop <= t_start && return 0.0
        win_ms = (t_stop - t_start).value
        cov_ms = 0
        idx = src.index
        for i in 1:length(idx.start_time)
            ov_start = max(t_start, idx.start_time[i])
            ov_stop  = min(t_stop,  idx.end_time[i])
            ov_stop > ov_start && (cov_ms += (ov_stop - ov_start).value)
        end
        return cov_ms / win_ms
    end

    cases = [
        # (label,                                      query_start,                       query_stop)
        ("entirely before all files",                  T0 - Millisecond(60_000),          T0 - Millisecond(1)),
        ("entirely after all files",                   ends[end] + Millisecond(1),        ends[end] + Millisecond(60_000)),
        ("window falls in gap between files 60–61",    ends[60],                          starts[61]),
        ("spans file-60/61 boundary",                  ends[60] - Millisecond(2_000),     starts[61] + Millisecond(2_000)),
        ("single file at start of deployment",         starts[1],                         ends[1]),
        ("single file at end of deployment",           starts[end],                       ends[end]),
        ("full deployment span",                       starts[1],                         ends[end]),
        ("partial overlap at leading edge of file 1",  starts[1] - Millisecond(5_000),   ends[1]),
        ("partial overlap at trailing edge of last",   starts[end],                       ends[end] + Millisecond(5_000)),
        ("mid-deployment spanning files 30–90",        starts[30],                        ends[90]),
        ("window spanning exactly one file boundary",  ends[50],                          ends[51]),
        ("degenerate window (stop == start)",          starts[1],                         starts[1]),
    ]

    for (label, q_start, q_stop) in cases
        expected = _linear_coverage(src, q_start, q_stop)
        @test EcoAcoustics.coverage_fraction(src, q_start, q_stop) ≈ expected atol=1e-10
    end
end
