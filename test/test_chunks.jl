using Test
using EcoAcoustics
using Dates
import DataFrames: DataFrame, nrow, names
import WAV

# ─── Shared fixtures ──────────────────────────────────────────────────────────
#
# Two source configurations, each built once and reused across testsets:
#
#   _C_SFS  — SingleFileSource wrapping one 3-second file at 1000 Hz.
#             Signal value = 0.5, start time = 2023-01-02 08:00:00.
#
#   _C_IFS_SRC — IndexedFileSource with two 1-second files and a 1-second gap:
#             file1: 2023-02-01 12:00:00 – 12:00:01, value =  0.5
#             gap:   2023-02-01 12:00:01 – 12:00:02
#             file2: 2023-02-01 12:00:02 – 12:00:03, value = -0.5
#
# SM3M-formatted filenames are required so that parse_filename extracts
# timestamps correctly.  The format is T1-C__0__yyyymmdd_HHMMSS.wav.

const _C_FS = 1000    # Hz
const _C_N  = 1000    # samples per 1-second window at _C_FS Hz

# ── SingleFileSource ─────────────────────────────────────────────────────────

const _C_SFS_DIR  = mktempdir()
const _C_SFS_T    = DateTime(2023, 1, 2, 8, 0, 0)
const _C_SFS_FILE = joinpath(_C_SFS_DIR, "T1-C__0__20230102_080000.wav")

# 3-second file: 3 × _C_N samples, all 0.5
WAV.wavwrite(fill(0.5, 3 * _C_N), _C_SFS_FILE; Fs = _C_FS)

const _C_SFS     = SingleFileSource(_C_SFS_FILE; recorder = "sm3m")
const _C_SFS_END = _C_SFS_T + Millisecond(3000)   # 08:00:03

# ── IndexedFileSource with gap ────────────────────────────────────────────────

const _C_IFS_DIR   = mktempdir()
const _C_IFS_T1    = DateTime(2023, 2, 1, 12, 0, 0)
const _C_IFS_T2    = DateTime(2023, 2, 1, 12, 0, 2)
const _C_IFS_FILE1 = joinpath(_C_IFS_DIR, "T1-C__0__20230201_120000.wav")
const _C_IFS_FILE2 = joinpath(_C_IFS_DIR, "T1-C__0__20230201_120002.wav")
const _C_IFS_ARROW = joinpath(_C_IFS_DIR, "index.arrow")

WAV.wavwrite(fill( 0.5, _C_N), _C_IFS_FILE1; Fs = _C_FS)
WAV.wavwrite(fill(-0.5, _C_N), _C_IFS_FILE2; Fs = _C_FS)

EcoAcoustics.build_index(_C_IFS_DIR; recorder = "sm3m", output_path = _C_IFS_ARROW)

const _C_IFS_SRC   = IndexedFileSource(_C_IFS_ARROW)
const _C_IFS_T1_END = _C_IFS_T1 + Millisecond(_C_N)   # 12:00:01
const _C_IFS_T2_END = _C_IFS_T2 + Millisecond(_C_N)   # 12:00:03

# ─── chunks — bad arguments ──────────────────────────────────────────────────

@testset "chunks: bad gap_handling throws" begin
    @test_throws ArgumentError chunks(_C_SFS; chunk_seconds = 1.0, gap_handling = :bad)
end

@testset "chunks: non-positive chunk_seconds throws" begin
    @test_throws ArgumentError chunks(_C_SFS; chunk_seconds = 0.0)
    @test_throws ArgumentError chunks(_C_SFS; chunk_seconds = -1.0)
end

@testset "chunks: non-positive stride_seconds throws" begin
    @test_throws ArgumentError chunks(_C_SFS; chunk_seconds = 1.0, stride_seconds = 0.0)
end

# ─── chunks — SingleFileSource, non-overlapping ───────────────────────────────

@testset "chunks: SingleFileSource non-overlapping" begin
    ch = collect(chunks(_C_SFS; chunk_seconds = 1.0))

    # Three 1-second windows: [0,1), [1,2), [2,3) — all within the 3-second file.
    @test length(ch) == 3

    @test ch[1].starttime == _C_SFS_T
    @test ch[2].starttime == _C_SFS_T + Millisecond(1000)
    @test ch[3].starttime == _C_SFS_T + Millisecond(2000)

    @test nsamples(ch[1]) == _C_N
    @test nsamples(ch[2]) == _C_N
    @test nsamples(ch[3]) == _C_N

    # Signal value is 0.5 throughout the file.
    @test all(isapprox.(ch[1].sig, 0.5; atol = 1e-10))
    @test all(isapprox.(ch[2].sig, 0.5; atol = 1e-10))
    @test all(isapprox.(ch[3].sig, 0.5; atol = 1e-10))
end

# ─── chunks — overlapping windows (stride < chunk) ───────────────────────────

@testset "chunks: overlapping windows" begin
    # chunk=2s, stride=1s → windows [0,2), [1,3), [2,4).
    # The third window [2,4) extends 1 s past the 3-second file; it is
    # emitted with 1 s of real data and 1 s of zero-fill.
    ch = collect(chunks(_C_SFS; chunk_seconds = 2.0, stride_seconds = 1.0))

    @test length(ch) == 3
    @test nsamples(ch[1]) == 2 * _C_N
    @test nsamples(ch[2]) == 2 * _C_N
    @test nsamples(ch[3]) == 2 * _C_N   # 1000 real + 1000 zeros

    # First half of ch[3] is real audio (0.5); second half is zero-fill.
    @test all(isapprox.(ch[3].sig[1:_C_N],          0.5; atol = 1e-10))
    @test all(ch[3].sig[_C_N + 1 : 2 * _C_N] .== 0.0)
end

# ─── chunks — IndexedFileSource, :skip gap handling ──────────────────────────

@testset "chunks: IndexedFileSource gap, :skip" begin
    # With chunk_seconds=1.0 and :skip:
    #   [12:00:00, 12:00:01] — coverage 1.0, file1 → emit
    #   [12:00:01, 12:00:02] — coverage 0.0, gap   → skip
    #   [12:00:02, 12:00:03] — coverage 1.0, file2 → emit
    ch = collect(chunks(_C_IFS_SRC; chunk_seconds = 1.0, gap_handling = :skip))

    @test length(ch) == 2
    @test ch[1].starttime == _C_IFS_T1
    @test ch[2].starttime == _C_IFS_T2

    @test all(isapprox.(ch[1].sig,  0.5; atol = 1e-10))
    @test all(isapprox.(ch[2].sig, -0.5; atol = 1e-10))
end

# ─── chunks — IndexedFileSource, :zero_fill ───────────────────────────────────

@testset "chunks: IndexedFileSource gap, :zero_fill" begin
    # All three 1-second windows are emitted; the gap window is zero-filled.
    ch = collect(chunks(_C_IFS_SRC; chunk_seconds = 1.0, gap_handling = :zero_fill))

    @test length(ch) == 3

    @test ch[1].starttime == _C_IFS_T1
    @test ch[2].starttime == _C_IFS_T1_END   # gap window starts at file1 end
    @test ch[3].starttime == _C_IFS_T2

    @test all(isapprox.(ch[1].sig,  0.5; atol = 1e-10))
    @test all(ch[2].sig .== 0.0)              # gap window is pure zeros
    @test all(isapprox.(ch[3].sig, -0.5; atol = 1e-10))
end

# ─── chunks — :error raises on gap ───────────────────────────────────────────

@testset "chunks: :error raises on partial coverage" begin
    # chunk_seconds=2.0: first window [12:00:00, 12:00:02] spans file1 and the
    # gap → coverage 0.5 < 1.0 → read_audio_range raises ArgumentError.
    it = chunks(_C_IFS_SRC; chunk_seconds = 2.0, gap_handling = :error)
    @test_throws ArgumentError first(it)
end

# ─── process_chunks — bad arguments ──────────────────────────────────────────

@testset "process_chunks: bad arguments throw" begin
    f = chunk -> (x = 1.0,)
    @test_throws ArgumentError process_chunks(_C_SFS, f; chunk_seconds=1.0,
                                              gap_handling=:bad, progress=false)
    @test_throws ArgumentError process_chunks(_C_SFS, f; chunk_seconds=1.0,
                                              on_error=:bad, progress=false)
    @test_throws ArgumentError process_chunks(_C_SFS, f; chunk_seconds=1.0,
                                              parallel=:bad, progress=false)
    @test_throws ArgumentError process_chunks(_C_SFS, f; chunk_seconds=1.0,
                                              device=:bad, progress=false)
end

# ─── process_chunks — basic schema and row count ─────────────────────────────

@testset "process_chunks: basic schema and row count" begin
    # Sum of squared samples = 1000 × 0.5² = 250. mean_sq = 250 / 1000 = 0.25.
    f  = chunk -> (mean_sq = sum(abs2, chunk.sig) / length(chunk.sig),)
    df = process_chunks(_C_SFS, f;
                        chunk_seconds = 1.0,
                        parallel      = :none,
                        progress      = false)

    @test df isa DataFrame
    @test "start_time"         ∈ names(df)
    @test "coverage_fraction"  ∈ names(df)
    @test "mean_sq"            ∈ names(df)
    @test nrow(df) == 3

    @test df[1, :start_time] == _C_SFS_T
    @test df[2, :start_time] == _C_SFS_T + Millisecond(1000)
    @test df[3, :start_time] == _C_SFS_T + Millisecond(2000)

    @test all(isapprox.(df[!, :coverage_fraction], 1.0; atol = 1e-10))
    @test all(isapprox.(df[!, :mean_sq],           0.25; atol = 1e-10))
end

# ─── process_chunks — gap source with :skip ───────────────────────────────────

@testset "process_chunks: gap source with :skip" begin
    f  = chunk -> (sum_abs = sum(abs, chunk.sig),)
    df = process_chunks(_C_IFS_SRC, f;
                        chunk_seconds = 1.0,
                        gap_handling  = :skip,
                        parallel      = :none,
                        progress      = false)

    # Two chunks emitted (gap window excluded).
    @test nrow(df) == 2
    @test df[1, :coverage_fraction] ≈ 1.0 atol = 1e-10
    @test df[2, :coverage_fraction] ≈ 1.0 atol = 1e-10
end

# ─── process_chunks — on_error=:skip drops failed chunk ──────────────────────

@testset "process_chunks: on_error=:skip drops failed chunk" begin
    attempt = Ref(0)
    f = chunk -> begin
        attempt[] += 1
        attempt[] == 1 && error("forced error on chunk 1")
        return (val = 1.0,)
    end

    df = process_chunks(_C_SFS, f;
                        chunk_seconds = 1.0,
                        parallel      = :none,
                        on_error      = :skip,
                        progress      = false)

    # 3 chunks, first one failed → 2 rows.
    @test nrow(df) == 2
    # Remaining rows have the metric column.
    @test "val" ∈ names(df)
    @test all(isapprox.(df[!, :val], 1.0; atol = 1e-10))
end

# ─── process_chunks — on_error=:record stores error row ──────────────────────

@testset "process_chunks: on_error=:record stores error row" begin
    attempt = Ref(0)
    f = chunk -> begin
        attempt[] += 1
        attempt[] == 1 && error("forced error on chunk 1")
        return (val = 1.0,)
    end

    df = process_chunks(_C_SFS, f;
                        chunk_seconds = 1.0,
                        parallel      = :none,
                        on_error      = :record,
                        progress      = false)

    # All 3 chunks present; first one is an error row.
    @test nrow(df) == 3
    @test "error" ∈ names(df)

    # Error row: has a non-missing error string; val is missing.
    @test !ismissing(df[1, :error])
    @test ismissing(df[1, :val])

    # Success rows: val is set; error is missing.
    @test ismissing(df[2, :error])
    @test !ismissing(df[2, :val])
end

# ─── process_chunks — parallel=:none and :threads give same results ───────────

@testset "process_chunks: parallel=:none and :threads agree" begin
    f = chunk -> (total = sum(chunk.sig),)

    df_seq = process_chunks(_C_SFS, f;
                            chunk_seconds = 1.0,
                            parallel      = :none,
                            progress      = false)
    df_par = process_chunks(_C_SFS, f;
                            chunk_seconds = 1.0,
                            parallel      = :threads,
                            progress      = false)

    @test nrow(df_seq) == nrow(df_par)

    # Sort both by start_time: thread execution order is deterministic for
    # pre-indexed slots, but sorting is cheap insurance.
    sort!(df_seq, :start_time)
    sort!(df_par, :start_time)

    @test df_seq[!, :start_time] == df_par[!, :start_time]
    @test all(isapprox.(df_seq[!, :total], df_par[!, :total]; atol = 1e-10))
end
