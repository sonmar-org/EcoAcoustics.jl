using Test
using EcoAcoustics
using Dates
using DataFrames

# Per-cycle band-metrics table tests. band_metrics/band_table are composition
# over compute_psd + compute_spl, so the assertions check (a) the flattening is
# faithful to compute_spl, (b) the table schema/typing is correct, and
# (c) partial-failure handling behaves.

const _FS = 1000.0

# One calibrated 5 s recording of a 100 Hz tone at the given amplitude.
_rec(amp; cal = ScalarCalibration(-150.0f0)) =
    Audiodata([amp * cos(2π * 100 * (i - 1) / _FS) for i in 1:5000],
              Float32(_FS), DateTime(2024, 1, 1); calibration = cal)

const _BANDS = Dict(:tol_100 => (89.1, 112.0), :b200_300 => (200.0, 300.0))

# ─── band_metrics: schema and faithfulness to compute_spl ─────────────────────

@testset "band_metrics: row schema and values" begin
    audio = _rec(1.0)
    row   = band_metrics(audio; bands = _BANDS, window_seconds = 1.0)

    # Meta columns.
    @test row[:start_time] == audio.starttime
    @test row[:duration_s] == 5.0
    # 5000 samples, 1 s window (1000), 50% overlap (hop 500): (5000-1000)/500+1 = 9.
    @test row[:n_frames] == 9

    # Every band × every metric key is present (2 bands × 11 metrics + 3 meta).
    @test length(row) == 2 * length(EcoAcoustics.METRIC_ORDER) + 3

    # Values match compute_spl on the same PSD, and :max is the per-frame maximum.
    psd = compute_psd(audio; window_seconds = 1.0)
    spl = compute_spl(psd; bands = _BANDS)
    for (k, b) in spl.bands
        @test row[Symbol(k, :_mean_dB)]   == b.mean_dB
        @test row[Symbol(k, :_median_dB)] == b.median_dB
        @test row[Symbol(k, :_L1_dB)]     == b.L1_dB
        @test row[Symbol(k, :_L99_dB)]    == b.L99_dB
        @test row[Symbol(k, :_max_dB)]    == maximum(b.spl_dB)
    end

    # n_frames equals the number of PSD frames the percentiles are taken over.
    @test row[:n_frames] == length(psd.time)
end

@testset "band_metrics: empty bands rejected" begin
    @test_throws AssertionError band_metrics(_rec(1.0);
        bands = Dict{Symbol,Tuple{Float64,Float64}}(), window_seconds = 1.0)
end

@testset "band_metrics: uncalibrated audio rejected" begin
    uncal = Audiodata([cos(2π * 100 * (i - 1) / _FS) for i in 1:5000],
                      Float32(_FS), DateTime(2024, 1, 1))   # NoCalibration, unknown recorder
    @test_throws AssertionError band_metrics(uncal; bands = _BANDS,
                                             window_seconds = 1.0,
                                             cal = NoCalibration())
end

# ─── band_table: assembly, ordering, typing ───────────────────────────────────

@testset "band_table: one row per item, typed columns, order" begin
    items = [0.5, 1.0, 2.0, 4.0]
    df = band_table(items; reader = _rec, bands = _BANDS, window_seconds = 1.0)

    @test size(df) == (4, 2 * length(EcoAcoustics.METRIC_ORDER) + 3)   # 4 × 25

    # Meta columns come first, in fixed order, then band metrics.
    @test names(df)[1:3] == ["start_time", "duration_s", "n_frames"]

    # Columns are concretely typed (not Any).
    @test eltype(df.start_time)       === DateTime
    @test eltype(df.duration_s)       === Float64
    @test eltype(df.n_frames)         === Int
    @test eltype(df.tol_100_mean_dB)  === Float64

    # Physical check: doubling amplitude raises the band level by 20·log10(2).
    step = 20 * log10(2)
    @test all(diff(df.tol_100_mean_dB) .≈ step)
end

@testset "band_table: column set derives from bands only (stable schema)" begin
    # Even if some items fail, the schema is fixed by `bands`.
    df = band_table([1.0]; reader = _rec, bands = _BANDS, window_seconds = 1.0)
    expected = vcat([:start_time, :duration_s, :n_frames],
                    EcoAcoustics._band_metric_colnames(_BANDS))
    @test Symbol.(names(df)) == expected
end

# ─── band_table: partial-failure handling ─────────────────────────────────────

@testset "band_table: on_error skip / fail" begin
    reader = x -> x === :bad ? error("corrupt") : _rec(Float64(x))

    # :skip (default) drops the bad item and keeps the good ones.
    df = band_table([1, :bad, 2]; reader = reader, bands = _BANDS,
                    window_seconds = 1.0)
    @test size(df, 1) == 2

    # :fail re-throws.
    @test_throws ErrorException band_table([:bad]; reader = reader,
                                           bands = _BANDS, window_seconds = 1.0,
                                           on_error = :fail)

    # invalid on_error value.
    @test_throws ArgumentError band_table([1.0]; reader = _rec, bands = _BANDS,
                                          window_seconds = 1.0, on_error = :nope)
end

@testset "band_table: parallel=:threads matches serial and validates arg" begin
    items = [0.5, 1.0, 2.0, 4.0, 0.25, 3.0]
    ser = band_table(items; reader = _rec, bands = _BANDS, window_seconds = 1.0,
                     parallel = :none)
    par = band_table(items; reader = _rec, bands = _BANDS, window_seconds = 1.0,
                     parallel = :threads)
    # Rows are assembled by item index, so threaded output is row-for-row identical
    # (deterministic order) regardless of thread count.
    @test names(par) == names(ser)
    @test par == ser

    # A pre-built shared plan gives the same result.
    plan = make_spectrogram_plan(_FS, 1.0)
    par2 = band_table(items; reader = _rec, bands = _BANDS, window_seconds = 1.0,
                      parallel = :threads, fft_plan = plan)
    @test par2 == ser

    # Bad parallel value.
    @test_throws ArgumentError band_table(items; reader = _rec, bands = _BANDS,
                                          window_seconds = 1.0, parallel = :nope)

    # Skips still work under threads (index-aligned; bad item dropped, order kept).
    reader = x -> x === :bad ? error("corrupt") : _rec(Float64(x))
    dfp = band_table([1, :bad, 2, 3]; reader = reader, bands = _BANDS,
                     window_seconds = 1.0, parallel = :threads)
    @test size(dfp, 1) == 3
end

@testset "band_table: empty input gives typed 0-row table" begin
    df = band_table(Audiodata[]; bands = _BANDS, window_seconds = 1.0)
    @test size(df, 1) == 0
    @test size(df, 2) == 2 * length(EcoAcoustics.METRIC_ORDER) + 3
    @test eltype(df.start_time)      === DateTime
    @test eltype(df.tol_100_mean_dB) === Float64
end
