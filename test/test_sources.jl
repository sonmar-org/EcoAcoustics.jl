using Test
using EcoAcoustics
using Dates

const WAV_FILE = joinpath(TEST_DIR, "test_real.wav")

@testset "SingleFileSource construction" begin
    src = EcoAcoustics.SingleFileSource(WAV_FILE)
    @test src.path == WAV_FILE
    @test src.audio isa EcoAcoustics.Audiodata
    @test EcoAcoustics.nsamples(src.audio) > 0
end

@testset "time_range" begin
    src = EcoAcoustics.SingleFileSource(WAV_FILE)
    t_start, t_stop = EcoAcoustics.time_range(src)
    @test t_start isa DateTime
    @test t_stop  isa DateTime
    @test t_stop > t_start
    # endtime must be consistent with signal length and fs
    expected_end = src.audio.starttime +
                   Nanosecond(round(Int, 1e9 * EcoAcoustics.nsamples(src.audio) /
                                         src.audio.fs))
    @test t_stop == expected_end
end

@testset "coverage_fraction" begin
    src = EcoAcoustics.SingleFileSource(WAV_FILE)
    t_start, t_stop = EcoAcoustics.time_range(src)
    dur = t_stop - t_start

    # Fully inside → 1.0
    @test EcoAcoustics.coverage_fraction(src, t_start, t_stop) ≈ 1.0

    # Half overlap at the head (window starts one duration before file start)
    @test EcoAcoustics.coverage_fraction(src, t_start - dur, t_stop) ≈ 0.5

    # Half overlap at the tail (window ends one duration after file end)
    @test EcoAcoustics.coverage_fraction(src, t_start, t_stop + dur) ≈ 0.5

    # No overlap (window entirely before file)
    @test EcoAcoustics.coverage_fraction(src, t_start - 2*dur, t_start - dur) ≈ 0.0

    # No overlap (window entirely after file)
    @test EcoAcoustics.coverage_fraction(src, t_stop + dur, t_stop + 2*dur) ≈ 0.0
end

@testset "read_audio_range: full file" begin
    src = EcoAcoustics.SingleFileSource(WAV_FILE)
    t_start, t_stop = EcoAcoustics.time_range(src)
    chunk = EcoAcoustics.read_audio_range(src, t_start, t_stop)

    @test chunk isa EcoAcoustics.Audiodata
    @test chunk.starttime == t_start
    @test EcoAcoustics.nsamples(chunk) == EcoAcoustics.nsamples(src.audio)
    @test chunk.sig ≈ src.audio.sig
    @test chunk.is_calibrated == false
end

@testset "read_audio_range: partial overlap — head" begin
    # Window starts one file-duration before the file; only the second half overlaps.
    src = EcoAcoustics.SingleFileSource(WAV_FILE)
    t_start, t_stop = EcoAcoustics.time_range(src)
    dur   = t_stop - t_start
    n_file = EcoAcoustics.nsamples(src.audio)

    chunk = EcoAcoustics.read_audio_range(src, t_start - dur, t_stop)

    # The last n_file samples must exactly match the file signal.
    @test chunk.sig[end - n_file + 1 : end] ≈ src.audio.sig atol=1e-10
    # Everything before that must be zeros.
    @test all(chunk.sig[1 : end - n_file] .== 0.0)
    @test chunk.starttime == t_start - dur
end

@testset "read_audio_range: partial overlap — tail" begin
    # Window ends one file-duration after the file; only the first half overlaps.
    src = EcoAcoustics.SingleFileSource(WAV_FILE)
    t_start, t_stop = EcoAcoustics.time_range(src)
    dur    = t_stop - t_start
    n_file = EcoAcoustics.nsamples(src.audio)

    chunk = EcoAcoustics.read_audio_range(src, t_start, t_stop + dur)

    # The first n_file samples must exactly match the file signal.
    @test chunk.sig[1 : n_file] ≈ src.audio.sig atol=1e-10
    # Everything after that must be zeros.
    @test all(chunk.sig[n_file + 1 : end] .== 0.0)
    @test chunk.starttime == t_start
end

@testset "read_audio_range: no overlap" begin
    src = EcoAcoustics.SingleFileSource(WAV_FILE)
    t_start, t_stop = EcoAcoustics.time_range(src)
    dur = t_stop - t_start

    # Window entirely before the file
    chunk = EcoAcoustics.read_audio_range(src, t_start - 2*dur, t_start - dur)
    @test EcoAcoustics.nsamples(chunk) == 0

    # Window entirely after the file
    chunk = EcoAcoustics.read_audio_range(src, t_stop + dur, t_stop + 2*dur)
    @test EcoAcoustics.nsamples(chunk) == 0
end

@testset "read_audio_range: gap_handling=:error" begin
    src = EcoAcoustics.SingleFileSource(WAV_FILE)
    t_start, t_stop = EcoAcoustics.time_range(src)
    dur = t_stop - t_start

    # Full coverage — should not throw
    @test_nowarn EcoAcoustics.read_audio_range(src, t_start, t_stop;
                                               gap_handling=:error)

    # Partial overlap — should throw
    @test_throws ArgumentError EcoAcoustics.read_audio_range(
        src, t_start - dur, t_stop; gap_handling=:error)

    # No overlap — should throw
    @test_throws ArgumentError EcoAcoustics.read_audio_range(
        src, t_stop + dur, t_stop + 2*dur; gap_handling=:error)
end

@testset "read_audio_range: metadata propagation" begin
    src = EcoAcoustics.SingleFileSource(WAV_FILE; recorder="sm3m",
                                        site_id="T1-C")
    t_start, t_stop = EcoAcoustics.time_range(src)
    chunk = EcoAcoustics.read_audio_range(src, t_start, t_stop)

    @test chunk.metadata.recorder == "sm3m"
    @test chunk.metadata.site_id  == "T1-C"
    @test chunk.calibration       == src.audio.calibration
    @test chunk.is_calibrated     == false
end
