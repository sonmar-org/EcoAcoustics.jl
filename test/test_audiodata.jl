using Test
using EcoAcoustics
using Dates

@testset "Audiodata basic construction" begin
    sig = randn(1_000)
    fs  = 48_000
    t0  = DateTime(2025, 1, 1, 0, 0, 0)

    aud = EcoAcoustics.Audiodata(
        sig,
        fs,
        t0;
        timezone    = "UTC",
        lat         = 38.5,
        lon         = -74.5,
        site_id     = "A5M",
        recorder    = "Rockhopper",
        recorder_id = "RH428",
        calibration = EcoAcoustics.NoCalibration(),
    )

    @test aud.sig isa Vector{Float64}
    @test length(aud.sig) == 1_000
    @test aud.fs == 48_000f0
    @test aud.starttime == t0
    @test aud.is_calibrated == false
    @test aud.calibration isa EcoAcoustics.NoCalibration

    @test aud.metadata.timezone    == "UTC"
    @test aud.metadata.lat         == 38.5
    @test aud.metadata.lon         == -74.5
    @test aud.metadata.site_id     == "A5M"
    @test aud.metadata.recorder    == "Rockhopper"
    @test aud.metadata.recorder_id == "RH428"
end

@testset "Audiodata endtime helper" begin
    # 48000 samples at 48000 Hz = exactly 1 second
    aud = EcoAcoustics.Audiodata(ones(48_000), 48_000, DateTime(2025, 1, 1))
    @test EcoAcoustics.endtime(aud) == DateTime(2025, 1, 1) + Second(1)

    @test EcoAcoustics.nsamples(aud) == 48_000
    @test EcoAcoustics.duration(aud) == Second(1)
end

@testset "Audiodata invariants" begin
    sig = randn(1_000)
    t0  = DateTime(2025, 1, 1)

    # is_calibrated defaults to false
    aud = EcoAcoustics.Audiodata(sig, 48_000, t0)
    @test aud.is_calibrated == false

    # fs must be positive
    @test_throws ArgumentError EcoAcoustics.Audiodata(sig, 0, t0)
    @test_throws ArgumentError EcoAcoustics.Audiodata(sig, -1, t0)

    # empty signal is now allowed (used as no-overlap sentinel by read_audio_range)
    @test EcoAcoustics.nsamples(EcoAcoustics.Audiodata(Float64[], 48_000, t0)) == 0
end
