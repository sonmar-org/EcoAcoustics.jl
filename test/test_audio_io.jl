using Test
using EcoAcoustics
using Dates

    # Helper function for audio metadata validation
    function test_audio_metadata(filename::String; expected_fs=48000, expected_time=nothing, lat=nothing, lon=nothing)
        path = joinpath(TEST_DIR, filename)
        @assert isfile(path) "Missing test file: $path"

#        aud = (lat === nothing || lon === nothing) ?
#              read_audio(path, "none") :
#              read_audio(path, "none"; lat=lat, lon=lon)

        aud = EcoAcoustics.read_audio(path)
        
        @test aud.fs == expected_fs
        if expected_time !== nothing
            @test aud.starttime == expected_time
        end
        # if lat !== nothing && lon !== nothing
        #     @test aud.lat !== missing
        #     @test aud.lon !== missing
        #     @test ≈(aud.lat, lat; atol=1e-5)
        #     @test ≈(aud.lon, lon; atol=1e-5)            
        # end
    end

# @testset "Load WAV file and parse timestamp" begin
#     test_audio_metadata("20250101_120000.wav";
#         expected_time=DateTime(2025, 1, 1, 12, 0, 0),
#         lat=0.0, lon=0.0)
# end

#     @testset "Load FLAC file and parse timestamp" begin
#           test_audio_metadata("test_32.30642_-122.61458_20250101_120000.flac"; expected_time=DateTime(2025, 1, 1, 12, 0, 0))
#     end

# @testset "Filename timestamp in middle of filename" begin
#     test_audio_metadata("example_20250101_120000_audio.wav";
#         expected_time=DateTime(2025, 1, 1, 12, 0, 0),
#         lat=0.0, lon=0.0)
# end

#     @testset "Parse location coordinates +lat,lon" begin
#         test_audio_metadata("test_32.30642_-122.61458_20250101_120000.wav";
#             expected_time=DateTime(2025, 1, 1, 12, 0, 0),
#             lat=32.30642, lon=-122.61458)
#     end

#     @testset "Parse location coordinates without +lat,lon" begin
#         test_audio_metadata("test_noLoc_20250101_120000.wav";
#             expected_time=DateTime(2025, 1, 1, 12, 0, 0),
#             lat=0.0, lon=0.0)
#     end

# @testset "child Audiodata preserves key metadata" begin
#     # path = joinpath(testdir, filename)
#     # @assert isfile(path) "Missing test file: $path"
#     a = read_audio(joinpath(testdir,"test_32.30642_-122.61458_20250101_120000.wav"), "DMON2"; site_id="A5M")
#     b = audiodata_from_parent(a; sig=a.sig[1:1_000], starttime=a.starttime, endtime=a.starttime + Dates.Second(1))
#     @test b.fs == a.fs
#     @test b.recorder == a.recorder
#     @test b.cal == a.cal
#     @test b.lat == a.lat && b.lon == a.lon
#     @test b.site_id == a.site_id
# end
