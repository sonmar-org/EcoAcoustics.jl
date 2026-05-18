using Test
using EcoAcoustics
using Dates
using WAV
using LibSndFile   # important: this registers FLAC, WAV, AIFF, etc. with FileIO

@testset "read_audio basic WAV (test_real.wav)" begin
    audio_file = joinpath(TEST_DIR, "test_real.wav")

    audio = EcoAcoustics.read_audio(audio_file)
 
    # Input File     : 'test_real.wav'
    # Channels       : 1
    # Sample Rate    : 48000
    # Precision      : 54-bit
    # Duration       : 00:00:10.00 = 480001 samples ~ 750.002 CDDA sectors
    # File Size      : 3.84M
    # Bit Rate       : 3.07M
    # Sample Encoding: 64-bit Floating Point PCM

    @test audio.fs == 48000
    @test length(audio.sig) == 480_001    # 10 s at 48 kHz (+1 sample)
    @test eltype(audio.sig) == Float64  # 64-bit float PCM
    @test ndims(audio.sig) == 1  # 1 Channel
    @test all(isfinite, audio.sig)

    # Sanity checks
    @test maximum(abs.(audio.sig)) > 0
    @test maximum(abs.(audio.sig)) < 10  # loosen this if needed
end

@testset "read_audio basic FLAC (test_real.flac)" begin
    audio_file = joinpath(TEST_DIR, "test_real.flac")

    audio = EcoAcoustics.read_audio(audio_file)
 
    # Input File     : 'test_real.flac'
    # Channels       : 1
    # Sample Rate    : 48000
    # Precision      : 24-bit
    # Duration       : 00:00:10.00 = 480001 samples ~ 750.002 CDDA sectors
    # File Size      : 298k
    # Bit Rate       : 238k
    # Sample Encoding: 24-bit FLAC
    # Comment        : 'Comment=Processed by SoX'

    @test audio.fs == 48000
    @test length(audio.sig) == 480_001    # 10 s at 48 kHz (+1 sample)
    @test eltype(audio.sig) == Float64  # 64-bit float PCM
    @test ndims(audio.sig) == 1  # 1 Channel
    @test all(isfinite, audio.sig)

    # Sanity checks
    @test maximum(abs.(audio.sig)) > 0
    @test maximum(abs.(audio.sig)) < 10  # loosen this if needed
end

#  This is for future AIF support by LibSndFile...maybe we do something
# @testset "read_audio basic AIFF (test_real.aiff)" begin
#     audio_file = joinpath(TEST_DIR, "test_real.aiff")

#     audio = read_audio(audio_file)
 
#     # Input File     : 'test_real.aiff'
#     # Channels       : 1
#     # Sample Rate    : 48000
#     # Precision      : 32-bit
#     # Duration       : 00:00:10.00 = 480001 samples ~ 750.002 CDDA sectors
#     # File Size      : 1.92M
#     # Bit Rate       : 1.54M
#     # Sample Encoding: 32-bit Signed Integer PCM
#     # Comment        : 'Processed by SoX'

#     @test audio.fs == 48000
#     @test length(audio.sig) == 480_001    # 10 s at 48 kHz (+1 sample)
#     @test eltype(audio.sig) == Float64  # 64-bit float PCM
#     @test ndims(audio.sig) == 1  # 1 Channel
#     @test all(isfinite, audio.sig)

#     # Sanity checks
#     @test maximum(abs.(audio.sig)) > 0
#     @test maximum(abs.(audio.sig)) < 10  # loosen this if needed
# end

@testset "read_audio WAV vs FLAC equivalence" begin
    wav_file  = joinpath(TEST_DIR, "test_real.wav")
    flac_file = joinpath(TEST_DIR, "test_real.flac")

    aud_wav  = EcoAcoustics.read_audio(wav_file)
    aud_flac = EcoAcoustics.read_audio(flac_file)

    # Basic consistency
    @test aud_wav.fs == aud_flac.fs
    @test length(aud_wav.sig) == length(aud_flac.sig)

    diff = maximum((aud_wav.sig .- aud_flac.sig))
    #    @info "max diff WAV vs FLAC" diff

    # FLAC decompresses float WAVs via integer PCM, producing ~6e-8 max differences.
    # 1e-7 is a strict but realistic threshold (~–140 dB), treating WAV and FLAC as equivalent.
    @test diff < 1e-7
end



# @testset "read_audio WAV roundtrip" begin
#     fs = 48_000                    # Hz
#     duration_s = 0.01              # 10 ms
#     n_samples = Int(round(fs * duration_s))

#     # Simple test signal: 1 kHz sine, mono
#     t = collect(0:(n_samples-1)) ./ fs
#     sig = sin.(2π * 1000 .* t)

#     # WAV.jl expects samples as N×C matrix (N samples, C channels)
#     data = reshape(sig, :, 1)      # N×1 mono

#     # Write to a temporary WAV file
#     path = joinpath(TEST_DIR, "test_read_audio.wav")
#     WAV.wavwrite(data, fs, path)

#     # Define a start time
#     starttime = DateTime(2025, 1, 1, 0, 0, 0)

#     # Call our function under test
#     aud = read_audio(path; starttime=starttime,
#                      recorder="test_recorder", cal=0.0f0)

#     # Basic checks
#     @test aud.fs == 48_000f0
#     @test length(aud.sig) == n_samples
#     @test aud.starttime == starttime

#     # Check that endtime - starttime is about 10 ms
#     duration = aud.endtime - aud.starttime
#     # duration is a Millisecond-based Period; convert to milliseconds
#     duration_ms = Dates.value(duration)  # integer milliseconds
#     @test duration_ms == round(Int, 1000 * duration_s)
# end
