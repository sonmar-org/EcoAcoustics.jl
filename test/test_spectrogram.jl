using Test
using EcoAcoustics

# All tests assert on known mathematical properties derived from first principles.
# No "run it and see what comes out" tests — every expected value is computed
# analytically before calling spectrogram.

# ─── 1. DC signal (rectangular window) ───────────────────────────────────────
#
# For signal = ones(N), rectangular window:
#   rfft(ones(N)) = [N, 0, 0, ..., 0]  exactly
#   → DC bin (stft[1]) = N + 0im, all other bins = 0

@testset "DC signal — rectangular window" begin
    N = 256
    signal = ones(Float64, N)
    r = spectrogram(signal;
                    fs               = 1000.0,
                    window_seconds   = N / 1000.0,
                    overlap_fraction = 0.0,
                    window           = :rectangular)

    @test size(r.stft, 2) == 1           # one frame
    @test real(r.stft[1, 1]) ≈ N        atol = 1e-10   # DC = sum(window) * amplitude
    @test imag(r.stft[1, 1]) ≈ 0.0      atol = 1e-10
    @test all(abs.(r.stft[2:end, 1]) .< 1e-10)         # no energy above DC
end

# ─── 2. Bin-aligned unit sine (rectangular window) ───────────────────────────
#
# For x[n] = sin(2π k n / N), rectangular window, nfft = N:
#
#   X[m] = Σ_{n=0}^{N-1} sin(2πkn/N) · exp(-j2πmn/N)
#         = (1/2j) [ Σ exp(j2π(k-m)n/N) - Σ exp(-j2π(k+m)n/N) ]
#
#   For m = k (and 0 < k < N/2, so 2k ≠ 0 mod N):
#     first sum = N (geometric series, exponent = 0)
#     second sum = 0 (exponent = -j4πkn/N, full N rotations ≠ 0)
#   → X[k] = N/(2j) = -j·N/2
#   → abs(X[k]) = N/2, real(X[k]) = 0, imag(X[k]) = -N/2
#
#   The 1-indexed Julia rfft bin for 0-indexed bin k is stft[k+1].
#   All other bins are exactly zero for bin-aligned frequencies with a
#   rectangular window (no spectral leakage).

@testset "Bin-aligned unit sine — rectangular window" begin
    N  = 256
    fs = 1000.0
    k  = 8                          # 0-indexed bin; frequency = k * fs / N = 31.25 Hz
    f0 = k * fs / N                 # exact bin alignment: f0 * N/fs = k (integer)
    signal = sin.(2π * f0 .* (0:N-1) ./ fs)

    r = spectrogram(signal;
                    fs               = fs,
                    window_seconds   = N / fs,
                    overlap_fraction = 0.0,
                    window           = :rectangular)

    @test abs(r.stft[k+1, 1]) ≈ N / 2  atol = 1e-9   # magnitude = N/2
    @test real(r.stft[k+1, 1]) ≈ 0.0   atol = 1e-9   # pure imaginary
    @test imag(r.stft[k+1, 1]) ≈ -N/2  atol = 1e-9   # -j·N/2

    # All other bins should be (near) zero: exact bin alignment with rectangular
    # window means no spectral leakage.
    other_bins = vcat(r.stft[1:k, 1], r.stft[k+2:end, 1])
    @test all(abs.(other_bins) .< 1e-9)
end

# ─── 3. Parseval's theorem (corrected single-sided formula) ──────────────────
#
# For a real signal x of even length N, the rfft returns X of length N/2+1.
# The two-sided DFT has the redundant upper half omitted, so the correct
# Parseval identity for the single-sided rfft spectrum is:
#
#   abs2(X[1]) + 2·Σ abs2(X[2:end-1]) + abs2(X[end]) = N · Σ abs2(x_windowed)
#       DC             interior bins        Nyquist
#
# The factor of 2 restores the energy from the dropped negative-frequency half.
# Without it (i.e., just sum(abs2, X)/N), the result is approximately half the
# time-domain energy — the single-sided correction is not yet applied here.
#
# This test uses a Hann window so the windowed frame ≠ the raw signal, and
# accesses _make_window via module prefix to reconstruct it.

@testset "Parseval's theorem — single-sided rfft (Hann window)" begin
    N      = 512
    signal = randn(Float64, N)
    r = spectrogram(signal;
                    fs               = 1000.0,
                    window_seconds   = N / 1000.0,
                    overlap_fraction = 0.0,
                    window           = :hann)

    X = r.stft[:, 1]

    # Reconstruct the windowed frame independently to form the RHS.
    w   = EcoAcoustics._make_window(:hann, N)
    x_w = w .* signal[1:N]

    parseval_freq = abs2(X[1]) + 2 * sum(abs2, X[2:end-1]) + abs2(X[end])
    parseval_time = N * sum(abs2, x_w)

    @test parseval_freq ≈ parseval_time  rtol = 1e-10
end

@testset "Parseval's theorem — single-sided rfft (rectangular window)" begin
    # Rectangular: windowed frame equals the raw signal, no need for _make_window.
    N      = 512
    signal = randn(Float64, N)
    r = spectrogram(signal;
                    fs               = 1000.0,
                    window_seconds   = N / 1000.0,
                    overlap_fraction = 0.0,
                    window           = :rectangular)

    X = r.stft[:, 1]
    parseval_freq = abs2(X[1]) + 2 * sum(abs2, X[2:end-1]) + abs2(X[end])
    parseval_time = N * sum(abs2, signal[1:N])

    @test parseval_freq ≈ parseval_time  rtol = 1e-10
end

# ─── 4. Frame count formula ───────────────────────────────────────────────────
#
# num_frames = div(length(signal) - window_length, hop) + 1
# Verified for several (N, window_length, hop) combinations.

@testset "Frame count formula" begin
    test_cases = [
        (1000, 256, 128),   # 50% overlap, clean division
        (1000, 256,  64),   # 75% overlap
        ( 512, 100,  50),   # non-power-of-2 window, 50% overlap
        ( 300, 100,  25),   # 75% overlap, short signal
    ]
    for (N, wl, hs) in test_cases
        signal   = randn(Float64, N)
        ol       = 1.0 - hs / wl
        r = spectrogram(signal;
                        fs               = 1000.0,
                        window_seconds   = wl / 1000.0,
                        overlap_fraction = ol,
                        window           = :rectangular)
        expected = div(N - wl, hs) + 1
        # Wrap in a named testset so failures report which (N, wl, hs) case failed.
        @testset "N=$N wl=$wl hs=$hs" begin
            @test size(r.stft, 2) == expected
            @test length(r.time)  == expected
        end
    end
end

# ─── 5. Partial trailing frames are dropped ───────────────────────────────────
#
# N=400, window_length=256, hop=128 (50% overlap):
#   Frame 1: samples 1..256
#   Frame 2: samples 129..384
#   Remaining: samples 385..400 = 16 samples < 256 → dropped
#
# div(400 - 256, 128) + 1 = div(144, 128) + 1 = 1 + 1 = 2 frames

@testset "Partial trailing frame dropped" begin
    N  = 400
    wl = 256
    hs = 128
    signal = randn(Float64, N)
    r = spectrogram(signal;
                    fs               = 1000.0,
                    window_seconds   = wl / 1000.0,
                    overlap_fraction = 0.5,
                    window           = :rectangular)

    @test size(r.stft, 2) == 2

    # Confirm there genuinely are leftover samples (the drop is non-trivial).
    leftover = N - ((size(r.stft, 2) - 1) * hs + wl)
    @test leftover == 16   # 400 - (1*128 + 256) = 16 samples remain, too few
end

# ─── 6. Window choice propagates + energy ordering ───────────────────────────
#
# Window energy ordering (for large N, from lowest to highest):
#   Blackman ≈ 0.305·N  < Hann ≈ 0.375·N  < Hamming ≈ 0.397·N  < Rectangular = N
#
# Derivation:
#   Rectangular: Σ 1² = N
#   Hann:        Σ (0.5 - 0.5·cos)² ≈ N·(0.25 + 0.25·0.5)           = 0.375N
#   Hamming:     Σ (0.54 - 0.46·cos)² ≈ N·(0.54² + 0.46²/2)         ≈ 0.397N
#   Blackman:    Σ (0.42 - 0.5·cos + 0.08·cos2)² ≈ N·(0.42² + 0.5²/2 + 0.08²/2) ≈ 0.305N
# (Cross-terms vanish by orthogonality for large N.)

@testset "Window choice propagates" begin
    N      = 512
    signal = randn(Float64, N)
    kw     = (fs=1000.0, window_seconds=N/1000.0, overlap_fraction=0.0)

    r_hann     = spectrogram(signal; kw..., window=:hann)
    r_hamming  = spectrogram(signal; kw..., window=:hamming)
    r_blackman = spectrogram(signal; kw..., window=:blackman)
    r_rect     = spectrogram(signal; kw..., window=:rectangular)

    # Stored window symbol matches
    @test r_hann.window     == :hann
    @test r_hamming.window  == :hamming
    @test r_blackman.window == :blackman
    @test r_rect.window     == :rectangular

    # Energy ordering: Blackman < Hann < Hamming < Rectangular
    @test r_blackman.window_energy < r_hann.window_energy
    @test r_hann.window_energy     < r_hamming.window_energy
    @test r_hamming.window_energy  < r_rect.window_energy

    # Rectangular is exactly N (sum of N ones squared)
    @test r_rect.window_energy ≈ Float64(N)  atol = 1e-12

    # Different windows produce different STFTs on the same signal
    @test r_hann.stft != r_rect.stft
end

# ─── 7. Plan reuse — identical results to internal plan construction ──────────
#
# make_spectrogram_plan returns a plan for the same nfft that spectrogram would
# build internally. FFTW is deterministic for the same algorithm and input, so
# results must be bit-for-bit identical (not just approximately equal).

@testset "Plan reuse — make_spectrogram_plan" begin
    fs             = 1000.0
    window_seconds = 0.256    # 256 samples at 1 kHz
    signal         = randn(Float64, 1000)

    plan = make_spectrogram_plan(fs, window_seconds)

    r_plan   = spectrogram(signal; fs=fs, window_seconds=window_seconds, fft_plan=plan)
    r_nopla  = spectrogram(signal; fs=fs, window_seconds=window_seconds)

    @test r_plan.stft         == r_nopla.stft
    @test r_plan.time         == r_nopla.time
    @test r_plan.freqs        == r_nopla.freqs
    @test r_plan.window_energy == r_nopla.window_energy
    @test r_plan.nfft         == r_nopla.nfft
    @test r_plan.hop          == r_nopla.hop

    # Plan can be used across multiple calls without mutation
    r_again = spectrogram(signal; fs=fs, window_seconds=window_seconds, fft_plan=plan)
    @test r_again.stft == r_plan.stft

    # Plan size mismatch: plan built for 256 samples, spectrogram expects 512
    wrong_plan = make_spectrogram_plan(fs, window_seconds)    # plan for 256
    @test_throws AssertionError spectrogram(signal;
                                            fs             = fs,
                                            window_seconds = 0.512,   # 512-sample window
                                            fft_plan       = wrong_plan)
end

# ─── 8. Float64 contract ──────────────────────────────────────────────────────
#
# The fallback method rejects non-Float64 inputs with an ArgumentError explaining
# why Float64 is required and how to fix it.

@testset "Float64 contract — non-Float64 input rejected" begin
    kw = (fs=1000.0, window_seconds=0.1)

    @test_throws ArgumentError spectrogram(ones(Float32, 200);  kw...)
    @test_throws ArgumentError spectrogram(ones(Float16, 200);  kw...)
    @test_throws ArgumentError spectrogram(ones(Int32,   200);  kw...)

    # Verify the message is informative (not a generic MethodError)
    err = try
        spectrogram(ones(Float32, 200); kw...)
    catch e
        e
    end
    @test contains(err.msg, "Float64")
    @test contains(err.msg, "Float32")
    @test contains(err.msg, "Float64.(signal)")
end

# ─── 9. Edge cases ────────────────────────────────────────────────────────────

@testset "Edge cases — assertion failures" begin
    # Empty signal
    @test_throws AssertionError spectrogram(Float64[];
                                            fs=1000.0, window_seconds=0.1)

    # Signal shorter than window_length
    @test_throws AssertionError spectrogram(randn(50);
                                            fs=1000.0, window_seconds=0.1)

    # nfft < window_length
    @test_throws AssertionError spectrogram(randn(1000);
                                            fs=1000.0, window_seconds=0.1, nfft=10)

    # overlap_fraction = 1.0 (upper bound excluded)
    @test_throws AssertionError spectrogram(randn(1000);
                                            fs=1000.0, window_seconds=0.1,
                                            overlap_fraction=1.0)

    # overlap_fraction negative
    @test_throws AssertionError spectrogram(randn(1000);
                                            fs=1000.0, window_seconds=0.1,
                                            overlap_fraction=-0.1)
end

# ─── 10. Time and frequency axis structure ─────────────────────────────────────
#
# Numerical checks on the axis vectors for a known configuration.
# fs=1000 Hz, N=256 samples → freq resolution = 1000/256 ≈ 3.90625 Hz/bin.
# Single frame (no overlap): time[1] = (0*hop + N/2) / fs = 128/1000 = 0.128 s.

@testset "Time and frequency axis values" begin
    N  = 256
    fs = 1000.0
    signal = randn(Float64, N)
    r = spectrogram(signal;
                    fs               = fs,
                    window_seconds   = N / fs,
                    overlap_fraction = 0.0,
                    window           = :rectangular)

    # Frequency axis
    @test r.freqs[1]   == 0.0                       # DC
    @test r.freqs[end] ≈ fs / 2                     # Nyquist = 500 Hz
    @test length(r.freqs) == N ÷ 2 + 1              # 129 bins
    @test r.freqs[2] - r.freqs[1] ≈ fs / N  atol=1e-9  # bin spacing = 3.90625 Hz

    # Time axis — PAMGuide convention: center at sample N/2 (not (N-1)/2)
    @test length(r.time) == 1
    @test r.time[1] ≈ (N / 2) / fs  atol = 1e-12   # = 0.128 s

    # Two-frame case: verify both centers.
    # N=256, overlap=0.5 → hop=128. Signal length 384 = 1*128 + 256 gives exactly 2 frames.
    # frame 1 center: (0*128 + 128) / 1000 = 0.128 s
    # frame 2 center: (1*128 + 128) / 1000 = 0.256 s
    signal2 = randn(Float64, N + N ÷ 2)   # 384 samples → 2 frames exactly
    r2 = spectrogram(signal2;
                     fs               = fs,
                     window_seconds   = N / fs,
                     overlap_fraction = 0.5,
                     window           = :rectangular)
    @test length(r2.time) == 2
    @test r2.time[1] ≈ 0.128  atol = 1e-12
    @test r2.time[2] ≈ 0.256  atol = 1e-12
end

# ─── 11. nfft zero-padding ────────────────────────────────────────────────────
#
# Zero-padding (nfft > window_length) increases frequency resolution without
# changing the analysis window. The result has more bins, the time axis is
# unchanged, and the window_energy is unchanged (energy depends only on the
# window, not the FFT size).

@testset "nfft zero-padding" begin
    N    = 256
    fs   = 1000.0
    signal = randn(Float64, N)
    kw   = (fs=fs, window_seconds=N/fs, overlap_fraction=0.0, window=:hann)

    r_no_pad  = spectrogram(signal; kw...)
    r_pad     = spectrogram(signal; kw..., nfft = 2 * N)

    # Zero-padding doubles the frequency bins
    @test size(r_pad.stft, 1) == 2 * N ÷ 2 + 1        # = N + 1 bins
    @test size(r_no_pad.stft, 1) == N ÷ 2 + 1         # = N/2 + 1 bins

    # Number of frames is unchanged (time resolution unchanged)
    @test size(r_pad.stft, 2) == size(r_no_pad.stft, 2)

    # Window energy is unchanged (the window hasn't changed)
    @test r_pad.window_energy == r_no_pad.window_energy

    # Nyquist frequency is still fs/2
    @test r_pad.freqs[end] ≈ fs / 2

    # nfft stored correctly
    @test r_pad.nfft == 2 * N
end
