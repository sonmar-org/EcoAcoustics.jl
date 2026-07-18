using Test
using EcoAcoustics
using Dates
using Statistics

# LTSA tests. Every assertion is derived analytically; no "run it and see".
#
# Conventions used throughout:
#   - fs is chosen so freqs land on integers (freqs = 0 : fs/nfft : fs/2).
#   - Bin-aligned cosines (integer cycles per FFT window) + a rectangular
#     window put all energy in one bin with zero leakage, so the peak bin is
#     unambiguous.
#   - Columns are non-overlapping consecutive spans of
#     samples_per_column = round(Int, average_span_seconds * fs) samples.
#     n_columns = div(nsamples, samples_per_column); the trailing partial
#     column is dropped (DD-28).
#   - cal = NoCalibration() is passed explicitly wherever the calibration path
#     is not under test, so no missing-calibration warning is emitted.

# Helper: build an Audiodata from a raw Float64 vector at sampling rate fs.
_audio(sig, fs; kwargs...) = Audiodata(sig, Float32(fs), DateTime(2024, 1, 1);
                                       kwargs...)

# ─── Test 1: time-varying spectral content lands in the right columns ─────────

@testset "LTSA: bands appear in the correct columns" begin
    fs   = 1000.0
    span = 1.0                       # 1-s columns → 1000 samples each
    nfft = 1000                      # window == span → 1 frame per column
    # df = fs / nfft = 1 Hz, so freqs = 0,1,...,500 and freq f is at index f+1.

    # Three columns, each a pure bin-aligned cosine: 100 Hz, 250 Hz, 100 Hz.
    # f * span is an integer number of cycles (100, 250, 100) → bin-aligned.
    freqs_per_col = [100.0, 250.0, 100.0]
    sig = Float64[]
    for f in freqs_per_col
        # Local index i = 1..1000; (i-1)/fs is the time of each sample.
        append!(sig, [cos(2π * f * (i - 1) / fs) for i in 1:1000])
    end

    a  = _audio(sig, fs)
    lt = compute_ltsa(a; average_span_seconds = span, fft_window_seconds = span,
                      window = :rectangular, cal = NoCalibration())

    @test size(lt.matrix, 2) == 3                    # three columns
    @test lt.freqs[end] == fs / 2                    # Nyquist at 500 Hz

    # The peak bin of each column is at that column's frequency.
    for (c, f) in enumerate(freqs_per_col)
        peak_freq = lt.freqs[argmax(lt.matrix[:, c])]
        @test peak_freq == f
    end

    # Cross-check: in column 1 (100 Hz), the 250 Hz bin is essentially zero
    # relative to the 100 Hz bin — no energy leaked into the other band.
    i_100 = 101    # freqs[101] = 100 Hz
    i_250 = 251    # freqs[251] = 250 Hz
    @test lt.matrix[i_250, 1] < lt.matrix[i_100, 1] * 1e-6
    # And symmetrically in column 2 (250 Hz): the 100 Hz bin is near zero.
    @test lt.matrix[i_100, 2] < lt.matrix[i_250, 2] * 1e-6
end

# ─── Test 2: each column equals compute_psd + average_psd over its span ───────

@testset "LTSA: column = energetic mean of that span (multi-frame)" begin
    fs      = 1000.0
    span    = 2.0                    # 2000 samples per column
    fftwin  = 0.5                    # 500-sample FFT window
    overlap = 0.5                    # hop 250 → 7 frames per column: averaging matters
    sig     = randn(Float64, 4000)  # two full columns
    a       = _audio(sig, fs)

    lt = compute_ltsa(a; average_span_seconds = span, fft_window_seconds = fftwin,
                      fft_overlap = overlap, window = :hann, cal = NoCalibration())
    @test size(lt.matrix, 2) == 2

    # Rebuild each column independently via the primitives and compare. Because
    # compute_ltsa uses exactly the same spectrogram → compute_psd → average_psd
    # chain on the same samples, the result is bit-identical (diff == 0).
    for c in 1:2
        seg   = _audio(sig[(c-1)*2000 + 1 : c*2000], fs)
        direct = average_psd(compute_psd(seg; window_seconds = fftwin,
                                         overlap_fraction = overlap,
                                         window = :hann, cal = NoCalibration()))
        @test lt.matrix[:, c] == direct
    end
end

# ─── Test 3: pre-built fft_plan is reused and gives identical output ──────────

@testset "LTSA: fft_plan reuse" begin
    fs     = 1000.0
    span   = 2.0
    fftwin = 0.5
    sig    = randn(Float64, 4000)      # two columns
    a      = _audio(sig, fs)

    # A plan built for the inner FFT size (window_seconds = fftwin) is threaded
    # through to every column. Passing it must not change the result at all —
    # the plan only selects FFTW's algorithm, not the numbers it computes.
    plan   = make_spectrogram_plan(fs, fftwin)
    lt_noplan = compute_ltsa(a; average_span_seconds = span,
                             fft_window_seconds = fftwin, window = :hann,
                             cal = NoCalibration())
    lt_plan   = compute_ltsa(a; average_span_seconds = span,
                             fft_window_seconds = fftwin, window = :hann,
                             fft_plan = plan, cal = NoCalibration())
    @test lt_plan.matrix == lt_noplan.matrix
    @test lt_plan.freqs  == lt_noplan.freqs

    # A plan built for a different FFT size is caught by spectrogram's size
    # assertion. Here the plan is for a 1 s window (1000 samples) but the columns
    # use a 0.5 s window (500 samples) → mismatch on the first column.
    wrong_plan = make_spectrogram_plan(fs, 1.0)
    @test_throws AssertionError compute_ltsa(a; average_span_seconds = span,
                                             fft_window_seconds = fftwin,
                                             fft_plan = wrong_plan,
                                             cal = NoCalibration())
end

# ─── Test 4: calibration propagation (flag, units, scaling) ───────────────────

@testset "LTSA: calibration propagation" begin
    fs  = 1000.0
    sig = randn(Float64, 3000)
    a   = _audio(sig, fs)

    # Uncalibrated (explicit NoCalibration): full-scale units, flag false.
    lt_un = compute_ltsa(a; average_span_seconds = 1.0, cal = NoCalibration())
    @test lt_un.is_calibrated == false
    @test ltsa_units(lt_un)   == :fullscale²_per_Hz
    @test lt_un.cal           isa NoCalibration

    # ScalarCalibration: physical units, flag true. A scalar cal multiplies the
    # PSD by 10^(−S/10) (Merchant/CLAUDE.md convention), so the calibrated
    # matrix is the uncalibrated one scaled by that constant factor.
    S      = -153.0f0
    lt_sc  = compute_ltsa(a; average_span_seconds = 1.0,
                          cal = ScalarCalibration(S))
    @test lt_sc.is_calibrated == true
    @test ltsa_units(lt_sc)   == :µPa²_per_Hz
    @test lt_sc.cal           isa ScalarCalibration
    factor = 10.0^(-Float64(S) / 10)
    @test lt_sc.matrix ≈ lt_un.matrix .* factor rtol = 1e-9

    # Pre-calibrated signal (is_calibrated=true, no PSD-layer cal): the matrix
    # equals the uncalibrated one (NoCalibration applied), but the flag is true
    # because the signal is already in physical units. cal = nothing lets the
    # cascade run; audio.is_calibrated short-circuits it without a warning.
    a_pre  = _audio(sig, fs; is_calibrated = true)
    lt_pre = compute_ltsa(a_pre; average_span_seconds = 1.0)
    @test lt_pre.is_calibrated == true
    @test ltsa_units(lt_pre)   == :µPa²_per_Hz
    @test lt_pre.matrix        == lt_un.matrix          # NoCalibration → no scaling
end

# ─── Test 4: column boundaries and column_times ───────────────────────────────

@testset "LTSA: column count, partial-column drop, column_times" begin
    fs = 1000.0

    # Exactly divisible: 5000 samples / (1 s = 1000 samples) → 5 columns.
    a5 = _audio(randn(Float64, 5000), fs)
    lt5 = compute_ltsa(a5; average_span_seconds = 1.0, cal = NoCalibration())
    @test size(lt5.matrix, 2) == 5
    @test lt5.column_times == [0.0, 1.0, 2.0, 3.0, 4.0]
    @test length(lt5.column_times) == size(lt5.matrix, 2)

    # Partial trailing column dropped: 5500 samples / 1000 → 5 columns (500
    # leftover samples are discarded), not 6.
    a55 = _audio(randn(Float64, 5500), fs)
    lt55 = compute_ltsa(a55; average_span_seconds = 1.0, cal = NoCalibration())
    @test size(lt55.matrix, 2) == 5

    # Fractional span: 0.5 s → samples_per_column = round(Int, 0.5*1000) = 500.
    # 3000 samples / 500 → 6 columns; times step by 500/fs = 0.5 s.
    # fft_window_seconds = 0.25 (250 samples) so a 500-sample column holds it.
    a3 = _audio(randn(Float64, 3000), fs)
    lt3 = compute_ltsa(a3; average_span_seconds = 0.5, fft_window_seconds = 0.25,
                       cal = NoCalibration())
    @test size(lt3.matrix, 2) == 6
    @test lt3.column_times ≈ [0.0, 0.5, 1.0, 1.5, 2.0, 2.5]

    # Field types and axis consistency.
    @test lt5.matrix       isa Matrix{Float64}
    @test lt5.column_times isa Vector{Float64}
    @test lt5.fs           === Float32(fs)
    @test lt5.average_span_seconds === 1.0
    @test length(lt5.freqs) == size(lt5.matrix, 1)      # freq axis matches rows
end

# ─── Test 5: edge cases ───────────────────────────────────────────────────────

@testset "LTSA: edge cases assert loudly" begin
    fs = 1000.0

    # Audio shorter than one column span → nothing to average.
    short = _audio(randn(Float64, 500), fs)          # 0.5 s
    @test_throws AssertionError compute_ltsa(short; average_span_seconds = 1.0,
                                             cal = NoCalibration())

    # average_span_seconds larger than the whole recording → same guard.
    one_sec = _audio(randn(Float64, 1000), fs)
    @test_throws AssertionError compute_ltsa(one_sec; average_span_seconds = 2.0,
                                             cal = NoCalibration())

    # Zero-length audio.
    empty = _audio(Float64[], fs)
    @test_throws AssertionError compute_ltsa(empty; average_span_seconds = 1.0,
                                             cal = NoCalibration())

    # Non-positive spans.
    a = _audio(randn(Float64, 3000), fs)
    @test_throws AssertionError compute_ltsa(a; average_span_seconds = 0.0,
                                             cal = NoCalibration())
    @test_throws AssertionError compute_ltsa(a; average_span_seconds = -1.0,
                                             cal = NoCalibration())
    @test_throws AssertionError compute_ltsa(a; average_span_seconds = 1.0,
                                             fft_window_seconds = 0.0,
                                             cal = NoCalibration())

    # Column narrower than the FFT window: span 0.1 s (100 samples) cannot hold
    # a 1 s (1000-sample) FFT window.
    @test_throws AssertionError compute_ltsa(a; average_span_seconds = 0.1,
                                             fft_window_seconds = 1.0,
                                             cal = NoCalibration())
end
