using Test
using EcoAcoustics
using Statistics

# All expected values are derived analytically from first principles.
# No "run it and see" tests — every assertion has a derivation comment.
#
# Key formulas used throughout:
#
#   rfft of length-N rectangular-windowed signal:
#     window_energy = N  (Σ 1² = N)
#     df            = fs / N
#
#   PSD normalization (Merchant 2015 Eq. 1, after single-sided correction):
#     psd[k] = correction_k × |X[k]|² / (fs × N)
#
#   Parseval identity for single-sided rfft:
#     |X[0]|² + 2·Σ|X[1:end-1]|² + |X[end]|² = N · Σ|x[n]|²
#
#   Therefore:  Σ_k psd[k] × df = mean(signal.^2)

# ─── Spectrogram even-nfft assertion (DD-07) ─────────────────────────────────

@testset "spectrogram: odd nfft rejected (DD-07)" begin
    # window_length = round(Int, 0.101 * 10000) = 1010 (even) — fine.
    # Force an odd nfft explicitly to trigger the assertion.
    signal = randn(Float64, 2000)
    @test_throws AssertionError spectrogram(signal; fs = 10000.0,
                                            window_seconds = 0.1,
                                            nfft = 1001)   # 1001 is odd

    # Even nfft passes (no assertion).
    r = spectrogram(signal; fs = 10000.0, window_seconds = 0.1, nfft = 1000)
    @test r.nfft == 1000
end

# ─── PSDResult struct ─────────────────────────────────────────────────────────

@testset "PSDResult: field types and sizes" begin
    N   = 1024
    fs  = 10000.0
    spec = spectrogram(randn(Float64, N * 4);
                       fs = fs, window_seconds = N / fs, window = :rectangular)
    result = compute_psd(spec)

    @test result.psd_linear    isa Matrix{Float64}
    @test size(result.psd_linear) == (spec.nfft ÷ 2 + 1, length(spec.time))
    @test result.freqs         === spec.freqs   # same object, no copy
    @test result.time          === spec.time
    @test result.fs            === spec.fs
    @test result.window_energy === spec.window_energy
    @test result.nfft          === spec.nfft
    @test result.cal           isa NoCalibration
    @test result.is_calibrated == false
end

@testset "PSDResult: is_calibrated flag" begin
    N    = 1024
    fs   = 10000.0
    spec = spectrogram(randn(Float64, N); fs = fs, window_seconds = N / fs,
                       window = :rectangular)

    @test !compute_psd(spec, NoCalibration()).is_calibrated
    @test  compute_psd(spec, ScalarCalibration(-153.0f0)).is_calibrated
    @test  compute_psd(spec,
                       TFCalibration([0.0, fs/2], [-100.0, -100.0],
                                     "test", :test, "test")).is_calibrated
end

# ─── psd_units dispatch ───────────────────────────────────────────────────────
#
# psd_units dispatches on PSDResult.is_calibrated (a Bool), not on cal type.
# This correctly handles the pre-calibrated-signal case: if a signal is
# calibrated in the time domain before spectrogram, the PSD wrapper sets
# cal=NoCalibration at the PSD layer but is_calibrated reflects the audio's
# is_calibrated flag, which would be true.  Dispatching on cal type would
# return :fullscale²_per_Hz in that case; dispatching on is_calibrated is
# semantically correct regardless of how calibration was applied.

@testset "psd_units: dispatch on is_calibrated flag" begin
    N    = 512
    fs   = 8000.0
    spec = spectrogram(randn(Float64, N); fs = fs, window_seconds = N / fs,
                       window = :rectangular)
    tf   = TFCalibration([0.0, fs/2], [-100.0, -100.0], "test", :test, "test")

    # NoCalibration → is_calibrated=false → :fullscale²_per_Hz
    @test psd_units(compute_psd(spec, NoCalibration()))             == :fullscale²_per_Hz
    # ScalarCalibration → is_calibrated=true → :µPa²_per_Hz
    @test psd_units(compute_psd(spec, ScalarCalibration(-153.0f0))) == :µPa²_per_Hz
    # TFCalibration → is_calibrated=true → :µPa²_per_Hz
    @test psd_units(compute_psd(spec, tf))                          == :µPa²_per_Hz
end

# ─── Parseval: sum(psd) × df ≈ mean(signal²) ─────────────────────────────────
#
# Derivation (rectangular window, single frame, nfft = N):
#   window_energy = N
#   single-sided corrected PSD: Σ_k psd[k]
#     = [|X[0]|² + 2Σ|X[1:end-1]|² + |X[end]|²] / (fs × N)
#     = N · Σ|x[n]|² / (fs × N)              [Parseval identity above]
#     = Σ|x[n]|² / fs
#   × df = Σ|x[n]|² / fs × fs/N = mean(x.^2)  ✓

@testset "Parseval: sum(psd) × df ≈ mean(signal²)" begin
    N      = 1024
    fs     = 10000.0
    signal = randn(Float64, N)
    spec   = spectrogram(signal;
                         fs = fs, window_seconds = N / fs,
                         overlap_fraction = 0.0, window = :rectangular)
    result = compute_psd(spec)
    df     = Float64(fs) / spec.nfft

    @test sum(result.psd_linear[:, 1]) * df ≈ mean(signal .^ 2)  rtol = 1e-10
end

# ─── DC signal: all power in bin 0 ───────────────────────────────────────────
#
# signal = ones(N), rectangular window:
#   rfft(ones(N)) = [N, 0, 0, …]
#   psd[1, 1] = |N|² / (fs × N) = N/fs   (DC: no ×2 correction)
#   psd[2:end, 1] = 0

@testset "DC signal: all power in DC bin" begin
    N      = 1024
    fs     = 10000.0
    signal = ones(Float64, N)
    spec   = spectrogram(signal;
                         fs = fs, window_seconds = N / fs,
                         overlap_fraction = 0.0, window = :rectangular)
    result = compute_psd(spec)

    @test result.psd_linear[1, 1] ≈ Float64(N) / fs  rtol = 1e-10
    @test maximum(abs, result.psd_linear[2:end, 1]) < 1e-20
end

# ─── Bin-aligned unit sine: verify single-sided correction and normalisation ──
#
# signal = sin(2π·k₀·n/N), rectangular window, nfft = N:
#   rfft: X[k₀+1] = −j·N/2  (0-indexed k₀, interior bin)
#   |X[k₀+1]|² = (N/2)²
#   single-sided correction: ×2  (interior bin)
#   psd[k₀+1] = 2·(N/2)² / (fs × N) = N/(2·fs)
#
# Consistency check via Parseval:
#   sum × df = N/(2·fs) × (fs/N) = 1/2 = mean(sin²(·)) ✓

@testset "Bin-aligned unit sine: PSD at sine bin, other bins zero" begin
    N   = 1024
    fs  = 10000.0
    k₀  = 20              # 0-indexed interior bin; freq = k₀ * fs/N = 195.3125 Hz
    f₀  = k₀ * fs / N
    signal = sin.(2π * f₀ .* (0:N-1) ./ fs)
    spec   = spectrogram(signal;
                         fs = fs, window_seconds = N / fs,
                         overlap_fraction = 0.0, window = :rectangular)
    result = compute_psd(spec)

    expected = Float64(N) / (2.0 * fs)
    @test result.psd_linear[k₀ + 1, 1] ≈ expected  rtol = 1e-10

    # All other bins are zero (no spectral leakage for bin-aligned sine with
    # rectangular window).  Zero the tested bin before checking.
    rest = copy(result.psd_linear[:, 1])
    rest[k₀ + 1] = 0.0
    @test maximum(abs, rest) < 1e-20
end

# ─── Single-sided correction: DC and Nyquist are NOT doubled ─────────────────
#
# Verify that DC and Nyquist rows are not multiplied by 2.
# Method: use a signal whose rfft has known non-zero values at DC and Nyquist
# and zero elsewhere, then check psd values against the analytic formula
# without the ×2 factor.
#
# signal = cos(π·n) for n=0..N-1: alternating +1/-1.
#   rfft: X[N/2+1] = N  (Nyquist bin), all others = 0.
#   psd[end] = N² / (fs × N) = N/fs  (no ×2)

@testset "Single-sided correction: Nyquist bin not doubled" begin
    N   = 1024
    fs  = 10000.0
    # Cosine at exactly Nyquist: cos(π·n) = [+1, −1, +1, −1, …]
    signal = [(-1.0)^n for n in 0:N-1]
    spec   = spectrogram(signal;
                         fs = fs, window_seconds = N / fs,
                         overlap_fraction = 0.0, window = :rectangular)
    result = compute_psd(spec)

    # Nyquist bin: no ×2, so psd = N² / (fs × N) = N/fs
    @test result.psd_linear[end, 1] ≈ Float64(N) / fs  rtol = 1e-10
    @test maximum(abs, result.psd_linear[1:end-1, 1]) < 1e-20
end

# ─── NoCalibration: psd_linear unchanged, is_calibrated = false ──────────────

@testset "NoCalibration: no-op, is_calibrated false" begin
    N    = 512
    fs   = 8000.0
    spec = spectrogram(randn(Float64, N);
                       fs = fs, window_seconds = N / fs, window = :hann)

    result = compute_psd(spec, NoCalibration())
    @test !result.is_calibrated
    @test result.cal isa NoCalibration

    # Round-trip: NoCalibration PSD equals default (no cal argument)
    result2 = compute_psd(spec)
    @test result.psd_linear == result2.psd_linear
end

# ─── ScalarCalibration: uniform power scale ───────────────────────────────────
#
# ScalarCalibration(S) multiplies the linear PSD by 10^(−S/10).
# For S = −100 dB: factor = 10^(100/10) = 10^10.

@testset "ScalarCalibration: PSD scaled by 10^(-S/10)" begin
    N    = 1024
    fs   = 10000.0
    signal = randn(Float64, N)
    spec   = spectrogram(signal; fs = fs, window_seconds = N / fs,
                         window = :rectangular)

    S      = -100.0f0
    factor = 10.0 ^ (100.0 / 10.0)   # 10^10
    result_raw = compute_psd(spec, NoCalibration())
    result_cal = compute_psd(spec, ScalarCalibration(S))

    @test result_cal.psd_linear ≈ result_raw.psd_linear .* factor  rtol = 1e-10
    @test result_cal.is_calibrated
end

# ─── TFCalibration flat: equivalent to ScalarCalibration ─────────────────────
#
# A flat TFCalibration at tf_dB = S is equivalent to ScalarCalibration(S):
#   TFCalibration divides by tf_lin = 10^(S/10) per bin.
#   ScalarCalibration multiplies by 10^(−S/10) = 1 / 10^(S/10).
# Both are identical for constant S across all frequencies.

@testset "TFCalibration flat: equivalent to ScalarCalibration" begin
    N    = 1024
    fs   = 10000.0
    signal = randn(Float64, N)
    spec   = spectrogram(signal; fs = fs, window_seconds = N / fs,
                         window = :rectangular)

    S      = -100.0
    cal_tf = TFCalibration([0.0, fs / 2.0], [S, S], "test", :test, "test")
    cal_sc = ScalarCalibration(Float32(S))
    result_tf = compute_psd(spec, cal_tf)
    result_sc = compute_psd(spec, cal_sc)

    # Tolerance 1e-10: the two paths differ by at most 1 ULP (~2e-16 relative).
    # x/a and x*(1/a) can differ by 1 ULP even when 1.0/a == b exactly, because
    # IEEE 754 division and multiplication are separate rounding operations.
    @test result_tf.psd_linear ≈ result_sc.psd_linear  rtol = 1e-10
end

# ─── TFCalibration ramp: per-bin scaling ─────────────────────────────────────
#
# A two-point ramp TF: tf_dB = 0 at DC, tf_dB = −20 at Nyquist.
# At exactly DC: tf_lin = 10^(0/10) = 1 → psd unchanged.
# At exactly Nyquist: tf_lin = 10^(−20/10) = 0.01 → psd ÷ 0.01 = psd × 100.
#
# Verify by comparing against manually scaled psd_raw at those two bins.

@testset "TFCalibration ramp: per-bin scaling at DC and Nyquist" begin

    N    = 1024
    fs   = 10000.0
    signal = randn(Float64, N)
    spec   = spectrogram(signal; fs = fs, window_seconds = N / fs,
                         window = :rectangular)

    nyq    = fs / 2.0
    cal_tf = TFCalibration([0.0, nyq], [0.0, -20.0], "test", :test, "test")

    result_raw = compute_psd(spec, NoCalibration())
    result_cal = compute_psd(spec, cal_tf)

    # DC bin: tf_dB = 0 → tf_lin = 1 → no change
    @test result_cal.psd_linear[1, :] ≈ result_raw.psd_linear[1, :]  rtol = 1e-10

    # Nyquist bin: tf_dB = −20 → tf_lin = 0.01 → ×100
    @test result_cal.psd_linear[end, :] ≈ result_raw.psd_linear[end, :] .* 100.0  rtol = 1e-8
end

# ═══════════════════════════════════════════════════════════════════════════════
# D4: Convenience wrappers, average_psd, to_dB
# ═══════════════════════════════════════════════════════════════════════════════

# ─── average_psd ──────────────────────────────────────────────────────────────

@testset "average_psd: single frame equals column 1" begin
    N      = 1024
    fs     = 10000.0
    signal = randn(Float64, N)
    spec   = spectrogram(signal; fs = fs, window_seconds = N / fs,
                         overlap_fraction = 0.0, window = :rectangular)
    result = compute_psd(spec)

    @test size(result.psd_linear, 2) == 1     # one frame
    avg = average_psd(result)
    @test avg isa Vector{Float64}
    @test length(avg) == size(result.psd_linear, 1)
    @test avg ≈ result.psd_linear[:, 1]  rtol = 1e-14
end

@testset "average_psd: multi-frame mean in linear power" begin
    # Four identical frames → average equals each frame exactly.
    N    = 512
    fs   = 8000.0
    # Signal length 4*N with 0% overlap gives exactly 4 frames.
    signal = randn(Float64, N * 4)
    spec   = spectrogram(signal; fs = fs, window_seconds = N / fs,
                         overlap_fraction = 0.0, window = :hann)
    result = compute_psd(spec)

    avg = average_psd(result)
    # mean of rows must equal Statistics.mean along dim=2
    expected = vec(mean(result.psd_linear; dims=2))
    @test avg ≈ expected  rtol = 1e-14
end

@testset "average_psd: Parseval consistency — sum(avg)*df == mean of per-frame Parseval sums" begin
    # Derivation:
    #   avg[k]           = mean over frames j of psd[k, j]
    #   sum(avg) * df    = (sum over k of mean_j(psd[k,j])) * df
    #                    = mean_j(sum_k(psd[k,j]) * df)       [sum and mean commute]
    #                    = mean(per_frame_parseval)
    # Each per_frame_parseval[j] = sum(psd[:,j]) * df ≈ mean(frame_j.^2)
    # (the per-frame Parseval identity for rectangular window, already tested above).
    N  = 512
    fs = 8000.0
    signal = randn(Float64, N * 8)
    spec   = spectrogram(signal; fs = fs, window_seconds = N / fs,
                         overlap_fraction = 0.0, window = :rectangular)
    result = compute_psd(spec)
    df     = Float64(fs) / spec.nfft

    avg = average_psd(result)
    per_frame_parseval = [sum(result.psd_linear[:, j]) * df
                          for j in axes(result.psd_linear, 2)]
    @test sum(avg) * df ≈ mean(per_frame_parseval)  rtol = 1e-12
end

# ─── to_dB ────────────────────────────────────────────────────────────────────

@testset "to_dB: PSDResult — returns Matrix{Float64} of correct shape" begin
    N    = 512
    fs   = 8000.0
    spec = spectrogram(randn(Float64, N); fs = fs, window_seconds = N / fs,
                       window = :rectangular)
    result = compute_psd(spec)
    db     = to_dB(result)

    @test db isa Matrix{Float64}
    @test size(db) == size(result.psd_linear)
    @test db ≈ 10.0 .* log10.(result.psd_linear)  rtol = 1e-12
end

@testset "to_dB: known values (PSDResult and Vector)" begin
    # 10*log10(1.0) = 0.0
    # 10*log10(10.0) = 10.0
    # 10*log10(100.0) = 20.0
    # 10*log10(0.01) = -20.0
    v = [1.0, 10.0, 100.0, 0.01]
    db = to_dB(v)
    @test db isa Vector{Float64}
    @test db[1] ≈  0.0  atol = 1e-12
    @test db[2] ≈ 10.0  atol = 1e-12
    @test db[3] ≈ 20.0  atol = 1e-12
    @test db[4] ≈ -20.0 atol = 1e-12
end

@testset "to_dB: composition with average_psd" begin
    N    = 512
    fs   = 8000.0
    spec = spectrogram(randn(Float64, N * 4); fs = fs, window_seconds = N / fs,
                       overlap_fraction = 0.0, window = :hann)
    result = compute_psd(spec)
    avg    = average_psd(result)
    db_avg = to_dB(avg)

    @test db_avg isa Vector{Float64}
    @test length(db_avg) == length(avg)
    @test db_avg ≈ 10.0 .* log10.(avg)  rtol = 1e-12
end

# ─── compute_psd(Audiodata): calibration resolution ──────────────────────────

@testset "compute_psd(Audiodata): NoCalibration — PSD matches manual route" begin
    N      = 1024
    fs     = 10000.0
    signal = randn(Float64, N)
    audio  = Audiodata(signal, Float32(fs), DateTime(2023, 1, 1))

    result_auto = compute_psd(audio; window_seconds = N / fs, window = :rectangular)

    spec          = spectrogram(signal; fs = Float64(fs), window_seconds = N / fs,
                                window = :rectangular)
    result_manual = compute_psd(spec, NoCalibration())

    @test result_auto.psd_linear ≈ result_manual.psd_linear  rtol = 1e-12
    @test !result_auto.is_calibrated
end

@testset "compute_psd(Audiodata): ScalarCalibration applied at PSD layer" begin
    N      = 1024
    fs     = 10000.0
    signal = randn(Float64, N)
    cal    = ScalarCalibration(-153.0f0)
    audio  = Audiodata(signal, Float32(fs), DateTime(2023, 1, 1);
                       calibration = cal)

    result_auto = compute_psd(audio; window_seconds = N / fs, window = :rectangular)

    spec          = spectrogram(signal; fs = Float64(fs), window_seconds = N / fs,
                                window = :rectangular)
    result_manual = compute_psd(spec, cal)

    @test result_auto.psd_linear ≈ result_manual.psd_linear  rtol = 1e-12
    @test result_auto.is_calibrated
    @test result_auto.cal === cal
end

@testset "compute_psd(Audiodata): Rockhopper auto-calibration via get_profile" begin
    # Rockhopper recorder: lookup_calibration returns NoCalibration (legacy path
    # removed), but compute_psd auto-resolves via get_profile(:rockhopper).
    N      = 4800
    fs     = 48000.0
    signal = randn(Float64, N)
    audio  = Audiodata(signal, Float32(fs), DateTime(2023, 1, 1);
                       recorder = "rockhopper")

    result_auto = compute_psd(audio; window_seconds = N / fs, window = :hann)

    # Manual route: explicit get_profile
    rh_tf         = get_profile(:rockhopper).tf
    spec          = spectrogram(signal; fs = fs, window_seconds = N / fs, window = :hann)
    result_manual = compute_psd(spec, rh_tf)

    @test result_auto.psd_linear ≈ result_manual.psd_linear  rtol = 1e-12
    @test result_auto.is_calibrated
    @test result_auto.cal isa TFCalibration
    @test psd_units(result_auto) == :µPa²_per_Hz
end

@testset "compute_psd(Audiodata): pre-calibrated audio skips PSD-layer cal" begin
    # apply_calibration! produces is_calibrated=true. The PSD wrapper must not
    # apply another calibration on top.
    N      = 1024
    fs     = 10000.0
    signal = randn(Float64, N)
    cal    = ScalarCalibration(-100.0f0)
    audio  = Audiodata(signal, Float32(fs), DateTime(2023, 1, 1);
                       calibration = cal)
    audio_precal = apply_calibration(audio)
    @test audio_precal.is_calibrated

    result = compute_psd(audio_precal; window_seconds = N / fs, window = :rectangular)

    # PSD-layer cal is NoCalibration (signal was pre-calibrated)
    @test result.cal isa NoCalibration
    # is_calibrated must be true: audio.is_calibrated propagates into PSDResult
    # even though resolved_cal is NoCalibration at the PSD layer.
    @test result.is_calibrated
    @test psd_units(result) == :µPa²_per_Hz
    # PSD values match spectrogram of pre-calibrated signal with no PSD-layer cal
    spec     = spectrogram(audio_precal.sig; fs = Float64(fs), window_seconds = N / fs,
                           window = :rectangular)
    expected = compute_psd(spec, NoCalibration())
    @test result.psd_linear ≈ expected.psd_linear  rtol = 1e-12
end

@testset "compute_psd(Audiodata): unknown recorder warns, returns uncalibrated" begin
    N      = 512
    fs     = 8000.0
    audio  = Audiodata(randn(Float64, N), Float32(fs), DateTime(2023, 1, 1);
                       recorder = "unknown")
    # Expect a warning about missing calibration.
    result = @test_logs (:warn, r"no calibration found") compute_psd(
        audio; window_seconds = N / fs, window = :rectangular)
    @test !result.is_calibrated
    @test result.cal isa NoCalibration
end

# ─── compute_psd: calibration cascade priority ────────────────────────────────
#
# _psd_calibration resolves calibration in priority order (DD-14). Each subtest
# triggers exactly one exit and asserts the distinguishing signal for that branch.
# Keeping all four exits together makes a future reordering immediately visible
# as a test failure rather than a silent change in output units.

@testset "compute_psd: calibration cascade priority" begin
    N  = 1024
    fs = 10000.0

    @testset "cascade step 1: pre-calibrated signal → NoCalibration" begin
        # audio.is_calibrated == true fires before any other check.
        # _psd_calibration returns NoCalibration() (no further PSD-layer cal),
        # but the Audiodata wrapper propagates audio.is_calibrated → is_cal = true.
        # Steps 2 and 3 are not reached; pre-condition assertions confirm it.
        cal          = ScalarCalibration(-100.0f0)
        audio        = Audiodata(randn(Float64, N), Float32(fs), DateTime(2023, 1, 1);
                                 calibration = cal)
        audio_precal = apply_calibration(audio)
        @test audio_precal.is_calibrated            # step 1 trigger: pre-condition

        result = compute_psd(audio_precal; window_seconds = N / fs, window = :rectangular)
        @test result.cal isa NoCalibration          # cascade returned NoCalibration at PSD layer
        @test result.is_calibrated                  # propagated from audio.is_calibrated
        @test psd_units(result) == :µPa²_per_Hz    # unit reflects physical state, not cal type
    end

    @testset "cascade step 2: explicit calibration attachment → use as-is" begin
        # audio.calibration isa !NoCalibration fires before get_profile lookup (step 3).
        # Identity check (===) proves result.cal is the exact object from audio.calibration,
        # not a TFCalibration from a profile.
        explicit = ScalarCalibration(-153.0f0)
        audio    = Audiodata(randn(Float64, N), Float32(fs), DateTime(2023, 1, 1);
                             calibration = explicit)
        @test !audio.is_calibrated                  # step 1 not triggered
        @test !(audio.calibration isa NoCalibration) # step 2 trigger: pre-condition

        result = compute_psd(audio; window_seconds = N / fs, window = :rectangular)
        @test result.cal === explicit               # same object: step 2 returned audio.calibration
        @test result.is_calibrated
    end

    @testset "cascade step 3: recognized recorder → shipped TF" begin
        # audio.calibration is NoCalibration and is_calibrated is false, so the
        # cascade reaches get_profile(:rockhopper) and returns its TFCalibration.
        # Steps 1 and 2 are not reached; pre-condition assertions confirm it.
        N_rh  = 4800
        fs_rh = 48000.0
        audio = Audiodata(randn(Float64, N_rh), Float32(fs_rh), DateTime(2023, 1, 1);
                          recorder = "rockhopper")
        @test !audio.is_calibrated                  # step 1 not triggered
        @test audio.calibration isa NoCalibration   # step 2 not triggered; step 3 trigger

        result = compute_psd(audio; window_seconds = N_rh / fs_rh, window = :hann)
        @test result.cal isa TFCalibration          # cascade resolved via get_profile
        @test result.is_calibrated
        @test psd_units(result) == :µPa²_per_Hz
    end

    @testset "cascade step 4: unknown recorder → warn + NoCalibration" begin
        # All three earlier steps fail: is_calibrated is false, calibration is
        # NoCalibration, and the recorder has no registered profile. The cascade
        # emits @warn and returns NoCalibration — the silent-failure exit.
        # @test_logs verifies the warning fires; without it, a future change that
        # silently drops the warning would go undetected.
        audio = Audiodata(randn(Float64, N), Float32(fs), DateTime(2023, 1, 1);
                          recorder = "unknown_recorder_xyz")
        @test !audio.is_calibrated                  # step 1 not triggered
        @test audio.calibration isa NoCalibration   # steps 2 and 3 not triggered

        result = @test_logs (:warn, r"no calibration found") compute_psd(
            audio; window_seconds = N / fs, window = :rectangular)
        @test result.cal isa NoCalibration
        @test !result.is_calibrated                 # distinguishes from step 1
        @test psd_units(result) == :fullscale²_per_Hz
    end
end

# ─── compute_psd(AbstractAudioSource, start, stop) ────────────────────────────

@testset "compute_psd(AbstractAudioSource): delegates to Audiodata wrapper" begin
    # Use a SingleFileSource over the real test WAV (48 kHz, 10 s, recorder=unknown).
    # Verify that compute_psd on the source produces the same result as
    # reading the audio manually and calling the Audiodata wrapper.
    wav_path = joinpath(@__DIR__, "test_files", "test_real.wav")
    src      = SingleFileSource(wav_path; recorder = "unknown")
    t_start, t_stop = time_range(src)
    # Trim to 1 second to keep the test fast.
    t_end = t_start + Dates.Millisecond(1000)

    result_src = @test_logs (:warn, r"no calibration found") compute_psd(
        src, t_start, t_end;
        gap_handling   = :zero_fill,
        window_seconds = 0.1,
        window         = :hann)

    audio          = read_audio_range(src, t_start, t_end; gap_handling = :zero_fill)
    result_manual  = @test_logs (:warn, r"no calibration found") compute_psd(
        audio; window_seconds = 0.1, window = :hann)

    @test result_src isa PSDResult
    @test result_src.psd_linear ≈ result_manual.psd_linear  rtol = 1e-12
    @test result_src.freqs      == result_manual.freqs
end

# ─── compute_psd(AbstractAudioSource) — full-source form ─────────────────────

@testset "compute_psd(src): round-trip equivalence with range form" begin
    # Full-source form must produce identical output to the range form over
    # the same window (time_range(src) start→stop).
    wav_path = joinpath(@__DIR__, "test_files", "test_real.wav")
    src      = SingleFileSource(wav_path; recorder = "unknown")
    t_start, t_stop = time_range(src)

    result_full  = @test_logs (:warn, r"no calibration found") compute_psd(
        src; window_seconds = 0.1, window = :hann)
    result_range = @test_logs (:warn, r"no calibration found") compute_psd(
        src, t_start, t_stop; window_seconds = 0.1, window = :hann)

    @test result_full.psd_linear == result_range.psd_linear   # identical: same code path
    @test result_full.freqs      == result_range.freqs
    @test result_full.nfft       == result_range.nfft
end

@testset "compute_psd(src): keyword forwarding — window_seconds and nfft" begin
    wav_path        = joinpath(@__DIR__, "test_files", "test_real.wav")
    src             = SingleFileSource(wav_path; recorder = "unknown")
    fs              = 48000.0
    t_start, t_stop = time_range(src)

    result = @test_logs (:warn, r"no calibration found") compute_psd(
        src; window_seconds = 2.0, overlap_fraction = 0.5, window = :hann)

    expected_nfft = round(Int, 2.0 * fs)
    @test result.nfft == expected_nfft
    @test size(result.psd_linear, 1) == expected_nfft ÷ 2 + 1

    # Read the audio to get the exact sample count so expected_frames mirrors
    # spectrogram()'s formula without hardcoding a value tied to this file's duration.
    # gap_handling=:zero_fill matches the default used by compute_psd(src; ...).
    n               = nsamples(read_audio_range(src, t_start, t_stop;
                                                gap_handling = :zero_fill))
    window_length   = round(Int, 2.0 * fs)
    hop             = max(1, round(Int, (1.0 - 0.5) * window_length))
    expected_frames = div(n - window_length, hop) + 1
    @test size(result.psd_linear, 2) == expected_frames
end

@testset "compute_psd(src): calibration kwarg forwarded, not dropped" begin
    # Pass an explicit cal; verify it appears in PSDResult.cal and that
    # is_calibrated reflects the provided calibration (not auto-resolution).
    wav_path = joinpath(@__DIR__, "test_files", "test_real.wav")
    src      = SingleFileSource(wav_path; recorder = "unknown")
    explicit = ScalarCalibration(-100.0f0)

    result = compute_psd(src; window_seconds = 0.1, cal = explicit)

    @test result.cal === explicit
    @test result.is_calibrated
end
