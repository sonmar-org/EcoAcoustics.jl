using Test
using EcoAcoustics
using Dates
using FFTW

# ─── helpers ──────────────────────────────────────────────────────────────────

# Build a minimal RecordingMetadata-equivalent NamedTuple for lookup_calibration.
const _cal_meta = (
    timestamp   = nothing,
    timezone    = nothing,
    lat         = missing,
    lon         = missing,
    site_id     = nothing,
    recorder_id = nothing,
)

# Write a small synthetic TF CSV to a tempfile; return the path.
# Caller is responsible for cleanup (rm).
function _write_test_tf_csv(rows::Vector{Tuple{Float64,Float64}})
    path = tempname() * ".csv"
    open(path, "w") do io
        for (f, db) in rows
            println(io, "$f,$db")
        end
    end
    return path
end

# Helper: construct a minimal TFCalibration with placeholder provenance fields.
# Used throughout these tests to avoid repeating the source/format/notes boilerplate.
function _test_tf(freqs::Vector{Float64}, tf_dB::Vector{Float64})
    TFCalibration(freqs, tf_dB, "test", :test, "test calibration — no gain conversion")
end

# ─── _load_tf_csv unit tests ──────────────────────────────────────────────────

@testset "_load_tf_csv" begin
    path = _write_test_tf_csv([(1.0, 72.27), (100.0, 50.22), (1000.0, 45.00)])
    try
        freqs, tf_db = EcoAcoustics._load_tf_csv(path)

        @test length(freqs) == 3
        @test length(tf_db) == 3

        # Frequencies preserved as-is
        @test freqs[1] ≈ 1.0f0
        @test freqs[2] ≈ 100.0f0
        @test freqs[3] ≈ 1000.0f0

        # Second column is negated on read (CSV positive → package negative)
        @test tf_db[1] ≈ -72.27f0  atol = 1e-3
        @test tf_db[2] ≈ -50.22f0  atol = 1e-3
        @test tf_db[3] ≈ -45.00f0  atol = 1e-3
    finally
        rm(path; force=true)
    end

    # Empty file
    empty_path = tempname() * ".csv"
    open(empty_path, "w") do io; end
    @test_throws AssertionError EcoAcoustics._load_tf_csv(empty_path)
    rm(empty_path; force=true)

    # Non-ascending frequencies
    bad_path = _write_test_tf_csv([(1000.0, 50.0), (100.0, 60.0)])
    @test_throws AssertionError EcoAcoustics._load_tf_csv(bad_path)
    rm(bad_path; force=true)
end

# ─── NoCalibration ────────────────────────────────────────────────────────────

@testset "NoCalibration" begin
    signal = [0.1, -0.2, 0.5, 0.0, -0.9]
    out    = similar(signal)

    apply_calibration!(out, signal, NoCalibration())

    @test out == signal                   # exact copy, not approximate
    @test out !== signal                  # different buffer (copyto!, not alias)

    # Empty signal is a no-op
    empty_signal = Float64[]
    empty_out    = Float64[]
    apply_calibration!(empty_out, empty_signal, NoCalibration())
    @test isempty(empty_out)

    # Audiodata-level: NoCalibration returns same object, is_calibrated stays false
    a = Audiodata(signal, 48000.0, DateTime(2023, 1, 1);
                  calibration = NoCalibration())
    a2 = apply_calibration(a)
    @test a2 === a                        # same object, not a copy
    @test !a2.is_calibrated
end

# ─── ScalarCalibration ────────────────────────────────────────────────────────

@testset "ScalarCalibration" begin
    # SM3M sensitivity: -153 dB → linear multiplier 10^(153/20) ≈ 4.4668e7
    cal     = ScalarCalibration(-153.0f0)
    signal  = [0.5, -0.5, 1.0, 0.0]
    out     = zeros(Float64, length(signal))
    expected_scale = 10^(153.0 / 20)   # ≈ 4.4668e7

    apply_calibration!(out, signal, cal)

    @test out ≈ signal .* expected_scale  rtol = 1e-6

    # is_calibrated flips via Audiodata wrapper
    a = Audiodata(signal, 48000.0, DateTime(2023, 1, 1);
                  calibration = cal)
    @test !a.is_calibrated
    a_cal = apply_calibration(a)
    @test a_cal.is_calibrated
    @test a_cal.sig ≈ signal .* expected_scale  rtol = 1e-6
    @test !a.is_calibrated                 # original unchanged

    # Double-application throws
    @test_throws ArgumentError apply_calibration(a_cal)

    # Out-of-place wrapper
    sig_cal = apply_calibration(signal, cal)
    @test sig_cal ≈ signal .* expected_scale  rtol = 1e-6
end

@testset "ScalarCalibration GPU" begin
    # Skip if no CUDA-capable GPU is available
    cuda_available = false
    try
        @eval using CUDA
        cuda_available = CUDA.functional()
    catch
    end

    if cuda_available
        cal    = ScalarCalibration(-153.0f0)
        signal = randn(Float64, 1024)
        expected_scale = 10^(153.0 / 20)

        out_cpu = zeros(Float64, 1024)
        apply_calibration!(out_cpu, signal, cal)

        sig_gpu = CuArray(signal)
        out_gpu = CUDA.zeros(Float64, 1024)
        apply_calibration!(out_gpu, sig_gpu, cal)

        @test Array(out_gpu) ≈ out_cpu  rtol = 1e-12
    else
        @test_skip "No CUDA GPU available"
    end
end

# ─── TFCalibration — time domain ─────────────────────────────────────────────

# All flat-TF tests are exact (not approximate): for a spectrally flat TF,
# FFT-multiply-IFFT is equivalent to scalar multiplication by linearity, so
# irfft(tf_mag * rfft(x), N) = tf_mag * x to within floating-point round-trip
# error (~1e-13 relative). No steady-state window or transient trimming needed.

@testset "TFCalibration time domain — flat 0 dB (identity)" begin
    # tf_dB = 0 → tf_mag = 10^(0/20) = 1 → output == input (FFT round-trip precision)
    cal    = _test_tf([1.0, 10000.0], [0.0, 0.0])
    signal = randn(2048)
    out    = similar(signal)
    apply_calibration!(out, signal, cal; fs = 10000.0)
    @test out ≈ signal  rtol = 1e-10
end

@testset "TFCalibration time domain — flat -20 dB (×10 scale)" begin
    # tf_dB = -20 → -tf_dB/20 = 1 → tf_mag = 10 → output == 10 × input
    cal    = _test_tf([1.0, 10000.0], [-20.0, -20.0])
    signal = randn(2048)
    out    = similar(signal)
    apply_calibration!(out, signal, cal; fs = 10000.0)
    @test out ≈ signal .* 10.0  rtol = 1e-10
end

@testset "TFCalibration time domain — short signal (100 samples)" begin
    # Flat TF is exact even for short signals: scalar multiplication is circular-
    # boundary-free. Verifies length is preserved and no out-of-bounds access.
    cal    = _test_tf([1.0, 10000.0], [-20.0, -20.0])
    signal = randn(100)
    out    = similar(signal)
    apply_calibration!(out, signal, cal; fs = 10000.0)
    @test length(out) == 100
    @test out ≈ signal .* 10.0  rtol = 1e-10
end

@testset "TFCalibration time domain — pre-computed plans" begin
    # Pre-computed plans must give the same result as the no-plans path.
    cal    = _test_tf([1.0, 10000.0], [-20.0, -20.0])
    signal = randn(2048)
    buf    = zeros(2048)
    fwd    = FFTW.plan_rfft(buf)
    plans  = (fwd, inv(fwd))
    out_plans   = similar(signal)
    out_noplans = similar(signal)
    apply_calibration!(out_plans,   signal, cal; fs = 10000.0, plans = plans)
    apply_calibration!(out_noplans, signal, cal; fs = 10000.0)
    @test out_plans ≈ out_noplans  rtol = 1e-14
end

# ─── TFCalibration — frequency domain (PSD path) ─────────────────────────────

@testset "TFCalibration frequency domain (apply_calibration_psd!)" begin
    cal   = _test_tf([1.0, 10000.0], [-100.0, -100.0])  # flat -100 dB
    freqs = [0.0, 500.0, 1000.0, 5000.0]
    psd   = [50.0, 60.0, 70.0, 80.0]           # arbitrary dBFS values
    out   = copy(psd)

    apply_calibration_psd!(out, freqs, cal)

    # Each value should increase by 100 (subtracting -100 adds 100)
    @test out ≈ psd .+ 100.0  atol = 1e-4

    # Length mismatch throws
    @test_throws AssertionError apply_calibration_psd!(out, freqs[1:3], cal)
end

# ─── lookup_calibration returns TFCalibration for Rockhopper (legacy path) ───
# NOTE: In deliverable 2, Rockhopper moves to get_profile(:rockhopper). Until
# then, the CALIBRATION_PROFILES entry with tf_path still routes through here.

@testset "lookup_calibration Rockhopper — legacy path removed (now uses get_profile)" begin
    # Rockhopper was removed from CALIBRATION_PROFILES in Deliverable 2.
    # lookup_calibration now returns NoCalibration() with a warning.
    # Use get_profile(:rockhopper).tf for Rockhopper calibration.
    cal = EcoAcoustics.lookup_calibration("dummy.flac", "rockhopper", _cal_meta;
                                          strict = false)
    @test cal isa NoCalibration
end
