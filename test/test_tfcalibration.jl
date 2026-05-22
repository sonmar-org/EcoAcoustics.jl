using Test
using EcoAcoustics

# Tests for the TFCalibration struct, inner constructor validation, and the
# two interpolation helpers (_interp_tf for dB-domain, _interp_linear_power
# for linear-power domain per DD-11).

# ─── Helper ───────────────────────────────────────────────────────────────────

function _tf(freqs, tf_dB)
    TFCalibration(Float64.(freqs), Float64.(tf_dB),
                  "test", :test, "test calibration")
end

# ─── Constructor validation ───────────────────────────────────────────────────

@testset "TFCalibration constructor — length < 2 rejected" begin
    @test_throws ArgumentError TFCalibration(
        [100.0], [-10.0], "test", :test, "")
end

@testset "TFCalibration constructor — length mismatch rejected" begin
    @test_throws ArgumentError TFCalibration(
        [10.0, 100.0, 1000.0], [-10.0, -12.0],   # 3 freqs vs 2 tf_dB
        "test", :test, "")
end

@testset "TFCalibration constructor — non-ascending frequencies rejected" begin
    # Descending order
    @test_throws ArgumentError TFCalibration(
        [1000.0, 100.0], [-10.0, -12.0], "test", :test, "")
    # Flat (duplicate)
    @test_throws ArgumentError TFCalibration(
        [100.0, 100.0], [-10.0, -12.0], "test", :test, "")
    # First two ascending, then backtrack
    @test_throws ArgumentError TFCalibration(
        [100.0, 200.0, 150.0], [-10.0, -11.0, -12.0], "test", :test, "")
end

@testset "TFCalibration constructor — valid input accepted" begin
    cal = _tf([1.0, 10.0, 100.0, 1000.0], [-229.0, -228.0, -226.0, -224.0])
    @test cal.frequency == [1.0, 10.0, 100.0, 1000.0]
    @test cal.tf_dB     == [-229.0, -228.0, -226.0, -224.0]
    @test cal.source    == "test"
    @test cal.format    === :test
    @test cal.conversion_notes == "test calibration"
end

# ─── _interp_tf (dB-domain interpolation) ────────────────────────────────────

@testset "_interp_tf — flat TF returns uniform value everywhere" begin
    # A flat TF of -50 dB should return -50 at any query frequency.
    cal = _tf([0.0, 1000.0], [-50.0, -50.0])
    @test EcoAcoustics._interp_tf(cal.frequency, cal.tf_dB, 0.0)    ≈ -50.0
    @test EcoAcoustics._interp_tf(cal.frequency, cal.tf_dB, 250.0)  ≈ -50.0
    @test EcoAcoustics._interp_tf(cal.frequency, cal.tf_dB, 500.0)  ≈ -50.0
    @test EcoAcoustics._interp_tf(cal.frequency, cal.tf_dB, 1000.0) ≈ -50.0
end

@testset "_interp_tf — linear ramp in dB (interior and edge)" begin
    # TF: -10 dB at 0 Hz, -20 dB at 100 Hz. Linear interpolation gives -15 at 50 Hz.
    cal = _tf([0.0, 100.0], [-10.0, -20.0])
    @test EcoAcoustics._interp_tf(cal.frequency, cal.tf_dB, 0.0)   ≈ -10.0  atol = 1e-10
    @test EcoAcoustics._interp_tf(cal.frequency, cal.tf_dB, 50.0)  ≈ -15.0  atol = 1e-10
    @test EcoAcoustics._interp_tf(cal.frequency, cal.tf_dB, 100.0) ≈ -20.0  atol = 1e-10
end

@testset "_interp_tf — clamping at boundaries" begin
    cal = _tf([10.0, 1000.0], [-30.0, -50.0])
    # Below lower bound → returns value at lower bound
    @test EcoAcoustics._interp_tf(cal.frequency, cal.tf_dB, 0.0)    ≈ -30.0  atol = 1e-10
    @test EcoAcoustics._interp_tf(cal.frequency, cal.tf_dB, 5.0)    ≈ -30.0  atol = 1e-10
    # Above upper bound → returns value at upper bound
    @test EcoAcoustics._interp_tf(cal.frequency, cal.tf_dB, 2000.0) ≈ -50.0  atol = 1e-10
end

# ─── _interp_linear_power (DD-11: linear-power-scale interpolation) ───────────
#
# DD-11 interpolation procedure (Raven Workbench convention):
#   1. Convert TF dB to linear power: tf_lin = 10^(tf_dB/10)
#   2. Linearly interpolate tf_lin at the target frequency
#   3. The caller divides the PSD by tf_lin (equivalently: psd *= 10^(-tf_dB/10))
#
# These tests verify step 2 by computing expected values by hand.

@testset "_interp_linear_power — flat TF returns same value at all freqs" begin
    cal    = _tf([0.0, 1000.0], [-10.0, -10.0])
    tf_lin = 10 .^ (cal.tf_dB ./ 10)    # [0.1, 0.1]
    @test EcoAcoustics._interp_linear_power(cal.frequency, tf_lin, 0.0)    ≈ 0.1  atol = 1e-12
    @test EcoAcoustics._interp_linear_power(cal.frequency, tf_lin, 500.0)  ≈ 0.1  atol = 1e-12
    @test EcoAcoustics._interp_linear_power(cal.frequency, tf_lin, 1000.0) ≈ 0.1  atol = 1e-12
end

@testset "_interp_linear_power — ramp: midpoint hand-computed" begin
    # TF: -10 dB at 0 Hz, -20 dB at 100 Hz.
    # Linear power values: 10^(-10/10) = 0.1, 10^(-20/10) = 0.01.
    # Midpoint (50 Hz) linear interpolation: (0.1 + 0.01) / 2 = 0.055.
    # Compare: dB-linear midpoint would be (-10 + -20)/2 = -15 dB = 10^(-15/10) ≈ 0.03162.
    # The linear-power result (0.055) differs from the dB-linear result (0.03162),
    # confirming these are distinct interpolation strategies.
    cal    = _tf([0.0, 100.0], [-10.0, -20.0])
    tf_lin = 10 .^ (cal.tf_dB ./ 10)    # [0.1, 0.01]

    result_50 = EcoAcoustics._interp_linear_power(cal.frequency, tf_lin, 50.0)
    @test result_50 ≈ 0.055  atol = 1e-12   # hand-computed: (0.1 + 0.01) / 2

    # Endpoints must be exact:
    @test EcoAcoustics._interp_linear_power(cal.frequency, tf_lin, 0.0)   ≈ 0.1   atol = 1e-12
    @test EcoAcoustics._interp_linear_power(cal.frequency, tf_lin, 100.0) ≈ 0.01  atol = 1e-12
end

@testset "_interp_linear_power — clamping outside range" begin
    cal    = _tf([10.0, 1000.0], [-10.0, -20.0])
    tf_lin = 10 .^ (cal.tf_dB ./ 10)    # [0.1, 0.01]
    # Below lower bound → value at lower bound
    @test EcoAcoustics._interp_linear_power(cal.frequency, tf_lin, 0.0)    ≈ 0.1   atol = 1e-12
    # Above upper bound → value at upper bound
    @test EcoAcoustics._interp_linear_power(cal.frequency, tf_lin, 2000.0) ≈ 0.01  atol = 1e-12
end

# ─── apply_calibration! — PSD matrix (linear power domain) ───────────────────

@testset "apply_calibration! PSD matrix — NoCalibration is no-op" begin
    psd   = [1.0 2.0; 3.0 4.0; 5.0 6.0]  # 3 freq bins × 2 frames
    freqs = [0.0, 100.0, 200.0]
    orig  = copy(psd)
    apply_calibration!(psd, freqs, NoCalibration())
    @test psd == orig  # unchanged, element-wise exact
end

@testset "apply_calibration! PSD matrix — ScalarCalibration uniform scale" begin
    # ScalarCalibration(-20.0f0): factor = 10^(20/10) = 100.0.
    # Every element of psd_linear multiplied by 100.
    cal   = ScalarCalibration(-20.0f0)
    psd   = [1.0 2.0; 3.0 4.0]
    freqs = [0.0, 100.0]
    apply_calibration!(psd, freqs, cal)
    @test psd ≈ [100.0 200.0; 300.0 400.0]  atol = 1e-10
end

@testset "apply_calibration! PSD matrix — TFCalibration flat -10 dB" begin
    # Flat TF of -10 dB: tf_lin = 0.1 everywhere.
    # Each PSD element divided by 0.1 → multiplied by 10.
    cal   = _tf([0.0, 1000.0], [-10.0, -10.0])
    psd   = [2.0 4.0; 6.0 8.0]
    freqs = [0.0, 500.0]
    orig  = copy(psd)
    apply_calibration!(psd, freqs, cal)
    @test psd ≈ orig .* 10.0  atol = 1e-10
end

@testset "apply_calibration! PSD matrix — TFCalibration ramp, per-bin check" begin
    # TF: -10 dB at 0 Hz (tf_lin = 0.1), -20 dB at 100 Hz (tf_lin = 0.01).
    # At 50 Hz (midpoint): linear power interpolation gives tf_lin = 0.055.
    # psd[1,:] at f=0 Hz:  value / 0.1   = value × 10.0
    # psd[2,:] at f=50 Hz: value / 0.055 ≈ value × 18.18...
    # psd[3,:] at f=100 Hz: value / 0.01  = value × 100.0
    cal    = _tf([0.0, 100.0], [-10.0, -20.0])
    freqs  = [0.0, 50.0, 100.0]
    psd    = ones(3, 2)   # all 1.0, size 3 bins × 2 frames
    apply_calibration!(psd, freqs, cal)

    @test psd[1, :] ≈ fill(1.0 / 0.1,   2)  atol = 1e-10   # 10.0
    @test psd[2, :] ≈ fill(1.0 / 0.055, 2)  atol = 1e-8    # ≈ 18.18
    @test psd[3, :] ≈ fill(1.0 / 0.01,  2)  atol = 1e-10   # 100.0
end

@testset "apply_calibration! PSD matrix — size mismatch asserts" begin
    cal   = _tf([0.0, 1000.0], [-10.0, -10.0])
    psd   = ones(3, 2)
    freqs_wrong = [0.0, 500.0]           # 2 freqs, psd has 3 rows
    @test_throws AssertionError apply_calibration!(psd, freqs_wrong, cal)
end
