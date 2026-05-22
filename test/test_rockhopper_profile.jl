using Test
using EcoAcoustics

# Tests for the Rockhopper calibration profile (Deliverable 2).
#
# Every test that validates a specific tf_dB value is computed analytically:
#   tf_dB = AnalogSensitivity_dB - 20*log10(5.0)
#   20*log10(5.0) = 13.9794... dB
# Values are checked against the shipped CSV, not against the legacy rockhopper_TF.csv.
# The format === :rockhopper_calcurves_csv assertion is the key guard that proves
# all tests route through load_tf_calcurves (the new path), not _load_tf_csv (the old path).

# ─── RockhopperProfile singleton ─────────────────────────────────────────────

@testset "RockhopperProfile singleton identity" begin
    p1 = RockhopperProfile()
    p2 = RockhopperProfile()
    # === checks object identity (same memory address), not just equality.
    # The const _RH_PROFILE guarantees this.
    @test p1 === p2
end

@testset "RockhopperProfile.Vmax_peak_V" begin
    @test RockhopperProfile().Vmax_peak_V === 5.0
end

@testset "RockhopperProfile.tf — format confirms new path, not legacy" begin
    # If this were routing through _load_tf_csv, format would be :legacy_two_column_csv.
    @test RockhopperProfile().tf.format === :rockhopper_calcurves_csv
end

@testset "RockhopperProfile.tf — 2798 data points" begin
    tf = RockhopperProfile().tf
    @test length(tf.frequency) == 2798
    @test length(tf.tf_dB)     == 2798
end

@testset "RockhopperProfile.tf — first point (f=0 Hz)" begin
    # CSV row 1 (after header): Frequency_Hz=0, AnalogSensitivity=-215.74
    # tf_dB[1] = -215.74 - 20*log10(5.0) = -215.74 - 13.9794 = -229.7194
    tf = RockhopperProfile().tf
    @test tf.frequency[1] ≈ 0.0     atol = 1e-10
    @test tf.tf_dB[1]     ≈ -229.7194  atol = 1e-3
end

@testset "RockhopperProfile.tf — last point (f≈99884 Hz)" begin
    # CSV last row: Frequency_Hz≈99884.94, AnalogSensitivity≈-206.6965
    # tf_dB[end] = -206.6965 - 13.9794 ≈ -220.676
    tf = RockhopperProfile().tf
    @test tf.frequency[end] ≈ 99884.937  atol = 1.0
    @test tf.tf_dB[end]     ≈ -220.676   atol = 1e-2
end

@testset "RockhopperProfile.tf — conversion_notes mentions 13.9794" begin
    notes = RockhopperProfile().tf.conversion_notes
    @test occursin("13.9794", notes)
end

@testset "RockhopperProfile.tf — source path points to shipped CSV" begin
    src = RockhopperProfile().tf.source
    @test occursin("rockhopper_tf_calibration.csv", src)
    @test isfile(src)   # the file must actually exist at the recorded path
end

# ─── get_profile dispatch ─────────────────────────────────────────────────────

@testset "get_profile(:rockhopper) === RockhopperProfile()" begin
    # Both the Symbol path and the direct call must return the same singleton.
    @test get_profile(:rockhopper) === RockhopperProfile()
end

@testset "get_profile — unrecognised recorder raises ArgumentError" begin
    err = @test_throws ArgumentError get_profile(:nonexistent_recorder)
    # Error message must name the unrecognised symbol.
    @test occursin("nonexistent_recorder", err.value.msg)
end

# ─── load_tf_calcurves ────────────────────────────────────────────────────────

@testset "load_tf_calcurves — shipped CSV matches get_profile singleton" begin
    # Load the shipped CSV directly and compare pointwise with the singleton.
    # This verifies that the singleton was built from the same file it claims to be.
    shipped_path = joinpath(pkgdir(EcoAcoustics), "src", "recorders",
                            "calibration_data", "rockhopper_tf_calibration.csv")
    direct_tf = load_tf_calcurves(shipped_path; vmax_peak_V = 5.0)
    singleton_tf = RockhopperProfile().tf

    @test maximum(abs, direct_tf.tf_dB .- singleton_tf.tf_dB) < 1e-10
    @test maximum(abs, direct_tf.frequency .- singleton_tf.frequency) < 1e-10
end

@testset "load_tf_calcurves — missing required column raises ArgumentError" begin
    # Write a CSV that has Frequency_Hz but not AnalogSensitivity.
    path = tempname() * ".csv"
    open(path, "w") do io
        println(io, "Frequency_Hz,SomeOtherColumn")
        println(io, "100.0,-50.0")
        println(io, "200.0,-51.0")
    end
    err = @test_throws ArgumentError load_tf_calcurves(path; vmax_peak_V = 5.0)
    # Error message must name the missing column and list required ones.
    @test occursin("AnalogSensitivity_dB_re_1VperRefPress", err.value.msg)
    @test occursin("Frequency_Hz",                          err.value.msg)
    rm(path; force = true)
end

@testset "load_tf_calcurves — extra columns are silently accepted" begin
    # A CSV with required columns plus an extra one should load without error.
    # The extra column must not affect the tf_dB values.
    path = tempname() * ".csv"
    open(path, "w") do io
        println(io, "Frequency_Hz,AnalogSensitivity_dB_re_1VperRefPress,ExtraColumn")
        println(io, "10.0,-200.0,999.9")
        println(io, "100.0,-198.0,888.8")
        println(io, "1000.0,-195.0,777.7")
    end
    tf = load_tf_calcurves(path; vmax_peak_V = 5.0)
    shift = 20.0 * log10(5.0)
    @test tf.tf_dB ≈ [-200.0, -198.0, -195.0] .- shift  atol = 1e-10
    @test tf.frequency ≈ [10.0, 100.0, 1000.0]           atol = 1e-10
    rm(path; force = true)
end

@testset "load_tf_calcurves — gain conversion (hand-computed)" begin
    # At vmax_peak_V = 10.0: shift = 20*log10(10) = 20.0 dB exactly.
    # tf_dB should equal AnalogSensitivity - 20.0 exactly.
    path = tempname() * ".csv"
    open(path, "w") do io
        println(io, "Frequency_Hz,AnalogSensitivity_dB_re_1VperRefPress")
        println(io, "0.0,-210.0")
        println(io, "1000.0,-205.0")
    end
    tf = load_tf_calcurves(path; vmax_peak_V = 10.0)
    @test tf.tf_dB ≈ [-230.0, -225.0]  atol = 1e-10
    rm(path; force = true)
end

@testset "load_tf_calcurves — non-ascending frequency rejected" begin
    path = tempname() * ".csv"
    open(path, "w") do io
        println(io, "Frequency_Hz,AnalogSensitivity_dB_re_1VperRefPress")
        println(io, "1000.0,-200.0")
        println(io, "100.0,-198.0")   # descending — invalid
    end
    @test_throws ArgumentError load_tf_calcurves(path; vmax_peak_V = 5.0)
    rm(path; force = true)
end
