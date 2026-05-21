using Test
using EcoAcoustics

@testset "Calibration lookup (SM3M)" begin
    # Minimal metadata NamedTuple matching parse_filename output shape
    meta = (
        timestamp   = nothing,
        timezone    = nothing,
        lat         = missing,
        lon         = missing,
        site_id     = nothing,
        recorder_id = nothing,
    )

    cal = EcoAcoustics.lookup_calibration("dummy.wav", "sm3m", meta; strict=false)

    @test cal isa EcoAcoustics.ScalarCalibration

    # From hydrophones.jl: -165 + 12 + 0 + 20*log10(1/1) = -153 dB
    @test cal.system_sensitivity_dB ≈ -153f0 atol = 1e-3
end

