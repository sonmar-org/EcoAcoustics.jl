using Test
using EcoAcoustics
using Logging
using Statistics

# All expected values are derived analytically from first principles.
# No "run it and see" tests — every assertion has a derivation comment.
#
# Shared synthetic PSD design:
#   fs = 1000 Hz  (Nyquist = 500 Hz)
#   nfft = 1000   (df = fs / nfft = 1.0 Hz, exactly)
#   n_freqs = 501 (rfft gives nfft/2 + 1 bins: 0, 1, 2, ..., 500 Hz)
#   n_frames = 5
#   psd_level = 1.0 µPa²/Hz, uniform across all bins and frames
#   is_calibrated = true  (so psd_units returns :µPa²_per_Hz)
#
# With this setup, for any band [f_lo, f_hi]:
#   n_bins     = searchsortedlast(freqs, f_hi) - searchsortedfirst(freqs, f_lo) + 1
#   band_power = n_bins × df × psd_level    (µPa²)
#   SPL_water  = 10 × log10(band_power)     (dB re 1 µPa, pref=1.0)
#   SPL_air    = 10 × log10(band_power / 400)  (dB re 20 µPa, pref=20.0)
#
# Because the PSD is flat and all frames are identical, every aggregate
# statistic (mean_dB, percentiles) equals the single-frame SPL value.

function _make_test_psd(;
        n_freqs      = 501,
        n_frames     = 5,
        fs           = Float32(1000.0),
        nfft         = 1000,
        psd_level    = 1.0,
        is_calibrated = true)
    # window_energy = nfft for a rectangular window; unused by compute_spl
    # but required by the PSDResult constructor.
    freqs   = collect(range(0.0, Float64(fs) / 2; length=n_freqs))
    time    = collect(0.0 : 1.0 : Float64(n_frames - 1))
    psd_mat = fill(psd_level, n_freqs, n_frames)
    return PSDResult(psd_mat, freqs, time, fs, Float64(nfft), nfft,
                     NoCalibration(), is_calibrated)
end

# ─── Default broadband band ───────────────────────────────────────────────────

@testset "compute_spl: default broadband band" begin
    psd    = _make_test_psd()
    result = compute_spl(psd)   # no bands argument → :broadband => (10.0, Nyquist)

    @test haskey(result.bands, :broadband)
    @test result.bands[:broadband].band == (10.0, Float64(psd.fs) / 2.0)

    # SPLResult metadata propagated from PSDResult.
    @test result.time === psd.time
    @test result.fs   === psd.fs
    @test result.units       === :dB_re_1µPa
    @test result.environment === :water
end

# ─── Flat PSD at known level → analytical SPL ────────────────────────────────

@testset "compute_spl: flat PSD, single band, known analytical value" begin
    psd = _make_test_psd()
    # Band [100, 200]:
    #   searchsortedfirst(0:500, 100) → index 101 (freqs[101] = 100.0)
    #   searchsortedlast(0:500, 200)  → index 201 (freqs[201] = 200.0)
    #   n_bins = 201 - 101 + 1 = 101
    #   band_power = 101 × 1.0 Hz × 1.0 µPa²/Hz = 101 µPa²
    #   SPL = 10 × log10(101) ≈ 20.04 dB re 1 µPa
    expected_power = 101 * 1.0 * 1.0   # n_bins × df × psd_level (µPa²)
    expected_spl   = 10.0 * log10(expected_power)

    result = compute_spl(psd; bands = Dict(:test_band => (100.0, 200.0)))
    b      = result.bands[:test_band]

    # Per-frame series: all 5 frames should equal expected_spl exactly.
    @test all(v -> v ≈ expected_spl, b.spl_dB)
    @test length(b.spl_dB) == 5

    # Aggregate statistics: flat series → all stats equal the one value.
    @test b.mean_dB   ≈ expected_spl  atol=1e-10
    @test b.median_dB ≈ expected_spl  atol=1e-10
    @test b.L1_dB     ≈ expected_spl  atol=1e-10
    @test b.L99_dB    ≈ expected_spl  atol=1e-10

    # Band edges stored verbatim.
    @test b.band == (100.0, 200.0)
end

# ─── Two non-overlapping bands ────────────────────────────────────────────────

@testset "compute_spl: two non-overlapping bands" begin
    psd = _make_test_psd()
    # Band :low [100, 200]: 101 bins → 101 µPa² → 10×log10(101)
    # Band :high [300, 350]: 51 bins → 51 µPa²  → 10×log10(51)
    #   searchsortedfirst(freqs, 300) = 301, searchsortedlast(freqs, 350) = 351
    #   n_bins = 351 - 301 + 1 = 51
    expected_low  = 10.0 * log10(101.0)
    expected_high = 10.0 * log10(51.0)

    result = compute_spl(psd; bands = Dict(
        :low  => (100.0, 200.0),
        :high => (300.0, 350.0)))

    @test haskey(result.bands, :low)
    @test haskey(result.bands, :high)
    @test result.bands[:low].mean_dB  ≈ expected_low  atol=1e-10
    @test result.bands[:high].mean_dB ≈ expected_high atol=1e-10
end

# ─── Overlapping bands are independent ────────────────────────────────────────

@testset "compute_spl: overlapping bands computed independently" begin
    psd = _make_test_psd()
    # :wide [100, 350]:  n_bins = 351-101+1 = 251 → 251 µPa²
    # :lower [100, 250]: n_bins = 251-101+1 = 151 → 151 µPa²
    # :upper [200, 350]: searchsortedfirst(freqs,200)=201; 351-201+1=151 → 151 µPa²
    # Each band's result is independent of the others.
    expected_wide  = 10.0 * log10(251.0)
    expected_lower = 10.0 * log10(151.0)
    expected_upper = 10.0 * log10(151.0)

    result = compute_spl(psd; bands = Dict(
        :wide  => (100.0, 350.0),
        :lower => (100.0, 250.0),
        :upper => (200.0, 350.0)))

    @test result.bands[:wide].mean_dB  ≈ expected_wide  atol=1e-10
    @test result.bands[:lower].mean_dB ≈ expected_lower atol=1e-10
    @test result.bands[:upper].mean_dB ≈ expected_upper atol=1e-10
end

# ─── Sub-10-Hz warning ────────────────────────────────────────────────────────

@testset "compute_spl: sub-10-Hz band triggers consolidated warning" begin
    psd = _make_test_psd()
    # :sub10 has low_Hz = 5.0 < 10 → warn; :normal does not → no warn for it.
    # The consolidated warning lists all offending band labels, so "sub10"
    # must appear in the message string.
    @test_logs (:warn, r"sub10") compute_spl(psd; bands = Dict(
        :sub10  => (5.0, 100.0),
        :normal => (100.0, 200.0)))

    # A band entirely above 10 Hz must not trigger any warning.
    @test_logs compute_spl(psd; bands = Dict(:normal => (100.0, 200.0)))
end

# ─── Band validation throws ───────────────────────────────────────────────────

@testset "compute_spl: band validation ArgumentError" begin
    psd = _make_test_psd()   # fs=1000, Nyquist=500

    # high_Hz > Nyquist (500 Hz).
    @test_throws ArgumentError compute_spl(psd;
        bands = Dict(:over_nyq => (100.0, 600.0)))

    # low_Hz > high_Hz (reversed edges).
    @test_throws ArgumentError compute_spl(psd;
        bands = Dict(:reversed => (500.0, 100.0)))

    # low_Hz == high_Hz (zero-width; low ≥ high fires the order check).
    @test_throws ArgumentError compute_spl(psd;
        bands = Dict(:zero_width => (100.0, 100.0)))

    # Multiple offending bands: all labels must appear in one error message.
    err = nothing
    try
        compute_spl(psd; bands = Dict(
            :bad_nyq => (100.0, 600.0),
            :ok      => (100.0, 200.0)))
    catch e
        err = e
    end
    @test err isa ArgumentError
    @test occursin("bad_nyq", err.msg)
end

# ─── Calibration assertion (DD-21) ───────────────────────────────────────────

@testset "compute_spl: uncalibrated PSD raises AssertionError (DD-21)" begin
    psd_uncal = _make_test_psd(; is_calibrated = false)

    # Must throw AssertionError (not ArgumentError or MethodError).
    @test_throws AssertionError compute_spl(psd_uncal)

    # Message must mention :µPa²_per_Hz so the user knows what was expected.
    err = try
        compute_spl(psd_uncal)
        nothing
    catch e
        e
    end
    @test err isa AssertionError
    @test occursin("µPa²_per_Hz", err.msg)
end

# ─── Air environment ──────────────────────────────────────────────────────────

@testset "compute_spl: air reference gives SPL 20×log10(20) ≈ 26 dB lower" begin
    psd   = _make_test_psd()
    bands = Dict(:b => (100.0, 200.0))

    result_water = compute_spl(psd; bands = bands, environment = :water)
    result_air   = compute_spl(psd; bands = bands, environment = :air)

    # Water SPL   = 10×log10(P / 1²)   = 10×log10(P)
    # Air SPL     = 10×log10(P / 20²)  = 10×log10(P) − 10×log10(400)
    # Difference  = 10×log10(400) = 20×log10(20) ≈ 26.021 dB
    expected_diff = 20.0 * log10(20.0)

    @test result_water.bands[:b].mean_dB - result_air.bands[:b].mean_dB ≈
          expected_diff atol=1e-10

    # All aggregate statistics shift by the same constant.
    @test result_water.bands[:b].L99_dB - result_air.bands[:b].L99_dB ≈
          expected_diff atol=1e-10
end

# ─── units and environment fields ─────────────────────────────────────────────

@testset "compute_spl: units and environment fields set correctly" begin
    psd = _make_test_psd()

    water = compute_spl(psd; environment = :water)
    @test water.units       === :dB_re_1µPa
    @test water.environment === :water

    air = compute_spl(psd; environment = :air)
    @test air.units       === :dB_re_20µPa
    @test air.environment === :air
end

# ─── D3: compute_spl(Audiodata) wrapper ──────────────────────────────────────

@testset "compute_spl(Audiodata): round-trip equivalence with primitive" begin
    # The wrapper must produce results identical to calling compute_psd then
    # compute_spl manually. Use ScalarCalibration so the calibration cascade
    # resolves deterministically at step 2 (audio.calibration isa !NoCalibration).
    N      = 4800
    fs     = 48000.0
    signal = randn(Float64, N)
    audio  = Audiodata(signal, Float32(fs), DateTime(2023, 1, 1);
                       calibration = ScalarCalibration(-153.0f0))
    bands  = Dict(:b => (100.0, 1000.0))

    # Wrapper path.
    result_wrapper = compute_spl(audio;
                                 bands = bands, window_seconds = 0.1)

    # Manual path: same audio, same parameters.
    psd            = compute_psd(audio; window_seconds = 0.1)
    result_direct  = compute_spl(psd; bands = bands)

    b_w = result_wrapper.bands[:b]
    b_d = result_direct.bands[:b]

    # spl_dB must be element-wise identical (same FFT plan size → same result).
    @test b_w.spl_dB   ≈ b_d.spl_dB   atol=1e-12
    @test b_w.mean_dB  ≈ b_d.mean_dB  atol=1e-12
    @test b_w.L99_dB   ≈ b_d.L99_dB   atol=1e-12
    @test result_wrapper.time  == result_direct.time
    @test result_wrapper.units === result_direct.units
end

@testset "compute_spl(Audiodata): Rockhopper auto-calibration" begin
    # compute_psd auto-resolves to get_profile(:rockhopper).tf (DD-13, DD-14
    # cascade step 3) when recorder="rockhopper" and no explicit calibration
    # is supplied. The resulting PSDResult is calibrated, so compute_spl's
    # assertion passes and units reflect physical µPa.
    N     = 4800
    fs    = 48000.0
    audio = Audiodata(randn(Float64, N), Float32(fs), DateTime(2023, 1, 1);
                      recorder = "rockhopper")

    result = compute_spl(audio; window_seconds = 0.1)

    @test result.units === :dB_re_1µPa
    @test result.environment === :water
    # Default broadband band spans 10 Hz to Nyquist.
    @test result.bands[:broadband].band == (10.0, Float64(fs) / 2.0)
end

@testset "compute_spl(Audiodata): custom bands forwarded correctly" begin
    N     = 4800
    fs    = 48000.0
    audio = Audiodata(randn(Float64, N), Float32(fs), DateTime(2023, 1, 1);
                      calibration = ScalarCalibration(-153.0f0))
    bands = Dict(:low => (100.0, 1000.0), :high => (5000.0, 10000.0))

    result = compute_spl(audio; bands = bands, window_seconds = 0.1)

    @test haskey(result.bands, :low)
    @test haskey(result.bands, :high)
    @test result.bands[:low].band  == (100.0,  1000.0)
    @test result.bands[:high].band == (5000.0, 10000.0)
    # No :broadband default — the explicit Dict was forwarded, not replaced.
    @test !haskey(result.bands, :broadband)
end

@testset "compute_spl(Audiodata): window_seconds changes spl_dB frame count" begin
    # Verify that window_seconds is forwarded to compute_psd → spectrogram.
    # With N=48000, fs=48000, overlap=0.5:
    #   window_seconds=0.1 → window_length=4800, hop=2400:
    #     n_frames = floor((48000-4800)/2400) + 1 = 19
    #   window_seconds=0.5 → window_length=24000, hop=12000:
    #     n_frames = floor((48000-24000)/12000) + 1 = 3
    N     = 48000
    fs    = 48000.0
    audio = Audiodata(randn(Float64, N), Float32(fs), DateTime(2023, 1, 1);
                      calibration = ScalarCalibration(-153.0f0))

    result_short = compute_spl(audio; window_seconds = 0.1)
    result_long  = compute_spl(audio; window_seconds = 0.5)

    @test length(result_short.bands[:broadband].spl_dB) == 19
    @test length(result_long.bands[:broadband].spl_dB)  == 3
end
