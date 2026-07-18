using Test
using EcoAcoustics

# Band-generator tests. All expected counts and key names are derived from the
# ANSI S1.6 (octave) and S1.11 (third-octave) preferred-frequency tables in
# bands.jl, and from the millidecade index formula f_c(n) = 10^(n/1000).
#
# Band edge formula:
#   octave:      (f_c / 2^(1/2), f_c × 2^(1/2))        base-2
#   third-octave (f_c / 10^(1/20), f_c × 10^(1/20))    base-10 decidecade (DD-30)
#                where f_c = 10^(n/10) is the exact decidecade center
#   millidecade: (10^((n-0.5)/1000), 10^((n+0.5)/1000))

# ─── octave_bands ─────────────────────────────────────────────────────────────

@testset "octave_bands: count and keys for [10, 10000] Hz" begin
    # ANSI S1.6 preferred centers in [10, 10000]:
    #   16, 31.5, 63, 125, 250, 500, 1000, 2000, 4000, 8000 → 10 bands
    # 16000 > 10000 is excluded; 8 < 10 is excluded.
    result = octave_bands(10.0, 10000.0)
    @test length(result) == 10

    # All expected keys present.
    for f in [16.0, 63.0, 125.0, 250.0, 500.0, 1000.0, 2000.0, 4000.0, 8000.0]
        @test haskey(result, Symbol("oct_$(round(Int, f))"))
    end
    # Non-integer preferred center: 31.5 → :oct_31_5
    @test haskey(result, :oct_31_5)

    # Boundary: 16000 > high_Hz = 10000, must be absent.
    @test !haskey(result, :oct_16000)
end

@testset "octave_bands: band edges for 1000 Hz" begin
    # f_c = 1000 Hz; edge factor = 2^(1/2) = √2 ≈ 1.41421
    # f_lo = 1000 / √2 ≈ 707.107   f_hi = 1000 × √2 ≈ 1414.214
    result = octave_bands(1000.0, 1000.0)
    @test length(result) == 1
    f_lo, f_hi = result[:oct_1000]
    @test f_lo ≈ 1000.0 / sqrt(2) atol=1e-10
    @test f_hi ≈ 1000.0 * sqrt(2) atol=1e-10
end

@testset "octave_bands: non-integer label :oct_31_5" begin
    result = octave_bands(20.0, 50.0)
    # Only 31.5 has its center in [20, 50].
    @test haskey(result, :oct_31_5)
    f_lo, f_hi = result[:oct_31_5]
    @test f_lo ≈ 31.5 / sqrt(2) atol=1e-10
    @test f_hi ≈ 31.5 * sqrt(2) atol=1e-10
end

@testset "octave_bands: empty result for impossible range" begin
    # No ANSI octave center falls in (15, 16) Hz — gap between 8 and 16 Hz bands.
    @test isempty(octave_bands(15.1, 15.9))
end

# ─── tol_bands ────────────────────────────────────────────────────────────────

@testset "tol_bands: count and keys for [10, 100] Hz" begin
    # ANSI S1.11 preferred centers in [10, 100]:
    #   10, 12.5, 16, 20, 25, 31.5, 40, 50, 63, 80, 100 → 11 bands
    result = tol_bands(10.0, 100.0)
    @test length(result) == 11

    @test haskey(result, :tol_10)
    @test haskey(result, :tol_12_5)   # non-integer: 12.5
    @test haskey(result, :tol_31_5)   # non-integer: 31.5
    @test haskey(result, :tol_100)

    # 125 > 100 must be absent.
    @test !haskey(result, :tol_125)
end

@testset "tol_bands: band edges for 1000 Hz (base-10 decidecade, DD-30)" begin
    # :tol_1000 → decidecade band n = round(10·log10(1000)) = 30 → exact center
    # f_c = 10^3 = 1000 Hz; base-10 half-bandwidth factor = 10^(1/20) ≈ 1.12202.
    # f_lo = 1000 / 10^(1/20) ≈ 891.25   f_hi = 1000 × 10^(1/20) ≈ 1122.02
    result = tol_bands(1000.0, 1000.0)
    @test length(result) == 1
    f_lo, f_hi = result[:tol_1000]
    @test f_lo ≈ 1000.0 / 10.0^(1/20) atol=1e-10
    @test f_hi ≈ 1000.0 * 10.0^(1/20) atol=1e-10
end

@testset "tol_bands: base-10 exact center for non-round nominal (63 Hz, DD-30)" begin
    # :tol_63 keeps the nominal label 63 but uses the EXACT base-10 center
    # 10^1.8 ≈ 63.096 Hz — this is what makes it match ISO/ADEON/PAMGuide.
    result = tol_bands(63.0, 63.0)
    f_lo, f_hi = result[:tol_63]
    fc = 10.0^(round(Int, 10*log10(63.0)) / 10)     # 10^1.8 ≈ 63.0957
    @test fc ≈ 63.0957344 atol=1e-4
    @test f_lo ≈ fc / 10.0^(1/20) atol=1e-10
    @test f_hi ≈ fc * 10.0^(1/20) atol=1e-10
end

@testset "tol_bands: upper-bound clipping" begin
    # Centers at 63000 and 80000 fall in [10, 80000]; 100000 > 80000 is excluded.
    result = tol_bands(10.0, 80000.0)
    @test haskey(result, :tol_80000)
    @test !haskey(result, :tol_100000)
end

# ─── decidecade_bands ─────────────────────────────────────────────────────────

@testset "decidecade_bands is identical to tol_bands" begin
    # By definition the two are the same band scheme (DD-20).
    @test decidecade_bands(10.0, 10000.0) == tol_bands(10.0, 10000.0)
    @test decidecade_bands(1.0, 200000.0) == tol_bands(1.0, 200000.0)
end

# ─── millidecade_bands ────────────────────────────────────────────────────────

@testset "millidecade_bands: count for [900, 1100] Hz" begin
    # n_lo = ceil(1000 × log10(900)) = ceil(2954.24) = 2955
    # n_hi = floor(1000 × log10(1100)) = floor(3041.39) = 3041
    # count = 3041 - 2955 + 1 = 87
    result = millidecade_bands(900.0, 1100.0)
    @test length(result) == 87

    # First and last band indices must be present.
    @test haskey(result, :mdec_2955)
    @test haskey(result, :mdec_3041)

    # :mdec_3000 has center 10^3 = 1000 Hz, which is in [900, 1100].
    @test haskey(result, :mdec_3000)
end

@testset "millidecade_bands: center and edges for :mdec_3000" begin
    # n = 3000 → f_c = 10^(3000/1000) = 10^3 = 1000 Hz
    # f_lo = 10^(2999.5/1000) ≈ 999.885 Hz
    # f_hi = 10^(3000.5/1000) ≈ 1000.115 Hz
    result = millidecade_bands(900.0, 1100.0)
    f_lo, f_hi = result[:mdec_3000]
    @test f_lo ≈ 10.0^(2999.5 / 1000.0) atol=1e-10
    @test f_hi ≈ 10.0^(3000.5 / 1000.0) atol=1e-10
    # Center frequency (geometric mean of edges) ≈ 1000 Hz.
    @test sqrt(f_lo * f_hi) ≈ 1000.0 atol=1e-6
end

@testset "millidecade_bands: MANTA index convention" begin
    # :mdec_4000 → n=4000, f_c = 10^4 = 10000 Hz
    # :mdec_2000 → n=2000, f_c = 10^2 = 100 Hz
    r10k = millidecade_bands(9900.0, 10100.0)
    @test haskey(r10k, :mdec_4000)
    r100 = millidecade_bands(99.0, 101.0)
    @test haskey(r100, :mdec_2000)
end

@testset "millidecade_bands: ArgumentError for non-positive low_Hz" begin
    @test_throws ArgumentError millidecade_bands(0.0, 1000.0)
    @test_throws ArgumentError millidecade_bands(-10.0, 1000.0)
end

@testset "millidecade_bands: empty result for impossible range" begin
    # high_Hz < low_Hz → n_hi < n_lo → empty loop.
    @test isempty(millidecade_bands(1100.0, 900.0))
end
