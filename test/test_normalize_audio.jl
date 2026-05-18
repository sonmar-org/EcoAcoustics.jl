using Test
using EcoAcoustics

@testset "_normalize_loaded_audio" begin
    data = randn(100, 1)       # fake mono matrix
    fs   = 48000.0

    sig, fs_out = EcoAcoustics._normalize_loaded_audio((data, fs))

    @test sig isa Vector{Float64}
    @test length(sig) == 100
    @test fs_out == 48000f0
end
