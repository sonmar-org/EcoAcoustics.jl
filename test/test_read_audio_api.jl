using Test
using EcoAcoustics
using Dates

@testset "read_audio API" begin
    @test hasmethod(EcoAcoustics.read_audio, Tuple{AbstractString})
end

