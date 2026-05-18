using Test
using EcoAcoustics
using Dates

@testset "Rockhopper filename parsing" begin
    fname = "139635MD01_197K_A6M_RH428_20231012_000654Z.flac"

    meta = EcoAcoustics.parse_filename(fname; recorder = "rockhopper")

    @test meta.site_id == "A6M"
    @test meta.recorder_id == "RH428"

    @test meta.timestamp isa Dates.DateTime
    @test meta.timestamp == DateTime(2023, 10, 12, 0, 6, 54)

    @test meta.timezone == "UTC"

    @test meta.lat === missing
    @test meta.lon === missing
end

@testset "SM3M filename parsing" begin
    fname = "T1-C__0__20170912_181500.wav"

    meta = EcoAcoustics.parse_filename(fname; recorder = "sm3m")

    @test meta.site_id == "T1-C"
    @test meta.recorder_id === nothing

    @test meta.timestamp == DateTime(2017, 9, 12, 18, 15, 0)
    @test meta.timezone === nothing
end

@testset "LS1X filename parsing" begin
    fname = "20210319T163400_2614231252441225_2.0dB_3.8V_ver2.00.wav"

    meta = EcoAcoustics.parse_filename(fname; recorder = "ls1x")

    @test meta.site_id === nothing
    @test meta.recorder_id == "2614231252441225"

    @test meta.timestamp == DateTime(2021, 3, 19, 16, 34, 0)
    @test meta.timezone === nothing
end

@testset "SNAP filename parsing (renamed)" begin
    fname = "T3C_180630_235000.wav"

    meta = EcoAcoustics.parse_filename(fname; recorder = "snap")

    @test meta.site_id == "T3C"
    @test meta.recorder_id === nothing
    
    @test meta.timestamp == DateTime(2018, 6, 30, 23, 50, 0)
    @test meta.timezone === nothing
end
