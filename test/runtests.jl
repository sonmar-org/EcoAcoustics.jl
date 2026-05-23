@info "LOAD_PATH in test process" Base.load_path()

using EcoAcoustics
using Test

const TEST_DIR = joinpath(@__DIR__, "test_files")

@testset "EcoAcoustics.jl" begin

    include("test_audiodata.jl")
    include("test_normalize_audio.jl")
    include("test_read_audio_api.jl")
    include("test_read_audio.jl")
    include("test_audio_io.jl")
    include("recorder_tests/test_filename_formats.jl")
    include("recorder_tests/test_calibration.jl")
    include("test_sources.jl")
    include("test_index.jl")
    include("test_indexed_file_source.jl")
    include("test_chunks.jl")
    include("test_apply_calibration.jl")
    include("test_tfcalibration.jl")
    include("test_rockhopper_profile.jl")
    include("test_spectrogram.jl")
    include("test_psd.jl")
    include("test_psd_manta_validation.jl")
end
