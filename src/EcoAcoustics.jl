module EcoAcoustics

using Dates
using DataFrames
using ProgressMeter
using FFTW
using KernelAbstractions

# ─── Public API ───────────────────────────────────────────────────────────────

export Calibration, NoCalibration, ScalarCalibration, TFCalibration
export Audiodata, RecordingMetadata
export endtime, duration, nsamples

export apply_calibration!, apply_calibration, apply_calibration_psd!

export RecorderProfile, CalibrationProfile
export RECORDER_PROFILES, CALIBRATION_PROFILES
export lookup_calibration

export parse_filename
export read_audio

export AbstractAudioSource, SingleFileSource, IndexedFileSource
export time_range, read_audio_range, coverage_fraction
export build_index, load_index
export chunks, process_chunks

export SpectrogramResult, spectrogram, make_spectrogram_plan

# ─── Includes ─────────────────────────────────────────────────────────────────

include("audio/dsp_helpers.jl")
include("audio/calibration.jl")
include("audio/Audiodata.jl")
include("gpu/apply_calibration_kernel.jl")

include("recorders/recorders.jl")
include("recorders/rockhopper.jl")
include("recorders/sm3m.jl")
include("recorders/ls1x.jl")
include("recorders/snap.jl")

include("audio/parse_filename.jl")
include("audio/read_audio.jl")

include("sources/AbstractAudioSource.jl")
include("sources/SingleFileSource.jl")
include("sources/index_builder.jl")
include("sources/IndexedFileSource.jl")
include("sources/chunks.jl")

include("soundscape/spectrogram.jl")

end
