module EcoAcoustics

using Dates
using Statistics
using DataFrames
using ProgressMeter
using FFTW
using KernelAbstractions

# ─── Public API ───────────────────────────────────────────────────────────────

export Calibration, NoCalibration, ScalarCalibration, TFCalibration
export Audiodata, RecordingMetadata
export endtime, duration, nsamples

export apply_calibration!, apply_calibration

export RecorderProfile, CalibrationProfile
export RECORDER_PROFILES, CALIBRATION_PROFILES
export lookup_calibration
export NO_CALIBRATION_WARN_PREFIX

export AbstractRecorderProfile, RockhopperProfile, get_profile
export load_tf_calcurves

export parse_filename
export read_audio

export AbstractAudioSource, SingleFileSource, IndexedFileSource
export time_range, read_audio_range, coverage_fraction
export build_index, load_index
export chunks, process_chunks

export SpectrogramResult, spectrogram, make_spectrogram_plan
export PSDResult, compute_psd, psd_units, average_psd, to_dB

export BandSPL, SPLResult, compute_spl
export octave_bands, tol_bands, decidecade_bands, millidecade_bands
export compute_octave, compute_tol, compute_decidecade, compute_millidecade

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
include("soundscape/psd.jl")
include("soundscape/bands.jl")
include("soundscape/spl.jl")

end
