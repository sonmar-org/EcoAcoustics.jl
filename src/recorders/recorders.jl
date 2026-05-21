"""
    RecorderProfile

Describes how a specific recorder family encodes metadata in its filenames.
One profile per recorder family. All profiles are registered in
`RECORDER_PROFILES` and looked up by `parse_filename`.

Fields
------
- `name::String`: Human-readable recorder name.
- `filename_style::Symbol`: `:tokenized` — basename is split on `split_char`.
- `split_char::Union{Char,Nothing}`: Token delimiter.
- `date_token`, `time_token`: 1-based indices of date and time tokens when
  they appear separately in the filename.
- `datetime_token`: 1-based index of a combined date-time token.
  Use either (`date_token` + `time_token`) or `datetime_token`, not both.
- `datetime_format::Union{String,Nothing}`: `Dates.DateFormat` string for
  parsing the combined date-time string.
- `site_token`, `recorder_id_token`: Token indices for site and recorder ID.
- `timezone_source::Symbol`: `:suffix_Z` if a trailing `Z` marks UTC;
  `:none` otherwise.
- `has_xml::Bool`: Whether this recorder may supply sidecar XML metadata.
- `fix_wav_header::Bool`: Whether WAV files from this recorder require
  in-memory RIFF header correction before reading. Set `true` for recorders
  that systematically write a data-chunk size larger than the actual file
  content (e.g. DMON2). When `true`, `read_audio` calls `_wavread_corrected`
  silently instead of `wavread`. Default `false` for all other recorders.
- `parser::Union{Nothing,Function}`: Optional function that overrides all
  generic token parsing. Required when filename rules need custom logic
  (e.g. 2-digit years, XML lookups).
"""
struct RecorderProfile
    name::String
    filename_style::Symbol
    split_char::Union{Char,Nothing}
    date_token::Union{Int,Nothing}
    time_token::Union{Int,Nothing}
    datetime_token::Union{Int,Nothing}
    datetime_format::Union{String,Nothing}
    site_token::Union{Int,Nothing}
    recorder_id_token::Union{Int,Nothing}
    timezone_source::Symbol
    has_xml::Bool
    fix_wav_header::Bool
    parser::Union{Nothing,Function}
end

# Keyword constructor so per-recorder files only name the fields they use.
# Unspecified fields default to nothing/:none/false.
function RecorderProfile(name::String;
                         filename_style::Symbol                 = :tokenized,
                         split_char::Union{Char,Nothing}        = nothing,
                         date_token::Union{Int,Nothing}         = nothing,
                         time_token::Union{Int,Nothing}         = nothing,
                         datetime_token::Union{Int,Nothing}     = nothing,
                         datetime_format::Union{String,Nothing} = nothing,
                         site_token::Union{Int,Nothing}         = nothing,
                         recorder_id_token::Union{Int,Nothing}  = nothing,
                         timezone_source::Symbol                = :none,
                         has_xml::Bool                          = false,
                         fix_wav_header::Bool                   = false,
                         parser::Union{Nothing,Function}        = nothing)
    RecorderProfile(name, filename_style, split_char,
                    date_token, time_token, datetime_token, datetime_format,
                    site_token, recorder_id_token,
                    timezone_source, has_xml, fix_wav_header, parser)
end

"""
    CalibrationProfile

Static calibration parameters for a recorder family.

The total scalar system sensitivity (dB re 1 V/µPa) is computed as:
    sensitivity + preamp_gain + board_gain + 20*log10(1/Vadc_0pk)

Fields
------
- `label::String`: Human-readable description.
- `sensitivity::Float64`: Hydrophone sensitivity [dB re 1 V/µPa].
- `preamp_gain::Float64`: Preamplifier gain [dB].
- `board_gain::Float64`: ADC board gain [dB].
- `Vadc_0pk::Float64`: ADC full-scale peak voltage [V].
- `tf_path::Union{String,Nothing}`: Path to a transfer-function file for
  frequency-dependent calibration (v2). `nothing` for scalar-only recorders.
"""
struct CalibrationProfile
    label::String
    sensitivity::Float64
    preamp_gain::Float64
    board_gain::Float64
    Vadc_0pk::Float64
    tf_path::Union{String,Nothing}
end

# tf_path defaults to nothing for scalar-only recorders.
function CalibrationProfile(label::String,
                            sensitivity::Float64,
                            preamp_gain::Float64,
                            board_gain::Float64,
                            Vadc_0pk::Float64;
                            tf_path::Union{String,Nothing} = nothing)
    CalibrationProfile(label, sensitivity, preamp_gain, board_gain, Vadc_0pk, tf_path)
end

"""
    RECORDER_PROFILES

Registry of filename-parsing profiles, keyed by lowercase recorder name
(e.g. `"sm3m"`, `"rockhopper"`). Each entry is a [`RecorderProfile`](@ref)
populated by the per-recorder source file under `src/recorders/`.

To add a new recorder, create `src/recorders/yourrecorder.jl`, add an entry
here, and include the file in `src/EcoAcoustics.jl`.
"""
const RECORDER_PROFILES = Dict{String,RecorderProfile}()

"""
    CALIBRATION_PROFILES

Registry of scalar calibration profiles, keyed by lowercase recorder name.
Each entry is a [`CalibrationProfile`](@ref) populated by the per-recorder
source file under `src/recorders/`.

Recorders without a registered profile return [`NoCalibration`](@ref) from
[`lookup_calibration`](@ref). Frequency-dependent (TF) calibration is handled
separately via [`TFCalibration`](@ref) and is not stored here.
"""
const CALIBRATION_PROFILES = Dict{String,CalibrationProfile}()

"""
    lookup_calibration(path, recorder, meta; strict=false)

Purpose:     Return a `Calibration` object for a recording based on its
             recorder family. Warns when no profile is found; raises when
             `strict=true`.

Arguments:
- `path::AbstractString`: File path. Currently unused; reserved for per-file
  TF calibration lookup in v2.
- `recorder::AbstractString`: Recorder family name (e.g. `"sm3m"`).
- `meta`: Named tuple from `parse_filename`. Currently unused; reserved for
  per-serial calibration overrides in v2.
- `strict::Bool = false`: If `true`, throw instead of returning `NoCalibration`.

Returns:     A `Calibration` subtype: `ScalarCalibration` when a profile
             exists, `NoCalibration` otherwise.

Constraints: Only scalar calibration is returned in v1. TF calibration
             (frequency-dependent) is constructed separately when needed.

Fails when:  `strict=true` and no calibration profile is registered for
             the given recorder.

Example:
```julia
cal = lookup_calibration("file.wav", "sm3m", meta)
```
"""
function lookup_calibration(path::AbstractString,
                            recorder::AbstractString,
                            meta;
                            strict::Bool = false)
    key = lowercase(String(recorder))

    if haskey(CALIBRATION_PROFILES, key)
        cp = CALIBRATION_PROFILES[key]
        if cp.tf_path !== nothing
            freqs, tf_db = _load_tf_csv(cp.tf_path)
            return TFCalibration(freqs, tf_db)
        end
        total_sens_db = cp.sensitivity + cp.preamp_gain + cp.board_gain +
                        20 * log10(1 / cp.Vadc_0pk)
        return ScalarCalibration(Float32(total_sens_db))
    end

    if strict
        error("lookup_calibration: no calibration profile for recorder " *
              "\"$recorder\". Add a CalibrationProfile to CALIBRATION_PROFILES " *
              "or pass strict=false to proceed with NoCalibration().")
    end

    @warn "No calibration entry for recorder; returning NoCalibration(). " *
          "Metric computations will be in raw ADC units, not physical units." recorder=recorder

    return NoCalibration()
end
