"""
    RecordingMetadata

Provenance fields for a single audio recording: where it was made, with what,
and any deployment identifiers. Kept separate from the physics fields on
`Audiodata` so the two concerns don't blur.

Fields
------
- `timezone::Union{String,Nothing}`:
    Timezone the recorder clock was set to, as an IANA name (e.g.
    `"America/New_York"`, `"UTC"`). `nothing` if unknown. When `timezone`
    is not `"UTC"`, `starttime` reflects local recorder time and has not
    been converted to UTC — the caller is responsible for conversion.

- `lat::Union{Float64,Missing}` / `lon::Union{Float64,Missing}`:
    Recorder location in decimal degrees. `missing` if unknown.

- `site_id::Union{String,Nothing}`:
    Deployment or site identifier (e.g. `"T1-C"`). `nothing` if unknown.

- `recorder::String`:
    Recorder family name (e.g. `"sm3m"`, `"rockhopper"`). Defaults to
    `"unknown"`.

- `recorder_id::Union{String,Nothing}`:
    Recorder serial number or unit ID. Distinct from the hydrophone serial
    number, which is tracked separately in v2. `nothing` if unknown.
"""
struct RecordingMetadata
    timezone::Union{String,Nothing}
    lat::Union{Float64,Missing}
    lon::Union{Float64,Missing}
    site_id::Union{String,Nothing}
    recorder::String
    recorder_id::Union{String,Nothing}
end


"""
    Audiodata

Container for a single mono audio recording and its associated metadata.

Fields
------
- `sig::Vector{Float64}`:
    Raw audio samples in linear amplitude. Always `Float64`, always
    uncalibrated. May be empty (length 0) when used as a no-overlap sentinel
    by `read_audio_range`.

- `fs::Float32`:
    Sampling rate in Hz. Must be > 0.

- `starttime::DateTime`:
    Recording start time as reported by the recorder clock. May reflect local
    time rather than UTC — check `metadata.timezone`.

- `is_calibrated::Bool`:
    `false` at construction. Set to `true` only by `apply_calibration!`.
    Metric functions (SPL, PSD, etc.) check this flag and warn when it is
    `false`.

- `calibration::Calibration`:
    Calibration object associated with this recording. Stored but not applied
    at I/O time. One of `NoCalibration`, `ScalarCalibration`, `TFCalibration`.

- `metadata::RecordingMetadata`:
    Provenance fields (location, recorder identity, site, timezone).

Notes
-----
- Use `endtime(a)`, `duration(a)`, and `nsamples(a)` helpers rather than
  computing from fields directly.
- Calibration is applied later via `apply_calibration!`, never at I/O.
"""
struct Audiodata
    sig::Vector{Float64}
    fs::Float32
    starttime::DateTime
    is_calibrated::Bool
    calibration::Calibration
    metadata::RecordingMetadata

    function Audiodata(sig::AbstractVector{<:Real},
                       fs::Real,
                       starttime::DateTime;
                       is_calibrated::Bool = false,
                       calibration::Calibration = NoCalibration(),
                       timezone::Union{String,Nothing} = nothing,
                       lat::Union{Float64,Missing} = missing,
                       lon::Union{Float64,Missing} = missing,
                       site_id::Union{String,Nothing} = nothing,
                       recorder::AbstractString = "unknown",
                       recorder_id::Union{String,Nothing} = nothing)

        fs <= 0 && throw(ArgumentError("Audiodata: fs must be > 0, got $fs."))

        meta = RecordingMetadata(timezone, lat, lon, site_id,
                                 String(recorder), recorder_id)
        new(Float64.(sig), Float32(fs), starttime, is_calibrated,
            calibration, meta)
    end
end

"""
    endtime(a::Audiodata) -> DateTime

Purpose:     Compute the recording end time from `starttime`, sample count,
             and sampling rate.

Arguments:
- `a::Audiodata`: The recording.

Returns:     `DateTime`. Computed using nanosecond arithmetic to minimise
             rounding error, then stored as a `DateTime` (millisecond
             resolution). Sub-millisecond remainders are truncated.

Constraints: Returns `starttime` unchanged for an empty recording
             (`nsamples = 0`). See `docs/src/sources.md` for why
             `read_audio_range` avoids round-tripping through `endtime`
             when slicing at the file boundary.

Fails when:  Never.

Example:
```julia
endtime(a)   # e.g. DateTime("2023-10-12T00:07:54.000")
```
"""
endtime(a::Audiodata) =
    a.starttime + Nanosecond(round(Int, 1e9 * length(a.sig) / a.fs))

"""
    duration(a::Audiodata) -> Millisecond

Purpose:     Compute the duration of the recording.

Arguments:
- `a::Audiodata`: The recording.

Returns:     `Millisecond`. Use `.value` to extract the integer count if
             needed.

Constraints: Returns `Millisecond(0)` for an empty recording.

Fails when:  Never.

Example:
```julia
duration(a)          # e.g. Millisecond(60000) for a 60-second recording
duration(a).value    # 60000 (integer milliseconds)
```
"""
duration(a::Audiodata) = endtime(a) - a.starttime

"""
    nsamples(a::Audiodata) -> Int

Purpose:     Return the number of audio samples in the recording.

Arguments:
- `a::Audiodata`: The recording.

Returns:     Non-negative `Int`. Returns `0` for an empty recording — the
             sentinel value returned by `read_audio_range` when a requested
             window falls entirely in a gap.

Fails when:  Never.

Example:
```julia
nsamples(a)   # e.g. 2_880_000 for 60 s at 48 000 Hz
```
"""
nsamples(a::Audiodata) = length(a.sig)

"""
    apply_calibration(a::Audiodata; plans=nothing) -> Audiodata

Purpose:     Apply the calibration stored in `a.calibration` to the recording's
             signal and return a new `Audiodata` with `is_calibrated = true`.
             The in-place method `apply_calibration!` is used internally; this
             wrapper handles allocation and struct reconstruction.

Arguments:
- `a::Audiodata`: The recording to calibrate. Must have `is_calibrated = false`.
- `plans::Union{Nothing, Tuple} = nothing`: Pre-computed FFTW plan tuple
  `(forward_plan, inverse_plan)` for the `TFCalibration` time-domain path.
  Ignored for `ScalarCalibration` and `NoCalibration`. See `apply_calibration!`
  for details on creating and re-using plans across chunks.

Returns:     A new `Audiodata` with:
             - `sig` containing the calibrated signal.
             - `is_calibrated = true`.
             - All other fields (`fs`, `starttime`, `calibration`, `metadata`)
               copied unchanged from `a`.

             Special case: if `a.calibration isa NoCalibration`, returns `a`
             unchanged with `is_calibrated` still `false` — no transform is
             possible without calibration data.

Constraints:
- `a.is_calibrated` must be `false`. Double-application throws.
- `a.fs` is passed automatically to the TFCalibration method; no caller action needed.

Fails when:  `a.is_calibrated` is already `true`.

Example:
```julia
a_cal = apply_calibration(a)
@assert a_cal.is_calibrated
@assert !a.is_calibrated   # original unchanged
```

Do not use when: `a.calibration isa NoCalibration` and calibrated output is
required — no transformation is possible. Check `a.calibration` first.
"""
function apply_calibration(a::Audiodata; plans::Union{Nothing, Tuple} = nothing)
    a.is_calibrated && throw(ArgumentError(
        "apply_calibration: recording is already calibrated. " *
        "Double-application of calibration is not permitted."))

    # NoCalibration: no transform is possible; return unchanged, is_calibrated stays false.
    a.calibration isa NoCalibration && return a

    out = similar(a.sig)
    if a.calibration isa TFCalibration
        apply_calibration!(out, a.sig, a.calibration; fs = a.fs, plans = plans)
    else
        apply_calibration!(out, a.sig, a.calibration)
    end

    return Audiodata(out, a.fs, a.starttime;
                     is_calibrated = true,
                     calibration   = a.calibration,
                     timezone      = a.metadata.timezone,
                     lat           = a.metadata.lat,
                     lon           = a.metadata.lon,
                     site_id       = a.metadata.site_id,
                     recorder      = a.metadata.recorder,
                     recorder_id   = a.metadata.recorder_id)
end
