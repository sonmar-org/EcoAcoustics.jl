"""
    SingleFileSource <: AbstractAudioSource

An audio source backed by a single WAV or FLAC file.

The entire file is read and cached in memory at construction time. All
subsequent `read_audio_range` calls slice the in-memory signal — no further
disk I/O occurs. This is appropriate for single-file, one-off analyses where
the file fits comfortably in RAM. For archive-scale work spanning many files,
use `IndexedFileSource` instead.

Fields
------
- `path::String`: Path to the audio file on disk.
- `audio::Audiodata`: Full recording, loaded at construction.

Do not use when: The file is too large to hold in RAM, or you are processing
                 a directory of files — use `IndexedFileSource` for those cases.
"""
struct SingleFileSource <: AbstractAudioSource
    path::String
    audio::Audiodata
end

"""
    SingleFileSource(path; kwargs...) -> SingleFileSource

Purpose:     Construct a `SingleFileSource` by reading the audio file at
             `path` into memory. All keyword arguments are passed through to
             `read_audio`.

Arguments:
- `path::AbstractString`: Path to a WAV or FLAC file.
- `kwargs...`: Any keyword accepted by `read_audio` — `recorder`,
  `starttime`, `lat`, `lon`, `site_id`, `strict`.

Returns:     A `SingleFileSource` with the full signal cached in `audio`.

Constraints: The file must be readable and mono. The entire signal is loaded
             into memory; files larger than available RAM will cause an
             out-of-memory error.

Fails when:  The file cannot be read, is not mono, or `read_audio` throws.

Example:
```julia
src = SingleFileSource("T1-C__0__20170912_181500.wav"; recorder="sm3m")
t1, t2 = time_range(src)
```
"""
function SingleFileSource(path::AbstractString; kwargs...)
    audio = read_audio(path; kwargs...)
    SingleFileSource(String(path), audio)
end

"""
    time_range(source::SingleFileSource) -> (DateTime, DateTime)

Purpose:     Return the start and end time of the audio file.

Arguments:
- `source::SingleFileSource`: The source to query.

Returns:     `(starttime, endtime)` as a 2-tuple of `DateTime`. `endtime` is
             derived from `starttime`, sample count, and sampling rate —
             it is not read from the file header.

Constraints: `endtime` is computed via nanosecond arithmetic but stored in a
             `DateTime`, which has millisecond resolution. The sub-millisecond
             remainder is truncated. See `read_audio_range` for how this is
             handled correctly when slicing.

Fails when:  Never — the values are cached at construction time.

Example:
```julia
t_start, t_stop = time_range(src)
```
"""
time_range(source::SingleFileSource) = (source.audio.starttime, endtime(source.audio))

"""
    coverage_fraction(source::SingleFileSource, start, stop) -> Float64

Purpose:     Return the fraction of `[start, stop]` covered by this file's
             audio. Used by `read_audio_range` and the future `chunks` iterator
             to annotate each chunk with its data quality.

Arguments:
- `source::SingleFileSource`: The source to query.
- `start::DateTime`, `stop::DateTime`: The time window of interest.

Returns:     `Float64` in `[0.0, 1.0]`. Computed as the millisecond length of
             the intersection of `[start, stop]` with the file's time range,
             divided by the millisecond length of `[start, stop]`.

Constraints: `start` must be before `stop`. Returns `0.0` for degenerate
             windows (`start >= stop`).

Fails when:  Never throws; degenerate windows return `0.0`.

Example:
```julia
frac = coverage_fraction(src, DateTime(2017,9,12,18,15), DateTime(2017,9,12,18,16))
```
"""
function coverage_fraction(source::SingleFileSource,
                           start::DateTime,
                           stop::DateTime)
    file_start, file_end = time_range(source)
    # DateTime subtraction returns Millisecond; .value extracts the integer ms count.
    requested_ms = (stop - start).value
    requested_ms <= 0 && return 0.0
    overlap_start = max(start, file_start)
    overlap_stop  = min(stop,  file_end)
    overlap_ms    = max(0, (overlap_stop - overlap_start).value)
    return overlap_ms / requested_ms
end

"""
    read_audio_range(source::SingleFileSource, start, stop; gap_handling=:zero_fill)

Purpose:     Return an `Audiodata` whose signal covers exactly `[start, stop]`,
             sliced from the cached file signal and zero-padded where the
             requested window extends beyond the file.

Arguments:
- `source::SingleFileSource`: The source to read from.
- `start::DateTime`: Start of the requested window (inclusive).
- `stop::DateTime`: End of the requested window (exclusive).
- `gap_handling::Symbol = :zero_fill`: Behaviour when the window extends
  beyond the file:
  - `:zero_fill` — pad out-of-file regions with zeros.
  - `:error` — throw if `coverage_fraction < 1.0`.

Returns:     `Audiodata` with `starttime = start`. Signal length is
             `leading_zeros + real_samples + trailing_zeros`, where each
             piece is computed independently from its own time span (see
             Constraints). Returns an empty `Audiodata` (zero-length signal)
             when the window has no overlap with the file at all.

Constraints: `start < stop` is required.

             **DateTime precision and sample counts:** `DateTime` has
             millisecond resolution. The file's true duration in samples
             rarely divides to an exact millisecond count, so converting
             the total window duration to samples and back introduces a
             rounding error of up to 1 sample. To avoid this, boundaries
             that coincide with the file edge use the exact sample count
             directly rather than deriving it from time arithmetic. Only
             gap regions (zero-fill outside the file) go through
             millisecond rounding, where a ±1 sample error is
             inconsequential.

             The returned `Audiodata` always has `is_calibrated = false`.
             Calibration and metadata are copied from the source file.

Fails when:
- `start >= stop`.
- `gap_handling = :error` and `coverage_fraction < 1.0`.
- `gap_handling` is not `:zero_fill` or `:error`.

Example:
```julia
using Dates
chunk = read_audio_range(src,
                         DateTime(2017,9,12,18,15,0),
                         DateTime(2017,9,12,18,16,0))
```

Do not use when: You need calibrated audio — call `apply_calibration!` on
                 the returned `Audiodata` before metric computation.
"""
function read_audio_range(source::SingleFileSource,
                          start::DateTime,
                          stop::DateTime;
                          gap_handling::Symbol = :zero_fill)

    start >= stop && throw(ArgumentError(
        "read_audio_range: start ($start) must be before stop ($stop)."))

    gap_handling in (:zero_fill, :error) || throw(ArgumentError(
        "read_audio_range: gap_handling must be :zero_fill or :error, " *
        "got :$gap_handling."))

    frac = coverage_fraction(source, start, stop)

    if gap_handling == :error && frac < 1.0
        throw(ArgumentError(
            "read_audio_range: window [$start, $stop] is only " *
            "$(round(100 * frac; digits=1))% covered by \"$(source.path)\". " *
            "Pass gap_handling=:zero_fill to allow partial coverage."))
    end

    # Collect metadata fields into a named tuple so every return path can
    # spread them into the Audiodata constructor with `meta_kwargs...`.
    # This avoids repeating all six fields at each return site.
    meta = source.audio.metadata
    meta_kwargs = (calibration = source.audio.calibration,
                   timezone    = meta.timezone,
                   lat         = meta.lat,
                   lon         = meta.lon,
                   site_id     = meta.site_id,
                   recorder    = meta.recorder,
                   recorder_id = meta.recorder_id)

    # No overlap — return empty Audiodata as the documented sentinel.
    frac == 0.0 && return Audiodata(Float64[], source.audio.fs, start; meta_kwargs...)

    file_start, file_end = time_range(source)
    fs     = source.audio.fs
    n_file = length(source.audio.sig)

    # --- Leading zeros and source start index ---
    # When the window starts before the file, zero-fill the gap and begin
    # reading from sample 0. When the window starts inside the file, compute
    # the start sample from the time offset. The src_offset = 0 branch avoids
    # any rounding — the file always starts at sample 0.
    if start < file_start
        leading    = round(Int, (file_start - start).value / 1000.0 * fs)
        src_offset = 0
    else
        leading    = 0
        src_offset = round(Int, (start - file_start).value / 1000.0 * fs)
    end

    # --- Trailing zeros and source stop index ---
    # When the window ends at or after the file end, use the exact sample count
    # (n_file) rather than a time-derived estimate. This is the key fix for the
    # DateTime millisecond-truncation issue: endtime() uses nanosecond arithmetic
    # but DateTime stores only milliseconds, so the round-trip can lose 1 sample.
    # Using n_file directly is always correct.
    if stop >= file_end
        trailing = round(Int, (stop - file_end).value / 1000.0 * fs)
        src_stop = n_file
    else
        trailing = 0
        src_stop = round(Int, (stop - file_start).value / 1000.0 * fs)
    end

    n_real = clamp(src_stop - src_offset, 0, n_file - src_offset)
    n_out  = leading + n_real + trailing

    n_out == 0 && return Audiodata(Float64[], fs, start; meta_kwargs...)

    sig_out = zeros(Float64, n_out)
    sig_out[leading + 1 : leading + n_real] .=
        source.audio.sig[src_offset + 1 : src_offset + n_real]

    return Audiodata(sig_out, fs, start; meta_kwargs...)
end
