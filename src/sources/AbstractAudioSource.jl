"""
    AbstractAudioSource

Abstract base type for all audio data sources.

A source is a descriptor that knows *where* audio data lives and how to
retrieve it by time. Concrete subtypes implement two required methods and
optionally override `coverage_fraction`.

Required methods
----------------
Every subtype must implement:

- [`time_range(source)`](@ref) — the start and end time of the source.
- [`read_audio_range(source, start, stop)`](@ref) — return audio for the
  requested time window.

Provided methods
----------------
A default implementation is provided for:

- [`coverage_fraction(source, start, stop)`](@ref) — returns `1.0` (assumes
  no gaps). Subtypes with gap information should override this.

Gap handling
------------
`read_audio_range` accepts a `gap_handling` keyword that controls behaviour
when the requested range extends beyond or between available audio:

- `:zero_fill` (default) — pad missing regions with zeros; returned signal
  covers exactly `[start, stop]`.
- `:error` — throw if coverage is less than 1.0.

See [`SingleFileSource`](@ref) for a concrete example.

Zarr readiness
--------------
All chunking and time-based access in the package flows through this
interface. Adding a new storage backend (e.g. Zarr, S3) means implementing
a new subtype here — nothing above this layer changes.
"""
abstract type AbstractAudioSource end

"""
    time_range(source::AbstractAudioSource) -> (DateTime, DateTime)

Purpose:     Return the start and end time of the audio source as a 2-tuple.

Arguments:
- `source::AbstractAudioSource`: Any concrete audio source.

Returns:     `(starttime::DateTime, endtime::DateTime)` covering all audio
             available from this source.

Constraints: `starttime ≤ endtime`. An empty source may return equal times.

Fails when:  Called on a subtype that has not implemented this method.

Example:
```julia
t_start, t_stop = time_range(source)
```
"""
function time_range(source::AbstractAudioSource)
    error("time_range not implemented for $(typeof(source)). " *
          "Every AbstractAudioSource subtype must define time_range(source).")
end

"""
    read_audio_range(source, start, stop; gap_handling=:zero_fill) -> Audiodata

Purpose:     Return an `Audiodata` covering exactly the time window
             `[start, stop]`, reading from `source`.

Arguments:
- `source::AbstractAudioSource`: The audio source to read from.
- `start::DateTime`: Start of the requested window (inclusive).
- `stop::DateTime`: End of the requested window (exclusive).
- `gap_handling::Symbol = :zero_fill`: How to handle regions where no audio
  is available. `:zero_fill` pads with zeros; `:error` raises if coverage
  is less than 1.0.

Returns:     `Audiodata` whose signal covers `[start, stop]`. When the
             requested window falls entirely outside the source, returns an
             empty `Audiodata` (zero-length signal). Check
             `coverage_fraction(source, start, stop)` to determine how much
             of the returned signal is real audio vs. zero-fill.

Constraints: `start < stop`. The returned `Audiodata` always has
             `is_calibrated = false`; calibration is applied separately.

Fails when:
- Called on a subtype that has not implemented this method.
- `gap_handling = :error` and coverage is less than 1.0.

Example:
```julia
using Dates
chunk = read_audio_range(source, DateTime(2023,10,12,0,0,0),
                                 DateTime(2023,10,12,0,1,0))
```

Do not use when: You need calibrated audio — call `apply_calibration!` on
                 the returned `Audiodata` before metric computation.
"""
function read_audio_range(source::AbstractAudioSource,
                          start::DateTime,
                          stop::DateTime;
                          gap_handling::Symbol = :zero_fill)
    error("read_audio_range not implemented for $(typeof(source)). " *
          "Every AbstractAudioSource subtype must define " *
          "read_audio_range(source, start, stop; gap_handling).")
end

"""
    coverage_fraction(source, start, stop) -> Float64

Purpose:     Return the fraction of the time window `[start, stop]` that
             contains real audio data (as opposed to silence or gaps).

Arguments:
- `source::AbstractAudioSource`: The audio source to query.
- `start::DateTime`: Start of the window.
- `stop::DateTime`: End of the window.

Returns:     `Float64` in `[0.0, 1.0]`. `1.0` means the window is fully
             covered by audio. `0.0` means no audio exists in the window.

Constraints: The default implementation always returns `1.0` — it does not
             validate `start < stop` or inspect the source at all. Subtypes
             that have gap information — such as `IndexedFileSource` — should
             override this method and enforce `start < stop`.

Example:
```julia
frac = coverage_fraction(source, DateTime(2023,10,12), DateTime(2023,10,13))
```
"""
function coverage_fraction(source::AbstractAudioSource,
                           start::DateTime,
                           stop::DateTime)
    return 1.0
end
