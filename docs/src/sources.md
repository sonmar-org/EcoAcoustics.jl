# Audio sources and time-based access

EcoAcoustics.jl separates the question of *where audio lives* from the
question of *how it is processed*. The bridge between the two is the
`AbstractAudioSource` interface.

---

## The source abstraction

An **audio source** is a Julia struct that knows how to retrieve audio for
any requested time window. Callers ask for audio by time:

```julia
chunk = read_audio_range(source, DateTime(2023,10,12,18,15),
                                 DateTime(2023,10,12,18,16))
```

and the source figures out which file(s) to read, how to assemble the
signal, and how to fill any gaps. Nothing above this layer opens files
directly.

`AbstractAudioSource` is defined in `src/sources/AbstractAudioSource.jl`.
Any concrete subtype must implement two methods:

| Method | Returns |
|--------|---------|
| `time_range(source)` | `(start::DateTime, stop::DateTime)` — the full span of available audio |
| `read_audio_range(source, start, stop)` | `Audiodata` covering `[start, stop]` exactly |

A third method, `coverage_fraction(source, start, stop)`, has a default
implementation (always returns `1.0`) that subtypes can override when they
have gap information.

### Why an abstract interface?

The same chunking code, metric computations, and analysis pipelines work
identically regardless of whether audio comes from one WAV file, an
indexed directory of thousands of files, or (in a future version) a cloud
Zarr archive. The source handles the I/O details; everything else sees only
`read_audio_range`.

Adding a new storage backend — for example, S3-hosted Zarr — means writing
one new subtype in `src/sources/`. No code above that layer changes.

---

## SingleFileSource

`SingleFileSource` is the simplest concrete source: it wraps one WAV or
FLAC file. At construction time, the entire file is read into memory and
cached in an `Audiodata` struct. All subsequent `read_audio_range` calls
slice the in-memory signal — no further disk I/O occurs.

```julia
src = SingleFileSource("T1-C__0__20170912_181500.wav"; recorder="sm3m")

t_start, t_stop = time_range(src)
chunk = read_audio_range(src, t_start, t_stop + Minute(1))
```

**When to use it:** For one-off analysis of a single file, unit tests, and
interactive exploration. The file must fit in RAM.

**When not to use it:** When processing an archive of many files — use
`IndexedFileSource` instead.

---

## Gap handling

`read_audio_range` accepts a `gap_handling` keyword that controls what
happens when the requested time window extends beyond the available audio:

| Value | Behaviour |
|-------|-----------|
| `:zero_fill` (default) | Regions outside the file are filled with zeros. The returned signal always spans exactly `[start, stop]`. |
| `:error` | An `ArgumentError` is thrown if any part of the window is uncovered. |

The fraction of the returned signal that is real audio (as opposed to
zero-fill) is available from `coverage_fraction`:

```julia
frac = coverage_fraction(src, start, stop)
# 1.0 = fully covered; 0.0 = entirely outside the file
```

Metric functions and the `chunks` iterator use `coverage_fraction` to
annotate each result row, so downstream analyses can distinguish high-
quality data from padded windows.

### Return value for no-overlap windows

If the requested window has no overlap with the file at all,
`read_audio_range` returns an `Audiodata` with a zero-length signal
(`nsamples(chunk) == 0`). Downstream code must handle this case. This
avoids the need for `nothing`-checks everywhere — the same type is returned
in all cases.

---

## The DateTime precision constraint

`DateTime` in Julia has **millisecond** resolution. This becomes important
at the boundary between time arithmetic and sample counts.

Consider a file recorded at 96 000 Hz with 10.5 seconds of audio:

```
nsamples = 1 008 000
duration = 1 008 000 / 96 000 = 10.500 000 s  (exact)
```

The `endtime` helper computes this as nanoseconds first and then stores in
a `DateTime`:

```julia
endtime(a) = a.starttime + Nanosecond(round(Int, 1e9 * nsamples / fs))
#          = starttime + Nanosecond(10_500_000_000)
#          = starttime + 10 500 ms   ✓ (exact in this example)
```

Now consider a file where the duration does *not* divide to an exact
millisecond:

```
nsamples = 1 000 001
duration = 1 000 001 / 96 000 ≈ 10.416 677 083 s
         = 10 416.677 083 ms
```

Adding `Nanosecond(10_416_677_083)` to a `DateTime` truncates to milliseconds:

```julia
starttime + Nanosecond(10_416_677_083)
# → starttime + 10 416 ms   (677 µs dropped)
```

If `read_audio_range` uses this truncated `endtime` to derive how many
samples to copy, it gets:

```
samples_to_copy = round(10_416 ms / 1000 * 96_000) = 1 000 000 ≠ 1 000 001
```

Off by one sample at the tail of every file that doesn't divide evenly.

### The fix: split-boundary arithmetic

`SingleFileSource.read_audio_range` avoids this by splitting sample-count
computation at the file boundary:

- When the window **ends at or after the file end** (`stop >= file_end`),
  use the exact `n_file` sample count for the real-audio portion. No time
  arithmetic, no rounding.
- When the window **starts at or before the file start** (`start <= file_start`),
  set `src_offset = 0`. The file always starts at sample 0.
- Only **gap regions** (zero-fill outside the file) are derived from time
  arithmetic, where a ±1 sample error in silence is inconsequential.

This is a deliberate design decision, documented in the source code.

### Click detectors require sample-level precision

The `DateTime` truncation issue is why click detector output timestamps
are stored as a `(file_start::DateTime, sample_offset::Int64)` pair, not
as a bare `DateTime`. The `file_start` anchors the recording, and
`sample_offset` gives the exact integer position within it. The wall-clock
time of a click is then:

```julia
event_time(file_start, sample_offset, fs) =
    file_start + Nanosecond(round(Int, 1e9 * sample_offset / fs))
```

This is lossless. Storing clicks as plain `DateTime` would lose up to
0.5 ms of precision at 96 kHz — about 48 samples. Click detector design
is planned for v2.

---

## Metadata propagation

`read_audio_range` copies calibration and all `RecordingMetadata` fields
from the source file to every returned `Audiodata`. The returned chunk
always has `is_calibrated = false`, regardless of the source file's state.
Calibration is applied later via `apply_calibration!`.

This means the chunk carries the same recorder identity, location, and
calibration object as the original file, ready for metric computation.

---

## Iterating over chunks

`chunks` is the low-level iteration primitive. It returns a lazy
`ChunksIterator` that yields one `Audiodata` per window on demand:

```julia
src = SingleFileSource("deployment.wav"; recorder="sm3m")
for chunk in chunks(src; chunk_seconds = 60.0)
    println(chunk.starttime, "  nsamples=", nsamples(chunk))
end
```

The iterator has `Base.IteratorSize = SizeUnknown` — `length` is not defined
because the number of emitted chunks depends on gap filtering and cannot be
known without iterating. Use `collect` to materialise all chunks, or process
them in a `for` loop.

### Gap handling in `chunks`

The `gap_handling` keyword has three modes:

| Value | Behaviour |
|-------|-----------|
| `:skip` (default) | Windows entirely in a gap (`coverage_fraction == 0`) are silently skipped. Windows with partial coverage are still emitted — their signal is zero-filled at the gap boundary. Matches how ecologists treat duty-cycled recordings: off-period windows are invisible; on-period windows near a boundary carry a `coverage_fraction < 1` annotation. |
| `:zero_fill` | All windows in the source time span are emitted. Gap regions are zero-filled. Use when a uniform time axis is required (e.g., LTSA). |
| `:error` | All windows are emitted. `read_audio_range` throws `ArgumentError` if any window has partial or zero coverage. |

### Overlapping windows

Set `stride_seconds < chunk_seconds` for overlapping windows:

```julia
for chunk in chunks(src; chunk_seconds = 2.0, stride_seconds = 1.0)
    # consecutive chunks share 1 second of signal
end
```

Windows that extend past the source end are still emitted with the
out-of-file region zero-filled (`:zero_fill` / `:skip` modes) or raising
(`:error` mode).

---

## `process_chunks` — batch metric computation

`process_chunks` is the primary entry point for archive-scale computation.
It applies a user function to every chunk and returns a `DataFrame`:

```julia
df = process_chunks(src,
                    chunk -> (rms = sqrt(sum(abs2, chunk.sig) / length(chunk.sig)),);
                    chunk_seconds = 60.0,
                    on_error      = :skip)
# columns: start_time, coverage_fraction, rms
```

The user function `f(chunk::Audiodata)` must return a `NamedTuple`. Its
fields become columns in the output. `start_time` and `coverage_fraction`
are always prepended; do not return those field names from `f`.

### Parallelism

`parallel=:threads` (default) dispatches over all available Julia threads
via `Threads.@threads`. Windows are pre-materialised into an indexed array so
each thread writes its result into a pre-allocated slot — row order in the
output DataFrame matches window order regardless of thread scheduling.

`f` must be thread-safe. Avoid shared mutable state in closures passed as `f`.

`parallel=:none` runs sequentially, which is useful for debugging and for `f`
functions with internal mutable state.

### Error handling

| `on_error` | Behaviour |
|------------|-----------|
| `:fail` (default) | Rethrows any exception immediately. Recommended during development. |
| `:skip` | Drops failed chunks. That window does not appear in the output. Use for long production runs over imperfect archives. |
| `:record` | Adds a row with an `error::String` column. Metric columns are `missing` for that row. |

### Progress

`progress=true` (default) displays a `ProgressMeter` bar on stdout.
The meter is thread-safe.

---

## Planned source types

| Type | Status | Description |
|------|--------|-------------|
| `SingleFileSource` | v1, complete | One WAV/FLAC file in memory |
| `IndexedFileSource` | v1, complete | Arrow timestamp index over a directory |
| `ZarrSource` | v2, deferred | Cloud-native chunked storage |

All three implement the same `AbstractAudioSource` interface. Analysis
pipelines written for `SingleFileSource` run unchanged on an
`IndexedFileSource` covering a multi-terabyte archive.
