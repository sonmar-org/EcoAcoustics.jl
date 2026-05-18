# index_builder.jl is always included before this file, so Arrow, Logging,
# _NoCalibrationWarnLogger, _EXPECTED_INDEX_COLUMNS, and load_index are
# already in scope from that file.

# Warn when zero-fill would produce more than this many seconds of silence.
# Five minutes is the threshold: gaps shorter than this are routine duty-cycle
# pauses; gaps longer than this are likely to surprise a caller expecting
# mostly real audio.
const _GAP_WARN_THRESHOLD_S = 300.0

# ─── Internal: schema validation ─────────────────────────────────────────────

function _validate_index_schema(tbl::Arrow.Table, context::String)
    present      = Set(propertynames(tbl))
    missing_cols = setdiff(_EXPECTED_INDEX_COLUMNS, present)
    isempty(missing_cols) || error(
        "$context: Arrow.Table is missing expected index columns: " *
        join(sort(string.(missing_cols)), ", ") * ".")
end

# ─── Struct ───────────────────────────────────────────────────────────────────

"""
    IndexedFileSource <: AbstractAudioSource

An audio source backed by a timestamp index over a directory of WAV and FLAC
files.

Where `SingleFileSource` loads one file into memory at construction time,
`IndexedFileSource` works from a pre-built index: a lightweight table (one row
per file) that records each file's path, timestamps, and sample count. At query
time, only the files that overlap the requested window are read from disk.

This is the standard source for archive-scale work. Build the index once with
`build_index`, save it with `output_path`, and construct `IndexedFileSource`
from the saved file on every subsequent run.

Fields
------
- `index::Arrow.Table`: The loaded timestamp index. One row per audio file,
  sorted ascending by `start_time`.
- `root::String`: Directory against which relative `file_path` values in the
  index are resolved when reading files.

Do not use when: You are working with a single file — use `SingleFileSource`
                 instead. Do not construct directly; use one of the two
                 constructor methods below.
"""
struct IndexedFileSource <: AbstractAudioSource
    index::Arrow.Table
    root::String
end

# ─── Constructors ────────────────────────────────────────────────────────────

"""
    IndexedFileSource(arrow_path; root=dirname(abspath(arrow_path))) -> IndexedFileSource

Purpose:     Construct an `IndexedFileSource` by loading a previously saved
             Arrow index file. This is the standard constructor for production
             use.

Arguments:
- `arrow_path::AbstractString`: Path to an Arrow IPC file written by
  `build_index`. Relative paths are resolved against the current directory.
- `root::AbstractString`: Directory that file paths in the index are relative
  to. Defaults to the directory containing the Arrow file, which is correct
  when the index was saved alongside the audio data.

Returns:     `IndexedFileSource` ready to serve time-range queries.

Constraints: The Arrow file must contain all eleven columns of the EcoAcoustics
             index schema (validated by `load_index`). The index must have at
             least one row. `root` should point to the directory that was passed
             to `build_index` when the index was created.

Fails when:
- `arrow_path` does not exist or is not a file.
- The Arrow file is missing required columns.
- The index contains zero rows.

Example:
```julia
src = IndexedFileSource("/data/deploy1/index.arrow")
t1, t2 = time_range(src)
```
"""
function IndexedFileSource(arrow_path::AbstractString;
                           root::AbstractString = dirname(abspath(arrow_path)))
    isfile(arrow_path) || throw(ArgumentError(
        "IndexedFileSource: file not found: \"$arrow_path\"."))
    # load_index validates the schema, so no further check is needed here.
    tbl = load_index(arrow_path)
    length(tbl.start_time) == 0 && throw(ArgumentError(
        "IndexedFileSource: index at \"$arrow_path\" contains zero rows. " *
        "Run build_index on a directory containing audio files."))
    IndexedFileSource(tbl, String(root))
end

"""
    IndexedFileSource(tbl::Arrow.Table; root) -> IndexedFileSource

Purpose:     Construct an `IndexedFileSource` from an in-memory `Arrow.Table`,
             typically the value returned directly by `build_index` when no
             `output_path` was given.

Arguments:
- `tbl::Arrow.Table`: A table with the EcoAcoustics index schema. Schema is
  validated before construction.
- `root::AbstractString`: Directory that file paths in the table are relative
  to. Required; there is no sensible default when no file path is available.

Returns:     `IndexedFileSource` ready to serve time-range queries.

Constraints: `tbl` must contain all eleven expected columns. The table must
             have at least one row.

Fails when:
- `tbl` is missing required columns.
- `tbl` contains zero rows.

Example:
```julia
tbl = build_index("/data/deploy1"; recorder="sm3m")
src = IndexedFileSource(tbl; root="/data/deploy1")
```
"""
function IndexedFileSource(tbl::Arrow.Table; root::AbstractString)
    _validate_index_schema(tbl, "IndexedFileSource")
    length(tbl.start_time) == 0 && throw(ArgumentError(
        "IndexedFileSource: table contains zero rows. " *
        "Run build_index on a directory containing audio files."))
    IndexedFileSource(tbl, String(root))
end

# ─── time_range ───────────────────────────────────────────────────────────────

"""
    time_range(source::IndexedFileSource) -> (DateTime, DateTime)

Purpose:     Return the earliest start time and latest end time across all
             files in the index.

Arguments:
- `source::IndexedFileSource`: The source to query.

Returns:     `(first_start, last_end)` as a 2-tuple of `DateTime`. This spans
             the entire deployment, including any gaps between files.

Constraints: The index must contain at least one row (enforced at construction).

Fails when:  Never — both values are scanned from the index at call time.

Example:
```julia
t_start, t_stop = time_range(src)
```
"""
function time_range(source::IndexedFileSource)
    return (minimum(source.index.start_time), maximum(source.index.end_time))
end

# ─── coverage_fraction ────────────────────────────────────────────────────────

"""
    coverage_fraction(source::IndexedFileSource, start, stop) -> Float64

Purpose:     Return the fraction of the window `[start, stop]` that is covered
             by audio files in the index. Gaps between files count as zero
             coverage.

Arguments:
- `source::IndexedFileSource`: The source to query.
- `start::DateTime`, `stop::DateTime`: The time window of interest.

Returns:     `Float64` in `[0.0, 1.0]`. `1.0` means all of `[start, stop]` is
             covered by audio. `0.0` means no files overlap the window.

Constraints: Assumes files in the index do not overlap each other — a properly
             built index always satisfies this. If files do overlap, coverage
             may be over-counted. Returns `0.0` for degenerate windows
             (`stop <= start`).

Fails when:  Never.

Example:
```julia
frac = coverage_fraction(src, DateTime(2023,1,1), DateTime(2023,1,2))
```
"""
function coverage_fraction(source::IndexedFileSource,
                           start::DateTime,
                           stop::DateTime)
    stop <= start && return 0.0
    window_ms  = (stop - start).value  # Millisecond count as Int64
    covered_ms = 0
    idx = source.index
    for i in 1:length(idx.start_time)
        overlap_start = max(start, idx.start_time[i])
        overlap_stop  = min(stop,  idx.end_time[i])
        if overlap_stop > overlap_start
            covered_ms += (overlap_stop - overlap_start).value
        end
    end
    return covered_ms / window_ms
end

# ─── read_audio_range ─────────────────────────────────────────────────────────

"""
    read_audio_range(source::IndexedFileSource, start, stop; gap_handling=:zero_fill)

Purpose:     Return an `Audiodata` covering exactly `[start, stop]`, reading
             from whichever files in the index overlap the requested window.
             Files are read from disk only as needed; files outside the window
             are never touched.

Arguments:
- `source::IndexedFileSource`: The source to read from.
- `start::DateTime`: Start of the requested window (inclusive).
- `stop::DateTime`: End of the requested window (exclusive).
- `gap_handling::Symbol = :zero_fill`: How to handle gaps between files:
  - `:zero_fill` — gap regions are filled with zeros in the output signal.
    `coverage_fraction` on the returned `Audiodata` reflects the fill.
  - `:error` — raise if coverage is less than 1.0.

Returns:     `Audiodata` whose signal covers exactly `[start, stop]`, with
             `starttime = start`. When the window falls entirely in a gap,
             returns an empty `Audiodata` (zero-length signal). Metadata
             (`recorder`, `site_id`, calibration, etc.) comes from the first
             file that contributes samples to the output — if the window spans
             files from different deployments, the first file's metadata is used
             for the whole result. This simplification is noted in Constraints.

Constraints: `start < stop` is required. Only `:zero_fill` and `:error` are
             valid for `gap_handling` — `:skip` raises with an explanation (use
             `chunks()` for gap-skipping iteration). All overlapping files must
             have the same sample rate. The returned `Audiodata` always has
             `is_calibrated = false`. Metadata is from the first contributing
             file — if the window spans files with different recorders or
             site_ids, the first file's metadata applies to the whole result.

Fails when:
- `start >= stop`.
- `gap_handling` is not `:zero_fill` or `:error`.
- Overlapping files have different sample rates.
- `gap_handling = :error` and coverage is less than 1.0.

Example:
```julia
using Dates
chunk = read_audio_range(src,
                         DateTime(2023,1,1,0,0,0),
                         DateTime(2023,1,1,0,1,0))
```

Do not use when: You need calibrated audio — call `apply_calibration!` on the
                 returned `Audiodata` before metric computation.
"""
function read_audio_range(source::IndexedFileSource,
                          start::DateTime,
                          stop::DateTime;
                          gap_handling::Symbol = :zero_fill)

    start >= stop && throw(ArgumentError(
        "read_audio_range: start ($start) must be before stop ($stop)."))

    gap_handling in (:zero_fill, :error) || throw(ArgumentError(
        "read_audio_range: gap_handling must be :zero_fill or :error, " *
        "got :$gap_handling. Use chunks() for gap-skipping iteration " *
        "over an IndexedFileSource."))

    idx = source.index

    # ── Find overlapping row indices ──────────────────────────────────────────
    # A row overlaps [start, stop) when its interval [start_time, end_time)
    # intersects: row.start_time < stop AND row.end_time > start.
    # The index is already sorted by start_time (guaranteed by build_index).
    overlap_indices = Int[]
    for i in 1:length(idx.start_time)
        idx.start_time[i] < stop && idx.end_time[i] > start || continue
        push!(overlap_indices, i)
    end

    # ── Empty result (window falls entirely in a gap) ─────────────────────────
    # Use the first row's fs and recorder as fallbacks for the Audiodata
    # constructor — the signal is length-0, so these values are placeholders.
    if isempty(overlap_indices)
        return Audiodata(Float64[], idx.fs[1], start;
                         recorder = String(idx.recorder[1]))
    end

    # ── Sample rate consistency check ─────────────────────────────────────────
    # All files contributing to the output must share a sample rate. Mixing
    # rates would require resampling, which is out of scope for v1.
    fs_first = idx.fs[overlap_indices[1]]
    for i in overlap_indices
        idx.fs[i] == fs_first || throw(ArgumentError(
            "read_audio_range: overlapping files have different sample rates " *
            "($(idx.fs[i]) Hz vs $fs_first Hz). Mixed-fs archives are not " *
            "supported in v1."))
    end
    fs = fs_first

    # ── Coverage check (used for :error and gap warning) ─────────────────────
    # Compute once here — before any file I/O — so :error fails fast and the
    # gap warning fires before the user waits for a large allocation.
    frac       = coverage_fraction(source, start, stop)
    window_s   = (stop - start).value / 1000.0
    gap_s      = (1.0 - frac) * window_s

    if gap_handling == :error && frac < 1.0
        throw(ArgumentError(
            "read_audio_range: window [$start, $stop] has coverage " *
            "$(round(100 * frac; digits=1))%. " *
            "Pass gap_handling=:zero_fill to allow gaps."))
    end

    if gap_handling == :zero_fill && gap_s > _GAP_WARN_THRESHOLD_S
        @warn "read_audio_range: $(round(Int, gap_s / 60)) minutes of gaps " *
              "will be zero-filled. Use chunks() with gap_handling=:skip to " *
              "iterate over recorded audio only." start=start stop=stop
    end

    # ── Allocate output buffer ────────────────────────────────────────────────
    # Pre-filled with zeros: any gap between files is already handled — no
    # explicit gap-filling code is needed.
    n_out   = round(Int, (stop - start).value / 1000.0 * fs)
    sig_out = zeros(Float64, n_out)

    first_meta = nothing
    first_cal  = nothing
    suppress   = _NoCalibrationWarnLogger(Logging.current_logger())

    # ── Copy each file's contribution into sig_out ────────────────────────────
    for i in overlap_indices
        full_path = joinpath(source.root, String(idx.file_path[i]))

        # Read the full file. Pass starttime from the index — this uses the
        # authoritative indexed timestamp rather than re-parsing the filename.
        audio = Logging.with_logger(suppress) do
            read_audio(full_path;
                       recorder  = String(idx.recorder[i]),
                       starttime = idx.start_time[i])
        end

        n_file = length(audio.sig)

        # The portion of this file that falls within [start, stop).
        overlap_start = max(idx.start_time[i], start)
        overlap_stop  = min(idx.end_time[i],   stop)

        # out_offset: how many samples from the start of sig_out before this
        # file's contribution begins.
        out_offset = round(Int, (overlap_start - start).value / 1000.0 * fs)

        # src_offset: how many samples into the file's signal before the
        # overlapping region begins.
        src_offset = round(Int, (overlap_start - idx.start_time[i]).value / 1000.0 * fs)

        # src_stop: the sample index within the file where the overlap ends.
        # When the overlap extends to the file's last sample, use n_file
        # directly rather than a time-derived value. This avoids the DateTime
        # millisecond-truncation off-by-one that SingleFileSource also guards
        # against: endtime() uses nanosecond arithmetic, but DateTime stores
        # only milliseconds, so the round-trip can lose 1 sample.
        src_stop =
            overlap_stop >= idx.end_time[i] ?
            n_file :
            round(Int, (overlap_stop - idx.start_time[i]).value / 1000.0 * fs)

        # clamp ensures we never read past the file's signal or write past
        # sig_out, even if rounding produces an off-by-one.
        n_copy = clamp(src_stop - src_offset,
                       0,
                       min(n_file - src_offset, n_out - out_offset))

        if n_copy > 0
            sig_out[out_offset + 1 : out_offset + n_copy] .=
                audio.sig[src_offset + 1 : src_offset + n_copy]
        end

        # Capture metadata from the first overlapping file.
        if first_meta === nothing
            first_meta = audio.metadata
            first_cal  = audio.calibration
        end
    end

    meta = first_meta
    return Audiodata(sig_out, fs, start;
                     calibration = first_cal,
                     timezone    = meta.timezone,
                     lat         = meta.lat,
                     lon         = meta.lon,
                     site_id     = meta.site_id,
                     recorder    = meta.recorder,
                     recorder_id = meta.recorder_id)
end
