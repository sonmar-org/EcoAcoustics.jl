import Logging
import Arrow

# ─── Internal: logger that silences the calibration warning during indexing ──
#
# read_audio calls lookup_calibration on every file. When no CalibrationProfile
# is registered for the recorder, lookup_calibration emits a @warn. That
# warning is correct and useful when reading a single file interactively, but
# becomes noise when indexing thousands of files — every file would warn.
# This logger passes everything through to the outer logger except that one
# specific message.

struct _NoCalibrationWarnLogger{L <: Logging.AbstractLogger} <: Logging.AbstractLogger
    inner::L
end

Logging.min_enabled_level(l::_NoCalibrationWarnLogger) =
    Logging.min_enabled_level(l.inner)

Logging.shouldlog(l::_NoCalibrationWarnLogger, level, mod, group, id) =
    Logging.shouldlog(l.inner, level, mod, group, id)

function Logging.handle_message(l::_NoCalibrationWarnLogger,
                                level, msg, mod, group, id, file, line; kwargs...)
    level == Logging.Warn &&
        startswith(string(msg), NO_CALIBRATION_WARN_PREFIX) && return
    Logging.handle_message(l.inner, level, msg, mod, group, id, file, line; kwargs...)
end

# ─── Supported formats ───────────────────────────────────────────────────────
#
# Only WAV and FLAC are indexed. WAV is read via WAV.jl (pure Julia). FLAC is
# read via FileIO/LibSndFile. Other formats (AIF, AIFF, OGG, ...) are not
# supported — convert them to FLAC or WAV first with a tool such as sox:
#
#   sox input.aif output.flac
#
# FLAC is lossless and smaller than WAV, so it is the recommended conversion
# target for archival data.

const _AUDIO_EXTENSIONS = (".wav", ".flac")

# ─── Schema ──────────────────────────────────────────────────────────────────

const _EXPECTED_INDEX_COLUMNS = Set{Symbol}([
    :file_path, :start_time, :end_time, :fs, :nsamples,
    :recorder, :recorder_id, :site_id,
    :hydrophone_id, :time_uncertainty_ms, :notes,
])

# ─── build_index ─────────────────────────────────────────────────────────────

"""
    build_index(directory; recorder, output_path=nothing) -> Arrow.Table

Purpose:     Walk a directory tree, read every WAV and FLAC file, and build a
             timestamp index as an `Arrow.Table` with one row per file. The
             index is the persistent record of what files exist, when they were
             recorded, and how long each one is. Once built and saved, it is
             loaded by `IndexedFileSource` for all subsequent time-range queries
             — no re-reading of audio files is required. Only `.wav` and `.flac`
             files are indexed; other formats (AIF, AIFF, OGG) are silently
             skipped. Convert them to FLAC or WAV with `sox` before indexing.

Arguments:
- `directory::AbstractString`: Root directory to scan. Walks recursively; all
  `.wav` and `.flac` files at any depth are included. Case-insensitive
  extension matching (`.WAV` is included).
- `recorder::AbstractString`: Recorder family name, e.g. `"sm3m"`,
  `"rockhopper"`. Every file in the directory is assumed to come from this
  recorder type. This drives filename parsing (to extract timestamps) and
  calibration lookup. Mixed-recorder directories are not supported in v1 —
  index each recorder's subdirectory separately.
- `output_path::Union{AbstractString,Nothing} = nothing`: If provided, the
  index is written to this path as an Arrow IPC file. Load it later with
  `load_index`. If `nothing`, the index is returned in memory only.

Returns:     `Arrow.Table` with one row per successfully indexed file, sorted
             ascending by `start_time`. Columns match the EcoAcoustics index
             schema exactly:
             `file_path` (String), `start_time` (DateTime), `end_time`
             (DateTime), `fs` (Float32), `nsamples` (Int64), `recorder`
             (String), `recorder_id` (Union{String,Missing}), `site_id`
             (Union{String,Missing}), `hydrophone_id` (Union{String,Missing}),
             `time_uncertainty_ms` (Union{Float64,Missing}), `notes`
             (Union{String,Missing}).

Constraints: Only WAV and FLAC files are indexed. Files with other extensions
             are silently skipped — this is intentional, not a bug. Convert
             AIF/AIFF or other formats with `sox` before indexing. All files
             must belong to the same recorder type. File paths in
             the returned table are relative to `directory` — `IndexedFileSource`
             resolves them against a root at read time. `end_time` is computed
             from the actual decoded sample count, not the header-claimed value,
             so truncated files record their true duration. `recorder_id` and
             `site_id` are populated only when the recorder's filename format
             encodes them. `hydrophone_id`, `time_uncertainty_ms`, and `notes`
             are always `missing` in v1 (reserved for future use). Calibration
             warnings are suppressed during indexing — calibration is irrelevant
             at index-build time.

Fails when:
- `directory` does not exist or is not a directory.
- Individual unreadable files are skipped with a `@warn` and counted; the
  builder does not abort.

Example:
```julia
tbl = build_index("/data/deploy1"; recorder="sm3m",
                  output_path="/data/deploy1/index.arrow")
tbl.start_time          # access the start_time column directly
load_index("/data/deploy1/index.arrow")  # reload from disk later
```

Do not use when: Files span more than one recorder type — build separate
                 indexes per recorder subdirectory and combine in v2.
"""
function build_index(directory::AbstractString;
                     recorder::AbstractString,
                     output_path::Union{AbstractString,Nothing} = nothing)

    isdir(directory) || throw(ArgumentError(
        "build_index: not a directory: \"$directory\"."))

    # ── Collect file paths ────────────────────────────────────────────────────
    # walkdir yields (root, subdirs, filenames) for every directory in the tree.
    # We join root + filename to get the absolute path for each file.
    audio_files = String[]
    for (root, _, files) in walkdir(directory)
        for file in files
            lowercase(splitext(file)[2]) in _AUDIO_EXTENSIONS || continue
            push!(audio_files, joinpath(root, file))
        end
    end

    n_total   = length(audio_files)
    n_indexed = 0
    n_skipped = 0

    # ── Per-column accumulators ───────────────────────────────────────────────
    # Typed explicitly so Arrow.write produces the correct schema.
    # Union{String,Missing} is required for optional fields — using Any[]
    # would lose type information and produce incorrect Arrow column types.
    col_file_path           = String[]
    col_start_time          = DateTime[]
    col_end_time            = DateTime[]
    col_fs                  = Float32[]
    col_nsamples            = Int64[]
    col_recorder            = String[]
    col_recorder_id         = Union{String,Missing}[]
    col_site_id             = Union{String,Missing}[]
    col_hydrophone_id       = Union{String,Missing}[]
    col_time_uncertainty_ms = Union{Float64,Missing}[]
    col_notes               = Union{String,Missing}[]

    # Build the suppression logger once and reuse it across all files.
    suppress = _NoCalibrationWarnLogger(Logging.current_logger())

    # ── Read each file ────────────────────────────────────────────────────────
    for path in audio_files
        try
            # read_audio gives us the actual decoded signal — length(sig) is
            # the true sample count. The calibration warning is suppressed
            # because calibration is irrelevant at index-build time.
            audio = Logging.with_logger(suppress) do
                read_audio(path; recorder = String(recorder))
            end

            n          = Int64(length(audio.sig))
            fs         = audio.fs
            start_time = audio.starttime
            # end_time arithmetic matches endtime(::Audiodata) exactly.
            end_time   = start_time + Nanosecond(round(Int, 1e9 * n / fs))
            # Store path relative to directory so the index is portable.
            rel_path   = relpath(path, directory)

            push!(col_file_path,           rel_path)
            push!(col_start_time,          start_time)
            push!(col_end_time,            end_time)
            push!(col_fs,                  fs)
            push!(col_nsamples,            n)
            push!(col_recorder,            String(recorder))
            # something(x, missing) returns x when x !== nothing, else missing.
            # This converts the Union{String,Nothing} metadata fields to the
            # Union{String,Missing} type the index schema requires.
            push!(col_recorder_id,         something(audio.metadata.recorder_id, missing))
            push!(col_site_id,             something(audio.metadata.site_id,     missing))
            push!(col_hydrophone_id,       missing)
            push!(col_time_uncertainty_ms, missing)
            push!(col_notes,               missing)

            n_indexed += 1

        catch e
            # Log-and-continue: corrupted or unreadable files become implicit
            # gaps in coverage rather than build failures. The stack trace is
            # attached so the cause is diagnosable.
            @warn "build_index: skipping unreadable file" path=path exception=(e, catch_backtrace())
            n_skipped += 1
        end
    end

    # ── Sort by start_time ────────────────────────────────────────────────────
    # sortperm returns a permutation vector: order[i] is the index of the i-th
    # smallest element. Indexing every column with the same permutation keeps
    # all columns aligned.
    order = sortperm(col_start_time)

    nt = (
        file_path           = col_file_path[order],
        start_time          = col_start_time[order],
        end_time            = col_end_time[order],
        fs                  = col_fs[order],
        nsamples            = col_nsamples[order],
        recorder            = col_recorder[order],
        recorder_id         = col_recorder_id[order],
        site_id             = col_site_id[order],
        hydrophone_id       = col_hydrophone_id[order],
        time_uncertainty_ms = col_time_uncertainty_ms[order],
        notes               = col_notes[order],
    )

    # ── Write Arrow file (optional) ───────────────────────────────────────────
    if output_path !== nothing
        Arrow.write(output_path, nt)
    end

    @info "build_index: complete" total=n_total indexed=n_indexed skipped=n_skipped

    # ── Return Arrow.Table ────────────────────────────────────────────────────
    # If we wrote to disk, load from there so the returned table is backed by
    # the file (memory-mappable for large indexes). Otherwise serialise to an
    # in-memory buffer and load from that.
    if output_path !== nothing
        return Arrow.Table(output_path)
    else
        buf = IOBuffer()
        Arrow.write(buf, nt)
        seekstart(buf)
        return Arrow.Table(buf)
    end
end

# ─── load_index ──────────────────────────────────────────────────────────────

"""
    load_index(path) -> Arrow.Table

Purpose:     Load a previously built index from an Arrow IPC file and return
             it as an `Arrow.Table`. Validates that all expected schema columns
             are present before returning.

Arguments:
- `path::AbstractString`: Path to an Arrow file written by `build_index`.

Returns:     `Arrow.Table` with the same column layout as `build_index`. Access
             columns as `tbl.start_time`, `tbl.file_path`, etc. The table is
             memory-mapped when possible, so loading a large index does not
             require reading the whole file into RAM.

Constraints: The file must have been written by `build_index` or otherwise
             contain all eleven expected columns. Extra columns are ignored.
             Column types are as written; no type coercion is performed.

Fails when:
- `path` does not exist or cannot be read.
- Any of the eleven required columns is absent.

Example:
```julia
tbl = load_index("/data/deploy1/index.arrow")
tbl.start_time
```
"""
function load_index(path::AbstractString)
    tbl = Arrow.Table(path)
    present     = Set(propertynames(tbl))
    missing_cols = setdiff(_EXPECTED_INDEX_COLUMNS, present)
    isempty(missing_cols) || error(
        "load_index: \"$path\" is missing expected columns: " *
        join(sort(string.(missing_cols)), ", ") *
        ". The file may have been built with a different version of EcoAcoustics.")
    return tbl
end
