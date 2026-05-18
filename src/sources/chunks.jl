"""
    ChunksIterator

Lazy iterator returned by `chunks`. Advances through an `AbstractAudioSource`
in fixed-duration windows. Use `chunks` to construct; do not build directly.

`Base.length` is not defined — the number of emitted chunks depends on gap
filtering and cannot be known without iterating. This is declared explicitly
via `Base.IteratorSize(::Type{ChunksIterator}) = Base.SizeUnknown()` so that
generic code (e.g., `collect`) never falls through to a wrong default.
"""
struct ChunksIterator
    source::AbstractAudioSource
    chunk_dur::Millisecond          # window width; rounded from chunk_seconds
    stride_dur::Millisecond         # step between window starts; rounded from stride_seconds
    t_start::DateTime               # earliest possible window start (= time_range(source)[1])
    t_end::DateTime                 # iteration stops when window start >= t_end
    gap_handling::Symbol            # :skip | :zero_fill | :error
end

Base.eltype(::Type{ChunksIterator}) = Audiodata
Base.IteratorSize(::Type{ChunksIterator}) = Base.SizeUnknown()

function Base.iterate(iter::ChunksIterator, t = iter.t_start)
    # In :skip mode, advance past any window whose requested span has zero
    # overlap with the source. Windows with partial coverage are still emitted —
    # their coverage_fraction in the signal reflects the fill.
    if iter.gap_handling == :skip
        while t < iter.t_end &&
              coverage_fraction(iter.source, t, t + iter.chunk_dur) == 0.0
            t += iter.stride_dur
        end
    end

    t >= iter.t_end && return nothing

    t_stop = t + iter.chunk_dur
    # :skip reuses :zero_fill for the actual read: gap filtering was done above,
    # so any remaining partial-coverage boundary window is intentionally included.
    gr    = iter.gap_handling == :skip ? :zero_fill : iter.gap_handling
    chunk = read_audio_range(iter.source, t, t_stop; gap_handling = gr)
    return chunk, t + iter.stride_dur
end

# ─── chunks ──────────────────────────────────────────────────────────────────

"""
    chunks(source; chunk_seconds, stride_seconds=chunk_seconds, gap_handling=:skip)

Purpose:     Return a lazy iterator of `Audiodata` windows covering the time
             span of `source`. Windows are generated on demand — no audio is
             read until the iterator is advanced. This is the low-level
             iteration primitive; `process_chunks` is the preferred entry point
             for batch metric computation.

Arguments:
- `source::AbstractAudioSource` — the audio source to iterate over.
- `chunk_seconds::Real` — duration of each window in seconds. Must be positive.
- `stride_seconds::Real` — step between successive window start times in
  seconds. Default equals `chunk_seconds` (non-overlapping). Values smaller
  than `chunk_seconds` produce overlapping windows.
- `gap_handling::Symbol` — controls behaviour at coverage gaps:
  * `:skip` (default) — windows where `coverage_fraction == 0.0` (entirely
    in a gap with no audio) are skipped. Windows with partial coverage are
    still emitted; their `coverage_fraction` annotates the fill fraction.
    Use this for duty-cycled recordings where off-period windows should be
    invisible to the analysis.
  * `:zero_fill` — all windows in the source time span are emitted, with gap
    regions zero-filled in the returned `Audiodata`.
  * `:error` — all windows are emitted; `read_audio_range` throws
    `ArgumentError` if any window has partial or zero coverage.

Returns:     A `ChunksIterator`, eltype `Audiodata`. Usable in `for` loops and
             with `collect`. `length` is not defined — the iterator has
             `Base.SizeUnknown()` because the chunk count depends on gap
             filtering and cannot be known without iterating.

Constraints: `chunk_seconds` and `stride_seconds` are rounded to millisecond
             precision (the resolution of `DateTime`). Sub-millisecond values
             are coerced silently.

Fails when:
- `gap_handling` is not `:skip`, `:zero_fill`, or `:error` — throws
  `ArgumentError`.
- `chunk_seconds` or `stride_seconds` is not positive — throws `ArgumentError`.

Example:
```julia
src = SingleFileSource("T1-C__0__20170912_181500.wav"; recorder="sm3m")
for chunk in chunks(src; chunk_seconds = 60.0)
    println(chunk.starttime, "  nsamples=", nsamples(chunk))
end
```

Do not use when: You need parallel execution and a DataFrame of results — use
             `process_chunks` instead. Use `read_audio_range` directly for
             one-off arbitrary time windows.
"""
function chunks(source::AbstractAudioSource;
                chunk_seconds::Real,
                stride_seconds::Real = chunk_seconds,
                gap_handling::Symbol = :skip)
    gap_handling ∈ (:skip, :zero_fill, :error) || throw(ArgumentError(
        "chunks: gap_handling must be :skip, :zero_fill, or :error; got :$gap_handling"))
    chunk_seconds  > 0 || throw(ArgumentError(
        "chunks: chunk_seconds must be positive; got $chunk_seconds"))
    stride_seconds > 0 || throw(ArgumentError(
        "chunks: stride_seconds must be positive; got $stride_seconds"))

    t_start, t_end = time_range(source)
    chunk_dur  = Millisecond(round(Int, chunk_seconds  * 1000))
    stride_dur = Millisecond(round(Int, stride_seconds * 1000))
    return ChunksIterator(source, chunk_dur, stride_dur, t_start, t_end, gap_handling)
end

# ─── _chunk_windows ──────────────────────────────────────────────────────────

# Purpose:     Materialise the (start, stop) window pairs that `chunks` would
#              emit into a pre-allocated Vector, applying the same gap-filtering
#              logic. Used by `process_chunks` to build an indexable work list
#              for parallel dispatch via `Threads.@threads`.
# Constraints: Millisecond rounding and :skip filtering are identical to
#              `ChunksIterator.iterate` — the two must remain in sync.
# Fails when:  Same conditions as `chunks`.
function _chunk_windows(source::AbstractAudioSource;
                        chunk_seconds::Real,
                        stride_seconds::Real = chunk_seconds,
                        gap_handling::Symbol = :skip)
    t_start, t_end = time_range(source)
    chunk_dur  = Millisecond(round(Int, chunk_seconds  * 1000))
    stride_dur = Millisecond(round(Int, stride_seconds * 1000))

    windows = Tuple{DateTime,DateTime}[]
    t = t_start
    while t < t_end
        t_stop = t + chunk_dur
        if gap_handling != :skip || coverage_fraction(source, t, t_stop) > 0.0
            push!(windows, (t, t_stop))
        end
        t += stride_dur
    end
    return windows
end

# ─── process_chunks ──────────────────────────────────────────────────────────

"""
    process_chunks(source, f; chunk_seconds, stride_seconds=chunk_seconds,
                   parallel=:threads, device=:cpu, gap_handling=:skip,
                   on_error=:fail, progress=true)

Purpose:     Apply a function to every audio chunk drawn from `source` and
             collect the results as a `DataFrame`. The primary entry point for
             archive-scale soundscape metric computation. Analysis pipelines
             written against `SingleFileSource` run unchanged on an
             `IndexedFileSource` spanning a multi-terabyte archive.

Arguments:
- `source::AbstractAudioSource` — the audio source to iterate.
- `f::Function` — applied to each `Audiodata` chunk. Must return a
  `NamedTuple`. Each field becomes a column in the output DataFrame. Do not
  return fields named `start_time` or `coverage_fraction`; those names are
  reserved and will be overwritten if present.
- `chunk_seconds::Real` — window duration in seconds (positive).
- `stride_seconds::Real` — step between window start times. Default equals
  `chunk_seconds` (non-overlapping).
- `parallel::Symbol`:
  * `:threads` (default) — dispatches over Julia threads via
    `Threads.@threads`. Thread-safe because each `read_audio_range` call
    opens and closes its own file independently. `f` must be thread-safe;
    avoid shared mutable state in closures.
  * `:none` — sequential processing. Useful for debugging and for `f`
    functions with internal mutable state.
  * `:distributed` — not implemented in v1; falls back to `:threads` with a
    warning.
- `device::Symbol` — `:cpu` (default, only active option in v1). `:gpu` and
  `:auto` fall back to `:cpu` with a warning.
- `gap_handling::Symbol` — `:skip` (default), `:zero_fill`, or `:error`.
  Same semantics as `chunks`.
- `on_error::Symbol`:
  * `:fail` (default) — rethrows any exception immediately. Recommended
    during development so bugs surface at once.
  * `:skip` — drops the failed chunk; that row is absent from the output.
    Use for long production runs where isolated corrupt files should not
    abort the job.
  * `:record` — adds a row with a string `error` column describing the
    exception. Metric columns are `missing` for that row.
- `progress::Bool` — `true` (default) displays a `ProgressMeter` progress
  bar on stdout. Thread-safe: `ProgressMeter.next!` uses internal locking.

Returns:     A `DataFrame` with columns `start_time::DateTime` and
             `coverage_fraction::Float64`, plus whatever columns `f` returns.
             Rows follow window order (pre-indexed slots preserve order even
             under concurrent execution). Returns an empty `DataFrame` if no
             windows are found or all chunks failed with `on_error=:skip`.

Constraints: `f` must be thread-safe when `parallel=:threads`. `chunk_seconds`
             and `stride_seconds` are rounded to millisecond precision.

Fails when:
- Any argument validation fails (same conditions as `chunks`).
- `on_error`, `parallel`, or `device` is not one of the listed symbols.
- Any chunk raises an exception and `on_error=:fail`.

Example:
```julia
src = IndexedFileSource(tbl, "/data/deployment/")
df  = process_chunks(src,
                     chunk -> (rms = sqrt(sum(abs2, chunk.sig) / length(chunk.sig)),);
                     chunk_seconds = 60.0,
                     on_error      = :skip)
# df columns: start_time, coverage_fraction, rms
```

Do not use when: You need custom control flow between chunks (early exit,
             stateful accumulators across chunks) — use `chunks` directly.
"""
function process_chunks(source::AbstractAudioSource, f;
                        chunk_seconds::Real,
                        stride_seconds::Real = chunk_seconds,
                        parallel::Symbol     = :threads,
                        device::Symbol       = :cpu,
                        gap_handling::Symbol = :skip,
                        on_error::Symbol     = :fail,
                        progress::Bool       = true)
    gap_handling ∈ (:skip, :zero_fill, :error) || throw(ArgumentError(
        "process_chunks: gap_handling must be :skip, :zero_fill, or :error"))
    on_error ∈ (:fail, :skip, :record) || throw(ArgumentError(
        "process_chunks: on_error must be :fail, :skip, or :record"))
    parallel ∈ (:none, :threads, :distributed) || throw(ArgumentError(
        "process_chunks: parallel must be :none, :threads, or :distributed"))
    device ∈ (:cpu, :gpu, :auto) || throw(ArgumentError(
        "process_chunks: device must be :cpu, :gpu, or :auto"))

    if parallel == :distributed
        @warn "process_chunks: parallel=:distributed is not implemented in v1; " *
              "falling back to :threads"
        parallel = :threads
    end
    if device != :cpu
        @warn "process_chunks: device=$device is not implemented in v1; " *
              "falling back to :cpu"
    end

    windows = _chunk_windows(source; chunk_seconds, stride_seconds, gap_handling)
    n = length(windows)
    n == 0 && return DataFrame()

    # :skip already removed zero-coverage windows in _chunk_windows. Use :zero_fill
    # for the actual reads so that partial-coverage boundary windows are returned
    # rather than raising; fully-in-gap windows were already excluded above.
    gr = gap_handling == :skip ? :zero_fill : gap_handling

    results = Vector{Any}(undef, n)
    meter   = progress ? Progress(n; desc = "Processing chunks: ") : nothing

    # Merge reserved fields with whatever f returns. start_time and
    # coverage_fraction are written second so they overwrite any same-named
    # fields from f (documented constraint: f should not return these names).
    function _process_one(i)
        t0, t1 = windows[i]
        chunk  = read_audio_range(source, t0, t1; gap_handling = gr)
        cov    = coverage_fraction(source, t0, t1)
        user   = f(chunk)
        return merge((start_time = t0, coverage_fraction = cov), user)
    end

    if parallel == :threads
        Threads.@threads for i in 1:n
            try
                results[i] = _process_one(i)
            catch e
                if on_error == :fail
                    rethrow()
                elseif on_error == :skip
                    results[i] = missing
                else  # :record
                    results[i] = (start_time       = windows[i][1],
                                  coverage_fraction = NaN,
                                  error             = sprint(showerror, e))
                end
            end
            isnothing(meter) || next!(meter)
        end
    else  # :none
        for i in 1:n
            try
                results[i] = _process_one(i)
            catch e
                if on_error == :fail
                    rethrow()
                elseif on_error == :skip
                    results[i] = missing
                else  # :record
                    results[i] = (start_time       = windows[i][1],
                                  coverage_fraction = NaN,
                                  error             = sprint(showerror, e))
                end
            end
            isnothing(meter) || next!(meter)
        end
    end

    valid = [r for r in results if !ismissing(r)]
    isempty(valid) && return DataFrame()

    # When on_error=:record and errors occurred, success rows and error rows have
    # different NamedTuple schemas (metric columns vs. an `error` string column).
    # vcat with cols=:union unifies them, filling absent columns with missing.
    # For the common case (all same schema), convert to a concrete typed Vector
    # so DataFrames builds columns without a row-by-row vcat.
    has_mixed = on_error == :record && any(r -> :error ∈ keys(r), valid)
    if has_mixed
        return reduce(vcat, [DataFrame([r]) for r in valid]; cols = :union)
    end

    T = typeof(first(valid))
    return DataFrame(convert(Vector{T}, valid))
end
