# ═══════════════════════════════════════════════════════════════════════════════
# Per-cycle band-metrics table
#
# Turns one recording (typically one duty-cycle file) into a single tidy row of
# candidate-band metrics, and a collection of recordings into a DataFrame with
# one row per cycle. This is the hand-off product for the R-side soundscape
# analysis (Ch. 1): each candidate band contributes a set of dB columns
# (energetic mean, median, the nine EA percentiles, and the max) computed over
# the file's PSD frames.
#
# Composition only: band_metrics = compute_psd → compute_spl, then flatten the
# BandSPL structs into named columns. No new signal math.
#
# NOTE (test-run parameterisation): the metric set, the FFT window, and the
# choice to take percentiles over PSD frames are provisional for the first
# full-pipeline run and may change for production. See CANDIDATE_ACOUSTIC_BANDS.org.
# ═══════════════════════════════════════════════════════════════════════════════

# Metrics emitted per band, in column order. `:max` is the maximum of the
# per-frame series; the rest map directly to BandSPL fields. EA percentile
# convention: L_n = n-th percentile (level BELOW which n% of frames fall), so
# L1/L5 are the quiet floor and L95/L99 the loud tail — the reverse of the
# acoustic "exceeded n%" convention (see docs/src/explanations/spl.md).
const METRIC_ORDER = (:mean, :median, :L1, :L5, :L10, :L25, :L75, :L90, :L95, :L99, :max)

# Purpose:     Return the dB value of one metric from a BandSPL. `:max` is
#              computed from the per-frame series; all others are stored fields.
# Constraints: `m` must be one of METRIC_ORDER.
# Fails when:  `m` is not a recognised metric name (throws).
function _metric_value(b::BandSPL, m::Symbol) :: Float64
    m === :mean   && return b.mean_dB
    m === :median && return b.median_dB
    m === :L1     && return b.L1_dB
    m === :L5     && return b.L5_dB
    m === :L10    && return b.L10_dB
    m === :L25    && return b.L25_dB
    m === :L75    && return b.L75_dB
    m === :L90    && return b.L90_dB
    m === :L95    && return b.L95_dB
    m === :L99    && return b.L99_dB
    m === :max    && return maximum(b.spl_dB)
    error("_metric_value: unknown metric :$m")
end

# Purpose:     Deterministic list of the band-metric column names for a given
#              band set: one `<bandkey>_<metric>_dB` symbol per (band, metric),
#              bands in sorted-key order, metrics in METRIC_ORDER. Used to build
#              typed DataFrame columns in a stable order.
# Constraints: `bands` non-empty for a non-trivial table.
# Fails when:  Never.
function _band_metric_colnames(bands::Dict{Symbol, Tuple{Float64, Float64}}) :: Vector{Symbol}
    cols = Symbol[]
    for k in sort(collect(keys(bands)))
        for m in METRIC_ORDER
            push!(cols, Symbol(k, :_, m, :_dB))
        end
    end
    return cols
end

"""
    band_metrics(audio::Audiodata; bands, window_seconds, overlap_fraction=0.5,
                 window=:hann, nfft=nothing, fft_plan=nothing, cal=nothing,
                 environment=:water) -> Dict{Symbol, Any}

Purpose:     Compute one tidy row of candidate-band metrics for a single
             recording (typically one duty-cycle file). Calls
             [`compute_psd`](@ref) then [`compute_spl`](@ref) and flattens the
             per-band [`BandSPL`](@ref) results into named columns. The
             percentiles are taken over the file's PSD frames.

Arguments:
- `audio::Audiodata`: The recording. Must resolve to a calibration (recorder
  profile, attached calibration, or `cal=`), or `compute_spl` will reject it.
- `bands::Dict{Symbol,Tuple{Float64,Float64}}`: Required (DD-27). Each entry maps
  a label to a `(low_Hz, high_Hz)` band.
- `window_seconds::Real`: FFT window duration in seconds (the per-frame
  timescale over which percentiles are taken). Required.
- `overlap_fraction`, `window`, `nfft`, `fft_plan`: Forwarded to `compute_psd`.
- `cal::Union{Calibration,Nothing}=nothing`: Override calibration, forwarded to
  `compute_psd`. Use for the SM3M per-era workaround (e.g.
  `cal = ScalarCalibration(-165f0)` for post-Apr-2018 files).
- `environment::Symbol=:water`: `:water` (reference 1 µPa) or `:air` (20 µPa).

Returns:     `Dict{Symbol,Any}` — one row. Always contains `:start_time`
             (`DateTime`), `:duration_s` (`Float64`), `:n_frames` (`Int`, the
             number of PSD frames the percentiles are computed over), plus, for
             every band `k` and metric `m` in
             `(:mean,:median,:L1,:L5,:L10,:L25,:L75,:L90,:L95,:L99,:max)`, a key
             `Symbol(k, :_, m, :_dB)` holding a dB value.

Constraints:
- Percentiles are over PSD frames, so `n_frames` (≈ file_duration ÷
  hop) sets their statistical resolution. A very short file yields few frames
  and noisy percentiles.
- EA percentile direction: `L1`/`L5` are the quiet floor, `L95`/`L99` the loud
  tail (see `docs/src/explanations/spl.md`).

Fails when:  `bands` is empty (`AssertionError`); the audio cannot be calibrated
             (`AssertionError` from `compute_spl`, DD-21); any `compute_psd` /
             `compute_spl` failure (short signal, band above Nyquist, …).

Example:
```julia
audio = read_audio("cycle_0001.wav"; recorder = "sm3m")
bands = Dict(:tol_63 => (56.2, 70.8), :b100_200 => (100.0, 200.0))
row   = band_metrics(audio; bands = bands, window_seconds = 1.0)
row[:b100_200_mean_dB]   # energetic-mean 100–200 Hz level for this cycle
```

Do not use when: You need a coarse per-column series over a long continuous
             recording — use [`compute_spl`](@ref) on an `LTSAResult` instead.
"""
function band_metrics(audio::Audiodata;
                      bands::Dict{Symbol, Tuple{Float64, Float64}},
                      window_seconds::Real,
                      overlap_fraction::Real          = 0.5,
                      window::Symbol                  = :hann,
                      nfft::Union{Int,Nothing}        = nothing,
                      fft_plan                        = nothing,
                      cal::Union{Calibration,Nothing} = nothing,
                      environment::Symbol             = :water) :: Dict{Symbol, Any}

    @assert !isempty(bands) "band_metrics: `bands` must not be empty"

    psd = compute_psd(audio; window_seconds = window_seconds,
                      overlap_fraction = overlap_fraction, window = window,
                      nfft = nfft, fft_plan = fft_plan, cal = cal)
    spl = compute_spl(psd; bands = bands, environment = environment)

    row = Dict{Symbol, Any}()
    row[:start_time] = audio.starttime
    row[:duration_s] = nsamples(audio) / Float64(audio.fs)
    row[:n_frames]   = length(psd.time)

    for (k, b) in spl.bands
        for m in METRIC_ORDER
            row[Symbol(k, :_, m, :_dB)] = _metric_value(b, m)
        end
    end
    return row
end

"""
    band_table(items; reader=identity, bands, window_seconds,
               overlap_fraction=0.5, window=:hann, nfft=nothing, cal=nothing,
               environment=:water, on_error=:skip) -> DataFrame

Purpose:     Build a one-row-per-cycle table of candidate-band metrics from a
             collection of recordings. Maps [`band_metrics`](@ref) over `items`
             and assembles the rows into a typed `DataFrame` — the hand-off
             product for the R-side soundscape analysis.

Arguments:
- `items`: Any iterable of things to turn into rows — `Audiodata` objects, or
  file paths (with a matching `reader`), or anything `reader` maps to
  `Audiodata`.
- `reader = identity`: Function `item -> Audiodata`. Default `identity` treats
  each item as an `Audiodata`. For files, pass e.g.
  `reader = f -> read_audio(f; recorder = "sm3m")`. Reading happens inside the
  per-item `try`, so unreadable/corrupt files are skipped under `on_error=:skip`.
- `bands`, `window_seconds`, `overlap_fraction`, `window`, `nfft`, `cal`,
  `environment`: Forwarded to `band_metrics` for every item.
- `fft_plan = nothing`: Optional pre-built FFTW plan shared across all items
  (and threads). When `nothing` and `parallel=:threads`, one plan is built from
  the first item's sampling rate and reused; when `nothing` and `parallel=:none`,
  each item builds its own. Build with [`make_spectrogram_plan`](@ref).
- `parallel::Symbol = :none`: `:none` processes items serially; `:threads` uses
  `Threads.@threads` over items (start Julia with `-t N`). Output is identical
  either way — rows are assembled by item index, so the threaded result is
  deterministic and matches the serial result row-for-row.
- `on_error::Symbol = :skip`: `:skip` logs and drops any item whose read or
  metric computation throws; `:fail` re-throws (wrapped in a task exception under
  `:threads`).
- `progress::Bool = false`: show a progress bar.

Returns:     `DataFrame`, one row per successfully processed item, **in item
             order**. Columns: `start_time`, `duration_s`, `n_frames`, then
             `<bandkey>_<metric>_dB` for every band (sorted-key order) × metric
             (`mean, median, L1, L5, L10, L25, L75, L90, L95, L99, max`).
             Columns are typed (`DateTime`, `Float64`, `Int`, `Float64`) even
             when the table is empty. Pass items already sorted by time (e.g.
             `sort(readdir(...))`) to get a time-ordered table.

Constraints:
- Column set is derived from `bands` alone, so every row has the same schema
  regardless of which items succeeded.
- Under `parallel=:threads`, the shared FFT plan assumes a uniform sampling rate
  across items (true for a single-recorder deployment); an item with a different
  `fs` fails the spectrogram size assertion and is skipped.
- `items` is materialised with `collect` for index-based assembly, so pass a
  vector of file paths, not a generator that reads audio eagerly.

Fails when:  `on_error` is not `:skip` or `:fail` (`ArgumentError`); under
             `on_error=:fail`, the first per-item error propagates.

Example:
```julia
files = readdir("deployment/"; join = true)          # 3 months of 5-min files
bands = Dict(:tol_63   => (56.2, 70.8),
             :b100_200 => (100.0, 200.0),
             :b1_4k    => (1000.0, 4000.0))
df = band_table(files;
                reader = f -> read_audio(f; recorder = "sm3m"),
                bands = bands, window_seconds = 1.0)
# df: one row per duty cycle; write to CSV for R with CSV.write("bands.csv", df)
```

Do not use when: You need gap-aware coverage over a continuous timeline — that
             is the source-backed archive-scale workflow, not this per-file map.
"""
function band_table(items;
                    reader                          = identity,
                    bands::Dict{Symbol, Tuple{Float64, Float64}},
                    window_seconds::Real,
                    overlap_fraction::Real          = 0.5,
                    window::Symbol                  = :hann,
                    nfft::Union{Int,Nothing}        = nothing,
                    fft_plan                        = nothing,
                    cal::Union{Calibration,Nothing} = nothing,
                    environment::Symbol             = :water,
                    parallel::Symbol                = :none,
                    on_error::Symbol                = :skip,
                    progress::Bool                  = false) :: DataFrame

    on_error in (:skip, :fail) ||
        throw(ArgumentError("band_table: on_error must be :skip or :fail, got :$on_error"))
    parallel in (:none, :threads) ||
        throw(ArgumentError("band_table: parallel must be :none or :threads, got :$parallel"))

    # Materialise the item list so items can be indexed — required for the
    # thread-safe, order-preserving assembly below (each item writes into its
    # own slot; no shared push!). For the file-path workflow this is a cheap
    # vector of strings; do NOT pass a generator that reads audio eagerly.
    itemvec = collect(items)
    n       = length(itemvec)
    # One slot per item; `nothing` marks an item that errored / was skipped.
    results = Vector{Union{Nothing, Dict{Symbol, Any}}}(nothing, n)

    # Shared FFT plan for the threaded path. Building ONE plan and reusing it
    # across threads avoids concurrent FFTW planning (not guaranteed safe) and
    # repeated per-file plan builds. The first item's sampling rate sizes the
    # plan; this assumes a uniform fs across items (true for a single-recorder
    # deployment). An item whose fs differs fails the spectrogram size assertion
    # and is skipped. If `fft_plan` is supplied, it is used directly.
    shared_plan = fft_plan
    if parallel === :threads && shared_plan === nothing && n > 0
        try
            probe       = reader(itemvec[1])
            shared_plan = make_spectrogram_plan(Float64(probe.fs), window_seconds, nfft)
        catch e
            on_error === :fail && rethrow(e)
            @warn "band_table: could not build shared FFT plan from the first item; " *
                  "each call will build its own plan" exception = e
        end
    end

    prog = progress ? Progress(n; desc = "band_table ", dt = 1.0) : nothing

    # Per-item work writing into a private index — safe to call from many threads.
    function _process!(i)
        try
            audio      = reader(itemvec[i])
            results[i] = band_metrics(audio; bands = bands,
                                      window_seconds = window_seconds,
                                      overlap_fraction = overlap_fraction,
                                      window = window, nfft = nfft,
                                      fft_plan = shared_plan, cal = cal,
                                      environment = environment)
        catch e
            on_error === :fail && rethrow(e)
            @warn "band_table: skipping item due to error" item = itemvec[i] exception = e
        end
        prog === nothing || next!(prog)
        return nothing
    end

    if parallel === :threads
        Threads.@threads for i in 1:n
            _process!(i)
        end
    else
        for i in 1:n
            _process!(i)
        end
    end

    # Collect surviving rows in item order (results is index-ordered, so the
    # threaded output is identical to the serial output — deterministic).
    rows      = Dict{Symbol, Any}[r for r in results if r !== nothing]
    n_skipped = n - length(rows)
    n_skipped > 0 && @info "band_table: skipped $n_skipped of $n item(s) with errors"

    # Build typed columns. Comprehensions read from Dict{Symbol,Any}, so the
    # element type is annotated explicitly (DateTime/Float64/Int) — this keeps
    # columns concretely typed even when `rows` is empty.
    df = DataFrame()
    df.start_time = DateTime[r[:start_time] for r in rows]
    df.duration_s = Float64[r[:duration_s] for r in rows]
    df.n_frames   = Int[r[:n_frames] for r in rows]
    for name in _band_metric_colnames(bands)
        df[!, name] = Float64[r[name] for r in rows]
    end
    return df
end
