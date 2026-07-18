# ═══════════════════════════════════════════════════════════════════════════════
# Long-Term Spectral Average (LTSA)
#
# An LTSA summarises how acoustic energy is distributed across frequency over long
# timescales (hours to days to a full deployment). The recording is divided into
# consecutive, non-overlapping time columns of `average_span_seconds` each; within
# each column a PSD is computed (via the spectrogram → compute_psd chain) and its
# FFT frames are averaged in linear power (the energetic mean). Each averaged PSD
# becomes one column of a frequency × time matrix.
#
# This matrix is the reusable pre-aggregation product for downstream soundscape
# work: arbitrary frequency bands (TOLs, octave bands, hand-picked exploratory
# bands) are integrated out of the same matrix without re-reading audio, and
# percentiles are taken down the time axis.
#
# Design decisions: energetic-mean-only averaging and non-overlapping columns are
# recorded in docs/design_decisions.md (DD-28). Composes compute_psd (DD-08..DD-16)
# and spectrogram (DD-01..DD-07).
# ═══════════════════════════════════════════════════════════════════════════════

"""
    LTSAResult

Purpose:     Container for a computed Long-Term Spectral Average — a
             frequency × time matrix of linear power spectral density, plus the
             axes and calibration provenance needed to interpret and plot it.

Fields:
- `matrix::Matrix{Float64}`: PSD values, shape `(n_freqs, n_columns)`. Row `k`
  is frequency `freqs[k]`; column `j` is time column `j`. Linear units:
  µPa²/Hz when `is_calibrated`, full-scale²/Hz otherwise. Convert to dB with
  [`to_dB`](@ref).
- `freqs::Vector{Float64}`: Frequency axis in Hz, length `n_freqs`. Identical
  for every column (same `fs`, window, and `nfft` across columns).
- `column_times::Vector{Float64}`: Start time of each column in seconds relative
  to the start of the analysed audio, length `n_columns`. Column `j` covers
  `column_times[j] .. column_times[j] + average_span_seconds`.
- `fs::Float32`: Sampling rate in Hz.
- `average_span_seconds::Float64`: Column width in seconds (the averaging span).
- `is_calibrated::Bool`: `true` when `matrix` is in physical units (µPa²/Hz).
- `cal::Calibration`: The calibration applied when building the columns.

Constraints: Every column is an energetic (linear-power) mean of the FFT frames
             falling inside its span; there is no dB-mean or median form.
             Columns never overlap and never cross each other's boundaries.

Units symbol: query with [`ltsa_units`](@ref).
"""
struct LTSAResult
    matrix::Matrix{Float64}
    freqs::Vector{Float64}
    column_times::Vector{Float64}
    fs::Float32
    average_span_seconds::Float64
    is_calibrated::Bool
    cal::Calibration
end

"""
    ltsa_units(result::LTSAResult) -> Symbol

Purpose:     Return the physical unit of `result.matrix` as a `Symbol`, mirroring
             [`psd_units`](@ref). Dispatches on the `is_calibrated` flag, not on
             the calibration type.

Returns:     `:µPa²_per_Hz` when `result.is_calibrated`, else
             `:fullscale²_per_Hz`.

Fails when:  Never.

Example:
```julia
lt = compute_ltsa(audio; average_span_seconds = 60.0)
ltsa_units(lt)   # :µPa²_per_Hz for a calibrated recording
```
"""
ltsa_units(result::LTSAResult) =
    result.is_calibrated ? :µPa²_per_Hz : :fullscale²_per_Hz

"""
    to_dB(result::LTSAResult) -> Matrix{Float64}

Purpose:     Convert the linear LTSA matrix to decibels: `10 × log10(matrix)`.
             Returns a new matrix the same shape as `result.matrix`; the result
             is not modified.

             Units after conversion:
             - `ltsa_units(result) == :µPa²_per_Hz` → dB re 1 µPa²/Hz
             - `ltsa_units(result) == :fullscale²_per_Hz` → dBFS/Hz

Returns:     `Matrix{Float64}`, shape `(n_freqs, n_columns)`, all values in dB.

Constraints: Matrix values must be strictly positive; `log10(0) = -Inf` and
             `log10(negative) = NaN` are left as-is (IEEE 754), consistent with
             `to_dB` for PSDResult.

Fails when:  Never.

Example:
```julia
lt    = compute_ltsa(audio; average_span_seconds = 60.0)
lt_dB = to_dB(lt)          # Matrix{Float64}, dB re 1 µPa²/Hz if calibrated
```
"""
to_dB(result::LTSAResult) = 10.0 .* log10.(result.matrix)

# Purpose:     Compute the averaged (energetic-mean) PSD for one LTSA column.
#              Slices `sig` to the `col_index`-th non-overlapping span of
#              `samples_per_column` samples, runs the STFT with the given FFT
#              parameters, converts to linear PSD, and averages over the frames.
# Constraints: `col_index` must be in `1:n_columns`. The slice must hold at least
#              one full FFT window — guaranteed by the
#              `samples_per_column >= window_length` assertion in `compute_ltsa`.
#              `cal` must be pre-resolved: this helper applies it directly and
#              never warns about missing calibration (the warning is emitted once
#              in `compute_ltsa`, not per column).
# Fails when:  The requested slice runs past the end of `sig` (a `col_index` too
#              large for the signal) — throws `BoundsError` via the `@view`.
function _ltsa_column(sig::AbstractVector{Float64},
                      col_index::Int,
                      samples_per_column::Int,
                      fs::Float64,
                      fft_window_seconds::Real,
                      fft_overlap::Real,
                      window::Symbol,
                      nfft::Union{Int,Nothing},
                      fft_plan,
                      cal::Calibration) :: PSDResult
    # 1-indexed first sample of this column; columns are back-to-back with no gap.
    col_start = (col_index - 1) * samples_per_column + 1

    # @view slices without copying; the STFT reads it directly. seg is an
    # AbstractVector{Float64}, so spectrogram dispatches to its Float64 method.
    seg = @view sig[col_start : col_start + samples_per_column - 1]

    spec = spectrogram(seg;
                       fs               = fs,
                       window_seconds   = fft_window_seconds,
                       overlap_fraction = fft_overlap,
                       window           = window,
                       nfft             = nfft,
                       fft_plan         = fft_plan)

    # compute_psd(spec, cal) returns the frame-resolved linear PSD (freq × frame);
    # average_psd collapses the frames to one spectrum per column.
    return compute_psd(spec, cal)
end

"""
    compute_ltsa(audio::Audiodata; average_span_seconds,
                 fft_window_seconds=1.0, fft_overlap=0.5, window=:hann,
                 nfft=nothing, fft_plan=nothing, cal=nothing) -> LTSAResult

Purpose:     Compute a Long-Term Spectral Average from an `Audiodata` recording.
             The signal is divided into consecutive non-overlapping columns of
             `average_span_seconds` each; within each column a PSD is computed
             and its FFT frames are averaged in linear power (energetic mean).
             Each averaged PSD is one column of the returned frequency × time
             matrix. A trailing partial column (fewer than `average_span_seconds`
             worth of samples) is dropped.

Arguments:
- `audio::Audiodata`: The recording to analyse. `audio.sig` is `Vector{Float64}`.
- `average_span_seconds::Real`: Column width in seconds — the averaging span
  (e.g. `60.0` for 1-minute columns, `3600.0` for hourly). Must be `> 0` and
  large enough to hold at least one FFT window
  (`average_span_seconds × fs ≥ fft_window_seconds × fs`).
- `fft_window_seconds::Real = 1.0`: Inner STFT window duration in seconds,
  forwarded to [`spectrogram`](@ref) as `window_seconds`.
- `fft_overlap::Real = 0.5`: Inner STFT frame overlap fraction in `[0, 1)`,
  forwarded to [`spectrogram`](@ref) as `overlap_fraction`.
- `window::Symbol = :hann`: Inner STFT window shape (`:hann`, `:hamming`,
  `:blackman`, `:rectangular`).
- `nfft::Union{Int,Nothing} = nothing`: FFT length. `nothing` uses the window
  length; a larger even value zero-pads.
- `fft_plan = nothing`: Pre-built FFTW plan reused across every column (all
  columns share one FFT size). Build with [`make_spectrogram_plan`](@ref); it
  must match the `nfft` the inner `spectrogram` chooses, or an `AssertionError`
  is thrown.
- `cal::Union{Calibration,Nothing} = nothing`: Override calibration. When
  `nothing` (default), calibration is auto-resolved once from `audio` via the
  same cascade as [`compute_psd`](@ref) (pre-calibrated signal → recorder
  profile → warn and fall back to full-scale). When supplied, `cal` is applied
  directly and the auto-resolution (and its possible warning) is skipped.

Returns:     [`LTSAResult`](@ref). `matrix` is `(n_freqs, n_columns)` in µPa²/Hz
             (calibrated) or full-scale²/Hz (uncalibrated); check
             [`ltsa_units`](@ref) or `result.is_calibrated`.
             `n_columns = floor(nsamples(audio) / (average_span_seconds × fs))`.

Constraints:
- Columns are non-overlapping and consecutive: column `j` covers samples
  `(j-1)·S+1 .. j·S`, where `S = round(average_span_seconds × fs)`. The inner
  FFT frames of one column never cross into an adjacent column.
- Column averaging is the energetic (linear-power) mean only, matching
  [`average_psd`](@ref) and the Merchant 2015 convention. There is no dB-mean
  or median form (DD-28).
- If the audio length is not an exact multiple of the column span, the final
  partial column is dropped — its samples do not appear in any column.
- Calibration is resolved once, before the column loop; a missing-calibration
  warning (if any) is emitted a single time, not per column.

Fails when:
- `average_span_seconds <= 0` — `AssertionError`.
- `fft_window_seconds <= 0` — `AssertionError`.
- `average_span_seconds × fs < fft_window_seconds × fs` (a column cannot hold one
  FFT window) — `AssertionError`.
- `nsamples(audio) < round(average_span_seconds × fs)` (audio shorter than one
  full column) — `AssertionError`; there is no column to emit.
- Any condition that makes the inner [`spectrogram`](@ref) or [`compute_psd`](@ref)
  fail (e.g. odd `nfft`, mismatched `fft_plan`).

Example:
```julia
audio = read_audio("deployment_hour.flac"; recorder = "sm3m")
lt    = compute_ltsa(audio; average_span_seconds = 60.0)   # 1-minute columns
size(lt.matrix)        # (n_freqs, n_minutes)
lt_dB = to_dB(lt)      # dB re 1 µPa²/Hz (SM3M scalar calibration auto-resolved)
```

Do not use when:
- The recording is shorter than one `average_span_seconds` column — there is
  nothing to average; call [`compute_psd`](@ref) directly instead.
- A single deployment-wide average spectrum is wanted rather than a time series
  — use [`compute_psd`](@ref) + [`average_psd`](@ref) over the whole signal.
"""
function compute_ltsa(audio::Audiodata;
                      average_span_seconds::Real,
                      fft_window_seconds::Real        = 1.0,
                      fft_overlap::Real               = 0.5,
                      window::Symbol                  = :hann,
                      nfft::Union{Int,Nothing}        = nothing,
                      fft_plan                        = nothing,
                      cal::Union{Calibration,Nothing} = nothing) :: LTSAResult

    fs = Float64(audio.fs)

    @assert average_span_seconds > 0 begin
        "compute_ltsa: average_span_seconds=$average_span_seconds must be > 0"
    end
    @assert fft_window_seconds > 0 begin
        "compute_ltsa: fft_window_seconds=$fft_window_seconds must be > 0"
    end

    # Column width and inner-window width in samples. round(Int, ·) converts the
    # real-valued seconds×Hz products to sample counts (same convention as
    # spectrogram uses for its own window length).
    samples_per_column = round(Int, average_span_seconds * fs)
    window_length      = round(Int, fft_window_seconds * fs)

    @assert samples_per_column >= window_length begin
        "compute_ltsa: a column of average_span_seconds=$average_span_seconds " *
        "holds $samples_per_column samples, fewer than the FFT window " *
        "(fft_window_seconds=$fft_window_seconds → $window_length samples). " *
        "Increase average_span_seconds or decrease fft_window_seconds."
    end

    # Number of whole columns; div is integer floor division, so a trailing
    # partial column is dropped (DD-28).
    total_samples = nsamples(audio)
    n_columns     = div(total_samples, samples_per_column)
    @assert n_columns >= 1 begin
        "compute_ltsa: audio has $total_samples samples, fewer than one column " *
        "span (average_span_seconds=$average_span_seconds → $samples_per_column " *
        "samples). Nothing to average — use compute_psd directly for short audio."
    end

    # Resolve calibration ONCE (not per column) so a missing-calibration warning
    # fires at most once. When audio.is_calibrated, _psd_calibration returns
    # NoCalibration (signal already in physical units); is_cal must still be true.
    resolved_cal = cal !== nothing ? cal : _psd_calibration(audio)
    is_cal       = audio.is_calibrated || !(resolved_cal isa NoCalibration)

    # Column 1 establishes the frequency axis and matrix size. Every column has
    # the identical frequency axis (same fs / window / nfft), so we take freqs
    # from the first column and reuse it.
    first_psd = _ltsa_column(audio.sig, 1, samples_per_column, fs,
                             fft_window_seconds, fft_overlap, window,
                             nfft, fft_plan, resolved_cal)
    freqs  = first_psd.freqs
    matrix = Matrix{Float64}(undef, length(freqs), n_columns)
    matrix[:, 1] = average_psd(first_psd)

    for j in 2:n_columns
        psd = _ltsa_column(audio.sig, j, samples_per_column, fs,
                           fft_window_seconds, fft_overlap, window,
                           nfft, fft_plan, resolved_cal)
        matrix[:, j] = average_psd(psd)
    end

    # Column start times in seconds, relative to the audio start. Derived from
    # the actual integer sample offsets so the axis reflects the rounding used
    # for samples_per_column rather than the requested seconds.
    column_times = [(j - 1) * samples_per_column / fs for j in 1:n_columns]

    return LTSAResult(matrix, freqs, column_times, audio.fs,
                      Float64(average_span_seconds), is_cal, resolved_cal)
end

"""
    compute_spl(ltsa::LTSAResult; bands, environment=:water) -> SPLResult

Purpose:     Integrate frequency bands out of an LTSA, producing one band-level
             time series per band. Each LTSA time column is treated as a single
             spectral estimate: for every column the power is summed across the
             band's bins, giving a per-column band SPL series plus the energetic
             mean and percentiles over columns. This is the band time series used
             for temporal / diel / seasonal soundscape plots. Reuses the same
             `_integrate_bands` core as `compute_spl(::PSDResult)` — the only
             difference is that a "column" here is an LTSA time column rather than
             an FFT frame.

Arguments:
- `ltsa::LTSAResult`: A calibrated LTSA (µPa²/Hz). Uncalibrated input fails.
- `bands::Dict{Symbol,Tuple{Float64,Float64}}`: Required (DD-27). Each entry maps
  a label to a `(low_Hz, high_Hz)` band; a bin is included when its centre lies
  in `[low_Hz, high_Hz]`. No default.
- `environment::Symbol = :water`: `:water` (reference 1 µPa) or `:air` (20 µPa).

Returns:     [`SPLResult`](@ref). Each `BandSPL.spl_dB` is the per-column series
             (length = number of LTSA columns); `SPLResult.time` is
             `ltsa.column_times` (seconds from the audio start). Units
             `:dB_re_1µPa` for water.

Constraints:
- `ltsa` must be calibrated (`ltsa_units(ltsa) === :µPa²_per_Hz`).
- Percentiles are taken over LTSA columns, so their statistical meaning depends
  on the column span (`ltsa.average_span_seconds`). A handful of coarse columns
  gives a coarse distribution.

Fails when:  Uncalibrated LTSA (`AssertionError`, DD-21); a band with
             `low_Hz ≥ high_Hz` or `high_Hz > Nyquist` (`ArgumentError`); no bins
             fall in a band (`ArgumentError`).

Example:
```julia
lt  = compute_ltsa(audio; average_span_seconds = 60.0)
spl = compute_spl(lt; bands = Dict(:b100_200 => (100.0, 200.0)))
spl.bands[:b100_200].spl_dB    # per-minute band level series (dB re 1 µPa)
spl.bands[:b100_200].mean_dB   # energetic mean over all columns
```

Do not use when: You want a fine within-file distribution for percentiles — a
             file's PSD frames (`compute_spl(compute_psd(audio); bands)`) give
             hundreds of samples per file, whereas an LTSA offers only as many
             samples as it has columns.
"""
function compute_spl(ltsa::LTSAResult;
                     bands::Dict{Symbol, Tuple{Float64, Float64}},
                     environment::Symbol = :water) :: SPLResult
    return _integrate_bands(ltsa.matrix, ltsa.freqs, ltsa.column_times, ltsa.fs,
                            ltsa_units(ltsa);
                            bands = bands, environment = environment)
end
