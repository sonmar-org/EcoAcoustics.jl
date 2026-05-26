# ─── Power Spectral Density ──────────────────────────────────────────────────
#
# Computes PSD from a SpectrogramResult using the Merchant et al. (2015)
# normalisation (Eq. 1). Calibration is applied here; the spectrogram layer is
# unit-agnostic. See docs/design_decisions.md DD-08 through DD-11 for the
# conventions embedded in this file.

"""
    PSDResult

Result of a power spectral density computation from a [`SpectrogramResult`](@ref).
Units are linear power: full-scale²/Hz when `is_calibrated = false`, µPa²/Hz
when `is_calibrated = true`.

For dB values, pass this object to `to_dB` (deliverable 4); for frame-averaged
PSD, pass to `average_psd`.

Fields
------
- `psd_linear::Matrix{Float64}`:
    PSD in linear power units, size `(n_freqs, n_frames)`. Row `k` is frequency
    bin `k`; column `j` is frame `j`. Units: full-scale²/Hz (`is_calibrated =
    false`) or µPa²/Hz (`is_calibrated = true`). Produced by squaring the STFT
    magnitudes, applying the single-sided correction (DD-10), dividing by
    `fs × window_energy` (Merchant 2015 Eq. 1), then calling
    `apply_calibration!` (DD-08).
- `freqs::Vector{Float64}`:
    Bin-centre frequencies in Hz, inherited from [`SpectrogramResult`](@ref).
    DC at index 1, Nyquist at the last index.
- `time::Vector{Float64}`:
    Frame-centre times in seconds from the signal start, inherited from
    [`SpectrogramResult`](@ref). PAMGuide convention (DD-03).
- `fs::Float32`:
    Sample rate in Hz, inherited from [`SpectrogramResult`](@ref).
- `window_energy::Float64`:
    Sum of squared window coefficients (Σwᵢ²), inherited from
    [`SpectrogramResult`](@ref). Carried for downstream band integration and
    chain-consistency checks — consumers can verify the PSD was computed with
    the expected window without access to the original `SpectrogramResult`.
- `nfft::Int`:
    FFT length used to produce this PSD, inherited from
    [`SpectrogramResult`](@ref). Used for frequency-band integration (bin width
    = `fs / nfft`) and zero-padding bookkeeping.
- `cal::Calibration`:
    Calibration object that was applied. One of `NoCalibration`,
    `ScalarCalibration`, or `TFCalibration`. `NoCalibration()` when no
    calibration was requested. Query units programmatically via `psd_units`.
- `is_calibrated::Bool`:
    `true` when `cal` is not `NoCalibration`. Provided for fast Boolean guards
    in metric functions; use `psd_units(result)` when the specific unit symbol
    is needed.

References
----------
Merchant et al. (2015) Measuring Acoustic Habitats. Methods in Ecology and
Evolution, 6, 257–265. Normalisation equation (Eq. 1) and single-sided
correction convention.
"""
struct PSDResult
    psd_linear::Matrix{Float64}
    freqs::Vector{Float64}
    time::Vector{Float64}
    fs::Float32
    window_energy::Float64
    nfft::Int
    cal::Calibration
    is_calibrated::Bool
end

"""
    psd_units(result::PSDResult) -> Symbol

Purpose:     Return a Symbol describing the physical units of a `PSDResult`'s
             `psd_linear` field.

Arguments:
- `result::PSDResult`: A computed PSD result.

Returns:
- `:µPa²_per_Hz` when `result.is_calibrated` is `true`. This covers both the
  common case (calibration applied at the PSD layer via `apply_calibration!`)
  and the time-domain pre-calibration case (signal already in µPa before
  spectrogram, PSD wrapper set `is_calibrated = true`).
- `:fullscale²_per_Hz` when `result.is_calibrated` is `false`; power is
  relative to ADC full-scale squared per Hz, not a physical unit.

Constraints: No-fail. Dispatches solely on the `Bool` field `is_calibrated`.

Example:
```julia
r = compute_psd(audio; window_seconds = 1.0)   # auto-calibrates
psd_units(r)   # :µPa²_per_Hz
```
"""
psd_units(p::PSDResult) = p.is_calibrated ? :µPa²_per_Hz : :fullscale²_per_Hz

"""
    compute_psd(spec, cal=NoCalibration()) -> PSDResult

Purpose:     Compute power spectral density from a [`SpectrogramResult`](@ref).
             Converts the complex STFT into a calibrated (or uncalibrated) PSD
             matrix using the Merchant et al. (2015) normalisation convention.
             This is the lowest-level PSD primitive; all higher-level routines
             (`average_psd`, `to_dB`, SPL, TOL, LTSA) consume a `PSDResult`.

             Four steps (DD-08 through DD-11):
             1. Square STFT magnitudes: `abs2.(spec.stft)` → full-scale²
             2. Single-sided correction (DD-10): interior frequency bins ×2;
                DC (row 1) and Nyquist (last row) unchanged.
             3. Normalise by `fs × window_energy` → full-scale²/Hz
                (Merchant 2015 Eq. 1, DD-09).
             4. Apply calibration via `apply_calibration!(psd_linear, freqs, cal)`
                → µPa²/Hz when calibration is available.

Arguments:
- `spec::SpectrogramResult`: Output of [`spectrogram`](@ref). Must have been
  produced with an even `nfft` (enforced by the assertion in `spectrogram`).
- `cal::Calibration = NoCalibration()`: Calibration to apply. Dispatches to the
  correct `apply_calibration!` method for `NoCalibration`, `ScalarCalibration`,
  or `TFCalibration`. When `NoCalibration()` (default), the PSD remains in
  full-scale²/Hz.

Returns:     [`PSDResult`](@ref) with `psd_linear` in µPa²/Hz (`is_calibrated =
             true`) or full-scale²/Hz (`is_calibrated = false`). `freqs` and
             `time` are the same `Vector` objects as in `spec` — no copy.

Constraints:
- `spec.nfft` must be even. This is guaranteed by `spectrogram` (DD-07) but is
  not re-checked here. Passing a `SpectrogramResult` constructed by any path
  other than `spectrogram` with an odd `nfft` will silently misapply the
  single-sided correction (the last row will be treated as Nyquist when it is
  not).
- `psd_linear` is always a freshly allocated matrix — mutations do not affect
  `spec`.
- Calibration is applied in-place on `psd_linear` after normalization. The
  original `spec.stft` is never modified.

Fails when:  Never — all failure modes are caught upstream in `spectrogram` and
             `apply_calibration!`.

Example:
```julia
signal = randn(Float64, 48000)
spec   = spectrogram(signal; fs = 48000.0, window_seconds = 1.0)
result = compute_psd(spec)                             # uncalibrated
result_cal = compute_psd(spec, RockhopperProfile().tf) # Rockhopper TF calibration
```

Do not use when:
- The dB form is needed directly — call `to_dB(compute_psd(spec, cal))`.
- Frame-averaged PSD is needed — call `average_psd(compute_psd(spec, cal))`.
- The input signal has not been calibrated in the time domain and a
  `ScalarCalibration` is available — both paths (time-domain `apply_calibration!`
  then uncalibrated `compute_psd`, and uncalibrated `compute_psd` then PSD-matrix
  `apply_calibration!`) give the same result, but the PSD-matrix path here is
  preferred for efficiency in chunked processing (one fewer FFT round-trip).

References:
Merchant et al. (2015) Measuring Acoustic Habitats. Methods in Ecology and
Evolution, 6, 257–265.

Ainslie MA, Miksis-Olds JL, Martin B, Heaney K, de Jong CAF, von
Benda-Beckmann AM, Lyons AP (2018) ADEON Underwater Soundscape and Modeling
Metadata Standard. Soundscape Specification deliverable v1.0.
Section 2.2.1 of the DPS defines the same normalisation using a
pre-normalised Hann window (explicit sqrt(8/3) factor visible in the formula);
that form is algebraically identical to the Merchant 2015 convention
implemented here: both yield `psd[k] = 2|X[k]|² / (fs × Σwᵢ²)` for interior
bins. Confirmed empirically: PAMGuide cross-validation Δmean < 0.02 dB.
"""
# ─── Private calibration resolver ────────────────────────────────────────────
#
# Purpose:  Determine which Calibration to apply at the PSD layer given an
#           Audiodata object. Three-step resolution, in priority order:
#           1. audio.is_calibrated == true  → signal already in physical units;
#              return NoCalibration (no further correction at PSD layer).
#           2. audio.calibration is not NoCalibration → an explicit calibration
#              was attached at I/O time (e.g. ScalarCalibration for SM3M); use it.
#           3. audio.calibration is NoCalibration → try get_profile(recorder) to
#              auto-resolve a typed recorder profile. If the profile carries a
#              TFCalibration in a field named :tf, return it.
#           4. No calibration found → warn and return NoCalibration.
#
# Constraints: Never throws. ArgumentError from get_profile is caught and treated
#              as "no profile registered for this recorder".
# Fails when:  Only re-throws non-ArgumentError exceptions from get_profile.
function _psd_calibration(audio::Audiodata) :: Calibration
    audio.is_calibrated && return NoCalibration()
    !(audio.calibration isa NoCalibration) && return audio.calibration

    # calibration is NoCalibration but signal is uncalibrated.
    # Try get_profile to auto-resolve via Val{recorder_id} dispatch (DD-13).
    recorder = audio.metadata.recorder
    try
        prof = get_profile(Symbol(recorder))
        if hasproperty(prof, :tf) && prof.tf isa TFCalibration
            return prof.tf
        end
    catch e
        e isa ArgumentError || rethrow(e)
    end

    @warn "compute_psd: no calibration found for recorder '$(recorder)'; " *
          "PSD will be in full-scale²/Hz. " *
          "Use get_profile(:$(recorder)).tf for recorders with typed profiles, " *
          "or supply calibration via apply_calibration!(audio) before calling " *
          "compute_psd."
    return NoCalibration()
end

function compute_psd(spec::SpectrogramResult,
                     cal::Calibration = NoCalibration()) :: PSDResult

    # Step 1 — squared magnitudes. abs2(z) = real(z)^2 + imag(z)^2, equivalent
    # to abs(z)^2 but avoids the intermediate sqrt. Shape: (n_freqs, n_frames).
    psd_linear = abs2.(spec.stft)

    # Step 2 — single-sided correction (DD-10). The rfft discards the negative-
    # frequency half of the DFT. Interior bins (neither DC nor Nyquist) each
    # represent two DFT bins of equal magnitude, so their power must be doubled.
    # DC (row 1, k=0) and Nyquist (last row, k=nfft÷2) are real-valued in the
    # two-sided DFT and appear only once; they are not doubled.
    # This is only correct for even nfft, enforced upstream by spectrogram (DD-07).
    n_freqs = size(psd_linear, 1)
    if n_freqs > 2
        # Indexing 2:end-1 in Julia is 1-based: skips row 1 (DC) and row n_freqs (Nyquist).
        psd_linear[2:end-1, :] .*= 2.0
    end

    # Step 3 — Merchant 2015 Eq. (1) normalisation.
    # Divides by fs (Hz) × window_energy (Σwᵢ², dimensionless but units of
    # samples²) to obtain units of power per Hz (full-scale²/Hz).
    # fs is Float32; Float64 conversion avoids precision loss in the product.
    psd_linear ./= Float64(spec.fs) * spec.window_energy

    # Step 4 — calibration. Dispatches on cal type:
    #   NoCalibration     → no-op; psd stays in full-scale²/Hz
    #   ScalarCalibration → uniform multiply by 10^(−S/10) → µPa²/Hz
    #   TFCalibration     → per-bin divide by interpolated tf_lin → µPa²/Hz
    apply_calibration!(psd_linear, spec.freqs, cal)

    return PSDResult(psd_linear, spec.freqs, spec.time, spec.fs,
                     spec.window_energy, spec.nfft, cal,
                     !(cal isa NoCalibration))
end

# ─── Convenience wrapper: Audiodata ──────────────────────────────────────────

"""
    compute_psd(audio::Audiodata; window_seconds, overlap_fraction=0.5,
                window=:hann, nfft=nothing, fft_plan=nothing) -> PSDResult

Purpose:     Convenience wrapper: compute PSD directly from an `Audiodata`
             object. Calls [`spectrogram`](@ref) internally and resolves
             calibration automatically via `_psd_calibration`:

             1. If `audio.is_calibrated` is `true`, no PSD-layer calibration
                is applied (`NoCalibration`). The PSD is in physical units
                because the signal was already calibrated in the time domain.
                `PSDResult.is_calibrated` is set to `true` to reflect this.
             2. If `audio.calibration` is not `NoCalibration`, it is applied
                at the PSD layer (e.g. `ScalarCalibration` for SM3M).
             3. Otherwise, `get_profile(Symbol(audio.metadata.recorder))` is
                tried. If the profile has a `:tf` field that is a
                `TFCalibration`, it is used (handles Rockhopper automatically).
             4. If no calibration is found, a warning is emitted and the PSD
                is returned in full-scale²/Hz.

Arguments:
- `audio::Audiodata`: The recording to analyse.
- `window_seconds::Real`: Analysis window duration. Forwarded to `spectrogram`.
- `overlap_fraction::Real = 0.5`: Frame overlap fraction. Forwarded to `spectrogram`.
- `window::Symbol = :hann`: Window function. Forwarded to `spectrogram`.
- `nfft::Union{Int,Nothing} = nothing`: FFT length. Forwarded to `spectrogram`.
- `fft_plan = nothing`: Pre-computed FFTW plan. When supplied, forwarded to
  `spectrogram` to avoid per-call plan construction overhead. Build with
  [`make_spectrogram_plan`](@ref). Must match the `nfft` that `spectrogram`
  would choose; a size mismatch throws `AssertionError` (DD-04).
- `cal::Union{Calibration,Nothing} = nothing`: Override calibration. When
  supplied, bypasses `_psd_calibration` and applies `cal` directly. When
  `nothing` (default), calibration is auto-resolved from `audio.metadata`.

Returns:     [`PSDResult`](@ref). Units depend on the resolved calibration —
             check `psd_units(result)` or `result.is_calibrated`.

Constraints:
- `audio.sig` must be non-empty and long enough for at least one frame.
- All other constraints from [`spectrogram`](@ref) apply.
- When `audio.is_calibrated = true`, `PSDResult.is_calibrated` is set to
  `true` even though `cal` is `NoCalibration` at the PSD layer. This ensures
  `psd_units` returns `:µPa²_per_Hz` for pre-calibrated signals.

Fails when:  Same conditions as [`spectrogram`](@ref).

Example:
```julia
audio = read_audio("recording.flac"; recorder="rockhopper")
result = compute_psd(audio; window_seconds = 1.0)   # auto-calibrates via Rockhopper profile
result_explicit = compute_psd(audio; window_seconds = 1.0, cal = NoCalibration())
```
"""
function compute_psd(audio::Audiodata;
                     window_seconds::Real,
                     overlap_fraction::Real          = 0.5,
                     window::Symbol                  = :hann,
                     nfft::Union{Int,Nothing}        = nothing,
                     fft_plan                        = nothing,
                     cal::Union{Calibration,Nothing} = nothing) :: PSDResult
    spec = spectrogram(audio.sig;
                       fs               = Float64(audio.fs),
                       window_seconds   = window_seconds,
                       overlap_fraction = overlap_fraction,
                       window           = window,
                       nfft             = nfft,
                       fft_plan         = fft_plan)
    resolved_cal = cal !== nothing ? cal : _psd_calibration(audio)
    result = compute_psd(spec, resolved_cal)
    # Propagate audio.is_calibrated: a pre-calibrated signal is in physical
    # units even if resolved_cal is NoCalibration (cascade step 1 above).
    # The primitive sets is_calibrated = !(cal isa NoCalibration), which misses
    # this case. Reconstruct only when the flag needs to change.
    is_cal = audio.is_calibrated || result.is_calibrated
    is_cal == result.is_calibrated && return result
    return PSDResult(result.psd_linear, result.freqs, result.time, result.fs,
                     result.window_energy, result.nfft, result.cal, is_cal)
end

# ─── Convenience wrapper: AbstractAudioSource ─────────────────────────────────

"""
    compute_psd(src, start, stop; gap_handling=:zero_fill,
                window_seconds, ...) -> PSDResult

Purpose:     Convenience wrapper: read a time window from an
             `AbstractAudioSource` and compute PSD. Equivalent to calling
             `read_audio_range(src, start, stop)` followed by the
             `Audiodata` overload of `compute_psd`. Calibration is resolved
             automatically from the resulting `Audiodata` (see
             `compute_psd(audio::Audiodata; ...)`).

Arguments:
- `src::AbstractAudioSource`: Source to read from.
- `start::DateTime`: Start of the time window (inclusive).
- `stop::DateTime`: End of the time window (exclusive).
- `gap_handling::Symbol = :zero_fill`: Gap handling passed to
  `read_audio_range`. `:zero_fill` inserts silence for gaps;
  `:error` raises on any gap.
- `window_seconds::Real`, `overlap_fraction`, `window`, `nfft`, `fft_plan`:
  Forwarded to `spectrogram` (see `compute_psd(audio::Audiodata; ...)`).

Returns:     [`PSDResult`](@ref) for the requested window.

Constraints:
- `stop > start` (enforced by `read_audio_range`).
- The window must contain at least one full analysis frame after gap handling.
- For large sources, prefer `process_chunks` with a `compute_psd` callback
  rather than loading the full time range at once.

Fails when:  `stop <= start`, or the resulting `Audiodata` is too short for one frame.

Example:
```julia
src    = IndexedFileSource(load_index("archive.arrow"), "archive/")
t0, t1 = time_range(src)
result = compute_psd(src, t0, t0 + Minute(5); window_seconds = 1.0)
```
"""
function compute_psd(src::AbstractAudioSource,
                     start::DateTime,
                     stop::DateTime;
                     gap_handling::Symbol            = :zero_fill,
                     window_seconds::Real,
                     overlap_fraction::Real          = 0.5,
                     window::Symbol                  = :hann,
                     nfft::Union{Int,Nothing}        = nothing,
                     fft_plan                        = nothing,
                     cal::Union{Calibration,Nothing} = nothing) :: PSDResult
    audio = read_audio_range(src, start, stop; gap_handling = gap_handling)
    return compute_psd(audio;
                       window_seconds   = window_seconds,
                       overlap_fraction = overlap_fraction,
                       window           = window,
                       nfft             = nfft,
                       fft_plan         = fft_plan,
                       cal              = cal)
end

"""
    compute_psd(src::AbstractAudioSource; window_seconds, ...) -> PSDResult

Purpose:     Compute PSD over the full time extent of `src`. Thin delegation
             to `compute_psd(src, start, stop; ...)` — see that method for all
             parameter semantics, calibration resolution, and constraints.

Arguments:   Same keyword arguments as the range form. `start` and `stop` are
             taken from `time_range(src)`.

Example:
```julia
src    = SingleFileSource("recording.flac"; recorder = "rockhopper")
result = compute_psd(src; window_seconds = 1.0)
```
"""
function compute_psd(src::AbstractAudioSource; kwargs...) :: PSDResult
    t_start, t_stop = time_range(src)
    return compute_psd(src, t_start, t_stop; kwargs...)
end

# ─── average_psd ──────────────────────────────────────────────────────────────

"""
    average_psd(result::PSDResult) -> Vector{Float64}

Purpose:     Average a multi-frame PSD matrix over time, returning a single
             spectral estimate. The mean is computed in linear power space —
             averaging in dB is not equivalent and is not done here. Convert
             to dB after averaging via `to_dB`.

Arguments:
- `result::PSDResult`: A computed PSD result, typically with multiple frames.

Returns:     `Vector{Float64}` of length `size(result.psd_linear, 1)`. Element
             `k` is the mean of `result.psd_linear[k, :]` across all frames.
             Units match `result` — check `psd_units(result)`.

Constraints:
- Averaging in linear power is the Merchant 2015 convention for LTSA and TOL
  band integration. Do not average in dB.
- For a single-frame PSDResult, the output equals `result.psd_linear[:, 1]`.

Fails when:  Never.

Example:
```julia
result = compute_psd(audio; window_seconds = 1.0)
avg    = average_psd(result)         # Vector{Float64}, length = n_freqs
db     = to_dB(avg)                  # dB re µPa²/Hz (if calibrated)
```
"""
function average_psd(result::PSDResult) :: Vector{Float64}
    # mean(m; dims=2) returns a (n_freqs, 1) Matrix; dropdims collapses it to
    # a (n_freqs,) Vector.
    return dropdims(mean(result.psd_linear; dims=2), dims=2)
end

# ─── to_dB ────────────────────────────────────────────────────────────────────

"""
    to_dB(result::PSDResult) -> Matrix{Float64}
    to_dB(v::AbstractVector{Float64}) -> Vector{Float64}

Purpose:     Convert linear PSD values to decibels: `10 × log10(x)`.
             Works on a full [`PSDResult`](@ref) (returns the same shape as
             `result.psd_linear`) or on a `Vector{Float64}` such as the
             output of `average_psd`.

             Units after conversion:
             - `psd_units(result) == :µPa²_per_Hz` → output in dB re 1 µPa²/Hz
             - `psd_units(result) == :fullscale²_per_Hz` → output in dBFS/Hz

Arguments:
- `result::PSDResult`: Source PSD to convert.
- `v::AbstractVector{Float64}`: Averaged PSD vector (from `average_psd`).

Returns:     `Matrix{Float64}` (same size as `psd_linear`) or `Vector{Float64}`
             (same length as `v`). All values in decibels.

Constraints:
- Input values must be strictly positive (> 0). `log10(0)` returns `-Inf`;
  negative input returns `NaN`. Both are left as-is — no clamping or warning.
- Do not call on dB values; `to_dB(to_dB(x)) ≠ x`.

Fails when:  Never (no assertions; behaviour on non-positive input is defined
             by IEEE 754 log10 semantics).

Example:
```julia
result  = compute_psd(audio; window_seconds = 1.0)
db_mat  = to_dB(result)               # Matrix{Float64}
db_mean = to_dB(average_psd(result))  # Vector{Float64}
```
"""
to_dB(result::PSDResult)              = 10.0 .* log10.(result.psd_linear)
to_dB(v::AbstractVector{Float64})     = 10.0 .* log10.(v)
