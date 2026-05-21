"""
    Calibration

Abstract base type for all calibration representations.

Three concrete subtypes form the calibration hierarchy:

- [`NoCalibration`](@ref) — no calibration data available.
- [`ScalarCalibration`](@ref) — a single system sensitivity in dB re 1 V/µPa,
  valid uniformly across all frequencies.
- [`TFCalibration`](@ref) — a frequency-dependent transfer function in dB re 1 V/µPa.

Calibration objects are stored in [`Audiodata`](@ref) at I/O time and applied
later via `apply_calibration!`. The `is_calibrated` flag on `Audiodata` records
whether calibration has been applied; metric functions check this flag and warn
when computing on uncalibrated data.
"""
abstract type Calibration end

"""
    NoCalibration <: Calibration

Marker type indicating that no calibration data is available for this recording.

`lookup_calibration` returns `NoCalibration()` when no [`CalibrationProfile`](@ref)
is registered for the recorder. Metric outputs computed from uncalibrated data
are in units of dBFS (decibels relative to full-scale digital counts) rather
than dB re 1 µPa. The `is_calibrated` flag on [`Audiodata`](@ref) remains `false`.
"""
struct NoCalibration <: Calibration end

"""
    ScalarCalibration <: Calibration

A single scalar system sensitivity, valid uniformly across all frequencies.

Fields
------
- `system_sensitivity_dB::Float32`:
    Total system sensitivity in dB re 1 V/µPa. Conventionally negative (e.g.
    −153 dB for an SM3M at default gain). Computed from the recorder's
    hydrophone sensitivity, preamplifier gain, board gain, and ADC full-scale
    peak voltage:

        system_sensitivity_dB = hydrophone_sensitivity + preamp_gain + board_gain
                                 + 20·log₁₀(1 / Vadc_0pk)

    SPL conversion: dB_SPL = 20·log₁₀(rms_linear) − system_sensitivity_dB

Returned by `lookup_calibration` when a [`CalibrationProfile`](@ref) is
registered for the recorder. Recorders with frequency-dependent response (e.g.
Rockhopper) use [`TFCalibration`](@ref) instead.
"""
struct ScalarCalibration <: Calibration
    system_sensitivity_dB::Float32
end

"""
    TFCalibration <: Calibration

Frequency-dependent system sensitivity, stored as a transfer function (TF).

Used for recorders whose sensitivity varies significantly across the frequency
range of interest — for example, the Rockhopper hydrophone, whose TF is supplied
per unit by Cornell.

Fields
------
- `freqs::Vector{Float32}`:
    Frequencies in Hz at which the sensitivity is defined. Must be sorted
    strictly ascending. Typically covers the recorder's usable bandwidth.
    Does not include 0 Hz; DC bin requests are clamped to `freqs[1]`.
- `tf_db::Vector{Float32}`:
    System sensitivity in dB re 1 V/µPa at each frequency in `freqs`. Length
    must equal `length(freqs)`. Conventionally negative. Loaded by
    `lookup_calibration` via `_load_tf_csv`, which negates the raw CSV values
    (stored positive by manufacturer convention) on read.

`TFCalibration` objects are self-contained: frequency and sensitivity arrays are
embedded at lookup time and no file path survives into the runtime object.
"""
struct TFCalibration <: Calibration
    freqs::Vector{Float32}
    tf_db::Vector{Float32}
end

# ─── Private helpers ──────────────────────────────────────────────────────────

# Purpose:  Linearly interpolate a TF (in dB) at frequency f.
#           Clamped at both ends — handles freq = 0 (DC bin always absent
#           from recorded TF data) by returning cal_db[1].
# Constraints: cal_freqs must be sorted ascending. cal_db same length.
# Fails when: never (clamping prevents out-of-range access).
function _interp_tf(cal_freqs::AbstractVector{Float32},
                    cal_db::AbstractVector{Float32},
                    f::Real)::Float64
    f64 = clamp(Float64(f), Float64(cal_freqs[1]), Float64(cal_freqs[end]))
    idx = searchsortedfirst(cal_freqs, Float32(f64))
    idx == 1                   && return Float64(cal_db[1])
    idx > length(cal_freqs)    && return Float64(cal_db[end])
    f1 = Float64(cal_freqs[idx - 1]);  f2 = Float64(cal_freqs[idx])
    d1 = Float64(cal_db[idx - 1]);     d2 = Float64(cal_db[idx])
    t  = (f64 - f1) / (f2 - f1)
    return d1 + t * (d2 - d1)
end

"""
    _load_tf_csv(path) -> (freqs::Vector{Float32}, tf_db::Vector{Float32})

Purpose:     Read a two-column, no-header CSV file (freq_hz, sensitivity_dB)
             and return the frequency and sensitivity arrays with the sign
             convention used throughout this package (negative dB).

Arguments:
- `path::AbstractString`: Absolute path to the CSV file.

Returns:     `(freqs, tf_db)` where `freqs` is in Hz (ascending) and `tf_db`
             is in dB re 1 V/µPa (negative by convention). The second column
             of the CSV is negated on read: CSV values are positive (manufacturer
             convention); in-memory values are negative (package convention).

Constraints:
- File must have exactly two comma-separated columns per non-empty row.
- Frequencies must be strictly ascending.
- Frequency 0 Hz need not be present; DC bin interpolation is handled by
  clamping in `_interp_tf`.

Fails when:
- File not found or unreadable.
- A row does not have exactly two columns.
- A value cannot be parsed as `Float32`.
- Frequencies are not strictly ascending.
- File contains no data rows.
"""
function _load_tf_csv(path::AbstractString)
    freqs = Float32[]
    tf    = Float32[]
    open(path, "r") do io
        for (lineno, line) in enumerate(eachline(io))
            stripped = strip(line)
            isempty(stripped) && continue
            parts = split(stripped, ',')
            @assert length(parts) == 2 begin
                "$(basename(path)) line $lineno: expected 2 columns, got $(length(parts))"
            end
            f  = tryparse(Float32, strip(parts[1]))
            db = tryparse(Float32, strip(parts[2]))
            @assert f  !== nothing "_load_tf_csv: cannot parse frequency on line $lineno of $(basename(path))"
            @assert db !== nothing "_load_tf_csv: cannot parse dB value on line $lineno of $(basename(path))"
            push!(freqs, f)
            push!(tf, db)
        end
    end
    @assert !isempty(freqs) "_load_tf_csv: $path contains no data rows"
    @assert all(diff(freqs) .> 0) begin
        "_load_tf_csv: frequencies in $path must be strictly ascending"
    end
    return freqs, -tf   # negate: CSV is positive (Cornell convention), package is negative
end

# ─── apply_calibration! ───────────────────────────────────────────────────────

"""
    apply_calibration!(out, signal, ::NoCalibration)

Purpose:     Copy `signal` into `out` unchanged. Called when no calibration
             profile is available; preserves the in-place API contract without
             applying a spurious transformation.

Arguments:
- `out::AbstractVector`: Output buffer. Same length as `signal`.
- `signal::AbstractVector`: Input signal in normalised ADC units [−1, 1].
- `::NoCalibration`: Calibration marker (no data).

Returns:     `out`.

Constraints: `length(out)` must equal `length(signal)` (enforced by `copyto!`).
             Output is in the same units as input — normalised ADC, not µPa.

Fails when:  `length(out) ≠ length(signal)` (thrown by `copyto!`).

Example:
```julia
apply_calibration!(out, signal, NoCalibration())
```
"""
function apply_calibration!(out::AbstractVector,
                            signal::AbstractVector,
                            ::NoCalibration)
    copyto!(out, signal)
    return out
end

"""
    apply_calibration!(out, signal, cal::TFCalibration; fs, plans=nothing)

Purpose:     Apply frequency-dependent (transfer-function) calibration to a
             time-domain signal in-place, converting raw ADC amplitudes to
             calibrated physical-unit amplitudes (µPa). Calibration is applied
             in the frequency domain: the signal is FFT'd, each bin's amplitude
             is multiplied by the interpolated TF magnitude, and the result is
             inverse-FFT'd. No FIR filter is constructed.

Arguments:
- `out::AbstractVector{Float64}`: Output buffer, same length as `signal`.
  Modified in-place.
- `signal::AbstractVector{Float64}`: Input in normalised ADC units [−1, 1].
- `cal::TFCalibration`: Calibration object with `freqs` and `tf_db` arrays.
- `fs::Real`: Sampling rate in Hz. Must be > 0.
- `plans::Union{Nothing, Tuple} = nothing`: Pre-computed FFT plan tuple
  `(forward_plan, inverse_plan)`. When `nothing` (default), plans are
  computed internally on each call — correct for one-off use. For chunked
  processing, pre-compute once per chunk length and pass here to avoid
  repeated overhead. See the example below.

Returns:     `out`, modified in-place.

Constraints:
- `length(out)` must equal `length(signal)`.
- `fs` must be > 0.
- If `plans` is provided, both plans must have been created for arrays of
  `length(signal)`. A size mismatch raises at FFT apply time.
- TF interpolation clamps at `cal.freqs[1]` and `cal.freqs[end]`; DC (0 Hz)
  returns the value at `cal.freqs[1]`.
- Boundary behaviour is circular (inherent to FFT convolution). For signals
  of 0.1 s or longer, wrap-around is negligible for smooth TF curves. For
  very short signals, zero-pad if boundary effects matter.

Fails when:
- `length(out) ≠ length(signal)`.
- `fs ≤ 0`.

Example:
```julia
# One-off use — plans computed automatically:
out = similar(a.sig)
apply_calibration!(out, a.sig, cal; fs = a.fs)

# Chunked use — pre-compute plans once for the chunk length, then re-use:
buf      = zeros(chunk_samples)
fwd      = FFTW.plan_rfft(buf)
plans    = (fwd, inv(fwd))
apply_calibration!(out, chunk.sig, cal; fs = chunk.fs, plans = plans)
```

Do not use when: Only PSD-based metrics are needed — use `apply_calibration_psd!`
instead (no FFT overhead, simpler).
"""
function apply_calibration!(out::AbstractVector{Float64},
                            signal::AbstractVector{Float64},
                            cal::TFCalibration;
                            fs::Real,
                            plans::Union{Nothing, Tuple} = nothing)
    @assert length(out) == length(signal) begin
        "apply_calibration!: out and signal must have the same length"
    end
    @assert fs > 0 "apply_calibration!: fs must be > 0, got $fs"

    # An FFT plan caches the algorithm FFTW chose for this signal length.
    # Re-using a pre-computed plan across many chunks avoids that selection
    # overhead on every call. fwd_plan transforms signal → spectrum;
    # inv_plan transforms spectrum → signal.
    fwd_plan, inv_plan = if isnothing(plans)
        P = FFTW.plan_rfft(signal)
        (P, inv(P))
    else
        plans
    end

    # rfft returns the single-sided complex spectrum: length(signal)÷2 + 1 bins,
    # covering 0 Hz (DC) through fs/2 (Nyquist). The range below generates exactly
    # those bin-centre frequencies.
    X     = fwd_plan * signal
    freqs = range(0.0, Float64(fs) / 2; length = length(X))

    # Interpolate the TF sensitivity (negative dB by convention) onto each FFT bin,
    # then convert to a linear amplitude multiplier. Because tf_db is negative
    # (e.g. −153 dB), negating it gives a large positive exponent, and tf_mag >> 1:
    # tf_db = −153 dB  →  −tf_db/20 = 7.65  →  tf_mag ≈ 4.47×10⁷  (ADC → µPa).
    # This matches ScalarCalibration: scale = 10^(−system_sensitivity_dB / 20).
    tf_db_interp = [_interp_tf(cal.freqs, cal.tf_db, f) for f in freqs]
    X .*= 10 .^ (-tf_db_interp ./ 20)

    # irfft requires the explicit output length because rfft output of size n÷2+1
    # does not uniquely identify whether the original signal had n or n+1 samples.
    out .= inv_plan * X
    return out
end

"""
    apply_calibration_psd!(psd, freqs, cal::TFCalibration)

Purpose:     Apply frequency-dependent calibration to a power spectral density
             vector in-place, converting from dBFS/Hz to dB re 1 µPa²/Hz.
             This is the metric-pipeline path; use the time-domain method when
             the calibrated waveform itself is needed.

Arguments:
- `psd::AbstractVector`: PSD values in dB (dBFS/Hz). Modified in-place.
- `freqs::AbstractVector`: Frequency axis in Hz. Same length as `psd`.
- `cal::TFCalibration`: Calibration TF. `tf_db` must follow the negative
  convention (conventionally negative dB re 1 V/µPa). Interpolated onto
  `freqs` with edge clamping.

Returns:     `psd`, modified in-place.

Constraints:
- `length(psd)` must equal `length(freqs)`.
- Calibration is subtracted in dB: `psd[i] -= tf_db(freqs[i])`. Because
  `tf_db` is negative, the PSD values increase after subtraction — correct,
  as this converts dBFS to physical units. Example: if tf_db = −153 dB at
  frequency f, then psd[i] increases by 153 dB.
- Frequencies outside `cal.freqs` bounds (including freq = 0) are clamped to
  the nearest edge value and noted in docs; no warning is issued per call.

Fails when:  `length(psd) ≠ length(freqs)`.

Example:
```julia
apply_calibration_psd!(psd_dbfs, freqs_hz, cal)   # now psd is in dB re 1 µPa²/Hz
```

Do not use when: The calibrated time-domain waveform is needed — use
`apply_calibration!` with `cal::TFCalibration` instead.
"""
function apply_calibration_psd!(psd::AbstractVector,
                                freqs::AbstractVector,
                                cal::TFCalibration)
    @assert length(psd) == length(freqs) begin
        "apply_calibration_psd!: psd and freqs must have the same length, " *
        "got $(length(psd)) and $(length(freqs))"
    end
    for i in eachindex(psd)
        psd[i] -= _interp_tf(cal.freqs, cal.tf_db, freqs[i])
    end
    return psd
end

# ─── Out-of-place wrappers ────────────────────────────────────────────────────

"""
    apply_calibration(signal, cal; kwargs...) -> Vector{Float64}

Purpose:     Out-of-place version of `apply_calibration!` for time-domain
             signals. Allocates a new output buffer and returns it.

Arguments:
- `signal::AbstractVector{Float64}`: Input in normalised ADC units.
- `cal::Calibration`: Calibration object.
- `kwargs...`: Forwarded to `apply_calibration!`. For `TFCalibration`,
  `fs` is required; `plans` is optional (see `apply_calibration!` for details).
  Ignored for `NoCalibration` (no kwargs apply).

Returns:     `Vector{Float64}` same length as `signal`.

Constraints: For `TFCalibration`, `fs` must be provided.

Fails when:  Same conditions as the matching `apply_calibration!` method.

Example:
```julia
sig_cal = apply_calibration(signal, cal; fs = Float64(a.fs))
```
"""
function apply_calibration(signal::AbstractVector{Float64}, cal::Calibration; kwargs...)
    out = similar(signal)
    if cal isa TFCalibration
        apply_calibration!(out, signal, cal; kwargs...)
    else
        apply_calibration!(out, signal, cal)
    end
    return out
end
