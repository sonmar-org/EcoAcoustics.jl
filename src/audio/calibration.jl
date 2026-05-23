"""
    Calibration

Abstract base type for all calibration representations.

Three concrete subtypes form the calibration hierarchy:

- [`NoCalibration`](@ref) — no calibration data available.
- [`ScalarCalibration`](@ref) — a single system sensitivity in dB re 1 V/µPa,
  valid uniformly across all frequencies.
- [`TFCalibration`](@ref) — a frequency-dependent transfer function in dB re
  full-scale per µPa (Raven Workbench canonical form).

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

Used for recorders whose sensitivity varies across the frequency range of
interest — for example, the Rockhopper hydrophone, whose TF is supplied per
unit by Cornell.

Fields
------
- `frequency::Vector{Float64}`:
    Frequencies in Hz at which the sensitivity is defined. Strictly ascending,
    no duplicate values, length ≥ 2.
- `tf_dB::Vector{Float64}`:
    System sensitivity in dB re full-scale per µPa at each frequency in
    `frequency`. This is the Raven Workbench canonical form (DD-12): values are
    in dB relative to the ADC full-scale amplitude per µPa. Conventionally
    negative (e.g. −229.7 dB at DC for the Rockhopper at 5 V full-scale peak).
    `length(tf_dB)` must equal `length(frequency)`.
- `source::String`:
    Path to the file from which this TF was loaded. Used for provenance and
    reproducibility; allows `test_psd_manta_validation.jl` to verify that the
    package-shipped TF matches the validation copy pointwise.
- `format::Symbol`:
    Identifies the file format the TF was loaded from.
    `:rockhopper_calcurves_csv` for files loaded via `load_tf_calcurves`.
    Other symbols defined per future recorder loaders.
- `conversion_notes::String`:
    Human-readable description of any conversion applied when loading raw file
    data into canonical form. For `load_tf_calcurves`, records the dB shift
    used to convert from AnalogSensitivity to dB re full-scale per µPa.

Calibration application
-----------------------
To convert full-scale PSD (fs²/Hz) to physical PSD (µPa²/Hz):

    PSD[µPa²/Hz] = PSD[fs²/Hz] × 10^(−tf_dB(f) / 10)

The linear-power sensitivity `tf_lin(f) = 10^(tf_dB(f) / 10)` is interpolated
at each target frequency using linear interpolation in linear-power space (DD-11;
see `_interp_linear_power`). Constant extrapolation applies outside `frequency`.

See [`load_tf_calcurves`](@ref) for the public loader that constructs this type
from a Cornell-format CSV. See DD-12 in `docs/design_decisions.md` for the
canonical-format rationale.

`TFCalibration` objects are self-contained: no file needs to remain accessible
at runtime.
"""
struct TFCalibration <: Calibration
    frequency::Vector{Float64}
    tf_dB::Vector{Float64}
    source::String
    format::Symbol
    conversion_notes::String

    function TFCalibration(frequency::Vector{Float64},
                           tf_dB::Vector{Float64},
                           source::String,
                           format::Symbol,
                           conversion_notes::String)
        length(frequency) >= 2 || throw(ArgumentError(
            "TFCalibration: frequency must have length ≥ 2, got $(length(frequency))"))
        length(frequency) == length(tf_dB) || throw(ArgumentError(
            "TFCalibration: frequency and tf_dB must have the same length; " *
            "got $(length(frequency)) and $(length(tf_dB))"))
        # diff(v) produces [v[2]-v[1], v[3]-v[2], ...]; all(> 0) requires every
        # consecutive pair to be strictly increasing — rules out both non-monotonic
        # sequences and duplicate values in one check.
        all(diff(frequency) .> 0) || throw(ArgumentError(
            "TFCalibration: frequency must be strictly ascending with no " *
            "duplicates. Found a non-ascending or repeated value."))
        new(frequency, tf_dB, source, format, conversion_notes)
    end
end

# ─── Private helpers ──────────────────────────────────────────────────────────

# Purpose:  Linearly interpolate an arbitrary y-series at x, given sorted x-series xs.
#           Clamps at both ends — queries below xs[1] return ys[1]; queries above
#           xs[end] return ys[end] (constant extrapolation).
# Constraints: xs must be sorted ascending. Same length as ys. Length ≥ 2.
# Fails when: never (clamping prevents out-of-range access).
function _interp(xs::AbstractVector{<:Real},
                 ys::AbstractVector{<:Real},
                 x::Real)::Float64
    x64 = clamp(Float64(x), Float64(xs[1]), Float64(xs[end]))
    # searchsortedfirst returns the first index i such that xs[i] >= x64.
    # After clamping, x64 ∈ [xs[1], xs[end]], so idx ∈ [1, length(xs)].
    idx = searchsortedfirst(xs, x64)
    idx == 1                 && return Float64(ys[1])
    idx > length(xs)         && return Float64(ys[end])
    x1 = Float64(xs[idx - 1]);  x2 = Float64(xs[idx])
    y1 = Float64(ys[idx - 1]);  y2 = Float64(ys[idx])
    t  = (x64 - x1) / (x2 - x1)   # fraction of the way from xs[idx-1] to xs[idx]
    return y1 + t * (y2 - y1)
end

# Purpose:  Interpolate TF sensitivity in dB at frequency f.
#           Used by the time-domain apply_calibration! path. Interpolates in
#           dB space (linear-in-dB), which is correct for the FFT-based
#           time-domain calibration where we want the dB value to convert to
#           a linear amplitude multiplier.
# Constraints: cal_freqs strictly ascending; cal_db same length.
# Fails when: never (clamping via _interp prevents out-of-range).
_interp_tf(cal_freqs::AbstractVector{<:Real},
           cal_db::AbstractVector{<:Real},
           f::Real)::Float64 = _interp(cal_freqs, cal_db, f)

# Purpose:  Interpolate TF sensitivity in linear power at frequency f.
#           Used by the PSD-layer apply_calibration! path. Interpolates in
#           linear power space per DD-11 (Raven Workbench convention): the TF
#           in dB is first converted to linear power (10^(tf_dB/10)), then
#           linearly interpolated, then the caller divides by the result.
#           The caller pre-computes tf_lin = 10 .^ (tf_dB ./ 10) once and
#           passes it here to avoid repeated exponentiation.
# Constraints: cal_freqs strictly ascending; tf_lin same length.
# Fails when: never (clamping via _interp prevents out-of-range).
_interp_linear_power(cal_freqs::AbstractVector{<:Real},
                     tf_lin::AbstractVector{<:Real},
                     f::Real)::Float64 = _interp(cal_freqs, tf_lin, f)

# Purpose:  Read a two-column, no-header CSV file (freq_hz, sensitivity_dB)
#           and return the frequency and sensitivity arrays with the sign
#           convention used throughout this package (negative dB).
#           Legacy helper; superseded by load_tf_calcurves for new recorders.
#
# Arguments:
# - `path::AbstractString`: Absolute path to the CSV file.
#
# Returns:     `(freqs, tf_db)` where `freqs` is in Hz (ascending, Float32) and
#              `tf_db` is in dB re 1 V/µPa (negative by convention, Float32).
#              The second column of the CSV is negated on read: CSV values are
#              positive (manufacturer convention); in-memory values are negative
#              (package convention).
#
# Constraints:
# - File must have exactly two comma-separated columns per non-empty row.
# - Frequencies must be strictly ascending.
# - File contains no header row (all rows are parsed as data).
#
# Fails when:
# - File not found or unreadable.
# - A row does not have exactly two columns.
# - A value cannot be parsed as Float32.
# - Frequencies are not strictly ascending.
# - File contains no data rows.
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

# ─── apply_calibration! — time domain ────────────────────────────────────────

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

             This is the time-domain path. For PSD-only calibration use
             `apply_calibration!(psd_linear, freqs, cal::TFCalibration)` instead.

Arguments:
- `out::AbstractVector{Float64}`: Output buffer, same length as `signal`.
  Modified in-place.
- `signal::AbstractVector{Float64}`: Input in normalised ADC units [−1, 1].
- `cal::TFCalibration`: Calibration object with `frequency` and `tf_dB` arrays.
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
- TF interpolation clamps at `cal.frequency[1]` and `cal.frequency[end]`;
  DC (0 Hz) returns the value at `cal.frequency[1]`.
- Boundary behaviour is circular (inherent to FFT convolution). For signals
  of 0.1 s or longer, wrap-around is negligible for smooth TF curves.

Fails when:
- `length(out) ≠ length(signal)`.
- `fs ≤ 0`.

Example:
```julia
# One-off use — plans computed automatically:
out = similar(a.sig)
apply_calibration!(out, a.sig, cal; fs = Float64(a.fs))

# Chunked use — pre-compute plans once for the chunk length, then re-use:
buf      = zeros(chunk_samples)
fwd      = FFTW.plan_rfft(buf)
plans    = (fwd, inv(fwd))
apply_calibration!(out, chunk.sig, cal; fs = Float64(chunk.fs), plans = plans)
```

Do not use when: Only PSD-based metrics are needed — use the PSD-matrix form
`apply_calibration!(psd_linear, freqs, cal)` instead (no FFT overhead, simpler).
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

    # GPU-blocker: this interpolation loop uses _interp_tf / searchsortedfirst
    # which are not available inside GPU device code. Before this
    # apply_calibration! method can be retargeted to operate on CuArray or
    # other GPU arrays, the loop must be replaced with a KernelAbstractions.jl
    # @kernel that implements sorted binary search on the device side.
    # The NoCalibration and ScalarCalibration methods do not have this
    # constraint.
    # Interpolate the TF sensitivity (in dB) onto each FFT bin using dB-linear
    # interpolation, then convert to a linear amplitude multiplier.
    # Because tf_dB is negative (Raven canonical form, e.g. −229 dB), negating
    # it gives a large positive exponent, and tf_mag >> 1:
    # tf_dB = −229 dB  →  −tf_dB/20 = 11.45  →  tf_mag ≈ 2.8×10¹¹  (ADC → µPa).
    tf_db_interp = [_interp_tf(cal.frequency, cal.tf_dB, f) for f in freqs]
    X .*= 10 .^ (-tf_db_interp ./ 20)

    # irfft requires the explicit output length because rfft output of size n÷2+1
    # does not uniquely identify whether the original signal had n or n+1 samples.
    out .= inv_plan * X
    return out
end

# ─── apply_calibration! — PSD matrix (linear power domain) ───────────────────
#
# These three methods handle the PSD-layer calibration path. First arg is a
# Matrix{Float64} of linear PSD values (full-scale²/Hz), not a time-domain
# signal. Julia dispatches the correct method based on the first argument type.
#
# These methods are called from compute_psd after the single-sided correction
# and window-energy normalisation have been applied. All three return the matrix
# with output units as documented.

"""
    apply_calibration!(psd_linear, freqs, ::NoCalibration)

Purpose:     No-op for the PSD matrix calibration path. `psd_linear` is
             returned unchanged. Units remain full-scale²/Hz, not µPa²/Hz.
             Called by `compute_psd` when no calibration is available.

Arguments:
- `psd_linear::AbstractMatrix{Float64}`: PSD in linear units. Not modified.
- `freqs::Vector{Float64}`: Bin-centre frequencies. Not used.
- `::NoCalibration`: Calibration marker (no data).

Returns:     `psd_linear` unchanged.

Constraints: None. Both arguments are accepted at any size.

Fails when:  Never.

Example:
```julia
apply_calibration!(psd_matrix, freqs, NoCalibration())   # no-op
```
"""
function apply_calibration!(psd_linear::AbstractMatrix{Float64},
                            freqs::Vector{Float64},
                            ::NoCalibration)
    return psd_linear
end

"""
    apply_calibration!(psd_linear, freqs, cal::ScalarCalibration)

Purpose:     Apply scalar (frequency-independent) calibration to a linear PSD
             matrix in-place, converting from full-scale²/Hz to µPa²/Hz.
             Every element is multiplied by `10^(−system_sensitivity_dB / 10)`.
             Because `system_sensitivity_dB` is conventionally negative (e.g.
             −153 dB), the multiplier is large (e.g. ≈ 2×10¹⁵ for SM3M).

Arguments:
- `psd_linear::AbstractMatrix{Float64}`: PSD in full-scale²/Hz. Modified in-place.
- `freqs::Vector{Float64}`: Bin-centre frequencies in Hz. Not used directly;
  accepted for dispatch consistency with the TFCalibration method.
- `cal::ScalarCalibration`: Holds `system_sensitivity_dB` in dB re full-scale
  per µPa (conventionally negative).

Returns:     `psd_linear`, modified in-place. Units after return: µPa²/Hz.

Constraints: The same scale factor is applied to all frequency bins (scalar).
             For frequency-dependent calibration use the TFCalibration method.

Fails when:  Never (no dimension constraints — scalar broadcast always succeeds).

Example:
```julia
cal = ScalarCalibration(-153.0f0)      # SM3M
apply_calibration!(psd_matrix, freqs, cal)
```
"""
function apply_calibration!(psd_linear::AbstractMatrix{Float64},
                            freqs::Vector{Float64},
                            cal::ScalarCalibration)
    # Power conversion: amplitude factor is 10^(-S/20) (as in the KA kernel),
    # so power factor is 10^(-S/10) = (10^(-S/20))^2.
    factor = 10.0 ^ (-Float64(cal.system_sensitivity_dB) / 10.0)
    psd_linear .*= factor
    return psd_linear
end

"""
    apply_calibration!(psd_linear, freqs, cal::TFCalibration)

Purpose:     Apply frequency-dependent TF calibration to a linear PSD matrix
             in-place, converting from full-scale²/Hz to µPa²/Hz. Each frequency
             bin (row) is divided by the TF's linear power sensitivity at that
             frequency, interpolated using linear-power-scale interpolation
             (DD-11). This is the Raven Workbench convention.

Arguments:
- `psd_linear::Matrix{Float64}`: PSD in full-scale²/Hz, size (n_freqs, n_frames).
  Modified in-place. Each row is one frequency bin; each column is one frame.
- `freqs::Vector{Float64}`: Bin-centre frequencies in Hz. Length must equal
  `size(psd_linear, 1)`. Typically `SpectrogramResult.freqs`.
- `cal::TFCalibration`: Calibration TF. `tf_dB` values in dB re full-scale per
  µPa (canonical form, DD-12). Interpolated at each bin frequency using
  linear-power-scale interpolation (DD-11); constant extrapolation at edges.

Returns:     `psd_linear`, modified in-place. Units after return: µPa²/Hz.

Constraints:
- `size(psd_linear, 1)` must equal `length(freqs)`.
- Operates on linear (not dB) PSD values. Passing dB PSD values will produce
  incorrect results without error.
- Extrapolation beyond `cal.frequency` bounds uses the nearest endpoint value.

Fails when:  `size(psd_linear, 1) ≠ length(freqs)`.

Example:
```julia
apply_calibration!(psd_matrix, spec.freqs, rockhopper_tf)
```

Do not use when: The PSD is already in dB. Work on linear PSD and convert to
dB after calibration via `to_dB(psd_result)`.
"""
function apply_calibration!(psd_linear::Matrix{Float64},
                            freqs::Vector{Float64},
                            cal::TFCalibration)
    @assert size(psd_linear, 1) == length(freqs) begin
        "apply_calibration!: psd_linear has $(size(psd_linear, 1)) rows " *
        "but freqs has $(length(freqs)) entries"
    end

    # Pre-compute TF in linear power once: tf_lin[i] = 10^(tf_dB[i] / 10).
    # This is the denominator we divide each bin's PSD by. Computing it here
    # avoids repeated exponentiation inside the per-bin loop.
    tf_lin = 10 .^ (cal.tf_dB ./ 10)

    # GPU-blocker: this interpolation loop uses _interp_linear_power / searchsortedfirst
    # which are not available inside GPU device code. Before this
    # apply_calibration! method can be retargeted to operate on CuArray or
    # other GPU arrays, the loop must be replaced with a KernelAbstractions.jl
    # @kernel that implements sorted binary search on the device side.
    # The NoCalibration and ScalarCalibration methods do not have this
    # constraint.
    # For each frequency bin (row), divide all frames by the interpolated TF.
    # @views ensures psd_linear[k, :] is a zero-copy slice, so ./ is in-place
    # on the matrix without allocating a temporary row vector.
    for k in axes(psd_linear, 1)
        factor = _interp_linear_power(cal.frequency, tf_lin, freqs[k])
        @views psd_linear[k, :] ./= factor
    end
    return psd_linear
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
