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
- `sens_db::Float32`:
    Total system sensitivity in dB re 1 V/µPa. Conventionally negative (e.g.
    −153 dB for an SM3M at default gain). Computed from the recorder's
    hydrophone sensitivity, preamplifier gain, board gain, and ADC full-scale
    peak voltage:

        sens_db = hydrophone_sensitivity + preamp_gain + board_gain
                  + 20·log₁₀(1 / Vadc_0pk)

    SPL conversion: dB_SPL = 20·log₁₀(rms_linear) − sens_db

Returned by `lookup_calibration` when a [`CalibrationProfile`](@ref) is
registered for the recorder. Recorders with frequency-dependent response (e.g.
Rockhopper) use [`TFCalibration`](@ref) instead.
"""
struct ScalarCalibration <: Calibration
    sens_db::Float32
end

"""
    TFCalibration <: Calibration

Frequency-dependent system sensitivity, stored as a transfer function (TF).

Used for recorders whose sensitivity varies significantly across the frequency
range of interest — for example, the Rockhopper hydrophone, whose TF is supplied
per unit by the manufacturer.

Fields
------
- `freqs::Vector{Float32}`:
    Frequencies in Hz at which the sensitivity is defined. Must be sorted
    ascending. Typically covers the recorder's usable bandwidth.
- `tf_db::Vector{Float32}`:
    System sensitivity in dB re 1 V/µPa at each frequency in `freqs`. Length
    must equal `length(freqs)`. Conventionally negative.

`TFCalibration` objects are self-contained: frequency and sensitivity arrays are
embedded at lookup time and no file path survives into the runtime object.
Frequency-dependent calibration is implemented in v2; this type is defined in v1
for architectural completeness.
"""
struct TFCalibration <: Calibration
    freqs::Vector{Float32}
    tf_db::Vector{Float32}
end
