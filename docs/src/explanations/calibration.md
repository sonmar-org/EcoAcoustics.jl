# Calibration

This page explains the calibration system in EcoAcoustics.jl: what system sensitivity means physically, how the three calibration types work, when to use each method, and what the `is_calibrated` flag guarantees.

---

## Why calibration matters

A hydrophone recording is a time series of normalised numbers in [−1, 1], where 1.0 represents the ADC full-scale voltage. These numbers have no physical units until you apply the system's end-to-end sensitivity.

**System sensitivity** (in dB re 1 V/µPa) is the conversion factor from voltage to underwater sound pressure. It accounts for:
- The hydrophone element's sensitivity (how many volts per pascal of pressure)
- Any preamplifier or board gain
- The ADC full-scale voltage

For the SM3M at default gain, the total is −153 dB re 1 V/µPa. That means a normalised sample of 1.0 corresponds to approximately 4.47 × 10⁷ µPa of pressure (a very loud sound — this is the upper limit of what the recorder can capture without clipping). The underwater acoustic reference level is 1 µPa (not 1 Pa as in air), which is why the numbers are large.

**The sign is conventionally negative.** The system sensitivity formula is:

```
system_sensitivity_dB = hydrophone_sensitivity + preamp_gain + board_gain
                        + 20 · log₁₀(1 / Vadc_0pk)
```

Because hydrophone sensitivity is typically around −165 dB and gain adds back less than that, the total remains negative.

**SPL from a calibrated signal:**

```
dB_SPL = 20 · log₁₀(rms(signal)) − system_sensitivity_dB
```

Because `system_sensitivity_dB` is negative, subtracting it adds a large positive number — converting from dBFS (relative to ADC full scale) to dB re 1 µPa.

---

## The three calibration types

### `NoCalibration`

No calibration data is available for this recorder. `apply_calibration!` copies the signal unchanged, and `is_calibrated` stays `false`. Metric outputs will be in dBFS rather than dB re 1 µPa, and metric functions will warn you.

Use `NoCalibration` when you are exploring data or doing relative comparisons where absolute physical units don't matter.

### `ScalarCalibration`

One number — `system_sensitivity_dB` — applies uniformly across all frequencies. This is appropriate for recorders whose frequency response is flat over the band of interest (e.g. SM3M, LS1X).

```julia
cal = ScalarCalibration(-153.0f0)   # SM3M at default gain
```

Applying it multiplies every sample by `10^(-system_sensitivity_dB / 20)`:

| Sensitivity | Linear multiplier |
|---|---|
| −153 dB | ≈ 4.47 × 10⁷ |
| −165 dB | ≈ 1.78 × 10⁸ |

The GPU-accelerated path for `ScalarCalibration` uses a KernelAbstractions kernel that runs on whatever device the array lives on — pass a `CuArray` to use the GPU, a plain `Array` to use the CPU.

### `TFCalibration`

A frequency-dependent sensitivity curve, stored as (`freqs`, `tf_db`) arrays. Used when the recorder's sensitivity varies across frequency — for example, the Rockhopper hydrophone from Cornell, which is supplied with a per-unit TF curve.

```julia
cal = lookup_calibration("file.flac", "rockhopper", meta)
# returns TFCalibration with 2797 frequency points, 1–98000 Hz
```

The TF data is loaded from `rockhopper_TF.csv` at the time `lookup_calibration` is called and embedded directly in the `TFCalibration` object. No file path survives into the runtime object.

---

## Time-domain vs frequency-domain application

`TFCalibration` has two application methods because it is used in two different pipeline contexts.

### Time-domain: `apply_calibration!`

Use this when you need the calibrated waveform — for example, before listening to a clip, computing a spectrogram that will be displayed, or extracting click waveforms for further analysis.

```julia
out = similar(a.sig)
apply_calibration!(out, a.sig, cal; fs = a.fs)
```

Internally, this constructs an FIR filter from the TF data and convolves it with the signal:

1. Interpolate the TF onto an FFT frequency grid (DC to Nyquist, `fir_length ÷ 2 + 1` bins).
2. Convert dB → linear magnitude.
3. Inverse FFT to get the impulse response.
4. Circularly shift the impulse response to the center of the array so the window peak aligns with the main lobe.
5. Apply a Hann window to reduce sidelobe energy.
6. Convolve with the signal using `filtfilt` (default, zero-phase) or `filt` (causal).

**FIR length** controls frequency resolution and transient duration. The default of 512 samples gives ≈ 188 Hz per bin at 96 kHz — sufficient for smooth hydrophone TF curves. Increase to 1024 if the TF has narrow spectral features.

**Phase modes:**

| `phase` | Method | Group delay | Transient |
|---|---|---|---|
| `:zero` (default) | `filtfilt` | None | ≈ `fir_length − 1` samples at **both** ends |
| `:causal` | `filt` | `(fir_length − 1) ÷ 2` samples | ≈ `fir_length − 1` samples at **start** only |

At fs = 196,000 Hz (Rockhopper) with `fir_length = 512`, the causal group delay is ≈ 1.3 ms. At fs = 48,000 Hz (SM3M), it is ≈ 5.3 ms. This matters for click timing in v2 — use `:zero` (filtfilt) for any analysis that needs sample-accurate timestamps.

**Note on filtfilt design:** `filtfilt` applies the filter twice (forward + backward), giving an effective magnitude response of |H(f)|². To achieve the correct calibration |H(f)|² = TF_linear, the FIR is designed from `sqrt.(TF_linear)` for the `:zero` path. The `:causal` path uses TF_linear directly. Both paths produce the same calibrated output; they differ only in phase and edge behaviour.

### Frequency-domain: `apply_calibration_psd!`

Use this in metric pipelines (SPL, LTSA) where you are working directly with PSD values in dB. It is much faster than the time-domain method: one interpolation per frequency bin, no FFT.

```julia
apply_calibration_psd!(psd_dbfs, freqs_hz, cal)
# psd_dbfs is now in dB re 1 µPa²/Hz
```

Internally, this subtracts the interpolated TF value at each frequency bin from the PSD:

```
psd_calibrated[i] = psd[i] − tf_db(freqs[i])
```

Because `tf_db` is negative (e.g. −153 dB), subtracting it adds a large positive number — the PSD values increase, which is correct when converting from dBFS to µPa².

Frequencies outside the TF range (including freq = 0, the DC bin) are clamped to the nearest edge value without warning. This is expected: the DC bin is always present in an FFT output but absent from typical manufacturer TF data.

---

## The `is_calibrated` flag

`Audiodata.is_calibrated` is `false` at construction and set to `true` only by `apply_calibration`. It is a hard guarantee:

- Metric functions that require physical units (SPL, PSD, TOL, LTSA) check this flag and **warn** when it is `false`, labelling output as dBFS.
- Applying calibration twice throws `ArgumentError`. There is no "re-calibrate" path — create a new `Audiodata` from the raw file if you need to start over.
- `NoCalibration` does not set the flag to `true`. A recording with no calibration profile cannot be marked as calibrated.

---

## Typical workflow

```julia
# 1. Read the file (calibration looked up but not applied)
a = read_audio("T1-C__0__20170912_181500.wav")
@assert !a.is_calibrated
@assert a.calibration isa ScalarCalibration

# 2. Apply calibration — returns a new Audiodata
a_cal = apply_calibration(a)
@assert a_cal.is_calibrated

# 3. Now compute metrics on the calibrated recording
# (metric functions will not warn about missing calibration)
```

For archive-scale processing via `process_chunks`, apply calibration inside the chunk function:

```julia
compute_spl(source) do chunk
    chunk_cal = apply_calibration(chunk)
    (; spl_db = rms_to_spl(chunk_cal.sig))
end
```
