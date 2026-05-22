# Rockhopper Calibration

The Cornell Rockhopper is a small underwater acoustic recorder designed for
long-term passive monitoring. It stores audio as 16-bit FLAC at sample rates
up to 200 kHz. EcoAcoustics.jl ships a fleet-level frequency-dependent
calibration for the Rockhopper and applies it automatically.

---

## What is the Rockhopper calibration?

The Rockhopper has a hydrophone, a preamplifier, and an ADC. The combined
system response varies with frequency: sensitivity is not constant across the
recording bandwidth. To convert ADC units (full-scale counts) to physical
pressure units (µPa), you need the **transfer function** — the system
sensitivity at every frequency in dB.

EcoAcoustics.jl uses the `AnalogSensitivity_dB_re_1VperRefPress` column from
the Cornell `RH_calCurves` spreadsheet, converted to the package canonical
form (dB re full-scale per µPa). The calibration is loaded once at startup and
reused across all Rockhopper PSD computations.

---

## The calibration pipeline

```
ADC full-scale       →    TF correction    →    µPa²/Hz
(dimensionless)           (DD-12, DD-13)        (physical)
```

In detail:

1. The raw signal is stored in full-scale units (`[-1, +1]`) in `Audiodata.sig`.
2. The spectrogram is computed from the raw signal — no calibration yet (DD-08).
3. `compute_psd(audio)` calls `_psd_calibration(audio)`, which auto-resolves
   the Rockhopper TF via `get_profile(:rockhopper).tf` (DD-14).
4. `apply_calibration!(psd_linear, freqs, tf)` divides each PSD frequency bin
   by the interpolated TF value in linear power:
   `psd_calibrated[k] = psd_raw[k] / tf_lin[k]`
   where `tf_lin[k] = 10^(tf_dB[k] / 10)` and `tf_dB[k]` is the
   interpolated sensitivity at frequency `freqs[k]`.

---

## Units and sign convention

The canonical TF form in EcoAcoustics.jl is:

```
tf_dB = AnalogSensitivity_dB − 20·log10(Vmax_peak_V)
```

where `Vmax_peak_V = 5.0` V for the Rockhopper. This shift (~13.98 dB)
converts from "dB re 1 V per µPa at the ADC input" to "dB re ADC full-scale
per µPa." After the TF correction, PSD values are in µPa²/Hz and `psd_units`
returns `:µPa²_per_Hz`.

The sign convention is: `tf_dB` is conventionally negative (e.g. −160 dB
at low frequencies), meaning the hydrophone + preamplifier + ADC system
attenuates the physical pressure signal before digitisation. Dividing by a
small number (10^(tf_dB/10) ≪ 1) scales the PSD back up to physical units.

---

## Automatic calibration resolution

For any `Audiodata` with `recorder = "rockhopper"`, `compute_psd` resolves
the calibration automatically:

```julia
audio  = read_audio("RH428_20230101_120000.flac"; recorder = "rockhopper")
result = compute_psd(audio; window_seconds = 1.0)
# result.is_calibrated == true
# psd_units(result) == :µPa²_per_Hz
```

This works because `_psd_calibration` (DD-14) checks
`get_profile(:rockhopper).tf`, which returns the shipped `TFCalibration`.

If you have a per-unit calibration file (not the fleet-level default), load
it explicitly and pass it to `compute_psd`:

```julia
my_tf  = load_tf_calcurves("/path/to/RH428_calCurves.csv"; vmax_peak_V = 5.0)
spec   = spectrogram(audio.sig; fs = Float64(audio.fs), window_seconds = 1.0)
result = compute_psd(spec, my_tf)
```

---

## The shipped calibration file

The calibration data is stored at:

```
src/recorders/calibration_data/rockhopper_tf_calibration.csv
```

It was derived from the `RH_calCurves` sheet of `Rockhopper_TF.xlsx`
(Cornell University). The CSV contains 2798 frequency points from 0 Hz to
approximately 99 884 Hz with headers:

```
Frequency_Hz,SensorSensitivity_dB_re_1V_perRefPress,PreampGain_dB,
RecorderGain_dB,AnalogSensitivity_dB_re_1VperRefPress
```

Only `Frequency_Hz` and `AnalogSensitivity_dB_re_1VperRefPress` are used by
`load_tf_calcurves`. The other columns are retained in the CSV for provenance.

---

## What NOT to use

The `RH_calCurves` spreadsheet also has a column called `TF` (column 13).
This column folds in a PeakToRMS conversion factor of
`20·log10(2√2) ≈ 9.03 dB`. That factor is appropriate for converting a
peak-calibrated sine to RMS in the time domain, but it is **not** appropriate
for PSD, which operates on instantaneous squared samples. Using the `TF` column
instead of `AnalogSensitivity` will produce PSD values that are 9.03 dB too
high. The `load_tf_calcurves` function expects the `AnalogSensitivity` column;
it cannot detect which column you extracted from the column name alone.

---

## Cross-validation

PSD output is cross-validated against MANTA (a PAMGuide-compatible tool).
The reference data and expected tolerance are described in
`test/validation/rockhopper/README.md`. The test file is
`test/test_psd_manta_validation.jl`; it skips gracefully when reference data
is absent.
