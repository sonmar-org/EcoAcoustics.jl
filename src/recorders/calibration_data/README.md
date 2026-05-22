# Shipped Calibration Data

This directory contains calibration files distributed with the EcoAcoustics.jl
package. Files here are version-controlled and pinned to the package version;
calibration updates are package version bumps (DD-13).

---

## rockhopper_tf_calibration.csv

**Provenance:** Derived from the Cornell University Rockhopper calibration
workbook `Rockhopper_TF.xlsx`, sheet `RH_calCurves`, columns 1–5:

| Column | Field |
|--------|-------|
| 1 | `Frequency_Hz` |
| 2 | `SensorSensitivity_dB_re_1V_perRefPress` |
| 3 | `PreampGain_dB` |
| 4 | `RecorderGain_dB` |
| 5 | `AnalogSensitivity_dB_re_1VperRefPress` |

The `AnalogSensitivity` column equals `SensorSensitivity + PreampGain +
RecorderGain` (columns 2+3+4). It is the total system sensitivity at the ADC
input in dB re 1 V per reference pressure, before ADC full-scale normalisation.

**Extraction:** One-time manual export from the Excel workbook. No numerical
transformation was applied during extraction. The CSV values are identical to
those in the spreadsheet.

**Canonical conversion:** `load_tf_calcurves` subtracts `20·log10(5.0) =
13.9794 dB` from `AnalogSensitivity` to convert to dB re full-scale per µPa
(the canonical TFCalibration form, DD-12). The `Vmax_peak_V = 5.0` V value is
the Rockhopper ADC full-scale peak voltage.

**Known deployments:** Deployment 139635MD01 (Mid-Atlantic Bight, 17 Sep 2023)
was empirically validated against MANTA published millidecade PSD output.
Results: max bin diff 0.05 dB, mean 0.001 dB, std 0.028 dB across 391 bins
(10–400 Hz) × 3 one-minute windows. See `test/test_psd_manta_validation.jl`
and `docs/src/calibration/rockhopper.md`.

**Validation test:** `test/test_psd_manta_validation.jl` asserts that loading
this file via `load_tf_calcurves` and loading the provenance copy at
`test/validation/rockhopper/Rockhopper_TF_calibration.csv` produce pointwise
identical `tf_dB` values (tolerance 1e-10). If the two files diverge, the test
fails, flagging a calibration drift that must be resolved before release.

**Do not modify** this file without a corresponding package version bump and an
update to the empirical validation test. Changing the calibration without
updating the version would silently invalidate any results computed with the
prior version.
