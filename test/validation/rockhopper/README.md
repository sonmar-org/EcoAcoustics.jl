# Rockhopper PSD Validation — Reference Data

This directory holds reference data for cross-validating EcoAcoustics.jl PSD
output against MANTA (a PAMGuide-compatible tool).

---

## Shipped audio clips

Two 120-second Rockhopper FLAC clips are checked into the repository:

| File | Duration | Sample rate |
|------|----------|-------------|
| `rh_clip_004600Z.flac` | 120 s | 197 368 Hz |
| `rh_clip_184620Z.flac` | 120 s | 197 368 Hz |

The validation test (`test/test_psd_manta_validation.jl`) uses
`rh_clip_004600Z.flac` as its audio source.

---

## Required file (not shipped)

The MANTA reference CSV must be generated locally and placed here:

```
test/validation/rockhopper/rh_clip_004600Z_manta_psd.csv
```

Format: two columns, one header row:

```
Frequency_Hz,PSD_dB_re_uPa2_per_Hz
0.0,-42.1
1.0,-38.7
...
```

---

## Generating the reference with MANTA

1. Open MANTA. Load `rh_clip_004600Z.flac`.
2. Set the Rockhopper system sensitivity. Use the frequency-dependent TF
   from the Cornell calibration sheet (`RH_calCurves`, column
   `AnalogSensitivity_dB_re_1VperRefPress`). The EcoAcoustics.jl equivalent
   is `get_profile(:rockhopper).tf`.
3. Set analysis parameters to match the test:

   | Parameter | Value |
   |-----------|-------|
   | Window function | Hann |
   | Window length | 1 second (197 368 samples at 197 368 Hz) |
   | Overlap | 50 % |
   | Calibration | Frequency-dependent TF (Rockhopper) |

4. Compute the time-averaged PSD (linear average over all frames, then
   convert to dB re 1 µPa²/Hz). Export as CSV with the header row and
   columns named exactly `Frequency_Hz` and `PSD_dB_re_uPa2_per_Hz`.
5. Place the file at `rh_clip_004600Z_manta_psd.csv` in this directory.
6. Run `julia --project=. test/runtests.jl` — the validation testsets will
   now execute instead of skipping.

---

## Tolerance

EcoAcoustics.jl is expected to agree with MANTA within:

- Peak error < 0.1 dB across all bins
- Mean error < 0.05 dB

---

## Known sources of small differences

- **Frame-count edge effects.** 120 s ÷ 0.5 s hop = 240 non-overlapping
  hops. With a 1 s window and 50% overlap, there are 239 full frames.
  If MANTA counts frames differently, the last frame's contribution to the
  average may differ slightly.

- **TF interpolation.** EcoAcoustics.jl interpolates the TF linearly in dB
  between the 2798 calibration points. MANTA uses the same convention;
  differences vs. cubic-spline interpolation tools are expected to be
  < 0.02 dB.
