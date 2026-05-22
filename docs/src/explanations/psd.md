# Power Spectral Density

The **power spectral density** (PSD) describes how acoustic power is
distributed across frequency. It is the primary output of the
EcoAcoustics.jl soundscape pipeline and the foundation for SPL, TOL, LTSA,
and soundscape indices. Every metric above PSD in the pipeline is a
transformation of a `PSDResult`.

---

## Position in the pipeline

```
spectrogram → complex STFT matrix  [SpectrogramResult]
      │
      ▼
compute_psd → power per Hz bin  [PSDResult, µPa²/Hz after calibration]
      │
      ├──→ SPL / TOL       (task 10)  →  broadband or third-octave level  [dB re 1 µPa]
      │
      └──→ LTSA            (task 11)  →  long-term spectral average  [dB re 1 µPa²/Hz]
                │
                └──→ Soundscape indices  (task 12)  →  ACI, entropy
```

`spectrogram` is unit-agnostic — it produces complex amplitudes in full-scale
units and applies no calibration. All physical unit conversion happens at the
PSD layer (DD-08).

---

## Normalization: Merchant et al. (2015) Equation 1

EcoAcoustics.jl uses the normalization defined in Merchant et al. (2015)
*Measuring Acoustic Habitats*, Methods in Ecology and Evolution, 6, 257–265.
This is also the convention used by PAMGuide, MANTA, and Triton.

For each frequency bin `k` and time frame `j`, the PSD before calibration is:

```
psd[k, j] = correction_k × |STFT[k, j]|² / (fs × W)
```

where:
- `|STFT[k, j]|²` is the squared magnitude of the complex STFT at bin `k`, frame `j`
  (`abs2` in Julia).
- `fs` is the sample rate in Hz.
- `W = Σwᵢ²` is the window energy — the sum of squared window coefficients
  (`SpectrogramResult.window_energy`). For a Hann window of length N,
  `W ≈ 3N/8`.
- `correction_k` is the single-sided correction factor (see below).

The result is in units of **full-scale² per Hz** before calibration, and
**µPa² per Hz** after calibration.

---

## Single-sided correction (DD-10)

`spectrogram` uses Julia's `rfft`, which produces only the non-negative
frequency bins of the DFT. For a real-valued signal of length N with even nfft:
- Bin 0 (DC, frequency 0 Hz) and bin N/2 (Nyquist) are **real-valued** in the
  two-sided DFT — they have no negative-frequency counterpart. Their power is
  complete as-is.
- Every interior bin `k = 1, …, N/2 − 1` represents two DFT bins of equal
  magnitude (the positive and negative frequency pair). The `rfft` discards the
  negative-frequency half, so the power must be **doubled** to account for the
  missing contribution.

In code:

```julia
psd_linear[2:end-1, :] .*= 2.0   # interior bins ×2
# DC (row 1) and Nyquist (last row) are not doubled
```

This requires `nfft` to be even. `spectrogram` enforces this with an
`@assert iseven(nfft)` (DD-07).

---

## The even-nfft constraint (DD-07)

The single-sided correction treats `psd_linear[end, :]` as the Nyquist bin.
This is only correct when `nfft` is even: for even nfft, `rfft` of length N
produces N/2+1 bins (DC + N/2−1 interior + Nyquist). For odd nfft, the last
bin is NOT Nyquist and would be incorrectly excluded from the ×2 correction.

`spectrogram` rejects odd nfft at the source:

```julia
@assert iseven(nfft_actual) "spectrogram: nfft must be even (DD-07)..."
```

If your window length is odd, pass an explicit even `nfft` (e.g.
`nfft = window_length + 1`).

---

## Calibration at the PSD layer

The spectrogram is calibration-agnostic. Calibration is applied at the PSD
layer because that is where we know the physical quantity (power per Hz). This
preserves the complex STFT for future click/whistle detectors, which need
phase information but not physical calibration.

`compute_psd` dispatches on the calibration type:

| Calibration type | Effect on `psd_linear` |
|---|---|
| `NoCalibration` | No change; output in full-scale²/Hz |
| `ScalarCalibration(S)` | Uniform multiply by `10^(-S/10)` → µPa²/Hz |
| `TFCalibration` | Per-bin divide by interpolated `tf_lin[k]` → µPa²/Hz |

`psd_units(result)` returns `:µPa²_per_Hz` when `result.is_calibrated` is
`true`, and `:fullscale²_per_Hz` otherwise (DD-15).

---

## Averaging and the dB step

**Always average in linear power, not in dB.** The Merchant 2015 LTSA
convention averages spectral estimates in linear power (µPa²/Hz), then
converts to dB for display or band integration.

```julia
result = compute_psd(audio; window_seconds = 1.0)
avg    = average_psd(result)       # Vector{Float64}, µPa²/Hz
db     = to_dB(avg)                # dB re 1 µPa²/Hz
```

Averaging in dB is mathematically equivalent only when the PSD has no temporal
variation. In practice, bioacoustic data always has variation, and dB averaging
underestimates power when the dynamic range is high — see DD-16.

`to_dB` is a thin wrapper: `10 .* log10.(x)`. It works on a `PSDResult`
(returns a `Matrix{Float64}`) or on a `Vector{Float64}` from `average_psd`.

---

## Pre-calibrated signals

A signal can be calibrated in the time domain before spectrogram computation:

```julia
audio_cal = apply_calibration(audio)    # audio_cal.is_calibrated == true
result    = compute_psd(audio_cal; window_seconds = 1.0)
```

In this case, `compute_psd` detects `audio_cal.is_calibrated == true` and
applies `NoCalibration()` at the PSD layer (DD-14). The PSD values are
physically correct (in µPa²/Hz) because the signal was already calibrated, but
`result.is_calibrated` will be `false` — it records whether PSD-layer
calibration was applied, not whether the signal is in physical units. This
known semantic gap is documented in the `compute_psd(Audiodata)` docstring.

---

## Chunked processing at archive scale

For archive-scale recordings, use `process_chunks` with a `compute_psd`
callback rather than loading the full time range into memory:

```julia
src    = IndexedFileSource(load_index("archive.arrow"), "archive/")
t0, t1 = time_range(src)

results = process_chunks(src, chunk -> compute_psd(chunk; window_seconds = 60.0);
                         chunk_seconds = 3600,    # 1-hour PSD per chunk
                         gap_handling  = :skip)

# results is a DataFrame with one row per chunk.
# Each row's PSDResult is in results.result; convert to dB for analysis.
```

Each chunk is a 1-hour `Audiodata`; `compute_psd` computes the PSD for that
chunk independently. The calibration is resolved per-chunk from the chunk's
metadata. This pattern is constant in memory regardless of archive size.

For a continuous time-averaged PSD over the full archive, average the per-chunk
results in linear power:

```julia
all_avg = [average_psd(r) for r in results.result]
grand_avg = mean(stack(all_avg, dims=2); dims=2) |> vec
db_grand_avg = to_dB(grand_avg)
```

---

## References

Merchant, N.D., et al. (2015). Measuring Acoustic Habitats. *Methods in
Ecology and Evolution*, 6, 257–265.
[DOI:10.1111/2041-210X.12330](https://doi.org/10.1111/2041-210X.12330)

Design decisions: DD-07 through DD-16 in `docs/design_decisions.md`.
