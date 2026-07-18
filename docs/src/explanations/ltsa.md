# Long-Term Spectral Average (LTSA)

A **long-term spectral average** (LTSA) shows how acoustic energy is distributed
across frequency over long timescales — hours, days, or a full deployment. It is
the standard first look at a passive-acoustic dataset: a single frequency × time
image in which diel cycles, vessel passages, weather events, and seasonal shifts
are all visible at a glance.

In EcoAcoustics.jl the LTSA is also the **reusable pre-aggregation product** for
downstream soundscape work. It is computed once per recording as a matrix of
linear power spectral density; arbitrary frequency bands and percentile
statistics are then derived from that matrix without touching the audio again.

---

## Position in the pipeline

```
spectrogram → complex STFT matrix           [SpectrogramResult]
      │
      ▼
compute_psd → power per Hz, per frame        [PSDResult, µPa²/Hz after calibration]
      │
      ▼
compute_ltsa → power per Hz, per time column [LTSAResult, µPa²/Hz]
      │
      ├──→ band time series   (integrate freq bands out of each column)
      │
      └──→ percentile statistics (down the time axis, per band)
```

`compute_ltsa` composes the layers below it — it adds no new spectral-estimation
or normalization math (DD-28). Each column is exactly
`average_psd(compute_psd(spectrogram(slice)))`.

---

## What an LTSA column is

The recording is divided into consecutive, **non-overlapping** time columns, each
`average_span_seconds` long (typically 60 s, 300 s, or 3600 s). Within one
column:

1. `spectrogram` computes the short-time Fourier transform using the inner FFT
   window (`fft_window_seconds`, default 1 s) and overlap (`fft_overlap`,
   default 0.5). A 60 s column with a 1 s / 50 % window contains ~119 FFT frames.
2. `compute_psd` converts each frame to power per Hz (µPa²/Hz once calibrated).
3. `average_psd` averages those frames **in linear power** — the energetic mean.

The averaged spectrum is one column of `LTSAResult.matrix`. There is no dB-mean
and no median mode: averaging in linear power is the physically correct operation
and matches the PSD and SPL layers (DD-16, DD-22, DD-28).

```
                 time columns  (average_span_seconds each)
              ┌────┬────┬────┬────┬────┐
   f  Nyquist │    │    │    │    │    │
   r          │    │    │    │    │    │   each cell = energetic mean of the
   e          │    │    │    │    │    │   FFT frames inside that column,
   q          │    │    │    │    │    │   at that frequency bin
   s        0 │    │    │    │    │    │
              └────┴────┴────┴────┴────┘
```

---

## Two timescales: the column span vs the FFT window

An LTSA has **two** independent time parameters, and they do different jobs:

| Parameter | Symbol | Controls | Typical value |
|---|---|---|---|
| `average_span_seconds` | column width | time resolution of the LTSA image | 60 s – 3600 s |
| `fft_window_seconds` | inner FFT window | frequency resolution of every column | 1 s |

- **`fft_window_seconds` sets frequency resolution.** Bin width is
  `df = fs / nfft ≈ 1 / fft_window_seconds`. A 1 s window gives 1 Hz bins; a
  0.1 s window gives 10 Hz bins. Longer windows resolve tonals (vessel lines,
  turbine tonals) more finely but smear short transients.
- **`average_span_seconds` sets time resolution.** It is how much the column
  averages over. Longer columns are smoother and cheaper to store/plot but hide
  short events; shorter columns show detail at the cost of a longer, noisier
  image.

The only hard coupling between them: a column must be at least as long as one FFT
window (`average_span_seconds × fs ≥ fft_window_seconds × fs`), or a column could
not contain a single frame. `compute_ltsa` asserts this.

---

## Non-overlapping columns and the dropped partial column (DD-28)

Column `j` covers samples `(j-1)·S + 1 … j·S`, where
`S = round(Int, average_span_seconds × fs)`. Columns never overlap, and the FFT
frames of one column never cross into the next.

The number of columns is `n_columns = div(nsamples, S)` — integer floor
division. If the recording length is not an exact multiple of the column span,
the final `< S` samples are **dropped**, not zero-padded or emitted as a short
column. This keeps every column a true `average_span_seconds` average, so
statistics taken across columns (means, percentiles) are not biased by a
short-and-therefore-different final column. See DD-28 for the full rationale.

`LTSAResult.column_times[j]` is the start time of column `j` in seconds, relative
to the start of the analysed audio, derived from the integer sample offset so it
reflects the actual rounding of `S`.

---

## Why the matrix is stored in linear units

`LTSAResult.matrix` holds **linear** PSD (µPa²/Hz when calibrated), not decibels.
This is deliberate: the LTSA is the input to two further operations that both
require linear power —

- **Band time series.** Integrating a frequency band (e.g. the 63 Hz
  one-third-octave, or a hand-picked 100–200 Hz band) is a sum of linear power
  across the band's bins, done per column. The candidate bands for the first
  soundscape analyses are listed in `CANDIDATE_ACOUSTIC_BANDS.org`; every one of
  them can be integrated out of the same LTSA matrix without re-reading audio.
- **Percentiles over time.** Exceedance levels (L10, L50, L90, …) and maxima are
  taken down the time axis of a band's linear series, then converted to dB.

Convert to dB only for display, with `to_dB(result)` — it returns a new matrix
and does not mutate the stored one.

---

## Calibration

Calibration is resolved **once**, before the column loop, using the same cascade
as `compute_psd` (pre-calibrated signal → attached calibration → recorder profile
→ warn and fall back to full-scale; DD-14). Resolving once means a
missing-calibration warning is emitted at most a single time rather than once per
column.

```julia
audio = read_audio("deployment_hour.wav"; recorder = "sm3m")
lt    = compute_ltsa(audio; average_span_seconds = 60.0)   # SM3M scalar cal auto-resolved
ltsa_units(lt)     # :µPa²_per_Hz
```

`LTSAResult.is_calibrated` follows the same rule as the PSD layer:
`audio.is_calibrated || !(resolved_cal isa NoCalibration)`. `ltsa_units`
dispatches on that flag (returning `:µPa²_per_Hz` or `:fullscale²_per_Hz`),
mirroring `psd_units` (DD-15). To override auto-resolution, pass `cal = ...`
explicitly.

---

## Worked example

```julia
using EcoAcoustics

audio = read_audio("deployment_hour.wav"; recorder = "sm3m")

# 1-minute LTSA columns, 1 s / 50 % inner FFT window (1 Hz bins).
lt = compute_ltsa(audio; average_span_seconds = 60.0,
                  fft_window_seconds = 1.0, fft_overlap = 0.5)

size(lt.matrix)        # (n_freqs, n_minutes)
lt.freqs               # frequency axis, Hz
lt.column_times        # [0.0, 60.0, 120.0, ...] seconds
ltsa_units(lt)         # :µPa²_per_Hz

lt_dB = to_dB(lt)      # Matrix{Float64}, dB re 1 µPa²/Hz — for plotting
```

### Reusing one FFT plan at archive scale

Every column uses the same inner FFT size, so a single pre-built FFTW plan can be
shared across all of them. For long recordings this avoids FFTW's per-call
plan-selection overhead:

```julia
plan = make_spectrogram_plan(audio.fs, 1.0)   # 1 s window at this fs
lt   = compute_ltsa(audio; average_span_seconds = 60.0,
                    fft_window_seconds = 1.0, fft_plan = plan)
```

The plan changes only *how* FFTW computes each transform, never the result — the
output is bit-identical to the no-plan call (verified in `test/test_ltsa.jl`).

---

## Parameter selection guidance

| Recording length | Suggested `average_span_seconds` | Notes |
|---|---|---|
| Minutes to ~1 hour | 1 – 10 s | Fine detail; individual vessel passes resolved |
| Hours to ~1 day | 60 s | The common default; diel structure clear |
| Days to a deployment | 300 – 3600 s | Smooth, storable; seasonal/weekly structure |

Choose `fft_window_seconds` from the *frequency* structure you need to resolve,
not the recording length: 1 s (1 Hz bins) is a good default for soundscape work;
drop to 0.1 s only when short transients matter more than tonal resolution.

---

## When not to use `compute_ltsa`

- **Recording shorter than one column span.** There is nothing to average across
  time — call [`compute_psd`](@ref) directly. `compute_ltsa` asserts rather than
  returning an empty result.
- **A single deployment-wide average spectrum is wanted** (one spectrum, not a
  time series) — use `average_psd(compute_psd(...))` over the whole signal.
- **Overlapping or dB-averaged columns are required.** Not supported by design
  (DD-28); both are non-breaking to add later if a concrete need appears.

---

## References

Merchant, N.D., et al. (2015). Measuring Acoustic Habitats. *Methods in
Ecology and Evolution*, 6, 257–265.
[DOI:10.1111/2041-210X.12330](https://doi.org/10.1111/2041-210X.12330)

Design decision: DD-28 in `docs/design_decisions.md`. Composes the PSD layer
(DD-08 through DD-16) and the spectrogram layer (DD-01 through DD-07).
