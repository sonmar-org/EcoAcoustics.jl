# Sound Pressure Level

**Sound pressure level** (SPL) is the primary soundscape metric produced by
EcoAcoustics.jl. It quantifies the acoustic power in one or more frequency
bands as a function of time, expressed in decibels relative to a reference
pressure. Every SPL value in this package is derived from an integrated
`PSDResult` — it is not computed directly from the waveform.

---

## Position in the pipeline

```
spectrogram → complex STFT matrix      [SpectrogramResult]
      │
      ▼
compute_psd → power per Hz bin          [PSDResult, µPa²/Hz]
      │
      ▼
compute_spl → band-integrated SPL       [SPLResult, dB re 1 µPa or 20 µPa]
   │   │
   │   ├── compute_tol / compute_octave / compute_millidecade
   │   │       (populate the bands Dict automatically)
   │   │
   │   └── compute_spl(audio::Audiodata; ...)
   │           (convenience wrapper; runs compute_psd internally)
   │
   └──→ LTSA, soundscape indices    (planned)
```

The PSD layer is mandatory. `compute_spl` does not re-run the FFT — it
integrates a `PSDResult` that already exists. This keeps the two operations
independently testable and avoids recomputing the spectrogram when multiple
band configurations are needed from the same recording.

---

## Integration formula

For each frequency band `[f_lo, f_hi]` and each time frame `j`, the band
power is:

```
band_power[j] = Σ_k  psd[k, j] × df
```

where:
- The sum runs over all PSD bins `k` with `f_lo ≤ freqs[k] ≤ f_hi`.
- `df = freqs[2] − freqs[1]` is the uniform bin width in Hz.
- `psd[k, j]` is in µPa²/Hz (calibrated).

The result is in **µPa²** — acoustic power integrated over the band.

This is a **Riemann sum**, not a trapezoid approximation. The choice is
deliberate: the PAMGuide convention treats each bin as a rectangle of width
`df` centered at `freqs[k]`. Using the trapezoidal rule would give a
different result at the band edges and break numerical agreement with PAMGuide
and MANTA. Do not change this without testing against PAMGuide reference
output.

### Bin selection: `searchsortedfirst` / `searchsortedlast`

Bins are selected with Julia's `searchsortedfirst(freqs, f_lo)` and
`searchsortedlast(freqs, f_hi)`. These return indices such that
`f_lo ≤ freqs[i_lo]` and `freqs[i_hi] ≤ f_hi` — inclusive at both edges.
A bin at exactly `f_lo` or exactly `f_hi` is included.

---

## dB conversion

```
SPL[j] = 10 × log10(band_power[j] / pref²)
```

| Environment  | `pref`    | `pref²`    | Output units     |
|:-------------|:----------|:-----------|:-----------------|
| `:water`     | 1 µPa     | 1 µPa²     | `:dB_re_1µPa`    |
| `:air`       | 20 µPa    | 400 µPa²   | `:dB_re_20µPa`   |

The reference pressure is squared because SPL is a power-domain quantity
(`10 × log10`, not `20 × log10`). A signal with `band_power = 1 µPa²` in
water gives `SPL = 10 × log10(1/1) = 0 dB re 1 µPa`. In air the same power
gives `10 × log10(1/400) = −26.02 dB re 20 µPa`.

---

## Aggregate statistics in BandSPL

Each `BandSPL` stores the per-frame `spl_dB` time series plus six aggregate
statistics computed from it.

### Energetic mean (DD-22)

```
mean_dB = 10 × log10(mean(10 .^ (spl_dB ./ 10)))
```

This converts to linear power, averages, then converts back. It is **not** the
arithmetic mean of `spl_dB`. The distinction matters whenever the SPL time
series has temporal variation — which it always does for real bioacoustic data:

```
mean(10·log10.(P)) ≠ 10·log10(mean(P))
```

The energetic mean is the value `L` such that a constant level of `L` dB
would produce the same total acoustic energy as the actual time-varying
signal. This matches the Merchant 2015 convention and MANTA. Arithmetic
averaging in dB underestimates the energetic mean when the signal contains
high transients (e.g., impulsive ship noise).

See DD-16 for the same argument applied to `average_psd`.

### Percentiles — ISO 18405 exceedance levels (DD-31)

`L_n` is the **exceedance level**: the SPL **exceeded n% of the frames**,
following the ISO 18405 terminology standard (also ADEON, OSPAR/JOMOPANS,
EU-MSFD). Equivalently, `L_n` is the **(100 − n)th statistical percentile**:

```
Ln = quantile(spl_dB, 1 - n/100)     # L5 = quantile(0.95), L95 = quantile(0.05)
```

So:

- **`L5` is the loud tail** (exceeded only 5% of the time);
- **`L95` is the quiet background** (exceeded 95% of the time — the standard
  ambient-noise indicator);
- `L1 ≥ L50 (median) ≥ L99`.

**This is the inverse of raw-percentile labelling.** On an SPD plot, a line
labelled `5%` is the *quiet* 5th percentile — that corresponds to EA's `L95`.
Keep table columns and SPD-plot labels on the same convention; do not mix `L_n`
(exceedance) with `%` (percentile).

Suggested methods wording: *"Percentile statistics are reported as exceedance
levels following ISO 18405: `L_n` is the band SPL exceeded n% of the averaging
period, so `L5` is the loud tail and `L95` the quiet background."*

### Single-frame case

When the PSD has only one frame, all aggregate statistics equal that frame's
value. `quantile` returns the single value for any percentile; `mean` of one
value is that value. This is correct behavior — there is no temporal
distribution to summarize. Do not add length > 1 guards around these
computations.

---

## Default band: broadband (DD-19)

Calling `compute_spl(psd)` with no `bands` argument uses a single broadband
band:

```julia
Dict{Symbol, Tuple{Float64, Float64}}(:broadband => (10.0, nyquist))
```

The 10 Hz lower floor matches the low-frequency response limit of most
calibrated hydrophones. The upper limit is the Nyquist frequency of the
recording (`psd.fs / 2`). Bands below 10 Hz are not rejected — they trigger
a `@warn` listing the offending band labels. The reasoning: some deployments
(infrasound, low-frequency baleen whale work) legitimately need sub-10-Hz
integration, and silently excluding those bands would be worse than warning.

---

## Band validation (DD-19)

Before any integration, `compute_spl` validates all bands in one pass. All
offending labels are collected before throwing, so the error message lists
every problem at once rather than forcing the user to fix them one at a time:

| Condition                  | Response                  |
|:---------------------------|:--------------------------|
| `low_Hz ≥ high_Hz`         | `ArgumentError` (all labels) |
| `high_Hz > Nyquist`        | `ArgumentError` (all labels) |
| `low_Hz < 10 Hz`           | `@warn` (all labels)      |
| No PSD bins in range       | `ArgumentError` (one label) |

The "no bins in range" case arises when the frequency resolution `df` is
coarser than the band width — for example, millidecade bands below ~435 Hz
are narrower than a 1 Hz bin spacing. The error message includes `df` so
the user knows the fix (use a longer window for finer frequency resolution).

---

## Calibration requirement (DD-21)

```julia
@assert psd_units(psd) === :µPa²_per_Hz "..."
```

`compute_spl` hard-asserts that the PSD is calibrated before touching any
values. An uncalibrated PSD is in full-scale²/Hz — integrating it produces a
number with no acoustic meaning, and silently returning it as "SPL" would be
incorrect. The assertion fires with a message pointing to the fix: use
`compute_psd` with a calibration, or use `compute_spl(audio::Audiodata; ...)`
which resolves calibration automatically via the DD-14 cascade.

Using `@assert` rather than `throw(ArgumentError(...))` is intentional: this
is a programming error, not a user input error. `@assert` can be disabled with
`--check-bounds=no` in performance-critical contexts (not recommended here, but
available).

---

## Archive-scale use: `fft_plan`

All `Audiodata`-accepting methods (`compute_spl`, `compute_tol`,
`compute_octave`, `compute_decidecade`, `compute_millidecade`) accept an
`fft_plan` keyword argument. It is forwarded without modification through the
call chain:

```
compute_spl(audio; fft_plan) → compute_psd(audio; fft_plan) → spectrogram(sig; fft_plan)
```

Build the plan once with `make_spectrogram_plan`, then reuse it across
`process_chunks` callbacks to avoid repeated FFTW wisdom lookups. For a
48 kHz recording with 1-second windows:

```julia
plan    = make_spectrogram_plan(48000.0, 1.0)    # built once
results = process_chunks(src, chunk ->
    compute_tol(chunk; window_seconds = 1.0, fft_plan = plan);
    chunk_seconds = 3600)
```

The plan's FFT length must match the resolved `nfft`. A mismatch throws
`AssertionError` (DD-04). Passing `fft_plan = nothing` (the default) causes
`spectrogram` to build its own plan per call — correct but slower for
archive-scale work.

`fft_plan` is not available on the `PSDResult` dispatch of any function —
the FFT has already been computed at that point, and the plan serves no
purpose.

---

## References

Merchant, N.D., et al. (2015). Measuring Acoustic Habitats. *Methods in
Ecology and Evolution*, 6, 257–265.
[DOI:10.1111/2041-210X.12330](https://doi.org/10.1111/2041-210X.12330)

Design decisions: DD-18 through DD-23 in `docs/design_decisions.md`.
