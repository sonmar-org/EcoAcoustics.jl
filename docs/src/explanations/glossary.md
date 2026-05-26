# Glossary

Terms used in EcoAcoustics.jl documentation, drawn from Merchant et al. (2015),
NOAA Passive Acoustic Monitoring Data Standards documentation, and standard DSP
texts. Direct access to ISO 18405:2017 was not available during writing; any
definition sourced from secondary literature is noted as such. Definitions marked
*[uncertain — secondary source only]* should be verified against the ISO
standard before use in regulatory or publication contexts.

**Planned additions (v1.1 / v1.2).** When click and whistle detectors enter the
package, additional terms and citations will be added. Anticipated references
include: Frasier et al. for click detection (KERNO-F); Roch & Silbido /
Baumgartner et al. for whistle contour tracing (Silbido remora); Martin et al.
(2021) and the MANTA group for millidecade spectra; Wiggins & Hildebrand (2007)
for the Triton / HARP architecture underlying LTSA design. These are noted here
so that users comparing EcoAcoustics.jl output against those tools know the
methodological lineage.

---

## STFT — Short-Time Fourier Transform

The **Short-Time Fourier Transform** is computed by dividing a signal into
overlapping frames and applying the Discrete Fourier Transform (DFT) to each
frame independently. For frame `i` of a signal `x`, with window `w` of length
N and hop size `h`:

```
X_i[k] = Σ_{n=0}^{N−1}  w[n] · x[i·h + n] · exp(−j2πkn/N)
```

The result is a matrix of complex numbers: column `i` is the spectrum of frame
`i`, row `k` is the time history of frequency bin `k`. In EcoAcoustics.jl, the
STFT is computed by `spectrogram` and stored in `SpectrogramResult.stft`.

ISO 18405:2017 defines frequency-domain analysis of underwater sound in the
context of PSD (§3.1.2); the STFT is the standard computational path to that
quantity. *[Reference via Merchant 2015 — ISO text not directly consulted.]*

---

## Single-sided spectrum

A real-valued signal's DFT is *conjugate-symmetric*: the bin at negative
frequency −k carries exactly the same information as the bin at positive
frequency +k. The **single-sided spectrum** retains only the bins from 0 Hz
(DC) through fs/2 (Nyquist), discarding the redundant negative-frequency half.

In FFTW.jl, `rfft(x)` for a real vector of length N (even) returns a complex
vector of length N/2 + 1:

- Index 1: DC bin (0 Hz)
- Indices 2 to N/2: positive-frequency interior bins
- Index N/2 + 1: Nyquist bin (fs/2)

The single-sided amplitude correction — multiplying interior bins by √2 (or
interior power bins by 2) to restore the energy lost by dropping the negative
half — is applied at the PSD layer, not in `spectrogram`. This separation
keeps the STFT layer independent of normalisation conventions.

---

## Frame, window, hop, overlap

**Frame**: one short segment of the signal submitted to a single FFT. Length is
`window_length = round(Int, window_seconds × fs)` samples. Frames are labelled
1 to `num_frames`.

**Window** (or *analysis window*, *tapering window*): a smooth bell-shaped
function multiplied element-wise into each frame before the FFT. Its purpose is
to reduce *spectral leakage*: without tapering, the abrupt cut at each frame
boundary acts as a rectangular pulse and spreads energy across all frequency
bins. Common choices and their trade-offs:

| Window      | Sidelobe level | Main lobe width | Typical use                     |
|:----------- | --------------:| ---------------:|:-------------------------------- |
| Hann        |      −31 dB    |     Moderate    | General acoustic analysis        |
| Hamming     |      −43 dB    |     Moderate    | When slightly lower leakage needed |
| Blackman    |      −58 dB    |     Wide        | Narrow tonal signals             |
| Rectangular |       0 dB     |     Narrowest   | Calibration checks only          |

**Hop** (also *step size*, *frame advance*): the number of samples between the
start of one frame and the start of the next.
`hop = round(Int, (1 − overlap_fraction) × window_length)`. A hop equal to
half the window length (overlap = 50%) means consecutive frames share half their
samples.

**Overlap**: the fraction of the window shared between consecutive frames.
`overlap_fraction = 1 − hop / window_length`. PAMGuide default is 50% overlap.
Higher overlap gives finer time resolution (more frames) at proportionally
higher computational cost. `overlap_fraction = 1.0` is undefined (infinite
frames); `spectrogram` rejects it with an `AssertionError`.

---

## Window energy

```
window_energy = Σ wᵢ²   (equivalently: ‖w‖²)
```

The sum of squared window coefficients. Required for converting raw STFT
magnitudes to power spectral density: the PSD at frequency bin `k` is
proportional to `|X[k]|² / (fs × window_energy)`.

Merchant et al. (2015) Section 2.1 calls this the normalisation factor and uses
the symbol W² for the sum. EcoAcoustics.jl stores it on `SpectrogramResult` so
the PSD layer has it without recomputing the window.

The window is applied to each frame **unmodified** (not pre-normalised to unit
energy). This is the PAMGuide / Merchant convention. Some other software
normalises the window so that `Σwᵢ² = 1` before applying it; that convention
absorbs the correction into the window itself and produces different absolute
power values. Always check which convention a tool uses when comparing results.

---

## Power spectral density (PSD, *P*_ss)

The **power spectral density** gives a signal's power per unit frequency. For
calibrated underwater acoustic data the units are µPa²/Hz. For a single-sided
STFT bin `k` from a calibrated signal:

```
P_ss[k] = |X[k]|²  /  (fs × window_energy)       [µPa²/Hz]
```

The subscript *ss* denotes "single-sided" (all energy concentrated in the
positive-frequency range, with interior bins doubled to account for the dropped
negative-frequency mirror). After taking 10 · log₁₀, the result is in
dB re 1 µPa²/Hz.

In EcoAcoustics.jl, PSD is computed at the metric layer (task 9) from a
`SpectrogramResult`. Notation follows Merchant et al. (2015).

*ISO 18405:2017 §3.1.2 — definition via Merchant 2015; ISO text not directly
consulted. [Uncertain — secondary source only.]*

---

## Sound pressure level (SPL, *L*_p)

The **sound pressure level** is the RMS sound pressure expressed in decibels
relative to the underwater acoustic reference of 1 µPa:

```
L_p = 20 · log₁₀(p_rms / p_ref)     p_ref = 1 µPa  (underwater)
```

In EcoAcoustics.jl, SPL is obtained by integrating the calibrated PSD over an
explicitly specified frequency band, then converting to dB. Third-octave-level
(TOL), decidecade, and arbitrary-band SPL all follow the same approach. There is
no function named "broadband SPL" — callers always pass an explicit frequency
range to `compute_spl`; see DD-27.

**Reference level note.** The underwater reference (1 µPa) differs from the
airborne reference (20 µPa). A value in dB re 1 µPa is not comparable to a
value in dB re 20 µPa; always check the reference when comparing across
literature sources or measurement systems.

*ISO 18405:2017 §3.1.3 — definition via Merchant 2015; ISO text not directly
consulted. [Uncertain — secondary source only.]*

---

## Temporal Observation Window (TOW) and Temporal Analysis Window (TAW)

Terms introduced in the ADEON Data Processing Specification (Ainslie et al.
2018, §2.1.1) to distinguish two nested time scales in soundscape reporting.

**Temporal Analysis Window (TAW):** The short window over which a single
spectral estimate (PSD frame) is computed. In EcoAcoustics.jl this is the
`window_seconds` parameter passed to `compute_psd` or `compute_spl`. The DPS
standard TAW durations are 1 s, 60 s, 1 h, 1 d, 1 month, and 1 year. The
package default of `window_seconds = 1.0` matches the DPS primary TAW.

**Temporal Observation Window (TOW):** The longer period over which multiple
TAW frames are aggregated into a single reported metric (mean, percentiles,
CDF). In EcoAcoustics.jl, all PSD frames computed from a single `Audiodata`
object or `compute_spl` call form one TOW; the aggregate statistics in
`BandSPL` (mean, median, percentiles) summarise that TOW. For archive-scale
reporting, the TOW is determined by how `process_chunks` groups chunks and
collects results — the package does not currently name a TOW parameter
explicitly.

---

## Decidecade

A **decidecade** band spans one-tenth of a decade on a logarithmic frequency
axis: a frequency ratio of `10^(1/10) ≈ 1.259`. This is numerically
indistinguishable from a **third-octave** band (ratio `2^(1/3) ≈ 1.260`);
the two terms refer to the same practical band width.

"Decidecade" is the ISO 18405:2017 preferred term for underwater acoustics.
"Third-octave" (ANSI S1.11) remains in common use in older literature and
cross-domain tools. EcoAcoustics.jl provides both `tol_bands` (legacy name)
and `decidecade_bands` (ISO name) as identical functions; output keys use the
`:tol_N` prefix for compatibility with PAMGuide, MANTA, and Merchant 2015.

The ISO 18405 symbol for decidecade-band sound pressure level is L_p,ddec.
EcoAcoustics.jl computes this via `compute_decidecade` or `compute_tol`.

**Edge formula note:** EcoAcoustics.jl uses ANSI S1.11 base-2 edges
(`f_c × 2^(±1/6)`). The ADEON DPS specifies base-10 edges (`f_c × 10^(±1/20)`).
The relative difference is ~0.04%; at 1 Hz PSD resolution the band boundaries
are identical to within one bin. See DD-20.

---

## Energetic mean

The **energetic mean** of a set of SPL values `{L₁, L₂, …, Lₙ}` is:

```
L_energetic = 10 × log₁₀(1/n × Σ 10^(Lᵢ/10))
```

This averages in the linear power domain and then converts back to dB. It is
the physically correct average when the quantity of interest is *mean power*
over time.

The arithmetic mean of dB values (`(L₁ + L₂ + … + Lₙ) / n`) does not
preserve total power and is incorrect for acoustic averaging. Merchant 2015
fig. 4 reports the energetic mean; the ADEON DPS calls it "the arithmetic mean
computed on linear values" (meaning: compute the mean on the linear quantities
`10^(Lᵢ/10)`, then express the result in dB).

In EcoAcoustics.jl, `BandSPL.mean_dB` is the energetic mean of the per-frame
SPL time series within the observation window.

---

## System-weighted vs frequency-flat broadband SPL

The ADEON DPS Figure 1 shows two parallel paths to broadband sound pressure
level (L_p,broadband):

**Frequency-flat path** (what EcoAcoustics.jl implements): Compute the
calibrated PSD using the full frequency-dependent calibration, then integrate
over the band: `L_p = 10 × log₁₀(Σ_k P(fₖ) × Δf)`. The result uses the
correct calibration at every frequency bin.

**System-weighted path** (not implemented in v1): Apply a single
representative sensitivity value — typically the hydrophone sensitivity at
250 Hz — directly to the time-domain signal, then take RMS. This is faster
(no FFT required for the broadband value alone) but uses only one point on
the calibration curve. For instruments with flat frequency response the two
paths agree exactly. For instruments with frequency-dependent response
(e.g. Rockhopper with a TF calibration curve), the frequency-flat path is the
physically correct choice.

See DD-26 for the rationale for implementing only the frequency-flat path in v1.

**No "broadband SPL" function.** EcoAcoustics.jl does not expose a
`compute_broadband_spl` function or a default broadband band. Callers always
pass an explicit frequency range: `compute_spl(psd; bands = Dict(:full =>
(10.0, 24000.0)))`. This makes the lower edge — which depends on the instrument's
calibrated response — explicit and auditable rather than silently assumed. See
DD-27.

---

## ISO level notation

ISO 18405:2017 defines standardised symbols for underwater acoustic levels.
The symbols used in ADEON DPS documentation and in EcoAcoustics.jl docstrings:

| Symbol | Name | EcoAcoustics.jl function / field |
|:--|:--|:--|
| L_p | Sound pressure level (broadband) | `compute_spl`, `BandSPL.mean_dB` |
| L_pk | Peak sound pressure level | not implemented (v1) |
| L_p,ddec | Decidecade-band sound pressure level | `compute_decidecade` / `compute_tol` |
| L_E | Sound exposure level | not implemented (v1) |
| L_E,ddec | Decidecade-band sound exposure level | not implemented (v1) |
| β | Kurtosis (impulsiveness metric) | not implemented (v1) |

All levels are relative to 1 µPa (underwater). In-air equivalents use 20 µPa;
set `environment = :air` in EcoAcoustics.jl to switch reference pressures.

*Definitions from ISO 18405:2017. [ISO text not directly consulted; sourced
via ADEON DPS (Ainslie et al. 2018) and Merchant 2015.]*

---

## References for this glossary

- Merchant, N.D., Fristrup, K.M., Johnson, M.P., Tyack, P.L., Witt, M.J.,
  Blondel, P. & Parks, S.E. (2015). Measuring acoustic habitats. *Methods in
  Ecology and Evolution*, **6**, 257–265.
- NOAA Passive Acoustic Monitoring Data Standards documentation. Exact edition
  not recorded at time of writing.
- ISO 18405:2017 *Underwater acoustics — Terminology*. Definitions noted as
  *[uncertain — secondary source only]* were sourced from the above secondary
  literature, not from the ISO text directly.
- Ainslie MA, Miksis-Olds JL, Martin B, Heaney K, de Jong CAF, von
  Benda-Beckmann AM, Lyons AP (2018) ADEON Underwater Soundscape and Modeling
  Metadata Standard. Soundscape Specification deliverable v1.0.
