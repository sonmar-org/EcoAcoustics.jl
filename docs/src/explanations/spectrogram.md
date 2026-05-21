# The Spectrogram

A **spectrogram** — more precisely, a Short-Time Fourier Transform (STFT) — is
the result of slicing a long acoustic recording into short, overlapping segments
called *frames*, computing the frequency content of each one, and arranging the
results into a matrix. Each column of that matrix is one moment in time; each row
is one frequency bin. The value in each cell is a complex number whose magnitude
tells you how much energy was present at that frequency during that frame, and
whose phase encodes the timing of the oscillation within the window. Plotting the
squared magnitudes as a colour image produces the visual spectrograms you see in
Raven, PAMGuard, and Triton — EcoAcoustics.jl computes the same underlying
object.

---

## The pipeline

Every acoustic metric in EcoAcoustics.jl is built from the spectrogram upward.
The structure is:

```
spectrogram (task 8)  →  complex STFT matrix  [SpectrogramResult]
      │
      ▼
PSD   (task 9)        →  power per Hz bin  [µPa²/Hz after calibration]
      │
      ├──→ SPL / TOL  (task 10)  →  broadband or third-octave level  [dB re 1 µPa]
      │
      └──→ LTSA       (task 11)  →  long-term spectral average  [dB re 1 µPa²/Hz]
                │
                └──→ Soundscape indices  (task 12)  →  ACI, temporal entropy,
                                                        permutation entropy
```

`spectrogram` computes the raw complex STFT and nothing else. Every
transformation above it — squaring magnitudes, correcting for window energy,
converting to dB, averaging over time, integrating over frequency bands —
happens at the next layer. This separation means you can always inspect the
intermediate STFT directly, and it means the same `SpectrogramResult` can feed
the metric pipeline and future click/whistle detectors without modification.

---

## Why complex output

The `stft` field of `SpectrogramResult` is a `Matrix{ComplexF64}` — each entry
has a real part and an imaginary part. Most acoustic metrics only need the
*magnitude* of each entry (written `abs(X)`, or `abs2(X)` for power). The PSD
layer squares the magnitudes and applies the window-energy correction. So why
keep the complex form here?

**1. Click detection (v1.1 and later) requires phase.**
A click is an impulsive transient: its energy arrives across many frequencies
at nearly the same moment. *Phase coherence* across frequency bins — the
relationship between the imaginary and real parts — distinguishes a genuine
impulsive event from broadband continuous noise. Magnitude alone cannot capture
this. Storing the raw complex STFT preserves the option without any extra
computational cost at the spectrogram stage.

**2. Squaring is irreversible.**
Once you compute `abs2.(stft)`, the phase information is gone and cannot be
recovered without rerunning the spectrogram. By keeping the complex form, any
downstream algorithm can compute what it needs from a single `SpectrogramResult`.

If you only need power spectra, compute `abs2.(result.stft)` at the call site —
that is exactly what the PSD module does.

---

## Window energy convention

Before taking the FFT of each frame, the signal samples are multiplied
element-wise by a *tapering window* — a smooth bell-shaped curve that equals
zero at the edges and peaks in the middle. This step is necessary to prevent
**spectral leakage**: at a frame boundary, the signal is abruptly cut off, and
that sharp edge spreads energy across all frequency bins like a smear. Tapering
the samples to zero before cutting eliminates the smear.

The cost of tapering is reduced energy. A Hann window, for example, downweights
the samples near the edges so heavily that the windowed frame carries only about
3/8 of the raw frame's energy. If the PSD layer did not correct for this,
every spectrum would read low by that factor.

The correction factor is:

```
window_energy = Σ wᵢ²
```

This is the sum of squared window coefficients — the L² norm of the window
vector squared. Merchant et al. (2015) Section 2.1 names this the normalisation
factor and gives the PSD formula as:

```
P_ss[k] = |X[k]|² / (fs × window_energy)     [µPa²/Hz]
```

Storing `window_energy` on `SpectrogramResult` means the PSD layer has
everything it needs without recomputing the window or knowing which window was
used.

**Important convention note.** Some software pre-normalises the window so that
`Σwᵢ² = 1` before applying it; this absorbs the correction into the window
itself. EcoAcoustics.jl does *not* do this. The window is applied unmodified,
and `window_energy` is the unnormalised squared norm. This matches the Merchant
2015 and PAMGuide convention. If you compare outputs against other tools, check
which convention they use.

Approximate window energies for a window of length N:

| Window       | `window_energy`   | Notes                              |
|:------------ | -----------------:|:---------------------------------- |
| Rectangular  |           N       | No tapering; maximum leakage       |
| Hamming      | ≈ 0.397 × N       | Moderate leakage suppression       |
| Hann         | ≈ 0.375 × N       | Good general-purpose choice        |
| Blackman     | ≈ 0.305 × N       | Lowest leakage; widest main lobe   |

---

## Worked numerical example

Suppose you have a 1-second recording of a 1 kHz pure tone at fs = 48 000 Hz,
and you call `spectrogram` with a 0.1-second Hann window and 50% overlap.

**Step 1 — derive the parameters:**

```
fs             = 48 000 Hz
window_seconds = 0.1 s
window_length  = round(Int, 0.1 × 48 000) = 4 800 samples
hop            = round(Int, 0.5 × 4 800)  = 2 400 samples   (50% overlap)
nfft           = 4 800   (default: equals window_length, no zero-padding)

num_frames  = div(48 000 − 4 800, 2 400) + 1
            = div(43 200, 2 400) + 1
            = 18 + 1
            = 19 frames
```

**Step 2 — the output matrix:**

```
result.stft  →  size (2401, 19)
                  │     └─ 19 frames (columns)
                  └─ nfft÷2+1 = 2401 frequency bins (rows)
```

**Step 3 — the frequency axis:**

Frequency resolution = fs / nfft = 48 000 / 4 800 = **10 Hz per bin**.

The 1 kHz tone lands on bin 100 (0-indexed), because 100 × 10 Hz = 1 000 Hz
exactly. In Julia's 1-indexed arrays, this is `result.stft[101, :]`.

**Step 4 — the time axis:**

The PAMGuide convention places each frame's center at sample
`(i−1)×hop + window_length/2`:

```
time[1]  = (0 × 2400 + 2400) / 48 000 = 2 400 / 48 000 = 0.050 s
time[2]  = (1 × 2400 + 2400) / 48 000 = 4 800 / 48 000 = 0.100 s
  ⋮
time[19] = (18 × 2400 + 2400) / 48 000 = 45 600 / 48 000 = 0.950 s
```

Frames are spaced `hop / fs = 2 400 / 48 000 = 0.05 s` apart. The first
center is at 0.05 s, not 0 s: the window covers samples 1–4 800, and the
center of those samples is at sample 2 400 = 0.05 s.

**Step 5 — expected magnitude at 1 kHz:**

For a unit-amplitude sine at a bin-aligned frequency with a Hann window of
length N, the STFT magnitude at that bin is approximately N/4 = 4 800/4 =
**1 200**. This follows from the Hann window's DFT: the window's energy is
distributed mainly at DC and at ±1 bin, so a signal that falls exactly on a bin
produces a peak of approximately N/4 rather than the rectangular window's N/2.

**Verifying in Julia:**

```julia
using EcoAcoustics

fs             = 48_000.0
window_seconds = 0.1
f0             = 1_000.0    # exactly on bin 100

signal = sin.(2π * f0 .* (0 : 47_999) ./ fs)
result = spectrogram(signal; fs = fs, window_seconds = window_seconds)

# Bin 100 (0-indexed) → Julia index 101
println(abs(result.stft[101, 1]))   # ≈ 1200
println(result.time[1])             # ≈ 0.05
println(result.freqs[101])          # ≈ 1000.0
```

---

## References

- Merchant, N.D., Fristrup, K.M., Johnson, M.P., Tyack, P.L., Witt, M.J.,
  Blondel, P. & Parks, S.E. (2015). Measuring acoustic habitats. *Methods in
  Ecology and Evolution*, **6**, 257–265.
  Window energy convention (Section 2.1), PSD normalisation, and frame-center
  definition.

See also the [Glossary](glossary.md) for definitions of STFT, single-sided
spectrum, window energy, PSD, and SPL.

**Planned citations (v1.1 / v1.2).** When click and whistle detectors enter
the package, additional references will be added: Frasier et al. for click
detection (KERNO-F methodology used for MAB validation), Roch & Silbido /
Baumgartner et al. for whistle contour tracing (Silbido remora), Martin et al.
(2021) / MANTA for millidecade spectra, and Wiggins & Hildebrand (2007) for
the Triton / HARP architecture that informed the LTSA design.
