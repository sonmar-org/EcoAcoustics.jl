# Design Decisions Log

Non-obvious choices made during implementation — what was decided, why, and
where the code lives. Read this before changing anything in calibration,
sources, or the index schema.

Decisions are listed in chronological session order. Add a new entry at the
bottom whenever a non-obvious choice is made; do not edit past entries unless
correcting a factual error.

---

## Prior sessions (pre-spectrogram)

Decisions from the foundation work sessions (Audiodata, calibration, sources,
index, chunks) were made before this log existed. The key commitments from those
sessions are captured in CLAUDE.md under "Architecture" and "Key Constraints".
If a decision from a prior session needs to be recorded here, add it when the
relevant code is next touched.

---

## Session: Spectrogram primitive (task 8)

### DD-01 — Complex STFT output, not power

**Decided:** `SpectrogramResult.stft` is `Matrix{ComplexF64}` — the raw complex
single-sided STFT. Magnitude squaring and single-sided amplitude correction are
not applied here.

**Why:** Click detection (v1.1) requires phase information. Phase coherence
across frequency bins distinguishes impulsive transients (clicks) from
continuous broadband noise. Squaring the magnitudes is irreversible; once done,
phase is unrecoverable without rerunning the spectrogram. Storing the complex
form costs nothing and preserves the option for all downstream consumers.

**Where:** `SpectrogramResult.stft` field type (`src/soundscape/spectrogram.jl`).
PSD computation (task 9) will call `abs2.(stft)` and apply the single-sided
correction at that layer.

---

### DD-02 — Window energy stored unnormalized on result

**Decided:** `window_energy = sum(w -> w^2, window_vec)` stores the raw L²-norm
squared of the window vector. The window is applied to each frame without
pre-normalisation (i.e., `Σwᵢ²` is not forced to 1.0).

**Why:** This matches the Merchant et al. (2015) and PAMGuide convention.
Pre-normalising the window (so `Σwᵢ² = 1`) is an alternative used by some
software — it absorbs the correction into the window itself and changes absolute
power values. Storing the unnormalised energy keeps the spectrogram layer
convention-agnostic; the PSD layer divides by `fs × window_energy` to get
µPa²/Hz per Merchant Eq. (1).

**Where:** `SpectrogramResult.window_energy` field; `spectrogram` function,
window-energy block (`src/soundscape/spectrogram.jl`). Documented in
`SpectrogramResult` docstring and `docs/src/explanations/spectrogram.md`.

---

### DD-03 — Frame-center time convention: window_length/2, not (window_length-1)/2

**Decided:** Frame center time for frame `i` (1-indexed) is
`((i-1)*hop + window_length/2) / fs`. For even `window_length`, this places the
center at sample index `window_length/2` (0-indexed), which is one sample past
the true midpoint at `(window_length-1)/2`.

**Why:** Matches the PAMGuide / Merchant 2015 convention. The true midpoint of
a symmetric window of even length N is at index (N-1)/2 (0-indexed), but
PAMGuide uses N/2 for consistency across even and odd N. Tools built against
PAMGuide (Triton, MANTA) use the same definition; deviating would break
cross-tool time-axis alignment.

**Where:** Time vector computation in `spectrogram` function; documented
explicitly in both the `time` field of the `SpectrogramResult` docstring and
the "Worked numerical example" section of `docs/src/explanations/spectrogram.md`.
Users comparing against non-PAMGuide tools are warned in both places.

---

### DD-04 — FFT plan size validation via FFTW.jl internal `.sz` field

**Decided:** When `fft_plan` is passed to `spectrogram`, its input size is
validated with `@assert first(fft_plan.sz) == nfft_actual`. A clear
`AssertionError` is thrown upfront rather than letting FFTW throw a
`DimensionMismatch` inside the loop.

**Why:** The alternative (letting FFTW throw) produces an error message that
points inside `spectrogram`'s inner loop rather than at the call-site argument.
For advanced users who pre-build plans and pass them into `process_chunks`
callbacks, the upfront message is materially more debuggable.

**Trade-off:** `.sz` is an internal FFTW.jl field (`NTuple{N,Int}` holding
input dimensions), not part of the public API. It has been stable across
FFTW.jl 1.x but could change in a major version bump. A code comment marks this
in the source.

**Where:** Plan validation block in `spectrogram`
(`src/soundscape/spectrogram.jl`). If FFTW.jl removes `.sz`, fall back to
wrapping the plan application in a try/catch for `DimensionMismatch`.

---

### DD-05 — Parseval test uses corrected single-sided formula

**Decided:** The Parseval validation test in `test_spectrogram.jl` uses:

```
abs2(X[1]) + 2·Σ abs2(X[2:end-1]) + abs2(X[end]) = N · Σ abs2(x_windowed)
```

rather than the simpler `sum(abs2, X) / N = sum(abs2, x_windowed)`.

**Why:** The simpler formula is wrong for a single-sided rfft spectrum. The
rfft output drops the negative-frequency half of the DFT. The factor of 2 on
interior bins restores the energy from that dropped half. Without the factor of
2, the left-hand side equals approximately half the right-hand side for a random
signal — the test would pass only if single-sided amplitude correction (×2) had
already been applied, which it has not at the spectrogram layer. The corrected
formula was identified during the pre-implementation design review; the original
brief had the wrong formula.

**Where:** Test 3 ("Parseval's theorem") in `test/test_spectrogram.jl`. The
comment block in that test derives the formula from first principles and
explains why the naïve `sum(abs2)/N` form fails.

---

### DD-06 — Float32 input rejected via fallback method, not type parameter alone

**Decided:** A second method `spectrogram(signal::AbstractVector; kwargs...)`
throws `ArgumentError` with an explanatory message when a non-Float64 signal is
passed.

**Why:** The typed primary method `spectrogram(signal::AbstractVector{Float64};
...)` causes Julia to emit `MethodError: no method matching
spectrogram(::Vector{Float32}, ...)` for Float32 input — correct but opaque to
a user who doesn't know why Float64 is required. The fallback method catches all
non-Float64 element types and produces a message that names the type received,
states the reason (numerical precision in long FFTs), and gives the fix
(`Float64.(signal)`).

**Where:** Fallback method immediately after the main `spectrogram` method
(`src/soundscape/spectrogram.jl`). Tested in `test/test_spectrogram.jl` test 8,
which also checks that the message string is informative.
