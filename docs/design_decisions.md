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

---

## Session: PSD layer (task 9)

### DD-07 — Even-nfft constraint enforced in `spectrogram`, not `compute_psd`

**Decided:** `spectrogram` asserts `iseven(nfft_actual)` immediately after
determining the FFT length. `compute_psd` does not re-check.

**Why:** The single-sided correction in `compute_psd` treats `psd_linear[end, :]`
as the Nyquist bin and does not double it. This is only correct when nfft is even:
for even nfft, `rfft` of a length-N signal produces N/2+1 bins — DC, N/2−1
interior bins, and one Nyquist bin. For odd nfft, the last rfft bin is NOT Nyquist
(there is no Nyquist for odd DFT length), so the correction would be silently
wrong. Placing the assertion in `spectrogram` catches the error at the source with
a message that explains the PSD dependency. Placing it in `compute_psd` would mean
the error is only caught when the user computes a PSD — too late if the
`SpectrogramResult` is stored or used for another purpose first.

**Where:** `spectrogram` function, after `nfft_actual` is determined
(`src/soundscape/spectrogram.jl`). Message cites DD-07 and gives the fix
(`nfft = window_length + 1` for odd window lengths). Tested in
`test/test_psd.jl` ("spectrogram: odd nfft rejected").

---

### DD-08 — Calibration applied at the PSD layer, not the spectrogram layer

**Decided:** `spectrogram` never applies calibration. `compute_psd` applies it
via `apply_calibration!(psd_linear, freqs, cal)` after normalization.

**Why:** The spectrogram produces complex amplitudes in full-scale units. It has
two downstream consumers: the PSD pipeline (needs physical calibration) and future
click/whistle detectors (need phase information but not physical units). Applying
calibration in the spectrogram would mix concerns and force click detection to work
with calibrated, magnitude-only data. Keeping the spectrogram layer calibration-free
preserves both options from the same `SpectrogramResult`.

**Where:** `compute_psd(spec, cal)` in `src/soundscape/psd.jl`; documented in
spectrogram and PSD docstrings and in `docs/src/explanations/psd.md`.

---

### DD-09 — Merchant 2015 Eq. (1) normalization: divide by `fs × window_energy`

**Decided:** The PSD normalization is:

```
psd[k, j] = correction_k × |STFT[k, j]|² / (fs × W)
```

where `W = Σwᵢ²` (`SpectrogramResult.window_energy`).

**Why:** This matches the PAMGuide / MANTA / Merchant 2015 convention. An
alternative convention normalizes the window so that `Σwᵢ² = 1` before applying it
to each frame, effectively absorbing the correction into the window. That convention
is used by some signal-processing textbooks but not by the marine bioacoustics
community. Using the Merchant form ensures that PSD values agree with MANTA to
within numerical tolerance (< 0.1 dB) without any post-hoc scaling.

**Where:** Step 3 of `compute_psd(spec, cal)` in `src/soundscape/psd.jl`.
`window_energy` is stored unnormalized in `SpectrogramResult` (DD-02) precisely so
it can be used here without re-computing it.

---

### DD-10 — Single-sided rfft correction: interior bins ×2, DC and Nyquist ×1

**Decided:** After computing `abs2.(stft)`, interior frequency bins
(`psd_linear[2:end-1, :]`) are multiplied by 2. DC (row 1) and Nyquist (last row)
are not modified.

**Why:** `rfft` drops the negative-frequency half of the DFT. For a real-valued
signal, every interior DFT bin `k` has a negative-frequency counterpart at bin
`N-k` with equal magnitude. The power in the positive-frequency bin therefore
represents only half of the total power at that frequency; doubling restores it.
DC (k=0) and Nyquist (k=N/2) are real-valued in the two-sided DFT — they are their
own complex conjugates and appear only once. Doubling them would incorrectly inflate
DC and Nyquist power by 3 dB. The correction is validated by three tests: the
bin-aligned sine test checks that an interior bin gets the expected power; the
Nyquist test checks that the Nyquist bin is NOT doubled; and the Parseval identity
test checks that the total power sums correctly across all bins.

**Where:** Step 2 of `compute_psd(spec, cal)` in `src/soundscape/psd.jl`. Only
valid for even nfft (DD-07).

---

### DD-11 — `psd_linear` freshly allocated; `spec.stft` never modified

**Decided:** `compute_psd` allocates a new matrix with `abs2.(spec.stft)`, then
modifies it in-place for the correction and normalization steps. The original
`spec.stft` is left unchanged.

**Why:** A `SpectrogramResult` may be consumed by multiple downstream computations:
PSD (task 9), future click detection, and potentially whistle detection. If
`compute_psd` modified `spec.stft` in-place, computing the PSD would destroy phase
information needed by click detectors. Allocating a new matrix costs one N_freqs ×
N_frames Float64 allocation; for typical parameters (512 bins, 100 frames) this is
~400 KB, negligible compared to the audio data.

**Where:** Step 1 of `compute_psd(spec, cal)` in `src/soundscape/psd.jl`.

---

### DD-12 — `load_tf_calcurves` normalizes to canonical form `dB re full-scale per µPa`

**Decided:** `load_tf_calcurves` reads the Cornell multi-column header CSV and
converts `AnalogSensitivity_dB_re_1VperRefPress` to the canonical TFCalibration
form by subtracting `20·log10(vmax_peak_V)`. The resulting `TFCalibration.format`
is `:rockhopper_calcurves_csv`.

**Why:** The Cornell calibration is in dB re 1 V/µPa at the ADC input. The package
canonical form is dB re ADC full-scale per µPa, which is what `apply_calibration!`
expects. The conversion is `tf_dB_canonical = AnalogSensitivity_dB − 20·log10(Vmax_peak_V)`.
For `Vmax_peak_V = 5.0 V`, the shift is `20·log10(5) ≈ 13.98 dB`. The `format` field
in `TFCalibration` (`:rockhopper_calcurves_csv`) distinguishes this multi-column header
format from other formats a user might supply (e.g., a bare 2-column Raven Expedition
CSV, which is already in canonical form and must NOT have the 13.98 dB shift applied
again). The shift amount and direction are recorded verbatim in
`TFCalibration.conversion_notes` for auditability.

**Where:** `load_tf_calcurves` in `src/recorders/rockhopper.jl`. The `format` and
`conversion_notes` fields are defined in `TFCalibration` in
`src/audio/calibration.jl`.

---

### DD-13 — `get_profile` uses `Val{recorder_id}` dispatch

**Decided:** `get_profile(id::Symbol)` creates a `Val{id}` value and dispatches to
`get_profile(::Val{:recorder_name})`. Adding a new recorder profile requires only
one new method in the recorder's source file.

**Why:** A dictionary-based registry (`PROFILE_REGISTRY[id]`) is the obvious
alternative but requires every recorder file to register itself at module load
time via a global mutation. `Val` dispatch achieves the same extensibility without
a shared mutable registry. The symbol-to-Val bridge in `recorders.jl` is the only
central code; each recorder file adds its own `get_profile(::Val{:name})` method.
The fallback method throws `ArgumentError` with the recorder name, making the
error message specific.

**Where:** `get_profile(id::Symbol)` bridge and `get_profile(::Val{T}) where T`
fallback in `src/recorders/recorders.jl`; `get_profile(::Val{:rockhopper})` method
in `src/recorders/rockhopper.jl`.

---

### DD-14 — `_psd_calibration` three-step cascade

**Decided:** When `compute_psd(audio::Audiodata)` is called, calibration is
resolved by `_psd_calibration(audio)` in priority order:
1. `audio.is_calibrated == true` → signal already in physical units; return
   `NoCalibration()` (no further correction at the PSD layer).
2. `audio.calibration isa !NoCalibration` → an explicit calibration was attached
   at I/O time (e.g. `ScalarCalibration` for SM3M); use it.
3. Try `get_profile(Symbol(recorder)).tf`. If the profile carries a `TFCalibration`
   in a field named `:tf` (checked via `hasproperty`), return it. Catch
   `ArgumentError` from unregistered recorders without re-throwing.
4. Nothing found → `@warn` and return `NoCalibration()`.

**Why:** This cascade makes Rockhopper recordings "just work" without the user
having to supply calibration explicitly: `compute_psd(audio; window_seconds=1.0)`
auto-resolves to the shipped TF calibration. Step 1 prevents double-calibration
when the user has already called `apply_calibration!` in the time domain. Step 2
handles SM3M and LS1X (scalar calibration attached by `read_audio`). Step 3 handles
recorders with typed profiles (currently only Rockhopper). `hasproperty` is used
rather than requiring `AbstractRecorderProfile` to define a `:tf` field, so adding
a recorder without a TF (e.g., one that only has a scalar sensitivity) doesn't
require interface changes.

**Where:** `_psd_calibration(audio::Audiodata)` in `src/soundscape/psd.jl`.

---

### DD-15 — `psd_units` dispatches on `is_calibrated::Bool`, not on `cal` type

**Decided:** `psd_units(p::PSDResult)` returns `:µPa²_per_Hz` when
`p.is_calibrated` is `true`, and `:fullscale²_per_Hz` when `false`. There are no
separate `psd_units(::NoCalibration)` etc. overloads.

**Why:** If `psd_units` dispatched on `p.cal`, a signal pre-calibrated in the time
domain (via `apply_calibration!`) and then processed with `compute_psd` (which
receives `cal = NoCalibration()` from `_psd_calibration` step 1) would incorrectly
report `:fullscale²_per_Hz` even though the PSD values are physically in µPa²/Hz.
Dispatching on `is_calibrated::Bool` is semantically correct for all cases where
`PSDResult.is_calibrated` accurately reflects the physical state.

**Implementation note:** The primitive sets `is_calibrated = !(cal isa NoCalibration)`,
which is correct for direct use of the primitive. The `Audiodata` wrapper corrects
the flag after the fact: `is_cal = audio.is_calibrated || result.is_calibrated`.
When both are false the wrapper returns the primitive's result directly (no
allocation); when `audio.is_calibrated` is true and `resolved_cal` is
`NoCalibration`, the wrapper reconstructs the `PSDResult` with `is_calibrated =
true`. This ensures `psd_units` returns `:µPa²_per_Hz` for pre-calibrated signals.

**Where:** `psd_units(p::PSDResult)` in `src/soundscape/psd.jl`. Tested in
`test/test_psd.jl` ("psd_units: dispatch on is_calibrated flag").

---

### DD-16 — `average_psd` averages in linear power, not in dB

**Decided:** `average_psd(result::PSDResult)` computes `mean(result.psd_linear;
dims=2)` in linear power (µPa²/Hz or full-scale²/Hz), then collapses the result to
a `Vector{Float64}`. Conversion to dB is done afterward by the caller via `to_dB`.

**Why:** Averaging in dB is NOT equivalent to averaging in linear power:
`mean(10·log10.(x)) ≠ 10·log10(mean(x))` except when all values are equal. For
bioacoustic data with high temporal variation (ship noise, biological choruses),
dB averaging underestimates the mean power — sometimes by several dB. The Merchant
2015 LTSA convention and all PAMGuide-family tools (MANTA, Triton) average in
linear power. EcoAcoustics.jl follows this convention to ensure numerical agreement
with reference tools.

**Where:** `average_psd(result::PSDResult)` in `src/soundscape/psd.jl`. Convention
explicitly stated in the function docstring and `docs/src/explanations/psd.md`.

---

## Session: Code review remediation (set 4)

### DD-17 — `strict` parameter convention across `read_audio` and `lookup_calibration`

**Decided:** Both functions follow the same rule: `strict=false` (default) emits
`@warn` on missing or unrecognised metadata and returns a sentinel value
(`NoCalibration()` for an unrecognised recorder, `DateTime(0)` for a missing
timestamp); `strict=true` throws an `ArgumentError`. This was the pre-existing
behaviour; DD-17 records it explicitly so future modifications to either function
maintain the consistency.

**Why:** A user who reads one function's docstring and calls the other should not
be surprised. Consistent semantics across `read_audio` and `lookup_calibration`
allow a single mental model: `strict=false` is exploration mode (forgiving, warns
to keep you informed); `strict=true` is production mode (fail fast). The asymmetry
that prompted the code review finding arose from an earlier version where
`lookup_calibration` used `error()` (throwing `ErrorException`) rather than
`throw(ArgumentError(...))`. The remediation session upgraded the exception type to
`ArgumentError` and sharpened both docstrings to state the rule explicitly.

**Where:** `lookup_calibration` in `src/recorders/recorders.jl` (docstring and
strict-mode throw); `read_audio` in `src/audio/read_audio.jl` (docstring `strict`
argument description). The `strict` kwarg is passed from `read_audio` into
`lookup_calibration` and `parse_filename` at the call sites inside `read_audio`.
