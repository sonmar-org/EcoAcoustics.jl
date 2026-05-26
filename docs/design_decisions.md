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

---

## Session: SPL layer (task 10)

### DD-18 — BandSPL struct stores per-frame series and six aggregate statistics

**Decided:** `BandSPL` holds `spl_dB::Vector{Float64}` (one value per PSD time
frame) plus ten pre-computed `Float64` aggregate scalars: `mean_dB`, `median_dB`,
`L1_dB`, `L5_dB`, `L10_dB`, `L25_dB`, `L75_dB`, `L90_dB`, `L95_dB`, `L99_dB`.
Aggregates are computed once at construction and stored on the struct; they are
not recomputed on access.

**Why:** The primary use pattern is: run `compute_spl` over an archive, then
query the statistics for visualization, reporting, or thresholding. If aggregates
were properties that recomputed on every access, archive-scale use would recompute
`quantile` millions of times. Storing them on the struct costs ten Float64 (80
bytes) per band per analysis window — negligible compared to the per-frame vector.

The ten statistics match the distribution shape reported in Merchant 2015 fig. 4
and the full percentile set specified in ADEON DPS Table C-1 (Ainslie et al. 2018):
L1/L99 (extreme values), L5/L95 (outer indicators), L10/L90 (inner indicators),
L25/L75 (quartiles), median (robust central tendency), and the energetic mean
(physically correct average). OSPAR and EU Marine Strategy Framework Directive
monitoring protocols typically report a subset of these.

**Where:** `struct BandSPL` in `src/soundscape/spl.jl`.

---

### DD-19 — Band validation rules: throws on Nyquist breach, warns on sub-10-Hz

**Decided:** `compute_spl` validates every entry in the caller-supplied `bands`
Dict before integration. Rules:

- `high_Hz > Nyquist` → `ArgumentError` (listing all offending labels)
- `low_Hz ≥ high_Hz` → `ArgumentError` (listing all offending labels)
- `low_Hz < 10 Hz` → `@warn` (listing all offending labels; computation proceeds)
- No PSD bins in band → `ArgumentError` (naming the band)

**Why:** The 10 Hz floor is the typical lower limit of calibrated hydrophone
response. Sub-10-Hz integration is not prohibited — infrasound and low-frequency
baleen whale work legitimately require it — but it warrants a warning because it
is almost always unintentional. Throwing on `high_Hz > Nyquist` is correct:
integrating above Nyquist is undefined and would silently include wrap-around
aliasing energy.

Collecting all offending labels before throwing (rather than stopping at the
first failure) is a usability choice: large band Dicts from band generators may
have multiple problems simultaneously, and iterative fix-run-fail cycles are
unnecessarily slow.

**Where:** Band validation block in `compute_spl(psd::PSDResult; ...)` in
`src/soundscape/spl.jl`. Tested in `test/test_spl.jl`.

---

### DD-20 — Band generators use ANSI preferred centers; decidecade is an alias for tol

**Decided:**
1. `octave_bands` and `tol_bands` use the tabulated ANSI preferred center
   frequencies (`OCTAVE_PREFERRED_HZ`, `TOL_PREFERRED_HZ`), not the values
   computed from the exact formula (`1000 × 2^(n−10)`, `1000 × 2^((n−30)/3)`).
2. `decidecade_bands` is a one-line alias for `tol_bands`. Its output is
   identical in every respect.

**Why (preferred centers):** The ANSI preferred values intentionally differ
from exact formula values at non-power-of-two centers (31.5 vs 31.25, 63 vs
62.5, 3.15 vs 3.155). PAMGuide, PAMGuard, MANTA, and Merchant 2015 all use the
preferred values. Using formula values would cause frequency-axis label mismatches
when comparing output against those tools, requiring a manual mapping step for
any cross-tool validation.

**Why (decidecade alias):** "Third-octave" (ANSI S1.11) and "decidecade"
(ISO 18405:2017) refer to the same 1/3-decade band scheme — the numerical
equivalence `log10(2^(1/3)) ≈ 0.1003` makes them functionally identical. The
term "decidecade" entered the underwater acoustics literature around 2018 and
is now preferred in that community. Providing both names under one
implementation avoids the maintenance burden of keeping two independent band
tables synchronized.

**Where:** `src/soundscape/bands.jl`. Documented in
`docs/src/explanations/spl_bands.md`.

---

### DD-21 — `compute_spl` asserts `psd_units(psd) === :µPa²_per_Hz` before integration

**Decided:** `compute_spl(psd::PSDResult; ...)` begins with:

```julia
@assert psd_units(psd) === :µPa²_per_Hz "compute_spl requires a calibrated PSD..."
```

An uncalibrated PSD results in `AssertionError`, not a silent wrong result.

**Why:** An uncalibrated PSD is in full-scale²/Hz. Integrating it and reporting
the result as "dB re 1 µPa" is physically meaningless — the number has no
acoustic interpretation. The failure is a programming error (forgot to calibrate),
not a user input error (bad argument at runtime). `@assert` communicates this
distinction: it is the correct Julia idiom for precondition violations rather
than runtime argument errors. The message points directly to the fix.

**Why not `ArgumentError`:** `ArgumentError` implies the caller passed a bad
value through the public API. An uncalibrated PSD is not a bad value; it is a
misconfigured pipeline. The semantic distinction matters for error handling:
production code that catches `ArgumentError` for band validation should not
accidentally swallow a precondition violation of this kind.

**Where:** First executable statement of `compute_spl(psd::PSDResult; ...)` in
`src/soundscape/spl.jl`. Tested in `test/test_spl.jl`
("compute_spl: uncalibrated PSD raises AssertionError").

---

### DD-22 — Energetic mean: average in linear power, then convert to dB

**Decided:** `BandSPL.mean_dB` is:

```julia
mean_dB = 10.0 * log10(mean(10.0 .^ (spl_dB ./ 10.0)))
```

not `mean(spl_dB)`.

**Why:** This is identical in motivation to DD-16 (which makes the same choice
for `average_psd`). For SPL: `mean(10·log10.(P)) ≠ 10·log10(mean(P))` unless
all frame powers are equal. The energetic mean is the physically correct average
— it is the constant level that would produce the same total acoustic energy as
the time-varying signal. For a 10-minute recording with 9.5 min of 90 dB ambient
noise and 30 s of a 120 dB ship passage, the arithmetic mean of dB gives ~90 dB;
the energetic mean gives ~102 dB. The energetic mean correctly captures the
ship's dominant contribution to the acoustic environment.

**Where:** `mean_dB` computation in the per-band integration loop in
`compute_spl(psd::PSDResult; ...)` in `src/soundscape/spl.jl`. The choice is
documented in the `BandSPL` docstring under "Energetic mean".

---

### DD-23 — `fft_plan` forwarded through `compute_spl(audio)` → `compute_psd` → `spectrogram`

**Decided:** All `Audiodata`-accepting SPL methods (`compute_spl`,
`compute_tol`, `compute_octave`, `compute_decidecade`, `compute_millidecade`)
accept `fft_plan = nothing` and forward it through `compute_psd` to
`spectrogram` without modification.

**Why:** Archive-scale SPL computation (e.g., hourly TOL bands over a year of
data) processes thousands of chunks at the same sample rate and window length.
Without plan reuse, `spectrogram` reconstructs an FFTW plan on every chunk via
wisdom lookup — this is fast per-call but measurable at scale. Building the plan
once with `make_spectrogram_plan` and forwarding it eliminates the per-chunk
overhead entirely.

The forwarding is a one-line add at each layer (`fft_plan = nothing` in the
kwarg list; `fft_plan` passed through to the next call). The alternative —
building the plan inside `compute_spl` from the inferred `nfft` — would require
`compute_spl` to know the FFT length before calling `compute_psd`, which
couples the two layers incorrectly.

**Plan size mismatch:** If the caller passes a plan built for a different FFT
length, DD-04's assertion in `spectrogram` fires with a clear error message.
This is tested via the plan-size-mismatch path in `test/test_spl.jl`.

**Where:** `fft_plan` kwarg in `compute_spl(audio::Audiodata; ...)` and all
four convenience wrappers in `src/soundscape/spl.jl`; forwarded into
`compute_psd(audio::Audiodata; fft_plan)` in `src/soundscape/psd.jl`.

---

### DD-24 — PSD cross-validation excludes bins below 30 Hz; symmetric vs periodic Hann window

**Decided:** `_compare_psd` in `test/test_psd_pamguide_validation.jl` accepts
a `freq_lo` keyword argument (default `0.0` Hz). All statistics — mean absolute
residual, percent exceeding threshold, and neighborhood anomaly sampling — are
restricted to frequency bins with `freq_hz >= freq_lo`. All PAMGuide PSD
testsets pass `freq_lo = 30.0`.

**Why — window convention mismatch at low frequencies:**

EcoAcoustics uses a **symmetric Hann window** (`_make_window` in
`src/audio/dsp_helpers.jl`):

```
w[n] = 0.5 * (1 − cos(2π·n / (N−1))),   n = 0, …, N−1
```

PAMGuide (`PG_DFT.m` line 70) uses a **periodic Hann window**:

```
w[n] = 0.5 − 0.5·cos(2π·n / N),   n = 1, …, N
```

The difference is in the denominator (N−1 vs N). For N=2000 (1-second window
at 2 kHz), the window energies differ by only 0.0022 dB — negligible for
mean statistics. However, the two windows have different spectral leakage
patterns: the periodic form is the DFT of a rectangular window convolved with
itself, while the symmetric form is not periodic over N samples. At very low
frequencies, where only 10–29 cycles fit inside a 1-second window, this leakage
difference makes per-frame PSD values essentially uncorrelated between the two
tools — differences of 5–12 dB per frame at a single 10 Hz bin are observed,
even though the time-averaged mean over all frames agrees within 0.04 dB.

The pct-exceeding-0.5 dB metric is sensitive to this: CallingPeriod (5,614
frames) had 1.45% of bins exceeding 0.5 dB before the `freq_lo` cut, driven
entirely by bins 10–28 Hz. VesselPassage (29,895 frames) showed the same
pattern at 10 Hz (~50% of its frames exceeded 0.5 dB there) but diluted the
overall percentage below the 1.0% threshold due to the larger frame count.

**Why 30 Hz:** A 1-second window at 2 kHz contains exactly 30 complete cycles
at 30 Hz. Below that, spectral leakage from neighboring bins contributes a
significant fraction of the measured bin energy, and the two window shapes
diverge in how much leakage they admit. Above 30 Hz, per-frame agreement
between EcoAcoustics and PAMGuide is within 0.5 dB across all three test
recordings. The time-averaged mean at 10–29 Hz is within tolerance and confirms
there is no systematic calibration error at those frequencies.

**What this is not:** This is not a calibration or normalization error. It is a
documentation of a known per-frame variance effect at very low frequencies.
The PAMGuide CSV calibration offset (~0.10 dB) that existed in the original
CSVs was corrected separately by regenerating the CSVs with the correct
calibration parameters (Mh=−203.0 dB, G=+33.2 dB, vADC=1.5 V).

**Do not change `_make_window` to periodic Hann** to match PAMGuide. The
symmetric (N−1) form is the standard DSP textbook definition and is what MANTA
uses. Switching would break the MANTA cross-validation test and change the
window energy slightly. The periodic form is a MATLAB convention, not a physics
requirement.

**Where:** `freq_lo` kwarg in `_compare_psd` in
`test/test_psd_pamguide_validation.jl`; documented in
`test/validation/pamguide/README.md`.

---

## Session: ADEON DPS alignment

### DD-25 — Partial-bin edge interpolation not implemented; hard bin boundaries used

**Decided:** `compute_spl` selects integration bins via
`searchsortedfirst(psd.freqs, f_lo)` and `searchsortedlast(psd.freqs, f_hi)`.
A bin is either fully included (its centre frequency falls within the band) or
fully excluded. No fractional power contribution is applied to bins that
straddle a band edge.

**The DPS requirement:** ADEON DPS §2.2.2 (Ainslie et al. 2018) explicitly
describes partial-bin interpolation for decidecade bands: a PSD bin straddling
a band edge contributes only the fraction of its 1-Hz width that falls within
the band. EA does not implement this.

**Why the error is bounded:** With a 1-second analysis window at any integer
sample rate, the frequency resolution is exactly 1 Hz and bin centres are at
integer Hz values (0, 1, 2, ..., fs/2). Decidecade band edges computed from
ANSI S1.11 preferred centres via `f_c × 2^(±1/6)` are irrational. No bin
centre coincides with a band edge. The missing partial-bin power equals the
single straddled bin's power times the fractional shortfall at each edge
(at most one bin per edge, at most one bin's worth of power). For broadband
integration (hundreds of bins) this is < 0.005 dB. For narrow decidecade bands
above 100 Hz (≥ 10 bins) the bias is < 0.05 dB — within the SPL validation
tolerance. At very low frequencies where decidecade bands span only 2–4 bins
(e.g. the 10 Hz band at 2 kHz), the fractional error can reach ~0.1 dB.

**Why not implemented in v1:**
1. PAMGuide cross-validation passes within ±0.05 dB without it, confirming
   the deviation is below the practical validation tolerance for real recordings.
2. The implementation requires sub-Hz interpolation logic at each band edge —
   non-trivial and has no other use in v1.
3. Partial-bin interpolation makes a difference only when comparing against
   software that implements it (not PAMGuide, which also uses hard boundaries).

**If sub-Hz-accurate decidecade SPL is required** (regulatory reporting against
a reference tool that implements interpolation), this decision should be
revisited.

**Where:** Integration loop in `compute_spl(psd::PSDResult; ...)` in
`src/soundscape/spl.jl`.

---

### DD-26 — System-weighted broadband SPL not implemented; frequency-flat path only

**Decided:** `compute_spl` implements the *frequency-flat* broadband SPL path
only: integrate the calibrated PSD over the band, `L_p = 10 × log₁₀(Σ P(f)Δf)`.
It does not implement the *system-weighted* path defined in ADEON DPS Figure 1
(Ainslie et al. 2018).

**What system-weighted broadband SPL is:** ADEON DPS Figure 1 shows two parallel
paths. The right (frequency-flat) path is what EA implements. The left
(system-weighted) path applies a single representative sensitivity at a fixed
reference frequency (typically 250 Hz) to the time-domain signal and takes RMS.
The system-weighted approach is common in instruments where a single scalar
sensitivity is the calibration product and frequency response is certified flat
within the band.

**Why frequency-flat is preferred:**
1. Strictly more general: uses the full TF calibration when available, correctly
   accounting for frequency-dependent sensitivity across the band.
2. Decomposable: broadband SPL is exactly the energetic sum of all per-band
   contributions. System-weighted broadband is not decomposable into decidecade
   bands from the same pipeline.
3. For instruments with flat frequency response (DMON2, SM3M), the two paths are
   numerically identical.
4. For instruments with significant frequency-dependent response (Rockhopper),
   the PSD-integration path with TFCalibration is the physically correct choice.
   A single-frequency sensitivity estimate would introduce a systematic bias that
   grows with the slope of the TF curve.

**Consequence for users:** EA's "broadband SPL" is always the frequency-flat
value. Comparisons against tools that report system-weighted broadband SPL for
instruments with non-flat TF curves will show a systematic offset. This is
documented in the `compute_spl` "Do not use when" section and in the glossary.

**Where:** `compute_spl` in `src/soundscape/spl.jl`; glossary entry in
`docs/src/explanations/glossary.md`.

---

### DD-27 — No `bands=nothing` default; `bands` is a required keyword argument

**Decided:** `compute_spl(psd; bands)` has no default for `bands`. Callers must
always supply an explicit `Dict{Symbol, Tuple{Float64,Float64}}`. There is no
function named "broadband SPL" and no hidden default that creates a broadband
band.

**Why:** A `bands=nothing` default that silently creates a `:broadband =>
(10.0, Nyquist)` band is a hidden assumption. The lower edge (10 Hz) and the
labelling (:broadband) are not derivable from the signal; they reflect a
measurement intent that the caller must state. Hiding them:

1. Makes the output key `:broadband` appear in results without any code-visible
   evidence of where it comes from, breaking traceability.
2. Tempts callers to treat "broadband" as a well-defined single number rather
   than a band-limited integral over an explicitly chosen range. Different
   deployments have different lower limits of hydrophone response; 10 Hz is
   a reasonable default for many instruments but is wrong for others.
3. Is inconsistent with `compute_tol`, `compute_octave`, and
   `compute_millidecade`, which all require an explicit frequency range.

Making `bands` required forces the caller to be explicit: `bands = Dict(:full =>
(10.0, 24000.0))` is unambiguous, auditable, and reproducible.

**What removed the default:** Task 7 in the ADEON DPS alignment session
(2026-05). The `bands=nothing` guard and `resolved_bands` local variable were
deleted; both `compute_spl` overloads now declare `bands` as a keyword with no
default.

**Consequence for existing callers:** Any call to `compute_spl(psd)` or
`compute_spl(audio)` without `bands=...` now fails at method dispatch with
`UndefKeywordError`. The fix is to add the explicit band dict. Convenience
wrappers (`compute_tol`, `compute_octave`, `compute_decidecade`,
`compute_millidecade`) are unaffected — they always pass an explicit `bands`
derived from the band generator.

**Where:** `compute_spl` signatures in `src/soundscape/spl.jl`. Validation tests
updated in `test/test_spl.jl` and `test/test_spl_pamguide_validation.jl`.
