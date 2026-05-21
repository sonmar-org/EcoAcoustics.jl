# ─── Spectrogram (STFT) ───────────────────────────────────────────────────────
#
# Lowest-level time-frequency primitive. All downstream metrics (PSD, SPL,
# LTSA, soundscape indices) consume a SpectrogramResult. No calibration,
# magnitude correction, or single-sided amplitude scaling is applied here —
# those belong to the PSD layer.

"""
    SpectrogramResult

Result of a Short-Time Fourier Transform (STFT) computation. Stores the complex
single-sided spectrum for every analysis frame, the time and frequency axes, and
the window metadata required for downstream power spectral density (PSD)
normalisation.

This is the lowest-level time-frequency object in EcoAcoustics.jl. All
downstream metrics — PSD, SPL, LTSA, and soundscape indices — consume a
`SpectrogramResult`. It does not contain calibrated units; calibration is applied
at the metric layer.

Fields
------
- `stft::Matrix{ComplexF64}`:
    Complex single-sided STFT. Size: `(nfft ÷ 2 + 1) × num_frames`. Rows are
    frequency bins (DC at row 1, Nyquist at the last row); columns are time
    frames. No magnitude correction, single-sided amplitude scaling, or
    calibration has been applied.
- `time::Vector{Float64}`:
    Frame center times in seconds from the start of the input signal. Center of
    frame `i` (1-indexed) is computed as `((i-1) * hop + window_length / 2) / fs`.
    This matches the PAMGuide convention; the true midpoint of an even-length
    window would be `(window_length - 1) / 2`, but PAMGuide and Merchant 2015
    use `window_length / 2` for consistency. Users comparing against other tools
    (Triton, MANTA) should verify their frame-center conventions.
- `freqs::Vector{Float64}`:
    Bin-centre frequencies in Hz. Single-sided: 0 Hz (DC) through `fs/2`
    (Nyquist). Length equals `nfft ÷ 2 + 1`. Bin `k` (0-indexed) has frequency
    `k * fs / nfft`.
- `fs::Float32`:
    Sample rate in Hz. Stored for reproducibility and downstream use.
- `window_energy::Float64`:
    Sum of squared window coefficients: `Σ wᵢ²`. Required for PSD normalisation
    per Merchant et al. (2015). The window is applied to each frame without
    pre-scaling; `window_energy` captures the resulting energy reduction so the
    PSD layer can correct for it.
- `window::Symbol`:
    Window function used: `:hann`, `:hamming`, `:blackman`, or `:rectangular`.
    Stored for reproducibility.
- `nfft::Int`:
    FFT length used. May exceed the window length when zero-padding is requested
    (via the `nfft` kwarg to [`spectrogram`](@ref)). Stored for reproducibility.
- `hop::Int`:
    Hop size in samples between consecutive frame starts. Equal to
    `round(Int, (1 - overlap_fraction) * window_length)`. Stored for
    reproducibility and for downstream time-axis reconstruction.

References
----------
Merchant et al. (2015) Measuring Acoustic Habitats. Methods in Ecology and
Evolution, 6, 257–265. Window energy convention (Section 2.1) and frame-center
definition.
"""
struct SpectrogramResult
    stft::Matrix{ComplexF64}
    time::Vector{Float64}
    freqs::Vector{Float64}
    fs::Float32
    window_energy::Float64
    window::Symbol
    nfft::Int
    hop::Int
end

"""
    spectrogram(signal; fs, window_seconds, overlap_fraction=0.5,
                window=:hann, nfft=nothing, fft_plan=nothing)
                -> SpectrogramResult

Purpose:     Compute the complex Short-Time Fourier Transform (STFT) of a mono
             Float64 signal. Returns a [`SpectrogramResult`](@ref) containing the
             complex single-sided spectrum for every analysis frame, the time and
             frequency axes, and the window metadata required for downstream PSD
             computation. This is the lowest-level time-frequency primitive in
             EcoAcoustics.jl; PSD, SPL, LTSA, and soundscape indices all consume
             its output.

Arguments:
- `signal::AbstractVector{Float64}`: Input signal in arbitrary linear units
  (calibrated µPa or raw ADC — this function is unit-agnostic). Must be Float64.
- `fs::Real`: Sampling rate in Hz. Strictly positive.
- `window_seconds::Real`: Analysis window duration in seconds.
  `window_length = round(Int, window_seconds * fs)` samples. Must yield
  `window_length ≥ 2`.
- `overlap_fraction::Real = 0.5`: Fraction of the window that overlaps between
  consecutive frames. Must be in `[0, 1)`. `0.0` = no overlap (frames are
  adjacent, no sample is used twice); `0.5` = 50% overlap (PAMGuide default).
  `1.0` is excluded — it would produce infinite frames.
- `window::Symbol = :hann`: Tapering window applied to each frame before the FFT.
  One of `:hann`, `:hamming`, `:blackman`, `:rectangular`. Applied via
  `_make_window`.
- `nfft::Union{Int,Nothing} = nothing`: FFT length. Must be ≥ `window_length`.
  When `nothing` (default), equals `window_length` (no zero-padding). Set larger
  for finer frequency resolution without changing the analysis window duration —
  e.g. `nfft = 2 * window_length` halves the bin spacing.
- `fft_plan`: Pre-computed FFTW plan for arrays of length `nfft`. When `nothing`
  (default), a plan is built internally. For batch use via `process_chunks`,
  pre-build once with [`make_spectrogram_plan`](@ref) and pass here to avoid
  FFTW's plan-selection overhead on every chunk call.

Returns:     [`SpectrogramResult`](@ref). `stft` is complex, single-sided, size
             `(nfft ÷ 2 + 1) × num_frames`. No magnitude correction or
             calibration has been applied.

Constraints:
- `signal` must be `Float64`. All arithmetic in EcoAcoustics.jl uses Float64 for
  numerical precision; Float32 input is rejected with an informative error.
- Partial trailing frames are dropped: a frame is only computed if the full
  `window_length` samples are available. This is the Merchant 2015 convention.
  `num_frames = div(length(signal) - window_length, hop) + 1`.
- When `nfft > window_length`, the zero-pad region is appended after the windowed
  frame inside the STFT loop. Padding is transparent to the caller.
- Calibration is not applied here. For physically meaningful dB values, pass a
  calibrated signal (after `apply_calibration!`) and apply calibration at the
  PSD layer.

Fails when:
- `signal` is empty.
- `length(signal) < window_length`.
- `window_seconds * fs < 2` (window too short for a meaningful FFT).
- `nfft < window_length`.
- `overlap_fraction` is not in `[0, 1)`.
- `fft_plan` was built for a different length than `nfft` (checked upfront via
  the plan's internal `.sz` field).
- `window` is not one of the four recognised symbols (thrown by `_make_window`).

Example:
```julia
# 1-second signal at 48 kHz, 0.1 s window, 50% overlap, Hann window
signal = randn(Float64, 48000)
result = spectrogram(signal; fs = 48000.0, window_seconds = 0.1)

# result.stft size: (2401, 19) — 2401 frequency bins, 19 frames
# result.freqs[end]    == 24000.0   (Nyquist = fs/2)
# result.time[1]       ≈ 0.05       (center of first frame: 4800/2 / 48000)
# result.window_energy ≈ 1800.0     (Hann window energy ≈ 3N/8 for N=4800)
```

Do not use when:
- The signal has not been calibrated and physically meaningful SPL is required.
  Apply `apply_calibration!` first.
- The signal is multi-channel. Multi-channel audio is not supported in v1; pass
  one channel at a time.

References:
Merchant et al. (2015) Measuring Acoustic Habitats. Methods in Ecology and
Evolution, 6, 257–265. Window energy convention and frame-center definition.
"""
function spectrogram(signal::AbstractVector{Float64};
                     fs::Real,
                     window_seconds::Real,
                     overlap_fraction::Real   = 0.5,
                     window::Symbol           = :hann,
                     nfft::Union{Int,Nothing} = nothing,
                     fft_plan                 = nothing)::SpectrogramResult

    # ── Input validation ──────────────────────────────────────────────────────

    @assert length(signal) > 0 "spectrogram: signal must not be empty"

    # round(Int, x) rounds to nearest integer — appropriate for sample counts
    # derived from real-valued seconds × Hz products.
    window_length = round(Int, window_seconds * Float64(fs))
    @assert window_length >= 2 begin
        "spectrogram: window_length=$window_length " *
        "(window_seconds=$window_seconds × fs=$fs) must be ≥ 2"
    end

    # isnothing is the idiomatic Julia check for nothing (type-stable).
    nfft_actual = isnothing(nfft) ? window_length : nfft
    @assert nfft_actual >= window_length begin
        "spectrogram: nfft=$nfft_actual must be ≥ window_length=$window_length"
    end

    @assert 0 <= overlap_fraction < 1 begin
        "spectrogram: overlap_fraction=$overlap_fraction must be in [0, 1)"
    end

    # max(1, ...) guards against the degenerate case overlap_fraction ≈ 1.0 due
    # to floating-point rounding producing hop = 0 before the assertion above fires.
    hop = max(1, round(Int, (1 - overlap_fraction) * window_length))

    @assert length(signal) >= window_length begin
        "spectrogram: signal length $(length(signal)) is shorter than " *
        "window_length=$window_length"
    end

    # ── FFT plan ──────────────────────────────────────────────────────────────

    # plan_rfft caches FFTW's optimal algorithm for a specific input length.
    # Reusing a pre-built plan across many chunks avoids per-call plan-selection
    # overhead. If the caller supplies a plan, validate its input size before the
    # loop to give a clear error rather than a cryptic FFTW dimension mismatch.
    #
    # NOTE: .sz is an internal FFTW.jl field (NTuple of input dimensions) — not
    # part of the public API. Revisit if FFTW.jl has a major version bump.
    #
    # Julia idiom: if/else/end is an expression that returns a value. The result
    # of the chosen branch is assigned to fft_plan_actual.
    fft_plan_actual = if isnothing(fft_plan)
        FFTW.plan_rfft(zeros(Float64, nfft_actual))
    else
        @assert first(fft_plan.sz) == nfft_actual begin
            "spectrogram: fft_plan was built for input length " *
            "$(first(fft_plan.sz)), but nfft=$nfft_actual. " *
            "Use make_spectrogram_plan(fs, window_seconds) to build a matching plan."
        end
        fft_plan
    end

    # ── Window and energy ─────────────────────────────────────────────────────

    window_vec = _make_window(window, window_length)
    # sum(w -> w^2, v) computes Σvᵢ² without allocating an intermediate vector.
    window_energy = sum(w -> w^2, window_vec)

    # ── Frame count ───────────────────────────────────────────────────────────

    # Partial trailing frames (where fewer than window_length samples remain)
    # are dropped. div(a, b) is integer floor division.
    num_frames = div(length(signal) - window_length, hop) + 1
    @assert num_frames >= 1 begin
        "spectrogram: signal length $(length(signal)) produces 0 frames " *
        "(window_length=$window_length, hop=$hop)"
    end

    # ── Pre-allocation ────────────────────────────────────────────────────────

    # rfft of a length-nfft real signal returns nfft÷2+1 complex bins
    # (DC at index 1 through Nyquist at the last index). undef is safe here
    # because every cell is written in the loop below.
    stft = Matrix{ComplexF64}(undef, nfft_actual ÷ 2 + 1, num_frames)

    # frame_buffer holds one windowed, zero-padded frame. Initialised to zero
    # so the zero-padding region (indices window_length+1 : nfft_actual) stays
    # zero across all iterations — only indices 1:window_length are overwritten
    # per frame.
    frame_buffer = zeros(Float64, nfft_actual)

    # ── STFT loop ─────────────────────────────────────────────────────────────

    for i in 1:num_frames
        # start: 1-indexed position of the first sample of frame i.
        start = (i - 1) * hop + 1

        # @view creates a zero-copy view of the signal slice; without @view,
        # signal[start:end] allocates a copy. The .= broadcasts the element-wise
        # product directly into frame_buffer without allocation.
        frame_seg = @view signal[start : start + window_length - 1]
        frame_buffer[1:window_length] .= window_vec .* frame_seg

        # fft_plan_actual * frame_buffer computes the rfft out-of-place and
        # returns a new Vector{ComplexF64} of length nfft_actual÷2+1. Capital X
        # is the conventional notation for the frequency-domain form of x.
        # The assignment writes into the pre-allocated matrix column.
        stft[:, i] = fft_plan_actual * frame_buffer
    end

    # ── Time axis ─────────────────────────────────────────────────────────────

    # PAMGuide convention: center of frame i (1-indexed) is at sample index
    # (i-1)*hop + window_length/2 (0-indexed within the full signal).
    # In Julia, window_length/2 with Int arguments returns Float64 (/ is true
    # division), so the bracket expression is Float64 without an explicit cast.
    # Julia idiom: [expr for i in 1:n] is an array comprehension — it evaluates
    # expr for each i and collects the results into a Vector. The ./ then divides
    # every element of that vector by fs_f64.
    fs_f64 = Float64(fs)
    time   = [(i - 1) * hop + window_length / 2 for i in 1:num_frames] ./ fs_f64

    # ── Frequency axis ────────────────────────────────────────────────────────

    # range(0, fs/2; length=N) produces N evenly-spaced values from 0 to fs/2.
    # For N = nfft÷2+1, the step is (fs/2)/(nfft/2) = fs/nfft — exactly the
    # bin spacing of a length-nfft DFT. collect converts the lazy StepRangeLen
    # to a Vector{Float64} as required by SpectrogramResult.
    freqs = collect(range(0.0, fs_f64 / 2; length = nfft_actual ÷ 2 + 1))

    return SpectrogramResult(stft, time, freqs, Float32(fs),
                             window_energy, window, nfft_actual, hop)
end

# Float32 (or any other non-Float64 element type): reject with an informative
# error. Without this method, Julia would emit a MethodError naming
# AbstractVector{Float64} — correct but unhelpful to a new user.
function spectrogram(signal::AbstractVector; kwargs...)
    throw(ArgumentError(
        "spectrogram requires a Float64 signal; got $(eltype(signal)). " *
        "EcoAcoustics.jl uses Float64 throughout for numerical precision in " *
        "long FFTs. Convert with Float64.(signal) before calling."))
end

"""
    make_spectrogram_plan(fs, window_seconds, nfft=nothing) -> plan

Purpose:     Build an FFTW plan for the FFT size that `spectrogram` will use with
             the same `fs`, `window_seconds`, and `nfft` arguments. Pre-building
             the plan once and passing it to every `spectrogram` call via the
             `fft_plan` kwarg avoids FFTW's per-call plan-selection overhead,
             which matters when processing thousands of chunks.

Arguments:
- `fs::Real`: Sampling rate in Hz. Must match the `fs` passed to `spectrogram`.
- `window_seconds::Real`: Window duration in seconds. Must match the
  `window_seconds` passed to `spectrogram`.
- `nfft::Union{Int,Nothing} = nothing`: FFT length. When `nothing` (default),
  equals `window_length = round(Int, window_seconds * fs)`. Must match the `nfft`
  passed to `spectrogram`.

Returns:     An FFTW plan sized for arrays of length `nfft`. Pass directly as
             the `fft_plan` kwarg to `spectrogram`.

Constraints:
- `fs`, `window_seconds`, and `nfft` must be identical to the arguments used in
  the corresponding `spectrogram` call. A mismatch will be caught by the size
  assertion in `spectrogram`.
- The formula `window_length = round(Int, window_seconds * fs)` and
  `nfft_actual = isnothing(nfft) ? window_length : nfft` must mirror those in
  `spectrogram` exactly. If either formula changes in `spectrogram`, update
  this function to match.

Fails when:  No explicit validation — `spectrogram` will assert on a size
             mismatch when the plan is used.

Example:
```julia
# Pre-build a plan for 48 kHz data with a 0.1 s window, then reuse it across
# many chunks fed through process_chunks:
plan = make_spectrogram_plan(48000.0, 0.1)

result1 = spectrogram(chunk1.sig; fs = 48000.0, window_seconds = 0.1, fft_plan = plan)
result2 = spectrogram(chunk2.sig; fs = 48000.0, window_seconds = 0.1, fft_plan = plan)
```
"""
function make_spectrogram_plan(fs::Real,
                               window_seconds::Real,
                               nfft::Union{Int,Nothing} = nothing)
    # These two lines must mirror the corresponding computation in spectrogram().
    # If either formula changes there, update here to match.
    window_length = round(Int, window_seconds * Float64(fs))
    nfft_actual   = isnothing(nfft) ? window_length : nfft
    return FFTW.plan_rfft(zeros(Float64, nfft_actual))
end
