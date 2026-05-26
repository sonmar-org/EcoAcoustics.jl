# ─── Sound Pressure Level ─────────────────────────────────────────────────────
#
# Band-integrated SPL computed from a calibrated PSDResult. Integration math,
# calibration assertion, and band-validation rules are in
# docs/design_decisions.md DD-18 through DD-23.

"""
    BandSPL

SPL statistics for one frequency band in an [`SPLResult`](@ref). Produced
by [`compute_spl`](@ref); not intended for direct construction.

Fields
------
- `band::Tuple{Float64, Float64}`:
    Integration band as `(low_Hz, high_Hz)` exactly as supplied by the caller.
    Edges are in Hz. The PSD bins included are those with centre frequency
    `low_Hz ≤ freq ≤ high_Hz`.
- `spl_dB::Vector{Float64}`:
    Per-frame SPL time series, one value per PSD time frame. Units are
    given by the enclosing `SPLResult.units` field (`:dB_re_1µPa` or
    `:dB_re_20µPa`).
- `mean_dB::Float64`:
    Energetic mean: `10 × log10(mean(10 .^ (spl_dB ./ 10)))`. Computes the
    mean in linear power, then converts to dB. This is NOT the arithmetic
    mean of `spl_dB` — see Merchant 2015 and DD-16 for why linear-domain
    averaging is required for physically correct results.
- `median_dB::Float64`:
    50th percentile of `spl_dB`.
- `L1_dB::Float64`:
    1st percentile of `spl_dB` — the level below which 1% of frames fall.
    Estimates the acoustic noise floor.
- `L5_dB::Float64`:
    5th percentile of `spl_dB` — the level below which 5% of frames fall.
    Low-ambient indicator.
- `L10_dB::Float64`:
    10th percentile of `spl_dB`.
- `L25_dB::Float64`:
    25th percentile of `spl_dB` — first quartile.
- `L75_dB::Float64`:
    75th percentile of `spl_dB` — third quartile.
- `L90_dB::Float64`:
    90th percentile of `spl_dB` — high-activity indicator.
- `L95_dB::Float64`:
    95th percentile of `spl_dB` — the level below which 95% of frames fall.
    High-transient indicator.
- `L99_dB::Float64`:
    99th percentile of `spl_dB` — the level below which 99% of frames fall.
    Extreme-transient proxy.

Percentile convention
---------------------
`L_n` is the **n-th percentile** of the SPL time series — the level *below*
which n% of frames fall. This is the standard statistical convention and
matches Merchant et al. (2015) fig. 4 and modern soundscape literature.

Note: some engineering standards use the inverse convention (L_n = level
*exceeded* n% of the time). EcoAcoustics.jl uses the statistical convention
throughout. When comparing output against other tools, verify which convention
they use — a reported L1 in one tool may equal L99 in another.

Single-frame note
-----------------
When the PSD has only one frame, all aggregate statistics (`mean_dB`,
`median_dB`, `L1_dB`, `L5_dB`, `L95_dB`, `L99_dB`) equal that single frame's
SPL value. There is no temporal distribution to summarize. This is correct
behavior, not a degenerate case to guard against.

References
----------
Merchant et al. (2015) Measuring Acoustic Habitats. Methods in Ecology and
Evolution, 6, 257–265. Percentile convention (fig. 4) and energetic mean.
"""
struct BandSPL
    band::Tuple{Float64, Float64}
    spl_dB::Vector{Float64}
    mean_dB::Float64
    median_dB::Float64
    L1_dB::Float64
    L5_dB::Float64
    L10_dB::Float64
    L25_dB::Float64
    L75_dB::Float64
    L90_dB::Float64
    L95_dB::Float64
    L99_dB::Float64
end

"""
    SPLResult

Result of a band-integrated SPL computation over one or more frequency bands.
Contains a per-frame SPL time series and aggregate statistics for each band.

Fields
------
- `bands::Dict{Symbol, BandSPL}`:
    One [`BandSPL`](@ref) entry per frequency band. Keys are the labels from
    the band specification Dict supplied to [`compute_spl`](@ref); values
    hold the per-frame time series and aggregate statistics for that band.
- `time::Vector{Float64}`:
    Frame-centre times in seconds from the signal start. Inherited from the
    [`PSDResult`](@ref) used in the computation (PAMGuide convention, DD-03).
- `fs::Float32`:
    Sample rate in Hz. Per the package convention (CLAUDE.md), `Float32`
    matches hardware precision.
- `units::Symbol`:
    Physical units of all SPL values in `bands`. One of:
    - `:dB_re_1µPa`  — underwater reference (1 µPa, set by `environment=:water`)
    - `:dB_re_20µPa` — in-air reference (20 µPa, set by `environment=:air`)
- `environment::Symbol`:
    `:water` or `:air`. Determines which reference pressure was used.

Constraints
-----------
- `units` must be `:dB_re_1µPa` or `:dB_re_20µPa`. Any other value throws
  `ArgumentError` at construction.
- Every `BandSPL` in `bands` must have `length(b.spl_dB) == length(time)`.
  If any band has a mismatched frame count, `ArgumentError` is thrown naming
  the offending band label.

Fails when:
- `units ∉ (:dB_re_1µPa, :dB_re_20µPa)` → `ArgumentError`
- Any `bands[label].spl_dB` length ≠ `length(time)` → `ArgumentError`

Example:
```julia
result = compute_spl(psd; bands = Dict(:full => (10.0, 24000.0)))
result.bands[:full].mean_dB    # energetic mean, dB re 1 µPa
result.bands[:full].spl_dB     # per-frame time series
result.time                         # frame-centre seconds from signal start
```
"""
struct SPLResult
    bands::Dict{Symbol, BandSPL}
    time::Vector{Float64}
    fs::Float32
    units::Symbol
    environment::Symbol

    function SPLResult(bands, time, fs, units, environment)
        # Validate units symbol before storing anything.
        units in (:dB_re_1µPa, :dB_re_20µPa) ||
            throw(ArgumentError(
                "SPLResult: units must be :dB_re_1µPa or :dB_re_20µPa; got :$units"))
        # Every band's spl_dB must have the same length as the time vector.
        # Checked here so callers get a clear error at construction rather
        # than a silent shape mismatch downstream.
        n = length(time)
        for (label, b) in bands
            length(b.spl_dB) == n ||
                throw(ArgumentError(
                    "SPLResult: band :$label has $(length(b.spl_dB)) SPL frames " *
                    "but time vector has $n entries"))
        end
        new(bands, time, fs, units, environment)
    end
end

# ─── compute_spl primitive ────────────────────────────────────────────────────

"""
    compute_spl(psd::PSDResult; bands, environment=:water) -> SPLResult

Purpose:     Compute band-integrated sound pressure level from a calibrated
             [`PSDResult`](@ref). For each requested frequency band, integrates
             the PSD over the band's bins, converts to dB SPL, and computes
             per-frame and aggregate statistics.

             Integration formula for each band [f_lo, f_hi]:
               band_power[j] = Σ_k psd[k,j] × df        (µPa²)
             where the sum runs over all bins k with f_lo ≤ freqs[k] ≤ f_hi
             and df = freqs[2] − freqs[1] (uniform bin width).

               SPL[j] = 10 × log10(band_power[j] / pref²)
             where pref = 1.0 µPa (water) or 20.0 µPa (air).

Arguments:
- `psd::PSDResult`: Calibrated PSD. Must satisfy
  `psd_units(psd) === :µPa²_per_Hz`. Use [`compute_psd`](@ref) with
  an appropriate calibration, or call `compute_spl(audio::Audiodata; ...)`
  which resolves calibration automatically.
- `bands::Dict{Symbol, Tuple{Float64,Float64}}`:
  Frequency bands to integrate. Each entry maps a label `Symbol` to a
  `(low_Hz, high_Hz)` tuple in Hz. There is no default — bands must
  always be supplied explicitly (DD-27). Overlapping bands are permitted;
  each is integrated independently. For standard band sets use the
  convenience wrappers: [`compute_tol`](@ref), [`compute_octave`](@ref),
  [`compute_millidecade`](@ref).
- `environment::Symbol = :water`:
  Acoustic medium. Determines the reference pressure:
  - `:water` → pref = 1 µPa, output units `:dB_re_1µPa`
  - `:air`   → pref = 20 µPa, output units `:dB_re_20µPa`

Returns:     [`SPLResult`](@ref) with one [`BandSPL`](@ref) per label.
             Each `BandSPL.spl_dB` has the same length as `psd.time`.

Constraints:
- `psd_units(psd)` must be `:µPa²_per_Hz`. Uncalibrated PSDs are in
  full-scale²/Hz — integrating them produces numbers with no acoustic
  meaning. The assertion fires with a message pointing to the fix (DD-21).
- PSD frequency bins must be uniformly spaced (guaranteed for PSD from
  `spectrogram`; asserted here to document the assumption).
- Band edges are validated before any integration:
  - `high_Hz > Nyquist` → `ArgumentError` listing all offending labels
  - `low_Hz ≥ high_Hz`  → `ArgumentError` listing all offending labels
  - `low_Hz < 10 Hz`    → consolidated `@warn` listing offending labels
    (below typical hydrophone response; may be intentional, so only warned)

Fails when:
- `psd_units(psd) ≠ :µPa²_per_Hz`         → `AssertionError` (DD-21)
- Any band has `high_Hz > psd.fs / 2`      → `ArgumentError`
- Any band has `low_Hz ≥ high_Hz`          → `ArgumentError`
- Any band has no PSD bins in range        → `ArgumentError` naming the band

Example:
```julia
psd    = compute_psd(audio; window_seconds = 1.0)
result = compute_spl(psd; bands = Dict(
    :full  => (10.0, Float64(psd.fs) / 2),
    :tonal => (18000.0, 22000.0)))
result.bands[:full].mean_dB    # energetic mean SPL, dB re 1 µPa
result.bands[:tonal].L99_dB   # 99th-percentile tonal SPL
```

Do not use when:
- `psd` is uncalibrated — the assertion will fire. Use the `Audiodata`
  wrapper or supply calibration explicitly to `compute_psd`.
- Standard band sets (octave, third-octave, millidecade) are needed —
  use [`compute_octave`](@ref), [`compute_tol`](@ref), or
  [`compute_millidecade`](@ref), which populate `bands` automatically.
- *System-weighted* broadband SPL (ADEON DPS Figure 1, left path) is
  required. This function implements the *frequency-flat* path only:
  the full calibrated PSD is integrated over the band. System-weighted
  SPL applies a single sensitivity value at a representative frequency
  (typically 250 Hz) directly to the time-domain signal and is not
  available in v1. For instruments with flat frequency response the two
  paths agree; for instruments with frequency-dependent response
  (e.g. Rockhopper TF calibration), the frequency-flat path is the
  physically correct choice. See DD-26.

References:
Merchant et al. (2015) Measuring Acoustic Habitats. Methods in Ecology and
Evolution, 6, 257–265.
"""
function compute_spl(psd::PSDResult;
                     bands::Dict{Symbol, Tuple{Float64, Float64}},
                     environment::Symbol = :water) :: SPLResult

    # DD-21: assert calibration before any computation. Uncalibrated PSDs are
    # in full-scale²/Hz; integrating them gives meaningless numbers.
    @assert psd_units(psd) === :µPa²_per_Hz (
        "compute_spl requires a calibrated PSD (DD-21). " *
        "psd_units(psd) must be :µPa²_per_Hz; got $(psd_units(psd)). " *
        "Call compute_psd with a calibration, or use " *
        "compute_spl(audio::Audiodata; ...) which resolves calibration automatically.")

    nyquist = Float64(psd.fs) / 2.0

    # ── Band validation (DD-19) ───────────────────────────────────────────────
    # Collect all offending labels before throwing, so the user sees every
    # problem in one error message rather than fixing them one at a time.
    invalid_order = Symbol[]
    above_nyquist = Symbol[]
    below_10hz    = Symbol[]

    for (label, (f_lo, f_hi)) in bands
        f_lo >= f_hi   && push!(invalid_order, label)
        f_hi > nyquist && push!(above_nyquist, label)
        f_lo < 10.0    && push!(below_10hz, label)
    end

    isempty(invalid_order) ||
        throw(ArgumentError(
            "compute_spl: band(s) with low_Hz ≥ high_Hz: " *
            join(sort(string.(invalid_order)), ", ")))
    isempty(above_nyquist) ||
        throw(ArgumentError(
            "compute_spl: band(s) with high_Hz > Nyquist ($nyquist Hz): " *
            join(sort(string.(above_nyquist)), ", ")))
    isempty(below_10hz) || @warn(
        "compute_spl: band(s) with low_Hz < 10 Hz — below typical hydrophone " *
        "response range; results may be unreliable. Affected bands: " *
        join(sort(string.(below_10hz)), ", "))

    # ── Uniform bin spacing ───────────────────────────────────────────────────
    # rfftfreq always produces uniform spacing. Asserted here to document the
    # assumption: compute_spl cannot be used with a non-uniform frequency axis.
    @assert length(psd.freqs) >= 2 "compute_spl: PSD has fewer than 2 frequency bins"
    df = psd.freqs[2] - psd.freqs[1]
    @assert all(d -> d ≈ df, diff(psd.freqs)) "compute_spl: non-uniform PSD bin spacing"

    # ── Reference pressure squared (µPa²) ────────────────────────────────────
    # water: pref = 1 µPa  → pref² = 1.0 µPa²
    # air:   pref = 20 µPa → pref² = 400.0 µPa²
    pref_sq = environment === :water ? 1.0 : 400.0

    # ── Per-band integration ──────────────────────────────────────────────────
    result_bands = Dict{Symbol, BandSPL}()
    for (label, (f_lo, f_hi)) in bands
        # First bin with centre ≥ f_lo; last bin with centre ≤ f_hi.
        i_lo = searchsortedfirst(psd.freqs, f_lo)
        i_hi = searchsortedlast(psd.freqs, f_hi)
        i_lo <= i_hi || throw(ArgumentError(
            "compute_spl: no PSD bins fall in band :$label " *
            "($f_lo Hz – $f_hi Hz; bin width = $df Hz)"))

        # Sum power over band bins and multiply by bin width.
        # @view avoids copying the row slice before summing — the SubArray is
        # read in-place. sum over dims=1 gives (1, n_frames); vec collapses to
        # Vector{Float64} of length n_frames. Units: µPa².
        power_per_frame = vec(sum(@view(psd.psd_linear[i_lo:i_hi, :]); dims=1)) .* df

        # Per-frame SPL.
        spl_dB = 10.0 .* log10.(power_per_frame ./ pref_sq)

        # Energetic mean: average in linear power domain, then convert to dB.
        # Arithmetic mean of dB is incorrect for signals with temporal variation
        # (DD-16); the energetic mean is physically correct.
        mean_dB = 10.0 * log10(mean(10.0 .^ (spl_dB ./ 10.0)))

        # Nine percentiles in a single quantile() call, matching DPS Table C-1.
        # qs[1]=L1, qs[2]=L5, qs[3]=L10, qs[4]=L25, qs[5]=L50(median),
        # qs[6]=L75, qs[7]=L90, qs[8]=L95, qs[9]=L99.
        qs = quantile(spl_dB, [0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99])

        result_bands[label] = BandSPL(
            (f_lo, f_hi), spl_dB, mean_dB,
            qs[5],   # median_dB  (L50)
            qs[1],   # L1_dB
            qs[2],   # L5_dB
            qs[3],   # L10_dB
            qs[4],   # L25_dB
            qs[6],   # L75_dB
            qs[7],   # L90_dB
            qs[8],   # L95_dB
            qs[9])   # L99_dB
    end

    units = environment === :water ? :dB_re_1µPa : :dB_re_20µPa
    return SPLResult(result_bands, psd.time, psd.fs, units, environment)
end

# ─── Convenience wrapper: Audiodata ──────────────────────────────────────────

"""
    compute_spl(audio::Audiodata; bands, environment=:water,
                window_seconds=1.0, overlap_fraction=0.5,
                window=:hann, nfft=nothing) -> SPLResult

Purpose:     Convenience wrapper: compute band-integrated SPL directly from
             an `Audiodata` object. Calls [`compute_psd`](@ref) internally,
             then delegates to the `compute_spl(psd::PSDResult; ...)` primitive.
             Calibration is resolved automatically via the three-step cascade
             in `compute_psd` (DD-14):
             1. `audio.is_calibrated == true` → no further calibration at PSD layer
             2. `audio.calibration isa !NoCalibration` → apply it at PSD layer
             3. `get_profile(Symbol(recorder))` → use recorder TF if available
                (handles Rockhopper automatically)
             4. No calibration found → `@warn`; `compute_spl` then asserts and fails

             This means an uncalibrated recording with no resolvable profile will
             fail at the `compute_spl` assertion (DD-21), not silently return
             meaningless SPL values.

Arguments:
- `audio::Audiodata`: Recording to analyse.
- `bands`: Forwarded to `compute_spl(psd; ...)`. Required — see that method
  for semantics (DD-27).
- `window_seconds::Real = 1.0`: Analysis window duration. Forwarded to
  [`compute_psd`](@ref) → [`spectrogram`](@ref).
- `overlap_fraction::Real = 0.5`: Frame overlap fraction. Forwarded to
  [`compute_psd`](@ref).
- `window::Symbol = :hann`: Window function. One of `:hann`, `:hamming`,
  `:blackman`, `:rectangular`. Forwarded to [`compute_psd`](@ref).
- `nfft::Union{Int,Nothing} = nothing`: FFT length. `nothing` uses
  `window_length` (no zero-padding). Forwarded to [`compute_psd`](@ref).
- `fft_plan = nothing`: Pre-built FFTW plan from [`make_spectrogram_plan`](@ref).
  Forwarded through `compute_psd` → `spectrogram`. Pass a plan when
  processing many chunks at the same sample rate and window size to avoid
  repeated FFTW wisdom lookups. The plan's FFT length must match the
  resolved `nfft`; a mismatch throws `AssertionError` (DD-04).

Returns:     [`SPLResult`](@ref). All fields identical to calling
             `compute_spl(compute_psd(audio; ...), bands=bands, environment=environment)`.

Constraints:
- The calibration cascade must resolve to a non-`NoCalibration` result,
  or `audio.is_calibrated` must be `true`. If neither holds, `compute_psd`
  emits a warning and returns an uncalibrated `PSDResult`, after which
  `compute_spl` asserts and throws `AssertionError` (DD-21).
- All constraints from [`spectrogram`](@ref) and [`compute_psd`](@ref) apply.

Fails when:  Same conditions as `compute_spl(psd::PSDResult; ...)` plus any
             failure modes of [`compute_psd`](@ref).

Example:
```julia
audio = read_audio("recording.flac"; recorder = "rockhopper")
result = compute_spl(audio;
                     window_seconds = 1.0,
                     bands = Dict(:full => (10.0, Float64(audio.fs) / 2)))
result.bands[:full].mean_dB   # energetic mean SPL, dB re 1 µPa
```

Do not use when:
- You already have a `PSDResult` — pass it directly to
  `compute_spl(psd; ...)` to avoid recomputing the spectrogram.
- Custom spectral parameters are needed alongside specific bands — build
  the `PSDResult` explicitly, inspect it, then call the primitive.
"""
function compute_spl(audio::Audiodata;
                     bands::Dict{Symbol, Tuple{Float64, Float64}},
                     environment::Symbol        = :water,
                     window_seconds::Real       = 1.0,
                     overlap_fraction::Real     = 0.5,
                     window::Symbol             = :hann,
                     nfft::Union{Int, Nothing}  = nothing,
                     fft_plan                   = nothing) :: SPLResult
    psd = compute_psd(audio; window_seconds, overlap_fraction, window, nfft, fft_plan)
    return compute_spl(psd; bands, environment)
end

# ─── compute_tol ──────────────────────────────────────────────────────────────

"""
    compute_tol(psd::PSDResult; low_Hz=10.0, high_Hz=fs/2, environment=:water) -> SPLResult
    compute_tol(audio::Audiodata; low_Hz=10.0, high_Hz=fs/2,
                environment=:water, window_seconds=1.0, overlap_fraction=0.5,
                window=:hann, nfft=nothing, fft_plan=nothing) -> SPLResult

Purpose:     Compute band-integrated SPL in ANSI S1.11 third-octave (decidecade)
             bands. Generates the band set via [`tol_bands`](@ref), then delegates
             to [`compute_spl`](@ref). The `PSDResult` method requires a
             pre-computed PSD; the `Audiodata` method computes the PSD internally.
             [`compute_decidecade`](@ref) is an alias for this function.

Arguments:
- `psd::PSDResult` or `audio::Audiodata`: Input data.
- `low_Hz::Real = 10.0`: Lower frequency bound. Bands with ANSI preferred
  center < `low_Hz` are excluded. Default 10 Hz matches the typical lower
  limit of hydrophone response.
- `high_Hz::Real = fs/2`: Upper frequency bound (Nyquist of the recording).
  Bands with ANSI preferred center > `high_Hz` are excluded. Defaults to
  Nyquist to clip at the recording's frequency limit automatically.
- `environment::Symbol = :water`: Reference pressure. `:water` → 1 µPa
  (`:dB_re_1µPa`); `:air` → 20 µPa (`:dB_re_20µPa`).
- `window_seconds`, `overlap_fraction`, `window`, `nfft`: Forwarded to
  [`compute_psd`](@ref) (Audiodata method only; not available on PSDResult method).
- `fft_plan`: Pre-built FFTW plan from [`make_spectrogram_plan`](@ref).
  Forwarded through `compute_psd` → `spectrogram` (Audiodata method only).
  Pass when processing many chunks at the same window size to avoid repeated
  FFTW wisdom lookups. The plan length must match the resolved `nfft` (DD-04).
  **Not available on the PSDResult method** — the FFT has already been computed.

Returns:     [`SPLResult`](@ref) with one [`BandSPL`](@ref) per ANSI S1.11
             band whose preferred center falls in `[low_Hz, high_Hz]`. Band
             label keys match [`tol_bands`](@ref) output (`:tol_N`, `:tol_12_5`,
             etc.). Returns a result with an empty `bands` Dict if no ANSI S1.11
             preferred centers fall in `[low_Hz, high_Hz]`.

Constraints:
- Same calibration requirement as [`compute_spl`](@ref): `psd_units(psd)` must
  be `:µPa²_per_Hz`. Uncalibrated input throws `AssertionError` (DD-21).
- Band edges validated by [`compute_spl`](@ref): `high_Hz > Nyquist` throws
  `ArgumentError`; `low_Hz < 10 Hz` triggers `@warn`.

Fails when:  Same conditions as [`compute_spl`](@ref).

Example:
```julia
psd    = compute_psd(audio; window_seconds = 1.0)
result = compute_tol(psd)
# Keys: :tol_10, :tol_12_5, :tol_16, ..., :tol_N (up to Nyquist)
result.bands[:tol_1000].mean_dB   # energetic mean 1-kHz third-octave SPL

# Or directly from Audiodata:
result = compute_tol(audio; window_seconds = 1.0)
```

Do not use when:
- Octave-band resolution suffices — use [`compute_octave`](@ref).
- Sub-Hz resolution is needed — use [`compute_millidecade`](@ref).
- Custom band boundaries are required — use [`compute_spl`](@ref) with an
  explicit `bands` Dict.

References:
ANSI S1.11-2004 (R2009). Specification for Octave-Band and Fractional-Octave-Band
Analog and Digital Filters. Acoustical Society of America.
ISO 18405:2017. Underwater Acoustics — Terminology. International Organization
for Standardization.
"""
function compute_tol(psd::PSDResult;
                     low_Hz::Real        = 10.0,
                     high_Hz::Real       = Float64(psd.fs) / 2.0,
                     environment::Symbol = :water) :: SPLResult
    return compute_spl(psd; bands = tol_bands(low_Hz, high_Hz), environment)
end

function compute_tol(audio::Audiodata;
                     low_Hz::Real             = 10.0,
                     high_Hz::Real            = Float64(audio.fs) / 2.0,
                     environment::Symbol      = :water,
                     window_seconds::Real     = 1.0,
                     overlap_fraction::Real   = 0.5,
                     window::Symbol           = :hann,
                     nfft::Union{Int,Nothing} = nothing,
                     fft_plan                 = nothing) :: SPLResult
    psd = compute_psd(audio; window_seconds, overlap_fraction, window, nfft, fft_plan)
    return compute_spl(psd; bands = tol_bands(low_Hz, high_Hz), environment)
end

# ─── compute_octave ───────────────────────────────────────────────────────────

"""
    compute_octave(psd::PSDResult; low_Hz=10.0, high_Hz=fs/2, environment=:water) -> SPLResult
    compute_octave(audio::Audiodata; low_Hz=10.0, high_Hz=fs/2,
                   environment=:water, window_seconds=1.0, overlap_fraction=0.5,
                   window=:hann, nfft=nothing, fft_plan=nothing) -> SPLResult

Purpose:     Compute band-integrated SPL in ANSI S1.6 octave bands. Generates
             the band set via [`octave_bands`](@ref), then delegates to
             [`compute_spl`](@ref). The `PSDResult` method requires a
             pre-computed PSD; the `Audiodata` method computes the PSD internally.

Arguments:
- `psd::PSDResult` or `audio::Audiodata`: Input data.
- `low_Hz::Real = 10.0`: Lower frequency bound. Bands with ANSI S1.6 preferred
  center < `low_Hz` are excluded. Default 10 Hz matches the typical lower
  limit of hydrophone response.
- `high_Hz::Real = fs/2`: Upper frequency bound (Nyquist). Bands with ANSI S1.6
  preferred center > `high_Hz` are excluded.
- `environment::Symbol = :water`: Reference pressure. `:water` → 1 µPa;
  `:air` → 20 µPa.
- `window_seconds`, `overlap_fraction`, `window`, `nfft`: Forwarded to
  [`compute_psd`](@ref) (Audiodata method only).
- `fft_plan`: Pre-built FFTW plan from [`make_spectrogram_plan`](@ref).
  Forwarded through `compute_psd` → `spectrogram` (Audiodata method only).
  **Not available on the PSDResult method** — the FFT has already been computed.

Returns:     [`SPLResult`](@ref) with one [`BandSPL`](@ref) per ANSI S1.6
             octave band whose preferred center falls in `[low_Hz, high_Hz]`.
             Band label keys match [`octave_bands`](@ref) output (`:oct_16`,
             `:oct_31_5`, `:oct_1000`, etc.).

Constraints:
- Same calibration requirement as [`compute_spl`](@ref): `psd_units(psd)`
  must be `:µPa²_per_Hz` (DD-21).
- `high_Hz` must not exceed Nyquist — passes through the band-edge validation
  in `compute_spl`, which throws `ArgumentError` if any band edge exceeds
  Nyquist.

Fails when:  Same conditions as [`compute_spl`](@ref).

Example:
```julia
psd    = compute_psd(audio; window_seconds = 1.0)
result = compute_octave(psd)
# Keys: :oct_16, :oct_31_5, :oct_63, ..., up to Nyquist
result.bands[:oct_1000].mean_dB   # energetic mean 1-kHz octave-band SPL

# Or directly from Audiodata:
result = compute_octave(audio; window_seconds = 1.0)
```

Do not use when:
- Third-octave or finer resolution is needed — use [`compute_tol`](@ref) or
  [`compute_millidecade`](@ref).
- Custom band boundaries are required — use [`compute_spl`](@ref) directly.

References:
ANSI S1.6-1984 (R2006). Preferred Frequencies, Frequency Levels, and Band
Numbers for Acoustical Measurements. Acoustical Society of America.
"""
function compute_octave(psd::PSDResult;
                        low_Hz::Real        = 10.0,
                        high_Hz::Real       = Float64(psd.fs) / 2.0,
                        environment::Symbol = :water) :: SPLResult
    return compute_spl(psd; bands = octave_bands(low_Hz, high_Hz), environment)
end

function compute_octave(audio::Audiodata;
                        low_Hz::Real             = 10.0,
                        high_Hz::Real            = Float64(audio.fs) / 2.0,
                        environment::Symbol      = :water,
                        window_seconds::Real     = 1.0,
                        overlap_fraction::Real   = 0.5,
                        window::Symbol           = :hann,
                        nfft::Union{Int,Nothing} = nothing,
                        fft_plan                 = nothing) :: SPLResult
    psd = compute_psd(audio; window_seconds, overlap_fraction, window, nfft, fft_plan)
    return compute_spl(psd; bands = octave_bands(low_Hz, high_Hz), environment)
end

# ─── compute_decidecade ───────────────────────────────────────────────────────

"""
    compute_decidecade(psd::PSDResult; kwargs...) -> SPLResult
    compute_decidecade(audio::Audiodata; kwargs...) -> SPLResult

Purpose:     Alias for [`compute_tol`](@ref). "Decidecade" (ISO 18405:2017) is
             the post-2018 underwater acoustics term for the same band scheme that
             ANSI S1.11 calls "third-octave". Output is identical to `compute_tol`
             with the same arguments.

Arguments:   Same as [`compute_tol`](@ref).
Returns:     Same as [`compute_tol`](@ref).
Constraints: Same as [`compute_tol`](@ref).
Fails when:  Same as [`compute_tol`](@ref).

Example:
```julia
compute_decidecade(psd) == compute_tol(psd)   # true (same band scheme)
```
"""
compute_decidecade(psd::PSDResult;   kw...) = compute_tol(psd;   kw...)
compute_decidecade(audio::Audiodata; kw...) = compute_tol(audio; kw...)

# ─── compute_millidecade ──────────────────────────────────────────────────────

"""
    compute_millidecade(psd::PSDResult; low_Hz=10.0, high_Hz=fs/2, environment=:water) -> SPLResult
    compute_millidecade(audio::Audiodata; low_Hz=10.0, high_Hz=fs/2,
                        environment=:water, window_seconds=1.0, overlap_fraction=0.5,
                        window=:hann, nfft=nothing, fft_plan=nothing) -> SPLResult

Purpose:     Compute band-integrated SPL in millidecade bands. Generates the
             band set via [`millidecade_bands`](@ref), then delegates to
             [`compute_spl`](@ref). The `PSDResult` method requires a pre-computed
             PSD; the `Audiodata` method computes the PSD internally.

             Millidecade bands are very fine resolution: near 1 kHz, each band
             spans roughly 2.3 Hz, producing ~87 bands in the 900–1100 Hz range
             versus 1 third-octave band. This resolution is the standard for
             NOAA NRS and MANTA long-term soundscape monitoring. For broadband
             surveys or comparisons with PAMGuide / Triton output, use
             [`compute_tol`](@ref) instead.

Arguments:
- `psd::PSDResult` or `audio::Audiodata`: Input data.
- `low_Hz::Real = 10.0`: Lower frequency bound. Must be > 0 (required by
  `millidecade_bands` for `log10`). Bands with center < `low_Hz` are excluded.
- `high_Hz::Real = fs/2`: Upper frequency bound (Nyquist). Bands with center
  > `high_Hz` are excluded.
- `environment::Symbol = :water`: Reference pressure. `:water` → 1 µPa;
  `:air` → 20 µPa.
- `window_seconds`, `overlap_fraction`, `window`, `nfft`: Forwarded to
  [`compute_psd`](@ref) (Audiodata method only).
- `fft_plan`: Pre-built FFTW plan from [`make_spectrogram_plan`](@ref).
  Forwarded through `compute_psd` → `spectrogram` (Audiodata method only).
  **Not available on the PSDResult method** — the FFT has already been computed.

Returns:     [`SPLResult`](@ref) with one [`BandSPL`](@ref) per millidecade
             band whose center falls in `[low_Hz, high_Hz]`. Band label keys
             match [`millidecade_bands`](@ref) output (`:mdec_N` where N is
             the integer index such that `f_c = 10^(N/1000)`). The number of
             bands scales with the frequency range: ~230 bands per decade, so
             a [10 Hz, 24 kHz] range yields ~3400 bands.

Constraints:
- Same calibration requirement as [`compute_spl`](@ref): `psd_units(psd)`
  must be `:µPa²_per_Hz` (DD-21).
- `low_Hz` must be > 0; `millidecade_bands` throws `ArgumentError` otherwise.
- For MANTA-compatible output, use the default `:mdec_N` label convention
  — do not rename band keys.

Fails when:
- Same conditions as [`compute_spl`](@ref).
- `low_Hz ≤ 0` → `ArgumentError` from [`millidecade_bands`](@ref).

Example:
```julia
psd    = compute_psd(audio; window_seconds = 1.0)
result = compute_millidecade(psd)
# Keys: :mdec_N for each millidecade band in [10 Hz, Nyquist]
result.bands[:mdec_3000].mean_dB   # energetic mean 1-kHz millidecade SPL

# Or directly from Audiodata:
result = compute_millidecade(audio; window_seconds = 1.0)
```

Do not use when:
- Third-octave resolution suffices — `compute_millidecade` produces ~14×
  more bands per decade than `compute_tol` and proportionally more output to
  store and process.
- Comparing against tools that use ANSI preferred-center third-octave bands —
  millidecade band edges do not align with ANSI S1.11 preferred centers.

References:
Miksis-Olds, J.L., et al. (2021). Ocean sound analysis software for making
ambient noise trends accessible (MANTA). Frontiers in Marine Science, 8.
Hatch, L.T., et al. (2016). Quantifying loss of acoustic communication space
for right whales in and around a U.S. national marine sanctuary. Conservation
Biology, 26(6), 983–994.
"""
function compute_millidecade(psd::PSDResult;
                              low_Hz::Real        = 10.0,
                              high_Hz::Real       = Float64(psd.fs) / 2.0,
                              environment::Symbol = :water) :: SPLResult
    return compute_spl(psd; bands = millidecade_bands(low_Hz, high_Hz), environment)
end

function compute_millidecade(audio::Audiodata;
                              low_Hz::Real             = 10.0,
                              high_Hz::Real            = Float64(audio.fs) / 2.0,
                              environment::Symbol      = :water,
                              window_seconds::Real     = 1.0,
                              overlap_fraction::Real   = 0.5,
                              window::Symbol           = :hann,
                              nfft::Union{Int,Nothing} = nothing,
                              fft_plan                 = nothing) :: SPLResult
    psd = compute_psd(audio; window_seconds, overlap_fraction, window, nfft, fft_plan)
    return compute_spl(psd; bands = millidecade_bands(low_Hz, high_Hz), environment)
end
