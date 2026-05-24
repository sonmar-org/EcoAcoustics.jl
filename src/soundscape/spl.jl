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
result = compute_spl(psd; bands = Dict(:broadband => (10.0, 24000.0)))
result.bands[:broadband].mean_dB    # energetic mean, dB re 1 µPa
result.bands[:broadband].spl_dB     # per-frame time series
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
    compute_spl(psd::PSDResult; bands=nothing, environment=:water) -> SPLResult

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
- `bands::Union{Nothing, Dict{Symbol, Tuple{Float64,Float64}}} = nothing`:
  Frequency bands to integrate. Each entry maps a label `Symbol` to a
  `(low_Hz, high_Hz)` tuple. When `nothing`, a single broadband band
  from 10 Hz to Nyquist is used: `Dict(:broadband => (10.0, fs/2))`.
  Overlapping bands are permitted; each is integrated independently.
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
    :broadband => (10.0, 24000.0),
    :tonal     => (18000.0, 22000.0)))
result.bands[:broadband].mean_dB   # energetic mean broadband SPL, dB re 1 µPa
result.bands[:tonal].L99_dB        # 99th-percentile tonal SPL
```

Do not use when:
- `psd` is uncalibrated — the assertion will fire. Use the `Audiodata`
  wrapper or supply calibration explicitly to `compute_psd`.
- Standard band sets (octave, third-octave, millidecade) are needed —
  use [`compute_octave`](@ref), [`compute_tol`](@ref), or
  [`compute_millidecade`](@ref), which populate `bands` automatically.

References:
Merchant et al. (2015) Measuring Acoustic Habitats. Methods in Ecology and
Evolution, 6, 257–265.
"""
function compute_spl(psd::PSDResult;
                     bands::Union{Nothing, Dict{Symbol, Tuple{Float64, Float64}}} = nothing,
                     environment::Symbol = :water) :: SPLResult

    # DD-21: assert calibration before any computation. Uncalibrated PSDs are
    # in full-scale²/Hz; integrating them gives meaningless numbers.
    @assert psd_units(psd) === :µPa²_per_Hz (
        "compute_spl requires a calibrated PSD (DD-21). " *
        "psd_units(psd) must be :µPa²_per_Hz; got $(psd_units(psd)). " *
        "Call compute_psd with a calibration, or use " *
        "compute_spl(audio::Audiodata; ...) which resolves calibration automatically.")

    # Default band: broadband from 10 Hz to Nyquist (DD-19).
    nyquist = Float64(psd.fs) / 2.0
    resolved_bands = bands === nothing ?
        Dict{Symbol, Tuple{Float64, Float64}}(:broadband => (10.0, nyquist)) : bands

    # ── Band validation (DD-19) ───────────────────────────────────────────────
    # Collect all offending labels before throwing, so the user sees every
    # problem in one error message rather than fixing them one at a time.
    invalid_order = Symbol[]
    above_nyquist = Symbol[]
    below_10hz    = Symbol[]

    for (label, (f_lo, f_hi)) in resolved_bands
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
    for (label, (f_lo, f_hi)) in resolved_bands
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

        # All five percentiles in a single quantile() call.
        # qs[1]=L1, qs[2]=L5, qs[3]=median(L50), qs[4]=L95, qs[5]=L99.
        qs = quantile(spl_dB, [0.01, 0.05, 0.50, 0.95, 0.99])

        result_bands[label] = BandSPL(
            (f_lo, f_hi), spl_dB, mean_dB,
            qs[3],   # median_dB
            qs[1],   # L1_dB
            qs[2],   # L5_dB
            qs[4],   # L95_dB
            qs[5])   # L99_dB
    end

    units = environment === :water ? :dB_re_1µPa : :dB_re_20µPa
    return SPLResult(result_bands, psd.time, psd.fs, units, environment)
end

# ─── Convenience wrapper: Audiodata ──────────────────────────────────────────

"""
    compute_spl(audio::Audiodata; bands=nothing, environment=:water,
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
- `bands`, `environment`: Forwarded to `compute_spl(psd; ...)`. See that method
  for semantics and defaults.
- `window_seconds::Real = 1.0`: Analysis window duration. Forwarded to
  [`compute_psd`](@ref) → [`spectrogram`](@ref).
- `overlap_fraction::Real = 0.5`: Frame overlap fraction. Forwarded to
  [`compute_psd`](@ref).
- `window::Symbol = :hann`: Window function. One of `:hann`, `:hamming`,
  `:blackman`, `:rectangular`. Forwarded to [`compute_psd`](@ref).
- `nfft::Union{Int,Nothing} = nothing`: FFT length. `nothing` uses
  `window_length` (no zero-padding). Forwarded to [`compute_psd`](@ref).

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
result = compute_spl(audio; window_seconds = 1.0)   # auto-calibrates
result.bands[:broadband].mean_dB   # energetic mean broadband SPL, dB re 1 µPa
```

Do not use when:
- You already have a `PSDResult` — pass it directly to
  `compute_spl(psd; ...)` to avoid recomputing the spectrogram.
- Custom spectral parameters are needed alongside specific bands — build
  the `PSDResult` explicitly, inspect it, then call the primitive.
"""
function compute_spl(audio::Audiodata;
                     bands::Union{Nothing, Dict{Symbol, Tuple{Float64, Float64}}} = nothing,
                     environment::Symbol        = :water,
                     window_seconds::Real       = 1.0,
                     overlap_fraction::Real     = 0.5,
                     window::Symbol             = :hann,
                     nfft::Union{Int, Nothing}  = nothing) :: SPLResult
    psd = compute_psd(audio; window_seconds, overlap_fraction, window, nfft)
    return compute_spl(psd; bands, environment)
end
