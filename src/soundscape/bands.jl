# ─── Band-generator helpers ───────────────────────────────────────────────────
#
# Returns Dict{Symbol, Tuple{Float64, Float64}} mapping band labels to
# (low_Hz, high_Hz) integration edges for standard acoustic frequency band sets.
# Design decisions: DD-20 in docs/design_decisions.md.

# ── ANSI S1.6 preferred octave-band center frequencies (Hz) ──────────────────
# Source: ANSI S1.6-1984 (R2006), Table 1.
# These tabulated preferred values are used by PAMGuide, PAMGuard, MANTA, and
# Merchant 2015. DO NOT replace with exact formula values (1000 × 2^(n−10)):
# the preferred values intentionally round exact powers of two to cleaner
# numbers (e.g. 31.5 rather than 31.25, 63 rather than 62.5) to match the
# published standard. Using formula values would cause frequency-axis
# mismatches when comparing output against reference tools.
const OCTAVE_PREFERRED_HZ = Float64[
    0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 31.5, 63.0, 125.0, 250.0, 500.0,
    1000.0, 2000.0, 4000.0, 8000.0, 16000.0, 31500.0, 63000.0, 125000.0,
    250000.0, 500000.0
]

# ── ANSI S1.11 preferred third-octave (decidecade) center frequencies (Hz) ───
# Source: ANSI S1.11-2004 (R2009), Table 1.
# Same rationale as OCTAVE_PREFERRED_HZ above: tabulated preferred values, not
# exact formula outputs (1000 × 2^((n−30)/3)). These match the conventions used
# by Merchant 2015, PAMGuide, PAMGuard, and MANTA.
const TOL_PREFERRED_HZ = Float64[
    1.0, 1.25, 1.6, 2.0, 2.5, 3.15, 4.0, 5.0, 6.3, 8.0,
    10.0, 12.5, 16.0, 20.0, 25.0, 31.5, 40.0, 50.0, 63.0, 80.0,
    100.0, 125.0, 160.0, 200.0, 250.0, 315.0, 400.0, 500.0, 630.0, 800.0,
    1000.0, 1250.0, 1600.0, 2000.0, 2500.0, 3150.0, 4000.0, 5000.0, 6300.0, 8000.0,
    10000.0, 12500.0, 16000.0, 20000.0, 25000.0, 31500.0, 40000.0, 50000.0,
    63000.0, 80000.0, 100000.0, 125000.0, 160000.0, 200000.0
]

# ── Band-edge multipliers ─────────────────────────────────────────────────────
# For a band of width 1/N octaves, edges = f_c / 2^(1/(2N)) and f_c × 2^(1/(2N)).
# N=1 (octave):       multiplier = 2^(1/2) = √2  ≈ 1.4142
# N=3 (third-octave): multiplier = 2^(1/6)       ≈ 1.1225
const _OCTAVE_EDGE_FACTOR = 2.0^(1/2)
const _TOL_EDGE_FACTOR    = 2.0^(1/6)

# ── Symbol helper ─────────────────────────────────────────────────────────────
# Purpose:  Build a band label Symbol from a prefix string and a preferred
#           center frequency. Integer-valued centers produce :prefix_1000;
#           non-integer centers replace the decimal point with an underscore
#           to produce valid Symbols: :prefix_12_5, :prefix_31_5.
# Constraints: f_c must be a finite Float64. Only called with values from
#              OCTAVE_PREFERRED_HZ or TOL_PREFERRED_HZ.
function _band_label(prefix::String, f_c::Float64) :: Symbol
    if isinteger(f_c)
        Symbol(prefix * "_" * string(round(Int, f_c)))
    else
        # Julia's string() uses shortest-roundtrip representation, so 12.5
        # prints as "12.5", 31.5 as "31.5", etc. Replace "." → "_".
        Symbol(prefix * "_" * replace(string(f_c), "." => "_"))
    end
end

# ─── octave_bands ─────────────────────────────────────────────────────────────

"""
    octave_bands(low_Hz, high_Hz) -> Dict{Symbol, Tuple{Float64, Float64}}

Purpose:     Return a Dict of ANSI S1.6 octave-band integration edges for every
             band whose preferred center frequency falls in `[low_Hz, high_Hz]`.
             Suitable for passing directly to [`compute_spl`](@ref) or
             [`compute_octave`](@ref).

Arguments:
- `low_Hz::Real`:  Lower bound in Hz. Bands with center < `low_Hz` are excluded.
- `high_Hz::Real`: Upper bound in Hz. Bands with center > `high_Hz` are excluded.
  Pass `fs/2` (Nyquist) to clip at the recording's frequency limit.

Returns:     `Dict{Symbol, Tuple{Float64, Float64}}`. Each key is a Symbol of the
             form `:oct_N` where `N` is the ANSI preferred center frequency as a
             bare integer (`:oct_1000`) or with a decimal underscore for
             non-integer preferred values (`:oct_31_5`). Each value is
             `(low_Hz, high_Hz)` computed as
             `(f_c / 2^(1/2), f_c × 2^(1/2))` — the ANSI S1.6 band edges.
             Returns an empty Dict if no bands fall in range.

Constraints:
- Center frequencies are the ANSI S1.6 preferred values, not exact formula
  outputs. This matches PAMGuide, PAMGuard, MANTA, and Merchant 2015.
- Band edges are derived from the preferred center, not from adjacent centers:
  edges may slightly overlap or gap between adjacent bands, as is standard.
- The caller is responsible for Nyquist clipping. Pass `high_Hz = fs/2` to
  exclude bands above the recording's frequency limit. No Nyquist check is
  applied inside this function.

Fails when:  Never. Returns empty Dict for ranges with no matching bands.

Example:
```julia
bands = octave_bands(10.0, 10000.0)   # 10 bands: 16 Hz through 8000 Hz
# keys: :oct_16, :oct_31_5, :oct_63, :oct_125, ..., :oct_8000
bands[:oct_1000]   # (707.1..., 1414.2...)  = (1000/√2, 1000×√2)
```

Do not use when:
- Third-octave resolution is required — use [`tol_bands`](@ref) or
  [`decidecade_bands`](@ref).
- Continuous frequency coverage at high resolution is required — use
  [`millidecade_bands`](@ref).

References:
ANSI S1.6-1984 (R2006). Preferred Frequencies, Frequency Levels, and Band
Numbers for Acoustical Measurements. Acoustical Society of America.
"""
function octave_bands(low_Hz::Real, high_Hz::Real) :: Dict{Symbol, Tuple{Float64, Float64}}
    result = Dict{Symbol, Tuple{Float64, Float64}}()
    for f_c in OCTAVE_PREFERRED_HZ
        low_Hz <= f_c <= high_Hz || continue
        result[_band_label("oct", f_c)] = (f_c / _OCTAVE_EDGE_FACTOR,
                                           f_c * _OCTAVE_EDGE_FACTOR)
    end
    return result
end

# ─── tol_bands ────────────────────────────────────────────────────────────────

"""
    tol_bands(low_Hz, high_Hz) -> Dict{Symbol, Tuple{Float64, Float64}}

Purpose:     Return a Dict of ANSI S1.11 third-octave (decidecade) band
             integration edges for every band whose preferred center frequency
             falls in `[low_Hz, high_Hz]`. Third-octave bands are the standard
             choice for soundscape reporting, biological call-level measurements,
             and most PAMGuide-family analyses.

Arguments:
- `low_Hz::Real`:  Lower bound in Hz. Bands with center < `low_Hz` are excluded.
- `high_Hz::Real`: Upper bound in Hz. Bands with center > `high_Hz` are excluded.
  Pass `fs/2` to clip at Nyquist.

Returns:     `Dict{Symbol, Tuple{Float64, Float64}}`. Keys are `:tol_N` Symbols
             using the ANSI preferred center frequency. Non-integer preferred
             values use an underscore decimal separator: `:tol_12_5`, `:tol_31_5`.
             Band edges are `(f_c / 2^(1/6), f_c × 2^(1/6))` — the ANSI S1.11
             third-octave edges derived from the preferred center.

Constraints:
- Center frequencies are the ANSI S1.11-2004 preferred values.
- `tol_bands` and [`decidecade_bands`](@ref) are identical. The term
  "decidecade" (ISO 18405:2017) superseded "third-octave" in the underwater
  acoustics literature around 2018; older literature (including Merchant 2015)
  uses "third-octave". Both refer to 1/3-decade frequency bands with the same
  preferred center frequencies and edge formula.
- **Band-edge formula differs from the ADEON DPS** (Ainslie et al. 2018):
  `tol_bands` uses ANSI S1.11 base-2 edges (`f_c × 2^(±1/6) ≈ f_c × 1.1225`);
  the DPS specifies base-10 edges (`f_c × 10^(±1/20) ≈ f_c × 1.1220`). The
  relative difference is ~0.04% (~0.04 Hz at 1 kHz), negligible for integer-Hz
  PSD resolution. This deviation is deliberate and documented in DD-20; using
  ANSI preferred centers ensures label compatibility with PAMGuide, MANTA, and
  Merchant 2015.
- Caller is responsible for Nyquist clipping.

Fails when:  Never. Returns empty Dict for ranges with no matching bands.

Example:
```julia
bands = tol_bands(10.0, 1000.0)
# 17 bands: :tol_10, :tol_12_5, :tol_16, ..., :tol_1000
bands[:tol_1000]   # (890.9..., 1122.5...)  = (1000/2^(1/6), 1000×2^(1/6))
```

Do not use when:
- Octave-band resolution is sufficient — use [`octave_bands`](@ref).
- High-resolution (sub-1-Hz bandwidth) SPL is needed — use
  [`millidecade_bands`](@ref).

References:
ANSI S1.11-2004 (R2009). Specification for Octave-Band and Fractional-Octave-Band
Analog and Digital Filters. Acoustical Society of America.
ISO 18405:2017. Underwater Acoustics — Terminology. International Organization
for Standardization.
"""
function tol_bands(low_Hz::Real, high_Hz::Real) :: Dict{Symbol, Tuple{Float64, Float64}}
    result = Dict{Symbol, Tuple{Float64, Float64}}()
    for f_c in TOL_PREFERRED_HZ
        low_Hz <= f_c <= high_Hz || continue
        result[_band_label("tol", f_c)] = (f_c / _TOL_EDGE_FACTOR,
                                           f_c * _TOL_EDGE_FACTOR)
    end
    return result
end

# ─── decidecade_bands ─────────────────────────────────────────────────────────

"""
    decidecade_bands(low_Hz, high_Hz) -> Dict{Symbol, Tuple{Float64, Float64}}

Purpose:     Alias for [`tol_bands`](@ref). Returns ANSI S1.11 third-octave
             (decidecade) bands with `:tol_N` label keys. The "decidecade"
             (ISO 18405:2017) and "third-octave" (ANSI S1.11) terms refer to the
             same band scheme; the former is preferred in the post-2018 underwater
             acoustics literature. Output is identical to `tol_bands`.

Arguments:   Same as [`tol_bands`](@ref).
Returns:     Same as [`tol_bands`](@ref).
Constraints: Same as [`tol_bands`](@ref), including the ADEON DPS base-10 edge
             deviation documented there and in DD-20.
Fails when:  Never.

Example:
```julia
decidecade_bands(10.0, 1000.0) == tol_bands(10.0, 1000.0)   # true
```
"""
decidecade_bands(low_Hz::Real, high_Hz::Real) = tol_bands(low_Hz, high_Hz)

# ─── millidecade_bands ────────────────────────────────────────────────────────

"""
    millidecade_bands(low_Hz, high_Hz) -> Dict{Symbol, Tuple{Float64, Float64}}

Purpose:     Return a Dict of millidecade band integration edges for every band
             whose center frequency falls in `[low_Hz, high_Hz]`. A millidecade
             is 1/1000 of a decade on a log-frequency axis; each band spans a
             frequency ratio of 10^(1/1000) ≈ 1.0023 (about 0.23%). This is the
             high-resolution band scheme used by NOAA's National Recording System
             (NRS) and MANTA for long-term soundscape monitoring.

             Resolution note: millidecade bands are very fine. Near 1 kHz, each
             band spans roughly 2.3 Hz. The range [900, 1100] Hz contains ~87
             bands, compared with 1 third-octave band. This resolution makes
             millidecade bands suitable for high-resolution PSD aggregation across
             long archives while retaining stable, integer-indexed band identities
             — the property that MANTA exploits.

Arguments:
- `low_Hz::Real`:  Lower bound in Hz. Must be > 0. Bands with center < `low_Hz`
  are excluded.
- `high_Hz::Real`: Upper bound in Hz. Bands with center > `high_Hz` are excluded.
  Pass `fs/2` to clip at Nyquist.

Returns:     `Dict{Symbol, Tuple{Float64, Float64}}`. Keys are `:mdec_N` where
             `N` is the integer band index such that `f_c = 10^(N/1000)`.
             `:mdec_3000` has center `10^3 = 1000 Hz`; `:mdec_4000` has center
             `10^4 = 10000 Hz`. This MANTA convention is unambiguous — adjacent
             bands never have the same center Hz when rounded to an integer.
             Band edges are `(10^((N−0.5)/1000), 10^((N+0.5)/1000))`.

Constraints:
- `low_Hz` must be > 0 (required for `log10`). Throws `ArgumentError` otherwise.
- `high_Hz` must be > 0. Throws `ArgumentError` otherwise.
- Returns empty Dict if no band centers fall in `[low_Hz, high_Hz]`.
- Caller is responsible for Nyquist clipping; pass `fs/2` as `high_Hz`.

Fails when:
- `low_Hz ≤ 0`  → `ArgumentError`
- `high_Hz ≤ 0` → `ArgumentError`

Example:
```julia
bands = millidecade_bands(900.0, 1100.0)   # 87 bands
bands[:mdec_3000]    # center = 10^3 = 1000 Hz; ≈ (999.9, 1000.1) Hz
length(bands)        # 87
```

Do not use when:
- Third-octave resolution suffices — `millidecade_bands` produces ~14× more
  bands per decade than `tol_bands` and proportionally more results to store.
- Comparing against tools that use ANSI preferred-center third-octave bands —
  millidecade band edges do not align with ANSI S1.11 preferred centers.

References:
Hatch, L.T., et al. (2016). Quantifying loss of acoustic communication space
for right whales in and around a U.S. national marine sanctuary. Conservation
Biology, 26(6), 983–994. Early use of millidecade bands.

Miksis-Olds, J.L., et al. (2021). Ocean sound analysis software for making
ambient noise trends accessible (MANTA). Frontiers in Marine Science, 8.
MANTA millidecade band convention and index scheme.
"""
function millidecade_bands(low_Hz::Real, high_Hz::Real) :: Dict{Symbol, Tuple{Float64, Float64}}
    low_Hz  > 0 || throw(ArgumentError(
        "millidecade_bands: low_Hz must be > 0; got $low_Hz"))
    high_Hz > 0 || throw(ArgumentError(
        "millidecade_bands: high_Hz must be > 0; got $high_Hz"))

    # Find the smallest and largest integer indices n such that
    # 10^(n/1000) ∈ [low_Hz, high_Hz].
    n_lo = ceil(Int,  1000.0 * log10(Float64(low_Hz)))
    n_hi = floor(Int, 1000.0 * log10(Float64(high_Hz)))

    result = Dict{Symbol, Tuple{Float64, Float64}}()
    for n in n_lo:n_hi
        result[Symbol("mdec_$n")] = (10.0^((n - 0.5) / 1000.0),
                                     10.0^((n + 0.5) / 1000.0))
    end
    return result
end
