# SPL Band Schemes

EcoAcoustics.jl provides four standard frequency band schemes for SPL
computation. Each scheme is a `Dict{Symbol, Tuple{Float64, Float64}}` mapping
a band label to `(low_Hz, high_Hz)` integration edges. The same Dict format
is accepted by `compute_spl` for custom band sets, so all four are
interchangeable at the API level.

---

## Overview: which scheme to use

| Scheme        | Function              | Bands in [10 Hz, 24 kHz] | Use case |
|:--------------|:----------------------|:-------------------------|:---------|
| Octave        | `octave_bands`        | ~10 bands                | Coarse overview; regulatory reporting; comparison with older literature |
| Third-octave  | `tol_bands`           | ~30 bands                | Standard soundscape reporting; PAMGuide / Triton output; Merchant 2015 |
| Decidecade    | `decidecade_bands`    | identical to above       | Same as third-octave; preferred term in post-2018 underwater acoustics |
| Millidecade   | `millidecade_bands`   | ~3400 bands              | High-resolution PSD aggregation; MANTA / NOAA NRS output |

When in doubt, use third-octave / decidecade. It is the most common output
format in published soundscape literature and the default unit for comparing
recordings across sites and instruments.

---

## ANSI preferred center frequencies (DD-20)

**Do not replace the tabulated preferred values with exact formula values.**
This is the most important constraint in this file.

Both octave and third-octave band sets use the ANSI tabulated preferred center
frequencies — not the values computed from the exact band-number formula. The
formulas are:

```
octave centers:        1000 × 2^(n − 10)     for n = 0, 1, 2, …
third-octave centers:  1000 × 2^((n − 30)/3) for n = 0, 1, 2, …
```

The ANSI preferred values round these exact results to cleaner numbers. For
example:
- Exact formula: `1000 × 2^(-5)` = 31.25 Hz → **preferred: 31.5 Hz**
- Exact formula: `1000 × 2^(-4)` = 62.5 Hz → **preferred: 63 Hz**
- Exact formula: `1000 × 2^(-5/3)` ≈ 3.155 Hz → **preferred: 3.15 Hz**

PAMGuide, PAMGuard, MANTA, and Merchant (2015) all use the preferred values.
If EcoAcoustics.jl used exact formula values instead, frequency-axis labels
would not match reference tools, and cross-tool comparisons would require a
manual mapping step.

The preferred values are stored in the module-level constants
`OCTAVE_PREFERRED_HZ` and `TOL_PREFERRED_HZ` in `src/soundscape/bands.jl`.

---

## Band edges

Integration edges are derived from the band center, not from adjacent centers.
Octave bands use the **base-2** convention; third-octave / decidecade bands use
the **base-10** convention with exact decidecade centers `f_c = 10^(n/10)`
(DD-30):

```
octave band:       f_lo = f_c / 2^(1/2),    f_hi = f_c × 2^(1/2)    (±√2, base-2)
decidecade band:   f_lo = f_c / 10^(1/20),  f_hi = f_c × 10^(1/20)  (≈ ±1.1220, base-10)
```

For decidecade bands the label uses the ANSI preferred (nominal) center for
readability (`:tol_63`), while the edges use the exact base-10 center
(63.096 Hz) — the ISO 18405 / ADEON / MANTA decidecade definition. Band levels
are computed by **frequency-domain spectral integration** of the PSD. This is the
same method family ADEON and MANTA use — verified against the primary sources
(Martin et al. 2021 hybrid millidecade; ADEON Data Processing Specification),
which compute decidecades by summing / integrating the PSD in the frequency
domain, not with a filter bank (DD-30). EA matches PAMGuide's PSD integrated into
the same bands to ~0.001 dB.

Two differences are worth knowing (both in DD-30):

1. **PAMGuide's dedicated `PG_TOL`** third-octave output uses the older
   **filter-bank** method (time-domain IIR octave filters), which differs from
   spectral integration by up to ~1.7 dB at low frequencies (63–315 Hz) where the
   filter skirts span a large fraction of the band.
2. **ADEON** splits edge bins that straddle two bands **fractionally** by percent
   overlap; EA uses **hard** bin edges (a bin is fully in or out — DD-25). So EA
   matches the ADEON *approach* but is not bit-identical at band edges; the
   difference is small and concentrated at low frequency.

Neither is an error — they are method / edge-handling choices.

Because edges are computed from centers rather than from adjacent centers,
adjacent bands may **slightly overlap or gap**. This is standard practice and is
present in PAMGuide and MANTA output. It is not a bug.

---

## Decidecade vs. third-octave (DD-20)

`decidecade_bands` is a pure alias for `tol_bands`. The two terms refer to
the same band scheme:
- **Third-octave** (ANSI S1.11): traditional term used in air acoustics and
  in Merchant 2015 and all PAMGuide-family literature.
- **Decidecade** (ISO 18405:2017): a 1/10-decade frequency band. This term
  entered the underwater acoustics literature around 2018 and is now preferred
  in that community. EA implements the base-10 decidecade (DD-30), so
  `tol_bands` and `decidecade_bands` are the same base-10 band — the near-match
  `log10(2^(1/3)) ≈ 0.1003` is *not* exact and the ~0.08% edge difference is
  material at low frequencies (see DD-30).

Both functions produce the same output. Choose whichever term matches your
literature. `compute_decidecade` delegates directly to `compute_tol` with no
overhead.

---

## Symbol naming for non-integer centers

Band label Symbols cannot contain dots. Non-integer preferred centers use an
underscore decimal separator:

| Preferred center | Symbol       |
|:-----------------|:-------------|
| 12.5 Hz          | `:tol_12_5`  |
| 31.5 Hz          | `:oct_31_5`  |
| 3.15 Hz          | `:tol_3_15`  |
| 6.3 Hz           | `:tol_6_3`   |
| 1000 Hz          | `:tol_1000`  |

The `_band_label` helper in `bands.jl` uses Julia's shortest-roundtrip
`string()` representation (e.g., `12.5` prints as `"12.5"`, not `"12.50"`),
then replaces `"."` with `"_"`. The key concern is round-trippability: the
string must reconstruct the original float without ambiguity.

---

## Millidecade bands: resolution and index convention

A **millidecade** is 1/1000 of a decade on a log-frequency axis. Each band
spans a frequency ratio of `10^(1/1000) ≈ 1.00230` — about 0.23%.

### Resolution is very fine

Near 1 kHz, each millidecade band spans approximately 2.3 Hz. The 900–1100 Hz
range contains **87 bands** compared with one third-octave band covering the
same range. Over a [10 Hz, 24 kHz] range — a typical 48 kHz recording — there
are roughly **3,400 millidecade bands**. This is the design intent: MANTA and
the NOAA NRS use millidecade bands specifically because the fine resolution
allows high-resolution PSD archives to be stored in a compact, stable,
integer-indexed format.

Before using `compute_millidecade`, consider:
- **Storage**: 3,400 `BandSPL` entries per time step, each with a 5-element
  vector of aggregate statistics, adds up quickly for long archives.
- **Frequency resolution requirement**: a millidecade band near 10 Hz spans
  only ~0.023 Hz. A PSD with `df = 1 Hz` (1-second window at 1 kHz) is far
  too coarse — many bands will contain zero bins and `compute_spl` will throw.
  For millidecade bands at 10 Hz you need `df < 0.02 Hz`, which requires a
  window of at least 50 seconds. In practice, millidecade analyses use long
  windows (typically 60 s) or restrict the low-frequency cutoff to where the
  frequency resolution is adequate.

### MANTA index convention

Millidecade bands use an integer index `n` such that the center frequency is:

```
f_c(n) = 10^(n / 1000)
```

Band labels are `:mdec_N` where `N` is this integer:
- `:mdec_3000` → `f_c = 10^3 = 1000 Hz`
- `:mdec_4000` → `f_c = 10^4 = 10000 Hz`
- `:mdec_2000` → `f_c = 10^2 = 100 Hz`

This MANTA convention is unambiguous: no two adjacent bands round to the same
center Hz when expressed as an integer. Using center Hz in the symbol
(`:mdec_1000`) would alias once Hz-rounded centers of adjacent bands coincide.

Band edges:
```
f_lo(n) = 10^((n − 0.5) / 1000)
f_hi(n) = 10^((n + 0.5) / 1000)
```

---

## Nyquist clipping

Band generators clip bands by **center frequency**, not by band edge. Passing
`high_Hz = fs/2` (Nyquist) to a generator may include a band centered exactly
at Nyquist, whose upper edge will exceed Nyquist. `compute_spl` will reject
that band with `ArgumentError`.

When using the convenience wrappers (`compute_tol`, `compute_octave`, etc.),
choose `high_Hz` such that the uppermost included center's upper edge stays
within Nyquist. For third-octave (decidecade) bands, the safe upper limit is
`Nyquist / 10^(1/20) ≈ 0.891 × Nyquist`. For a 48 kHz recording
(Nyquist = 24 kHz), this is ~21.4 kHz — so the 20 kHz band is the topmost
safe third-octave band.

The wrappers default to `high_Hz = Float64(psd.fs) / 2.0`, which is correct
when Nyquist does not coincide with an ANSI preferred center. For the common
case (fs = 48, 44.1, 96, 192 kHz), no ANSI preferred center sits exactly at
Nyquist and the default works without adjustment.

---

## References

ANSI S1.6-1984 (R2006). Preferred Frequencies, Frequency Levels, and Band
Numbers for Acoustical Measurements. Acoustical Society of America.

ANSI S1.11-2004 (R2009). Specification for Octave-Band and Fractional-Octave-Band
Analog and Digital Filters. Acoustical Society of America.

ISO 18405:2017. Underwater Acoustics — Terminology. International Organization
for Standardization.

Miksis-Olds, J.L., et al. (2021). Ocean sound analysis software for making
ambient noise trends accessible (MANTA). *Frontiers in Marine Science*, 8.
[DOI:10.3389/fmars.2021.703650](https://doi.org/10.3389/fmars.2021.703650)

Design decisions: DD-20 in `docs/design_decisions.md`.
