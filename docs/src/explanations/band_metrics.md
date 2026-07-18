# Band metrics and the per-cycle table

Soundscape analysis rarely uses the full spectrum directly. Instead it reduces
each recording to a handful of **frequency-band levels** — a 63 Hz shipping
indicator, a broadband 100–200 Hz summary, and so on — and tracks how those
levels behave over time. EcoAcoustics.jl provides this in two composable steps:

1. **Band integration** — sum power across a band's bins to get a band SPL series
   (`compute_spl` on a `PSDResult` or an `LTSAResult`).
2. **Per-cycle table** — flatten those band levels into one tidy row per
   recording and assemble a `DataFrame` for statistical analysis in R
   (`band_metrics`, `band_table`).

Both steps are pure composition over [`compute_psd`](@ref) and
[`compute_spl`](@ref) — no new signal processing (DD-29).

---

## Band integration: `compute_spl` on a PSD or an LTSA

`compute_spl` integrates each requested band out of a calibrated
frequency × time matrix. Both callers run the identical band-integration core
(`_integrate_bands`); they differ only in what a "column" of the series means:

| Input | A column is… | Series length | Use for |
|---|---|---|---|
| `PSDResult` | one FFT frame (~1 s) | hundreds per file | within-file distributions, percentiles |
| `LTSAResult` | one LTSA time column (e.g. 60 s) | one per column | long-recording band time series |

```julia
# Band time series across a long recording (one value per LTSA column):
lt  = compute_ltsa(audio; average_span_seconds = 60.0)
spl = compute_spl(lt; bands = Dict(:b100_200 => (100.0, 200.0)))
spl.bands[:b100_200].spl_dB     # per-minute 100–200 Hz level
```

Each band becomes a [`BandSPL`](@ref) holding the per-column series (`spl_dB`),
the energetic mean (`mean_dB`), and nine percentiles. The bands are supplied as a
`Dict{Symbol,Tuple{Float64,Float64}}` — there is no default band set (DD-27); you
always state the frequency ranges explicitly.

### Percentile direction (read this before interpreting Ln)

`L_n` is the **exceedance level** — the level *exceeded* n% of the columns —
following the ISO 18405 standard (ADEON, OSPAR/JOMOPANS; DD-31). Equivalently
`L_n` is the (100−n)th percentile:

- `L1`, `L5` → the **loud tail** (high levels; brief vessel passages).
- `L95`, `L99` → the **quiet background** (low levels; `L95` is the standard
  ambient indicator).

Note this is the **inverse of raw-percentile SPD-plot labels**, where a `5%`
line is the quiet 5th percentile (= EA's `L95`). See the [SPL page](spl.md) for
the full discussion and the suggested methods wording. When you want a
loud-event indicator, use `L1`, `L5`, or `max` — not `L95`.

---

## The per-cycle table

For archive-scale exploratory work — e.g. a three-month, duty-cycled deployment
where each cycle is one self-contained 5-minute file — the natural unit is
**one row per cycle**, with a set of band metrics as columns. That table is the
hand-off to the R-side statistical modelling.

### `band_metrics` — one row

```julia
audio = read_audio("cycle_0001.wav"; recorder = "sm3m")
bands = Dict(:tol_63 => (56.1, 70.7), :b100_200 => (100.0, 200.0))
row   = band_metrics(audio; bands = bands, window_seconds = 1.0)
```

`band_metrics` computes the file's PSD, integrates the bands, and flattens the
result into a `Dict{Symbol,Any}`:

- **Meta columns:** `start_time` (`DateTime`), `duration_s`, `n_frames` — the
  number of PSD frames the percentiles are taken over (a quality-control signal:
  a short or partial file shows up as few frames).
- **Band metric columns:** for every band `k` and metric `m`, a key
  `Symbol(k, :_, m, :_dB)` — e.g. `:tol_63_mean_dB`. The metrics are
  `mean, median, L1, L5, L10, L25, L75, L90, L95, L99, max`.

Percentiles are taken over the file's **PSD frames** (~1 s each), giving a rich
within-cycle distribution — fine enough that a single vessel pass inside a cycle
moves `L5`/`max` (the loud-tail exceedance levels).

### `band_table` — the DataFrame

```julia
files = readdir("deployment/"; join = true)
df = band_table(files;
                reader = f -> read_audio(f; recorder = "sm3m"),
                bands  = bands,
                window_seconds = 1.0)
```

`band_table` maps `band_metrics` over `items` and assembles a typed `DataFrame`,
one row per cycle:

- **`reader`** maps each item to an `Audiodata`. The read happens *inside* the
  per-item error guard, so corrupt/unreadable files are skipped
  (`on_error = :skip`, the default) rather than aborting the whole run. Pass
  `on_error = :fail` to stop on the first error.
- **Schema is fixed by `bands`**, so every row has the same columns regardless of
  which files succeeded. Columns are concretely typed (`DateTime`, `Float64`,
  `Int`) even when the table is empty.

Write it out for R with `CSV.write("bands.csv", df)`.

### SM3M per-era calibration

Absolute band levels are only comparable once the correct per-era hydrophone
sensitivity is applied. The SM3M profile currently encodes a single era
(−153 dB, pre-Apr-2018); for post-Apr-2018 files pass the correct calibration
explicitly via the `cal` keyword until per-era resolution is implemented:

```julia
band_metrics(audio; bands = bands, window_seconds = 1.0,
             cal = ScalarCalibration(-165.0f0))   # post-Apr-2018 era
```

See `CANDIDATE_ACOUSTIC_BANDS.org` and `ROADMAP.org` for the full data context,
and `examples/band_table_candidate_bands.jl` for the complete Chapter 1 workflow.

---

## Provisional parameters

The choices in the Chapter 1 pipeline — the candidate band set, the 1 s FFT
window, and taking percentiles over PSD frames — are the **test-run**
configuration, not fixed conventions. They are expected to change for production
runs once the exploratory analysis identifies which bands and timescales carry
signal (DD-29).

---

## References

Design decisions: DD-27 (explicit bands), DD-29 (band-integration core and the
per-cycle table) in `docs/design_decisions.md`. Composes the PSD, SPL, and LTSA
layers.
