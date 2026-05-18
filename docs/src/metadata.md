# Filename parsing, recorder profiles, and calibration

EcoAcoustics.jl keeps three concerns separate:

1. **Reading audio samples** — `read_audio` handles file I/O.
2. **Understanding filenames** — `parse_filename` + `RecorderProfile` handle
   recorder-specific naming conventions.
3. **Calibration** — `CalibrationProfile` + `lookup_calibration` supply the
   physical sensitivity constants. Calibration is stored in `Audiodata` but
   never applied at I/O time.

---

## Recorder profiles: one file per recorder

Each supported recorder family has its own file under `src/recorders/`:

```
src/recorders/
├── recorders.jl      # RecorderProfile and CalibrationProfile structs,
│                     # RECORDER_PROFILES and CALIBRATION_PROFILES dicts,
│                     # lookup_calibration function
├── rockhopper.jl     # Rockhopper (Embedded Ocean Systems)
├── sm3m.jl           # SM3M (Wildlife Acoustics)
├── ls1x.jl           # LS1X (Loggerhead Instruments)
└── snap.jl           # SNAP (NOAA/OSU, retired)
```

Adding a new recorder means adding one file and one `include` line in
`src/EcoAcoustics.jl`. Nothing else changes.

---

## Using built-in recorder profiles

Pass the recorder name to `read_audio`:

```julia
using EcoAcoustics

ad = read_audio("139635MD01_197K_A6M_RH428_20231012_000654Z.flac";
                recorder = "rockhopper")
```

Internally, `read_audio` calls `parse_filename` which looks up
`RECORDER_PROFILES["rockhopper"]` and parses the filename tokens:

```julia
meta = parse_filename("139635MD01_197K_A6M_RH428_20231012_000654Z.flac";
                      recorder = "rockhopper")
# (timestamp=DateTime(2023,10,12,0,6,54), timezone="UTC",
#  lat=missing, lon=missing, site_id="A6M", recorder_id="RH428")
```

The metadata is then merged with calibration and any user overrides to
produce an `Audiodata` struct.

---

## Providing metadata manually

When filenames have no structure, supply metadata directly:

```julia
using Dates, EcoAcoustics

ad = read_audio("mystery_001.wav";
                recorder  = "unknown",
                starttime = DateTime(2021, 5, 1, 12, 0, 0),
                lat       = 38.5,
                lon       = -74.5,
                site_id   = "MAB01")
```

`parse_filename` returns all-`nothing`/`missing` for unknown recorders and
emits a warning. Your provided metadata populates the final `Audiodata`.

---

## The Audiodata struct

`read_audio` returns an `Audiodata` containing the signal and two
sub-structures:

```julia
struct Audiodata
    sig::Vector{Float64}          # raw uncalibrated samples
    fs::Float32                   # sampling rate [Hz]
    starttime::DateTime           # recorder clock time (may be local, not UTC)
    is_calibrated::Bool           # false until apply_calibration! is called
    calibration::Calibration      # stored but not applied at I/O
    metadata::RecordingMetadata   # provenance fields
end
```

`endtime` is a derived helper — not a stored field — so it is always
consistent with the actual signal length:

```julia
endtime(ad)   # starttime + signal_length / fs (nanosecond precision)
nsamples(ad)  # length(ad.sig)
duration(ad)  # endtime(ad) - ad.starttime
```

Provenance fields live in `RecordingMetadata`:

```julia
struct RecordingMetadata
    timezone::Union{String,Nothing}    # IANA name; nothing if unknown
    lat::Union{Float64,Missing}
    lon::Union{Float64,Missing}
    site_id::Union{String,Nothing}
    recorder::String
    recorder_id::Union{String,Nothing}
end
```

Access them via `ad.metadata.site_id`, `ad.metadata.recorder`, etc.

!!! note "Timezone and UTC"
    `starttime` reflects the recorder's clock. When `timezone` is not
    `"UTC"`, the timestamp is local time and has **not** been converted to
    UTC by this package. Convert to UTC using the IANA timezone name before
    computing absolute time differences across deployments.

---

## Calibration and `is_calibrated`

`Audiodata.calibration` stores the calibration object for the recording.
It is **never applied at I/O time**. The `is_calibrated` flag starts as
`false` and is set to `true` only by `apply_calibration!` (v1, in
progress). Metric functions (SPL, PSD, LTSA) check the flag and warn when
computing on uncalibrated data, labelling output as dBFS.

Three calibration types exist:

| Type | Description |
|------|-------------|
| `NoCalibration()` | No calibration available |
| `ScalarCalibration(sens_db)` | Single total sensitivity in dB re 1 V/µPa |
| `TFCalibration(freqs, tf_db)` | Frequency-dependent transfer function |

**Scalar calibration** is returned by `lookup_calibration` when a
`CalibrationProfile` is registered for the recorder. The total sensitivity
is:

```
sens_db = hydrophone_sensitivity + preamp_gain + board_gain
          + 20·log10(1 / Vadc_0pk)
```

**Transfer-function (TF) calibration** applies to recorders like the
Rockhopper, where sensitivity is frequency-dependent. The TF is stored as a
two-column file (frequency in Hz, sensitivity in dB) supplied by the
manufacturer. At calibration lookup time the file is read once and the
arrays are embedded directly in the `TFCalibration` struct, which then
travels with `Audiodata`. No path reference survives into the runtime
object. TF calibration is implemented in v2.

---

## Adding a new recorder

**Checklist:**

1. Create `src/recorders/yourrecorder.jl` following the examples below.
2. Add `include("recorders/yourrecorder.jl")` to `src/EcoAcoustics.jl`
   after the other recorder includes.
3. Add a `@testset` block to `test/recorder_tests/test_filename_formats.jl`
   with at least one representative filename.
4. If calibration values are known, add a test to
   `test/recorder_tests/test_calibration.jl`.

### Tokenized filenames

Most recorders encode metadata as underscore-separated tokens. Suppose
filenames follow:

```
ST01_20250101_120000_ID42.wav
```

Only specify the fields that apply; everything else defaults to `nothing`.

```julia
# src/recorders/soundtrap.jl
# SoundTrap (Ocean Instruments) — example

RECORDER_PROFILES["soundtrap"] = RecorderProfile(
    "SoundTrap";
    split_char        = '_',
    date_token        = 2,               # "20250101"
    time_token        = 3,               # "120000"
    datetime_format   = "yyyymmdd_HHMMSS",
    site_token        = 1,               # "ST01"
    recorder_id_token = 4,               # "ID42"
)

CALIBRATION_PROFILES["soundtrap"] = CalibrationProfile(
    "SoundTrap 300 HF default",
    -177.0,   # hydrophone sensitivity [dB re 1 V/µPa]
    +16.0,    # preamp gain [dB]
      0.0,    # board gain [dB]
      1.0,    # ADC full-scale peak voltage [V]
)
```

Test:

```julia
@testset "SoundTrap filename parsing" begin
    meta = EcoAcoustics.parse_filename("ST01_20250101_120000_ID42.wav";
                                       recorder = "soundtrap")
    @test meta.timestamp == DateTime(2025, 1, 1, 12, 0, 0)
    @test meta.site_id   == "ST01"
    @test meta.recorder_id == "ID42"
end
```

### Custom parser functions (e.g., SNAP)

When token positions alone are insufficient, supply a `parser` function that
receives the filename basename (without extension) and returns the standard
metadata named tuple. See `src/recorders/snap.jl` for a worked example —
SNAP uses a 2-digit year (`YYMMDD`) that requires explicit expansion before
parsing.

The custom `parser` overrides all generic token logic. Only the function
output matters; the token-index fields on the profile are ignored when a
parser is present.

### Recorders with XML sidecars (e.g., DMON)

Some recorders (e.g., DMON, PAMS) store metadata in a per-file XML sidecar
alongside each WAV. The `has_xml` flag on `RecorderProfile` marks this:

```julia
RECORDER_PROFILES["dmon"] = RecorderProfile(
    ...,
    has_xml = true,
    parser  = _parse_dmon,   # reads both the wav path and the XML sidecar
)
```

The custom `parser` function receives the file path, locates the sidecar
(e.g., same basename with `.xml` extension), parses it, and returns the
standard metadata named tuple. All XML-parsing logic is isolated to
`src/recorders/dmon.jl` — nothing else in the package changes. DMON support
is planned for a future release.

---

## No generic guessing

Calling `parse_filename` with an unrecognised recorder name:

```julia
parse_filename("foo_20200101_120000.wav"; recorder="unknownbox")
```

returns all `nothing`/`missing` fields and emits a warning. There is no
generic heuristic fallback. Silent mis-parsing in scientific workflows is
worse than an explicit warning.

You must either define a `RecorderProfile` or supply metadata manually via
`read_audio`.

---

## How `read_audio` assembles the result

1. `parse_filename(path; recorder)` — extract what the filename knows.
2. `lookup_calibration(path, recorder, meta)` — retrieve calibration.
3. Merge in order of precedence (highest wins):
   - User-supplied kwargs (`starttime`, `lat`, `lon`, `site_id`)
   - Filename-parsed metadata
   - Recorder defaults
4. Construct `Audiodata` with `is_calibrated = false`.

This merge order means you can always override a bad filename parse without
modifying the profile.
