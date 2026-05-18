```@meta
CurrentModule = EcoAcoustics
```

# EcoAcoustics.jl

EcoAcoustics.jl is a Julia package for passive acoustic monitoring analysis.
It is designed for recordings from underwater hydrophone deployments — its
initial target is cetacean acoustics in the Mid-Atlantic Bight — but the
architecture is general-purpose and not species-specific.

**Core workflow:**

1. **Read audio** — `read_audio` loads a WAV or FLAC file and returns an
   `Audiodata` struct containing the signal, sample rate, timestamps, location,
   and calibration parameters. Filename metadata (timestamps, recorder identity,
   deployment site) is parsed automatically for supported recorder families.

2. **Build an index** — for archive-scale work, `build_index` scans a directory
   of audio files once and writes a lightweight Arrow index (one row per file).
   `IndexedFileSource` then loads audio by time window on demand, reading only
   the files that overlap the requested range.

3. **Apply calibration** — `apply_calibration!` converts raw ADC counts to
   physical units (dB re 1 µPa). Calibration is stored with the recording but
   never applied at I/O time, so raw and calibrated data cannot be confused.

4. **Compute metrics** — soundscape indices (SPL, PSD, ACI, entropy, LTSA) are
   computed via chunked streaming over any `AbstractAudioSource`. Implementations
   are in progress for v1.

**Supported recorders:** Rockhopper (Embedded Ocean Systems), SM3M (Wildlife
Acoustics), LS1X (Loggerhead Instruments), SNAP (NOAA/OSU). Adding a new
recorder requires one file and one line in `src/EcoAcoustics.jl`.

**Supported formats:** WAV and FLAC. Files in other formats (AIF, OGG) should
be converted to FLAC with `sox` before indexing.

---

## API Reference

```@index
```

```@autodocs
Modules = [EcoAcoustics]
Private = false
```
