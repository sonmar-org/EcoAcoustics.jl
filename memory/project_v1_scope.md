---
name: v1 scope
description: Full feature list for EcoAcoustics.jl v1.0
type: project
---

v1.0 feature list (as of 2026-05-15):

**Audio I/O**
- Read wav and flac files
- Parse dates per recorder profile
- Proper metadata architecture
- Zarr format support (or alternative decision still pending)
- Read from directories of files
- Handle chunking for HPC resources

**Calibration**
- Rockhopper, SM3M, LS1X

**Power calculations**
- With timestamp tracking

**Spectral metrics**
- SPL (including TOLs, arbitrary bands) and PSD
- Spectral averages for LTSA
- Entropy (temporal, permutation)
- ACI (Acoustic Complexity Index)

**Plotting**
- SPL, PSD, metrics, LTSA

**Why:** Dissertation research tool for passive acoustic monitoring of bottlenose dolphin in Mid-Atlantic Bight. Also designed for general-purpose use by others.

**How to apply:** Scope decisions should stay within this list. Do not add features outside it without discussion.
