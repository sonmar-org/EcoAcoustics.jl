---
name: Calibration decisions
description: Rockhopper TF policy and serial number handling
type: project
---

Rockhopper calibration uses a single TF for all units (only one TF observed across all deployed Rockhoppers).

Serial number (recorder_id) stays in Audiodata metadata for traceability but is NOT used for calibration lookup in v1. All Rockhoppers share the same CalibrationProfile.

TFCalibration struct is already defined with freqs::Vector{Float32} and tf_db::Vector{Float32}. Rockhopper profile must populate this from a TF file (path in CalibrationProfile.tf_path).

**Why:** No per-unit variation observed in field deployments. May need revisiting if new units show deviation.

**How to apply:** Do not add per-serial calibration lookup logic for Rockhopper in v1. Keep tf_path in CalibrationProfile for loading the TF file.
