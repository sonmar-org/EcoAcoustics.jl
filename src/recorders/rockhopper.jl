# Rockhopper underwater recorder (Cornell University)
#
# Example filename: 139635MD01_197K_A6M_RH428_20231012_000654Z.flac
# Tokens (underscore-separated):
#   1: deployment prefix ("139635MD01")
#   2: unknown field ("197K")
#   3: site_id ("A6M")
#   4: recorder_id ("RH428")
#   5: date ("20231012")
#   6: time with UTC suffix ("000654Z")  — trailing Z confirms UTC
#
# Calibration: frequency-dependent (TF). The package ships the Cornell fleet-level
# calibration derived from the RH_calCurves sheet of Rockhopper_TF.xlsx, columns
# 1–5 (AnalogSensitivity). Loaded once at module init into _ROCKHOPPER_TF.
# See src/recorders/calibration_data/README.md and docs/src/calibration/rockhopper.md.

RECORDER_PROFILES["rockhopper"] = RecorderProfile(
    "Rockhopper";
    split_char        = '_',
    date_token        = 5,
    time_token        = 6,
    datetime_format   = "yyyymmdd_HHMMSS",    # Z suffix stripped before parsing
    site_token        = 3,
    recorder_id_token = 4,
    timezone_source   = :suffix_Z,            # trailing Z on time token means UTC
)

# ─── load_tf_calcurves ────────────────────────────────────────────────────────

"""
    load_tf_calcurves(path; vmax_peak_V) -> TFCalibration

Purpose:     Load a Cornell-format Rockhopper calibration CSV and return a
             `TFCalibration` in the package canonical form (dB re full-scale per
             µPa, DD-12). The primary caller is this file's module-init block;
             the function is exported so advanced users can load a custom or
             updated calibration file.

Arguments:
- `path::AbstractString`: Path to a CSV file with a named-column header row.
  Required columns: `Frequency_Hz` and `AnalogSensitivity_dB_re_1VperRefPress`.
  Optional columns (`SensorSensitivity_dB_re_1V_perRefPress`, `PreampGain_dB`,
  `RecorderGain_dB`) are present in the Cornell file for provenance but are not
  used arithmetically by this loader. Additional unknown columns are silently
  ignored. Column order does not matter.
- `vmax_peak_V::Float64`: ADC full-scale peak voltage in volts. `5.0` for the
  Rockhopper. Used to convert AnalogSensitivity to the canonical form:
  `tf_dB = AnalogSensitivity_dB − 20·log10(vmax_peak_V)`.

Returns:     [`TFCalibration`](@ref) with:
             - `frequency`: strictly ascending Hz values from the `Frequency_Hz` column.
             - `tf_dB`: canonical sensitivity in dB re full-scale per µPa.
             - `source`: the absolute path supplied.
             - `format`: `:rockhopper_calcurves_csv`.
             - `conversion_notes`: records the exact dB shift applied.

Constraints:
- `Frequency_Hz` must be strictly ascending (no duplicates).
- `vmax_peak_V` must be > 0.
- The resulting `TFCalibration` is validated by the inner constructor (≥ 2 points,
  strictly ascending, matching lengths).

Fails when:
- `path` is not found or unreadable.
- A required column is missing (error names the missing column and lists required ones).
- `Frequency_Hz` is not strictly ascending.
- Any row has the wrong number of fields.

**Legacy file formats not supported — read before using a custom file:**

1. **Cornell `TF` column (column 13 of RH_calCurves).** The Excel sheet also
   contains a column called `TF` that folds in a PeakToRMS factor of
   20·log10(2√2) ≈ 9.03 dB. This factor is appropriate for converting a
   peak-calibrated sine to RMS but is *not* appropriate for PSD, which operates
   on instantaneous squared samples. A CSV built from the `TF` column rather than
   the `AnalogSensitivity` column (column 5) will produce PSD values that are
   9.03 dB too high. This loader cannot detect the error from column names alone
   if the CSV was constructed incorrectly; verify you have extracted column 5.

2. **Bare 2-column Raven Expedition CSV (no header row, columns: frequency, dB re
   full-scale per µPa).** Those values are already in canonical form; passing such
   a file to *this* loader would subtract 20·log10(vmax_peak_V) ≈ 13.98 dB a
   second time, producing a calibration that is 13.98 dB too negative and PSD
   values that are 13.98 dB too low. Use `_load_tf_csv` (private, headerless-CSV
   reader) and construct `TFCalibration` manually if you have a bare 2-column file.

Example:
```julia
# Package-internal use (module init):
tf = load_tf_calcurves(
    joinpath(@__DIR__, "calibration_data", "rockhopper_tf_calibration.csv");
    vmax_peak_V = 5.0)

# Advanced user with a custom calibration file:
my_tf = load_tf_calcurves("/path/to/my_RH_calCurves_extract.csv"; vmax_peak_V = 5.0)
psd   = compute_psd(audio; cal = my_tf)
```

Do not use when:
- You have a bare 2-column CSV (no header). See the legacy warning above.
- You have a CSV derived from the Cornell `TF` column (column 13). See above.
"""
function load_tf_calcurves(path::AbstractString; vmax_peak_V::Float64)::TFCalibration
    vmax_peak_V > 0 || throw(ArgumentError(
        "load_tf_calcurves: vmax_peak_V must be > 0, got $(vmax_peak_V)"))

    # Read and filter blank lines. readlines returns one String per line.
    all_lines = filter(!isempty, strip.(readlines(path)))
    isempty(all_lines) && throw(ArgumentError(
        "load_tf_calcurves: file has no content: \"$(basename(path))\""))

    # First line is the header. strip.(split(..., ',')) trims whitespace from
    # each field name so column matching is not whitespace-sensitive.
    headers = String.(strip.(split(all_lines[1], ',')))

    # Verify required columns are present. Check both before throwing so the
    # error message is complete (not just the first missing one).
    required = ("Frequency_Hz", "AnalogSensitivity_dB_re_1VperRefPress")
    missing_cols = filter(c -> c ∉ headers, required)
    if !isempty(missing_cols)
        throw(ArgumentError(
            "load_tf_calcurves: required column(s) not found in " *
            "\"$(basename(path))\": $(join(missing_cols, ", ")). " *
            "Required: $(join(required, ", ")). " *
            "Found: $(join(headers, ", "))."))
    end

    freq_idx   = findfirst(==("Frequency_Hz"),                          headers)
    analog_idx = findfirst(==("AnalogSensitivity_dB_re_1VperRefPress"), headers)

    # Parse data rows. Each line must have the same number of comma-separated
    # fields as the header; parse only the two columns we need.
    frequency = Float64[]
    analog_dB = Float64[]
    for (i, line) in enumerate(all_lines[2:end])
        parts = String.(strip.(split(line, ',')))
        length(parts) == length(headers) || throw(ArgumentError(
            "load_tf_calcurves: data row $(i) has $(length(parts)) fields but " *
            "header has $(length(headers)) columns in \"$(basename(path))\""))
        push!(frequency, parse(Float64, parts[freq_idx]))
        push!(analog_dB, parse(Float64, parts[analog_idx]))
    end

    isempty(frequency) && throw(ArgumentError(
        "load_tf_calcurves: no data rows in \"$(basename(path))\""))

    # Monotonicity is also enforced by TFCalibration's inner constructor, but
    # we check here to produce a file-specific error message.
    all(diff(frequency) .> 0) || throw(ArgumentError(
        "load_tf_calcurves: Frequency_Hz column must be strictly ascending " *
        "(no duplicates). Found non-ascending or repeated values in " *
        "\"$(basename(path))\"."))

    # Convert AnalogSensitivity to canonical form (dB re full-scale per µPa).
    # AnalogSensitivity is in dB re 1V per reference pressure at the ADC input.
    # Subtracting 20·log10(vmax_peak_V) normalises to full-scale (DD-12).
    shift_dB = 20.0 * log10(vmax_peak_V)
    tf_dB    = analog_dB .- shift_dB

    shift_str = string(round(shift_dB; digits = 4))
    notes = "AnalogSensitivity_dB_re_1VperRefPress reduced by " *
            "20·log10($(vmax_peak_V)) = $(shift_str) dB to convert to " *
            "dB re full-scale per µPa (canonical form, DD-12)."

    return TFCalibration(frequency, tf_dB, path, :rockhopper_calcurves_csv, notes)
end

# ─── RockhopperProfile singleton ─────────────────────────────────────────────

"""
    RockhopperProfile <: AbstractRecorderProfile

Calibration profile for the Cornell Rockhopper underwater recorder.

The shipped calibration (derived from the Cornell RH_calCurves AnalogSensitivity
column) is loaded once at module init and bound to `_RH_PROFILE`. Every call to
`RockhopperProfile()` returns the same singleton object.

Fields
------
- `Vmax_peak_V::Float64`: ADC full-scale peak voltage. `5.0` V for the Rockhopper.
- `tf::TFCalibration`: The shipped transfer-function calibration in canonical form.
"""
struct RockhopperProfile <: AbstractRecorderProfile
    Vmax_peak_V::Float64
    tf::TFCalibration
end

# Load the shipped calibration CSV once at module init. @__DIR__ expands to the
# directory of this source file at compile time, so the path is always relative
# to the package installation, not the working directory.
const _ROCKHOPPER_TF = let
    path = joinpath(@__DIR__, "calibration_data", "rockhopper_tf_calibration.csv")
    load_tf_calcurves(path; vmax_peak_V = 5.0)
end

const _RH_PROFILE = RockhopperProfile(5.0, _ROCKHOPPER_TF)

# Zero-argument constructor returns the pre-built singleton. Every call is
# identical: === identity is guaranteed because _RH_PROFILE is a const.
RockhopperProfile() = _RH_PROFILE

# Val dispatch: get_profile(:rockhopper) → RockhopperProfile().
# The symbol-to-Val bridge (get_profile(id::Symbol)) and the fallback error
# are defined in recorders.jl.
get_profile(::Val{:rockhopper}) = RockhopperProfile()
