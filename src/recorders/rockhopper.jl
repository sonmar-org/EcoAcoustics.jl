# Rockhopper underwater recorder (Cornell University)
#
# Example filename: 139635MD01_197K_A6M_RH428_20231012_000654Z.flac
# Tokens (underscore-separated):
#   1: recorder serial prefix ("139635MD01")
#   2: unknown field ("197K")
#   3: site_id ("A6M")
#   4: recorder_id ("RH428")
#   5: date ("20231012")
#   6: time with UTC suffix ("000654Z")  — trailing Z confirms UTC
#
# Calibration: frequency-dependent (TF), supplied per unit by Cornell.
# The TF is per-unit; the CSV at tf_path encodes one representative curve.
# lookup_calibration reads the CSV and constructs a TFCalibration object.
# CSV values are positive (Cornell's source convention); _load_tf_csv negates
# on read so TFCalibration.tf_db is negative (package convention).

RECORDER_PROFILES["rockhopper"] = RecorderProfile(
    "Rockhopper";
    split_char      = '_',
    date_token      = 5,
    time_token      = 6,
    datetime_format = "yyyymmdd_HHMMSS",    # Z suffix is stripped before parsing
    site_token      = 3,
    recorder_id_token = 4,
    timezone_source = :suffix_Z,            # trailing Z on time token means UTC
)

CALIBRATION_PROFILES["rockhopper"] = CalibrationProfile(
    "Rockhopper (Cornell University, per-unit TF)",
    0.0, 0.0, 0.0, 1.0;                    # gains embedded in the TF; scalar fields unused
    tf_path = joinpath(@__DIR__, "rockhopper_TF.csv"),
)
