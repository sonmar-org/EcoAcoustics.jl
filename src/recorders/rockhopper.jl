# Rockhopper underwater recorder (Embedded Ocean Systems)
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
# Calibration: frequency-dependent (TF), supplied per unit by the manufacturer.
# No single scalar applies across all units. TF calibration is v2; until then
# lookup_calibration returns NoCalibration for Rockhopper.

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
