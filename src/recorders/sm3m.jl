# SM3M underwater recorder (Wildlife Acoustics)
#
# Example filename: T1-C__0__20170912_181500.wav
# Tokens (underscore-separated):
#   1: site_id ("T1-C")
#   2: empty (double-underscore separator)
#   3: channel index ("0")
#   4: empty
#   5: date ("20170912")
#   6: time ("181500")
#
# No recorder_id in the filename; unit identity is tracked via deployment logs.
# Timestamps reflect the recorder's local clock; no UTC convention.
#
# Calibration: scalar. Values below are for the SM3M internal hydrophone at
# standard gain. Verify against the deployment configuration for your unit.
# Total sensitivity: -165 + 12 + 0 + 20·log10(1/1.0) = -153 dB re 1 V/µPa.

RECORDER_PROFILES["sm3m"] = RecorderProfile(
    "SM3M";
    split_char      = '_',
    date_token      = 5,
    time_token      = 6,
    datetime_format = "yyyymmdd_HHMMSS",
    site_token      = 1,
)

CALIBRATION_PROFILES["sm3m"] = CalibrationProfile(
    "SM3M default",
    -165.0,   # hydrophone sensitivity [dB re 1 V/µPa]
    +12.0,    # preamp gain [dB]
     0.0,     # board gain [dB]
     1.0,     # ADC full-scale peak voltage [V]
)
