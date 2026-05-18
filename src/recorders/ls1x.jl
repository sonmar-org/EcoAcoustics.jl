# LS1X underwater recorder (Loggerhead Instruments)
#
# Example filename: 20210319T163400_2614231252441225_2.0dB_3.8V_ver2.00.wav
# Tokens (underscore-separated):
#   1: combined datetime ("20210319T163400")
#   2: recorder serial number ("2614231252441225")
#   3: gain setting ("2.0dB")       — not currently parsed
#   4: battery voltage ("3.8V")     — not currently parsed
#   5: firmware version ("ver2.00")
#
# Gain and battery voltage are encoded in the filename but not extracted yet.
# A future version could use these to flag per-file calibration adjustments.
# Timestamps reflect the recorder's local clock; no UTC convention.
#
# Calibration: scalar, not yet registered. Add a CalibrationProfile entry
# here once sensitivity, preamp gain, board gain, and Vadc_0pk are confirmed.

RECORDER_PROFILES["ls1x"] = RecorderProfile(
    "LS1X";
    split_char      = '_',
    datetime_token  = 1,                  # combined date-time token
    datetime_format = "yyyymmddTHHMMSS",
    recorder_id_token = 2,
)
