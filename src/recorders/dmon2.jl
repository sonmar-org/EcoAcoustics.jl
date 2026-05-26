# DMON2 underwater recorder (WHOI / Woods Hole Group)
#
# Filename convention:
#   DMON2 units write WAV files with an attached XML sidecar that carries
#   timestamp, recorder serial number, and deployment metadata. The WAV
#   filenames themselves have no standardized structure. Filename-based
#   timestamp parsing is therefore NOT implemented here.
#
#   Options under consideration (see CLAUDE.md):
#     A) Parse original DMON2 names + XML sidecar at index-build time.
#     B) Run a one-time rename script (e.g. YYMMDD_HHMMSS_Label.wav)
#        and then treat files like any other recorder.
#   Neither approach is implemented in v1. Timestamps must be supplied
#   explicitly (e.g. via read_audio kwargs or the index builder once a
#   strategy is chosen).
#
# Calibration (DMON2 default configuration):
#   Hydrophone sensitivity : -203.0  dB re 1 V/µPa
#   Preamplifier gain      :  +20.0  dB
#   ADC board gain         :  +13.2  dB
#   ADC full-scale peak    :    1.5  V
#   ─────────────────────────────────────────────────────────────────
#   System sensitivity S   = -203.0 + 20.0 + 13.2 + 20·log10(1/1.5)
#                          ≈ -173.32  dB re 1 V/µPa
#
# WAV header note:
#   DMON2 units systematically write a data-chunk size larger than the
#   actual file content. fix_wav_header = true causes read_audio to use
#   _wavread_corrected instead of wavread for these files.
#
# XML sidecar note:
#   has_xml = true flags that every WAV file may carry a same-stem .xml
#   sidecar. XML parsing is deferred to v2 (see CLAUDE.md).

RECORDER_PROFILES["dmon2"] = RecorderProfile(
    "DMON2";
    fix_wav_header = true,
    has_xml        = true,
    # No filename parser — DMON2 filenames carry no parseable metadata.
    # Timestamps and recorder IDs must come from the XML sidecar (v2)
    # or be supplied explicitly via read_audio kwargs.
)

CALIBRATION_PROFILES["dmon2"] = CalibrationProfile(
    "DMON2 default",
    -203.0,   # hydrophone sensitivity [dB re 1 V/µPa]
     +20.0,   # preamplifier gain [dB]
     +13.2,   # ADC board gain [dB]
       1.5,   # ADC full-scale peak voltage [V]
              # → system sensitivity ≈ -173.32 dB re 1 V/µPa
)
