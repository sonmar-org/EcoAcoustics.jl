# SNAP recorder (NOAA/OSU — retired; data archives still in use)
#
# Files are renamed post-deployment as: SITE_YYMMDD_HHMMSS.wav
# Example: T3C_180630_235000.wav
#
# The YYMMDD date uses a 2-digit year. We assume all SNAP data are from
# 2000–2099 and expand YY → 20YY. If you have pre-2000 recordings, this
# parser must be revised.
#
# No recorder_id in the filename; unit identity is tracked via deployment logs.
#
# Calibration: not registered. SNAP units were deployed with varied hydrophones
# and gain settings; provide calibration per-deployment via read_audio kwargs.

function _parse_snap_filename(basename_noext::AbstractString)
    tokens = split(basename_noext, '_')

    # Need at least SITE, YYMMDD, HHMMSS
    if length(tokens) < 3
        return (timestamp=nothing, timezone=nothing,
                lat=missing, lon=missing,
                site_id=nothing, recorder_id=nothing)
    end

    site_id = String(tokens[1])
    date6   = tokens[2]   # "YYMMDD"
    time6   = tokens[3]   # "HHMMSS"

    # Reject if date token isn't exactly 6 characters
    if length(date6) != 6
        return (timestamp=nothing, timezone=nothing,
                lat=missing, lon=missing,
                site_id=site_id, recorder_id=nothing)
    end

    # Expand 2-digit year: "180630" → "20180630_235000"
    dtstr = "20" * date6 * "_" * time6
    ts = try
        DateTime(dtstr, DateFormat("yyyymmdd_HHMMSS"))
    catch
        nothing
    end

    return (timestamp=ts, timezone=nothing,
            lat=missing, lon=missing,
            site_id=site_id, recorder_id=nothing)
end

RECORDER_PROFILES["snap"] = RecorderProfile(
    "SNAP (renamed)";
    split_char = '_',
    parser     = _parse_snap_filename,   # custom parser handles the 2-digit year
)
