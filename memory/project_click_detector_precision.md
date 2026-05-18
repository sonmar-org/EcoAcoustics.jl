---
name: project-click-detector-precision
description: Click detector must anchor timestamps to file_start, not chunk boundary — read_audio_range boundaries are millisecond-precise
metadata:
  type: project
---

When the click detector is implemented (v2), timestamps must be stored as
`(file_start::DateTime, sample_offset::Int64)`, not as a bare `DateTime`.

**Why:** `read_audio_range` chunk boundaries are millisecond-precise (Julia
`DateTime` resolution). At 192 kHz, 0.5 ms of rounding is up to 96 samples.
If a click timestamp is computed relative to the chunk boundary rather than
the file anchor, that error propagates into detection accuracy.

**How to apply:** The click detector receives a chunk from `read_audio_range`
whose `starttime` reflects the requested window start (millisecond-precise).
The detector finds click sample index `k` within the chunk, then records:

    file_start    = source.audio.starttime   # the file's own anchor
    sample_offset = src_offset_of_chunk + k  # exact integer, from file start

Wall-clock time is recovered as:
    event_time = file_start + Nanosecond(round(Int, 1e9 * sample_offset / fs))

This is consistent with the CLAUDE.md architecture spec and with how
`endtime(a::Audiodata)` already works.

Cross-reference: CLAUDE.md "Audio Data Model" section, `event_time` helper.
