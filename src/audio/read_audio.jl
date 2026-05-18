using Dates
using FLAC: FLACDecoder
using WAV: wavread

# ─── Internal helpers for RIFF header reading / writing ──────────────────────

# Purpose:     Read a little-endian unsigned 32-bit integer from a byte vector.
# Constraints: `offset` is 1-based. Reads bytes[offset .. offset+3].
#              Caller must ensure offset + 3 <= length(bytes).
# Fails when:  BoundsError if the slice falls outside `bytes`.
function _le_uint32(bytes::Vector{UInt8}, offset::Int)
    UInt32(bytes[offset])           |
    UInt32(bytes[offset + 1]) << 8  |
    UInt32(bytes[offset + 2]) << 16 |
    UInt32(bytes[offset + 3]) << 24
end

# Purpose:     Write a little-endian unsigned 32-bit integer into a byte vector
#              at a 1-based offset. Mutates `bytes` in place.
# Constraints: `offset` is 1-based. Writes to bytes[offset .. offset+3].
#              Caller must ensure offset + 3 <= length(bytes).
# Fails when:  BoundsError if the slice falls outside `bytes`.
function _write_le_uint32!(bytes::Vector{UInt8}, offset::Int, value::UInt32)
    bytes[offset]     =  value        & 0xff
    bytes[offset + 1] = (value >> 8)  & 0xff
    bytes[offset + 2] = (value >> 16) & 0xff
    bytes[offset + 3] = (value >> 24) & 0xff
end

# ─── Internal: WAV reading with in-memory RIFF header correction ─────────────
#
# Purpose:    Read a WAV file whose RIFF data-chunk size field claims more
#             bytes than the file actually contains. Corrects the field
#             in-memory and passes the patched buffer to wavread.
#
# Constraints: Only corrects the case where actual < claimed (truncation).
#              If the file is not recognisable RIFF/WAVE, or if wavread still
#              fails after correction, the original exception is re-thrown.
#
# Fails when:  File cannot be read, is not a RIFF/WAVE file, or wavread fails
#              even after correction.
function _wavread_corrected(path::AbstractString;
                            warn::Bool,
                            original_error::Union{Exception,Nothing} = nothing)

    _rethrow(e) = original_error !== nothing ? throw(original_error) : throw(e)

    bytes = try
        read(path)
    catch e
        _rethrow(e)
    end

    # Verify RIFF/WAVE signature: bytes 1–4 = "RIFF", bytes 9–12 = "WAVE".
    if length(bytes) < 12 ||
       String(bytes[1:4]) != "RIFF" ||
       String(bytes[9:12]) != "WAVE"
        original_error !== nothing ? throw(original_error) :
            error("_wavread_corrected: not a RIFF/WAVE file: \"$path\"")
    end

    # Scan forward through RIFF sub-chunks to find the "data" chunk.
    # Each chunk: 4-byte ID + 4-byte little-endian size + `size` bytes of data.
    # RIFF pads chunks to even byte boundaries, so advance by size + size % 2.
    i = 13  # 1-based index of first sub-chunk (right after the 12-byte header)
    while i + 7 <= length(bytes)
        chunk_id   = String(bytes[i : i + 3])
        chunk_size = Int(_le_uint32(bytes, i + 4))

        if chunk_id == "data"
            # Bytes actually present after the 8-byte chunk header (4 ID + 4 size).
            actual = length(bytes) - (i + 8) + 1

            if actual < chunk_size
                # Truncated file: rewrite the data-chunk size and the top-level
                # RIFF chunk size to match the actual file content.
                _write_le_uint32!(bytes, i + 4, UInt32(actual))
                _write_le_uint32!(bytes, 5,     UInt32(length(bytes) - 8))

                if warn
                    @warn "read_audio: WAV header corrected — file appears " *
                          "truncated. Data chunk claimed $chunk_size bytes " *
                          "but file contains $actual. " *
                          "Likely a battery-dropout file." path=path
                end
            end
            break
        end

        # Advance past this chunk (8-byte header + data + alignment padding).
        i += 8 + chunk_size + chunk_size % 2
    end

    # Feed the (possibly corrected) bytes to wavread via an in-memory buffer.
    # If wavread still fails, re-throw the original error so the caller sees
    # a consistent exception — not a secondary one from the corrected read.
    try
        return wavread(IOBuffer(bytes))
    catch e
        _rethrow(e)
    end
end

# ─── Internal: normalize audio loader outputs into (sig::Vector{Float64}, fs::Float32)
#
# Purpose:     Convert the raw return value from WAV.jl or FLAC.jl into a
#              plain (Vector{Float64}, Float32) pair for the Audiodata constructor.
# Constraints: Input must be a Tuple whose first element is either a
#              Vector{<:Real} (mono) or a Matrix{<:Real} with exactly one
#              column (mono stored as N×1). Any other shape throws immediately.
# Fails when:  Input is not a Tuple, or audio has more than one channel.
#
# Supported shapes:
#   WAV.jl  → (Matrix{Float64}(N,1), fs, nbits, chunks)  — N×1 matrix, Float64
#   FLAC.jl → (Matrix{Float32}(N,1), Int32(fs))          — N×1 matrix, Float32
#             or (Vector{Float32}(N), Int32(fs))          — 1D for truncated reads
function _normalize_loaded_audio(x)
    if !(x isa Tuple)
        error("_normalize_loaded_audio: expected a Tuple from WAV.jl or FLAC.jl, " *
              "got $(typeof(x))")
    end

    data = x[1]   # audio samples: Vector or Matrix, Float32 or Float64
    fs   = x[2]   # sample rate (Int32 from FLAC.jl, Float64 from WAV.jl)

    # FLAC.jl returns an N×1 Matrix{Float32} (normal) or a Vector{Float32} (truncated
    # read, where the stream ended early). WAV.jl returns an N×1 Matrix{Float64}.
    # All three cases reduce to the same mono Vector{Float64} here.
    # `Float64.(...)` converts each element. `vec(data[:, 1])` extracts column 1
    # as a 1D vector (`data[:, 1]` = all rows, first column).
    sig =
        data isa AbstractVector                          ? Float64.(data) :
        (data isa AbstractMatrix && size(data, 2) == 1) ? Float64.(vec(data[:, 1])) :
        error("Only mono audio supported. Got array of size $(size(data)).")

    return sig, Float32(fs)
end


"""
    read_audio(path; starttime, recorder, lat, lon, site_id, strict)

Purpose:     Read a mono audio file and return an `Audiodata` struct with signal,
             sample rate, timestamps, and recorder metadata.

Arguments:
- `path::AbstractString`:
    Path to the audio file. Supported formats: `.wav`, `.flac`.
- `starttime::Union{DateTime,Nothing} = nothing`:
    Override for recording start time. When provided, takes precedence over any
    timestamp parsed from the filename — including `DateTime(0)` if you explicitly
    want year zero. When `nothing` (default), the filename is parsed first; if that
    also fails, falls back to `DateTime(0)` with a warning (if `strict=true`).
- `recorder::AbstractString = "unknown"`:
    Recorder family name (e.g. `"rockhopper"`, `"sm3m"`, `"ls1x"`, `"snap"`).
    Used to select filename parsing rules and calibration profile.
- `lat::Union{Float64,Missing} = missing`:
    Recorder latitude in decimal degrees (positive = North). Overrides any value
    parsed from the filename. Leave as `missing` when location is unknown.
- `lon::Union{Float64,Missing} = missing`:
    Recorder longitude in decimal degrees (positive = East). Overrides any value
    parsed from the filename. Leave as `missing` when location is unknown.
- `site_id::Union{String,Nothing} = nothing`:
    Deployment site identifier (e.g. `"T1-C"`). Override for filename-parsed value.
- `strict::Bool = false`:
    If `true`, emit warnings when calibration or metadata fields are missing.

Returns:     `Audiodata` with signal in raw linear amplitude (uncalibrated), sample
             rate as `Float32`, timestamps as `DateTime` (UTC where known), and
             calibration stored but not applied.

Constraints: Only mono audio is supported. Multi-channel files will throw an error.
             Keyword arguments act as overrides: filename-parsed values are used when
             a keyword argument is not supplied (i.e. still at its default).

Fails when:  File cannot be read or is corrupt; audio is not mono; signal is
             empty after decoding; file extension is not `.wav` or `.flac`
             (throws `ArgumentError` — convert AIF with `sox input.aif output.flac`).

Example:
```julia
a = read_audio("T1-C__0__20170912_181500.wav"; recorder="sm3m")
a = read_audio("data.wav"; starttime=DateTime(2021,3,19,16,34,0), recorder="ls1x")
```

Do not use when: Audio has multiple channels (stereo/hydrophone arrays) — mono is
                 enforced.
"""
function read_audio(path::AbstractString;
                    starttime::Union{DateTime,Nothing} = nothing,
                    recorder::AbstractString = "unknown",
                    lat::Union{Float64,Missing} = missing,
                    lon::Union{Float64,Missing} = missing,
                    site_id::Union{String,Nothing} = nothing,
                    strict::Bool = false)

    ext = lowercase(splitext(path)[2])

    # Parse metadata from filename, then look up calibration.
    filename_meta = parse_filename(path; recorder=recorder, strict=strict)
    cal           = lookup_calibration(path, recorder, filename_meta; strict=strict)

    # Merge starttime: caller override takes precedence; fall back to filename; warn if neither.
    merged_starttime =
        starttime !== nothing               ? starttime :
        filename_meta.timestamp !== nothing ? filename_meta.timestamp :
        (strict && @warn("read_audio: no timestamp found in filename or arguments." *
                         " Defaulting to DateTime(0).", path=path); DateTime(0))

    # Merge site_id: keyword argument overrides filename.
    merged_site_id =
        site_id !== nothing             ? site_id :
        filename_meta.site_id

    # Merge timezone: only source is the filename parser for now.
    merged_timezone = filename_meta.timezone

    # Dispatch on file extension. In Julia, if/elseif/else is an expression —
    # `raw` receives whatever value the matching branch returns (the audio data
    # as a tuple), rather than executing the branch for side-effects only.
    #
    # .wav  — WAV.jl (wavread). Recorders with fix_wav_header=true (e.g. DMON2)
    #         always use the in-memory RIFF correction path silently; all others
    #         try wavread directly and fall back to correction with a warning on
    #         any exception (battery-dropout safety net).
    #
    # .flac — FLAC.jl direct (FLACDecoder). FileIO/LibSndFile is intentionally
    #         bypassed; see CLAUDE.md Key Constraints. Returns (Matrix{Float32},
    #         Int32) which _normalize_loaded_audio handles via the Tuple branch.
    #
    # other — unsupported in v1. AIF files should be converted with:
    #         sox input.aif output.flac
    profile = get(RECORDER_PROFILES, String(recorder), nothing)
    raw = if ext == ".wav"
        if profile !== nothing && profile.fix_wav_header
            _wavread_corrected(path; warn=false)
        else
            try
                wavread(path)
            catch e
                _wavread_corrected(path; warn=true, original_error=e)
            end
        end

    elseif ext == ".flac"
        # `local f` is required because Julia scopes variables to the block in
        # which they are first assigned. Without it, `f` would be invisible outside
        # the try block, and the lines below the catch could not use it.
        local f
        try
            f = FLACDecoder(String(path))  # FLACDecoder requires String, not AbstractString
        catch e
            # sprint(showerror, e) converts the exception object to a human-readable
            # string so it can be embedded in our own error message.
            throw(ErrorException(
                "read_audio: cannot read FLAC file \"$path\": $(sprint(showerror, e))"
            ))
        end
        # channels == 0 means no STREAMINFO block was parsed — empty or wholly corrupt file.
        if f.metadata.channels == 0
            error("read_audio: \"$path\" is not a valid FLAC file or contains no audio " *
                  "(STREAMINFO metadata missing)")
        end
        # read and length here dispatch on FLACDecoder — FLAC.jl extends the standard
        # Base.read and Base.length to work on its decoder type. Returns
        # Array{Float32,2}(nsamples, channels); may be fewer rows than length(f)
        # if the file is truncated (battery dropout).
        data = read(f, length(f))
        (data, Int32(f.metadata.samplerate))

    else
        throw(ArgumentError(
            "read_audio: unsupported audio format \"$ext\" in \"$path\". " *
            "EcoAcoustics.jl v1 supports .wav and .flac only. " *
            "Convert AIF files with: sox input.aif output.flac"
        ))
    end

    sig, fs = _normalize_loaded_audio(raw)

    return Audiodata(sig, fs, merged_starttime;
                     timezone    = merged_timezone,
                     lat         = lat,
                     lon         = lon,
                     site_id     = merged_site_id,
                     recorder    = String(recorder),
                     recorder_id = filename_meta.recorder_id,
                     calibration = cal)
end

