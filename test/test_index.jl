using Test
using EcoAcoustics
using Dates
using Logging
import Arrow
import WAV

# Helper: write a small synthetic mono WAV to `path`.
# Uses Float64 samples so the round-trip through wavread is exact.
function _write_test_wav(path; fs = 1000, n = 1000, value = 0.5)
    WAV.wavwrite(fill(value, n), path; Fs = fs)
end

# ─── build_index ─────────────────────────────────────────────────────────────

@testset "build_index: bad directory throws" begin
    @test_throws ArgumentError EcoAcoustics.build_index("/no/such/path"; recorder = "sm3m")
end

@testset "build_index: schema and row count" begin
    mktempdir() do dir
        _write_test_wav(joinpath(dir, "T1-C__0__20230101_120000.wav"))
        tbl = EcoAcoustics.build_index(dir; recorder = "sm3m")
        @test tbl isa Arrow.Table
        @test length(tbl.start_time) == 1
        # All eleven schema columns must be present.
        @test Set(propertynames(tbl)) == EcoAcoustics._EXPECTED_INDEX_COLUMNS
    end
end

@testset "build_index: skips non-audio files" begin
    mktempdir() do dir
        _write_test_wav(joinpath(dir, "T1-C__0__20230101_120000.wav"))
        write(joinpath(dir, "notes.txt"), "ignore me")
        tbl = EcoAcoustics.build_index(dir; recorder = "sm3m")
        @test length(tbl.start_time) == 1
    end
end

@testset "build_index: recursive directory walk" begin
    mktempdir() do dir
        subdir = mkdir(joinpath(dir, "sub"))
        _write_test_wav(joinpath(dir,    "T1-C__0__20230101_120000.wav"))
        _write_test_wav(joinpath(subdir, "T1-C__0__20230101_120010.wav"))
        tbl = EcoAcoustics.build_index(dir; recorder = "sm3m")
        @test length(tbl.start_time) == 2
    end
end

@testset "build_index: sorted ascending by start_time" begin
    mktempdir() do dir
        # Write the later-timestamped file first to verify the sort is not
        # filesystem-order-dependent.
        _write_test_wav(joinpath(dir, "T1-C__0__20230101_120010.wav"))
        _write_test_wav(joinpath(dir, "T1-C__0__20230101_120000.wav"))
        tbl = EcoAcoustics.build_index(dir; recorder = "sm3m")
        @test length(tbl.start_time) == 2
        @test tbl.start_time[1] < tbl.start_time[2]
        @test tbl.start_time[1] == DateTime(2023, 1, 1, 12, 0,  0)
        @test tbl.start_time[2] == DateTime(2023, 1, 1, 12, 0, 10)
    end
end

@testset "build_index: relative file paths" begin
    mktempdir() do dir
        _write_test_wav(joinpath(dir, "T1-C__0__20230101_120000.wav"))
        tbl = EcoAcoustics.build_index(dir; recorder = "sm3m")
        # Path is relative to the indexed directory, not absolute.
        @test tbl.file_path[1] == "T1-C__0__20230101_120000.wav"
    end
end

@testset "build_index: nsamples is actual decoded count" begin
    mktempdir() do dir
        # Use an odd count so coincidental alignment cannot mask an error.
        _write_test_wav(joinpath(dir, "T1-C__0__20230101_120000.wav"); n = 1234)
        tbl = EcoAcoustics.build_index(dir; recorder = "sm3m")
        @test tbl.nsamples[1] == 1234
    end
end

@testset "build_index: end_time matches nsamples and fs" begin
    mktempdir() do dir
        _write_test_wav(joinpath(dir, "T1-C__0__20230101_120000.wav"); n = 1000, fs = 1000)
        tbl = EcoAcoustics.build_index(dir; recorder = "sm3m")
        expected_end = tbl.start_time[1] +
                       Nanosecond(round(Int, 1e9 * tbl.nsamples[1] / tbl.fs[1]))
        @test tbl.end_time[1] == expected_end
    end
end

@testset "build_index: output_path and round-trip" begin
    mktempdir() do dir
        _write_test_wav(joinpath(dir, "T1-C__0__20230101_120000.wav"))
        arrow_path = joinpath(dir, "index.arrow")
        tbl = EcoAcoustics.build_index(dir; recorder = "sm3m", output_path = arrow_path)
        @test isfile(arrow_path)

        tbl2 = EcoAcoustics.load_index(arrow_path)
        @test tbl2 isa Arrow.Table
        @test length(tbl2.start_time) == 1
        @test tbl2.start_time[1]  == tbl.start_time[1]
        @test tbl2.nsamples[1]    == tbl.nsamples[1]
        @test tbl2.file_path[1]   == tbl.file_path[1]
    end
end

@testset "build_index: unreadable file is skipped" begin
    # A file with a .wav extension but garbage content cannot be parsed.
    # build_index logs a warning and continues; only the valid file is indexed.
    # Two warnings fire in order: no-timestamp from parse_filename on corrupt.wav,
    # then skipping from build_index's catch block.
    mktempdir() do dir
        _write_test_wav(joinpath(dir, "T1-C__0__20230101_120000.wav"))
        write(joinpath(dir, "corrupt.wav"), "not a wav file")
        @test_logs (:warn, r"no timestamp") (:warn, r"skipping unreadable") min_level=Logging.Warn begin
            tbl = EcoAcoustics.build_index(dir; recorder = "sm3m")
            @test length(tbl.start_time) == 1
        end
    end
end

@testset "build_index: WAV header correction" begin
    # Write a valid WAV, then inflate its data-chunk size field so that the
    # header claims more bytes than the file actually contains. This simulates
    # a truncated file (e.g. battery dropout mid-write). build_index must still
    # index the file and report the actual, not the claimed, sample count.
    mktempdir() do dir
        path = joinpath(dir, "T1-C__0__20230101_120000.wav")
        _write_test_wav(path; n = 1000)

        # Scan the RIFF sub-chunks to find "data" and inflate its size field.
        bytes = read(path)
        i = 13  # first sub-chunk starts right after the 12-byte RIFF header
        while i + 7 <= length(bytes)
            if String(bytes[i : i + 3]) == "data"
                current = EcoAcoustics._le_uint32(bytes, i + 4)
                EcoAcoustics._write_le_uint32!(bytes, i + 4, current + UInt32(10_000))
                break
            end
            sz = Int(EcoAcoustics._le_uint32(bytes, i + 4))
            i += 8 + sz + sz % 2
        end
        write(path, bytes)

        # build_index should recover via _wavread_corrected and record the
        # actual sample count, not the inflated header claim. The truncated-file
        # correction emits a warning, which is the expected behavior under test.
        @test_logs (:warn, r"WAV header corrected") min_level=Logging.Warn begin
            tbl = EcoAcoustics.build_index(dir; recorder = "sm3m")
            @test length(tbl.start_time) == 1
            @test tbl.nsamples[1] == 1000
        end
    end
end

# ─── load_index ──────────────────────────────────────────────────────────────

@testset "load_index: missing columns throws" begin
    mktempdir() do dir
        bad_path = joinpath(dir, "bad.arrow")
        Arrow.write(bad_path, (x = [1, 2, 3],))
        @test_throws ErrorException EcoAcoustics.load_index(bad_path)
    end
end
