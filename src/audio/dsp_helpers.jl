# ─── DSP building blocks ──────────────────────────────────────────────────────
#
# Low-level signal-processing helpers shared across calibration, spectrogram,
# and future soundscape primitives. Nothing in this file opens audio files or
# references domain types (Audiodata, Calibration, etc.).

# Purpose:  Construct a symmetric window vector of length N.
#           Used by the spectrogram and any other signal-processing primitive
#           that needs a smooth tapering window. Lives here rather than in
#           calibration.jl so future modules can use it without importing
#           calibration types.
# Constraints: N ≥ 2. All windows peak at n = N÷2 (1-indexed: index N÷2+1).
#              Hann and Blackman are zero at both endpoints.
# Fails when: window is not one of the four recognised symbols.
#
# Sidelobe properties:
#   :hann        first sidelobe −31 dB, rolloff −18 dB/octave
#   :hamming     first sidelobe −43 dB, rolloff  −6 dB/octave
#   :blackman    first sidelobe −58 dB, rolloff −18 dB/octave
#   :rectangular no sidelobe attenuation (use only for very smooth TF data)
function _make_window(window::Symbol, N::Int)::Vector{Float64}
    n = 0:(N - 1)
    if window === :hann
        return @. 0.5 * (1.0 - cos(2π * n / (N - 1)))
    elseif window === :hamming
        return @. 0.54 - 0.46 * cos(2π * n / (N - 1))
    elseif window === :blackman
        return @. 0.42 - 0.5 * cos(2π * n / (N - 1)) + 0.08 * cos(4π * n / (N - 1))
    elseif window === :rectangular
        return ones(Float64, N)
    else
        throw(ArgumentError(
            "_make_window: window must be :hann, :hamming, :blackman, or " *
            ":rectangular, got :$window"))
    end
end
