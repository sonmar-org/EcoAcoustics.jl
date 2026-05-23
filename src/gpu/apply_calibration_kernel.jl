# KernelAbstractions kernel for ScalarCalibration.
#
# This file is the architectural template for GPU work in EcoAcoustics.jl.
# The pattern demonstrated here — @kernel + get_backend dispatch — is the
# one all future GPU primitives (spectrogram, SPL, etc.) will follow.
#
# How it works:
#   1. @kernel defines a device-agnostic kernel. KA compiles it for whichever
#      backend the output array belongs to.
#   2. apply_calibration! inspects the array with get_backend, which returns
#      CPU() for Array{Float64} or CUDABackend() for CuArray, etc.
#   3. The kernel is instantiated with the backend, then launched with ndrange.
#   4. synchronize() blocks until the kernel completes before returning.
#
# Users opt into GPU by passing GPU arrays. No :device kwarg exists; the array
# type IS the device declaration.
#
# GPU requirements (not hard dependencies — KA is the only required dep):
#   NVIDIA: using CUDA   (registers CUDABackend)
#   AMD:    using AMDGPU  (registers ROCBackend)
#   Apple:  using Metal   (registers MetalBackend)

@kernel function _scalar_cal_kernel!(out, @Const(signal), scale)
    i = @index(Global, Linear)
    @inbounds out[i] = signal[i] * scale
end

"""
    apply_calibration!(out, signal, cal::ScalarCalibration)

Purpose:     Apply scalar (frequency-independent) calibration to `signal`
             in-place, multiplying every sample by the linear equivalent of
             `cal.system_sensitivity_dB`. Dispatches to a
             KernelAbstractions kernel whose backend is determined by the array
             type of `out` at runtime.

Arguments:
- `out::AbstractArray`: Output buffer, same shape as `signal`. For CPU execution
  pass `Array{Float64}`; for GPU pass the corresponding GPU array type (e.g.
  `CuArray{Float64}` after `using CUDA`).
- `signal::AbstractArray`: Input in normalised ADC units [−1, 1]. Must have the
  same element type and length as `out`.
- `cal::ScalarCalibration`: Calibration holding `system_sensitivity_dB` in dB
  re 1 V/µPa (conventionally negative). The linear multiplier is
  10^(−system_sensitivity_dB / 20). For SM3M at −153 dB this is ≈ 4.47 × 10⁷.

Returns:     `out`, modified in-place.

Constraints:
- `length(out)` must equal `length(signal)`.
- `out` and `signal` must reside on the same device (both CPU arrays or both
  GPU arrays of the same backend). Mixing backends causes a KA error.
- GPU backends are opt-in: `using CUDA` (or AMDGPU / Metal) must be called
  before passing GPU arrays. KernelAbstractions is the only hard dependency.

Fails when:
- `length(out) ≠ length(signal)`.
- `out` and `signal` are on different devices.

Example:
```julia
# CPU (always works)
apply_calibration!(out, signal, cal)

# GPU (requires using CUDA)
using CUDA
out_gpu    = CuArray{Float64}(undef, length(signal))
signal_gpu = CuArray(signal)
apply_calibration!(out_gpu, signal_gpu, cal)
@assert Array(out_gpu) ≈ Array(out)   # CPU and GPU results identical
```
"""
function apply_calibration!(out::AbstractArray,
                            signal::AbstractArray,
                            cal::ScalarCalibration)
    @assert length(out) == length(signal) begin
        "apply_calibration!: out and signal must have the same length, " *
        "got $(length(out)) and $(length(signal))"
    end
    scale   = 10.0 ^ (-Float64(cal.system_sensitivity_dB) / 20.0)
    backend = KernelAbstractions.get_backend(out)
    kernel! = _scalar_cal_kernel!(backend)
    kernel!(out, signal, scale; ndrange = length(signal))
    KernelAbstractions.synchronize(backend)
    return out
end
