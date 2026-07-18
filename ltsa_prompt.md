EcoAcoustics.jl LTSA implementation.

Working directory
~/dissertation/EcoAcoustics.jl/. Confirm with `pwd` before starting.

Session goal
Add Long-Term Spectral Average (LTSA) computation to EcoAcoustics.jl.
LTSA shows energy distribution across frequency over long timescales
(hours to days or even years) by averaging PSDs within configurable column spans
(typically 1 minute, 5 minutes, or 1 hour each). First v1.0.0 roadmap
item after the foundation layer.

Required reading before starting
- CLAUDE.md
- ROADMAP.org (top of v1.0.0 section confirms LTSA is the active item)
- docs/design_decisions.md (DD-01 through DD-24)
- src/soundscape/psd.jl — the PSD primitive being composed
- src/soundscape/spl.jl — for the pattern of composing PSDs into
  higher-level results (band-integrated SPL is structurally similar
  to time-averaged LTSA)
- test/test_psd_manta_validation.jl — for the validation pattern

Architecture decisions already settled
- Column averaging: energetic mean only (consistent with average_psd
  primitive and Merchant convention). No dB-mean or median modes in
  this initial implementation; can be added later if needed.
- Column boundaries: non-overlapping consecutive spans. Column 1
  covers 0..span seconds, column 2 covers span..2*span, etc. Final
  partial column dropped if audio doesn't divide evenly.
- Inner FFT: composes compute_psd with its default keyword behavior.
  fft_window and fft_overlap pass through to compute_psd.
- Calibration: same propagation pattern as compute_psd. is_calibrated
  flag, cal provenance, units symbol dispatched on calibration state.

Deliverables — pause for review between each

Deliverable 1 — LTSAResult struct and compute_ltsa function
- Add src/soundscape/ltsa.jl
- LTSAResult struct: matrix::Matrix{Float64} (freq × time, linear
  PSD units), freqs::Vector{Float64}, column_times::Vector{Float64}
  (start time of each column in seconds), fs::Float32,
  average_span_seconds::Float64, is_calibrated::Bool,
  cal::AbstractCalibration
- compute_ltsa(audio::Audiodata; average_span_seconds,
  fft_window_seconds=1.0, fft_overlap=0.5, cal=nothing) → LTSAResult
- Internal loop: walk audio in average_span_seconds chunks, call
  compute_psd on each chunk, take energetic mean across the column's
  FFT frames using existing average_psd primitive, store resulting
  PSD as one column of the matrix
- Add to_dB method for LTSAResult (Matrix{Float64} freq×time in dB)
- Add ltsa_units(result) dispatch matching psd_units pattern
- Export compute_ltsa, LTSAResult, ltsa_units, to_dB (if not already)
- Include in src/EcoAcoustics.jl

Full docstring per CLAUDE.md standard: Purpose, Arguments, Returns,
Constraints, Fails when, Example, Do not use when. Be explicit about
what happens when audio length does not divide evenly by
average_span_seconds (partial column dropped).

Pause for review.

Deliverable 2 — Tests
- Add test/test_ltsa.jl
- Test 1: synthetic signal with time-varying spectral content
  (e.g., bandlimited noise in band A for first half, band B for
  second half). Confirm LTSA shows the bands in the right columns
  with right relative levels.
- Test 2: energetic-mean equivalence — compare LTSA column for a
  known stationary signal to compute_psd with average_psd over the
  same span. Should match exactly (within floating-point precision).
- Test 3: calibration propagation — calibrated audio produces
  calibrated LTSA with correct is_calibrated flag and units symbol.
  NoCalibration audio produces uncalibrated LTSA.
- Test 4: column boundaries — confirm correct number of columns for
  audio of known length, confirm partial-column handling (dropped),
  confirm column_times vector matches.
- Test 5: edge cases — what happens with audio shorter than one
  column span? With zero-length audio? With average_span_seconds
  larger than audio? Document and assert expected behavior.
- Include in test/runtests.jl.

Pause for review.

Deliverable 3 — Documentation
- Add docs/src/explanations/ltsa.md
- Explain what LTSA is, when to use it (multi-hour to multi-day
  recordings, soundscape characterization), how to interpret it,
  parameter selection guidance (column span vs FFT window tradeoffs)
- Reference relevant DD entries
- Add to docs/make.jl pages list
- Add a DD entry (DD-25) documenting the energetic-mean and
  non-overlapping-column decisions

Pause for review.

Out of scope for this session
- Plot recipe — deferred to a separate small session after the
  computational layer is solid
- dB-mean or median column averaging — deferred unless needed
- Overlapping column boundaries — deferred unless needed

End of session
- Run full test suite, confirm no regressions
- Commit deliverables to feature/ltsa branch with descriptive
  commit messages per deliverable
- Push to origin/feature/ltsa
- Do NOT merge to dev — Robert reviews the full diff before merging
