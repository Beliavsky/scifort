# Complex sparse exponential validation — 2026-11-06

Test environment: GNU Fortran (checked `-O0 -fcheck=all`, optimized `-O3`),
NumPy and SciPy 1.17.0.

- Native new-feature regression: `test_sparse_complex_expm` passes with
  bounds checking and in a fully optimized `-O3` build. Covers CSR/CSC actions,
  matrix-free callbacks, a 4,096-element matrix-free stress case, multiple right-hand
  sides, exact Jordan-block exponentials, diagonal matrices, positive and
  negative times, sampled intervals (including endpoint exclusion), zero-time
  identities, failure propagation, invalid tolerance/bounds and nonfinite data.
- Full bounds-checked standalone sparse test runner: **53/53 passed**.
- Randomized reference seeds 20261106 (180 cases per mode) and 90506 (250
  cases per mode): **430/430 per compiler mode**. Each case tests six action
  paths: single-time CSR, CSC, and matrix-free, plus interval CSR, CSC, and
  matrix-free. This yields **2,580 comparisons per mode, 5,160 total**.
- Maximum observed relative error versus independent dense exponentials:
  at most 5.8e-14 across the listed runs, with magnitudes varying widely.
- SciPy `expm_multiply` is also evaluated as a cross-check; the independent
  dense exponential is the primary reference when SciPy's interval stepping
  loses accuracy on highly nonnormal matrices.

- Full optimized standalone native sparse suite: **not completed**. It
  exceeded the execution limit while compiling the pre-existing large
  `src/sparse_linalg.f90`, before producing test results. The new complex
  exponential source and native test were independently built at `-O3`,
  and the SciPy comparisons compiled them at `-O3` as well.

The project-wide FPM suite was not run in this environment.
