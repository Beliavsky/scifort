# Reuse review: Jason Christopherson's signal and dynamics libraries

Reviewed October 19, 2026. SciFort's objective is SciPy algorithm coverage,
not reproducing entire independent Fortran libraries or vendoring dependencies.

| Upstream project | Repository license | Dependencies listed by its FPM manifest | Preferred action |
|---|---|---|---|
| `peaks` | MIT, copyright 2023 Jason Christopherson | standalone peak detection | Add an independently written excursion detector, with credit to Christopherson and Eli Billauer |
| `spectrum` | MIT, copyright 2023-2026 Jason Christopherson | `linalg` | Reuse *algorithms* selectively; avoid duplicating SciFort's existing spectrum routines |
| `dynamics` | MIT, copyright 2022-2026 Jason Christopherson | `linalg`, `nonlin`, `diffeq`, `spectrum`, `fstats`, `peaks`, with BLAS/LAPACK link requirements | Prefer a few relevant FRF/system-identification/modal-analysis routines over the whole package |

Original repositories and manifests:

- https://github.com/jchristopherson/peaks
- https://github.com/jchristopherson/spectrum
- https://github.com/jchristopherson/spectrum/blob/main/fpm.toml
- https://github.com/jchristopherson/dynamics
- https://github.com/jchristopherson/dynamics/blob/main/fpm.toml

## Implemented in this SciFort update

- `argrelextrema(x, comparator, order, mode, status)` for one-dimensional
  `real(dp)` signals. The comparator is a pure scalar logical procedure.
- `argrelmax` and `argrelmin`, using one-based index arrays.
- `peak_detect(x, delta, imx, mx, imn, mn, status)` as a **SciFort extension**,
  using the excursion-based max/min state machine described by Eli Billauer
  and used in `peaks`. It is not a SciPy-compatible synonym for `find_peaks`.
  It requires a strictly positive `delta` and finite samples, returns only
  confirmed swings, and omits the last unconfirmed one.

No source files were copied from the upstream libraries. The attribution and
original MIT license notice are recorded in CODE_PROVENANCE.md and
THIRD_PARTY_LICENSES.md. No extra FPM dependencies have been introduced.

### Verification

- All 44 signal-related Fortran source modules compile under GNU Fortran 14.2.0.
- 76/76 standalone signal tests pass with `-O1 -fcheck=all`.
- 76/76 standalone signal tests pass with `-O3 -fcheck=all`.
- The two targeted new tests pass at `-O3` with floating-point traps.
- 450 randomized cases provide 1,350 direct SciPy 1.17 comparisons for
  `argrelextrema`, `argrelmax`, and `argrelmin` plus 450 separate comparisons
  to a Python Billauer-style state-machine reference for `peak_detect`.
- Project-wide `fpm test` was not run; FPM is not installed in the execution
  environment. Do not substitute standalone signal testing for a full FPM run.

## Next sensible gaps

1. `scipy.signal.find_peaks_cwt`, for wavelet-based peak finding. This requires
   different underlying logic from both SciFort's current `find_peaks` and the
   Billauer-style excursion method.
2. `scipy.signal.place_poles`, if not already provided under another name.
   It would complement SciFort's existing state-space and LTI response code.
3. Expand spectral analysis by axis/channel and cover more SciPy-supported
   window and averaging choices. SciFort already credits `spectrum` for its
   original 1D Welch and cross-spectrum algorithm design.
4. Consider specific `dynamics` modal-fitting or frequency-response kernels
   only after matching their functionality to SciPy's public APIs and checking
   dependencies on a per-file basis.

## Subsequent algorithm coverage

SciFort next implemented `find_peaks_cwt` and analytic `place_poles` without
incorporating `peaks`, `spectrum`, or `dynamics` source code or adding a new
FPM dependency. Wavelet ridge-detection conventions are adapted from SciPy
1.17 (BSD-3-Clause), with explicit attribution in CODE_PROVENANCE.md and
THIRD_PARTY_LICENSES.md. The advanced partially actuated MIMO case remains
on the roadmap.
