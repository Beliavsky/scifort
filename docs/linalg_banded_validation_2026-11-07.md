# Independent validation of the dense banded solvers — 2026-11-07

## Implemented capabilities

The new `src/linalg_banded.f90` implements real/complex binary64 analogues of
`scipy.linalg.solve_banded`, `solveh_banded`, `cholesky_banded` and
`cho_solve_banded`, supporting vector and multiple right-hand sides, SciPy's
one-based general packed-band layout, and packed upper/lower positive-definite
layouts. General elimination uses *within-band* partial pivoting and tracks the
extra superdiagonals, rather than allocating the full coefficient matrix.

## Reproducible tests

Tests compiled with GNU Fortran 14.2 and compared numerically against the
SciPy 1.17.0 installed in the test environment. No SciPy library is linked
into the Fortran implementation.

| Test | Checked (-O0 -fcheck=all) | Optimized (-O3) |
| --- | ---: | ---: |
| Random seed 109, order 1–70, 750 candidates | 744 / 744 passed | 744 / 744 passed |
| Random seed 991, order 1–120, 650 candidates | 645 / 645 passed | 645 / 645 passed |
| Random seed 84, order 1–50, 200 candidates | 198 / 198 passed | 198 / 198 passed |
| Native real/complex pivoting, Cholesky and shape errors | Passed | Passed |
| Native tridiagonal order 50,000 | Passed | Passed |

Some random candidates were excluded before testing because their generated
coefficient matrix exceeded a condition number of `1e12` or SciPy could not
compute a reference solve; counts above are *actual comparisons*, not
candidate counts. The Python generator includes randomly varied bandwidths
(including full, diagonal, and triangular matrices), real/complex cases,
upper/lower Hermitian packed formats, factor reuse, 1–3 right-hand sides,
forced pivots, and ill-conditioned systems. The native tests also check
singular elimination, complex non-Hermitian diagonal rejection, incorrect
packed-band shapes, and factor values against an analytically chosen L.

- Seed 109: worst relative error against SciPy `6.97e-14`; worst normwise
  backward error `8.99e-17`.
- Seed 991: worst relative error against SciPy `1.48e-11`; worst normwise
  backward error `1.04e-16`. The largest forward discrepancy occurred on an
  intentionally nearly singular pivoting problem and is not evidence of
  equally strong forward accuracy for arbitrary ill-conditioned inputs.
- Seed 84: worst relative error `9.92e-15`; backward error `8.00e-17`.

The normalized backward error is
`||A*X-B||_F / (||A||_F*||X||_F + ||B||_F)`, whereas the
relative forward comparison uses `||X-X_scipy||_F / ||X_scipy||_F`.
These are not a proof of numerical accuracy in all cases.

The native order-50,000 tridiagonal problem verifies that the new module
stores only the active band. Its `3 x 50000` input has no dense `50000 x
50000` counterpart. It verifies the constant solution to approximately
`2e-13` absolute tolerance with bounds checking enabled.

Reproduce with:

```sh
python tools/linalg_reference/run_banded_checks.py --mode both --cases 750 --max-order 70 --seed 109 --stress
python tools/linalg_reference/run_banded_checks.py --mode both --cases 650 --max-order 120 --seed 991 --stress
python tools/linalg_reference/run_banded_checks.py --mode both --cases 200 --max-order 50 --seed 84 --stress
```

Additional pre-existing matrix-function regression comparisons were run in
checked and optimized modes: `funm` including adaptive Taylor, `sqrtm`/
`signm`, `logm`/fractional powers, the six trigonometric/hyperbolic matrix
functions, and `expm_cond`. They passed on focused test subsets.

**Build limitation:** The full FPM test suite was not run because FPM is
unavailable here. The standalone GNU Fortran build, using only `kinds.f90`,
`linalg.f90` and `linalg_banded.f90`, succeeds under `-std=f2018`.

## Implementation provenance and limitations

The kernel is an independently written realization of standard banded LU
elimination with partial pivoting and banded Cholesky; no SciPy or LAPACK
source was copied. Source and tests are MIT licensed. There is no dense
matrix reconstruction, runtime BLAS/LAPACK requirement, batching of
leading dimensions, iterative refinement or condition-number estimate.
Extremely ill-conditioned problems require extra care even when status is
zero. See [public API and examples](linalg_banded.md).
