# Complex sparse LSMR validation — 2026-11-06

## Implementation and independent reference

- `src/sparse_linalg_complex_lsmr.f90`: canonical complex(dp) CSR/CSC LSMR.
- `test/test_sparse_complex_lsmr.f90`: rectangular/square exact systems,
  inconsistent/rank-deficient systems, nonzero initial guesses, centered
  Tikhonov damping, invalid input, and 6,000-by-120 sparse stress problem.
- `tools/sparse_reference/complex_lsmr_driver.f90` and
  `compare_sparse_complex_lsmr_scipy.py`: streaming independent comparisons
  with SciPy 1.17.0 and NumPy 2.3.5.

## Randomized numerical comparisons

| Seed | Cases per build | Checked `-fcheck=all` | Optimized `-O3` | Max difference from SciPy |
| --- | ---: | --- | --- | ---: |
| 20261106 | 360 | 360 passed | 360 passed | 4.259e-11 |
| 246802 | 220 | 220 passed | 220 passed | 4.437e-11 |
| 123456 | 256 | 256 passed | 256 passed | 5.414e-11 |

**Total:** 836 random cases per mode, 1,672 across both modes. The final
seed includes zero-sized dimensions. Each case checks the computed solution
against SciPy's LSMR, the minimum-norm dense SVD solution of the centered
regularized system, and direct computation of the returned true residual and
normal-equation residual. The maximum difference from SciPy was 5.414e-11.

The randomized cases include CSR and CSC matrices, overdetermined and
underdetermined systems, consistent and inconsistent right-hand sides,
explicitly rank-deficient cases, zero rows/columns, zero matrices, damped
systems and arbitrary complex nonzero starting guesses.

## Native sparse regression suite

- **48/48** native programs passed with `-O0 -g -fcheck=all -fbacktrace`.
- **48/48** native programs passed with **all library modules compiled at
  `-O3`** and native test drivers compiled at `-O0`.
- New LSMR focused test compiled and passed separately with **both** the
  LSMR and test drivers at `-O3`.
- The focused test's large 6,000-by-120 problem is built from 18,000 COO
  triplets and does not construct the dense matrix.

## Reproduction

```sh
python tools/sparse_reference/compare_sparse_complex_lsmr_scipy.py --mode both --cases 360
python tools/sparse_reference/compare_sparse_complex_lsmr_scipy.py --mode both --cases 220 --seed 246802
python tools/sparse_reference/compare_sparse_complex_lsmr_scipy.py --mode both --cases 256 --seed 123456
python tools/sparse_reference/run_native_sparse.py --mode checked
python tools/sparse_reference/run_native_sparse.py --mode optimized --fast-tests
```

GNU Fortran 14.2.0 was used. `fpm test` for the entire SciFort project was
not run because FPM is not available in this environment. The standalone
native runner covers sparse modules only; the SciFort modules outside that
suite are not claimed to have been regression-tested here.

The LSMR recurrence inherits BSD-3-Clause provenance from the prior real
implementation; it is not represented as new MIT-only code.
