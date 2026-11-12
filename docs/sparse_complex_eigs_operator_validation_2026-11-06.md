# Validation — complex matrix-free Arnoldi and reusable shift-invert LU

Validated on 2026-11-06 with GNU Fortran 14.2.0, Python 3.13,
NumPy 2.3.5, and SciPy 1.17.0.

## Native Fortran tests

- **52/52 registered standalone sparse regression programs passed** under
  `gfortran -std=f2018 -O0 -g -fcheck=all -fbacktrace`, via
  `python tools/sparse_reference/run_native_sparse.py --mode checked`.
- The new `test_sparse_complex_eigs_operator.f90` and the existing
  `test_sparse_complex_eigs.f90` and `test_sparse_complex_eigs_shift.f90`
  also passed with **both the library and test programs compiled at `-O3`**.
- The new documentation example compiled and ran at `-O3`; it computed
  `lambda=(3.0,0.3)` with residual approximately `3.937e-16`.
- The native test includes a **2,500-dimensional matrix-free** complex
  upper-bidiagonal operator, finding two interior eigenvalues near a shift.
  Its operator is never assembled as a CSR or dense matrix in this test.
- The native suite checks all six selection modes, provided initial guesses,
  real/imaginary true residuals, LU shift-cache reuse across CSC/CSR,
  invalidating a cache after a matrix or shift change, singular shifts,
  an unavailable inverse callback, invalid tolerances and `ncv`, and
  failed/nonfinite user callbacks.

The full `run_native_sparse.py --mode optimized --fast-tests` did **not**
complete: compilation of the pre-existing large `src/sparse_linalg.f90`
module exceeded the execution time limit before any optimized suite results
were produced. Therefore the result **does not** claim a complete optimized
regression pass. All modified eigensolver modules and their three focused
programs *were* separately compiled and tested at `-O3`.

## Independent randomized numerical comparisons

The script `tools/sparse_reference/compare_sparse_complex_eigs_operator_numpy.py`
constructs randomized nonsymmetric **complex upper-bidiagonal** operators.
Their complete spectra are known independently from their diagonals; no
Fortran eigensolver output is used in generating the reference eigenvalues.
It compares both the matrix-free callback and the reusable cached-CSR path,
including original-system residual norms and SciPy shift-invert checks.

| Seed | Cases per build | Checked | Optimized |
|---|---:|---|---|
| `20261106` | 120 | Passed, max error `7.371e-14` | Passed, max error `8.614e-14` |
| `20261107` | 240 | Passed, max error `8.050e-14` | Passed, max error `6.442e-14` |
| `20261108` | 90 | Passed, max error `6.962e-14` | Passed, max error `4.494e-14` |

**Total:** 450 randomized cases per compiler mode (900 across the two),
with two independent solve paths verified in each case. In the third batch,
**36 direct SciPy `scipy.sparse.linalg.eigs` comparisons per build** also
passed (72 across checked and optimized modes). The largest eigenvalue
error against the independently computed reference spectrum was below
`8.7e-14` across all tests.

Reproduce the SciPy batch with:

```sh
python tools/sparse_reference/compare_sparse_complex_eigs_operator_numpy.py \
    --mode both --cases 90 --seed 20261108 --scipy-samples 36
```

## Limitations

This is an explicitly restarted complex Arnoldi/Rayleigh–Ritz algorithm,
not ARPACK implicit restarting. `SM` requires `ncv=n`; an inexact user
inverse may harm spectral selection even when returned original-system
residuals are small. Matrix-free `norm_estimate` influences the acceptance
scale, so supplying an unreasonably large estimate can loosen convergence.
Factor caching duplicates one CSR matrix for strict matching. No
project-wide FPM suite was run because FPM is not installed in this
execution environment.
