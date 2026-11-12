# Banded LU factorization and factor reuse: independent validation (2026-11-07)

## What's new

`scifort_linalg_banded` now exposes two generics:

- `lu_factor_banded`: partial-pivoted real/complex LU factorization of a packed
  general banded matrix, in expanded `upper+2*lower+1` storage with one-based
  integer pivot indices.
- `lu_solve_banded`: solves normal, transpose, or conjugate-transpose systems
  against a previously computed factor, for one or many RHS columns.

`solve_banded` has been refactored to call these new kernels. The public API
for the pre-existing one-step solver is unchanged. The algorithm reuses the
same row-exchange sequence on the right-hand side that is used in the LU
factorization, and reverses that sequence for transposed systems.

## Reproducible tests

The following commands were run with GNU Fortran 14.2.0 and SciPy 1.17.0.
Native tests and independent SciPy comparisons were run for both checked
`-O0 -g -fcheck=all -fbacktrace` and optimized `-O3` builds:

```sh
OPENBLAS_NUM_THREADS=1 python tools/linalg_reference/run_banded_lu_checks.py --mode both --cases 600 --max-order 85 --seed 171 --stress
OPENBLAS_NUM_THREADS=1 python tools/linalg_reference/run_banded_lu_checks.py --mode both --cases 550 --max-order 100 --seed 809 --stress
OPENBLAS_NUM_THREADS=1 python tools/linalg_reference/run_banded_checks.py --mode both --cases 650 --max-order 120 --seed 991 --stress
OPENBLAS_NUM_THREADS=1 python tools/linalg_reference/run_banded_checks.py --mode both --cases 750 --max-order 70 --seed 109 --stress
```

| Test | Checked | Optimized |
| --- | ---: | ---: |
| Seed 171: 600 LU factorizations, 6 comparisons per factor | 3,600 / 3,600 | 3,600 / 3,600 |
| Seed 809: 550 LU factorizations, 6 comparisons per factor | 3,300 / 3,300 | 3,300 / 3,300 |
| Previous banded test seed 991 | 645 / 645 | 645 / 645 |
| Previous banded test seed 109 | 744 / 744 | 744 / 744 |
| Native LU tests (forced pivots, N/T/C, failures) | Passed | Passed |
| Native LU and original one-step order-50,000 systems | Passed | Passed |

**13,800** new SciPy comparisons passed across both builds, plus **2,778**
regressions from the previous banded-suite cases. Each LU case tests real or
complex arithmetic, a random band, 1–3 RHSs, both of two independent RHS
batches using a single factorization, and each of the `N`, `T` and `C`
operations. Test cases cover diagonal, one-sided triangular, narrow, wider,
and sometimes full bands. Pivot-forcing cases include zero-near diagonal
entries and later pivot-row swaps. The stress generator includes matrices
of order 80, 120, 200, and 350; native stress tests additionally factor an
order-50,000 tridiagonal system with `O(n)` memory.

New LU comparison results:

- Seed 171: maximum relative difference from SciPy `1.247e-11`, maximum
  normwise backward residual `1.579e-16`.
- Seed 809: maximum relative difference `7.928e-11`, maximum backward
  residual `1.726e-16`.
- Both builds gave equivalent accuracy on these cases.

A test comparison's relative difference is `||X-X_scipy||_F / ||X_scipy||_F`.
The backward residual is `||op(A)*X-B||_F / (||op(A)||_F*||X||_F+||B||_F)`.
The SciPy reference for `N`/`T`/`C` uses independent dense `scipy.linalg.solve`
on the original matrix and appropriate (conjugate) transpose, not the
SciFort LU factorization. No dense matrix is created *inside* the Fortran
routines. These tests cannot guarantee small forward errors for extremely
ill-conditioned matrices. The generator bounds the condition number of
small/medium random reference matrices by `1e7` to make forward agreement a
meaningful check; large matrices are deliberately made diagonally dominant
(except for a controlled initial pivot case).

Other dense matrix-function tests also passed in checked and optimized builds:
`run_matrix_functions_checks.py` with 80 cases and `run_funm_checks.py` with
60 general and 104 adaptive-Taylor cases per mode. The Fortran
`example/linalg_banded_example.f90` compiled and ran successfully.

## Limitations

- Exact-zero pivot detection, without a reciprocal condition estimate;
  nearly singular matrices can yield large forward errors.
- Partial pivoting occurs only within `lower` rows below the current pivot,
  as required for banded Gaussian elimination.
- The same LU factors support `A`, `A^T` and `A^H`, but no optional iterative
  refinement, equilibration, or reliable error bound is provided.
- Factor matrices and pivot arrays are caller-managed and must have the
  correct expanded shapes. They are zeroed if factorization fails.
- The project-wide `fpm test` suite was not run here because FPM is not
  installed. Standalone native builds passed under `-std=f2018`.

Source, tests, documentation and harness are MIT-licensed original SciFort
work. SciPy provides only numerical reference outputs; no external
BLAS/LAPACK code or additional runtime dependency is introduced.
