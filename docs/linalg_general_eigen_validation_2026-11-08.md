# General nonsymmetric eigenproblem validation — 2026-11-08

## Environment

- GNU Fortran (`gfortran`) with Fortran 2018 source, runtime checks and
  `-ffpe-trap=invalid,zero,overflow` in checked mode; `-O3` optimized mode.
- Local SciPy 1.17.0 / NumPy 2.3.5, with one BLAS thread.
- Core source dependencies: `src/kinds.f90`, `src/linalg.f90`,
  `src/linalg_matrix_functions.f90`, `src/linalg_general_eigen.f90`.
- No LAPACK or external numerical libraries linked into the Fortran tests.

## Tests completed

The independent `run_general_eigen_checks.py` driver compares SciFort and SciPy
on real/complex nonsymmetric matrices with 0..40 rows, including zero,
tridiagonal, triangular, diagonal, Hermitian, skew-Hermitian, strictly upper
triangular, clustered, repeated, and wildly scaled cases.

| Input seed | Cases in checked build | Cases at `-O3` | Largest scaled eigenvalue difference | Largest normalized eigenvector residual |
|---|---:|---:|---:|---:|
| 99181 | 540 | 540 | 1.704e-14 | 9.299e-15 |
| 20261108 | 540 | 540 | 7.941e-15 | 7.105e-15 |
| **Total** | **1080** | **1080** | **1.704e-14** | **9.299e-15** |

Additional development runs, including a lower-triangular high-nonnormality
regression that originally revealed an avoidable spectral perturbation,
passed after a triangular fast path was added.

`test/test_linalg_general_eigen.f90` passed in both modes. It verifies
analytic complex-conjugate eigenvalues, normalized left/right eigenvector
residuals, complex nonsymmetric matrices, zero matrices, zero dimensions,
invalid shapes, the 40-order cap, and nonfinite inputs.

Neighboring regressions completed in both modes:

- `run_hessenberg_checks.py --cases 65 --max-order 25 --seed 431`
- `run_dense_eigen_checks.py --cases 65 --max-order 24 --seed 525`

The checked `funm` native test and a 35-case SciPy batch passed; a larger
combined `funm` regression did not finish before its time limit and is **not**
claimed as passing. Full `fpm test` is untested because FPM is unavailable.

## Numerical interpretation

The spectrum is unordered: comparisons match each computed eigenvalue to a
SciPy eigenvalue by minimum-cost bipartite assignment. Eigenvectors are checked
through their normalized left and right residuals, not componentwise
comparison, because scaling and basis choices are nonunique. No guarantee is
made for individual eigenvalues of severely defective or very nonnormal
matrices, even when their eigenpair backward residuals are small.

## Reproduction

```bash
python tools/linalg_reference/run_general_eigen_checks.py \
    --mode both --cases 540 --max-order 40 --seed 99181
python tools/linalg_reference/run_general_eigen_checks.py \
    --mode both --cases 540 --max-order 40 --seed 20261108
```
