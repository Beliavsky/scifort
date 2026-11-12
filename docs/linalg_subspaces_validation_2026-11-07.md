# Dense subspaces: SciPy comparison and native validation (2026-11-07)

## Tested environment

- GNU Fortran 14.2.0 with `-std=f2018`, separately compiled with
  `-O0 -g -fcheck=all -fbacktrace` and `-O3`.
- Python 3.13.5, NumPy 2.3.5, and SciPy 1.17.0.
- Standalone driver compilation uses, in order, `src/kinds.f90`,
  `src/linalg.f90`, `src/linalg_svd.f90`, and
  `src/linalg_subspaces.f90`; no external BLAS/LAPACK link flags.

## Randomized comparisons

Three commands passed in **both** checked and optimized configurations:

```sh
python tools/linalg_reference/run_subspaces_checks.py --mode both --cases 300 --max-order 32 --seed 1029 --stress
python tools/linalg_reference/run_subspaces_checks.py --mode both --cases 280 --max-order 64 --seed 2092 --stress
python tools/linalg_reference/run_subspaces_checks.py --mode both --cases 180 --max-order 84 --seed 3847 --stress
```

The **760 distinct cases in each compiler mode** give **1,520 executed
SciPy comparison cases**. Input families include real/complex and mixed-type
matrices; full-rank, exact-rank-deficient, duplicated-column, all-zero, and
geometrically graded spectra; tall, wide, square and empty shapes; contained
subspaces, near-parallel (1e-9/1e-5 radians) and near-orthogonal subspaces;
`1e-100`/`1e100` scaling; and both the default and explicit `rcond` cutoff.

Since the orthonormal basis itself is nonunique, tests compare **orthogonal
projectors**, not basis coefficients. The Fortran range projector is compared
with `scipy.linalg.orth(A) @ orth(A).conj().T`; the kernel projector is compared
similarly using `scipy.linalg.null_space(A)`. Further checks verify
`Q^H Q = I`, `Z^H Z = I`, numerical kernel residuals, returned dimensions,
and descending principal angles compared with SciPy on bases computed using
identical rank cutoffs.

Maximum errors across the three completed commands and two compiler modes:

| Statistic | Maximum observed |
| --- | ---: |
| Norm-relative range-projector discrepancy | `1.69e-10` |
| Norm-relative null-projector discrepancy | `2.00e-10` |
| `||A Z||_F / ||A||_F` | `1.22e-15` |
| Dimension-normalized orthogonality discrepancy | `3.18e-15` |
| Absolute principal-angle difference (radians) | `2.11e-8` |

The principal-angle maximum is an absolute rather than relative error;
this matters for zero or nearly zero angles. Comparisons use an absolute
`2e-7` tolerance for angles and norm-relative projector differences.

## Native Fortran tests and existing regressions

`test/test_linalg_subspaces.f90` passed standalone in both the checked and
optimized builds. It includes analytic wide rank-deficient real and complex
examples, null-space projector and orthonormality checks, exact zero and
empty dimensions, custom rank cutoffs, invalid options, nonfinite input,
row-count mismatches, mixed real/complex arguments, and angles within
`1e-9` radians of zero and of 90 degrees.

`example/linalg_subspaces_example.f90` compiled and ran with checked GNU
Fortran flags. It returned rank 2, nullity 3, a maximum null residual
approximately `1.1e-16`, and the expected 30-degree principal angle.

The following existing randomized regression commands passed in both modes:

```sh
python tools/linalg_reference/run_svd_checks.py --mode both --cases 120 --max-order 32 --seed 710
python tools/linalg_reference/run_qr_checks.py --mode both --cases 70 --max-order 30 --seed 431
python tools/linalg_reference/run_dense_checks.py --mode both --cases 25 --max-order 24 --seed 218
python tools/linalg_reference/run_dense_cholesky_checks.py --mode both --cases 20 --max-order 24 --seed 390
```

These constitute 240 SVD, 140 QR, 50 dense-LU, and 40 dense-Cholesky
SciPy comparisons (each counting both compiler modes). The existing test
runners also execute their corresponding native tests.

The full repository-wide `fpm test` was **not** run: FPM is not installed in
this environment. It also depends on dependencies managed by the project.

## Numerical limits

The standalone one-sided Jacobi SVD and the two-pass Gram-Schmidt completion
are designed for small-to-medium dense matrices, rather than large sparse
problems. Matrices with singular values near the numerical rank threshold
may have different effective ranks under SciPy's LAPACK-backed method. A
basis has arbitrary phases/signs and should be tested via its projector.
