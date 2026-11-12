# Generalized Schur validation — November 8, 2026

## Build and comparison method

- GNU Fortran: checked `-O0 -fcheck=all -ffpe-trap=invalid,zero,overflow`
  and optimized `-O3`, Fortran 2018.
- Independent numerical reference: NumPy 2.3.5 and SciPy 1.17.0,
  `scipy.linalg.qz(..., output='complex')` and `eigvals`.
- Two independent seeds: `31209` with 680 input pairs and `88402` with
  700 input pairs in **each** build mode (2,760 pair/build comparisons).
- Order 0–40; real, complex, and mixed coefficient pairs; diagonal,
  triangular, dense, rank-deficient, singular, and extreme-scale inputs.
- Native tests for rectangular inputs, undersized output, NaN input,
  singular pencils, singular-B and both-singular regular pairs.

## Completed measurements

| Metric | Observed maximum |
|---|---:|
| Normalized reconstruction of A or B | 1.265e-10 |
| Unitarity residual (Frobenius) | 4.33e-14 |
| Projective eigenvalue disagreement for highly nonnormal cases | 2.82e-3 |
| Generalized eigenvalue backward error in those difficult cases | 1.6e-17 |

Most matrices have substantially smaller reconstruction and eigenvalue
comparison differences. In a few severely nonnormal cases, two backward-stable
eigensolvers disagree noticeably about the eigenvalue locations. These
cases were accepted **only after** checking their eigenvalue backward errors
against the original homogeneous pencil, not by relaxing every tolerance.
Homogeneous/projective comparisons are used for infinite eigenvalues and
other singular-B cases; direct `alpha/beta` comparison is not used when
`beta` vanishes.

Reproduce in the project root:

```bash
python tools/linalg_reference/run_qz_checks.py --mode both --cases 680 --max-order 40 --seed 31209
python tools/linalg_reference/run_qz_checks.py --mode both --cases 700 --max-order 40 --seed 88402
```

This is an inverse-based generalized Schur construction. It is **not** an
implementation of LAPACK's iterative QZ algorithm, and the tests do not
establish stability for arbitrary nearly singular regular pencils.
The complete `fpm test` suite was not run (FPM not installed).
