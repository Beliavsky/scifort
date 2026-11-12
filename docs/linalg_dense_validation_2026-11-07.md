# Dense triangular / LU numerical validation (2026-11-07)

This update adds `scifort_linalg_dense` and reference harnesses for the public
SciPy-style `solve_triangular`, `lu_factor`, and `lu_solve` generics. The new
Fortran implementation uses independently coded forward/back substitution
and Gaussian elimination with partial pivoting.

## Commands and results

GNU Fortran 14.2.0, NumPy 2.3.5, and SciPy 1.17.0 were available. With
`fpm` not installed, the following independent standalone GNU Fortran runs
were used. Each mode compiles and runs the native Fortran assertions and a
Fortran text-I/O comparison driver, then compares all solutions to SciPy:

```sh
python tools/linalg_reference/run_dense_checks.py --mode both --cases 450 --seed 217 --max-order 60
python tools/linalg_reference/run_dense_checks.py --mode both --cases 600 --seed 218 --max-order 80
python tools/linalg_reference/run_dense_checks.py --mode both --cases 70 --seed 902 --max-order 60 --stress
python tools/linalg_reference/run_dense_checks.py --mode both --cases 220 --seed 221 --max-order 40
```

| Run | Checked | `-O3` | Worst relative solution difference | Worst backward residual |
| --- | ---: | ---: | ---: | ---: |
| Seed 217 | 450/450 | 450/450 | 1.651e-14 | 8.580e-17 |
| Seed 218 | 600/600 | 600/600 | 5.723e-15 | 1.346e-16 |
| Stress seed 902 | 70/70 | 70/70 | 1.723e-15 | 7.103e-17 |
| Seed 221, additional failure paths | 220/220 | 220/220 | 5.011e-14 | 1.088e-16 |
| **Total** | **1,340/1,340** | **1,340/1,340** | | |

These 2,680 independent comparisons cover real and complex matrices,
`n=1` through stress-order 256, one through four right-hand sides, forced
pivoting, and `'N'/'T'/'C'` solves. Scalar complex conjugation is tested
independently of ordinary transposition. Triangular systems test both upper
and lower storage, implicit unit diagonals, and the ignored opposite triangle.
The LU driver also solves a *second* RHS with the same factors in every case.

The comparison error is
`||X_Fortran-X_SciPy||_F / max(1, ||X_SciPy||_F)`; normwise backward residual
is `||op(A)*X-B||_F / (||op(A)||_F*||X||_F + ||B||_F)`.

The native test includes invalid output shapes, malformed pivot arrays,
unsupported transpose codes, exact singular matrices and nonfinite inputs.
All tested builds passed. The test suite cannot establish forward stability
for very ill-conditioned matrices, and it does not provide a formal proof
of correctness or performance advantage over optimized BLAS/LAPACK routines.

The project-wide `fpm test` was not run because `fpm` is unavailable in
this environment. No new compiled runtime dependency is required.
