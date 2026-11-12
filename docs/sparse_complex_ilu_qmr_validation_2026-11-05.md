# Complex ILU(0), BiCG and QMR validation (2026-11-05)

Environment: GNU Fortran 14.2.0, NumPy 2.3.5, SciPy 1.17.0.
The Fortran source requires no external numerical libraries. SciPy and NumPy
are used by the optional Python reference test only.

## Native regression

The complete standalone sparse runner:

```bash
python tools/sparse_reference/run_native_sparse.py --mode checked
```

passed **44/44 native sparse regression programs** with
`-O0 -g -fcheck=all -fbacktrace`. This includes all 43 prior programs and
`test/test_sparse_complex_ilu_bicg_qmr.f90`, which tests ordinary and
Hermitian-adjoint ILU solves, 2-D CSR/CSC BiCG/QMR, preconditioned BiCGSTAB
and GMRES, invalid pivots/diagonals, diagonal shifts, zero-sized systems and
invalid preconditioner combinations.

An attempt at the standalone runner with `--mode optimized` exceeded the
execution time limit **during compilation**, before native optimized
suite results were available. No full `-O3` sparse-suite result is claimed.

The focused independent checker compiles **all of its nine supporting
Fortran modules** in the selected compiler mode, then runs the new native
regression and the external randomized references:

```bash
python tools/sparse_reference/run_complex_ilu_qmr_checks.py --mode both --per-solver 180
```

| Build | Native focused test | SciPy/NumPy randomized comparisons |
|---|---|---:|
| `-O0 -g -fcheck=all` | Passed | 1,800/1,800 |
| `-O3` | Passed | 1,800/1,800 |

Both compiler modes passed 900 cases for seed 1728 and 900 for seed 9894.
The maximum norm of the solution error in each mode was approximately
`8.32e-10`. A separate 220-case nonsymmetric indefinite-system stress
comparison also passed with the updated complex solvers compiled at `-O3`.

The reference driver accepts streamed, unsorted, sometimes duplicated COO
coordinates with explicit zeros. The Python reference script compares
ordinary and adjoint ILU(0) solves against an independently assembled
Python zero-fill LU factor, all four Krylov solvers against dense NumPy
solutions, and non-ILU solver cases against SciPy's `bicg`, `qmr`,
`bicgstab` and `gmres`. CSR/CSC wrappers, Jacobi/ILU/identity choices,
nonzero initial guesses, complex coefficients, and true residual reporting
are exercised.

## Scope and limitations

- Complex values use `complex(dp)` (binary64 real and imaginary parts);
  CSR/CSC indices use default integers.
- ILU(0) is unpivoted and uses the existing sparsity pattern. Missing
  diagonal entries, very small pivots, or nonfinite results may fail.
  SciPy's SuperLU `spilu` is a different, more general algorithm.
- BiCG and QMR can break down on general nonsymmetric systems; no
  look-ahead or automatic pivoting is implemented.
- QMR uses a single left ILU/Jacobi preconditioner, not SciPy's two separate
  left/right `M1`/`M2` operators.
- No full project-wide `fpm test` was run because FPM is not installed.
