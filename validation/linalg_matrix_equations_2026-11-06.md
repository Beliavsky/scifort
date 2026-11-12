# Validation: dense Sylvester and Lyapunov matrix equations

Date: 2026-11-06. Reference: NumPy 2.3.5 and SciPy 1.17.0.
Compiler: GNU Fortran 14.2.0, checked (`-O0 -g -fcheck=all`) and optimized (`-O3`).

Added `scifort_linalg_matrix_equations` and tested `solve_sylvester`,
`solve_continuous_lyapunov`, and `solve_discrete_lyapunov`, each with real and
complex binary64 matrices.

## Tests actually performed

- `python tools/linalg_reference/run_native_linalg.py --mode both`:
  passed both checked and optimized native regression tests. Includes 2x3
  Sylvester solutions, nontrivial real/complex Lyapunov solutions, Hermitian
  transpose, forced row pivoting, singular equations, invalid shapes,
  oversized matrices, and NaN input.
- `compare_matrix_equations_scipy.py --mode both --cases 360 --seed 20261106`:
  360/360 per build. Maximum relative difference from SciPy `3.843e-13`,
  maximum true-equation relative residual `1.287e-12`.
- `--mode both --cases 1200 --seed 20261106`:
  1200/1200 per build. Maximum relative difference from SciPy `3.843e-13`,
  maximum relative residual `1.287e-12`. (Includes the first 360 cases;
  those 360 should not be counted twice in unique-case totals.)
- `--mode both --cases 1200 --seed 20261107`:
  1200/1200 per build, max difference `2.617e-15`, residual `1.789e-15`.
- `--mode both --cases 160 --seed 20261108 --max-dim 10`:
  160/160 per build, max difference `6.426e-15`, residual `1.735e-15`.
- `example/linalg_matrix_equations_example.f90` compiled and ran:
  maximum residual `6.9389e-18`.
- Standalone module compiled with `gfortran -std=f2018 -pedantic -fcheck=all`
  and no extended line-length option; native test passed.

Total *distinct* randomized scenarios: **2,560 per compiler mode**, or
**5,120 checked/optimized executions**. The 360-case seed-20261106 run is a
subset of the 1200-case run and is not double-counted. The reference script
checks both numerical agreement with SciPy and the true matrix-equation
residual. The documented tolerance is `1e-10` for both relative comparisons.

## Scope and limitations

The full project-wide `fpm test` and the pre-existing standalone sparse
regression runner were **not** run for this extension. The new source does not
modify any existing numerical module. The method is a pivoted dense
Kronecker direct solve restricted to at most 400 scalar unknowns, not a
scalable Schur-based solver. This is particularly important for large or
nearly singular matrix equations.
