# Schur matrix-equation validation — November 7, 2026

## What changed

The existing real/complex `scifort_linalg_matrix_equations` generics gained an
optional `method='auto'|'direct'|'schur'` argument. `auto` retains the original
Kronecker system for at most 400 scalar unknowns; beyond that it uses a
Bartels–Stewart style triangular reduction. Schur reductions are restricted to
coefficient matrix orders no greater than 40.

The three supported equations are:

- `solve_sylvester`: `A X + X B = Q` (two Schur decompositions)
- `solve_continuous_lyapunov`: `A X + X A^H = Q` (one Schur decomposition)
- `solve_discrete_lyapunov`: `X - A X A^H = Q` (one Schur decomposition)

All Fortran implementations are dependency-free SciFort code. SciPy 1.17.0
provides reference results only, never runtime source or linkage.

## Focused SciPy comparisons

The reproducible command:

```bash
python tools/linalg_reference/run_matrix_equations_schur_checks.py --mode both --cases 60 --max-order 24 --seed 404
python tools/linalg_reference/run_matrix_equations_schur_checks.py --mode both --cases 60 --max-order 40 --large --seed 711
python tools/linalg_reference/run_matrix_equations_schur_checks.py --mode both --cases 72 --max-order 32 --seed 902
python tools/linalg_reference/run_matrix_equations_schur_checks.py --mode both --cases 40 --max-order 40 --large --seed 933
python tools/linalg_reference/run_matrix_equations_schur_checks.py --mode both --cases 30 --max-order 32 --seed 1999
```

| Problem group | Checked | Optimized | Maximum relative solution difference (across modes) |
|---|---:|---:|---:|
| Orders 1–24, seed 404 | 60/60 | 60/60 | 8.31e-15 |
| Orders 21–40, seed 711 | 60/60 | 60/60 | 9.97e-15 |
| Orders 1–32, seed 902 | 72/72 | 72/72 | 8.50e-15 |
| Orders 21–40, seed 933 | 40/40 | 40/40 | 9.37e-15 |
| Orders 1–32, seed 1999 | 30/30 | 30/30 | 7.82e-15 |
| **Total** | **262/262** | **262/262** | **9.97e-15** |

The Python comparison harness runs small SciPy references with a single
BLAS thread to avoid host-dependent thread oversubscription. The randomized
checked comparison harness uses `-O1 -g -fcheck=all` (with
runtime checking); optimized uses `-O3`. Native tests additionally use
`-O0 -g -fcheck=all` and the example was run with floating-point traps.

The maximum normalized backward residual from these tests was
approximately `2.76e-15` (different normalizations apply in each equation).
The cases covered real and complex data, shifted random matrices, normal
matrices, upper-triangular matrices with repeated diagonal eigenvalues,
modestly nonnormal and bidiagonal Jordan-like examples, nonsquare Sylvester
solutions, and all three equations.

## Legacy compatibility and other regressions

```bash
python tools/linalg_reference/run_native_linalg.py --mode both
python tools/linalg_reference/compare_matrix_equations_scipy.py --mode both --cases 180 --max-dim 8 --seed 101
python tools/linalg_reference/compare_riccati_scipy.py --mode both --cases 45 --max-n 4 --seed 103
```

- Native matrix-equation, new Schur matrix-equation, and Riccati native tests:
  all passed in checked and optimized builds.
- Original matrix-equation SciPy comparisons: 180/180 each build,
  maximum relative difference `2.44e-15`.
- Riccati SciPy comparisons: 45/45 each build,
  maximum relative difference `5.32e-12`.
- New Schur example compiled under `-std=f2018`, checked with floating-point
  traps and with `-O3`, and ran with maximum errors at machine precision.
- Native failure paths cover true scalar resonances, unknown method strings,
  shape mismatches, empty matrices, nonfinite entries, the preserved direct
  size limit, and the new Schur order limit.

## Interpretation and limitations

The near-machine-precision results above apply to the tested, moderately
scaled families. They are **not** error bounds for arbitrary matrices.
The Schur implementation is intrinsic-only and uses an iterative QR algorithm
with a matrix-order cap of 40, not LAPACK. It may be slow or fail for difficult
nonnormal, clustered or nearly resonant cases; the diagonal denominators are
screened and the returned result undergoes a backward-residual check.
A small backward residual is not a certificate of small forward error for
an ill-conditioned Sylvester or Lyapunov operator. Full FPM integration tests
were not run because FPM was not available in this environment.
