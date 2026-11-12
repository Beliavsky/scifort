# Validation: bounded linear least squares (2026-11-10)

- GNU Fortran 14.2.0, checked `-O0 -fcheck=all` and optimized `-O3 -fcheck=all`.
- Each build passed 1,200 reproducibly seeded real NumPy/SciPy comparison
  cases (2,400 total). Cases include nonnegative, finite box and mixed effectively
  unbounded boxes; tall, wide, rank-deficient, and zero matrices.
- All 2,400 native solves returned feasible, first-order-optimal results with
  objective values no worse than the SciPy 1.17 reference within the required
  test tolerances. No SciPy reference exceptions occurred.
- Worst relative objective gap above the SciPy reference: `3.315e-14`.
- Largest scaled first-order KKT violation: `6.598e-15`.
- Standalone native assertion suite passed 32 checks in each build.
- Previous differential-evolution native regressions: 104 plus 42 checks per
  build, passed in checked and optimized builds.
- SLSQP BVLS bridge: compiled and ran against a source-compatible test
  double only; **actual pinned SLSQP runtime validation remains pending**.
- Full `fpm test` not executed (FPM unavailable in this environment).

Comparison driver and reproducibility script:
`tools/optimize_linear_least_squares_compare.f90` and
`tools/compare_optimize_linear_least_squares_scipy.py`.
