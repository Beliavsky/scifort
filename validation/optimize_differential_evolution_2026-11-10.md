# Differential evolution validation — 2026-11-10

The new `scifort_optimize_differential_evolution` module was compiled and
run with **GNU Fortran 14.2.0**, **SciPy 1.17.0**, and **NumPy 2.3.5**.

## Commands

From the SciFort root directory:

```text
python tools/run_differential_evolution_tests.py
```

This runner compiles the native Fortran test and the Fortran reference driver
from source in temporary directories, with separate checked and optimized
builds. The checked build uses `-O0 -fcheck=all -fbacktrace`, and the optimized
build uses `-O3 -march=native`; both use `-std=f2018 -Wall -Wextra -Werror`
for native tests.

## Results

| Test | Checked | Optimized |
|---|---:|---:|
| Native Fortran assertions | 104 passed | 104 passed |
| Independent SciPy reference cases | 72 passed | 72 passed |
| Near-global objective `< 1e-3` | 72/72 | 72/72 |
| Equivalent SciPy runs near global optimum | 71/72 | 71/72 |
| Fortran outcomes as good or better than SciPy, tolerance `1e-3` | 72/72 | 72/72 |
| Maximum absolute error in recomputed objective | `3.55e-15` | `3.55e-15` |

Independent comparisons test all 12 standard mutation/crossover strategies
on shifted quadratic, Rosenbrock, Rastrigin, Ackley, boundary-constrained,
and fixed-coordinate quadratic objectives. Scipy's pseudorandom generator
differs, so the comparisons verify independent objective values, bound
satisfaction, known optimality, and solution quality rather than requiring
identical stochastic samples, output vectors, or iteration counts.

Native tests additionally cover fixed/all-fixed bounds, single variables,
Latin hypercube stratum coverage, deterministic seeds, accepted/rejected
starting vectors, objective NaNs, invalid parameters, forced evaluation
budgets including exhaustion immediately after initialization, callback
termination, population shape, and objective/evaluation counter accounting.

Example: `example/optimize_differential_evolution_example.f90` builds as a
standalone program against `src/kinds.f90` and the new module, and reached
objective 0.0 on a two-dimensional Rastrigin instance with seed 2026.

## Limitations

The complete FPM dependency-linked test suite was **not** run because FPM
is unavailable in this environment. Facade export through `src/optimize.f90`
was code-reviewed, but not compiled through FPM here; the new module and its
independent example and tests were compiled directly. No new FPM dependencies
were added. The implementation lacks constraints, integer variables,
polishing, QMC initializations, vectorized evaluation, and parallel workers.
