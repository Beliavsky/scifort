# Deterministic SHGO subset: implementation and validation (2026-11-10)

## Scope

- A new independent `scifort_optimize_shgo` core generates nested dyadic
  Kuhn-Freudenthal simplex graphs, identifies candidate local minima as births
  of sublevel-set connected components, and caches all prior mesh evaluations.
- A new public `shgo` adapter connects the global core to SciFort's existing
  bounded Nelder-Mead and Powell solvers. Both procedures export result and
  callback types through `scifort_optimize` without changing existing APIs.
- Strict objective-call and mesh-size budgets, optional target stopping,
  deterministic results, duplicate-minimum consolidation and callback-based
  termination are supported.
- No third-party source was imported. Only existing project dependencies are
  needed by the public local-search adapter.

## Reproducible tests

```sh
bash scripts/validate_optimize_shgo.sh
```

| Test | Checked GNU Fortran | Optimized GNU Fortran |
|---|---:|---:|
| Standalone native assertions | 61 passed | 61 passed |
| Independent SciPy SHGO cases | 33 passed | 33 passed |
| Level 2/3 SciPy objective comparisons | 22 passed | 22 passed |
| Maximum normalized objective recomputation error | 0 | 0 |
| Existing DIRECT native assertions | 109 passed | 109 passed |
| Existing DIRECT SciPy comparisons | 60 passed | 60 passed |

**Compiler:** GNU Fortran 14.2.0. Checked flags: `-std=f2008 -O0 -g
-Wall -Wextra -fcheck=all -ffpe-trap=invalid,zero,overflow`. Optimized
flags: `-std=f2008 -O3 -Wall -Wextra`. SciPy reference: 1.17.0.

The 33 SciPy comparisons use 11 functions at three mesh depths: shifted
quadratics (1-D to 3-D), rotated quadratic, Rastrigin, Ackley, Himmelblau,
tilted double-well, oscillatory multiwell, absolute-value minimum, and a
piecewise two-basin objective. They verify returned points and independently
recomputed objective values, compare level-2/3 objective quality with SciPy,
check bounds, budget accounting, and known benchmark targets. The standalone
reference driver injects an independent bounded coordinate-pattern local
solver to test the global mesh engine without the optional SciFort local
solver dependencies.

Native tests include nonfinite functions, invalid and overflowing bounds,
strict budgets and local-solver budget overruns, deterministic repeatability,
H_0 candidate counts, duplicate removal, custom local-search callbacks, failed
local-minimizer returns, target stopping, and early callback termination.

For additional regression coverage, `scripts/validate_optimize_direct.sh`
was run successfully in checked and optimized modes.

The `src/optimize/shgo_minimize.f90` adapter also compiled and passed a
minimal Fortran wiring smoke test against a temporary stub matching the
existing SciFort optimizer interfaces. **This is a syntax/interface check,
not a numerical validation of the real Nelder-Mead/Powell adapter.**
The required Fortran Package Manager (`fpm`) is not installed in this
execution environment, and its external Git dependencies were not fetched.
Consequently, a full `fpm test`/public-adapter numerical build was **not**
run. It remains a required integration test on an FPM-equipped system.

## Limitations and next priorities

- Structured *uniform* rather than adaptive triangulation; mesh growth is
  exponential in dimension. The public API is bounded to 1-8 dimensions.
- Only H_0 sublevel-set component births are used; no homology persistence
  computation, rank-differential stopping, or published finite-convergence
  guarantee is asserted for this restricted solver.
- No linear/nonlinear constraints, arbitrary sampling methods, or parallel
  objective execution.
- `n` sets a minimum **first full mesh** size, not an exact SciPy sample
  count; evaluation histories and local minima may differ from SciPy.
- Next useful additions: adaptive triangulation based on edge indicators,
  constrained feasible-complex generation, and homology-rank stopping,
  followed by full FPM regression on Windows and Linux.
