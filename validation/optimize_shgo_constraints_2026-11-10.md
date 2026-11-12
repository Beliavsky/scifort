# Constrained SHGO verification — 2026-11-10

## Scope

The optional `linear_a/linear_lb/linear_ub` and
`nonlinear/nonlinear_lb/nonlinear_ub` blocks accept exact equality components,
lower/upper inequalities and mixed forms. Every returned solution is
independently checked against the bound constraints; infeasible local
minimizer outputs are rejected even when the callback reports success.

The core uses feasibility-first graph rankings so that strict equality
surfaces without sampled feasible vertices can still generate local-polishing
starts. Its `nfev` counter counts only actual objective calls, whereas
`ncev` counts constraint-vector checks (including local SLSQP evaluations
through the public adapter). Infeasible mesh vertices do not call the user's
objective. The returned `fun` and `funl` contain **untransformed** objective
values, never internal graph scores.

## Reproduction

```bash
bash scripts/validate_optimize_shgo.sh
```

With GNU Fortran 14.2.0, in both checked mode (`-O0`, runtime bounds checks,
IEEE invalid/zero/overflow traps) and optimized mode (`-O3`):

| Suite | Passed in each mode |
|---|---:|
| Original SHGO | 61 assertions |
| Adaptive SHGO | 46 assertions |
| Homology-growth stopping | 36 assertions |
| Conforming 2D | 47 assertions |
| Conforming 3D | 47 assertions |
| **New constraint tests** | **58 assertions** |
| Previous SciPy SHGO quality comparisons | 53 cases |
| **New constrained analytic/reference problems** | **6 cases** |
| **New successful SciPy constrained SHGO comparisons** | **6 cases** |
| Public adapter signature smoke test | Passed |

The constrained tests include no-feasible-sample equalities, feasible
restoration using a custom constrained polisher, mixed linear and nonlinear
bounds, impossible constraints, inaccurate, fabricated-objective, nonfinite, or infeasible minimizer results,
strict objective-call caps, plateau/connectivity interactions and malformed
constraint input including NaN bounds.

## Independent numerical reference cases

The six new comparisons cover:

1. A quadratic with a one-sided linear inequality; optimum 0.18.
2. A quadratic with an offset linear equality; optimum 0.
3. A quadratic with a nonlinear disk inequality; analytic boundary optimum.
4. A quadratic with a nonlinear circle equality; analytic boundary optimum.
5. Combined equality and nonlinear upper inequality; optimum 0.
6. Three-dimensional conforming sampling with a linear equality; optimum 0.

All six SciPy 1.17.0 `scipy.optimize.shgo` runs succeeded and agreed with
SciFort's independently verified constrained results to within the stated
tolerances. The Fortran comparison driver supplies **known-solution test
polishers** to isolate objective/constraint accounting and candidate
acceptance. It does **not** validate convergence of SciFort's actual SLSQP
adapter, and is not evidence of identical numerical SHGO algorithms.

## Integration boundary

`src/optimize/shgo_minimize.f90` forwards constrained problems to
`scifort_optimize_minimize::slsqp`. The adapter compiles against a
**test-only stand-in** that invokes equality/inequality callbacks with the
same interface. The actual SLSQP dependency was not installed here, and
`fpm test` against the complete SciFort dependency graph could not be run.
This is the remaining integration validation item; the tests should be
repeated on a machine with the declared FPM dependencies installed.

## Limitations

The method is not a complete reproduction of SciPy SHGO's constrained
simplicial homology algorithm, nor a proof of global optimality. Its graph
ranking under constraints is a heuristic and cannot be interpreted as the
ordinary objective-sublevel H0 filtration. Equality manifolds that are
very thin, curved, disconnected or outside local solver basins may still
fail. `constraint_tol` may permit a small residual, which is visible in
`maxcv`; `shgo_infeasible` means **no feasible point found by this run**,
not a mathematical infeasibility certificate. Exact vertex evaluation
budget preflights remain conservative for constrained runs.
