# SciFort `linprog` validation — 2026-11-10

## Environment and repeatability

- GNU Fortran 14.2.0, Python 3, NumPy 2.3.5, SciPy 1.17.0.
- Checked compilation: `-std=f2008 -Wall -Wextra -O0 -g -fcheck=all -fbacktrace`.
- Optimized compilation: `-std=f2008 -Wall -Wextra -O3`.
- Reproduce both builds with `python3 scripts/check_optimize_linprog.py`.

## Verified results

| Suite | Checked build | Optimized build |
|---|---:|---:|
| Native Fortran assertions | 27/27 | 27/27 |
| Independent SciPy/HiGHS LP cases | 1,450/1,450 | 1,450/1,450 |
| Example (`x=[3,1]`, objective `-5`) | Pass | Pass |

The independent reference suite includes 1,181 optimal, 169 unbounded,
and 100 infeasible cases per build. It exercises mixed and unbounded bounds,
free variables, overdetermined/redundant equalities, degenerate vertices,
fixed variables, and independently scaled constraints/objectives. In the
successful SciPy cases, the largest relative objective difference was
`4.34984e-13` and the largest normalized primal constraint violation was
`1.95399e-14`. The native tests include invalid inputs, a zero-variable
LP, and large-RHS feasibility/infeasibility checks.

The reference harness sets `presolve=False` in HiGHS because with presolve
*enabled* one generated problem that is mathematically feasible and unbounded
was classified infeasible; HiGHS with presolve disabled and SciFort both
classified it unbounded. This discrepancy is not counted as a solver failure;
the case and corresponding settings are reproducible in the Python harness.

## Scope and limitations

This is a first independent two-phase **tableau** simplex implementation,
primarily aimed at small/moderate dense continuous LPs. It lacks sparse
matrix storage, an industrial revised simplex, presolve, warm starts,
`integrality`, sensitivity marginals and advanced numerical scaling.
It is not established that it is robust or fast for large/ill-conditioned
industrial LPs. The complete project-wide `fpm test` suite and full
third-party dependency graph could not be executed because FPM is not
installed in this environment. Only the explicitly listed standalone
native and SciPy comparisons were run. No prior numerical source modules
were altered; the public optimization facade gained new exports.
