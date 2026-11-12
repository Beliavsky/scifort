# `scipy.optimize.linprog`-style dense linear programming

`src/optimize/linprog.f90` implements a **self-contained, continuous,
finite-dimensional, dense** linear-programming solver. It is a first
algorithmic implementation, not a drop-in replacement for HiGHS.

## Problem and API

Minimize `dot_product(c,x)` subject to

- `matmul(a_ub,x) <= b_ub`
- `matmul(a_eq,x) == b_eq`
- `lower <= x <= upper`

Call `result = linprog(c, a_ub, b_ub, a_eq, b_eq, lower, upper, tol, maxiter)`
from either `scifort_optimize_linprog` or the `scifort_optimize` facade.
The constraint matrices and corresponding right-hand sides must be supplied
in pairs, but both kinds of constraints are optional. Unspecified bounds
are **zero below** and **unbounded above**, as in SciPy; use IEEE negative
infinity for an unrestricted lower bound and IEEE positive infinity for an
unrestricted upper bound. Arrays use `real(dp)`.

`type(linprog_result)` has fields:

- `x`: solution or last feasible iterate for valid input (not allocated if
  Phase I diagnoses infeasibility or validation fails)
- `fun`: `dot_product(c, x)` when `x` is allocated
- `slack`: inequality residual `b_ub - A_ub*x`
- `con`: equality residual `b_eq - A_eq*x`
- `success`, `nit`, `status`, and `message`

Status codes follow SciPy: `linprog_optimal=0`, `linprog_iteration_limit=1`,
`linprog_infeasible=2`, `linprog_unbounded=3`, `linprog_numerical=4`.
SciFort adds `linprog_invalid=-1` for malformed inputs. Check `success`
before using the returned point as an optimizer. Pivots (`nit`) need not
match SciPy/HiGHS.

## Algorithm and limitations

This is an independent dense **two-phase primal tableau simplex** solver.
The transform `x=x0+T*y`, `y>=0`, handles finite bounds, upper-only variables,
and unrestricted variables. Upper bounds become extra inequalities.
Artificial variables supply a Phase I feasible basis for equality and
reversed inequality rows. Redundant equalities are removed, and Bland's
smallest-index entering/leaving rule avoids combinatorial cycling.

It is intended for **small to moderate dense LPs**. It does not implement
HiGHS's revised simplex, interior-point method, presolve, sparse constraints,
integer programming (`integrality`), dual marginals, or solver warm starts.
`tol` controls tableau pivot and feasibility decisions, not HiGHS's separate
primal and dual feasibility tolerances. Large or ill-conditioned problems
may require a robust revised-simplex/interior-point backend; results should
be independently checked in important applications.

The implementation does **not** use or vendor third-party LP code. We
considered [DualSimplex](https://github.com/thchang/DualSimplex), MIT licensed,
but it specializes in the asymmetric dual of dense `A*x<=b` problems and
requires external LAPACK/BLAS operations and basis handling. Incorporating
it directly would not provide the complete `linprog` input formulation without
substantial adaptation. It remains an optional future backend candidate.

## Example

```fortran
use scifort_kinds, only : dp
use scifort_optimize, only : linprog, linprog_result
! ...
type(linprog_result) :: res
real(dp) :: a(2,2)
a(1,:) = [1.0_dp,1.0_dp]
a(2,:) = [1.0_dp,3.0_dp]
res = linprog([-1.0_dp,-2.0_dp],a_ub=a,b_ub=[4.0_dp,6.0_dp])
! res%x=[3,1]; res%fun=-5
```

See `example/optimize_linprog_example.f90` for a complete runnable example.

## Reproducible tests

From the SciFort root, run

```
python3 scripts/check_optimize_linprog.py
```

This runs GNU Fortran checked and optimized builds, native tests, an example,
and an independent SciPy/HiGHS comparison harness. It does not run the full
`fpm test` suite; that must be performed separately in an environment with
FPM and the project's pinned dependencies available.
