# Dense mixed-integer linear programming (`scipy.optimize.milp` subset)

SciFort now provides `milp` in `scifort_optimize_milp`, re-exported by the
`scifort_optimize` facade. The implementation is an **independent dense
branch-and-bound solver**, **not a HiGHS interface**, intended for small and
moderate problems. It reuses SciFort's two-phase simplex `linprog` module.

## Problem and Fortran API

Minimize `dot_product(c, x)` subject to

```text
b_l <= matmul(a, x) <= b_u
lower <= x <= upper
```

with the selected variables constrained to integer or semi-continuous domains.
Use a negated coefficient vector for maximization.

```fortran
res = milp(c, integrality=itype, a=a, b_l=bl, b_u=bu, &
    lower=lower, upper=upper, node_limit=10000, &
    mip_rel_gap=0.0_dp, tol=1.0e-9_dp)
```

Inputs are `real(dp)` except `integrality` and the integer limits. The
`integrality(n)` values follow SciPy's convention:

| Value | Variable domain |
|---|---|
| 0 | Continuous within its bounds (default) |
| 1 | Integer within its bounds |
| 2 | Either zero **or** any value in its original `[lower,upper]` range |
| 3 | Either zero **or** an integer in its original `[lower,upper]` range |

Unspecified variable lower bounds are zero and upper bounds are `+infinity`.
Use `ieee_value(1.0_dp, ieee_positive_inf)` for infinity, or its negative.
The `a(m,n)`, `b_l(m)`, and `b_u(m)` arrays must be supplied *together*, or
all omitted. Finite equal row bounds indicate equality constraints; infinite
endpoints leave one side unconstrained. A zero-row matrix is permitted.
Individual nonlinear constraints, sparse matrices, and multiple independent
constraint objects are not included in this interface.

## Result fields

`type(milp_result)` contains:

- `x(n)` and `fun`: best known **feasible** incumbent and its objective. `x`
  is unallocated if no feasible point has been found. `fun` is `HUGE` in that
  case rather than a meaningful objective.
- `status`: 0 optimal or requested relative gap achieved; 1 node, simplex, or
  CPU-time budget reached; 2 infeasible; 3 proven continuous unboundedness;
  4 numerical failure or unresolved integer-relaxation unboundedness;
  -1 invalid input.
- `success`: true for status 0 only.
- `mip_node_count`: number of actual LP relaxations solved.
- `nit`: total simplex pivots in all relaxations.
- `mip_dual_bound`: inherited valid lower bound on any unresolved branch,
  or the primal objective once the entire search is proven optimal.
- `mip_gap`: `(fun - dual_bound) / max(1,abs(fun))` when both primal and dual
  bounds are available; otherwise `HUGE`.
- `message`: explanatory termination message.

Options: `node_limit` defaults to 10000 LP relaxations, `mip_rel_gap` to
`1.0e-4`, and `tol` to `1.0e-9`. `lp_maxiter` sets a per-node simplex pivot
budget. `time_limit` is measured in **CPU seconds** and tested *between*
node solves, unlike SciPy's wall-clock HiGHS time limit. For an exact search,
set `mip_rel_gap=0.0_dp`. Node and time limits can leave a valid incumbent
without proving it globally optimal; inspect `success`, `mip_gap`, and
`mip_dual_bound`.

## Method and limitations

The root problem drops integrality and relaxes the disjoint semi-variable
domains into their convex intervals. Each node solves a continuous LP using
the existing two-phase simplex implementation (with a corrected large-RHS
minimum-ratio pivot tie check). Nonintegral variables split
at floor/ceiling; semi-variables lying outside `{0} U [lower,upper]` split
into a zero branch and a nonzero interval branch. The solver chooses the
unexplored node with the best inherited LP lower bound, verifies primal
feasibility before recording an incumbent, and uses certified parent bounds
to prune. Relative gap and evaluation-budget stopping are optional.

This is a **small-problem correctness-oriented baseline**, not a replacement
for SciPy's mature sparse HiGHS MILP solver. It has no presolve, cutting
planes, advanced branching heuristics, primal rounding heuristics, sparse
matrix support, or parallelism. The number of LP solves may grow
exponentially. Pure-integer unbounded LP relaxations yield a conservative
`milp_numerical` status because LP unboundedness alone does not prove integer
unboundedness. Floating-point feasibility checks are tolerance-based and
extremely ill-conditioned problems may be unresolved.

## Example and reproducible tests

The runnable Fortran example is `example/optimize_milp_example.f90`. Using
the published SciPy example, it returns the integer solution `x=[1,2]`
with `fun=-2`. Build the standalone native test using the sources
`src/kinds.f90`, `src/optimize/linprog.f90`, and `src/optimize/milp.f90`.

```shell
python validation/compare_optimize_milp_scipy.py
python validation/compare_optimize_milp_scipy.py --optimize
python validation/compare_optimize_milp_bruteforce.py
```

The independent SciPy comparison covers randomized continuous, integer,
semi-continuous, and semi-integer problems and validates the returned solution
against every bound and row constraint. The exhaustive reference checks
small bounded integer problems by enumerating all integer assignments.
Neither validation script requires FPM. The full repository-wide FPM build
and integration suite remain to be run with the actual dependency set.

SciPy reference: <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.milp.html>.

The MILP implementation also motivated a targeted numerical bug fix in
`linprog`: the Phase I leaving-row ratio comparison no longer treats
order-one differences between order-`1e12` ratios as ties. A dedicated
large-RHS equality regression was added to `test_optimize_linprog.f90`.
