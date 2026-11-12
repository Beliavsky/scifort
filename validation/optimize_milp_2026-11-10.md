# Dense MILP validation — 2026-11-10

Implementation: `src/optimize/milp.f90` (best-bound branch-and-bound using
SciFort `linprog` LP relaxations, MIT). Full FPM integration not run because
`fpm` is unavailable in the validation environment.

## Standalone GNU Fortran 14.2.0 tests

- `test/test_optimize_milp.f90`: **78 assertions** in a runtime-checked build,
  including SciPy's published two-integer example, both mixed and pure
  integer objectives, variable types 0/1/2/3, equality and inequality rows,
  exact node budgets, requested gaps, invalid inputs, and large RHS values.
- The same native test passed under `-O3 -march=native`.
- `validation/compare_optimize_milp_scipy.py --count 320`: **320 problems,
  640 checks** passed against SciPy 1.17.0/HiGHS. The same number passed with
  `--optimize`. Models mix continuous/integer/semi-continuous/semi-integer
  variables and random one-/two-sided and equality constraints. Fortran
  results undergo separate variable-domain, row-feasibility, and objective
  recomputation checks.
- `validation/compare_optimize_milp_bruteforce.py`: **220 independent pure
  integer models, 348 assertions** passed. All bounded assignments were
  enumerated to verify optima or infeasibility independently of SciPy.
- `example/optimize_milp_example.f90`: published SciPy solution
  `x=[1,2], fun=-2`, and semi-continuous solution `x=3` reproduced.

## Reused simplex fix and unchanged optimization regressions

A previously existing LP Phase I ratio comparison treated values near
`1e12` that differ by one as ties because it scaled the user tolerance
by the RHS ratio. This selected an incorrect leaving row and misclassified
feasible equality-constrained LPs. `src/optimize/linprog.f90` now uses a
machine-precision-scale ratio tie threshold, with a regression test.

- `scripts/check_optimize_linprog.py`: **29 native checks** and **1,450
  SciPy comparison cases** passed in each of checked and optimized builds.
- Existing `quadratic_assignment` native suite: **469 assertions** passed in
  checked and optimized builds, and **5,156** optimized cross-check
  assertions passed separately against SciPy.
- Existing standalone SHGO quality driver: **33 SciPy comparisons** passed.

## Scope and cautions

The Fortran solver is dense and repeatedly solves an LP tableau at every
branch-and-bound node. It has no sparse support, presolve, cuts, parallelism,
or primal heuristics. General mixed-integer LP relaxation unboundedness is
reported as unresolved rather than proven mixed-integer unboundedness.
CPU time limits are checked between nodes. No performance parity with
SciPy's production HiGHS engine is claimed. Full FPM integration and
cross-platform/compiler testing remain future work.
