# SHGO analytic constraint Jacobians and feasibility restoration — 2026-11-10

## Change scope

This update extends SciFort's original bounded simplicial SHGO variant with:

- An optional analytic `nonlinear_jac(x,a)` (shape `(m,nvar)`) accepted by
  `shgo_core` and the public `shgo` adapter. Linear rows use the exact matrix
  supplied by the caller; combined SLSQP rows are organized as equalities,
  lower inequalities, and negated upper inequalities. When a nonlinear block
  lacks an analytic Jacobian, SLSQP uses its existing numerical derivative.
- **Opt-in** damped projected Gauss-Newton restoration of sampled starts to
  linear or nonlinear equality/inequality manifolds. For nonlinear components
  it uses the user analytic Jacobian or bounded numerical differences.
  A graph-minimum/vertex seed list helps avoid zero-Jacobian plateau starts.
- Strict independent objective checks before accepting restored solutions.
  Invalid/nonfinite evaluations and impossible constraints are not reported
  as feasible. Objective evaluations contribute to `nfev`; the two separate
  constraint work counters are `ncev` and `ncjev`.
- New `nrestore` (attempted seeds) and `nrestored` (successful finite feasible
  evaluations) fields. The original untransformed objective is always returned.

All new options have unchanged defaults (`restore_feasibility=.false.`);
existing constrained and unconstrained SHGO calls preserve their prior path.

## Reproduction

```bash
bash scripts/validate_optimize_shgo.sh
bash scripts/validate_optimize_direct.sh
```

Validated using **GNU Fortran 14.2.0**, with `-std=f2008 -O0 -g -Wall
-Wextra -fcheck=all -ffpe-trap=invalid,zero,overflow` and separately
`-std=f2008 -O3 -Wall -Wextra`. SciPy **1.17.0** was used for comparisons.
Both scripts completed successfully in both build modes.

| Suite | Checked | Optimized |
|---|---:|---:|
| SHGO original native | 61 | 61 |
| SHGO adaptive native | 46 | 46 |
| SHGO homology-growth native | 36 | 36 |
| SHGO conforming 2D native | 47 | 47 |
| SHGO conforming 3D native | 47 | 47 |
| SHGO constrained native | 58 | 58 |
| **New SHGO restoration native** | **45** | **45** |
| Prior SciPy SHGO comparisons | 59 | 59 |
| **New SciPy SLSQP feasibility comparisons** | **5** | **5** |
| Public adapter test (test-only SLSQP stub) | pass | pass |
| DIRECT native | 109 | 109 |
| DIRECT SciPy comparisons | 60 | 60 |

The new comparisons test a linear equality, nonlinear sum equality, nonlinear
squared-radius equality with analytic and numeric derivatives, and a product
equality. The driver independently verifies that every returned objective is
actually evaluated at its returned feasible point, compares it with SciPy
SLSQP numerical constrained references as an *objective lower bound*, and
rejects an impossible equality. This **does not** claim that projection alone
is an optimizer or reproduces SciPy SHGO objective minima.

The derivative-checking local-optimizer stub calls SciFort's canonical
constraint values and Jacobian and compares the rows against independent
finite differences. It checks mixed linear equality and nonlinear two-sided
constraints, including the sign of upper-bound derivatives. It intentionally
does not emulate SLSQP convergence.

## Known limitations and next verification

- A Gauss-Newton projection is a local feasibility heuristic, **not** a
  certificate of feasibility/infeasibility, global minimum, or connected
  component discovery. Rank-deficient, nearly-singular, disconnected, or
  extremely narrow regions can defeat it. A feasible projected point is not
  necessarily a stationary local minimum.
- `constraint_tol` is an absolute per-component tolerance. Near-zero residual
  roundoff must be interpreted relative to that tolerance; very large or
  badly scaled constraint components can require user rescaling.
- The existing full uniform-grid budget preflight remains conservative even
  though infeasible grid vertices do not consume objective evaluations.
- **Full `fpm test` and integration against the actual FPM-provided SLSQP
  module were not possible in this environment: `fpm` is not installed.**
  The public SHGO adapter was compiled against the repository's explicitly
  labeled derivative-checking test stub only. This verifies interface shape,
  argument routing and Jacobian-row conversion, but not real SLSQP linking
  or numerical convergence.
