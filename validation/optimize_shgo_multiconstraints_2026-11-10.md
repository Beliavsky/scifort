# SHGO independent nonlinear constraint objects — 2026-11-10

## Scope

This update adds `constraints(:)` to the standalone `shgo_core` and public
`shgo` interfaces, with type `shgo_nonlinear_constraint` containing independent
procedure pointers (`fun`, optional `jac`) and separate allocatable `lb` and
`ub` arrays. It supports mixed equality/inequality components, multiple
vector-valued callbacks, mixed analytic and numerical Jacobians, optional
linear constraints, and the existing opt-in feasibility restoration.

The legacy `nonlinear=`/`nonlinear_jac=` interface remains unchanged and
cannot be mixed with `constraints=`. No new solver dependency was introduced.

## Reproduction

```bash
bash scripts/validate_optimize_shgo.sh
bash scripts/validate_optimize_direct.sh
```

Validated with GNU Fortran 14.2.0 under both
`-std=f2008 -O0 -g -Wall -Wextra -fcheck=all
-ffpe-trap=invalid,zero,overflow` and `-std=f2008 -O3 -Wall -Wextra`.
SciPy 1.17.0 was used for independent numerical comparisons.

| Suite | Checked | Optimized |
|---|---:|---:|
| Original SHGO native assertions | 61 | 61 |
| Adaptive SHGO native assertions | 46 | 46 |
| Homology-growth native assertions | 36 | 36 |
| Conforming 2D native assertions | 47 | 47 |
| Conforming 3D native assertions | 47 | 47 |
| Previous constrained SHGO assertions | 58 | 58 |
| Feasibility-restoration native assertions | 45 | 45 |
| **New independent-block assertions** | **40** | **40** |
| Existing SciPy comparison cases | 64 | 64 |
| **New SciPy/SLSQP numerical comparison cases** | **5** | **5** |
| Public adapter / derivative-layout smoke test | passed | passed |
| Example (mixed exact/FD derivatives) | passed | passed |
| DIRECT native checks | 109 | 109 |
| DIRECT SciPy comparison cases | 60 | 60 |

Five new problems independently verify objective values at each reported
point and re-evaluate all constraint components. Two have mixed analytic and
finite-difference derivatives; one uses a vector-valued third block with
one-sided bounds; one is three-dimensional. Constraints are checked to within
`1.1e-8` absolute violation. Four cases use SciPy's SHGO as an objective-value
reference and one (a curved equality) uses SciPy SLSQP because SHGO may report
nominal success for a point violating that equality. A 3D case uses Sobol
sampling on the SciPy side for the same reason. The comparator independently
checks SciPy feasibility rather than trusting `success` alone.

The native assertions cover invalid callback pointers, unallocated or
mismatched bound arrays, NaNs, wrong-sign infinities, mixed old/new API
rejection, contradictory equalities, objective-call limits, combined linear
and multiple nonlinear constraints, cached 3D conforming meshes, and separate
analytic vs. finite-difference Jacobian work.

The public adapter is checked against a deliberately restricted *test-only*
SLSQP interface that validates equality/inequality ordering, sign
transformations, and mixed analytic/numerical Jacobian entries via finite
differences. This is **not an SLSQP numerical implementation**. `fpm` is not
installed in the test environment, so the real SciFort SLSQP dependency was
not linked and full `fpm test` has not been run. Actual constrained local
polishing against that dependency remains to be verified.

## Limitations

- The `shgo_nonlinear_constraint` records borrow callback procedure targets;
  they must remain valid until the optimization call returns.
- User-supplied Jacobian dimensions and derivative correctness are not
  automatically verified at runtime in release mode; the adapter smoke test
  exercises the canonical layout on representative examples.
- Constraint callbacks are sequential. Results use the existing `ncev` / 
  `ncjev` accounting; `ncev` counts batches and extra finite-difference
  evaluations and is not necessarily the sum of all per-block callback calls.
- Multiple independently supplied blocks do not change the underlying SHGO
  sampling policy, and the local feasibility projection still cannot prove
  mathematical feasibility/infeasibility or global optimality.
- This is not the entire SciPy `NonlinearConstraint` class; Hessians, sparse
  Jacobians, `keep_feasible`, alternative local methods, and custom
  difference-step options are not exposed.
