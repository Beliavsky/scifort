# SHGO conforming adaptive 2D mesh validation (2026-11-10)

## Scope

The opt-in `conforming=.true.` API with `adaptive=.true.` is supported for two
variables. A leaf is triangulated from its sampled center to its complete
boundary chain, including all hanging dyadic nodes. Validation verifies:

- Every triangle has strictly positive area in exact integer lattice arithmetic.
- Sum of doubled triangle areas equals twice the normalized domain area.
- Interior mesh edges appear in two triangles; edges on the domain boundary
  appear once; no segment is used more than twice.
- Strict `maxfun` and `max_mesh_points` include interior samples and forbid
  partially committed splits.
- Reproducibility, nested refinement, flat and disconnected plateau minima,
  local polishing, target/callback stopping, nonfinite functions, and errors.
- Existing uniform, nonconforming adaptive and HGR mode behavior unchanged.

The geometric tests do not prove global optimization convergence. Mesh
conformity is not available for three or more variables.

## Reproduction

```bash
bash scripts/validate_optimize_shgo.sh
```

This tests the core under GNU Fortran 14.2.0 with `-fcheck=all` and `-O3`.
The suite additionally compares five 2D objectives with SciPy 1.17.0:
shifted quadratic, Himmelblau, Rastrigin, Ackley, and six-hump camel.
Results are quality comparisons against SciPy's locally refined solution;
they do not assert identical samples, graph structure, or stopping rules.

## Test results

Both checked and optimized builds:

| Test group | Passing checks per build |
|---|---:|
| Original SHGO core | 61 native assertions |
| Nonconforming adaptive SHGO | 46 native assertions |
| Homology-growth stopping | 36 native assertions |
| **Conforming adaptive SHGO** | **47 native assertions** |
| Original SciPy quality cases | 33 comparisons |
| Nonconforming adaptive SciPy quality cases | 6 comparisons |
| HGR SciPy quality cases | 4 comparisons |
| **Conforming SciPy quality cases** | **5 comparisons** |
| Public adapter smoke tests | Passed with test-only local solver stub |

Example conforming-only cases:

| Objective | Distinct vertices | Leaf cells | Triangles | Sampled best |
|---|---:|---:|---:|---:|
| Shifted quadratic | 36 | 13 | 59 | 0.02125 |
| Himmelblau | 228 | 97 | 439 | 0.1304016 |
| Rastrigin | 64 | 25 | 110 | 0.0 |
| Ackley | 36 | 13 | 58 | ~0.0 |
| Six-hump camel | 234 | 97 | 452 | -1.0264681 |

All five conforming meshes passed exact edge/area checks. Native assertions
also verify the small 7-leaf/30-triangle nonuniform mesh, full refinement
with 16 leaves/64 triangles, and deeper asymmetric refinement.

## Limits

A complete end-to-end build against the actual SciFort FPM dependencies
remains unverified because FPM was not available in the test environment.
The adapter's signature/option forwarding and example were checked against
a test-only stand-in for the local minimizer module.

A partially completed round may have valid individual cell splits, but its
triangulation has not been rebuilt or validated; `conforming_valid` is false
and `ntriangles` is reset to zero in that case.
