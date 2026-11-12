# SHGO conforming adaptive 3D tetrahedral mesh validation (2026-11-10)

## Method and contract

In three dimensions, `adaptive=.true., conforming=.true.` builds a mesh of
independent dyadic leaf cubes. Each leaf center and every distinct leaf-face
center is sampled using a common exact integer-coordinate key cache. A dyadic
face subdivision is used wherever all its four edge midpoints and the face
center have been sampled. All boundary hanging vertices are included in the
triangulation; a patch is fanned about its sampled face center or triangulated
by exact-integer ear clipping. Each oriented face triangle is coned to its
leaf center to form a nondegenerate tetrahedron. Geometry checks verify:

- Each tetrahedron has positive **exact 64-bit integer lattice** determinant.
- The signed two-dimensional triangle areas sum to each cube face's area.
- Each triangle face in the tetrahedral complex occurs exactly twice, except
  outer-boundary triangles, which occur once **and lie on the true box face**.
- All samples (including face centers) are cached and obey atomic per-cell
  `maxfun` / `max_mesh_points` preflights. Incomplete selected rounds invalidate
  topology diagnostics rather than advertising the last valid triangulation.
- The simplex count has a safety cap of `min(300000,24*max_mesh_points)` to avoid
  excessive memory use when collecting global incidence information.

Geometric conformity does **not** establish global optimization convergence,
persistent homology or SciPy's `minhgrd` behavior. Dimensions >=4 still use
the original nonconforming adaptive graph unless `conforming` is requested,
in which case the arguments are rejected.

## Reproduction

```sh
bash scripts/validate_optimize_shgo.sh
```

The script compiles the standalone core in both checked and `-O3` GNU Fortran
modes, runs **47 new 3D native assertions per mode**, compares five independent
benchmarks with SciPy SHGO 1.17.0, and checks the public wrapper using a test-only
minimizer stub. It also reruns all previous SHGO core suites. `fpm` and the
complete repository dependencies were not available for an integrated build.

## Quality comparisons (no local polishing in SciFort driver)

| Function | Vertices | Leaf cubes | Tetrahedra | SciFort sampled best | SciPy local-polished best |
|---|---:|---:|---:|---:|---:|
| Shifted quadratic | 1,145 | 176 | 5,208 | 0.0009375 | ~0 |
| Rastrigin (3D) | 457 | 64 | 1,788 | 0 | 1.98992 |
| Ackley (3D) | 457 | 64 | 1,788 | ~0 | 2.95828 |
| Double well | 1,199 | 176 | 5,460 | 0.00625 | ~0 |
| Rosenbrock (3D) | 1,145 | 176 | 5,160 | 0.694458 | ~0 |

Both suites passed in checked and optimized builds. These values are sampled
minima without SciFort local polishing, while SciPy SHGO uses Sobol sampling
and local minimization; thus the comparison checks objective quality within
case-specific tolerances, **not** identical algorithm behavior or fair equal
function-evaluation budgets. A SciPy run can itself become trapped in local
minima (as in the Rastrigin and Ackley rows).

## Regression highlights

- Initial `3 x 3 x 3` structured grid: 71 cached vertices, 8 leaves,
  192 validated tetrahedra.
- One selected split: 127 vertices, 15 leaves, 414 tetrahedra.
- Full one-level split of all eight leaves: 429 vertices, 64 leaves,
  1,536 tetrahedra.
- Seven deeply uneven rounds: 6 splits, 50 leaves; all faces validate.
- A split requiring 56 new evaluations is rejected at `maxfun=126`, accepted
  at `maxfun=127`. A partially completed *batch* clears `conforming_valid`
  and `ntetrahedra` until the next complete validation.
- Flat objective gives one minimum plateau; two disconnected low-value regions
  are correctly recognized as two minima. Early stopping and callback tests pass.

## Limitations

Only 3D octree cells are supported by this tetrahedralization; general
polyhedral meshes, non-box constraints, 4D+ conforming complexes, persistence
and certified global convergence are outside its scope. The safety cap on
tetrahedra may reject a problem before its point budget is exhausted.
