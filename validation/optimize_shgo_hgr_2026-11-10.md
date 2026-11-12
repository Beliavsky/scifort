# SHGO sampled homology-growth stabilization — validation (2026-11-10)

## Scope

An opt-in `hgr_stability_iters` option stops when the number of finite sampled
simplicial-graph minima stays unchanged for a requested number of successive
**completed rounds**. A new `hgr_min_iters` guard defaults to 3 completed
rounds; a nonnegative `hgr_f_tol` guard defaults to `1e-8` normalized objective
improvement. An all-nonfinite graph never qualifies. The result exposes
`hgr_estimate` and `hgr_streak`, and `shgo_hgr_stable = 5` is a distinct
termination status. Existing defaults are unchanged.

Adaptive equal-value plateaus are grouped by union-find over current leaf
edges before checking downhill neighbors. Thus constant surfaces produce a
single connected plateau rather than artificial minima caused by removed
coarse-grid edges and historical insertion indices. Distinct flat basins
remain distinct. The adaptive graph is still nonconforming at hanging nodes.

## Test command and results

```sh
bash scripts/validate_optimize_shgo.sh
```

Run on GNU Fortran 14.2.0 and SciPy 1.17.0. The two compile configurations
are `-std=f2008 -O0 -g -Wall -Wextra -fcheck=all
-ffpe-trap=invalid,zero,overflow` and `-std=f2008 -O3 -Wall -Wextra`.

| Check | Checked | Optimized |
|---|---:|---:|
| Existing SHGO native assertions | 61 | 61 |
| Existing adaptive native assertions | 46 | 46 |
| New HGR stopping and plateau assertions | 36 | 36 |
| Existing independent SciPy benchmark cases | 33 | 33 |
| Existing adaptive comparisons with SciPy | 6 | 6 |
| New HGR-stopped objective comparisons with SciPy | 4 | 4 |
| Public local-adapter smoke test (with test-only minimizer stub) | passed | passed |

The new native tests cover: opt-in behavior, the minimum observation count,
longer stability windows, 1D two-well discovery, material improvement resets,
full and selective adaptive meshes, connected flat plateaus and multiple
separate flat basins, deferred local solver calls, callbacks, strict budgets,
all-nonfinite inputs, target priority, reproducibility, and invalid options.

Four additional Fortran-driver/SciPy quality comparisons independently
recompute the returned objective for a double well, a shifted quadratic, a
2D sphere, and a 2D double basin. The SciPy comparison is against its
unmodified `shgo` search, not its `minhgrd` stopping rule; only objective
quality, feasibility and HGR early-termination bookkeeping are checked.

| Benchmark | HGR evaluations | Five-level full grid | HGR objective | SciPy objective |
|---|---:|---:|---:|---:|
| Double well, 1D | 17 | 33 | 0 | 0 |
| Shifted quadratic, 1D | 17 | 33 | 0.0025 | approximately 0 |
| Sphere, 2D | 81 | 1,089 | 0 | 0 |
| Two basins, 2D | 289 | 1,089 | 0 | 0 |

These illustrate that HGR early termination can save evaluations but may
leave a worse sampled objective than SciPy's locally polished solution.

### Example: constant 1D objective, 5 nominal levels

- Without HGR stopping: all 5 levels are attempted (33 mesh evaluations).
- With `hgr_stability_iters=2` and default minimum rounds: stop after level
  3 (9 evaluations) with `hgr_estimate=1`, `hgr_streak=2`.
- These counts are expected from the dyadic mesh construction and verified
  by the native test; they do **not** establish performance on other functions.

### Important limitations

This is a discrete graph-local-minima **count** (births in an H_0 sublevel
filtration), not the final H_0 rank or a persistent-homology computation.
It is not an exact reproduction of SciPy's `minhgrd` criterion, and equal
counts can conceal different basins or missed minima. The additional best-
value condition reduces but does not eliminate premature-stop risk. The
adaptive triangulation is nonconforming and has no proven SHGO guarantee.

The full `fpm test` run was **not available** in this environment. The SHGO
adapter was compiled and exercised against the repository's test-only stub,
not the real Nelder-Mead/Powell modules. This establishes interface forwarding
and does not establish real FPM integration. No new external dependencies were
introduced and no third-party implementation source was copied.
