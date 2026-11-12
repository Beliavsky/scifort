# Adaptive SHGO extension — implementation and validation (2026-11-10)

## Changes

- Opt-in `adaptive=.true.` and `adaptive_fraction=0.25_dp` support, appended
  to the `shgo_core` and public `shgo` argument lists for source compatibility.
- Deterministic leaf-cell scoring (87.5% best finite corner, 12.5% average),
  exact integer-lattice vertex hashes, and selective dyadic subdivision.
- The candidate graph uses Kuhn/Freudenthal edges on each active leaf.
  In contrast to uniform mode, the *adaptive graph is not conforming* at
  hanging vertices. No topological convergence guarantee is claimed.
- Per-cell atomic preflight against both the objective-call and mesh-vertex
  budget. A partially completed *round* retains completely subdivided cells,
  but does not increment `nit`. Diagnostic `nmesh`, `nleaf`, `nrefined` report
  the actual committed mesh. Existing budget-guarded local polishing,
  objective accounting, callbacks, target stopping and deduplication remain.
- No third-party algorithms or source code were copied, and no new runtime
  dependencies were added.

## Reproducible validation

```sh
bash scripts/validate_optimize_shgo.sh
```

| Check per build | GNU Fortran `-O0 -fcheck=all` | GNU Fortran `-O3` |
|---|---:|---:|
| Existing SHGO native assertions | 61 passed | 61 passed |
| **New adaptive native assertions** | **46 passed** | **46 passed** |
| Existing independent SciPy SHGO cases | 33 passed | 33 passed |
| Existing level 2/3 objective checks | 22 passed | 22 passed |
| **New adaptive SciPy quality comparisons** | **6 passed** | **6 passed** |
| Public SHGO adapter signature smoke test | passed | passed |

Compiler: GNU Fortran 14.2.0. SciPy reference: 1.17.0. Checked flags:
`-std=f2008 -O0 -g -Wall -Wextra -fcheck=all
-ffpe-trap=invalid,zero,overflow`; optimized flags:
`-std=f2008 -O3 -Wall -Wextra`.

The new native tests cover 1D/2D/3D graphs; `adaptive_fraction=1` equivalence
with uniform grids; exact reuse of mesh vertices; meaningful reduction of
objective evaluations; determinism; invalid parameters; finite/infinite
objective values; initial/incremental mesh and objective budgets; local search;
callback cancellation; and immediate target stopping. The independent SciPy
script recomputes values at Fortran-returned points and compares objective
quality against SciPy SHGO refinements. It does **not** claim the same sample
order or numerical convergence guarantee.

### Sampling cost and objective quality (without local polishing)

| Benchmark | Adaptive evaluations | Complete grid of same nominal depth | Adaptive objective | SciPy locally refined objective |
|---|---:|---:|---:|---:|
| Shifted quadratic, 2D | 23 | 81 | 0.02125 | ~0 |
| Himmelblau, 2D | 131 | 1,089 | 1.62134 | ~0 |
| Rastrigin, 2D | 39 | 289 | 0 | 0 |
| Ackley, 2D | 23 | 81 | ~0 | ~0 |
| Tilted double-well, 1D | 5 | 9 | -0.1 | -0.10062 |
| Shifted quadratic, 3D | 321 | 4,913 | 0.021875 | ~0 |

The comparison is illustrative, not a matched-effort performance benchmark:
SciPy applies local minimization by default, whereas the figures above are
**sample-only** results. The adaptive heuristic is not guaranteed to give
lower objective values than uniform sampling at every budget. For instance,
it may underexplore an isolated basin whose initial corners look poor.

## Integration limitation

`src/optimize/shgo_minimize.f90` compiles and passes a public API smoke test
against the **test-only** stub module in
`validation/optimize_shgo_test_minimize_stub.f90`. This verifies Fortran
signature compatibility and keyword forwarding, not real numerical behavior
of the Nelder-Mead or Powell backends. The actual solver implementation and
full repository test suite still require an FPM installation with the
SciFort dependencies:

```sh
fpm test
```

Future work: conforming adaptive simplicial refinement; rigorous
persistent-homology stopping beyond the subsequently added graph-minimum
stabilization heuristic; linear and nonlinear constraints; better exploration
versus exploitation; and FPM-backed end-to-end testing.
