# Deterministic simplicial global search (`scipy.optimize.shgo` subset)

SciFort provides `shgo` for bounded global minimization in 1-8 dimensions,
with distinct local-minimum discovery. The implementation is original and
independent. It follows the general simplicial-topology strategy of SHGO,
**but does not reproduce the full SciPy algorithm or its theoretical
convergence guarantees**.

## Example

```fortran
use scifort_kinds, only : dp
use scifort_optimize, only : shgo, shgo_result
implicit none
type(shgo_result) :: res
real(dp) :: bounds(2,2)
bounds(:,1) = [-6.0_dp, -6.0_dp]
bounds(:,2) = [ 6.0_dp,  6.0_dp]
res = shgo(himmelblau, bounds, iters=3, maxfun=15000)
! res%x / res%fun: best point / scalar objective
! res%xl(:,i) / res%funl(i): i-th local minimum, sorted by value
```

The full executable program with the Himmelblau callback is at
[`example/optimize_shgo_example.f90`](../example/optimize_shgo_example.f90).

## API

```fortran
res = shgo(fun, bounds, n, iters, maxfun, max_mesh_points, method, &
           local_maxiter, local_iter, minimize_every_iter, f_min, &
           f_tol, merge_tol, minimizer, callback, adaptive, adaptive_fraction, &
           hgr_stability_iters, hgr_min_iters, hgr_f_tol, conforming, &
           linear_a, linear_lb, linear_ub, nonlinear, nonlinear_lb, &
           nonlinear_ub, constraint_tol, nonlinear_jac, restore_feasibility, &
           restoration_starts, restoration_maxiter, constraints)
```

- `fun`: scalar `real(dp)` function of a real vector. `bounds` is `(dim,2)`
  with finite, strictly increasing bounds. Optional linear and nonlinear
  constraint blocks follow the bounds convention of SciFort's differential
  evolution solver, described below.
- `iters`: number of successively finer, nested simplicial meshes (default 1).
  Mesh width along each axis at level `k` is `2**k+1` vertices. The
  default first mesh has `3**dim` vertices; two levels require `5**dim`
  unique mesh evaluations, three levels `9**dim`, etc. Coarse vertices
  are cached and never re-evaluated.
- `n`: an optional **lower bound** on the number of vertices in the first
  complete dyadic mesh. It is *not* the exact sample count of SciPy's `n`.
  The smallest appropriate level at least 1 is selected.
- `maxfun`: strict maximum actual user objective calls (default 20,000).
  A full new grid is sampled only if enough calls remain to complete it;
  otherwise the existing best point and previously found minima are returned.
  The center is sampled first, even if the first full mesh cannot fit.
- `max_mesh_points`: maximum vertices in any one full uniform mesh, or
  maximum *distinct sampled vertices* in adaptive mode (default 100,000).
  Adaptive cells also have a bounded leaf-cell store.
- `adaptive`: enable selective local dyadic cell refinement (default
  `.false.`; existing uniform-mesh behavior is unchanged).
- `adaptive_fraction`: fraction of active leaf cells to split in each round,
  in the interval `(0,1]` (default `0.25`). At least one leaf is split per
  round; `1.0` refines every leaf, reproducing the complete uniform vertex
  grid. Each split replaces one cell by `2**dim` children. Candidate cells
  are ranked deterministically by a blend of their best and average finite
  corner values (87.5% minimum, 12.5% mean). This is an *exploitation*
  heuristic rather than SciPy's adaptive mesh algorithm.
- `conforming`: optional **SciFort extension** for `adaptive=.true.` in **2D or 3D**,
  default `.false.`. In 2D, each leaf has a fan triangulation from its sampled
  center to all boundary vertices. In 3D, each dyadic cube has a sampled center
  and six sampled face centers (shared across leaves); its faces are tiled by
  matching triangles, then coned to the cell center to make positively oriented
  tetrahedra. Exact lattice-coordinate checks validate nonzero volumes,
  coverage of all cell faces, and interior/exterior face multiplicities.
  Additional objective evaluations are counted within the same strict budget.
  Unsupported dimensions or `conforming=.true.` without `adaptive=.true.`
  return `shgo_invalid`. These complexes are geometrically conforming, but the
  adaptive basin-selection strategy still has no global convergence guarantee.
- `hgr_stability_iters`: **opt-in heuristic stopping**, default `0` (off).
  A positive value requires that many consecutive *completed* refinement
  rounds in which the count of graph local minima is unchanged **and** the
  best objective has not improved materially. For example, `2` requires two
  comparisons with the preceding round, therefore at least three completed
  rounds. A zero-count graph never qualifies.
- `hgr_min_iters`: minimum number of completed refinement rounds before HGR
  stopping can fire, default `3`, required `>= 2`. The initial observation
  never counts as a stable comparison.
- `hgr_f_tol`: nonnegative relative best-value improvement tolerance for the
  stabilization check, default `1e-8`. Improvement is considered small if
  `(old_best - new_best) / max(1, abs(old_best), abs(new_best))` does not exceed
  this tolerance. The normalization avoids overflow with large objectives.
  Set `0` to require no observed decrease, up to floating-point arithmetic.
- `method`: local-polishing algorithm, `nelder-mead` (default without
  constraints), `powell`, or **`slsqp` (default with constraints)**, all
  reused from SciFort. Constrained runs reject `nelder-mead` and `powell`
  unless a custom constrained `minimizer` is provided. `local_maxiter`
  defaults to 200.
- `local_iter`: at most this many promising vertices are passed to local
  minimization per mesh level, smallest sampled values first (default all).
  Zero disables all local minimization.
- `minimize_every_iter`: perform local minimization at each level (default
  `.true.`), or only at the final requested level. If a mesh/budget cap
  stops refinement early, the deferred final polishing does not run. When an
  HGR early-stop is about to fire, pending local polishing runs first; a
  material improvement cancels that tentative early stop.
- `f_min`: optional known global objective; reaching `f_min` within `f_tol`
  times `max(1,abs(f_min))` stops early. `f_tol` defaults to `1e-8`.
- `merge_tol`: maximum normalized coordinate difference for merging repeated
  local minima (default `1e-4`), applied per coordinate relative to its bound
  width.
- `minimizer`: optional custom local solver in place of SciFort's built-in
  solver. Its callback receives `(fun,lower,upper,start,maxfev,x,value,nfev,
  success)`; the passed objective is budget-guarded and checks bounds.
- `callback`: optional `(x,value,stop)` callback after each complete mesh
  or adaptive refinement round; set `stop=.true.` to return `shgo_callback_stop`.

The standalone `shgo_core` has the same options except `method` and
`local_maxiter`, and skips polishing when `minimizer` is omitted. It depends
only on `scifort_kinds` and Fortran standard intrinsic modules, so the core
can be compiled and tested without the optional library dependencies.

### Result

`shgo_result` contains `x(dim)`, `fun`, sorted `xl(dim,nminima)`,
`funl(nminima)`, `nfev`, `nlfev` (local-only evaluations), `ncev`
(constraint-vector evaluations), `maxcv` (largest constraint violation at
returned point), `feasible` (all constraints satisfied), `nit`
(complete levels/rounds), `nlocal` (polishing calls), `nminima`, `nmesh`
(last full uniform mesh size or number of cached adaptive vertices), `nleaf`
(active adaptive leaves), `nrefined` (cumulative cell subdivisions),
`ntriangles` (triangles in the latest fully validated conforming 2D mesh),
`ntetrahedra` (tetrahedra in the latest fully validated conforming 3D mesh),
`conforming_valid` (the current 2D/3D complex passed its exact lattice
geometry and simplex-incidence checks; false for incomplete rounds), `status`,
`success`, and `message`, plus `hgr_estimate` (number of graph minima on the
last completed refinement where candidates were examined; `-1` if none) and
`hgr_streak` (consecutive qualifying stability comparisons; `0` when the
feature is disabled). `nfev` counts *actual*
user objective calls, not reported counts from external local solvers.
`success` reports discovery of a finite feasible objective, not proof that
this is the global optimum.

Termination statuses:

| Constant | Value | Meaning |
|---|---:|---|
| `shgo_completed` | 0 | Requested levels completed |
| `shgo_maxfun` | 1 | Insufficient budget or excessive local calls |
| `shgo_target` | 2 | Known objective target reached |
| `shgo_callback_stop` | 3 | Callback requested termination |
| `shgo_mesh_limit` | 4 | Next grid/cell split exceeds mesh/leaf limit |
| `shgo_hgr_stable` | 5 | Sampled graph minima and best objective stabilized |
| `shgo_infeasible` | 6 | No feasible point was found within completed sampling |
| `shgo_invalid` | -1 | Invalid dimensions, bounds, or options |
| `shgo_numerical` | -2 | No finite feasible objective sample |

## Linear and nonlinear constraints

The bounded solver now accepts the same vector-constraint convention as
`differential_evolution`:

- `linear_a(m,nvar)`, `linear_lb(m)`, `linear_ub(m)` enforce
  `linear_lb <= matmul(linear_a,x) <= linear_ub`.
- `nonlinear(x,c)` is a user subroutine with
  `real(dp), intent(in) :: x(:)` and `real(dp), intent(out) :: c(:)`.
  Both `nonlinear_lb(m)` and `nonlinear_ub(m)` must accompany it.
- `constraints(:)` is an alternative to the single-vector `nonlinear=` API.
  It accepts one or more `shgo_nonlinear_constraint` records, each with its
  own callback, component-wise bound arrays, and optional Jacobian. See the
  dedicated example below. Do **not** pass `constraints=` together with any
  legacy `nonlinear`, `nonlinear_lb`, `nonlinear_ub`, or `nonlinear_jac`
  argument: this returns `shgo_invalid`.
- Bound components with **identical finite lower/upper limits** are
  equalities. Finite one-sided and two-sided bounds are inequalities.
  Use IEEE `-Inf`/`+Inf` for unused sides, *not* `-huge`/`huge`.
- `nonlinear_jac(x,a)` supplies analytic derivatives of the **untransformed**
  nonlinear constraint components. The callback receives `x(nvar)` and fills
  `a(m,nvar)`, with **one constraint component per row**. In the public
  `shgo` adapter, its rows are reordered and negated for upper bounds to
  match SciFort SLSQP's *equalities first, inequalities >= 0* convention.
  Linear constraint Jacobians use `linear_a` exactly. With nonlinear
  constraints and no `nonlinear_jac`, SLSQP retains its existing finite-
  difference Jacobian behavior. Supplying `nonlinear_jac` without a matching
  `nonlinear` function is invalid. `ncjev` counts user nonlinear Jacobian
  evaluations, separately from the constraint-vector count `ncev`.
- `restore_feasibility=.true.` enables a SciFort-only deterministic,
  **opt-in** feasibility projection from sampled vertices; default `.false.`
  preserves all prior sampling and polishing behavior. The core minimizes
  squared *active* equality/inequality residuals using damped Gauss-Newton
  steps in normalized coordinates and backtracking within the variable box.
  It uses `nonlinear_jac` if provided, otherwise bound-aware finite differences.
- `restoration_starts` (default 3) is the maximum projection starts per mesh
  level, including simplicial graph minima and additional cached vertices.
  Set 0 to disable the projection starts even if `restore_feasibility` is true.
  `restoration_maxiter` (default 25, positive) limits Gauss-Newton iterations
  per start. `nrestore` counts attempted starts; `nrestored` counts feasible
  points found with a verified finite **original objective** evaluation.
  A restored point is added to the result pool without claiming it is a
  locally optimal point; use SLSQP polishing for local refinement.
- Projection operations evaluate constraints but **never call the user
  objective at infeasible points**. Feasible trial points respect `maxfun`;
  new objective evaluations are included in `nfev`. Finite-difference
  constraint-vector calls increase `ncev`, not `nfev`. The same conservative
  uniform-grid budget preflight still applies when projection is enabled.
- `constraint_tol` (default `1e-8`) is the absolute per-component bound
  violation tolerated by SciFort when accepting a point as feasible.
  `maxcv` reports the actual largest violation, not a scaled score.

```fortran
use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf
real(dp) :: a(1,2), bounds(2,2), inf
inf = ieee_value(0.0_dp, ieee_positive_inf)
bounds(:,1) = 0.0_dp
bounds(:,2) = 1.0_dp
a(1,:) = [1.0_dp,1.0_dp]
! Require x(1)+x(2) >= 0.9. SLSQP is chosen automatically.
res = shgo(f, bounds, iters=2, maxfun=500, &
           linear_a=a, linear_lb=[0.9_dp], linear_ub=[inf])
if (res%success) print *, res%x, res%fun, res%maxcv
```

See [`example/optimize_shgo_constraints_example.f90`](../example/optimize_shgo_constraints_example.f90)
for an executable example mixing linear and nonlinear inequalities.

**Feasibility-first sampling.** Vertices are ranked by their maximum
constraint violation until a feasible vertex is found. A feasible vertex
receives a monotone, bounded objective score that always outranks an
infeasible vertex. Graph construction therefore continues to provide local
restoration starts even for equality surfaces that do not pass through any
sampled grid vertex. These scores are used **only** for mesh ranking;
`res%fun`, `res%funl`, and `res%xl` always contain original objectives and
feasible points. The solver skips expensive user objective calls at
infeasible mesh points, but constraint callbacks are counted in `ncev`.

The public `shgo` adapter sends equality and inequality components to
SciFort's existing SLSQP local optimizer. `shgo_core` remains independent
of FPM dependencies and can instead accept a custom constraint-aware
`minimizer`. Its callback receives the **original objective**, including
when the supplied point violates a nonlinear constraint, and must solve
or restore feasibility itself. The core independently verifies a custom
solver's returned point and rejects infeasible local results.

**Limitations:** Feasibility restoration is a local least-squares heuristic,
not a general constraint solver. It can fail on rank-deficient or zero-Jacobian
manifolds, incompatible constraints, narrow/disconnected feasible sets, poor
initial seeds, or finite-difference roundoff. It does not identify *all*
connected components, certify feasibility/infeasibility, or guarantee an
optimum. The selected starts are deterministic and may revisit a basin;
`nrestored` is not a count of topologically distinct components. The objective
at each feasible restored point is independently evaluated. Finding a feasible
point is not guaranteed, especially for thin/disconnected constraints or
problems with poor local convergence.
A nonzero `maxcv` up to `constraint_tol` is accepted. The original complete-
mesh budget preflight remains conservative even when infeasible samples
skip objective calls. An infeasible exit is not a mathematical proof that
no feasible point exists. The heuristic graph-minimum count with
constraints is based on a feasibility-first ranking and is no longer a
pure objective sublevel-set H0 birth count; HGR stabilization should be
used cautiously. SLSQP can evaluate the user objective at infeasible
points during local restoration. SciFort's geometrically conforming 2D/3D
meshes work with constraints but do not gain a convergence guarantee.

The native core, example interface, and public adapter with a test-only
SLSQP stand-in are exercised in the standalone validation script. Full
linking against SciFort's real SLSQP dependency still requires `fpm test`.

## Multiple independent nonlinear constraint objects

Use the public `shgo_nonlinear_constraint` type from `scifort_optimize`
(or directly from `scifort_optimize_shgo`). The type has four fields:

```fortran
type(shgo_nonlinear_constraint) :: blocks(2)
blocks(1)%fun => sum_constraint
blocks(1)%jac => sum_jacobian      ! Optional exact Jacobian
blocks(1)%lb = [0.7_dp]
blocks(1)%ub = [0.7_dp]
blocks(2)%fun => difference_constraint
! blocks(2)%jac deliberately absent: only its rows use finite differences.
blocks(2)%lb = [0.1_dp]
blocks(2)%ub = [0.1_dp]
res = shgo(f,bounds,iters=1,constraints=blocks, &
           restore_feasibility=.true.)
```

Each `fun(x,c)` is a user subroutine, and each `jac(x,a)` (when associated)
returns a Jacobian of shape `(size(lb),size(x))`, with one component per row.
The `lb` and `ub` fields are allocatable `real(dp)` arrays and must be
allocated and have equal lengths. Array assignment such as
`blocks(1)%lb = [0.7_dp]` allocates them automatically.

- Procedure components are **pointers** to user callbacks. Use `=>` to
  associate them, and ensure their procedure targets remain valid until the
  `shgo` call returns. A missing `fun`, invalid/NaN bound, wrong-sign
  infinite limit, reversed bound pair, or unallocated bounds is an input
  error. `shgo_constraint_blocks_valid(blocks)` checks record metadata;
  it cannot verify the sizes or correctness of user callback output.
- Component arrays are flattened in the order of the input blocks.
  SLSQP converts the resulting components to its canonical ordering:
  all finite **equalities first**, followed by each finite lower inequality
  and negated upper inequality. Linear rows appear before nonlinear rows
  within each section. A vector block may freely mix equalities and one-
  or two-sided inequalities.
- Each block's analytic `jac` is evaluated only when needed. A block without
  `jac` uses bound-aware one-sided finite differences; analytic rows from
  other blocks are **retained**, rather than discarded. Restoration and
  the public local-polisher adapter support this mixed derivative mode.
  `ncjev` counts nonlinear Jacobian callback invocations. `ncev` counts
  constraint checks and additional finite-difference constraint evaluations
  (not necessarily individual callbacks in an aggregated check).
- `constraints=` may be empty (no nonlinear constraints), and may be used
  alongside `linear_a` / `linear_lb` / `linear_ub`. The original single-
  callback `nonlinear=` and `nonlinear_jac=` arguments still work, with their
  previous defaults. A block with zero components is permitted but is
  skipped during evaluations.

For a complete program with a two-component vector block and a mixture of
analytic and numerical Jacobians, see
[`example/optimize_shgo_multiconstraints_example.f90`](../example/optimize_shgo_multiconstraints_example.f90).

**Scope:** These constraint objects represent SciPy's general lower/upper
bound semantics but do not implement the complete `NonlinearConstraint`
class (`hess`, `keep_feasible`, sparse Jacobians and method configuration are
not supported). The core's feasibility projection is a local heuristic and
cannot prove global constrained optimality. The public SLSQP adapter has
been checked with a derivative-validating test stub; real SLSQP integration
requires a full FPM build with external numerical dependencies.

## Numerical method and limitations

On each dyadic rectangular grid, adjacent cells are triangulated by the
Kuhn-Freudenthal monotone simplex construction. Two vertices share a graph
edge when they differ by one grid step on a nonempty subset of axes, all in
one coordinatewise direction. We order all vertices by objective and use
a deterministic tie break. A graph local minimum creates a new component
in the 0-dimensional sublevel-set filtration. These component births provide
candidates for bounded local minimization; distinct local solutions are
retained and sorted. The global best sampled or locally evaluated value is returned.
The best incumbent is included in the returned candidate-minimum pool, even
if a custom local minimizer visits it but later reports failure.

Uniform grids remain the default. With `adaptive=.true.`, the solver starts
from the same complete grid, then selectively bisects promising *leaf cells*.
Each leaf cell is represented by the edge graph of its Kuhn simplices;
vertex coordinates are keyed exactly on a common integer dyadic lattice,
so each distinct point is sampled only once. Each split is committed only if
its entire `(3**dim)` subgrid fits in the remaining objective/mesh budgets.
A round may stop after some but not all selected cells have been refined:
`nmesh` and `nleaf` reflect all *committed* splits, whereas `nit` counts only
fully completed rounds. Local polishing uses the same strict objective budget
and can prevent later refinement.

In adaptive mode, vertices with exactly equal objective values are grouped
into connected **plateaus** using a disjoint-set (union-find) pass over the
active leaf graph. A plateau supplies just one candidate, and only when none
of its vertices has a strictly downhill neighbor. This avoids spurious
minima caused by insertion-order tie breaking as coarse edges disappear.
Nearly equal values are not merged.

### Homology-growth stopping (SciFort extension)

The `hgr_estimate` is the number of local-minimum births in the sampled
zero-dimensional sublevel-set graph. It is **not** the final zeroth Betti
number, not a persistence diagram, and not the exact homology group rank
estimated by SciPy's `minhgrd` option. This implementation does **not**
reproduce that SciPy stopping rule. Requiring a stable count for successive
rounds, a minimum number of observations, and no material improvement is a
deliberately conservative *heuristic*, not proof of global optimality.

```fortran
res = shgo(f, bounds, iters=6, hgr_stability_iters=2, &
           hgr_min_iters=3, hgr_f_tol=1.0e-8_dp)
if (res%status == shgo_hgr_stable) then
    print '(a,i0)', 'stable sampled graph-minimum count: ', res%hgr_estimate
end if
```

The count may remain unchanged even while a new unsampled basin exists.
In nonconforming adaptive mode, missing inter-cell graph edges add further
uncertainty. With `conforming=.true.`, these hanging-node gaps are repaired,
but the sampled-minimum count is still only a heuristic H_0 birth count. For difficult multimodal functions, increase the stability
and minimum-round requirements or disable the stopping option altogether.
Material sampling or local improvement resets the streak. Both modes remain
deterministic.

**Adaptive topology:** the default adaptive corner-edge graph has hanging
vertices at coarse-fine boundaries and is **not** conforming. In 2D,
`conforming=.true.` repairs the topology through explicit fan triangles.
Every cell center is sampled once, and each cell boundary is subdivided at
all cached hanging vertices. Triangles have positive exact lattice area;
total area and all interior/exterior edge incidences must match.

In **3D**, each leaf cube has its center and all six face-center samples in
the exact vertex cache. Each face is subdivided into dyadic subfaces wherever
all requisite edge midpoints are present. Subfaces are triangulated by a fan
about a sampled face center when available, or exact-integer convex polygon
ear clipping otherwise. Tetrahedra are formed by coning each face triangle
to the owning cube center. The implementation validates positive exact
lattice-volume determinants, exact face-area coverage, and triangular face
incidences: each interior triangle belongs to exactly two tetrahedra,
while a single-incidence face must be on an exterior domain boundary.
These checks detect coarse-fine face cracks and inverted tetrahedra. Failed
checks yield `shgo_numerical`; a safety limit of 300,000 tetrahedra (also
bounded by 24 times `max_mesh_points`) yields `shgo_mesh_limit`.

This provides genuine geometric simplicial complexes **in 2D and 3D**, but
it does not make adaptive cell selection globally exhaustive or establish the
convergence properties of SciPy SHGO. In particular, low-scoring-cell selection
can neglect unsampled basins. For wider exploration, raise
`adaptive_fraction` or use uniform sampling. This is not SciPy's adaptive
SHGO algorithm. The exact checks validate mesh connectivity and coverage,
not objective smoothness or optimizer convergence.

The solver does not implement conforming adaptive triangulation in 4D or higher,
higher-dimensional triangulation of scattered points, persistent homology,
SciPy-equivalent homology-rank stopping, analytic constraint Jacobians,
independent lists of constraint objects, symmetry reduction, or Sobol/Halton sampling. Nor does it implement SciPy's
`minimizer_kwargs`, `options`, `workers`, or arbitrary `sampling_method`.
Nonfinite objective values are treated as unusable; a nonfinite local result
is discarded. Mesh construction still grows exponentially in dimension (for
uniform mode, roughly `O((2**level+1)**dim * 2**dim)`), so low-dimensional
problems are the main use case. Neither mode certifies global optimality on a
finite sample mesh.

## Adaptive example

```fortran
res = shgo(himmelblau, bounds, iters=5, maxfun=5000, &
           adaptive=.true., adaptive_fraction=0.25_dp, &
           local_iter=3, method='nelder-mead')
! res%nleaf: current number of leaf cells
! res%nrefined: total number of cell subdivisions
```

See [`example/optimize_shgo_adaptive_example.f90`](../example/optimize_shgo_adaptive_example.f90)
for an executable uniform-versus-adaptive comparison. `iters` counts the
initial full grid as round 1 in both modes. The new adaptive controls are
added at the *end* of the existing argument list for source compatibility.
The three HGR controls are appended after those adaptive controls, likewise
preserving the existing calling conventions. The `conforming` option is
appended last for backward source compatibility.

## Conforming 2D example

```fortran
res = shgo(himmelblau, bounds, iters=4, adaptive=.true., &
           conforming=.true., adaptive_fraction=0.25_dp, &
           maxfun=5000, method='nelder-mead')
print *, res%ntriangles, res%conforming_valid
```

See [`example/optimize_shgo_conforming_example.f90`](../example/optimize_shgo_conforming_example.f90).
When all leaves are refined at each round, the 2D conforming mode uses
additional center evaluations: a 3x3 initial grid plus four cell centers
requires 13 calls, not nine. At the next uniform-equivalent level, 41
vertices generate 64 triangles. All samples and each complete split obey
the same strict `maxfun` and `max_mesh_points` caps.

## Conforming 3D example

```fortran
res = shgo(f, bounds3, iters=4, adaptive=.true., &
           conforming=.true., adaptive_fraction=0.25_dp, &
           local_iter=0)
print *, res%ntetrahedra, res%conforming_valid
```

See [`example/optimize_shgo_conforming_3d_example.f90`](../example/optimize_shgo_conforming_3d_example.f90).
With a `3 x 3 x 3` starting grid, the initial conforming complex has **71**
unique objective samples (27 corners, 36 shared face centers, 8 cell centers)
and **192** tetrahedra. A single selected octree split adds 56 previously
uncached samples and produces a validated 15-leaf/414-tetrahedron complex.
Sampling is **not** the SciPy SHGO sampling method; the 3D construction is a
SciFort extension. Face centers and tetrahedra add cost compared with the
nonconforming corner-only adaptive graph. Neither local minima detection nor
stabilization is a persistent-homology certificate.

## Provenance and validation

Based on the independent mathematical description in S. C. Endres,
C. Sandrock, and W. W. Focke (2018), *A simplicial homology algorithm for
Lipschitz optimisation*, Journal of Global Optimization 72, 181-217,
[doi:10.1007/s10898-018-0645-y](https://doi.org/10.1007/s10898-018-0645-y).
SciPy reference API: <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.shgo.html>.
No upstream SHGO/SciPy implementation code was copied. See
[`CODE_PROVENANCE.md`](../CODE_PROVENANCE.md) for details.

Run stand-alone native and SciPy comparisons with:

```sh
bash scripts/validate_optimize_shgo.sh
```

The full SciFort test suite, including the bounded Nelder-Mead/Powell
adapter, requires FPM and the project's existing dependencies:

```sh
fpm test
```

Detailed results: [`validation/optimize_shgo_2026-11-10.md`](../validation/optimize_shgo_2026-11-10.md),
[`validation/optimize_shgo_adaptive_2026-11-10.md`](../validation/optimize_shgo_adaptive_2026-11-10.md),
[`validation/optimize_shgo_hgr_2026-11-10.md`](../validation/optimize_shgo_hgr_2026-11-10.md),
[`validation/optimize_shgo_conforming_2026-11-10.md`](../validation/optimize_shgo_conforming_2026-11-10.md),
[`validation/optimize_shgo_conforming_3d_2026-11-10.md`](../validation/optimize_shgo_conforming_3d_2026-11-10.md),
and [`validation/optimize_shgo_constraints_2026-11-10.md`](../validation/optimize_shgo_constraints_2026-11-10.md).

## Additional verification: analytic Jacobians and feasibility restoration

Run `bash scripts/validate_optimize_shgo.sh`. The new test
`test/optimize_shgo_restoration_reference.f90` checks exact and numerical
constraint Jacobians, equality manifolds absent from the mesh, nonlinear
curves/products, inconsistent bounds, invalid callback results, the 3D
conforming mesh, and objective budgets. The test-only SLSQP stand-in now
checks the canonical Jacobian against independent finite differences. See
[the separate validation record](../validation/optimize_shgo_restoration_2026-11-10.md)
and the [standalone example](../example/optimize_shgo_restoration_example.f90).

**Full FPM integration against the actual SciFort SLSQP solver has not been
run in this environment.** The standalone SLSQP stub verifies interfaces and
derivative layout, not optimizer convergence or external dependency linkage.
