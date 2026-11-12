# Clough–Tocher C¹ interpolation on scattered 2D data

SciFort provides `scifort_interpolate_clough_tocher` and exports the constructors
from `scifort_interpolate`. The new routines implement the SciPy 1.17-compatible
`CloughTocher2DInterpolator` surface, plus `griddata(method='cubic')` through
the public `scifort_interpolate` facade.

```fortran
use scifort_kinds, only: dp
use scifort_interpolate, only: clough_tocher_2d_interpolator, &
    clough_tocher_2d_interpolator_type, griddata

real(dp) :: sites(4,2), values(4), query(2,2)
real(dp), allocatable :: answer(:,:), result(:)
type(clough_tocher_2d_interpolator_type) :: ct
integer :: status

sites(1,:) = [0.0_dp, 0.0_dp]
sites(2,:) = [1.0_dp, 0.0_dp]
sites(3,:) = [0.0_dp, 1.0_dp]
sites(4,:) = [1.0_dp, 1.0_dp]
values = [0.0_dp, 1.0_dp, 1.0_dp, 2.0_dp]
query(1,:) = [0.25_dp, 0.40_dp]
query(2,:) = [0.75_dp, 0.20_dp]

ct = clough_tocher_2d_interpolator(sites, values, &
    tol=1.0e-10_dp, maxiter=800, status=status)
if (status /= 0) error stop "Clough-Tocher initialization failed"
answer = ct%evaluate(query)
result = griddata(sites, values, query, method='cubic')
```

`clough_tocher_2d_interpolator(x,y, rescale,fill_value,tol,maxiter,status)`
accepts `x(nsites,2)`, `y(nsites)` or `y(nsites,nchannels)`, and returns a
**real** or **complex** interpolation object selected by the ordinate type.
`%evaluate(x)` returns an allocatable vector of channel values; `%evaluate(xi)`
returns `(nquery,nchannels)`. `%value(x,channel=1)` returns a scalar.
`%simplices()` returns the 1-based triangulation indices; `%valid()`,
`%converged()`, and `%channels()` expose object state.

`griddata` supports the existing `nearest` and 2D `linear` methods unchanged;
`cubic` is now available for 2D coordinates. The low-level
`scifort_interpolate_scattered:griddata` interface continues to provide
nearest/linear only; import the public facade for cubic.

## Numerical method

1. Bowyer–Watson constructs the triangular mesh, optionally using a wider
   floating-point kind for circumcircle predicates. A boundary-notch repair
   and local Delaunay edge flips complete/regularize the mesh. An independent
   convex-hull area check rejects any unresolved gaps. When the compiler
   cannot represent a 33-digit kind, the predicates revert to `real(dp)`;
   the geometry checks remain enabled.
2. A global Gauss–Seidel solve estimates a gradient at each vertex from the
   minimum-curvature network, one independent system per real or complex
   output channel. Convergence is measured by a relative/absolute tolerance.
3. The interpolant uses four subtriangles per triangle with affine-invariant
   cubic Bernstein–Bézier controls. Neighboring triangle centroids determine
   the common cross-edge derivative direction, producing a C¹ surface.

The formula and algorithm follow the original SciPy Cython implementation,
with attribution and required BSD-3-Clause notices retained at
[`LICENSES/scipy-interpnd-BSD-3-Clause.txt`](../LICENSES/scipy-interpnd-BSD-3-Clause.txt).
SciPy uses Qhull, whereas SciFort uses the project's native independent
Bowyer–Watson mesh; triangle choices can differ for cocircular sites.

**Limits:** exactly two spatial dimensions and unique finite sites are
required. A strongly elongated or nearly collinear cloud can fail native
triangulation; `rescale=.true.` often helps. A failed mesh is rejected, not
silently treated as a full convex hull. Robust adaptive geometric predicates
and higher-dimensional cubic scattered interpolation are not implemented.
Construction currently uses O(n²) Bowyer–Watson work. **Cubic Bernstein–Bézier
controls are cached** for every triangle and channel; queries no longer
recompute the vertex-gradient or neighboring-centroid formulas. The cache
stores 19 complex control values per channel per triangle.

A reusable 2D uniform-bin spatial index maps grid cells to candidate triangle
bounding boxes. It stores candidate lists in compressed sparse row (CSR)
format, with triangle IDs in their original order so that shared-edge
lookup is deterministic. The number of bins grows with the number of
triangles, and a memory cap of 32 references per triangle triggers an exact
linear-search fallback for elongated/pathological meshes. Lookup is typically
much faster than scanning every triangle, but worst-case search is still
O(ntri). No external spatial-index dependency is required.

To preserve queries on shared edges when the physical coordinates have large
offsets relative to their ranges, bounding-box expansion and barycentric
containment tolerances account for floating-point roundoff introduced by
coordinate normalization. This cannot eliminate the inherent conditioning
problems of extremely thin or degenerate triangles.

Extremely ill-conditioned Delaunay predicates may still fail or choose a
slightly different mesh from Qhull. The optional extended precision does not
require any additional declared FPM dependency, but its compiler-specific
support and linkage should be verified on Windows/MinGW-w64 and Intel Fortran.

The `status` return codes are 0 success, 1 invalid arguments/duplicates,
2 geometrically degenerate or incomplete triangulation, 3 gradient estimate
did not converge in `maxiter` iterations. With code 3, the object remains
valid and can still evaluate, but `%converged()` is false; the caller should
increase `maxiter` or `tol` before relying on its results. The default
tolerance is `1.0e-6` and maximum iterations 400, matching SciPy.
Outside-hull queries return NaNs unless `fill_value` is specified.

## Reproducible testing

```bash
python tools/interpolate_reference/compare_clough_tocher_scipy.py --cases 160
python tools/interpolate_reference/compare_clough_tocher_scipy.py --cases 160 --seed 8675309
python tools/interpolate_reference/run_all_native_interpolation.py
```

The reference generates SciPy Qhull interior queries plus site and outside-hull
queries, and compares complex multichannel values for real/complex datasets,
including rescaled anisotropic clouds. Native regressions cover exact affine
reproduction, shared-edge C¹ continuity, invalid inputs, fill behavior, and
public `griddata` dispatch. The test suite cannot establish Qhull-like
robustness for nearly cocircular/degenerate point configurations.

No LAPACK/BLAS, Qhull, network access, or mandatory new FPM dependency is
required for the native code. To use full SciFort FPM builds, the existing
`fpm.toml` dependencies are unchanged.

## Indexed evaluation benchmark (2026-10-31)

The benchmark compiles the same deterministic Fortran driver against the
previous ZIP and the updated source, without downloading anything:

```bash
python tools/interpolate_reference/benchmark_clough_tocher.py \
    --baseline-zip /path/to/scifort-interpolate-clough-tocher-cubic-2026-10-31.zip
```

Using GNU Fortran 14.2 `-O3` on one Linux environment, for 450 sites,
883 triangles, four channels, and 9,000 queries repeated four times:

| Phase | Previous | Cached/indexed |
|---|---:|---:|
| Construction | 0.1794 s | 0.1862 s |
| Mean evaluation time | 0.1800 s | 0.00374 s |

The measured evaluation speedup was about **48.2×**; the numerical checksums
matched. These are single-run CPU timings, **not guaranteed speedups** and
exclude construction, for which caching introduces some additional work.

New native regression `test/test_interpolate_clough_tocher_cache.f90` tests
210 anisotropically distributed sites, 403 triangles, 1,824 queries, real
and complex multichannel values, every site's ordinate, triangle centroids,
shared-edge midpoints, outside-hull fills, and a nonconverged-gradient
object. The existing SciPy comparison suite continues to exercise general
randomized data. Both builds (checked and `-O3`) are tested with
`tools/interpolate_reference/run_all_native_interpolation.py`.
