# Scattered-data nearest-neighbor and dimension-generic simplex interpolation

The new `scifort_interpolate_scattered` module provides native Fortran
implementations of three `scipy.interpolate` interfaces:

| Fortran constructor | SciPy counterpart | Supported dimensions |
|---|---|---|
| `nearest_nd_interpolator` | `NearestNDInterpolator` | Any positive number of coordinate dimensions |
| `linear_nd_interpolator` | `LinearNDInterpolator` | 2D–5D native Delaunay; 2D or higher with supplied simplices |
| `griddata` | `griddata(..., method='nearest'/'linear')` | Any positive dimension for `nearest`; 2D–5D for `linear` |

The API is exported both from `scifort_interpolate` and the module directly.
No third-party library, system BLAS/LAPACK, Qhull, or extra FPM dependency
is required.

## Basic usage

```fortran
use scifort_kinds, only: dp
use scifort_interpolate, only: linear_nd_interpolator, &
    linear_nd_interpolator_type, nearest_nd_interpolator, &
    nearest_nd_interpolator_type, griddata

real(dp) :: points(4,2), values(4), query(2,2)
real(dp), allocatable :: answer(:), batch(:,:)
type(linear_nd_interpolator_type) :: linear
integer :: status

points(1,:) = [0.0_dp, 0.0_dp]
points(2,:) = [1.0_dp, 0.0_dp]
points(3,:) = [0.0_dp, 1.0_dp]
points(4,:) = [1.0_dp, 1.0_dp]
values = [0.0_dp, 1.0_dp, 1.0_dp, 2.0_dp]
query(1,:) = [0.25_dp, 0.4_dp]
query(2,:) = [0.75_dp, 0.2_dp]

linear = linear_nd_interpolator(points, values, &
    fill_value=-1.0_dp, status=status)
answer = linear%evaluate(query(1,:))  ! one value per output channel
batch = linear%evaluate(query)        ! (nquery, nchannel)

! A one-shot API for scalar data produces a rank-one array:
answer = griddata(points, values, query, method='linear')
```

Input `points(nsites,ndim)` consists of observation coordinates.
`values(nsites)` or `values(nsites,nchannels)` stores real or complex
observations, with the site dimension first (SciPy's axis-0 convention).
The constructor returns the corresponding real or complex interpolator type.
`%value(q,channel=1)` returns a scalar, and `%evaluate(q)` returns a vector
for a single query, or a matrix for multiple queries.

All constructors accept optional `rescale=.true.`, which scales each
coordinate axis by its range before finding nearest sites or triangulating.
The linear constructors also accept a real or complex `fill_value`, defaulting
to NaN, for points outside the convex hull. The nearest-neighbor method
ignores `fill_value`, as in SciPy.

Inspect the triangulation with `linear%simplices()`, returning
`(ndim+1,nsimplices)` one-based site indices (three per triangle in 2D,
four per tetrahedron in 3D, five in 4D, six in 5D).
Simplex ordering is implementation-dependent and need not match Qhull's ordering.
When `simplices=` is supplied, `%simplices()` returns the supplied array
unchanged, of shape `(ndim+1,nsimplex)`.

### Supplied triangulations in arbitrary dimensions

The `linear_nd_interpolator` constructor also accepts an optional
`simplices=` integer array with shape `(ndim+1, nsimplex)`. Indices are **1-based**,
matching Fortran. When provided, SciFort skips Delaunay construction and
uses a dimension-independent, extended-precision, partially pivoted
barycentric solver with precomputed affine transforms and a neighbor-walk
locator. A dimension-generic bounding-volume hierarchy provides a fallback. This is supported for real and complex values, scalar
or multichannel ordinates, and any positive number of spatial dimensions
**at least two**, provided the input simplices are nonsingular.

For example, a single four-dimensional simplex can be evaluated without a
separate 4D-specific Fortran interpolation routine:

```fortran
use scifort_kinds, only: dp
use scifort_interpolate, only: linear_nd_interpolator, linear_nd_interpolator_type
implicit none
real(dp) :: nodes(5,4), values(5), query(4)
integer :: simplices(5,1), i
type(linear_nd_interpolator_type) :: model
nodes=0.0_dp
do i=1,4
    nodes(i+1,i)=1.0_dp
end do
values=1.0_dp+matmul(nodes,[2.0_dp,-1.0_dp,0.5_dp,3.0_dp])
simplices(:,1)=[1,2,3,4,5]
model=linear_nd_interpolator(nodes,values,simplices=simplices)
query=[0.1_dp,0.2_dp,0.15_dp,0.05_dp]
print *, model%value(query)   ! 1.225
```

To reuse a SciPy/Qhull tessellation, pass the integer simplex IDs from
`scipy.spatial.Delaunay(points).simplices` **plus one** (Python indices start
at zero). The site order must be unchanged. SciFort checks simplex indices,
repeated vertices and affine degeneracy; it does **not** check whether the
supplied simplices form a complete, consistent, nonoverlapping tessellation
or whether they satisfy the Delaunay condition.
Evaluation first attempts a cached neighbor walk starting from the nearest
input site, and falls back to a dimension-generic bounding-volume hierarchy
when the walk cannot locate a containing simplex. For a supplied nonmanifold,
disconnected, or overlapping tessellation, the BVH fallback preserves the
previous ability to find any containing simplex, but the caller remains
responsible for tessellation correctness. In high dimensions the number of
simplices, memory, and time can increase rapidly. The one-shot public
`griddata(method='linear')` constructs native Delaunay triangulations in
2D–5D; higher dimensions still require supplied simplices through the
reusable `linear_nd_interpolator` constructor.

#### Why assumed-rank is not necessary for spatial dimension

In Fortran, `points(nsites,ndim)` has **rank two** for every `ndim`;
`real(dp) :: points(:,:)` is already dimension-independent. `assumed-rank`
(`points(..)`) describes a variable **number of array indices**, not a
variable number of coordinate columns. It would require `select rank` to
recover an ordinary rank-two coordinate view and would not generalize the
Delaunay predicates. `assumed-rank` may be useful later for ordinate arrays
whose output-channel shape has more than one trailing dimension; those
arrays can be flattened at the API boundary. Internally the interpolator
needs only `(nsites,ndim)`, `(ndim+1,nsimplex)` and the barycentric solver.

## Native 4D/5D Delaunay construction and indexed simplex lookup

A single rank-two, dimension-generic Fortran backend implements incremental
four- and five-dimensional Delaunay construction:

```fortran
use scifort_interpolate, only: delaunay_nd, linear_nd_interpolator
integer, allocatable :: simplices(:,:)
! points(nsite,ndim) and values(nsite) are defined by the caller.
call delaunay_nd(points, simplices, status)
! In native 4D, simplices has shape (5,nsimplices), and in native 5D (6,nsimplices).
interpolator = linear_nd_interpolator(points, values, status=status)
```

`delaunay_nd` also dispatches to the existing native 2D and 3D builders.
The direct routine returns `scattered_status_invalid_input` for unsupported
dimensions (including one-dimensional input and dimensions over five),
nonfinite coordinates, repeated sites, and more than **160 input sites** in
4D/5D. With invalid geometry it returns `scattered_status_degenerate`.
The high-dimensional cap bounds scratch-memory allocation and avoids claiming
a production-scale substitute for Qhull.

The new builder uses a greedy, well-conditioned initial simplex, generalized
circumsphere caching, incremental cavity and visible-hull construction, and
face cancellation with a hash. Predicates use `selected_real_kind(33,4931)`
when available, falling back to `real64`. Unlike Qhull, predicates are **not
adaptive exact predicates**; extremely thin, degenerate, or cospherical data
may need the externally supplied SciPy/Qhull simplices route. Complexity can
grow very rapidly with both the number of observations and dimension.

The query path is dimension-independent and uses **three stages**: (1) a KD-tree
nearest-site starting simplex; (2) cached inverse affine transformations and a
walk across shared facets; and (3) a balanced bounding-volume hierarchy
fallback for disconnected supplied meshes, boundary ambiguity, or failed
walks. No assumed-rank coordinate arguments or per-dimension
barycentric code is used. Previously supported `simplices=` inputs in
2D–10D continue to work.

A reproducible `-O3` benchmark, using the *same externally supplied simplices
and 7,500 query evaluations for both versions*, measured these one-run CPU
query times (excluding construction and tessellation):

| Data | Previous linear scan | Cached neighbor walk/BVH | Speedup |
|---|---:|---:|---:|
| 4D, 65 sites, 884 simplices | 0.4267 s | 0.0637 s | 6.7× |
| 5D, 50 sites, 1,810 simplices | 3.0511 s | 0.1096 s | 27.8× |

These are single-run timings on one GNU Fortran build, not performance
guarantees. Run `tools/interpolate_reference/benchmark_simplex_bvh.py` with
`--baseline-zip /path/to/previous_scifort.zip` to replicate the comparison.

For independent SciPy reference comparisons in both compiler configurations:

```bash
python tools/interpolate_reference/compare_auto_nd_scipy.py
python tools/interpolate_reference/compare_auto_nd_scipy.py --seed 8675309
```

A new native test, `test/test_interpolate_delaunay_nd.f90`, verifies affine
polynomial reproduction, real and complex channels, axis rescaling,
out-of-hull fill, and invalid input. The independent Python test checks exact
simplex sets against centered Qhull and complex interpolation against SciPy.

### Three-dimensional example

```fortran
use scifort_kinds, only: dp
use scifort_interpolate, only: linear_nd_interpolator, &
    linear_nd_interpolator_type
implicit none
real(dp) :: sites(5,3), values(5), query(3)
type(linear_nd_interpolator_type) :: interpolator
sites(1,:) = [0.0_dp, 0.0_dp, 0.0_dp]
sites(2,:) = [1.0_dp, 0.0_dp, 0.0_dp]
sites(3,:) = [0.0_dp, 1.0_dp, 0.0_dp]
sites(4,:) = [0.0_dp, 0.0_dp, 1.0_dp]
sites(5,:) = [0.2_dp, 0.2_dp, 0.2_dp]
values = 2.0_dp*sites(:,1) - sites(:,2) + 3.0_dp*sites(:,3)
interpolator = linear_nd_interpolator(sites,values)
query = [0.25_dp, 0.2_dp, 0.1_dp]
print *, interpolator%value(query)   ! equals 0.6 for this affine function
```


## Numerical methods and limits

Nearest-neighbor interpolation builds a balanced KD-tree on the original
or rescaled coordinates and uses bounding-box pruning. Ties select the
lowest original input index, including duplicate input sites. For typical
low-dimensional data, queries inspect far fewer than all sites, although
worst-case time remains linear in the input size. The independent index
also supports reusable `k`-nearest and radius queries; see
[`interpolate_kdtree.md`](interpolate_kdtree.md).

Linear interpolation incrementally constructs Delaunay simplices in 2D–5D
using the Bowyer–Watson cavity algorithm. The 2D/3D specialized builders are
retained; 4D and 5D use the same dimension-independent implementation. It then
searches simplices for a containing triangle/tetrahedron and computes
barycentric weights. The 3D version begins with a nondegenerate tetrahedron,
updates the visible convex hull with each new site, and avoids the potential
hull gaps that can result from removing finite supertetrahedron vertices.
It uses an open-addressed facet hash for finding exterior faces and extended
precision for orientation and in-sphere predicates where supported.
The 3D builder recomputes an exterior-facet hash at each insertion and is
intended for moderate-sized clouds rather than millions of observations.
Construction now caches each tetrahedron's circumsphere center and squared
radius; a conservative distance check rejects most non-cavity tetrahedra
before calling the extended-precision in-sphere determinant. Points near a
sphere boundary retain the original predicate, and the cached spheres are
updated on tetrahedron compaction/insertion. This reduces construction time
but does not change the worst-case quadratic search across tetrahedra.
After tetrahedralization, the 3D constructor builds a **balanced bounding-volume
hierarchy (BVH)** over tetrahedron bounding boxes. Each query traverses this
index, checks individual tetrahedron bounding boxes, then tests the remaining
candidates with the original extended-precision barycentric routine. Leaf size
is eight tetrahedra. The algorithm typically examines only a small fraction of
the mesh, rather than scanning every tetrahedron; worst-case lookup is still
linear in the number of tetrahedra. Building and storing the BVH requires
additional `O(ntetrahedra)` space and approximately `O(ntetrahedra log ntetrahedra)`
work under typical median-partition conditions. The 2D implementation uses
approximately quadratic construction work and scans triangles on evaluation.
The 2D in-circle predicate uses an extended-precision kind when available (with
fallback to double precision); hull-area verification, conservative boundary
notch repair, and local edge flips improve reliability. This covers
well-conditioned scattered inputs but is not a substitute for Qhull's
robust predicates and scaling behavior for strongly degenerate or huge
point clouds.

All native linear constructors reject duplicate input sites and geometrically
degenerate clouds. A native triangulation in `ndim` dimensions requires at
least `ndim+1` affinely independent input sites. The 3D routine centers and uniformly
scales coordinates internally for more stable Delaunay decisions; requesting
`rescale=.true.` first applies SciPy-style independent axis rescaling.
Nonfinite site coordinates and malformed dimensions return an error status. Outside-hull queries are
successful evaluations returning the requested fill value.

`griddata` imported from the public `scifort_interpolate` facade also supports
`method='cubic'` in 2D through the new Clough–Tocher module; see
[`interpolate_clough_tocher.md`](interpolate_clough_tocher.md). Importing
`griddata` directly from the lower-level `scifort_interpolate_scattered`
module still provides only `nearest` and `linear`.
Native Delaunay construction above 5D and direct consumption of Python
`scipy.spatial.Delaunay` **objects** are not implemented; its **integer simplex
array** can be passed through `simplices=` after converting to 1-based indices.
3D simplex location
uses the BVH described above; the 2D linear path has not yet been indexed.
This native implementation does not guarantee the
same degeneracy handling as Qhull's adaptive predicates; near-coplanar or
cospherical clouds can be rejected or choose a different valid triangulation.
In 2D, incomplete convex-hull meshes are detected and rejected explicitly.
In 3D, mesh construction is tested against Qhull's complete hull volume and
simplex topology on randomized, anisotropic and large-offset clouds.

## 3D tetrahedron lookup: numerical behavior and performance

The BVH affects only `LinearNDInterpolator` and `griddata(method='linear')`
on **three-dimensional** site arrays. The Fortran public API, tetrahedra,
interpolation coefficients, convex-hull fills, and extended-precision
barycentric calculation are unchanged. BVH leaves retain original tetrahedron
IDs and choose the lowest-numbered containing tetrahedron at shared faces,
matching the previous full scan. Axis-wise bounding-box margins allow the
same small boundary tolerances as before. Exact-edge queries on highly
ill-conditioned tetrahedra can still be limited by query-coordinate rounding;
the spatial index deliberately does not relax the pre-existing containment
tolerance.

On one GNU Fortran `-O3` run with 160 random sites, 865 tetrahedra, 12,000
three-dimensional queries, and three output channels, evaluation took:

| Algorithm | Construction CPU seconds | Evaluation CPU seconds |
|---|---:|---:|
| Previous full scan | 0.188 | 12.150 |
| BVH (with per-tetra AABB pruning) | 0.187 | 0.135 |

The measured evaluation speedup was **89.9×** (single-run CPU timing, not a
portable performance guarantee). Both paths produced the same checksum and
out-of-hull count. The BVH itself does not accelerate Delaunay construction; a separate
circumsphere cache now accelerates that earlier phase.

Run the indexed benchmark alone, or compare against a locally available
previous-version ZIP:

```bash
python tools/interpolate_reference/benchmark_tetra_bvh.py
python tools/interpolate_reference/benchmark_tetra_bvh.py \
    --baseline-zip /path/to/scifort-interpolate-linearnd-3d-delaunay-2026-10-31.zip
```

The corresponding native regression `test/test_interpolate_tetra_bvh.f90`
checks 447 tetrahedra in a reproducible 90-site cloud, including interior,
near-face, near-edge, vertex, real/complex and multichannel evaluations.

## 3D Delaunay construction benchmark (cached circumspheres)

For a reproducible comparison against the immediately preceding project ZIP:

```bash
python tools/interpolate_reference/benchmark_delaunay_3d_cached.py \
    --baseline-zip /path/to/scifort-interpolate-tetra-bvh-2026-10-31.zip
```

Four-run median GNU Fortran `-O3` timings on one machine, including 3D
construction and BVH build but excluding queries:

| Observation sites | Previous construction | Cached construction | Speedup |
|---:|---:|---:|---:|
| 160 | 0.1886 s | 0.0285 s | 6.6× |
| 320 | 0.8474 s | 0.0800 s | 10.6× |
| 440 | 1.6007 s | 0.1285 s | 12.5× |

For each workload, the tetrahedron count, out-of-hull count and interpolation
checksum matched the previous version. These are specific benchmark
measurements, not portable performance guarantees.

## Higher-dimensional reference verification

Run `tools/interpolate_reference/compare_simplex_nd_scipy.py` to compare
Fortran interpolation given SciPy/Qhull simplices in 2D–7D. The native test
`test/test_interpolate_simplex_nd.f90` additionally tests single-simplex
affine reproduction and input validation in dimensions 2–10.
Native triangulation is independently compared against SciPy/Qhull in 4D–5D.

## RBF nearest-neighbor improvement

The existing `src/interpolate_rbf.f90` now constructs one balanced KD-tree
when local fitting (`neighbors=k`) is requested, replacing the full
`O(n*ndim)` distance scan and `O(n log k)` heap selection at **each** query.
The KD-tree uses exact Euclidean queries with stable input-index tie handling.
Worst-case query cost is still linear in n; the per-query dense local system
solve remains approximately `O(k^3)` without subproblem caching. Global RBF
fits are unchanged.

## Verification

```bash
python tools/interpolate_reference/compare_scattered_scipy.py
python tools/interpolate_reference/compare_scattered_3d_scipy.py
python tools/interpolate_reference/compare_scattered_3d_scipy.py --seed 8675309
python tools/interpolate_reference/check_delaunay_3d_scipy.py
python tools/interpolate_reference/compare_scattered_scipy.py --seed 8675309 --cases 600
python tools/interpolate_reference/compare_rbf_scipy.py
python tools/interpolate_reference/run_all_native_interpolation.py
python tools/interpolate_reference/benchmark_tetra_bvh.py
```

The independent reference suite compares SciFort to installed SciPy 1.17's
`NearestNDInterpolator`, `LinearNDInterpolator`, and `griddata` across
1–5 output channels, real/complex inputs, varying coordinate scales,
rescaling, extrapolation fill, original site queries, and interior points.
The native suite also checks exact affine reproduction, invalid inputs,
degenerate/duplicate sites, and direct public-facade calls. The 3D numerical
suite covers both `LinearNDInterpolator` and `griddata(method='linear')`, real
and complex 1–4-channel data, per-axis rescaling, sample locations, and fill
values. A separate 3D geometry reference verifies volumes and exact simplex
sets against *centered* Qhull input to avoid translation-induced Qhull
roundoff differences at coordinate offsets near 1e6.

The implementation is independently written under SciFort's MIT license.
SciPy and Qhull were used only for black-box reference comparisons; no
third-party source has been copied into SciFort.
