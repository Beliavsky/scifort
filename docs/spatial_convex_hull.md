# 2D and 3D convex hulls (`scipy.spatial.ConvexHull` subset)

`scifort_spatial` exports `convex_hull(points, status=...)` and
`convex_hull_type`. The implementation constructs a full-dimensional convex
polygon or convex polyhedron from an arbitrary set of real(dp) points, with
interior points and exact duplicate points allowed. It uses only native Fortran,
with no Qhull, BLAS, or LAPACK dependencies.

## Example

```fortran
use scifort_kinds, only : dp
use scifort_spatial, only : convex_hull, convex_hull_type
implicit none
type(convex_hull_type) :: hull
real(dp) :: points(5,2)
integer :: status

points(1,:) = [0.0_dp,0.0_dp]
points(2,:) = [1.0_dp,0.0_dp]
points(3,:) = [1.0_dp,1.0_dp]
points(4,:) = [0.0_dp,1.0_dp]
points(5,:) = [0.5_dp,0.5_dp]
hull = convex_hull(points,status=status)
if (hull%valid()) then
    print *, 'perimeter:', hull%area
    print *, 'enclosed area:', hull%volume
    print *, 'indices of polygon vertices:', hull%vertices
end if
```

The constructor accepts `points(n,d)`, where `d=2` or `3` and `n>=d+1`.
Input points are given in **rows**, following SciPy. All indices in results
are **one-based**.

## Public components

| Component | Shape | Description |
| --- | --- | --- |
| `points` | `(n,d)` | Original input coordinates |
| `vertices` | `(nvertices)` | Original indices of extreme points; CCW around a polygon (2D), in ascending input order (3D) |
| `simplices` | `(d,nfacets)` | Boundary segments (2D) or triangular facets (3D) |
| `neighbors` | `(d,nfacets)` | Adjacent facet opposite each corresponding local vertex |
| `equations` | `(d+1,nfacets)` | Outward unit normal in entries `1:d`, followed by plane offset |
| `area` | Scalar | Polygon perimeter (2D), polyhedron surface area (3D) |
| `volume` | Scalar | Polygon area (2D), polyhedron volume (3D) |
| `status` | Scalar | `spatial_distance_success=0` or `spatial_distance_invalid_input=1` |

For every hull facet `f` and every point `x` inside the hull, the equation
satisfies `dot_product(equations(1:d,f),x) + equations(d+1,f) <= 0` within
floating-point tolerance. Each `neighbors(k,f)` is the neighboring facet
sharing all vertices except `simplices(k,f)`. The arrays are **transposed
relative to SciPy's row-major conceptual shapes**: Fortran columns represent
facets, rather than rows.

The `valid()` method returns false if construction fails (for example,
collinear 2D sites or coplanar 3D sites). Results for invalid geometry are
not meaningful. The routine centers and uniformly rescales input internally
to improve the conditioning of geometric orientation tests. It does not modify
the original coordinates or permit nonsensical nonfinite values.

## Algorithms and limits

The 2D method uses a lexicographically sorted monotone chain, removing
non-extreme collinear boundary points. Sorting currently uses insertion
sort, so this implementation is intended for moderate data sets. The 3D
method starts with a well-spread tetrahedron and incrementally replaces
visible hull facets with a consistently oriented horizon mesh. It also
constructs a watertight per-facet neighbor table and obtains the volume
by summing signed tetrahedral volumes around an interior reference point.

The 3D constructor now independently validates Euler topology and the
supporting plane of every triangular facet **against every original site**.
If the incremental mesh is incomplete or non-manifold, it retries up to 24
deterministic point orderings, remaps all indices to the input array, and
revalidates the result. A final fallback enumerates globally supporting
planes, projects their approximately coplanar sites to 2D and triangulates
each polygonal face. Only validated hulls are returned as successful.
The fallback is approximately O(n^4) in the worst case, so it is intended
for moderate data sizes and is invoked only when the faster construction
cannot pass its checks. Its plane classification uses an internal tolerance
on uniformly normalized coordinates, followed by a stricter geometric
containment validation. Facets closer than floating-point resolution may
be represented by different vertices than Qhull even if the geometric
region, surface area, and volume agree.

This implementation does **not** provide SciPy/Qhull's general `ndim>3`
hulls, `qhull_options`, incremental `add_points`, `coplanar`, or `good`
facet flags. Extremely close-to-coplanar inputs and extremely large point
clouds would benefit from exact/adaptive orientation predicates and more
scalable hull algorithms. The previously reported 26-point polar-hull
failure is now covered by a native regression and independent SciPy oracle,
but the implementation cannot guarantee recovery of every arbitrarily
ill-conditioned full-dimensional input. The separate
`halfspace_intersection` implementation continues to enumerate plane
intersections directly. Coplanar facets may be triangulated differently
from Qhull, and identical duplicate coordinates may have different
representative input indices.

## Reproducible tests

```sh
python tools/spatial_reference/run_native_spatial.py --mode both
python tools/spatial_reference/compare_convex_hull_scipy.py --mode both --seed 6919 --count 250
python tools/spatial_reference/compare_convex_hull_scipy.py --mode both --seed 20731 --count 250
python tools/spatial_reference/compare_convex_hull_nearcoplanar_scipy.py --mode both
```

The independent oracle uses SciPy `ConvexHull` for 2D/3D geometry only;
no SciPy or Qhull source code is included. Facet numbering and triangulation
are not assumed to match. Each comparison tests the hull extreme-point
positions, area, volume, outward facet unit normals, point containment,
neighbor topology, and 2D CCW ordering.

The new deterministic 26-point fixture is stored in
`tools/spatial_reference/convex_hull_nearcoplanar_26.txt` and embedded in
`test/test_spatial_convex_hull_nearcoplanar.f90`, which exercises twelve
input permutations. The separate near-coplanar SciPy script tests 100
synthetic 3D cube-face clouds with coordinate perturbations from zero to
`1e-8`. Their geometry is compared numerically, but their exact extreme
indices are not asserted: at roundoff-scale perturbations those indices
need not agree between valid triangulations.
