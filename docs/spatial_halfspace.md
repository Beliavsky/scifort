# Bounded halfspace intersections (`scipy.spatial.HalfspaceIntersection` subset)

SciFort exports `halfspace_intersection(halfspaces, interior_point, status=...)`
and `halfspace_intersection_type` from `scifort_spatial`. It finds the unique
vertices of a bounded full-dimensional convex polygon (2D) or polyhedron (3D)
represented as the intersection of affine halfspaces. No Qhull or external
linear algebra dependency is required.

## Example

```fortran
use scifort_kinds, only : dp
use scifort_spatial, only : halfspace_intersection, halfspace_intersection_type
implicit none
real(dp) :: inequalities(4,3)
type(halfspace_intersection_type) :: region
integer :: status

! -1 <= x <= 1 and -2 <= y <= 2, as a*x + b <= 0.
inequalities(1,:) = [ 1.0_dp, 0.0_dp, -1.0_dp]
inequalities(2,:) = [-1.0_dp, 0.0_dp, -1.0_dp]
inequalities(3,:) = [ 0.0_dp, 1.0_dp, -2.0_dp]
inequalities(4,:) = [ 0.0_dp,-1.0_dp, -2.0_dp]
region = halfspace_intersection(inequalities, [0.0_dp,0.0_dp], status=status)
if (region%valid()) then
    print *, 'number of vertices:', size(region%intersections,1)
    print *, 'vertices (one per row):'
    print *, region%intersections
end if
```

### Inputs

- `halfspaces(n,d+1)`: each row is an inequality
  `sum(halfspaces(i,1:d)*x) + halfspaces(i,d+1) <= 0`.
- `interior_point(d)`: **strictly interior** reference point. It must satisfy
  all inequalities with strict `< 0`; SciFort does not compute this point.
- `status` (optional): `spatial_distance_success=0` or
  `spatial_distance_invalid_input=1`.

Only dimensions **2 and 3** are supported. Inputs must be finite. Redundant
and duplicate inequalities are supported when their strict-interior condition
is met. A zero normal with negative offset is redundant; a zero normal with
nonnegative offset is not strictly feasible.

### Output type

| Component | Shape | Meaning |
| --- | --- | --- |
| `intersections` | `(nv,d)` | Unique primal vertices, one per row, unordered |
| `dual_points` | `(n,d)` | `-A_i / (A_i . interior_point + b_i)` |
| `dual_facets` | `(d,nf)` | One-based indices of triangular (3D) or line (2D) polar facets |
| `dual_vertices` | `(nactive)` | One-based input constraint indices on the polar hull |
| `dual_equations` | `(d+1,nf)` | Unit outward polar normals and affine offsets, column by facet |
| `dual_area` | Scalar | Perimeter of 2D polar hull, surface area of 3D polar hull |
| `dual_volume` | Scalar | Area of 2D polar hull, volume of 3D polar hull |
| `status` | Scalar | Zero on success, one on failure |

The `halfspaces` and `interior_point` components retain the original inputs.
`valid()` reports whether the result is usable. Unlike SciPy, SciFort stores
facet-index and facet-equation arrays in **Fortran column-major orientation**:
column `j` describes facet `j`. `intersections` and `dual_points` follow
SciPy's row-of-points convention. `dual_facets` triangulates non-simplicial
polar faces in 3D and may use different triangulations or vertex orders than
SciPy. Returned primal vertices are unique; SciPy/Qhull can include multiple
nearly identical intersections when a dual polygonal face is triangulated.

### Algorithm and limitations

The implementation scales coordinates around the supplied interior point,
computes normalized polar constraints, then enumerates intersections of all
pairs of boundary lines (2D) or triples of planes (3D). Every candidate must
satisfy **every** inequality; nearby duplicate vertices are coalesced. A
recession-direction test rejects unbounded sets. Each primal vertex determines
a polar supporting face, and the existing native 2D convex hull is used to
triangulate that face and compute polar area and volume. This avoids relying
on the incremental 3D convex hull builder for halfspace vertex completeness.

Worst-case work is approximately **O(n^3 * n)** for enumerating and checking
3D candidates, where `n` is the number of inequalities. The algorithm is
intended for moderate `n`, not very large polyhedral systems. It does not
implement SciPy's arbitrary-dimensional Qhull backend, `incremental`,
`add_halfspaces`, `qhull_options`, or calculation of an interior point.
Extremely thin/near-degenerate polytopes may be rejected or have vertices
merged within the documented floating-point tolerance. Invalid, unbounded,
and numerically unsupported inputs share status code 1; no partial result
should be used in that case.

### Reproducible tests

```sh
python tools/spatial_reference/run_native_spatial.py --mode both
python tools/spatial_reference/compare_halfspace_scipy.py --mode both --count 200 --seed 20261102
python tools/spatial_reference/compare_halfspace_scipy.py --mode both --count 220 --seed 20261111
```

The Python comparator uses SciPy 1.17 `HalfspaceIntersection` solely as a
numerical oracle. Geometrically coincident SciPy intersection points are
coalesced for comparison, since their multiplicity depends on Qhull's polar
facet triangulation. Both native status tests and independent randomized
geometry checks are included in the archive.
