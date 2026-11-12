# SciFort spherical Voronoi diagrams (2026-11-02)

`scifort_spatial` now includes the SciPy-compatible **`spherical_voronoi`**
constructor for **two- or three-dimensional** generators constrained to the
surface of a common circle/sphere. This is independently implemented,
MIT-licensed Fortran with no external convex-hull, BLAS, or LAPACK dependency.
The constructor supports arbitrary generator order, a non-origin center,
positive radius, and configurable geometric validation tolerance.

## API

```fortran
use scifort_kinds, only : dp
use scifort_spatial, only : spherical_voronoi, spherical_voronoi_type
implicit none
type(spherical_voronoi_type) :: diagram
real(dp) :: generators(6,3)
real(dp), allocatable :: areas(:)
integer, allocatable :: region(:)
integer :: status

! Populate the six 3-D generator rows before calling.
diagram = spherical_voronoi(generators, radius=1.0_dp, &
                            center=[0.0_dp,0.0_dp,0.0_dp], &
                            status=status)
if (status == 0) then
    areas = diagram%calculate_areas()
    region = diagram%region(1)
    ! diagram%vertices contains the (vertices, 3) coordinates.
    ! region contains one-based indices into diagram%vertices.
end if
```

`diagram` contains these public components:

- `points(n,d)`: original input coordinates.
- `vertices(m,d)`: Voronoi vertices on the sphere/circle.
- `simplices(d,m)`: outward hull triangular facets in 3D and cyclic
  generator-pair hull edges in 2D. All indices are **one-based**.
- `region_offsets(n+1)` and `region_vertices(:)`: compressed sparse row
  storage of the per-generator Voronoi region vertex lists. For generator
  `i`, use `region_vertices(region_offsets(i):region_offsets(i+1)-1)` or
  the convenience method `diagram%region(i)`.
- `radius`, `center`, `status`, and `diagram%valid()`.

The 3D region indices are in consistent cyclic order, obtained by walking
adjacent triangular hull faces rather than sorting floating-point angles.
Their starting positions and the global numbering of facets **may differ from
SciPy**, but the geometric vertices and corresponding per-generator regions
are equivalent. 2D regions contain the two endpoints of each generator's
Voronoi arc.

### Region measures: important SciPy 2D distinction

- **3D:** `calculate_areas()` returns the area of each *generator's* spherical
  Voronoi polygon. The values sum to `4*pi*radius**2`.
- **2D:** To match the implementation in SciPy 1.17.0, `calculate_areas()`
  returns the original circular *hull-edge* arc lengths (one per Voronoi
  vertex, indexed in the same order as `simplices`), **not** the arc length
  of each generator's Voronoi region. These sum to `2*pi*radius`.
- **Both dimensions:** The SciFort extension `region_measures()` returns one
  measure **per generator**: a Voronoi arc length in 2D or a spherical polygon
  area in 3D.

The separate 2D method avoids silently conflating SciPy's hull-edge arc
length convention with per-generator Voronoi regions. For equally spaced
2D points the two sets of lengths happen to coincide.

### Validation and limitations

`points` has shape `(n,d)`, with `d` equal to 2 or 3, `n>=d+1`, and all
coordinates finite. `radius` defaults to 1.0, `center` defaults to the origin,
and `threshold` defaults to `1e-6` as a relative duplicate and on-sphere
tolerance. Inputs must lie within that relative radial discrepancy of the
specified sphere and must not have generator pairs closer than `threshold`
when normalized by the radius. A 3D generator set without full affine rank
is rejected. On failure, `status=spatial_distance_invalid_input=1` and
`diagram%valid()` is false; its vertex/region arrays are not populated.

3D hull construction uses an incremental, outward-oriented triangle mesh,
with a wide initial tetrahedron, visible-face removal and horizon replacement.
The resulting hull is traversed via edge adjacency to produce ordered regions.
Region areas sum spherical triangle solid angles evaluated with a stable
`atan2` formulation. It uses `O(n)` scratch capacity and is intended for
modest numbers of generators, not as a replacement for Qhull for very large
or deliberately near-degenerate data sets. Coplanar hull facets and nearly
coincident points can require more robust exact predicates and Qhull-style
symbolic perturbation, which are **not** implemented here.

## Reproducible tests

```sh
python tools/spatial_reference/run_native_spatial.py --mode both
python tools/spatial_reference/compare_voronoi_scipy.py --mode both
python tools/spatial_reference/compare_voronoi_scipy.py --mode both --seed 87915 --count 210
```

Randomized tests compare 2D/3D geometry, vertex clouds and generator-region
membership (independent of indexing), 3D per-generator areas, 2D hull-edge
arcs, translated/scaled spheres, and clustered/hemispherical generators
against SciPy 1.17.0. Both bounds-checked and optimized GNU Fortran builds
are covered. No third-party source was imported.

Other spatial features: [distance functions](spatial_distance.md),
[Hausdorff and SLERP](spatial_geometry.md), [Procrustes](spatial_procrustes.md).
