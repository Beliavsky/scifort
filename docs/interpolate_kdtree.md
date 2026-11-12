# Reusable balanced KD-tree for nearest-neighbor interpolation

The independent MIT-licensed `scifort_interpolate_kdtree` module implements a
balanced spatial index (Euclidean by default) and integrates it with the existing
`NearestNDInterpolator` and local `RBFInterpolator` code. Both also work
without a BLAS/LAPACK library or additional FPM dependency. SciPy 1.17's
`scipy.spatial.cKDTree` is used only for independent numerical validation.

## Fortran API

```fortran
use scifort_kinds, only: dp
use scifort_interpolate, only: kd_tree_type, kd_tree

type(kd_tree_type) :: tree
real(dp) :: sites(4,2)
integer, allocatable :: neighbors(:), within_radius(:)
real(dp), allocatable :: distance2(:)
integer :: best, status, visited

sites(1,:) = [0.0_dp, 0.0_dp]
sites(2,:) = [1.0_dp, 0.0_dp]
sites(3,:) = [0.0_dp, 1.0_dp]
sites(4,:) = [1.0_dp, 1.0_dp]
tree = kd_tree(sites,status=status)
best = tree%nearest([0.2_dp, 0.1_dp])
call tree%query([0.2_dp, 0.1_dp], neighbors, distance2=distance2, &
    k=3, visited=visited)
call tree%query_ball_point([0.2_dp, 0.1_dp], 0.9_dp, within_radius)
```

`sites(npoints,ndim)` must have at least one point and one dimension, with
finite real coordinates. The query dimension must agree and its coordinates
must be finite. The public interface uses **one-based** site indices;
`query()` reports **squared Euclidean distances** by default, avoiding an
unnecessary square root for internal distance comparisons. Indices are returned in
ascending `(squared_distance, original_input_index)` order. Duplicate sites
are permitted; exact distance ties favor the first input site. Radius
queries include the boundary and return matching original indices sorted in
ascending index order. `status` is zero on success, one on invalid input.

`nearest()` returns zero and sets an error status for an invalid query;
`query()` and `query_ball_point()` return empty arrays on invalid requests.
`%valid()`, `%dimensions()`, and `%size()` inspect a constructed tree.
`visited=` reports the number of exact point-to-query distances calculated,
excluding bounding-box tests; it is useful for testing pruning, but it is
not a hardware-independent runtime predictor.

This is a focused subset of `scipy.spatial.KDTree` functionality. **Exact
`query_pairs`, `query_ball_tree`, `count_neighbors`, and `sparse_distance_matrix`
are now implemented**, including pruning and one-based CSR/COO results; see
[spatial_kdtree.md](spatial_kdtree.md). The constructor additionally supports
periodic `boxsize` (scalar or per-dimension vector), and every search accepts
Minkowski `p>=1` or `p=infinity`; default Euclidean behavior is unchanged.
Use `query(distance=...)` for ordinary, unsquared Minkowski distances. The
legacy `distance2` output remains the square of the selected metric.
`count_neighbors` now also supports multiple radii and cumulative or
noncumulative outputs; `weighted_count_neighbors` handles weighted pair
sums (see [spatial_kdtree.md](spatial_kdtree.md)). Approximate searches
(`eps`) remain unsupported. SciFort intentionally
uses one-based indices. The implementation uses binary64 data.

## Algorithms and complexity

The builder uses deterministic median partitioning at each level, cyclic
coordinate-axis selection, and an axis-aligned bounding box per subtree.
Construction averages `O(n log n)` comparisons and uses `O(nd)` space.
The exact query maintains a max-heap of the best `k` candidates and prunes
subtrees whose bounding boxes are strictly farther than its worst current
candidate. The final heap sort costs `O(k log k)`. Typical low-dimensional
nearest queries visit a small fraction of the input; **worst-case** query
cost is still `O(nd + n log k)`. Performance can degrade in high dimensions.
The radius query prunes by the subtree's minimum distance from the ball and
sorts `m` returned matches in `O(m log m)` time.

NearestND builds its tree on the same rescaled coordinates that its
previous full scan used. The 2D Delaunay/LinearND code is unchanged. Local
RBF fitting builds its tree once on raw observation sites; the per-query
`k`-nearest selection is accelerated, but each local RBF system is still
solved from scratch with dense Gaussian elimination (roughly `O(k^3)`).
There is no local coefficient cache, GPU implementation, or shared thread
mutable state in the index.

## Verification and benchmarking

```bash
python tools/interpolate_reference/compare_kdtree_scipy.py
python tools/interpolate_reference/compare_kdtree_scipy.py --seed 8675309 --cases 220
python tools/interpolate_reference/compare_scattered_scipy.py
python tools/interpolate_reference/compare_rbf_scipy.py
python tools/interpolate_reference/run_all_native_interpolation.py
python tools/interpolate_reference/run_kdtree_benchmark.py 12000 1500 2
```

The native regression also checks exact tie ordering (including duplicate
sites), `k=1` through `k=n`, radius queries against an independently computed
brute-force reference, invalid inputs, one through seven dimensions, and
subtree-pruning visit counts. The randomized SciPy comparison checks
`cKDTree.query` and `cKDTree.query_ball_point` in one through seven dimensions.
SciPy's ordering is not defined for exact duplicate-site ties; those are
validated against a local brute-force oracle rather than its tie behavior.

A single GNU Fortran `-O3` run on the local test system using 12,000 2D
sites and 1,500 queries gave approximately 0.00134 s for the KD-tree nearest
queries after a 0.00299 s build, versus 0.03745 s for brute force.
Average exact site-distance evaluations were 15.12 / 12,000 (0.126%).
These are exploratory single-run CPU timings and not reproducible guarantees
or a substitute for benchmarking on target systems.

## Provenance

All tree implementation and validation code was independently authored as
SciFort MIT-licensed code. No C++/Python SciPy source, `nanoflann`, Qhull,
or third-party Fortran spatial-index source was copied. Existing FPM
requirements remain unchanged.
