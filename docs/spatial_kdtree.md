# `scipy.spatial` — native KD-tree nearest, pair, and periodic searches

The MIT-licensed `scifort_interpolate_kdtree` module extends its existing
balanced point index with four exact operations corresponding to
`scipy.spatial.cKDTree`: `query_pairs`, `query_ball_tree`, `count_neighbors`,
and `sparse_distance_matrix`. `count_neighbors` also supports radius arrays
with cumulative or histogram-bin outputs, and `weighted_count_neighbors`
returns real-valued pair-product sums. These support Minkowski orders `p>=1` and
periodic boxes, including selected nonperiodic axes. The same public type is accessible via both
`scifort_spatial` and `scifort_interpolate`.  These methods use only native
Fortran and do not require BLAS, LAPACK, C++, or external dependencies.

## Fortran API

```fortran
use, intrinsic :: iso_fortran_env, only: int64
use scifort_kinds, only: dp
use scifort_spatial, only: kd_tree, kd_tree_type

type(kd_tree_type) :: a, b
integer, allocatable :: pairs(:,:), row_ptr(:), indices(:), row(:), column(:)
real(dp), allocatable :: distances(:)
integer(int64) :: npairs

a = kd_tree(reshape([0.0_dp, 1.0_dp, 3.0_dp], [3,1]))
b = kd_tree(reshape([0.0_dp, 2.0_dp], [2,1]))

call a%query_pairs(1.0_dp, pairs)
! pairs(1,:) == [1,2]; the distance from site 2 to site 3 is 2,
! so there is only one pair within radius 1.

call a%query_ball_tree(b, 1.0_dp, row_ptr, indices)
! neighbors for source point i: indices(row_ptr(i):row_ptr(i+1)-1)

npairs = a%count_neighbors(b, 1.0_dp)
! Cumulative pair counts at three cutoffs (any radius order is accepted).
block
    integer(int64), allocatable :: cumulative_counts(:), bin_counts(:)
    real(dp), allocatable :: weight_sums(:)
    cumulative_counts = a%count_neighbors(b, [1.0_dp, 0.0_dp, 2.0_dp])
    bin_counts = a%count_neighbors(b, [0.0_dp, 1.0_dp, 2.0_dp], cumulative=.false.)
    weight_sums = a%weighted_count_neighbors(b, [0.0_dp, 1.0_dp, 2.0_dp], &
        weights_self=[2.0_dp, 1.0_dp, 3.0_dp], weights_other=[1.0_dp, 4.0_dp])
end block
call a%sparse_distance_matrix(b, 1.0_dp, row, column, distances)
```

**API contracts:**

- `query_pairs(radius, pairs, status, visited)` returns an allocated integer
  array `pairs(npairs,2)` of **one-based**, strictly increasing `(i,j)`
  indices with `i<j`.  Results are in lexicographic order and include exact
  duplicate-site pairs when the radius is zero.
- `query_ball_tree(other, radius, row_ptr, indices, status, visited)` returns
  **one-based CSR**.  `size(row_ptr)=self%size()+1`, `row_ptr(1)=1`, and
  `row_ptr(size(row_ptr))=size(indices)+1`.  Destination indices are sorted
  ascending within each row.  Repeated and coincident input sites retain
  their separate indices.
- `count_neighbors(other, radius, status, visited, p)` returns an
  `integer(int64)` count of **ordered** `(source,destination)` pairs.  For a
  tree compared with itself, self-matches and both `(i,j)` and `(j,i)` are
  included.  This differs intentionally from `query_pairs`, which returns
  only `i<j` pairs.
- `count_neighbors(other, radii, status, visited, p, cumulative)` accepts a
  rank-one array of thresholds and returns an allocated `integer(int64)`
  vector. `cumulative` defaults to true, counting distances `<= radii(i)`;
  radii need not be sorted and duplicates are allowed. For
  `cumulative=.false.`, radii must be nondecreasing, and the bins are
  `distance <= radii(1)` and `radii(i-1) < distance <= radii(i)` thereafter.
  Equal adjacent thresholds yield empty bins. An empty `radii` vector is valid.
- `weighted_count_neighbors(other, radius, weights_self, weights_other,
  status, visited, p)` returns one `real(dp)` pair-weight sum; passing `radii`
  rather than `radius` instead returns an allocatable `real(dp)` vector and
  permits optional `cumulative`. Each qualifying `(i,j)` contributes
  `weights_self(i)*weights_other(j)`. Either weight array may be omitted
  (unit weights on that side); signed and zero weights are supported. The
  weights must be finite and have lengths equal to their respective tree
  sizes. Weighted results are not restricted to integers.
- `sparse_distance_matrix(other, max_distance, row, column, distance,
  status, visited, p)` returns COO-style, one-based row and column arrays and
  **unsquared** Minkowski distances, sorted by `(row,column)`.  Explicit
  zero distances are retained for identical points and self-matches.  There
  is no required sparse matrix package.

All radii must be finite, nonnegative numbers.  Inputs must be valid trees
with matching positive coordinate dimensions.  Every search uses a **closed
ball** (`distance <= radius`).  The optional `status` is
`kdtree_success=0` or `kdtree_invalid_input=1`.  On invalid input, array
results are empty and scalar `count_neighbors` returns zero (with an error
status).
`visited` reports actual site-pair distance evaluations, *not* bounding-box
checks; it may be zero for wholly included subtree pairs. Unsorted radii in
noncumulative mode, negative/nonfinite radii, and invalid weight arrays cause
`kdtree_invalid_input`; vector-returning methods then return empty arrays.  Initial tree
construction, `nearest`, `query`, and `query_ball_point` retain their
existing contracts documented in [interpolate_kdtree.md](interpolate_kdtree.md).

## Periodic boxes and Minkowski distance

`kd_tree(points, boxsize=...)` accepts either a **scalar** box length, which
applies to all coordinates, or a **rank-one** array of length `size(points,2)`.
Each strictly positive length `L(j)` makes coordinate `j` periodic. A zero
length makes that coordinate nonperiodic; omitting `boxsize` also means all
coordinates are nonperiodic. All box lengths must be finite and nonnegative.
For periodic coordinates, input sites must lie in `[0,L(j))`; otherwise the
constructor returns `kdtree_invalid_input`. Queries may lie outside the box
and are wrapped to the fundamental domain. Cross-tree operations require
identical `boxsize` arrays in both trees.

All search operations accept an **optional trailing `p`** argument:

- `query(..., p=p, distance=dist)`, `nearest(..., p=p)`, and
  `query_ball_point(..., p=p)`;
- `query_pairs(..., p=p)`, `query_ball_tree(..., p=p)`,
  `count_neighbors(..., p=p)`, `weighted_count_neighbors(..., p=p)`,
  and `sparse_distance_matrix(..., p=p)`.

`p` defaults to 2 (Euclidean). Valid values are real `p>=1` or positive
IEEE infinity (Chebyshev). If a periodic coordinate has box length `L`, its
absolute coordinate gap is `min(abs(dx), L-abs(dx))`, after wrapping both
coordinates; all gaps are combined into the chosen Minkowski norm. Searches
remain **exact and inclusive**. Invalid `p`, including NaN, and incompatible
periodic topologies return `kdtree_invalid_input`.

```fortran
use scifort_kinds, only: dp
use scifort_spatial, only: kd_tree, kd_tree_type
type(kd_tree_type) :: tree
integer, allocatable :: ids(:), pairs(:,:)
real(dp), allocatable :: distances(:)

tree = kd_tree(reshape([0.05_dp, 0.95_dp], [2,1]), boxsize=1.0_dp)
call tree%query([1.02_dp], ids, k=2, p=1.0_dp, distance=distances)
! ids are [1,2], distances approximately [0.03,0.07]
call tree%query_pairs(0.101_dp, pairs, p=1.0_dp)
! pairs(1,:) is [1,2], crossing the periodic boundary
```

For backward compatibility, `query(distance2=...)` returns the **square of
the chosen metric** (ordinary squared Euclidean distance by default), even
when `p/=2`. Use the new `query(distance=...)` optional argument for true,
unsquared Minkowski distances. One-based indices and sorted outputs are
unchanged. Note that SciPy uses zero-based indices and returns unsquared
distances in its usual `query` result.

## Algorithm

The balanced median-partitioned tree stores axis-aligned bounding boxes and
subtree cardinalities.  Every radius search prunes bounding boxes whose
minimum possible distance exceeds the radius.  `count_neighbors` uses a
**dual-tree traversal**: if the maximum possible pairwise distance between
two subtrees is within the radius, their full cardinality product is added
without visiting any individual sites. For an array of radii, a single
traversal resolves all thresholds together; it sorts cumulative query radii
internally and records interval contributions to a difference array. This
avoids an independent traversal per threshold. Weighted queries precompute
subtree sums of site weights and add their products when an entire pair of
subtrees qualifies.  Otherwise it splits whichever side
contains more sites, recursing through the point and child subtrees.  Pair
and CSR enumeration use one query-tree search per source point and grow
output storage geometrically.  Only actual matches are retained in the
result; worst-case output size is necessarily quadratic for dense searches.

Closed-ball comparisons use a scaled Minkowski norm, avoiding spurious
underflow from raising approximately `1e-200` differences to a power or
overflow from raising approximately `1e200` differences to a power.
For periodic boxes, minimum bounding-box distances consider the wrapped
axis intervals, while safe upper bounds use at most half the periodic length
per coordinate. Pruning remains correct, though periodic wide boxes can
yield looser bounds than nonperiodic boxes.

## Deliberate differences from SciPy

Only exact-radius searches are supported: approximate `eps` queries
are not yet exposed. Scalar/vector radius counting and noncumulative binning
are available, but scalar-radius **enumeration** (`query_pairs`,
`query_ball_tree`, `query_ball_point`, and `sparse_distance_matrix`) is
unchanged. The return type of existing scalar `count_neighbors` remains
`integer(int64)`, so weighted counting uses a separate bound method. Negative `boxsize` lengths are rejected (zero means
nonperiodic). Different box topologies on two trees cannot be mixed. Indices
are one-based and outputs are Fortran arrays or CSR/COO triplets rather than
Python sets, lists, or SciPy sparse matrices.  Tree construction currently
rejects empty point sets, which SciPy's KDTree can represent.

## Independent validation

```bash
python tools/spatial_reference/compare_kdtree_pairs_scipy.py
python tools/spatial_reference/compare_kdtree_pairs_scipy.py --seed 8675309
python tools/spatial_reference/compare_kdtree_periodic_scipy.py
python tools/spatial_reference/compare_kdtree_periodic_scipy.py --seed 8675309
python tools/spatial_reference/compare_kdtree_weighted_count_scipy.py
python tools/spatial_reference/compare_kdtree_weighted_count_scipy.py --seed 8675309
python tools/spatial_reference/run_native_spatial.py
python tools/interpolate_reference/compare_kdtree_scipy.py
python tools/interpolate_reference/run_all_native_interpolation.py
```

The comparison script builds checked and optimized GNU Fortran executables
and independently verifies SciPy 1.17's `cKDTree` operations.  Native
regressions test one through seven dimensions, mixed periodic axes, scalar
and vector `boxsize`, Minkowski `p=1,1.5,2,3,infinity`, duplicate sites, exact
closed-ball boundaries, zero-radius matches, invalid inputs, self-counting,
fully covered subtree pruning, vector-radius histograms, signed and
one-sided weights, and extreme coordinate magnitudes.

## Provenance

All new KD-tree procedures, periodic bounding-box calculations, stable
Minkowski norms, CSR/COO handling, weighted and vector-radius dual-tree
count traversal, tests,
and this documentation are independently written MIT-licensed SciFort code.
No code has been copied from SciPy, FLANN, Qhull, or other spatial libraries.
