# Hierarchical agglomerative clustering

The public module `scifort_cluster_hierarchy` implements the core of
`scipy.cluster.hierarchy` with an idiomatic native Fortran interface.
The numerical code is independently implemented from the published
Lance-Williams update equations, without BLAS, LAPACK, or additional FPM
packages.

## APIs

```fortran
use scifort_kinds, only : dp
use scifort_cluster_hierarchy, only : linkage, single, complete, average, &
    weighted, centroid, median, ward, fcluster, fclusterdata, &
    fcluster_distance, fcluster_maxclust, cophenetic_distances, &
    cophenetic_correlation, leaves_list, &
    inconsistent, maxdists, maxrstat, maxinconsts, &
    fcluster_inconsistent, fcluster_monocrit, fcluster_maxclust_monocrit, &
    cut_tree, is_monotonic, num_obs_linkage, &
    optimal_leaf_ordering, to_tree, leaders, correspond, cluster_tree_node, &
    is_isomorphic, is_valid_linkage, is_valid_im, &
    to_mlab_linkage, from_mlab_linkage, dendrogram, cluster_dendrogram, &
    cluster_success, cluster_invalid_input

real(dp) :: observations(5,2)
real(dp), allocatable :: z(:,:), zordered(:,:), cophenetic(:), stats(:,:), crit(:)
integer, allocatable :: labels(:), order(:), cuts(:,:), leader_pairs(:,:)
type(cluster_tree_node), allocatable :: nodes(:)
type(cluster_dendrogram) :: layout
integer :: status

! Initialize observations with five two-dimensional data points.
z = linkage(observations, method='ward', status=status)
if (status /= cluster_success) error stop 'invalid linkage input'

labels = fcluster_maxclust(z, 3, status=status)
labels = fcluster(z, 3.0_dp, criterion='maxclust', status=status)
labels = fclusterdata(observations, 1.5_dp, criterion='distance', &
    method='ward', status=status)
cophenetic = cophenetic_distances(z, status=status)
! Correlation is computed against the original condensed distances,
! for example those returned by scifort_spatial_distance.pdist.
order = leaves_list(z, status=status)

stats = inconsistent(z, depth=2, status=status)
crit = maxinconsts(z, stats, status=status)
labels = fcluster_inconsistent(z, 0.8_dp, status=status, r=stats)
labels = fcluster_monocrit(z, crit, 1.0_dp, status=status)
labels = fcluster_maxclust_monocrit(z, crit, 3, status=status)
cuts = cut_tree(z, n_clusters=[2,3,4], status=status)

! Original merge distances and memberships remain unchanged.
zordered = optimal_leaf_ordering(z, observations, status=status)
order = leaves_list(zordered)
nodes = to_tree(z, status=status)
leader_pairs = leaders(z, labels, status=status)
layout = dendrogram(z, status=status)
! layout%icoord(k,:) and layout%dcoord(k,:) hold the U-link vertices.
! One-based root node: nodes(2*size(observations,1)-1)
```

`linkage` is a generic function with two overloads:

- `linkage(x(:,:), method, metric, status, optimal_ordering)` computes pairwise distances
  between the rows of `x` using SciFort's `pdist` implementation.
- `linkage(y(:), method, status, optimal_ordering)` accepts the condensed upper-triangle vector
  in SciPy `pdist` order, with length `n*(n-1)/2`.

Both return an allocatable `real(dp)` array `z(n-1,4)`. Each row records
`[left_child, right_child, merge_distance, leaf_count]`. **Node indices are
one-based**, unlike SciPy's zero-based linkage matrices: leaf nodes are `1:n`,
and the merged node on row `k` has identifier `n+k`.

`method` defaults to `'single'` and accepts the seven standard rules:
`'single'`, `'complete'`, `'average'`, `'weighted'`, `'centroid'`, `'median'`,
and `'ward'`. The seven identically named convenience generics
`single`, `complete`, `average`, `weighted`, `centroid`, `median`, and `ward`
accept the same condensed or observation inputs and return exactly the
corresponding `linkage(..., method='...')` result, including the optional
`optimal_ordering` setting. The observation overload also accepts `metric`.
The observation overload's `metric` defaults to `'euclidean'`
and supports the metrics exposed by `scifort_spatial_distance.pdist`.
For nondefault metric parameters (such as Minkowski `p` or feature weights),
call `pdist` explicitly and pass its condensed result to `linkage`.
Both overloads accept `optimal_ordering=.true.` to reorder each pair of
children after constructing the tree, minimizing the original distances
between consecutive leaves. The option defaults to `.false.`.
Centroid, median, and Ward linkage have their standard clustering meaning
only for Euclidean distances.

`fcluster_distance(z, threshold, status)` yields one-based positive labels
for flat clusters whose **maximum descendant merge height** is at most
`threshold`. The maximum-descendant rule is important for nonmonotone centroid
and median dendrograms. `fcluster_maxclust(z, max_clusters, status)` chooses
the smallest such threshold yielding at most the requested number of clusters.
Cluster labels are assigned in first-leaf encounter order, not necessarily
SciPy's numeric label ordering; partitions agree with the documented criteria.

`cophenetic_distances(z, status)` returns an `n*(n-1)/2` vector of merge
heights for leaf pairs, in condensed ordering. `leaves_list(z, status)`
returns the one-based leaf order obtained by traversing left children before
right children.

## End-to-end clustering and cophenetic correlation

`fcluster(z, threshold, criterion, depth, r, monocrit, status)` dispatches
among all five SciPy flat-clustering rules:

- `inconsistent` (default): calculate (or reuse `r`) and cut using the
  maximum descendant inconsistency coefficient. `depth=2` by default.
- `distance`: threshold the maximum descendant merge distance.
- `maxclust`: find the smallest distance threshold producing no more than
  the requested number of groups.
- `monocrit`: threshold a supplied monotone per-merge `monocrit` vector.
- `maxclust_monocrit`: find the minimum monotone-criterion threshold
  yielding no more than the requested cluster count.

The two routines accept `real(dp)` or ordinary `integer` thresholds through
explicit generic overloads. Positive fractional real values
for the maximum-cluster criteria are truncated toward zero (like SciPy).
A threshold that truncates below 1, or nonfinite thresholds, are rejected
rather than returning every sample as a singleton; use the typed
`fcluster_maxclust(z, integer_count)` routines when possible. The
`monocrit` rules require the named vector and validate its tree
monotonicity. The result is a one-based integer label vector. Positive label
*values* may differ from SciPy, but cluster membership agrees. For invalid
input, the result has length zero and `status=cluster_invalid_input`.

`fclusterdata(x, threshold, criterion, metric, depth, method, r, status)`
implements SciPy's convenience workflow by computing pairwise distances
from `x(n,features)`, constructing a linkage tree, and cutting it.
`criterion='inconsistent'`, `metric='euclidean'`, `depth=2`, and
`method='single'` are the defaults. It accepts `'inconsistent'`,
`'distance'`, or `'maxclust'` as documented for SciPy `fclusterdata`;
for custom monotone criteria, build `z` and call `fcluster` directly.
`r` may contain precomputed inconsistency statistics for the *same*
linkage tree. The routine supports metrics implemented by SciFort's `pdist`
without its additional optional parameterization (for weighted Minkowski,
custom `p`, etc., use `pdist` and `linkage` explicitly). It requires at
least two observations and returns one-based labels.

`cophenetic_correlation(z, y, status)` returns the **Pearson correlation**
between original condensed distances `y` and `cophenetic_distances(z)`.
It is the scalar coefficient returned by SciPy's `cophenet(z, Y)`;
`cophenetic_distances(z)` already provides the vector output.
The implementation normalizes each vector before computing covariance
and variance to prevent squaring large distances. It returns quiet NaN
with **success status** if either vector has zero variance (undefined
correlation), matching SciPy; malformed linkage, distance lengths, or
nonfinite/negative distances return invalid-input status and quiet NaN.

Example:

```fortran
use scifort_kinds, only : dp
use scifort_spatial_distance, only : pdist
use scifort_cluster_hierarchy, only : linkage, fclusterdata, &
    cophenetic_correlation, cluster_success
implicit none
real(dp) :: x(4,1), correlation
real(dp), allocatable :: z(:,:), y(:)
integer, allocatable :: labels(:)
integer :: status

x(:,1) = [0.0_dp, 1.0_dp, 4.0_dp, 10.0_dp]
labels = fclusterdata(x, 3.0_dp, criterion='distance', status=status)
if (status /= cluster_success) error stop 'fclusterdata failed'
z = linkage(x, method='single', status=status)
y = pdist(x)
correlation = cophenetic_correlation(z, y, status=status)
if (status /= cluster_success) error stop 'cophenetic correlation failed'
```

## Optimal leaf ordering and tree utilities

`optimal_leaf_ordering` is overloaded for condensed distances and observation
matrices:

- `optimal_leaf_ordering(z, y(:), status)` takes a full condensed `pdist`
  vector containing the original pairwise dissimilarities.
- `optimal_leaf_ordering(z, x(:,:), status, metric)` takes observation rows;
  `metric` is optional and defaults to `'euclidean'`.

Both return a new, valid linkage matrix whose merge distances, counts, and
unordered child pairs are **identical to the input**. Only the order of left
and right children can change. They minimize
`sum(distance(leaves(k),leaves(k+1)), k=1,...,n-1)` over *all* leaf orders
consistent with the fixed tree. The solver uses endpoint-pair dynamic
programming and traceback of child orientations, not a greedy heuristic;
it handles tied distances and trees with inverted merge heights.

The dynamic program is useful for modest or medium-size dendrograms. Its
runtime is at most **O(n^4)** with the direct endpoint-pair recurrence, and
its working storage is **O(n^2)**: each possible pair of leaves creates a
bridge decision only at its lowest common ancestor, and intermediate cost
matrices are discarded once their parent has been processed. SciPy's compiled
implementation may be faster on larger datasets. Input distances are scaled
by their largest value before optimization to avoid overflowing accumulated
path costs; results can be affected by floating-point ties.

`to_tree(z, status)` returns an allocated array of `cluster_tree_node` entries,
indexed directly by one-based node ID. Nodes have integer `id`, `left`,
`right`, `parent`, `count` and real `distance` fields. Leaves have `count=1`,
`distance=0`, and `left=right=0`. The root is `nodes(2*n-1)`, and its
`parent=0`. This indexed tree avoids pointer ownership/allocation and makes
iterative tree algorithms straightforward.

`leaders(z, labels, status)` returns an integer array `(n_clusters,2)`:
`[leader_node_id, flat_cluster_label]` in each row. It checks that every
label corresponds to exactly one complete subtree, accepting arbitrary
integer label IDs (not only 1,...,k); scattered/disconnected assignments
are rejected with an empty result and `cluster_invalid_input`. As with all
node IDs here, leaders are **one-based**, unlike SciPy. Row ordering is by
node ID, not necessarily SciPy's leader-array order.

`correspond(z, y, status)` returns a logical indicating whether linkage `z`
and condensed distance vector `y` encode the **same observation count**.
As in SciPy, it does not verify that merge distances match the supplied
pairwise distances. A count mismatch returns `.false.` with success status;
invalid linkage or nonfinite/negative pairwise distances produce an
invalid-input status.

## Validation, partition isomorphism and MATLAB interoperability

The hierarchy module also offers the standard SciPy-style predicates and
linkage conversions:

- `is_valid_linkage(z)` validates a **one-based** linkage matrix with shape
  `(n-1,4)`. It checks integer, previously created child IDs; use of each
  child exactly once; finite nonnegative linkage distances; and exact
  cardinalities of merged clusters.
- `is_valid_im(r)` accepts a finite, nonempty `(n-1,4)` inconsistency matrix
  with nonnegative means, standard deviations, and link counts. It rejects
  nonfinite entries, a stricter policy than SciPy 1.17's internal validator.
- `is_isomorphic(t1, t2, status)` returns true exactly when the two integer
  label vectors encode the same set partition, irrespective of relabeling.
  The labels may be negative, zero, or nonconsecutive. Its sorting algorithm
  runs in **O(n log n)** time and **O(n)** additional memory. Different
  lengths produce `.false.` and `cluster_invalid_input`; empty vectors
  represent the same partition, returning `.true.`. During randomized
  comparison, SciPy 1.17 sometimes returned *false positives* for a repeated
  cluster label appearing nonadjacently after sorting by the other assignment;
  SciFort deliberately checks both directions to avoid this error.
- `to_mlab_linkage(z, status)` returns the MATLAB-compatible three-column
  linkage matrix. Because SciFort already uses one-based node identifiers,
  no ID shift is needed; the fourth column is simply removed.
- `from_mlab_linkage(m, status)` accepts a valid MATLAB `(n-1,3)` linkage
  matrix and reconstructs its exact `(n-1,4)` SciFort representation,
  including node cardinalities. Both conversions reject invalid topology,
  and return zero-row arrays with `cluster_invalid_input` on failure.

## Dendrogram geometry without plotting

`dendrogram(z, status, count_sort, distance_sort)` provides an allocatable
`cluster_dendrogram` structure with these fields:

- `leaves(n)` — the left-to-right sequence of **one-based** leaf IDs.
- `icoord(n-1,4)` — the horizontal-axis coordinates of the four vertices
  for each U-shaped branch (standard SciPy spacing: 5, 15, 25, ... for leaves).
- `dcoord(n-1,4)` — the corresponding vertical coordinates, with linkage
  heights on the horizontal bars.

Branches are returned in **depth-first postorder**, matching
`scipy.cluster.hierarchy.dendrogram(z, no_plot=True)` with no truncation.
The routine does not create an image, require graphics libraries, calculate
colors, or modify `z`. Optional `count_sort` and `distance_sort` each accept
`'none'` (default), `'ascending'`, or `'descending'` and follow SciPy's
ordering/tie rules. They cannot both be non-`'none'` in one call. Incorrect
arguments or invalid linkage return empty arrays and an invalid-input status.

```fortran
use scifort_cluster_hierarchy, only : linkage, dendrogram, cluster_dendrogram
use scifort_kinds, only : dp
implicit none
real(dp) :: x(4,1)
real(dp), allocatable :: z(:,:)
type(cluster_dendrogram) :: layout
x(:,1) = [0.0_dp, 1.0_dp, 4.0_dp, 10.0_dp]
z = linkage(x, method='average')
layout = dendrogram(z, count_sort='ascending')
! Each row of layout%icoord and layout%dcoord specifies four vertices.
```

This implements the numeric output, not SciPy's plotting and truncation
options (`truncate_mode`, `p`, `labels`, colors and orientation).

## Inconsistency statistics and extended flat clustering

`inconsistent(z, depth, status)` returns the `(n-1,4)` per-merge
statistics of SciPy's `inconsistent`: the mean linkage height, **sample**
standard deviation, number of inspected links, and inconsistency coefficient.
The default `depth=2` includes the current merge and one level of descendant
merges. `depth=0` acts like `depth=1`, matching SciPy. A zero standard
deviation gives an inconsistency coefficient of zero. The calculation scales
merge heights before accumulating the sample variance to reduce intermediate
overflow/cancellation; numerical values can therefore differ slightly from
SciPy's finite-precision accumulation.

`maxdists(z)` returns the maximum merge height over each non-singleton
node and its descendants, a monotonically nondecreasing *tree* criterion
(even if the linkage row heights are nonmonotone). `maxrstat(z,r,column)`
returns the maximum value of the selected statistic over each subtree;
**`column` is 1 through 4**, not Python's 0 through 3.
`maxinconsts(z,r)` is the `column=4` specialization. The supplied
statistic array must be finite and have shape `(n-1,4)`.

Additional flat-clustering operations:

- `fcluster_inconsistent(z, threshold, status, depth, r)` groups nodes whose
  maximum descendant inconsistency coefficient does not exceed `threshold`.
  Supply an optional precomputed `r` matrix to avoid recalculating it.
- `fcluster_monocrit(z, criterion, threshold, status)` groups clusters under a
  nondecreasing *tree* criterion of length `n-1`.
- `fcluster_maxclust_monocrit(z, criterion, max_clusters, status)` chooses the
  smallest criterion threshold producing at most the requested count.

A criterion is accepted only if its value on each parent merge is at least
that of its non-singleton children. An invalid criterion, depth, count, or
statistic matrix produces `cluster_invalid_input` and an empty result.

## Full and selected cut trees

`cut_tree(z, n_clusters, height, status)` returns a matrix whose rows are the
original observations and whose columns are selected partitions. It has
three mutually exclusive modes:

- Neither optional selector: returns all `n` columns for `n, n-1, ..., 1`
  clusters, following linkage-table row order.
- `n_clusters=[k1,k2,...]`: returns exactly those requested cluster counts,
  each in its own column, allowing arbitrary order and repeated counts.
- `height=[h1,h2,...]`: returns each height-based cut in its own column.
  **Only monotone linkage trees are accepted.** As in SciPy's `cut_tree`,
  only merges with height strictly **less than** `h` are applied.

Unlike threshold-based `fcluster_maxclust`, `cut_tree` with `n_clusters`
returns the exact number even when several merge heights are tied.
Every output uses **one-based labels**, in first-leaf encounter order.
The two selectors cannot be combined. `is_monotonic(z)` checks successive
linkage heights and `num_obs_linkage(z)` returns the number of leaves (or
zero for invalid input).

**Compatibility note:** SciPy's `cut_tree` sorts merge events by their
heights, including on nonmonotone centroid/median trees. SciFort retains
linkage row order when making exact-count cuts on such trees, preserving the
original agglomeration sequence; partitions may differ in this special case.
Its multiple-count/multiple-height requests also avoid a SciPy 1.17 issue in
which unsorted selections may leave some columns with default zeros. The
numerical comparison script checks these cases against separate scalar
SciPy calls.

All routines accept an optional integer `status` output. On invalid input the
status is `cluster_invalid_input` (1) and an empty array is returned. Invalid
linkage trees (incorrect node references, repeated children, invalid counts,
nonfinite heights) are rejected; library procedures do not abort execution.

## Numerical and performance notes

- The algorithm maintains a full pairwise cluster distance matrix and
  exhaustively selects the nearest active pair. It uses **O(n^2) memory** and
  **O(n^3) time**, suitable for moderate-size datasets. SciPy's optimized
  minimum-spanning-tree and nearest-neighbor-chain algorithms are faster for
  the methods that admit them.
- Centroid, median and Ward linkages use scale-normalized squared-distance
  updates to reduce overflow/underflow in the Lance-Williams formulas.
- Merge ties are broken deterministically by ascending one-based child IDs;
  SciPy can choose a different but equally valid tied dendrogram. Optimal
  leaf ordering is supported with `optimal_ordering=.true.` or the standalone
  `optimal_leaf_ordering` routine. Dendrogram plotting remains out of scope.
- SciFort's `maxclust` computes the mathematical minimum threshold yielding at
  most the requested number of clusters. SciPy's internal `maxclust` algorithm
  has implementation-specific behavior on some very small trees. Labels are
  not assumed numerically identical between implementations.
- Only real-valued input is supported. `dendrogram` plotting and several
  validation/conversion convenience routines remain gaps.
  No new third-party dependencies are required.

## Verification

`test/test_cluster_hierarchy_flatdata.f90` tests `fcluster` criterion
selection, `fclusterdata` default/overridden options and precomputed
statistics, SciPy's documented 12-observation example, reference and
scale-invariant cophenetic correlations, degenerate two-leaf/constant
correlations, and invalid inputs. `compare_hierarchy_flatdata_scipy.py`
compares randomized observation-to-labels workflows against SciPy 1.17
for seven linkage methods, six metrics, and three flat-cut criteria,
including direct `fcluster` equivalence and correlation coefficients.

`test/test_cluster_hierarchy.f90` tests exact reference trees, condensed
and observation inputs, all linkage methods, cophenetic distances, flat
clusters, leaf orders, duplicated points, nonmonotone centroid heights,
large coordinate scales and invalid inputs. It also tests the new inconsistent
statistics, subtree maxima, monotone criteria, exact-count and strict-height
cuts, duplicate merge heights, and invalid-input handling.

`tools/cluster_reference/compare_hierarchy_scipy.py` compares against
`scipy.cluster.hierarchy.linkage`, `cophenet`, `leaves_list` and `fcluster`
with independent random seeds. The new
`compare_hierarchy_diagnostics_scipy.py` compares statistics and additional
partition cuts using SciPy-generated linkage matrices; it uses a relative
tolerance for inconsistency coefficients to account for different accumulation
rounding. An additional native regression, `test/test_cluster_hierarchy_ordering.f90`,
checks the optimal-ordering overloads, `linkage(optimal_ordering=.true.)`,
node/parent relationships, arbitrary leader labels, inverted trees,
zero/tied and extremely large distances, and rejection of invalid inputs.
The independent `compare_hierarchy_ordering_scipy.py` compares 190 randomized
linkage trees in **each** build mode, and explicitly enumerates every
legal child-flip orientation for 68 small random trees per mode. The algorithm
reaches the global minimum in these tests; several SciPy 1.17 orientations
have slightly higher path costs, so the comparisons require the new
ordering to be no worse than SciPy rather than match its exact permutation.
`compare_hierarchy_leaders_scipy.py` compares randomized flat-cluster leader
mappings, including arbitrary integer IDs, against the SciPy public API.
These scripts use SciPy as a numerical oracle, not as an implementation source.
The scripts compare partitions rather than numeric label values. Example:

```bash
python tools/cluster_reference/run_native_cluster.py --mode both
python tools/cluster_reference/compare_hierarchy_flatdata_scipy.py --mode both --cases 504
python tools/cluster_reference/compare_hierarchy_scipy.py --mode both --cases 300
python tools/cluster_reference/compare_hierarchy_diagnostics_scipy.py --mode both --cases 180
python tools/cluster_reference/compare_hierarchy_ordering_scipy.py --mode both --cases 190
python tools/cluster_reference/compare_hierarchy_leaders_scipy.py --mode both --cases 200
```

## Algorithm references and provenance

The algorithms follow the standard hierarchical clustering definitions and
Lance-Williams recurrence documented in:

- SciPy 1.17.0 `linkage` documentation:
  https://docs.scipy.org/doc/scipy-1.17.0/reference/generated/scipy.cluster.hierarchy.linkage.html
- SciPy `fcluster` and `fclusterdata` documentation:
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.fcluster.html
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.fclusterdata.html
- SciPy `cophenet` documentation:
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.cophenet.html
- SciPy `optimal_leaf_ordering`, `to_tree`, and `leaders` documentation:
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.optimal_leaf_ordering.html
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.to_tree.html
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.leaders.html
- SciPy `inconsistent`, `maxRstat` and `cut_tree` documentation:
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.inconsistent.html
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.maxRstat.html
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.cut_tree.html
- F. Murtagh and P. Contreras, *Algorithms for hierarchical clustering: an
  overview*, WIREs Data Mining and Knowledge Discovery 2 (2012), 86-97.

No source code was imported or copied from SciPy or an external implementation.

## Additional regression and SciPy comparison

```sh
python tools/cluster_reference/run_native_cluster.py --mode both
python tools/cluster_reference/compare_cluster_hierarchy_layout_scipy.py --cases 250 --optimization=-O0
python tools/cluster_reference/compare_cluster_hierarchy_layout_scipy.py --cases 250 --optimization=-O3
```

The comparison uses 250 seeded randomized hierarchy examples for each build,
covering all seven linkage types and five sorting modes. It compares exact
leaf sequences and branch coordinates against SciPy 1.17. It also checks
MATLAB linkage roundtrips and mathematical partition isomorphism, and reports
any divergences from SciPy's predicate separately.
