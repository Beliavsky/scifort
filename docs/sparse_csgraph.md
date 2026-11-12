# Sparse graph algorithms (`scipy.sparse.csgraph` subset)

The independent MIT-licensed `scifort_sparse_csgraph` module provides a
one-based **compressed sparse row (CSR)** graph representation, traversal,
shortest paths (including K shortest simple paths), spanning forests, dense and sparse graph Laplacians,
maximum flows, bipartite matchings (including minimum-cost full matchings),
and bandwidth-reducing reorderings. The implementation is
pure modern Fortran with no additional dependencies or borrowed third-party
source. It has deliberately explicit status codes in place of Python errors.

## Example

```fortran
use scifort_kinds, only : dp
use scifort_sparse_csgraph, only : csr_graph, csgraph_from_edges, &
    shortest_path, johnson, laplacian, depth_first_order, &
    path_from_predecessors, csgraph_success
implicit none

type(csr_graph) :: graph
real(dp), allocatable :: costs(:,:), lap(:,:), diagonal(:)
integer, allocatable :: predecessors(:,:), order(:), parents(:), path(:)
integer :: status

! Weighted directed arcs 1->2 (2), 2->3 (-1), 1->3 (10).
graph = csgraph_from_edges(3, [1,2,1], [2,3,3], &
    [2.0_dp,-1.0_dp,10.0_dp], status)
if (status /= csgraph_success) error stop 'invalid graph'

call shortest_path(graph, costs, predecessors, status, method='J')
! costs(1,:) = [0,2,1]; predecessors(1,:) = [0,1,2].
call path_from_predecessors(1, 3, predecessors(1,:), path, status)
! path = [1,2,3].

call depth_first_order(graph, 1, order, parents, status)
call laplacian(graph, lap, diagonal, status, use_out_degree=.true.)
```

## Graph representation

`type(csr_graph)` exposes `n`, `indptr(n+1)`, `indices(:)`, and `data(:)`.
Vertices and CSR offsets are **one-based**: `indptr(1)=1` and
`indptr(n+1)=size(indices)+1`. All edge weights must be finite `real(dp)`.
The zero-vertex graph is supported.

- `csgraph_from_edges(n, source, target, weight, status)` creates a CSR graph
  from parallel arrays (stored zero-weight and duplicate edges allowed).
- `csgraph_from_dense(a, status)` takes a finite square adjacency matrix;
  **zero means no edge**. Use the edge-list constructor to preserve a
  zero-weight edge.
- `csgraph_to_dense(graph, status)` returns a dense adjacency matrix.
  It retains the final edge when parallel arcs target the same cell: use
  explicit aggregation if duplicates must be summed.
- `valid_csgraph(graph)` verifies dimensions, offsets, vertex IDs, and finiteness.

### Implemented algorithms

Let `n` be the number of vertices and `m` the number of stored arcs.

| Routine | Operation | Typical complexity |
|---|---|---|
| `connected_components` | Weak/strong components | O(n+m) |
| `breadth_first_order`, `depth_first_order` | Iterative BFS/DFS orders and predecessors | O(n+m) |
| `breadth_first_tree`, `depth_first_tree` | Weighted, oriented CSR search trees | O(n+m) |
| `dijkstra` | Nonnegative single-source shortest paths; binary heap | O((n+m) log(m+1)) |
| `bellman_ford` | Single-source signed-weight shortest paths and reachable negative-cycle detection | O(nm) |
| `floyd_warshall` | All-pairs shortest paths; signed weights, negative-cycle detection | O(n³), O(n²) storage |
| `johnson` | All-pairs / selected-source shortest paths; signed weights | O(nm + n(m+n) log(m+1)) for all sources |
| `shortest_path` | Dispatcher selecting Floyd–Warshall, Dijkstra, Bellman–Ford, or Johnson | Depends on method |
| `yen` | K shortest simple paths, using a candidate heap and Johnson reweighting for signed costs | Approximately O(nm + K n (m log(m+1) + K n²)) worst-case including duplicate scans |
| `minimum_spanning_tree` | Kruskal minimum spanning forest | O(m log(m+1)) |
| `reverse_cuthill_mckee` | Structural bandwidth-reducing vertex ordering; degree-sorted BFS from pseudo-peripheral roots | O((n+m) log(n+1)) plus root search |
| `maximum_flow` | Dinic integral maximum flow and antisymmetric net-flow CSR | O(n² m) worst case |
| `laplacian` | Dense in-/out-degree, normalized, or symmetrized graph Laplacian | O(n²+m), O(n²) storage |
| `laplacian_sparse` | Sparse CSR Laplacian with the same degree/normalization options | O(n+m), O(n+m) storage |
| `minimum_weight_full_bipartite_matching` | Sparse rectangular min/max-cost assignment | O(k[m log(m+1)+n_rows+n_cols]) worst-case, where k=min(n_rows,n_cols) |
| `reconstruct_path` | Construct a weighted oriented CSR tree from predecessors | O(n+m) |
| `path_from_predecessors` | **SciFort extension:** recover a source-to-target vertex list | O(n) |

The graph traversal/tree routines follow the stored edge order and preserve the
original weight of the selected parent-to-child edge. With
`directed=.false.`, a reversed stored edge is also traversable and can supply
the tree edge weight. Multiple valid traversals may have different orders.

### K shortest simple paths (Yen's algorithm)

`call yen(graph, source, sink, k, distances, predecessors, status, &
 directed, unweighted)` finds up to **K distinct, loopless vertex paths** in
nondecreasing cost order. The input is the existing weighted CSR graph, with
**one-based** source, sink and stored CSR offsets; `k >= 0`. `directed` defaults
to `.true.`; `unweighted=.true.` makes every stored arc cost one. Explicit
zero-cost edges are traversable, and parallel arcs are treated as alternative
costs for the **same vertex-to-vertex connection**, not distinct paths.

`distances(:)` has length `M <= k` (or fewer if the target is unreachable).
`predecessors(:,:)` has shape `(M, graph%n)`, and stores the one-based parent
vertex for each visited non-root node. The sentinel is **zero** for source and
vertices outside that path, rather than SciPy's `-9999`. To turn a row into
its ordered vertices, use `path_from_predecessors(source, sink, &
 predecessors(i,:), path, status)`.

Signed weights are supported via global Bellman–Ford potentials, followed by
nonnegative reduced-cost Dijkstra for each spur. A **negative cycle anywhere**
in the graph returns `csgraph_negative_cycle` (2) and no paths, consistent with
Johnson's global cycle checks. In an undirected graph any negative-weight
edge creates a negative round-trip cycle, unless `unweighted=.true.`.
`source == sink` yields one length-zero route with cost zero if `k > 0` and no
negative cycle. `k=0` returns empty results without performing a cycle check.
Input errors or unrepresentable nonfinite arithmetic return status 1.

```fortran
use scifort_kinds, only : dp
use scifort_sparse_csgraph, only : csr_graph, csgraph_from_edges, &
    yen, path_from_predecessors, csgraph_success
implicit none
type(csr_graph) :: graph
real(dp), allocatable :: distances(:)
integer, allocatable :: predecessors(:,:), path(:)
integer :: status

graph = csgraph_from_edges(4, [1,1,2,2,3], [2,3,4,3,4], &
    [1.0_dp,2.0_dp,1.0_dp,1.0_dp,3.0_dp])
call yen(graph, 1, 4, 3, distances, predecessors, status)
if (status /= csgraph_success) error stop 'yen failed'
! distances = [2.0_dp,5.0_dp,5.0_dp]
call path_from_predecessors(1, 4, predecessors(1,:), path, status)
! path = [1,2,4]
```

Path tie-breaking can differ from SciPy. Yen's algorithm may require O(K n)
shortest-path searches in the worst case; unlike all-pairs algorithms, it
holds a heap of candidate paths rather than constructing an `n x n` matrix.
The implementation uses O(n+m) scratch memory per search, plus stored
accepted/candidate paths. Duplicates are checked by vertex sequence, and
thus the duplicate-checking cost can dominate when K is very large.

### Reverse Cuthill–McKee ordering

`call reverse_cuthill_mckee(graph, permutation, status, symmetric_mode)`
returns a **one-based** new-to-old permutation containing every vertex exactly
once, including isolated vertices. The default `symmetric_mode=.false.`
builds the structural pattern of `A + transpose(A)`; set it true only when
the stored pattern is already symmetric. Edge weights are ignored, while
explicit zero-weight arcs still count as structural edges. Duplicate arcs
and diagonal entries do not inflate vertex degree. Root selection uses
pseudo-peripheral BFS refinement and ascending-degree neighbor traversal;
valid permutations and bandwidth are compared to SciPy, but ties and exact
permutations need not match SciPy.

### Integral maximum flow

`call maximum_flow(graph, source, sink, flow_value, flow, status)` implements
Dinic's level-graph algorithm with a bounded, iterative augmenting-path stack.
The graph is **directed**; `source` and `sink` are distinct one-based vertices.
The input `csr_graph` still stores `real(dp)` data, but capacities must be
finite, nonnegative, **exactly integral**, and at most `2**53 - 1` per arc.
The output `flow_value` is `integer(int64)`, while `flow` is a coalesced CSR
**antisymmetric net-flow matrix**: `flow(i,j)=-flow(j,i)`. For original
capacities `C`, the net flow satisfies `-C(j,i) <= flow(i,j) <= C(i,j)`.
The algorithm supports zero-capacity arcs, parallel and antiparallel arcs,
self-loops, and disconnected networks. The CSR result contains only nonzero
net flows and can safely be converted with `csgraph_to_dense`. It rejects
flows that would overflow signed 64-bit accumulation or require a net CSR
entry not exactly representable in `real(dp)`.

Unlike SciPy's result-object API, SciFort returns separate `flow_value`,
`flow`, and `status` outputs. Its available algorithm is Dinic; an
Edmonds–Karp alternative may be added later.

### Rectangular bipartite matching and structural rank

Bipartite graphs require two separately sized partitions, so SciFort adds
`type(csr_bipartite_graph)` with `n_rows`, `n_cols`, `indptr(:)`, and
`indices(:)`. It is **structural**: the values of connections are irrelevant.
Use `bipartite_from_edges(n_rows,n_cols,rows,columns,status)` for any stored
edge (including explicit zero-valued edges from another representation), or
`bipartite_from_dense(a,status)` for nonzero entries of a rectangular matrix.
`bipartite_from_csgraph(square,status)` converts the existing square CSR.
`valid_bipartite_graph` checks dimensions and one-based CSR offsets.

`call maximum_bipartite_matching(bigraph, matching, status, perm_type)` uses
**iterative Hopcroft–Karp** to find a maximum-cardinality matching. The
`perm_type='row'` default returns a length-`n_cols` vector, with each column's
matched **one-based row**; `'column'` returns a length-`n_rows` vector
containing each row's matched column. The unmatched sentinel is **zero**,
whereas SciPy uses `-1`. Distinct valid maximum matchings may be returned
because traversal tie-breaking is not standardized. Explicitly stored zeros
and repeated structural edges count as connections.

`call structural_rank(bigraph, rank, status)` computes maximum matching
cardinality as the structural rank of the rectangular sparse pattern. Both
empty partitions and non-perfect matchings are supported. The matching
algorithm uses O(m+n_rows+n_cols) storage and O(m sqrt(n_rows+n_cols))
worst-case time for a simple bipartite graph.

### Minimum-weight full bipartite matching

The new **weighted** rectangular CSR representation is
`type(csr_weighted_bipartite_graph)` with `n_rows`, `n_cols`,
`indptr(:)`, `indices(:)`, and `data(:)`. Unlike the separate structural
`csr_bipartite_graph`, it stores finite `real(dp)` edge costs. Use:

- `weighted_bipartite_from_edges(n_rows,n_cols,rows,columns,weights,status)`
  to preserve explicitly stored **zero-weight** connections, negative costs,
  and duplicate edges. When duplicates connect the same vertices, the solver
  chooses the cheapest (or, with `maximize=.true.`, most expensive) edge.
- `weighted_bipartite_from_dense(costs,status)` to interpret **nonzero**
  entries as existing edges. To represent a genuine zero-cost edge use the
  edge-list constructor instead.
- `weighted_bipartite_from_csgraph(graph,status)` to convert square weighted
  CSR graphs with their stored zero-weight edges intact.
- `valid_weighted_bipartite_graph(graph)` to validate the weighted pattern.

`call minimum_weight_full_bipartite_matching(graph,row_ind,col_ind,&
 total_weight,status,maximize)` returns two equally sized vectors of
**one-based** endpoints and a `real(dp)` sum of the original selected costs.
`row_ind` is sorted in ascending order. By default it minimizes the sum;
`maximize=.true.` maximizes it. Every vertex on the *smaller* side is
matched, with no repeated row or column. Tie-breaking is unspecified.

On success, `status=csgraph_success` (0). If no such full matching exists,
`status=csgraph_no_full_matching` (3), and both output index vectors are empty.
Invalid input or numerical overflow gives `csgraph_invalid_input` (1).
A graph with an empty partition has an empty matching with zero total cost.
The implementation uses successive shortest augmenting paths with a binary
heap and feasible dual potentials, avoiding a dense cost matrix. Its auxiliary
storage is O(m+n_rows+n_cols). This is a separate, independently written
algorithm from the SciPy LAPJVsp implementation; it may have different
runtime and tie-breaking characteristics.

```fortran
use scifort_kinds, only : dp
use scifort_sparse_csgraph, only : csr_weighted_bipartite_graph, &
    weighted_bipartite_from_edges, minimum_weight_full_bipartite_matching
implicit none
type(csr_weighted_bipartite_graph) :: costs
integer, allocatable :: row_ind(:), col_ind(:)
integer :: status
real(dp) :: total
costs = weighted_bipartite_from_edges(2,3,[1,1,2,2],[1,2,1,3], &
    [0.0_dp,4.0_dp,3.0_dp,-1.0_dp],status)
call minimum_weight_full_bipartite_matching(costs,row_ind,col_ind,total,status)
! row_ind=[1,2], col_ind=[1,3], total=-1.0_dp
```

Example:

```fortran
use, intrinsic :: iso_fortran_env, only : int64
use scifort_kinds, only : dp
use scifort_sparse_csgraph, only : csr_graph, csr_bipartite_graph, &
    csgraph_from_edges, bipartite_from_edges, maximum_flow, &
    maximum_bipartite_matching, reverse_cuthill_mckee
implicit none
type(csr_graph) :: network, net_flow
type(csr_bipartite_graph) :: pattern
integer, allocatable :: matches(:), perm(:)
integer(int64) :: value
integer :: status
network = csgraph_from_edges(3,[1,1,2],[2,3,3],[4.0_dp,1.0_dp,3.0_dp])
call maximum_flow(network,1,3,value,net_flow,status) ! value=4
call reverse_cuthill_mckee(network,perm,status)
pattern = bipartite_from_edges(2,3,[1,2,2],[3,1,2])
call maximum_bipartite_matching(pattern,matches,status,'row') ! 2 matches
```

### All-pairs shortest paths

`floyd_warshall(graph, distances, predecessors, status, directed, unweighted)`
returns full `n x n` arrays. The algorithm is suitable for small or dense
graphs; it allows signed edges, but rejects negative cycles. Its costs are
symmetric in undirected mode, where reciprocal arcs are compared using the
smaller weight. `unweighted=.true.` sets each stored arc's cost to one,
regardless of the stored weight.

`johnson(graph, distances, predecessors, status, directed, unweighted,
indices)` runs Bellman–Ford from an implicit super-source followed by
binary-heap Dijkstra with reduced costs. It supports negative edges but
rejects **negative cycles anywhere in the graph**, including disconnected
components. `indices=[3,1]` selects these *one-based* sources, in that order,
returning a `2 x n` distance matrix; absent `indices` means all `n` sources.
An empty `indices` array returns a `0 x n` matrix.

Both routines return `+Inf` for unreachable pairs and `predecessors=0` for
unreachable or same-vertex pairs (SciPy uses the sentinel `-9999`). Other
predecessors are one-based vertex IDs. On invalid input or a negative cycle,
the arrays are allocated with shape `(0,0)` and the status identifies why.

`shortest_path(graph, distances, predecessors, status, method, directed,
unweighted, indices)` is the common dispatcher. `method` accepts `'auto'`
(default), `'FW'`, `'D'`, `'BF'`, or `'J'`. The automatic method uses
Floyd–Warshall on dense graphs with `n <= 256`, at least 25% stored density,
and all-source requests; otherwise it uses Johnson. Method `'D'` rejects
negative-weight graphs unless `unweighted=.true.`. Unlike SciPy, SciFort
supports `indices` with explicit `'FW'` by computing all pairs then selecting
rows. Method `'BF'` detects cycles reachable from the selected sources;
`'J'` and `'FW'` check the entire graph.

**Compatibility note:** SciFort treats a negative-weight self-loop as a
negative cycle in both Floyd–Warshall and Johnson. SciPy 1.17's Floyd–Warshall
ignores negative diagonal loops, whereas its Johnson implementation reports
the negative cycle.

`dijkstra` and `bellman_ford` retain their earlier single-source signatures:

```fortran
call dijkstra(graph, source, distances, predecessors, status, &
              directed=.true., limit=10.0_dp)
call bellman_ford(graph, source, distances, predecessors, status, &
                  directed=.true.)
```

Dijkstra rejects any negative edge and an invalid or negative limit.
Bellman–Ford detects negative cycles **reachable from its source**. For
undirected traversals, reciprocal edges with unequal weights are interpreted
as alternative arcs, so their smaller cost is available in either direction.
SciPy documents undefined behavior for some undirected Dijkstra searches
with asymmetric reciprocal weights.

### Laplacian

`laplacian(graph, lap, diagonal, status, normed, use_out_degree,
symmetrized)` constructs a **dense** `real(dp)` Laplacian. By default it
uses the incoming (column-sum) degree, consistent with SciPy. Set
`use_out_degree=.true.` for outgoing row sums, or `symmetrized=.true.` to
construct the Laplacian of `A + transpose(A)` (not divided by two).
The diagonal of the adjacency matrix is ignored and parallel weights are
summed, rather than treated as separate degree contributions.

With `normed=.false.`, `lap = diag(degree) - A_offdiag` and `diagonal`
contains the weighted degree. With `normed=.true.`, the off-diagonal entries
are divided by `sqrt(degree(i)*degree(j))` and diagonal entries are one for
non-isolated vertices or zero for zero-degree vertices. The returned
`diagonal` holds `sqrt(degree)`, with zero degrees represented by one.
Negative weighted degree is rejected for normalization (status 1 rather than
SciPy's NaN). Unnormalized Laplacians allow signed degrees. `laplacian_sparse(graph, lap, diagonal, status, normed, use_out_degree,
 symmetrized)` implements the **same options** but returns a `csr_graph`
result in one-based sparse format. Repeated off-diagonal edges are summed,
exact zero cancellations are dropped, and stored self-loops do not contribute
to degrees. It stores only nonzero matrix elements: at most `n+m` entries in
the ordinary case or `n+2m` with `symmetrized=.true.`, and uses O(n+m)
working memory. CSR entries are not promised to be sorted by column.
The result can be converted to dense for small graphs using
`csgraph_to_dense(lap,status)`. On invalid input or unrepresentable results,
`status=csgraph_invalid_input` and the returned `lap` should not be used.

```fortran
use scifort_kinds, only : dp
use scifort_sparse_csgraph, only : csr_graph, csgraph_from_edges, laplacian_sparse
implicit none
type(csr_graph) :: graph, graph_lap
real(dp), allocatable :: diagonal(:)
integer :: status
graph = csgraph_from_edges(3,[1,2],[2,3],[1.0_dp,2.0_dp],status)
call laplacian_sparse(graph,graph_lap,diagonal,status, &
    symmetrized=.true.,normed=.true.)
```

### Trees and path reconstruction

`breadth_first_tree(graph, start, status, directed)` and
`depth_first_tree(graph, start, status, directed)` return an oriented CSR
forest with only the edges discovered from the specified root. Their weights
are taken from the corresponding graph arcs; a reverse arc may be used in
undirected mode.

`reconstruct_path(graph, predecessors, status, directed)` matches SciPy's
*tree-construction* API (with one-based predecessors and zero for roots).
It validates the entire parent forest, rejects parent cycles, and requires
all referenced arcs to exist in the given graph.

`call path_from_predecessors(source, target, predecessors, path, status)`
is a separate SciFort extension yielding a one-dimensional array of vertices
from source through target. It returns an empty `path` for an unreachable
target; an invalid predecessor cycle returns status 1.

## Status codes

| Name | Value | Meaning |
|---|---:|---|
| `csgraph_success` | 0 | Successful computation |
| `csgraph_invalid_input` | 1 | Invalid argument, CSR, index, method, or normalized degree |
| `csgraph_negative_cycle` | 2 | Negative-weight cycle detected |
| `csgraph_no_full_matching` | 3 | Full bipartite assignment is impossible |

## Testing

No new libraries, Fortran packages, or vendored sources are required.
The native and reference tests need `gfortran`; the randomized comparison
scripts additionally require NumPy and SciPy.

```sh
python tools/cluster_reference/run_native_csgraph.py --mode both
python tools/csgraph_reference/compare_csgraph_scipy.py --mode both
python tools/csgraph_reference/compare_csgraph_paths_scipy.py --mode both
python tools/csgraph_reference/compare_csgraph_extended_scipy.py --mode both
python tools/csgraph_reference/compare_csgraph_dispatch_scipy.py --mode both
python tools/csgraph_reference/compare_csgraph_advanced_scipy.py --mode both
python tools/csgraph_reference/compare_csgraph_weighted_scipy.py --mode both
python tools/csgraph_reference/compare_csgraph_yen_scipy.py --mode both
```

The comparison scripts check signed and unsigned shortest-path distances,
K-shortest simple routes (also checked by exhaustive enumeration),
negative cycles, predecessor paths, normalized/unnormalized graph Laplacians,
DFS reachability, dispatcher methods, and BFS/DFS weighted trees against
SciPy 1.17. Existing native tests also cover empty graphs, invalid CSR inputs,
zero-weight and parallel arcs, negative diagonal self-loops,
matching cardinality/structural rank, maximum-flow conservation, and ordering
validity. The advanced randomized suite additionally uses exhaustive min-cut
and matching checks on small graphs.

Not yet implemented: Edmonds–Karp maximum flow or the full set of SciPy's
CSR matrix conversion helpers. The CSR object is a numerical graph container,
not a clone of SciPy's general-purpose sparse matrix classes.

### Testing

The standalone native regression runner builds the graph, weighted-matching,
and hierarchy shortcut tests with both checked and optimized GNU Fortran:

```shell
python tools/cluster_reference/run_native_csgraph.py --mode both
```

Randomized SciPy and exhaustive small-assignment reference checks:

```shell
python tools/csgraph_reference/compare_csgraph_weighted_scipy.py --mode both
```

The reference script requires Python, NumPy, SciPy and `gfortran` only for
validation; the Fortran module itself adds no runtime dependencies.
