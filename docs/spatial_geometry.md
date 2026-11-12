# SciFort `scipy.spatial` geometry (2026-11-02)

The public `scifort_spatial` module now includes **directed Hausdorff
point-cloud distances** and **geometric spherical interpolation**. This is
original MIT-licensed Fortran with no new dependencies. The library uses
`real(dp)` (binary64) coordinates, with point clouds stored as
`(observations, dimensions)` arrays, matching the convention used by
`pdist`, `cdist`, and SciPy.

## Hausdorff distances

```fortran
use scifort_kinds, only : dp
use scifort_spatial, only : directed_hausdorff, hausdorff_distance
implicit none
real(dp) :: a(2,2), b(2,2), d, symmetric
integer :: index_a, index_b, status
logical :: from_a

a(1,:) = [0.0_dp, 0.0_dp]
a(2,:) = [3.0_dp, 0.0_dp]
b(1,:) = [0.0_dp, 0.0_dp]
b(2,:) = [1.0_dp, 0.0_dp]

d = directed_hausdorff(a, b, index_u=index_a, index_v=index_b, &
                       status=status)  ! d=2, pair=(2,2)
symmetric = hausdorff_distance(a, b, index_u=index_a, &
    index_v=index_b, from_u=from_a, status=status)  ! symmetric=2
```

The directed distance is the exact max-min Euclidean distance:

\[
h(U,V)=\max_{u\in U}\min_{v\in V}\|u-v\|_2.
\]

`index_u` and `index_v` are optional output indices, **one-based** and
identifying a pair attaining the directed distance. Input point clouds must
contain at least one observation and one coordinate, must have the same
feature dimension, and must be finite. An invalid call returns NaN,
`spatial_distance_invalid_input=1`, and zero indices. Valid calls return
`spatial_distance_success=0`.

The implementation uses deterministic original-row ordering and early
abandonment when a source row cannot improve the incumbent max-min distance.
Pass `early_break=.false.` to disable this optimization for verification.
Unlike SciPy, SciFort does **not shuffle point order or use random state**.
When multiple pairs attain the maximum, witness indices may differ from
SciPy's default shuffled witnesses, though the distance is the same.

### SciFort extension: symmetric Hausdorff distance

`hausdorff_distance(u, v, index_u=, index_v=, from_u=, status=)` computes
`max(h(u,v),h(v,u))`. This convenience procedure is a **SciFort extension**;
SciPy users typically call `directed_hausdorff` twice. The optional
`from_u` output is `.true.` when the maximizing direction is from `u` to
`v`; in a directional tie it prefers `u` to `v`. Indices always refer to
the original `u` and `v` arrays, regardless of direction.

Distances are evaluated using scaled Euclidean norms to avoid premature
underflow or overflow when the true binary64 result is representable.
Overflow of the true distance may still yield positive infinity.
The algorithm has worst-case `O(n_source*n_target*n_dimensions)` time,
`O(n_dimensions)` scratch space, and can terminate many inner loops early.

**Documented edge difference:** SciPy 1.17.0 may produce unspecified
witness indices for empty point arrays; SciFort rejects empty inputs. Inputs
containing NaN or infinity are also rejected explicitly, instead of
attempting to reproduce undefined or version-dependent results.

## `geometric_slerp`

```fortran
use scifort_kinds, only : dp
use scifort_spatial, only : geometric_slerp
implicit none
real(dp), allocatable :: point(:), path(:,:)
integer :: status

point = geometric_slerp([1.0_dp,0.0_dp,0.0_dp], &
    [0.0_dp,1.0_dp,0.0_dp], 0.5_dp, status=status)
! point = [1/sqrt(2), 1/sqrt(2), 0]

path = geometric_slerp([1.0_dp,0.0_dp,0.0_dp], &
    [0.0_dp,1.0_dp,0.0_dp], [0.0_dp,0.5_dp,1.0_dp], status=status)
! path shape: (3 times, 3 dimensions)
```

`geometric_slerp(start, finish, t[, tol, status])` is rank-overloaded for
scalar and rank-one `t`. It traces the **shorter great-circle arc** between
unit vectors in two or more dimensions. Fraction values may be unsorted,
repeated, or outside `[0,1]` (extrapolation, consistent with SciPy 1.17+).
It calculates a unit tangent at the initial point and evaluates the arc
with trigonometric functions, improving stability over a raw
`sin(theta)`-denominator formula when the points are close together.

Inputs must be finite, must have the same dimension (at least two), and
must have norms within `1e-7` of unity. Accepted vectors are normalized
internally. The optional `tol` is a **nonnegative antipodal chord-distance
threshold** (default `1e-7`). If the two endpoints are exact or close
antipodes according to `tol`, the shortest arc is ambiguous, so SciFort
returns `spatial_distance_invalid_input=1` and an empty `(0)` or `(0,0)`
result. This intentionally differs from SciPy 1.17+, which may emit a warning
and select an arbitrary path for antipodal inputs. An empty array of times
is otherwise accepted and returns shape `(0,n_dimensions)`.

## Reproducible checks

```sh
python tools/spatial_reference/run_native_spatial.py --mode both
python tools/spatial_reference/compare_geometry_scipy.py
python tools/spatial_reference/compare_geometry_scipy.py --seed 84321
```

The independent numerical oracle is SciPy 1.17.0. The two comparison seeds
cover 840 Hausdorff point-cloud scenarios per build configuration, including
repeated points, ties, one-dimensional clouds, and widely scaled coordinates,
and 720 geometric-slerp scenarios with scalar and vector parameters, near-
coincident endpoints, near-antipodal (but unambiguous) arcs, and extrapolation.

Further geometry: [Procrustes analysis](spatial_procrustes.md).

Other spatial APIs: [distance matrices and metrics](spatial_distance.md),
[KD-trees and spatial queries](spatial_kdtree.md).
