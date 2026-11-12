# Rank-two rigid-transform grids (SE(3))

`scifort_spatial_rigid_transform_grid` adds two-batch-axis arrays of rigid
transformations, corresponding to the rank-two subset of
[`scipy.spatial.transform.RigidTransform`](https://docs.scipy.org/doc/scipy-1.17.0/reference/generated/scipy.spatial.transform.RigidTransform.html).
It reuses SciFort's native rotation and SE(3) implementations, and does not
introduce any runtime dependencies. All array dimensions are in Fortran order;
component dimensions are last, and indices are **one-based**.

```fortran
use scifort_kinds, only : dp
use scifort_spatial_rigid_transform_grid, only : &
    rigid_transform_grid_type, rigid_grid_from_exp_coords
use scifort_spatial_rigid_transform, only : rigid_transform_stack_type

type(rigid_transform_grid_type) :: poses, inverse_poses
type(rigid_transform_stack_type) :: mean_over_first_axis
real(dp) :: twists(2,3,6)
real(dp), allocatable :: transformed(:,:,:), matrices(:,:,:,:)

! Initialize twists(:,:,1:3) as rotation vectors in radians,
! and twists(:,:,4:6) as body-fixed translational exponential coordinates.
poses = rigid_grid_from_exp_coords(twists)
matrices = poses%as_matrix()                ! shape (2,3,4,4)
transformed = poses%apply([1.0_dp,0.0_dp,0.0_dp])  ! (2,3,3)
inverse_poses = poses%inverse()
mean_over_first_axis = poses%mean_axis(1) ! stack of length 3
```

## Factories

| Factory | Input dimensions | Notes |
|---|---|---|
| `rigid_grid_identity(dims, status)` | `dims = [n,m]` | Positive dimensions required |
| `rigid_grid_from_transforms(poses, status)` | `(n,m)` `rigid_transform_type` | Explicit collection |
| `rigid_grid_from_translation(t, status)` | `(n,m,3)` | Identity rotation |
| `rigid_grid_from_rotation(rgrid, status)` | `rotation_grid_type` | Zero translation |
| `rigid_grid_from_components(t, rotations, status)` | `t(n,m,3)`, rotation array `(n,m)` or `rotation_grid_type` | Independent axis-wise broadcasting |
| `rigid_grid_from_matrix(matrices, status)` | `(n,m,4,4)` | Proper homogeneous matrices |
| `rigid_grid_from_exp_coords(twists, status)` | `(n,m,6)` | Rotation vector followed by translational twist |
| `rigid_grid_from_dual_quat(dq, scalar_first, status)` | `(n,m,8)` | Two quaternion components, optionally scalar-first |

All factories accept an optional integer `status` (0 for success, 1 for
invalid input). An invalid grid has `%valid()==.false.`.

## Methods

- `%valid()`, `%shape()` return validity and `(n,m)` dimensions.
- `%at(i,j)` selects a scalar rigid transform; `%row(i)` and `%column(j)`
  return a one-dimensional `rigid_transform_stack_type`.
- `%translation()`, `%as_matrix()`, `%as_exp_coords()`, and
  `%as_dual_quat(scalar_first=)` return component-last real arrays.
- `%inverse()` and `%power(exponent, status=)` act elementwise.
- `%compose(other)` accepts a scalar transform, rank-one transform stack,
  or another rank-two grid. Composition is **right-to-left**: the other
  transform is applied first. Each of the two batch dimensions can broadcast
  when its size is one. Rank-one stacks align with the *second* grid axis.
- `%apply(points, inverse=)` accepts a common point `(3)`, a rank-one set
  `(k,3)` aligned with the second batch axis, or a rank-two array `(a,b,3)`.
  The transform and point axes broadcast independently.
- `%mean(weights=, status=)` returns a single SE(3) transform, averaging all
  grid elements, while `%mean_axis(axis, weights=, status=)` returns a
  one-dimensional transform stack: axis 1 yields length `m`, axis 2 yields
  length `n`. The translation mean is arithmetic, while the rotation mean is
  the weighted **chordal L2** mean inherited from `scifort_spatial_transform`.
  Nonnegative finite weight matrices can broadcast independently on either
  batch axis. Every reduced slice must contain a strictly positive weight.

## SciPy compatibility notes

This implements rank-two broadcasting from SciPy 1.17, not arbitrary NumPy
rank-N transform arrays or general-purpose Python indexing. Shape mismatch
or invalid parameters return undefined grids or zero-length outputs instead of
raising Python exceptions. All indexing is one-based; `axis=1` here is
`axis=0` in NumPy/SciPy, and `axis=2` here is `axis=1` in NumPy/SciPy.

Rigid-transform matrix construction **requires** a proper rotation matrix;
it does not project an approximately orthogonal input onto SO(3) as SciPy may.
Fractional powers use the exact mathematical exp/log contract, which can
differ slightly from SciPy's general `RigidTransform.__pow__` implementation.

## Tests

```text
python tools/spatial_reference/run_native_spatial.py --mode both
python tools/spatial_reference/compare_rigid_grid_scipy.py --mode both
```

The native regression checks one-based indexing, unit and incompatible
broadcasts, round trips through all representations, inverse and fractional
powers, per-axis/global means, invalid weights and shapes, and transform-stack
extraction. The reference driver independently checks randomized translations,
rotations, point arrays, covariance-independent weighted averaging, and
broadcast combinations against SciPy 1.17.0.

All implementation, test, and documentation code is independent MIT-licensed
SciFort work. No SciPy implementation source was copied.
