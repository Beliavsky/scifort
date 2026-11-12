# Arbitrary-rank rotation and rigid-transform batches

The module `scifort_spatial_transform_nd` extends SciFort's scalar, stack, and
rank-two rotation/rigid-transform interfaces with **runtime-rank batch shapes**.
The input component dimension is **last**, as in `scipy.spatial.transform`.
The module is written independently in Fortran and adds no dependencies.

Fortran array ranks are fixed at compile time, so values returned by the new
API use **flattened component arrays** plus separately available batch shapes:

- `%shape()` returns the batch shape, with `[]` denoting a scalar transform.
- `%as_quat()` returns `(product(shape),4)`; `%as_rotvec()` returns `(n,3)`;
  `%as_matrix()` returns `(n,3,3)` for rotations or `(n,4,4)` for rigid transforms.
- `%apply(points, result_shape=...)` returns `(n,3)` and optionally the
  *broadcasted* batch shape through `result_shape` (allocatable integer array).
- `%at(indices)` uses **one-based** multi-indices; flat data are ordered with
  the first Fortran batch dimension varying fastest.
- `status=0` means success; invalid inputs return status `1` and an invalid
  object (or a zero-row numerical result for `%apply`).

## Example: three batch axes

```fortran
use scifort_kinds, only : dp
use scifort_spatial_transform_nd, only : &
    rotation_nd_type, rigid_nd_type, &
    rotation_nd_from_quat, rigid_nd_from_components

type(rotation_nd_type) :: a, b, combined
type(rigid_nd_type) :: poses
real(dp) :: qa(2,1,3,4), qb(1,4,1,4), translations(2,1,3,3)
real(dp) :: vectors(1,4,1,3)
real(dp), allocatable :: points(:,:), quats(:,:)
integer, allocatable :: batch_shape(:)

! Initialize qa/qb with finite nonzero quaternions and translations/vectors.
a = rotation_nd_from_quat(qa)        ! batch shape (2,1,3)
b = rotation_nd_from_quat(qb)        ! batch shape (1,4,1)
combined = a%compose(b)             ! broadcast shape (2,4,3)
quats = combined%as_quat()          ! shape (24,4)
points = combined%apply(vectors, result_shape=batch_shape) ! (24,3)
poses = rigid_nd_from_components(translations,a)
```

`reshape(quats, [2,4,3,4])` reconstructs a component-last quaternion tensor.
Similarly, `reshape(points, [batch_shape, 3])` reconstructs the vector array.

## Public constructors

| Constructor | Input | Notes |
|---|---|---|
| `rotation_nd_from_quat(q, scalar_first, status)` | `(*batch,4)` | Normalizes quaternions |
| `rotation_nd_from_rotvec(v, degrees, status)` | `(*batch,3)` | Radians by default |
| `rotation_nd_from_matrix(m, status)` | `(*batch,3,3)` | Proper rotation matrices |
| `rotation_nd_from_rotations(r, batch_shape, status)` | Flat scalar rotations and shape | Interop with existing API |
| `rotation_nd_identity(batch_shape, status)` | Integer batch shape | Identity rotations |
| `rigid_nd_from_matrix(m, status)` | `(*batch,4,4)` | Proper homogeneous transforms |
| `rigid_nd_from_exp_coords(x, status)` | `(*batch,6)` | SE(3) exponential coordinates |
| `rigid_nd_from_dual_quat(q, scalar_first, status)` | `(*batch,8)` | Real/dual quaternion pairs |
| `rigid_nd_from_components(t, r, status)` | `t(*batch,3)` + rotation ND | Broadcasts translation and rotation batches |
| `rigid_nd_from_transforms(tf, batch_shape, status)` | Flat scalar transforms and shape | Interop with existing API |
| `rigid_nd_identity(batch_shape, status)` | Integer batch shape | Identity transforms |

Constructor arrays are **assumed-rank** and use `SELECT RANK` for ranks 1–15
(vector-like data) or 2–15 (matrices), corresponding to at most 14 or 13
batch axes respectively. The data are copied into a rank-independent
representation. Every *non-scalar* batch dimension must be positive.

`rotation_nd_type` provides `valid`, `size`, `shape`, `at`, `as_quat`,
`as_rotvec`, `as_matrix`, `inverse`, `compose`, and `apply`.
`rigid_nd_type` provides `valid`, `size`, `shape`, `at`, `translation`,
`as_matrix`, `as_exp_coords`, `as_dual_quat`, `inverse`, `compose`, `apply`,
and `power` (uniform scalar exponent).

## Broadcasting rules

Shape matching follows NumPy's **right-aligned** rules: dimensions agree,
or one operand has extent one. A missing leading dimension behaves as one;
array access, flattening, and indices still follow **Fortran order**.
Broadcasting applies to rotation and SE(3) composition, point/vector
application, and construction of rigid transforms from rotation/translation
arrays. No hidden expansion is performed when shapes are incompatible.

## Differences and limitations

The scalar, stack, and rank-two grid APIs remain unchanged. This separate
module is for arbitrary batch ranks; it does not return dynamic-rank Fortran
arrays, or implement NumPy-style slicing, advanced indexing, stack reduction,
or empty batch dimensions. For these operations use the existing fixed-rank
interfaces or reshape the flattened results explicitly. As before, matrix
constructors require proper orthogonal rotation blocks rather than projecting
imperfect matrices onto SO(3) by SVD. The implementation supports the usual
scalar-first option for quaternion and dual-quaternion I/O.

## Validation

Run native tests using
`python tools/spatial_reference/run_native_spatial.py --mode checked`.
The independent SciPy numerical comparisons use
`python tools/spatial_reference/compare_rotation_nd_scipy.py --mode both`.
The reference driver compares composition, point application, and homogeneous
matrices under randomized **four-axis** broadcasting, including singleton
axes on either operand; the native suite additionally tests rank-15 inputs,
constructor round trips, inverse/powers, and invalid input handling.
