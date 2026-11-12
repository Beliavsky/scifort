# SciFort `scipy.spatial.transform`: 3D rotations and Slerp

SciFort provides an independent, native Fortran subset of the
[`scipy.spatial.transform`](https://docs.scipy.org/doc/scipy-1.17.0/reference/spatial.transform.html)
package in the module `scifort_spatial_transform`. There are no additional
Fortran library dependencies. `real(dp)` is binary64.

## Rotations

`rotation_type` represents an **active, right-handed 3D rotation**. Internally,
rotations are stored as normalized quaternions in scalar-last order
`[x,y,z,w]`, the default ordering for SciPy `Rotation`.

```fortran
use scifort_kinds, only : dp
use scifort_spatial_transform, only : rotation_type, &
    rotation_from_rotvec, rotation_from_quat, rotation_from_matrix, &
    rotation_from_euler, rotation_from_mrp, rotation_from_davenport, &
    rotation_identity
implicit none
type(rotation_type) :: rz, rx, combined, restored
real(dp) :: v(3), q(4), matrix(3,3)
integer :: status

rz = rotation_from_rotvec([0.0_dp,0.0_dp,90.0_dp], &
                          degrees=.true., status=status)
rx = rotation_from_quat([1.0_dp,0.0_dp,0.0_dp,1.0_dp])
combined = rz%compose(rx)             ! Apply rx, then rz.
v = combined%apply([1.0_dp,0.0_dp,0.0_dp])
matrix = combined%as_matrix()
q = combined%as_quat(canonical=.true.)
restored = rotation_from_matrix(matrix,status=status)
```

| Public procedure/method | Meaning |
|---|---|
| `rotation_identity()` | Identity rotation |
| `rotation_from_quat(quat,scalar_first=,status=)` | Normalize scalar-last quaternion, or scalar-first `[w,x,y,z]` |
| `rotation_from_rotvec(rotvec,degrees=,status=)` | Axis times rotation angle (radians by default) |
| `rotation_from_matrix(matrix,status=)` | Convert a proper orthogonal 3x3 matrix |
| `rotation_from_euler(seq,angles,degrees=,status=)` | Active Euler rotation, one to three intrinsic or extrinsic axes |
| `rotation_from_mrp(mrp,status=)` | Rotation from three modified Rodrigues parameters |
| `rotation_from_davenport(axes,order,angles,degrees=,status=)` | Generalized Euler rotations about 1–3 user-specified Cartesian axes |
| `r%as_quat(canonical=,scalar_first=)` | Return quaternion, optionally canonicalized |
| `r%as_matrix()` | Return 3x3 active rotation matrix |
| `r%as_rotvec(degrees=)` | Principal axis-angle vector, magnitude in [0,pi] radians |
| `r%as_euler(seq,degrees=,gimbal_lock=,status=)` | Extract 3-axis Euler angles; optional gimbal-lock flag |
| `r%as_mrp()` | Principal modified Rodrigues vector |
| `r%as_davenport(axes,order,degrees=,gimbal_lock=,status=)` | Decompose a rotation about 3 user-specified Davenport axes |
| `r%apply(vector,inverse=)` | Rotate a single 3D Cartesian vector |
| `r%inverse()` | Inverse rotation |
| `r%compose(other)` | Apply `other` first, then `r` (SciPy `r * other`) |
| `r%magnitude()` | Principal rotation angle in radians |
| `r%valid()` | Check whether a constructor succeeded |

Quaternions are normalized with scaled arithmetic before any rotation is
created, so otherwise finite quaternions on the order of `1e-300` or `1e300`
work. `q` and `-q` represent the same rotation; `canonical=.true.` chooses
one representative. Input errors set `spatial_rotation_invalid_input=1`,
return an invalid `rotation_type`, and yield NaNs from real-valued accessors.
Success is `spatial_rotation_success=0`. Status arguments are optional.

`rotation_from_matrix` deliberately requires an already orthogonal matrix
(with determinant +1 to within `1e-8`); **unlike SciPy**, it does not project
arbitrary slightly nonorthogonal matrices onto SO(3). Higher-than-rank-two rotation-array broadcasting is **not yet implemented**.
One-dimensional stacks and two-dimensional rotation grids are documented
separately.


## Euler-angle conventions

`rotation_from_euler(seq, angles, degrees=, status=)` supports **one, two, or
three axes** with no adjacent duplicates. Lowercase axes denote **extrinsic**
rotations about fixed axes; uppercase axes denote **intrinsic** rotations about
moving axes. The cases cannot be mixed. For example, `xyz` applies x first,
y second, z third; `XYZ` applies x first in a moving frame, then y, then z.
These produce different orientations for nonzero angles. The returned type is
an active right-handed rotation, consistent with SciPy's `Rotation` class.

`r%as_euler(seq, degrees=, gimbal_lock=, status=)` requires a **three-character**
valid sequence. It supports the six Tait–Bryan and six proper Euler sequences
in each of the two frames, for **24 three-axis conventions**. The two outer
angles lie in `[-pi,pi]` radians. The middle angle lies in `[-pi/2,pi/2]`
for three different axes, or `[0,pi]` when the first and third axes are equal.

At gimbal lock, the third angle is set to zero, and optional logical output
`gimbal_lock` is true. The first two angles still describe the orientation.
SciFort uses a `1e-7` threshold in the relevant matrix norm, matching the
practical singularity tolerance in the SciPy reference. Unlike SciPy's Python
method, which issues a warning, SciFort reports the condition via the flag.
For invalid sequences or rotations, `status` is 1 and the returned angles
are NaNs.

```fortran
use scifort_spatial_transform, only : rotation_type, rotation_from_euler
! ... declare r, angles, locked, status ...
r = rotation_from_euler('ZYX', [30.0_dp,-20.0_dp,10.0_dp], &
                        degrees=.true., status=status)
angles = r%as_euler('ZYX', degrees=.true., gimbal_lock=locked, status=status)
```

## Davenport angles (generalized Euler rotations)

`rotation_from_davenport(axes,order,angles,degrees=,status=)` follows
[`Rotation.from_davenport`](https://docs.scipy.org/doc/scipy-1.17.0/reference/generated/scipy.spatial.transform.Rotation.from_davenport.html).
Supply a real array `axes(k,3)` containing one axis per row, where `k=1..3`,
and a vector `angles(k)`. Each axis may have any finite nonzero magnitude.
Adjacent axes **must be perpendicular**, but the first and third axes may be
oblique, equal, or antiparallel. `order` is `'e'`/`'extrinsic'` (fixed-frame)
or `'i'`/`'intrinsic'` (moving-frame). The overloaded single-axis form also
accepts a Cartesian vector and one scalar angle.

`r%as_davenport(axes,order,degrees=,gimbal_lock=,status=)` requires exactly
three axes. It returns three principal angles; at gimbal lock, the third is
zero and the logical `gimbal_lock` flag is true. SciPy emits a warning instead
of returning this flag. Just like Euler angles, Davenport representations
may be nonunique; especially when the first/third axes are parallel, angle
triples may differ while representing exactly the same orientation. The
reconstructed rotation, not raw angle equality, is the appropriate check.

```fortran
use scifort_kinds, only : dp
use scifort_spatial_transform, only : rotation_type, rotation_from_davenport
real(dp) :: axes(3,3), theta(3)
type(rotation_type) :: r
integer :: status
logical :: locked
axes(1,:) = [1.0_dp,0.0_dp,0.0_dp]
axes(2,:) = [0.0_dp,1.0_dp,0.0_dp]
axes(3,:) = [0.5_dp,0.0_dp,sqrt(3.0_dp)/2.0_dp]
r = rotation_from_davenport(axes,'extrinsic',[0.2_dp,0.4_dp,-0.1_dp], &
                            status=status)
theta = r%as_davenport(axes,'extrinsic',gimbal_lock=locked,status=status)
```

Batched variants are provided by `rotation_stack_from_davenport`,
`rotation_stack%as_davenport`, `rotation_grid_from_davenport`, and
`rotation_grid%as_davenport`. Input angles are shaped `(n,k)` for stacks,
`(n,m,k)` for grids, with one logical lock flag per result. Input errors yield
status 1, invalid results, and zero-sized output arrays for batched accessors.

## Modified Rodrigues parameters

MRPs are defined by `mrp = axis*tan(angle/4)`. SciFort provides
`rotation_from_mrp(mrp, status=)` and `r%as_mrp()`. Construction accepts any
finite three-vector, including its non-principal *shadow representation*
(norm larger than 1); extraction always returns the principal MRP with norm
at most 1, describing a rotation through at most 180 degrees. For finite
MRPs with extremely large magnitudes, scale normalization avoids overflowing
squared components (similarly for tiny MRPs). As with quaternions, many
parameter values can represent the same physical rotation.

```fortran
use scifort_spatial_transform, only : rotation_from_mrp
! ... declare r, mrp ...
r = rotation_from_mrp([0.0_dp,0.0_dp,0.414213562373095_dp])
mrp = r%as_mrp()  ! 90-degree rotation about z
```


## Weighted alignment and rotation averaging

`rotation_align_vectors(a,b,weights=,rssd=,sensitivity=,status=)` solves
SciPy's `Rotation.align_vectors` weighted 3D Wahba/Kabsch problem. `a` and
`b` are both `(n,3)` arrays with **one vector per row**, and the returned
rotation acts on `b` to approximate `a` without changing vector lengths.
Weights must be positive and finite, except that **one** may be positive
infinity. Omitted weights default to one. An optional `rssd` returns the root
weighted sum of squared residual distances. Requesting `sensitivity(3,3)`
returns SciPy's normalized rotation-error sensitivity matrix when it is
well-defined. Sensitivity is undefined for a single vector, infinite primary
weight, or a singular information matrix: those requests return status 1.

```fortran
use scifort_kinds, only : dp
use scifort_spatial_transform, only : rotation_type, &
    rotation_align_vectors, rotation_mean
! Given real(dp) :: a(n,3), b(n,3), w(n), rssd, sensitivity(3,3)
!       integer :: status; type(rotation_type) :: aligned, average
aligned = rotation_align_vectors(a,b,weights=w, &
    rssd=rssd,sensitivity=sensitivity,status=status)
! For previously constructed rotations(:) of type(rotation_type):
average = rotation_mean(rotations,status=status)
```

For one nonzero vector pair, alignment uses the shortest SO(3) rotation;
near-antipodal pairs use a cancellation-resistant half-angle construction.
For an infinite primary weight, SciFort aligns that pair exactly in direction
and then optimizes the remaining rotation around the primary axis using the
finite-weight secondary pairs. The primary contribution to `rssd` is excluded,
matching SciPy. The output for exactly antipodal or unidentifiable configurations
is not unique and can differ from SciPy without changing the fitting objective.
A zero source or destination vector in the single-vector case is rejected.
For other unidentifiable configurations, SciFort returns status 1 rather than
choosing an arbitrary quaternion. Small residuals are evaluated directly
rather than from differences of large nearly equal squared norms.

`rotation_mean(rotations,weights=,status=)` implements SciPy's **chordal L2**
rotation mean, not an intrinsic/geodesic Karcher mean. It accepts a one-dimensional
array of `rotation_type` and optional finite nonnegative weights. At least one
weight must be positive. Mean computations are invariant to quaternion sign and
normalize large/small weights to avoid overflow. With nonunique minimizers,
the returned mean is deterministic but need not have SciPy's quaternion sign
or exact choice of representative orientation. Like SciPy, these algorithms
operate on orientation data only; vector magnitudes are not normalized for
`rotation_align_vectors`, so they contribute implicit statistical weight.

Internally, both estimators use original, fixed-size symmetric 4x4 Jacobi
iterations: Davenport's quaternion objective for alignment and Markley's
quaternion outer-product objective for the chordal mean. No BLAS, LAPACK, or
other external source code is required. The sensitivity uses the inverse
symmetric Hessian of the weighted alignment objective and the arithmetic mean
of the observation weights, matching SciPy's normalization.

## Spherical interpolation of orientations

`rotation_slerp(start,finish,fraction,status=)` interpolates on the shortest
path between two orientations. The fraction is any finite number; values
outside `[0,1]` extrapolate as a SciFort convenience extension.

```fortran
use scifort_spatial_transform, only : rotation_type, &
    rotation_from_rotvec, rotation_slerp, rotation_slerp_path
! ... declare and initialize variables ...
r_mid = rotation_slerp(r1,r2,0.5_dp)
call rotation_slerp_path(key_times, key_rotations, &
    sample_times, interpolated, status)
```

`rotation_slerp_path` follows SciPy's `Slerp` convention:

- `key_times(:)` are finite, strictly increasing, and contain at least 2 points.
- `key_rotations(:)` are `rotation_type`, one per keyframe.
- `sample_times(:)` may be unsorted or repeated but must all be in the
  inclusive interval `[key_times(1),key_times(n)]`.
- `interpolated(:)` is an allocatable array of `rotation_type`, allocated by
  the procedure. Invalid calls yield status 1 and an empty array.
- The interpolation is piecewise geodesic with constant angular velocity
  within each segment, not a globally smooth rotation spline.

For `slerp`, the quaternion sign of the second endpoint is flipped where
necessary to take the shortest orientation path. Nearly identical rotations
use a normalized linear blend to reduce cancellation. The definition follows
standard quaternion/SO(3) geometry; see Ken Shoemake,
*Animating rotation with quaternion curves*, SIGGRAPH 1985,
https://doi.org/10.1145/325334.325242 .

This API is **different** from the existing `scifort_spatial%geometric_slerp`,
which interpolates between points on a unit sphere in arbitrary dimensions
rather than between rigid-body orientations.

## Tests and known limitations

```text
python tools/spatial_reference/run_native_spatial.py --mode both
python tools/spatial_reference/compare_rotation_scipy.py --mode both
python tools/spatial_reference/compare_rotation_euler_mrp_scipy.py --mode both
python tools/spatial_reference/compare_rotation_estimation_scipy.py --mode both
```

`compare_rotation_scipy.py` checks quaternion and rotation-vector conversions,
matrices, vector application and inverse application, composition, principal
angle, two-keyframe Slerp, and nonuniform multi-keyframe Slerp against SciPy.
It compares quaternions modulo sign. The reference normalizes extremely
scaled quaternions first, because the SciPy 1.17.0 implementation may overflow
or underflow when normalizing unscaled inputs.

The scalar `rotation_type%apply` handles one 3-vector at a time. The separate
[`scifort_spatial_rotation_stack`](spatial_rotation_stack.md) now supports
one-dimensional batched rotations, vector broadcasting, composition, and
stack conversions; full NumPy higher-rank broadcasting remains unsupported. Trigonometric computations at extraordinarily
large rotation-vector angles can be limited by floating-point argument
reduction; inputs must be finite. There is no SO(3) projection of noisy
matrices in the scalar constructor; cubic `RotationSpline` is provided in a
separate module (see below).

## Smooth angular-rate and acceleration interpolation

The separate [`scifort_spatial_rotation_spline`](spatial_rotation_spline.md)
module provides `rotation_spline` and `rotation_spline_type`, implementing
`scipy.spatial.transform.RotationSpline`-style cubic SO(3) interpolation with
continuous angular rate and acceleration. It builds on the rotation types
described here, without changing the existing Slerp API.

## Batched rotation-stack interface

For one-dimensional collections, use the separate
[`scifort_spatial_rotation_stack`](spatial_rotation_stack.md) module. It builds
on this scalar quaternion API without changing existing call signatures.
