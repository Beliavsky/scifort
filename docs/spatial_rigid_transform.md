# Rigid transformations in SE(3)

`scifort_spatial_rigid_transform` implements an independent modern Fortran
subset of [`scipy.spatial.transform.RigidTransform`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.transform.RigidTransform.html),
using the pre-existing SciFort quaternion `rotation_type`.

A transform has an active rotation `R` and a translation `t`, and maps a
Cartesian point `x` to `R*x + t`. The order of composition is the same as
SciPy: `A%compose(B)` applies `B` **first**, then `A`.

## Scalar API

```fortran
use scifort_kinds, only : dp
use scifort_spatial_transform, only : rotation_type, rotation_from_euler
use scifort_spatial_rigid_transform, only : &
    rigid_transform_type, rigid_from_components, rigid_from_matrix, &
    rigid_from_exp_coords, rigid_from_dual_quat, rigid_identity

type(rotation_type) :: r
type(rigid_transform_type) :: transform, middle, reconstructed
real(dp) :: point(3), moved(3), homogeneous(4,4), coordinates(6), dual(8)
integer :: status

r = rotation_from_euler('z', [90.0_dp], degrees=.true.)
transform = rigid_from_components([2.0_dp, 0.0_dp, 1.0_dp], r, status)
if (status /= 0) error stop 'Invalid transform'

point = [1.0_dp, 0.0_dp, 0.0_dp]
moved = transform%apply(point) ! approximately [2,1,1]
homogeneous = transform%as_matrix()
coordinates = transform%as_exp_coords()
dual = transform%as_dual_quat()
reconstructed = rigid_from_dual_quat(dual)
middle = transform%power(0.5_dp)
```

Factories (`status` is optional, `0` on success and `1` on error):

| Factory | Input |
|---|---|
| `rigid_identity()` | SE(3) identity |
| `rigid_from_translation(t, status)` | Length-3 Cartesian translation |
| `rigid_from_rotation(r, status)` | `rotation_type`, zero translation |
| `rigid_from_components(t, r, status)` | Apply `r` first, then add `t` |
| `rigid_from_matrix(m, status)` | Proper 4x4 homogeneous matrix |
| `rigid_from_exp_coords(c, status)` | Twist `(rx,ry,rz,vx,vy,vz)` |
| `rigid_from_dual_quat(q, scalar_first, status)` | Unit dual quaternion, real + dual parts |

Methods:

- `%valid()`, `%translation()`, `%rotation()`
- `%as_matrix()`, `%as_exp_coords()`, `%as_dual_quat(scalar_first=)`
- `%apply(vector, inverse=)` for a length-3 vector, or
  `%apply(points, inverse=)` for points of shape `(n,3)`
- `%inverse()`, `%compose(other)`, `%power(exponent, status=)`.

Exponential-coordinate translation uses the **left SO(3) Jacobian**, with
Taylor series for angles below `1e-4` radians. The inverse log map uses
its inverse Jacobian, avoiding catastrophic cancellation near zero. The
rotation logarithm uses the principal representative with angle in
`[0,pi]`. Consequently, a screw transform raised to a fractional power
follows the principal (shortest) rotational branch.

Dual quaternions are stored in scalar-last order `(qx,qy,qz,qw, dx,dy,dz,dw)`;
`scalar_first=.true.` instead returns `(qw,qx,qy,qz,dw,dx,dy,dz)`.
The input real quaternion is normalized, and the orthogonality component
of the dual part is projected away before recovering the translation.

## One-dimensional transform stacks

`rigid_transform_stack_type` provides a one-dimensional batch interface.

| Factory | Shape |
|---|---|
| `rigid_stack_from_components(translations, rotations, status)` | `(n,3)` translations, `rotation_type(:)` rotations; length `n` or `1` |
| `rigid_stack_from_matrix(matrices, status)` | `(n,4,4)` |
| `rigid_stack_from_exp_coords(coordinates, status)` | `(n,6)` |
| `rigid_stack_from_dual_quat(duals, scalar_first, status)` | `(n,8)` |

Methods include `%valid()`, `%size()`, `%at(i)` (one-based), `%as_matrix()`,
`%translation()`, `%as_exp_coords()`, `%as_dual_quat()`, `%inverse()`,
`%compose(other)` (stack or scalar, with size-one broadcasting),
`%power(real_exponent)`, `%apply(vector)` (result `(n,3)`),
`%apply(points)` (point count `n` or `1`, or a singleton stack), and
`%mean(weights=, status=)`.

The mean uses the **weighted arithmetic translation mean** and the existing
**weighted chordal L2 rotation mean**. Weights must be finite and
nonnegative, with at least one positive entry.

## Limitations and safety

- The 4x4 matrix factory **validates** an already orthogonal proper rotation
  (`rotation_from_matrix` tolerance `1e-8`) and a homogeneous bottom row
  (tolerance `1e-10`). It does **not** project arbitrary 3x3 matrices onto
  SO(3) by singular-value decomposition as SciPy does.
- Scalar transforms and **one-dimensional** stacks are supported; full
  NumPy-style arbitrary-rank broadcasting, indexing and mean reduction
  axes are not yet supported.
- Invalid factories return `%valid() == .false.` and nonzero status. Invalid
  vector/matrix conversions produce NaN values. Invalid stack operations
  return undefined stacks. Empty stacks are unsupported.
- Very large translation components may overflow intermediate SE(3)
  operations; the library checks most such results for finite values.
- SciPy 1.17's backend implementation of `RigidTransform.__pow__` can differ
  at roughly `1e-7` from the documented exp/log identity for noninteger
  exponents. SciFort directly uses the mathematical
  `from_exp_coords(as_exp_coords()*n)` contract and verifies this identity
  against SciPy's explicit exp/log methods. This behavior is intentional.

## Test commands

```text
python tools/spatial_reference/run_native_spatial.py --mode both
python tools/spatial_reference/compare_rigid_transform_scipy.py --mode both
```

The new native regression exercises identity, round trips through each
representation, integer/fractional powers, exact compositions, batched point
application, stack broadcasting, weighted means, and validation failures.
The independent comparator checks two randomized seeds against SciPy 1.17.0
in checked and optimized GNU Fortran modes, including small-angle and
near-180-degree rotations. A second independent driver tests 600 arbitrary
SE(3) twists, including rotation angles near multiples of `2*pi`.

All code was written independently for SciFort, under the repository MIT
license. No third-party implementation or new dependency is incorporated.
