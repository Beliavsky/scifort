# SciFort batched 3D rotations and rotation stacks

`scifort_spatial_rotation_stack` provides a one-dimensional stack of 3D active
rotations analogous to SciPy's vectorized `scipy.spatial.transform.Rotation`.
It builds on the scalar [`rotation_type`](spatial_transform.md) rather than
changing that API. The stack is a derived type with allocatable rotation storage,
so it does not require a fixed compile-time size or BLAS/LAPACK.

## Construction and representation

| Constructor | Input shape | Result |
| --- | --- | --- |
| `rotation_stack_from_quat(quats, scalar_first=, status=)` | `(n,4)` | Quaternion stack |
| `rotation_stack_from_rotvec(rotvecs, degrees=, status=)` | `(n,3)` | Rotation-vector stack |
| `rotation_stack_from_mrp(mrps, status=)` | `(n,3)` | Modified Rodrigues parameter stack |
| `rotation_stack_from_euler(seq, angles, degrees=, status=)` | `(n,len(seq))` | Euler rotation stack |
| `rotation_stack_from_matrix(matrices, status=)` | `(n,3,3)` | Matrix rotation stack |
| `rotation_stack_from_rotations(rotations, status=)` | `rotation_type(n)` | Stack existing rotations |
| `rotation_stack_identity(n, status=)` | `n >= 0` | `n` identity rotations |
| `rotation_stack_random(n, rng=, status=)` | `n >= 0` | Haar-uniform random SO(3) rotations |

Every constructor validates the complete collection. Invalid input returns an
undefined stack (`stack%valid() == .false.`) and status 1; success returns 0.
There is no partially initialized result. Constructors other than `identity`
require at least one rotation.

On a valid stack, call `as_quat(canonical=,scalar_first=)`, `as_rotvec(degrees=)`,
`as_mrp()`, `as_matrix()`, or `as_euler(seq,degrees=,gimbal_lock=,status=)`.
The output shapes are `(n,4)`, `(n,3)`, `(n,3)`, `(n,3,3)`, and `(n,3)`.
Euler gimbal-lock flags are allocatable `logical(:)` and are optionally returned.
Unlike NumPy's flexible multidimensional broadcasting, this module represents
one-dimensional stacks only. Indices are **one-based**. Arrays are represented
in the same *logical* row orientation as SciPy, although Fortran uses
column-major physical storage.

## Batched operations and broadcasting

The following methods preserve scalar-rotation conventions:

- `stack%apply(v, inverse=, status=)` with `v(3)` applies each of the `n`
  rotations to one vector and returns `(n,3)`.
- `stack%apply(vectors, inverse=, status=)` with `vectors(m,3)` returns
  `(max(n,m),3)` if `n==m`, `n==1`, or `m==1`; otherwise status 1 and `(0,3)`.
- `stack%compose(other, status=)` performs **left-after-right** composition,
  broadcasting a `rotation_stack_type` or a single `rotation_type` on the right.
  Reversed composition is generally a different rotation.
- `stack%inverse()`, `stack%power(exponent)`, and `stack%magnitude()` return
  inverse rotations, quaternion-based rotation powers, and principal angles.
- `stack%mean(weights=,status=)` uses the existing weighted chordal L2 mean.
- `stack%approx_equal(other, atol=, degrees=, status=)` returns per-rotation
  logical comparisons, with the same size-one broadcasting, using geodesic
  angular tolerance (`atol` defaults to `1e-8` radians).
- `stack%at(i,status=)` extracts one scalar rotation; `stack%slice(indices,status=)`
  extracts/reorders arbitrary one-based indices (including duplicates).
- `rotation_stack_concatenate([stack_a,stack_b],status=)` joins stacks without
  performing rotation composition.

`power` applies the exponent to the principal rotation vector and may
choose a different quaternion sign than SciPy while representing the same SO(3)
rotation. Nonfinite exponents and components, incompatible broadcast lengths,
and invalid indices are rejected. As a convenient Fortran extension,
`rotation_stack_identity(0)` and empty slices are valid zero-length stacks;
nonempty pairwise operations do not accept zero-length stacks.

## Example

```fortran
use scifort_kinds, only : dp
use scifort_spatial_rotation_stack, only : &
    rotation_stack_type, rotation_stack_from_euler

type(rotation_stack_type) :: orientations, inverse_orientations
real(dp), allocatable :: directions(:,:), quaternions(:,:)
integer :: status

orientations = rotation_stack_from_euler('ZYX', &
    reshape([0.0_dp,0.0_dp,0.0_dp, &
             45.0_dp,0.0_dp,0.0_dp], [2,3], order=[2,1]), &
    degrees=.true., status=status)
if (status /= 0) error stop 'Invalid orientations'

directions = orientations%apply([1.0_dp,0.0_dp,0.0_dp])
quaternions = orientations%as_quat(canonical=.true.)
inverse_orientations = orientations%inverse()
```

## Regression tests

```text
python tools/spatial_reference/run_native_spatial.py --mode both
python tools/spatial_reference/compare_rotation_stack_scipy.py --mode both
```

`test/test_spatial_rotation_stack.f90` checks construction, conversions,
rotation application, both broadcasting directions, order of composition,
inversion, repeated/unsorted indices, concatenation, powers, means, and invalid
shapes. The independent Python comparison tests randomized quaternion scales,
near-180-degree rotations, quaternion signs, vector broadcasting,
representations, composition, rotation powers, and weighted means against SciPy.

The implementation performs no implicit orientation repair for noisy input
matrices; like scalar `rotation_from_matrix`, it requires an approximately
proper orthogonal 3x3 rotation. The companion [`scifort_spatial_rotation_grid`](spatial_rotation_grid.md)
module supports two batch axes with independent broadcasting and adds
`rotation_random`, `rotation_stack_random`, and `rotation_grid_random` with
caller-owned reproducible streams. Full arbitrary-rank NumPy rotation arrays
and every SciPy broadcasting combination remain future extensions.
Finite SO(3) symmetry groups and minimum-angle reduction are available in
[`scifort_spatial_rotation_group`](spatial_rotation_group.md).

## Davenport-angle conversions

`rotation_stack_from_davenport(axes,order,angles,degrees=,status=)` constructs
rotations from generalized Euler angles. `axes(k,3)` contains one Cartesian
axis per row, with consecutive axes perpendicular. The angle array has
shape `(n,k)` for stacks or `(n,m,k)` for grids.

`stack%as_davenport(axes,order,...)` returns `(n,3)` angles and an
optional `(n)` gimbal-lock flag array.

See [Davenport angle documentation](spatial_transform.md#davenport-angles-generalized-euler-rotations).
