# SciFort Haar-uniform random rotations and rank-two rotation arrays

The MIT-licensed `scifort_spatial_rotation_stack` and
`scifort_spatial_rotation_grid` modules extend SciFort's implementation of
`scipy.spatial.transform.Rotation`. These APIs add random rotations and a
**two-dimensional batch shape**, while leaving existing scalar rotations and
one-dimensional rotation stacks intact. They require no BLAS, LAPACK, or
third-party runtime dependency.

## Haar-uniform random rotations

| Constructor | Result | Details |
| --- | --- | --- |
| `rotation_rng(seed)` | `rotation_rng_type` | Reproducible caller-owned stream; `seed` has kind `int64` |
| `rotation_random(rng=,status=)` | `rotation_type` | One random 3D orientation |
| `rotation_stack_random(n,rng=,status=)` | `rotation_stack_type` | `n >= 0` random orientations |
| `rotation_grid_random([n,m],rng=,status=)` | `rotation_grid_type` | Rank-two grid of random rotations; `n,m >= 1` |

The Shoemake quaternion transform maps three uniform variates to the invariant
(Haar) probability measure on SO(3). It does **not** produce uniform Euler
angles or independent uniformly distributed quaternion coordinates.

When `rng` is omitted, generation uses the compiler's intrinsic
`random_number` stream. For reproducible, independent streams, create a local
`rotation_rng_type` and pass it as an `intent(inout)` argument; the same stream
can be advanced by scalar, stack and grid calls without modifying the
process-wide intrinsic random state. A grid traverses its **first** Fortran
batch axis fastest. Success status is 0; invalid inputs return status 1 and an
invalid rotation or grid.

The explicit stream uses a Park–Miller 31-bit integer generator and combines
two updates for each floating-point uniform variate. Its integer arithmetic is
overflow-free on signed 64-bit integers, making the sequence independent of
compiler `random_seed` implementations. The generator is intended for
**repeatable simulation**, not cryptography or extremely large Monte Carlo
runs; its underlying period is `2147483646`. In particular, SciFort and
NumPy/SciPy **do not produce identical sequences from the same integer seed**.

### Example

```fortran
use, intrinsic :: iso_fortran_env, only : int64
use scifort_spatial_rotation_stack, only : rotation_rng, rotation_rng_type, &
    rotation_stack_type, rotation_stack_random
use scifort_spatial_rotation_grid, only : rotation_grid_type, rotation_grid_random

type(rotation_rng_type) :: rng
type(rotation_stack_type) :: sample
type(rotation_grid_type) :: grid
real(dp), allocatable :: rotated(:,:,:)

rng = rotation_rng(20261102_int64)
sample = rotation_stack_random(100, rng)
grid = rotation_grid_random([2,3], rng)
rotated = grid%apply([1.0_dp,0.0_dp,0.0_dp])
```

In complete code, also `use scifort_kinds, only : dp`.

## Rank-two rotation grids

`rotation_grid_type` stores an `(n,m)` array of SO(3) rotations. This
corresponds to the two-leading-batch-axis subset of SciPy >=1.17's
N-dimensional `Rotation`. All arrays in the Fortran API use **one-based
indexing** and component dimensions last.

| API | Input or result shape |
| --- | --- |
| `rotation_grid_from_quat(q, scalar_first=, status=)` | `q(n,m,4)` |
| `rotation_grid_from_rotvec(v, degrees=, status=)` | `v(n,m,3)` |
| `rotation_grid_from_mrp(v, status=)` | `v(n,m,3)` |
| `rotation_grid_from_euler(seq, angles, degrees=, status=)` | `angles(n,m,len(seq))`, 1–3 axes |
| `rotation_grid_from_matrix(matrices, status=)` | `matrices(n,m,3,3)` |
| `rotation_grid_from_rotations(rotations, status=)` | `rotation_type(n,m)` |
| `grid%shape()` | Integer `[n,m]` |
| `grid%at(i,j,status=)` | One `rotation_type` |
| `grid%as_quat(canonical=,scalar_first=)` | `(n,m,4)` |
| `grid%as_rotvec(degrees=)` | `(n,m,3)` |
| `grid%as_mrp()` | `(n,m,3)` |
| `grid%as_matrix()` | `(n,m,3,3)` |
| `grid%as_euler(seq,degrees=,gimbal_lock=,status=)` | `(n,m,3)`, optional logical lock grid `(n,m)` |
| `grid%magnitude()` | Principal angles `(n,m)` |
| `grid%inverse()` | `(n,m)` rotation grid |
| `grid%mean(axis,weights=,status=)` | 1D `rotation_stack_type` along remaining axis |

`grid%apply(vector, inverse=,status=)` broadcasts a single three-component
vector and returns `(n,m,3)`; alternatively,
`grid%apply(vectors, inverse=,status=)` broadcasts vectors of shape
`(k,l,3)` to the result shape `(max(n,k),max(m,l),3)`.
`grid%compose(other,status=)` uses active left-after-right composition and
broadcasts a second grid by comparing each batch axis independently.
`compose` also accepts a single scalar `rotation_type` on the right.

In both grid-grid composition and grid-vector application, each pair of axis
lengths must be equal or one of the two must be 1. Invalid input or
incompatible shape returns status 1 and an invalid grid or an empty array;
constructors do not return partially initialized results. Zero-sized grid
batch axes are not yet supported.

Weighted `grid%mean(axis=1,weights=w)` reduces the first Fortran axis and
returns `m` rotations; `axis=2` reduces the second and returns `n` rotations.
Only weights shared across slices (a 1D vector of length equal to the reduced
axis) are currently supported. Because Fortran axes are one-based, these
correspond to SciPy `axis=0` and `axis=1`, respectively. With SciPy's
N-dimensional weighted mean, the reference comparison broadcasts the weights
across the unreduced axis (e.g., `weights[:,None]` for `axis=0`).

## Example: broadcasting

```fortran
use scifort_kinds, only : dp
use scifort_spatial_rotation_grid, only : rotation_grid_from_rotvec, rotation_grid_type

type(rotation_grid_type) :: orientations, composed
real(dp) :: rotvecs(2,3,3), vectors(2,1,3)
real(dp), allocatable :: rotated(:,:,:)
integer :: status

! Assume rotvecs and vectors have been initialized.
orientations = rotation_grid_from_rotvec(rotvecs,status=status)
if (status /= 0) error stop 'Invalid rotations'
rotated = orientations%apply(vectors,status=status)
if (status /= 0) error stop 'Invalid broadcast shape'
composed = orientations%compose(orientations%inverse())
```

## Numerical comparisons and limitations

```text
python tools/spatial_reference/run_native_spatial.py --mode both
python tools/spatial_reference/compare_rotation_grid_scipy.py --mode both
```

The independent reference driver checks quaternion, matrix, MRP, Euler and
rotation-vector representations, conversions, 2D broadcasting, inversion,
composition, angular magnitudes, axis-wise weighted means, and seeded random
streams. The Python reference independently reproduces the documented
Park–Miller/Shoemake sequence and statistically tests isotropic rotated vectors
against uniform-sphere expectations and SciPy's `Rotation.random` output.

This implementation covers **rank-two** rotation arrays, not the full
arbitrary-rank NumPy/SciPy API. It does not yet provide arbitrary slicing of
2D rotation grids, elementwise grid powers or all NumPy broadcasting cases.
`rotation_grid_from_matrix` shares the scalar constructor's restriction to
approximately proper 3x3 rotation matrices; it does not project an arbitrary
noisy matrix onto SO(3). New code is original MIT-licensed Fortran without
third-party implementation imports.

## Davenport-angle conversions

`rotation_grid_from_davenport(axes,order,angles,degrees=,status=)` constructs
rotations from generalized Euler angles. `axes(k,3)` contains one Cartesian
axis per row, with consecutive axes perpendicular. The angle array has
shape `(n,k)` for stacks or `(n,m,k)` for grids.

`grid%as_davenport(axes,order,...)` returns `(n,m,3)` angles and an
optional `(n,m)` gimbal-lock flag array.

See [Davenport angle documentation](spatial_transform.md#davenport-angles-generalized-euler-rotations).
