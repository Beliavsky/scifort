# Smooth cubic rotation splines (`scipy.spatial.transform.RotationSpline`)

The Fortran module `scifort_spatial_rotation_spline` provides SO(3) orientation
interpolation through keyframes with **continuous angular velocity and
continuous angular acceleration**. This is distinct from the existing
piecewise constant-angular-velocity `rotation_slerp_path`.

No new runtime libraries are needed: the angular-rate continuity equations are
solved using a small, pivoted, **3x3 block-tridiagonal** solver, whose storage
and arithmetic scale linearly with the number of keyframes. The algorithm
follows the mathematical formulation used in SciPy 1.17's `RotationSpline`.
The Fortran implementation and its test driver are written for SciFort.

## Public API

| Fortran symbol | Description |
|---|---|
| `rotation_spline(times, rotations, status=)` | Construct spline from at least two keyframes |
| `spline%valid()` | Check whether construction succeeded |
| `spline%at(time, status=)` | Return orientation (`rotation_type`) at time |
| `spline%angular_rate(time, status=)` | Return 3-vector angular velocity (rad per time unit) |
| `spline%angular_acceleration(time, status=)` | Return 3-vector angular acceleration (rad per time unit squared) |
| `spline%evaluate_path(times, rotations, angular_rates=, angular_accelerations=, status=)` | Batched orientations, with optional `(n,3)` rate and acceleration arrays |

All timestamps and rotations are finite; timestamps must be **strictly
increasing**, and the number of rotations must equal the number of timestamps.
At least two keyframes are needed. Output orientations and numerical rates
are valid at *unordered* evaluation times, including at keyframes, before
and after the time interval. Extrapolation is permitted, consistent with
SciPy's `RotationSpline.__call__`.

For invalid input the factory returns an object with `valid()==.false.` and
sets `status=spatial_rotation_invalid_input` (1). For extreme near-antipodal
paths where the nonlinear solve yields catastrophically large coefficients
that no longer reproduce the original orientations, it returns
`status=spatial_rotation_no_convergence` (2), defined in
`scifort_spatial_rotation_spline`, instead of an inaccurate interpolant.
Invalid scalar queries return an invalid `rotation_type`, or NaNs for real
derivatives; an invalid batch returns empty arrays. Success is
`spatial_rotation_success=0` (defined alongside the invalid-input constant
by `scifort_spatial_transform`).

```fortran
use scifort_kinds, only : dp
use scifort_spatial_transform, only : rotation_type, rotation_from_euler
use scifort_spatial_rotation_spline, only : rotation_spline_type, rotation_spline

! Declare and initialize:
real(dp) :: times(4), omega(3), alpha(3)
type(rotation_type) :: keys(4), attitude
type(rotation_spline_type) :: spline
integer :: status

! ... initialize times with increasing values and keys with valid rotations ...
spline = rotation_spline(times, keys, status=status)
if (status /= 0) error stop 'invalid attitude spline'
attitude = spline%at(2.5_dp)
omega = spline%angular_rate(2.5_dp)
alpha = spline%angular_acceleration(2.5_dp)
```

## Mathematical formulation

For each time segment, define the relative rotation vector
`v_i = log(R_i^{-1} R_{i+1})`. The interpolation is
`R(t) = R_i exp([v(t)]_x)` with a cubic polynomial in elapsed segment time.
The coefficients enforce the keyframe orientations exactly. The angular
velocity and acceleration include the noncommutative SO(3) Jacobian terms,
not just the ordinary first and second derivatives of the rotation-vector
polynomial. The unknown angular rates at the internal knots are solved by a
nonlinear iteration. Each linearized step has a block-tridiagonal system of
3x3 blocks. The **two-keyframe** special case reduces to a constant-rate
geodesic, as in Slerp.

Near the identity, the SO(3) Jacobians and nonlinear acceleration correction
use Taylor expansions to avoid cancellation. Endpoints use the same
one-sided derivative and extrapolation convention as SciPy. Rotation
quaternion signs are not significant: `q` and `-q` represent the same
orientation, so reference tests compare relative rotation angles.

## Verification

```text
python tools/spatial_reference/run_native_spatial.py --mode both
python tools/spatial_reference/compare_rotation_spline_scipy.py --mode both --seed 41 --cases 240
python tools/spatial_reference/compare_rotation_spline_scipy.py --mode both --seed 1729 --cases 240
```

The randomized tests include up to 160 keyframes. In a two-seed run with
480 configurations, two divergent near-antipodal paths were correctly
rejected with status 2 and 478 valid paths matched SciPy. These unstable
paths are *not* counted as successful numerical comparisons. In the first
rejected case SciPy itself fails to reproduce some specified keyframes,
consistent with severe coefficient cancellation; both examples are retained
as regression checks.

The comparison driver is a standalone GNU Fortran program that reads keyframe
quaternions and sample times, constructs SciFort splines, and writes results.
The Python comparison independently constructs SciPy `RotationSpline` objects
and checks rotations, angular rates and accelerations at keyframes, interior
points, and outside the keyframe interval. The native regression checks
keyframe interpolation, continuity, invalid input and batched queries.

## Limitations

* Only 3D SO(3) rotations (`rotation_type`) are supported, with one spline
  instance at a time; no higher-dimensional rigid transformations or
  broadcasting of batches of different splines.
* The constructor follows SciPy's ten-iteration angular-rate scheme. Some
  iterates may not satisfy the strict nonlinear stopping criterion, as in the
  SciPy reference. Catastrophic coefficient cancellation is detected by
  checking reconstructed keyframe rotations, and returns status 2 rather
  than an invalid trajectory; convergence is otherwise not guaranteed.
* The Jacobian evaluation assumes relative rotation vectors have principal
  magnitudes no larger than pi; exactly pi and nonunique paths can have
  convention-dependent derivative results.
* Extremely large extrapolation intervals can overflow the cubic polynomial.
  Such numerical results should not be interpreted as physically meaningful
  attitude trajectories.

See [SciPy RotationSpline](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.transform.RotationSpline.html)
and the [Smooth Attitude Interpolation](https://github.com/scipy/scipy/files/2932755/attitude_interpolation.pdf)
reference cited by SciPy. See `CODE_PROVENANCE.md` and
`THIRD_PARTY_LICENSES.md` for attribution and licensing context.
