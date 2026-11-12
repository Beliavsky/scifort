# CubicSpline interpolation in SciFort

`scifort_interpolate` now supports `scipy.interpolate.CubicSpline`-style
interpolation for real arrays, **without additional runtime dependencies**.

```fortran
use scifort_kinds, only: dp
use scifort_interpolate, only: cubic_spline, piecewise_cubic_interpolator

type(piecewise_cubic_interpolator) :: spline
real(dp) :: x(5), y(5)
integer :: status

x = [0.0_dp, 0.5_dp, 1.0_dp, 2.0_dp, 3.0_dp]
y = sin(x)
spline = cubic_spline(x, y, bc_type='natural', status=status)
if (status /= 0) error stop 'invalid cubic spline samples'
print *, spline%evaluate(1.3_dp)
print *, spline%derivative(1.3_dp)
print *, spline%integrate(0.5_dp, 2.0_dp)
```

## Boundary conditions

`cubic_spline(x, y, bc_type, status, extrapolate, bc_left, bc_right,
left_value, right_value)` supports:

- `'not-a-knot'` (default), making adjacent outer cubic pieces identical in
  their third derivative. Two-knot special case is linear; the three-knot
  double-not-a-knot case is the unique quadratic.
- `'natural'`: both endpoints have zero second derivative.
- `'clamped'`: both endpoints have zero first derivative by default; to supply
  nonzero values, use `bc_left='first'`, `bc_right='first'`, and the appropriate
  `left_value`/`right_value`.
- `'periodic'`: matching endpoint ordinates are required. The default
  extrapolation and definite integration extend periodically, including
  across more than one full period. As in SciPy, explicitly setting
  `extrapolate=.true.` instead extends outermost cubics, and
  `extrapolate=.false.` yields NaN outside the original domain.

`bc_left` and `bc_right` independently override the two endpoints with
`'not-a-knot'`, `'natural'`, `'clamped'`, `'first'`, or `'second'`. Optional
`left_value` and `right_value` give the prescribed first or second derivative
values (default 0), but are not accepted for periodic splines. A periodic
boundary must be set at both ends.

## Multi-output interpolation

The `cubic_spline` generic also accepts `y(n, m)`, with `n=size(x)` and
`m` independent output channels. It returns
`type(piecewise_cubic_array_interpolator)` with the same method names:

```fortran
use scifort_kinds, only: dp
use scifort_interpolate, only: cubic_spline, piecewise_cubic_array_interpolator

type(piecewise_cubic_array_interpolator) :: s
real(dp) :: x(4), y(4,2)
real(dp), allocatable :: values(:,:), derivatives(:,:), integrals(:)

x = [0.0_dp, 0.5_dp, 1.0_dp, 2.0_dp]
y(:,1) = sin(x)
y(:,2) = x*x
s = cubic_spline(x, y, bc_type='natural')
values = s%evaluate([0.2_dp, 0.7_dp])  ! shape (2, 2)
derivatives = s%derivative([0.2_dp,0.7_dp],order=1)
integrals = s%integrate(0.0_dp, 2.0_dp) ! shape (2)
```

`%coefficients()` returns coefficients **in ascending order**, with
shape `(4, n-1)` for scalar ordinates or `(4, n-1, m)` for multichannel
ordinates; SciPy's `PPoly.c` instead stores descending powers. The matrix
constructor takes `left_values(m)` and `right_values(m)` instead of the
scalar constructor's `left_value`/`right_value`.

Both objects provide scalar and rank-one query arrays, derivatives of any
nonnegative order (zero after the cubic order), exact piecewise-polynomial
integration, explicit status codes and `%valid()` checks. Only rank-one knot
coordinates, real-valued data and output channels along the second axis are
supported in this increment. Arbitrary higher ordinate ranks and general
N-dimensional coordinate grids remain future work.

## Verification

```sh
python tools/interpolate_reference/compare_cubic_spline_scipy.py
python tools/interpolate_reference/compare_cubic_spline_multi_scipy.py
```

Both scripts compile separate GNU Fortran checked and optimized executables
and compare coefficients, value/derivative orders 0–4, extrapolation and exact
integrals to SciPy 1.17.0. Default seeds produce 847 scalar and 480 matrix
scenarios per compiler configuration, covering 197,351 and 214,704 values.
Second seeds can be selected with `--seed NUMBER`. Also run
`test/test_interpolate_cubic_spline.f90` for native API and rejection tests.

Implementation is original MIT-licensed code based on the standard cubic
spline continuity and end-boundary equations. No upstream numerical library
source was copied. See [upstream reuse review](interpolation_upstream_reuse.md)
for licensed future additions to SciFort interpolation.
