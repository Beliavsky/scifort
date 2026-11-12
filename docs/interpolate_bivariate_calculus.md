# Bivariate spline differentiation and exact knot refinement (2026-11-01)

SciFort's `fitpack_bivariate_tck_type` now supports **constructing a new
bivariate B-spline for a mixed partial derivative**, corresponding to
SciPy's `BivariateSpline.partial_derivative(dx,dy)`, as well as exact
shape-preserving knot refinement on either axis. Both operations preserve the
FITPACK TCK coefficient convention and work with surfaces from
`rect_bivariate_spline` or the optional `perazz/fitpack` fitter.

The implementation uses SciFort's native B-spline machinery and has **no
external FITPACK, BLAS, LAPACK or other runtime dependency**. No Delaunay code
is changed.

## Public API

```fortran
use scifort_kinds, only : dp
use scifort_interpolate, only : fitpack_bivariate_tck_type, &
    rect_bivariate_spline, bisplev, bispline_success

real(dp) :: x(6), y(7), z(6,7), q, qx, qy
integer :: i, j, status
type(fitpack_bivariate_tck_type) :: surface, diff, refined

x = [0.0_dp,0.1_dp,0.25_dp,0.4_dp,0.8_dp,1.0_dp]
y = [0.0_dp,0.15_dp,0.3_dp,0.5_dp,0.65_dp,0.8_dp,1.0_dp]
do j = 1, size(y)
    do i = 1, size(x)
        z(i,j) = 1.0_dp + x(i)**2*y(j) + sin(x(i)+y(j))
    end do
end do
surface = rect_bivariate_spline(x,y,z,status=status)
if (status /= bispline_success) error stop 'surface fit failed'

! Mixed partial derivative is a reusable B-spline, not a sampled grid.
diff = surface%partial_derivative(1,1,status=status)
if (status /= bispline_success) error stop 'derivative failed'
qx = 0.35_dp
qy = 0.52_dp
q = bisplev(qx,qy,diff)

! Insert the same interior x knot twice, then one y knot.
! The represented surface remains mathematically identical.
refined = surface%insert_knot(axis=1,x=0.35_dp,m=2,status=status)
if (status /= bispline_success) error stop 'x knot insertion failed'
refined = refined%insert_knot(axis=2,x=0.52_dp,status=status)
if (status /= bispline_success) error stop 'y knot insertion failed'
if (abs(bisplev(qx,qy,refined)-bisplev(qx,qy,surface)) > 1.0e-10_dp) &
    error stop 'refinement changed values'

! Integrals are available for derivative and refined splines alike.
q = diff%integral(0.0_dp,1.0_dp,0.0_dp,1.0_dp)
```

Both operations also have functional aliases:

```fortran
! use scifort_interpolate, only : bispline_partial_derivative, bispline_insert_knot
diff = bispline_partial_derivative(surface,dx=1,dy=1,status=status)
refined = bispline_insert_knot(surface,axis=1,x=0.35_dp,m=2,status=status)
```

### `%partial_derivative(dx,dy,status)`

- `dx`, `dy`: explicit nonnegative derivative orders satisfying `dx < kx`
  and `dy < ky`, as required by SciPy's `BivariateSpline` derivative object.
  Zero/zero returns a valid copy with only active coefficients.
- The output degrees are `kx-dx`, `ky-dy`, with one outer knot removed per
  differentiation order along the corresponding axis. Coefficients are
  transformed by exact B-spline differentiation identities (no fitting).
- `status` is optional. On an invalid TCK, invalid order, or division by a
  zero-width knot interval, it returns `bispline_invalid_input` (1) and an
  invalid result, rather than dividing by zero.
- Subsequent `bisplev`, `bisplev_points`, and `%integral` calls work on the
  transformed TCK. For the same query coordinates, an original TCK's
  `bisplev(dx=...,dy=...)` agrees with evaluation of the derivative TCK.

The coefficient recurrence for x differentiation is

\[
C^{(x)}_{i,j} = \frac{k_x}{t^x_{i+k_x+1}-t^x_{i+1}}
    (C_{i+1,j}-C_{i,j}),
\]

and analogously in y. The original FITPACK coefficients are indexed as
`C(i,j) = c((i-1)*ny+j)` (y index fastest), so the code transforms the two
axes using different strides and updates the knots/degrees at each pass.

### `%insert_knot(axis,x,m,status)`

- `axis=1` selects the x knot sequence, `axis=2` the y sequence.
- `x` must be finite and within the selected axis's active spline domain.
- `m` is the **nonnegative** number of insertions (default 1). Zero returns
  the same surface in canonical, unpadded coefficient form.
- Existing interior knots may be repeated up to total multiplicity `degree+1`.
  Inserting again after this limit is rejected. Both axes preserve their
  original polynomial degrees; only the chosen knot sequence and the number
  of coefficient fibers change.
- A fully repeated interior knot (`degree+1` copies) is valid for evaluating
  the spline but can make its *derivative representation* singular. If a
  derivative object will be needed afterwards, avoid creating a full break.
- Errors return `bispline_invalid_input` (1) and an invalid result. On success,
  all values and supported derivatives and integrals remain unchanged to
  floating-point precision. Knot insertion is **not smoothing or refitting**.

Boehm's standard one-dimensional knot-insertion recurrence is applied to
all fibers of the selected tensor axis, using the existing tested SciFort
`BSpline%insert_knot` algorithm. The output coefficients are stored in the
same y-fastest flattened FITPACK order. Extra padded coefficients sometimes
returned by FITPACK are ignored by either transformation.

## Validation and limitations

`tools/interpolate_reference/compare_bivariate_calculus_scipy.py` uses SciPy
1.17.0's bivariate FITPACK spline construction and its 1D `BSpline` operators
as independent references for the resulting knots, coefficients, evaluations,
mixed derivatives and integrals. The direct SciPy
`BivariateSpline.partial_derivative` is also compared wherever its wrapper
accepts the requested degree pair. (The SciPy 1.17.0 FITPACK wrapper rejects
some otherwise valid asymmetric-degree derivative orders; the 1D `BSpline`
reference covers those cases.)

At the default seed, **180/180 randomized scenarios pass in each** GNU
Fortran 14.2.0 bounds-checked and optimized `-O3` configurations. The cases
include 1–5 x/y degrees (including asymmetric), irregular and nonsquare
rectangular grids, repeated knots, zero to several insertions, derivative
orders of different sizes, out-of-range evaluation coordinates, and signed
rectangle integrals. Maximum normalized knot/coefficients/evaluation errors
are approximately `0` / `9.6e-12` / `6.8e-12`; integrals agree to roughly
`5e-14` in the default-seed test set. A second seed of 200 independent cases
per build also passed, with maximum normalized derivative-integral error
below `7e-13`.

The native regression includes analytically differentiable bicubic surfaces,
exact integration of the transformed surface, invariance under two-axis knot
insertion, optional-status validation, pre-existing knots, fully repeated
knot handling and padded FITPACK coefficient arrays. **29/29 native public
interpolation tests pass** with bounds checking.

To reproduce:

```sh
python tools/interpolate_reference/compare_bivariate_calculus_scipy.py \
  --cases 180 --mode both
python tools/interpolate_reference/compare_bivariate_calculus_scipy.py \
  --cases 200 --seed 20261102 --mode both
python tools/interpolate_reference/run_all_native_interpolation.py --only checked
python tools/interpolate_reference/run_all_native_interpolation.py --only optimized
```

These are real, finite, nonperiodic tensor-product B-splines. They are **not**
NURBS, rational splines, adaptive fitters, or a replacement for FITPACK's
weighted smoothing. Insertion is repeated over coefficient fibers, which may
not be the fastest implementation for millions of coefficients. For large
regular grids, banded tensor-product interpolation and optimized batch
contraction remain separate performance opportunities.
