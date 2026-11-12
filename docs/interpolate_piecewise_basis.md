# General piecewise polynomials: PPoly and BPoly

`scifort_interpolate` now exports original Fortran implementations of two
representations in `scipy.interpolate`:

- `ppoly_type` and `ppoly`: the **local power basis**, with coefficient
  `c(k+1-j, interval, channel)` multiplying `(z-x(interval))**j`.
- `bpoly_type` and `bpoly`: the **Bernstein basis**, with coefficients
  `c(j+1,interval,channel)` multiplying
  `binomial(k,j)*t**j*(1-t)**(k-j)`, where
  `t=(z-x(interval))/(x(interval+1)-x(interval))`.

Both are independently implemented in `src/interpolate_piecewise.f90`.
The same `ppoly()` and `bpoly()` constructors also accept complex coefficients
and return `complex_ppoly_type` and `complex_bpoly_type`, respectively.
Their conventions were checked against SciPy 1.17.0; no SciPy
implementation code or third-party interpolation source was copied.

## Example

```fortran
use scifort_kinds, only : dp
use scifort_interpolate, only : ppoly_type, bpoly_type, ppoly, &
    bpoly_from_power, ppoly_from_bernstein

type(ppoly_type) :: pp, recovered
type(bpoly_type) :: bp

! f(x) = 3*x + 2 on both [0,1] and [1,2].
pp = ppoly(reshape([3.0_dp, 2.0_dp, 3.0_dp, 5.0_dp], [2,2]), &
    [0.0_dp, 1.0_dp, 2.0_dp])
bp = bpoly_from_power(pp)
recovered = ppoly_from_bernstein(bp)

print *, pp%value(1.25_dp)     ! 5.75
print *, bp%value(1.25_dp)     ! 5.75
print *, pp%integral(0.0_dp, 2.0_dp)  ! 10.0
```

Construction accepts rank-2 `(degree+1, intervals)` coefficients for a
single channel or rank-3 `(degree+1, intervals, channels)` coefficients.
All public coefficient arrays use **the same dimension order** as the
standard SciPy `axis=0` convention. `knots()` and `coefficients()` return
allocatable copies. `valid()` and optional constructor `status=` permit
nonfatal input validation.

`%evaluate(z)` returns a vector with one entry per output channel for
scalar `z`, and a matrix of shape `(size(z),channels)` for vector `z`.
For single-output cases, `%value(z)` returns a scalar and `%values(z)`
returns a vector. The optional `nu=` argument requests the derivative order
at evaluation time.

Additional type-bound methods are `%derivative(nu)`,
`%antiderivative(nu)`, `%integrate(a,b)` (a vector of channel integrals),
`%integral(a,b)` (a scalar single-channel convenience method),
`%channels()`, `%knots()`, and `%coefficients()`.
Negative orders to `%derivative()` and `%antiderivative()` perform the
opposite operation. Repeated antiderivatives are continuous at every knot.

Breakpoints may be nondecreasing or nonincreasing. Exactly repeated knots
are supported in the power evaluator; a zero-width Bernstein interval has
undefined nonzero derivatives and is treated conservatively. At internal
knots the evaluator uses the interval to the right in the supplied knot
order, matching SciPy. `extrapolate=.false.` returns quiet NaNs outside
its finite range. `periodic=.true.` extends the interval periodically;
periodic integration handles multiple periods and reversed limits.

## Extending PPoly and BPoly piecewise polynomials

`call poly%extend(new_coefficients, new_breakpoints, status)` appends or prepends
segments to a valid existing power-basis or Bernstein-basis polynomial.
`new_coefficients` may have either rank 2 `(degree+1, new_intervals)` for a
single channel or rank 3 `(degree+1, new_intervals, channels)`. The number of
new breakpoints must equal the number of new segments: the existing first or
last breakpoint is **not** repeated in the input. Both the existing and new
breakpoints can increase or decrease in the same direction. Mixed directions,
overlapping domains, inconsistent channels, and nonfinite coefficients are
rejected (`status=piecewise_status_invalid_input`) without changing the
original object.

For `PPoly`, old or new lower-degree sections are promoted by *leading*
zero-padding, maintaining descending-power ordering. For `BPoly`, degree
elevation is required to preserve polynomial values; SciFort uses successive
Bernstein degree-elevation recurrences, including complex coefficients.

```fortran
use scifort_kinds, only : dp
use scifort_interpolate, only : ppoly, ppoly_type

type(ppoly_type) :: pp
integer :: status

pp = ppoly(reshape([1.0_dp,0.0_dp],[2,1]),[0.0_dp,1.0_dp])
! Add (x-1)**2 + 1 on [1,2] (second-degree segment).
call pp%extend(reshape([1.0_dp,0.0_dp,1.0_dp],[3,1]), &
    [2.0_dp],status=status)
print *, status, pp%value(1.5_dp)  ! 0, 1.25
```

SciPy 1.17.0 has a bug in some *prepending* calls to `PPoly.extend` with
mismatched polynomial degrees. The randomized tests compare SciFort against
the mathematically correct explicit piecewise coefficient concatenation for
those cases, and directly against SciPy for the remaining cases.

## General Hermite interpolation (`BPoly.from_derivatives`)

The `bpoly_from_derivatives` constructor accepts `x`, `y`, and `counts`:

- `x(:)` contains strictly *increasing* knots.
- `y(i,j,k)` gives the derivative of order `j-1` at knot `i` for channel
  `k`. The rank-2 scalar-channel overload uses `y(i,j)`.
- `counts(i)` indicates how many derivatives, from order zero upward, are
  available at knot `i`; unused slots in `y` are ignored.

Optional `orders=[k_1,...,k_(n-1)]` requests a local **polynomial degree**
on each interval. As a Fortran convenience, `order=k` sets the same degree
on all intervals; these two optional arguments are mutually exclusive. The
local interpolant retains as many derivative constraints as required from
both ends, taking an extra right-end constraint when the required number is
odd. If one endpoint lacks constraints, the other supplies the balance.
Invalid degrees or impossible constraints return `status=1` and an invalid
polynomial without modifying the input arrays.

```fortran
! y(i,j) contains the derivative of order j-1 at knot i.
! Two intervals use degree 3 and degree 1, respectively.
bp = bpoly_from_derivatives([0.0_dp,1.0_dp,2.0_dp], &
    derivatives, [2,1,2], orders=[2,1])
```

SciPy 1.17.0 has an asymmetric-endpoint bookkeeping defect in its order
selection. The additional tests compare those cases to an independently
solved Hermite linear system and use the SciPy public API directly where it
agrees with that independent reference.

Each interval is the unique Hermite interpolant satisfying the specified
endpoint derivatives. Intervals of lower degree are converted to the global
maximum Bernstein degree by exact degree elevation, as in SciPy. This is a
broader construction than `CubicHermiteSpline` because it supports any
number of supplied derivatives.

Example:

```fortran
use scifort_kinds, only : dp
use scifort_interpolate, only : bpoly_type, bpoly_from_derivatives

type(bpoly_type) :: bp
real(dp) :: derivatives(3,2)

derivatives = 0.0_dp
derivatives(:,1) = [0.0_dp, 1.0_dp, 0.0_dp]
derivatives(1,2) = 1.0_dp  ! f'(0)
derivatives(3,2) = -1.0_dp ! f'(2)
bp = bpoly_from_derivatives([0.0_dp, 1.0_dp, 2.0_dp], &
    derivatives, [2,1,2])
print *, bp%value(0.5_dp)
```

`ppoly_from_cubic()` converts an existing SciFort PCHIP, Akima, Makima,
CubicHermiteSpline, or CubicSpline object (scalar or multichannel) into a
`PPoly` without refitting. Its coefficient order is reversed to the SciPy
power-basis convention while preserving the original extrapolation policy.

## Real polynomial roots and equation solving

`ppoly_type` provides `%roots(discontinuity,extrapolate,channel,status)` and
`%solve(y,discontinuity,extrapolate,channel,status)` for **real-valued**
polynomials. `y` defaults to zero, `discontinuity` to true, and `channel`
to the first channel (1-based). Use the channel selector to obtain roots of
each column of a multichannel polynomial. A zero-polynomial interval is
represented by its starting knot followed by a quiet NaN, matching SciPy.
Degenerate zero-width intervals contribute the point only, without NaN.
Optional `status` returns `piecewise_status_invalid_input` for an invalid
channel, uninitialized object, or nonfinite right-hand side. An invalid
request yields an empty array instead of raising an exception.

The root finder recursively isolates the real zeros of each derivative and
bisects sign-changing monotonic subintervals. It therefore detects even-
multiplicity tangencies, unlike sign-change-only solvers, and needs no
system BLAS/LAPACK. Roots are emitted in interval-processing order and may
differ in ordering from SciPy, especially when extrapolation finds several
roots outside the knot domain. Their numerical locations and sentinel counts
are independently verified against SciPy. Extremely ill-conditioned,
clustered or high-degree roots can still challenge double precision.

```fortran
use scifort_kinds, only : dp
use scifort_interpolate, only : ppoly, ppoly_type

type(ppoly_type) :: p
real(dp), allocatable :: zeros(:)

! (x-0.25)**2 on the first interval; a positive polynomial on the second.
p = ppoly(reshape([1.0_dp,-0.5_dp,0.0625_dp, &
    1.0_dp,0.0_dp,1.0_dp],[3,2]),[0.0_dp,1.0_dp,2.0_dp], &
    extrapolate=.false.)
zeros = p%roots()                ! [0.25]
```

## B-spline conversion (`PPoly.from_spline`)

`ppoly_from_spline(spline,status=...)` accepts either a
`bspline_interpolator_type` (scalar output) or
`bspline_array_interpolator_type` (multiple channels). It forms the local
power basis by evaluating each spline derivative at every original knot and
dividing by its factorial. As in SciPy, **the full knot vector is retained**,
including repeated end knots and zero-width sections. The default PPoly
extrapolation policy is true, like SciPy's `PPoly.from_spline`; optional
`extrapolate` and `periodic` flags override this.

```fortran
use scifort_kinds, only : dp
use scifort_interpolate, only : bspline, ppoly_from_spline, &
    bspline_interpolator_type, ppoly_type

type(bspline_interpolator_type) :: spline
type(ppoly_type) :: powers

spline = bspline([0.0_dp,0.0_dp,0.0_dp,1.0_dp,2.0_dp,2.0_dp,2.0_dp], &
    [1.0_dp,2.0_dp,3.0_dp,4.0_dp],2)
powers = ppoly_from_spline(spline)
```

## Complex coefficients

`complex_ppoly_type` and `complex_bpoly_type` are created by passing
`complex(dp)` coefficients to `ppoly()` and `bpoly()`. Both provide the same
`%evaluate`, `%value`, `%values`, `%derivative`, `%antiderivative`,
`%integrate`, and `%integral` methods as the real-valued types, including
periodic extension, descending knots and multichannel outputs. Output values
and coefficients are `complex(dp)`; breakpoints and query coordinates remain
real `real(dp)` arrays.

Use `complex_bpoly_from_power()` and `complex_ppoly_from_bernstein()` to
change basis; these preserve coefficient channels and extrapolation options.
Complex root solving is not provided, consistent with SciPy's
`PPoly.solve` real-polynomial restriction.

## Scope and known numerical limits

- One real interpolation axis, with an optional trailing channel dimension,
  is supported. Arbitrary-rank ordinate arrays and a general coefficient
  axis argument are not yet implemented.
- Both `%extend` and the variable-degree `orders` argument to
  `bpoly_from_derivatives` are implemented. The coefficient layout remains
  `(degree+1, intervals, channels)`; SciFort supports a single trailing
  channel dimension rather than arbitrary-rank output arrays.
- The implementation rejects nonfinite input coefficients and completely
  constant knot vectors. SciPy's validation has some differences.
- For Bernstein derivatives above the polynomial degree, SciFort returns the
  mathematically correct zero polynomial (with the usual bounds policy).
  SciPy 1.17 sometimes raises `ValueError` when the requested order is more
  than degree+1.
- SciPy 1.17 has inconsistent integration behavior for *descending* knots
  with periodic extrapolation. Those cases are compared independently against
  SciPy's nonperiodic antiderivative plus an explicit periodic decomposition.
- Converting high-degree polynomials to the power basis may be ill-conditioned;
  strong extrapolation can amplify ordinary double-precision roundoff. The
  checked numerical references use tight tolerances except when high-degree
  converted/antiderivative polynomials are evaluated outside the data span.
  Degree 20–30 or higher may be numerically unstable even in SciPy.

## Reproduce verification

```bash
python tools/interpolate_reference/compare_piecewise_scipy.py
python tools/interpolate_reference/compare_piecewise_scipy.py --seed 8675309
python tools/interpolate_reference/compare_bernstein_hermite_scipy.py
python tools/interpolate_reference/compare_bernstein_hermite_scipy.py --seed 8675309
python tools/interpolate_reference/compare_piecewise_roots_scipy.py
python tools/interpolate_reference/compare_piecewise_roots_scipy.py --seed 8675309
python tools/interpolate_reference/compare_complex_piecewise_scipy.py
python tools/interpolate_reference/compare_complex_piecewise_scipy.py --seed 8675309
python tools/interpolate_reference/compare_piecewise_extend_scipy.py
python tools/interpolate_reference/compare_piecewise_extend_scipy.py --seed 8675309
python tools/interpolate_reference/compare_bernstein_orders_scipy.py
python tools/interpolate_reference/compare_bernstein_orders_scipy.py --seed 8675309
```

All listed reference scripts compile a fresh GNU Fortran driver with checked
(`-O0 -fcheck=all`) and optimized (`-O3`) builds before comparing numerical
results with SciPy. Native `test/test_interpolate_piecewise.f90` and
`test/test_interpolate_piecewise_roots_complex.f90` and
`test/test_interpolate_piecewise_extend_orders.f90` exercise constructors,
calculus, conversion, descending and periodic breakpoints, invalid input, scalar/matrix Hermite interpolation and the
existing CubicSpline bridge. FPM auto-discovers the test file when installed.
