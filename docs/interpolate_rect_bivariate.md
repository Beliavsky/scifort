# Rectangular-grid bivariate splines and exact surface integration (2026-11-01)

SciFort now provides an independent, dependency-free implementation of the
**exact-interpolation** case of SciPy's `RectBivariateSpline(x, y, z, s=0)`.
It builds a tensor-product B-spline for data sampled on a rectangular grid,
using the existing native 1D not-a-knot collocation solver. Its result is the
same `fitpack_bivariate_tck_type` already accepted by `bisplev`,
`bisplev_points`, and the optional `perazz/fitpack` backend.

## Public API

```fortran
use scifort_kinds, only : dp
use scifort_interpolate, only : fitpack_bivariate_tck_type, &
    rect_bivariate_spline, bisplev, bispline_integral, bispline_success

real(dp) :: x(5), y(6), z(5,6), answer, area
integer :: i, j, status
type(fitpack_bivariate_tck_type) :: surface

x = [0.0_dp,0.2_dp,0.5_dp,0.8_dp,1.0_dp]
y = [0.0_dp,0.1_dp,0.3_dp,0.6_dp,0.8_dp,1.0_dp]
do j = 1,size(y)
    do i = 1,size(x)
        z(i,j) = sin(x(i)) + x(i)*y(j) + y(j)**2
    end do
end do
surface = rect_bivariate_spline(x,y,z,kx=3,ky=3,status=status)
if (status /= bispline_success) error stop 'Cannot interpolate grid'
answer = bisplev(0.4_dp,0.7_dp,surface)
area = surface%integral(0.0_dp,1.0_dp,0.0_dp,1.0_dp)
! Identical to bispline_integral(0.0_dp,1.0_dp,0.0_dp,1.0_dp,surface).
```

- `x(:)`, `y(:)` must be strictly increasing and finite, with at least
  `kx+1` and `ky+1` sites respectively. Irregular spacing is supported.
- `z(size(x),size(y))` is a complete **real** rectangular grid, x varying
  along the first Fortran array dimension. All observations must be finite.
- `kx`, `ky` are independent B-spline degrees from 1 to 5, default 3.
- The result has open (clamped) not-a-knot knot vectors `tx`, `ty` and
  coefficients `c((ix-1)*ny+iy)`, where the **y coefficient varies fastest**.
- `status` is optional. `bispline_success=0`, `bispline_invalid_input=1`,
  and `bispline_singular=2` (collocation solver failure). The returned TCK is
  invalid on an error and `surface%valid()` returns false.

The returned TCK also supports constructing derivative-spline objects and
exact shape-preserving x/y knot insertion; see
[bivariate spline calculus and refinement](interpolate_bivariate_calculus.md).

For evaluation, `bisplev(x,y,surface)` takes scalar coordinates or two
coordinate vectors for a tensor-product grid; `bisplev_points` evaluates
paired x and y query vectors. Evaluations outside the support are clamped to
the nearest boundary, matching FITPACK `bispev` conventions. The supported
partial derivative orders must be **strictly less than** `kx` and `ky`.

The constructor supports exact interpolation (`s=0`) only. It does **not**
implement `RectBivariateSpline`'s `s>0` FITPACK smoothing, custom bounding
boxes, missing grid values, or per-observation weights. For scattered
observation smoothing or fixed-knot least-squares fitting, use SciFort's
separately documented optional `perazz/fitpack` backend. The implementation
uses 1D dense collocation solves repeatedly for each coordinate direction;
large grids may benefit from shared banded factorization in a future change.

## Exact rectangular integration

Both the TCK method

```fortran
v = surface%integral(xa,xb,ya,yb,status=status)
```

and the free function

```fortran
v = bispline_integral(xa,xb,ya,yb,surface,status=status)
```

perform **analytical tensor-product integration**, using the B-spline
antiderivative along y and then along x. They are available for any valid
bivariate FITPACK TCK, including surfaces fitted by the optional backend.
Integration clips each axis to the spline's active domain (there is no
out-of-bounds extrapolation), returns zero for an empty intersection, and
preserves the orientation: reversing one pair of bounds negates the result;
reversing both preserves it. A nonfinite limit or an invalid TCK produces a
quiet NaN with status 1. No numerical quadrature or external dependency is
required.

## Method and validation

The mathematical construction is tensor-product B-spline interpolation:
first collocate independent y-columns on the x axis, then collocate the
transposed coefficient matrix on the y axis. This relies on existing
SciFort `make_interp_spline` without importing any SciPy source.
The exact rectangular integral follows from the tensor product of the
standard 1D B-spline integral identities.

Independent SciPy 1.17.0 numerical comparisons exercise degrees 1 through 5,
asymmetric knot sequences, rectangular (nonsquare) grids, noisy surface
values, out-of-domain queries, mixed partial derivatives, and signed/clipped
integrals. **120/120** scenarios passed in each of GNU Fortran 14.2.0's
bounds-checked and `-O3` configurations. Maximum relative-to-unit-scale
errors in the test sample:

| Quantity | Maximum observed error |
| --- | ---: |
| x/y knot vectors | 0 |
| coefficient array | 9.52e-14 |
| scalar values and derivatives | 3.18e-13 |
| exact integrals | 3.20e-15 |

A separate native regression verifies grid-node interpolation, analytic
bicubic polynomial values and derivatives, an analytic rectangular integral,
all 25 pairs of degrees on their minimum-sized grids, tensor and paired
queries, zero-width/outside rectangles, input errors, and NaN behavior.
All 28 native interpolation regression programs passed in the checked build.

To reproduce, with Python, NumPy, SciPy and GNU Fortran available:

```sh
python tools/interpolate_reference/compare_rect_bivariate_scipy.py --cases 120 --mode both
python tools/interpolate_reference/run_all_native_interpolation.py --only checked
python tools/interpolate_reference/run_all_native_interpolation.py --only optimized
```

The reference driver is `tools/interpolate_reference/rect_bivariate_driver.f90`.
The core is independent MIT-licensed SciFort code. SciPy is used as an
independent reference only; there is no FITPACK vendoring or new FPM
runtime dependency.
