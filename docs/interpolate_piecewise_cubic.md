# 1-D shape-preserving and local cubic interpolation

The native `scifort_interpolate` module introduces SciPy-compatible **numerical
algorithms** for three kinds of cubic Hermite interpolation:

| SciPy concept | Fortran factory | Behavior |
| --- | --- | --- |
| `scipy.interpolate.PchipInterpolator` | `pchip_interpolator(x, y, status, extrapolate)` | Nonovershooting monotone-piecewise data, default extrapolation **true** |
| `scipy.interpolate.Akima1DInterpolator` | `akima1d_interpolator(x, y, method, status, extrapolate)` | Local weighted slopes, default method `'akima'`, default extrapolation **false** |
| `Akima1DInterpolator(method='makima')` | `akima1d_interpolator(x, y, 'makima', ...)` | Modified Akima weights, suppresses many overshoots |
| `scipy.interpolate.CubicHermiteSpline` | `cubic_hermite_spline(x, y, dydx, status, extrapolate)` | Prescribed first derivatives at knots, default extrapolation **true** |
| `scipy.interpolate.pchip_interpolate` | `pchip_interpolate(x, y, z, order, status, extrapolate)` | Convenience scalar or vector evaluation/derivative |

Each factory returns a `type(piecewise_cubic_interpolator)` object, exposing
`%evaluate(z)`, `%derivative(z, order)`, `%integrate(a,b)`, `%valid()`,
`%knots()` and `%coefficients()`. Scalar or rank-one real query arrays are
accepted. `evaluate()` defaults to value order zero, and `derivative()` to
first order; orders four or greater return zero on the supported domain.
`extrapolate` can be overridden on an individual evaluation, derivative or
integral. The default constructor value is retained otherwise.

```fortran
program example_cubic
    use scifort_kinds, only : dp
    use scifort_interpolate, only : piecewise_cubic_interpolator, &
        pchip_interpolator, akima1d_interpolator, cubic_hermite_spline
    implicit none
    type(piecewise_cubic_interpolator) :: p, a, h
    real(dp) :: x(4), y(4)
    integer :: status

    x = [0.0_dp, 1.0_dp, 2.0_dp, 4.0_dp]
    y = [0.0_dp, 1.0_dp, 1.2_dp, 3.0_dp]
    p = pchip_interpolator(x, y, status)
    if (status /= 0) error stop 'invalid samples'
    print *, 'pchip(1.5)=', p%evaluate(1.5_dp)
    print *, 'pchip derivative=', p%derivative(1.5_dp)
    print *, 'integral [0,4]=', p%integrate(0.0_dp, 4.0_dp)

    a = akima1d_interpolator(x, y, method='makima')
    print *, 'makima at knots=', a%evaluate(x)

    h = cubic_hermite_spline(x, y, [1.0_dp,0.5_dp,0.5_dp,1.0_dp])
    print *, 'prescribed-slope Hermite=', h%evaluate(1.5_dp)
end program
```

## Mathematical definition

For each interval `[x(i),x(i+1)]`, let `h=x(i+1)-x(i)`,
`s=(y(i+1)-y(i))/h`, and `d0,d1` be the supplied or calculated knot
derivatives. The function is stored as an ordinary polynomial in
`t=z-x(i)`:

```text
p(t) = y(i) + d0*t + ((3*s-2*d0-d1)/h)*t^2
       + ((d0+d1-2*s)/h^2)*t^3
```

PCHIP uses the Fritsch–Butland weighted harmonic mean at inner knots, a
one-sided endpoint estimate, and zero slope across extrema and plateaus.
Akima/Makima use four neighboring interval slopes, extended linearly at the
ends. When the Akima weight sum is effectively zero, the same outer-slope
tie convention as the SciPy reference is used. Integration is **exact** for
each cubic piece; cumulative integrals are cached at knots. The binary
interval search handles arbitrary nonuniform knot spacing.

## Validation and limitations

- `x` must contain at least two **strictly increasing finite** numbers;
  `y` and, for Hermite, `dydx` must match the knot count and be finite.
- Invalid input returns an uninitialized object (`%valid()==.false.`) with
  status `interpolate_status_invalid_input=1`; evaluation then returns a
  quiet NaN. Invalid or nonfinite queries also return NaN.
- Unlike SciPy's `PPoly.c`, `%coefficients()` has rows in **ascending**
  polynomial power (constant, linear, quadratic, cubic).
- Current coverage is **real one-dimensional `y`**, with real scalar/rank-one
  queries. Multichannel `y`, arbitrary interpolation axis, periodic
  extrapolation, complex ordinates, root-finding and inverse solutions are not
  implemented. No `CubicSpline` boundary-condition or smoothing APIs are
  claimed. If finite input values are so large that intermediate polynomial
  coefficients overflow binary64, fitting reports an invalid status when
  those values are detected; use appropriately scaled data.
- This is not an interface to SciPy/Python at runtime. The production
  implementation requires only `scifort_kinds` and intrinsic Fortran modules.

Tests: `test/test_interpolate_cubic.f90` (native assertions) and
`tools/interpolate_reference/compare_pchip_akima_scipy.py` (480 randomized
SciPy comparisons per mode, including derivatives through order four, exact
integrals, extrapolation, and polynomial coefficients). Run from any directory:

```sh
python tools/interpolate_reference/compare_pchip_akima_scipy.py
```

The Python test compiles fresh checked and optimized GNU Fortran executables
in temporary directories. It requires NumPy and SciPy for independent
reference results; the Fortran library does not.

## Mathematical and API references

- F. N. Fritsch and J. Butland, *A Method for Constructing Local Monotone
  Piecewise Cubic Interpolants*, SIAM J. Sci. Stat. Comput. **5** (1984),
  300–304. DOI: `10.1137/0905021`.
- H. Akima, *A New Method of Interpolation and Smooth Curve Fitting Based on
  Local Procedures*, J. ACM **17** (1970), 589–602.
  DOI: `10.1145/321607.321609`.
- C. Moler and C. Ionita, *Makima Piecewise Cubic Interpolation* (2019),
  https://blogs.mathworks.com/cleve/2019/04/29/makima-piecewise-cubic-interpolation/
- SciPy's publicly documented algorithm conventions (reference only):
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.PchipInterpolator.html
  and
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.Akima1DInterpolator.html

All implementation code is original Fortran based on these published
mathematical formulas; no third-party implementation was copied.
