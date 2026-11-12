# `make_lsq_spline`: weighted B-spline fitting with complete knots

**Added 2026-11-01.** SciFort implements a real-valued subset of
`scipy.interpolate.make_lsq_spline` with independent, MIT-licensed Fortran.
It returns the existing `bspline_interpolator_type` or
`bspline_array_interpolator_type`. SciPy is needed only for optional numerical
validation, **not** at runtime. No BLAS or LAPACK dependency is introduced.

## Usage

```fortran
use scifort_kinds, only: dp
use scifort_interpolate, only: make_lsq_spline, bspline_interpolator_type
real(dp) :: x(8), y(8), w(8), interior(2)
real(dp), allocatable :: knots(:)
type(bspline_interpolator_type) :: fitted
integer :: status

x = [0.0_dp,0.3_dp,0.5_dp,0.8_dp,1.1_dp,1.5_dp,2.1_dp,3.0_dp]
y = sin(x)
w = 1.0_dp
interior = [0.7_dp,1.8_dp]
knots = [spread(x(1),1,4),interior,spread(x(8),1,4)]
fitted = make_lsq_spline(x,y,knots,k=3,w=w,status=status,method='qr')
if (status /= 0) error stop 'least-squares fit failed'
print *, fitted%evaluate(1.2_dp)
print *, fitted%coefficients()
```

## API and differences from FITPACK

`make_lsq_spline(x, y, t, k, w, status, extrapolate, method)` accepts a **full
knot vector** `t`, including the endpoint repetitions, rather than the
interior-only `t` parameter expected by `lsq_univariate_spline`.
`k` defaults to 3 and supports 0 through 15. Observations `x` must be finite
and nondecreasing (equal neighboring positions are allowed); `y` must be
finite with shape `(m)` or `(m, channels)`. The first axis always corresponds
to the observations (SciPy's `axis=0`). Knots must be finite, nondecreasing,
and span every observation. Interior knots may repeat. The spline must have
at least `k+1` data points and no more coefficients than observations.

Optional `w(m)` consists of finite, **nonnegative** residual weights. A zero
weight excludes the corresponding observation; the remaining nonzero-weight
rows must still give full column rank. As in SciPy, the criterion is

\[
    \min_{c} \sum_{i=1}^{m} \left[w_i
      \left(y_i-\sum_j c_j B_{j,k}(x_i)\right)\right]^2.
\]

The solution has shared knots but independent coefficients per channel.
`%evaluate`, `%derivative`, `%integrate` and `%coefficients` reuse the existing
B-spline object methods. `extrapolate` defaults to `.true.`;
`extrapolate=.false.` returns NaN for evaluations outside the base interval.

`method='qr'` is the only implemented method; `method='norm-eq'` is explicitly
rejected rather than silently using a different solver. Other SciPy options,
including non-leading interpolation axes, complex-valued ordinates and
`check_finite=.false.`, are not yet exposed. Unlike the FITPACK wrapper,
`make_lsq_spline` does not require strictly increasing `x` and accepts zero
weights and higher polynomial degrees.

All constructors return an invalid, non-evaluable object on failure:

- `bspline_status_success` (0): fit succeeded.
- `bspline_status_invalid_input` (1): invalid shapes, knot order/range,
  weights, degree, or solver method.
- `bspline_status_singular` (2): rank-deficient collocation matrix or
  numerically singular fit. A sparse sample region or excessive internal knot
  multiplicity can trigger this result.

## Method and complexity

SciFort uses online Givens rotations to factor the weighted local B-spline
collocation rows directly. Only `k+1` basis functions can be nonzero at a
sample site, so the upper-triangular factor is stored as a compact band.
The algorithm avoids forming normal equations and reduces memory to
`O((k+1) * ncoeff + ncoeff * channels)`; computational effort is
`O(m * k * (k + channels) + ncoeff * k * channels)` (with the natural
`k+1` interpretation for degree zero). Multiple output channels reuse the
same factorization. This implementation shares its QR kernel with
`lsq_univariate_spline`, so the legacy API retains its FITPACK-style input
restrictions without duplicating numerical code.

## Validation

The FPM test `test/test_interpolate_make_lsq_spline.f90` contains independent
SciPy 1.17 values for degrees 0, 1, 2, 3, 5, 9, and 15, plus multichannel
consistency, stationarity, repeated positions/knots, zero weights, rank loss,
invalid inputs, extrapolation, and compatibility with
`lsq_univariate_spline`. It can also be run standalone:

```sh
gfortran -std=f2018 -fcheck=all src/kinds.f90 src/interpolate_bspline.f90 \
    test/test_interpolate_make_lsq_spline.f90 -o test_make_lsq
./test_make_lsq
```

The optional NumPy/SciPy oracle builds a temporary GNU Fortran executable,
feeds it generated weighted least-squares problems, and compares the entire
coefficient matrix to `scipy.interpolate.make_lsq_spline(method='qr')`:

```sh
python validation/compare_make_lsq_spline.py --seed 20261101 --count 160
python validation/compare_make_lsq_spline.py --seed 8675309 --count 320
```

The two seeds passed with maximum relative coefficient discrepancy below
`5e-12` in the tests run for this change. Rank-deficient cases are handled
by the Fortran status code instead of trusting SciPy's occasional `inf`
coefficient output for unsupported knot configurations.
