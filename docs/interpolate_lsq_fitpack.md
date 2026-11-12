# FITPACK-style weighted least-squares and interpolating splines

SciFort implements two FITPACK-style `scipy.interpolate` constructors. For
SciPy's complete-knot factory, see
[`make_lsq_spline`](interpolate_make_lsq_spline.md), which shares the QR
solver but supports degrees 0–15, repeated interior knots, and zero weights.

- `lsq_univariate_spline` corresponds to `scipy.interpolate.LSQUnivariateSpline`.
  It fits a spline of specified degree to data using the **user-provided
  interior knots** and a weighted least-squares objective.
- `interpolated_univariate_spline` corresponds to
  `scipy.interpolate.InterpolatedUnivariateSpline`, for an exact fit through
  strictly increasing data sites using not-a-knot endpoint conditions.

Unlike `make_smoothing_spline`, these constructors use explicit knot or
interpolation conditions. There is no GCV parameter search or roughness
penalty. The older FITPACK `UnivariateSpline(s>0)` automatic-knot smoothing
algorithm is **not implemented** by these interfaces.

## API

```fortran
use scifort_kinds, only : dp
use scifort_interpolate, only : lsq_univariate_spline, &
    interpolated_univariate_spline, bspline_interpolator_type

real(dp) :: x(8), y(8), interior(2), weights(8)
type(bspline_interpolator_type) :: fit, exact
integer :: status

x = [0.0_dp, 0.3_dp, 0.5_dp, 0.8_dp, 1.1_dp, 1.5_dp, 2.1_dp, 3.0_dp]
y = sin(x)
interior = [0.7_dp, 1.8_dp]
weights = 1.0_dp

fit = lsq_univariate_spline(x, y, interior, w=weights, k=3, status=status)
if (status /= 0) error stop 'spline fit failed'
print *, fit%evaluate(1.2_dp)
print *, fit%derivative(1.2_dp)
print *, fit%integrate(0.5_dp, 2.0_dp)

exact = interpolated_univariate_spline(x, y, k=3, status=status)
if (status /= 0) error stop 'exact interpolation failed'
```

For matrix-valued data pass `y(n,channels)` rather than a vector and receive
`bspline_array_interpolator_type`, whose `%evaluate`, `%derivative` and
`%integrate` methods return one result per channel. Data are laid out with
samples in the **first dimension** (equivalent to SciPy's `axis=0`).

For `lsq_univariate_spline`:

- `x` must be strictly increasing and finite, with at least `k+1` samples;
  `y` must be finite and have matching sample count.
- `t` is an ordered vector of **interior knots**, without endpoint
  multiplicities. It may be empty (global polynomial fit). Each interior knot
  must lie strictly between the domain boundaries; knots may not repeat.
- `w(n)` is optional, positive and finite. The objective is
  `sum((w_i*(y_i-f(x_i)))**2)`, **so observation weights are squared**.
  This is FITPACK's convention; it differs from SciFort's
  `make_smoothing_spline`, whose objective has `w_i*(y_i-f(x_i))**2`.
- `k` is the spline degree, from 1 through 5 (default 3).
- `bbox=[left,right]` optionally extends the spline support beyond the first
  and last observations; it must contain the complete data range.
- `extrapolate=.false.` causes evaluations outside the support interval to
  return NaN. Otherwise the endpoint polynomial pieces are extrapolated.
- `status=0` means success, `status=1` invalid input and `status=2` a singular
  or numerically rank-deficient least-squares system. When status is nonzero,
  the returned spline is invalid and should not be evaluated.

For `interpolated_univariate_spline`, the supported arguments are `x`, `y`,
`k`, `status`, and `extrapolate`. It accepts degrees 1–5 and requires distinct
increasing sites. It reuses SciFort's existing `make_interp_spline` code.

The returned objects use SciFort's existing `BSpline` representation, with
complete knot vectors and non-padded coefficient arrays accessible through
`%knots()` and `%coefficients()`.

## Numerical method

At sample point `x_i`, at most `k+1` B-spline basis functions are nonzero.
SciFort constructs those local basis entries using the existing de Boor
recursion. For each observation, an **online Givens QR** update eliminates
those entries directly into an upper-triangular banded factor `R` with
half-bandwidth `k`. The factorization solves

\[
\min_c \sum_i w_i^2\left(y_i-\sum_j c_j B_{j,k}(x_i)\right)^2.
\]

`R` uses `O((k+1)*ncoeff)` storage, and the cost of factorization is
`O(n_samples * k * (k + n_channels))`. Multiple channels share the same QR
factor and have independent right-hand sides. Back substitution is banded.

This avoids the squared condition number from constructing normal equations
and does not allocate an `n_samples x ncoeff` dense design matrix. A
rank-deficient basis is rejected with status 2; for example, too many knots
in a region without observations. It does **not** automatically adjust or
remove invalid user knots.

## Differences from the full legacy SciPy API

- The FITPACK `ext` modes 1 (zero outside), 2 (raise) and 3 (constant outside)
  are not represented; use `extrapolate=.false.` for NaN outside instead.
- `bbox` is supported for the LSQ constructor, but not the exact constructor.
- Only strictly increasing distinct `x` are accepted. The SciPy FITPACK
  wrappers may accept certain non-strict data sequences in their least-squares
  routines.
- The interface supports real outputs; complex-valued FITPACK spline fitting
  has not been implemented. Matrix-valued real data is an added convenience.
- No automatic-knot `UnivariateSpline(s>0)` or `LSQUnivariateSpline` residual
  accessor is provided. Residuals can be evaluated using `%evaluate(x)`.

## Validation

Independent SciPy 1.17 reference tests cover degrees 1–5, zero or several
interior knots, positive nonuniform weights, optional extended bounds,
real multichannel data, exact interpolating splines, coefficients,
extrapolation, derivatives and integrals. Both bounds-checked and `-O3`
GNU Fortran builds are exercised. Native tests additionally cover singular
and invalid inputs, affine reproduction, weighted residual orthogonality,
and an 8,000-point/340-interior-knot fit.

```bash
python tools/interpolate_reference/compare_fitpack_lsq_scipy.py
python tools/interpolate_reference/compare_fitpack_lsq_scipy.py --seed 8675309
python tools/interpolate_reference/run_all_native_interpolation.py
```

The implementation is independent and introduces no third-party Fortran
source, FITPACK, BLAS, or LAPACK dependency. SciPy is used solely for
numerical comparisons.

## FITPACK tuple and positive-s smoothing extensions (2026-11-01)

See [the FITPACK interoperability guide](interpolate_fitpack_tck_backend.md)
for the new dependency-free tuple API (`splrep(s=0)` and fixed-knot LSQ)
and the optional `perazz/fitpack` backend for automatic knots and `s>0`.
