# Two-dimensional spline interpolation and filtering

The `scifort_signal` module exports three procedures corresponding to
`scipy.signal.cspline2d`, `qspline2d`, and `spline_filter`.

```fortran
use scifort_signal, only : dp, cspline2d, qspline2d, spline_filter, &
    spline_status_success
real(dp), allocatable :: image(:,:), coeff(:,:), smooth(:,:)
integer :: status

! image is an already populated rank-two real(dp) array.
coeff = cspline2d(image, status=status)
if (status /= spline_status_success) error stop 'cspline2d failed'

coeff = cspline2d(image, lamb=0.5_dp, precision=1.0e-9_dp)
coeff = qspline2d(image)
smooth = spline_filter(image, lmbda=5.0_dp)
```

All interfaces accept `real(dp)` or `complex(dp)` rank-two arrays and return
allocatable arrays of the same type. The spline-coefficient routines accept
optional `lamb`, `precision`, and `status` arguments; `spline_filter` accepts
optional `lmbda` and `status`. Default spline coefficient `lamb` is zero.
Default `spline_filter` `lmbda` is five. `precision` defaults to `1e-6`,
matching SciPy's double-precision boundary-series truncation. A negative or
at-least-one precision value selects that default, as in SciPy.

The implementation applies two separable symmetric IIR passes. Ordinary
cubic and quadratic interpolation use stable first-order sections;
smoothed cubic coefficients use second-order sections. `spline_filter`
performs cubic coefficient smoothing followed by the separable FIR kernel
`[1, 4, 1]/6` along both dimensions. With `lamb <= 1/144`, `cspline2d`
uses the ordinary cubic interpolation branch, matching SciPy.

## Evaluation at arbitrary coordinates

The coefficient arrays returned by `cspline2d` and `qspline2d` can be
sampled using `cspline2d_eval` and `qspline2d_eval`, respectively.
The evaluators use tensor-product cardinal B-splines with mirror-symmetric
reflections at image boundaries, corresponding to SciPy's
`scipy.ndimage.map_coordinates(coefficients, coordinates, order=2 or 3,
mode="mirror", prefilter=False)`.
Unlike SciPy, **coordinates are one-based**: `(1,1)` is the first pixel.
Pass *coefficient arrays*, not unfiltered image samples, unless the latter
are already spline coefficients.

```fortran
use scifort_signal, only : dp, cspline2d, cspline2d_eval, &
    cspline2d_grid, spline_status_success
real(dp), allocatable :: image(:,:), coeff(:,:), points(:), grid(:,:)
real(dp) :: row(3), col(3), point
integer :: status

! image has previously been initialized with finite sampled data.
coeff = cspline2d(image, status=status)
if (status /= spline_status_success) error stop 'coefficient failure'

row = [1.0_dp, 2.25_dp, 10.5_dp]
col = [1.0_dp, 7.50_dp, -3.0_dp]
points = cspline2d_eval(coeff, row, col, status=status)
point = cspline2d_eval(coeff, 2.25_dp, 7.50_dp)

! Evaluate on all combinations of row and column coordinates.
grid = cspline2d_grid(coeff, row, col)
```

Each `*_eval` generic supports scalar coordinate pairs, equal-length rank-one
arrays of paired coordinates, and equal-shape rank-two coordinate arrays. The
`*_grid` routines take independent rank-one row/column coordinate vectors and
return a rank-two grid with shape `(size(row), size(col))`. All overloads accept
real or complex double-precision coefficient matrices. Noncontiguous array
sections and positions outside the sample range are supported. Both singleton
axes are handled without division by zero.

An invalid coefficient array, nonfinite coordinate, or mismatched paired-query
shape returns `spline_status_invalid_input` (when `status` is supplied) and an
empty array; scalar overloads return zero with an error status. Empty coordinate
lists are allowed and result in empty outputs with a success status. These
conventions apply to both cubic and quadratic interpolation.

The evaluators are independent SciFort implementations of the compact-support
B-spline basis and mirror-index arithmetic. They introduce no runtime
libraries and use SciPy solely for independent numerical comparison.
`test/test_signal_spline2d_eval.f90` tests the public API, with
`tools/signal_reference/compare_signal_spline2d_eval_scipy.py` checking 256
seeded real and complex cases (17,920 coordinates). SciPy 1.17 agrees within
`9.16e-16` maximum absolute error on these cases.

## Differences and limits

- Both real and complex `dp` input are supported throughout. For complex
  smoothed splines, SciFort applies the same real-coefficient linear operator
  to both complex components. SciPy 1.17 rejects the corresponding complex
  `cspline2d` case and historically converts complex128 `spline_filter`
  input to complex64; SciFort does neither reduction.
- Negative smoothing parameters, smoothing quadratic coefficients, NaN/Inf
  samples, and empty image dimensions return `spline_status_invalid_input`
  with an empty `(0,0)` result. For a one-pixel interpolation axis, SciFort
  preserves the pixel. SciPy can fail to converge for very short axes.
- Tolerances for symmetric IIR initial conditions are finite, so narrow
  images or aggressive user-specified `precision` can differ from a
  mathematically exact infinite-series solution.
- The routines are intended for stable, finite inputs and do not use any
  external LAPACK, BLAS, or FPM package.

## Verification

`test/test_signal_spline2d.f90` checks fixed SciPy outputs, real/complex
linearity, selected noncontiguous sections, and invalid arguments.
`tools/signal_reference/compare_signal_spline2d_scipy.py` compiles the
standalone Fortran driver and runs 44 seeded SciPy 1.17 reference cases
covering cubic and quadratic coefficients, positive cubic smoothing,
`precision` choices, and separable filtering. Complex smoothing/filtering
that SciPy cannot directly represent is checked against the separate real
and imaginary SciPy calculations. The reference script uses a tolerance of
`5e-7` for coefficient arrays and `2e-5` for `spline_filter` outputs, the
latter reflecting SciPy's internal single-precision FIR coefficients.

The numerical recurrences are adapted from the SciPy spline filtering
algorithms (BSD-3-Clause), with attribution in `CODE_PROVENANCE.md` and
`THIRD_PARTY_LICENSES.md`.
