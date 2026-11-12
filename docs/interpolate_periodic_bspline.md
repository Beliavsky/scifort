# Periodic interpolating B-splines (`make_interp_spline`)

SciFort implements `scipy.interpolate.make_interp_spline(..., bc_type='periodic')`
for **real scalar and multichannel ordinates**, using the existing
`scifort_interpolate` and `scifort_interpolate_bspline` public modules.
The implementation needs no FITPACK, external BLAS/LAPACK or other dependencies.

## Example

```fortran
use scifort_kinds, only : dp
use scifort_interpolate, only : bspline_interpolator_type, make_interp_spline
implicit none
real(dp), parameter :: pi = acos(-1.0_dp)
real(dp) :: x(17), y(17), period
integer :: i, status
type(bspline_interpolator_type) :: spline

do i = 1, 17
    x(i) = real(i-1,dp)*2.0_dp*pi/16.0_dp
end do
y = sin(x)
y(17) = y(1) ! periodic endpoint values must match
spline = make_interp_spline(x,y,k=3,bc_type='periodic',status=status)
if (status /= 0) error stop 'periodic spline construction failed'
period = x(17)-x(1)
print *, spline%evaluate(0.4_dp)
print *, spline%evaluate(0.4_dp+2.0_dp*period) ! same value
print *, spline%derivative(0.4_dp,1)
print *, spline%integrate(x(1),x(17))
```

For multiple ordinates, pass `y(n,m)` with samples in the first axis.
The returned `bspline_array_interpolator_type` uses a common knot vector
and factors **one collocation system** for all `m` columns.

## Algorithm

For degree `k` (2–15), the periodic knot vector contains `n + 2*k`
knots, and the interpolant has `n + k - 1` B-spline coefficients per
channel. A square collocation system imposes `n` interpolation values
and equality of endpoint derivatives of orders `1` through `k-1`.
The knots extend sample gaps cyclically, with the standard half-step
interior shift for even degrees. The period is `x(n)-x(1)`.

By default, evaluation and derivatives wrap queries onto the base interval.
Setting `extrapolate=.false.` disables wrapping and returns NaN for
out-of-domain evaluations. `%integrate()` respects periodic extension
and supports intervals crossing one or more periods. `%insert_knot()`
retains exact spline values when applied inside the base interval.

For compatibility with SciPy, degrees 0 and 1 use their usual piecewise
constant/linear interpolation construction instead of the periodic
collocation system; these low-order splines do **not** automatically wrap.

Input rules:

- `x` must be finite, strictly increasing and contain at least two samples.
- Observations must be finite, with matching first/last values in every
  channel (SciFort tolerance: `1e-12 * max(1,maxval(abs(y)))`).
- Explicit `knots=` is unsupported with periodic boundary conditions, as
  in SciPy. Requested degree must be in 0–15.
- `status=0` indicates success, `1` an invalid input and `2` a singular
  collocation system. An unsuccessful constructor returns an invalid spline.

The current implementation uses dense partial-pivoting Gaussian
elimination, with cubic factorization cost in the number of coefficients.
It is intended for moderate-sized one-dimensional datasets; a future
banded/cyclic solver would improve scaling for large datasets.

## Independent validation

```sh
python tools/interpolate_reference/compare_periodic_bspline_scipy.py \
    --cases 144 --mode both
python tools/interpolate_reference/run_all_native_interpolation.py
```

The reference comparison independently constructs SciPy 1.17.0 splines,
and checks knots, coefficients, periodic evaluation, derivatives and exact
integrals for irregular nodes, multiple channels and degrees 1–15.
Two independent seeds with 150 and 120 randomized cases each were checked
under both `-fcheck=all` and `-O3` builds. The maximum scale-normalized
coefficient difference was less than `8e-11` and the maximum derivative
difference less than `2e-12` on these cases.
The native test checks endpoint derivative continuity, multichannel
agreement, knot refinement, integration across multiple periods, and input
validation. The implementation does not import or link against SciPy.

The standard knot extension and periodic collocation convention were
adapted from SciPy 1.17.0's BSD-3-Clause `scipy/interpolate/_bsplines.py`;
see `CODE_PROVENANCE.md` and `THIRD_PARTY_LICENSES.md` for notices.
