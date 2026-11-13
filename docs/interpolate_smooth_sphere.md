# Scattered-data spherical smoothing

`smooth_sphere_bivariate_spline` is a native binary64 counterpart to SciPy's
`SmoothSphereBivariateSpline`. It fits irregularly located, weighted real
observations using FITPACK's adaptive bicubic spherical spline solver.
Longitude is periodic, pole values are single-valued, and the adjacent
coefficient rows satisfy FITPACK's first-derivative pole constraints.
No optional backend or mandatory dependency is needed.

```fortran
program example
    use scifort_kinds, only : dp
    use scifort_constants, only : scifort_pi
    use scifort_interpolate, only : smooth_sphere_bivariate_spline, &
        sphere_bivariate_spline_type
    implicit none
    type(sphere_bivariate_spline_type) :: fit
    real(dp) :: theta(120), phi(120), r(120)
    integer :: i, status
    do i = 1, size(theta)
        theta(i) = scifort_pi*modulo(real(i,dp)*0.618033988749895_dp,1.0_dp)
        phi(i) = 2*scifort_pi*modulo(real(i,dp)*0.414213562373095_dp,1.0_dp)
        r(i) = cos(theta(i))+0.2_dp*sin(theta(i))*cos(phi(i)) &
            +0.02_dp*sin(7*theta(i)+3*phi(i))
    end do
    fit = smooth_sphere_bivariate_spline(theta,phi,r,s=0.1_dp,status=status)
    if (fit%valid()) then
        print *, fit%value(0.8_dp,0.3_dp)
        print *, fit%residual(), fit%rank(), status
    end if
end program example
```

The constructor and result type are exported by `scifort_interpolate` and
`scifort_interpolate_sphere`. For user-selected knots, use
[fixed-knot least squares](interpolate_sphere.md); for rectangular data and
explicit pole options, use [rectangular spherical splines](interpolate_rect_sphere.md).

## Construction

```text
fit = smooth_sphere_bivariate_spline(theta,phi,r,w,s,eps,ntest,npest,status)
```

- `theta`, `phi`, and `r` are same-length finite arrays with at least two
  observations. Colatitudes are in `[0,pi]`, longitudes in `[0,2*pi]`, and
  coordinates are in radians. Input ordering is arbitrary and repeated sites
  are allowed. These are scattered pairs, rather than coordinate axes.
- Optional finite weights `w > 0` default to one. The objective is
  `fp = sum((w*(r-f(theta,phi)))**2)`. Inverse measurement standard deviations
  are appropriate weights; these are not variance weights.
- Optional finite `s >= 0` defaults to zero, matching the SciPy high-level
  class. It is a weighted squared-residual target in original units, not a
  penalty multiplier. Automatic knot selection and rational smoothing-parameter
  iteration normally achieve `abs(fp-s) <= 0.001*s`. A sufficiently large
  target returns a two-parameter polynomial with a smaller residual.
- Optional `eps` is finite and in `(0,1)`, default `1e-16`. The banded solver
  treats triangular diagonals at or below `sqrt(eps)*max_diagonal` as deficient.
  This is FITPACK's diagonal test, rather than the SVD test used by SciFort's
  separate fixed-knot constructor.
- `ntest` and `npest` are **SciFort extension** controls for the maximum full
  latitude/longitude knot counts, including padding. Both default to
  `8+int(sqrt(real(size(theta)/2,dp)))`, with integer division inside the root,
  matching SciPy 1.15.3's default capacities. Both must be at least eight.
  A nonpolynomial fit needs at least nine latitude and eleven longitude knots.
  Larger capacities can help when a small target exhausts the defaults.

The solver starts from a polynomial with flat poles, then adds knots according
to residual-weighted coordinate locations. Its positive smoothing-parameter
iteration uses a relative residual tolerance of 0.001 and at most 20 steps.
Each constructor call starts fresh; warm starts, fixed-knot mode, and custom
iteration/tolerance controls are not exposed here.

`s=0` requests interpolation but does not guarantee it for arbitrary data or
the default capacities. Conflicting repeated observations cannot be interpolated,
and very small targets may hit knot, conditioning, or iteration limits.
Failures return an invalid object; partially converged fits are not accepted.

## Evaluation, diagnostics, and statuses

The shared `sphere_bivariate_spline_type` provides scalar `%value`, paired
`%values`, Cartesian `%grid`, `%get_knots`, `%coefficients`, and `%residual`.
Their signatures and derivative/query behavior follow the
[spherical evaluation API](interpolate_sphere.md#evaluation-and-metadata).
Derivative orders independently range from zero through two. Finite query
longitudes wrap modulo `2*pi`, a **SciFort extension**; training longitudes
still require the range above.

`%residual()` reports the actual weighted squared residual in original
data/weight units. It can return positive infinity if the physical squared
residual is outside binary64 range; this does not invalidate representable
spline values. `%rank()` and `%degrees_of_freedom()` are **SciFort extension**
diagnostics. The latter is the reduced model dimension, not statistical
residual degrees of freedom: two in the initial polynomial case, otherwise
`6+(size(full_theta_knots)-8)*(size(full_phi_knots)-7)`.

`%rank()` returns -1 for the initial polynomial, whose rank is not estimated.
Otherwise it reports the final triangular system's estimated rank, including
penalty rows when smoothing is active. It is not necessarily the rank of the
unpenalized observation matrix. A deficient fit is minimum norm in FITPACK's
reduced coefficient coordinates, which does not establish identifiable
coefficients or unique predictions away from observed sites.

| Status | Meaning |
| --- | --- |
| 0 | Valid fit |
| 1 | Invalid input |
| 2 | Valid fit with estimated rank below model dimension |
| 3 | Numerical/range failure, unattainable target, overlapping candidate knots, or failed iteration |
| 4 | Knot capacity or workspace guard exceeded |

SciPy 1.15.3's high-level class rejects its rank-warning return codes. SciFort
retains valid deficient fits with status 2, consistent with its existing
fixed-knot API. Check `%valid()` when omitting `status`.

## Numerical and compatibility limits

The constructor caps selected real/integer workspace, coefficient/input
storage, and knot/panel arrays at 2,000,000 entries before allocating solver
arrays. This is a resource guard, not a limit on total process memory. The
banded solver has no separate 256-coefficient cap. Extremely large capacity
requests are checked with bounded integer arithmetic before workspace products.

Observations and weights are normalized, and the target is transformed in
logarithms. Positive normalized weights that underflow to zero, positive scaled
targets below `tiny(dp)`, nonfinite solver results, and unrepresentable physical
coefficients return status 3. Very large scaled targets select the polynomial
limit without overflowing. These guards do not guarantee accuracy for every
clustered, weakly constrained, or ill-conditioned data set.

Full-rank reference fits agree closely with SciPy. Rank thresholds, finite
precision, and adaptive knot choices can affect results; bitwise or general
rank-deficient SciPy equivalence is not claimed. The adapted minimum-norm path
preserves the transformed observation vector and recomputes its actual residual
before accepting a knot set. An independent dense SVD verifies a rank-four
equatorial-data fixture, including predictions away from the observations.
See the [validation record](../validation/smooth_sphere_2026-11-13.md).

## Provenance

The adaptive solver and rank-reduction helpers are adapted from Federico
Perini's [perazz/FITPACK](https://github.com/perazz/fitpack), based on Paul
Dierckx's FITPACK. Existing native FITPACK numerical helpers and spline
evaluators are reused. BSD-3-Clause notices and the complete
[license](../LICENSES/perazz-fitpack-BSD-3-Clause.txt) are retained. The source
pin and significant corrections are recorded in [provenance](../CODE_PROVENANCE.md).
No C ABI is added.
