# Fixed-knot spherical least-squares splines

`lsq_sphere_bivariate_spline` provides a native Fortran counterpart to SciPy's
`LSQSphereBivariateSpline`. It fits real scattered observations on a sphere
using a constrained bicubic tensor B-spline, with user-supplied knots and
positive residual weights. No optional backend or new dependency is needed.

For gridded interpolation or adaptive smoothing, see
[rectangular spherical splines](interpolate_rect_sphere.md), which use
the same result type with different fitting diagnostics and pole options.

```fortran
program example
    use scifort_kinds, only : dp
    use scifort_constants, only : scifort_pi
    use scifort_interpolate, only : lsq_sphere_bivariate_spline, sphere_bivariate_spline_type
    implicit none
    type(sphere_bivariate_spline_type) :: fit
    real(dp) :: theta(60), phi(60), r(60), tt(1), tp(3)
    integer :: k, status
    do k = 1, 60
        theta(k) = scifort_pi*real(mod(k-1,5),dp)/4.0_dp
        phi(k) = 2.0_dp*scifort_pi*real((k-1)/5,dp)/12.0_dp
        r(k) = cos(theta(k))
    end do
    tt = [scifort_pi/2.0_dp]
    tp = scifort_pi*[0.5_dp,1.0_dp,1.5_dp]
    fit = lsq_sphere_bivariate_spline(theta,phi,r,tt,tp,status=status)
    if (fit%valid()) then
        print *, fit%value(scifort_pi/3.0_dp,0.7_dp)
        print *, fit%rank(), fit%residual()
    end if
end program example
```

The fit is available from `scifort_interpolate` or `scifort_interpolate_sphere`.
The example approximates cosine; a finite cubic latitude spline does not
represent cosine exactly.

## Construction

```text
fit = lsq_sphere_bivariate_spline(theta,phi,r,tt,tp,w,eps,status)
```

- `theta`, `phi`, and `r` are same-length finite arrays, with at least two
  observations. Input ordering is arbitrary and repeated observations are
  allowed. Colatitudes satisfy `0 <= theta <= pi`, and training longitudes
  satisfy `0 <= phi <= 2*pi`. Coordinates and knots are in radians.
- `tt` contains strictly increasing interior colatitude knots in `(0,pi)`;
  it may be empty. `tp` contains at least three strictly increasing interior
  longitude knots in `(0,2*pi)`. Unlike SciPy's more general documented input,
  fewer longitude knots are not supported by this initial native interface.
- Optional `w` contains finite, strictly positive weights; the default is one.
  The objective is `sum((w*(r-f(theta,phi)))**2)`, so inverse measurement standard
  deviations are appropriate weights. These are not variance weights.
- Optional `eps` must be finite and in `(0,1)`, with default `1e-16`. The
  relative singular-value cutoff is `max(sqrt(eps), max(m,n)*epsilon(dp))`,
  where `m` is the observation count and `n` is the reduced coefficient count.

The reduced dimension is `n = 6 + size(tt)*(size(tp)+1)`. The current dense
backend requires `n <= 256` and `size(tp)+1 <= 256`. It also limits the summed
entries of the observation, coordinate-basis, and folded longitude design
matrices to 2,000,000. This is a workspace guard rather than a bound on total
process memory; the SVD also allocates temporary arrays. Large problems need
a future sparse/banded backend. With these bounds, fitting still has dense
SVD cost and should be reserved for modest knot counts.

Construction returns these optional statuses:

| Status | Constant | Result |
| --- | --- | --- |
| 0 | `sphere_spline_success` | Valid full-rank fit |
| 1 | `sphere_spline_invalid_input` | Invalid fit |
| 2 | `sphere_spline_rank_deficient` | Valid minimum-norm fit, with a rank warning |
| 3 | `sphere_spline_numerical_failure` | Invalid fit: solver/range failure |
| 4 | `sphere_spline_size_limit` | Invalid fit: dense size guard exceeded |

Check `%valid()` if the status is omitted. `%rank()` and
`%degrees_of_freedom()` are **SciFort extension** diagnostics; the latter
returns the reduced model dimension, not the residual statistical degrees of
freedom. A rank-deficient fit is not evidence of identifiable coefficients.

## Evaluation and metadata

```text
y = fit%value(theta,phi,dtheta,dphi,status)
y = fit%values(theta_array,phi_array,dtheta,dphi,status)
z = fit%grid(theta_array,phi_array,dtheta,dphi,status)
call fit%get_knots(full_theta_knots,full_phi_knots)
c = fit%coefficients()
fp = fit%residual()
```

`%values` evaluates paired arrays of equal length. `%grid` evaluates their
Cartesian product, returning shape `(size(theta_array),size(phi_array))`;
the coordinates need not be sorted. `dtheta` and `dphi` are optional partial
derivative orders, each in `0..2`, default zero. Derivatives are with respect
to angular coordinates in radians, without surface-gradient metric factors.

Query colatitudes must stay in `[0,pi]`. Finite query longitudes wrap modulo
`2*pi`, including negative longitudes: automatic wrapping is a **SciFort
extension**. Very large angles have the phase resolution of their binary64
input. Invalid queries or orders return NaN and status 1; an unrepresentable
rescaled value returns NaN and status 3. Array queries retain valid entries
alongside NaNs and return a nonzero status if any query fails. With multiple
different failures, the last failing query determines the array status.

`%get_knots` returns full padded knot vectors. `%coefficients` returns the
standard FITPACK tensor coefficients, with longitude varying fastest;
these may be used with a bicubic tensor spline evaluator. Both return empty
arrays on an invalid fit. `%residual` returns the weighted sum of squared
residuals in original data/weight units, or NaN for an invalid fit. It can
return positive infinity when the physical squared residual is outside the
binary64 range; this does not invalidate an otherwise representable fit.

## Constraints, compatibility, and numerical limits

The longitude coefficients repeat cyclically, giving a periodic cubic spline.
At each pole the coefficient row is constant. The adjacent row is restricted
to the span of a constant and the periodic cubic interpolants of cosine and
sine at the longitude knots, following FITPACK's spherical formulation.
This is not an unconstrained rectangular fit with its coordinates renamed.

SciFort uses a dense Jacobi SVD in the reduced coefficient coordinates.
FITPACK uses banded Givens reduction and a different rank test. Full-rank,
well-conditioned fits agree closely in the reference cases; near a rank
threshold, the estimated rank, minimum-norm coefficients, and predictions
away from observed points can differ. No bitwise, large-problem, or general
rank-deficient FITPACK equivalence is claimed.

Observations and weights are scaled internally. Construction rejects weight
ranges whose normalized positive weights underflow to zero, unrepresentable
physical coefficients, singular longitude collocation, and solver failures.
Stored scaled coefficients reduce overflow during ordinary evaluation; this
does not provide a uniform accuracy guarantee for extremely clustered knots
or every high-order derivative. The existing B-spline and SVD kernels retain
their numerical limits. See the [validation record](../validation/sphere_spline_2026-11-13.md).

Adaptive scattered-data spherical smoothing remains future work. The separate
rectangular interface supports interpolation and smoothing. Both fit real binary64
data and add no C ABI.

## Provenance

The periodic and pole coefficient constraints are adapted from `fpsphe` and
`fprpsp` in [perazz/FITPACK](https://github.com/perazz/fitpack), a modern Fortran
implementation by Federico Perini based on Paul Dierckx's FITPACK.
The adapted source is BSD-3-Clause, with the complete
[license terms](../LICENSES/perazz-fitpack-BSD-3-Clause.txt) retained.
SciFort reuses its existing B-spline evaluators and SVD solver.
See [provenance](../CODE_PROVENANCE.md) for the exact upstream commit and changes.
