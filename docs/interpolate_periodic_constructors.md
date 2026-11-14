# Periodic native spline constructors

`make_splrep` and `make_splprep` in `scifort_interpolate` accept
`bc_type='periodic'` for degrees 1-5. They return the existing scalar and
array B-spline objects, respectively, and support interpolation (`s=0`)
and automatic- or supplied-knot smoothing (`s>0`) without an optional dependency.

```fortran
use scifort_kinds, only: dp
use scifort_interpolate, only: make_splrep, make_splprep, &
    bspline_interpolator_type, bspline_array_interpolator_type
type(bspline_interpolator_type) :: fit
type(bspline_array_interpolator_type) :: curve
real(dp), allocatable :: parameters(:)
integer :: status
! Initialize increasing sites u(:), closed y(:), and closed x(samples,channels).
fit = make_splrep(u, y, s=0.03_dp, bc_type='periodic', status=status)
curve = make_splprep(x, u=u, s=0.05_dp, bc_type='periodic', &
    parameters=parameters, status=status)
! Check status == 0 and the object's valid() method before evaluation.
```

## Closed input and interval

The first and last ordinates must match. For a parametric curve, every
coordinate component must match. SciFort uses the componentwise tolerance
`abs(first-last) <= 1e-12*max(1,abs(first),abs(last))`, evaluated with scaling
to avoid overflow. All coordinates, including the final row, must be finite.
This closure tolerance is a documented compatibility difference from the
current SciPy constructor's endpoint validation.

Parameter sites must be finite and strictly increasing, with more than `k`
samples. Omitted parametric `u` uses normalized Euclidean chord lengths
including the supplied closing row. The constructor does not append a point
or replace an unmatched final point. Coincident consecutive points require
supplied increasing parameters. Optional scalar `xb,xe` or parametric `ub,ue`
must equal the first and last sites. An enclosing larger interval is invalid.

These modern constructor requirements are separate from the legacy tuple
APIs: `splrep(per=.true.)` and `splprep(per=.true.)` retain their existing
ignored-final-observation and automatic-closure behavior.

## Interpolation and smoothing

For `s=0`, `w`, `t`, and `nest` remain invalid. Native periodic interpolation
matches endpoint values and derivatives through order `k-1`. Parameters and
all coordinates are normalized before the existing periodic collocation
solve; the physical knots and coefficients are checked before returning.
The result has unpadded coefficients, including the usual cyclic redundancy.

For `s>0`, weights must be finite and positive at **every** site, including
the closing site. The residual and knot selection use the first `m-1` sites;
the last point represents the same periodic location as the first. The final
weight is validated but does not enter the fit. The scalar residual is
`sum((w(:m-1)*(fit(u(:m-1))-y(:m-1)))**2)`. Parametric fitting sums that
residual across all coordinates with one common smoothing constraint.
The FITPACK relative residual tolerance is `0.001`. A sufficiently large
target produces a weighted constant or constant vector.

Optional `nest` is a hard capacity from `2*(k+1)` through `m+2*k` for
periodic smoothing. An infeasible capacity returns numerical status 3.
For positive `s`, optional `t` supplies the full periodic knot vector and
excludes `nest`. No knots are added or removed, including at the constant
limit. Accepted knots are finite and nondecreasing, with
`2*(k+1) <= size(t) <= m+2*k`, `t(k+1)` equal to the first site, and
`t(size(t)-k)` equal to the closing site. The base interval knots must be
strictly increasing. Repeated interiors are invalid. Exterior knots may be
cyclic, noncyclic, partially repeated, or clamped.

The solver minimizes squared jumps of the kth derivative at interior knots
under the same residual constraint. A sufficiently large target gives the
weighted constant on those knots. A target below the attainable residual
beyond its tolerance returns status 3. Noncyclic exterior vectors and
rank-deficient cases use the [rank-aware dense fallback](interpolate_fixed_general.md),
which enforces endpoint derivatives explicitly and has documented size,
range, and conditioning limits. This does not enable legacy periodic `task=-1`.

## Evaluation and numerical limits

As with the underlying native and SciPy interpolation APIs, interpolation
objects of degree greater than one wrap queries periodically. Degree-one
interpolation uses ordinary endpoint-polynomial extrapolation. Smoothing
objects also use ordinary endpoint-polynomial extrapolation, matching the
inspected modern SciPy smoothing result. Boundary continuity does not imply
that these smoothing objects wrap out-of-domain queries. Scalar and array
evaluation, derivatives, integrals, and knot insertion retain the existing
object methods and policies.

Periodic fitting needs headroom for extended exterior knots. The existing
periodic smoothing guards apply: site magnitudes at most
`huge/(4*(k+1))`, physical span and distinct knot gaps at least `tiny`,
normalized sample and supplied-knot gaps at least `tiny**(1/(4*(k+1)))`, and relative active
weights at least `tiny**0.125`. Normalization uses one coordinate scale,
preserving the Euclidean geometry. Unrepresentable targets, unresolved rank,
allocation failures, or a failed physical residual check return status 3.

Periodic interpolation uses the same site-magnitude and normalized-gap
guards, checks physical knot gaps, and verifies active sample values within
`512*epsilon` after common coordinate normalization. The parametric
constructor additionally retains its existing conservative interpolation
range guards, described in its guide. Dense periodic collocation targets
moderate sample counts; no large-problem performance parity is claimed.
These checks do not guarantee good conditioning for all accepted problems.

The supplied-knot path is also checked against an independent cyclic
B-spline basis and augmented least-squares constrained solve; see its
[validation record](../validation/periodic_fixed_2026-11-14.md).

## Compatibility and references

The inspected SciPy 1.18.0 public `make_splprep` interpolation shortcut omits
forwarding `bc_type`. SciFort applies periodic boundaries also for parametric
`s=0`; interpolation comparisons therefore use SciPy's explicit
`make_interp_spline(...,bc_type='periodic')` mathematical counterpart.

Automatic smoothing retains the attributed native legacy FITPACK path.
Exact modern knot-selection and iteration parity remain gaps, as do SciPy's
soft capacity behavior and repeated interior knots. Numerical validation
checks unwrapped endpoint derivatives, independent B-spline evaluation,
weighted residuals, and cyclic normal equations with a shared nonnegative
multiplier. See the [validation record](../validation/periodic_constructors_2026-11-14.md).

See the [scalar constructor](interpolate_make_splrep.md),
[parametric constructor](interpolate_make_splprep.md),
[TCK guide](interpolate_fitpack_tck_backend.md), and the official SciPy
[scalar](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.make_splrep.html)
and [parametric](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.make_splprep.html)
constructor documentation.
