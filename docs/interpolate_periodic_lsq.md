# Periodic fixed-knot least squares

Native `splrep` and `splprep` support `per=.true., task=-1`, degrees 1-5,
without an optional FITPACK dependency. The result minimizes the weighted
squared residual on supplied interior knots, without a derivative-jump penalty.
This differs from positive-`s` smoothing on supplied knots in `make_splrep`
and `make_splprep`.

## Inputs and tuple conventions

Scalar `splrep` takes **interior knots only**:

```fortran
fit = splrep(u, y, w=w, k=3, task=-1, t=[0.25_dp,0.5_dp,0.75_dp], &
             per=.true., status=status)
```

Parametric `splprep` takes a **full knot vector**, as in its existing
nonperiodic interface:

```fortran
curve = splprep(points, u=u, w=w, k=3, task=-1, &
    t=[(0.0_dp,i=1,4),0.25_dp,0.5_dp,0.75_dp,(1.0_dp,i=1,4)], &
    per=.true., status=status)
```

These examples assume the base interval is `[0,1]`; initialize inputs and
check `status==0` and the returned tuple's `valid()` method before evaluation.
Parametric input has shape `(samples,coordinates)`. As before, explicitly
supply `task=-1`; supplying `t` with default `task=0` is invalid in SciFort.

Interior knots must be finite, strictly increasing, strictly between the
period endpoints, and nonempty. Their maximum count is `m-2`, where `m`
includes the closing site. The parametric full vector therefore has length
`2*k+3` through `m+2*k`. Empty-interior fixed fits are rejected, matching the
inspected FITPACK task validation. Constant data can still be fitted on a
nonempty supplied knot set.

For parametric input, entries `t(k+2:size(t)-k-1)` select the interior knots.
The full input must be finite, but endpoint and exterior template entries
are ignored and rebuilt, matching FITPACK's periodic contract. They need not
be clamped or cyclic. Input arrays are not mutated. Output endpoints equal
the first and last parameter sites, and exterior knots are canonical cyclic
extensions with the same period. The interior physical values are retained.

Scalar tuples pad coefficients with `k+1` trailing zeros. Parametric tuples
store unpadded coefficients by coordinate and return the supplied or generated
parameter sites. Cyclic redundancy is included in both coefficient layouts.
Existing `splev` extrapolation and tuple calculus behavior is preserved;
periodic boundary continuity does not enable periodic query wrapping.

## Periodic observations and parameters

Only the first `m-1` observations and weights enter the objective:
`sum((w(:m-1)*(fit(u(:m-1))-y(:m-1)))**2)` for a scalar, summed over coordinates
for a curve. The final observation and weight are ignored, including nonfinite
values, as in the existing legacy periodic fitting paths. Active observations
must be finite and active weights positive and finite.

Explicit parameter sites must all be finite and strictly increasing, including
the closing site. Without parametric `u`, existing normalized Euclidean chord
lengths include the closing leg from the last active row back to the first.
The supplied final row does not affect chord lengths. Coincident consecutive
active points or a zero closing chord still require explicit parameters.
Optional scalar `bbox` or parametric `ub,ue` must match the parameter endpoints.

`s` retains its existing finite/nonnegative validation but does not constrain
a `task=-1` fit. Changing a valid `s` leaves the least-squares result unchanged.
The result need not satisfy a smoothing residual target.

## Rank, range, and size limits

An independent cyclic B-spline design is solved with the existing native
pivoted QR and SVD least-squares routines. Data and weights use common scaling
before the solve. Endpoint values and derivatives through order `k-1` follow
from the cyclic knot/coefficient representation.

**SciFort extension:** rank-deficient observation matrices use the SVD
minimum-norm solution for the independent cyclic coefficients. Successful
fitting does not imply those coefficients are uniquely determined by the
data. SciPy comparisons are strict for the tested full-rank cases; no matching
SciPy failure behavior is claimed for rank-deficient cases.

The dense solver supports up to **256 independent cyclic coefficients**,
**10,000 sites**, and **256 coordinate channels**. With `q` interior knots,
there are `q+1` independent coefficients and `q+k+1` unpadded coefficients.
Invalid inputs return status 1; size, allocation, numerical, or range failures
return status 3. `task=1` continuation remains unavailable.

Existing periodic range guards apply: site magnitude at most
`huge/(4*(k+1))`, physical span and knot gaps at least `tiny`, normalized
sample/knot gaps at least `tiny**(1/(4*(k+1)))`, and active weights relative
to their maximum at least `tiny**0.125`. Common coordinate scaling must
preserve nonzero data and returned coefficients must be representable.
These checks do not guarantee good conditioning for every accepted problem.

See the [validation record](../validation/periodic_lsq_2026-11-14.md),
[TCK guide](interpolate_fitpack_tck_backend.md),
[modern periodic constructors](interpolate_periodic_constructors.md), and
[provenance](../CODE_PROVENANCE.md).
The public compatibility references are SciPy's
[`splrep`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.splrep.html)
and [`splprep`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.splprep.html)
documentation.
