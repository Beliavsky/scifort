# Modern automatic selection in spline constructors

For positive `s` with no supplied `t`, `make_splrep` and `make_splprep` now
use native [`generate_knots`](interpolate_generate_knots.md) for up to 10,000
sites and 256 channels. The constructors consume the sequence, retain its
last full knot vector, then run the existing derivative-jump smoothing solve
on those knots. This covers scalar/shared-channel, nonperiodic/periodic fits
and degrees 1-5. No public signature or C ABI is added.

Interpolation (`s=0`) retains the existing interpolation path. Supplied `t`
retains the existing fixed-knot path, including its accepted exterior knots,
rank-aware fallback, polynomial/constant limits and invalid-input behavior.
The legacy tuple `splrep` and `splprep` algorithms also retain their existing
selection paths. Automatic constructor results can now differ from those
legacy tuple fits because they can choose different knots.

## Soft knot-count targets

`nest` remains mutually exclusive with `t` and invalid for `s=0`. Its accepted
range remains `2*(k+1)` through `m+k+1`, or through `m+2*k` for periodic fits.
Within the modern selection limits it is a soft target rather than a storage
capacity: the generator can return more knots than `nest`.

For example, if a polynomial fit fails the residual bound but adding one
interior knot makes the bound feasible, a nonperiodic call with
`nest=2*(k+1)` can succeed with `2*(k+1)+1` knots. The analogous periodic
constant-to-first-interior refinement can also exceed the supplied count.
The tests cover both cases for degrees 1-5, scalar and shared-channel data.

This does not guarantee that any small `nest` will permit the requested
residual. If the final selected knots cannot meet `s` within tolerance, the
constructor returns status 3, an invalid object, and no `parameters` output.
SciFort does not turn an infeasible bounded-error request into a successful
approximation with an unmet residual. SciPy can instead return a spline
after a root-search warning in such cases; failure behavior is not identical.

Above 10,000 sites or 256 channels, the established banded automatic path
remains available to preserve prior coverage. That path retains its earlier hard-capacity
`nest` behavior and legacy knot selection. Modern selection and soft-count
coverage for those larger problems remain gaps. No fallback to legacy
selection hides a range or numerical failure for an input within the modern
size limits.

## Numerical and compatibility details

The weighted objective and residual tolerance are unchanged: weights multiply
residuals before squaring, one target is shared across all channels, and the
relative residual tolerance is `0.001`. A sufficiently loose target returns
a least-squares polynomial, or a periodic constant, without forcing equality
with `s`. Periodic input closure, ignored closing residual/weight, parameter
generation, coefficient layouts and extrapolation policies retain their
documented contracts.

At the interpolation knot-count limit, the generator uses the actual sample
endpoints, even when positive-target nonperiodic bounds enclose a larger
interval. The smoothing solve uses those selected base endpoints. Queries
outside the resulting base interval use the ordinary extrapolation policy.

All generator range guards apply, followed by the existing fixed-knot
smoothing checks. A successful knot selection is not a guarantee that the
smoothing solve will succeed at an unrepresentably small target or on an
ill-conditioned problem. No returned fit silently violates the residual bound.

Periodic vectors are canonicalized after every insertion. This intentionally
differs from SciPy v1.18.0's stale exteriors within refinement batches and at
some capacity stops, so exact periodic selection-path parity is not claimed.
Nonperiodic selected knots match the independent SciPy cases in validation.

Even with matching knots, coefficients and higher derivatives need not match
SciPy to roundoff. Existing scalar/parametric smoothing root searches can
stop at different residuals within the permitted band. The validation checks
the requested target, native B-spline representation, common-multiplier normal
equations, and a residual-dependent coefficient-difference bound. Raw SciPy
value/derivative differences are reported separately; large relative
differences in high derivatives or extrapolation are not a claim of strict
numerical parity. See the
[validation record](../validation/modern_selection_2026-11-14.md).

See [banded large-input selection](interpolate_large_selection.md) for storage,
rank fallback and checks through 10,000 sites.

Repeated interior knots, modern selection above those limits, broader ill-
conditioned behavior, exact root-iteration parity, C ABI and wrappers remain
outside this addition.
