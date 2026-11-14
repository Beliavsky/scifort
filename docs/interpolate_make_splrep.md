# Scalar smoothing constructor: `make_splrep`

`scifort_interpolate` exports `make_splrep`, returning the existing native
`bspline_interpolator_type`. It supports interpolation and bounded-residual
smoothing for degrees 1-5 without a new dependency. The result provides
`evaluate`, `derivative`, `derivative_spline`, `antiderivative`, `integrate`,
`insert_knot`, `knots`, `coefficients`, `degree`, and `valid` methods.

```fortran
use scifort_kinds, only: dp
use scifort_interpolate, only: make_splrep, bspline_interpolator_type
type(bspline_interpolator_type) :: fit
integer :: status
! Initialize strictly increasing x(:), finite y(:), and positive w(:).
fit = make_splrep(x, y, w=w, s=0.15_dp, status=status)
! Check status == 0 and fit%valid() before evaluation.
```

## Arguments and behavior

| Argument | Meaning and default |
|---|---|
| `x(:), y(:)` | Finite scalar samples; `size(x) > k`; x strictly increasing |
| `w(:)` | Positive finite residual multipliers; default all ones |
| `xb, xe` | Finite enclosing interval; default first/last x |
| `k` | Polynomial degree 1-5; default 3 |
| `s` | Nonnegative total weighted squared-residual target; default zero |
| `t(:)` | Full nondecreasing knots; interior multiplicity up to `k+1` (SciFort extension for repeats); positive `s` only |
| `nest` | Soft modern knot-count target; upper bound `size(x)+k+1`, or `size(x)+2*k` for periodic fits; large-input fallback retains a hard capacity |
| `status` | Existing FITPACK status convention: 0 success, 1 invalid, 2 unsupported, 3 numerical failure |
| `bc_type` | Default or explicit `'not-a-knot'`, or `'periodic'` |

For `s=0`, the constructor calls the native interpolator with the selected boundary. As in
SciPy 1.15.3, supplying `w`, `t`, or `nest` with `s=0` is invalid, even when
the supplied weights are ones. Bounds must equal the sample endpoints for
interpolation; unlike SciPy's interpolation shortcut, they are validated.

For nonperiodic `s>0`, the fit minimizes squared jumps of the kth derivative subject
to `sum((w*(fit(x)-y))**2) <= s`, with FITPACK relative tolerance `0.001`.
Automatic fitting uses native `generate_knots` within its size limits, then
the existing fixed-knot smoothing solver. Larger inputs retain legacy banded selection.
Coefficients are unpadded: `size(c) = size(knots)-k-1`.

When `t` is supplied, its base interval endpoints equal the bounds,
its exterior knots are nondecreasing, and interior multiplicities are at most
`k+1`. Its length must be between `2*(k+1)`
and `size(x)+k+1`. The constructor retains that full knot vector and solves
the smoothing problem on it. This is **not** a call to `make_lsq_spline`:
the latter minimizes residual without penalizing derivative jumps. If the
target exceeds the zero-penalty residual, the zero-penalty fit is represented
on the supplied knots. With simple interiors this is a polynomial. If the
least-squares residual already exceeds the target,
status 3 and an invalid object are returned. Nonclamped exterior knots and
rank-deficient regularized fits use a [documented dense fallback](interpolate_fixed_general.md).
Repeated interiors are a **SciFort extension** with a larger zero-penalty
space; see [the repeated-knot contract](interpolate_repeated_knots.md).

`t` and `nest` are mutually exclusive. Modern selection can exceed `nest`;
if the last selected vector still cannot meet the residual bound, status 3
is returned. Large-input legacy fallback retains its earlier hard capacity.
See [modern selection and soft counts](interpolate_modern_selection.md).

The solver shifts coordinates before scaling, normalizes ordinates and
weights, and verifies the final physical representation. The normalized
sample and supplied-knot gaps must be at least
`tiny(1.0_dp)**(1.0_dp/(4*(k+1)))`, relative weights must be at least
`tiny(1.0_dp)**0.125`, and distinct physical knot gaps at least `tiny(1.0_dp)`.
Unrepresentable targets/spans, unresolved rank, allocation failures, and
unmet residuals return status 3. These guards do not guarantee good
conditioning for all accepted problems.

## Compatibility and remaining gaps

Automatic selection uses the modern generator for at most 10,000 sites and 256 channels.
See [large-input selection](interpolate_large_selection.md) for banded storage
and bounded rank fallback.
Nonperiodic selected knots match the independent SciPy comparison cases.
Periodic canonicalization and root-stopping differences remain documented;
matching knots do not promise coefficients or higher derivatives that match
to roundoff. See [modern selection validation](../validation/modern_selection_2026-11-14.md)
and the earlier [constructor validation](../validation/make_splrep_2026-11-13.md).

Periodic interpolation and automatic/supplied-knot smoothing are supported; see the
[periodic constructor guide](interpolate_periodic_constructors.md) for closed
input, residual, capacity, numerical, and extrapolation contracts and accepted full periodic knot forms. The native parametric counterpart
[`make_splprep`](interpolate_make_splprep.md) is also available. Repeated
interior knots and modern selection/soft counts for larger inputs remain gaps. See
[SciPy's constructor documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.make_splrep.html)
and the [native TCK guide](interpolate_fitpack_tck_backend.md).
