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
| `t(:)` | Optional **full clamped** knot vector with simple interior knots |
| `nest` | Optional hard knot capacity, from `2*(k+1)` through `size(x)+k+1` |
| `status` | Existing FITPACK status convention: 0 success, 1 invalid, 2 unsupported, 3 numerical failure |
| `bc_type` | Default or explicit `'not-a-knot'`; `'periodic'` returns unsupported status 2 |

For `s=0`, the constructor calls the native not-a-knot interpolator. As in
SciPy 1.15.3, supplying `w`, `t`, or `nest` with `s=0` is invalid, even when
the supplied weights are ones. Bounds must equal the sample endpoints for
interpolation; unlike SciPy's interpolation shortcut, they are validated.

For `s>0`, the fit minimizes squared jumps of the kth derivative subject
to `sum((w*(fit(x)-y))**2) <= s`, with FITPACK relative tolerance `0.001`.
Automatic fitting chooses knots through the licensed native FITPACK solver.
Coefficients are unpadded: `size(c) = size(knots)-k-1`.

When `t` is supplied, it includes exactly `k+1` copies of each bound and
strictly increasing interior knots. Its length must be between `2*(k+1)`
and `size(x)+k+1`. The constructor retains that full knot vector and solves
the smoothing problem on it. This is **not** a call to `make_lsq_spline`:
the latter minimizes residual without penalizing derivative jumps. If the
target exceeds the polynomial residual, a polynomial is represented on the
supplied knots. If the least-squares residual already exceeds the target,
status 3 and an invalid object are returned. Nonclamped knots and repeated
interior knots are outside this constructor's supported knot forms.

`t` and `nest` are mutually exclusive. `nest` is a hard capacity in SciFort;
if it prevents meeting the target, status 3 is returned. SciPy instead
treats `nest` as a target that may be exceeded, so this is an intentional
compatibility limitation. No successful fit silently violates the residual
contract because of insufficient workspace.

The solver shifts coordinates before scaling, normalizes ordinates and
weights, and verifies the final physical representation. The normalized
sample and supplied-knot gaps must be at least
`tiny(1.0_dp)**(1.0_dp/(4*(k+1)))`, relative weights must be at least
`tiny(1.0_dp)**0.125`, and distinct physical knot gaps at least `tiny(1.0_dp)`.
Unrepresentable targets/spans, unresolved rank, allocation failures, and
unmet residuals return status 3. These guards do not guarantee good
conditioning for all accepted problems.

## Compatibility and remaining gaps

Automatic selection follows the native legacy FITPACK path. Modern SciPy's
`generate_knots` path can select different knots and produce different
coefficients and higher derivatives while satisfying the same mathematical
contract. Exact modern automatic-knot parity is not claimed. Supplied-knot
and interpolation behavior are independently checked against modern SciPy;
automatic fits are checked against legacy FITPACK and the smoothing equations.
See the [validation record](../validation/make_splrep_2026-11-13.md).

The current SciPy interface also supports periodic boundary selection;
this constructor does not yet implement it. Native periodic fitting remains
available through `splrep(per=.true.)`. Modern `make_splprep`, broader knot
forms, and SciPy's soft `nest` behavior remain gaps. See
[SciPy's constructor documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.make_splrep.html)
and the [native TCK guide](interpolate_fitpack_tck_backend.md).
