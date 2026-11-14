# Parametric smoothing constructor: `make_splprep`

`scifort_interpolate` exports `make_splprep`, returning the existing native
`bspline_array_interpolator_type`. Coordinates have shape `(samples,channels)`;
all channels share one knot vector and one total weighted residual target.
Degrees 1-5, interpolation, automatic-knot smoothing, and supplied-knot
smoothing require no optional FITPACK dependency.

```fortran
use scifort_kinds, only: dp
use scifort_interpolate, only: make_splprep, bspline_array_interpolator_type
type(bspline_array_interpolator_type) :: curve
real(dp), allocatable :: parameters(:)
integer :: status
! Initialize finite coordinates x(samples,channels).
curve = make_splprep(x, s=0.15_dp, parameters=parameters, status=status)
! Check status == 0 and curve%valid() before evaluation.
```

## Arguments and result

| Argument | Meaning and default |
|---|---|
| `x(:,:)` | Finite coordinates, at least `k+1` samples and one channel |
| `w(:)` | Positive finite residual multipliers, one per sample; default ones |
| `u(:)` | Finite strictly increasing parameters; default normalized Euclidean chord lengths |
| `ub, ue` | Finite enclosing parameter bounds; default first/last parameter |
| `k` | Degree 1-5; default 3 |
| `s` | Nonnegative total weighted squared-residual target; default zero |
| `t(:)` | Full clamped knots with simple interiors; nonperiodic positive `s` only |
| `nest` | Hard knot capacity; upper bound `size(x,1)+k+1`, or `size(x,1)+2*k` for periodic fits |
| `status` | 0 success, 1 invalid input, 2 unsupported option, 3 numerical/allocation failure |
| `parameters(:)` | Optional allocatable output, allocated only on success; receives supplied or calculated sites |
| `bc_type` | Default or `'not-a-knot'`, or `'periodic'` |

The `parameters` output is the Fortran equivalent of SciPy's second return
value. It is separate from the read-only optional `u` input. Inputs are never
modified. Invalid calls return an invalid object and leave `parameters`
unallocated, including when it was allocated on entry.

The result provides `valid`, `channels`, `evaluate`, `derivative`, `integrate`,
`insert_knot`, `knots`, and `coefficients`. Scalar queries return a channel
vector; array queries return `(queries,channels)`. Coefficients are unpadded,
with shape `(size(knots)-k-1,channels)`. Supporting more than ten channels
preserves the existing SciFort extension to legacy FITPACK's channel limit.

## Interpolation and smoothing

For `s=0`, the constructor uses native interpolation with the selected boundary. As in SciPy
1.15.3, supplying `w`, `t`, or `nest` is invalid in this mode. Interpolation
bounds must equal the parameter endpoints; SciPy's interpolation shortcut
ignores supplied bounds. Consecutive coincident samples cannot supply an
automatic parameterization, but are allowed with strictly increasing `u`.
Automatic chord lengths use common coordinate scaling and a stable norm.

For nonperiodic `s>0`, the solver minimizes the sum of squared jumps of the kth
derivative over all channels, subject to
`sum((spread(w,2,channels)*(curve(u)-x))**2) <= s`. The FITPACK relative
residual tolerance is `0.001`. This is one shared constraint, not a separate
target for each channel. A sufficiently large target yields a degree-k
least-squares polynomial in each coordinate.

Supplied `t` must include exactly `k+1` copies of each bound and strictly
increasing simple interior knots; its length lies between `2*(k+1)` and
`size(x,1)+k+1`. It requests derivative-jump smoothing on those knots,
not unconstrained least squares. A polynomial limit retains the supplied
full knots. A target below the achievable least-squares residual, unresolved
rank, or insufficient `nest` returns status 3 and an invalid object.

Smoothing shifts parameters before scaling and normalizes all coordinates
by one common factor, preserving the Euclidean geometry. Weights are
normalized separately. The returned physical knots and coefficients must
meet the residual contract. Normalized parameter and supplied-knot gaps
must be at least `tiny(1.0_dp)**(1.0_dp/(4*(k+1)))`, relative weights at
least `tiny(1.0_dp)**0.125`, and distinct physical knot gaps at least
`tiny(1.0_dp)`. Unrepresentable spans/targets and failed numerical checks
return status 3. These guards do not guarantee good conditioning.

Interpolation additionally rejects parameter magnitudes above `huge/4`,
spans below `tiny`, unresolved normalized gaps, and coordinate magnitudes
above `huge/(4*samples*(k+1))`, reserving range for the existing interpolator.
Adjacent interpolation parameters must also be separated by at least `tiny`;
this avoids overflowing intermediate basis arithmetic on subnormal intervals.

## Compatibility limits

Automatic smoothing uses the existing licensed FITPACK parametric solver.
Modern SciPy's knot generator can choose different knots. Its root-search
initialization also differs from the legacy parametric path; the new
fixed-knot mode uses reciprocal-mean initialization as in the scalar solver.
Residual tolerances can still produce differences in coefficients and
especially higher derivatives. Exact modern numerical parity is not claimed.
Validation checks interpolation against SciPy and verifies independent
B-spline evaluation, residuals, and shared-multiplier smoothing equations.
See the [validation record](../validation/make_splprep_2026-11-14.md).

`nest` is a hard capacity in SciFort; SciPy treats it as a target that may
be exceeded. Periodic interpolation and automatic smoothing are supported;
see the [periodic constructor guide](interpolate_periodic_constructors.md)
for closed input, residual, capacity, numerical, and extrapolation contracts.
Supplied periodic knots, broader knot forms, exact modern automatic selection,
and soft capacities remain gaps. Legacy periodic fitting remains available
through `splprep(per=.true.)`.

See [SciPy's constructor documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.make_splprep.html),
the [scalar constructor](interpolate_make_splrep.md), and the
[native TCK guide](interpolate_fitpack_tck_backend.md).
