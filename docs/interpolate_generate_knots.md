# Native adaptive knot generation

`generate_knots` and `knot_generator_type` are exported by
`scifort_interpolate`. They generate successive full B-spline knot vectors
for scalar or shared-channel weighted least squares, without an optional
FITPACK dependency. The degree is 1-5 (default 3).

```fortran
use scifort_kinds, only: dp
use scifort_interpolate, only: generate_knots, knot_generator_type
implicit none
type(knot_generator_type) :: gen
real(dp) :: x(6), y(6)
real(dp), allocatable :: t(:)
logical :: done
integer :: status, i

x = [(real(i-1,dp), i=1,6)]
y = sin(x)
gen = generate_knots(x,y,s=0.01_dp,status=status)
if (status == 0) then
    do
        call gen%next(t,done,status)
        if (status /= 0 .or. done) exit
        ! Use t as full knots for make_lsq_spline, or inspect the refinement.
        write(*,'(a,i0)') 'Knot count: ',size(t)
    end do
end if
```

Initialization copies the inputs. Each `next` returns one allocated vector
and advances only that generator. Assignment makes an independent copy of
the state. `done=.true., status=0` means normal exhaustion, with an allocated
zero-length vector. Repeated calls after exhaustion have the same result.
Initialization and iteration must both be checked: later solves or midpoint
rounding can fail even after successful initialization. A default-constructed
or failed state returns `done=.true., status=1`. An iteration failure returns
`done=.true., status=3` and invalidates the state. Failure output is empty.

## Arguments and residual convention

```text
gen = generate_knots(x,y,w=...,xb=...,xe=...,k=...,s=...,nest=...,
                     bc_type=...,status=...)
call gen%next(t,done,status)
```

- `x` is finite and strictly increasing, with more than `k` sites.
- `y` is finite, either `(sites)` or `(sites,channels)`.
- `w` contains finite positive residual weights, default one. The target is
  `sum((w(i)*(fit(i,j)-y(i,j)))**2)` over samples and channels. Weights are
  squared in this objective. This differs from `make_smoothing_spline`.
- `s` is finite and nonnegative, default zero. Positive targets accept an
  LSQ residual below `s`, or within `0.001*s` above it.
- `xb` and `xe` default to the first and last sites. For positive-target
  nonperiodic selection they may enclose a larger interval. Zero-target
  interpolation uses the actual site endpoints, as SciPy does; supplied
  bounds are still validated as finite and enclosing.
- `bc_type` is `'not-a-knot'` (default) or `'periodic'`.
- `nest` is an optional soft stopping count, at least `2*(k+1)`. If the
  residual still exceeds the target at the initial minimum count, refinement
  can return a vector exceeding `nest`. It does not promise a satisfied
  residual target when stopping for this reason.

With `s=0`, weights and `nest` are rejected. Exactly one interpolation knot
vector is returned. Odd nonperiodic degrees use data-site interior knots;
even degrees use adjacent-site midpoints with the end pairs omitted according
to the not-a-knot convention. This is a knot generator, not a spline fit.

With `s>0`, nonperiodic selection starts with no interior knots. It repeatedly
fits weighted least squares, splits residual sums among knot intervals, and
inserts a middle data site from the largest eligible interval. Several knots
may be inserted between outputs. Each yielded vector is longer than its
predecessor. At the interpolation count, the final vector switches to the
interpolation convention, including midpoint knots for even degrees.

## Periodic selection

Periodic input requires matching first and last observations, component by
component within `1e-12*max(1,abs(first),abs(last))`. Both bounds must equal the
actual site endpoints. The final observation and weight are validated but
excluded from the positive-target residual. They are not independent samples.

Selection first tests a constant periodic fit. If that suffices, it returns
the `2*(k+1)` knot vector for that fit and stops. Otherwise it starts with one
interior data-site knot and refines a cyclic coefficient basis. Positive-target
vectors always have canonical exterior knots, including at a `nest` stop.
Their interpolation limit is `size(x)+2*k` knots. Even degrees shift the
interior interpolation sites to the periodic midpoint convention.

For `s=0, k=1`, SciPy's special clamped linear knot convention is retained.
For other periodic zero-target degrees, the single vector is cyclic. Ordinary
`make_lsq_spline` does not impose periodic coefficient constraints: periodic
generated knots must be paired with a periodic fitter. Legacy periodic
`splrep(task=-1)` takes only the interior slice
`t(k+2:size(t)-k-1)`, while `splprep(task=-1)` takes full knots; their constant
case with no interior knots is not accepted. See
[periodic fixed-knot least squares](interpolate_periodic_lsq.md).

## Compatibility and limits

The selection control flow and knot conventions adapt the exact SciPy
v1.18.0 implementation under BSD-3-Clause. Native least-squares evaluation
uses existing nonperiodic banded QR and the attributed periodic FITPACK QR
for larger fits, with bounded native QR/SVD fallbacks. Nonperiodic full
sequences are independently compared with SciPy 1.15.3; periodic checks use pinned v1.18.0 selection control flow
and independent NumPy or sparse LSQR cyclic residual evaluators.

Intentional differences and limits:

- Strictly increasing sites, positive weights, degrees 1-5, at most 10,000
  sites and 256 channels. SciPy also admits some zero-weight and repeated-site
  cases. Those are not supported here.
- Periodic closure uses the explicit tolerance above, rather than NumPy's
  `allclose` defaults.
- Cyclic exterior knots are refreshed after every insertion and before every
  residual solve. SciPy v1.18.0's Python loop leaves them stale within a batch
  and can return stale exteriors at a capacity stop. Periodic selection paths
  can therefore differ; exact upstream periodic sequence parity is not claimed.
- A minimum-norm QR/SVD fallback for rank-deficient intermediate fits is a
  **SciFort extension**. It is retained for at most 256 independent
  coefficients and one million active sample-by-coefficient entries.
  Larger rank failures return status 3 without an unbounded dense allocation.
  Matching SciPy failure behavior is not promised.
- Sites, data, weights and the target are normalized for solving. A physical
  interval must be at least `tiny(dp)`; bounds need headroom of
  `huge(dp)/(4*(k+1))`. Normalized adjacent sites must differ by at least
  `tiny(dp)**(1/(4*(k+1)))`. Active relative weights must be at least
  `tiny(dp)**0.125`; normalization must preserve nonzero data. A positive
  normalized target below `tiny(dp)` is rejected. Normalized coefficients
  exceeding `1e100` are rejected. These are range guards, not conditioning or
  forward-error guarantees. Unrepresentable/collapsed output knot intervals
  fail explicitly.
- Floating-point residual ties or near-threshold fits can change selection
  paths. A yielded vector is a proposal; it does not certify attainable
  accuracy below floating-point resolution.

Status `1` denotes invalid input; `3` denotes a size, range, representation,
or numerical failure. Status `0` denotes successful initialization, a returned
vector, or normal exhaustion as distinguished by `done`.

`make_splrep` and `make_splprep` now use this generator for automatic smoothing
within its size limits, followed by their existing fixed-knot smoothing solve.
See [constructor integration and soft counts](interpolate_modern_selection.md)
for larger-input fallback and infeasible-target behavior. See also
[banded large-input selection](interpolate_large_selection.md). C ABI and
wrappers remain deferred. See the [validation record](../validation/generate_knots_2026-11-14.md).
