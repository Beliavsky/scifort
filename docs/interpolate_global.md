# Krogh and Floater–Hormann interpolators (2026-10-31)

SciFort now implements a **self-contained Fortran** subset of
`scipy.interpolate.KroghInterpolator`, `krogh_interpolate`, and
`FloaterHormannInterpolator` in `src/interpolate_global.f90`. The symbols are
re-exported by `scifort_interpolate`. The original source is MIT licensed;
SciPy was an external numerical oracle, not a source of copied implementation.

## Public API

```fortran
use scifort_kinds, only: dp
use scifort_interpolate, only: krogh_interpolator, krogh_interpolator_type, &
    floater_hormann_interpolator, floater_hormann_interpolator_type

type(krogh_interpolator_type) :: polynomial
type(floater_hormann_interpolator_type) :: rational
integer :: status

! Hermite conditions: f(0)=0, f'(0)=2, f(1)=0.
polynomial = krogh_interpolator([0.0_dp,0.0_dp,1.0_dp], &
    [0.0_dp,2.0_dp,0.0_dp], status=status)
print *, polynomial%value(0.5_dp)                  ! 0.5
print *, polynomial%derivative(0.5_dp)             ! 0.0
print *, polynomial%derivatives([0.5_dp],der=3)    ! (order,query,channel)

! Floater-Hormann degree 1 rational interpolation.
rational = floater_hormann_interpolator([0.0_dp,0.5_dp,1.0_dp,1.5_dp], &
    [0.0_dp,1.0_dp,0.0_dp,-1.0_dp], d=1,status=status)
print *, rational%value(0.7_dp)
```

Construction generics also accept complex `y` values and rank-two ordinate
arrays shaped `(n_nodes,n_channels)`. For complex `y`, the returned types are
`complex_krogh_interpolator_type` and
`complex_floater_hormann_interpolator_type`. Node and query coordinates are
real `real(dp)`; nodes for Krogh may be in arbitrary group order, provided
all repeated occurrences of a node are contiguous. Floater–Hormann requires
a strictly monotonic set of distinct real nodes (ascending or descending).

The Krogh object's `%value(q,channel,nu)` evaluates a selected channel and
nonnegative derivative order, `%derivative(q,nu,channel)` defaults to order 1,
`%values(q(:),nu)` returns `(n_queries,n_channels)`, and
`%derivatives(q(:),der)` returns all orders `0:der-1` with shape
`(der,n_queries,n_channels)`. Omitted `der` defaults to the node count.
The `krogh_interpolate(x,y,q,nu,status)` generic provides a scalar-ordinate
convenience interface. `%nodes`, `%channels`, and `%valid` provide
introspection; constructors use `global_interp_success` (0) or
`global_interp_invalid_input` (1).

The Floater–Hormann object supplies `%value(q,channel)`, `%values(q(:))`,
`%weights()`, `%degree()`, `%nodes()`, `%channels()`, and `%valid()`.
The `d` argument defaults to **3**, just as in SciPy; pass `d=0`, `d=1`, or
`d=2` explicitly when there are fewer than four nodes. `0 <= d < n` is
required. The normalized barycentric weights are defined only up to a
common factor and may not have the same absolute normalization as SciPy.

## Algorithms

Krogh interpolation constructs Newton coefficients from **confluent divided
differences**, using repeated knots to encode successive derivatives and
factorial-divided differences. Values and derivatives are computed with a
nested differentiated Horner recurrence, without expanding high-degree
monomial coefficients. Complexity is `O(n² * channels)` to construct and
`O(n * der * channels)` to evaluate derivatives at one query.

Floater–Hormann weights are assembled from local inverse-distance products
using log-sum-exp with a shared normalization and alternating signs. Values
are evaluated through the barycentric quotient, returning stored node values
exactly at matching knots. The naive weight construction is
`O(n * (d+1)²)`; evaluations use `O(n * channels)` per query. The special
case `d=n-1` reproduces polynomial interpolation.

## Numerical scope and limits

- Floater–Hormann real nodes only; complex interpolation *values* are
  supported. Complex nodes, `poles()`, `roots()`, and `residues()` from
  SciPy's rational class are not implemented.
- Krogh interpolation can encode arbitrarily repeated Hermite nodes but
  does not accept arbitrary-rank ordinates or a non-leading interpolation
  axis. Additional output dimensions can be flattened into channels.
- The global polynomial representation may become ill-conditioned for many
  nodes. SciPy likewise cautions about polynomial interpolation beyond
  roughly 30 points. High derivatives of tightly spaced data may differ
  significantly in floating-point arithmetic even when the functions are
  mathematically identical.
- The methods do not expose an error-status argument on *evaluation*;
  construction returns a status. Invalid/nonfinite queries and negative
  derivative orders currently produce zero, rather than a SciPy exception.
- Neither algorithm introduces a BLAS/LAPACK requirement or requires
  `polyroots-fortran`. The existing opt-in root backend remains separate.

## Independent reference tests

```bash
python tools/interpolate_reference/compare_global_interpolators_scipy.py
python tools/interpolate_reference/compare_global_interpolators_scipy.py --seed 8675309
python tools/interpolate_reference/run_all_native_interpolation.py
```

The comparison suite generates Hermite conditions from known polynomials and
compares to both the generating polynomial and SciPy 1.17. It also tests
random Floater–Hormann data, extrapolation and exact nodes, four derivative
orders, real/complex multichannel arrays, both compiler optimization modes,
and the Krogh derivative-array interface. It avoids judging extrapolation
of high-degree global polynomials against a single ill-conditioned reference.
The native tests cover special cases and reject invalid inputs.
