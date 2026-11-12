# Barycentric polynomial interpolation (2026-10-31)

SciFort implements `scipy.interpolate.BarycentricInterpolator` and the
one-dimensional `scipy.interpolate.barycentric_interpolate` convenience
function in `src/interpolate_barycentric.f90`. They are re-exported by
`scifort_interpolate`. The implementation is original MIT-licensed SciFort
code: no SciPy or other upstream source has been copied.

## API

```fortran
use scifort_kinds, only : dp
use scifort_interpolate, only : barycentric_interpolator, &
    barycentric_interpolator_type, barycentric_interpolate

type(barycentric_interpolator_type) :: p
real(dp), allocatable :: result(:,:), values(:)
integer :: status
p = barycentric_interpolator([-1.0_dp,0.0_dp,1.0_dp], &
    [1.0_dp,0.0_dp,1.0_dp],status=status)
result = p%values([0.25_dp,0.5_dp]) ! (queries, channels)
print *, p%value(0.25_dp)           ! 0.0625
print *, p%derivative(0.25_dp)      ! 0.5
values = barycentric_interpolate([-1.0_dp,0.0_dp,1.0_dp], &
    [1.0_dp,0.0_dp,1.0_dp],[0.25_dp,0.5_dp])
```

The constructor generic also accepts rank-two ordinates `(n_nodes,n_channels)`
and complex `y` values. Complex inputs produce a
`complex_barycentric_interpolator_type` object; interpolation nodes and query
coordinates must be real. Arbitrarily ordered, irregular nodes are supported.

Methods on real and complex objects:

- `%value(x,channel,nu)` evaluates one channel at one real coordinate.
- `%values(query,nu)` returns shape `(size(query), n_channels)`; `nu=0`
  evaluates the polynomial and `nu>=1` gives derivatives.
- `%derivative(x,nu,channel)` is equivalent to `%value(x,channel,nu)` with
  `nu=1` by default.
- `%set_values(y,status)` updates all channels without changing nodes.
- `%add_nodes(xnew,ynew,status)` appends a batch of distinct nodes and values,
  rebuilding weights. Invalid updates are rejected without mutating the object.
- `%nodes()`, `%weights()`, `%channels()`, `%valid()` provide introspection.

An invalid constructor returns an object with `%valid() == .false.` and
`status = barycentric_invalid_input`. Node and ordinate values must be
finite, node coordinates must be distinct, and insertion must include the
same channel count. A negative derivative order currently returns zeros;
callers should supply nonnegative orders. Derivatives at or above the
polynomial degree (`nu>=size(nodes)`) return exact zeros.

## Algorithms and limitations

The weights use a log-product computation with a common normalization, to
avoid underflow or overflow from products of knot distances. Evaluation uses
scaled barycentric quotients and returns exact stored ordinates when the
query coincides with a node. Derivative values are formed using the
barycentric spectral differentiation matrix at nodes and interpolated by
the same barycentric formula. The cost to build weights is `O(n^2)`, to
compute the order-`k` node derivatives is `O(k*n^2*n_channels)`, and to
interpolate one query is `O(n*n_channels)`.

High-order global polynomial interpolation remains susceptible to the
Runge phenomenon and to ill-conditioning for many or clustered nodes.
The normalized weights may underflow in extreme cases; constructor
validation rejects that configuration. Derivatives, especially of high
order, can magnify rounding errors. This implementation does not currently
accept user-provided barycentric weights or generate SciPy's randomized
`rng` node permutation.

## Independent reference tests

```bash
python tools/interpolate_reference/compare_barycentric_scipy.py
python tools/interpolate_reference/compare_barycentric_scipy.py --seed 8675309
```

The driver compiles both checked and optimized GNU Fortran versions and
compares real/complex scalar/multichannel values and derivatives (orders
0–3), with and without incremental node insertion, against SciPy 1.17.
The native test checks exact polynomial reproduction, transactional node
insertion, invalid nodes, channel outputs and the convenience generic.

**SciPy 1.17 discrepancy**: `BarycentricInterpolator.add_xi(xnew,ynew)` was
observed to produce incorrect barycentric weights after insertion of
**an even number of nodes together**, while rebuilding a fresh SciPy
`BarycentricInterpolator` using the concatenated data returned the
mathematically correct polynomial. The test therefore checks the new
SciFort results against a freshly rebuilt SciPy object for insertion
cases, reports how many SciPy `add_xi` discrepancies it found, and also
checks against exact low-degree generating polynomials. The finding is
specific to the tested SciPy version, not a claim about other versions.
