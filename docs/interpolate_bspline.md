# SciPy-style B-splines (`scifort_interpolate`)

**Date:** 2026-10-31. The independent MIT-licensed implementation in
`src/interpolate_bspline.f90` is re-exported by `scifort_interpolate`.
No new FPM, LAPACK, BLAS, or Python runtime dependencies are required.

## Public API

```fortran
use scifort_kinds, only: dp
use scifort_interpolate, only: bspline, make_interp_spline, &
    bspline_interpolator_type, bspline_array_interpolator_type, &
    bspline_design_matrix, bspline_status_success

type(bspline_interpolator_type) :: s, ds, as, inserted
real(dp) :: x(6), y(6), query(3)
integer :: code

x = [0.0_dp,0.2_dp,0.6_dp,1.2_dp,2.0_dp,3.0_dp]
y = sin(x)
query = [0.1_dp,0.7_dp,1.9_dp]
s = make_interp_spline(x,y,k=3,status=code)
if (code /= bspline_status_success) error stop 'invalid spline'
print *, s%evaluate(query)
print *, s%derivative(query,nu=2)
print *, s%integrate(0.0_dp,2.0_dp)
ds = s%derivative_spline()
as = s%antiderivative()
inserted = s%insert_knot(0.9_dp,m=2)
print *, s%coefficients()
print *, s%knots()
```

### Constructors

- `bspline(t,c,k,status,extrapolate,periodic)` takes an explicit, finite,
  **nondecreasing** knot vector, coefficients and *polynomial degree* `k`.
  It supports degrees `k >= 0`, repeated interior knots, and extra trailing
  coefficient elements (ignored, matching SciPy). `periodic=.true.` wraps
  evaluation coordinates modulo the base span but does **not** solve a
  periodic interpolation problem.
- `make_interp_spline(x,y,k,bc_type,status,extrapolate,knots)` interpolates
  strictly increasing real coordinates. Default `k=3`, `bc_type='not-a-knot'`.
  The default knot convention matches SciPy for degrees 0–15, including
  half-integer positions for even orders. Optional user-supplied `knots`
  must satisfy the collocation and degree/constraint count. `natural` and
  `clamped` endpoint conditions are supported for **cubic splines only**.
  `periodic` boundary conditions are now supported for degrees 2–15 (and
  SciPy-compatible degree-0/1 aliases); see
  [periodic B-splines](interpolate_periodic_bspline.md). Arbitrary
  derivative-constraint tuples remain unimplemented.
- Both constructors are overloaded for `c(:)` / `y(:)` and
  **matrix-valued** `c(:,:)` / `y(:,:)`, where the first dimension is the
  coefficient/sample axis and the second dimension denotes output channels,
  corresponding to SciPy's `axis=0` convention. Their returned derived types
  are `bspline_interpolator_type` and `bspline_array_interpolator_type`.
- The interpolation constructor uses a **dense pivoted collocation solver**.
  It is intended for moderate 1D interpolation systems; it is **not** an
  efficient large-n or multidimensional tensor-product spline backend.

### Methods and numerical behavior

Scalar splines provide `%evaluate(x)` and `%derivative(x,nu)`, each with scalar
and vector overloads. `%derivative_spline(nu=1)` and
`%antiderivative(nu=1)` return spline representations. `%integrate(a,b)`
evaluates the *exact analytic* piecewise-polynomial definite integral using
an antiderivative coefficient recurrence, with periodic wrapping if requested.

`%insert_knot(x,m=1)` uses Boehm's shape-preserving knot-insertion recurrence;
the returned object has the same function values. For nonclamped splines,
SciPy's FITPACK-based insertion may also change exterior knots; SciFort
preserves those original exterior knots rather than reproducing that
implementation detail. Knot insertion at a knot already having multiplicity
`k+1` is rejected. No source objects are mutated by these operations.

`%coefficients()`, `%knots()`, `%degree()`, and `%valid()` provide introspection.
Multichannel types additionally offer `%channels()` and channel-shaped scalar
and batch evaluation, derivatives and integrals, coefficient extraction and
knot insertion.

`bspline_design_matrix(x,t,k,status,extrapolate)` returns a **dense** real
matrix (rows: queries; columns: basis functions), unlike SciPy's CSR sparse
array. With the default `extrapolate=.false.`, out-of-range queries return
`bspline_status_invalid_input` and a 0x0 array instead of raising an exception.
Use `matmul(matrix,coefficients)` to reproduce `BSpline` evaluation.

Invalid constructors return an uninitialized object with `%valid()==.false.`
and a status of 1. A numerically singular collocation system returns status 2.
Out-of-range scalar evaluation returns NaN when `extrapolate=.false.`; integral
bounds are clipped to the base interval in that mode. The default behavior is
endpoint-polynomial extrapolation. Derivative spline construction reports an
invalid-input status for a knot multiplicity that makes a derivative object
nonrepresentable, in line with SciPy.

## Current limitations

- **Real-valued** data, rank-one or rank-two `axis=0` ordinates only;
  complex data and arbitrary ordinate axes are not implemented yet.
- `make_interp_spline` lacks periodic boundary-condition fitting and custom
  derivative lists, and its dense solver is not suitable for very large grids.
- The design matrix is dense rather than sparse; this standalone B-spline
  module does not export an explicit multidimensional knot/coefficient object.
  **Tensor-product** 1–6D interpolation of sampled regular-grid data (degrees
  1, 3, 5) is now available through `scifort_interpolate_regular_grid`.
- SciPy 1.17 has a degree-zero antiderivative/integrate endpoint behavior at
  the final knot that does not match the mathematically continuous integral.
  SciFort preserves the continuous integral convention; the independent
  comparisons exclude this single compatibility exception.
- Fortran status-code conventions replace Python exceptions. Exact
  floating-point equivalence is not expected for ill-conditioned high-degree
  extrapolation; compare with relative/scale-aware tolerances.

## Tests

```sh
python tools/interpolate_reference/compare_bspline_scipy.py
python tools/interpolate_reference/compare_bspline_scipy.py --seed 8123
```

The verifier compiles native Fortran in both checked and optimized modes and
compares 1,390 randomized/deterministic cases per seed and compiler build,
covering more than 1.07 million values per build. It checks construction
coefficients, evaluation, derivatives, exact definite
integrals, derivative/antiderivative objects, knot insertion and dense basis
matrices with SciPy 1.17. Its test corpus includes degrees 0–7, nonuniform
and repeated knots, nonclamped knots, real multichannel outputs, periodic
extrapolation, and mixed endpoint conditions. The independent SciPy reference
requires NumPy and SciPy only for testing. Native tests are in
`test/test_interpolate_bspline.f90`.

## Upstream attribution and reuse

[Jacob Williams' bspline-fortran](https://github.com/jacobwilliams/bspline-fortran)
provides BSD-style licensed tensor-product splines in 1–6D with both
subroutine and object interfaces, and includes its own basic BLAS routines.
Its APIs primarily interpolate sampled grid data, with spline **order**
(`degree + 1`) as an argument, whereas SciPy's `BSpline` uses explicit knot,
coefficient, and **degree** arrays. This new module was written independently
from the standard de Boor, Cox–de Boor, and Boehm mathematical recurrences;
**no upstream source was copied or added as a dependency**. For future
higher-dimensional spline support, evaluate a revision-pinned FPM dependency
(e.g. upstream's documented `tag='7.3.0'`) and retain its full BSD notices if
its code is ever copied. See [`interpolation_upstream_reuse.md`](interpolation_upstream_reuse.md)
and [`CODE_PROVENANCE.md`](../CODE_PROVENANCE.md).

## Periodic construction (2026-11-01)

For `make_interp_spline(x,y,k=...,bc_type='periodic')`, SciFort now fits
cyclic B-splines with matching endpoint derivatives, using shared
collocation for multichannel data. See
[periodic B-spline construction and validation](interpolate_periodic_bspline.md)
for examples, input rules and limitations.
