# SciFort rectilinear-grid interpolation (2026-10-31)

The native module `scifort_interpolate_regular_grid` is re-exported from
`scifort_interpolate`. It implements nearest, multilinear, recursive monotone
PCHIP, and degree 1, 3, and 5 **tensor-product interpolating B-splines**
compatible with the `'pchip'`, `'slinear'`, `'cubic'`, and `'quintic'` methods of
[`scipy.interpolate.RegularGridInterpolator`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.RegularGridInterpolator.html).
An `interpn`-style one-shot convenience function uses these methods as well.
Real and complex data and one through six independent coordinate dimensions
are supported without introducing an external numerical dependency.

This is an original implementation based on standard tensor-product linear
interpolation and cell location. SciPy is used for independent numerical tests,
not as an implementation source. The independent Fortran package
[`finterp`](https://github.com/jacobwilliams/finterp) is a promising optional
back-end for a later increment; its BSD license permits reuse with attribution.
**No `finterp` or other upstream interpolation source has been copied here.**

## Basic use

```fortran
program example
    use scifort_kinds, only : dp
    use scifort_interpolate, only : grid_axis, regular_grid_interpolator, &
        regular_grid_interpolator_type, interpn
    implicit none
    type(grid_axis) :: axes(2)
    type(regular_grid_interpolator_type) :: rgi
    real(dp) :: values(3,2), queries(2,2)

    axes(1) = grid_axis([0.0_dp,1.0_dp,3.0_dp])
    axes(2) = grid_axis([5.0_dp,2.0_dp]) ! descending is supported
    values(:,1) = [3.0_dp,5.0_dp,9.0_dp]
    values(:,2) = [-3.0_dp,-1.0_dp,3.0_dp]

    ! Function f(x,y) = 2*x + 2*y - 7, exactly multilinear.
    rgi = regular_grid_interpolator(axes, values)
    print *, rgi%evaluate([1.5_dp,3.5_dp])
    queries(1,:) = [0.25_dp,3.0_dp]
    queries(2,:) = [2.5_dp,4.0_dp]
    print *, rgi%evaluate_values(queries)

    ! The one-shot method returns an (nqueries,1) array for scalar data.
    print *, interpn(axes, values, queries)
end program example
```

The above values are an affine function `2*x + 2*y - 7`; the multilinear
interpolator returns that function exactly up to roundoff. Coordinate arrays
may be unevenly spaced, ascending or descending in each dimension. Singleton
axes are supported for nearest and multilinear interpolation, but not for
B-spline methods. Each axis is passed as a `grid_axis` with an allocatable
`points(:)` member; the convenience structure constructor `grid_axis([..])`
allocates it automatically.

## Interfaces and shapes

| Entry point | Purpose |
|---|---|
| `regular_grid_interpolator(axes, values, ...)` | Construct reusable real-valued grid object |
| `complex_regular_grid_interpolator(axes, values, ...)` | Construct reusable complex-valued grid object |
| `grid%evaluate(x(:), status, method, nu)` | Evaluate one point; scalar output requires one channel |
| `grid%evaluate_channels(x(:), status, method, nu)` | Return all output channels at one point |
| `grid%evaluate(x(:,:), status, method, nu)` | Return `(nqueries,nchannels)` array; **queries are rows** |
| `grid%evaluate_values(x(:,:), status, method, nu)` | Return length-`nqueries` vector for one channel |
| `grid%valid()`, `%dimensions()`, `%channels()` | Interpolator properties and validation |
| `interpn(axes,values,x,...)` | One-shot real-valued constructor/evaluation |
| `complex_interpn(axes,values,x,...)` | One-shot complex-valued constructor/evaluation |

The `values` argument is an assumed-rank array of rank `N`, with leading
axis lengths `(n1,...,nN)`, or rank `N+1` with shape `(n1,...,nN,nchannels)`.
`N` must be 1 through 6 and `nchannels` must be positive. Fortran array
storage is column-major, so the first coordinate varies fastest in the native
flattened storage. The output-channel index is always the final input dimension.
The supported scalar data types are `real(dp)` and `complex(dp)` (separate
constructors); integral inputs should be converted to `real(dp)` explicitly.

Optional constructor parameters:

- `method='linear'` (default), `'nearest'`, `'pchip'`, `'slinear'`,
  `'cubic'` or `'quintic'` (also overrideable per evaluation). PCHIP requires
  at least **4** knots per axis, and the B-spline methods require 2, 4, and
  6 knots respectively, matching SciPy's `RegularGridInterpolator` minimums.
  Complex data are supported for all methods **except PCHIP** (SciPy 1.17
  rejects complex-valued PCHIP).
- `bounds_error=.true.` (default): out-of-domain samples return quiet NaN and
  `grid_status_out_of_bounds=2` through the optional `status` argument.
  SciPy raises `ValueError` for this case; SciFort reports an error code.
- `bounds_error=.false.`: out-of-domain samples use `fill_value` (default NaN).
- `extrapolate=.true.` together with `bounds_error=.false.`: use the first/last
  linear cell outside each axis, or the nearest endpoint for nearest-neighbor
  interpolation. B-splines use the polynomial extension of the boundary
  spline span (this can grow rapidly for quintic splines). Singleton axes
  have a constant extrapolation for linear/nearest methods.

Nearest-neighbor ties select the lower **physical** coordinate, including for
descending axes, as in SciPy. Query coordinates containing NaN yield NaN
rather than raising a bounds error. Complex default fill and NaN query results
are `NaN+0i`; a caller-supplied complex fill value is preserved.

The status constants are `grid_status_success=0`,
`grid_status_invalid_input=1`, `grid_status_out_of_bounds=2`, and
`grid_status_singular=3` (ill-conditioned interpolation system).
A failed constructor returns an invalid object whose `%valid()` is false.

## Tensor-product B-splines and partial derivatives

With `method='slinear'`, `'cubic'` or `'quintic'`, SciFort constructs a
not-a-knot B-spline basis for each coordinate and solves the separable
collocation systems successively using a **banded LU factorization without
pivoting**. The B-spline collocation matrices are totally nonnegative at
valid Schoenberg-Whitney points. Each factor has a half-bandwidth at most the
spline degree `k` (1, 3, or 5). Factorization requires `O(n_d*k²)` operations
and `O(n_d*k)` memory per axis; compared with the earlier dense LU, this
eliminates `O(n_d³)` work and `O(n_d²)` storage. A pivot-threshold check
reports `grid_status_singular` for numerically degenerate grids rather than
silently accepting an unstable solve. Local Cox–de Boor support evaluates
without summing over all grid coefficients. No upstream code or external
BLAS/LAPACK dependency was incorporated.

```fortran
use scifort_kinds, only: dp
use scifort_interpolate, only: grid_axis, regular_grid_interpolator, &
    regular_grid_interpolator_type
implicit none
type(grid_axis) :: axes(2)
type(regular_grid_interpolator_type) :: rgi
real(dp) :: values(4,4), x(4), y(4)
integer :: i, j

x = [-1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp]
y = [3.0_dp, 2.0_dp, 1.0_dp, 0.0_dp]   ! descending coordinate axis
axes(1) = grid_axis(x)
axes(2) = grid_axis(y)
do j = 1, 4
    do i = 1, 4
        values(i,j) = x(i)**3 + x(i)*y(j)**2
    end do
end do
rgi = regular_grid_interpolator(axes, values, method='cubic')
print *, rgi%evaluate([0.2_dp, 1.5_dp])
print *, rgi%evaluate([0.2_dp, 1.5_dp], nu=[1, 0])  ! d/dx
print *, rgi%evaluate([0.2_dp, 1.5_dp], nu=[1, 1])  ! d²/dxdy
```

For a degree-`k` spline, `nu=[r1,...,rN]` supplies nonnegative per-axis
partial derivative orders, including mixed derivatives; any order larger
than `k` gives zero. The optional `nu` keyword works with scalar and batch
`%evaluate`, `%evaluate_channels` and `%evaluate_values`, including complex
interpolators. This is analogous to SciPy's `RegularGridInterpolator.__call__`
`nu` keyword. It is **not** available for `'linear'` or `'nearest'` methods.
At a knot where a high-order derivative is discontinuous, SciFort selects the
same physical right-hand spline segment as SciPy, including descending axes.

B-spline evaluation also supports an explicit `method=` override from a
previously constructed grid; switching methods creates a temporary spline fit.
Batch overrides reuse that fit for all points. For repeat queries, construct
with the final method to avoid extra work.

Run the independent SciPy regression, which uses SciPy's direct sparse solver
for accurate high-degree reference values, with both compiler configurations
and two random seeds:

```sh
python tools/interpolate_reference/compare_tensor_grid_scipy.py
python tools/interpolate_reference/compare_tensor_grid_scipy.py --seed 8675309
```

Both seeds run 189 cases and compare 2,997 output values per compiler build,
covering one to four dimensions, different orientations/spacings, complex and
multichannel data, fills/bounds/extrapolation, and mixed derivatives. The
native `test/test_interpolate_tensor_grid.f90` also verifies a six-dimensional
cubic and quintic interpolation grid, exact polynomials, and invalid inputs.
Strongly extrapolated quintic derivatives can be numerically ill-conditioned;
use scale-aware comparison tolerances in these cases.

## Multidimensional shape-preserving PCHIP

The `method='pchip'` interpolator constructs a fresh 1D Fritsch–Butland
monotone cubic interpolation along each data fiber. It collapses dimensions
**from the last axis to the first**, reproducing the recursive order used by
SciPy's `RegularGridInterpolator`. As PCHIP slopes depend nonlinearly on the
data, a different evaluation order is generally not equivalent. Every axis
needs at least four knots. Nonuniform and descending axes, extrapolation,
multiple output channels and `interpn` are supported for real `dp` data.

```fortran
use scifort_kinds, only: dp
use scifort_interpolate, only: grid_axis, regular_grid_interpolator, &
    regular_grid_interpolator_type
implicit none
type(grid_axis) :: axes(2)
type(regular_grid_interpolator_type) :: rgi
real(dp) :: samples(4,4)
integer :: i, j
axes(1) = grid_axis([0.0_dp,0.2_dp,1.0_dp,3.0_dp])
axes(2) = grid_axis([5.0_dp,3.0_dp,0.0_dp,-2.0_dp])
do j = 1,4
    do i = 1,4
        samples(i,j) = axes(1)%points(i)**2 + axes(2)%points(j)
    end do
end do
rgi = regular_grid_interpolator(axes,samples,method='pchip', &
    bounds_error=.false.,extrapolate=.true.)
print *, rgi%evaluate([0.6_dp,2.0_dp])
```

`nu=` partial derivatives are **not** available for PCHIP (matching SciPy).
PCHIP is evaluated on-demand for each query; unlike a fitted tensor-product
B-spline, its sequential slopes must be updated as the axes are collapsed.
Complex-valued ordinates are rejected for `method='pchip'`.

The independent SciPy comparison exercises 1–6 dimensions, randomized
irregular axes (including reversed order), nonlinear and flat data,
multichannel output, bounds/fill/extrapolation, and NaN queries:

```sh
python tools/interpolate_reference/compare_tensor_pchip_scipy.py
python tools/interpolate_reference/compare_tensor_pchip_scipy.py --seed 8675309
```

Each seed comprises 232 cases and compares 3,528 returned values across
both compiler configurations; the largest observed scale-normalized
error over both seeds was approximately `4.0e-14`. Native coverage resides
in `test/test_interpolate_pchip_grid.f90`.

## Banded collocation performance

The self-contained benchmark
`tools/interpolate_reference/benchmark_tensor_grid_solver.f90` reports
CPU seconds for fitting a one-dimensional quintic spline and maximum
interpolation error at its input knots. It can also be compiled against the
previous archived implementation to compare dense and banded algorithms.
On one Linux GNU Fortran `-O3` build, the following *single-run* CPU times
were observed (not a cross-machine performance guarantee):

| Grid size | Previous dense LU | New banded LU | Maximum new nodal error |
|---|---:|---:|---:|
| 128 | 0.002612 s | 0.000034 s | `4.44e-16` |
| 256 | 0.042759 s | 0.000067 s | `4.44e-16` |
| 512 | 0.769019 s | 0.000104 s | `4.44e-16` |
| 4096 | Not benchmarked | 0.000746 s | `4.44e-16` |

For larger problems, the multi-dimensional coefficient array itself may
still consume considerable memory; the banded solver reduces only the
additional *per-axis* collocation-factor storage.

## Validation and limitations

To run the checked and optimized independent reference comparison:

```sh
python tools/interpolate_reference/compare_regular_grid_scipy.py
python tools/interpolate_reference/compare_regular_grid_scipy.py --seed 8675309
```

The native FPM smoke test is `test/test_interpolate_regular_grid.f90`.
The reference suites cover dimensions 1--6, complex/real samples, one to three
channels, singleton axes, irregular and reversed coordinates, nearest/linear
methods, fills, extrapolation, NaN queries, tie-breaking, bounds errors and
function values (540 randomized scenarios in each seed).

Limitations compared with full `scipy.interpolate.RegularGridInterpolator`:

- The `'slinear_legacy'`, `'cubic_legacy'`, and `'quintic_legacy'`
  variants are not implemented. `nu` applies only to `'slinear'`, `'cubic'`,
  and `'quintic'` methods. PCHIP requires finite real-valued ordinates.
  B-spline methods also reject input containing nonfinite ordinates,
  rather than propagating NaN through collocation.
- The spline solver avoids forming the enormous global tensor-product
  collocation matrix, but spline evaluation still scales with the number of
  locally supported tensor basis functions `(k+1)**N` per query, and total
  stored coefficients scale with the complete grid size. Banded LU avoids
  partial pivoting for the totally nonnegative collocation systems and
  reports near-zero pivots as singular. Method overrides between different
  B-spline degrees incur a refit; construct a reusable object with the
  desired method for large query batches.
- Input values support a **single optional trailing channel dimension**;
  Python-style multiple trailing batch axes must be flattened into channels.
- Arbitrary rank query arrays are not accepted: `x` must be `(nqueries,N)`
  or one rank-one point; a single-channel batch may use `%evaluate_values`.
- The native `status` mechanism replaces Python exceptions; a failing strict
  batch returns status 2 and NaNs only at out-of-range points (other rows are
  still evaluated).
- No system BLAS/LAPACK or new FPM package is needed for this implementation.

The separate BSD-licensed `finterp` package would be valuable for testing and
possibly optimizing the same 1--6D numerical kernel; this increment keeps the
new SciPy-style API self-contained and independently testable, without
vendoring source or changing dependency versions. An FPM dependency on
`bspline-fortran` remains a good candidate for expanding tensor-product
higher-order splines.
