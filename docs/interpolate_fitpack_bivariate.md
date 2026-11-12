# FITPACK bivariate smoothing and fixed-knot surface fitting (2026-11-01)

This update adds SciPy-compatible **tensor-product bivariate spline** support,
without touching Delaunay interpolation.

For data on a **complete rectangular grid**, SciFort now also supports
`RectBivariateSpline(s=0)`-style exact interpolation and exact rectangular
integration for any bivariate TCK; see
[rectangular-grid bivariate splines](interpolate_rect_bivariate.md).

## Dependency-free TCK evaluation

`src/interpolate_fitpack_bivariate.f90` defines:

- `fitpack_bivariate_tck_type`: knots `tx(:)`, `ty(:)`; flattened coefficients
  `c(:)`; and degrees `kx`, `ky` (1–5).
- `bisplev(x, y, tck, dx, dy, status)`: scalar evaluation and mixed partial
  derivatives.
- `bisplev(x(:), y(:), tck, ...)`: an `size(x) × size(y)` **tensor-product
  grid** (x index first), matching `scipy.interpolate.bisplev`.
- `bisplev_points(x(:),y(:),tck,...)`: **paired** observation queries. The
  arrays must have the same length.

Any valid fitted TCK can also be transformed into a reusable mixed-partial
B-spline with `%partial_derivative(dx,dy)` or refined by exact knot insertion
with `%insert_knot(axis,x,m)`. See
[bivariate spline calculus](interpolate_bivariate_calculus.md).

Coefficients use FITPACK/SciPy ordering `c((ix-1)*ny+iy)` where `ny` is the
number of basis functions along y. All procedures use native Fortran
floating-point values and SciFort's existing 1D B-spline engine. This
module has **no FITPACK or BLAS/LAPACK dependency** and is exported from
`scifort_interpolate`.

SciPy's `bisplev` convention is respected: evaluation outside a knot domain
is **clamped to the nearest domain boundary**, not extrapolated. Mixed
partial-derivative orders are nonnegative and strictly below the respective
spline degrees, as in SciPy's FITPACK evaluator. Invalid scalar evaluations
return NaN and a nonzero status.

Example (a TCK produced by an optional FITPACK fitting call):

```fortran
use scifort_kinds, only: dp
use scifort_interpolate, only: fitpack_bivariate_tck_type, bisplev

type(fitpack_bivariate_tck_type) :: tck
! Populate tck%tx, %ty, %c, %kx, %ky before use.
if (.not. tck%valid()) error stop 'Invalid TCK'
print *, bisplev(0.3_dp,0.6_dp,tck)
print *, bisplev(0.3_dp,0.6_dp,tck,dx=1,dy=1)
```

## Optional FITPACK 2.0.0 fitting

`optional/fitpack_backend/src/scifort_optional_fitpack.f90` adds:

- `bivariate_smoothing_spline(x,y,z,w,k,s,status,residual)`: FITPACK `surfit`
  with automatic knots and a target weighted residual `s` (SciPy's
  `SmoothBivariateSpline` / `bisplrep(task=0)`), with a common degree `k`
  along the two axes, default 3.
- `lsq_bivariate_spline(x,y,z,tx,ty,w,kx,ky,status,residual)`: fixed interior
  knots with independent degrees `kx` and `ky` (SciPy's
  `LSQBivariateSpline` / `bisplrep(task=-1)`). `tx`, `ty` contain *only*
  interior knots, not repeated boundary knots. Default degrees are 3.

Both accept real scattered observations and optional positive FITPACK residual
weights; `residual` is `sum((w*(z - fitted))**2)` at the observation sites.
Return statuses: 0 = success, 1 = invalid input, 2 = upstream fitting error.
No approximate replacement is returned on failure.

The adapter computes correct `surfit` workspace requirements for degrees 1–5,
including asymmetric degrees, and invokes **one fixed-knot least-squares
pass** rather than allowing an upstream convenience method to change the
provided knots. Upstream FITPACK code remains separate and unchanged.

This is **not a full emulation of every SciPy option**: bounding boxes are
currently inferred from the input ranges; automatic-knot fitting only offers
common x/y degrees, and the adapter handles a single real ordinate channel
per fit. Positive-`s` fits can choose different knots/coefficients from
SciPy's bundled FITPACK even when both satisfy the residual target.

The upstream `fitpack_surface` convenience constructor allocates substantial
workspaces and may use much more memory than simple rectangular-grid spline
fitting. Benchmark memory and convergence for large scattered datasets.

## Verified tests

With the **user-supplied FITPACK 2.0.0 source** and GNU Fortran 14.2.0:

- Core: **130/130** randomized SciPy bivariate TCK evaluations per build,
  checked and `-O3`; degrees 1–5, different degrees per axis, arbitrary
  interior knots, out-of-range queries, and mixed derivatives. The largest
  absolute discrepancy is `2.79e-9`, in high-order derivatives; all tested
  values pass the scale-aware SciPy tolerances.
- Genuine optional backend: **100/100** randomized checked/optimized `-O3`
  reference scenarios on the main seed, with 50 adaptive and 50 fixed-knot
  cases in each build. Maximum observed absolute TCK-evaluation difference
  in the optimized build: `3.41e-13`; maximum weighted-residual discrepancy:
  `4.26e-14`. All fixed-knot results agreed with SciPy's
  `LSQBivariateSpline` within the specified floating-point tolerances.
- Independent seed: **40/40** genuine upstream cases in checked and optimized
  builds, with maximum absolute TCK-evaluation difference `5.15e-14`.
- Public SciFort native interpolation suite: **27/27** in checked and
  optimized GNU Fortran builds, including the new bivariate TCK test.
- The optional bivariate native test checks degrees 1–5, residuals, mixed
  partial derivatives and input validation in both builds.
- A separate synthetic FITPACK-interface compile test passes, but this is
  *not* substituted for the genuine upstream numerical validation.

## Reproduction

Without FITPACK (only GNU Fortran, Python, SciPy and NumPy):

```sh
python tools/interpolate_reference/compare_fitpack_bivariate_tck_scipy.py
python tools/interpolate_reference/run_all_native_interpolation.py
```

With the **original** FITPACK source ZIP (offline, no FPM needed):

```sh
python tools/interpolate_reference/check_optional_fitpack_local.py \
  /path/to/fitpack-main.zip --cases 80
```

The script now also compiles, links and runs the new bivariate native test
and SciPy comparison. It retains its original 1D/periodic/parametric tests.

To use the supplied source with FPM locally:

```sh
python tools/interpolate_reference/prepare_optional_fitpack_fpm.py \
  /path/to/fitpack-main.zip
cd optional/fitpack_backend
fpm test
cd ../..
python tools/interpolate_reference/prepare_optional_fitpack_fpm.py --restore
```

The optional manifest still has its existing **1.0.0 Git tag**; the new
bivariate wrappers are verified against **2.0.0 from the attached ZIP**,
not the old remote pin. The FPM build itself and Windows/MinGW-w64 remain
unverified; do not infer those results from the native GNU Fortran tests.

## Licensing and provenance

Only SciFort wrapper and evaluation code was added, released under SciFort's
MIT license. **No FITPACK source was copied or vendored.** The optional library
is [`perazz/fitpack`](https://github.com/perazz/fitpack), by Federico Perini,
modernizing Paul Dierckx's FITPACK. Its BSD-3-Clause license and notice are
preserved in `LICENSES/perazz-fitpack-BSD-3-Clause.txt`; see
`CODE_PROVENANCE.md`.
