# Periodic and parametric FITPACK smoothing (optional backend)

SciFort now wraps **genuine `perazz/fitpack`** for scalar periodic splines and
adaptive-knot parametric smoothing in the separate
`optional/fitpack_backend` FPM subproject. The main SciFort library does **not**
depend on FITPACK, BLAS, or LAPACK. Neither FITPACK source nor built libraries
are vendored in this repository.

## Fortran API

```fortran
use scifort_kinds, only: dp
use scifort_interpolate_fitpack_compat, only: &
    fitpack_tck_type, fitpack_parametric_tck_type, splev
use scifort_optional_fitpack, only: &
    periodic_univariate_spline, parametric_smoothing_spline

real(dp) :: x(100), y(100), points(100,2), weights(100)
integer :: status
real(dp) :: residual
! Initialize x, y, points, weights before calling.

type(fitpack_tck_type) :: periodic_curve
type(fitpack_parametric_tck_type) :: smooth_curve

periodic_curve = periodic_univariate_spline( &
    x, y, w=weights, k=3, s=0.2_dp, &
    status=status, residual=residual)
! Check status == 0 before evaluating.

smooth_curve = parametric_smoothing_spline( &
    points, u=x, w=weights, k=3, s=0.2_dp, per=.true., &
    status=status, residual=residual)
```

This is illustrative code: the arrays must be populated, and variable
declarations must precede executable statements.

- `periodic_univariate_spline(x,y,...)` implements the periodic FITPACK scalar
  smoothing algorithm underlying `scipy.interpolate.splrep(per=1)`. The **last
  ordinate is ignored**, matching the periodic FITPACK convention. `s=0`
  interpolates the first `n-1` ordinates.
- `parametric_smoothing_spline(points,...)` is equivalent in purpose to
  `scipy.interpolate.splprep` with positive `s`, using FITPACK's adaptive
  knot-selection algorithm. Coordinates have shape `(npoints, ndim)` where
  `1 <= ndim <= 10`. The TCK includes fitted parameter values `u` and is
  evaluated by `splev`.
- `per=.true.` enables closed parametric FITPACK (`clocur`). It treats the
  final observation as equal to the first during fitting, even when the
  original input values do not match exactly. This reproduces SciPy's
  `splprep(per=1)` final-ordinate convention without modifying the caller's
  data. Both function values and derivatives are periodic at the endpoints.
- `s` defaults to the number of observations, degree `k` defaults to 3 and
  must be between 1 and 5; weights must be strictly positive. The weighted
  objective uses `sum((w_i * residual_i)**2)` over **all coordinates**,
  omitting the last observation for periodic fits. `s=0` requests exact
  interpolation, with normal floating-point error.
- Returned `status`: 0 success, 1 invalid input, 2 upstream fit failure.
  No approximate fallback is silently substituted if FITPACK fails.
- The returned TCK can be evaluated and differentiated with SciFort's
  `splev`; scalar TCK also supports `splint`, `splder`, `splantider` and
  `spalde`. Outside the knot domain, default `splev` extrapolates the end
  polynomial pieces, **not periodically wrapped values**. Restrict queries
  to the fit domain when strict periodic extension is required.

## Offline reproduction from the supplied 2.0.0 ZIP

With GNU Fortran, GNU `ar`, Python 3.11+, NumPy and SciPy installed:

```sh
python tools/interpolate_reference/check_optional_fitpack_local.py \
    /path/to/fitpack-main.zip --cases 120

python tools/interpolate_reference/check_optional_fitpack_local.py \
    /path/to/fitpack-main.zip --cases 120 --optimized --no-cpp \
    --work build/fitpack-optimized
```

The checked command compiles **all 14 upstream core and 26 C-API modules**,
executes the genuine upstream Fortran/C++ test executable, and runs both the
old and new optional SciFort adapter tests plus SciPy-reference comparisons.
The optimized command recompiles the original upstream numerical source and
both SciFort adapters with `-O3` (omitting the slow upstream C++ test suite).
The `--work` argument preserves executables and detailed logs.

To build with FPM from the **same local ZIP**, without downloading a Git
revision or manually editing the manifest:

```sh
python tools/interpolate_reference/prepare_optional_fitpack_fpm.py \
    /path/to/fitpack-main.zip
cd optional/fitpack_backend
fpm test --compiler gfortran
cd ../..
python tools/interpolate_reference/prepare_optional_fitpack_fpm.py --restore
```

The helper extracts the supplied FITPACK 2.x source into an ignored local
folder, saves the optional subproject's original manifest, and temporarily
uses an FPM `path` dependency. `--restore` returns to the original pinned
manifest and removes only the generated local folder. **The FPM command itself
has not been tested here**, because FPM is not installed. The current remote
manifest still points to FITPACK tag `1.0.0`; these new features were verified
against the supplied **2.0.0 `main` source**. A Git tag or source commit for
the supplied ZIP was not independently verified.

## Numerical validation and SciPy differences

Using GNU Fortran 14.2.0, FITPACK 2.0.0 as supplied, and SciPy 1.17.0:

- Genuine upstream checked Fortran/C++ tests: **66 passed, 0 failed**.
- Both checked and optimized SciFort native adapter tests: **passed**.
- The new randomized reference: **120 cases per seed and build**, with 40
  periodic scalar, 40 open parametric, and 40 closed parametric fits.
  **Two seeds and two builds: 480 cases total.**
- Independent `scipy.interpolate.BSpline` reconstruction from FITPACK's
  returned knots/coefficients: maximum absolute discrepancy across those
  runs **1.46e-11** (including exterior extrapolation).
- Independent first-derivative reconstruction: maximum **6.99e-10**.
- Independent weighted-residual verification: maximum **6.67e-16**.
- SciPy fits with matching knots and **`s=0`** agree to better than
  **8e-14** in evaluated values.

**Important limitation:** For `s>0`, FITPACK 2.0.0 and SciPy can produce
*materially different smoothing curves even when they choose identical knot
locations*. In one second-seed case, the maximum discrepancy within the fit
domain was **0.037**; differences in polynomial **extrapolation** beyond the
fitted domain reached **117**, despite each curve meeting its smoothing
criterion. This is **not** a claim of numerical identity with SciPy's
positive-`s` smoothing solution. The adapter accurately represents the
provided FITPACK fit, but users requiring identical SciPy results should
compare their data and parameters explicitly.

There are also cases where the adaptive knots differ. The comparison script
counts these separately rather than treating different fits as a successful
pointwise match. The upstream FITPACK algorithm controls the positive-`s`
solution; SciFort does not rewrite its smoothing or knot selection.

This integration is limited to real-valued coordinates and ordinates.
Windows/MinGW-w64 and the full SciFort FPM suite are not yet verified.

## Copyright and attribution

This optional adapter is original MIT-licensed SciFort code. It *uses* the
BSD-3-Clause `perazz/fitpack` library by Federico Perini, based on Paul
Dierckx's FITPACK algorithms. See the existing
`LICENSES/perazz-fitpack-BSD-3-Clause.txt` and `CODE_PROVENANCE.md`. If you
redistribute upstream FITPACK, retain its original copyright and BSD notices.
