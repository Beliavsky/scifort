# FITPACK TCK compatibility and optional adaptive-knot smoothing

SciFort provides a new **dependency-free** `scifort_interpolate` API for
FITPACK-style spline tuples `(t,c,k)` and a separate **optional** backend for
the adaptive-knot FITPACK smoothing algorithm maintained by Federico Perini
at <https://github.com/perazz/fitpack>.

## Licensing and attribution

- Upstream: Federico Perini, `perazz/fitpack`, a modernization of Paul
  Dierckx's FITPACK, based partly on the FITPACK distribution in SciPy.
- License: BSD-3-Clause; see upstream `LICENSE.md`, reproduced verbatim in
  `LICENSES/perazz-fitpack-BSD-3-Clause.txt` for clarity.
- The copyright notice, license conditions and disclaimer must accompany
  redistributed source/binaries as specified. Do not imply author endorsement.
- The optional FPM package fetches the full upstream implementation separately.
  The native [spherical least-squares interface](interpolate_sphere.md) now
  adapts its fixed-knot periodic/pole constraints under BSD-3-Clause. Native
  [rectangular spherical interpolation](interpolate_rect_sphere.md) also
  adapts the banded grid solver for `s=0`. Neither adds a mandatory dependency.
- Source references: <https://github.com/perazz/fitpack/blob/main/LICENSE.md>,
  <https://github.com/perazz/fitpack/blob/main/fpm.toml>,
  <https://perazz.github.io/fitpack/tutorial_curve.html>.

## Dependency-free `splrep` / `splev` / `splint` / `splder` / `splantider` / `spalde`

The main public module `scifort_interpolate` exports `fitpack_tck_type` and:

| Routine | Supported behavior |
|---|---|
| `splrep(x,y,k=...,s=0)` | Exact interpolation with FITPACK-compatible padded tuple |
| `splrep(x,y,task=-1,t=interior,w=...,bbox=...)` | Weighted least-squares fit with user-specified interior knots |
| `splev(x,tck,der=...,ext=...)` | Scalar/vector evaluation and derivatives; `ext=0,1,2,3` |
| `splint(a,b,tck)` | Analytic integral clipped to spline domain, like FITPACK |
| `splder(tck,n=...)`, `splantider(tck,n=...)` | Derivative/antiderivative splines |
| `spalde(x,tck)` | All derivative orders at an interior coordinate |

The `fitpack_tck_type` holds public allocatable `t(:)`, `c(:)` and integer `k`.
For `splrep`, `c` is padded with `k+1` trailing zeros, as in SciPy's FITPACK
`(t,c,k)` representation. The routine accepts degrees 1–5 and validates
input. FITPACK residual weights `w` enter as `sum((w*(y-f(x)))**2)`.

**Limitations:** `splrep(task=0,s>0)` does *not* implement FITPACK's automatic
knot placement. It returns `fitpack_tck_unsupported_smoothing` (`status=2`)
and an invalid result rather than silently returning a different smoother.
`task=1`, advanced `per` periodic smoothing, and `quiet` are not included.
The automatic-knot case is handled by the optional package described below.

For `splev`, `ext=0` extrapolates endpoint segments, `ext=1` returns zero,
`ext=2` returns status 4 and NaN for an out-of-range query, and `ext=3`
returns the endpoint value. Integral limits outside the knot domain are
clipped, matching SciPy/FITPACK's `splint` behavior.

```fortran
use scifort_kinds, only : dp
use scifort_interpolate, only : fitpack_tck_type, splrep, splev, splint

type(fitpack_tck_type) :: tck
integer :: status

tck = splrep(x, y, k=3, s=0.0_dp, status=status)
if (status /= 0) error stop 'splrep failed'
print *, splev(0.5_dp,tck)
print *, splint(x(1), x(size(x)), tck)
```

## Dependency-free parametric curves: `splprep` and `splev`

The public `scifort_interpolate` module also exports
`fitpack_parametric_tck_type` and `splprep`. This implements SciPy's
**nonperiodic parametric** `splprep` for exact interpolation (`s=0`) and for
fixed-knot weighted least squares (`task=-1`), with degrees 1–5 and an
arbitrary positive number of coordinate channels. This is new SciFort code
using existing B-spline solvers, not copied from FITPACK or SciPy.

```fortran
use scifort_kinds, only : dp
use scifort_interpolate, only : fitpack_parametric_tck_type, splprep, splev

type(fitpack_parametric_tck_type) :: curve
integer :: status
real(dp), allocatable :: values(:,:)

! points(sample, coordinate), unlike scipy.splprep's list of coordinates.
curve = splprep(points, k=3, s=0.0_dp, status=status)
if (status /= 0) error stop 'splprep failed'
values = splev([0.25_dp, 0.50_dp, 0.75_dp], curve)
! values(query,coordinate); scalar splev returns a coordinate vector.
```

`curve%t(:)` is the shared knot vector; `curve%c(:,coordinate)` stores **unpadded**
coefficients; `curve%u(:)` stores the original or automatically selected
parameter values; `curve%k` holds the degree. The automatic parameter is
normalized cumulative Euclidean chord length between successive points.
Call `curve%valid()` to check initialization, or `curve%channels()` to get the
coordinate count. `splev` accepts scalar and vector queries, derivative
orders, and the four FITPACK extrapolation policies.

For `task=-1`, pass **the full knot vector** `t`, including `k+1` repetitions
of each endpoint; this is SciPy `splprep`'s convention and is *different*
from SciPy `splrep(task=-1)`, where `t` contains only interior knots:

```fortran
curve = splprep(points, u=u, w=w, k=3, task=-1, &
                t=[spread(u(1),1,4), interior_knots, &
                   spread(u(size(u)),1,4)], status=status)
```

Weighted fits minimize the sum over observation sites and coordinate
components of `(w(i) * residual(i,j))**2`. For `s=0`, weights do not change
an exact interpolation. Explicit `u` must be finite and strictly increasing.
`ub`/`ue` may enlarge the domain for a fixed-knot least-squares fit; exact
interpolation currently requires the domain endpoints to coincide with
`u(1)` and `u(n)`.

**Not yet supported by the dependency-free core:** automatic knot selection for `s>0`, periodic
`per=.true.` interpolation, or a FITPACK task=1 continuation. Unsupported
features return status `fitpack_tck_unsupported_smoothing` (2), rather than
silently returning a different fit. The optional upstream backend remains
the intended extension for these cases, now that the user-supplied FITPACK 2.0.0 adapter has been compiled and tested offline.

The independent regression compares `splprep` parameters, knots,
coefficients, fitted points and first through third derivatives against
SciPy's FITPACK implementation. Run:

```bash
python tools/interpolate_reference/compare_fitpack_parametric_scipy.py
python tools/interpolate_reference/compare_fitpack_parametric_scipy.py --seed 8675309
python tools/interpolate_reference/run_all_native_interpolation.py
```

## Optional `perazz/fitpack` backend (automatic knot selection)

The `optional/fitpack_backend` directory is a *separate* FPM subproject. Its
manifest pins the **verified Git tag `1.0.0`**. Do not confuse its tag with
its internal metadata version, which has changed across upstream revisions.
Run, where FPM and network access are available:

```bash
cd optional/fitpack_backend
fpm test --compiler gfortran
```

This dependency provides a `fitpack_curve` for SciPy's `UnivariateSpline`-style
adaptive knot selection and constrained weighted residual target `s`.
The adapter provides:

```fortran
use scifort_optional_fitpack, only : univariate_spline
use scifort_interpolate_fitpack_compat, only : fitpack_tck_type,splev

type(fitpack_tck_type) :: curve
integer :: status
real(dp) :: residual

curve = univariate_spline(x,y,s=5.0_dp,k=3, &
                          status=status,residual=residual)
if (status /= 0) error stop 'fit failed'
print *, splev(0.5_dp,curve), residual
```

`s` defaults to `size(x)` and may be zero. Weights must be positive; knots
are selected by the pinned upstream FITPACK fitting procedure. The adapter
returns a valid SciFort TCK tuple so that the existing `splev`, `splint`,
`splder`, `splantider`, and `spalde` routines can be used unchanged.
Only scalar real-valued `y` is currently supported by this adapter.

### Validation status — 2026-11-01

**The user-supplied `fitpack-main.zip` identifies itself as FITPACK 2.0.0.**
This specific source archive was compiled using GNU Fortran 14.2.0, without
FPM or any network download. The original optional manifest still pins
`1.0.0`: the successful *local 2.0.0 build does not verify that Git tag* or
establish the exact Git commit from which the user ZIP was downloaded.

- All **14 core** and **26 C-API** FITPACK source modules compile.
- FITPACK's complete upstream Fortran/C++ test executable passed:
  **66 tests passed, 0 failed**, with checked compilation.
- The genuine `scifort_optional_fitpack` adapter compiled and passed
  `test_fitpack_smoothing.f90` in checked and optimized (`-O3`) builds.
- **140/140 randomized** FITPACK/SciFort TCK and weighted-residual reference
  cases passed in each build. Evaluated values match independent
  `scipy.interpolate.BSpline` with exactly the same returned knots and
  coefficients, to at most `1.7e-15`; reported weighted residuals match
  independently evaluated residuals to at most `9.8e-15`.
- FITPACK 2.0.0 and SciPy 1.17.0 choose *different knots* in **24/140**
  randomized cases (particularly small positive `s`); the maximum value
  discrepancy in these different-knot comparisons was approximately `1.1`.
  The residual constraint and the SciFort TCK conversion remain correct.
  Do **not** infer strict fit identity merely from matching smoothing targets.
- The full upstream `-O3` test executable did **not** complete within the
  test execution limit; the optimized adapter native and numerical checks
  **did** pass. The entire SciFort `fpm test` suite and Windows/MinGW-w64
  remain unverified.

Build and verify directly from the original ZIP, without vendoring source
into SciFort:

```bash
python tools/interpolate_reference/check_optional_fitpack_local.py \
    /path/to/fitpack-main.zip
```

The script builds a static `libfitpack.a`, runs the upstream test executable
(including C++ tests), compiles SciFort's optional adapter, and compares its
output against SciPy. Run with `--optimized --no-cpp` to concentrate on the
optimized adapter without rerunning the longer upstream C++ suite. Use
`--work /path/to/output` to retain objects, module files, and test logs.
Requires GNU Fortran, `ar`, C++17/G++, Python 3.11+, NumPy and SciPy (or
`--no-cpp` / `--no-numerical` to skip the respective optional tests).
On Windows, use matching MinGW-w64 `gfortran`, `g++`, and `ar` binaries;
**Windows compilation has not yet been tested here**.

Other reference paths:

- `compare_fitpack_tck_scipy.py` checks SciFort's dependency-free TCK API,
  without using `perazz/fitpack`.
- `check_optional_fitpack_interface.py` is only the older synthetic mock test;
  it does **not** replace the genuine upstream build now supported above.
- `compare_optional_fitpack_scipy.py` accepts `--executable` to test a locally
  compiled genuine FITPACK driver without FPM. It reports different-knot
  SciPy results separately instead of falsely requiring equivalent
  adaptive-knot sequences.

### Optional periodic and parametric adaptive smoothing (FITPACK 2.0.0)

The optional backend now implements `periodic_univariate_spline` and
`parametric_smoothing_spline`, including `per=.true.` for closed curves.
See [`docs/interpolate_fitpack_periodic_parametric.md`](interpolate_fitpack_periodic_parametric.md)
for API details, verified numerical results, FPM offline setup, and important
positive-`s` differences from SciPy. These routines do not extend the
core `splprep` automatically; call the optional module explicitly.
