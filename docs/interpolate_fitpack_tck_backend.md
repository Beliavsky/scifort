# FITPACK TCK compatibility and optional adaptive-knot smoothing

SciFort provides a new **dependency-free** `scifort_interpolate` API for
FITPACK-style spline tuples `(t,c,k)` and scalar/parametric automatic-knot smoothing,
plus a separate **optional** backend for additional FITPACK fitting interfaces maintained by Federico Perini
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
  [rectangular spherical splines](interpolate_rect_sphere.md) also adapt
  the banded grid solver and adaptive smoothing. Native
  [scattered spherical smoothing](interpolate_smooth_sphere.md) adapts the
  scattered solver and rank reduction. These add no mandatory dependency.
- Native scalar `splrep(s>0)` also adapts the pinned upstream `fpcurf`
  automatic-knot and smoothing solver; the optional package is not required.
- Native nonperiodic `splprep(s>0)` adapts `fppara` and its vector Givens
  helpers from the same licensed source.
- Native periodic `splrep(per=.true.)` adapts `fpperi`, its interpolation
  initializer, wrapped triangular solve, and two-part Givens helper.
- Native periodic `splprep(per=.true.)` adapts `fpclos`, its interpolation
  initializer, and the vector two-part Givens helper, reusing those factors.
- Source references: <https://github.com/perazz/fitpack/blob/main/LICENSE.md>,
  <https://github.com/perazz/fitpack/blob/main/fpm.toml>,
  <https://perazz.github.io/fitpack/tutorial_curve.html>.

## Dependency-free `splrep` / `splev` / `splint` / `splder` / `splantider` / `spalde`

The main public module `scifort_interpolate` exports `fitpack_tck_type` and:

| Routine | Supported behavior |
|---|---|
| `splrep(x,y,k=...,s=0)` | Exact interpolation with FITPACK-compatible padded tuple |
| `splrep(x,y,w=...,k=...,s=positive)` | Native automatic-knot scalar smoothing |
| `splrep(x,y,w=...,k=...,s=...,per=.true.)` | Native periodic interpolation or smoothing |
| `splrep(x,y,task=-1,t=interior,w=...,bbox=...)` | Weighted least-squares fit with user-specified interior knots |
| `splev(x,tck,der=...,ext=...)` | Scalar/vector evaluation and derivatives; `ext=0,1,2,3` |
| `splint(a,b,tck)` | Analytic integral clipped to spline domain, like FITPACK |
| `splder(tck,n=...)`, `splantider(tck,n=...)` | Derivative/antiderivative splines |
| `spalde(x,tck)` | All derivative orders at an interior coordinate |

The `fitpack_tck_type` holds public allocatable `t(:)`, `c(:)` and integer `k`.
For `splrep`, `c` is padded with `k+1` trailing zeros, as in SciPy's FITPACK
`(t,c,k)` representation. The routine accepts degrees 1–5 and validates
input. FITPACK residual weights `w` enter as `sum((w*(y-f(x)))**2)`.

For nonperiodic `task=0, s>0`, `splrep` selects knots and solves FITPACK's derivative-jump
smoothing problem with weighted residual `sum((w*(y-fit))**2)`. The residual
target has relative tolerance `0.001`; a sufficiently large target returns
the weighted least-squares polynomial of degree `k`. Positive smoothing may
use an enclosing `bbox`; exact interpolation still requires the sample bounds.
The existing default `s=0` is preserved even when weights are supplied, unlike
SciPy's weight-dependent default in its legacy `splrep` interface.

The internal solver normalizes coordinates, observations, and weights, then
returns physical knots and padded coefficients. Numerical or allocation
failure returns `fitpack_tck_singular` (3) and an invalid tuple. This includes
unrepresentable scaled targets, targets below achieved rounding error, and
conservative exponent-range guards: normalized adjacent sample gaps must be
at least `tiny(1.0_dp)**(1.0_dp/(4*(k+1)))`, and weights divided by their maximum
must be at least `tiny(1.0_dp)**0.125`. Normalization must preserve distinct
sites and nonzero samples. These guards do not guarantee good conditioning
for every accepted data set. Working storage is linear in sample count for
the supported degrees; no new dependency or C ABI is added.

**Remaining limitations:** `task=1`, periodic fixed-knot fits, `quiet`, and
the modern `make_splrep` interface are not included. Nonperiodic parametric
smoothing is native through `splprep`, as described below. See the
[native smoothing validation](../validation/splrep_smoothing_2026-11-13.md),
including a pinned upstream regression that fails in the older SciPy 1.15.3
reference but succeeds in the adopted Fortran implementation.

### Native periodic scalar fitting

`splrep(x,y,w=...,k=...,s=...,per=.true.)` fits a periodic scalar curve of
degree 1-5 over `[x(1),x(m)]`, with `m > k`. Omitted `per` is false; omitted
`s` remains zero. `per` follows `status` in the argument list, preserving
existing positional calls. `s=0` interpolates the first `m-1` ordinates;
`s>0` selects knots and smooths with relative residual tolerance `0.001`.
The last ordinate and last weight are ignored, including their finiteness
and positivity checks. All coordinates must be finite and strictly increasing;
the first `m-1` ordinates must be finite and active weights finite and positive.
The period is `x(m)-x(1)`. A supplied `bbox` must equal these endpoints.

The fit matches endpoint values and derivatives through order `k-1`.
Knots extend beyond both endpoints and the active coefficients obey periodic
wrapping; `c` still has `k+1` trailing zeros. Large targets yield the weighted
least-squares constant. Periodicity constrains the fit, but does not change
`splev` extrapolation or `splint` clipping: queries outside the period are
not automatically wrapped. Existing derivative and integral operations work
within the base period.

Periodic `task=-1` with supplied knots returns status 2 and an invalid tuple;
`task=1` remains invalid (status 1). Numerical/allocation failures return 3.
The normalized gap and weight guards above also apply. Additionally,
`maxval(abs(x)) <= huge(1.0_dp)/(4*(k+1))`, the physical period and every
physical knot gap must be at least `tiny(1.0_dp)`, and physical knots must
remain distinct after rescaling. These conservative guards reserve range
for the exterior knots. The returned physical representation is checked
against the residual target; interpolation allows normalized rounding error
of `m*(128*epsilon(1.0_dp))**2` in squared weighted residual.

Exact SciPy knot selection, coefficients, and higher derivatives are not
guaranteed. The pinned modernized FITPACK core can follow a different
adaptive path while satisfying the same smoothing equations and residual
contract. See the [periodic validation record](../validation/splrep_periodic_2026-11-13.md).

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

tck = splrep(x, y, w=w, k=3, s=0.5_dp, status=status)
if (status /= 0) error stop 'smoothing failed'
```

## Dependency-free parametric curves: `splprep` and `splev`

The public `scifort_interpolate` module also exports
`fitpack_parametric_tck_type` and `splprep`. This implements SciPy's
**nonperiodic parametric** `splprep` for exact interpolation (`s=0`),
automatic-knot smoothing (`s>0`), and fixed-knot weighted least squares
(`task=-1`), with degrees 1-5 and an arbitrary positive number of coordinate
channels within storage limits. Interpolation and fixed-knot paths use
original SciFort B-spline code; the positive-s path adapts the licensed
perazz/FITPACK `fppara` solver. Support for more than ten channels is a
**SciFort extension** to the legacy SciPy `splprep` interface.

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

curve = splprep(points, w=w, k=3, s=0.25_dp, status=status)
if (status /= 0) error stop 'parametric smoothing failed'
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
For nonperiodic fits, `ub`/`ue` may enlarge the fixed-knot or positive-s domain; exact
interpolation currently requires the domain endpoints to coincide with
`u(1)` and `u(n)`.

For nonperiodic `task=0,s>0`, one shared knot vector and one smoothing parameter serve all
channels. The total target is `sum_i sum_j (w(i)*residual(i,j))**2 <= s`,
with FITPACK relative target tolerance `0.001`. A sufficiently large target
returns the weighted least-squares polynomial curve. The existing default
`s=0` is retained; SciPy's legacy routine has a different omitted-s default.
Coefficients remain unpadded, with shape `(size(t)-k-1,channels)`.

Positive-s chord lengths subtract coordinates first when the difference is
representable, preserving spacing at large offsets. They scale first when
subtraction would overflow. One common coordinate scale preserves the joint
geometry across channels. The solver normalizes the
parameter domain, coordinates, and weights, then restores physical knots
and coefficients. Parameters are shifted before scaling to preserve spacing
at large offsets; the final tuple's residual is checked in physical coordinates.
The physical parameter span must be representable, and distinct physical
knot intervals must be at least `tiny(1.0_dp)`. It uses the scalar solver's
normalized parameter-gap/relative-weight
exponent guards and returns status 3 with an invalid tuple for numerical or
allocation failure, unrepresentable targets, or unmet residuals. Storage is
linear in sample count and channel count for degrees 1-5. There is no new
linked dependency or C ABI.

The nonperiodic perazz and SciPy smoothing iterations use different initial guesses;
normalization can also change the iteration path. Coefficients and derivatives,
especially higher derivatives of noisy curves, can differ substantially even
when both fits meet the target. Residual tolerance is not a derivative-error
bound. Independent validation checks the weighted residual and penalized
normal equations, as well as agreement with the pinned original core and
evaluation of the native tuple by SciPy. See the
[parametric smoothing validation](../validation/splprep_smoothing_2026-11-13.md).

### Native periodic parametric fitting

`splprep(points,u=...,w=...,k=...,s=...,per=.true.)` fits a closed curve with
one shared knot vector and one total weighted residual over all channels.
Degrees 1-5, `m > k` sites, and any positive channel count within storage
limits are supported. The existing signature, unpadded coefficient layout,
and default `s=0` are preserved. `s=0` interpolates the first `m-1` rows;
positive `s` smooths to relative residual tolerance `0.001`. A large target
returns a constant vector whose components are weighted means. Support for
more than ten channels remains a **SciFort extension** to legacy `splprep`.

The final coordinate row and weight are ignored, even if nonfinite or
nonpositive. The input arrays are not modified. The final point is treated
as a copy of the first when computing automatic chord lengths, so the
closing leg runs from row `m-1` to row 1. Zero-length successive or closing
legs make automatic parameters invalid; supplied strictly increasing `u`
allows repeated coordinate points. With supplied `u`, all `m` parameter
values, including the final one, must be finite and strictly increasing.
The period is `u(m)-u(1)`; automatic parameters span `[0,1]`.

The fitted values and derivatives through order `k-1` match at the seam,
subject to floating-point conditioning. The knots extend beyond both
endpoints, and the channel coefficients wrap periodically. This constrains
the fit; `splev` keeps its existing extrapolation policies and does not
automatically wrap out-of-period queries.

If supplied, `ub` and `ue` must equal the parameter endpoints. This deliberate
validation differs from SciPy's periodic path, which ignores these bounds.
Periodic `task=-1` with supplied knots returns status 2 and an invalid tuple;
`task=1` and malformed input return status 1. Numerical/range/allocation
failures return status 3. The normalized gap and relative-weight guards
described above apply. As with periodic scalar fits, parameter magnitudes
reserve headroom for exterior knots: `maxval(abs(u)) <= huge(1.0_dp)/(4*(k+1))`.
The physical span and all physical knot gaps must be at least `tiny(1.0_dp)`;
rounding must preserve distinct knots. The actual physical representation's
residual is verified before returning success. Interpolation admits a
normalized squared-residual rounding allowance of
`m*channels*(128*epsilon(1.0_dp))**2`.

See the [closed-curve validation](../validation/splprep_closed_2026-11-13.md)
for original-core and SciPy comparisons, independent cyclic normal equations,
and the remaining conditioning limitations. Exact adaptive paths and
high-derivative accuracy are not guaranteed for every accepted problem.

**Not yet supported for parametric `splprep`:** periodic fixed-knot fits,
FITPACK `task=1` continuation, and the modern `make_splprep` constructor.

The independent regression compares `splprep` parameters, knots,
coefficients, fitted points and first through third derivatives against
SciPy's FITPACK implementation. Run:

```bash
python tools/interpolate_reference/compare_fitpack_parametric_scipy.py
python tools/interpolate_reference/compare_fitpack_parametric_scipy.py --seed 8675309
python tools/interpolate_reference/run_all_native_interpolation.py
```

## Optional `perazz/fitpack` backend (automatic knot selection)

Scalar and parametric nonperiodic `splrep(s>0)` / `splprep(s>0)` are native.
Periodic scalar and parametric interpolation/smoothing are also native via
`per=.true.`; the optional package is not required for these paths.
The optional backend remains
useful for its separate spline-object API and broader fitting interfaces.

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
