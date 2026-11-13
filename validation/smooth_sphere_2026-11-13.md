# Scattered spherical spline validation

Date: 2026-11-13. Compiler: GNU Fortran (GCC) 17.0.0 20260510
(experimental), Windows. References: SciPy 1.15.3 and independent NumPy dense
least squares/SVD.

## Scope and outcomes

Seven focused programs passed at `-O0` and `-O3`, both with `-fcheck=all`,
`-ffpe-trap=invalid,zero,overflow`, and `-Werror=implicit-interface`:

- `test_interpolate_smooth_sphere`: 119 assertions for exact flat-pole cubic
  data, adaptive weighted targets, derivatives, periodic wrapping, paired/grid
  agreement, metadata, data/weight scaling by `1e150` and `1e-150`, zero-data
  interpolation, huge targets, invalid inputs, capacity/workspace guards,
  relative underflow, conflicting observations, and rank-deficient fitting.
- `test_interpolate_sphere`: 126 existing scattered fixed-knot checks, including
  dense SVD behavior also used as a reference in the new test.
- `test_interpolate_rect_sphere`: 79 rectangular interpolation checks.
- `test_interpolate_rect_sphere_smoothing`: 90 rectangular smoothing checks.
- `test_interpolate_bivariate_fitpack`: shared tensor evaluation dependency.
- `test_interpolate_bivariate_calculus`: tensor derivative/calculus integration.
- `test_interpolate_rect_bivariate`: neighboring interpolation API regression.

The implementation adds a separate scattered solver and constructor; existing
grid helpers are reused without numerical changes. The interpolation umbrella
exports the new constructor. Scope covers these interfaces and their tensor
evaluation paths. No shared numerical algorithm, dependency, or build
configuration changed. The full SciFort suite was not run, following the
repository's impact-based testing policy.

## Independent references

The checked-in `check_smooth_sphere.py` passed in both builds. It covers
27 fits: three observation counts (60, 120, 240), three data families (exact
flat-pole cubic, smooth harmonics, noisy harmonics), and three residual
targets. Sites are unsorted and include poles and repeated seam locations;
weights are nonuniform and positive. Each fit is evaluated at every training
pair, both poles, and 25 random query pairs, for values and five partial
derivatives. This gives 27,054 SciPy value/derivative comparisons per build,
plus full knot, physical coefficient, and weighted residual comparisons.

| Comparison | Acceptance bound | Maximum across both builds |
| --- | --- | --- |
| Pointwise `abs(native-reference)/(1+abs(reference))` | `3e-7` | `4.00361e-14` |
| Maximum coefficient error / `(1+max(abs(reference coefficients)))` | `3e-8` | `4.92980e-15` |
| Residual error / `(1+reference residual)` | `3e-9` | `1.56692e-15` |
| Full knot vectors | `3e-13` absolute, identical lengths | Passed |
| Residual recomputed from evaluated observations | `3e-10*(1+fp)` absolute | Passed |
| Target upper bound | `fp <= 1.0011*s` | Passed |

All 15 polynomial-limit fits are also checked independently. Their two basis
functions are `1-h` and `h`, where `h=(theta/pi)**2*(3-2*theta/pi)`. Weighted
NumPy least squares gives endpoint coefficients; predictions agree within
`3e-11` absolute. These fits report dimension two and rank -1 (not estimated).

## Rank-deficient correction

The rank fixture has 40 equatorial observations, uniform longitudes, unit
weights, and values `2+sin(phi)`. Initial knots are latitude `pi/2` and
longitude `pi/2,pi,3*pi/2`. The reduced coefficient space has dimension ten
but the observation design has rank four. The reference checker explicitly
constructs the constrained tensor map using cubic B-splines and periodic
cosine/sine collocation, then solves its observation matrix with NumPy SVD.
It does not use FITPACK's Givens or minimum-norm implementation.

SciFort returns a valid fit, status 2, rank four, dimension ten, and weighted
residual `0.003103635964350815`. All 35 physical coefficients and 258 values/
derivatives (including predictions away from the equator) agree with the
dense minimum-norm solution. Maximum absolute value/derivative error is
`2.66454e-15` in both builds, against a `3e-10` bound. The native dense SVD
constructor supplies an additional minimum-norm/residual agreement check.

The inspected modern FITPACK source clears the transformed observations
before the first deficient solve. SciFort preserves them, as the inspected
legacy routine does, and corrects read-before-write dummy intents in `fprank`.
The accumulated residual correction was also inaccurate on this fixture:
approximately `0.00156929`, while direct evaluation gives `0.003103635964350815`.
SciFort recomputes the actual objective before accepting or refining a
rank-deficient knot set. This correction matters to the choice between a
least-squares fit and further smoothing.

This fixture does not establish general rank-deficient SciPy parity. Rank
thresholds, coordinate scaling, and reduced-coordinate norms can affect weakly
determined fits. SciPy 1.15.3's high-level class rejects its rank-warning codes;
SciFort retains valid deficient fits with an explicit warning status.

## Upstream fixture

The upstream `mnsphe` example in
[`test/fitpack_tests.f90`](https://github.com/perazz/fitpack/blob/1907c7bad9f5c4557970e906bdd4909eb52b2d4c/test/fitpack_tests.f90)
uses 192 scattered coordinate pairs from
[`test/fitpack_test_data.f90`](https://github.com/perazz/fitpack/blob/1907c7bad9f5c4557970e906bdd4909eb52b2d4c/test/fitpack_test_data.f90).
A temporary checker reads the exact coordinates, converts degrees to radians,
and evaluates the upstream `testsp` formula, preserving its single-precision
literal 0.2 when promoted to binary64. Fresh native fits use its knot
capacities 15/19, `eps=1e-6`, and smoothing targets 500, 135, and 15. The
upstream continuation and fixed-knot modes are outside this constructor.

All three fits pass at `-O0` and `-O3`, returning the polynomial, rank-14, and
rank-46 fits respectively. Their actual residuals meet the requested budgets.
Full knot arrays match SciPy's fits, and 1,458 values/derivatives on the
upstream 9-by-9 query grid agree with maximum normalized error below `5e-14`.
No upstream source or fixture is redistributed. The full upstream package
test suite was not run.

## Commands and documentation checks

Changed modules were compiled ahead of the existing FPM archive for unchanged
modules. Directories were created before compilation. Debug commands below
use CMD syntax:

```cmd
gfortran -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -Jbuild/smooth-sphere-debug -Ibuild/smooth-sphere-debug -Ibuild/gfortran_2806A042DAA8BBF7 -c src/interpolate_sphere_grid_core.f90 -o build/smooth-sphere-debug/interpolate_sphere_grid_core.o
gfortran -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -Jbuild/smooth-sphere-debug -Ibuild/smooth-sphere-debug -Ibuild/gfortran_2806A042DAA8BBF7 -c src/interpolate_sphere_scattered_core.f90 -o build/smooth-sphere-debug/interpolate_sphere_scattered_core.o
gfortran -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -Jbuild/smooth-sphere-debug -Ibuild/smooth-sphere-debug -Ibuild/gfortran_2806A042DAA8BBF7 -c src/interpolate_sphere.f90 -o build/smooth-sphere-debug/interpolate_sphere.o
gfortran -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -Jbuild/smooth-sphere-debug -Ibuild/smooth-sphere-debug -Ibuild/gfortran_2806A042DAA8BBF7 -c src/interpolate.f90 -o build/smooth-sphere-debug/interpolate.o
gfortran -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -Ibuild/smooth-sphere-debug -Ibuild/gfortran_2806A042DAA8BBF7 test/test_interpolate_smooth_sphere.f90 build/smooth-sphere-debug/interpolate.o build/smooth-sphere-debug/interpolate_sphere.o build/smooth-sphere-debug/interpolate_sphere_scattered_core.o build/smooth-sphere-debug/interpolate_sphere_grid_core.o build/gfortran_3C159EB802C22FB9/scifort/libscifort.a -o build/smooth-sphere-debug/test_interpolate_smooth_sphere.exe
build\smooth-sphere-debug\test_interpolate_smooth_sphere.exe
```

The same link/run command was used for the other six programs. Optimized
commands used `-O3` in place of `-O0 -g -fbacktrace` and directory
`build/smooth-sphere-release`, retaining runtime checks and traps. The reference
driver substituted `validation/smooth_sphere_driver.f90` as source and linked
the spherical objects and archive, producing `reference.exe` in each directory.
It accepts batches of cases through programmatic stdin; FPM tests require no
keyboard input.

```cmd
python validation/check_smooth_sphere.py build/smooth-sphere-debug/reference.exe
python validation/check_smooth_sphere.py build/smooth-sphere-release/reference.exe
```

The [guide example](../docs/interpolate_smooth_sphere.md) compiled and ran with
debug flags/objects. It returned value `0.79714879478787759`, residual
`0.099936704223794789`, rank ten, and status zero. Source ASCII, new/changed
line lengths, dummy comments, affected local links, diff whitespace, and
retained license terms were checked. The upstream source/license hashes are
recorded in [provenance](../CODE_PROVENANCE.md); retained license text agrees
with the inspected license after whitespace normalization.

## Remaining limitations

The 2,000,000-entry resource guard, default knot capacities, diagonal rank
test, polynomial rank sentinel, fixed iteration/tolerance settings, target/
weight underflow guards, and conditioning limits are documented in the guide.
Interpolation at `s=0` is not guaranteed for arbitrary observations and default
capacities. C ABI expansion remains deferred.
