# Native closed parametric spline validation (2026-11-13)

Native `splprep(per=.true.)` supports closed-curve interpolation and
automatic-knot smoothing for degrees 1-5. The existing signature, default
`s=0`, shared knots, unpadded coefficient array, and support for more than
ten channels are preserved. Periodic fixed knots and continuation remain
unsupported. See the [API contract](../docs/interpolate_fitpack_tck_backend.md#native-periodic-parametric-fitting).

## References and method

Adapted kernels and the retained BSD-3-Clause license are recorded in
[provenance](../CODE_PROVENANCE.md). The upstream comparison independently
links `splprep_closed_upstream_driver.f90` against the unchanged original
perazz/FITPACK core at commit `1907c7bad9f5c4557970e906bdd4909eb52b2d4c`.
The reference driver replaces the final coordinate row with the first,
then calls `clocur`, `splev`, and `splder` directly. Original upstream
sources and compiled binaries are not committed.

`check_splprep_closed.py` uses deterministic seed 20261115 to generate
180 fits spanning degrees 1-5, 1-5 coordinate channels, nonuniform supplied
parameters or automatic closed chord lengths, nonuniform weights, minimal
sample counts, interpolation, positive targets, and constant limits. Final
coordinate rows deliberately differ from the first. Seven additional
fresh-fit cases read coordinates from the pinned upstream `test_closed_fit`
fixture at run time, preserving default-real literal rounding before
promotion to binary64. They cover cubic targets 900, 10, .1, .5, and zero,
and quintic targets .5 and zero. This is not the complete upstream
object-oriented test or its stateful continuation sequence.

For each fit, checks include parameter and knot agreement, coefficients,
fitted samples, and query derivatives through order three. The final
coordinate row and weight do not enter the weighted residual. Queries are
also evaluated independently from the native tuple by SciPy `BSpline`.
Endpoint derivatives through order `k-1` are checked using a roundoff bound
scaled by the derivative-coefficient norm, accommodating cancellation near
zero in fourth derivatives.

An independent cyclic basis folds the final `k` coefficient rows onto their
periodic counterparts. For positive adaptive fits, design and kth-derivative
jump matrices verify the penalized normal equations. One nonnegative penalty
multiplier is estimated jointly across all channels; maximum allowed scaled
stationarity error is `2e-8`. On common knots, squared weighted prediction
distance is bounded by the absolute residual difference plus
`1e-10*max(1,sum((w*x)**2))` roundoff allowance.

## Results

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
Reference: Python SciPy 1.15.3 and the original pinned Fortran core.
Both checked `-O0` and optimized `-O3` builds use:

```text
-g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface
```

- New public-API regression: 134 checks passed.
- Full `fpm test`: passed, exit code 0.
- Existing periodic scalar (115), scalar smoothing (108), and parametric
  smoothing (125) checks pass, along with parametric interpolation, TCK
  calculus, and fixed-knot spline regressions.
- Both references: 187 fits, 124 adaptive, 132,912 compared numbers each.
  All knot sets agree within the `3e-13` knot tolerance.
- Maximum normalized difference from original upstream: `2.63437e-11`.
- Maximum normalized difference from SciPy: `2.00694e-8`.
  Both pass the `3e-8` coefficient/value/derivative comparison limit.
- Maximum adaptive relative residual error: `0.000994893`, below `0.0011`.
- Maximum independent scaled stationarity error: `5.92517e-14`.
- Maximum independent tuple-evaluation normalized error: `3.56311e-12`.

The native regression also covers ignored NaN/infinite closing coordinates
and weights, unchanged input coordinates, correct closing-leg geometry,
invalid active inputs, finite/increasing final parameters, invalid bounds,
unsupported modes, zero-length chords, constant data, scalar/array agreement,
the scalar periodic identity for one channel, 12 channels, extreme balanced
coordinate/weight scales, subnormal weights, large parameter offsets/scales,
and explicit numerical failure for unrepresentable scaled targets.

These finite tests do not guarantee identical adaptive paths or high-order
derivative accuracy for every accepted curve. Residual tolerance does not
bound derivative errors, especially for noisy, tightly fitted high-degree
curves. Only the compiler and references above were tested. Periodic `ub`/`ue`
must equal the parameter endpoints, whereas SciPy's periodic path ignores
those bounds. `splev` retains its existing extrapolation behavior and does
not wrap out-of-period queries automatically.

## Reproduction

Compile the native driver against the current core, compatibility module,
public interpolation module, and SciFort support archive. Separately compile
the original upstream core and link the upstream driver against it. Then run:

```text
python validation/check_splprep_closed.py build/closed-debug/driver.exe --upstream-tests build/splrep-upstream/fitpack_curve_tests.f90
python validation/check_splprep_closed.py build/closed-debug/driver.exe --upstream-tests build/splrep-upstream/fitpack_curve_tests.f90 --upstream-executable build/closed-debug/upstream.exe
python validation/check_splprep_closed.py build/closed-release/driver.exe --upstream-tests build/splrep-upstream/fitpack_curve_tests.f90
python validation/check_splprep_closed.py build/closed-release/driver.exe --upstream-tests build/splrep-upstream/fitpack_curve_tests.f90 --upstream-executable build/closed-debug/upstream.exe
fpm test
```

Omitting `--upstream-tests` runs the 180 generated cases alone. The fixture
file comes from the pinned source, with its hash recorded in provenance.
Standalone builds compile fresh `src/interpolate_curfit_core.f90`,
`src/interpolate_fitpack_compat.f90`, and `src/interpolate.f90`, using the
optimization/debug flags above. The public new regression compiles through
`scifort_interpolate`; other unchanged support modules come from the archive.

The exact native compilation arguments used are reproduced below in CMD
form. Cache directory names identify this local FPM build and may differ
on another checkout. Repeat with `-O3` replacing `-O0` and `closed-release`
replacing `closed-debug` for the optimized build. In a batch file, double
the percent signs on loop variables.

```bat
set "fitflags=-O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -Ibuild/closed-debug -Ibuild/gfortran_2806A042DAA8BBF7 -Jbuild/closed-debug"
set "fitobjects=build/closed-debug/interpolate_curfit_core.o build/closed-debug/interpolate_fitpack_compat.o build/closed-debug/interpolate.o build/gfortran_3C159EB802C22FB9/scifort/libscifort.a"
gfortran %fitflags% -c src/interpolate_curfit_core.f90 -o build/closed-debug/interpolate_curfit_core.o
gfortran %fitflags% -c src/interpolate_fitpack_compat.f90 -o build/closed-debug/interpolate_fitpack_compat.o
gfortran %fitflags% -c src/interpolate.f90 -o build/closed-debug/interpolate.o
for %t in (closed periodic smoothing parametric parametric_smoothing tck lsq) do gfortran %fitflags% test/test_interpolate_fitpack_%t.f90 %fitobjects% -o build/closed-debug/test_interpolate_fitpack_%t.exe
for %t in (closed periodic smoothing parametric parametric_smoothing tck lsq) do build\closed-debug\test_interpolate_fitpack_%t.exe
gfortran %fitflags% validation/splprep_closed_driver.f90 %fitobjects% -o build/closed-debug/driver.exe
gfortran -O0 -g -fcheck=all -Ibuild/splrep-upstream-compiled validation/splprep_closed_upstream_driver.f90 build/splrep-upstream-compiled/core.o -o build/closed-debug/upstream.exe
```
