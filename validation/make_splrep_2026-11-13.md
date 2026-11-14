# Native scalar smoothing constructor validation (2026-11-13)

`make_splrep` returns a native `bspline_interpolator_type`. Its supported
scope is interpolation, automatic smoothing, and supplied full clamped
knots with simple interiors, degrees 1-5. The [constructor guide](../docs/interpolate_make_splrep.md)
documents hard capacities, boundary restrictions, and compatibility limits.

## Compiler and checks

GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows. Checked `-O0`
and optimized `-O3` builds use these additional flags:

```text
-g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface
```

- New public-API regression: 117 checks passed in each build.
- Full `fpm test`: passed, exit code 0.
- Existing closed parametric (134), periodic scalar (115), scalar smoothing
  (108), and parametric smoothing (125) checks pass, along with parametric
  interpolation, TCK calculus, and fixed-knot spline regressions.
- Independent Python reference: SciPy 1.15.3, seed 20261116.
- 165 cases: 40 interpolations, 80 automatic fits (with omitted or full
  explicit capacity), 40 supplied-knot fits, and five nonminimal supplied-knot
  polynomial limits. There are 85 active nonpolynomial smoothing constraints.
- 39,080 strict coefficient/value/derivative/integral comparisons, using
  modern SciPy for interpolation and supplied knots and legacy SciPy FITPACK
  for automatic fits. Maximum normalized error: `1.28873e-11`, below `3e-8`.
- Maximum adaptive relative residual error: `0.000967687`, below `0.0011`.
- Maximum independent scaled stationarity error: `1.1957e-14`, below `2e-8`.

The checker forms an independent B-spline design matrix and kth-derivative
jump matrix using SciPy `BSpline`. It verifies the penalized normal equations
with a nonnegative penalty multiplier. Polynomial limits instead require
negligible derivative jumps and a feasible residual. The five nonminimal
polynomial representations are independently constructed by fitting the
polynomial with SciPy and inserting the supplied interior knots, without
using the new native fixed-knot solver. Native tuple evaluations are also
independently checked against SciPy.

Native tests cover default interpolation, weights and bounds, unpadded object
coefficients, supplied-knot smoothing versus unpenalized LSQ, retention of
supplied knots in the polynomial limit, full/minimum/insufficient capacities,
invalid and nonfinite inputs, mutually exclusive options, infeasible fixed
targets, balanced scales, large coordinate offsets, and unrepresentable targets.
They exercise the constructor through the public `scifort_interpolate` module.

## Modern automatic-knot compatibility

The native constructor deliberately reuses the existing licensed FITPACK
automatic-selection path. Modern SciPy uses a different generator path.
In 16 of the 80 automatic cases, modern SciPy selects different knots.
The maximum normalized difference among the modern value/derivative/integral
diagnostics is `119.576`; these are diagnostics, not tight comparisons.
The native fits independently meet the residual and smoothing equations.
Exact modern automatic knot, coefficient, and derivative parity is not claimed.
Residual tolerance provides no general high-derivative error bound.

SciFort's `nest` is a hard limit and returns status 3 if a residual target
cannot be met. SciPy's limit can be exceeded. Only full clamped supplied
knots with simple interiors are supported. Interpolation rejects supplied
weights/knots/capacity, matching SciPy 1.15.3, but additionally validates
that interpolation bounds coincide with the sample endpoints. Periodic
constructor boundaries are explicitly unsupported; native periodic legacy
interfaces remain available. Only the compiler and reference above were tested.

## Reproduction

The exact standalone compilation arguments are shown in CMD form below.
Cache names identify this local FPM build and can differ on another checkout.
Create `build/modern-debug` first. For the optimized build, replace `-O0`
with `-O3` and `modern-debug` with `modern-release`. Double loop-variable
percent signs when placing these commands in a batch file.

```bat
set "fitflags=-O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -Ibuild/modern-debug -Ibuild/gfortran_2806A042DAA8BBF7 -Jbuild/modern-debug"
set "fitobjects=build/modern-debug/interpolate_curfit_core.o build/modern-debug/interpolate_fitpack_compat.o build/modern-debug/interpolate.o build/gfortran_3C159EB802C22FB9/scifort/libscifort.a"
gfortran %fitflags% -c src/interpolate_curfit_core.f90 -o build/modern-debug/interpolate_curfit_core.o
gfortran %fitflags% -c src/interpolate_fitpack_compat.f90 -o build/modern-debug/interpolate_fitpack_compat.o
gfortran %fitflags% -c src/interpolate.f90 -o build/modern-debug/interpolate.o
gfortran %fitflags% test/test_interpolate_make_splrep.f90 %fitobjects% -o build/modern-debug/test_interpolate_make_splrep.exe
build\modern-debug\test_interpolate_make_splrep.exe
for %t in (closed periodic smoothing parametric parametric_smoothing tck lsq) do gfortran %fitflags% test/test_interpolate_fitpack_%t.f90 %fitobjects% -o build/modern-debug/test_interpolate_fitpack_%t.exe
for %t in (closed periodic smoothing parametric parametric_smoothing tck lsq) do build\modern-debug\test_interpolate_fitpack_%t.exe
gfortran %fitflags% validation/make_splrep_driver.f90 %fitobjects% -o build/modern-debug/driver.exe
python validation/check_make_splrep.py build/modern-debug/driver.exe
fpm test
```

The native build compiles fresh core, compatibility, and public interpolation
modules, linking the existing archive for unchanged support code. No new
dependency, downloaded source, or generated binary is committed. The
[provenance record](../CODE_PROVENANCE.md) identifies the retained BSD-licensed
kernel and the original wrapper/control-flow extensions.
