# Rectangular spherical smoothing validation

Date: 2026-11-13. Compiler: GNU Fortran (GCC) 17.0.0 20260510
(experimental), Windows. Numerical references: SciPy 1.15.3, NumPy constrained
least squares, and pinned modern perazz/FITPACK.

## Scope and outcomes

Six focused programs passed at `-O0` and `-O3`, both with bounds/runtime
checks, floating-point exception traps, and implicit-interface errors:

- `test_interpolate_rect_sphere_smoothing`: 90 assertions covering active
  residual targets, the polynomial limit, exact poles, periodicity, mixed
  derivatives, metadata, scaling data by `1e150` and `1e-150` while scaling
  the target by their squares, zero/constant grids, huge targets, scaled-target
  underflow, and numerical failure at an unattainably small target.
- `test_interpolate_rect_sphere`: 79 interpolation/input-validation assertions.
- `test_interpolate_sphere`: 126 scattered spherical regression assertions.
- `test_interpolate_bivariate_fitpack`: shared tensor evaluation dependency.
- `test_interpolate_bivariate_calculus`: shared tensor derivatives/calculus.
- `test_interpolate_rect_bivariate`: neighboring rectangular interpolation API.

The existing interpolation reference checker also passed all 48 grids and
24,576 value/derivative comparisons in each build, with unchanged maximum
normalized error `4.49027e-12`. Its independent missing-pole estimate also
remains correct. The interpolation branch is kept separate from smoothing.

The smoothing checker covers 96 fits: three grids (6-by-8, 9-by-12, 12-by-16),
eight symmetric/asymmetric pole configurations, and four targets (0.01, 0.2,
2, 1000). Inputs include noisy harmonics, shifted longitude periods, nonuniform
sites, estimated/exact/soft pole observations, and flat poles. Every fit is
checked for validity, grid-residual agreement, and the upper residual bound
`fp <= 1.0011*s`. Queries include all input nodes, poles, and random pairs,
for values and five coordinate partial derivatives.

Strict SciPy comparison covers 12 C0 fits. For some C1-pole cases, SciPy 1.15.3
returned nonfinite coefficients/values or different outputs depending on
prior calls. These are not treated as usable correctness references, and no
general C1 numerical parity with that version is claimed. A separately
compiled full modern FITPACK core, with the two corrections below, matches
all 96 native fits in the debug comparison. Together with the 12 SciPy fits,
this gives 92,664 value/derivative reference comparisons, plus full knot,
coefficient, and residual comparisons.

| Comparison | Acceptance bound | Debug maximum, modern FITPACK and SciPy |
| --- | --- | --- |
| Pointwise `abs(native-ref)/(1+abs(ref))` | `3e-7` | `2.12252e-13` |
| Maximum coefficient error / `(1+max(abs(ref coefficients)))` | `3e-8` | `3.05992e-15` |
| Residual error / `(1+reference residual)` | `3e-9` | `4.29286e-15` |
| Full knot arrays | `3e-15` absolute, identical lengths | Passed |

All 24 polynomial-limit cases are independently checked with NumPy least
squares. A clamped cubic latitude B-spline has four coefficients; C1 poles
identify the adjacent endpoint coefficient with the endpoint, exact poles
fix it, and soft pole observations add a unit-weight equation. The grid
observations determine the remaining coefficients. Values and first/second
latitude derivatives agree within `3e-10` absolute; longitude and mixed
derivatives are below `3e-11`. This provides an independent check of pole
semantics without depending on the banded solver.

The upstream `mnspgr` fixture was also checked for its four smoothing cases:
the 11-by-14 `daspgr` data in
[`test/fitpack_test_data.f90`](https://github.com/perazz/fitpack/blob/1907c7bad9f5c4557970e906bdd4909eb52b2d4c/test/fitpack_test_data.f90),
with options from
[`test/fitpack_tests.f90`](https://github.com/perazz/fitpack/blob/1907c7bad9f5c4557970e906bdd4909eb52b2d4c/test/fitpack_tests.f90).
Fresh fits use `s=60`/C0, `s=0.05`/C0, `s=0.05`/C1 with exact poles, and
`s=0.05`/C1 with exact flat poles. Native residuals are respectively
`47.470346641007005`, `0.05000207899497225`, `0.050004520888705044`, and
`0.049992770542566005`. Across 3,696 additional value/derivative comparisons
with the modern reference, maximum normalized error is `6.45920e-14`.
The fixture and source are not redistributed; a temporary checker read the
downloaded data. Warm-start and fixed-knot modes in the upstream test are
outside this public interface. The full upstream package suite was not run.

The full SciFort suite was not run. The changed behavior is confined to the
spherical grid core, its wrapper, and the shared spherical object metadata;
unchanged numerical kernels, dependencies, and build configuration are covered
by the focused regression scope under the repository testing policy.

## Native commands

Changed modules were compiled ahead of the existing FPM archive for unchanged
modules. Debug commands, shown in CMD syntax:

```cmd
gfortran -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -Jbuild/rect-sphere-smoothing-debug -Ibuild/rect-sphere-smoothing-debug -Ibuild/gfortran_2806A042DAA8BBF7 -c src/interpolate_sphere_grid_core.f90 -o build/rect-sphere-smoothing-debug/interpolate_sphere_grid_core.o
gfortran -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -Jbuild/rect-sphere-smoothing-debug -Ibuild/rect-sphere-smoothing-debug -Ibuild/gfortran_2806A042DAA8BBF7 -c src/interpolate_sphere.f90 -o build/rect-sphere-smoothing-debug/interpolate_sphere.o
gfortran -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -Jbuild/rect-sphere-smoothing-debug -Ibuild/rect-sphere-smoothing-debug -Ibuild/gfortran_2806A042DAA8BBF7 -c src/interpolate.f90 -o build/rect-sphere-smoothing-debug/interpolate.o
gfortran -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -Jbuild/rect-sphere-smoothing-debug -Ibuild/rect-sphere-smoothing-debug -Ibuild/gfortran_2806A042DAA8BBF7 test/test_interpolate_rect_sphere_smoothing.f90 build/rect-sphere-smoothing-debug/interpolate.o build/rect-sphere-smoothing-debug/interpolate_sphere.o build/rect-sphere-smoothing-debug/interpolate_sphere_grid_core.o build/gfortran_3C159EB802C22FB9/scifort/libscifort.a -o build/rect-sphere-smoothing-debug/test_interpolate_rect_sphere_smoothing.exe
build\rect-sphere-smoothing-debug\test_interpolate_rect_sphere_smoothing.exe
```

The same link/run command was used for each of the other five test programs.
Optimized commands used `-O3` instead of `-O0 -g -fbacktrace` and directory
`build/rect-sphere-smoothing-release`, retaining checks and traps. Directories
were created before compilation. The reference executable substituted
`validation/rect_sphere_driver.f90` as source, used the corresponding spherical
objects and archive, and wrote `reference.exe` in the build directory.

```cmd
python validation/check_rect_sphere.py build/rect-sphere-smoothing-debug/reference.exe
python validation/check_rect_sphere.py build/rect-sphere-smoothing-release/reference.exe
python validation/check_rect_sphere_smoothing.py build/rect-sphere-debug/reference.exe --upstream build/rect-sphere-upstream/fixed_reference.exe
python validation/check_rect_sphere_smoothing.py build/rect-sphere-smoothing-release/reference.exe
```

The extended debug comparison used the earlier `build/rect-sphere-debug`
objects with identical numerical implementation; subsequent changes only
adjusted comments/indentation. The optimized checker runs all 96 native
invariant checks, 12 SciPy comparisons (10,296 value/derivative references),
and 24 independent polynomial checks. Its maximum normalized errors are
`6.63703e-14` for values/derivatives, `1.13123e-15` for coefficients, and
`2.73884e-15` for residuals.

## Modern FITPACK comparison

The exact inspected source is
[`src/fitpack_core.F90`](https://github.com/perazz/fitpack/blob/1907c7bad9f5c4557970e906bdd4909eb52b2d4c/src/fitpack_core.F90),
commit `1907c7bad9f5c4557970e906bdd4909eb52b2d4c`, SHA-256
`02e413063f7df855c45eed6119130f7ab7a07cfa3e2f6ac0e42b6f496df414ec`.
A temporary copy made only these two numerical corrections, also present
in SciFort:

1. In `fpgrsp`, initialize cosine/sine collocation when either pole requests
   continuity, rather than only the north pole.
2. In `new_knot_dimension`, alternate axis choices on equal proposed knot
   counts, following the spherical FITPACK rule. Longitude is selected when
   available and its proposed count is smaller, latitude is already capped,
   no axis has yet been refined, or counts tie after a latitude refinement.

The full pinned core and reference driver were compiled separately:

```cmd
gfortran -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Jbuild/rect-sphere-upstream -c build/rect-sphere-upstream/fixed_core.F90 -o build/rect-sphere-upstream/fixed_core.o
gfortran -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Ibuild/rect-sphere-upstream -Ibuild/gfortran_2806A042DAA8BBF7 validation/rect_sphere_upstream_driver.f90 build/rect-sphere-upstream/fixed_core.o build/gfortran_3C159EB802C22FB9/scifort/libscifort.a -o build/rect-sphere-upstream/fixed_reference.exe
```

This reference validates the source extraction and wrapper against the
established solver, rather than independently validating every algorithmic
choice. The driver shares SciFort's already-tested tensor evaluator; the SciPy
and analytic checks provide independent evaluation references. No downloaded
core, build artifact, or third-party fixture is included in this change.

## Remaining limitations

At least four sites per axis, coordinate-separation and workspace guards,
binary64 conditioning limits, a fixed tolerance/iteration cap, and the
scaled-target lower bound remain documented in the
[API guide](../docs/interpolate_rect_sphere.md). Scattered-data spherical
smoothing and C ABI expansion remain deferred. Source ASCII, line lengths,
local links, whitespace, retained license text, and the guide example were
checked after documentation updates.
