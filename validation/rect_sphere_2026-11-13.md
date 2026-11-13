# Rectangular spherical interpolation validation

Date: 2026-11-13. Compiler: GNU Fortran (GCC) 17.0.0 20260510
(experimental), Windows. Reference: Python SciPy 1.15.3 and NumPy.

## Scope and results

Five focused programs passed in both debug and optimized configurations:

- `test_interpolate_rect_sphere`: 79 assertions, including exact cubic latitude
  interpolation, derivatives, poles, periodic wrapping, array/grid agreement,
  metadata, missing-pole options, scaling by `1e200` and `1e-200`, invalid
  inputs, unsupported smoothing, small coordinate gaps, and the size guard.
- `test_interpolate_sphere`: 126 existing scattered spherical spline checks.
- `test_interpolate_bivariate_fitpack`: shared tensor evaluator.
- `test_interpolate_bivariate_calculus`: tensor derivatives/calculus.
- `test_interpolate_rect_bivariate`: neighboring rectangular interpolation API.

Both builds enabled bounds/runtime checks, floating-point exception traps,
and implicit-interface errors. The scope covers the new grid core, changes
to the shared spherical object, its evaluator, and umbrella integration.
Unchanged SVD and basis kernels were not modified. The full SciFort suite and
the full upstream FITPACK package suite were not run; this follows the local
impact-based testing policy. No build configuration or dependency changed.

`check_rect_sphere.py` passed 48 grids in each build: three grid sizes/origins,
two data families, and eight symmetric/asymmetric pole configurations. It
includes uniform/nonuniform grids, constants, smooth harmonics, supplied and
estimated poles, flat poles, exactness flags, and shifted longitude periods.
All nodes, both poles, and random query pairs are compared for value,
`dtheta`, `dphi`, `dtheta=2`, `dphi=2`, and mixed first derivatives, totaling
24,576 reference values per build. Full knots and physical coefficients are
also compared. Trigonometric data test interpolation and agreement, rather
than exact representation of the underlying trigonometric function.

| Comparison | Acceptance bound | Maximum observed, both builds |
| --- | --- | --- |
| Pointwise `abs(native-ref)/(1+abs(ref))` | `5e-8` | `4.49027e-12` |
| Maximum coefficient error / `(1+max(abs(ref coefficients)))` | `3e-9` | `3.36245e-14` |
| Absolute input-node error | `3e-11*(1+max(abs(data)))` | `2.22045e-15` |
| Full knot vectors | `3e-15` absolute | Passed |

The upstream `mnspgr` interpolation case (case 3 in
[`test/fitpack_tests.f90`](https://github.com/perazz/fitpack/blob/1907c7bad9f5c4557970e906bdd4909eb52b2d4c/test/fitpack_tests.f90))
was also checked using its exact 11-by-14 `daspgr_u`, `daspgr_v`, and
`daspgr_r` fixture from
[`test/fitpack_test_data.f90`](https://github.com/perazz/fitpack/blob/1907c7bad9f5c4557970e906bdd4909eb52b2d4c/test/fitpack_test_data.f90).
With C0 poles and no supplied pole observations, both builds interpolate the
fixture with maximum absolute node error `1.33227e-15`. Comparing all six
value/derivative outputs at the 154 nodes with SciPy gives 924 additional
comparisons per build and maximum normalized error `1.63561e-13`.
This check used a temporary script reading the downloaded fixture and calling
the checked-in reference driver; no upstream fixture is redistributed.

An additional constant-grid fixture deliberately supplies a different south
pole: four latitudes from 0.15 to `pi-0.2`, four uniform longitudes, grid 2.5,
south value -0.8, continuity `(false,true)`, and a flat south pole. SciPy
1.15.3 leaves the unspecified north pole at 2.5. Modern pole-range
initialization yields `1.609286231270378`. Independently, the preliminary
latitude space is a cubic with no interior knots; fixing its last two
Bernstein coefficients to -0.8 and solving the remaining two by NumPy least
squares gives the same north pole within `3e-12`. The native fit also
interpolates the input grid and south pole. This is a documented compatibility
difference, not a reference tolerance exception hidden in the grid comparison.

## Commands

The builds reused the existing FPM archive for unchanged modules, compiling
the changed modules ahead of it. Debug flags were:

```text
-O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface
```

Optimized flags used `-O3` instead of `-O0 -g -fbacktrace`, retaining runtime
checks and traps. Commands below use CMD syntax; the build directories were
created before compilation.

```cmd
gfortran -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -Jbuild/rect-sphere-debug -Ibuild/rect-sphere-debug -Ibuild/gfortran_2806A042DAA8BBF7 -c src/interpolate_sphere_grid_core.f90 -o build/rect-sphere-debug/interpolate_sphere_grid_core.o
gfortran -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -Jbuild/rect-sphere-debug -Ibuild/rect-sphere-debug -Ibuild/gfortran_2806A042DAA8BBF7 -c src/interpolate_sphere.f90 -o build/rect-sphere-debug/interpolate_sphere.o
gfortran -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -Jbuild/rect-sphere-debug -Ibuild/rect-sphere-debug -Ibuild/gfortran_2806A042DAA8BBF7 -c src/interpolate.f90 -o build/rect-sphere-debug/interpolate.o
gfortran -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -Jbuild/rect-sphere-debug -Ibuild/rect-sphere-debug -Ibuild/gfortran_2806A042DAA8BBF7 test/test_interpolate_rect_sphere.f90 build/rect-sphere-debug/interpolate.o build/rect-sphere-debug/interpolate_sphere.o build/rect-sphere-debug/interpolate_sphere_grid_core.o build/gfortran_3C159EB802C22FB9/scifort/libscifort.a -o build/rect-sphere-debug/test_interpolate_rect_sphere.exe
build\rect-sphere-debug\test_interpolate_rect_sphere.exe
```

The same link/run command was used for each of the other four tests. Release
commands used `build/rect-sphere-release` and the optimized flags. The reference
driver used the same flags and objects, substituting
`validation/rect_sphere_driver.f90` as source and `reference.exe` as output.

```cmd
python validation/check_rect_sphere.py build/rect-sphere-debug/reference.exe
python validation/check_rect_sphere.py build/rect-sphere-release/reference.exe
```

The [guide example](../docs/interpolate_rect_sphere.md) was compiled and run
with the debug objects and flags. Source ASCII, line lengths, documentation
links, diff whitespace, and retained license terms were checked.

## Remaining limits

Only `s=0` is supported, with at least four sites on each axis. The constructor
guards near-coincident coordinates and aggregate workspace; no general
ill-conditioned-grid guarantee or rank estimate is provided. Pole arrays,
longitude wrapping, rank diagnostics, and the legacy pole-estimation difference
are described in the guide. No C ABI or adaptive smoothing was added.
