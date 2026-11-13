# Fixed-knot spherical spline validation

Date: 2026-11-13. Compiler: GNU Fortran (GCC) 17.0.0 20260510
(experimental), Windows. Independent numerical reference: SciPy 1.15.3.

## Scope and outcomes

The new `scifort_interpolate_sphere` module fits the reduced fixed-knot
FITPACK spherical spline space using existing native B-spline and SVD kernels.
The `scifort_interpolate` umbrella re-exports its public API. Existing kernels,
dependencies, build configuration, and C ABI are unchanged.

Six focused programs passed at both `-O0` and `-O3`, with `-fcheck=all`,
`-ffpe-trap=invalid,zero,overflow`, and `-Werror=implicit-interface`:

- `test_interpolate_sphere`: 126 assertions for exact cubic latitude data,
  values, coordinate derivatives, poles, periodic wrapping, paired/grid
  queries, metadata, weights, data ordering, rank warnings, and invalid inputs.
  Scaling tests use data and weights multiplied by `1e200` and `1e-200`.
  Guard tests cover unrepresentable relative weights and physical coefficients.
- `test_interpolate_bspline`: basis/collocation dependency.
- `test_interpolate_fitpack_tck`: existing public TCK interface.
- `test_interpolate_bivariate_fitpack`: tensor evaluator dependency.
- `test_interpolate_bivariate_calculus`: bivariate derivatives and related
  integration/insertion APIs exported through the same umbrella.
- `test_linalg_svd`: SVD/minimum-norm solver dependency.

The rank-deficient analytic case has two observations at the north pole,
values `[2,4]` and residual weights `[1,2]`. It returns rank one, a valid fit,
status 2, fitted pole value 3.6, and weighted squared residual 3.2. General
rank-deficient equivalence with FITPACK is not claimed.

The reference checker passed 20 full-rank fits in each configuration,
covering zero, one, two, three, and four interior latitude knots; uniform
and nonuniform longitude/latitude knots; nonuniform positive weights;
unsorted scattered data; duplicated pole/seam observations; constants;
an exactly representable cubic latitude function; smooth angular data;
and noisy data. The smooth trigonometric functions test fitting and API
agreement, not exact representation by a finite spline space.

Each fit was compared at 37 query pairs for the value and five partial
derivatives: `(1,0)`, `(0,1)`, `(2,0)`, `(0,2)`, and `(1,1)`. This gives
4,440 independent value/derivative comparisons per build. Full knot vectors,
all physical tensor coefficients, and weighted residuals were also compared.
All reference fits returned full model rank and status zero.

| Comparison | Acceptance bound | Maximum observed in both builds |
| --- | --- | --- |
| Pointwise `abs(value-reference)/(1+abs(reference))` | `3e-9` | `1.53870e-12` |
| Max coefficient error divided by `1+max(abs(native coefficient))` | `2e-10` | `1.39571e-14` |
| Residual error divided by `1+reference residual` | `2e-10` | `9.64313e-15` |
| Exact constant/cubic latitude values | `1e-11` absolute | Passed |
| Full knot vectors | `2e-15` absolute | Passed |

The example in `docs/interpolate_sphere.md` compiled and ran with debug
checks. It returned approximately `0.48588862514538034`, rank 10, and
weighted squared residual `0.0062650199895133065` for its cosine data.
This example is an approximation; the exact cosine at its query is 0.5.

The full repository numerical suite was not run. Tests cover the added module,
the umbrella integration, and the existing numerical dependencies, following
the repository's impact-based test policy. No unrelated modules were retested.

## Commands

Both builds reused the existing FPM archive for unchanged modules and compiled
the new module and interpolation umbrella ahead of it. Debug flags were:

```text
-O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface
```

Optimized flags replaced `-O0 -g -fbacktrace` with `-O3`. Runtime checks and
traps remained enabled. Debug commands, shown in CMD-compatible syntax:

```cmd
gfortran -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -Jbuild/sphere-debug -Ibuild/sphere-debug -Ibuild/gfortran_2806A042DAA8BBF7 -c src/interpolate_sphere.f90 -o build/sphere-debug/sphere.o
gfortran -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -Jbuild/sphere-debug -Ibuild/sphere-debug -Ibuild/gfortran_2806A042DAA8BBF7 -c src/interpolate.f90 -o build/sphere-debug/interpolate.o
gfortran -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -Jbuild/sphere-debug -Ibuild/sphere-debug -Ibuild/gfortran_2806A042DAA8BBF7 test/test_interpolate_sphere.f90 build/sphere-debug/interpolate.o build/sphere-debug/sphere.o build/gfortran_3C159EB802C22FB9/scifort/libscifort.a -o build/sphere-debug/test_interpolate_sphere.exe
build\sphere-debug\test_interpolate_sphere.exe
```

The same link/run command was used for the other five test programs. Optimized
commands used `build/sphere-release` for module outputs, changed objects, and
executables. The reference driver linked `validation/sphere_spline_driver.f90`
with the corresponding `sphere.o` and the existing archive using the same flags.

```cmd
python validation/check_sphere_spline.py build/sphere-debug/reference.exe
python validation/check_sphere_spline.py build/sphere-release/reference.exe
```

The guide program used the same debug link command with its source extracted
to `build/sphere-debug/guide_example.f90`. Diff whitespace, new Fortran ASCII
and line lengths, local documentation links, and retained license terms were
checked. The pinned upstream license and the retained license have identical
text after whitespace normalization. No downloaded raw source or build artifact
is included in the change.

## Compatibility limits

The native interface supports real binary64 fixed-knot fitting, at least three
interior longitude knots, at most 256 reduced coefficients and 256 independent
longitude coefficients, and a 2,000,000-entry aggregate design-workspace guard.
It uses SVD singular-value truncation rather than FITPACK's banded rank test,
so rank decisions and deficient solutions may differ. Periodically wrapping
query longitudes and the rank/model-dimension accessors are documented SciFort
extensions. Automatic smoothing and rectangular-grid spherical fitting remain
unimplemented. See the guide for the complete numerical contract.
