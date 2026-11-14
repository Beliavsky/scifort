# Native scalar automatic-knot smoothing validation

Date: 2026-11-13. Compiler: GNU Fortran (GCC) 17.0.0 20260510
(experimental), Windows. SciPy reference: 1.15.3.

## Scope and behavior

The existing `splrep(task=0,s>0)` now calls an internal perazz/FITPACK-derived
solver for nonperiodic scalar smoothing, degrees 1-5, positive residual weights,
and optionally enlarged support bounds. Interpolation, fixed-knot fitting,
parametric fitting, tuple padding, evaluation/calculus, and existing status
names remain compatible. Default `s=0` is preserved. A numerical failure
returns status 3 and an invalid tuple. The optional package is not needed for
this scalar path. No new linked dependency or C ABI is introduced.

The normalized band solver allocates storage linear in sample count at these
degrees. The wrapper handles scaling and checks the returned representation's
weighted residual. FITPACK's `0.001` relative target tolerance and 20-iteration
smoothing limit are retained. A sufficiently large target returns the weighted
least-squares polynomial. Conservative exponent-range guards, unrepresentable
targets, and requested residuals below achieved roundoff are documented in the
[TCK guide](../docs/interpolate_fitpack_tck_backend.md).

The core is isolated; existing shared numerical implementations are unchanged.
The 24 standalone validation scripts that compile the TCK module directly now
include its new source dependency. The package manifest is unchanged.

## Native and independent reference checks

All four focused native programs passed at both `-O0` and `-O3`, with bounds/
runtime checks, invalid/divide-by-zero/overflow traps, and implicit-interface
errors enabled:

- `test_interpolate_fitpack_smoothing`: 108 checks for degrees 1-5, weighted
  residuals, automatic knots, minimal sample counts, polynomial limits, near
  interpolation, enlarged bounds, padded tuples, scalar/array agreement,
  derivative/integral identities, scale changes, invalid inputs, NaN/infinity,
  and explicit numerical failures for unsupported numerical ranges.
- `test_interpolate_fitpack_tck`: existing tuple/evaluation/calculus cases,
  updated to expect successful positive-s polynomial smoothing.
- `test_interpolate_fitpack_lsq`: existing interpolation and fixed-knot cases.
- `test_interpolate_fitpack_parametric`: existing parametric behavior.

`check_splrep_smoothing.py` passed in both builds:

| Reference | Fits | Adaptive / polynomial | Compared numbers | Maximum normalized difference |
| --- | --- | --- | --- | --- |
| SciPy 1.15.3 | 150 | 109 / 41 | 38,990 | `1.42249e-11` |
| Original pinned perazz core, including crash fixture | 153 | 112 / 41 | 39,934 | `1.42249e-11` |

Comparisons include selected knots, coefficients, sample/query values,
derivatives through order three (or degree), and integrals. Normalized numerical
error is `abs(native-reference)/(1+abs(reference))`, bounded by `3e-8`.
Knots/counts also have independent agreement assertions. Adaptive residuals
must meet the requested target to within `0.0011` relative tolerance; the
observed maximum relative target error was `0.000987891` in either build.

The reference core is compiled directly from inspected upstream source, not
from the adaptation. The optional upstream fixture reads the original
`test_fpknot_crash` data from the pinned test file, retaining its default-real
literal rounding, and uses explicit targets 0.1, 1, and 100. All three native
fits pass and agree with the original Fortran core. SciPy 1.15.3 succeeds for
the first two but reports iteration failure with NaN residual and unordered
knots at `s=100`. This older-reference defect is not copied into SciFort.
The upstream spline-object test and full upstream suite were not run unchanged.
No third-party test data are added to the repository.

Additional standalone regressions passed:

- `compare_fitpack_tck_scipy.py`: 90 cases per checked/optimized build;
  maximum absolute difference `5.96e-8` under its value/derivative tolerances.
- `compare_fitpack_parametric_scipy.py`: 90 cases per build, 12,591 samples
  times channels; maximum absolute error `2.62e-7`, scaled error `6.79e-12`.
- `check_optional_fitpack_interface.py`: synthetic adapter interface check.
  This checks the mock interface, not an upstream optional-package build.
- Changed Python scripts parse successfully.

The required full `fpm test` run passed with exit status 0, including the
new smoothing test and existing interpolation integration paths. The debug
run emitted array-temporary warnings but no test failures. Source ASCII,
new source line lengths, affected local documentation links, Python parsing,
retained license terms, and diff whitespace were checked. Downloaded upstream
source, build products, and scratch scripts remain outside the change.

## Commands and build inputs

Focused compiler flags:

```text
-O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface
-O3 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface
```

Changed objects precede the existing FPM archive for unchanged support
routines. Expanded debug commands are shown in CMD syntax; output directories
were created first. The same compiler flags and object list were used to link
and run each of the four focused test programs listed above.

```cmd
gfortran -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -Jbuild/splrep-debug -Ibuild/splrep-debug -Ibuild/gfortran_2806A042DAA8BBF7 -c src/interpolate_curfit_core.f90 -o build/splrep-debug/core.o
gfortran -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -Jbuild/splrep-debug -Ibuild/splrep-debug -Ibuild/gfortran_2806A042DAA8BBF7 -c src/interpolate_fitpack_compat.f90 -o build/splrep-debug/compat.o
gfortran -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -Jbuild/splrep-debug -Ibuild/splrep-debug -Ibuild/gfortran_2806A042DAA8BBF7 test/test_interpolate_fitpack_smoothing.f90 build/splrep-debug/compat.o build/splrep-debug/core.o build/gfortran_3C159EB802C22FB9/scifort/libscifort.a -o build/splrep-debug/test_interpolate_fitpack_smoothing.exe
build\splrep-debug\test_interpolate_fitpack_smoothing.exe
gfortran -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -Jbuild/splrep-debug -Ibuild/splrep-debug -Ibuild/gfortran_2806A042DAA8BBF7 validation/splrep_smoothing_driver.f90 build/splrep-debug/compat.o build/splrep-debug/core.o build/gfortran_3C159EB802C22FB9/scifort/libscifort.a -o build/splrep-debug/reference.exe
python validation/check_splrep_smoothing.py build/splrep-debug/reference.exe
```

Optimized commands replace `-O0` with `-O3` and the output/module directory
with `build/splrep-release`. The original upstream core and test file were
downloaded into ignored `build/splrep-upstream/` from the pinned URLs in
[provenance](../CODE_PROVENANCE.md). Its standalone reference was built/run as:

```cmd
gfortran -O0 -g -fcheck=all -Jbuild/splrep-upstream-compiled -c build/splrep-upstream/fitpack_core.F90 -o build/splrep-upstream-compiled/core.o
gfortran -O0 -g -fcheck=all -Ibuild/splrep-upstream-compiled validation/splrep_upstream_driver.f90 build/splrep-upstream-compiled/core.o -o build/splrep-upstream-compiled/reference.exe
python validation/check_splrep_smoothing.py build/splrep-debug/reference.exe --upstream-tests build/splrep-upstream/fitpack_curve_tests.f90 --upstream-executable build/splrep-upstream-compiled/reference.exe
python validation/check_splrep_smoothing.py build/splrep-release/reference.exe --upstream-tests build/splrep-upstream/fitpack_curve_tests.f90 --upstream-executable build/splrep-upstream-compiled/reference.exe
python tools/interpolate_reference/compare_fitpack_tck_scipy.py
python tools/interpolate_reference/compare_fitpack_parametric_scipy.py
python tools/interpolate_reference/check_optional_fitpack_interface.py
fpm test
```

## Remaining coverage

Parametric automatic smoothing, periodic scalar smoothing, continuation,
modern `make_splrep`/`make_splprep` constructors, and optional backend/wrapper
expansion remain separate work. Binary64 range and conditioning limits still
apply. The library returns an explicit numerical failure when its guarded
solver cannot provide an acceptable residual.
