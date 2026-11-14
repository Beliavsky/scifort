# Native parametric automatic-knot smoothing validation

Date: 2026-11-13. Compiler: GNU Fortran (GCC) 17.0.0 20260510
(experimental), Windows. SciPy: 1.15.3. Upstream core: perazz/FITPACK commit
`1907c7bad9f5c4557970e906bdd4909eb52b2d4c`.

## Scope

`scifort_interpolate:splprep(points,s=positive)` now fits nonperiodic curves
with a shared automatically selected knot vector and a single total weighted
residual target across coordinate channels. Degrees 1-5, positive observation
weights, explicit or chord-length parameters, and enclosing parameter bounds
are supported. Existing interpolation, fixed-knot fitting, default `s=0`,
unpadded coefficient storage, higher channel counts, and status conventions
remain. Periodic smoothing, continuation, modern constructors, and C ABI work
remain separate tasks.

The change adds `fppara` and its vector rotation helpers beside the scalar
solver in the existing core module. Existing shared numerical helper bodies,
manifest, external dependencies, and standalone build source lists are
unchanged. Full provenance and retained BSD terms are recorded in
[CODE_PROVENANCE.md](../CODE_PROVENANCE.md).

The wrapper normalizes all coordinates by one common scale, not one scale per
channel, preserving the joint metric. Parameters are shifted before scaling;
chord differences are subtracted before scaling when representable and scaled
before subtraction only to avoid overflow. The final physical tuple is checked
against the residual target after coefficient and knot rescaling. Numerical or
allocation failure returns status 3 with an invalid tuple. The original
relative target tolerance `0.001` and 20 smoothing iterations are retained.

## Native checks

Five focused programs passed in checked and optimized builds with bounds/
runtime checks, invalid/divide-by-zero/overflow traps, and implicit-interface
errors enabled:

- `test_interpolate_fitpack_parametric_smoothing`: 125 checks covering all
  degrees, shared knots, joint weighted residuals, unpadded shapes, explicit/
  chord parameters, enlarged bounds, polynomial limits, minimal sample counts,
  scalar/vector evaluation, derivatives, scale changes, zero/constant channels,
  twelve channels, NaN/infinity and invalid inputs, periodic status preservation,
  and explicit numerical-range failures.
- `test_interpolate_fitpack_parametric`: previous interpolation/fixed-knot and
  extrapolation cases, updated to expect positive-s polynomial smoothing.
- `test_interpolate_fitpack_smoothing`: all 108 scalar smoothing checks.
- `test_interpolate_fitpack_tck`: existing tuple/evaluation/calculus checks.
- `test_interpolate_fitpack_lsq`: existing interpolation/fixed-knot checks.

The new regressions include alternating coordinates near `1e308`, large common
coordinate offsets, parameter ranges near `1e200`, unrepresentable physical
spans/intervals, and parameters `1e16+2*(i-1)`. The latter exposed cancellation
in a divide-before-shift parameter transformation: target 0.35 initially gave
physical residual `0.44774006338062444`. Origin-first normalization and physical
verification give `0.35006620710120917`, within the requested tolerance.

The standalone `compare_fitpack_parametric_scipy.py` regression also passed
90 interpolation/fixed-knot cases in each build, covering 12,591 samples times
channels. Maximum absolute error was `2.62e-7`; scaled error was `6.79e-12`.

## Independent smoothing checks

`check_splprep_smoothing.py` exercises 150 cases across degrees 1-5 and one
through five channels, irregular explicit parameters or automatic parameters,
nonuniform weights, enlarged bounds, multiple smoothing levels, and polynomial
limits. An optional upstream test file adds six fresh fits drawn from
`test_parametric_fit`, preserving the default-real rounding of its literals.
These cover explicit parameters with targets 100, 1, 0.05 and 0.25, plus cubic
and quintic chord-parameter fits at 0.25. The stateful upstream continuation
sequence and full upstream suite are not run unchanged. No upstream data are
vendored.

All checks below passed for both native `-O0` and `-O3` builds.

| Check | Coverage | Acceptance | Observed maximum |
| --- | --- | --- | --- |
| Separately compiled original upstream core on equivalent normalized inputs | 156 fits: 96 adaptive, 60 polynomial; 93,431 coefficient/value/derivative numbers | normalized error < `3e-8` | `5.20149e-12` |
| Total adaptive weighted residual | 150 ordinary and 6 upstream-data fits | relative target error < `0.0011` | `0.000991072` |
| Independent penalized normal equations | Every native fitted representation, all channels jointly | scaled stationarity error < `2e-9` | `2.27145e-13` |
| Independent SciPy evaluation of native tuples | Sample/query values and derivatives through min(3,k) | normalized error < `3e-8` | `3.56652e-12` |
| Raw SciPy smoothing-path prediction bound | 150 fits, shared-knot/parameter agreement | squared distance / budget <= 1 | `0.139471` |
| Raw upstream smoothing-path prediction bound | 156 fits, shared-knot/parameter agreement | squared distance / budget <= 1 | `0.127897` |

Normalized numerical error is `abs(native-reference)/(1+abs(reference))`.
Knot and parameter comparisons also assert agreement within `3e-13` absolute/
relative tolerances. The original upstream executable is linked against its
unmodified core, separately compiled from the inspected pinned source. Reference
normalization is a mathematically equivalent affine transformation of the
parameter domain and common scaling of coordinates/weights/target.

For the independent stationarity check, SciPy constructs the design matrix
`B` and the jumps `D` of the kth derivative of each basis function at interior
knots. The checker finds one nonnegative multiplier `lambda` for all channels
and verifies

```text
B.T @ W**2 @ (B @ C - X) + lambda * D.T @ D @ C = 0.
```

The residual is scaled by the data normal-equation norm and the magnitudes of
the two gradient terms. Polynomial fits check the weighted least-squares normal
equations directly. This checks the mathematical smoothing problem independently
of the translated band solver and requires a common multiplier across channels.

## Iteration-path differences and limits

The pinned perazz core initializes its smoothing guess to the mean observation-
factor diagonal; SciPy 1.15.3 uses its reciprocal. Global weight normalization
can also change the iteration path. Both stop on a residual tolerance, so
coefficient or derivative identity across these paths is not promised.

Raw-reference comparisons therefore assert shared knots/parameters, the target
contract, and a smoothing-path prediction bound rather than a blanket small
coefficient/derivative error. For fixed knots the residual components in a
penalty eigenbasis vary monotonically, giving

```text
||W*(prediction_native-prediction_reference)||**2
  <= abs(residual_native-residual_reference) + roundoff_budget.
```

The checker uses `1e-10*max(1,sum((W*X)**2))` for the roundoff budget. Maximum
weighted sample distance divided by `sqrt(target)` was `0.015988` against raw
SciPy and `0.0148116` against raw upstream. These are diagnostics of different
accepted points on the same smoothing path, not roundoff errors.

High derivatives of noisy quintic curves can be especially sensitive. Across
91,387 SciPy diagnostic numbers, the largest normalized raw-reference difference
was `4.75739`, in a third derivative (maximum absolute difference in that case
about `1.12656e6`). Raw upstream diagnostics reached `7.71683`, also in a third
derivative (case maximum absolute difference about `15093.9`). These large
high-derivative differences are explicitly not claimed to pass the tight
normalized-upstream threshold; that threshold applies to equivalently normalized
inputs in the first table row. The independent tuple evaluation and stationarity
checks still pass. Residual tolerance is not a bound on derivative accuracy.

Binary64 range and conditioning limits remain. The physical parameter span must
be representable and distinct physical knots at least `tiny(1.0_dp)` apart.
Normalized gaps and relative weights use the documented conservative exponent
guards. Targets below achieved rounding error and failed rescaling return an
explicit numerical failure. More than ten channels remain a SciFort extension
to the legacy SciPy interface; a twelve-channel regression covers that path.

## Commands and build inputs

Focused flags (GNU Fortran version above):

```text
-O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface
-O3 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface
```

Changed objects precede the existing FPM archive for unchanged support modules.
Expanded debug commands in CMD syntax follow; directories were created first.
Each of the five test programs above used the same link/run pattern. Optimized
commands replace `-O0` with `-O3` and the output directory with
`build/splprep-release`.

```cmd
gfortran -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -Jbuild/splprep-debug -Ibuild/splprep-debug -Ibuild/gfortran_2806A042DAA8BBF7 -c src/interpolate_curfit_core.f90 -o build/splprep-debug/core.o
gfortran -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -Jbuild/splprep-debug -Ibuild/splprep-debug -Ibuild/gfortran_2806A042DAA8BBF7 -c src/interpolate_fitpack_compat.f90 -o build/splprep-debug/compat.o
gfortran -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -Jbuild/splprep-debug -Ibuild/splprep-debug -Ibuild/gfortran_2806A042DAA8BBF7 test/test_interpolate_fitpack_parametric_smoothing.f90 build/splprep-debug/compat.o build/splprep-debug/core.o build/gfortran_3C159EB802C22FB9/scifort/libscifort.a -o build/splprep-debug/test_interpolate_fitpack_parametric_smoothing.exe
build\splprep-debug\test_interpolate_fitpack_parametric_smoothing.exe
gfortran -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -Jbuild/splprep-debug -Ibuild/splprep-debug -Ibuild/gfortran_2806A042DAA8BBF7 validation/splprep_smoothing_driver.f90 build/splprep-debug/compat.o build/splprep-debug/core.o build/gfortran_3C159EB802C22FB9/scifort/libscifort.a -o build/splprep-debug/reference.exe
python validation/check_splprep_smoothing.py build/splprep-debug/reference.exe
python validation/check_splprep_smoothing.py build/splprep-release/reference.exe
```

The upstream core object is reused from the previous scalar validation, built
with the command below. Its source and license hashes were verified before
adaptation; the full retained terms match modulo whitespace. Only the new
standalone reference driver is added to the repository.

```cmd
gfortran -O0 -g -fcheck=all -Jbuild/splrep-upstream-compiled -c build/splrep-upstream/fitpack_core.F90 -o build/splrep-upstream-compiled/core.o
gfortran -O0 -g -fcheck=all -Ibuild/splrep-upstream-compiled validation/splprep_upstream_driver.f90 build/splrep-upstream-compiled/core.o -o build/splprep-upstream/upstream.exe
python validation/check_splprep_smoothing.py build/splprep-debug/reference.exe --upstream-executable build/splprep-upstream/upstream.exe --upstream-tests build/splrep-upstream/fitpack_curve_tests.f90
python validation/check_splprep_smoothing.py build/splprep-release/reference.exe --upstream-executable build/splprep-upstream/upstream.exe --upstream-tests build/splrep-upstream/fitpack_curve_tests.f90
python tools/interpolate_reference/compare_fitpack_parametric_scipy.py
fpm test
```

The required full `fpm test` run on the final numerical code passed with exit
status 0, including the new 125-check test and existing interpolation integration
paths. The run emitted array-temporary runtime warnings but no test failures.
Source ASCII, new/changed source widths, local documentation links, diff
whitespace, Python syntax, pinned source hash, and retained license terms passed
their checks. Upstream downloads, scratch scripts, and build products remain
outside the change.
