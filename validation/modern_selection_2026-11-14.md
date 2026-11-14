# Modern constructor selection validation, 2026-11-14

Environment: Windows, GNU Fortran (GCC) **17.0.0 20260510 (experimental)**,
Python **3.13**, SciPy **1.15.3**, binary64.

## Native scope

Both checked `-O0` and `-O3` builds passed:

- **123 new checks** in `test_interpolate_modern_selection`: degrees 1-5,
  feasible scalar/shared-channel soft-count overshoots, periodic constant-to-
  first-interior refinement, middle-site selection, shared weighted residuals,
  unwrapped seams, scalar/shared agreement, small feasible targets and returned
  parameters. Larger-input compatibility checks use 301 sites, active
  nonperiodic fitting and periodic constant limits against the existing tuple
  paths; they do not establish broader large-input modern coverage.
- Updated scalar constructor **112**, parametric constructor **179** and
  periodic constructor **196** checks. Assertions that automatic modern fits
  must equal legacy automatic fits were removed because knot selection now
  intentionally differs. Independent modern knot and equation checks replace
  that old compatibility assumption.
- Existing general fixed-knot **206**, supplied periodic **1,548**, periodic
  tuple LSQ **1,450**, knot generator **758**, closed tuple **134**, periodic
  tuple **115**, scalar smoothing **108**, parametric smoothing **125**, plus
  parametric interpolation, TCK and LSQ regressions.

The new adapter and the interpolation sources listed below were rebuilt with
bounds/runtime checks, backtraces, implicit-interface errors and floating-
point invalid/divide-by-zero/overflow traps. Supporting modules in the existing
archive were not rebuilt separately at both optimization settings.

## Independent results

Both builds produced the same printed metrics:

| Checker | Coverage / maximum error |
| --- | --- |
| Combined constructor selection | 100 successful fits, 20 explicit infeasible failures, 20 feasible soft-count overshoots, 2,685 selected knots compared |
| Combined residual relative error | 0.000864821 |
| Combined scaled common-multiplier stationarity | 1.21557e-11 |
| Combined representation difference | 3.08747e-12 |
| Scalar constructor | 165 fits, 26,737 strict SciPy comparisons, 55 residual-aware automatic comparisons, zero modern knot differences |
| Scalar scaled stationarity | 1.28397e-14 |
| Scalar residual relative error | 0.000967687 |
| Scalar maximum raw normalized SciPy difference | 1.48599e-7 |
| Scalar coefficient error / stopping bound | 0.00566473 |
| Parametric constructor | 165 fits, 51,779 strict interpolation/supplied-knot comparisons, 80 residual-aware automatic comparisons, zero modern knot differences |
| Parametric scaled stationarity | 5.59034e-15 |
| Parametric residual relative error | 0.000971529 |
| Parametric strict interpolation / supplied-knot errors | 3.12418e-11 / 2.03845e-11 |
| Parametric maximum raw normalized SciPy diagnostic difference | 1.37358 |
| Parametric coefficient error / stopping bound | 0.996875 |
| Periodic constructors | 200 fits, 23,020 strict SciPy interpolation comparisons, 100 active constraints |
| Periodic scaled stationarity / seam mismatch | 6.61376e-15 / 8.19602e-15 |
| Periodic residual relative error | 0.000984125 |
| Periodic legacy smoothing diagnostics | 67 knot differences, maximum raw normalized difference 488.121 |

The combined checker uses seed **20261119**, nonuniform sites/weights,
one/three channels, scalar/parametric and nonperiodic/periodic cases. It covers
zero targets, polynomial/constant limits, active smoothing, feasible soft
stops, infeasible small capacities, automatic chord parameters, enlarged
nonperiodic intervals and near-interpolation targets. It checks selected
knots with `rtol=2e-14, atol=2e-14`, parameter sites with `rtol=2e-14,
atol=2e-15`, every derivative through degree `k` against the returned B-spline,
and unwrapped periodic seams. Successful active residuals must agree with
the target within `0.0011*s`; scaled stationarity must be below `2e-8`.
Infeasible cases must return status 3, an invalid object and no parameters.

Nonperiodic selected vectors compare with installed SciPy's public
`generate_knots`. Periodic vectors compare with hash-verified v1.18.0 Python
selection control flow and an independent NumPy cyclic residual evaluator
that refreshes exteriors. Final capacity vectors are canonicalized explicitly.
This models SciFort's documented periodic correction, not an unmodified
SciPy 1.18 binary. Upstream source versions, hashes and full licenses are in
[provenance](../CODE_PROVENANCE.md) and the
[generator record](generate_knots_2026-11-14.md).

## Residual-aware comparisons

Matching selected knots do not force matching coefficients to roundoff:
the scalar/parametric root searches can stop at different residuals within
the accepted `0.001*s` band. The raw diagnostic differences above include
high derivatives and out-of-domain queries; they are not strict parity
metrics. The periodic legacy diagnostic also includes deliberately different
selected knots. No small blanket tolerance was substituted for those effects.

In addition to target and KKT checks, `spline_stop_checks.py` compares
coefficient differences using the convexity bound for exact minimizers of
the same derivative-jump objective at different residual constraints:

```text
norm(A*(c1-c2))**2 <= abs(fp1-fp2)
norm(c1-c2) <= sqrt(abs(fp1-fp2)) / sigma_min(A)
```

For example, if `fp2 >= fp1`, the second solution has no larger jump penalty.
Convexity gives `dot(D*c1,D*(c2-c1)) <= 0`. Stationarity of the first fit
makes `dot(A*c1-b,A*(c2-c1)) >= 0`; expanding the residual difference yields
the first inequality. The second follows for full-column-rank `A`. The helper
checks that rank and adds `256*epsilon*(norm(b)+norm(A)*max(norm(c)))**2`
to the squared-residual difference as an evaluation roundoff budget. The
reported error/bound ratio must not exceed one. The native normalized normal
equations are separately required to pass; the bound alone is not a proof
that either fit is correct.

Interpolation/supplied-knot strict comparisons retain their prior `3e-8`
normalized forward-error limit. Native derivative evaluation is also checked
against the returned SciPy B-spline representation, independently of the
upstream smoothing stop. No unmodified upstream suite run or exact root-
iteration parity is claimed.

## Commands and full-suite result

The local focused runner compiled these sources in dependency order:

```text
src/interpolate_curfit_core.f90
src/interpolate_fixed_smoothing.f90
src/interpolate_periodic_lsq.f90
src/interpolate_generate_knots.f90
src/interpolate_fitpack_compat.f90
src/interpolate.f90
```

For example, the adapter debug compile command was:

```bat
gfortran -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -Ibuild\modern-selection-debug -Ibuild\gfortran_2806A042DAA8BBF7 -Jbuild\modern-selection-debug -c src\interpolate_fitpack_compat.f90 -o build\modern-selection-debug\interpolate_fitpack_compat.o
```

Tests and drivers were compiled/linked with the same flags and fresh objects,
then `build\gfortran_3C159EB802C22FB9\scifort\libscifort.a`. Optimized runs
used `-O3` and `build\modern-selection-release`, retaining all checks/traps.
The independent commands in each build were:

```bat
python validation\check_make_splrep.py build\modern-selection-debug\make_splrep_driver.exe
python validation\check_make_splprep.py build\modern-selection-debug\make_splprep_driver.exe
python validation\check_periodic_constructors.py build\modern-selection-debug\periodic_constructors_driver.exe
python validation\check_modern_selection.py build\modern-selection-debug\modern_selection_driver.exe
```

The optimized run used the corresponding `modern-selection-release` driver
paths. Native test executables were run directly after linking, including
`build\modern-selection-debug\test_interpolate_modern_selection.exe` and
its optimized counterpart.

Full-suite command:

```bat
fpm test
```

Full-suite result: **passed, native exit code 0**, including all 123 new modern
selection checks and the updated constructor and existing interpolation,
distribution, optimization, signal, sparse and special-function tests. Native
stdout/stderr and the exit code were captured with Python `subprocess.run`,
avoiding shell stderr-formatting artifacts. Existing array-temporary runtime
warnings were present; no test failure was reported.

## Limits

Modern selection and soft counts apply only within the existing generator
limits. Larger inputs retain legacy banded selection and hard capacities.
Infeasible selected-vector targets fail explicitly. Repeated interiors,
broader ill-conditioned inputs, exact periodic selection/root-iteration
parity, C ABI and wrappers remain gaps. See the
[integration guide](../docs/interpolate_modern_selection.md).
