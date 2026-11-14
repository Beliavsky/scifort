# Native `make_splprep` validation, 2026-11-14

Environment: Windows, GNU Fortran (GCC) **17.0.0 20260510 (experimental)**,
Python **3.13**, SciPy **1.15.3**, binary64.

## Scope and results

The native constructor supports degrees 1-5, supplied/chord-length parameters,
interpolation, automatic/fixed-knot smoothing, shared knots, one weighted
residual target, polynomial limits, and a hard capacity. Periodic constructor
boundaries and nonclamped/repeated-interior knot forms are unsupported.

Both `-O0` and `-O3`, with bounds/runtime checks and invalid/zero/overflow
traps, passed:

- `test_interpolate_make_splprep`: **184 checks**, including all degrees,
  supplied/automatic parameters, shared residuals, polynomial limits retaining
  knots, scalar/array evaluation, reversed integration, capacity failures,
  rank failures, NaN/infinities, malformed inputs, large offsets/scales,
  unresolved and subnormal gaps, and twelve channels.
- Existing `make_splrep` (117), closed parametric (134), periodic scalar (115),
  scalar smoothing (108), parametric smoothing (125), parametric interpolation,
  TCK, and fixed-knot LSQ regression programs.
- The new independent checker: **165 fits**, seed **20261117**, 40 interpolation,
  80 automatic smoothing, 40 supplied-knot smoothing, and five nonminimal
  supplied-knot polynomial-limit cases; 96 active residual constraints.
- The existing scalar independent checker: 165 fits, preserving its recorded
  scalar results and compatibility limits.

The new checker makes **51,779 strict SciPy interpolation/supplied-knot
comparisons**, including coefficients, values, derivatives through order
`min(3,k)`, and integrals. Maximum normalized error `abs(a-b)/(1+abs(b))`:
interpolation **3.12418e-11**, supplied knots **2.03845e-11**; limit **3e-8**.
Automatic parameter sites also agree with SciPy within the checker's
`rtol=2e-14, atol=2e-15`; supplied sites and supplied knots are retained exactly.
Polynomial cases use an independent SciPy least-squares polynomial followed
by exact knot insertion, rather than relying on SciPy's inactive-constraint
root search.

All fits meet the residual upper bound `1.0011*s` (interpolation roundoff
allowance `1e-22`). Maximum active relative residual error is **0.000963567**.
Independent B-spline basis matrices and kth-derivative jump matrices verify
normal equations with **one nonnegative multiplier shared by all channels**:
maximum scaled stationarity **5.59034e-15**, limit **2e-9**. Inactive cases
check zero derivative jumps. Independent SciPy `BSpline` evaluation of the
native tuple has maximum normalized error **1.32643e-12**, limit **3e-8**.

Modern automatic numerical parity is deliberately not asserted: **17 of 80**
automatic fits choose different knots from modern SciPy. The maximum modern
value/derivative/integral normalized diagnostic is **230.31**, including
sensitive higher derivatives. These diagnostics are reported, not counted
as passing strict modern comparisons. Both knot selection and legacy root
initialization can differ. The new fixed-knot mode uses reciprocal-mean
initialization; legacy automatic initialization remains unchanged.

The final adapted core also passed the existing **156-fit** parametric
checker against the separately compiled, unchanged pinned original
perazz/FITPACK core, using six fresh-fit cases from its
`test_parametric_fit` data. All **93,431** normalized-original diagnostics
are below `3e-8`; maximum **5.20149e-12**. Maximum active residual relative
error **0.000991072**, independent stationarity **2.27145e-13**.
Raw, unnormalized upstream comparisons retain the previously documented
iteration differences (maximum normalized diagnostic **7.71683**); they
are not strict matches. This exercises fresh fits, not the full upstream
object/continuation suite.

Static checks passed: changed Fortran is ASCII, added lines respect the
free-form limit, new test/driver lines are at most 100 columns, Python parses,
local Markdown file targets resolve, and `git diff --check` is clean.
Cached upstream source/license/test hashes match the provenance records.

The full `fpm test` build completed and ran through the final test program;
all programs reported passes, including the 184 new constructor checks.
The log contains no FPM compilation or execution error. Existing runtime
array-temporary warnings were emitted by unchanged tests.

The PowerShell host wrapper returned **1**, rather than providing a clean
native exit status: redirecting native stderr with `2>&1` converted the first
runtime warning to a `NativeCommandError`. The native FPM exit code was not
captured separately, so an exit-zero result is not claimed. This host behavior
was reproduced with a Python program explicitly exiting zero: plain redirected
invocation returned wrapper status 1, while explicitly capturing and returning
`$LASTEXITCODE` returned 0. No numerical-suite rerun was needed for that
reporting artifact. Direct CMD use does not use this PowerShell redirection.

## Commands

Full repository suite:

```bat
fpm test
```

For checked builds the modified core, compatibility module, and public
umbrella module were compiled freshly, then linked before the existing FPM
support archive. Unchanged support objects were reused. The following CMD
commands reproduce the debug run from this workspace; use a separate output
directory and replace `-O0` with `-O3` for the optimized checked run:

```bat
set "out=build\modern-debug"
if not exist "%out%" mkdir "%out%"
set "fitflags=-O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -I%out% -Ibuild\gfortran_2806A042DAA8BBF7 -J%out%"
gfortran %fitflags% -c src\interpolate_curfit_core.f90 -o %out%\interpolate_curfit_core.o
gfortran %fitflags% -c src\interpolate_fitpack_compat.f90 -o %out%\interpolate_fitpack_compat.o
gfortran %fitflags% -c src\interpolate.f90 -o %out%\interpolate.o
set "fitobjects=%out%\interpolate_curfit_core.o %out%\interpolate_fitpack_compat.o %out%\interpolate.o build\gfortran_3C159EB802C22FB9\scifort\libscifort.a"
gfortran %fitflags% test\test_interpolate_make_splprep.f90 %fitobjects% -o %out%\test_interpolate_make_splprep.exe
%out%\test_interpolate_make_splprep.exe
gfortran %fitflags% validation\make_splprep_driver.f90 %fitobjects% -o %out%\make_splprep_driver.exe
C:\Programs\Python313\python.exe validation\check_make_splprep.py %out%\make_splprep_driver.exe
```

The same compile/link/run pattern was used for the dependent regression
programs listed above and `validation/make_splrep_driver.f90` with
`validation/check_make_splrep.py`. Workspace FPM module/archive directory
hashes may change after configuration changes. A CMD batch file uses doubled
loop variables if these commands are wrapped in a `for` loop.

For the existing original-core comparison:

```bat
gfortran %fitflags% validation\splprep_smoothing_driver.f90 %fitobjects% -o %out%\splprep_regression.exe
C:\Programs\Python313\python.exe validation\check_splprep_smoothing.py %out%\splprep_regression.exe --upstream-executable build\splprep-upstream\upstream.exe --upstream-tests build\splrep-upstream\fitpack_curve_tests.f90
```

The cached original executable is unchanged from the earlier
[parametric validation](splprep_smoothing_2026-11-13.md); it links the exact
pinned upstream core, not the SciFort adaptation. Cached source/executables
and orchestration scripts remain ignored build artifacts. See
[provenance](../CODE_PROVENANCE.md) for the pinned source and retained license.
