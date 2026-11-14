# Periodic native constructors: validation, 2026-11-14

Environment: Windows, GNU Fortran (GCC) **17.0.0 20260510 (experimental)**,
Python **3.13**, SciPy **1.15.3**, binary64.

## Scope and results

Both `make_splrep` and `make_splprep` now accept `bc_type='periodic'`
for closed-input interpolation and automatic smoothing, degrees 1-5.
Periodic supplied knots remain unsupported. Positive smoothing targets use
the first `m-1` sites, with one shared constraint for parametric coordinates.
Modern constructors validate all final observations and weights; the legacy
tuple APIs' ignored-final-observation behavior is preserved.

Both `-O0` and `-O3` builds with bounds/runtime checks and floating-point
invalid/zero/overflow traps passed:

- `test_interpolate_make_periodic`: **206 checks**, including all degrees,
  unwrapped endpoint derivatives, shared residuals, supplied/chord parameters,
  constant limits, minimal and full capacities, explicit capacity failures,
  closing-input validation, malformed bounds, interpolation/smoothing
  extrapolation, subnormal rejection, very small/large parameter intervals,
  common coordinate/weight scaling, and twelve channels.
- Existing scalar constructor (117), parametric constructor (184), closed
  tuple fitting (134), periodic scalar tuple fitting (115), scalar smoothing
  (108), parametric smoothing (125), parametric interpolation, TCK, and LSQ
  regressions. Old unsupported-periodic assertions now check unmatched
  endpoint rejection; the old check counts are preserved.
- Existing scalar and parametric constructor independent checkers (165 fits
  each) retain their previously recorded results and compatibility limits.

`validation/check_periodic_constructors.py`, seed **20261118**, checks
**200 fits**: 50 interpolations and 150 smoothing fits, including 100 active
constraints. Scalar and multichannel cases include supplied/chord parameters,
minimum sample counts, positive weights, explicit capacities, and queries
outside the base interval. Both checked builds produce the same recorded
metrics:

| Independent check | Maximum error | Required limit |
|---|---:|---:|
| SciPy periodic interpolation, normalized coefficient/value/derivative/integral error | 2.01482e-11 | 3e-8 |
| Native tuple evaluated independently with SciPy `BSpline` | 2.1104e-12 | 3e-8 |
| Unwrapped endpoint derivative mismatch, scaled by derivative coefficients | 3.65557e-15 | `2048*epsilon` |
| Active relative residual-target error | 0.000927215 | 0.0011 |
| Cyclic shared-multiplier normal equations | 5.87451e-15 | 2e-8 |

There are **23,020 strict SciPy interpolation comparisons**. Automatic
parameter sites agree within `rtol=2e-14, atol=2e-15`. Independent cyclic
basis/jump matrices check stationarity with one nonnegative multiplier
across coordinates. Returned tuples independently reproduce values,
derivatives through order `min(3,k)`, and integrals. Seam checks evaluate
unwrapped tuples so periodic query wrapping cannot conceal a mismatch.

Smoothing comparisons with SciPy's legacy periodic `splrep`/`splprep` are
diagnostic rather than strict parity assertions: **nine** knot sets differ;
the maximum normalized value/derivative/integral diagnostic is **137.999**,
including sensitive higher derivatives outside the base interval. Native
residuals and cyclic equations are checked independently. Exact modern
automatic-knot and iteration-path parity is not claimed. The core allocation
change leaves default legacy solver settings unchanged; existing tuple
regressions pass. This task does not rerun the original upstream full
object-oriented continuation suite.

Current SciPy **v1.18.0** public source and license were inspected separately
for periodic constructor semantics. No new SciPy implementation is translated
or installed. Installed 1.15.3 provides the periodic interpolation and legacy
fitting references. The v1.18.0 public parametric `s=0` shortcut omits forwarding
`bc_type`; this checker uses explicit periodic `make_interp_spline` to check
the periodic mathematical contract implemented by SciFort.

The inspected v1.18.0 `_fitpack_repro.py` has SHA-256
`82449ae292bb361675439b422cec6c084f9fbfa6519eabec06b1009c51a818bd`;
its root `LICENSE.txt` has SHA-256
`221e59f5e910fd7f94e44f0dac77436a11338c285c6346232e4a850a50da0e94`.
The cached pinned perazz core hash matches the provenance record. Static
checks pass for changed Fortran ASCII and free-form limits, new test/driver
100-column lines, Python syntax, local Markdown file targets, and
`git diff --check`. Cached inspection sources are not tracked.

## Commands

Full repository suite:

```bat
fpm test
```

The full run is launched through Python `subprocess.run`, with stdout/stderr
written directly to the ignored `build/periodic-constructors/full-test.log`
and the native return code recorded separately. This avoids PowerShell's
native-stderr redirection artifact from the previous validation.
The final full run **passed with native exit code 0**, including all 206
new periodic checks and the existing 184/117 constructor checks. Existing
runtime array-temporary warnings from unchanged tests are informational.

For checked builds, the changed core and compatibility module and the public
interpolation module were compiled freshly, then linked before the existing
FPM support archive. The following CMD commands reproduce the debug checks
in this workspace; use another output directory and `-O3` for the optimized
checked build:

```bat
set "out=build\modern-debug"
if not exist "%out%" mkdir "%out%"
set "fitflags=-O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -I%out% -Ibuild\gfortran_2806A042DAA8BBF7 -J%out%"
gfortran %fitflags% -c src\interpolate_curfit_core.f90 -o %out%\interpolate_curfit_core.o
gfortran %fitflags% -c src\interpolate_fitpack_compat.f90 -o %out%\interpolate_fitpack_compat.o
gfortran %fitflags% -c src\interpolate.f90 -o %out%\interpolate.o
set "fitobjects=%out%\interpolate_curfit_core.o %out%\interpolate_fitpack_compat.o %out%\interpolate.o build\gfortran_3C159EB802C22FB9\scifort\libscifort.a"
gfortran %fitflags% test\test_interpolate_make_periodic.f90 %fitobjects% -o %out%\test_interpolate_make_periodic.exe
%out%\test_interpolate_make_periodic.exe
gfortran %fitflags% validation\periodic_constructors_driver.f90 %fitobjects% -o %out%\periodic_constructors_driver.exe
C:\Programs\Python313\python.exe validation\check_periodic_constructors.py %out%\periodic_constructors_driver.exe
```

The listed dependent regressions and two existing independent checkers use
the same compile/link/run pattern. Unchanged support objects are reused;
FPM module/archive directory hashes may differ after configuration changes.
Cached sources, executables, logs, and orchestration scripts remain ignored
build artifacts. See [provenance](../CODE_PROVENANCE.md), the
[periodic guide](../docs/interpolate_periodic_constructors.md), and earlier
[scalar](splrep_periodic_2026-11-13.md) and
[parametric](splprep_closed_2026-11-13.md) periodic validation records.
