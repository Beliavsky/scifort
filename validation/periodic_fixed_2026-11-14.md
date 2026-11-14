# Supplied periodic constructor knots: validation, 2026-11-14

Environment: Windows, GNU Fortran (GCC) **17.0.0 20260510 (experimental)**,
Python **3.13**, SciPy **1.15.3**, binary64.

## Scope and checked results

Native `make_splrep` and `make_splprep`, degrees 1-5, accept full strict
periodic knots with positive `s`. The supplied vector is retained exactly,
including constant limits. The closing observation is validated and omitted
from the scalar/shared parametric residual, preserving the established
periodic constructor contract. Legacy periodic `task=-1` remains unsupported.

Both `-O0` and `-O3` builds with runtime/bounds checks and floating-point
invalid/zero/overflow traps passed:

- `test_interpolate_periodic_fixed`: **1,548 checks** covering all degrees,
  1-12 independent periodic coefficients, constant fits retaining knots,
  active residuals, unwrapped endpoint derivatives, 12 shared channels,
  invalid extensions/endpoints/repetitions/nonfinite knots, `t`/`nest`
  exclusion, interpolation rejection of `t`, infeasible residuals, rank
  failures, common coordinate/weight scaling, and tiny/huge parameter scales.
- Existing periodic constructors (206), scalar constructor (117), parametric
  constructor (184), closed tuple fitting (134), periodic scalar tuple fitting
  (115), scalar smoothing (108), parametric smoothing (125), parametric
  interpolation, TCK, and LSQ regressions.
- Existing independent scalar/parametric constructor checkers (165 fits each)
  and periodic constructor checker (200 fits), retaining their prior recorded
  numerical results and automatic-selection compatibility limits.

`check_periodic_fixed.py`, seed **20261114**, checks **130 cases**: **57 active
constraints**, **70 constant fits**, and **3 expected rank failures**. Cases
include degrees 1-5, 1/2/3/4/7/13/23 independent coefficients, random nonuniform
sample sites and clustered knots, positive nonuniform weights, scalar and
multichannel coordinates, derivatives, integrals, and extrapolation.

The independent checker constructs a cyclic basis from SciPy `BSpline`,
forms derivative jumps at interior knots, and solves augmented least-squares
systems with one shared nonnegative multiplier. A scalar root search matches
the achieved native residual, separating coefficient correctness from the
FITPACK `0.001` relative target tolerance. It does not call a FITPACK smoothing
routine as its coefficient oracle. Rank-deficient observation factors are
checked against numerical failure rather than assuming an inverse exists.

| Independent check | Maximum error | Required limit |
|---|---:|---:|
| Normalized constrained coefficient difference | 2.67859e-10 | 2e-7 |
| Independent tuple value/derivative/integral representation | 6.98492e-10 | 2e-8 |
| Unwrapped seam, scaled by derivative coefficients and knot-gap derivative bound | 2.04928e-16 | 2e-12 |
| Active relative residual-target error | 0.000896621 | 0.0011 |

Both checked builds produce these same printed metrics. Supplied knots and
parameter sites compare exactly. Small independent coefficient counts test
sequential cyclic extension; clustered knots exercise observation intervals
containing more than one knot. No query wrapping is used for seam checks.

## Provenance and limitations

The change extends the existing attributed periodic kernels from
[perazz/FITPACK](https://github.com/perazz/fitpack/blob/1907c7bad9f5c4557970e906bdd4909eb52b2d4c/src/fitpack_core.F90),
commit `1907c7bad9f5c4557970e906bdd4909eb52b2d4c`. Its exact cached source and
license were inspected again; source SHA-256 is
`02e413063f7df855c45eed6119130f7ab7a07cfa3e2f6ac0e42b6f496df414ec`.
The original BSD notices and full license remain. No new implementation or
linked dependency is imported, and no new upstream full object-oriented
continuation-suite run is claimed.

SciPy v1.18.0 `_fitpack_repro.py` and `_bsplines.py` (including `fpcheck`)
were inspected for interface and knot contracts, without translating their
implementation. The inspected `_bsplines.py` SHA-256 is
`b099333d1dc57d0d4ad5300c04cd2767350048cac2b9dee464e40ab469510caa`.
The prior [constructor record](periodic_constructors_2026-11-14.md) gives the
reproduction-source and SciPy license hashes. Installed SciPy 1.15.3 supplies
the independent basis and root/least-squares tools; no strict modern periodic
smoothing implementation parity is claimed.

Accepted knots are finite, strict periodic extensions, with exact base
endpoints and normalized extension agreement within `64*epsilon`. Broader
forms and rank-deficient regularized solutions remain gaps. Existing scaling,
rank, normalized-gap, and physical-residual guards apply. These checks do not
guarantee conditioning for all accepted inputs. Modern automatic selection,
soft capacities, continuation, C ABI, and wrappers are outside this change.
See the [periodic guide](../docs/interpolate_periodic_constructors.md) and
[provenance](../CODE_PROVENANCE.md).

## Commands

Full repository suite:

```bat
fpm test
```

The full suite is launched with Python `subprocess.run`, writing stdout and
stderr directly to `build/periodic-fixed/full-test.log` and saving the native
return code in `full-test-exit.txt` in that directory. This avoids PowerShell's
native-stderr redirection artifact. The full suite **passed with native exit
code 0**, including all 1,548 new checks. Existing runtime array-temporary
warnings from unchanged tests remain informational.

The checked builds compile the changed core/compatibility modules and the
public interpolation module freshly, linking these objects before the
existing FPM support archive. These CMD commands reproduce the debug checks
in this workspace; use `build\fixed-release` and `-O3` for the optimized
checked build. FPM directory hashes can differ after configuration changes.

```bat
set "out=build\fixed-debug"
if not exist "%out%" mkdir "%out%"
set "fitflags=-O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -I%out% -Ibuild\gfortran_2806A042DAA8BBF7 -J%out%"
gfortran %fitflags% -c src\interpolate_curfit_core.f90 -o %out%\interpolate_curfit_core.o
gfortran %fitflags% -c src\interpolate_fitpack_compat.f90 -o %out%\interpolate_fitpack_compat.o
gfortran %fitflags% -c src\interpolate.f90 -o %out%\interpolate.o
set "fitobjects=%out%\interpolate_curfit_core.o %out%\interpolate_fitpack_compat.o %out%\interpolate.o build\gfortran_3C159EB802C22FB9\scifort\libscifort.a"
gfortran %fitflags% test\test_interpolate_periodic_fixed.f90 %fitobjects% -o %out%\test_interpolate_periodic_fixed.exe
%out%\test_interpolate_periodic_fixed.exe
gfortran %fitflags% validation\periodic_fixed_driver.f90 %fitobjects% -o %out%\periodic_fixed_driver.exe
C:\Programs\Python313\python.exe validation\check_periodic_fixed.py %out%\periodic_fixed_driver.exe
```

Dependent regressions and existing independent checkers use the same
compile/link/run pattern. Build logs, cached inspection sources, orchestration
scripts, and executables remain ignored artifacts.
