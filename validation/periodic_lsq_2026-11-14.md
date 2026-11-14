# Periodic fixed-knot tuple least squares: validation, 2026-11-14

Environment: Windows, GNU Fortran (GCC) **17.0.0 20260510 (experimental)**,
Python **3.13**, SciPy **1.15.3**, binary64.

## Scope and results

Native `splrep` and `splprep` now support `per=.true., task=-1`, degrees 1-5.
Scalar input supplies interior knots; parametric input supplies a full vector,
whose endpoint/exterior template is canonicalized. The weighted least-squares
objective uses the first `m-1` observations/weights. Existing tuple layouts,
closing-leg chord parameters, and evaluation policies are retained.

Both `-O0` and `-O3` builds with bounds/runtime checks and floating-point
invalid/zero/overflow traps passed:

- `test_interpolate_periodic_lsq`: **1,450 checks**, covering degrees 1-5,
  2-12 independent coefficients, recovery of known cyclic coefficients,
  unwrapped endpoint derivatives, padded scalar/unpadded parametric layouts,
  twelve coordinate channels, full-vector canonicalization, input immutability,
  ignored nonfinite final observations/weights, invariance to valid `s`,
  rank-deficient minimum-norm solutions, invalid knots/tasks, failure cleanup,
  supplied/chord parameters, balanced coordinate/weight scaling, tiny/huge
  physical intervals, and explicit coefficient/channel limits.
- Existing general supplied-knot (206), periodic supplied-knot (1,548), periodic
  constructor (206), scalar constructor (117), parametric constructor (184),
  closed tuple (134), periodic scalar tuple (115), scalar smoothing (108),
  parametric smoothing (125), parametric interpolation, TCK, and LSQ regressions.
- All five earlier independent checkers (165 scalar constructor fits, 165
  parametric constructor fits, 200 periodic constructor fits, 130 periodic
  fixed-knot smoothing fits, 160 general fixed-knot smoothing fits), preserving
  their printed metrics and compatibility limits.

The two new guide examples also compiled and ran in the checked debug build.

`check_periodic_lsq.py`, seed **20261115**, checks **120 fits**: 60 scalar and
60 parametric, including 20 automatic chord-parameter cases and **three
rank-deficient** observation matrices. Fits cover degrees 1-5, nonuniform
sites and weights, 2/3/7/15 independent coefficients, noisy data, clustered
knots, endpoint derivatives, and queries outside the base interval.

For all cases the independent reference constructs a cyclic B-spline basis
with SciPy `BSpline` and solves weighted least squares using NumPy. The 117
full-rank cases additionally compare with SciPy's legacy periodic `splrep` and
`splprep`, for **56,946 strict coefficient/value/derivative/integral comparisons**.
Rank-deficient cases use the independent minimum-norm solve rather than
assuming SciPy accepts singular input.

| Independent check | Maximum in both checked builds | Required limit |
|---|---:|---:|
| Normalized cyclic coefficient error | 1.08996e-14 | 2e-7 |
| Normalized value/derivative/integral error | 5.38056e-11 | 2e-7 |
| Scaled least-squares normal equations | 2.47986e-16 | 2e-12 |
| Scaled unwrapped endpoint derivative mismatch | 4.93175e-16 | 2e-12 |

Both builds produce the same printed new-checker metrics. SciPy comparisons
use `rtol=2e-7, atol=2e-7`. Parameters and knots compare with `rtol=2e-14,
atol=2e-15`. The normal-equation residual is divided by
`norm(A)*(norm(b)+norm(A)*norm(c))`. Seam errors use a coefficient/knot-gap
derivative bound; no periodic query wrapping can conceal a mismatch.
The scalar integral uses `splint`; parametric integration checks the returned
tuple through the existing native array B-spline object's integral method.

## Contract and provenance

Fixed-knot fitting minimizes residual and does not enforce a smoothing target.
Valid `s` is ignored after the existing finite/nonnegative validation. Empty
interior sets are rejected, matching the inspected FITPACK task validation.
An explicit `task=-1` remains necessary in SciFort; no implicit task selection
is introduced. More than ten coordinates and rank-deficient minimum-norm
fitting are documented SciFort extensions.

The dense path is limited to **256 independent cyclic coefficients**,
**10,000 sites**, and **256 coordinate channels**. Finite-range, normalized
sample/knot-gap, relative-weight, physical-gap, and coefficient representability
checks apply. These guards do not guarantee good conditioning for all accepted
inputs. `task=1`, modern `generate_knots` selection, and soft capacities remain
separate work. C ABI and wrapper expansion remain deferred.

The new module composes existing SciFort B-spline, pivoted QR, and Jacobi SVD
routines; no third-party implementation or dependency is imported. The
mathematical specification follows Dierckx (1993) and public SciPy
[`splrep`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.splrep.html)
and [`splprep`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.splprep.html)
contracts. The pinned perazz/FITPACK `percur`/`clocur` task validation and
periodic canonicalization blocks and BSD-3-Clause license were inspected:
commit `1907c7bad9f5c4557970e906bdd4909eb52b2d4c`, source SHA-256
`02e413063f7df855c45eed6119130f7ab7a07cfa3e2f6ac0e42b6f496df414ec`.
That source is a behavioral reference; no new original-upstream full suite
run or imported fixed-knot solver is claimed. Existing attributed kernels
and notices are unchanged. See [provenance](../CODE_PROVENANCE.md) and the
[guide](../docs/interpolate_periodic_lsq.md).

## Commands

Full repository suite:

```bat
fpm test
```

The full run uses Python `subprocess.run`, writing stdout/stderr directly to
`build/periodic-lsq/full-test.log` and recording the native return code in
`full-test-exit.txt` in that directory. The full suite **passed with native
exit code 0**, including all 1,450 new checks and existing dependent regressions.
Existing runtime array-temporary warnings remain informational.

The checked builds compile the existing core and general fixed-smoothing
module, the new periodic least-squares module, compatibility module, and
public interpolation module freshly, linking those objects before the FPM
support archive. These CMD commands reproduce the new debug checks in this
workspace; use `build\periodic-lsq-release` and `-O3` for the optimized checked
build. FPM directory hashes may differ after configuration changes.

```bat
set "out=build\periodic-lsq-debug"
if not exist "%out%" mkdir "%out%"
set "fitflags=-O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface -I%out% -Ibuild\gfortran_2806A042DAA8BBF7 -J%out%"
gfortran %fitflags% -c src\interpolate_curfit_core.f90 -o %out%\interpolate_curfit_core.o
gfortran %fitflags% -c src\interpolate_fixed_smoothing.f90 -o %out%\interpolate_fixed_smoothing.o
gfortran %fitflags% -c src\interpolate_periodic_lsq.f90 -o %out%\interpolate_periodic_lsq.o
gfortran %fitflags% -c src\interpolate_fitpack_compat.f90 -o %out%\interpolate_fitpack_compat.o
gfortran %fitflags% -c src\interpolate.f90 -o %out%\interpolate.o
set "fitobjects=%out%\interpolate_curfit_core.o %out%\interpolate_fixed_smoothing.o %out%\interpolate_periodic_lsq.o %out%\interpolate_fitpack_compat.o %out%\interpolate.o build\gfortran_3C159EB802C22FB9\scifort\libscifort.a"
gfortran %fitflags% test\test_interpolate_periodic_lsq.f90 %fitobjects% -o %out%\test_interpolate_periodic_lsq.exe
%out%\test_interpolate_periodic_lsq.exe
gfortran %fitflags% validation\periodic_lsq_driver.f90 %fitobjects% -o %out%\periodic_lsq_driver.exe
C:\Programs\Python313\python.exe validation\check_periodic_lsq.py %out%\periodic_lsq_driver.exe
```

Dependent regressions and earlier independent checkers use the same
compile/link/run pattern. Cached sources, orchestration scripts, executables,
logs, and objects remain ignored build artifacts.
