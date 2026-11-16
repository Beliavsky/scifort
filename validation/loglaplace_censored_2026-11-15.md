# Log-Laplace censoring validation, 2026-11-15

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
Reference libraries: SciPy 1.15.3 and mpmath 1.3.0.

## Native checks

All 286 checks in `test/test_censored_loglaplace.f90` passed at both `-O0`
and `-O3`, each with `-g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace`.
Checks cover scalar/array agreement, Laplace transformation with exact-density
Jacobian, defaults, mixed observations, negative likelihood, free/fixed
parameters, fit objective agreement, optional controls, analytic shape MLE,
SciPy mixed fits, invalid/empty samples, endpoint densities, infinite interval
endpoints, and equal endpoints as exact observations.

Stress cases cover adjacent intervals at the join, the smallest positive
binary64 shape, nearly unit interval probability, underflowed tail probabilities
with finite logs, overflowing raw displacement, and unrepresentable log tails.
The near-unit probability check uses a ratio to reject a falsely rounded zero.

The mixed sample has exact `[0.5,1.2,2.5]`, left `[0.3]`, right `[3.5]`, and
intervals `[1.1,1.4]`, `[4.0,4.4]`. At shape 1.7, location 0, scale 1.3,
its mpmath log-likelihood is `-15.530991142547814`. At fixed location zero,
SciPy returns shape `0.9300659117552635`, scale `1.240969675551833`.
With shape fixed at 1.7 as well, its fitted scale is `1.2410156250000006`.
For exact `[1,2,4]`, location zero, scale two, the analytic shape MLE is
`3/(2*log(2))`.

Ordinary checks allow absolute errors of `2e-14` to `5e-14`; extreme log tails
and widths allow up to `3e-13`. The near-unit ratio tolerance is `3e-14`.
Different optimizers' SciPy fitted parameters allow `4e-4`; the analytic
shape MLE allows `4e-6`. The guide example compiled and ran with checked flags.

Native checked command (CMD, with matching local FPM dependencies):

```cmd
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Jbuild\loglaplace-native -Ibuild\loglaplace-native -Ibuild\gfortran_2806A042DAA8BBF7 src\special\log_gamma.f90 src\empirical_distribution.f90 src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored_loglaplace.f90 build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\loglaplace-native\test.exe
build\loglaplace-native\test.exe
```

The optimized run used `-O3` and separate `build\loglaplace-optimized` module
and executable paths. Changed modules were compiled explicitly before linking
existing FPM dependencies. Directory hashes are local to this build.

## Independent reference and integration checks

All 2,276 reference comparisons passed in both checked and optimized builds.
The 110 new cases comprise 100 randomized mixed log-Laplace samples and ten
interval stress cases, including adjacent join endpoints at large shape,
underflowed/overflowed standardized distances, clipping at the support endpoint,
overflowing raw displacement, smallest positive shape, and a nearly unit mass.
The mpmath oracle uses direct density/CDF/SF definitions and log probability
differences, at up to 800 digits; the harness seed is 20261112.

The new family's maximum scaled error was `4.65e-15`; the overall maximum
remained `1.08e-12` in an existing gamma stress case. The metric is
`abs(actual-reference)/max(1,abs(reference))`, with threshold `2e-11`.
Configurations are `-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace`
and `-O3 -fcheck=bounds`.

```cmd
python tools\generate_censored_reference.py
python validation\check_censored.py --module-dir build\gfortran_2806A042DAA8BBF7 --library build\gfortran_3C159EB802C22FB9\scifort\libscifort.a
fpm test
python tools\check_test_inputs.py
```

The full `fpm test` suite passed with exit code 0 in 1442.293 seconds
(24 minutes 2 seconds). No numerical changes followed that run.
The standalone reference generator exited 0 and reproduced the mixed reference
and fitted parameters above. The unattended-console-input check, source ASCII
and line-length checks, Python syntax, inventory counts, changed local links,
and `git diff --check` passed.

## Limits

Fits are bounded local searches without global-optimality or uniqueness
guarantees. A free location with shape below one can approach an unbounded
endpoint density; nonfinite objectives are rejected by the existing optimizer.
The existing scalar log-Laplace routines are unchanged and do not gain the new
censored kernels' extreme-range guarantees. No C ABI or dependency is added.
