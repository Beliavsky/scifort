# Gibrat censoring validation, 2026-11-15

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
Reference libraries: SciPy 1.15.3 and mpmath 1.3.0.

## Native checks

All 154 checks in `test/test_censored_gibrat.f90` passed at both `-O0` and
`-O3`, each with `-g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace`.
They cover mixed observations, defaults/keywords, negative likelihood,
lognormal shape-one identities, scalar and array agreement, support endpoints,
clipped and infinite-endpoint intervals, equal endpoints, invalid parameters,
empty data, extreme log tails, adjacent intervals, underflowed standardized
intervals, and overflowing raw displacements.

Fit tests cover the independent geometric-mean exact-data scale MLE, mixed
SciPy fit comparisons, objective agreement, lognormal parent-fit agreement,
original array-fit compatibility, optional controls, fixed scale, rejected free
location, invalid bounds, empty samples, and zero-density exact endpoints.
Bounds and returned parameters have two entries, `[loc,scale]`.

The independent mixed sample has exact `[0.5,1.2,2.5]`, left `[0.3]`, right
`[3.5]`, and intervals `[1.1,1.4]`, `[4.0,4.4]`. At location 0, scale 1.3,
its 100-digit mpmath log-likelihood is
`-14.60381704491309970298104829462452042778`. SciPy's Gibrat fit at fixed
location gives scale `1.3704101562500006`. The exact sample `[1,2,4]` has
analytic scale MLE 2 at location zero, independent of lognormal fitting code.

Ordinary scalar/interval comparisons allow `2e-14` absolute error, the mixed
reference `5e-14`, adjacent intervals `3e-14`, and extreme log tails and huge
raw displacements `3e-13`. Parent likelihood/fit identities are exact. The
SciPy fitted scale allows `2e-4` for different optimizers; the analytic MLE
allows `3e-6`.

The new guide example compiled and ran with the checked flags. Build products
were placed only under ignored `build/` directories. The native checked command
from CMD, with the matching local FPM module/archive paths, was:

```cmd
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Jbuild\gibrat-native -Ibuild\gibrat-native -Ibuild\gfortran_2806A042DAA8BBF7 src\special\log_gamma.f90 src\empirical_distribution.f90 src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored_gibrat.f90 build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\gibrat-native\test.exe
build\gibrat-native\test.exe
```

The optimized run used `-O3` with a separate `build\gibrat-optimized` module
and executable directory. The changed modules were compiled explicitly before
linking the existing FPM dependency archive. Directory hashes are local to this
build.

## Reference and integration checks

The reference harness adds 100 mixed Gibrat samples and four interval stress
cases, using independent mpmath evaluation of the lognormal mathematical
specification at shape 1. Stress cases include adjacent endpoints at 1 and
`exp(40)`, underflowed standardized endpoints, and an overflowing raw
displacement. High-precision calculations use up to 800 decimal digits.

```cmd
python validation\check_censored.py --module-dir build\gfortran_2806A042DAA8BBF7 --library build\gfortran_3C159EB802C22FB9\scifort\libscifort.a
fpm test
python tools\check_test_inputs.py
```

All 1,825 reference comparisons passed in each checked and optimized build.
The new Gibrat maximum scaled log-likelihood error was `1.24e-14`; the overall
maximum remained `1.08e-12` in an existing gamma stress case. The harness uses
`abs(actual-reference)/max(1,abs(reference)) <= 2e-11`, with seed 20261112.
The two configurations use `-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow
-fbacktrace` and `-O3 -fcheck=bounds`. The unattended-console-input check,
source/Python syntax, inventory counts, and local documentation links passed.
The full `fpm test` suite passed with exit code 0 in 1337.840 seconds
(22 minutes 18 seconds). No numerical changes followed that run.

## Limits

The new wrappers add no numerical approximation or runtime dependency. They
inherit existing lognormal-kernel limitations. Censored fitting requires fixed
location and is a bounded local scale search; shape is always 1. Existing
scalar/array APIs and C ABI remain unchanged.
