# Half-normal censoring validation, 2026-11-15

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
References: SciPy 1.15.3 and mpmath 1.3.0.

## Native checks

All 138 checks in `test/test_censored_halfnorm.f90` passed at both `-O0` and
`-O3`, each with `-g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace`.
They cover mixed observations, defaults and keywords, negative likelihoods,
scalar and array agreement, the independent RMS exact-data scale MLE, SciPy
mixed fits, free/fixed location and scale, objective agreement, optional
controls, empty/invalid inputs, impossible fits, and endpoint observations.
Exact support endpoints have finite positive density and are allowed in fits.

Numerical checks include extreme log survival, nearly unit CDF and SF, almost
full-support intervals retaining both excluded tails, clipped and infinite
intervals, equal endpoints as exact observations, adjacent far-tail intervals,
subnormal widths, underflowed standardized widths, overflowing raw displacements,
and unrepresentable log densities. Near-one values use relative ratio checks
so an erroneous zero cannot pass an absolute log-likelihood tolerance.

The independent mixed sample has exact `[0.5,1.2,2.5]`, left `[0.3]`, right
`[3.5]`, and intervals `[1.1,1.4]`, `[4.0,4.4]`. At location 0 and scale 1.3,
its 100-digit mpmath log-likelihood is
`-19.20396573028761760766855468296300866984`. SciPy's fixed-location fit gives
scale `2.7028686618572966`. At location zero, exact observations `[1,2,3]`
have analytic scale MLE `sqrt(14/3)`, independently of any optimization library.

Ordinary scalar/interval comparisons allow `2e-14` absolute error and the
mixed reference `7e-14`. The endpoint log density allows `2e-16`; log-scale
extremes allow up to `5e-13`, and the adjacent far-tail density-integral check
allows `2e-12`. Near-one ratios allow `5e-14` for CDF and `3e-13` for SF and
interval mass. Array-fit objectives allow `2e-14`. Fitted scale allows `4e-4`
for the SciPy comparison and `4e-6` for the analytic RMS MLE.

The documentation example compiled and ran with the checked flags. Products
were placed under ignored `build/` directories. Native checked command (CMD,
with the matching local FPM module and dependency archive paths):

```cmd
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Jbuild\halfnorm-native -Ibuild\halfnorm-native -Ibuild\gfortran_2806A042DAA8BBF7 src\special\log_gamma.f90 src\empirical_distribution.f90 src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored_halfnorm.f90 build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\halfnorm-native\test.exe
build\halfnorm-native\test.exe
```

The optimized run used `-O3` and a separate `build\halfnorm-optimized` module
and executable directory. Changed modules were explicitly compiled before
linking existing FPM dependencies. Directory hashes depend on the local build.

## Reference and integration checks

The reference harness adds 100 mixed half-normal samples and 12 interval stress
cases. Independent mpmath calculations use erf/erfc probabilities and their
logarithmic differences at up to 800 decimal digits. Cases include support
clipping, subnormal and underflowed widths, adjacent tail endpoints, nearly
full-support intervals, wide raw intervals and overflowing raw displacements,
and a scale whose conversion by a square-root-of-two factor could overflow.

```cmd
python validation\check_censored.py --module-dir build\gfortran_2806A042DAA8BBF7 --library build\gfortran_3C159EB802C22FB9\scifort\libscifort.a
fpm test
python tools\check_test_inputs.py
```

All 1,937 reference comparisons passed in each checked and optimized build.
The new half-normal maximum scaled log-likelihood error was `1.09e-14`; the
overall maximum remained `1.08e-12` in an existing gamma stress case. The harness
uses `abs(actual-reference)/max(1,abs(reference)) <= 2e-11`, with seed 20261112.
The two configurations use `-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow
-fbacktrace` and `-O3 -fcheck=bounds`. The console-input, source/Python syntax,
coverage inventory, and local documentation-link checks passed.
The full `fpm test` suite passed with exit code 0 in 1127.275 seconds
(18 minutes 47 seconds). No numerical changes followed that run.

## Limits

The new likelihood uses the existing three-point narrow normal-interval
quadrature and normal log-tail kernels, retaining their limitations. Fits are
bounded local searches with free or fixed location and scale, without global
optimality or uniqueness guarantees. Scalar/array APIs, C ABI, and runtime
dependencies are unchanged. The new censored kernels have their own overflow
and tiny-probability handling; this change does not extend those guarantees to
existing scalar routines.
