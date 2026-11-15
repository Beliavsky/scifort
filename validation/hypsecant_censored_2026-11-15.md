# Hyperbolic-secant censoring validation, 2026-11-15

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
Reference libraries: SciPy 1.15.3 and mpmath 1.3.0.

## Native checks

All 216 checks in `test/test_censored_hypsecant.f90` passed at both `-O0` and
`-O3`, each with `-g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace`.
Checks cover mixed observations, defaults and keywords, negative likelihood,
reflection, scalar/array agreement, free/fixed fitting parameters, an independent
symmetric-sample location MLE, SciPy mixed fits, objective agreement, optional
controls, invalid/empty data, infinite interval endpoints, and equal endpoints
as exact observations.

Stress checks cover extreme log CDF/SF, a nearly unit CDF, wide crossing intervals
retaining both omitted tails, positive and negative adjacent tail intervals,
subnormal crossing widths, underflowed standardized same-side and crossing
widths, overflowing raw displacements, and unrepresentable log densities and
interval log masses. Tiny near-one logs use relative ratio checks so a spurious
zero cannot pass an absolute log-likelihood tolerance.

The independent mixed sample has exact `[0.5,1.2,2.5]`, left `[0.3]`, right
`[3.5]`, and intervals `[1.1,1.4]`, `[4.0,4.4]`. At location 0.4 and scale 1.3,
its 100-digit mpmath log-likelihood is
`-16.31028768324200973331054425185042559084`. SciPy's unrestricted fit gives
location `1.7600566713787145`, scale `1.3906362832127406`. The symmetric exact
sample `[-1,0,1]` has location MLE zero at fixed scale one: its score vanishes
by symmetry, and the log-likelihood is strictly concave in location.

Ordinary scalar/interval comparisons allow `2e-14` absolute error, the mixed
reference `7e-14`, and the center log density `3e-16`. Log-scale extremes allow
up to `5e-13`, and adjacent far-tail density-integral checks allow `3e-12`.
Near-one ratios allow `3e-14`. SciPy fitted parameters allow `4e-4` for different
optimizers; the analytic symmetric-sample location allows `3e-6`.

The new guide example compiled and ran with the checked flags. All products
were placed under ignored `build/` directories. Native checked command (CMD,
using the matching local FPM module directory and dependency archive):

```cmd
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Jbuild\hypsecant-native -Ibuild\hypsecant-native -Ibuild\gfortran_2806A042DAA8BBF7 src\special\log_gamma.f90 src\empirical_distribution.f90 src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored_hypsecant.f90 build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\hypsecant-native\test.exe
build\hypsecant-native\test.exe
```

The optimized run used `-O3` and a separate `build\hypsecant-optimized` module
and executable directory. Changed modules were compiled explicitly before
linking existing FPM dependencies. Directory hashes are local to this build.

## Reference and integration checks

All 2,049 reference comparisons passed in both checked and optimized builds.
The 112 new cases comprise 100 mixed hyperbolic-secant samples, including
positive/negative tails and intervals crossing the center, and 12 stress cases.
Independent mpmath calculations use direct densities and logarithmic CDF/SF
differences, with up to 800 decimal digits and seed 20261112. They cover
subnormal/underflowed widths, adjacent tails, near-full-support intervals,
wide raw intervals, overflowing displacements, and large scales.

The maximum scaled error for the new family was `1.92e-15`; the overall maximum
remained `1.08e-12` in an existing gamma stress case. The harness uses
`abs(actual-reference)/max(1,abs(reference)) <= 2e-11`. Configurations are
`-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace` and
`-O3 -fcheck=bounds`.

```cmd
python tools\generate_censored_reference.py
python validation\check_censored.py --module-dir build\gfortran_2806A042DAA8BBF7 --library build\gfortran_3C159EB802C22FB9\scifort\libscifort.a
fpm test
python tools\check_test_inputs.py
```

The standalone reference generator exited 0 and reproduced the likelihood
and fit above. SciPy emits noninteger-shape warnings while fitting the existing
Erlang comparison; these do not indicate a reference or library failure.
Console-input, source/Python syntax, coverage counts, and local links passed.
The full `fpm test` suite passed with exit code 0 in 1687.101 seconds
(28 minutes 7 seconds). No numerical changes followed that run.

## Limits

The new likelihood reuses the existing log-arctangent and exponential-mass
kernels, retaining their numerical limitations. Fitting is a bounded local
search with free or fixed location and scale, without uniqueness or global
optimality guarantees. Scalar/array APIs, C ABI, and dependencies are unchanged.
The new censored kernels have their own overflow and tiny-probability handling;
this change does not extend those guarantees to existing scalar routines.
