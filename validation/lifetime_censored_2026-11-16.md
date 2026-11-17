# Exponential-power, exponentiated Weibull, and Mielke validation, 2026-11-16

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
Reference libraries: SciPy 1.15.3 and mpmath 1.3.0.

## Native numerical and API checks

All 1,305 checks in `test/test_censored_lifetime.f90` passed at both `-O0`
and `-O3`, each with
`-g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace`.
They cover mixed references and NNLF, scalar density/CDF/SF agreement,
CDF monotonicity, PPF round trips, original array-fitting interfaces,
fixed/partly free/fully free fits, constrained SciPy shape fits, invalid
parameters and data, empty samples, singular and zero endpoint densities,
support clipping, equal and infinite interval endpoints, fit iteration
limits, and impossible or singular fixed objectives. Both four-parameter
families have tests with all four parameters free.

Identities include unit-shape exponential-power/Gompertz, unit-exponent
Weibull, and equal-shape Mielke/Fisk. Native ratio checks retain small
nonzero logs, including exponential-power survival near `-1e-240` and
Mielke CDF near `-6.907755278982137e-298`. Opposing exponentiated-Weibull
shapes `a=1e-300,c=1e308` retain a finite log CDF near `-6.907755278982137e10`
at `x=1e-300`, although `c*log(x)` overflows. Other checks include an
underflowing Mielke shape ratio, large-shape Weibull survival at `x=1000`,
and overflowing endpoint shape products.

The scalar exponentiated-Weibull upper-tail correction is checked through
ordinary scalar/censored agreement and the exact unit-shape identity
`logsf(80,1,1)=-80`. The old log-CDF formula lost precision when its
inner probability approached one; the complementary `log1p` branch fixes
that observed cancellation.

The mixed sample uses exact `[0.5,1.2,2.5]`, left `[0.3]`, right `[3.5]`,
intervals `[1.1,1.4]`, `[4,4.4]`, location zero, and scale 1.3:

| Family and shapes | Independent log likelihood | SciPy fit with other parameters fixed |
| --- | ---: | ---: |
| Exponential-power, b=1.7 | -1097.7750366596047 | b=0.3581054687499994 |
| Exponentiated Weibull, a=1.7, c=1.2 | -15.776394603077286 | a=1.3926757812500012 |
| Mielke, k=1.7, s=1.2 | -14.709016801278233 | k=1.302636718750001 |

Ordinary scalar log tolerances are `7e-14`, mixed references `3e-13`,
probability round trips `8e-14`, and related-family identities `1e-13`.
Extreme log-value ratios allow `3e-13`; constrained shape comparisons
allow `2e-4` for the different optimizers.

Checked native commands (CMD, with matching local FPM dependencies):

```cmd
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Jbuild\lifetime-native -Ibuild\lifetime-native -Ibuild\gfortran_2806A042DAA8BBF7 src\stats\exponweib.f90 src\special\log_gamma.f90 src\empirical_distribution.f90 src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored_lifetime.f90 build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\lifetime-native\test.exe
build\lifetime-native\test.exe
```

The optimized run used `-O3` and separate `build\lifetime-optimized`
module/executable paths. Changed numerical and integration modules were
compiled explicitly before linking the archive. These standalone compiler
commands require the build directories to exist.

## Independent references

An isolated run of 440 new cases passed in checked and optimized builds,
with maximum scaled log-likelihood error `2.86e-13`. They include 300
mixed-censoring samples and 140 finite interval/point stresses. Ordinary
references use 100 decimal digits; new stresses use 750 digits. Cases
include adjacent and subnormal gaps, shape ratios beyond the binary64
range, opposing shapes, overflowing transformed logs, nearly full
intervals, and scales from `1e-300` to `1e300`.

The combined existing/new reference command is:

```cmd
python validation\check_censored.py --module-dir build\gfortran_2806A042DAA8BBF7 --library build\gfortran_3C159EB802C22FB9\scifort\libscifort.a
```

The metric is `abs(actual-reference)/max(1,abs(reference))`, with acceptance
threshold `2e-11`. Native ratio checks supplement it for tiny nonzero logs.
The oracle avoids nested exponentials beyond the representable log range;
nonfinite objectives are tested natively, while finite references are
required for the randomized comparison. The combined run passed all 5,909
comparisons per build at checked `-O0` and optimized `-O3`. The maximum
new-family error was `7.13e-14`; the overall maximum was `1.08e-12` in
an existing gamma stress case. The isolated run starts the RNG directly
at the configured seed, while the combined run preserves earlier samples.

## Integration and documentation

The full `fpm test` run passed with exit code 0 in 2,282.833 seconds
(38 minutes 3 seconds), including all 1,305 new native checks. Shared
likelihood/fitting/statistics integration and the scalar numerical correction
required the full suite. No executable changes followed the successful run;
final validation-report edits were documentation only.

The guide example compiled and ran with checked compiler flags. Source
ASCII/line-length checks, Python syntax, changed Markdown links, inventory
and API-list consistency (50 supported plus 59 remaining),
`python tools\check_test_inputs.py`, and `git diff --check` passed.
The complete `python tools\generate_censored_reference.py` run completed
successfully and reproduced the mixed references and SciPy fit values above.

## Limits

Binary64 is the public precision. Endpoint classification for exponentiated
Weibull follows the rounded binary64 product `a*c`; overflow gives zero
endpoint density. Other scalar extreme-parameter limitations are unchanged.
The censored callbacks have separate validated log-domain paths.

The bounded optimizer is local; free shapes, location, and scale can need
more than the default 1,000 sweeps. Free-fit tests allow 10,000 sweeps;
four-parameter tests use tolerance `1e-6`. Endpoint singularities can make
unrestricted-location likelihoods unbounded. Success does not certify a
global MLE or parameter identifiability. Fully fixed singular or impossible
objectives return `fit_status_no_finite_objective`.

Native Fortran coverage takes priority; C ABI and wrapper work remain deferred.
