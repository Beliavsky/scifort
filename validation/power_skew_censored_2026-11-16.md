# Power-normal, power-lognormal, and skew-Cauchy validation, 2026-11-16

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
Reference libraries: SciPy 1.15.3 and mpmath 1.3.0.

## Native API and numerical checks

All 1,279 checks in `test/test_censored_power_skew.f90` passed at both
`-O0` and `-O3`, each with
`-g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace`.
They cover mixed references and NNLF, scalar density/CDF/SF agreement,
monotonicity, PPF round trips, original array-fit compatibility, fixed and
partly/fully free fits, SciPy constrained shape fits, fit status handling,
invalid data and parameters, empty samples, support endpoints, and
equal/infinite interval endpoints. The power-lognormal test frees all four
parameters. Unit-power normal/lognormal and zero-skew Cauchy identities
are checked separately.

Extreme native checks include a tiny power rescuing a normal log density
and survival near `-5e299`, overflowing negative standardized coordinates,
and skewness at the nearest interior binary64 values to both endpoints.
A ratio check verifies that power-normal log survival at `x=-40,c=1e300`
retains `-3.65589354091503e-50` rather than rounding to zero.

The mixed sample has exact `[0.5,1.2,2.5]`, left `[0.3]`, right `[3.5]`,
and intervals `[1.1,1.4]`, `[4,4.4]`, with location zero and scale 1.3:

| Family and shapes | Independent log likelihood | SciPy fit, fixed other parameters |
| --- | ---: | ---: |
| Power-normal, c=1.7 | -34.00443518723973 | c=0.2776367187499993 |
| Power-lognormal, c=1.7, s=0.8 | -17.71829372436072 | c=0.7646484374999998 |
| Skew-Cauchy, a=0.4 | -15.408647115964277 | a=0.8061875000000008 |

Ordinary scalar log tolerances are `7e-14`, mixed references `3e-13`,
probability round trips `8e-14`, identities `1e-13`, and extreme log-value
ratios `3e-13`. Shape comparisons allow `2e-4` for different optimizers.

The checked native commands (CMD, using matching local FPM dependencies) were:

```cmd
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Jbuild\power_skew-native -Ibuild\power_skew-native -Ibuild\gfortran_2806A042DAA8BBF7 src\special\log_gamma.f90 src\empirical_distribution.f90 src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored_power_skew.f90 build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\power_skew-native\test.exe
build\power_skew-native\test.exe
```

The optimized run used `-O3` and separate `build\power_skew-optimized`
module/executable paths. Changed modules were compiled explicitly before
linking the existing archive. The build directories must exist before
running these standalone compiler commands.

## Independent high-precision comparisons

The command below passed 5,469 comparisons per build in checked `-O0`
and optimized `-O3` configurations. There are 409 new cases: 300 seeded
mixed-censoring samples and 109 interval/point stress cases. Normal power
references use high-precision erfc and the mathematical survival/density
identities; skew-Cauchy references use high-precision arctangents. Ordinary
references use 100 decimal digits; new stress cases use 750.

```cmd
python validation\check_censored.py --module-dir build\gfortran_2806A042DAA8BBF7 --library build\gfortran_3C159EB802C22FB9\scifort\libscifort.a
```

The metric is `abs(actual-reference)/max(1,abs(reference))`, with acceptance
threshold `2e-11`. The maximum new-family error was `1.24e-13`; the overall
maximum was `1.08e-12` in an existing gamma stress case. Native ratio checks
supplement this metric for tiny nonzero log probabilities. An additional
run of only the 409 new cases, starting the RNG at the same seed, passed both configurations
with maximum error `1.61e-13`.

Stress cases include adjacent and subnormal endpoint gaps, powers from
`1e-300` through `1e300`, lognormal shape `1e-300`, normalized distances
outside the binary64 range, tiny/huge scales, nearly full intervals, and
skewness nearest to +/-1. Nonfinite mathematical log objectives are checked
in the native program; the randomized comparison requires finite references.

## Integration and documentation checks

The full `fpm test` run passed with exit code 0 in 1,689.931 seconds
(28 minutes 10 seconds), including the new 1,279-check program. It was
required because the public censored-likelihood, fitting dispatch, and
statistics integration modules are shared. No executable changes followed
the successful run; final validation-report edits were documentation only.

The guide's new example compiled and ran with the same checked compiler
flags. `python tools\generate_censored_reference.py` completed successfully
and reproduced the mixed likelihoods and SciPy fit values above.
`python tools\check_test_inputs.py`, source ASCII/line-length checks,
Python syntax checks, changed Markdown links, inventory/API-list consistency
(47 supported plus 62 remaining), and `git diff --check` passed.

## Limits

Binary64 is the public precision. The bounded optimizer is local; jointly
free power, location, and scale can require more than the default 1,000
sweeps. The native free-fit checks use up to 10,000 sweeps; the four-parameter
power-lognormal case uses tolerance `1e-6`. Success does not certify a global
optimum or parameter identifiability. Impossible fixed objectives return
`fit_status_no_finite_objective`.

The new likelihoods reuse existing log-distance and lognormal-coordinate
kernels. Their extreme behavior is tested independently; existing scalar
distribution implementations are unchanged and can have different extreme
parameter limits. C ABI and wrapper additions remain deferred.
