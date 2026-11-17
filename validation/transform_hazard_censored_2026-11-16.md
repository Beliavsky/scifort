# Fatigue-life, Johnson SU, and generalized exponential validation, 2026-11-16

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
Reference libraries: SciPy 1.15.3 and mpmath 1.3.0.

## Native API and numerical checks

All 1,307 checks in `test/test_censored_transform_hazard.f90` passed at
both `-O0` and `-O3`, each with
`-g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace`.
They cover mixed references and NNLF, ordinary scalar density/CDF/SF
agreement, monotonicity, PPF round trips, original array-fit compatibility,
fixed/partly free/fully free fits, constrained SciPy shape fits, invalid
parameters/data, empty samples, support endpoints, equal/infinite interval
endpoints, fit iteration limits, and impossible fixed objectives. Johnson
SU is tested with negative, zero, and positive first shape. Free-fit checks
include all four Johnson parameters and all five generalized-exponential
parameters.

Identities include fatigue-life reciprocal tails, Johnson reflection, and
the generalized-exponential limit as the second rate tends to zero.
Extreme checks include a fatigue-life shape `1e-300` at the scale point,
underflowed standardized distances, a Johnson scale shape `1e300` at its
center, a subnormal Johnson interval, an overflowing Johnson product
rescued by an opposing first shape, and generalized-exponential rates
whose sum overflows. With `a=1e-300,b=1e300,c=1e-300` at `x=1`, the latter
family retains log density `-0.5` and its adjacent-interval mass.

The mixed sample has exact `[0.5,1.2,2.5]`, left `[0.3]`, right `[3.5]`,
intervals `[1.1,1.4]`, `[4,4.4]`, location zero, and scale 1.3:

| Family and shapes | Independent log likelihood | SciPy fit with other parameters fixed |
| --- | ---: | ---: |
| Fatigue-life, c=1.7 | -14.670431355890935 | c=1.2929687500000007 |
| Johnson SU, a=1.7, b=1.2 | -42.437977289016104 | a=-1.2559570312500017 |
| Generalized exponential, a=1.7, b=1.2, c=0.8 | -23.52587064286926 | a=0.15458984374999926 |

Native ordinary scalar log tolerances are `7e-14`, mixed references
`3e-13`, probability round trips `8e-14`, related-family identities
`1e-13` to `2e-13`, and extreme log-value/ratio checks `3e-13`.
Constrained shape fits allow `2e-4` for different optimizers.

Checked native commands (CMD, using matching local FPM dependencies):

```cmd
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Jbuild\transform_hazard-native -Ibuild\transform_hazard-native -Ibuild\gfortran_2806A042DAA8BBF7 src\special\log_gamma.f90 src\empirical_distribution.f90 src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored_transform_hazard.f90 build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\transform_hazard-native\test.exe
build\transform_hazard-native\test.exe
```

The optimized run used `-O3` and separate `build\transform_hazard-optimized`
module/executable paths. Changed modules were compiled explicitly before
linking the existing archive. Standalone compiler commands require these
build directories to exist.

## Independent references

An isolated 409-case run passed in checked and optimized builds, with
maximum scaled log-likelihood error `4.23e-14`. It includes 300 seeded
mixed-censoring samples and 109 finite interval/point stresses. Ordinary
references use 100 decimal digits; new stress cases use 750. They include
adjacent/subnormal gaps, overflowing or underflowing standardized ratios,
compensating extreme rates, tiny/huge shapes, support crossings, and
nearly full intervals.

The oracle computes normal transforms directly at high precision, with
independent complementary erfc tails. The generalized-exponential small
remainder is evaluated through a confluent-hypergeometric identity rather
than the Fortran series. Nonfinite log objectives are tested natively;
the comparison requires finite references. The combined command is:

```cmd
python validation\check_censored.py --module-dir build\gfortran_2806A042DAA8BBF7 --library build\gfortran_3C159EB802C22FB9\scifort\libscifort.a
```

The error metric is `abs(actual-reference)/max(1,abs(reference))`, with
threshold `2e-11`. Native checks supplement this metric for endpoint and
underflow behavior. The combined run passed all 6,318 comparisons per build
at checked `-O0` and optimized `-O3`. Its new-family maximum was `4.23e-14`;
the overall maximum was `1.08e-12` in an existing gamma stress case. The
isolated run starts the RNG directly at the configured seed; the combined
run preserves all earlier samples before adding the new families.

## Integration and documentation checks

The full `fpm test` run passed with exit code 0 in 3,434.581 seconds
(57 minutes 15 seconds), including all 1,307 new native checks. Shared
censored-likelihood, fitting dispatch, and statistics exports required
the full integration suite. No executable changes followed the successful
run; final validation-report edits were documentation only.

The complete `python tools\generate_censored_reference.py` run reproduced
the mixed likelihoods and SciPy fits above. The new guide example compiled
and ran with checked compiler flags. `python tools\check_test_inputs.py`,
source ASCII/line-length checks, Python syntax, changed Markdown links,
inventory/API-list consistency (53 supported plus 56 remaining), and
`git diff --check` passed.

## Limits

Binary64 is the public precision. Very large opposing Johnson shape terms
make `a+b*asinh(z)` ill-conditioned; rounding in the inverse hyperbolic
sine can greatly amplify the error near cancellation. The tests establish
safe arithmetic and the listed extreme cases, not uniform relative accuracy
through those cancellation regimes. Normal-transform tails share the
existing normal kernel's accuracy limits.

The bounded optimizer is local, and joint shape/location/scale fitting can
be weakly identifiable. Free-fit tests allow up to 10,000 sweeps; fully
free four/five-parameter cases use tolerance `1e-6`. A successful fit does
not certify a global optimum or identifiability. Impossible fixed objectives
return `fit_status_no_finite_objective`.

Existing scalar distribution implementations are unchanged and can have
different extreme-parameter limits. C ABI and wrapper additions remain deferred.
