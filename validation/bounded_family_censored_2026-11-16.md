# Uniform, power-function, and arcsine censoring validation, 2026-11-16

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
Reference libraries: SciPy 1.15.3 and mpmath 1.3.0.

## Native numerical and API checks

All 1,091 checks in `test/test_censored_bounded_family.f90` passed at both
`-O0` and `-O3`, each with
`-g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace`.
Checks cover mixed references and NNLF, ordinary scalar agreement at small,
unit and larger power-function shapes, monotonicity/range, CDF/PPF round trips,
original array likelihoods/fits, free/fixed support parameters, invalid/empty
inputs, optimizer controls, and impossible objectives. Support clipping,
infinite interval endpoints and equal endpoints as exact data are covered.

Identities include uniform as power-function shape one, power-function as
beta `(a,1)`, and arcsine as beta `(1/2,1/2)`. Native checks also cover lower
and upper density/CDF/SF limits, the scalar power-function lower-endpoint fix,
adjacent interval widths, the smallest positive shape, underflowed standardized
density arguments, a huge shape near a compensated upper endpoint, and an underflowed upper
distance compensated by a large shape.
Nearly unit interval masses use ratio checks: a location offset of `1e-100`
with interval `[2e-100,1]` and width one excludes two tiny endpoint regions.
This rejects a zero log mass and checks that both excluded tails are retained.

The mixed sample has exact `[0.1,0.3,0.6]`, left `[0.05]`, right `[0.8]`, and
intervals `[0.2,0.25]`, `[0.7,0.75]`, at location zero and width one:

| Family | Independent mpmath log likelihood |
| --- | ---: |
| Uniform | -10.596634733096073 |
| Power-function, a=0.6 | -10.086476648213639 |
| Arcsine | -10.496569053424526 |

The fixed-support SciPy power-function fit returns `a=0.6711914062499998`.
For exact `[0.25,0.5,1]` at location zero, width one, the analytic shape MLE
is `-3/sum(log(x)) = 1/log(2)`. For uniform exact `[2,3,5]` and fixed location
two, the analytic width MLE is three. The free uniform fit uses the existing
local optimizer: on this sample it can stop at location 1.5 and width about
3.5, since changing location alone is flat and decreasing width alone would
exclude data. The test checks finite objective agreement for that free fit,
and the analytic MLE for the fixed-location case. It does not claim global
optimality for the free fit. Arcsine and power-function free-support fits use
bounded boxes away from singular exact endpoints.

Ordinary scalar/identity/reference absolute log tolerances are `3e-14` to
`6e-14`; extreme small-shape/density checks allow up to `3e-13`.
The huge-shape log CDF near a compensated edge is about `-1e8` and allows
`1e-7` absolute error. Near-unit uniform/arcsine ratio tolerances are `5e-14`
and `2e-13`. The SciPy shape fit allows `2e-4` for different optimizers;
analytic shape/width MLEs allow `4e-6`.

Native checked command (CMD, with matching local FPM dependencies):

```cmd
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Jbuild\bounded_family-native -Ibuild\bounded_family-native -Ibuild\gfortran_2806A042DAA8BBF7 src\stats\powerlaw.f90 src\special\log_gamma.f90 src\empirical_distribution.f90 src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored_bounded_family.f90 build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\bounded_family-native\test.exe
build\bounded_family-native\test.exe
```

The optimized run used `-O3` and separate `build\bounded_family-optimized`
module and executable paths. Changed modules were compiled explicitly before
linking existing FPM dependencies. Directory hashes are local to this build.
The combined guide example compiled and ran with checked flags.

## Independent references and integration

All 3,180 reference comparisons passed in checked and optimized builds.
The 338 new cases comprise 100 mixed samples, nine interval stresses, and
three density/tail extremes per family, plus two compensated upper-edge
power-function tails. Independent mpmath definitions use up to 1,000 digits and the harness seed is 20261112. Cases include smallest
positive and adjacent widths, support clipping, a compensated nearly full
interval, tiny first shape, extreme normalized distances and displacements,
and finite upper edges with tiny offsets. The new family's maximum scaled
error was `4.25e-14`; the overall maximum remained `1.08e-12` in an existing
gamma stress case. The metric is
`abs(actual-reference)/max(1,abs(reference))`, with threshold `2e-11`.
Native ratio checks supplement this absolute scale for tiny log probabilities.

The harness configurations are
`-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace`
and `-O3 -fcheck=bounds`.

```cmd
python tools\generate_censored_reference.py
python validation\check_censored.py --module-dir build\gfortran_2806A042DAA8BBF7 --library build\gfortran_3C159EB802C22FB9\scifort\libscifort.a
fpm test
python tools\check_test_inputs.py
```

The standalone generator exited 0 and reproduced the likelihoods and fixed-
support SciPy power-function fit above. The full `fpm test` suite passed
with exit code 0 in 1351.364 seconds (22 minutes 31 seconds). No
numerical changes followed that successful run. The unattended-console-input
check, source ASCII and line lengths, Python syntax, inventory counts, changed
local links, and `git diff --check` passed.

## Limits

All three interfaces require a finite rounded upper support edge, inheriting
the beta helper contract. Nonfinite endpoint likelihoods are preserved and
rejected as fit objectives. Free arcsine/small-shape power-function endpoints
can make MLEs ill-posed. Fits retain the existing bounded local optimizer's
limits, including flat uniform location directions; success does not imply a
global optimum. Scalar power-function endpoint behavior is corrected, while
other scalar extreme-range behavior remains unchanged. C ABI and wrappers
remain deferred, and no dependency is added.
