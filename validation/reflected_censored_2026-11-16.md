# Double-gamma, double-Weibull, and log-gamma censoring validation, 2026-11-16

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
Reference libraries: SciPy 1.15.3 and mpmath 1.3.0.

## Native numerical and API checks

All 1,085 checks in `test/test_censored_reflected.f90` passed at both `-O0`
and `-O3`, each with
`-g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace`. The test covers mixed likelihoods and NNLF,
ordinary scalar density/CDF/SF agreement on both sides of the location,
monotonicity/range, CDF/PPF round trips, original array fit interfaces,
free/partly fixed/fully fixed bounded fits, SciPy constrained shape fits,
empty/invalid fit objectives, invalid shapes/data, infinite/equal interval
endpoints, center singularities, and the Laplace shape-one identities.

Extreme checks cover reflected intervals with normalized radii near `1e-600`,
the smallest positive interval width, both omitted tails of nearly full
intervals (using ratios to reject a spuriously zero result), overflowing
unscaled reflected distances, and tiny log-gamma shapes that rescue an
overflowing standardized coordinate. Large-shape log-gamma checks cover the
mode and finite density/survival logarithms when `exp(z)` overflows.
The latter reference at `c=1e308,z=710` is `-4.302034083277817e307`, from
100-digit evaluation of the analytic density and the DLMF 8.11(iii) endpoint
term; its first relative correction is below `2e-308`. The large-argument
Stirling remainder is also checked directly against `1/(12*c)`.

The mixed sample has exact `[0.5,1.2,2.5]`, left `[0.3]`, right `[3.5]`,
and intervals `[1.1,1.4]`, `[4,4.4]`, at shape 1.7, location zero, scale 1.3:

| Family | Independent mpmath log likelihood | SciPy shape fit, fixed location and scale |
| --- | ---: | ---: |
| Double-gamma | -16.34465703580217 | 1.8682617187500021 |
| Double-Weibull | -22.75474570177836 | 0.8721679687499999 |
| Log-gamma | -41.125493018751435 | 4.708984375000009 |

Native absolute log tolerances are `6e-14` for ordinary scalar agreement,
`3e-13` for mixed/extreme logs, and `2e-4` for fitted shapes with different
optimizers. Extreme ratios allow `5e-14` to `5e-13` as appropriate.

Native checked command (CMD, with matching local FPM dependencies):

```cmd
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Jbuild\reflected-native -Ibuild\reflected-native -Ibuild\gfortran_2806A042DAA8BBF7 src\special\log_gamma.f90 src\empirical_distribution.f90 src\stats\dweibull.f90 src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored_reflected.f90 build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\reflected-native\test.exe
build\reflected-native\test.exe
```

The optimized run uses `-O3` and separate `build\reflected-optimized`
module and executable paths. Changed modules are compiled explicitly before
linking FPM dependencies; directory hashes are local to this build.

## Independent references and integration

The standalone generator exited 0 and reproduced the likelihoods and
constrained SciPy fits above. The guide example compiled and ran with
checked flags. All 4,667 independent reference comparisons passed in checked
and optimized builds, including 390 new cases: 100 mixed samples per family
and 90 explicit interval stresses. The harness seed is 20261112 and reference
precision is up to 1,000 digits. The new families' maximum scaled log error
was `1.13e-13`; the overall maximum remained `1.08e-12` in an existing gamma
stress case. The metric is `abs(actual-reference)/max(1,abs(reference))`,
with acceptance threshold `2e-11`. Native ratio checks supplement this
metric for tiny log probabilities.

The oracle integrates gamma/Weibull radial intervals directly: subtracting
two reflected CDFs near one half would erase masses as small as
`exp(-138450)` even at 1,000 digits. The harness configurations are
`-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace`
and `-O3 -fcheck=bounds`. The full `fpm test` suite passed with exit code 0
in 1647.870 seconds (27 minutes 28 seconds), including all 1,085 new native
checks. No numerical changes followed that successful run. The unattended
console-input check, source ASCII/line lengths, Python syntax, inventory/API
counts, changed local links, and `git diff --check` passed.

```cmd
python tools\generate_censored_reference.py
python validation\check_censored.py --module-dir build\gfortran_2806A042DAA8BBF7 --library build\gfortran_3C159EB802C22FB9\scifort\libscifort.a
fpm test
python tools\check_test_inputs.py
```

The full suite was run because the Stirling remainder is shared numerical
code and to check integration of the new overloads with existing Fortran
interfaces. Subsequent documentation-only edits did not require another
numerical run.

## Limits

These APIs use binary64 and the existing incomplete-gamma kernel. Extremely
large shapes produce sharply concentrated distributions; coordinate and
`log(shape)` rounding can limit near-mode accuracy. Subnormal results retain
binary64 resolution. For reflected families, free-location fits below unit
shape can have unbounded likelihoods at exact observations. A finite bounded
optimizer result is not a global-MLE guarantee; fully fixed singular or
zero-density objectives return `fit_status_no_finite_objective`. Original
scalar/array interfaces remain available. No C ABI or wrappers are added.
