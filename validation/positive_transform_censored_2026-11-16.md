# Inverse-gamma, beta-prime, and Maxwell censoring validation, 2026-11-16

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
Reference libraries: SciPy 1.15.3 and mpmath 1.3.0.

## Native numerical and API checks

All 1,093 checks in `test/test_censored_positive_transform.f90` passed at
both `-O0` and `-O3`, each with
`-g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace`.
Coverage includes mixed references and NNLF, scalar density/CDF/SF agreement,
CDF monotonicity and PPF round trips, original array-fit interfaces,
free/partly fixed/fully fixed fits, constrained SciPy shape/scale fits,
empty/invalid fit objectives, invalid parameters/data, support endpoints
and clipping, equal/infinite interval endpoints, and beta-prime singularity. Shape-sum boundary checks verify symmetric
acceptance of the largest finite shape paired with a tiny positive shape,
and its finite endpoint density with first shape one.

Identities include inverse-gamma/inverse-Weibull at unit shape, beta-prime
with first shape one as Lomax, and the Maxwell/gamma CDF transformation.
Extremes include adjacent intervals, Maxwell radii near `1e-600`, inverse-gamma
and beta-prime far tails, both omitted tails of nearly full intervals, and
beta-prime log likelihoods near the largest representable magnitude.
The large beta-prime reference `-1.373112378206389e308` is independently
computed at 100 digits from the density specification with
`a=1e305,b=1.7e308,x=1e-300,scale=1e300`; its log CDF has the same rounded
binary64 value. Ratio checks distinguish tiny negative logs from zero.

The mixed sample has exact `[0.5,1.2,2.5]`, left `[0.3]`, right `[3.5]`,
and intervals `[1.1,1.4]`, `[4,4.4]`, at location zero, scale 1.3:

| Family and shapes | Independent log likelihood | SciPy constrained fit |
| --- | ---: | --- |
| Inverse-gamma, a=1.7 | -15.94191834999753 | a=1.3391601562500006, fixed location/scale |
| Beta-prime, a=1.7 and b=2.3 | -14.960817730539578 | a=2.2651367187500027, fixed b/location/scale |
| Maxwell | -19.55213612016153 | scale=1.4474879497966218, fixed location |

Native absolute log tolerances are `6e-14` for ordinary scalar agreement,
`3e-13` for mixed/ordinary narrow intervals, and `2e-12` for large extreme
logs. Relative log-value ratios allow `3e-13` to `5e-13`; constrained shapes
allow `2e-4` and Maxwell scale `2e-5` for different optimizers.

Native checked command (CMD, with matching local FPM dependencies):

```cmd
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Jbuild\positive_transform-native -Ibuild\positive_transform-native -Ibuild\gfortran_2806A042DAA8BBF7 src\special\log_gamma.f90 src\empirical_distribution.f90 src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored_positive_transform.f90 build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\positive_transform-native\test.exe
build\positive_transform-native\test.exe
```

The optimized run uses `-O3` and separate `build\positive_transform-optimized`
module/executable paths. Changed modules were compiled explicitly before
linking FPM dependencies; directory hashes are local to this build.

## Independent references and integration

The standalone generator exited 0 and reproduced the likelihoods and
constrained fits above. The guide example compiled and ran with checked
flags. All 5,060 independent reference comparisons passed in checked and
optimized builds. The 393 new cases comprise 100 mixed samples per family,
89 interval stresses, two complementary beta-prime tails with normalized
coordinates near `1e-600` or `1e600`, and two largest-finite-shape boundary
tails with the other shape fixed at one. The harness seed is 20261112. New cases
use 100-digit mpmath integrals; existing stress cases retain up to 1,000 digits.
Gamma references choose the smaller tail before subtraction and complement
both omitted tails for nearly full intervals. Beta-prime point references
retain both transformed coordinates before complementing the smaller tail.

The new families' maximum scaled log error was `5.23e-13`; the overall maximum
remained `1.08e-12` in an existing gamma stress case. The metric is
`abs(actual-reference)/max(1,abs(reference))`, with threshold `2e-11`.
Native ratio checks supplement this scale for tiny negative log probabilities.
Reference configurations are
`-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace`
and `-O3 -fcheck=bounds`. The full `fpm test` suite passed with exit code 0
in 2068.864 seconds (34 minutes 29 seconds), including all 1,093 new native
checks. No numerical changes followed that successful run. The unattended
console-input check, source ASCII/line lengths, Python syntax, inventory/API
counts, changed local links, and `git diff --check` passed.

```cmd
python tools\generate_censored_reference.py
python validation\check_censored.py --module-dir build\gfortran_2806A042DAA8BBF7 --library build\gfortran_3C159EB802C22FB9\scifort\libscifort.a
fpm test
python tools\check_test_inputs.py
```

The full suite was run to check the shared gamma-product refinement and
integration of the new overloads with the existing native API. Subsequent
documentation-only edits did not require another numerical run.

## Limits

The beta-prime censored API requires a finite rounded shape sum `a+b`;
otherwise it returns quiet NaN, including for empty samples. Incomplete-gamma
and beta kernel accuracy limits carry through to these transformations.
At very large shapes, coordinate rounding limits near-mode accuracy in
sharply concentrated distributions. Subnormal results retain binary64
resolution. Beta-prime free-location fits with first shape below one can
have unbounded likelihoods at exact observations. A fully fixed singular
or zero-density objective returns `fit_status_no_finite_objective`; bounded
optimizer success does not certify a global MLE. No C ABI or wrappers are added.
