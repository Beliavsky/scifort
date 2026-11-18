# ARGUS, Pearson III, and kappa4 censoring validation, 2026-11-17

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
Independent references: SciPy 1.15.3 and mpmath 1.3.0.

## Native regression

The [native test](../test/test_censored_gamma_kappa.f90) passes all 2,139
checks at both `-O0` and `-O3`, with runtime/bounds checks and
invalid/zero/overflow traps enabled. Coverage includes scalar log-density
and tail agreement, CDF monotonicity, CDF/PPF and SF/ISF round trips,
ordinary/adjacent intervals, affine scales 1e-200/1/1e200, invalid and empty
inputs, NaN/infinity behavior, one-sided infinite intervals, equal endpoints
as exact observations, support endpoints, fixed/one-parameter/all-free fits,
existing real-array fit compatibility, and NNLF.

Identities include Pearson's normal and shifted exponential limits and
reflection, and kappa4's GEV, logistic, exponential, and uniform cases.
Extreme regressions retain tiny ARGUS CDF/log-SF, underflowed kappa SF,
overflowed normal coordinates, and finite log densities/interval masses
near binary64's largest magnitude. Endpoint density classifications include
zero, finite, and infinite cases. Unsupported compressed kappa endpoints
and documented shape gaps are explicitly rejected.

CMD commands, with existing output directories and FPM dependencies:

```cmd
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Jbuild\gamma_kappa-native -Ibuild\gamma_kappa-native -Ibuild\gfortran_2806A042DAA8BBF7 src\special\log_gamma.f90 src\empirical_distribution.f90 src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored_gamma_kappa.f90 build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\gamma_kappa-native\test.exe
build\gamma_kappa-native\test.exe
```

The optimized command substitutes `-O3` and paths under
`build\gamma_kappa-optimized`, retaining all checks and traps. Changed
modules are explicitly compiled before the existing archive is linked.
Ordinary scalar, adjacent, and affine comparisons allow 3e-12 scaled error,
ordinary intervals 3e-11, and quantile round trips 5e-11. Most identities
allow 2e-12; the reflection identity allows 2e-13. The tiny ARGUS log-SF
is divided by its expected magnitude before comparison, with tolerance
2e-12. The error measure is abs(actual-reference)/max(1,abs(reference)).
Invalid/infinite classifications and exact support values use no tolerance.
Largest-magnitude checks compare with their leading asymptotic terms at
2e-13 scaled tolerance, not an absolute log-probability accuracy claim.

## Independent references

All 431 comparisons pass in checked and optimized builds: 300 seeded
mixed samples and 131 stresses. Maximum scaled error is `6.6229e-13`.

```cmd
python validation\check_censored_gamma_kappa.py --module-dir build\gfortran_2806A042DAA8BBF7 --library build\gfortran_3C159EB802C22FB9\scifort\libscifort.a
```

The [checker](check_censored_gamma_kappa.py) independently evaluates gamma
integrals, normal integrals, and elementary kappa antiderivatives using
mpmath. Ordinary samples and most stresses use 100 decimal digits;
ARGUS stresses use 300 digits to retain squared-coordinate increments
as small as 1e-200. Seed: 20261117. Finite comparisons require scaled
error <=2e-11. Narrow gamma and density integrals use independent mpmath
Gauss-Legendre quadrature with maximum refinement degree 3.

Stresses include ARGUS chi 1e-6/0.1/1/10/100 and points 1e-100/0.5/the
last binary64 point below one; Pearson skew -2000/-2/-0.002/0/1e-6/
0.002/2/2000, both exponential tails at magnitude 1e100, and adjacent
intervals; 32 supported kappa shape combinations from h in
[-10,-1,-1e-6,0,1e-6,1,10] and k in [-10,-1,0,1,10], each with an
adjacent interval; extreme affine scales 1e-300/1e300; and tiny uniform
and far logistic kappa limits. Unsupported compressed combinations are
excluded by the same explicit mathematical parameter contract as the API.

Kappa stress points use mathematical inverse formulas for positive h
and the interior point z=0.01 for nonpositive h. This avoids depending on
the installed SciPy inverse for negative h and avoids quantiles that round
onto an endpoint. Positive-h tiny CDF references use expm1 before taking
the log, retaining the small outer support factor independently.

Checked driver flags are `-O0 -g -fcheck=all
-ffpe-trap=invalid,zero,overflow`; optimized flags are `-O3 -fcheck=bounds`.

Independent SciPy bounded fits use exact [0.2,0.4,0.6], left [0.1],
right [0.8], and intervals [0.3,0.35] and [0.5,0.55]. Only scale varies
in [0.3,3], with loc=0:

| Family | Fixed shapes | Fitted scale | NNLF |
| --- | --- | ---: | ---: |
| argus | chi=2 | 0.8676833205360297 | 13.678940903065282 |
| pearson3 | skew=1 | 0.56885024604535128 | 13.368238839719634 |
| kappa4 | h=0.3, k=0.2 | 0.47467771562661354 | 10.672577136502522 |

Native comparisons allow 3e-5 scaled parameter error and 2e-8 scaled
NNLF error. All-free fits require finite improving candidates without
claiming global optima or identifiability.

## Full suite and limitations

The shared censoring/fit/statistics integration changed, so the full suite
is run in addition to focused checked and optimized validation:

```cmd
fpm test
```

The full suite passed with exit code 0 in 1203.263 seconds (20 minutes
3 seconds), including all 2,139 new checks and the existing censored,
scalar-distribution, fitting, C API, and integration regressions. Existing
array-temporary runtime warnings appeared and did not cause failures.
No numerical code changed after this successful run began; subsequent
edits only completed documentation and validation metadata.

The [API guide](../docs/censored_fitting.md) documents the numerical shape
limits: ARGUS chi in [1e-6,100]; Pearson's existing normal approximation
for abs(skew)<1.6e-5, otherwise 0.002<=abs(skew)<=2000; kappa absolute
shapes at most 10, each zero or magnitude at least 1e-6, with
k*log(h)<=log(10) when h>0. Wider shape ranges, the nonnormal small-skew
gap, and compressed kappa endpoints remain compatibility gaps.
Very close standardized/affine endpoints remain ill-conditioned relative
to binary64 input rounding. Scalar distribution implementations and C ABI
wrappers are unchanged.
