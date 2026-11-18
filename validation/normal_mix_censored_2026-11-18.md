# Normal-mixture censoring validation, 2026-11-18

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
Independent references: SciPy 1.15.3 and mpmath 1.3.0.

## Native regression

The [native test](../test/test_censored_normal_mix.f90) passes all 2,156
checks at both checked `-O0` and optimized `-O3`. It checks exact,
left/right, and interval-censored likelihoods, NNLFs, CDF monotonicity,
CDF/PPF and SF/ISF round trips, ordinary and adjacent-value intervals,
affine scales 1e-200/1/1e200, invalid and empty parameters/data,
NaN/infinity behavior, support boundaries, equal endpoints, and one-sided
infinite intervals. Fixed and one-parameter fits, improving all-free
candidates, and existing real-array overloads are exercised.

Additional checks cover skew-normal's normal limit and its closed-form
origin CDF for integer shapes -10 through 10; exponentially modified
normal tails at z=1000 and adjacent intervals there; tiny lower Gaussian
tails; finite log densities near binary64's largest Gaussian coordinate;
overflowed normal coordinates; double Pareto-lognormal reciprocal tails,
u/outer-scale confounding, subnormal inputs/intervals, and positive physical
distances that would overflow before scaling.

CMD commands, with existing output directories and FPM dependencies:

```cmd
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Jbuild\normal_mix-native -Ibuild\normal_mix-native -Ibuild\gfortran_2806A042DAA8BBF7 src\special\log_gamma.f90 src\empirical_distribution.f90 src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored_normal_mix.f90 build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\normal_mix-native\test.exe
build\normal_mix-native\test.exe
```

The optimized command substitutes `-O3` and paths under
`build\normal_mix-optimized`, retaining all checks and traps. Changed
modules are compiled before linking the existing archive. The error
measure is abs(actual-reference)/max(1,abs(reference)). Ordinary scalar
and adjacent/affine checks allow 3e-12, ordinary intervals 3e-11, quantile
round trips 5e-11, and identities 2e-12 or 3e-12. Extreme log tails,
subnormal cases, and origin-CDF checks allow 2e-13 or 3e-13. Tiny negative log-CDF/log-SF
checks divide by the nonzero reference and require relative error <=3e-11. Classifications
and exact support values use zero tolerance. Scaled log errors do not
establish uniform relative probability errors after exponentiation.

## Independent references

The [checker](check_censored_normal_mix.py) passes all 207 comparisons
in checked and optimized builds: 60 seeded mixed samples and 147 stresses.
Maximum scaled log error is `5.74777e-14`; maximum relative error for
the six tiny-log checks is `1.23679e-13`. It uses 100 decimal digits and
seed 2026111801. EMG references use direct erfc/exponential formulas;
skew-normal uses full semi-infinite multiprecision density quadrature,
independent of the Fortran finite GL16 rule; double Pareto-lognormal uses
the density and its integrated signed CDF identities rather than positive
EMG mixtures.
Finite interval references factor a midpoint density and physical width;
steep intervals use the better high-precision tail difference. Density-only
checks avoid unnecessary CDF integrations. Finite comparisons require
scaled error <=3e-11; six tiny-complement cases instead require relative
error <=3e-11 on the nonzero log probability.

```cmd
python validation\check_censored_normal_mix.py --module-dir build\gfortran_2806A042DAA8BBF7 --library build\gfortran_3C159EB802C22FB9\scifort\libscifort.a
```

Checked driver flags are `-O0 -g -fcheck=all
-ffpe-trap=invalid,zero,overflow`; optimized flags are `-O3 -fcheck=bounds`.
The optional `--reuse` flag uses ignored reference files under `build/`.
The checker covers seeded mixed samples, EMG k=0.01/0.1/1/100,
skew-normal shapes -10/-1/0/1/10, coordinates through magnitude 1000,
adjacent intervals, double Pareto-lognormal shape limits and inputs
1e-300/1/1e300, affine scales 1e-300/1e300, and subnormal positive
Pareto coordinates before scaling. The six relative-error cases check
EMG CDF/SF at 700/-35, skew-normal CDF/SF at 35/-35 with shapes 2/-2,
and Pareto CDF/SF at 1e100/1e-100, with log magnitudes near 1e-300.

Independent SciPy bounded fits use exact [0.2,0.4,0.6], left [0.1],
right [0.8], and intervals [0.3,0.35] and [0.5,0.55]. Only scale varies
in [0.3,3], with loc=0:

| Family | Fixed shapes | Fitted scale | NNLF |
| --- | --- | ---: | ---: |
| exponnorm | k=1.7 | 0.3000000072924005 | 10.081147658879386 |
| skewnorm | a=1 | 0.5118115365121739 | 10.105222635536382 |
| dpareto_lognorm | u=0.1, s=0.7, a=2, b=3 | 0.3000000072924005 | 10.008137885354882 |

The first and third optima are on the lower bound; SciPy's bounded scalar
optimizer returns nearby interior approximations. Native comparisons
allow 3e-5 scaled parameter error and 2e-8 scaled NNLF error. All-free
checks require finite improving candidates, not global optima.

## Full suite and limitations

Shared censoring/fit/statistics integration requires the full suite in
addition to the checked and optimized focused validation:

```cmd
fpm test
```

The full suite passed with exit code 0 in 2057.629 seconds (34 minutes
18 seconds), including all 2,156 new checks and the existing censored,
scalar-distribution, fitting, C API, and remaining repository tests.
Existing array-temporary warnings appeared in unrelated FFT, interpolation,
and spatial tests and did not cause failures. No executable code changed
after this successful full run began; subsequent edits completed only
validation metadata and documentation.

The [API guide](../docs/censored_fitting.md) documents shape ranges,
coordinate-overflow behavior, log-range limits, steep-interval sensitivity
to rounding, and u/scale nonidentifiability. Exact double Pareto-lognormal
support endpoints retain the existing scalar API's zero-density convention.
Wider shape ranges and improved EMG standardization for unrepresentable
coordinates remain numerical coverage gaps. Native coverage is the current
scope; C ABI expansion remains deferred.
