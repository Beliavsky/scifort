# Inverse-Gaussian family censored validation, 2026-11-17

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
Independent references: SciPy 1.15.3 and mpmath 1.3.0.

## Native checks

The [native regression](../test/test_censored_inverse_gaussian.f90) covers
inverse Gaussian, Wald, and reciprocal inverse Gaussian: scalar density/CDF/SF
agreement, monotonicity, CDF/PPF and SF/ISF round trips, ordinary and adjacent
intervals, scales 1e-300/1/1e300, invalid inputs, support endpoints, empty
data, infinite one-sided intervals, equal interval endpoints, fixed and
all-free censored fits, original array-fit compatibility, NNLF, Wald/shape-one
and reciprocal-tail identities, huge/tiny tails, subnormal shapes, and a
tiny-shape point one ULP from the center. A relative check retains a
complementary log probability of order 1e-150.

All 2,184 checks pass at both `-O0` and `-O3`, with all runtime checks and
invalid/zero/overflow traps enabled. One-parameter native fits agree with
the independent SciPy optima below within 3e-5 scaled parameter error and
2e-8 scaled NNLF error. The steep rising/falling adjacent-interval
regressions retain finite log masses of order -1e199 when rounded endpoint
log tails coincide.

Checked command (CMD, with existing output directory and FPM dependencies):

```cmd
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Jbuild\inverse_gaussian-native -Ibuild\inverse_gaussian-native -Ibuild\gfortran_2806A042DAA8BBF7 src\special\log_gamma.f90 src\empirical_distribution.f90 src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored_inverse_gaussian.f90 build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\inverse_gaussian-native\test.exe
build\inverse_gaussian-native\test.exe
```

The optimized run substitutes `-O3` and paths under
`build\inverse_gaussian-optimized`, retaining all runtime and trap flags.
Changed modules are compiled explicitly before the existing archive is linked.

Most scalar and round-trip log comparisons allow 2e-12/3e-12 scaled error;
ordinary interval differences allow 2e-11; adjacent/affine logs allow 3e-13.
The tiny-shape near-center stress allows 6e-13 relative error after logarithmic
rescaling of the very large normal penalty. Tiny log probabilities are compared
relative to their own magnitude. Invalid/infinite classifications and exact
special-case identities use no tolerance.

## Independent references

All 427 comparisons pass in both checked and optimized builds: 300 mixed
samples and 127 stresses. Maximum scaled error is `2.362e-12`.
The extra center stresses include shapes 1e-100/1e-300, one-ULP offsets,
affine center cancellation, and the minimum subnormal shape.

```cmd
python validation\check_censored_inverse_gaussian.py --module-dir build\gfortran_2806A042DAA8BBF7 --library build\gfortran_3C159EB802C22FB9\scifort\libscifort.a
```

The checker uses independent normal/erfc formulas and rescaled density
integration with mpmath at 100 digits for seeded mixed samples and 750 digits
for extreme cases. Seed: 20261117. Shape stresses range from 1e-100 to 1e300,
tail points and scales from 1e-300 to 1e300, with adjacent intervals and nearly
coincident normal coordinates. It compares the full mixed log likelihood,
requiring `abs(actual-reference)/max(1,abs(reference)) <= 2e-11`.

Narrow smooth interval references integrate a rescaled density with mpmath's
Gauss-Legendre rule, maximum refinement degree 3, under a 0.001 log-density
variation bound. This independently evaluates the integral with more nodes
than the native three-point rule. At the minimum subnormal center, mpmath
1.3.0's erfc series selector overflows while converting the squared argument
to a float; the oracle instead evaluates the defining Mills integral
`integral_0^infinity exp(-u-u*u/(2*B*B)) du / B` at 100 digits. All other
stress normal-tail arithmetic remains at 750 digits.

Checked driver flags are `-O0 -g -fcheck=all
-ffpe-trap=invalid,zero,overflow`; optimized flags are `-O3 -fcheck=bounds`.
The checker also prints independent one-parameter SciPy bounded-fit benchmarks
for the mixed sample: exact `[0.2,0.4,0.6]`, left `[0.1]`, right `[0.7]`,
intervals `[0.25,0.35]` and `[0.5,0.55]`. For inverse Gaussian/reciprocal,
loc=0 and scale=1.3 are fixed, and shape varies on `[0.3,5]`; for Wald,
loc=0 is fixed and scale varies on `[0.3,5]`.

| Family | Free parameter | NNLF |
| --- | ---: | ---: |
| Inverse Gaussian shape | 0.3349859947937545 | 10.735705684804694 |
| Wald scale | 0.47823955115344424 | 9.221866539671925 |
| Reciprocal inverse Gaussian shape | 4.999999874519365 | 11.23530622645976 |

The reciprocal optimum lies at the selected upper bound; this is a
bounded-sample benchmark, not evidence of an unconstrained maximum.

## Scope and limits

The shared censoring and fitting integration changed, so the full repository
suite is required in addition to focused checked and optimized runs:

```cmd
fpm test
```

The full suite passed with exit code 0 in 2165.310 seconds (36 minutes
5 seconds), including all 2,184 new checks and the existing scalar,
censored, fitting, and integration regressions. Existing array-temporary
warnings appeared in scattered interpolation and signal tests; they did
not cause failures. No numerical code changed after this run began.

The implementation adds private callbacks and fit dispatch without changing
the scalar distribution routines, optimizer, C ABI, dependencies, or global
state. Extremely concentrated affine models can be ill-conditioned relative
to input rounding. Log probabilities beyond binary64 range return negative
infinity, and finite logarithms can represent masses that would underflow
if exponentiated. These are documented separately from named-family coverage.
