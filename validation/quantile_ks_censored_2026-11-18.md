# Tukey lambda, KSone and fixed-branch von Mises validation

These additions provide native censored likelihoods, NNLFs and bounded-fit
overloads. Fixed-branch von Mises censoring is a **SciFort extension**;
periodic circular CDF values are not used as ordinary probabilities.
See [the API contract](../docs/censored_fitting.md#tukey-lambda-censoring)
for shape, affine-support, endpoint-density and circular-branch limits.

## Environment and commands

GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
Independent references use Python 3.13, mpmath 1.3.0 and SciPy 1.15.3;
these are validation dependencies only.

The native test `test/test_censored_quantile_ks.f90` was manually compiled
at `-O0` and `-O3`, both with `-g -fcheck=all
-ffpe-trap=invalid,zero,overflow -fbacktrace`. Sources, in order, were
`src/math.f90`, `src/special/log_gamma.f90`, `src/empirical_distribution.f90`,
`src/censored.f90`, `src/fit.f90`, `src/stats.f90`, and the test, linked
against the existing FPM archive with matching module include directories.

The independent checker's CMD command is:

```bat
python validation\check_censored_quantile_ks.py
```

It compiles its input-driven validation driver at `-O0 -g -fcheck=all
-ffpe-trap=invalid,zero,overflow` and `-O3 -fcheck=bounds`.
`--reuse` reuses reference fixtures under ignored
`build/quantile_ks-reference`. The driver is outside FPM's unattended tests
and does not add keyboard input to `fpm test`.

The required full-suite command is:

```bat
fpm test
```

## Numerical coverage

Native checks cover density and complementary tails, monotonicity,
quantile round trips, ordinary/adjacent/affine intervals, endpoints,
NaN/infinite/invalid parameters and data, empty-input behavior, fixed and
free fits, integer-constrained KSone fitting against enumeration,
existing real-array fit overloads, and likelihood/NNLF identities.
Special-case checks include Tukey lambda zero/logistic, lambda 1 and 2
uniform laws, KSone n=1/uniform, and fixed-branch von Mises/line agreement.
The two sides of the first KSone density knot and subnormal intervals
have explicit regressions. SciPy supplies independent Tukey-lambda and
line von Mises scale-fit references.

Independent references use 100-digit mpmath calculations: direct inversion
of the mathematical Tukey quantile, exact finite Smirnov sums and their
derivatives, and scaled positive von Mises density integrals. Adjacent
intervals use high-precision midpoint densities only when smooth relative
variation is negligible; wider short intervals use a separate Simpson
rule. Intervals split at the first KSone density jump; intervals reaching
a support endpoint use the corresponding exact high-precision tail.

The grid includes shapes at numerical limits, positive/negative lambda
near zero, n=1/2/3/10/100, kappa=0/1/100, standardized magnitudes through
1e300, scales 1e-200/1/1e200, positive/negative location-to-scale ratios
1e12 for unbounded Tukey lambda, adjacent finite endpoints, subnormal
widths, and tiny complementary logs. Binary64 computed affine endpoints
follow the documented support convention.

Some SciPy values are unsuitable as sole references: the first KSone
knot can use a different representative density, small nonzero Tukey
lambda can be approximated as zero, and concentrated von Mises tails
can have substantial error. These cases use mathematical multiprecision
references. Seeded direct SciPy comparisons cover ordinary values away
from these regimes.

## Results

- Checked and optimized native runs: 468 assertions each, zero failures.
- The gamma/helper regression also passes at `-O0` and `-O3` with
  `-g -fcheck=all -fbacktrace`, including 104 new exact small/subnormal
  helper checks. It links `src/math.f90`, `test/test_checks.f90`,
  `test/gamma_reference.f90`, and `test/test_gamma.f90` against the FPM archive.
- Both independent driver profiles passed all 355 multiprecision cases and
  210 seeded ordinary SciPy comparisons.
- Maximum scaled log error was 4.747929967498059e-12, using
  abs(actual-reference)/max(1,abs(reference)). Maximum relative error in
  tiny complementary logs was 1.1368683772161603e-13.
- A 2376-point adjacent-neighbor scan around representation switches and
  KS endpoint regions found at most 1.936228954946273e-13 downward log-CDF
  variation near the first KSone knot. This is within the numerical log
  tolerance; strict last-bit monotonicity is not guaranteed by the inherited
  finite-sum tail helper.
- Full `fpm test`: passed, exit code 0, in 2039.7613854999654 seconds
  (33 minutes 59.8 seconds), including the 468 new censored checks and
  the gamma/helper regressions. This run covers the final numerical code
  and tests. Subsequent edits only recorded validation results in documentation.
- Final source ASCII/line-limit, Python syntax, local-link, whitespace,
  and 104-supported/5-remaining inventory checks passed. Generated fixtures,
  executables, compiler modules and logs remain under ignored `build/`.


The gamma executable was also attempted with invalid/zero/overflow traps.
Its existing NaN inverse-argument test stops in the unchanged
`scifort_incomplete_gamma::valid_inverse_arguments` ordered comparison.
The gamma regression therefore ran with bounds/runtime checks and without
those traps. The new censored test and independent checked driver both
pass with those traps enabled.


## Limits

Tukey censoring currently supports abs(lam)<=5, with zero or abs(lam)>=1e-6;
positive shapes also have the documented affine restrictions. KSone
requires integer 1<=n<=100 and the documented finite affine support.
Fixed-branch von Mises inherits the line method's limits, including
0<=kappa<=100. Circular crossing arcs remain unsupported. The shared `log1p_safe` fix preserves tiny inputs exactly for
abs(x)<=epsilon/4, avoiding halving underflow for odd subnormal values.
These methods
do not establish complete SciPy parity or expand the C ABI/wrappers.

## Reproducing the focused native builds

These CMD commands show the flags and source order used in this checkout.
The compiler-specific FPM cache paths can differ in another checkout.
Separate output directories must already exist.

```bat
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace ^
  -Jbuild\quantile_ks-native -Ibuild\quantile_ks-native ^
  -Ibuild\gfortran_2806A042DAA8BBF7 ^
  src\math.f90 src\special\log_gamma.f90 src\empirical_distribution.f90 ^
  src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored_quantile_ks.f90 ^
  build\gfortran_3C159EB802C22FB9\scifort\libscifort.a ^
  -o build\quantile_ks-native\test.exe
build\quantile_ks-native\test.exe

gfortran -O3 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace ^
  -Jbuild\quantile_ks-optimized -Ibuild\quantile_ks-optimized ^
  -Ibuild\gfortran_2806A042DAA8BBF7 ^
  src\math.f90 src\special\log_gamma.f90 src\empirical_distribution.f90 ^
  src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored_quantile_ks.f90 ^
  build\gfortran_3C159EB802C22FB9\scifort\libscifort.a ^
  -o build\quantile_ks-optimized\test.exe
build\quantile_ks-optimized\test.exe
```

The two gamma/helper builds used these commands, followed by executing each
output file:

```bat
gfortran -O0 -g -fcheck=all -fbacktrace ^
  -Jbuild\quantile_ks-native -Ibuild\quantile_ks-native ^
  -Ibuild\gfortran_2806A042DAA8BBF7 ^
  src\math.f90 test\test_checks.f90 test\gamma_reference.f90 test\test_gamma.f90 ^
  build\gfortran_3C159EB802C22FB9\scifort\libscifort.a ^
  -o build\quantile_ks-native\gamma.exe
build\quantile_ks-native\gamma.exe

gfortran -O3 -g -fcheck=all -fbacktrace ^
  -Jbuild\quantile_ks-optimized -Ibuild\quantile_ks-optimized ^
  -Ibuild\gfortran_2806A042DAA8BBF7 ^
  src\math.f90 test\test_checks.f90 test\gamma_reference.f90 test\test_gamma.f90 ^
  build\gfortran_3C159EB802C22FB9\scifort\libscifort.a ^
  -o build\quantile_ks-optimized\gamma.exe
build\quantile_ks-optimized\gamma.exe
```
