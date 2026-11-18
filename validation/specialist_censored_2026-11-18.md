# Censored noncentral Student t and Landau validation

The [API contract](../docs/censored_fitting.md#noncentral-student-t-and-landau-censoring)
specifies shape limits, affine behavior, tail and interval methods.

## Environment and commands

GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
Independent references use Python 3.13, mpmath 1.3.0 and SciPy 1.15.3.
These are validation dependencies only.

The native test `test/test_censored_specialist.f90` is compiled at `-O0`
and `-O3`, each with `-g -fcheck=all -ffpe-trap=invalid,zero,overflow
-fbacktrace`. Sources in order are `src/empirical_distribution.f90`,
`src/censored.f90`, `src/fit.f90`, `src/stats.f90`, and the test; the executable
links the existing FPM archive with its matching module include directory.

Exact native commands, after creating the two output directories, are:

```bat
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Jbuild\specialist-native -Ibuild\specialist-native -Ibuild\gfortran_2806A042DAA8BBF7 src\empirical_distribution.f90 src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored_specialist.f90 build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\specialist-native\test.exe
build\specialist-native\test.exe
gfortran -O3 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Jbuild\specialist-optimized -Ibuild\specialist-optimized -Ibuild\gfortran_2806A042DAA8BBF7 src\empirical_distribution.f90 src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored_specialist.f90 build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\specialist-optimized\test.exe
build\specialist-optimized\test.exe
```

The independent checker's CMD command is:

```bat
python validation\check_censored_specialist.py
python validation\check_censored_specialist.py --reuse
```

It compiles the input-driven driver at `-O0` and `-O3`, both with
`-g -fcheck=all -ffpe-trap=invalid,zero,overflow`. `--reuse` reuses cached
references under ignored `build/specialist-reference`. Its compiler source
list is `src/censored.f90` and `validation/censored_reference_driver.f90`,
linked against the matching existing FPM archive. No test requires keyboard
input under FPM.

The full-suite command is:

```bat
fpm test
```

## Independent checks

100-digit mpmath calculations integrate positive normal moments for the
noncentral-t density. They integrate the normal numerator against gamma
CDFs for tails, reversing the conditioning used by the Fortran kernel.
Subnormal noncentral-t coordinates use the analytic first-order CDF limit
and zero-coordinate density; discarded corrections are negligible compared
with the binary64 target. The Landau oracle uses adaptive complex contour inversion, characteristic
function integration, and right-tail Laplace integration. No SciPy/Boost
implementation code is used. Ordinary seeded cases and scale fits use
public SciPy calls as additional checks.

Checks cover all observation types, ordinary and narrow/adjacent/subnormal
intervals, support and probability limits, invalid inputs, empty data,
NaN and infinities, shape bounds, extreme and tiny complementary logs,
large location-to-scale ratios, scales 1e-200/1/1e200, fully and partially overflowed
standardized interval endpoints, related central-t identities, retained real-array
fit overloads, bounded scale fitting, and likelihood/NNLF signs.

Some public SciPy 1.15.3 noncentral-t tails lose relative accuracy for
noncentralities near +/-10. Those cases use defining-integral multiprecision
references, including both reflections, rather than SciPy as the sole oracle.
The existing scalar Landau left-tail switch also has larger error near -2.7;
the new censored contour has separate reference checks there.

## Results

- Checked and optimized native runs: 405 assertions each, zero failures.
- Both independent driver profiles pass 279 multiprecision cases and 180
  seeded ordinary SciPy comparisons. Maximum measured error is
  4.667377595524158e-13, using scaled absolute log error and relative error
  for tiny complementary logs, against a 3e-11 multiprecision target.
  Ordinary SciPy comparisons use 3e-10 scaled log tolerance. Native checks
  against subtracted scalar CDFs allow 2e-9 because that reference loses
  precision in narrow intervals; independent multiprecision interval checks
  retain the tighter 3e-11 target.
- Independent scale-fit references: noncentral t with df=4, nc=1 gives
  scale=0.31430866732103147 and NNLF=9.610148768581032. Landau's constrained
  optimum is the lower bound scale=0.3, NNLF=12.241979042425669. Native
  regressions check both fits. SciPy's interior bounded scalar optimizer
  approaches the Landau boundary from above.
- An additional 32-interval stress scan covers near-adjacent and ordinary
  widths in remote, central, and transition regimes. All log masses are
  nonpositive and none are NaN.
- Final full `fpm test`: exit 0, elapsed 1671.974 seconds
  (27 minutes 52.0 seconds). It includes all 405 new native
  assertions on the final source state, including the 22 overflowed-interval
  regressions. The earlier interrupted run is not counted as validation.

The small-error target applies to log values. It does not guarantee small
relative probability error for very large negative log values, or strict
last-bit monotonicity. Infinite df and wider noncentral-t shapes remain
unsupported by the new censored method. C ABI and wrapper work is deferred.

No other compiler/platform combination was exercised for this batch.
