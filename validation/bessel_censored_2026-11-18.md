# Bessel-family censored likelihood validation

The native additions cover generalized inverse Gaussian (`geninvgauss`),
normal-inverse Gaussian (`norminvgauss`), and generalized hyperbolic
(`genhyperbolic`). See [the numerical contract](../docs/censored_fitting.md#generalized-inverse-gaussian-and-hyperbolic-censoring)
for shape limits and remote-tail approximations.

## Compiler and commands

GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
The native focused test is `test/test_censored_bessel.f90`.
Checked and optimized manual builds use respectively `-O0` and `-O3`, both
with `-g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace`.
Sources are `src/special/log_gamma.f90`, `src/empirical_distribution.f90`,
`src/censored.f90`, `src/fit.f90`, `src/stats.f90`, and the new test,
linked against the existing FPM `libscifort.a` with the matching module
include directory. Separate build directories hold each profile's modules.

The independent checker builds its driver at `-O0 -g -fcheck=all
-ffpe-trap=invalid,zero,overflow` and `-O3 -fcheck=bounds`.
Its CMD invocation in this checkout is:

```bat
python validation\check_censored_bessel.py
```

It used mpmath 1.3.0 and SciPy 1.15.3 and requires those packages for validation only; neither becomes a library
or FPM test dependency. `--reuse` reuses the generated reference fixture in
ignored `build/bessel-reference`.

The required full-suite command is:

```bat
fpm test
```

## Coverage and limitations

Native checks cover scalar density/CDF/SF agreement, monotonicity, quantile
round trips, ordinary/adjacent/affine intervals, empty and invalid input,
infinite interval conventions, support endpoints, bounded censored fits,
preservation of array-fit overloads, GIG reciprocal symmetry, GH reflection,
and the exact NIG/GH identity. Subnormal intervals and large-location
intervals are included. The GIG scale fit is checked for convergence and
improvement over its starting objective. A concentrated lower-tail regression
uses an independently evaluated 65-digit reference.

Independent references use high-precision mathematical PDFs and adaptive
mpmath quadrature. Hyperbolic ordinary tail references instead integrate
the independent normal variance-mean mixture over its generalized inverse
Gaussian mixing density. Remote tails use a scaled density integral;
working precision increases with endpoint magnitude. Narrow intervals
retain exact binary64 endpoints as multiprecision values. A separate composite
Simpson rule checks widths below 0.001. At widths below
1e-12, the midpoint density approximation has error well below the selected
log-error tolerance for these fixtures.

The stress cases include shape boundaries, scales 1e-200/1/1e200,
positive/negative location-to-scale ratios 1e8/1e12, adjacent endpoints,
very small complements, and standardized magnitudes through 1e100.
Ordinary seeded density comparisons use SciPy independently. Generic SciPy
quadrature is not an oracle for deep tails: some finite results saturate
near log(machine epsilon) or have substantial relative error.

The censored methods deliberately expose narrower shape ranges than the
existing scalar distributions. Remote endpoint approximations are assessed
using scaled log errors, not relative probability error or correct rounding.
Hyperbolic standardized-coordinate overflow uses infinite-coordinate limits,
even if scaling its exponential decay rate could yield a finite log tail;
that extension is not covered here. No C ABI/wrapper expansion or new
dependency is included.

## Results

- Checked and optimized native runs: 263 assertions each, zero failures.
- Checked and optimized independent-driver runs: 120 multiprecision cases
  and 60 seeded SciPy density comparisons each, all passed.
- Maximum reported reference error: 4.9093193123542804e-12. The measure is
  absolute log error divided by max(1,abs(reference log)); tiny complementary
  logs instead use relative error. The worst cases were large-location
  adjacent intervals.
- Regenerating references with the final checker reproduced passing results
  in both profiles. Maximum relative error in tiny complementary logs was
  1.3211653993039363e-13.
- Full `fpm test`: passed, exit code 0, in 1301.4249148999806 seconds
  (21 minutes 41.4 seconds), including all 263 new native checks. This run
  covers the final numerical implementation and tests. Subsequent edits only
  recorded validation results in documentation; the suite was not repeated.
- Final checks passed for source ASCII/line limits, Python syntax, local
  documentation links, whitespace, and the 101-supported/8-remaining inventory.
  Build products and generated reference fixtures remain under ignored `build/`.



## Focused native build commands

The following are the commands used, with CMD continuations and relative
paths for the same checkout. The FPM cache names identify this compiler/build
configuration and can differ in another checkout.

```bat
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace ^
  -Jbuild\bessel-native -Ibuild\bessel-native ^
  -Ibuild\gfortran_2806A042DAA8BBF7 ^
  src\special\log_gamma.f90 src\empirical_distribution.f90 src\censored.f90 ^
  src\fit.f90 src\stats.f90 test\test_censored_bessel.f90 ^
  build\gfortran_3C159EB802C22FB9\scifort\libscifort.a ^
  -o build\bessel-native\test.exe
build\bessel-native\test.exe

gfortran -O3 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace ^
  -Jbuild\bessel-optimized -Ibuild\bessel-optimized ^
  -Ibuild\gfortran_2806A042DAA8BBF7 ^
  src\special\log_gamma.f90 src\empirical_distribution.f90 src\censored.f90 ^
  src\fit.f90 src\stats.f90 test\test_censored_bessel.f90 ^
  build\gfortran_3C159EB802C22FB9\scifort\libscifort.a ^
  -o build\bessel-optimized\test.exe
build\bessel-optimized\test.exe
```
