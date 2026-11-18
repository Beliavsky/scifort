# Kolmogorov, angular, and resonance censoring validation, 2026-11-18

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
Independent references: SciPy 1.15.3 and mpmath 1.3.0.

## Native regression

The [native test](../test/test_censored_three_tail.f90) passes all 2,126
checks at checked -O0 and optimized -O3. It checks exact/left/right/interval
contributions, monotonicity, CDF/PPF and SF/ISF round trips, ordinary and
adjacent intervals, affine scales 1e-200/1/1e200, invalid and empty inputs,
NaN/infinity behavior, support boundaries, equal endpoints, and one-sided
infinite intervals. Fixed and scale-only fits, improving all-free
candidates, and existing real-array overloads are exercised.

Extra regressions cover extreme Kolmogorov log tails, tiny complementary
logs, subnormal angular/resonance intervals, angular support clipping and
midpoint probability, invalid crossing arcs, finite endpoint densities,
NNLF signs, shape/affine bounds, narrow intervals at large location offsets,
and physical distances that would overflow
before scaling. Independent scale-fit parameters allow scaled error 3e-5
and NNLFs 2e-8. Ordinary density/tail/adjacent checks allow 3e-12, wider
intervals 3e-11, and quantile round trips 5e-11. Extreme-tail/subnormal
regressions allow 3e-13; tiny-log regressions compare relative error to
3e-12. Exact classifications use zero tolerance.

CMD commands, with existing output directories and FPM dependencies:

```cmd
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Jbuild\three_tail-native -Ibuild\three_tail-native -Ibuild\gfortran_2806A042DAA8BBF7 src\special\log_gamma.f90 src\empirical_distribution.f90 src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored_three_tail.f90 build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\three_tail-native\test.exe
build\three_tail-native\test.exe
```

The optimized command substitutes -O3 and paths under
`build\three_tail-optimized`, retaining runtime checks and traps.
Modules are compiled before linking the existing archive.

## Independent references

The [checker](check_censored_three_tail.py) passes all 204 high-precision
cases in both profiles: 30 seeded mixed observations and 174 stress cases.
It uses 110 decimal digits and seed 2026111804. Maximum scaled log error
is 7.59556e-13; maximum relative tiny-log error is 1.72529e-13.
The error measure for ordinary finite logs is
abs(actual-reference)/max(1,abs(reference)), limited to 3e-11; five tiny
complementary-log cases instead use relative error <=3e-11. Infinite
reference logs require matching classification. These measures do not
establish uniform relative probability errors after exponentiation.

References use 24-term multiprecision theta series and differentiated
Kolmogorov densities, exact wrapped-Cauchy half-angle identities with the
computed binary64 angular endpoint convention, and direct positive
Breit-Wigner density quadrature. Upper resonance tails use reciprocal
quadrature, independently of the Fortran coefficient expansion. Ordinary
finite intervals integrate relative densities after factoring midpoint
density and physical width; steep Kolmogorov intervals use high-precision
log-tail differences.

Coverage includes Kolmogorov coordinates 1e-150 through 1e150, both sides
of the 0.82 series switch, angular concentrations 1e-300/.001/.4/.99,
upper angular endpoints and minimum-positive subnormal intervals,
resonance shapes .1/1/100, both sides of the reciprocal-series cutoff,
resonance coordinates through 1e300, adjacent intervals, and affine scales
1e-300/1e300. Another 54 interval cases cover positive/negative
location-to-scale ratios 1e8/1e10/1e12 with scales 1e-200/1/1e200.
Five relative-error cases cover both Kolmogorov tails,
small-angle wrapped-Cauchy survival, and both resonance tails.

```cmd
python validation\check_censored_three_tail.py --module-dir build\gfortran_2806A042DAA8BBF7 --library build\gfortran_3C159EB802C22FB9\scifort\libscifort.a
```

Checked driver flags are -O0 -g -fcheck=all
-ffpe-trap=invalid,zero,overflow; optimized flags are -O3 -fcheck=bounds.
The optional --reuse flag reads ignored high-precision fixtures under build/.
A second seeded stress check (2026111805) samples 200 shape/coordinate sets
per family. All 1,058 finite, nontrivial SciPy comparisons pass in both
profiles. Discrepancies above scaled error 3e-10 are independently
adjudicated against 110-digit formulas at tolerance 3e-11; 38 SciPy
1.15.3 discrepancies were resolved this way, agreeing with SciFort.
This avoids treating finite but inaccurate reference-library tails as
unquestioned expected values.

Independent SciPy scale fits use exact [0.2,0.4,0.6], left [0.1], right
[0.8], intervals [0.3,0.35] and [0.5,0.55], loc=0, and scale in [0.3,3]:

| Family | Shape | Fitted scale | NNLF |
| --- | --- | ---: | ---: |
| kstwobign | none | 0.36696009974403854 | 28.498269735035024 |
| wrapcauchy | c=0.4 | 0.3000000072924005 | 12.550326038133834 |
| rel_breitwigner | rho=1 | 0.46747001313866393 | 9.354988140535767 |

The wrapped-Cauchy optimum is at the lower bound; SciPy's bounded solver
returns a nearby interior point. All-free checks require finite improving
candidates, not a proven global optimum.

## Full suite and limitations

The final full-suite command was:

```cmd
fpm test
```

The full suite passed with exit code 0 in 2348.545 seconds (39 minutes
9 seconds), including all 2,126 new checks and existing censoring,
scalar-distribution, fitting, C API, and other repository tests. Existing
array-temporary warnings appeared in unrelated spatial and signal tests
and did not cause failures. No executable code changed after this final
successful run began; subsequent edits completed only validation metadata
and documentation.

The [API guide](../docs/censored_fitting.md) documents the fixed angular
branch, computed endpoint convention, affine/shape caps, and binary64 log
range. Mathematical series remainder bounds are in
[provenance](../CODE_PROVENANCE.md); they do not certify inherited
antiderivative or floating-point rounding errors. The ordinary resonance
tails reuse the existing BSD-derived antiderivative with its retained
notices. Crossing-boundary angular arcs, broader numerical ranges, and
C ABI/wrapper additions remain outside this batch.
