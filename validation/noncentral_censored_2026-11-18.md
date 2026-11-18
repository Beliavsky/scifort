# Noncentral censoring validation, 2026-11-18

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
Independent references: SciPy 1.15.3 and mpmath 1.3.0.

## Native regression

The [native test](../test/test_censored_noncentral.f90) passes all 2,127 checks
at both checked -O0 and optimized -O3. It exercises ordinary
exact/left/right/interval contributions, monotonicity, CDF/PPF and SF/ISF
round trips, adjacent intervals, affine scales 1e-200/1/1e200, invalid and
empty parameters/data, NaN/infinity conventions, exact support endpoints,
equal endpoints, and one-sided infinite intervals. It checks the three
central-limit identities, NNLFs, fixed/scale-only fits, improving all-free
fit candidates, existing real-array overloads, a rejected unconverged
series, and overflowed physical distances. Subnormal F density/CDF/interval
regressions retain physical widths. Multiprecision Rice tails and an F
density guard against known discrepancies with the installed SciPy.

CMD commands, with existing output directories and FPM dependencies:

```cmd
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Jbuild\noncentral-native -Ibuild\noncentral-native -Ibuild\gfortran_2806A042DAA8BBF7 src\special\log_gamma.f90 src\empirical_distribution.f90 src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored_noncentral.f90 build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\noncentral-native\test.exe
build\noncentral-native\test.exe
```

The optimized command substitutes -O3 and paths under
`build\noncentral-optimized`, retaining checks and traps. Modules are
compiled before linking the existing archive. Scaled log error means
abs(actual-reference)/max(1,abs(reference)); ordinary densities/tails and
adjacent intervals allow 3e-12, finite intervals 3e-11, quantile round trips
5e-11, and multiprecision endpoint/tail/density regressions 3e-13.
Tiny complementary logs use relative error <=3e-12. Exact classifications
use zero tolerance. These criteria do not establish uniform relative
probability error after exponentiation.

## Independent reference checks

Both profiles pass 144 high-precision cases (30 mixed and 114 stress
cases), with maximum scaled log error 9.80546e-15 and maximum relative
tiny-log error 2.27374e-13. The additional seeded comparison passes all
653 finite/nontrivial SciPy checks in both profiles, resolving five
SciPy discrepancies against multiprecision formulas.

The [checker](check_censored_noncentral.py) uses 110 decimal digits,
seed 2026111802, and 10 mixed samples per family by default. Density
references use Bessel I for Rice/chi-squared and hypergeometric 1F1 for F.
Tail references use multiprecision regularized gamma/beta mixtures;
relative-density finite quadrature supplies interval masses, independently
of the Fortran interval callbacks. Cases include shape limits, coordinates
1e-300 through 100, adjacent endpoints, affine scales 1e-300/1e300,
minimum positive subnormals, and six tiny complementary logs.

```cmd
python validation\check_censored_noncentral.py --module-dir build\gfortran_2806A042DAA8BBF7 --library build\gfortran_3C159EB802C22FB9\scifort\libscifort.a
```

Checked driver flags are -O0 -g -fcheck=all
-ffpe-trap=invalid,zero,overflow; optimized flags are -O3 -fcheck=bounds.
The optional --reuse flag reads ignored reference files under build/.
Finite references require scaled error <=3e-11; tiny complementary logs
instead require relative error <=3e-11 on the nonzero logarithm.

A second seeded check (2026111803) samples 200 shape/coordinate sets per
family, then checks finite nontrivial SciPy density/CDF/SF values. Cases
with discrepancies above 3e-10 are independently adjudicated with the
110-digit formulas and must agree with SciFort within 3e-11. This avoids
using known inaccurate SciPy 1.15.3 Rice tails and large-coordinate F
densities as unquestioned references. Both compiled profiles run this check.

Independent bounded fits use exact [0.2,0.4,0.6], left [0.1], right [0.8],
and intervals [0.3,0.35] and [0.5,0.55], loc=0, scale in [0.3,3]:

| Family | Shapes | SciPy fitted scale | NNLF |
| --- | --- | ---: | ---: |
| ncx2 | df=2, nc=1 | 0.3000000072924005 | 10.974656410900508 |
| rice | b=1 | 0.3000000072924005 | 10.055647712310694 |
| ncf | dfn=2, dfd=3, nc=1 | 0.3000000072924005 | 10.612152223011206 |

All optima lie at the lower bound; SciPy's bounded solver returns nearby
interior points. Native comparisons allow 3e-5 scaled parameter error and
2e-8 scaled NNLF error. All-free checks require finite improving candidates,
not convergence to a proven global optimum.

## Full suite and limitations

The final full-suite command was:

```cmd
fpm test
```

The full suite passed with exit code 0 in 1456.663 seconds (24 minutes
17 seconds), including all 2,127 new checks and the existing censoring,
scalar-distribution, fitting, C API, and other repository tests. Existing
array-temporary warnings appeared in unrelated spatial and signal tests
and did not cause failures. No executable code changed after this full
run began; subsequent edits completed only documentation and validation
metadata.

See the [API guide](../docs/censored_fitting.md) for finite shape ranges,
noncentral F transformed-scale requirements, support densities, and
convergence failures. The remainder bound controls series truncation only;
rounding and inherited special-function errors remain. Noncentral radial
tails can exceed the 10,001-term cap even with representable logarithms.
No C ABI or wrapper additions are included.
