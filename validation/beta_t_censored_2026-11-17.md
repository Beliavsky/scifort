# Student t, F, and Jones-Faddy censored validation, 2026-11-17

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
References: SciPy 1.15.3 and mpmath 1.3.0.

## Native and dependent checks

All 1,846 checks in [the native test](../test/test_censored_beta_t.f90)
passed at `-O0` and `-O3`, each with
`-g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace`.
The affected [positive-transform test](../test/test_censored_positive_transform.f90)
also passed all 1,093 checks in both builds. F uses the beta-prime private
callbacks with an optional log rate; the dependent test covers the existing
four-parameter beta-prime path, inverse gamma, and Maxwell.

New checks cover scalar density/CDF/SF agreement, monotonicity, CDF/PPF
and SF/ISF round trips, adjacent intervals, scales 1e-200/1/1e200, invalid
inputs and parameter ranges, empty data, infinite interval endpoints,
equal-endpoint exact observations, support endpoints, minimum normal beta
shapes, maximum shape sums, fixed/shape-only fits, all-free candidate
improvement, NNLF, and original array-fit compatibility. Identities cover
Student t/Cauchy, the normal limit, symmetric Jones-Faddy/Student t, and
F/beta-prime. Relative tests retain tiny complementary log probabilities
of order 1e-280 for extremely unequal shapes.
Reflection tests cover far negative densities, tails, and ordinary/tiny-shape
intervals against positive intervals with the shape parameters exchanged.

Ordinary scalar logs allow 2e-12 scaled error, round trips 3e-12, most
adjacent/extreme logs 3e-13, and maximum-shape center comparisons 2e-11.
Measures use a unit absolute floor; tiny relative log comparisons first
divide by their expected magnitude. Exact classifications use no tolerance.

Checked command (CMD, with existing output directory and FPM dependencies):

```cmd
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Jbuild\beta_t-native -Ibuild\beta_t-native -Ibuild\gfortran_2806A042DAA8BBF7 src\special\log_gamma.f90 src\empirical_distribution.f90 src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored_beta_t.f90 build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\beta_t-native\test.exe
build\beta_t-native\test.exe
```

The optimized run substitutes `-O3` and paths under `build\beta_t-optimized`,
retaining all runtime/trap flags. The dependent runs substitute
`test\test_censored_positive_transform.f90` with separate module/executable
paths under `build\beta_t-dependent-native` and `build\beta_t-dependent-optimized`.
Changed modules are compiled explicitly before linking the existing archive.

## Independent references

All 407 comparisons passed in both checked and optimized builds: 300 mixed
samples and 107 stresses. Maximum scaled error is `1.931e-12`, in a
Jones-Faddy mixed sample. Stresses cover tiny/unequal shapes, minimum normal
beta parameters, adjacent intervals, 1e-300/1e300 scales and tail points,
underflowed scaled widths, center crossings, and maximum symmetric/asymmetric
shape sums. Sixteen maximum-shape comparisons use direct density integrals
or exact beta identities at 90 digits.

```cmd
python validation\check_censored_beta_t.py --module-dir build\gfortran_2806A042DAA8BBF7 --library build\gfortran_3C159EB802C22FB9\scifort\libscifort.a
```

The checker uses mpmath beta integrals at 100 digits for 300 seeded mixed
samples, 1100 digits for stresses, and 90 digits for maximum-shape cases.
Direct density integration supplies adjacent and large symmetric near-center
references; shapes equal to one use exact beta identities. Reference tail
coordinates are formed independently, including when their complements
round to one. Finite comparisons require
`abs(actual-reference)/max(1,abs(reference)) <= 2e-11`.

Checked flags are `-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow`;
optimized flags are `-O3 -fcheck=bounds`. Seed: `20261117`.

Independent SciPy bounded optimization uses exact `[0.2,0.4,0.6]`, left
`[0.1]`, right `[0.7]`, intervals `[0.25,0.35]`, `[0.5,0.55]`, fixed loc=0,
scale=1.3, first-shape bounds `[0.3,5]`, and fixed second shape 1.7 where used:

| Family | First shape | NNLF |
| --- | ---: | ---: |
| Student t | 4.999999881819905 | 13.595649458888033 |
| F | 1.1741350701452165 | 12.162483724958465 |
| Jones-Faddy | 2.2316074118819835 | 13.427914058395741 |

Student t's optimum is at the upper bound; SciFort reaches 5 exactly.
Native tolerances are 3e-5 in shape and 2e-8 in scaled objective.

## Integration

The full `fpm test` run passed with exit code 0 in 1832.629 seconds
(30 minutes 33 seconds), using the compiler identified above. It includes
all 1,846 new checks and the existing beta-prime, fitting, statistics,
and C API integration tests. The full run covers the small private
beta-prime callback extension as well as the new code and satisfies the
supplied repository instructions. It used the final numerical source and
test state; subsequent edits only finalized documentation.

Existing array-temporary runtime warnings appeared in unchanged code.
The console-input scan (`python tools/check_test_inputs.py`), `git diff --check`,
and audits of source ASCII/line lengths, Python syntax, changed documentation
links, and the 74 implemented/35 remaining family inventory passed. No C ABI
or wrapper surface was expanded.

## Numerical and fitting limits

Require normal positive beta shapes (at least tiny(dp)) and sum <=1e6.
For Student t these shapes are df/2,df/2; for F dfn/2,dfd/2; for Jones-Faddy
a,b. Student t additionally accepts positive-infinite df in likelihoods
and NNLF as the normal limit. Subnormal beta shapes and larger finite sums
remain explicit numerical gaps. Existing beta-kernel accuracy and
representability limits apply; no new large-shape asymptotic is introduced.
All-free fits may reach the iteration limit with finite improving candidates.
Successful local fits do not establish global optimality or identifiability.
See the [API guide](../docs/censored_fitting.md) and
[provenance](../CODE_PROVENANCE.md). C ABI/wrapper expansion is deferred.
