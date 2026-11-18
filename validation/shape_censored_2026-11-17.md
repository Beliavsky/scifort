# Alpha, GEV, and Crystal Ball censoring validation, 2026-11-17

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
Independent references: SciPy 1.15.3 and mpmath 1.3.0.

## Native regression

The [native regression](../test/test_censored_shape.f90) contains 2,174 checks:
ordinary scalar PDF/CDF/SF agreement, CDF monotonicity, quantile round trips,
adjacent and ordinary intervals, Crystal Ball intervals on both sides of its
transition, affine scales 1e-200/1/1e200, invalid and empty inputs, NaN/infinity
classification, one-sided infinite intervals, equal endpoints as exact data,
support endpoints, fixed/one-parameter/all-free fits, original array-fit
generic compatibility, and NNLF. GEV at c=0 is checked against Gumbel.

Extreme regressions retain alpha's tiny SF and tiny complementary log-CDF,
steep adjacent alpha intervals, tiny Crystal Ball normalization complements,
overflowed standardized coordinates, and finite quadratic log tails near
the largest binary64 magnitude. The last cases compare with the leading
quadratic term: omitted linear, logarithmic, and normalization corrections
are below 1e-150 relative to that magnitude. This is not an absolute
tail-accuracy claim.

CMD commands, with existing output directories and FPM dependencies:

```cmd
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Jbuild\shape-native -Ibuild\shape-native -Ibuild\gfortran_2806A042DAA8BBF7 src\special\log_gamma.f90 src\empirical_distribution.f90 src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored_shape.f90 build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\shape-native\test.exe
build\shape-native\test.exe
```

The optimized command substitutes `-O3` and paths under
`build\shape-optimized`, retaining all checks and traps. Changed modules
are explicitly compiled before the existing archive is linked.
Both builds pass all 2,174 checks.

Ordinary scalar, adjacent, and affine log comparisons allow 3e-12 scaled
error; ordinary intervals 3e-11, quantile round trips 5e-11, Gumbel identity
2e-12. Tiny complementary logs are divided by their expected magnitude
before comparison, with tolerance 2e-13. The scaled error measure is
abs(actual-reference)/max(1,abs(reference)). Invalid/infinite classifications
and exact endpoint values use no tolerance.

## Independent references

All 498 comparisons pass in checked and optimized builds: 300 seeded
mixed samples, 100 adjacent affine transition intervals, and 98 stresses.
Maximum scaled error is `2.2606e-14`.

```cmd
python validation\check_censored_shape.py --module-dir build\gfortran_2806A042DAA8BBF7 --library build\gfortran_3C159EB802C22FB9\scifort\libscifort.a
```

The [checker](check_censored_shape.py) independently evaluates normal
integrals, GEV exponential/power antiderivatives, and the two integrated
Crystal Ball pieces using mpmath. Ordinary samples use 100 decimal digits;
stresses use 1100. Seed: 20261117. Finite comparisons require scaled error
<=2e-11. The final guard changes were also rechecked against all cached
high-precision inputs.

Stresses cover alpha shapes 1e-6/1/100 and points 1e-100/1e100/1e300;
GEV shapes -100/-1/-1e-6/0/1e-6/1/100; Gumbel tails through x=1000;
Crystal Ball beta 1e-6/2/100, m immediately above 1/3/1e6, Gaussian tails
at x=40, and power tails down to x=-1e300. Each has adjacent intervals.
Affine scales 1e-300/1e300 and overflowed standardized coordinates are
included. Checked flags are `-O0 -g -fcheck=all
-ffpe-trap=invalid,zero,overflow`; optimized flags are `-O3 -fcheck=bounds`.

The transition regressions retain the total physical width across the join;
independently standardized endpoints can otherwise change adjacent interval
masses substantially. For very narrow crossing intervals the implementation
uses the shared Gaussian density and slope, with a relative omitted correction
below 1e-12 under the documented width criterion. Native tests include 50
deterministic affine joins; the reference checker includes 100 seeded joins
using exact mpmath antiderivatives of both pieces.

Independent SciPy bounded fits use exact [0.5,0.8,1.1], left [0.4],
right [1.7], and intervals [0.6,0.7] and [0.9,1.0]. Only scale varies in
[0.3,3], with loc=0:

| Family | Fixed shapes | Fitted scale | NNLF |
| --- | --- | ---: | ---: |
| alpha | a=1.7 | 1.2307260337339929 | 10.399301854809206 |
| genextreme | c=0.2 | 0.90425515523373112 | 13.013831218213928 |
| crystalball | beta=1.5, m=2.7 | 1.1798418087757312 | 14.678207227692237 |

Native benchmarks allow 3e-5 scaled parameter error and 2e-8 scaled NNLF
error. All-free checks require finite improving candidates, without claiming
global optima or parameter identifiability.

## Full suite and limitations

The shared censoring/fit/statistics integration changed, so the full suite
is required in addition to focused checked and optimized validation:

```cmd
fpm test
```

The full suite passed with exit code 0 in 1152.242 seconds (19 minutes
12 seconds), including all 2,174 new checks and the existing censored,
scalar-distribution, fitting, C API, and integration regressions. Existing
array-temporary runtime warnings appeared and did not cause failures.
No numerical code changed after this successful run began; subsequent
edits only completed documentation and validation metadata.

The [API guide](../docs/censored_fitting.md) documents finite shape limits:
alpha a in [1e-6,100]; GEV abs(c)<=100 and c=0 or abs(c)>=1e-6;
Crystal Ball beta in [1e-6,100], m in (1,1e6]. Wider ranges, notably
very small nonzero GEV shapes, remain compatibility gaps. Extremely close
shape-dependent support endpoints and concentrated affine inputs remain
ill-conditioned relative to binary64 input rounding. These additions do
not change scalar distribution implementations or add C ABI wrappers.
