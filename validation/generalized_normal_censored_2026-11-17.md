# Generalized normal and R censored validation, 2026-11-17

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
References: SciPy 1.15.3 and mpmath 1.3.0.

## Native numerical and fitting checks

All 1,867 checks in [the native test](../test/test_censored_generalized_normal.f90)
passed at `-O0` and `-O3`, each with
`-g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace`.
They cover ordinary scalar density/CDF/SF agreement, CDF monotonicity,
CDF/PPF round trips, NNLF, original array-fit overloads, fixed/shape-only/
three-parameter fits, SciPy benchmarks, empty/invalid inputs, support
endpoints, equal/infinite interval endpoints, and adjacent intervals.

Identities cover generalized normal/Laplace, half generalized normal/
exponential, shape-two normal/half-normal with scale divided by sqrt(2),
and R shapes 1/2/3 versus arcsine/uniform/semicircular on matching support.
Stresses include an unrepresentable reciprocal or half-shape, the R shape
cap, overflowing log-gamma normalization with valid empty/nonempty data,
subnormal beta normalization, maximum supported R shape, subnormal
intervals, tiny center-crossing intervals, and reflected intervals.
The nearly uniform half-law endpoint log density `5.7721566490153286e-301`
at beta=1e300 is checked by a ratio against the independent Euler-constant
limit. A regression checks a steep powered-coordinate boundary layer
at beta=1e6, with independently computed log mass
`-13.630472120271284593` over the stated binary64 endpoints.
At the largest finite power, tests retain finite CDFs despite an overflowing
intermediate powered distance, and check the subnormal endpoint log CDF
`-1.2203636434998277e-309` using a relative comparison.

Ordinary scalar logs use `2e-12`, round trips `3e-12`, and most
adjacent/extreme logs and identities `3e-13`. Independent fits allow
`3e-5` in shape and `2e-8` in scaled objective, including a search approaching
a boundary. Exact support and invalid-input assertions use no tolerance.

Checked command (CMD, matching local FPM dependencies and an existing
output directory):

```cmd
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Jbuild\generalized_normal-native -Ibuild\generalized_normal-native -Ibuild\gfortran_2806A042DAA8BBF7 src\special\log_gamma.f90 src\empirical_distribution.f90 src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored_generalized_normal.f90 build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\generalized_normal-native\test.exe
build\generalized_normal-native\test.exe
```

The optimized run substitutes `-O3` and separate module/executable paths
under `build\generalized_normal-optimized`, retaining all runtime/trap flags.
Changed modules are compiled explicitly before linking the existing archive.

## Independent references

The [checker](check_censored_generalized_normal.py) passed all 383 cases in
each checked and optimized build: 300 seeded mixed samples and 83 stress
cases. Maximum scaled log-likelihood error was `1.855e-13`, in a generalized
normal mixed sample. Mixed samples use 100 decimal digits, scales from
`1e-100` to `1e100`, and interval widths down to adjacent floats. Stresses
use up to 1100 digits, powers from 1e-300 through the largest finite binary64
value, amplified scale-point
displacements, nearly full support, subnormal R shapes/widths, R shapes
through 1e6, and steep thin boundary layers. Eight upper-cap R tail and
adjacent-interval references use direct density quadrature at 90 digits,
including standardized coordinates +/-0.001 and +/-0.01.

The oracle independently integrates the mathematical densities with
incomplete-gamma and beta functions. A positive lower-gamma series avoids
unnecessary huge-shape upper-gamma work. Symmetric interval references use
the half-law logs before adding the two branches, retaining masses too
small to recover from probabilities rounded to one half. Memoization keys
include precision. The R oracle follows the documented rounded affine and
half-shape law. Finite comparisons require
`abs(actual-reference)/max(1,abs(reference)) <= 2e-11`; native tests separately
cover nonfinite results and tiny relative log-density errors.

```cmd
python validation\check_censored_generalized_normal.py --module-dir build\gfortran_2806A042DAA8BBF7 --library build\gfortran_3C159EB802C22FB9\scifort\libscifort.a
```

The checker prints exact compiler commands. Checked flags are
`-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow`; optimized flags are
`-O3 -fcheck=bounds`. Default seed: `20261117`.

Independent SciPy bounded optimization uses exact `[0.2,0.4,0.6]`, left
`[0.1]`, right `[0.7]`, intervals `[0.25,0.35]`, `[0.5,0.55]`, fixed loc=0,
scale=1.3, and shape bounds `[0.3,5]`:

| Family | SciPy shape | Negative log likelihood |
| --- | ---: | ---: |
| Generalized normal | 4.097050810785287 | 11.851573580953866 |
| Half generalized normal | 3.813252881858734 | 9.546996838141601 |
| R | 4.999999859344038 | 11.137246801642787 |

The R optimum is at the upper search bound 5; SciFort reaches it exactly.

## Integration and documentation

The full `fpm test` run passed with exit code 0 in 1846.995 seconds
(30 minutes 47 seconds), using the compiler identified above. This includes
all 1,867 new generalized-normal/R checks. The run used the final numerical
source and test state; subsequent changes only finalized documentation.
Existing array-temporary runtime warnings appeared in scattered interpolation.

The API-guide example compiled and ran successfully, returning half-generalized-
normal shape `3.8132529526948926` and NNLF `9.5469968381415988`.
The console-input scan (`python tools/check_test_inputs.py`), `git diff --check`,
and audits of source ASCII/line lengths, Python syntax, changed documentation
links, and the 68 implemented/41 remaining family inventory all passed.

## Numerical and fitting limits

Generalized-normal beta must have a finite binary64 reciprocal. R requires
`0<c<=1e6`, positive representable half-shape, and a finite noncollapsed
rounded affine support. Larger R shapes and reciprocal-overflow powers
remain explicit numerical coverage gaps, tracked separately from the
68 implemented/41 remaining named-family count.

R uses the rounded lower edge and half-shape; subnormal or shifted supports
can differ from the original scalar API. Fits use local bounded searches;
movable endpoints or singular limits can make maximization unbounded, and
success does not establish global optimality or identifiability. See the
[API guide](../docs/censored_fitting.md) and
[provenance](../CODE_PROVENANCE.md). No C ABI/wrapper expansion is included.
