# Folded distributions and generalized Pareto validation, 2026-11-17

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
References: SciPy 1.15.3 and mpmath 1.3.0.

## Native numerical and fitting checks

All 1,847 checks in [the native test](../test/test_censored_fold_pareto.f90)
passed at both `-O0` and `-O3`, each with
`-g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace`.
They cover scalar density/CDF/SF agreement, CDF monotonicity, CDF/PPF round
trips, original array-fit overloads, NNLF, fixed/shape-only/three-parameter
fits, empty and invalid inputs, support endpoints, equal/infinite interval
endpoints, adjacent intervals, and support clipping.

Independent identities cover folded normal/half-normal, folded
Cauchy/half-Cauchy, generalized Pareto/exponential, shape -1/uniform,
and shape 1/Lomax. Stresses include the smallest positive and negative
Pareto shapes with the smallest scale, normal squared-distance overflow,
finite Cauchy log density beyond standardized-distance overflow, tiny
Cauchy CDF with huge folding offset, subnormal intervals with huge scale,
and zero/finite/infinite bounded upper-endpoint densities.

Ordinary scalar log tolerances are `2e-12`, round trips `3e-12`, and
adjacent/extreme logs and identities `3e-13` where applicable. Independent
fit benchmarks allow `3e-5` in shape and `2e-8` in scaled objective,
including an optimizer approaching a boundary rather than reaching it.

Exact checked command (CMD, matching local FPM dependencies and an
existing output directory):

```cmd
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Jbuild\fold_pareto-native -Ibuild\fold_pareto-native -Ibuild\gfortran_2806A042DAA8BBF7 src\special\log_gamma.f90 src\empirical_distribution.f90 src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored_fold_pareto.f90 build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\fold_pareto-native\test.exe
build\fold_pareto-native\test.exe
```

The optimized run substitutes `-O3` and separate module/executable paths
under `build\fold_pareto-optimized`, retaining all runtime/trap flags.
Changed modules are compiled explicitly before linking the existing archive.

## Independent high-precision references

The [reference checker](check_censored_fold_pareto.py) passed all 390 cases
in each checked and optimized build: 300 seeded mixed samples and 90 finite
stress cases. Maximum scaled log-likelihood error was `8.582e-14`, in a
generalized Pareto shape stress. Ordinary cases use 100 decimal digits and
scales from `1e-100` to `1e100`; stresses use 1100 digits.

The oracle independently integrates the public densities: arbitrary-
precision normal erfc differences, positive Cauchy arctangent differences,
and generalized Pareto cumulative hazards. The bounded Pareto oracle
follows the documented rounded-width law. Stresses include zero, tiny,
huge, and negative shapes, huge folding centers, extreme normal tails,
adjacent intervals, nearly full support, subnormal widths, and overflowing
standardized Cauchy distances. It requires finite actual/reference values
with `abs(actual-reference)/max(1,abs(reference)) <= 2e-11`; native tests
separately check nonfinite endpoints and invalid-input contracts.

```cmd
python validation\check_censored_fold_pareto.py --module-dir build\gfortran_2806A042DAA8BBF7 --library build\gfortran_3C159EB802C22FB9\scifort\libscifort.a
```

The checker prints its exact compiler commands. Checked flags are
`-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow`; optimized flags are
`-O3 -fcheck=bounds`. Default seed: `20261117`.

Independent SciPy bounded optimization uses exact `[0.2,0.4,0.6]`,
left `[0.1]`, right `[0.7]`, and intervals `[0.25,0.35]`, `[0.5,0.55]`,
with fixed location zero and scale 1.3:

| Family | SciPy numerical shape | Negative log likelihood | Bounds |
| --- | ---: | ---: | --- |
| Folded normal | 1.760720062091216e-8 | 11.332155282192696 | `[0,3]` |
| Folded Cauchy | 7.616843834759234e-9 | 12.774437362573709 | `[0,3]` |
| Generalized Pareto | -0.9999999717159224 | 9.948277979270488 | `[-1,3]` |

SciFort reaches the exact boundaries 0, 0, and -1. At generalized Pareto
shape -1 the objective is `9.948277934580510`.

## Integration and documentation

The full integration command passed with exit code 0:

```cmd
fpm test
```

Elapsed time was 1,122.875 seconds (18 minutes 43 seconds), including
compilation and execution. The full run includes all 1,847 new native
checks and covers the shared likelihood module, fitting dispatch,
statistics exports, and their dependents. No executable changes followed
this successful run; subsequent validation-record edits are documentation only.

The new guide example compiled and ran with checked flags, reproducing the
generalized Pareto boundary fit. The console-input scan passed. Source
ASCII/line lengths, Python syntax, changed Markdown links, and inventory/API
consistency checks and `git diff --check` passed (65 supported plus 44 remaining). No new
dependency, generated artifact, or third-party implementation is tracked.

## Limits

Binary64 is the public precision. Subtraction near a large folding offset
is ill-conditioned; uniform relative accuracy through such cancellation is
not promised. Logarithms of normal densities/tails beyond the binary64
magnitude range return negative infinity.

Negative generalized Pareto shapes define a normalized law on a rounded
physical width; underflowed, overflowed, or collapsed supports return NaN.
Compensated gaps retain translated-edge residuals. Subnormal width rounding
can differ from the unchanged scalar APIs. Free endpoints and singular
limits can make maximization unbounded; bounded local fit success does not
guarantee a global maximum or unique identifiable parameters. See the
[API guide](../docs/censored_fitting.md) and
[provenance](../CODE_PROVENANCE.md). No C ABI/wrapper changes are included.
