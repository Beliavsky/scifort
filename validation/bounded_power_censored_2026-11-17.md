# Log-uniform, truncated Pareto, and generalized half-logistic validation

Date: 2026-11-17. Compiler: GNU Fortran (GCC) 17.0.0 20260510
(experimental), Windows. References: SciPy 1.15.3, mpmath 1.3.0.

## Native checks

The native program is [test_censored_bounded_power.f90](../test/test_censored_bounded_power.f90).
All 1,858 checks passed at both `-O0` and `-O3`, each with bounds/runtime
checks and traps for invalid operations, division by zero, and overflow.
They cover ordinary scalar density/CDF/SF agreement, CDF monotonicity,
CDF/PPF round trips, NNLF, original array-fit overloads, fixed/free fits,
SciPy shape benchmarks, invalid inputs, empty samples, endpoints, equal
and infinite interval endpoints, whole-support and adjacent intervals.
Both four-parameter fits and a negative-shape Pareto fit are exercised.

Independent identities cover the Pareto limit as either signed shape tends
to zero and generalized half-logistic's rational CDF/density at unit shape.
Stresses include the smallest positive shape/scale, underflowed hazards,
adjacent cutoffs, subnormal single-rounded products, collapsed/overflowed
supports, and zero/finite/infinite upper endpoint density. A regression
checks a cutoff product admitted by an upward-rounded quotient guard that
would otherwise overflow, together with its neighboring finite product.

Ordinary scalar log tolerances are `2e-12`, round trips `3e-12`, adjacent
intervals/extreme logs `3e-13`, and exact identities `1e-14` where applicable.
Fit benchmarks allow `3e-5` in shape and `3e-10` in scaled objective.

Checked command, using matching local FPM dependencies and an existing
output directory (CMD):

```cmd
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Jbuild\bounded_power-native -Ibuild\bounded_power-native -Ibuild\gfortran_2806A042DAA8BBF7 src\special\log_gamma.f90 src\empirical_distribution.f90 src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored_bounded_power.f90 build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\bounded_power-native\test.exe
build\bounded_power-native\test.exe
```

The optimized run substitutes `-O3` and separate module/executable paths
under `build\bounded_power-optimized`, retaining all runtime/trap flags.
Changed modules are compiled explicitly before linking the existing archive.

## Independent references

The [reference checker](check_censored_bounded_power.py) evaluates the
mathematical integrals directly with arbitrary precision. It follows the
documented rounded physical supports. There are 360 comparisons per build:
300 seeded mixed samples and 60 stress cases. Mixed samples use 100 decimal
digits, scales from `1e-100` to `1e100`, both Pareto shape signs, and intervals
as narrow as adjacent floating-point values. Stresses use 1100 digits for
tiny/huge shapes, endpoints, adjacent cutoffs, nearly full support,
subnormal widths, and rounded subnormal products.

Finite reference comparisons require
`abs(actual-reference)/max(1,abs(reference)) <= 2e-11`.
All 360 comparisons passed in each build. Maximum scaled error was
`1.6231e-13`, in a generalized half-logistic mixed sample. Native tests
separately cover nonfinite values and exact invalid-input contracts. The checker uses `-O0 -g -fcheck=all
-ffpe-trap=invalid,zero,overflow` and `-O3 -fcheck=bounds`, printing the
exact compiler commands:

```cmd
python validation\check_censored_bounded_power.py --module-dir build\gfortran_2806A042DAA8BBF7 --library build\gfortran_3C159EB802C22FB9\scifort\libscifort.a
```

Default seed: `20261117`. SciPy independently optimizes the first shape
with exact `[1.2,1.4,1.6]`, left `[1.1]`, right `[1.7]`, intervals
`[1.25,1.35]`, `[1.5,1.55]`, fixed location zero and scale one.
Log-uniform's upper shape and Pareto's cutoff are fixed at 3:

| Family | Shape | Negative log likelihood | Search bounds |
| --- | ---: | ---: | --- |
| Log-uniform | a=0.9306175899682245 | 10.427146663355687 | `[0.1,1.09]` |
| Truncated Pareto | b=2.112500505038408 | 9.220560683106488 | `[0.01,5]` |
| Generalized half-logistic | c=0.3800976249438511 | 11.999191342035351 | `[0.05,0.58]` |

## Integration

The full integration command passed with exit code 0:

```cmd
fpm test
```

Elapsed time was 1,798.869 seconds (29 minutes 59 seconds), including
compilation and execution. The full run includes all 1,858 new native
checks and covers the shared likelihood module, fitting dispatch,
statistics exports, and their dependents. No executable changes followed
this successful run; subsequent validation-record edits are documentation only.

The guide example compiled and ran with checked flags, reproducing the
truncated Pareto benchmark. Final checks passed for source ASCII and line
lengths, Python syntax, changed Markdown links, inventory/API consistency
(62 supported plus 47 remaining), `git diff --check`, and
`python tools\check_test_inputs.py`. No generated artifact, new dependency,
or third-party implementation was added to the tracked changes.

## Limits

Binary64 is the public precision. Rounded physical cutoffs/widths define
normalized bounded laws; underflowed, overflowed, or collapsed supports
return NaN. Addition residuals remain in endpoint gaps. This can differ
from the unchanged scalar APIs, especially for subnormal products.

Truncated Pareto excludes zero shape. Search its negative and positive
branches separately: local searches spanning zero can settle on the
poorer branch. Shape/scale redundancy affects log-uniform identifiability.
Movable endpoints and singularities can produce unbounded likelihoods.
Bounded local fit success does not establish a global maximum or unique
identifiable parameters. See [the API guide](../docs/censored_fitting.md)
and [provenance](../CODE_PROVENANCE.md). No C ABI/wrapper changes are included.
