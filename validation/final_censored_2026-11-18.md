# Censored stable and studentized-range validation

Date: 2026-11-18.

See the [API contracts](../docs/censored_fitting.md#levy-stable-and-studentized-range-censoring)
and [provenance](../CODE_PROVENANCE.md). All 109 named continuous families now
have censored entry points; numerical ranges remain explicitly restricted.

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
Python 3.13.3, SciPy 1.15.3, mpmath 1.3.0.

```cmd
C:\Programs\Python313\python.exe validation/check_censored_final.py
fpm test
```

The independent checker prints exact compiler/link commands and builds both
-O0 and -O3 with -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace.
It compiles empirical_distribution, censored, fit, and stats, links the native
regression program and comparison driver against the existing FPM archive,
and runs both. Optional --reuse reuses ignored reference fixtures while
rebuilding both profiles.

References use mpmath defining characteristic-function integrals, folded-t
incomplete-beta identities (including extreme tails), lower-end range moment
coefficients, positive adaptive normal-range/chi integrals, and ordinary SciPy
cross-checks. Small-probability adaptive references force a tight absolute
integration tolerance: SciPy's default studentized-range integration can lose
relative accuracy there. The independent stable central kernel avoids the
scalar stable evaluator's near-mode/near-one rounding; scalar-stable routines
are not tight density oracles for this new kernel. An existing scalar stable
PDF evaluation at zero triggered SIGFPE under traps during test development;
that unrelated implementation is unchanged. Independent Fourier references
cover zero instead.

The native tests cover mixed censoring, supported special cases, symmetry and
complementarity, endpoint behavior, adjacent intervals, affine scales, shape
restrictions, NaN/invalid data, fixed/free fits and retained array overloads.
Independent fit objectives use Cauchy and the folded-t identity.

Both -O0 and -O3 builds passed all 176 native checks and 277 reference
comparisons, with runtime checks and floating-point traps enabled. The maximum
comparison error across all references was 5.80e-10; each case met its own
criterion: 3e-10 scaled log error for independent references, 3e-8 for ordinary
SciPy comparisons, and relative error for tiny complementary multiprecision
logs. Support endpoints and invalid inputs use exact semantic checks.
The tests also retain finite densities, survival logs, and adjacent-interval
logs when x=1e300 and scale=1e-200 would overflow standardization.

The independent free-scale fits used this mixed sample: exact [.5,1,2], left
[.8], right [3], interval [1.2,1.4]. The stable Cauchy special-case fit gives
scale 1.3339167533337817 and NNLF 12.176776394885268. The k=2, df=10 folded-t
fit gives scale 1.2559117716600137 and NNLF 9.256777181772247. Native fit tests
compare parameters with tolerance 3e-6 and scaled objectives with 3e-10.

The checker uses the local FPM module/archive cache paths printed in its
compiler commands; build that cache first. Different compiler/cache layouts
require adjusting those paths. No generated fixtures or compiled artifacts
are tracked.

The full `fpm test` run passed with exit code 0 in 1611.175 seconds
(26 minutes 51 seconds), including all 176 new native checks and existing
distribution/fit regressions. Existing array-temporary warnings from scattered
interpolation were emitted. The full suite used the normal FPM configuration;
the focused -O0/-O3 builds above supplied bounds/runtime checks and traps.

Final audit passed: changed Fortran is ASCII and within the standard free-form
line limit; Python syntax and changed relative Markdown links are valid;
109 distinct supported family names match 109 documented likelihood entries,
with zero remaining named families. `git diff --check` passed. No generated
fixtures, binary files, build artifacts, or license-notice removals are included.
