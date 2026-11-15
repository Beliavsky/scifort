# Censored Cauchy validation, 2026-11-15

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
Reference tools: SciPy 1.15.3 and mpmath 1.3.0.

All 248 new native checks passed at both `-O0` and `-O3`, retaining
`-fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace` in both builds.
The change adds Cauchy-specific callbacks, one fit dispatch branch and a fit
overload. Existing accumulation, standardization, optimization, scalar Cauchy,
and array likelihood/fit algorithms are unchanged.

The regression checks mixed censoring, negative log-likelihood, free/fixed
fits, objective agreement and optimizer option forwarding, reflection,
ordinary array likelihood/fit and scalar-kernel agreement, and invalid data,
parameters, bounds and tolerances. A symmetric three-point exact sample
`[-1,0,1]` independently checks the fitted location 0 and scale `1/sqrt(3)`.
One-sided censoring checks bounded boundary convergence. Interval tests cover
ordinary mpmath references, nesting, median endpoints and crossings, branch
boundaries, subnormal widths, underflowed standardized widths, adjacent far
tails, overflowing endpoint products and raw displacements, standardized
values beyond binary64 range, equal/infinite endpoints, and empty data.

Near-one CDF, SF and wide-interval log probabilities are checked with relative
ratios, including a subnormal wide-interval log probability. This prevents
a rounded zero from passing solely because of an absolute tolerance. Other
checks use regime-specific absolute tolerances: ordinary interval references
allow `2e-15`, logarithmic extreme-scale cases allow up to `8e-13`, and
mixed fitted parameters allow `4e-4` for different local optimizers.

## Independent references

The mixed sample has exact values `[0.5,1.2,2.5]`, left limit `[0.3]`, right
limit `[3.5]`, and intervals `[1.1,1.4]` and `[4.0,4.4]`. At location 0.4 and
scale 1.3, the independent log-likelihood is `-16.18709661147815`. SciPy's
censored fit returns location `1.3816244599074707` and scale
`1.1371319252147074`.

`validation/check_censored.py` passes all 1,173 comparisons at both `-O0`
and `-O3`, including 119 Cauchy comparisons: 100 seeded mixed samples and
19 stress intervals. Cauchy references use 700-digit arithmetic and direct
high-precision arctangent CDF/SF differences, independently of the production
same-side interval identity and crossing branches. Random samples include
same-side and forced median-crossing intervals. The seed is 20261112.

The maximum scaled Cauchy log-likelihood error is `3.1e-15`; the overall
maximum remains `1.08e-12` in an existing gamma stress case. The reference
harness uses error scale `max(1,abs(reference))` and threshold `2e-11`.
The native relative checks above additionally constrain tiny near-one log
probabilities that the harness's absolute error measure would not resolve.

The Cauchy interval stress cases cover the smallest binary64 width, very
small and large scales, standardized widths near `1e-632`, adjacent
endpoints near +/-`1e300`, an interval near `1e308` with the smallest scale,
overflowing raw width/displacement, nearly full-support intervals, one
endpoint at the median, and the crossing branch boundary.

Commands (CMD; FPM directory hashes are local artifacts):

```cmd
python tools\generate_censored_reference.py
python validation\check_censored.py --module-dir build\gfortran_2806A042DAA8BBF7 --library build\gfortran_3C159EB802C22FB9\scifort\libscifort.a
```

The reference harness explicitly compiles the changed censored module and
its reference driver before linking unchanged archive dependencies. It uses
`-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace` and
`-O3 -fcheck=bounds`.

## Native builds and guide example

Create `build\cauchy-native` before compiling:

```cmd
mkdir build\cauchy-native
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Jbuild\cauchy-native -Ibuild\cauchy-native -Ibuild\gfortran_2806A042DAA8BBF7 src\special\log_gamma.f90 src\empirical_distribution.f90 src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored_cauchy.f90 build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\cauchy-native\test.exe
build\cauchy-native\test.exe
```

The optimized build uses `-O3` and `build\cauchy-optimized`, retaining all
runtime checks and arithmetic traps. The new guide example is wrapped in a
program, compiled with the checked flags above against the changed modules,
and run successfully. `python tools\check_test_inputs.py` also passes.

## Limits

Fits reuse the existing bounded local coordinate/pattern search. Cauchy
likelihoods can have multiple stationary points; convergence does not
establish global optimality, identifiability, or existence of an unconstrained
finite MLE. Use meaningful bounds and consider multiple guesses. Individual
Cauchy log densities and small-tail logs remain finite for finite observations
and valid finite parameters, even when the standardized distance overflows;
accumulated sums must still fit in binary64. Scalar/array routines retain
their existing numerical range. No observation weights, dependencies, or
C ABI symbols are added.

## Full-suite verification and final audit

The full `fpm test` run passed with exit code 0 in 1271.509 seconds
(about 21.2 minutes), including all 248 new Cauchy checks. The run
used the compiler above and FPM's default debug profile. The equivalent CMD
command for the captured invocation is:

```cmd
fpm test > build\cauchy-full-test.log 2>&1
```

The actual runner used Python `subprocess.run(['fpm','test'])` with stdout
and stderr redirected to that log, preserving the native exit status.
Existing array-temporary warnings were emitted; no test failures occurred.
Changed Fortran ASCII/132-column limits, Python syntax, newly added local
Markdown links, the censoring inventory counts, and `git diff --check`
passed. No build artifacts or unrelated edits are included. Subsequent edits
only record these results.
