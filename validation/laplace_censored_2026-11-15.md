# Censored Laplace validation, 2026-11-15

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
Reference tools: SciPy 1.15.3 and mpmath 1.3.0.

All 228 new native checks passed at both `-O0` and `-O3`, retaining
`-fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace` in both builds.
The existing 373 censored and 150 Gumbel checks are also exercised against
the changed modules. These cover the generic accumulator and fitter paths
reused by Laplace; no existing optimizer or distribution algorithm is changed.

The native Laplace regression checks mixed censoring, negative likelihood,
free/fixed fits and option forwarding, reflection, array likelihood/fit
agreement, exact-sample median and mean-absolute-deviation estimates, scalar
kernel agreement, ordinary/nested/cusp-crossing intervals, subnormal widths,
underflowed standardized widths, adjacent intervals in both extreme tails,
overflowing raw widths/displacements with representable standardized values,
tiny scales, infinite/equal endpoints, invalid/empty data and parameters,
fit error statuses, and all-right-censored boundary convergence.

The mixed sample uses exact values `[0.5,1.2,2.5]`, left limit `[0.3]`, right
limit `[3.5]`, and intervals `[1.1,1.4]` and `[4.0,4.4]`. At location 0.4 and
scale 1.3 the independent log-likelihood is `-16.62410935200545`. SciPy's
censored fit returns location `1.2500073446335311` and scale
`1.6532044386643419`; parameter comparisons allow `4e-4` for different local
optimizers. Exact-sample tests independently check the median and mean
absolute deviation rather than relying exclusively on SciPy fitting.

## Independent likelihood checks

`validation/check_censored.py` passes all 1,054 comparisons at both `-O0`
and `-O3`. There are 109 Laplace comparisons: 100 seeded mixed samples and
nine stress intervals. The maximum scaled Laplace log-likelihood error is
`9.15e-16`; the maximum across all families is `1.08e-12` in an existing
gamma stress case. The scale is `max(1,abs(reference))`, with acceptance
threshold `2e-11`. Native tests additionally use regime-specific absolute
tolerances, down to exact agreement where appropriate.

Laplace references use 700-digit arithmetic and independently subtract
high-precision log tails. Random samples include intervals on either side
of the center plus a forced center-crossing interval, with seed 20261112.
Stress cases include the smallest positive binary64 width, a subnormal
center-crossing interval, underflowed standardized widths on one side and
across the center, adjacent intervals near both +/-1000, overflowing raw
width/displacement, and a narrow scale. The reference checker compiles the
changed likelihood modules explicitly before linking unchanged dependencies.

Commands (CMD; local FPM hashes may differ):

```cmd
python tools\generate_censored_reference.py
python validation\check_censored.py --module-dir build\gfortran_2806A042DAA8BBF7 --library build\gfortran_3C159EB802C22FB9\scifort\libscifort.a
```

## Native builds and documentation example

Create `build\laplace-native` before compiling:

```cmd
mkdir build\laplace-native
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Jbuild\laplace-native -Ibuild\laplace-native -Ibuild\gfortran_2806A042DAA8BBF7 src\special\log_gamma.f90 src\empirical_distribution.f90 src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored_laplace.f90 build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\laplace-native\test.exe
build\laplace-native\test.exe
```

The optimized build uses `-O3` and `build\laplace-optimized`, retaining all
runtime checks and arithmetic traps. The existing censored/Gumbel regressions
use the same commands with their respective test source and output names.
The new guide example is wrapped in a program, compiled with the checked
flags above against the changed modules, and run successfully.

## Limits

Fits reuse the existing bounded local coordinate/pattern search. Convergence
does not establish global optimality, identifiability, or existence of an
unconstrained finite MLE. There are no observation weights or C ABI additions.
Individual unrepresentable log densities/tails return negative infinity;
finite sums must still fit in binary64. Existing array APIs are preserved.

## Full-suite and final checks

The full `fpm test` run passed with exit code 0 in 1379.196 seconds
(about 23 minutes), including all 228 new Laplace checks. The run used the
compiler above and the default FPM debug profile. The equivalent CMD command
for the captured invocation is:

```cmd
fpm test > build\laplace-full-test.log 2>&1
```

The actual runner used Python `subprocess.run(['fpm','test'])` with stdout
and stderr redirected to that log, preserving the native exit code. Existing
array-temporary runtime warnings were emitted; there were no test failures.
`python tools\check_test_inputs.py` also passed. Changed Fortran was checked
for ASCII and the 132-column limit; changed Python parsed successfully;
added local Markdown links and `git diff --check` passed. The guide example
compiled and ran. No build artifacts, new dependencies, or unrelated edits
are included. Subsequent edits only record these results.
