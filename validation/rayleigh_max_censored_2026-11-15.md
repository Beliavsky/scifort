# Rayleigh and Weibull maximum censoring validation, 2026-11-15

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
References: SciPy 1.15.3 and mpmath 1.3.0.

## Native tests

All 245 new checks passed at both `-O0` and `-O3`, using
`-g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace` in both builds.
The checks cover exact/left/right/interval observations, negative likelihoods,
free/fixed fitting parameters and iteration limits, objective agreement,
scalar/array API agreement, reflection, and shape-one/shape-two identities.
Rayleigh exact/right censoring independently checks the analytic fixed-location
scale MLE `sqrt(sum(exact**2)+sum(right**2))/(sqrt(2*n_exact))` at location zero.

Support tests cover zero Rayleigh endpoint density, zero/unit endpoint tails,
Weibull maximum endpoint densities for shapes below/equal/above one, sums of
opposite infinite likelihood terms returning NaN, clipped/impossible/full-support
intervals, invalid parameters and observations, empty data, and fit error
statuses for nonfinite objectives. The smallest and largest finite Weibull
shapes are checked at unit displacement without arithmetic overflow.

Stress tests cover subnormal widths, underflowed standardized widths, adjacent
far-tail intervals, huge raw widths/displacements with representable scaled
values, Rayleigh scales whose conversion to Weibull scale would overflow,
and unrepresentable log densities. Tiny near-one logs are checked by ratios,
including almost-full-support intervals with excluded mass at both ends.

Ordinary scalar agreement uses `2e-14` or `3e-14` absolute tolerances; mixed
likelihoods use `6e-14` or `7e-14`. Log-scale extremes allow up to `8e-13`,
and adjacent-tail density-integral approximations allow up to `2e-12`.
Near-one log ratios allow `3e-12` for exponential sensitivity to rounded
squared endpoints/hazards, with `5e-13` for a tiny Rayleigh hazard. Fitted
parameters allow `4e-4` for differing optimizers. The independent analytic
Rayleigh fit allows `3e-6`. No blanket tolerance is used for all regimes.

## Independent reference comparisons

`validation/check_censored.py` passed all 1,621 comparisons at checked `-O0`
and optimized `-O3`. There are 218 new comparisons: 100 mixed samples per
family, ten Rayleigh interval stress cases, and eight Weibull maximum cases.
They use 800-digit arithmetic, direct CDF/SF differences, and seed 20261112,
independently of the production hazard-increment interval calculation.
The maximum scaled error for the new families is `1.08e-14`; the overall
maximum remains `1.08e-12` in an existing gamma stress case. The harness
scales by `max(1,abs(reference))` with threshold `2e-11`. Native relative
near-one checks additionally resolve tiny logs hidden by that absolute scale.

The mixed Rayleigh sample has exact `[0.5,1.2,2.5]`, left `[0.3]`, right `[3.5]`,
and intervals `[1.1,1.4]`, `[4.0,4.4]`. At location 0 and scale 1.3, its
log-likelihood is `-17.949503441709364`. SciPy's fixed-location censored fit
returns scale `1.8108567877921697`.

The Weibull maximum sample reflects those observations, exchanges left/right
censoring, and reverses interval endpoints. At shape 1.7, location 0, scale 1.3,
its log-likelihood is `-20.5129071280305`. SciPy's fixed-location fit returns
shape `1.0192387961160727` and scale `2.22130765498472`. Native tests check
reflection against minimum-Weibull at five shapes, shape one against reflected
exponential, and shape two against reflected Rayleigh with the correct scale.
The reference generator also checks all three endpoint density regimes.

Commands (CMD; matching local FPM directory hashes are required):

```cmd
python tools\generate_censored_reference.py
python validation\check_censored.py --module-dir build\gfortran_2806A042DAA8BBF7 --library build\gfortran_3C159EB802C22FB9\scifort\libscifort.a
```

The reference harness compiles the changed censored module and driver explicitly
before linking unchanged dependencies. Its configurations are
`-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace` and
`-O3 -fcheck=bounds`.

## Native builds and documentation example

Create `build\hazard-native` before compiling:

```cmd
mkdir build\hazard-native
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Jbuild\hazard-native -Ibuild\hazard-native -Ibuild\gfortran_2806A042DAA8BBF7 src\special\log_gamma.f90 src\empirical_distribution.f90 src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored_rayleigh_max.f90 build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\hazard-native\test.exe
build\hazard-native\test.exe
```

The optimized build uses `-O3` and `build\hazard-optimized`, retaining all
runtime checks and traps. The guide example is wrapped in a program, compiled
with the checked configuration against the changed modules, and run.

## Limits

Existing scalar/array kernels and the bounded optimizer are unchanged. Fits
remain local box-constrained searches; success does not guarantee a global
finite MLE or identifiability. Free support locations and Weibull shapes below
one can create infinite endpoint likelihoods. Likelihoods preserve those
infinities; the optimizer rejects nonfinite candidates, and fully fixed
nonfinite objectives return `fit_status_no_finite_objective`. Sums must fit
in binary64. Weibull maximum narrow intervals inherit the minimum-Weibull
increment kernel's accuracy/range limits. No new runtime dependency, observation
weights, imported implementation code, or C ABI is introduced.

## Full-suite verification and final audit

The full `fpm test` run passed with exit code 0 in 1800.074 seconds
(about 30.0 minutes), including all 245 new native checks. The run
used the compiler above and FPM's default debug profile. The equivalent CMD
command for the captured invocation is:

```cmd
fpm test > build\hazard-full-test.log 2>&1
```

The actual runner used Python `subprocess.run(['fpm','test'])`, redirecting
stdout and stderr to that log and preserving the native exit code. Existing
array-temporary warnings were emitted; there were no test failures. The run
covers the final numerical code and tests.

`python tools\check_test_inputs.py` passed. Changed Fortran ASCII/132-column
limits, Python syntax, newly added local documentation links, the inventory
counts (15 supported families and 94 remaining), and `git diff --check`
passed. The guide example compiled and ran. No build artifacts, imported
implementation code, new dependencies, or unrelated edits are included.
Subsequent edits only record these validation results.
