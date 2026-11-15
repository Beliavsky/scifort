# Erlang censoring validation, 2026-11-15

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
References: SciPy 1.15.3 and mpmath 1.3.0.

## Native checks

All 175 checks in `test/test_censored_erlang.f90` passed at both `-O0` and
`-O3`, each with `-g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace`.
They cover mixed observations, default and keyword arguments, negative
likelihoods, positive noninteger shapes, gamma and exponential identities,
scalar and array agreement, shape-two analytic survival and interval masses,
adjacent and subnormal intervals, empty/invalid inputs, support densities,
and intervals with infinite endpoints. Fits cover fixed shape, continuous
shape search, gamma parameter/objective/diagnostic agreement, array-call
compatibility, optional iteration controls, and rejected free locations,
empty samples, and exact support endpoints.

The independent mixed sample has exact `[0.5,1.2,2.5]`, left `[0.3]`, right
`[3.5]`, and intervals `[1.1,1.4]`, `[4.0,4.4]`. At shape 2.3, location 0,
scale 1.3, its 100-digit mpmath log-likelihood is
`-16.30517655979195826580375891409484261436`. SciPy's gamma fit at fixed
location gives shape `0.9927593350620091` and scale `2.232420387730313`,
checking the Erlang/gamma continuous-shape convention. An independent
shape-one exact-data MLE is scale 2 for observations `[1,2,3]`.

Tolerances are absolute `2e-14` for ordinary scalar/interval identities,
`6e-14` for the mixed reference, `2e-13` for the far-tail survival,
and `5e-13` for the subnormal interval log mass. Parent-kernel and fit
diagnostic agreement are exact. Fitted parameters allow `4e-4` for different
optimizers; the analytic scale MLE allows `3e-6`.

The new Erlang documentation example compiled and ran with the checked flags.
All products were placed under ignored `build/` directories. To reproduce
native checks from CMD with the local FPM dependency paths:

```cmd
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Jbuild\erlang-native -Ibuild\erlang-native -Ibuild\gfortran_2806A042DAA8BBF7 src\special\log_gamma.f90 src\empirical_distribution.f90 src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored_erlang.f90 build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\erlang-native\test.exe
build\erlang-native\test.exe
```

Repeat with `-O3` and a separate output/module directory for the optimized
configuration. The command explicitly compiles changed modules before linking
existing FPM dependencies; directory hashes depend on the local build.

## Reference and integration checks

The reference harness adds 100 mixed Erlang samples with positive real shapes,
varied scales/locations, tail intervals, and adjacent endpoints. Erlang shares
an independently evaluated mpmath gamma specification, rather than calling
production kernels. The driver exercises the new Erlang entry point.

```cmd
python validation\check_censored.py --module-dir build\gfortran_2806A042DAA8BBF7 --library build\gfortran_3C159EB802C22FB9\scifort\libscifort.a
fpm test
python tools\check_test_inputs.py
```

All 1,721 reference comparisons passed in each checked and optimized build.
The new Erlang maximum scaled log-likelihood error was `7.51e-14`; the overall
maximum was `1.08e-12` in an existing gamma stress case. The harness uses
`abs(actual-reference)/max(1,abs(reference)) <= 2e-11` and seed 20261112.
Its flags are `-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace`
and `-O3 -fcheck=bounds`. The console-input and source/inventory/link audits
also passed. The full `fpm test` suite passed with exit code 0 in 3021.693 seconds
(50 minutes 22 seconds). No numerical changes followed this full run.

## Limits

The wrappers add no new numerical approximation or optimizer. They inherit
existing gamma kernel limitations. Censored fits require fixed location and
exclude exact support endpoints; likelihoods retain their mathematical endpoint
values. Shape fitting is continuous, not integer constrained. Fits are bounded
local optima, without a guarantee of a unique or global MLE. Existing array
APIs, C ABI and runtime dependencies are unchanged.
