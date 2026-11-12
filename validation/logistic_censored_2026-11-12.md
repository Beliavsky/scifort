# Censored logistic validation, 2026-11-12

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
Reference tools: SciPy 1.15.3 and mpmath 1.3.0.

All 373 native checks passed at both `-O0` and `-O3`, with bounds/runtime
checks and traps for invalid arithmetic, division by zero, and overflow.
All 729 high-precision reference comparisons passed in each reference build;
the maximum scaled log-likelihood error was `1.08e-12` (a gamma stress case).

The native regression covers all seven censored families. Logistic additions
cover a mixed-data fit compared with SciPy, fixed parameters, objective
agreement, array API compatibility, reflection and symmetric-mass identities,
monotonic interval masses, subnormal and underflowed standardized widths,
adjacent extreme-tail intervals, overflowing raw widths and displacements,
rounded subtraction-overflow thresholds, maximum finite standardized values,
and zero displacement at the smallest positive scale,
unrepresentable log tails, infinite/equal endpoints, empty data, invalid
parameters/data, and fit error statuses. NaN parameters are rejected before
ordered comparisons, including with arithmetic traps enabled.

Checked native command (create `build\logistic-native` first):

```cmd
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Jbuild\logistic-native -Ibuild\logistic-native -Ibuild\gfortran_2806A042DAA8BBF7 src\special\log_gamma.f90 src\empirical_distribution.f90 src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored.f90 build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\logistic-native\test_censored.exe
build\logistic-native\test_censored.exe
```

The optimized command replaces `-O0` with `-O3` and `logistic-native` with
`logistic-optimized`, retaining runtime checks and arithmetic traps.
FPM directory hashes are local artifacts; supply matching module and archive
paths for other builds. The modules listed above are compiled explicitly;
the archive supplies unchanged dependencies.

Reference commands:

```cmd
python tools\generate_censored_reference.py
python validation\check_censored.py --module-dir build\gfortran_2806A042DAA8BBF7 --library build\gfortran_3C159EB802C22FB9\scifort\libscifort.a
```

The harness uses seed 20261112 and 100 mixed cases per family, plus one
Weibull, three lognormal, nine gamma, nine beta, and seven logistic stress
cases (729 total). Logistic stress cases cover subnormal widths, scaled
underflow, both far tails, overflowing raw width, overflowing displacement
with a finite standardized value, and tiny scale. Random references use
100 digits; beta and logistic stress references use 700 digits.
Logistic reference intervals use high-precision tail subtraction independently
of the product identity used in Fortran. Reference builds use
`-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace` and
`-O3 -fcheck=bounds`.

The mixed-data SciPy fit gives location `1.847711699978416` and scale
`1.12816063694188`; parameter comparisons allow `3e-4` for optimizer differences.
The independent mixed log-likelihood at location 0.4 and scale 1.3 is
`-15.999635359473038`.

## Full-suite verification

The full `fpm test` run passed with exit code 0, including all 373 censored
checks. The exact shell invocation used to capture compiler/runtime output
and preserve the FPM exit status was:

```cmd
cmd /d /c "fpm test > build\logistic-full-test.log 2>&1"
```

The existing interpolation array-temporary runtime warnings were emitted;
there were no test failures. `python tools\check_test_inputs.py` passed.
Final whitespace, changed Fortran ASCII/line lengths, selected documentation
links, and the diff's file list were checked. No generated build artifact or
unrelated file is included in the change.

## Limits

Fits use the existing bounded local search. Convergence does not establish
a global optimum, identifiability, or existence of an unconstrained finite
MLE. No weights, new optimizer, dependency, or C ABI symbols are added.
Censored logistic log magnitudes outside binary64 range return negative
infinity; the existing likelihood accumulation range limits still apply.
