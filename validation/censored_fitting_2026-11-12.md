# Censored fitting validation, 2026-11-12

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
Reference tools: SciPy 1.15.3 and mpmath 1.3.0.

## Focused numerical checks

- 91 native checks passed with `-O0` and with `-O3`, each with runtime checks
  and traps for invalid arithmetic, division by zero, and overflow.
- The tests include analytic fixed-location exponential estimates, mixed
  censored fits compared with SciPy, invalid samples and parameters, support
  boundaries, infinite interval endpoints, underflow/overflow-prone estimates,
  adjacent-float intervals, and compatibility with uncensored calls.
- 301 independent likelihood comparisons passed per build in the randomized
  reference harness: 100 mixed-data cases per family and a large-shape Weibull
  stress case. The maximum error divided by `max(1,abs(reference))` was
  approximately `1.34e-14`, against 100-digit mpmath calculations.
- Fixed native reference values use 100 to 700 decimal digits. SciPy fit
  comparisons use parameter-specific tolerances and do not require matching
  iteration counts or optimizer trajectories.

Commands used for the checked native test (from the repository root):

```cmd
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Jbuild\censored_checks -Ibuild\censored_checks -Ibuild\gfortran_2806A042DAA8BBF7 src\empirical_distribution.f90 src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored.f90 build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\censored_checks\test_censored.exe
build\censored_checks\test_censored.exe
```

The optimized native check used the same command with `-O3` replacing `-O0`
and `build\censored_optimized` replacing `build\censored_checks`.
Both directories were created before compilation. The archive supplies the
unchanged library dependencies; the changed modules are recompiled explicitly
at the stated optimization level into isolated module directories.

```cmd
python tools\generate_censored_reference.py
python validation\check_censored.py --module-dir build\gfortran_2806A042DAA8BBF7 --library build\gfortran_3C159EB802C22FB9\scifort\libscifort.a
```

The randomized harness compiles its changed modules and input-driven driver
with `-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace` and with
`-O3 -fcheck=bounds`. The recorded FPM directory hashes are local build
artifacts; supply the matching module directory and archive for a new build.

## Repository checks

The complete `fpm test` build and run passed with exit code 0, including the
new 91-check test and existing fitting, likelihood, ECDF, logrank, and C ABI
regressions. Command used:

```cmd
fpm test > build\censored-full-test.log 2>&1
```

`python tools/check_test_inputs.py` passed. Changed Fortran sources were
checked for ASCII and the 132-column free-form limit, documentation local
links were checked, and `git diff --check` passed.

## Limits

These checks establish behavior for the tested cases, not complete equivalence
to SciPy or a proof of global fitting convergence. Censoring support is limited
to the documented families, and observation weights are not implemented.
No additional compiler or platform was tested for this change.
