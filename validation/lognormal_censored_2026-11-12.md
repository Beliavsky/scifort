# Censored lognormal validation, 2026-11-12

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
Reference tools: SciPy 1.15.3 and mpmath 1.3.0.

All 145 native checks passed in both checked (`-O0`) and optimized (`-O3`)
builds. All 404 independent likelihood comparisons passed in each reference
build (`-O0` with runtime traps, and `-O3 -fcheck=bounds`). The maximum error
divided by `max(1,abs(reference))` was `1.32e-14`.

The native regression test covers all four censored families.
Lognormal additions cover mixed exact/left/right/interval observations,
fixed-location fits compared with SciPy, exact-only analytic MLEs, agreement
with the existing uncensored API, support boundaries, infinite and equal
interval endpoints, adjacent-float and extreme-tail intervals, empty samples,
NaN/infinite parameters, invalid/free location bounds, and tiny/subnormal shape.

Checked native command (create the output directory first):

```cmd
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Jbuild\lognormal-native -Ibuild\lognormal-native -Ibuild\gfortran_2806A042DAA8BBF7 src\empirical_distribution.f90 src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored.f90 build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\lognormal-native\test_censored.exe
build\lognormal-native\test_censored.exe
```

The optimized native command replaces `-O0` with `-O3` and the output/module
directory with `build\lognormal-optimized`, retaining all runtime checks.
Both compile the changed modules explicitly and link unchanged dependencies
from the existing FPM archive. Directory hashes are local artifacts; other
builds must supply their matching module directory and library.

Independent reference commands:

```cmd
python tools\generate_censored_reference.py
python validation\check_censored.py --module-dir build\gfortran_2806A042DAA8BBF7 --library build\gfortran_3C159EB802C22FB9\scifort\libscifort.a
```

The reference harness uses seed 20261112 and 100 mixed cases per family, plus
one large-shape Weibull and three lognormal stress cases. The lognormal stress
cases include adjacent tiny endpoints with huge scale, adjacent huge endpoints,
and large shifted distances. Defining formulas are evaluated at 100 digits;
fixed native constants use up to 700 digits. The checker selects the better
conditioned reference tail and rejects nonfinite expected values.

Limitations: only this compiler/platform was checked. Bounded local convergence
does not prove a unique/global MLE or identifiability. Censored lognormal fitting
requires fixed location; no observation weights or new C ABI symbols are added.

## Full suite and repository checks

The complete `fpm test` build and execution passed. A final run through CMD
returned exit code 0, including the new 145-check regression and existing
likelihood, fitting, Krylov refresh, and C ABI tests. Command used from CMD:

```cmd
fpm test > build\lognormal-full-test-cmd.log 2>&1
```

The first full run was launched with PowerShell redirection, which treated
an existing `Integration successful.` stderr diagnostic as a native command
error. Its output reported passes, but its shell status was ambiguous. The
CMD run above verifies the actual successful exit status.

`python tools/check_test_inputs.py` passed. Changed Fortran files passed ASCII
and 132-column checks, affected documentation local links resolved, and
`git diff --check` passed. No build artifacts or downloaded sources were added.
