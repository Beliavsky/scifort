# Censored beta validation, 2026-11-12

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
Reference tools: SciPy 1.15.3 and mpmath 1.3.0.

All 307 native checks passed in checked (`-O0`) and optimized (`-O3`) builds
with runtime checks and traps for invalid arithmetic, division by zero, and
overflow. All 622 reference comparisons passed in each reference build; the
maximum scaled log-likelihood error was `1.08e-12` (a gamma stress case).

The native regression covers all six censored families.
Beta additions cover a fixed-support fit compared with SciPy, ordinary array
API agreement, the uniform and reflection identities, both singular support
endpoints, undefined sums, fixed-support validation, empty/invalid samples,
infinite/equal interval bounds, monotonicity, adjacent and edge-touching
intervals, independent endpoint normalization, compensated upper support sums,
underflowed coordinates, tiny/huge shape combinations, and large shapes at the mode.

Checked native command (create the output directory first):

```cmd
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Jbuild\beta-native -Ibuild\beta-native -Ibuild\gfortran_2806A042DAA8BBF7 src\special\log_gamma.f90 src\empirical_distribution.f90 src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored.f90 build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\beta-native\test_censored.exe
build\beta-native\test_censored.exe
```

The optimized command replaces `-O0` with `-O3` and the output/module directory
with `build\beta-optimized`, retaining runtime checks and arithmetic traps.
All changed numerical modules are compiled explicitly, including `log_gamma`;
the archive supplies unchanged dependencies. FPM directory hashes are local
artifacts and must match the module directory and archive on other builds.

Reference commands:

```cmd
python tools\generate_censored_reference.py
python validation\check_censored.py --module-dir build\gfortran_2806A042DAA8BBF7 --library build\gfortran_3C159EB802C22FB9\scifort\libscifort.a
```

The harness uses seed 20261112 and 100 mixed cases per family, plus one
Weibull, three lognormal, nine gamma, and nine beta stress cases (622 total).
Beta stress cases cover both edges, adjacent endpoints, scaled underflow,
a tiny shape, large shapes, and the compensated upper-support edge. Random
references use 100 digits; beta stress and fixed extreme constants use 700.
The tiny-shape/huge-shape SF reference uses the leading defining-integral term
with omitted relative correction below 1e-291. Reference builds use
`-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace` and
`-O3 -fcheck=bounds`. The changed shared gamma module is included in both.

## Shared-helper regression

The extreme-shape checked run exposed intermediate overflow in the compiler's
integer-power implementation for `x**negative_integer` inside the shared
Stirling remainder difference. GDB located the trap in `__powidf2` called
from `src/special/log_gamma.f90`. The one-line fix forms positive powers of
`1/x` instead. Coefficients and the mathematical series are unchanged.
The beta case with shapes 1e-20 and 1e6 under arithmetic traps exercises this
helper. A separate 1e-20/1e308 case exercises the bounded Stirling limit.

## Limits

Only this compiler/platform was checked. Existing incomplete-beta, log-beta,
and gamma-helper shape/range limits still apply. The initial censored fit
requires fixed location/scale and excludes exact support-endpoint samples;
local convergence does not prove a unique/global MLE or identifiability.
The rounded upper-support sum must be finite. Underflowed coordinates return
NaN if the far-shape/coordinate bound does not permit a controlled leading-term
approximation. This failure path is tested. Extreme-tail results can improve
on the existing scalar beta API because both endpoint distances are preserved.
No new dependencies, observation weights, or C ABI symbols are introduced.

## Full suite and repository checks

The complete FPM build and test execution passed through CMD with exit code 0.
The censored-fitting test reported all 307 checks passed in the full run.
Command used from CMD:

```cmd
fpm test > build\beta-full-test-cmd.log 2>&1
```

`python tools/check_test_inputs.py` passed. Changed Fortran sources passed ASCII
and 132-column checks, affected local documentation links resolved, and
`git diff --check` passed. Build outputs, diagnostic reproducers, and reference
JSON remain under ignored `build/`; no generated binaries or downloaded sources
were added to the change.
