# Censored gamma validation, 2026-11-12

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
Reference tools: SciPy 1.15.3 and mpmath 1.3.0.

All 216 native checks passed in both `-O0` and `-O3` builds with bounds/runtime
checks and traps for invalid arithmetic, division by zero, and overflow.

All 513 independent comparisons passed in each reference build, with maximum
error divided by `max(1,abs(reference))` of approximately `1.08e-12`.

The native regression covers all five censored families.
Gamma additions cover mixed observations, a fixed-location fit compared with
SciPy, an analytic fixed-shape scale estimate, array API agreement, support
endpoints, undefined sums of opposing infinities, narrow and wider intervals,
left/right tail and interval monotonicity, invalid inputs, infinite/equal
interval bounds, tiny scaled arguments, and huge-shape underflow saturation.
Gamma with shape one uses exactly the exponential censored likelihood.

Checked native command (create the output directory first):

```cmd
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Jbuild\gamma-native -Ibuild\gamma-native -Ibuild\gfortran_2806A042DAA8BBF7 src\empirical_distribution.f90 src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored.f90 build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\gamma-native\test_censored.exe
build\gamma-native\test_censored.exe
```

The optimized command replaces `-O0` with `-O3` and the output/module directory
with `build\gamma-optimized`, retaining runtime checks and arithmetic traps.
Changed modules are compiled explicitly; the archive supplies unchanged
dependencies. FPM directory hashes are local artifacts and must be replaced
with the matching module directory and archive on other builds.

Independent reference commands:

```cmd
python tools\generate_censored_reference.py
python validation\check_censored.py --module-dir build\gfortran_2806A042DAA8BBF7 --library build\gfortran_3C159EB802C22FB9\scifort\libscifort.a
```

The randomized harness uses seed 20261112 and 100 mixed cases per family,
plus one Weibull, three lognormal, and nine gamma stress cases (513 total).
Gamma stress cases include scaled underflow, adjacent endpoints, extreme
upper tails, shapes 1000 and 1e6 near their modes, and tiny shape 1e-20.
The defining formulas use 100-digit mpmath calculations; native fixed constants
use up to 700 digits. Reference checks run with `-O0 -g -fcheck=all
-ffpe-trap=invalid,zero,overflow -fbacktrace` and `-O3 -fcheck=bounds`.

Limitations: this validates tested cases, not all gamma parameter regimes or
SciPy parity. Only this compiler/platform was checked. Existing incomplete-gamma
kernel limitations, especially very large shapes, remain applicable. Censored
gamma fits require fixed location and exclude exact support-endpoint samples;
local convergence does not establish identifiability or global optimality.
No new dependencies, observation weights, or C ABI symbols are introduced.

## Full suite and repository checks

The complete FPM build and test execution passed through CMD with exit code 0.
The new censored-fitting test reported all 216 checks passed in the full run.
Command used from CMD:

```cmd
fpm test > build\gamma-full-test-cmd.log 2>&1
```

`python tools/check_test_inputs.py` passed. Changed Fortran sources passed ASCII
and 132-column checks, affected local documentation links resolved, and
`git diff --check` passed. Build outputs and reference JSON remain under ignored
`build/`; no generated binaries or downloaded sources were added.
