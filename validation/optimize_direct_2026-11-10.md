# DIRECT implementation and numerical validation — 2026-11-10

## Added

- `src/optimize/direct.f90`: independent MIT implementation of deterministic
  rectangle-based global optimization, with both original and locally biased
  rectangle-size selections.
- Exported `direct`, `direct_result`, callback interfaces, and status constants
  through `src/optimize.f90`.
- API guide, executable example, native Fortran assertions, and SciPy reference
  comparison runner.

## Reproducible validation

```sh
bash scripts/validate_optimize_direct.sh
```

| Check | GNU Fortran checked | GNU Fortran `-O3` |
|---|---:|---:|
| Native assertions | 109 passed | 109 passed |
| Independent SciPy/NumPy comparisons | 60 passed | 60 passed |
| High-budget known-optimum targets | 18/18 | 18/18 |
| Maximum objective recomputation difference | 0 | 0 |
| Rastrigin executable example | Passed | Passed |
| Previous basin-hopping regression script | Passed | Passed |
| Previous dual-annealing regression script | Passed | Passed |

Checked compilation: GNU Fortran 14.2.0 with `-std=f2008 -O0 -g
-Wall -Wextra -fcheck=all -ffpe-trap=invalid,zero,overflow`.
Optimized compilation: GNU Fortran 14.2.0 with `-std=f2008 -O3
-Wall -Wextra`. Independent comparisons: SciPy 1.17.0 and NumPy.

Reference cases cover 1-D through 4-D problems; Rastrigin, Ackley,
Rosenbrock, Himmelblau, quadratic, nonsmooth absolute deviations, and a
periodic multimodal function. Each case uses 120, 480, and 1,600 evaluation
budgets and both rectangle-size conventions; the test independently recalculates
the returned objective and checks variable bounds and strict function-call
limits. The SciPy reference implementation has different rectangle-selection
and approximate evaluation-budget semantics, so solutions are compared by
mathematical accuracy and not by their exact sampling histories.

Native tests cover reproducibility, shifted optima, invalid/nonfinite inputs,
callback termination, fixed evaluation budgets, and all standard stopping paths.

**Limitations:** The implementation is standalone (no SciPy backend code was
copied), and the locally biased metric is DIRECT-L-style rather than a bitwise
port of SciPy's implementation. No sparse/integrality/general constraints
are supported. FPM and its Git dependencies could not be downloaded in this
offline environment, so this release's full FPM integration remains untested.
