# SciFort constrained differential evolution — validation

## Scope

- Integer variables with fractional box endpoints and rejection of empty integer lattices.
- Linear inequality/equality bounds and nonlinear constraint callback values.
- Combined mixed-integer/linear/nonlinear feasible problems and wholly infeasible systems.
- Feasibility-first selection, separate objective/constraint counts, local seed reproducibility.
- Optional polishing procedure with immutable integer coordinates.

## Build and numerical tests

```sh
python tools/run_differential_evolution_tests.py --mode both
```

The complete runner builds with `gfortran -std=f2018 -Wall -Wextra -Werror`,
first with `-O0 -fcheck=all -fbacktrace` and then with `-O3 -march=native`.
Native tests include the previous unconstrained DE tests, and the new extended
constraint tests.

The comparison harness compiles a separate Fortran driver and compares its
results against SciPy `differential_evolution(polish=False)` for 5 benchmark
problems x 4 strategies x 4 seeds = 80 cases **per build**. Objective quality
is checked against independently known optima, and constraint feasibility is
computed independently, not merely accepted from SciFort's status flags.
Infeasible test cases are required to return an infeasible result without
calling the objective. Exact stochastic trajectory agreement is not expected.

In the initial run with GNU Fortran 14.2.0 and SciPy 1.17.0:

- Legacy 104 native checks per build: passed.
- New 42 native checks per build: passed.
- 72 legacy SciPy benchmark cases per build: passed.
- 80 constrained/mixed-integer SciPy cases per build: passed.
- Worst objective difference from SciPy among new cases: checked 3.38e-10,
  optimized 4.32e-10. The cases are distinct-seed stochastic comparisons.
- Maximum excess over known optimum among new cases: checked 2.48e-11,
  optimized 3.86e-11.

The small `de_polisher` callback is tested by the native test, and the separate
L-BFGS-B bridge was **compile-checked with a test double** for the optimizer
interface. The bridge has **not** been exercised against the actual pinned
L-BFGS-B dependency in this environment (FPM is unavailable).

## Limitations

- The feasibility criterion is an absolute `constraint_tol` (default 1e-10),
  not SciPy's exact acceptance rule; multiple infeasible constraints are ranked
  by sum of component violations.
- Exact equality constraints have a tiny feasible region and may require many
  generations, or may not be found.
- Constrained local polishing and `trust-constr` integration are not implemented.
- No Sobol/Halton DE population initialization, parallel evaluation, or
  arbitrary batched/vectorized objective calls.
- Full `fpm test` not run because FPM is unavailable.
