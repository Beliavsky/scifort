# Basin-hopping validation — 2026-11-10

Modules: `scifort_optimize_basinhopping` (engine) and
`scifort_optimize_basinhopping_minimize` (SciFort minimizer adapter).

- Standalone GNU Fortran checked build: `-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow`.
- Optimized GNU Fortran build: `-O3 -march=native`.
- Native checks cover strict uphill rejection at zero temperature, forced accept/reject, reproducible seeds, local evaluation reporting, no-improvement early stopping, callback stopping, custom steps, adaptation and invalid inputs.
- Independent numerical reference: 6 functions (1-D tilted double well, 2-D Rastrigin, 2-D Himmelblau, 2-D Ackley, 1-D periodic multimodal, 2-D tilted double well), 8 independent seeds, 120 moves each, against known optima and separate SciPy `basinhopping` runs using Nelder–Mead.
- Fortran and SciPy random streams differ, so objective accuracy is compared rather than exact trajectories.
- The wrapper's *compiler contract* is tested separately against a clearly labeled mock `scifort_optimize_minimize` module. This does not establish numerical integration with the real Nelder-Mead, Powell, and BFGS code. That requires FPM with its external dependencies.

Use `scripts/validate_optimize_basinhopping.sh` to reproduce the counts.
