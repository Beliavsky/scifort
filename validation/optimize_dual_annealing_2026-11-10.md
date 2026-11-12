# Dual annealing validation — November 10, 2026

- Core: `src/optimize/dual_annealing.f90`, independent generalized simulated
  annealing with Tsallis visiting/acceptance, cooling, restarts and callbacks.
- Adapter: `src/optimize/dual_annealing_minimize.f90`, SciFort bounded L-BFGS-B
  and Powell numerical methods (interfaces built with a test-only mock).
- Checked mode: `gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow`.
- Optimized mode: `gfortran -O3 -march=native`.
- Independent reference: SciPy 1.17.0, NumPy 2.3.5; eight objective functions,
  six independent seeds, 360 steps per case, no local minimization in either
  code. Benchmark functions: shifted quadratic (2/3-D), Ackley, Rastrigin,
  Rosenbrock, Griewank, Himmelblau and a 1-D periodic multimodal objective.
- Core native tests include reproducible seeds, ordinary and limiting visiting
  parameters, best tracking, callbacks, function budget, restart accounting,
  user-supplied minimization, nonfinite objectives, and invalid inputs.
- Local adapter tests verify public Fortran signatures and default bounded
  local-method selection. Mock testing is **not** a validation of actual SciFort
  L-BFGS-B or Powell numerical integration.
- No FPM installation in the standalone test environment, so the full FPM
  regression suite and actual pinned L-BFGS-B dependency remain untested.

Reproduce with `./scripts/validate_optimize_dual_annealing.sh`.
