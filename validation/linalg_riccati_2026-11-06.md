# Validation: real/complex continuous/discrete algebraic Riccati solvers

Date: 2026-11-06. gfortran 14.2.0, NumPy 2.3.5, SciPy 1.17.0.

## Performed checks

- Native `python tools/linalg_reference/run_native_linalg.py --mode both`:
  **2/2 native test programs** passed in each mode (new Riccati regression
  and prior Sylvester/Lyapunov regression). The new native test checks both
  standard and generalized real/complex cases, exact scalar CARE/DARE
  solutions (stable and unstable open loop), Hermitian output, singular
  state cost, nonstabilizing zero-cost discrete models, zero/nonpositive
  iteration budgets, shape errors, indefinite `R`, singular `E`, non-Hermitian
  `Q`, nonfinite input, and a 20-state boundary case.
- Independent `compare_riccati_scipy.py --mode both --cases 320
  --seed 20261106 --max-n 6`: **320/320** passed in each mode;
  maximum relative solution error `2.336e-10` and maximum normalized
  residual `1.175e-11` across the two modes.
- Independent `--mode both --cases 160 --seed 20261108 --max-n 10`:
  **160/160** passed in each mode; maximum relative solution error
  `4.051e-11`, maximum normalized residual `1.318e-11`.
- Independent `--mode both --cases 160 --seed 20261109 --max-n 6`:
  **160/160** passed in each mode; maximum relative solution error
  `5.853e-12`, maximum normalized residual `6.331e-14`.
- Total: **640 distinct randomized scenarios in each compiler mode**
  (1,280 checked/optimized executions), in addition to the native tests.
  Every scenario uses a direct SciPy QZ result plus an independently evaluated
  generalized Riccati residual. Tests cover real and complex solutions,
  stable/unstable open-loop dynamics, 1–10 states, multiple controls, and
  optional descriptor `E` and cross cost `S`.
- Example program compiles and runs under `-O3` and reports small residuals.

The checked build used `-O0 -g -fcheck=all -fbacktrace`, optimized build
`-O3`. No SciPy sources are included in the Fortran implementation.

## Reference accuracy and limitations

In preliminary tests SciPy 1.17's *default* `balanced=True` sometimes gave
non-negligible residuals for random generalized equations even though
`balanced=False` and SciFort satisfied the equation to near roundoff.
For those comparisons the harness uses SciPy `balanced=False` and
independently checks every reference's equation residual before using it.
The three final runs above reported **zero discarded SciPy references**.
These tests do not establish general superiority of one method over SciPy.

The implementation uses a matrix-sign method for CARE and structured doubling
for DARE, not SciPy's Schur/QZ method. The state dimension is limited to 20;
poorly conditioned, nearly marginally stable, unobservable or uncontrollable
cases may fail. The test uses a finite relative tolerance and numerical
stabilizability/detectability assumptions, and is not a proof that every
stabilizing solution exists or was selected in a degenerate problem.

FPM is unavailable here, so the full project-wide `fpm test` suite was **not**
run. No pre-existing sparse modules or algorithms were modified.
