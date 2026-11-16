# Lomax censoring validation, 2026-11-15

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
Reference libraries: SciPy 1.15.3 and mpmath 1.3.0.

## Native checks

All 355 checks in `test/test_censored_lomax.f90` passed at both `-O0` and
`-O3`, each with `-g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace`.
Checks cover defaults, scalar/array agreement, the shifted-Pareto density
identity, exponential transformation with exact-density Jacobian, mixed
observations, NNLF, free/fixed fits, objective agreement, optimizer controls,
empty/invalid inputs, support endpoints, clipped/impossible intervals,
infinite endpoints, and equal endpoints as exact observations.

Stress checks cover adjacent tail intervals, underflowed cumulative hazards,
nearly unit interval probabilities, overflowing standardized ratios and raw
displacements, unrepresentable log tails, and compensating huge shapes with
underflowed or subnormal displacements. Ratio checks reject a falsely rounded
zero and detect relative errors in the subnormal/huge-shape case.

The mixed sample has exact `[0.5,1.2,2.5]`, left `[0.3]`, right `[3.5]`, and
intervals `[1.1,1.4]`, `[4.0,4.4]`. At shape 1.7, location zero, scale 1.3,
its mpmath log-likelihood is
`-15.45809456699385922829247285612573042901`.
With location fixed at zero and scale at 1.3, SciPy's fitted shape is
`1.0771484375000002`. A free-scale SciPy fit approached the exponential limit
with shape and scale growing very large, so the fixed-scale fit is used as
the parameter reference. For exact `[0,1,3]`, location zero and scale one,
the analytic shape MLE is `3/sum(log1p(x)) = 1/log(2)`.

Ordinary scalar, identity and mixed-reference absolute tolerances range from
`2e-14` to `6e-14`; extreme log values allow up to `4e-13`. The near-unit
probability ratio allows `3e-14`; compensated/subnormal hazard ratios allow
`3e-13`. SciPy fitted shape allows `2e-4` for different optimizers, and the
analytic MLE allows `4e-6`. The guide example compiled and ran with checked flags.

Native checked command (CMD, with matching local FPM dependencies):

```cmd
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Jbuild\lomax-native -Ibuild\lomax-native -Ibuild\gfortran_2806A042DAA8BBF7 src\special\log_gamma.f90 src\empirical_distribution.f90 src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored_lomax.f90 build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\lomax-native\test.exe
build\lomax-native\test.exe
```

The optimized run used `-O3` and separate `build\lomax-optimized` module and
executable paths. Changed modules were compiled explicitly before linking
existing FPM dependencies. Directory hashes are local to this build.

## Reference and integration checks

All 2,503 reference comparisons passed in checked and optimized builds.
The independent harness adds 115 Lomax cases: 100 randomized mixed samples,
ten interval stress cases, and five exact/tail extremes. The mpmath oracle
uses direct density/CDF/SF definitions and log probability differences at up
to 800 digits, with seed 20261112. Cases include subnormal/underflowed widths,
smallest positive shape, nearly unit mass, support clipping, overflowing
ratios/displacements, and huge shape compensating tiny displacement.

The new family's maximum scaled error was `2.39e-15`; the overall maximum
remained `1.08e-12` in an existing gamma stress case.
The harness uses `abs(actual-reference)/max(1,abs(reference)) <= 2e-11`.
Configurations are `-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace`
and `-O3 -fcheck=bounds`. Native ratio checks supplement the harness's
absolute tolerance for tiny log probabilities.

```cmd
python tools\generate_censored_reference.py
python validation\check_censored.py --module-dir build\gfortran_2806A042DAA8BBF7 --library build\gfortran_3C159EB802C22FB9\scifort\libscifort.a
fpm test
python tools\check_test_inputs.py
```

The standalone generator exited 0 and reproduced the mixed reference and
fixed-scale SciPy shape fit. The full `fpm test` suite passed with exit code 0
in 1407.566 seconds (23 minutes 28 seconds). No numerical changes
followed that successful run. The unattended-console-input check, source ASCII
and line lengths, Python syntax, inventory counts, changed local links, and
`git diff --check` passed.

## Limits

Fits use the existing bounded local search without uniqueness or global-MLE
guarantees. Shape/scale can approach an exponential limit; free location can
approach data endpoints. Choose scientifically justified bounds or fixed
parameters. Scalar/array routines and C ABI are unchanged; the censored
kernels' extreme-range behavior does not extend existing scalar routines.
No dependency is added.
