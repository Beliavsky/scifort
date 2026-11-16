# Pareto censoring validation, 2026-11-15

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
Reference libraries: SciPy 1.15.3 and mpmath 1.3.0.

## Native checks

All 354 checks in `test/test_censored_pareto.f90` passed at both `-O0` and
`-O3`, each with `-g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace`.
Checks cover defaults, scalar/array agreement, exponential transformation with
exact-density Jacobian, the log-Laplace density identity, mixed observations,
NNLF, free/fixed fit parameters, objective agreement, optimizer controls,
empty/invalid inputs, support endpoints, clipped/impossible intervals,
infinite endpoints, and equal endpoints as exact observations.

Stress checks include adjacent tail intervals, underflowed cumulative hazard
with the smallest positive binary64 shape, a nearly unit interval mass,
overflowing standardized ratios and raw displacements, and unrepresentable
log tails. The nearly unit probability check uses a ratio to reject zero;
it caught and led to removal of a redundant hazard log/exp round trip.

The mixed sample has exact `[0.5,1.2,2.5]`, left `[0.3]`, right `[3.5]`, and
intervals `[1.1,1.4]`, `[4.0,4.4]`. At shape 1.7, location 0, scale 0.2,
its 100-digit mpmath log-likelihood is
`-24.25414363441046155039288912405825430402`.
With location fixed at zero and scale at 0.2, SciPy's fitted shape is
`0.45595703124999953`. The default SciPy fit with free scale returned an
infeasible candidate that excluded exact observations; this is not used as a
correctness reference. SciFort preserves its existing policy of rejecting
impossible observations rather than applying SciPy's finite fitting penalty.
For exact `[1,2,4]`, fixed location zero and scale one, the analytic shape
MLE is `3/sum(log(x)) = 1/log(2)`. With shape fixed at two and location zero,
the fitted scale is the sample minimum, one (the included support endpoint).

Ordinary scalar, identity and mixed-reference absolute tolerances are
`2e-14` to `6e-14`; extreme log values allow up to `4e-13`.
The near-unit probability ratio allows `3e-14`. The SciPy fitted shape allows
`2e-4` for different optimizers; analytic shape/scale MLEs allow `4e-6`/`3e-6`.

Native checked command (CMD, with matching local FPM dependencies):

```cmd
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Jbuild\pareto-native -Ibuild\pareto-native -Ibuild\gfortran_2806A042DAA8BBF7 src\special\log_gamma.f90 src\empirical_distribution.f90 src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored_pareto.f90 build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\pareto-native\test.exe
build\pareto-native\test.exe
```

The optimized run used `-O3` and separate `build\pareto-optimized` module
and executable paths. Changed modules were compiled explicitly before linking
existing FPM dependencies. Directory hashes are local to this build.

## Reference and integration checks

All 2,388 reference comparisons passed in checked and optimized builds.
The 112 new Pareto cases comprise 100 randomized mixed samples, nine interval
stress cases, and three exact/tail extremes. Independent mpmath density/CDF/SF
definitions and log probability differences use up to 800 digits, with seed
20261112. Stress cases include adjacent/support-crossing intervals at large
shape, extreme ratios and displacements, support clipping, nearly unit mass,
smallest positive shape, and underflowed cumulative hazard.

The new family's maximum scaled error was `4.93e-15`; the overall maximum
remained `1.08e-12` in an existing gamma stress case. The metric is
`abs(actual-reference)/max(1,abs(reference))`, with threshold `2e-11`.
Configurations are `-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace`
and `-O3 -fcheck=bounds`.

```cmd
python tools\generate_censored_reference.py
python validation\check_censored.py --module-dir build\gfortran_2806A042DAA8BBF7 --library build\gfortran_3C159EB802C22FB9\scifort\libscifort.a
fpm test
python tools\check_test_inputs.py
```

The standalone generator exited 0 and reproduced the mixed reference and
fixed-support SciPy shape fit. The guide example compiled and ran with the
checked flags. The full `fpm test` suite passed with exit code 0 in
2002.195 seconds (33 minutes 22 seconds). No numerical changes
followed that run. The unattended-console-input check, source ASCII and line
lengths, Python syntax, inventory counts, changed local links, and
`git diff --check` passed.

## Limits

Fits use the existing bounded local search without uniqueness or global-MLE
guarantees. Support-changing location/scale bounds need care. Scalar/array
routines and C ABI are unchanged; the new censored kernels' extreme-range
behavior does not extend the existing scalar routines. No dependency is added.
