# Moyal, Levy, and reflected Levy censoring validation, 2026-11-16

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
Reference libraries: SciPy 1.15.3 and mpmath 1.3.0.

## Native numerical and API checks

All 1,068 checks in `test/test_censored_normal_transform.f90` passed at both
`-O0` and `-O3`, each with
`-g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace`.
Checks cover mixed references and NNLF, ordinary scalar density/CDF/SF
agreement at three scales, CDF monotonicity/range, CDF/PPF round trips,
original array fit interfaces, free/fixed bounded fits, SciPy fixed-location
scale estimates, optimizer controls, invalid/empty inputs and bounds,
impossible objectives, support clipping and endpoint densities/tails,
infinite interval endpoints, and equal endpoints as exact measurements.

Transformation identities cover inverse-square and negative-log half-normal
intervals, a density Jacobian, and mixed Levy/reflected-Levy censoring with
reflected location. The mixed reflection check allows ordinary summation
roundoff because left and right contribution groups exchange order.
Extreme tests cover the smallest positive interval width, overflowing
standardized ratios, underflowed transformed coordinates, finite Moyal log
density/CDF beyond direct exponential overflow, a finite half coordinate
after the full coordinate overflows, extreme Levy lower tails, and log
magnitudes genuinely outside binary64 range. Near-unit Levy CDFs use a
ratio check to reject a spuriously zero log probability.

The positive mixed sample has exact `[0.5,1.2,2.5]`, left `[0.3]`, right
`[3.5]`, and intervals `[1.1,1.4]`, `[4,4.4]`. Reflected Levy uses negated
exact measurements, exchanged/negated left/right limits, and reflected
interval endpoints. Location is zero and scale is 1.3:

| Family | Independent mpmath log likelihood | SciPy fitted scale, fixed loc=0 |
| --- | ---: | ---: |
| Moyal | -15.130277471440076 | 1.1208090746233839 |
| Levy | -15.99466380610435 | 0.6926757812499997 |
| Reflected Levy | -15.99466380610435 | 0.6926757812499997 |

Ordinary absolute log tolerances range from `3e-14` to `2e-13`; finite
extreme tail/density logs allow `3e-13`, and extreme ratio checks allow
`2e-13` (the finite half-coordinate check uses `3e-15`). SciPy fitted scales
allow `2e-4` for different optimizers. All-free checks establish bounded
convergence and finite objective consistency, without claiming global MLEs.

Native checked command (CMD, with matching local FPM dependencies):

```cmd
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Jbuild\normal_transform-native -Ibuild\normal_transform-native -Ibuild\gfortran_2806A042DAA8BBF7 src\special\log_gamma.f90 src\empirical_distribution.f90 src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored_normal_transform.f90 build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\normal_transform-native\test.exe
build\normal_transform-native\test.exe
```

The optimized run used `-O3` and separate `build\normal_transform-optimized`
module and executable paths. Changed modules were compiled explicitly before
linking existing FPM dependencies; directory hashes are local to this build.

## Independent references and integration

All 3,827 reference comparisons passed in checked and optimized builds.
The 326 new cases comprise 100 mixed samples per family and 26 explicit
interval stresses. Independent definitions use up to 1,000 mpmath digits;
the harness seed is 20261112. New stresses cover adjacent/subnormal widths,
overflowing displacements, extreme positive/negative Moyal tails, tiny
standardized Levy distances, and rounded coincident transformed endpoints.
The new families' maximum scaled error was `4.14e-14`; the overall maximum
remained `1.08e-12` in an existing gamma stress case. The metric is
`abs(actual-reference)/max(1,abs(reference))`, with threshold `2e-11`.
Native ratio checks supplement this absolute scale for tiny log probabilities.

Some huge erfc arguments overflow a range check inside mpmath 1.3.0.
The independent reference therefore uses the leading asymptotic term in
[DLMF 7.12.1](https://dlmf.nist.gov/7.12.E1) for arguments above `1e50`.
The first omitted term relative to the leading term is below `5e-101`, far
below the comparison tolerance. This approximation is confined to Python
reference generation; library calculations reuse existing normal kernels.

The harness configurations are
`-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace`
and `-O3 -fcheck=bounds`.

```cmd
python tools\generate_censored_reference.py
python validation\check_censored.py --module-dir build\gfortran_2806A042DAA8BBF7 --library build\gfortran_3C159EB802C22FB9\scifort\libscifort.a
fpm test
python tools\check_test_inputs.py
```

The standalone generator exited 0 and reproduced the likelihoods and SciPy
scale fits above. The guide example compiled and ran with the checked flags.
The full `fpm test` suite passed with exit code 0 in 1863.002 seconds
(31 minutes 3 seconds). No numerical changes followed that successful run.
The unattended-console-input check, source ASCII and line lengths, Python
syntax, inventory/API counts, changed local links, and `git diff --check` passed.

## Limits

The new APIs return quiet NaN for invalid data or parameters and retain
nonfinite mathematical endpoint/objective values. True log magnitudes
outside binary64 range return negative infinity; subnormal values retain
binary64 resolution limits. Normal-kernel accuracy limits carry through
the transformations. Fits use the existing bounded local optimizer, which
can encounter boundary or ill-posed optima; success does not establish a
global optimum. No scalar distribution source or dependency is changed.
C ABI and wrapper work remain deferred.
