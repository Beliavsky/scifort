# Generalized logistic, Gompertz, and inverse Weibull censoring, 2026-11-16

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
Reference libraries: SciPy 1.15.3 and mpmath 1.3.0.

## Native numerical and API checks

All 1,089 checks in `test/test_censored_hazard_family.f90` passed at both
`-O0` and `-O3`, each with
`-g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace`.
Checks cover mixed samples and NNLF, scalar density/CDF/SF agreement at
small, unit and larger shapes, CDF monotonicity/range, CDF/PPF round trips,
the original array fit interfaces, free and fixed bounded fits, SciPy
fixed-support shape estimates, invalid/empty inputs, support clipping,
impossible objectives, infinite interval endpoints, and equal endpoints
as exact data. Logistic shape-one and reciprocal-Weibull identities include
the exact-observation Jacobian and transformed interval endpoints.

Extreme checks retain the smallest positive interval width and shape,
underflowed powers, huge standardized coordinates compensated by tiny shape,
and tiny softplus/expm1 values compensated by huge shape. Ratio checks reject
zero log probabilities for nearly full intervals, nearly unit CDFs, and
small Gompertz survival hazards. A scalar inverse Weibull small-shape
overflow guard failed under traps before the two-line correction; ordinary
scalar agreement now passes with traps enabled.

The mixed sample has exact `[0.5,1.2,2.5]`, left `[0.3]`, right `[3.5]`, and
intervals `[1.1,1.4]`, `[4,4.4]`, at shape 1.7, location zero and scale 1.3:

| Family | Independent mpmath log likelihood | SciPy fitted c, fixed loc=0 and scale=1.3 |
| --- | ---: | ---: |
| Generalized logistic | -15.095845974501499 | 3.0464843750000044 |
| Gompertz | -71.5847386509558 | 0.12675781249999923 |
| Inverse Weibull | -24.434523984882752 | 0.7249999999999989 |

Ordinary absolute log tolerances range from `3e-14` to `3e-13`; extreme
ratio checks allow `2e-13`, finite tail logs allow `3e-13`, and fitted shapes
allow `2e-4` for the different optimizers. These shape comparisons use free
`c` only; the all-free checks establish finite objective consistency and
bounded convergence, without claiming global optimality.

Native checked command (CMD, with matching local FPM dependencies):

```cmd
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Jbuild\hazard_family-native -Ibuild\hazard_family-native -Ibuild\gfortran_2806A042DAA8BBF7 src\stats\invweibull.f90 src\special\log_gamma.f90 src\empirical_distribution.f90 src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored_hazard_family.f90 build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\hazard_family-native\test.exe
build\hazard_family-native\test.exe
```

The optimized run used `-O3` and separate `build\hazard_family-optimized`
module and executable paths. Changed modules were compiled explicitly before
linking existing FPM dependencies; directory hashes are local to this build.
The new guide example compiled and ran with the checked flags.

## Independent references and integration

All 3,501 reference comparisons passed in checked and optimized builds.
The 321 new cases comprise 100 mixed samples per family and 21 explicit
interval stresses. Independent mpmath definitions use up to 1,000 digits;
the harness seed is 20261112. New stresses cover adjacent and subnormal
widths, extreme tails, nearly full intervals, tiny shapes, overflowing
standardized coordinates or exponentials compensated by shape, and
underflowed standardized distances. The new families' maximum scaled error
was `2.36e-14`; the overall maximum remained `1.08e-12` in an existing gamma
stress case. The metric is `abs(actual-reference)/max(1,abs(reference))`,
with threshold `2e-11`. Native ratio checks supplement this absolute scale
for tiny log probabilities.

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
shape fits above. The full `fpm test` suite passed with exit code 0 in
1112.616 seconds (18 minutes 33 seconds). No numerical changes followed
that successful run. The unattended-console-input check, source ASCII and
line lengths, Python syntax, inventory/API counts, changed local links,
and `git diff --check` passed.

## Limits

The new APIs preserve nonfinite mathematical endpoint/objective values and
return quiet NaN for invalid data or parameters. True log magnitudes beyond
binary64 range return negative infinity; subnormal values retain binary64
resolution limits. Bounded fits use the existing local optimizer and may
reach boundary optima or ill-posed MLEs; success does not establish a global
optimum. The scalar inverse Weibull guard is corrected, while other scalar
extreme-range behavior remains unchanged. No dependency is added. C ABI
and wrappers remain deferred.
