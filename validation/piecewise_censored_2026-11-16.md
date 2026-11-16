# Bradford, triangular, and trapezoidal censoring validation, 2026-11-16

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
Reference libraries: SciPy 1.15.3 and mpmath 1.3.0.

## Native numerical and API checks

All 1,092 checks in `test/test_censored_piecewise.f90` passed at both
`-O0` and `-O3`, each with
`-g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace`.
Checks cover mixed references and NNLF, ordinary scalar density/CDF/SF
agreement, monotonicity/range, CDF/PPF round trips, original array fit
interfaces, free/fixed bounded fits, SciPy constrained shape fits,
optimizer controls, empty/impossible fit objectives, invalid shapes/data,
support clipping, infinite/equal interval endpoints, endpoint modes,
and the finite rounded upper-edge contract.

Identities include triangular as equal-edge trapezoidal, uniform as
trapezoidal `(0,1)`, endpoint-mode triangles as beta `(1,2)` and `(2,1)`,
and the uniform limit of the smallest positive Bradford shape. Extreme
checks cover the smallest interval width, underflowed compensated upper
distances, and an interval crossing the smallest positive trapezoid
breakpoint. Its CDF beyond that breakpoint is also checked independently
using `F(z)=2*z-c+O(z*z)`, with negligible quadratic correction there.
Nearly full triangular/Bradford intervals use ratio checks that reject
a spuriously zero log probability and retain both omitted tails.

The mixed sample has exact `[0.1,0.3,0.6]`, left `[0.05]`, right `[0.8]`,
and intervals `[0.2,0.25]`, `[0.7,0.75]`, at location zero and width one:

| Family and reference shapes | Independent mpmath log likelihood | SciPy constrained fit |
| --- | ---: | --- |
| Bradford, c=1.7 | -10.201408354319534 | c=2.8134765625000036, fixed support |
| Triangular, c=0.3 | -13.055640654451915 | c=0.000048828124999111865, fixed support |
| Trapezoidal, c=0.2 and d=0.7 | -12.541072613831128 | c=3.885780586188048e-16, fixed d=0.7 and support |

The triangular and fixed-d trapezoidal optima lie at `c=0`; the SciPy
estimates above reflect optimizer tolerances. SciPy's default both-shape
trapezoidal fit on this sample instead returned `(1,1)`, with log likelihood
`-15.36946646093937`, worse than the initial `(0.2,0.7)` reference and
the `(0,0)` triangular case (`-10.99670999256986`). That default result is
not treated as an MLE reference. The constrained comparison supplies
initial shapes `(0.2,0.7)` and fixes the second shape at 0.7. Native all-free
checks establish finite bounded convergence/objective consistency, without
claiming global optimality.

Ordinary absolute log tolerances range from `4e-14` to `2e-13`;
extreme finite logs allow `4e-13`, near-unit ratios allow `3e-13`, and
SciPy fitted first shapes allow `2e-4` for different optimizers.

Native checked command (CMD, with matching local FPM dependencies):

```cmd
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Jbuild\piecewise-native -Ibuild\piecewise-native -Ibuild\gfortran_2806A042DAA8BBF7 src\special\log_gamma.f90 src\empirical_distribution.f90 src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored_piecewise.f90 build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\piecewise-native\test.exe
build\piecewise-native\test.exe
```

The optimized run used `-O3` and separate `build\piecewise-optimized`
module and executable paths. Changed modules were compiled explicitly before
linking existing FPM dependencies; directory hashes are local to this build.

## Independent references and integration

All 4,277 reference comparisons passed in checked and optimized builds.
The 450 new cases comprise 100 mixed samples per family and 150 explicit
interval stresses across limiting and ordinary shapes. Independent
mpmath density and integrated piecewise-polynomial definitions use up to
1,000 digits; the harness seed is 20261112. Stresses include corner-crossing,
adjacent/subnormal widths, tiny/huge Bradford shapes, subnormal trapezoid
breakpoints, compensated nearly full intervals, and support clipping.
The new families' maximum scaled error was `3.62e-15`; the overall maximum
remained `1.08e-12` in an existing gamma stress case. The metric is
`abs(actual-reference)/max(1,abs(reference))`, with threshold `2e-11`.
Native ratio checks supplement that scale for tiny log probabilities.

The harness configurations are
`-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace`
and `-O3 -fcheck=bounds`.

```cmd
python tools\generate_censored_reference.py
python validation\check_censored.py --module-dir build\gfortran_2806A042DAA8BBF7 --library build\gfortran_3C159EB802C22FB9\scifort\libscifort.a
fpm test
python tools\check_test_inputs.py
```

The standalone generator exited 0 and reproduced the likelihoods and
constrained SciPy fits above. The guide example compiled and ran with
checked flags. The full `fpm test` suite passed with exit code 0 in
1522.026 seconds (25 minutes 22 seconds). No numerical changes followed
that successful run. The unattended-console-input check, source ASCII and
line lengths, Python syntax, inventory/API counts, changed local links,
and `git diff --check` passed.

## Limits

All three new APIs require a finite rounded upper support edge, following
the existing beta endpoint helpers. Invalid data/parameters return quiet
NaN, while zero-probability observations retain negative infinity. Binary64
subnormal resolution still applies. Bradford uses its uniform limit below
`epsilon*exp(-2)`, with relative corrections below roundoff. Fits use the
existing bounded local optimizer; flat or boundary optima and ill-posed
MLEs remain possible, and success does not establish a global optimum.
Other scalar kernels' extreme-range behavior is unchanged. No dependency
or scalar distribution source is changed. C ABI and wrappers remain deferred.
