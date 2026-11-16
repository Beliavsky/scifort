# Fisk and Burr III/XII censoring validation, 2026-11-16

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
Reference libraries: SciPy 1.15.3 and mpmath 1.3.0.

## Native numerical and API checks

All 1,190 checks in `test/test_censored_burr_family.f90` passed at both `-O0`
and `-O3`, each with `-g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace`.
Checks cover three mixed references, NNLF, ordinary scalar agreement for small,
unit, and larger first shapes, CDF monotonicity and range, CDF/PPF round trips,
original array likelihoods and fits, endpoint densities, free/fixed fit
parameters, invalid/empty inputs, optimizer controls, and nonfinite objectives.
Infinite interval endpoints and equal endpoints as exact data are exercised.

Identities include Fisk as `d=1` in either Burr family; Fisk as a logarithmic
logistic transform with exact-density Jacobian; Burr XII at `c=1` as Lomax;
and reciprocal Burr III/Burr XII mixed samples with the reciprocal Jacobian.
Analytic second-shape MLEs use exact `[0.5,1,2]` at `c=1`, location zero,
scale one: `d=3/sum(log(1+x))` for XII, and `d=3/sum(log(1+1/x))` for III.
The reciprocal sample is identical, so the two analytic values agree.

The mixed sample has exact `[0.5,1.2,2.5]`, left `[0.3]`, right `[3.5]`,
and intervals `[1.1,1.4]`, `[4.0,4.4]`. At `c=1.7`, location zero and scale
1.3, with `d=0.8` for Burr families, the independent log likelihoods are:

| Family | mpmath log likelihood | SciPy fit with location/scale fixed |
| --- | ---: | --- |
| Fisk | -14.64543894349949 | c = 1.3590820312500007 |
| Burr III | -14.710130480007274 | d = 0.9417968749999999, c fixed at 1.7 |
| Burr XII | -14.609935975674242 | d = 0.8731445312499999, c fixed at 1.7 |

Scalar, identity and mixed-reference tolerances range from `2e-14` to `7e-14`
in absolute log error; extreme log tails allow up to `3e-12`.
An adjacent interval check uses density times width with `4e-14` tolerance.
Nearly unit interval probabilities use a ratio with `5e-14` tolerance to
reject an erroneous zero. Different optimizers' fitted parameters allow
`3e-4`; analytic second-shape MLEs allow `4e-6`. The harder all-free Burr III
case uses a smaller parameter box, 10,000 maximum sweeps and `1e-6`
parameter tolerance; the other native fits use the normal defaults.

The checked comparisons exposed existing scalar Fisk/Burr overflow guards
that evaluated `huge()/shape` for shapes below one. Burr also relied on
short-circuit evaluation, which Fortran does not guarantee. Explicit small-
multiplier branches now avoid those divisions and preserve the public formulas.
The small-shape scalar comparisons run with traps enabled in both builds.

Native checked command (CMD, with matching local FPM dependencies):

```cmd
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Jbuild\burr_family-native -Ibuild\burr_family-native -Ibuild\gfortran_2806A042DAA8BBF7 src\stats\burr.f90 src\stats\fisk.f90 src\special\log_gamma.f90 src\empirical_distribution.f90 src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored_burr_family.f90 build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\burr_family-native\test.exe
build\burr_family-native\test.exe
```

The optimized run used `-O3` and separate `build\burr_family-optimized` module
and executable paths. Changed modules, including the scalar guards, were
compiled explicitly before linking existing FPM dependencies. Directory hashes
are local to this build. The combined guide example compiled and ran with
checked flags and explicit compilation of the changed scalar modules.

## Independent references and integration

All 2,842 reference comparisons passed in checked and optimized builds.
The 339 new cases comprise 100 mixed samples and 11 interval stresses per
family, plus three compensated-shape exact/CDF/SF cases for each Burr family.
Independent mpmath definitions and log probability differences use up to
800 digits, with seed 20261112. Cases cover smallest positive widths and first
shape, underflowed/overflowing standardized ratios, overflowing displacements,
wide intervals, adjacent points at the join with huge shape, and same-side
or crossing softplus intervals. Opposing shapes `c=1e308`, `d=1e-308` exercise
representable hazards when the intermediate power logarithm overflows.

The batch maximum scaled error was `5.87e-14`; the overall maximum remained
`1.08e-12` in an existing gamma stress case. The metric is
`abs(actual-reference)/max(1,abs(reference))`, with threshold `2e-11`.
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
fixed-support SciPy fits above. The full `fpm test` suite passed with exit
code 0 in 1106.945 seconds (18 minutes 27 seconds). No numerical
changes followed that run. The unattended-console-input check, source ASCII
and line lengths, Python syntax, inventory counts, changed local links, and
`git diff --check` passed.

## Limits

Fits use the existing bounded local search without uniqueness or global-MLE
guarantees. Singular endpoint densities and weak shape/location/scale
identification can make free-location fits ill-posed or slow. Native Fortran
coverage takes priority; C ABI and wrappers remain deferred. Scalar routines
receive the small-shape guard fixes but do not otherwise gain the censored
kernels' extreme-range guarantees. No dependency is added.
