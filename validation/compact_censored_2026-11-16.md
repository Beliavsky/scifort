# Cosine, anglit, and semicircular censoring validation, 2026-11-16

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
Reference libraries: SciPy 1.15.3 and mpmath 1.3.0.

## Native checks

All 1,516 checks in `test/test_censored_compact.f90` passed at both `-O0`
and `-O3`, each with
`-g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace`.
They cover scalar density/CDF/SF agreement, monotonicity, PPF round trips,
reflection identities, NNLF agreement, original array-fit overloads,
fixed and free location/scale fitting, independent SciPy scale benchmarks,
invalid parameters and samples, empty samples, support endpoints, impossible
fixed objectives, equal/infinite interval endpoints, and whole-support or
outside-support intervals. Tail and interval tests include subnormal gaps,
extreme scales, the smallest positive scale, and adjacent observations near
support endpoints. Narrow integration uses gap coordinates so quadrature
nodes cannot collapse to rounded physical endpoints.

Ordinary scalar log comparisons use scaled absolute tolerance `2e-12`;
round trips use `1e-12`, reflection `2e-14`, and extreme or adjacent-interval
references `3e-13`. Different fitting algorithms allow `2e-6` for scale.

Checked native commands (CMD, with matching local FPM dependencies):

```cmd
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Jbuild\compact-native -Ibuild\compact-native -Ibuild\gfortran_2806A042DAA8BBF7 src\special\log_gamma.f90 src\empirical_distribution.f90 src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored_compact.f90 build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\compact-native\test.exe
build\compact-native\test.exe
```

The optimized run replaces `-O0` with `-O3` and uses separate module and
executable paths under `build\compact-optimized`. These standalone compiler
commands require the build directories to exist. Changed modules are
compiled explicitly before linking the existing FPM archive.

## Independent numerical references

The new independent reference checker passed 396 comparisons in each
checked and optimized build: 300 seeded mixed-censoring samples and 96
endpoint, adjacent-interval, nearly-full-interval, subnormal, extreme-scale,
and center-crossing cases. Maximum scaled log-likelihood error was
`1.63e-14`, in a semicircular mixed sample. Ordinary samples use 100 decimal
digits and scales between `1e-100` and `1e100`; stresses include `1e-300`
and `1e300`. Subnormal endpoint gaps use 800 digits.

The independent formulas integrate the public density specifications using
mpmath hypergeometric functions for cosine, the sine-square integral for
anglit, and regularized incomplete beta for semicircular. Extremely narrow
central intervals integrate a locally constant positive density, with
relative variation below `1e-30`. The checker follows the documented rounded
physical half-width contract. It requires finite reference/actual results
and `abs(actual-reference)/max(1,abs(reference)) <= 2e-11`; native checks
separately cover nonfinite endpoint behavior.

```cmd
python validation\check_censored_compact.py --module-dir build\gfortran_2806A042DAA8BBF7 --library build\gfortran_3C159EB802C22FB9\scifort\libscifort.a
```

This command prints the exact compiler invocations. The checked profile
uses `-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow`; the optimized
profile uses `-O3 -fcheck=bounds`. The default seed is `20261116`.

For exact `[-0.2,0.1,0.35]`, left `[-0.4]`, right `[0.45]`, and intervals
`[-0.3,-0.1]`, `[0.2,0.25]`, independent SciPy bounded optimization with
location fixed at zero and scale in `[0.7,2]` gives:

| Family | Scale | Negative log likelihood |
| --- | ---: | ---: |
| Cosine | 0.7 (lower bound) | 11.008968106170082 |
| Anglit | 1.0646672682233491 | 9.326188349907682 |
| Semicircular | 0.7078037543218011 | 9.339707758713789 |

SciPy's numerical bounded search returned cosine scale
`0.7000000136831084`, just above the boundary.

## Integration and documentation

The existing reference regression check passed all 6,318 comparisons in
each checked/optimized build, with maximum scaled error `1.08e-12` in an
existing gamma stress case. Together with the new checker, 6,714 references
passed per build. The regression command was:

```cmd
python validation\check_censored.py --module-dir build\gfortran_2806A042DAA8BBF7 --library build\gfortran_3C159EB802C22FB9\scifort\libscifort.a
```

The full `fpm test` run passed with exit code 0 in 1,365.013 seconds
(22 minutes 45 seconds), including all 1,516 new native checks. The full
run covers shared censored-likelihood, fitting dispatch, statistics exports,
and their integration paths. No executable changes followed that successful
run; final validation-record edits were documentation only.

The new guide example compiled and ran with checked flags. Source ASCII
and line lengths, Python syntax, changed Markdown links, coverage/API-list
consistency (56 supported plus 53 remaining), `git diff --check`, and
`python tools\check_test_inputs.py` passed.

## Numerical limits

The native precision is binary64. These likelihoods use the rounded physical
half-width with finite rounded support edges strictly straddling the center;
unrepresentable or collapsed supports return NaN. Compensated distances
retain the addition residual, so a rounded edge need not be the mathematical
endpoint. Subnormal half-width rounding can cause larger differences from
the original scalar API, whose implementation is unchanged. See
[the API guide](../docs/censored_fitting.md) for the complete support contract
and [provenance](../CODE_PROVENANCE.md) for the mathematical specifications.

Fits are bounded local searches and are not guaranteed global optima.
No C ABI or wrapper expansion is included.
