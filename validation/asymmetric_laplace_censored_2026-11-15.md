# Asymmetric Laplace censoring validation, 2026-11-15

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
Reference libraries: SciPy 1.15.3 and mpmath 1.3.0.

## Native checks

All 282 checks in `test/test_censored_laplace_asymmetric.f90` passed at both
`-O0` and `-O3`, each with
`-g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace`.
Checks cover mixed observations, defaults/keywords, negative likelihood,
symmetric Laplace and reciprocal-shape reflection identities, scalar/array
agreement, free/fixed shape/location/scale, independent analytic scale MLE,
SciPy mixed fits, objective agreement, optional controls, empty/invalid input,
infinite interval endpoints, and equal endpoints as exact cusp observations.

Stress checks cover extreme log tails, nearly unit CDF and interval masses,
including a same-side interval whose mass is nearly one, adjacent far-tail
intervals, subnormal crossing widths, underflowed standardized intervals,
overflowing raw displacements, the smallest/largest positive binary64 shapes,
underflowed CDF/SF probabilities with finite logs, and compensating shape/scale.
Near-one logs use ratio checks so an erroneous zero cannot pass an absolute
log-likelihood tolerance. Extreme cases run with floating-point traps enabled.

The independent mixed sample has exact `[0.5,1.2,2.5]`, left `[0.3]`, right
`[3.5]`, and intervals `[1.1,1.4]`, `[4.0,4.4]`. At shape 1.7, location 0.4,
scale 1.3, its 100-digit mpmath log-likelihood is
`-23.31171335066417074412690936788162160893`. At fixed location 0.4, SciPy
returns shape `0.4577391395316225`, scale `1.0270210720704078`. At fixed
shape 1.7 and location 0.4, its fitted scale is `3.6634935386148726`.
For exact `[-2,-1,1,3]` at shape 2 and location zero, scale MLE is 2.375:
the weighted absolute deviations sum to `2*(1+3)+(2+1)/2=9.5`, divided by 4.

Ordinary scalar/array comparisons allow `2e-14` absolute error, the mixed
reference `9e-14`, and symmetry identities `5e-14`. Extreme log tails allow
up to `3e-12`, log-width extremes up to `5e-13`, and the adjacent-tail
approximation `4e-12`. Near-one ratios allow `3e-13`. SciPy fitted parameters
allow `4e-4` for different optimizers; the independent analytic MLE allows
`4e-6`. Exact identities are used for the cusp normalization with matching
smallest shape/scale and extreme cusp log tails.

The new guide example compiled and ran with checked flags. All products were
placed under ignored `build/` directories. Native checked command (CMD,
using the matching local FPM module directory and dependency archive):

```cmd
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Jbuild\asymmetric-native -Ibuild\asymmetric-native -Ibuild\gfortran_2806A042DAA8BBF7 src\special\log_gamma.f90 src\empirical_distribution.f90 src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored_laplace_asymmetric.f90 build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\asymmetric-native\test.exe
build\asymmetric-native\test.exe
```

The optimized run used `-O3` and separate `build\asymmetric-optimized` module
and executable paths. Changed modules were compiled explicitly before linking
existing FPM dependencies. Directory hashes are local to this build.

## Reference and integration checks

All 2,166 reference comparisons passed in both checked and optimized builds.
The 117 new comparisons comprise 100 mixed asymmetric-Laplace samples, 13
interval stress cases, and four extreme-shape/scale exact or censored values.
Independent mpmath calculations use direct density/CDF/SF definitions and
logarithmic probability differences at up to 800 digits, with seed 20261112.

The maximum scaled error for the new family was `5.77e-14`; the overall maximum
remained `1.08e-12` in an existing gamma stress case. The harness uses
`abs(actual-reference)/max(1,abs(reference)) <= 2e-11`. Configurations are
`-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace` and
`-O3 -fcheck=bounds`.

```cmd
python tools\generate_censored_reference.py
python validation\check_censored.py --module-dir build\gfortran_2806A042DAA8BBF7 --library build\gfortran_3C159EB802C22FB9\scifort\libscifort.a
fpm test
python tools\check_test_inputs.py
```

The standalone reference generator exited 0 and reproduced the likelihood
and both SciPy fits above. Existing SciPy Erlang noninteger-shape warnings
are emitted during that separate comparison. The unattended-console-input,
source/Python syntax, inventory-count, and local-link checks passed.
The full `fpm test` suite passed with exit code 0 in 1099.481 seconds
(18 minutes 19 seconds). No numerical changes followed that run.

## Limits

The new likelihood reuses the existing softplus, log-hazard mass, and guarded
negative-hazard helpers, retaining their limitations. Fits are bounded local
searches without uniqueness or global-optimality guarantees. Scalar/array APIs,
C ABI, and dependencies are unchanged. The new censored kernels handle extreme
shapes and tiny probabilities separately; this change does not extend those
guarantees to existing scalar routines.
