# Truncated exponential, Johnson SB, and kappa3 validation, 2026-11-17

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
Reference libraries: SciPy 1.15.3 and mpmath 1.3.0.

## Native API and numerical checks

All 1,843 checks in `test/test_censored_bounded_shape.f90` passed at both
`-O0` and `-O3`, each with
`-g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace`.
They cover ordinary scalar density/CDF/SF agreement, monotonicity,
CDF/PPF round trips, NNLF agreement, original array-fit overloads,
fixed/partly free/free fits, independent SciPy shape benchmarks,
invalid parameters and samples, empty samples, support endpoints,
equal/infinite interval endpoints, and adjacent floating-point intervals.
The Johnson SB checks include a fit with all four parameters free.

Identities and stresses include kappa3/Lomax at unit shape, the tiny-shape
kappa3 CDF limit `exp(-1)`, overflowing endpoint log density, large shapes
amplifying a scale-point displacement, an unrepresentable `loc+scale`,
large Johnson shapes amplifying a midpoint displacement, and a subnormal
interval with a huge support scale. Truncated exponential retains the
nearly uniform log-density correction `2.5e-301` in a ratio-based test.
Both bounded families have explicit tests for odd subnormal widths;
these catch the half-step error from rounding a midpoint.

Ordinary scalar log tolerances are `2e-12`, CDF/PPF round trips `3e-12`,
adjacent/extreme log checks `3e-13`, and fixed-array-fit comparisons
`2e-13`. The tiny log-density test compares the ratio to its independent
reference at `3e-13`. Different fitting algorithms allow `2e-6` in shape.

Checked native commands (CMD, using matching local FPM dependencies):

```cmd
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Jbuild\bounded_shape-native -Ibuild\bounded_shape-native -Ibuild\gfortran_2806A042DAA8BBF7 src\special\log_gamma.f90 src\empirical_distribution.f90 src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored_bounded_shape.f90 build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\bounded_shape-native\test.exe
build\bounded_shape-native\test.exe
```

The optimized run substitutes `-O3` and separate module/executable paths
under `build\bounded_shape-optimized`. Changed modules are compiled
explicitly before linking the existing archive. These standalone commands
require the build directories to exist.

## Independent references

The new independent reference checker passed all 387 comparisons in each
checked and optimized build: 300 seeded mixed samples and 87 finite point,
tail, and interval stresses. Maximum scaled log-likelihood error was
`4.13e-14`, in a Johnson SB mixed sample. Ordinary mixed samples use 100
decimal digits and scales from `1e-100` to `1e100`; extreme cases use 800
digits. Stresses include tiny/huge shapes, opposing shape/scale magnitudes,
subnormal and adjacent gaps, nearly full intervals, amplified midpoint or
scale-point displacements, odd subnormal widths, and an overflowing
location-plus-scale point that remains a valid kappa3 parameterization.

The oracle evaluates the mathematical density integrals directly, using
arbitrary-precision exponential increments for truncated exponential,
complementary normal erfc tails for Johnson SB, and elementary reversed
hazards for kappa3. It follows the documented rounded-width law for
truncated exponential. It requires finite actual and reference log values,
with `abs(actual-reference)/max(1,abs(reference)) <= 2e-11`; native tests
separately check nonfinite endpoints and tiny relative log-density errors.

```cmd
python validation\check_censored_bounded_shape.py --module-dir build\gfortran_2806A042DAA8BBF7 --library build\gfortran_3C159EB802C22FB9\scifort\libscifort.a
```

The checker prints its exact compiler commands. Checked flags are
`-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow`; optimized flags
are `-O3 -fcheck=bounds`. The default seed is `20261117`.

Independent SciPy bounded optimization uses exact `[0.2,0.4,0.6]`,
left `[0.1]`, right `[0.7]`, and intervals `[0.25,0.35]`, `[0.5,0.55]`.
Location is fixed at zero and scale at 1.3; Johnson b is fixed at 1.2:

| Family | Shape | Negative log likelihood |
| --- | ---: | ---: |
| Truncated exponential | b=0.6650398044582356 | 8.384192417965537 |
| Johnson SB | a=1.120500197116846 | 9.120521933545074 |
| Kappa3 | a=3 (upper bound) | 11.921900024027959 |

Search bounds are `[0.6,3]`, `[-3,3]`, and `[0.1,3]`, respectively.
SciPy's numerical bounded search returned kappa3 shape
`2.9999999131353752`, just below the boundary.

## Integration and documentation checks

The full `fpm test` run passed with exit code 0 in 1,605.599 seconds
(26 minutes 46 seconds), including all 1,843 new native checks. This
integration run covers shared likelihood, fitting-dispatch, statistics
exports, and their dependents. No executable changes followed that
successful run; final validation-record edits were documentation only.

The new guide example compiled and ran with checked flags, reproducing the
Johnson shape benchmark. Source ASCII/line lengths, Python syntax, changed
Markdown links, inventory/API-list consistency (59 supported plus 50
remaining), `git diff --check`, and `python tools\check_test_inputs.py`
passed. No new dependency, generated artifact, or third-party implementation
was added to the tracked changes.

## Limits

Binary64 is the public precision. Truncated exponential uses a normalized
rounded physical support width; underflowed, overflowed, and collapsed
supports are rejected. Johnson SB requires a finite rounded upper edge
strictly above the lower edge, with compensated addition residuals retained.
Subnormal width rounding can produce differences from original scalar APIs,
which remain unchanged. Large opposing Johnson shape terms make its normal
coordinate ill-conditioned; accuracy is not uniform through cancellation.

Kappa3 does not require a finite location-plus-scale point. An endpoint
log density can overflow for a valid tiny shape, and free-location fitting
can have an unbounded likelihood in singular parameter limits. Fits are
bounded local searches, without a global-optimum or identifiability
guarantee. See [the API guide](../docs/censored_fitting.md) for complete
contracts and [provenance](../CODE_PROVENANCE.md) for mathematical sources.
No C ABI or wrapper expansion is included.
