# Bounded-family censoring validation, 2026-11-18

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
Independent references: SciPy 1.15.3 and mpmath 1.3.0.

## Native regression

The [native test](../test/test_censored_bounded_final.f90) passes 2,128
checks in both checked `-O0` and optimized `-O3` builds. Both enable
runtime/bounds checks and invalid/zero/overflow traps. Coverage includes
scalar density/tail agreement, monotonic CDFs, CDF/PPF and SF/ISF round
trips, ordinary and adjacent intervals, affine scales 1e-200/1/1e200,
invalid/empty data and parameters, NaN/infinity behavior, equal endpoints,
infinite one-sided intervals, fixed and one-parameter fits, improving
all-free candidates, and compatibility with existing array-fit overloads.

Identity checks cover the beta case of Gauss hypergeometric, Irwin-Hall's
n=1 uniform case, and the uniform line von Mises limit. Endpoint tests
include zero, finite, and singular densities; the last representable
intervals below the Gauss hypergeometric and n=32 supports; and Irwin-Hall
probabilities with log magnitudes near 24,000. The new family-module kernel
integral is checked for normalization, equal endpoints, and invalid bounds.
An integer-constrained Irwin-Hall fit verifies its discrete n search.

CMD commands, with existing output directories and FPM dependencies:

```cmd
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Jbuild\bounded_final-native -Ibuild\bounded_final-native -Ibuild\gfortran_2806A042DAA8BBF7 src\stats\gausshyper.f90 src\special\log_gamma.f90 src\empirical_distribution.f90 src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored_bounded_final.f90 build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\bounded_final-native\test.exe
build\bounded_final-native\test.exe
```

The optimized command substitutes `-O3` and paths under
`build\bounded_final-optimized`, retaining all checks and traps. Changed
modules are explicitly compiled before linking the existing archive.
The scaled error measure is abs(actual-reference)/max(1,abs(reference)).
Ordinary scalar and adjacent/affine comparisons allow 3e-12, ordinary
intervals 3e-11, and quantile round trips 5e-11. Tiny Irwin-Hall and kernel
normalization checks allow 3e-13; identities allow 2e-12. Classifications
and exact support values use no tolerance. These log-domain tolerances
are not claims of a uniform relative probability error after exponentiation.

## Independent references

The [checker](check_censored_bounded_final.py) evaluates 167 high-precision
comparisons: 60 seeded mixed samples and 107 stresses, in checked and
optimized builds. All references use 180 decimal digits and seed 20261118.
All 167 comparisons pass in both builds, with maximum scaled error
`1.97074e-14`. Finite comparisons require scaled error <=3e-11.

```cmd
python validation\check_censored_bounded_final.py --module-dir build\gfortran_2806A042DAA8BBF7 --library build\gfortran_3C159EB802C22FB9\scifort\libscifort.a
```

Gauss hypergeometric references use the independent beta*2F1 normalization
and mpmath density quadrature. Endpoint-power substitutions remove
singularities before integration; interior intervals factor out a reference
density and physical width so tiny masses do not defeat absolute quadrature
stopping rules. Irwin-Hall references use inclusion-exclusion polynomials,
independent of the Fortran positive B-spline recurrence. Line von Mises
uses Bessel I0 normalization and direct density quadrature, with the
documented binary64 support endpoints.

Stresses cover the Gauss hypergeometric shape-bound combinations, tiny
left tails at 1e-100, adjacent interior and last-upper-endpoint intervals;
Irwin-Hall n=1/2/3/16/32 and points down to 1e-200, adjacent intervals,
subnormal positive inputs and intervals, the last n=32 upper support
interval, and 15 interior polynomial-piece comparisons including n=16/32; line von Mises concentration
0/1e-10/1/100, near-boundary tails and adjacent intervals; and affine
scales 1e-300/1e300. Narrow physical widths are retained in reference
coordinates using high precision before subtraction and division.

Checked driver flags are `-O0 -g -fcheck=all
-ffpe-trap=invalid,zero,overflow`; optimized flags are `-O3 -fcheck=bounds`.
The optional `--reuse` flag uses ignored reference files under `build/`.

Independent SciPy bounded fits use exact [0.2,0.4,0.6], left [0.1],
right [0.8], and intervals [0.3,0.35] and [0.5,0.55]. Only scale varies
in [0.3,3], with loc=0:

| Family | Fixed shapes | Fitted scale | NNLF |
| --- | --- | ---: | ---: |
| gausshyper | a=1.7, b=2.7, c=0.5, z=1 | 1.2233448777774485 | 9.704537797653995 |
| irwinhall | n=3 | 0.33700479821178475 | 12.807592741295426 |
| vonmises_line | kappa=2 | 0.7485579168664227 | 12.366078083348853 |

Native comparisons allow 3e-5 scaled parameter error and 2e-8 scaled
NNLF error. All-free checks require finite improving candidates without
claiming global optima or identifiable parameters.

## Full suite and limitations

The shared censoring/fit/statistics integration changed, so the full suite
is required in addition to focused checked and optimized validation:

```cmd
fpm test
```

The full suite passed with exit code 0 in 2019.890 seconds (33 minutes
40 seconds). This includes all 2,128 new checks, existing censored and
scalar-distribution regressions, fitting and C API tests, and the remaining
repository tests. Existing array-temporary runtime warnings appeared in
unrelated interpolation/spatial tests and did not cause failures.
No numerical library or native test code changed after this successful
run began. Later edits extended the independent reference cases and
completed documentation and validation metadata.

The [API guide](../docs/censored_fitting.md) documents the validated shape
and affine limits. Wider ranges and improved endpoint conditioning remain
numerical compatibility gaps even though the bounded-family name gap is
closed. Standardization that rounds or underflows onto a support endpoint
uses endpoint behavior. Line von Mises uses a fixed real interval, with no
wrapping or arc semantics. Additional precision for ill-conditioned affine endpoints and
unrepresentable quantiles is not provided. Native coverage is the current
scope; C ABI expansion remains deferred.
