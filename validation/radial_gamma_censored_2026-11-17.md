# Chi-squared, chi, and Nakagami censored validation, 2026-11-17

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.

## Native checks

All 1,924 checks in [the native test](../test/test_censored_radial_gamma.f90)
passed at both `-O0` and `-O3`, each with
`-g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace`.
They cover scalar density/CDF/SF comparisons, CDF monotonicity, CDF/PPF and
SF/ISF round trips, adjacent intervals, affine scales 1e-200 through 1e200,
support and singular endpoints, invalid data/parameters, empty data,
infinite interval endpoints, exact equal-endpoint conventions, minimum
representable shapes, fixed and shape-only fits, all-free finite candidate
improvement, original array-fit overloads, and NNLF agreement.

Independent identities compare chi-squared with gamma/exponential, chi
with half-normal/Rayleigh/Maxwell, and Nakagami with half-normal/Rayleigh.
Regressions retain subnormal half-normal interval probabilities and
validate a Nakagami far tail at x=1e200, nu=1e-300. The scalar-tail guard
now runs without overflowing for shapes below one.

Ordinary scalar logs allow 2e-12 scaled error, inverse round trips 3e-12,
and adjacent/extreme logs and identities 3e-13. Exact classifications use
no tolerance. These measures have an absolute floor of one; tiny relative
tail comparisons divide by their expected magnitude before comparison.

Checked command (CMD, with an existing output directory and FPM cache):

```cmd
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Jbuild\radial_gamma-native -Ibuild\radial_gamma-native -Ibuild\gfortran_2806A042DAA8BBF7 src\special\log_gamma.f90 src\stats\nakagami.f90 src\empirical_distribution.f90 src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored_radial_gamma.f90 build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\radial_gamma-native\test.exe
build\radial_gamma-native\test.exe
```

The optimized run substitutes `-O3` and separate module/executable paths
under `build\radial_gamma-optimized`, retaining all runtime/trap flags.
Changed modules are compiled explicitly before the archive is linked.

## Independent references

The [independent checker](check_censored_radial_gamma.py) passed all 369
cases in each checked and optimized build: 300 seeded mixed samples and
69 stresses. Maximum scaled log-likelihood error was `6.389e-13`, in a
Nakagami mixed sample. Random samples use 100 decimal digits, scales from
1e-200 through 1e200, and interval widths down to adjacent floats. Stresses
use 1100 digits, with scales 1e-300/1e300, minimum representable shapes,
subnormal intervals, overflowing physical differences, nearly full support,
and underflowed powered distances. Six maximum-shape center/adjacent
references use 90 digits.

The oracle independently evaluates gamma integrals, retaining small-shape
upper probabilities with the alternating expansion from DLMF 8.7. It uses
the documented rounded half-shape. Finite log comparisons require
`abs(actual-reference)/max(1,abs(reference)) <= 2e-11`; native tests separately
check nonfinite classifications and relative far-tail behavior.

```cmd
python validation\check_censored_radial_gamma.py --module-dir build\gfortran_2806A042DAA8BBF7 --library build\gfortran_3C159EB802C22FB9\scifort\libscifort.a
```

The checker prints exact compiler commands. Checked flags are
`-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow`; optimized flags are
`-O3 -fcheck=bounds`. Seed: `20261117`. References use SciPy 1.15.3 and
mpmath 1.3.0. A separate SciPy bounded optimization uses exact
`[0.2,0.4,0.6]`, left `[0.1]`, right `[0.7]`, intervals `[0.25,0.35]` and
`[0.5,0.55]`, fixed loc=0/scale=1.3, and shape bounds `[0.3,5]`:

| Family | Shape | NNLF |
| --- | ---: | ---: |
| Chi-squared | 0.9099424378451767 | 11.137830726064342 |
| Chi | 0.5690179321507661 | 10.012294312027775 |
| Nakagami | 0.3000000078454339 | 10.649162735292245 |

Nakagami's optimum is at the lower search bound; SciFort reaches 0.3 exactly.
Native fit tolerances are 3e-5 in shape and 2e-8 in scaled objective.

## Integration

The full `fpm test` run passed with exit code 0 in 2480.727 seconds
(41 minutes 21 seconds), using the compiler identified above. It includes
all 1,924 new native checks and all dependent fitting/statistics integration
tests. This full run satisfies the supplied repository instructions and
covers the small scalar Nakagami change as well as the new private
censoring callbacks. It used the final numerical source/test state;
subsequent edits only finalized documentation.

Existing array-temporary runtime warnings appeared in unchanged interpolation
and signal tests. The console-input scan (`python tools/check_test_inputs.py`),
`git diff --check`, and audits of source ASCII/line lengths, Python syntax,
changed documentation links, and the 71 implemented/38 remaining inventory
passed. No C ABI or wrapper surface was expanded.

## Limits

Chi and chi-squared require positive representable fl(df/2), which defines
the numerical law for subnormal degrees of freedom. The current shape
limits are df<=2e6 and nu<=1e6; larger shapes return NaN. All three reuse the
existing gamma kernels, including their accuracy and representability
limits; no new large-shape asymptotic is introduced. Fitting uses local
bounded searches. An all-free search can reach the iteration limit with
a finite improving candidate. Endpoint/singular parameter limits can make
maximization unbounded. Success does not prove global optimality or
identifiability. See the [API guide](../docs/censored_fitting.md) and
[provenance](../CODE_PROVENANCE.md). C ABI/wrapper expansion is deferred.
