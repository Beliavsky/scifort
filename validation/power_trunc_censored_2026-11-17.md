# Generalized gamma and truncated-law censoring validation, 2026-11-17

Compiler: GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows.
Independent references: SciPy 1.15.3 and mpmath 1.3.0.

## Native regression

All 2,125 checks pass at both `-O0` and `-O3`, with bounds/runtime checks
and invalid/zero/overflow traps enabled. This includes explicit rejection
of subnormal or poorly resolved affine truncation bounds.

The [native test](../test/test_censored_power_trunc.f90) covers generalized
gamma, truncated normal, and truncated Weibull: scalar PDF/CDF/SF agreement,
monotonicity, CDF/PPF and SF/ISF round trips, adjacent/ordinary intervals,
scales 1e-200/1/1e200, both gamma-power signs, invalid/empty inputs,
support endpoints, infinite one-sided intervals, equal-endpoint exact
observations, fixed/one-parameter/all-free fits, original array-fit generic
compatibility, and NNLF. Identities cover inverse gamma and Weibull.
Extreme regressions cover two-ULP windows, tiny Weibull normalizers,
minimum normal gamma shape with a relative tiny-log-CDF comparison,
and the minimum/maximum signed powers.

Checked command (CMD, with existing output directory and FPM dependencies):

```cmd
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Jbuild\power_trunc-native -Ibuild\power_trunc-native -Ibuild\gfortran_2806A042DAA8BBF7 src\special\log_gamma.f90 src\empirical_distribution.f90 src\censored.f90 src\fit.f90 src\stats.f90 test\test_censored_power_trunc.f90 build\gfortran_3C159EB802C22FB9\scifort\libscifort.a -o build\power_trunc-native\test.exe
build\power_trunc-native\test.exe
```

The optimized run substitutes `-O3` and paths under
`build\power_trunc-optimized`, retaining all checks and traps. Changed
modules are compiled explicitly before linking the existing archive.
Ordinary scalar/adjacent/affine log comparisons allow 3e-12 scaled error,
ordinary intervals 3e-11, quantile round trips 5e-11, most identities
2e-12, and the relative minimum-shape comparison 1e-12. The measure is
`abs(actual-reference)/max(1,abs(reference))`, except the tiny log-CDF
is first divided by its expected magnitude. Invalid/infinite classifications
and exact support masses use no tolerance.

## Independent references

All 363 comparisons pass in both checked and optimized builds: 300 mixed
samples and 63 stresses. Maximum scaled error is `2.6691e-14`. Cached
reference inputs were rechecked against the final affine-guard changes.

```cmd
python validation\check_censored_power_trunc.py --module-dir build\gfortran_2806A042DAA8BBF7 --library build\gfortran_3C159EB802C22FB9\scifort\libscifort.a
```

The checker evaluates independent gamma/normal integrals and shifted
exponential identities with mpmath: 100 digits for 300 seeded mixed samples
and 1100 digits for extreme cases. Seed: 20261117. Narrow smooth intervals
use independently rescaled density quadrature with mpmath Gauss-Legendre,
maximum refinement degree 3. For a=1 generalized gamma, interval references
use the exact exponential identity at both signed powers, including very
steep adjacent tails. Finite comparisons require scaled error <=2e-11.

Checked driver flags are `-O0 -g -fcheck=all
-ffpe-trap=invalid,zero,overflow`; optimized flags are `-O3 -fcheck=bounds`.
Stresses include gamma shapes down to 1e-300, both power signs, 1e-100/1e100
tail points, scales 1e-300/1e300, normal truncation near -100/-81/40/80,
two-ULP truncation windows, huge Weibull hazards, and underflowed Weibull
normalizers with power 1e6 and upper bounds 1e-100/1e-300.

Independent SciPy bounded fits use exact `[0.5,0.8,1.1]`, left `[0.4]`,
right `[1.7]`, intervals `[0.6,0.7]` and `[0.9,1.0]`:

| Family | Free parameter | Fitted value | NNLF |
| --- | --- | ---: | ---: |
| Generalized gamma | a in [0.3,5], c=-1.3, loc=0, scale=1 | 1.7716401277384342 | 10.216323674112758 |
| Truncated normal | loc in [-0.1,0.1], a=-1, b=3, scale=1 | 0.09999999778083914 | 13.024419603423851 |
| Truncated Weibull | c in [0.3,5], a=0.1, b=3, loc=0, scale=1 | 1.6519866276629518 | 10.531775107244417 |

Native one-parameter fits allow 3e-5 scaled parameter error and 2e-8 scaled
NNLF error. The normal optimum is at the selected location bound; this
does not establish an unconstrained optimum. All-free checks require a
finite improving candidate without claiming global convergence.

## Full suite and limits

The shared censoring/fit/statistics integration changed, so the full suite
is run in addition to focused checked and optimized validation:

```cmd
fpm test
```

The full suite passed with exit code 0 in 1545.617 seconds (25 minutes
46 seconds), including all 2,125 new checks and the existing gamma-power,
censored-family, scalar-distribution, fitting, and integration regressions.
Existing array-temporary warnings appeared in interpolation and signal
tests and did not cause failures. No numerical code changed after this
successful full run began; subsequent edits only record validation results.

The [API guide](../docs/censored_fitting.md) records current shape/power,
finite-bound, affine-support, and narrow-window restrictions. Generalized
gamma requires normal positive a<=1e6 and 1e-6<=abs(c)<=1e6; normal
truncation bounds are finite and within [-100,100]; truncated Weibull
requires 1e-6<=c<=1e6. Very narrow standardized windows currently require
zero loc and a power-of-two scale. All truncation configurations require
normal nonzero scaled/physical endpoints and abs(loc)/(scale*(b-a))<=1e4.
General affine narrow windows and wider parameter ranges remain gaps.
Ill-conditioned affine models and log probabilities beyond binary64
range retain the documented limitations. Existing scalar routines,
optimizer, C ABI, and dependencies are unchanged.
