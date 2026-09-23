# Starter validation record

Validation date: 2026-08-20

## Completed checks

The starter source was checked with GNU Fortran 14.2.0 on Linux using Fortran
2018 mode, warnings, implicit-interface errors, bounds and runtime checks, and
a backtrace-enabled debug build.

The following checks passed:

- strict compilation of every Fortran module;
- execution of `test/test_scifort.f90`;
- optimized compilation and execution of the distribution example;
- construction of a shared library;
- compilation and execution of a C program against `include/scifort.h`;
- exported C-symbol inspection;
- TOML parsing of `fpm.toml`;
- repository ASCII-policy check;
- selected comparison of 948 results with SciPy 1.17.0.

The sampled SciPy comparison covered PDF, log-PDF, CDF, survival,
log-CDF, log-survival, PPF, and inverse-survival values. The largest sampled
absolute difference for a quantile was approximately `7.11e-15`. The largest
sampled relative difference for a normal probability or log-probability was
approximately `1.03e-13` in a far tail.

These checks are smoke and regression checks, not a complete accuracy proof or
a final accuracy specification.

## FPM validation status

The execution environment used to create this archive did not contain an FPM
binary and did not permit downloading one. Therefore, the exact `fpm build`
and `fpm test` commands were not executed locally. The source dependency order
and equivalent GNU Fortran builds were exercised directly.

The repository includes a GitHub Actions workflow that installs FPM 0.13.0 and
runs debug build, tests, release build, and the example on Linux, macOS, and
Windows. Run the following after unpacking to complete local FPM validation:

```text
fpm build --profile debug
fpm test --profile debug
fpm build --profile release
fpm run --example distributions_demo
```

## Incomplete gamma, gamma, and chi-square validation

Validation date: 2026-09-19

Environment: Windows 11, GNU Fortran 17.0.0 20260510 (experimental), FPM
0.12.0 alpha, Python 3.13 with mpmath 1.3.0 and SciPy 1.15.3.

Commands run, all passing (`test_gamma` and `test_scifort`):

```text
fpm test
fpm test --profile debug
fpm test --profile release
fpm test --flag "-O0 -g -fcheck=all -ffpe-trap=zero,overflow -fbacktrace -Wall -Wextra -std=f2018"
fpm run --example distributions_demo
python tools/check_ascii.py
```

The trapping build produced no compiler warnings and no division-by-zero or
overflow traps. It found one division by an underflowed Newton slope in the
inverse solver, which was fixed before this record.

Additional sweeps compiled the modules directly with `gfortran -O2` and
compared them with mpmath at 50 digits. The sweep scripts were not committed.

- `log1p_safe` and `expm1_safe`, 5000 arguments on `[-0.99, 5]` and
  magnitudes down to `1e-20`: maximum relative errors `5.6e-16` and
  `2.2e-16`. The 0.1.0 versions had `5.0e-13` and `3.0e-12`.
- `gammainc`, `gammaincc`, and their logarithms, 3000 random `(a, x)` with
  `1e-8 <= a <= 1e6`, concentrated near `x = a` and spread over
  `1e-10 <= x <= 1e3`: with `V` the tail and `x f(x) / V` its condition
  number, the maximum of
  `|error| / (V (1 + |log V| + x f(x) / V))` was `6.4e-16`.
- `gammaincinv` and `gammainccinv`, 1500 random `(a, p)` with
  `1e-6 <= a <= 1e5` and `p` down to `1e-300`: the relative error in `x`
  implied by the mpmath residual, divided by `1 + V / (x f(x))`, was at most
  `4.1e-14`. The worst case was `p = 7.7e-300`, where the logarithmic
  residual limits accuracy to about `eps |log p| / a`. All 740 results equal
  to zero were confirmed to have exact roots below the smallest positive
  normal number.
- Against SciPy 1.15.3 on 20000 random `(a, x)`, the median relative
  difference was `2.3e-16` for `P` and `1.1e-16` for `Q`. The maximum was
  `2.0e-12`. In each of the eight largest differences checked with mpmath,
  SciPy had the larger error (up to `2.0e-12`) and SciFort's error was at
  most `3.7e-14`.



Known limitations:

- Evaluation cost near `x = a` grows like `sqrt(a)`. No uniform asymptotic
  expansion for large `a` is implemented yet.
- Tail probabilities below about `1e-300` carry relative error of order
  `eps |log V|` from forming the prefactor in logarithmic form.
- Inverse values in the subnormal range are returned with reduced precision,
  and roots below them are returned as zero.
- Only GNU Fortran was tested. Intel `ifx` and LLVM Flang were not run.
- The C ABI exposes gamma and chi-square PDF, CDF, and PPF routines; the
  incomplete-gamma functions and the remaining distribution tails/logarithms
  are still native-Fortran-only.

## Incomplete beta, beta, Student t, and F validation

Validation date: 2026-09-19. Same environment as the previous section.

Commands run on the final code, all passing (`test_beta`, `test_gamma`, and
`test_scifort`):

```text
fpm test
fpm test --profile debug
fpm test --profile release
fpm test --flag "-O0 -g -fcheck=all -ffpe-trap=zero,overflow -fbacktrace -Wall -Wextra -Wno-integer-division -std=f2018"
fpm run --example distributions_demo
python tools/check_ascii.py
```

The trapping build produced no compiler warnings and no traps.
`-Wno-integer-division` silences only the notice for the constant
`huge(n) / 4` in the iteration limits.

Sweeps on the final code, with the modules compiled directly by
`gfortran -O2` and compared with mpmath (the sweep scripts were not
committed):

- `betainc` and `betaincc` with their logarithms, 1000 random `(a, b, x)`
  with parameters from `1e-6` to `1e5`, concentrated near the mean and in both
  tails: with `V` the tail, `f` the density, and `z = min(x, 1 - x)`, the
  maximum of `|error| / (V (1 + |log V| + z f / V))` was `3.8e-16`. For 21
  cases mpmath produced no reference that was stable across precisions.
- 150 cases with one parameter from `1e-7` to `0.3` and the other from `1e2`
  to `1e5` (the region that the gamma expansion now covers): the same
  normalized error was at most `2.4e-16`.
- Against SciPy 1.15.3 on 20000 random t, F, and beta cases, the median
  relative difference in `cdf` and `sf` was between `1.3e-16` and `3.8e-15`.
  The largest differences were checked with mpmath. In the worst t and F
  survival-function cases SciPy had the larger error (for F, `4.8e-11`
  against `5.8e-17`). The largest beta `cdf` difference, `6.3e-3`, was a
  SciPy error at a probability near `5.8e-271`; SciFort's error there was
  `2.7e-13`, consistent with `eps |log V|`.

The following sweeps were started on the final code but did not complete;
Claude Code stopped the background job because the system was low on
memory. Their results are therefore not claimed here:

- a second 1000-case forward sweep of the incomplete beta function;
- a 600-case sweep of `betaincinv` and `betainccinv`;
- a 900-case mpmath sweep of the t, F, and beta distribution functions.

Earlier runs of those sweeps on intermediate versions of the code found the
defects that led to the small-parameter complement, the gamma expansion, and
the F log-density rearrangement. The committed tests cover the inverses and
distribution functions against mpmath reference values with
conditioning-based tolerances.

Known limitations:

- For `|t|` beyond about `1e150 sqrt(df)`, and for the corresponding
  quantiles, the Student t tail uses the exact leading term; SciPy 1.15.3
  returns `0` and `-1e100` for `t.sf(1e200, 1)` and `t.ppf(1e-300, 1)`.
- Infinite F degrees of freedom return NaN.
- Continued-fraction evaluation near the transition point with both
  parameters large loses accuracy in proportion to `sqrt(a + b)`, which is
  within the conditioning of the problem.
- Only GNU Fortran was tested. The C ABI exposes beta, Student t, and F PDF,
  CDF, and PPF routines; their remaining tails and logarithmic routines are
  still native-Fortran-only.

## Poisson and binomial validation

Validation date: 2026-09-19. Same environment as the previous sections.

Commands run on the final code, all passing (`test_discrete`, `test_beta`,
`test_gamma`, and `test_scifort`):

```text
fpm test
fpm test --profile debug
fpm test --flag "-O0 -g -fcheck=all -ffpe-trap=zero,overflow -fbacktrace -Wall -Wextra -Wno-integer-division -std=f2018"
python tools/check_ascii.py
```

The trapping build produced no compiler warnings and no traps.

`test/test_discrete.f90` compares against mpmath reference values whose
tolerances combine the accuracy of the underlying special function with the
sensitivity of each probability to the rounding of `mu` or `p`. It also
checks that the mass function sums to the cumulative function and to one,
that each quantile is the smallest count reaching its probability, the
identities with the incomplete gamma and beta functions, the reflection
identity of the binomial distribution at `p = 1/4`, the degenerate cases
`mu = 0`, `p = 0`, `p = 1`, and `n = 0`, and that the integer and real
interfaces agree.

Comparison with SciPy 1.15.3 on 4000 random cases, with `mu` from `1e-3` to
`1e6`, `n` up to `1e6`, and `p` from `1e-5` to `1`, covering both the bulk
and the tails:

| quantity | median difference | largest difference |
| --- | --- | --- |
| Poisson pmf | `3.6e-15` | `2.3e-9` |
| Poisson cdf | `0` | `3.8e-12` |
| Poisson sf | `1.9e-16` | `3.9e-6` |
| binomial pmf | `1.4e-15` | `2.9e-11` |
| binomial cdf | `0` | `2.9e-11` |
| binomial sf | `4.4e-16` | `2.9e-11` |

Each of those six largest differences was checked against mpmath at 60
digits. In every case SciFort was the more accurate of the two, with errors
from `1.4e-15` to `1.7e-13` against SciPy's `2.3e-9` to `3.9e-6`. The largest
SciPy error was in `poisson.sf(981157, 976226)`.

Known limitations:

- The quantile search evaluates the cumulative function about
  `log2(range)` times; no closed-form starting bound is used beyond the
  normal approximation.
- Counts are held in binary64, so a count above `2**53` cannot be
  distinguished from its neighbors.
- Only GNU Fortran was tested. The C ABI exposes Poisson and binomial PMF,
  CDF, and PPF routines; the remaining tails and logarithmic routines are
  still native-Fortran-only.

## Bernoulli, geometric, and negative-binomial validation

Validation date: 2026-09-20.

Environment: Linux x86-64, GNU Fortran 14.2.0, Python with SciPy 1.17.0,
NumPy 2.3.5, and mpmath 1.3.0. FPM was not installed in this validation
environment, so the final tests were compiled and linked directly with
`gfortran` using the same module dependency order as the FPM manifest.

The complete distribution test set was rebuilt from a clean directory with
runtime checking, floating-point traps, and warnings promoted to errors. The
following options were used:

```text
-std=f2018 -O0 -g -fcheck=all -ffpe-trap=zero,overflow -fbacktrace \
-Wall -Wextra -Werror -Wno-compare-reals -Wno-integer-division
```

All of these programs passed on the final code:

```text
test_scifort
test_gamma
test_beta
test_discrete
test_discrete_extended
test_continuous
```

`python tools/check_ascii.py` also passed. The generated extended-discrete
reference file was reproduced byte-for-byte by
`tools/generate_discrete_extended_reference.py`.

`test/test_discrete_extended.f90` checks Bernoulli against binomial with one
trial, geometric against negative binomial with `n = 1`, shifted supports,
integer and real count interfaces, elemental array evaluation, parameter and
probability errors, degenerate `p = 1` cases, PMF sums, CDF identities, and
that each PPF/ISF result is the smallest integer satisfying its defining tail
inequality. Geometric and negative-binomial reference values are generated
with mpmath at 70 decimal digits.

A separate 2000-case random comparison of geometric and negative-binomial
forward probabilities against SciPy 1.17.0 gave median relative differences
of about `3e-15` for PMF and `1.5e-16` for CDF/SF. The largest absolute CDF or
SF difference was below `8.8e-14`; the largest absolute PMF difference was
below `5.6e-16`.

For quantiles, 1995 of 2000 PPF values and 1958 of 2000 ISF values were
identical to SciPy. Representative disagreements in the far tails were
recomputed independently with mpmath at 60 digits. In those checked cases the
SciFort result was the mathematically correct smallest count, while SciPy
1.17.0 differed by one or more counts. For example, with geometric
`p = 6.925563431893058e-7` and survival probability
`3.762754918489501e-12`, the high-precision result and SciFort both give
`37983712`, while SciPy 1.17.0 gives `37983706`.

Known limitations:

- Negative-binomial quantiles use a normal starting approximation followed by
  a doubling bracket and integer bisection, so their cost is logarithmic in
  the count range rather than constant time.
- Counts are represented in binary64 in the real-valued interfaces; above
  `2**53`, adjacent integer counts cannot all be represented distinctly.
- Only GNU Fortran was tested for this addition. Intel `ifx` and LLVM Flang
  were not run in this environment.
- The C ABI exposes Bernoulli, geometric, and negative-binomial PMF, CDF,
  and PPF routines; the remaining tails and logarithmic routines are still
  native-Fortran-only.

## Expanded scalar C ABI validation

Validation date: 2026-09-20.

Environment: Linux x86-64, GNU Fortran 14.2.0 and GCC 14.2.0. FPM was not
installed in this environment, so the ABI layer was rebuilt directly with
GNU Fortran in module dependency order.

The C ABI now exposes scalar PDF or PMF, CDF, and PPF entry points for all 20
implemented distributions. `test/test_c_api.f90` calls every scalar entry
point and compares it with the corresponding native Fortran result, including
non-default locations/scales and representative shape parameters. It also
checks invalid gamma and binomial parameters return NaN.

The source and ABI regression test were compiled with:

```text
-std=f2018 -O0 -g -fcheck=all -ffpe-trap=zero,overflow -fbacktrace \
-Wall -Wextra -Werror -Wno-compare-reals -Wno-integer-division
```

A separate C11 smoke program included `include/scifort.h`, called every new
scalar symbol, and linked against the Fortran objects. The C translation unit
was compiled with:

```text
gcc -std=c11 -Wall -Wextra -Werror -pedantic -Iinclude
```

Both the Fortran ABI regression test and the C smoke program passed. The C
header also passed a standalone syntax-only compilation. No C ABI entry points
for SF, log-PDF/log-PMF, log-CDF, log-SF, ISF, or the special functions have
been added yet.

## Explicit-state RNG and random-variate validation

Validation date: 2026-09-20.

Environment: Linux x86-64, GNU Fortran 14.2.0. FPM was not installed in this
environment, so the RNG and distribution samplers were compiled directly in
the same dependency order as the rest of the library.

`test/test_random.f90` fixes the first 12 open-interval binary64 uniforms from
the default six-component state, then checks the exact six-component state
after those draws. An explicit scalar seed of `12345` is also required to
reproduce the default state stream. Default-integer and `int64` seed overloads,
bulk uniform filling, state checkpoint/restore, invalid-state rejection, and
non-advancement on invalid distribution parameters are all regression tested.

A deterministic 100000-draw smoke test requires every uniform to lie strictly
inside `(0, 1)`, the sample mean to be within `0.005` of `0.5`, and the sample
variance to be within `0.003` of `1/12`. This is a regression and gross-error
check, not a substitute for a dedicated TestU01-style statistical assessment.

Every one of the 20 current distributions is exercised through both its
scalar `*_rvs` function and array `*_rvs_array` subroutine. Synchronized RNG
states are used to require exact equality between each sampler result and the
corresponding PPF evaluated at the same generated uniform. Invalid gamma and
binomial parameter cases return NaNs without consuming generator state.

The generator and sampler tests were compiled with:

```text
-std=f2018 -O0 -g -fcheck=all -ffpe-trap=zero,overflow -fbacktrace \
-Wall -Wextra -Werror -Wno-compare-reals -Wno-integer-division
```

Known limitations:

- The current random-variate layer uses inverse transforms for every family.
  This provides a small and consistent first implementation but is slower than
  specialized sampling algorithms, especially when a PPF performs iterative
  inversion or discrete quantile search.
- The project-local 52-bit conversion from MRG32k3a recurrence outputs is
  regression tested and uses equal-width rejection buckets, but no large
  empirical RNG test battery was run in this environment.
- Random-state and random-variate entry points are not yet exposed through the
  C ABI.
- Only GNU Fortran was run locally for this addition; Intel `ifx` and LLVM
  Flang were not available in this environment.

## Likelihood aggregation validation

Validation date: 2026-09-20.

`test/test_likelihood.f90` exercises `*_loglikelihood` for all 20 implemented
distributions and verifies exact agreement with an explicit sum of the
corresponding elemental `logpdf` or `logpmf` values. Representative `*_nnlf`
results are checked as the exact negative of the log-likelihood. All five
discrete families are tested through both integer and `real(dp)` observation
arrays, with the binomial integer interface also using an integer trial count.

The test additionally checks the empty-sample convention (log-likelihood and
negative log-likelihood both compare equal to zero) and NaN propagation for
invalid normal and gamma parameters. It is built with the same strict GNU
Fortran flags shown in the RNG validation section.

Known limitations:

- The current likelihood layer is an unweighted sum for independent,
  uncensored observations. It does not yet implement observation weights or
  censored-data likelihoods.
- Likelihood routines are native-Fortran-only and are not yet part of the C
  ABI.

## Analytic score validation

Validation date: 2026-09-20.

`test/test_score.f90` exercises the analytic likelihood gradients for all 20
implemented distribution families. It compares every regular score component
with a central finite-difference derivative of the independent
`*_loglikelihood` implementation. The comparisons include the shape
derivatives of gamma, beta, Student t, F, Weibull, Pareto, and negative
binomial distributions as well as location and scale derivatives.

The test separately checks five positive-argument digamma values against
high-precision reference constants, integer and `real(dp)` discrete-score
overload equivalence, and NaN behavior at selected nonregular support
boundaries. The positive-argument digamma helper uses recurrence followed by
the DLMF asymptotic series rather than a finite-difference approximation.

Known limitations:

- Scores are for independent, uncensored, equally weighted observations.
- Discrete scores intentionally omit integer support shifts and the binomial
  trial count. The negative-binomial score differentiates the positive-real
  extension of `n`, even though the fitting interface constrains `n` to an
  integer.
- Moving-support likelihoods can be nonregular at support boundaries. The
  affected components return NaN there.
- Score routines are native-Fortran-only and are not yet in the C ABI.

## Bounded maximum-likelihood fitting validation

Validation date: 2026-09-20.

`test/test_fit.f90` exercises all 20 family-specific fit wrappers. A free
normal fit is compared with the closed-form normal MLE. Fixed-location
exponential, Bernoulli, Poisson, geometric, binomial, and negative-binomial
fits are compared with their corresponding closed-form parameter estimates.
The test also verifies integer projection when binomial `n` is free, equal
lower/upper bounds as fixed parameters, integer and `real(dp)` observation
overloads for all discrete families, fixed-parameter smoke fits for every
continuous family, and `fit_status_no_finite_objective` for an invalid fixed
normal scale.

The optimizer is a deterministic projected coordinate/pattern search with no
hidden global state or procedure-pointer callbacks. A final object-file check
with `readelf` confirms that the GNU Fortran `fit.o` does not request an
executable stack.

Known limitations:

- The first optimizer is local and derivative-free; it does not guarantee a
  global optimum. Finite bounds are required for every parameter.
- The fitter does not yet support observation weights, censored likelihoods,
  covariance/Hessian estimates, or alternative optimizers.
- Fitting routines are native-Fortran-only and are not yet in the C ABI.

## 2026-09-20 combined checkpoint rebuild

After the expanded C ABI, explicit RNG, random-variate, likelihood, score, and
fitting changes were combined, the source was rebuilt from clean directories
in both strict and optimized configurations. These eleven programs passed in
each build:

```text
test_scifort
test_gamma
test_beta
test_discrete
test_discrete_extended
test_continuous
test_c_api
test_random
test_likelihood
test_score
test_fit
```

The strict build used the runtime-checking/warnings-as-errors flags documented
above. The optimized build used `-std=f2018 -O2`. A separate C11 smoke program
including `include/scifort.h` compiled with warnings as errors, linked to the
optimized Fortran objects, and passed. `python tools/check_ascii.py` also
passed. No object, module, or executable build artifacts were present in the
source tree. FPM was not installed in this local validation environment, so
the score/fitting additions in this checkpoint were not independently run
through `fpm test` here.

## Descriptive statistics and rank validation

Validation date: 2026-09-20. GNU Fortran 14.2.0 on Linux.

`test/test_descriptive.f90` compares mean, variance, standard deviation,
central moments through order four, linear quantiles, covariance, Pearson
correlation, and all five `rankdata` methods against values reproduced by
`tools/generate_statistics_reference.py` with NumPy 2.3.5 and SciPy 1.17.0.
The test also covers empty arrays, invalid `ddof`, invalid quantiles and rank
methods, NaN propagation, infinities in order statistics, and data near the
binary64 overflow limit.

The implementation scales finite data before centered powers and products.
The strict floating-point-trap build verifies that a sample containing
`+1e308` and `-1e308` can return a finite standard deviation and an infinite
variance without triggering an intermediate overflow exception.

Known limitations:

- NaN omission policies and masked data are not implemented; NaNs propagate.
- `quantile` currently implements only linear / Hyndman-Fan type 7
  interpolation.
- Descriptive reductions are one-dimensional and unweighted.

## Classical hypothesis-test validation

Validation date: 2026-09-20. Same compiler and platform.

`test/test_hypothesis.f90` compares one-sample, pooled independent, Welch,
and paired t-tests; Pearson and Spearman correlation tests; and asymptotic
Mann-Whitney U tests against SciPy 1.17.0. The comparison includes
`two-sided`, `less`, and `greater` alternatives, tied ranks, continuity
correction on/off for Mann-Whitney, Welch-Satterthwaite degrees of freedom,
constant samples, invalid alternatives, and mismatched sample sizes.

Pearson p-values use the exact beta null distribution for `n > 2` and the
`n = 2` special case. Spearman p-values use the usual Student-t
transformation. Mann-Whitney uses the tie-corrected normal approximation and
reports the U statistic of the first sample.

Known limitations:

- Exact small-sample Mann-Whitney p-values are not implemented.
- Spearman p-values are asymptotic rather than exact/permutation p-values.
- NaN omission policies are not implemented.

## Resampling validation

Validation date: 2026-09-22. GNU Fortran 14.2.0 direct-build environment.

`test/test_resampling.f90` retains the original bootstrap/RNG regression and
now verifies SciPy-style automatic switching to an exact independent-sample
permutation test when all 20 partitions of the 3-by-3 example can be
enumerated.

`test/test_resampling_inference.f90` covers the expanded framework:

- one-sample percentile, basic, and BCa bootstrap intervals and standard errors;
- SciPy-compatible one-sided bootstrap endpoint conventions;
- exact `independent`, `samples`, and `pairings` permutation distributions;
- preservation of RNG state for exact tests;
- one- and two-sample Monte Carlo tests with the plus-one correction;
- one- and two-sample power simulation and significance-threshold behavior;
- invalid-input handling without RNG advancement.

The exact permutation, deterministic Monte Carlo, and deterministic power
reference values are generated by SciPy 1.17.0 with
`tools/generate_resampling_inference_reference.py`. The generated reference
module is byte-for-byte reproducible. Random bootstrap streams intentionally use
SciFort's MRG32k3a generator rather than NumPy's generator, so bootstrap draw
sequences are validated by deterministic SciFort replay. The percentile/basic/
BCa interval arithmetic was additionally cross-checked against SciPy 1.17 using
the exact SciFort-generated bootstrap distribution; the interval formulas agree.

Bootstrap and random-permutation indices are formed by rejection from the RNG's
52-bit uniform grid, making each selected index equiprobable on that grid.
Randomized permutation and Monte Carlo p-values use the conservative plus-one
correction. Exact permutation tests omit this adjustment.

Current native-Fortran scope is one-dimensional real samples with scalar
callbacks. SciPy's array-axis/vectorized/batch broadcasting API and reuse of a
previous `bootstrap_result` to append additional draws are not implemented.

## 2026-09-20 Milestone 5 checkpoint rebuild

The final Milestone 5 candidate was compiled directly with GNU Fortran 14.2.0
because FPM is not installed in this validation environment. The user reported
that `fpm build` and `fpm test` succeeded for the immediately preceding
scores/fitting checkpoint; that result is not claimed for the new Milestone 5
files.

A strict build compiled all source files and all fourteen test programs with:

```text
-std=f2018 -O0 -g -fcheck=all -ffpe-trap=zero,overflow -fbacktrace
-Wall -Wextra -Werror -Wno-compare-reals -Wno-integer-division
```

The fourteen passing programs were:

```text
test_scifort
test_gamma
test_beta
test_discrete
test_discrete_extended
test_continuous
test_c_api
test_random
test_likelihood
test_score
test_fit
test_descriptive
test_hypothesis
test_resampling
```

A second clean build with `-std=f2018 -O2` and the same warning settings also
passed all fourteen programs. A standalone C11 smoke translation unit compiled
with `-Wall -Wextra -Werror -pedantic`, linked against the optimized Fortran
objects, and successfully called the normal CDF and Poisson PMF ABI entries.

`python tools/check_ascii.py` and `python tools/generate_statistics_reference.py`
were also run. No generated object, module, executable, or temporary reference
files are included in the source archive.

## Milestone 6: core special functions

The Milestone 6 checkpoint adds `test/test_special_elementary.f90` and
`tools/generate_special_reference.py`. Reference values were generated with
SciPy 1.17.0 and NumPy 2.3.5. The committed tests cover error and inverse-error
functions; log-gamma, log-beta, and digamma; normal-CDF and logistic transforms;
entropy helpers; cancellation-resistant elementary functions; Box-Cox forward
and inverse transforms; and one-dimensional log-sum-exp/softmax reductions.

One deliberately asymmetric beta case, `betaln(1e-3, 2000)`, is checked against
an 80-digit mpmath value rather than the SciPy result. SciFort returns
`6.899578232695083` in binary64, consistent with the high-precision value
`6.899578232695082464...`; SciPy 1.17.0 returned
`6.899578232692875` in this environment. This is documented as a validation
difference, not normalized away for parity.

A separate 1000-point random sweep over real arguments in `[-50, 50]`, excluding
points within `1e-4` of negative integer digamma/gamma poles, compared the final
code with SciPy 1.17.0. The largest absolute `gammaln` difference was
`5.7e-14`; the largest digamma difference scaled by `max(1, |psi(x)|)` was
`1.4e-14`, occurring close to a pole. Box-Cox, Box-Cox-1p, and `exprel` stayed
within a few units of binary64 roundoff on the same sweep. The sweep harness was
not committed.

The complete source tree was rebuilt from clean directories with the strict GNU
Fortran flags used for Milestone 5 and separately with `-std=f2018 -O2`. All
fifteen test programs passed in both configurations. The shipped application
and distribution example also compiled and ran against the optimized objects,
and an independent C11 program including `include/scifort.h` compiled with
warnings as errors, linked to the optimized Fortran objects, and ran
successfully. `python tools/check_ascii.py .`, the special-reference generator,
the 132-column source check, and the new-source FORD dummy-declaration check
also passed. FPM and `fprettify` are not installed in this validation
environment; the user reported that `fpm build` and `fpm test` passed on the
immediately preceding Milestone 5 checkpoint.

Known limitations:

- public `betaln` currently promises only finite positive shape arguments;
- inverse error functions inherit the current bisection-based normal-quantile
  speed and far-tail range;
- `logsumexp`, `softmax`, and `log_softmax` are one-dimensional and unweighted;
- complex-valued SciPy special functions are outside the current scope.


## Milestone 7 compact extension: arcsine and half-normal

Validation date: 2026-09-20. This deliberately small checkpoint adds only the
arcsine and half-normal families. `test/test_more_continuous.f90` checks SciPy
1.17.0 PDF/CDF/SF reference values, support/endpoints, quantile round trips,
likelihood sums, analytic scores against centered differences, deterministic
RNG mapping, fixed-parameter fit dispatch, and scalar C-ABI equality.

A clean direct GNU Fortran 14.2.0 build used the strict warning/runtime flags
shown for Milestone 5. All 17 test programs passed. FPM is not installed in
this validation environment; the user reported that `fpm build` and `fpm test`
passed on the immediately preceding Milestone 7 checkpoint.


## Milestone 7 compact extension: half-Cauchy and Lomax

Validation date: 2026-09-20. This checkpoint deliberately adds only two
continuous families. `test/test_additional_continuous.f90` checks SciPy 1.17.0
PDF/CDF reference values, support/endpoints, quantile and survival round trips,
likelihood reductions, analytic scores against centered finite differences,
deterministic RNG mapping, fixed-parameter fit dispatch, and scalar C-ABI
equality.

A clean direct GNU Fortran 14.2.0 strict build used the same warning/runtime
flags as the preceding Milestone 7 checkpoint. All 18 test programs passed.
The user reported that `fpm build` and `fpm test` passed for the preceding
arcsine/half-normal checkpoint; FPM is not installed in this validation
environment, so FPM success is not claimed for these two newly added families.


## Milestone 7 compact extension: chi and Maxwell

Validation date: 2026-09-20. This checkpoint deliberately adds only the chi
and Maxwell continuous families. `test/test_chi_maxwell.f90` checks SciPy
1.17.0 PDF/CDF reference values, support/endpoints, quantile and survival
round trips, the Maxwell = chi(df=3) identity, likelihood reductions,
analytic scores against centered finite differences, deterministic RNG
mapping, fixed-parameter fit dispatch, and scalar C-ABI equality.

A clean direct GNU Fortran 14.2.0 strict build used the same warning/runtime
flags as the preceding compact checkpoints. All 19 test programs passed. The
user reported that `fpm build` and `fpm test` passed for the preceding
half-Cauchy/Lomax checkpoint; FPM is not installed in this validation
environment, so FPM success is not claimed for the newly added chi/Maxwell
code.

## Milestone 7 compact extension: cosine and semicircular

The bounded cosine and semicircular families are covered by
`test/test_cosine_semicircular.f90`. The test checks SciPy 1.17.0 PDF/CDF
references, CDF/SF and PPF/ISF round trips, likelihood reductions, analytic
score finite differences, deterministic inverse-CDF RNG mapping, fixed-bound
fitting dispatch, and scalar C-ABI equality. The complete direct GNU Fortran
regression suite contains 20 test programs after this extension.


## Milestone 7 compact extension: Anglit and Moyal

Validation date: 2026-09-20. This checkpoint deliberately adds only the Anglit
and Moyal loc/scale families. `test/test_anglit_moyal.f90` checks SciPy 1.17.0
PDF/CDF reference values, CDF/SF and PPF/ISF round trips, likelihood reductions,
analytic scores against centered finite differences, deterministic inverse-CDF
RNG mapping, fixed-parameter fit dispatch, and scalar C-ABI equality.

A clean direct GNU Fortran 14.2.0 strict build used the same warning/runtime
flags as the preceding compact checkpoints, and all 21 test programs passed.
All 21 programs also passed in a clean `-O2` build. A standalone C11 smoke
program compiled with warnings as errors, linked against the optimized Fortran
objects, and called the new Anglit/Moyal ABI entries successfully. FPM is not
installed in this validation environment, so FPM success is not claimed for
this checkpoint.


## Milestone 7 compact extension: hyperbolic secant and half logistic

Validation date: 2026-09-20. This checkpoint deliberately adds only the
hyperbolic-secant and half-logistic loc/scale families.
`test/test_hypsecant_halflogistic.f90` checks SciPy 1.17.0 PDF/CDF reference
values, CDF/SF and PPF/ISF round trips, deep logarithmic tails, likelihood
reductions, analytic scores against centered finite differences, deterministic
inverse-CDF RNG mapping, fixed-parameter fit dispatch, and scalar C-ABI
equality.

A clean direct GNU Fortran 14.2.0 strict build uses `-std=f2018`, `-Wall`,
`-Wextra`, warnings as errors, `-fcheck=all`, and floating-point zero/overflow
traps. All 22 test programs pass. The new hyperbolic-secant/half-logistic regression
also passes in an `-O2` build, and a standalone C11 smoke program compiles
against `scifort.h`, links to the optimized Fortran objects, and calls both new
ABI families successfully. FPM is not installed in this validation environment,
so FPM success is not claimed for this checkpoint.


## Milestone 7 larger positive-support extension

Validation date: 2026-09-21. This checkpoint adds inverse-gamma,
inverse-Gaussian, Levy, and log-Laplace families.
`test/test_inverse_positive.f90` checks SciPy 1.17.0 PDF/CDF/PPF reference
values, tail and quantile round trips, likelihood reductions, analytic scores
against centered finite differences, deterministic inverse-CDF RNG mapping,
fixed-parameter fit dispatch, and scalar C-ABI equality.

A clean direct GNU Fortran 14.2.0 strict build used `-std=f2018`, `-Wall`,
`-Wextra`, warnings as errors, `-fcheck=all`, and floating-point zero/overflow
traps. All 23 test programs passed. The user reported that `fpm build` and
`fpm test` passed for the immediately preceding 37-family checkpoint; FPM is
not installed in this validation environment, so FPM success is not claimed
for these four newly added families.

## Bradford, truncated exponential, Fisk, and double-Weibull validation

Validation date: 2026-09-21.

Environment: Linux x86-64, GNU Fortran 14.2.0, Python with SciPy 1.17.0 and
NumPy 2.3.5. FPM was not installed in this validation environment, so the
final source and tests were compiled directly with GNU Fortran in module
dependency order.

The final checked source passed all 24 test programs with Fortran 2018 mode,
`-fcheck=all`, floating-point zero/overflow traps, backtraces, warnings, and
warnings-as-errors. `test/test_closed_form_continuous.f90` covers SciPy 1.17.0
PDF, CDF, survival, logarithmic-tail, PPF, and ISF reference values for the four
new families, plus tail/quantile round trips, analytic-score finite differences,
deterministic inverse-transform RNG checks, fixed-parameter fit dispatch, and
Fortran-to-C ABI equality.

The four new families were also rebuilt at `-O2`; the combined regression test
passed, and a standalone C11 program compiled against `include/scifort.h`,
linked to the optimized Fortran objects, and successfully called all four new
C-ABI families. `python tools/check_ascii.py` and the 132-column Fortran source
check passed. The numerical reference table is regenerated by
`tools/generate_closed_form_reference.py` using SciPy 1.17.0.

Known limitations:

- Bradford and truncated-exponential likelihood derivatives are nonregular when
  an observation lies exactly on a support boundary; their score routines
  return NaN for such boundary samples.
- Sampling uses inverse transforms for consistency with the other current
  families rather than specialized direct generators.
- This validation did not run FPM locally; run `fpm build` and `fpm test` after
  unpacking to confirm the new 45-family checkpoint in an FPM environment.

## Alpha, fatigue-life, generalized-logistic, and generalized-normal validation

Validation date: 2026-09-21.

This checkpoint adds four one-shape continuous families: Alpha,
Birnbaum-Saunders/fatigue-life, generalized logistic, and generalized normal.
`test/test_shape_family_continuous.f90` checks SciPy 1.17.0 PDF, CDF,
survival, logarithmic-tail, PPF, and ISF reference values, tail/quantile round
trips, likelihood reductions, analytic scores against centered finite
differences, deterministic inverse-transform RNG mapping, fixed-parameter fit
dispatch, and scalar C-ABI equality. Numerical references are regenerated by
`tools/generate_shape_family_reference.py`.

A clean direct GNU Fortran 14.2.0 checked build used Fortran 2018 mode,
`-fcheck=all`, floating-point zero/overflow traps, backtraces, warnings, and
warnings-as-errors (with the repository's existing real-comparison and
constant-integer-division warnings left nonfatal). All 25 test programs passed. The combined four-family regression also passes
in an `-O2` build. A standalone C11 program compiles against
`include/scifort.h`, links to the optimized Fortran objects, and calls the
Alpha, fatigue-life, generalized-logistic, and generalized-normal C ABI
functions successfully. Regenerating `test/shape_family_reference.f90` with
`tools/generate_shape_family_reference.py` reproduces the committed file
byte-for-byte. ASCII, 132-column Fortran-line, and clean-source-tree checks
also pass.

FPM is not installed in this validation environment, so the four newly added
families still require external `fpm build` and `fpm test` confirmation.

Known limitations:

- Random sampling uses inverse transforms for consistency with the current
  SciFort RNG layer rather than specialized direct generators.
- The fitting layer remains a bounded local projected coordinate/pattern
  search and requires finite bounds.
- Generalized-normal location scores are nonregular at an observation exactly
  equal to `loc` when `beta <= 1`; the score routine returns NaN in that case.

## Nakagami, power-normal, log-gamma, and Wald validation

Validation date: 2026-09-21.

This checkpoint adds Nakagami, power-normal, log-gamma, and Wald continuous
families and integrates them with the distribution, likelihood, analytic-score,
bounded-fit, explicit-state RNG, `scifort_stats`, and scalar C-ABI layers.
`test/test_transform_family_continuous.f90` checks SciPy 1.17.0 PDF, CDF,
survival, logarithmic-tail, PPF, and ISF reference values, tail/quantile round
trips, analytic scores against centered finite differences, deterministic
inverse-transform RNG mapping, fixed-parameter fit dispatch, the Wald identity
with inverse Gaussian shape one, and Fortran-to-C ABI equality.

A fresh exact-source GNU Fortran 14.2.0 strict build compiled all 73 library
sources in Fortran 2018 mode with `-fcheck=all`, floating-point zero/overflow
traps, warnings, and warnings-as-errors (with the repository's existing
real-comparison and constant-integer-division warnings left nonfatal). All
26 test programs passed. A separate `-O2` build of all 73 library sources also
passed the combined four-family regression. The shipped application and example
compiled and ran against that optimized build.

A standalone C11 program compiled against `include/scifort.h`, linked to the
optimized Fortran objects, and successfully called the Nakagami, power-normal,
log-gamma, and Wald C ABI entries. Regenerating
`test/transform_family_reference.f90` with
`tools/generate_transform_family_reference.py` reproduces the committed file
byte-for-byte using SciPy 1.17.0 reference outputs.

FPM is not installed in this validation environment, so this 53-family
checkpoint still requires external `fpm build` and `fpm test` confirmation.

Known limitations:

- Random sampling uses inverse transforms for consistency with the current
  SciFort RNG layer rather than specialized direct generators.
- The fitting layer remains a bounded local projected coordinate/pattern search
  and requires finite bounds.
- Wald deliberately exposes only `loc` and `scale`; it is the shape-one special
  case of the already supported inverse-Gaussian family.

## Gompertz, inverse-Weibull, beta-prime, and Burr XII validation

Validation date: 2026-09-21.

This checkpoint adds Gompertz, inverse-Weibull, beta-prime, and Burr XII
continuous families and integrates them with the distribution, likelihood,
analytic-score, bounded-fit, explicit-state RNG, `scifort_stats`, and scalar
C-ABI layers. `test/test_positive_shape_continuous.f90` checks SciPy 1.17.0
PDF, CDF, survival, logarithmic-tail, PPF, and ISF reference values,
tail/quantile round trips, likelihood reductions, analytic scores against
centered finite differences, deterministic inverse-transform RNG mapping,
fixed-parameter fit dispatch, and Fortran-to-C ABI equality.

A fresh exact-source GNU Fortran 14.2.0 strict build compiled all 77 library
sources in Fortran 2018 mode with `-fcheck=all`, floating-point zero/overflow
traps, backtraces, warnings, and warnings-as-errors (with the repository's
existing real-comparison and constant-integer-division warnings left
nonfatal). All 27 test programs passed. A separate `-O2` build of all 77
library sources also passed all 27 test programs.

A standalone C11 program compiled against `include/scifort.h`, linked to the
optimized Fortran objects, and successfully called the Gompertz,
inverse-Weibull, beta-prime, and Burr XII C ABI entries. Regenerating
`test/positive_shape_reference.f90` with
`tools/generate_positive_shape_reference.py` reproduces the committed file
byte-for-byte using SciPy 1.17.0 numerical outputs. The shipped application and
example also compile and run against the optimized objects. ASCII, 132-column,
new-source one-statement-per-line, FORD dummy-comment, and clean-source-tree
checks pass.

The strict build exposed an overflow trap in an initial Burr XII guard whose
compound `.and.` expression assumed short-circuit evaluation. Fortran does not
guarantee short-circuit evaluation, so the final implementation uses nested
branches before potentially overflowing `huge()/x` comparisons. The same
pattern was hardened in the Gompertz cumulative-hazard guard.

FPM is not installed in this validation environment, so this 57-family
checkpoint still requires external `fpm build` and `fpm test` confirmation.

Known limitations:

- Random sampling uses inverse transforms for consistency with the current
  SciFort RNG layer rather than specialized direct generators.
- The fitting layer remains a bounded local projected coordinate/pattern search
  and requires finite bounds.
- For beta-prime and Burr XII, location-score derivatives are nonregular when
  observations lie exactly on the moving lower support boundary; the score
  interface follows the project's existing boundary-NaN policy.


## Generalized half-logistic, exponential-power, exponentiated-Weibull, and power-lognormal validation

Validation date: 2026-09-21.

This checkpoint adds generalized half-logistic, exponential-power,
exponentiated-Weibull, and power-lognormal continuous families and integrates
them with the distribution, likelihood, analytic-score, bounded-fit,
explicit-state RNG, `scifort_stats`, and scalar C-ABI layers.
`test/test_transform_shape_four.f90` checks SciPy 1.17.0 PDF, CDF, and PPF
reference values, analytic scores against centered finite differences,
deterministic inverse-transform RNG mapping, fixed-parameter fit dispatch, and
Fortran-to-C ABI equality.

A fresh exact-source GNU Fortran 14.2.0 build compiled all 81 library sources
with `-std=f2008 -fcheck=all -ffpe-trap=zero,overflow` and ran all 28 test
programs successfully. A separate `-O2` build of the same 81 library sources
also ran all 28 test programs successfully.

A standalone C11 program compiled against `include/scifort.h`, linked to the
optimized Fortran objects, and successfully called one entry point from each
of the four new C-ABI families. The shipped application and example also
compiled and ran against the optimized objects. ASCII, 132-column,
new-source one-statement-per-line, and clean-source-tree checks pass.

FPM is not installed in this validation environment, so this 61-family
checkpoint still requires external `fpm build` and `fpm test` confirmation.

Known limitations:

- Random sampling uses inverse transforms for consistency with the current
  SciFort RNG layer rather than specialized direct generators.
- The fitting layer remains a bounded local projected coordinate/pattern search
  and requires finite bounds.
- The generalized half-logistic support endpoint depends on the shape
  parameter, so the ordinary likelihood derivative is nonregular when an
  observation lies exactly on that moving endpoint; the score follows the
  project's existing boundary-NaN convention in such cases.

## Left-Levy, Weibull maximum, R-distribution, and skew-Cauchy validation

Validation date: 2026-09-21.

This checkpoint adds left-Levy, Weibull maximum, R-distribution, and
skew-Cauchy continuous families and integrates them with distribution,
likelihood, analytic-score, bounded-fit, explicit-state RNG, `scifort_stats`,
and scalar C-ABI layers. `test/test_reflected_shape_four.f90` checks SciPy
1.17.0 PDF, CDF, PPF, and selected logarithmic-tail reference values,
quantile/tail identities, likelihood reductions, analytic scores against
centered finite differences, deterministic inverse-transform RNG mapping,
fixed-parameter fit dispatch, and Fortran-to-C ABI equality.

A fresh exact-source GNU Fortran 14.2.0 strict build compiled all 85 library
sources in Fortran 2018 mode with `-fcheck=all`, floating-point zero/overflow
traps, backtraces, warnings, and warnings-as-errors (with the repository's
existing real-comparison and constant-integer-division warnings left
nonfatal). All 29 test programs passed. A separate `-O2` build of all 85
library sources also passed all 29 test programs.

A standalone C11 program compiled against `include/scifort.h`, linked to the
optimized Fortran objects, and successfully called one entry point from each
of the four new C-ABI families. The shipped application and example also
compile and run against the optimized objects. ASCII, 132-column,
one-statement-per-line, and clean-source-tree checks pass.

At probability `1e-6`, the left-Levy quantile is about `-6.37e11`; SciFort and
SciPy 1.17.0 differ there by about `3.3e-10` relative because the two libraries
use independent normal-quantile kernels. The regression therefore uses a
slightly wider relative tolerance only for that extreme-tail value while
retaining tight tolerances elsewhere.

FPM is not installed in this validation environment, so this 65-family
checkpoint still requires external `fpm build` and `fpm test` confirmation.

Known limitations:

- Random sampling uses inverse transforms for consistency with the current
  SciFort RNG layer rather than specialized direct generators.
- The fitting layer remains a bounded local projected coordinate/pattern search
  and requires finite bounds.
- Weibull maximum and the R-distribution have parameter-dependent support
  boundaries; ordinary two-sided likelihood derivatives are nonregular when
  observations lie exactly on a moving boundary, and the score follows the
  project's existing boundary-NaN convention there.

## Double-gamma, asymmetric-Laplace, truncated-normal, and log-uniform validation

Validation date: 2026-09-21.

This checkpoint adds double-gamma, asymmetric Laplace, finite-bound truncated
normal, and log-uniform/reciprocal continuous families with distribution,
likelihood, analytic-score, bounded-fit, explicit-state RNG, `scifort_stats`,
and scalar C-ABI integration. `test/test_scipy_transform_four.f90` checks
SciPy 1.17.0 PDF/CDF/PPF and selected logarithmic-tail reference values,
analytic scores against centered finite differences, deterministic
inverse-transform RNG mapping, fixed-parameter fit dispatch, and C-ABI
equality.

A fresh exact-source GNU Fortran 14.2.0 strict build compiled all 89 library
sources in Fortran 2018 mode with `-fcheck=all`, floating-point zero/overflow
traps, backtraces, warnings, and warnings-as-errors, with the repository's
existing real-comparison and constant-integer-division warnings left nonfatal.
All 30 test programs passed. A separate `-O2` build of the same source also
passed all 30 test programs. A standalone C11 program compiled against
`include/scifort.h`, linked to the optimized Fortran objects, and successfully
called all four new C-ABI families.

`tools/generate_scipy_transform_reference.py` under SciPy 1.17.0 regenerates
`test/scipy_transform_reference.f90` byte-for-byte. The shipped application and
example compile and run against the optimized objects. ASCII, 132-column,
one-statement-per-line, and clean-source-tree checks pass.

FPM is not installed in this validation environment, so this 69-family
checkpoint still requires external `fpm build` and `fpm test` confirmation.

Known limitations:

- `truncnorm` currently accepts finite standardized `a` and `b` only; SciPy's
  infinite-endpoint cases remain future compatibility work.
- `loguniform` is overidentified if `a`, `b`, `loc`, and `scale` are all fitted
  simultaneously. Fix `scale=1` to reproduce SciPy's default fit convention.
- Truncated-normal and log-uniform support boundaries move with fitted
  parameters; ordinary two-sided scores are nonregular at exact boundaries and
  follow SciFort's boundary-NaN policy.
- Random sampling uses inverse transforms and the fitter remains a bounded local
  projected coordinate/pattern search requiring finite bounds.

## Folded normal, folded Cauchy, reciprocal inverse-Gaussian, and truncated Pareto validation

Validation date: 2026-09-21.

This checkpoint adds folded normal, folded Cauchy, reciprocal inverse-Gaussian,
and truncated Pareto continuous families with distribution, likelihood,
analytic-score, bounded-fit, explicit-state RNG, `scifort_stats`, and scalar
C-ABI integration. `test/test_folded_reciprocal_four.f90` checks SciPy 1.17.0
PDF, CDF, PPF, and selected logarithmic-tail reference values, likelihood
reductions, analytic scores against centered finite differences, deterministic
inverse-transform RNG mapping, fixed-parameter fit dispatch, and C-ABI
equality.

A fresh exact-source GNU Fortran 14.2.0 strict build compiled all 93 library
sources in Fortran 2018 mode with `-fcheck=all`, floating-point zero/overflow
traps, backtraces, warnings, and warnings-as-errors, with the repository's
existing real-comparison and constant-integer-division warnings left nonfatal.
All 31 test programs passed.

A separate `-O2` build of all 93 library sources compiled successfully, and
`test_folded_reciprocal_four` passed against those optimized objects. A
standalone C11 program compiled against `include/scifort.h`, linked to the
optimized Fortran objects, and successfully called all four new C-ABI
families. The attempted all-tests optimized batch exceeded the local execution
window before returning a result, so this checkpoint does not claim a complete
31/31 optimized run.

`tools/generate_folded_reciprocal_reference.py` under SciPy 1.17.0 regenerates
`test/folded_reciprocal_reference.f90` byte-for-byte. The shipped application
and example compile and run against the optimized objects. ASCII, 132-column,
one-statement-per-line, and clean-source-tree checks pass.

FPM is not installed in this validation environment, so this 73-family
checkpoint still requires external `fpm build` and `fpm test` confirmation.

Known limitations:

- Folded-normal and folded-Cauchy quantiles use safeguarded bisection rather
  than specialized rational approximations.
- Random sampling uses inverse transforms for consistency with the current
  SciFort RNG layer rather than specialized direct generators.
- The fitting layer remains a bounded local projected coordinate/pattern search
  and requires finite bounds.
- Truncated-Pareto support moves with `c`, `loc`, and `scale`; ordinary
  two-sided score derivatives are nonregular for observations exactly on a
  moving support boundary and follow SciFort's boundary-NaN convention.

## Exponentially modified normal, Johnson SB/SU, and trapezoid validation

Validation date: 2026-09-21.

This checkpoint adds exponentially modified normal, Johnson SB, Johnson SU,
and trapezoid continuous families with distribution, likelihood,
analytic-score, bounded-fit, explicit-state RNG, `scifort_stats`, and scalar
C-ABI integration. `test/test_normal_transform_four.f90` checks SciPy 1.17.0
PDF, CDF, PPF, and selected logarithmic-tail reference values, analytic scores
against centered finite differences, deterministic inverse-transform RNG
mapping, fixed-parameter fit dispatch, and C-ABI equality.

A fresh exact-source GNU Fortran 14.2.0 strict build compiled all 97 library
sources in Fortran 2018 mode with `-fcheck=all`, floating-point zero/overflow
traps, warnings, and warnings-as-errors, with the repository's established
real-comparison and constant-integer-division warnings left nonfatal. All 32
test programs passed.

`tools/generate_normal_transform_reference.py` under SciPy 1.17.0 regenerates
`test/normal_transform_reference.f90` byte-for-byte. Optimized validation for
this checkpoint covers the new four-family regression and standalone C ABI
smoke test rather than claiming a complete optimized all-tests run. The shipped
application and example also compile and run against the exact final source.
ASCII, 132-column, one-statement-per-line, and clean-source-tree checks pass.

FPM is not installed in this validation environment, so this 77-family
checkpoint still requires external `fpm build` and `fpm test` confirmation.

Known limitations:

- `exponnorm` quantiles use safeguarded bisection rather than a specialized
  closed-form or rational inverse.
- Trapezoid support and plateau boundaries move with fitted parameters;
  ordinary two-sided score derivatives are nonregular at exact boundaries and
  follow SciFort's boundary-NaN convention.
- Random sampling uses inverse transforms and the fitter remains a bounded
  local projected coordinate/pattern search requiring finite bounds.

## Burr, Mielke, Gibrat, and wrapped-Cauchy validation

Validation date: 2026-09-21.

This checkpoint adds Burr Type III, Mielke beta-kappa/Dagum, Gibrat, and
wrapped Cauchy continuous families with distribution, likelihood,
analytic-score, bounded-fit, explicit-state RNG, `scifort_stats`, and scalar
C-ABI integration. `test/test_burr_mielke_gibrat_wrapcauchy.f90` checks all
eight scalar distribution methods against generated SciPy 1.17.0 references,
CDF/PPF and SF/ISF tail round trips, support endpoints, invalid parameters, elemental
array behavior, analytic scores against centered finite differences,
deterministic inverse-transform RNG mapping and invalid-call state preservation,
fixed-parameter fit dispatch, and C-ABI equality. It also checks the Mielke/Burr
and Gibrat/lognormal mathematical identities.

GNU Fortran 14.2.0 compiled all 101 library sources in Fortran 2018 mode with
`-Wall -Wextra -Wimplicit-interface -fcheck=all -ffree-line-length-none`.
The complete 33-program test suite then passed under the same runtime-checking
configuration. A C11 header smoke program compiled with `-Wall -Wextra -Werror`,
called all 12 new scalar C entry points, linked, and ran successfully. The shipped
application and example also compile and run against the updated objects. The new
four-family test also passes with
`-ffpe-trap=invalid,zero,overflow`; five older IEEE-input tests intentionally
construct or propagate NaNs and therefore are not compatible with globally
trapping invalid floating-point operations.

`tools/generate_burr_mielke_gibrat_wrapcauchy_reference.py` under SciPy 1.17.0
regenerates `test/burr_mielke_gibrat_wrapcauchy_reference.f90` byte-for-byte.

FPM is not installed in this validation environment, so this 81-family
checkpoint still requires external `fpm build` and `fpm test` confirmation.

Known limitations:

- Random sampling uses inverse transforms for consistency with the current
  SciFort RNG layer rather than specialized direct generators.
- The fitting layer remains a bounded local projected coordinate/pattern search
  and requires finite bounds.
- Wrapped-Cauchy location and scale move both support endpoints; ordinary
  two-sided score derivatives are nonregular for observations exactly on a
  moving endpoint and follow SciFort's boundary-NaN convention.

## Generalized-gamma, half-generalized-normal, ARGUS, and Erlang validation

Validation date: 2026-09-21.

This checkpoint adds generalized-gamma, half-generalized-normal, ARGUS, and
Erlang continuous families with complete eight-function distribution APIs,
likelihood/NNLF, analytic scores, bounded/fixed fitting, explicit-state RNG,
`scifort_stats` exports, and scalar C-ABI integration.
`test/test_gengamma_halfgennorm_argus_erlang.f90` checks all eight scalar
methods against generated SciPy 1.17.0 references, positive- and negative-shape
generalized-gamma tail round trips, support endpoints, invalid parameters,
elemental array behavior, analytic scores against centered finite differences,
deterministic inverse-transform RNG mapping and invalid-call state preservation,
fixed-parameter fit dispatch, and C-ABI equality. It also checks generalized
gamma against gamma and Weibull, half-generalized-normal against exponential and
half-normal, and Erlang against gamma, including SciPy-compatible noninteger
shape evaluation.

A fresh exact-source GNU Fortran 14.2.0 build compiled all 109 library sources
in Fortran 2018 mode with `-Wall -Wextra -Wimplicit-interface -fcheck=all`.
The complete 35-program test suite passed under runtime checking. The new
four-family regression also passed when rebuilt with
`-ffpe-trap=invalid,zero,overflow`. A C11 smoke program compiled with
`-Wall -Wextra -Werror`, called all 12 new scalar C entry points, linked, and
ran successfully. The shipped application and example also compiled and ran
against the exact final static archive.

`tools/generate_gengamma_halfgennorm_argus_erlang_reference.py` under SciPy
1.17.0 regenerates `test/gengamma_halfgennorm_argus_erlang_reference.f90`
byte-for-byte. The generalized-gamma endpoint tests additionally cover the
zero, finite, and singular density limits at the lower support boundary when
`c>0`, plus the zero-density limit for `c<0`.

FPM and fprettify are not installed in this validation environment, so this
89-family checkpoint still requires external `fpm build`, `fpm test`, and any
repository-preferred fprettify confirmation.

Known limitations:

- ARGUS quantiles use safeguarded bisection on `[0,1]` rather than a specialized
  inverse approximation.
- Random sampling uses inverse transforms for consistency with the current
  SciFort RNG layer rather than specialized direct generators.
- The fitting layer remains a bounded local projected coordinate/pattern search
  and requires finite bounds.
- Erlang accepts every finite positive shape value numerically, matching the
  computation performed by SciPy, but SciFort does not reproduce SciPy's
  user-facing warning for noninteger shape values.

## Crystal Ball, Jones-Faddy skew-t, Pearson III, and relativistic Breit-Wigner validation

Validation date: 2026-09-21.

This checkpoint adds Crystal Ball, Jones-Faddy skew-t, Pearson III, and
relativistic Breit-Wigner continuous families with complete eight-function
distribution APIs, likelihood/NNLF, analytic scores, bounded/fixed fitting,
explicit-state RNG, `scifort_stats` exports, and scalar C-ABI integration.
`test/test_crystalball_jf_pearson_breitwigner.f90` checks all eight scalar
methods against generated SciPy 1.17.0 references, CDF/PPF and SF/ISF tail
round trips, invalid parameters and support endpoints, all new score components
against centered finite differences, deterministic inverse-transform RNG and
invalid-call state preservation, fixed-bound fit dispatch, and all 12 new C
entry points. Identity checks cover Jones-Faddy with equal shapes against
Student t, Pearson III at zero skew against the normal law, and the positive/
negative-skew Pearson III reflection relation.

A fresh exact-source GNU Fortran 14.2.0 build compiled all 113 library sources
in Fortran 2018 mode with `-Wall -Wextra -Wimplicit-interface -fcheck=all`.
The complete 36-program test suite passed under runtime checking. The new
four-family regression also passed with
`-ffpe-trap=invalid,zero,overflow`. A C11 smoke program compiled with
`-Wall -Wextra -Werror`, called all 12 new scalar C entry points, linked, and
ran successfully. The shipped application and distribution example both
compiled and ran against the exact final static archive.

`tools/generate_crystalball_jf_pearson_breitwigner_reference.py` under SciPy
1.17.0 regenerates `test/crystalball_jf_pearson_breitwigner_reference.f90`
byte-for-byte. A separate multi-parameter grid check covered Crystal Ball,
Jones-Faddy skew-t, positive/negative/near-zero Pearson III, and relativistic
Breit-Wigner values. PDF and CDF results agreed with SciPy at roughly binary64
roundoff over that grid. Quantile differences were small; for the extreme
relativistic Breit-Wigner upper tail, SciFort deliberately solves its direct
survival function rather than a rounded `1-CDF` value.

The Jones-Faddy beta-coordinate transform uses a scaled hypot-style evaluation
rather than forming `z*z`. This preserves finite logarithmic tails for very
large finite standardized observations; the regression explicitly exercises
`abs(z)=1e150`, where the naive squared transform overflows binary64.

The complex CDF antiderivative in `src/stats/rel_breitwigner.f90` is adapted
from SciPy 1.17.0 under BSD-3-Clause. The source notice, exact upstream path,
and adaptation details are recorded in `CODE_PROVENANCE.md`, and the complete
applicable license text is retained in `THIRD_PARTY_LICENSES.md`.

FPM and fprettify are not installed in this validation environment, so this
93-family checkpoint still requires external `fpm build`, `fpm test`, and any
repository-preferred fprettify confirmation.

Known limitations:

- Relativistic Breit-Wigner PPF and ISF use safeguarded bisection rather than a
  specialized analytic or rational inverse.
- Random sampling uses inverse transforms for consistency with the current
  SciFort RNG layer rather than specialized direct generators.
- The fitting layer remains a bounded local projected coordinate/pattern search
  and requires finite bounds.
- Pearson III intentionally uses SciPy's small-skew normal transition for
  numerical stability rather than evaluating an increasingly ill-conditioned
  gamma transform arbitrarily close to zero skew.

## Generalized exponential, skew-normal, Tukey lambda, and Rice validation

Validation date: 2026-09-21.

This checkpoint adds generalized exponential, skew-normal, Tukey lambda, and
Rice continuous families with complete eight-function distribution APIs,
likelihood/NNLF, analytic scores, bounded/fixed fitting, explicit-state RNG,
`scifort_stats` exports, and scalar C-ABI integration.
`test/test_genexpon_skewnorm_tukeylambda_rice.f90` checks all eight scalar
methods against generated SciPy 1.17.0 references, CDF/PPF and SF/ISF tail
round trips, support endpoints and invalid parameters, all new score components
against centered finite differences, deterministic inverse-transform RNG and
invalid-call state preservation, fixed-bound fit dispatch, and all 12 new C
entry points. Identity checks cover skew-normal at zero shape against normal,
Tukey lambda at zero against logistic and at one against uniform, and Rice at
zero shape against Rayleigh.

A fresh exact-source GNU Fortran 14.2.0 build compiled all 117 library sources
in Fortran 2018 mode with `-Wall -Wextra -Wimplicit-interface -fcheck=all`. The
complete 37-program test suite passed under runtime checking. The new regression
also passed with `-ffpe-trap=invalid,zero,overflow`. A C11 smoke program compiled
with `-Wall -Wextra -Werror`, called all 12 new scalar C entry points, linked,
and ran successfully.

`tools/generate_genexpon_skewnorm_tukeylambda_rice_reference.py` under SciPy
1.17.0 regenerates `test/genexpon_skewnorm_tukeylambda_rice_reference.f90`
byte-for-byte. Skew-normal validation includes a strongly suppressed left-tail
case, and Rice checks direct survival/inverse-survival paths rather than forming
upper tails as `1-CDF`. A separate multi-parameter spot grid covered positive
and negative skew-normal shapes, negative/zero/positive Tukey-lambda shapes,
widely different generalized-exponential rates, and Rice shapes from zero to
20. The largest relative discrepancy against SciPy on that grid was about
`2.7e-11`, occurring for a skew-normal CDF of order `1e-10`; ordinary values
were generally at binary64 roundoff.

The generalized-exponential quantile uses a dependency-free real principal
Lambert-W solver. Skew-normal CDF evaluation uses an independently implemented
Owen-T quadrature plus a transformed direct lower-tail integral where ordinary
subtraction would lose precision. Tukey-lambda inversion is based on its native
quantile representation. Rice tails use a Poisson mixture of direct regularized
incomplete-gamma tails, while its log-density and score use independently
implemented scaled modified-Bessel evaluations.

FPM and fprettify are not installed in this validation environment, so this
97-family checkpoint still requires external `fpm build`, `fpm test`, and any
repository-preferred fprettify confirmation.

Known limitations:

- Skew-normal quantiles use safeguarded bisection rather than a specialized
  rational approximation.
- Tukey-lambda CDF and density require inversion of the quantile function and
  therefore prioritize correctness over throughput.
- Rice CDF/SF use a Poisson mixture rather than a dedicated Marcum-Q/noncentral
  chi-square asymptotic algorithm, so very large shape parameters can be more
  expensive than SciPy's specialized implementation.
- Random sampling uses inverse transforms for consistency with the current
  SciFort RNG layer rather than specialized direct generators.
- The fitting layer remains a bounded local projected coordinate/pattern search
  and requires finite bounds.


## Double-Pareto lognormal, von Mises, and asymptotic Kolmogorov validation

Validation date: 2026-09-21.

This checkpoint adds double-Pareto lognormal, circular von Mises, von Mises on
a line, and the asymptotic two-sided Kolmogorov distribution with complete
eight-function APIs, likelihood/NNLF, analytic scores, bounded/fixed fitting,
explicit-state RNG, `scifort_stats` exports, and scalar C-ABI integration. The
Rice family's scaled modified-Bessel kernels were also promoted into
`scifort_special` as shared `i0e` and `i1e` functions and reused by von Mises.

`test/test_dpareto_vonmises_kstwobign.f90` checks all eight scalar methods
against generated SciPy 1.17.0 references, CDF/PPF and SF/ISF round trips,
support endpoints and invalid parameters, all new analytic score components
against centered finite differences, deterministic inverse-transform RNG and
invalid-call state preservation, fixed-bound fit dispatch, and all 12 new C
entry points. Structural checks cover reciprocal symmetry of a symmetric
log-scale double-Pareto lognormal law, equivalence of circular and line von
Mises inside the central period, circular CDF periodicity, the `kappa=0`
uniform limit, and direct Kolmogorov complements on both sides of the theta-
series switch.

`tools/generate_dpareto_vonmises_kstwobign_reference.py` requires SciPy 1.17.0
and regenerates `test/dpareto_vonmises_kstwobign_reference.f90`; reference
regeneration is expected to be byte-for-byte reproducible. The DPLN quantile
solver works on the logarithm of the positive standardized variate, von Mises
uses independently implemented adaptive quadrature rather than SciPy's CDF
kernel, and `kstwobign` switches between mathematically equivalent theta and
alternating Kolmogorov representations for tail stability.


A clean GNU Fortran 14.2.0 build compiled all 121 library sources with
Fortran 2018, warnings, implicit-interface diagnostics, and `-fcheck=all`; all
38 test programs passed. The new regression also passed with
`-ffpe-trap=invalid,zero,overflow`. A strict C11 smoke program compiled with
`-Wall -Wextra -Werror`, linked through the Fortran archive, and called all 12
new scalar ABI entry points. Both shipped demo programs compiled and ran.
Reference regeneration was byte-for-byte identical. FPM and fprettify are not
installed in this validation environment, so external `fpm build`, `fpm test`,
and repository-preferred formatting confirmation remain required.

A 151-value supplementary grid covered three substantially different DPLN
shape sets, von Mises concentrations through 20, both sides of the Kolmogorov
series switch, extreme quantiles, and shared `i0e`/`i1e` values. The largest
relative SciPy difference was about `1.8e-6` for a von Mises CDF of only
`3.9e-9` (about `7e-15` absolute). At `kstwobign x=0.83`, SciFort's direct
Kolmogorov density series differs from SciPy 1.17.0 by about `2.8e-8`; an
80-digit evaluation of the defining series agrees with SciFort, so this is
retained as a deliberate accuracy difference rather than fitted to SciPy's
approximation.

Known limitations:

- Circular von Mises CDF evaluation prioritizes direct quadrature accuracy and
  may differ slightly from SciPy's internal approximation for very large
  concentration parameters even when the underlying integral is more accurate.
- Von Mises and double-Pareto lognormal sampling currently use inverse
  transforms rather than specialized direct generators.
- The bounded fitting layer remains the project's local projected
  coordinate/pattern search and requires finite parameter bounds.

## Irwin-Hall, one-sided KS, noncentral chi-square, and noncentral F validation

Validation date: 2026-09-21.

This checkpoint adds Irwin-Hall, one-sided finite-sample Kolmogorov-Smirnov
(`ksone`), noncentral chi-square, and noncentral F distributions with complete
eight-function APIs, likelihood/NNLF, analytic continuous-parameter scores,
bounded/fixed fitting, explicit-state RNG, `scifort_stats` exports, and scalar
C-ABI integration. Irwin-Hall and `ksone` retain their sample-count/number-of-
uniforms shape as an integer-constrained fit parameter rather than inventing a
continuous derivative for an intrinsically discrete shape.

`test/test_irwinhall_ksone_ncx2_ncf.f90` checks all eight scalar methods against
generated SciPy 1.17.0 references, CDF/PPF and SF/ISF round trips, support and
invalid-parameter behavior, every exposed analytic score component against
centered finite differences, deterministic inverse-transform RNG and
invalid-call state preservation, fixed-bound fit dispatch, and all 12 new C
entry points. Structural checks cover the `n=1` uniform identities for
Irwin-Hall and `ksone` and the `nc=0` reductions of `ncx2` and `ncf` to central
chi-square and F distributions. Additional generated stress references cover a
`ksone` quantile at probability `1e-10`, an `ncx2(df=0.5)` quantile at `1e-9`,
an `ncx2` density of order `5e-72`, and a noncentral-F survival probability of
order `1e-9` at a variate near `3.8e8`.

`tools/generate_irwinhall_ksone_ncx2_ncf_reference.py` requires SciPy 1.17.0
and regenerates `test/irwinhall_ksone_ncx2_ncf_reference.f90` byte-for-byte.
Irwin-Hall uses a cardinal-B-spline recurrence rather than an alternating
power sum. `ksone` uses a logarithmically scaled exact Smirnov sum and a direct
small-`x` CDF branch. Noncentral chi-square/F use centered Poisson mixtures of
direct incomplete-gamma/beta lower and upper tails rather than forming one
tail by subtraction. Density mixtures keep Poisson weights in log space and
terminate on the scaled density contribution, which preserves far-tail terms
whose raw Poisson weight is tiny but whose conditional density is material.

A clean GNU Fortran 14.2.0 validation build compiled all 125 library sources
with Fortran 2018, warnings, implicit-interface diagnostics, and `-fcheck=all`;
all 39 test programs passed. The new regression also passed with
`-ffpe-trap=invalid,zero,overflow`. A strict C11 smoke program compiled with
`-Wall -Wextra -Werror`, linked through the Fortran archive, and called all 12
new scalar ABI entry points. Reference regeneration was byte-for-byte
identical.

Known limitations:

- The cardinal-B-spline Irwin-Hall implementation is correctness-first and can
  be slower than specialized asymptotic algorithms for very large `n`.
- The noncentral chi-square and noncentral F implementations use explicit
  Poisson mixtures rather than specialized large-noncentrality asymptotics, so
  extreme noncentrality can be more expensive than SciPy's specialized
  backends.
- `ksone` density derivatives are undefined at the finite-sample summation
  knots; the density itself uses a symmetric limiting average there.
- Random sampling uses inverse transforms for consistency with the current RNG
  layer rather than specialized direct generators.
- The bounded fitting layer remains the project's local projected
  coordinate/pattern search and requires finite parameter bounds.

## Randint, Planck, discrete Laplace, and logarithmic-series parity

`test/test_randint_planck_dlaplace_logser.f90` checks PMF, log-PMF, CDF, SF,
log-CDF, log-SF, PPF, and ISF against generated SciPy 1.17.0 references. It
also verifies discrete-uniform normalization and finite support endpoints, the
exact Planck/geometric identity, discrete-Laplace reflection symmetry, and the
logarithmic-series mass-plus-tail identity. Planck, discrete-Laplace, and
logarithmic-series analytic shape scores are checked against centered finite
differences; `randint` has no differentiable shape score because its two shape
parameters are integral bounds. The regression covers deterministic
inverse-transform RNG, invalid-call state preservation, fixed-bound fitting,
and all 12 new scalar C ABI entry points.

`tools/generate_randint_planck_dlaplace_logser_reference.py` requires SciPy
1.17.0 and regenerates
`test/randint_planck_dlaplace_logser_reference.f90` byte-for-byte. Additional
stress references exercise Planck and discrete-Laplace ISFs at `1e-12`, a
logarithmic-series ISF at `1e-10` with `p=0.99`, and a direct log-series upper
tail at count 1000. A supplementary shifted-lattice grid using `loc=0.5`
matched SciPy CDF/SF/PPF/ISF values with worst absolute difference about
`1.9e-15`.

A clean GNU Fortran 14.2.0 validation build compiled all 129 library sources
with Fortran 2018, warnings, implicit-interface diagnostics, and `-fcheck=all`;
all 40 test programs passed. The new regression also passed with
`-ffpe-trap=invalid,zero,overflow`. A strict C11 smoke program compiled with
`-Wall -Wextra -Werror`, linked through the Fortran archive, and called all 12
new scalar ABI entry points. Both shipped Fortran demos compiled and ran, and
reference regeneration was byte-for-byte identical. `fpm` and `fprettify` were
not installed in this validation environment.

## Beta-binomial, hypergeometric, negative-hypergeometric, and Boltzmann parity

`test/test_betabinom_hypergeom_nhypergeom_boltzmann.f90` checks PMF, log-PMF,
CDF, SF, log-CDF, log-SF, PPF, and ISF against generated SciPy 1.17.0
references. Structural checks cover the beta-binomial/Bernoulli reduction at
`n=1`, beta-binomial reflection when its beta shapes are equal, hypergeometric
population-complement symmetry, the standard negative-hypergeometric relation
to a hypergeometric draw, Boltzmann as a truncated Planck/geometric law, exact
finite-support normalization, endpoint conventions, and invalid parameters.
Analytic beta-binomial `a,b` scores and the Boltzmann rate score are checked
against centered finite differences; hypergeometric and negative-hypergeometric
have only integral shape parameters and therefore zero-length continuous score
vectors. The regression also covers deterministic inverse-transform RNG,
invalid-call state preservation, fixed-bound fitting, and all 12 new scalar C
ABI entry points.

`tools/generate_betabinom_hypergeom_nhypergeom_boltzmann_reference.py` requires
SciPy 1.17.0 and regenerates
`test/betabinom_hypergeom_nhypergeom_boltzmann_reference.f90` byte-for-byte.
Ordinary reference points come directly from SciPy. Two deliberately difficult
stress tails use direct high-precision defining formulas instead: for
`nhypergeom.sf(60,500,120,100)` the high-precision value is approximately
`5.972895432242913e-7`, while SciPy 1.17.0 returns
`5.972894242534110e-7`; SciFort returns approximately
`5.972895432245056e-7`. For `boltzmann.sf(99990,1e-8,100000)` the high-precision
value is approximately `8.995501154797533e-5`, whereas SciPy 1.17.0 returns
`8.995501147857077e-5`; SciFort agrees with the high-precision value to binary64
roundoff. These differences are retained rather than deliberately reproducing
SciPy's loss of precision.


A clean GNU Fortran 14.2.0 validation build for this slice compiled all 133
library sources with Fortran 2018, warnings, implicit-interface diagnostics,
and `-fcheck=all`; all 41 test programs passed. The new regression also passed
with `-ffpe-trap=invalid,zero,overflow`. A strict C11 smoke program compiled
with `-Wall -Wextra -Werror`, linked through the Fortran archive, and called all
12 new scalar ABI entry points. Both shipped Fortran demos compiled and ran,
and reference regeneration was byte-for-byte identical. `fpm` and `fprettify`
were not installed in this validation environment.


## Beta-negative-binomial, Yule-Simon, Zipf, and Zipfian slice

`test/test_betanbinom_yulesimon_zipf_zipfian.f90` compares PMF, log-PMF, CDF,
SF, log-CDF, log-SF, PPF, and ISF against values generated with SciPy 1.17.0.
It also checks Yule-Simon's `alpha=1` identity, the uniform `zipfian(a=0)`
limit, direct-tail complement identities, shifted lattices, invalid parameters,
analytic shape scores against centered finite differences, explicit-state RNG
behavior, fixed-bound fit dispatch, and all twelve scalar C ABI entry points.

The shared Hurwitz-zeta kernel is independently checked at `zeta(2,1)` and its
exponent derivative. `tools/generate_betanbinom_yulesimon_zipf_zipfian_reference.py`
regenerates the numerical reference module byte-for-byte with SciPy 1.17.0.

## Generalized inverse Gaussian, normal inverse Gaussian, and Skellam slice

Validation date: 2026-09-21.

`test/test_gig_nig_skellam.f90` checks all eight scalar distribution methods
against generated SciPy 1.17.0 references, complement identities, Skellam
normalization/reflection/shift behavior, invalid parameters, every exposed
analytic score component against centered finite differences, deterministic
inverse-transform RNG and invalid-call state preservation, fixed-bound fitting,
and all nine new scalar C ABI entry points.

The shared modified-Bessel-K kernel receives independent high-precision stress
checks at order 100 and argument 1000 for `log(K_nu)`, its order derivative,
and its argument derivative. Adaptive quadrature panel widths were added after
a wide-order stress sweep exposed a roughly `5.9e-5` error in the original
fixed-panel asymptotic/integral transition; the corrected value agrees with an
80-digit defining-integral evaluation to binary64 roundoff.

GIG and NIG tails are evaluated directly in transformed coordinates and in the
log domain. A wide-tail sweep exposed loss of endpoint resolution in the first
fixed-width NIG quadrature at `x=-200`; localizing monotone tails and adapting
panel width reduced deep-log-tail discrepancies to a few `1e-9` in absolute
log probability through the tested range. Independent high-precision stress
constants in the regression include GIG log-CDF/log-SF values near magnitude
`1150` and NIG log-CDF/log-SF values through `x=+/-300`.

Two difficult inverse-tail cases intentionally retain SciFort values rather
than imitating SciPy 1.17.0 approximation failures. For
`geninvgauss.isf(1e-6,p=-0.5,b=0.01)`, SciFort returns approximately
`1048.9614348032378`; direct high-precision integration at that point gives a
log survival matching `log(1e-6)` to about `9e-15`, while SciPy 1.17.0 returns
approximately `3083.614691`. For
`norminvgauss.isf(1e-6,a=1,b=0.99)`, SciFort returns approximately
`753.9034692822579`; direct high-precision integration matches the requested
log survival to about `5e-14`, while SciPy 1.17.0 is lower by roughly `0.005`.

The release validation for this slice compiled all 142 library sources with GNU
Fortran 14.2.0 in Fortran 2018 mode, warnings-as-errors, implicit-interface
diagnostics, and `-fcheck=all`; all 43 test programs passed. The new regression
also passed with `-ffpe-trap=invalid,zero,overflow`. A strict C11 smoke program
compiled with `-Wall -Wextra -Werror`, linked through the Fortran archive, and
called all nine new scalar ABI entry points. The shipped Fortran distribution
demo compiled and ran, and SciPy 1.17.0 reference regeneration was byte-for-byte
identical. The source tree passed ASCII, Fortran/public-header line-length,
new-slice FORD dummy-comment, and build/cache-artifact audits. `fpm` and
`fprettify` are not installed in this validation environment, so their
repository-preferred checks must still be run on a system where they are
available.

## Generalized hyperbolic and Fisher noncentral-hypergeometric slice

Validation date: 2026-09-21.

`test/test_genhyperbolic_fisher.f90` checks all eight scalar distribution
methods against generated SciPy 1.17.0 references, the `p=-0.5` generalized-
hyperbolic/normal-inverse-Gaussian identity, the `odds=1` Fisher/hypergeometric
identity, support shifts, normalization/complement identities, every exposed
analytic score against centered finite differences, explicit-state RNG and
invalid-call state preservation, fixed-bound fitting, and all six new scalar C
ABI entry points. The scaled-log Bessel-K routine is also checked against
`scipy.special.kve`.

A wider GH stress sweep used `p` from -2 through 2, skew parameters near both
open boundaries, and observations from -20 through 20. Away from the exact
closed `p<0, abs(b)=a` boundary (where SciPy 1.17.0 fails during CDF/quantile
integration), PDF/CDF differences were near binary64 roundoff; the largest
observed absolute log-tail difference was about `4.6e-10` at a log-survival
near -198.3. Fisher PMFs across several population/odds configurations agreed
with SciPy to about `1.5e-15` absolute error.

The release validation for this slice compiled all 144 library sources with
GNU Fortran 14.2.0 in Fortran 2018 mode with runtime/interface checks; all 44
test programs passed. The new regression also passed with
`-ffpe-trap=invalid,zero,overflow`. Reference regeneration was byte-for-byte
identical. `fpm` and `fprettify` were not installed in this validation
environment.


## Noncentral Student t and Gauss hypergeometric slice

Validation date: 2026-09-22.

`test/test_nct_gausshyper.f90` checks all eight scalar methods against generated
SciPy 1.17.0 references, the `nc=0` noncentral-t/Student-t identity, the `c=0`
and `z=0` Gauss-hypergeometric/beta identities, direct tail complements,
invalid parameters, every exposed analytic score component against centered
finite differences, deterministic inverse-transform RNG and invalid-call state
preservation, fixed-bound fitting, and all six new scalar C ABI entry points.

A 210-point supplementary noncentral-t sweep found that the initial common
quadrature spacing was insufficient for extremely small degrees of freedom.
The released implementation uses progressively finer double-exponential grids
for `df<1` and `df<0.3`. After that refinement, the largest CDF difference from
SciPy in the sweep was about `1.6e-14`; the largest relative PDF difference was
about `2.4e-10` only for a density of order `1e-10` (absolute error about
`6e-20`). An independent defining-integral regression at `df=0.2`, `nc=6`,
`x=10` checks CDF `0.2504117582414005` and PDF `0.01449116494956138`.

The regression also retains independent direct-integral tail constants where
SciPy 1.17.0 loses accuracy. In particular, `nct_logsf(10,30,-5)` is checked
against `-57.214203664073054`. Gauss-hypergeometric stress cases include a
log-survival near `-63.30`, a concentrated-parameter log-survival near `-55.37`,
and a hard log-density near `-96.02`; SciFort evaluates normalization and tails
through a beta-quantile transform instead of integrating endpoint-singular beta
factors directly.

The release validation for this slice compiled all 146 library sources with GNU
Fortran 14.2.0 in Fortran 2018 mode with warnings, implicit-interface checks, and
`-fcheck=all`; all 45 test programs passed. The project-wide run was split into
two execution batches to stay within the harness timeout. `fpm` and `fprettify`
were not installed in this validation environment.


## Landau slice

Validation date: 2026-09-22.

`test/test_landau.f90` checks all eight scalar methods against generated SciPy
1.17.0 references from `x=-5` through `x=100`, quantiles from `1e-5` through
`1-1e-5`, direct CDF/SF complement identities, loc/scale transformation, deep
log-tail inverse round trips (`1e-100` lower tail and `1e-12` upper tail),
invalid inputs, the analytic loc/scale likelihood score against centered finite
differences, deterministic inverse-transform RNG with invalid-call state
preservation, fixed-bound fitting, and all three scalar C ABI entry points.

The released implementation uses the classical left-tail asymptotic expansion
for the double-exponential tail, an exponentially damped characteristic-function
quadrature centrally, and a scaled nonoscillatory Landau integral for the heavy
right tail. The complete checked build compiled all 147 library sources and all
46 Fortran test programs passed.


## Wallenius and Poisson-binomial slice

Validation date: 2026-09-22.

`test/test_wallenius_poisson_binom.f90` checks all eight scalar methods against
generated SciPy 1.17.0 references, Wallenius `odds=1` reduction to ordinary
hypergeometric, Poisson-binomial reduction to binomial for equal probabilities,
normalization and shifted supports, every exposed analytic score against
centered finite differences, exact explicit-state RNG replay, invalid-call
state preservation, fixed-bound fitting, and all six new C ABI calls.

Wallenius probabilities are propagated by the exact finite-state sequential
biased-urn recursion, with the odds sensitivity propagated alongside the state
probabilities. Poisson-binomial masses and direct tails use Bernoulli-convolution
dynamic programming; its vector score uses leave-one-out convolutions. Ordinary
reference cases agree with SciPy 1.17.0 at essentially binary64 roundoff.

The checked integrated build compiled all 149 library sources with GNU Fortran
14.2.0 in Fortran 2018 mode with warnings, implicit-interface diagnostics, and
`-fcheck=all`; all 47 Fortran test programs passed. The new regression also
passed with `-ffpe-trap=invalid,zero,overflow`. A strict C11 smoke program built
with `-Wall -Wextra -Werror`, linked through the Fortran archive, and called all
six new ABI functions, including the pointer-plus-`size_t` Poisson-binomial
interface. SciPy 1.17.0 reference regeneration was byte-for-byte identical.
`fpm` and `fprettify` were not installed in this validation environment.


## Finite-sample two-sided Kolmogorov-Smirnov slice

Validation date: 2026-09-22.

`test/test_kstwo.f90` checks all eight methods against generated SciPy 1.17.0
references, the exact `n=1` distribution, the `2*ksone` upper-tail identity,
support endpoints and loc/scale transformation, a `1e-8` direct upper-tail
ISF/log-SF round trip, loc/scale score finite differences, explicit-state RNG
with invalid-call state preservation, fixed-bound fitting, and all three new C
ABI calls.

A supplementary 288-point sweep covered sample sizes from 1 through 100000.
The largest observed CDF difference from SciPy was about `4.8e-15` through
`n=140`; in the large-`n` Pelz-Good regime the largest SF difference was below
`9.4e-13`. The largest PDF difference was about `4.1e-9` absolute, reflecting
the finite-difference density evaluation. Direct SciFort tail inversion is more
self-consistent than SciPy 1.17.0's approximate inverse in some very small
upper tails: for `n=200`, SciPy's nominal `isf(1e-8)` has SF about
`1.0114351e-8`, whereas SciFort's direct-SF bisection returns an argument whose
SF is `1e-8` to binary64 rounding.

The exact release-source build compiled all 150 library sources with GNU Fortran
14.2.0 in Fortran 2018 mode with warnings, implicit-interface diagnostics, and
`-fcheck=all`; all 48 Fortran test programs passed in two batches. The new test
also passed after rebuilding the full library with
`-ffpe-trap=invalid,zero,overflow`. A strict C11 smoke program compiled with
`-Wall -Wextra -Werror`, linked through the Fortran archive, and called all
three new ABI functions. Both shipped Fortran demos compiled and ran, and SciPy
1.17.0 reference regeneration was byte-for-byte identical. `fpm` and
`fprettify` were not installed in this validation environment.


## Levy-stable and studentized-range slice

Validation date: 2026-09-22.

`test/test_levy_stable_studentized_range.f90` checks all eight scalar methods
against generated SciPy 1.17.0 references, Levy-stable reductions to normal,
Cauchy, and Levy laws, S1 `alpha=1` scale shifting and symmetry, studentized
range support and asymptotic behavior, loc/scale transformations, direct-tail
quantile round trips, numerical likelihood scores, deterministic inverse-CDF
RNG, fixed-bound fits, and all six scalar C ABI entry points.

On a representative Levy-stable grid away from SciPy's deliberately rounded
`alpha`-near-one saturation cases, the largest observed relative PDF difference
was about `5e-14` and the largest CDF absolute difference about `3e-15`. For
`alpha=1.4, beta=0.35`, PPF/ISF comparisons over probabilities from `1e-4` to
`0.9999` were within about `1.5e-10` absolute, with central quantiles near
binary64 roundoff. Extremely close to the one-sided support boundary for
`alpha<1, |beta|=1`, the angular quadrature has a practical probability floor;
subnormal-scale SciPy probabilities there should not be expected to agree
relatively.

For studentized range, representative finite- and infinite-df CDF values agreed
with SciPy to about `2e-12` absolute. PDF agreement was typically much tighter;
the largest observed relative discrepancy in the initial grid was about
`3.3e-5` for a very small far-tail density. For `k=5, df=10`, PPF/ISF values
from `1e-4` to `0.9999` agreed to about `1.4e-11` absolute or better.

The complete listed source graph contains 152 library sources and compiles with
GNU Fortran 14.2.0 in Fortran 2018 mode using `-Wall -Wextra`, implicit-interface
diagnostics promoted to errors, and `-fcheck=all`. The new integrated regression
passes under that checked build and also with
`-ffpe-trap=invalid,zero,overflow`; the previous `test_kstwo` regression still
passes. `fpm` and `fprettify` were not installed in this validation
environment.

## Multivariate-distribution slice

Validation date: 2026-09-22.

`test/test_multivariate.f90` checks the dense Cholesky factor and triangular
solve by reconstruction, the symmetric eigensolver by spectral reconstruction,
and numerical-rank handling on a singular PSD matrix. Multivariate-normal
PDF/log-PDF and entropy are checked against generated SciPy 1.17.0 references
in two and three dimensions, including SciPy's lower-triangle covariance
semantics and a rank-one covariance both on and off its affine support. MLE
fitting is compared with `scipy.stats.multivariate_normal.fit`, including a
fixed-mean case. Normal CDF/log-CDF validation includes ordinary lower-tail
probabilities, finite lower/upper boxes, singular covariance matrices, and
two-, three-, and four-dimensional cases; the four-dimensional cases exercise
the pivoted conditional factorization path. Marginal extraction is checked in
nontrivial requested index order.

Multivariate Student-t PDF/log-PDF and entropy are compared directly with SciPy
1.17.0, including singular shape semantics and one-dimensional reductions to the
univariate Student t. CDF tests cover ordinary lower tails, finite boxes, and a
four-dimensional pivot-path case. Marginals and deterministic explicit-state RNG
replay are also checked. Representative normal CDF discrepancies are around
`1e-6` or smaller at the configured test budgets; representative multivariate-t
CDF discrepancies are within a few `1e-6`, with test tolerances set to `2e-5`
for the randomized low-discrepancy integration path.

Dirichlet PDF/log-PDF, mean, variance, covariance, entropy, and the `K-1`
coordinate convention are compared with SciPy. Multinomial PMF/log-PMF, mean,
covariance, entropy, binomial reduction, and SciPy 1.17's final-probability
adjustment are checked. Dirichlet-multinomial PMF/log-PMF and moments are also
checked. Explicit-state random variates are tested for deterministic replay and
support invariants.

The complete source graph contains 159 library sources. All 50 Fortran test
programs passed with GNU Fortran 14.2.0 after this slice was integrated, in
separate execution batches to keep long numerical regressions within the harness
time limit. A checked rebuild uses Fortran 2018, `-Wall -Wextra`,
implicit-interface diagnostics promoted to errors, and `-fcheck=all`. The new
multivariate regression also passes with
`-ffpe-trap=invalid,zero,overflow`. `fpm` and `fprettify` are not installed in
this validation environment.

## Multivariate hypergeometric and matrix-distribution slice

Validation date: 2026-09-22.

`test/test_matrix_distributions.f90` compares `multivariate_hypergeom`,
`normal_inverse_gamma`, `matrix_normal`, `wishart`, `invwishart`, and `matrix_t`
against generated SciPy 1.17.0 numerical references. The checks include PMF/PDF
and log forms, means, variances/covariances, matrix moments and entropy, the
matrix-t example values published with SciPy 1.17, finite-support behavior,
invalid-parameter results, SciPy-style lower-triangle SPD semantics, and the
zero-/one-object multivariate-hypergeometric edge cases. Every newly exposed
random-variate routine is checked for deterministic replay from an explicit
SciFort RNG state and for its structural support invariants.

The matrix-normal, Wishart, inverse-Wishart, and matrix-t density reference cases
agree with SciPy 1.17.0 to binary64 roundoff. Wishart sampling uses the Bartlett
decomposition. Inverse-Wishart sampling uses the mathematically equivalent
Wishart-inversion construction rather than SciPy's direct Axen sampler, and the
matrix-t sampler uses its inverse-Wishart/matrix-normal mixture, so RNG validation
checks deterministic SciFort replay and distribution support rather than exact
SciPy sample values.

The integrated graph contains 166 Fortran library sources. All 51 Fortran test
programs passed with GNU Fortran 14.2.0 in Fortran 2018 mode after this slice was
integrated. A checked build used `-Wall -Wextra`, implicit-interface diagnostics
promoted to errors, and `-fcheck=all`; the warnings are the project's existing
intentional exact-real comparisons and one compile-time integer-division note.
The new matrix-distribution regression also passes with
`-ffpe-trap=invalid,zero,overflow`. The SciPy 1.17.0 reference module regenerates
byte-for-byte identically, and the repository passes the ASCII-only source/text
check. `fpm` and `fprettify` are not installed in this validation environment.

## Directional and random-matrix distribution slice

Validation date: 2026-09-22.

`test/test_directional_random_matrix.f90` validates `uniform_direction`,
`vonmises_fisher`, `ortho_group`, `special_ortho_group`, `unitary_group`,
`random_correlation`, and `random_table`, together with the real-order scaled
modified-Bessel-I kernels used by von Mises-Fisher. PDF/log-PDF, entropy, fitted
direction/concentration, and fixed-margin table mass/mean values are compared with
generated SciPy 1.17.0 references. The Bessel regression exercises both the
positive-term series and the large-argument asymptotic branch; at
`nu=50, x=1000` the scaled-log result is within `3e-11` and the adjacent-order
ratio within `8e-13` of `scipy.special.ive` references.

Every RNG API is checked for deterministic replay from explicit SciFort state. The
2D, 3D, and higher-dimensional von Mises-Fisher sampling branches are exercised
and checked for unit norm. Orthogonal/special-orthogonal/unitary draws are checked
for their defining Gram identities and determinants where applicable. Random
correlation draws are checked for symmetry, exact unit diagonal, prescribed
eigenvalues including a rank-deficient spectrum, and invalid-trace rejection.
Random tables are checked for exact row/column margins, nonnegativity, support
semantics, and SciPy-compatible zero-total mean behavior.

After this slice, the integrated graph contains 175 Fortran library sources and all
52 Fortran test programs pass with GNU Fortran 14.2.0 in Fortran 2018 mode. The
full library compiles with `-Wall -Wextra`, `-Wimplicit-interface`,
`-Werror=implicit-interface`, and `-fcheck=all`; the new regression additionally
passes with `-ffpe-trap=invalid,zero,overflow`. The generated SciPy 1.17.0 reference
module regenerates byte-for-byte identically. `fpm` and `fprettify` are not
installed in this validation environment.


## Gaussian KDE and quasi-Monte Carlo

Validation date: 2026-09-22.

`test/test_gaussian_kde.f90` compares unweighted and weighted covariance,
effective sample size, Scott/Silverman/constant bandwidths, PDF/log-PDF values,
Gaussian overlap, KDE-product overlap, one-dimensional box integration, inverse
covariance, and marginal densities against generated SciPy 1.17.0 references.
A two-dimensional finite-box integral is checked at a `3e-5` absolute tolerance
because it uses randomized multivariate-normal integration. Resampling is checked
for exact replay under equal SciFort RNG states. Singular data, invalid bandwidths,
and one-observation input are covered explicitly.

`test/test_qmc.f90` compares all four discrepancy formulas (CD, WD, MD, and
L2-star), iterative centered discrepancy and `update_discrepancy`, Euclidean
minimum-distance/MST and city-block geometric criteria, scaling, Van der Corput,
and the unscrambled Halton and Sobol sequences directly with SciPy 1.17.0 numerical
outputs. Sobol coverage includes exact first/base-two blocks, reset/fast-forward,
a 7D sequence after skipping 12,345 points, 64-bit direction numbers at dimension
21,201, cumulative base-two balance rules, and bit exhaustion. Scrambled Halton and
Sobol paths use SciFort's explicit RNG and are tested for deterministic replay and
unit-cube invariants. Latin-hypercube tests cover ordinary stratification and every
pair of coarse symbols in the strength-two OA construction.

The same QMC regression compares deterministic `MultinomialQMC` and
`MultivariateNormalQMC` outputs with SciPy 1.17.0, including a rank-one singular
covariance. The Box-Muller normal path is checked for exact SciFort-state replay and
finite output. Poisson-disk tests cover replay, minimum separation, unit and negative
rectangular bounds, surface proposals, fill-space exhaustion, and invalid parameter
paths. SciFort uses exact accepted-set neighbor checks, so the minimum-radius
invariant does not depend on SciPy's grid indexing implementation.

After this slice, the complete source-order graph contains 178 Fortran library
sources and all 54 Fortran test programs pass with GNU Fortran 14.2.0 in Fortran
2018 mode. A clean checked build uses `-Wall -Wextra`, `-Wimplicit-interface`,
`-Werror=implicit-interface`, and `-fcheck=all`; the QMC regression also passes with
`-ffpe-trap=invalid,zero,overflow`. The KDE and QMC reference modules regenerate
byte-for-byte identically. `fpm` and `fprettify` are not installed in this validation
environment.

## Extended hypothesis and contingency tests

The Fisher lower-tail case for `[[6, 2], [1, 4]]` also has an exact rational
reference: its hypergeometric masses for counts 2 through 7 are
`[7, 70, 175, 140, 35, 2] / 429`, so the lower tail at 6 is `427/429`.
These were independently checked with integer binomial coefficients and Python
standard-library `fractions.Fraction`, using the distribution specified in the
[SciPy Fisher exact documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.fisher_exact.html).
The near-one lower-tail assertion allows `32 * epsilon(1.0_dp)` absolute error
for rounding in log-gamma evaluations and the five-term logarithmic sum.
The smaller Fisher p-values retain their existing tolerances.

Portability recheck on 2026-09-23: GNU Fortran (GCC) 17.0.0 20260510
(experimental), Windows. The full `fpm test` suite passes. The affected test
also passes with runtime checks:
`fpm test --target test_hypothesis_extended --profile debug --flag "-fcheck=all -fbacktrace"`.

Validation date: 2026-09-22. GNU Fortran 14.2.0 on Linux.

`test/test_hypothesis_extended.f90` validates the public `scifort_stats` exports
against SciPy 1.17.0 reference values for Wilcoxon rank-sum, Kruskal-Wallis,
Friedman, signed-rank Wilcoxon, Cressie-Read power divergence and Pearson
chi-square, Fisher exact, contingency-table chi-square with and without Yates'
correction, classical and Welch one-way ANOVA, Bartlett, Levene/Brown-Forsythe,
and Fligner-Killeen tests. Signed-rank coverage includes exact no-tie cases,
tied/zero observations, Pratt and zero-split conventions, one-sided alternatives,
and asymptotic z statistics. Invalid all-identical Kruskal data and invalid named
power-divergence choices are also exercised.

After this slice, the complete source-order graph contains 179 Fortran library
sources and all 55 Fortran test programs pass. A clean library build succeeds
with Fortran 2018, `-Wall -Wextra`, `-Wimplicit-interface`,
`-Werror=implicit-interface`, and `-fcheck=all`; the extended-hypothesis regression
also passes with `-ffpe-trap=invalid,zero,overflow`. The warnings in the new module
are limited to deliberate exact real comparisons for zero/tie/special-case
semantics and gfortran allocatable-assignment diagnostics. `fpm` and `fprettify`
are not installed in this validation environment.

## Goodness-of-fit and normality tests

Validation date: 2026-09-22. GNU Fortran 14.2.0 on Linux.

`test/test_goodness_of_fit.f90` validates public `scifort_stats` exports against
values generated by SciPy 1.17.0 for D'Agostino skew and kurtosis tests, the
D'Agostino-Pearson omnibus test, Jarque-Bera, exact/asymptotic one-sample KS,
exact one- and two-sided two-sample KS, finite-sample and asymptotic
Cramer-von Mises tests, Shapiro-Wilk, and Anderson-Darling tests for normal,
exponential, logistic, left-Gumbel, and right-Gumbel models. The regression also
checks KS extremum locations/signs and invalid/constant-data behavior.

After this slice the source-order graph contains 180 Fortran library sources and
all 56 Fortran test programs pass. The library builds with Fortran 2018,
`-Wall -Wextra`, `-Wimplicit-interface`, `-Werror=implicit-interface`, and
`-fcheck=all`. The goodness-of-fit regression also passes with
`-ffpe-trap=invalid,zero,overflow`. The generated SciPy 1.17.0 reference module
regenerates byte-for-byte identically.

## Additional nonparametric and multi-sample tests

Validation date: 2026-09-22. GNU Fortran 14.2.0 on Linux.

`test/test_nonparametric_extended.f90` validates public `scifort_stats` exports
against SciPy 1.17.0 reference values for the Ansari-Bradley and Mood scale
tests, Epps-Singleton two-sample test, k-sample Anderson-Darling, and Mood median
test. Ansari coverage includes exact untied probabilities, tied asymptotic
probabilities, and all three alternatives. Mood covers tied and untied cases
for all alternatives. Epps-Singleton is checked with the default and a custom
three-frequency grid. Anderson k-sample coverage includes midrank, right-side,
and continuous variants, tied samples, standard critical values, and capped or
interpolated p-values. Median-test coverage includes below/above/ignore tie
handling, contingency tables, and the log-likelihood divergence.

After this slice the complete source-order graph contains 181 Fortran library
sources and all 57 Fortran test programs pass. The complete library builds in
Fortran 2018 mode with `-Wall -Wextra`, `-Wimplicit-interface`,
`-Werror=implicit-interface`, and `-fcheck=all`; the new nonparametric regression
also passes with `-ffpe-trap=invalid,zero,overflow`. Its SciPy 1.17.0 reference
module regenerates byte-for-byte identically. Exact real comparisons in the
new module are intentional for tie grouping and duplicate detection.

## Association, robust regression, and Page trend tests

Validation date: 2026-09-22. GNU Fortran 14.2.0 on Linux.

`test/test_association_extended.f90` validates the public `scifort_stats`
exports against SciPy 1.17.0 generated references. Coverage includes tied
Kendall tau-b/tau-c, exact and asymptotic untied Kendall probabilities and all
three alternatives, point-biserial correlation, all `linregress` alternatives,
Theil-Sen/Siegel repeated-x cases and intercept conventions, Brunner-Munzel
with Student-t and normal references, and exact/asymptotic Page trend tests with
pre-ranked/reordered conditions.

After this slice the source-order graph contains 182 Fortran library sources and
all 58 Fortran test programs pass. The library builds in Fortran 2018 mode with
`-Wall -Wextra`, `-Wimplicit-interface`, `-Werror=implicit-interface`, and
`-fcheck=all`. The association regression also passes with
`-ffpe-trap=invalid,zero,overflow`, and its SciPy 1.17.0 reference module
regenerates byte-for-byte identically.

## Exact contingency, effect-measure, and p-value utility tests

Validation date: 2026-09-22. GNU Fortran 14.2.0 on Linux.

`test/test_contingency_meta.f90` validates the public `scifort_stats` exports
against generated SciPy 1.17.0 references. Coverage includes all three exact
binomial alternatives; exact, Wilson, and continuity-corrected Wilson intervals;
pooled/unpooled Barnard tests; all Boschloo alternatives; conditional and sample
odds ratios and confidence intervals; relative risk; Cramer/Tschuprow/Pearson
association; 2D margins; all five p-value-combination methods; and BH/BY false-
discovery control. Edge coverage includes degenerate binomial probabilities,
zero-column unconditional tests, odds-ratio support endpoints, zero risks, and
infinite/NaN p-value-combination endpoint semantics.

After this slice the source-order graph contains 183 Fortran library sources and
all 59 Fortran test programs pass. The library builds in Fortran 2018 mode with
`-Wall -Wextra`, `-Wimplicit-interface`, `-Werror=implicit-interface`, and
`-fcheck=all`. The new regression also passes with
`-ffpe-trap=invalid,zero,overflow`. The SciPy 1.17.0 reference module regenerates
byte-for-byte identically. Barnard/Boschloo use SciFort's deterministic global-
grid plus golden-section nuisance maximizer rather than SciPy's SHGO machinery.


## Multiple-comparison and Poisson-means tests

Validation date: 2026-09-22. GNU Fortran 14.2.0 on Linux.

`test/test_multiple_comparisons.f90` validates `alexandergovern`, pooled Tukey
HSD, unequal-size Tukey-Kramer, heteroscedastic Games-Howell, all three Dunnett
alternatives with simultaneous confidence intervals, and the Poisson means
E-test with zero and nonzero null differences against generated SciPy 1.17.0
references. Tukey/Games-Howell values agree to the numerical accuracy of the
studentized-range implementation. Dunnett uses randomized-QMC-style
multivariate-t integration in SciPy and deterministic SciFort integration, so
its p-values and interval endpoints are checked with integration-appropriate
tolerances rather than bitwise equality.

After this slice the source-order graph contains 184 Fortran library sources and
all 60 Fortran test programs pass. The library builds in Fortran 2018 mode with
`-Wall -Wextra`, `-Wimplicit-interface`, `-Werror=implicit-interface`, and
`-fcheck=all`. The new regression also passes with
`-ffpe-trap=invalid,zero,overflow`, and its SciPy 1.17.0 reference module
regenerates byte-for-byte identically.

## Merge review numerical regressions

`test/test_merge_regressions.f90` was introduced before the implementation fixes.
The initial run failed 50 checks, covering Zipf NaN/infinity behavior, large
positive/negative observations in hypergeometric and discrete-uniform tails,
wide discrete-uniform quantiles, log-softmax normalization and translation
invariance, and scalar/independent multivariate-normal log-tail identities.
Additional cases cover signed/empty normal intervals, nonunit diagonal covariance
with means and ignored upper-triangle entries, singular deterministic coordinates,
invalid inputs, and log probabilities outside the finite binary64 range.
The latter case caught a negative-infinity subtraction in the first version of
the interval helper; it now returns negative infinity rather than NaN.

Correlated multivariate-normal log-CDF calculations continue to use the existing
probability-domain QMC integrator; the new exact log-domain path covers scalar
and diagonal covariance only. No general correlated extreme-tail accuracy claim
is made by these regressions.

Validation on 2026-09-23, Windows, GNU Fortran (GCC) 17.0.0 20260510
(experimental): both commands completed successfully:

```text
fpm test
fpm test --target test_merge_regressions --profile debug --flag "-fcheck=all -fbacktrace"
```

The full suite includes the expanded regression program. ASCII checks and
`git diff --check` also pass. The ASCII checker exempts only the explicitly
listed verbatim CC BY license text; Fortran sources remain subject to the check.
