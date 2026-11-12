## General piecewise extension and degree-selected Hermite interpolation (2026-10-31)

**GNU Fortran F2018**, with both checked (`-O0 -g -fcheck=all -fbacktrace`)
and optimized (`-O3`) compilation of the five-module public interpolation
closure:

- 9/9 native Fortran interpolation regression programs passed per build,
  including a new `test_interpolate_piecewise_extend_orders.f90`.
- The new `compare_piecewise_extend_scipy.py` passed **480 default-seed
  scenarios per build** (64,476 values per build) plus **120 second-seed
  scenarios per build**. It covers real and complex PPoly/BPoly, appended and
  prepended intervals, descending knots, degree changes and multichannel
  outputs. It uses explicit mathematical assembly for prepending mixed-degree
  power polynomials where SciPy 1.17 has an implementation defect.
- The new `compare_bernstein_orders_scipy.py` passed **320 scenarios for each
  of two seeds per build**, plus additional third-seed cases. Each scenario
  was compared against an independent Hermite linear-system solve.
  Many asymmetric cases fail or disagree with SciPy 1.17's public
  `BPoly.from_derivatives` because of its endpoint-count selection defect;
  that is reported by the Python test, not hidden or counted as direct SciPy
  agreement. For the default and second seed, 373/640 cases agreed directly
  with SciPy per build; all 640 agreed with the independent mathematical oracle.
- Selected existing PPoly/BPoly, Hermite, complex piecewise, PPoly roots,
  and B-spline conversion regression scripts also passed. The entire FPM
  test inventory was not run because FPM was unavailable.

Commands:

```bash
python tools/interpolate_reference/compare_piecewise_extend_scipy.py
python tools/interpolate_reference/compare_piecewise_extend_scipy.py --seed 8675309 --cases 120
python tools/interpolate_reference/compare_bernstein_orders_scipy.py
python tools/interpolate_reference/compare_bernstein_orders_scipy.py --seed 8675309
```

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

## Finite-difference nonlinear constraint Hessians (2026-10-08)

`nonlinear_constraint%hess_fd_scheme` now selects a two- or three-point dense
approximation of the weighted nonlinear-constraint Hessian using an analytic
`%jac`. The two-point stencil allows signed relative steps, while three-point
differences use centered or bound-aware one-sided stencils. The numerical
derivative calls are charged to the selected block's `constr_njev`, while
`constr_nhev` remains zero without a user-supplied Hessian callback.

`tools/generate_optimize_trust_constr_constraint_hess_fd_reference.py`
reproduces a two-variable nonlinear equality from SciPy 1.17.0 using
`NonlinearConstraint(..., hess='2-point')` and `hess='3-point'`. Both SciPy
versions converge to approximately `[0.9218519654, 0.1507593917]`; their
constraint Jacobian evaluation counts are respectively 33 and 60 (SciFort
does not claim identical counts because it assembles a dense Hessian rather
than SciPy's Hessian LinearOperator). The generated Fortran constants and
`test/test_optimize_trust_constr_constraint_hess_fd.f90` check the optimum,
constraint feasibility, per-block accounting, signed steps, bounds, invalid
Jacobian/Hessian combinations, and generic `minimize` dispatch. A second
reference is a lower-inequality constraint with retained box bounds and a
BFGS objective Hessian. SciPy converges near `[0.8, 0.36831]`; the Fortran
test also exercises the interior-point method with the requested numerical
constraint Hessian.

This code was compiled and run in isolated GNU Fortran 14.2 builds with
temporary interface stubs for unrelated L-BFGS-B/SLSQP dependencies. Those
test-only stubs are not part of the distribution. Full `fpm test` coverage
requires installing the declared external dependencies.


## Mixed trust-constr constraint-Hessian validation (2026-10-08)

The `test_optimize_trust_constr_mixed_hessian.f90` regression covers blocks
with an exact weighted Hessian product alongside blocks with no Hessian
callback, both with analytic and numerical three-point Jacobians. It also
checks objective Hessian products, the exact-Hessian override,
independent per-block evaluation counts,
and second-order one-sided Hessian probing at a retained variable bound.
The constrained three-variable quadratic reference minimizer, generated by
SciPy 1.17.0 `minimize(method='trust-constr')`, is approximately
`[1.29147583, -0.46790981, 0.48106041]`.

With external FPM dependencies unavailable in the current sandbox, the
`scifort_optimize_minimize` source was compiled with GNU Fortran using
**test-only, never-invoked placeholders** for the L-BFGS-B and SLSQP
interfaces (no such placeholders are included in this repository). The
new test and all eight pre-existing `test_optimize_trust_constr*.f90`
programs passed with `-O0 -g -fcheck=all`. This checks the relevant
trust-constr code but is **not** a complete `fpm test` run.

## Scalar Chandrupatla minimum refinement (2026-10-09)

`scifort_optimize_minimize_scalar` adds scalar `find_minimum` with SciPy
1.17-style Chandrupatla three-point refinement and an optional stopping
callback. `tools/generate_optimize_find_minimum_reference.py` produces 15
SciPy reference cases for quadratic, quartic, kinked, exponential, and flat
objectives; reversed input triples, custom tolerances, the iteration limit,
and invalid bracketing. `test/test_optimize_find_minimum.f90` also verifies
callback termination, finite-difference-free work counts, invalid settings
before objective evaluation, nonfinite values both initially and during
iteration, and direct use of `bracket_minimum` output.

Six scalar optimization/root-finding regression programs pass in both GNU
Fortran 14.2 debug and optimized builds with `-std=f2018`, bounds checks,
IEEE invalid/division-by-zero/overflow traps, and
`-Werror=implicit-interface`. Tests compile the public `scifort_optimize`
exports using temporary, test-only declarations of unrelated external
solver interfaces. Those declarations are not shipped. Full FPM testing
still requires its declared external dependencies and an installed FPM.

## Rank-one elementwise optimization (2026-10-09)

`test/test_optimize_elementwise_arrays.f90` tests the rank-one overloads of
`bracket_root`, `find_root`, `bracket_minimum`, and `find_minimum` against their
scalar counterparts and generated SciPy 1.17 batch results from
`tools/generate_optimize_elementwise_arrays_reference.py`. Tests check
independent statuses, iteration/evaluation counts, successful/failed problems,
per-element bounds, scalar and rank-one argument broadcasting, callbacks,
invalid batch dimensions, and empty arrays. These tests use the same scalar
callbacks as the non-batched routines.

The seven pre-existing scalar tests and the new rank-one test passed in both
debug and optimized GNU Fortran builds with `-std=f2018`,
`-Werror=implicit-interface`, `-fcheck=all`, and IEEE invalid/zero/overflow
traps. The public umbrella module was compiled using test-only stub modules
for unavailable external solver packages. The full FPM suite is still unverified.

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

Validation dates: 2026-09-20 and 2026-10-01 through 2026-10-02. GNU Fortran 14.2.0 on Linux.

`test/test_descriptive.f90` compares mean, variance, standard deviation,
central moments through order four, linear quantiles, covariance, Pearson
correlation, and all five `rankdata` methods against values reproduced by
`tools/generate_statistics_reference.py` with NumPy 2.3.5 and SciPy 1.17.0.
The test also covers empty arrays, invalid `ddof`, invalid quantiles and rank
methods, NaN propagation, infinities in order statistics, and data near the
binary64 overflow limit.

`test/test_descriptive_scipy_compat.f90` adds SciPy 1.17.0 reference checks for
weighted/unweighted `gmean`, `hmean`, and `pmean`; `sem` and `variation` with
multiple `ddof` values; IQR percentile ranges/scaling; median absolute
deviation scaling; moments about the mean and an explicit center; biased and
bias-corrected skewness; and Fisher/Pearson kurtosis. Domain and shape tests
cover zero/negative data, constant samples, empty arrays, invalid weight
vectors, and invalid percentile ranges. Extreme power-mean checks reproduce
SciPy's overflow-to-infinity and underflow-to-zero behavior.

For the 2026-10-01 through 2026-10-02 descriptive extensions, the 212-source
dependency closure required by
`scifort_stats` compiled with GNU Fortran 14.2.0 using Fortran 2018,
`-O0 -g -fcheck=all -Wall -Wextra -Wimplicit-interface`. The original descriptive
regression, SciPy-compatibility regression, rank-2 axis/keepdims regression, and
`test/test_descriptive_nan_policy.f90`, `test/test_descriptive_rank3.f90`,
`test/test_descriptive_rank3_weights.f90`, `test/test_descriptive_rank4.f90`,
`test/test_descriptive_highrank_means.f90`,
`test/test_descriptive_highrank_weights.f90`,
`test/test_descriptive_highrank_stats.f90`, and
`test/test_descriptive_named_scale.f90`, and
`test/test_descriptive_standardize.f90` all link through the public
`scifort_stats` facade and pass. The unweighted rank-3 regression compares all
ten helpers on every supported axis, whole-array IQR,
fixed-rank keepdims, invalid axes, NaN propagation, and NaN omission against
NumPy 2.3.5 / SciPy 1.17.0 references. The weighted rank-3 regression checks
rank-1/rank-2/rank-3 trailing weights, joint singleton expansion of data and
weights, fixed-rank keepdims, invalid broadcasts, NaN-weight omission, and
propagation against the same SciPy release. The rank-4 regression checks all ten
helpers across representative axes 1 through 4, whole-array IQR, fixed-rank
keepdims, invalid axes, and rank-4 NaN propagation/omission against SciPy 1.17.0.
The high-rank mean regression exercises `gmean`, `hmean`, and `pmean` at every
rank from 5 through 15, with rank-15 SciPy references for non-default axes,
keepdims, NaN propagation/omission, invalid axes, and zero extents. The
high-rank weighted regression exercises all three weighted means at every rank
from 4 through 15, every assumed-rank weight branch through full rank, scalar
weights, joint singleton expansion, invalid broadcasts/excess weight rank,
weighted keepdims, and NaN-weight propagation/omission against SciPy 1.17.0. The
high-rank descriptive-statistics regression exercises `sem`, `variation`, `iqr`,
`median_abs_deviation`, `moment`, `skew`, and `kurtosis` at every rank from 5
through 15. Rank-15 SciPy 1.17.0 comparisons cover non-default axes, fixed-rank
keepdims, whole-array IQR, NaN propagation/omission, invalid axes, and nonlinear
third-moment, unbiased-skewness, and SEM references. The same test passes with
`-ffpe-trap=invalid,zero,overflow`. The NaN-policy regression checks
vector/rank-2
`propagate` and `omit` behavior, all-NaN samples, weighted observations with NaN
weights, trailing broadcast weights, matrix `axis=None` IQR, and fixed-rank
keepdims. Its reference values are generated with NumPy 2.3.5 / SciPy 1.17.0 by
`tools/generate_statistics_reference.py`. The regression also passes with
`-ffpe-trap=invalid,zero,overflow`, verifying that NaN propagation/omission and
power-mean overflow/underflow do not depend on trapped intermediate operations.
Earlier checkpoints also passed an independent optimized closure build. The
current direct GNU Fortran integration build links and runs all 22
`test/test_descriptive*.f90` programs through `scifort_stats` with
`-ffpe-trap=invalid,zero,overflow`.

`test/test_descriptive_named_scale.f90` exercises the named normal scale for
`iqr` and `median_abs_deviation` at every rank from 1 through 15, including
whole-array IQR, fixed-rank keepdims, a non-default percentile range,
case-insensitive `normal`, invalid names, and `nan_policy='propagate'` /
`'omit'`. Reference values are generated with SciPy 1.17.0 by
`tools/generate_statistics_reference.py`.

`test/test_descriptive_frequency.f90` validates `percentileofscore`,
`scoreatpercentile`, `cumfreq`, and `relfreq` against SciPy 1.17.0 numerical
references. Coverage includes all four percentile-rank tie conventions, scalar and
rank-15 score arrays, whole-array and axis percentile reductions, scalar/vector
percentiles, interpolation methods, closed value limits, NaN ordering and omission,
weighted histograms, explicit/default real limits, empty/constant samples, and the
current SciPy multidimensional `relfreq` normalization convention.
`test/test_descriptive_frequency_errors.f90` uses child processes to verify invalid
policies, kinds, percentile ranges/methods/limits, histogram bin counts, nonfinite
automatic ranges, and weight-shape mismatches. Both pass with floating-point traps.

`test/test_descriptive_rankdata.f90` validates multidimensional `rankdata` and
`tiecorrect` against SciPy 1.17.0 references. Coverage includes all five tie methods,
whole-array Fortran-order flattening, both matrix axes, NaN propagation/omission, empty
input, legacy invalid-method behavior, tie factors with NaNs, and rank-15 axis/all
paths. `test/test_descriptive_rankdata_errors.f90` uses child processes to verify
`nan_policy='raise'` and invalid NaN-policy termination. Both pass with
`-ffpe-trap=invalid,zero,overflow`.

`test/test_descriptive_summary.f90` validates `mode`, `mode_all`, `mode_keepdims`,
`describe`, and `describe_all` against SciPy 1.17.0 references. Coverage includes modal
ties, NaN-as-a-category propagation, omission, empty mode input, matrix/keepdims and
rank-15 reductions, ordinary and bias-corrected description, arbitrary `ddof`, and SciPy's
masked-invalid omit path. `test/test_descriptive_summary_errors.f90` uses child processes
to verify invalid axes/policies, NaN raising, and empty `describe` termination. Both pass
with `-ffpe-trap=invalid,zero,overflow`.

`test/test_descriptive_standardize.f90` exercises `zscore`, `gzscore`, `gstd`,
and `gstd_keepdims` through rank 15. SciPy 1.17.0 comparisons cover both matrix
axes, `ddof=0`, `ddof=1`, negative `ddof`, NaN propagation/omission, geometric
log-domain standardization, invalid axes/policies, constants, infinities, and
rank-15 zero extents. The same regression passes with
`-ffpe-trap=invalid,zero,overflow`. `test/test_descriptive_nan_raise.f90`
checks the `'raise'` policy with child-process invocations so expected
`error stop` termination is verified without making the parent FPM test fail.
It covers vector and high-rank descriptive reductions, weighted NaNs,
standardization, whole-array reductions, tuple-axis `gstd`, and `zmap`'s
comparison-only raise semantics.

The implementation scales finite data before centered powers and products.
The strict floating-point-trap build verifies that a sample containing
`+1e308` and `-1e308` can return a finite standard deviation and an infinite
variance without triggering an intermediate overflow exception.

Known limitations:

- The SciPy-style helper slice supports `nan_policy='propagate'`, `'omit'`, and
  `'raise'` through rank 15 for all ten descriptive helpers. Masked-array
  semantics are not implemented.

- `quantile` currently implements only linear / Hyndman-Fan type 7
  interpolation and does not itself expose a NaN-omission policy.
- `gmean`, `hmean`, and `pmean` support broadcast weights through rank 15.
  Rank-3 calls accept rank-1/rank-2/rank-3 trailing weights; rank-4 through
  rank-15 calls additionally accept scalar or trailing-rank weights up to the
  input rank with NumPy-style joint singleton expansion. The other descriptive
  helpers in this compatibility slice are unweighted.
- `iqr` and `median_abs_deviation` expose SciPy-compatible named
  `scale='normal'` through rank 15. Dtype selection, custom MAD center callbacks,
  and alternate percentile interpolation modes are not yet exposed.

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
- higher-rank reductions and broadcasting that expands the input shape remain future reduction work;
- complex-valued SciPy special-function coverage is family-specific; remaining unsupported overloads are tracked in `ROADMAP.md`.


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
pair of coarse symbols in the strength-two OA construction. The public
`qmc_random_cd` optimizer is checked for deterministic replay, strict centered-
discrepancy improvement on the fixed six-point reference design, exact preservation
of each coordinate column, and bounded proposal/work counters. Engine-level
`optimization='random-cd'` coverage verifies non-increasing centered discrepancy for
LatinHypercube, Halton, Sobol, and unit-cube PoissonDisk draws. Optimized LHS output
is also checked to retain one point in every one-dimensional stratum. The public
`qmc_lloyd` routine is checked against a fixed one-iteration finite-Voronoi
reference generated with SciPy 1.17.0/Qhull using `Qbb Qc Qz` (no `QJ` joggle),
with deterministic replay and unit-cube invariants. Engine-level
`optimization='lloyd'` coverage exercises LatinHypercube, Halton, Sobol, and
unit-cube PoissonDisk draws. The fixed reference intentionally isolates the
mathematical Lloyd relaxation from SciPy's Qhull `QJ` perturbation; degenerate
geometries are documented as potentially different in SciFort.

The same QMC regression compares deterministic `MultinomialQMC` and
`MultivariateNormalQMC` outputs with SciPy 1.17.0, including a rank-one singular
covariance. The Box-Muller normal path is checked for exact SciFort-state replay and
finite output. Poisson-disk tests cover replay, minimum separation, unit and negative
rectangular bounds, surface proposals, fill-space exhaustion, and invalid parameter
paths. SciFort uses exact accepted-set neighbor checks, so the minimum-radius
invariant does not depend on SciPy's grid indexing implementation. Random-CD is
tested separately from that geometric invariant because post-sampling coordinate
exchange can intentionally destroy the original Poisson-disk separation; custom
PoissonDisk bounds with random-CD are rejected because centered discrepancy is a
unit-hypercube criterion.

For the random-CD and Lloyd extensions, the updated QMC regression passes with
GNU Fortran 14.2.0 both under `-fcheck=all -ffpe-trap=invalid,zero,overflow` and
independently under `-O2`. A clean checked source-order build of all 183 modules
required by the public `scifort_stats` umbrella also links and runs `test_qmc`
successfully with the Lloyd-enabled QMC modules. The QMC reference generator
reproduces the fixed unjoggled-Qhull Lloyd reference deterministically.

After the original QMC slice, the complete source-order graph contained 178 Fortran library
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

BWS extension validation date: 2026-10-03. GNU Fortran 14.2.0 on Linux.
The same regression program now checks `bws_test` against SciPy 1.17.0 exact
permutation results for the canonical Neuhaeuser example and a tied-rank
example under all three alternatives. It additionally exercises the randomized-
permutation path with explicit `rng_state` and invalid argument handling. The
230-source `scifort_stats` dependency closure compiles in Fortran 2018 mode with
`-fcheck=all -ffpe-trap=invalid,zero,overflow`, and the public-facade regression
passes under those flags. The SciPy reference module regenerates byte-for-byte
identically.

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

## Fast common Bessel special functions

Validation date: 2026-09-23. GNU Fortran 14.2.0 on Linux, Python 3.13.5,
SciPy 1.17.0, NumPy 2.3.5, and mpmath 1.3.0. `fpm` and `fprettify` were not
installed in the validation environment.

This slice adds `j0`, `j1`, `y0`, `y1`, `i0`, `i1`, `k0`, `k0e`, `k1`, and
`k1e` to `scifort_special`. The new regression covers ordinary SciPy 1.17.0
reference values, J/I parity, Y/K domains and zero endpoints, NaN/infinity
behavior, elemental array evaluation, scaled K behavior, and large-I values.
The `i0(713)` and `i1(713)` references were evaluated independently with
80-digit mpmath 1.3.0 values because SciPy 1.17.0's fast Cephes wrappers
already overflow there while the binary64 mathematical results remain finite.

Because FPM was unavailable, the required `fpm test` command could not be run.
A manual source-order build compiled all 185 library sources and 42 test-support
modules with:

```text
gfortran -std=f2018 -O0 -g -fcheck=all -Wall -Wextra \
  -Wimplicit-interface -Werror=implicit-interface -ffree-line-length-none
```

All 63 Fortran test programs then compiled, linked, and passed. The new and
pre-existing special-function regressions also passed when the same build added:

```text
-ffpe-trap=invalid,zero,overflow -fbacktrace
```

`python tools/check_ascii.py` passed. `python tools/generate_special_reference.py`
ran successfully with SciPy 1.17.0/NumPy 2.3.5/mpmath 1.3.0 and reproduced the
committed numerical constants used by the new regression. No third-party Bessel
implementation was imported, so no third-party license text was added.

On a system with FPM installed, run at minimum:

```text
fpm test
fpm test --flag "-O0 -g -fcheck=all -ffpe-trap=zero,overflow -fbacktrace -Wall -Wextra -std=f2018"
```

## Ordered cylindrical and modified Bessel special functions

Validation date: 2026-09-23. GNU Fortran 14.2.0 on Linux, Python 3.13.5,
SciPy 1.17.0, NumPy 2.3.5, and mpmath 1.3.0. `fpm` and `fprettify` were not
installed in the validation environment.

This slice adds integer-order `yn`/`kn` and real-order `iv`/`ive`/`kv`/`kve`
to `scifort_special`. `test/test_special_bessel_order.f90` covers positive and
negative integer orders, positive and negative real I orders, negative real
arguments at integer I order, the negative-noninteger I connection formula and
a sign-changing zero crossing, scaled/unscaled behavior, zero/infinity/NaN
domains, and elemental array evaluation. Reference values are from SciPy 1.17.0.

A separate development grid compared 115 ordered-function points with SciPy over
orders from -5.5 through 20 and positive arguments from `1e-6` through `100`.
The largest relative discrepancies observed were approximately `4.6e-14` for
`iv`, `2.7e-14` for `ive`, `9.5e-15` for `kv`/`kve`, `3.4e-14` for `yn`, and
`8.4e-15` for `kn`. This grid is a validation aid rather than a committed
reference fixture.

Because FPM was unavailable, the required `fpm test` command could not be run.
A manual source-order build compiled all 186 library sources and 42 test-support
modules with:

```text
gfortran -std=f2018 -O0 -g -fcheck=all -Wall -Wextra \
  -Wimplicit-interface -Werror=implicit-interface -ffree-line-length-none
```

All 64 Fortran test programs then compiled, linked, and passed. The new ordered
Bessel regression, the preceding classic-Bessel regression, and the elementary
special-function regression also passed when the same build added:

```text
-ffpe-trap=invalid,zero,overflow -fbacktrace
```

`python tools/check_ascii.py` passed. `python tools/generate_special_reference.py`
ran successfully with SciPy 1.17.0/NumPy 2.3.5/mpmath 1.3.0 and emitted the
ordered-Bessel reference values used by the regression. No third-party Bessel
implementation was imported, so no third-party license text or dependency was
added.

On a system with FPM installed, run at minimum:

```text
fpm test
fpm test --flag "-O0 -g -fcheck=all -ffpe-trap=zero,overflow -fbacktrace -Wall -Wextra -std=f2018"
```

## Real spherical Bessel special functions

Validation date: 2026-09-23. GNU Fortran 14.2.0 on Linux, Python 3.13.5,
SciPy 1.17.0, NumPy 2.3.5, and mpmath 1.3.0. `fpm` and `fprettify` were not
installed in the validation environment.

This slice adds real `spherical_jn`, `spherical_yn`, `spherical_in`, and
`spherical_kn` to `scifort_special`, with a SciPy-style optional logical
`derivative` argument. `test/test_special_spherical_bessel.f90` covers ordinary
SciPy 1.17.0 reference values, positive and negative arguments, derivative
parity, zero/infinity/NaN behavior, elemental arrays, tiny high-order values,
and large finite modified-K results near binary64 overflow.

A deterministic development grid compared 799 positive-argument points with
SciPy 1.17.0 over orders 0 through 30 plus 40, 50, 75, 100, 150, 200, 300,
and 500, with arguments from `1e-12` through 700 and extra points around each
order. Away from SciPy's tiny-argument derivative anomalies, the largest
relative discrepancies observed were approximately `1.8e-13` for
`spherical_jn`, `3.0e-14` for `spherical_yn`, `5.6e-12` for `spherical_in`,
and `4.1e-12` for `spherical_kn`/its derivative. A separate 132-point negative-
argument grid through order 100 agreed to about `1.6e-13` relative or better
when the reference magnitude exceeded `1e-12`.

Tiny derivatives and near-overflow K values were checked independently with
80- to 800-digit mpmath evaluations. For example, the mathematical limits give
`spherical_jn(1,1e-300,derivative=.true.)` and
`spherical_in(1,1e-200,derivative=.true.)` as one third, and mpmath gives
`spherical_in(22,1e-12,derivative=.true.) =
8.6703637264923585e-280`. The independent high-precision values
`spherical_kn(150,1) = 5.8857785702867412e306` and
`spherical_kn(150,0.99) = 2.6847508422762158e307` remain finite in binary64.
SciFort's direct differentiated series and logarithmic terminating-K formula
preserve these cases; SciPy 1.17.0 has visible loss in some of these extreme
reference evaluations, so they are not used as SciPy-regression targets.

Because FPM was unavailable, the required `fpm test` command could not be run.
A manual source-order build compiled all 187 library sources and 42 test-support
modules with:

```text
gfortran -std=f2018 -O0 -g -fcheck=all -Wall -Wextra \
  -Wimplicit-interface -Werror=implicit-interface -ffree-line-length-none
```

All 65 Fortran test programs compiled, linked, and passed. The new spherical
Bessel regression and the three preceding special-function regressions also
passed when the same build added:

```text
-ffpe-trap=invalid,zero,overflow -fbacktrace
```

The IEEE-trap run found and corrected two overflow-guard expressions that had
relied on non-guaranteed short-circuit evaluation; the final implementation
uses nested guards whose tests are themselves overflow-safe.

`python tools/check_ascii.py` passed. `python tools/generate_special_reference.py`
ran successfully with SciPy 1.17.0/NumPy 2.3.5/mpmath 1.3.0 and emitted the
spherical reference values used by the regression. No third-party spherical-
Bessel implementation was imported, so `THIRD_PARTY_LICENSES.md` and project
dependencies are unchanged.

On a system with FPM installed, run at minimum:

```text
fpm test
fpm test --flag "-O0 -g -fcheck=all -ffpe-trap=zero,overflow -fbacktrace -Wall -Wextra -std=f2018"
```

## Real-order cylindrical Bessel J/Y special functions

Validation date: 2026-09-23. GNU Fortran 14.2.0 on Linux, Python 3.13.5,
SciPy 1.17.0, NumPy 2.3.5, and mpmath 1.3.0. `fpm` and `fprettify` were not
installed in the validation environment.

This slice adds real `jv`, `jve`, `yv`, and `yve`, plus real-order `jn` as the
SciPy-compatible alias of `jv`. The positive-order numerical core is a modern
adaptation of the permissively licensed MATH77 `DBESJN`, `DBESYN`, and `DBESPQ`
routines. `test/test_special_bessel_jy.f90` covers small, middle-recurrence, and
large-asymptotic regions; negative and near-integer orders; valid and invalid
negative arguments; scaled real-axis aliases; zero/infinity/NaN behavior;
reflection overflow/underflow; and elemental arrays.

Development comparisons against SciPy 1.17.0 covered positive and negative
orders (including near integers and half/quarter orders), arguments from
`1e-12` through `1e4`, turning-point neighborhoods near `x=v`, and orders into
the thousands. On ordinary finite values, discrepancies were generally at the
few-ulp to low-`1e-12` relative level; an observed relative outlier near a Y
zero was about `1.0e-9` while its absolute error was about `6e-15`. A dedicated
near-integer grid around orders 0, 1, and 5 agreed to about `4.5e-14` relative
for J and `6.0e-15` for Y. High-order cases through order 5000.25 were also
checked at representative small, turning-point, and large arguments.

Several extreme reflection cases were independently checked with high-precision
mpmath because SciPy 1.17.0's AMOS-backed path can overflow or underflow before
the real mathematical result requires it. For example, for negative
noninteger orders near -95.586 and -71.377 at very small positive arguments,
mpmath confirms the negative sign of an overflowing J result that SciPy reports
as positive infinity; SciFort preserves the identity's sign. Likewise,
`yv(-100.5, 1e-12)` reduces exactly to `jv(100.5, 1e-12)` and underflows to
zero, whereas the corresponding SciPy call reports an overflow artifact.
SciFort's real-only `jve`/`yve` deliberately use the mathematical real-axis
scale factor of one rather than reproduce such AMOS-specific extreme artifacts.

Because FPM was unavailable, the required `fpm test` command could not be run.
A manual source-order build compiled all 188 library sources and 42 test-support
modules with:

```text
gfortran -std=f2018 -O0 -g -fcheck=all -Wall -Wextra \
  -Wimplicit-interface -Werror=implicit-interface -ffree-line-length-none
```

All 66 Fortran test programs compiled, linked, and passed. The new real-order
J/Y regression and the four neighboring special-function regressions also
passed when the same build added:

```text
-ffpe-trap=invalid,zero,overflow -fbacktrace
```

`python tools/check_ascii.py` passed. `python tools/generate_special_reference.py`
ran successfully with SciPy 1.17.0/NumPy 2.3.5/mpmath 1.3.0 and emitted the
ordinary `jv`/`jve`/`yv`/`yve` reference values used by the regression. The
MATH77 license and exact Netlib source URLs are recorded in
`CODE_PROVENANCE.md` and `THIRD_PARTY_LICENSES.md`, and the complete MATH77
notice/disclaimer is retained in the adapted source.

On a system with FPM installed, run at minimum:

```text
fpm test
fpm test --flag "-O0 -g -fcheck=all -ffpe-trap=zero,overflow -fbacktrace -Wall -Wextra -std=f2018"
```

## Real-axis Bessel derivatives

Validation date: 2026-09-23. GNU Fortran 14.2.0 on Linux, Python 3.13.5,
SciPy 1.17.0, NumPy 2.3.5, and mpmath 1.3.0. `fpm` and `fprettify` were not
installed in the validation environment.

This slice adds `jvp`, `yvp`, `ivp`, and `kvp` to `scifort_special`. Each is a
pure elemental real-axis function with optional nonnegative integer derivative
order `n` (default 1); `n=0` delegates exactly to the corresponding base Bessel
function. The implementation follows DLMF 10.6.7 for J/Y and 10.29.5 for I/K.
Normalized binomial weights avoid the raw-coefficient overflow of a literal
finite-sum implementation, and finite terms are accumulated relative to their
largest magnitude so large intermediate partial sums do not overflow
unnecessarily. A negative derivative order returns NaN rather than raising the
`ValueError` used by SciPy's Python interface.

A deterministic development grid compared all four derivatives with SciPy
1.17.0 for derivative orders 0 through 8, eight real Bessel orders from -5.5
through 25.5, and nine positive arguments from `1e-6` through `100` (2592
function evaluations total). The largest finite relative discrepancies were
approximately `1.3e-13` for `jvp`, `6.3e-13` for `yvp`, `5.9e-14` for `ivp`,
and `1.1e-13` for `kvp`; no point exceeded a mixed `2e-10` relative/absolute
screen. The committed regression also checks a 1000th J derivative at `x=1`,
where the normalized-weight implementation agrees with SciPy 1.17.0.

The IEEE-trap validation found and corrected an initial overflow guard whose
`huge/max_term` test could itself overflow for tiny `max_term`. The final guard
only divides by a magnitude greater than one, so its test is overflow-safe.
The new derivative regression plus the five neighboring special-function tests
pass with `-ffpe-trap=invalid,zero,overflow`.

A 2026-09-25 Windows gfortran regression report exposed a separate portability
issue in the base integer-order J path: the platform `BESSEL_JN` intrinsic could
overflow internally for very large order with `x << n`, causing the existing
`jvp(0,1,1000)` regression to return NaN. SciFort now selects its existing
DLMF defining J power series when `x <= sqrt((n+1)/2)` for exact nonnegative
integer order. In this region the term ratio is at most one eighth initially and
then decreases, so the series is rapidly convergent and the logarithmic leading
factor can return zero directly when the result is below binary64 range. The
focused derivative regression passes both ordinary runtime checking and
`-ffpe-trap=invalid,zero,overflow` with GNU Fortran 14.2.0.

Because FPM was unavailable, the required `fpm test` command could not be run.
A manual source-order build compiled all 189 library sources and 42 test-support
modules with:

```text
gfortran -std=f2018 -O0 -g -fcheck=all -Wall -Wextra \
  -Wimplicit-interface -Werror=implicit-interface -ffree-line-length-none
```

All 67 Fortran test programs compiled, linked, and passed. The new derivative
regression plus the five neighboring special-function regressions also passed
when the same build added:

```text
-ffpe-trap=invalid,zero,overflow -fbacktrace
```

`python tools/check_ascii.py` passed and
`python tools/generate_special_reference.py` completed successfully with SciPy
1.17.0/NumPy 2.3.5/mpmath 1.3.0. No third-party implementation code or new
license text was added in this derivative slice.

On a system with FPM installed, run at minimum:

```text
fpm test
fpm test --flag "-O0 -g -fcheck=all -ffpe-trap=zero,overflow -fbacktrace -Wall -Wextra -std=f2018"
```

## Real Airy special functions

Validation date: 2026-09-23. GNU Fortran 14.2.0 on Linux, Python 3.13.5,
SciPy 1.17.0, NumPy 2.3.5, and mpmath 1.3.0. `fpm` and `fprettify` were not
installed in the validation environment.

This slice adds the pure elemental subroutines `airy` and `airye` to
`scifort_special`. Each returns Ai, Ai', Bi, and Bi' through scalar output
dummies; ordinary elemental invocation therefore also supports conformable
arrays. The implementation is original MIT-licensed SciFort code derived from
the Airy differential equation and NIST DLMF 9.6.2--9.6.9 and 9.7. A power
series is used for `abs(x) <= 1`; the outer real branches reuse SciFort's
existing real-order Bessel kernels. Positive scaled values are evaluated
directly through `kve`/`ive`, and a leading DLMF asymptotic branch covers
positive arguments so large that the Airy exponential argument would overflow.
No SciPy, AMOS, Cephes, or Fullerton implementation code was copied or adapted.

A 316-point development grid compared both four-component interfaces with
SciPy 1.17.0 across the series/Bessel transition and logarithmically spaced
positive and negative real arguments through magnitude 1000. On the positive
real subset through `x=100`, the largest finite relative discrepancies were
about `8.9e-13` for Ai, `9.0e-13` for Ai', `8.9e-14` for Bi, `1.5e-13` for
Bi', `8.9e-13` for scaled Ai, `9.0e-13` for scaled Ai', `4.8e-14` for scaled
Bi, and `2.6e-14` for scaled Bi'. Across the wider oscillatory negative grid,
relative error grows near zeros where absolute values are tiny; representative
absolute errors remained around `1e-12` or smaller through the sampled range.

Very large negative arguments were also checked independently with 80-digit
mpmath values rather than assuming SciPy's real-axis phase reduction is exact.
At `x=-1e6`, for example, SciFort Ai differs from the high-precision value by
about `2.9e-10` absolute, while SciPy 1.17.0 differs by about `7.0e-10` at the
same point. Other components have comparable large-phase sensitivity. This
reflects the phase-reduction accuracy of the existing real-order `jv` kernel,
not a separate Airy recurrence failure.

SciFort intentionally extends the positive real range in two ways. When
unscaled Bi/Bi' exceed binary64 range, `airy` returns positive infinity rather
than the NaN that SciPy 1.17.0's AMOS-backed path can produce after overflow.
`airye` remains finite and, when the Airy exponential argument itself would
overflow, uses the leading scaled asymptotic terms; at that magnitude omitted
corrections are below binary64 resolution. For negative real input, the real
SciPy convention is retained: scaled Ai/Ai' components are NaN and scaled
Bi/Bi' components equal their unscaled values.

The focused Airy regression was built and run with:

```text
python3 /mnt/data/manual_build_scifort.py \
  /mnt/data/scifort_airy_work/scifort-main \
  /mnt/data/build_airy_focus --run-tests --only test_special_airy
```

using compiler flags:

```text
gfortran -std=f2018 -O0 -g -fcheck=all -Wall -Wextra \
  -Wimplicit-interface -Werror=implicit-interface -ffree-line-length-none
```

The focused result was 1/1 passed. A second build added
`-ffpe-trap=invalid,zero,overflow -fbacktrace` and ran the Airy regression plus
six neighboring special-function regressions; all 7/7 passed.

Because FPM was unavailable, the required `fpm test` command could not be run.
The same strict manual source-order build compiled all 190 library sources and
42 test-support/reference modules. All 68 Fortran test programs compiled,
linked, and passed in four returned chunks: 20/20, 20/20, 15/15, and 13/13.
The long merge-regression program, aggregate `test_scifort`, QMC/resampling,
and all direct special-function regressions are included in that total.

`python tools/check_ascii.py` passed. `python tools/generate_special_reference.py`
completed successfully with SciPy 1.17.0/NumPy 2.3.5/mpmath 1.3.0 and emitted
the ordinary/scaled Airy reference values used for development and regression
checks.

The initial Airy batch deliberately deferred `itairy` and Airy zero tables until
dedicated large-argument and oscillatory algorithms were available. The
completion batch below adds those real APIs. Complex-valued Airy functions and
complex Bi zeros remain future work.

On a system with FPM installed, run at minimum:

```text
fpm test
fpm test --flag "-O0 -g -fcheck=all -ffpe-trap=zero,overflow -fbacktrace -Wall -Wextra -std=f2018"
```


### Airy completion: integrals and real zeros

This slice adds `itairy`, `ai_zeros`, and `bi_zeros` to `scifort_special`.
`itairy` is original SciFort code derived from the Airy ODE and DLMF 9.10. It
integrates the local ODE basis series term by term through the small-argument
region, uses composite 16-point Gauss-Legendre quadrature only as a bridge over
moderate oscillatory arguments, and switches to ODE-derived asymptotic
antiderivative expansions in the tails. The known DLMF limits
`integral_0^inf Ai = 1/3`, `integral_0^inf Ai(-t) dt = 2/3`, and
`integral_0^inf Bi(-t) dt = 0` are used explicitly at infinite limits.

The zero routines use the DLMF 9.9(iv) leading zero phases as starting values
for Newton refinement through SciFort's existing `airy` kernel. The first five
Ai, Ai', Bi, and Bi' roots and their cross-values were checked against both
SciPy 1.17.0 and 80-digit mpmath values. The SciFort roots agree with the
mpmath values to binary64 roundoff; for some fourth/fifth roots SciPy's legacy
routine differs by several ulps, so the committed references use the
high-precision values rather than copying the lower-accuracy SciPy outputs.

High-precision `itairy` development checks used mpmath 1.3.0 at 80 digits for
`x = 1, 5, 8, 10, 12, 20`. Across those points, all four SciFort integrals
agreed to about `1e-13` absolute or better, while the exponentially growing
positive-Bi integral retained near-binary64 relative accuracy. This independent
check is important because SciPy 1.17.0's Zhang/Jin-backed `itairy` loses
substantial accuracy in part of the moderate range; at `x=8`, for example,
the high-precision Ai integral is approximately `0.33333331724248357`, whereas
SciPy 1.17.0 returns approximately `0.23791459295898676`. SciFort returns
approximately `0.33333331724238491`.

A focused strict build executed:

```text
python3 /mnt/data/manual_build_scifort.py   /mnt/data/scifort_airy_complete_work/scifort-main   /mnt/data/build_airy_complete_focus3   --run-tests --only test_special_airy --only test_special_airy_complete
```

with GNU Fortran 14.2.0 and flags:

```text
-std=f2018 -O0 -g -fcheck=all -Wall -Wextra
-Wimplicit-interface -Werror=implicit-interface
```

Result: 2/2 focused Airy tests passed.

A separate IEEE-trap build added:

```text
-ffpe-trap=invalid,zero,overflow
```

and ran `test_special_airy`, `test_special_airy_complete`,
`test_special_bessel_classic`, `test_special_bessel_derivatives`,
`test_special_bessel_jy`, `test_special_bessel_order`,
`test_special_elementary`, and `test_special_spherical_bessel`.
Result: 8/8 passed.

The strict full suite compiled 190 library sources and 42 test-support/reference
sources. A monolithic test invocation exceeded the execution window without a
test failure, so the exact compiled tree was rerun in returned chunks. All
69/69 Fortran test programs passed. The returned chunks covered tests 1--30,
31--42, 43--50, 51--60, and 61--69; `test_merge_regressions`, aggregate
`test_scifort`, both Airy tests, and every neighboring special-function test are
included in that total.

`python tools/check_ascii.py` passed. `python
tools/generate_special_reference.py` completed with Python 3.13.5, NumPy 2.3.5,
SciPy 1.17.0, and mpmath 1.3.0, emitting SciPy API references plus independent
80-digit mpmath Airy-integral and zero references.

FPM and fprettify were unavailable in the validation environment. On a system
with FPM installed, run at minimum:

```text
fpm test
fpm test --flag "-O0 -g -fcheck=all -ffpe-trap=zero,overflow -fbacktrace -Wall -Wextra -std=f2018"
```

## Carlson symmetric elliptic integrals

Validation date: 2026-09-24. GNU Fortran 14.2.0 on Linux, Python 3.13.5,
SciPy 1.17.0, NumPy 2.3.5, and mpmath 1.3.0. Git was 2.47.3. `fpm` and
`fprettify` were not installed in the validation environment.

This slice adds the pure elemental real binary64 functions `elliprc`, `elliprd`,
`elliprf`, `elliprg`, and `elliprj` to `scifort_special`. The RC/RD/RF/RJ
fifth-order duplication kernels adapt Jacob Williams'
`carlson-elliptic-integrals` release 1.0.0, which refactors the public-domain
SLATEC Carlson routines and is distributed under BSD-3-Clause. SciFort removes
diagnostics/error-status side effects, normalizes arguments by homogeneity,
adds IEEE/domain wrappers, and implements the real Cauchy-principal-value RC
and RJ continuations. RG is formed from the standard Carlson RF/RD identity.
The exact upstream source and license URLs and full retained notice are recorded
in `CODE_PROVENANCE.md` and `THIRD_PARTY_LICENSES.md`.

A deterministic development comparison used 1,000 positive random quadruples
with arguments log-uniform on `1e-12` through `1e12` plus 250 real
principal-value RJ cases with the symmetric arguments and pole magnitude
log-uniform on `1e-8` through `1e8`. Against SciPy 1.17.0, the largest finite
relative discrepancies were approximately `5.5e-16` for RC, `6.2e-16` for RF,
`8.3e-16` for RD, `2.4e-12` for RJ, and `1.8e-14` for RG. The RJ maximum occurs
in a strongly cancelling principal-value case with a very small result. No case
exceeded both `5e-12` relative and `5e-13` absolute error.

Uniform extreme-scale values were checked independently with mpmath 1.3.0 at
80 decimal digits. For example, scaling `(1,2,3,4)` by `1e200` gives RF about
`7.2694593546890819854e-101`, RD about `2.9046028102899064423e-301`, RJ about
`2.3984809974956776218e-301`, and RG about
`1.4018470999908950994e100`; SciFort agrees to binary64 roundoff. The inverse
`1e-200` scaling is likewise retained without avoidable intermediate overflow.
Extremely *unbalanced* finite arguments spanning hundreds of orders of
magnitude can lose more digits because the inherited fifth-order duplication
kernels eventually see normalized arguments near underflow. This is documented
as a numerical limitation rather than hidden by relaxed regressions.

The focused IEEE-trap run used:

```text
python3 /mnt/data/manual_build_scifort.py \
  /mnt/data/scifort_next_work/scifort-main \
  /mnt/data/build_carlson_trap_final \
  --run-tests --trap \
  --only test_special_airy \
  --only test_special_airy_complete \
  --only test_special_bessel_classic \
  --only test_special_bessel_derivatives \
  --only test_special_bessel_jy \
  --only test_special_bessel_order \
  --only test_special_carlson_elliptic \
  --only test_special_elementary \
  --only test_special_spherical_bessel
```

with flags:

```text
-std=f2018 -O0 -g -fcheck=all -Wall -Wextra
-Wimplicit-interface -Werror=implicit-interface
-ffpe-trap=invalid,zero,overflow
```

Result: 9/9 special-function programs passed. The Carlson regression includes
ordinary SciPy references, RC/RJ principal values, symmetry and degeneracy
identities, elemental array invocation, NaN/infinity/domain endpoints, singular
zero cases, `1e+/-200` homogeneity, and signed overflow protection at still
smaller scales.

Because FPM was unavailable, the exact current tree was then compiled manually
with the same strict flags except the IEEE traps. This compiled 191 library
sources and 42 test-support/reference sources. All 70 Fortran test programs
compiled, linked, and passed. The returned chunks were tests 1, 2--15, 16--30,
31--40, 41, 42--48, 49--52, 53--55, and 56--70. Two broader chunk invocations
exceeded the execution window and were deliberately not counted; their ranges
were rerun completely in the smaller returned chunks above. The final total
therefore consists only of explicit returned passes. `test_merge_regressions`,
aggregate `test_scifort`, QMC/resampling, and every direct special-function
regression are included in the 70/70 result.

`python3 tools/check_ascii.py` passed. `python3
tools/generate_special_reference.py` completed successfully and emitted SciPy
1.17.0 ordinary/principal-value references plus 80-digit mpmath extreme-scale
Carlson references.

On a system with FPM installed, run at minimum:

```text
fpm test
fpm test --flag "-O0 -g -fcheck=all -ffpe-trap=zero,overflow -fbacktrace -Wall -Wextra -std=f2018"
```

Complex Carlson arguments remain future work. SciPy supports complex RC/RD/RF/
RG/RJ values; the current SciFort API intentionally covers the real branch only.

## Legendre elliptic validation (2026-09-24)

This batch adds the real Legendre complete/incomplete elliptic interfaces
`ellipk`, `ellipkm1`, `ellipe`, `ellipkinc`, and `ellipeinc`. The implementation
is original SciFort code based on the NIST DLMF 19.25 Carlson connection formulas
and DLMF 19.12.2 near-singularity expansion, reusing SciFort's existing Carlson
RF/RD/RG kernels. No new third-party implementation code or license was added.

The validation environment was:

```text
GNU Fortran (Debian 14.2.0-19) 14.2.0
Python 3.13.5
NumPy 2.3.5
SciPy 1.17.0
mpmath 1.3.0
```

`fpm` and `fprettify` were not installed in the execution environment. The
repository checks were therefore exercised with direct GNU Fortran compilation,
and the FPM commands at the end of this section remain required local checks.

A deterministic development grid compared approximately 2,300 complete and
incomplete evaluations with SciPy 1.17.0. The parameter set included ordinary
values in `[-100,1)`, complementary parameters down to about `1e-300`, negative
parameters with magnitudes up to `1e300`, and amplitudes up to magnitude 1000.
Additional amplitude checks extended to about `1e300`. The largest finite
relative discrepancies in the main grid were approximately:

```text
ellipk     7.3e-16
ellipkm1   7.3e-16
ellipe     7.6e-16
ellipkinc  8.9e-16
ellipeinc  9.2e-16
```

No finite/nonfinite mismatches occurred on that ordinary real-domain grid.
`ellipkm1` was also checked at the minimum positive binary64 subnormal and at
`1e-300` against mpmath at 100 decimal digits. For the exact minimum subnormal
input, the high-precision reference is approximately
`373.60633032181052177588811346595717`; SciFort preserves the corresponding
binary64 value. Extreme negative-parameter complete-E values were independently
checked with mpmath as well. The near-`m=1` DLMF 19.12.2 series avoids the loss
of precision that would result from forming complete E only through a highly
unbalanced RG call.

The focused IEEE-trap check used:

```text
python3 /mnt/data/manual_build_scifort.py \
  /mnt/data/scifort_legendre_work/scifort-main \
  /mnt/data/build_legendre_trap_final \
  --run-tests --trap \
  --only test_special_airy \
  --only test_special_airy_complete \
  --only test_special_bessel_classic \
  --only test_special_bessel_derivatives \
  --only test_special_bessel_jy \
  --only test_special_bessel_order \
  --only test_special_carlson_elliptic \
  --only test_special_elementary \
  --only test_special_legendre_elliptic \
  --only test_special_spherical_bessel
```

with compiler flags:

```text
-std=f2018 -O0 -g -fcheck=all -Wall -Wextra
-Wimplicit-interface -Werror=implicit-interface
-ffpe-trap=invalid,zero,overflow -fbacktrace
```

Result: 10/10 direct special-function programs passed. This run exposed and
then regression-tested a Fortran-specific issue in the first draft: logical
`.or.` does not guarantee short-circuit evaluation, so NaN checks are now kept
separate from relational comparisons that could signal `IEEE_INVALID`.

The exact final numerical tree was then compiled manually with the same strict
flags except the IEEE traps. This compiled 192 library source files and 42 test
support/reference source files. All 71 Fortran test programs linked and passed.
The explicit returned chunks were tests 1--15, 16--30, 31--40, 41 alone,
42--44, 45--50, 51--55, and 56--71. A broader 42--55 invocation exceeded the
execution window after explicitly returning passes only through test 44; the
unreturned range was rerun completely in the smaller chunks above and the
partial invocation was not counted. `test_merge_regressions`, aggregate
`test_scifort`, QMC/resampling, and every direct special-function regression are
included in the 71/71 result.

`python3 tools/check_ascii.py` passed. `python3
tools/generate_special_reference.py` completed successfully and emitted SciPy
1.17.0 complete/incomplete Legendre references plus 100-digit mpmath references
for tiny complementary parameters and extreme negative parameters.

The current API intentionally covers the real branch only. Inputs with `m > 1`
return NaN rather than complex continuations. At `m = 1`, incomplete-K paths
that cross the logarithmic singularity return infinity with the orientation
sign; SciPy 1.17.0 returns positive infinity for negative crossings. For
`m = -infinity`, SciFort uses the mathematical real limiting signs, which can
differ from SciPy edge artifacts. Jacobi `ellipj` is validated separately in
the following section.

On a system with FPM and fprettify installed, run at minimum:

```text
fprettify -r src test
fpm test
fpm test --flag "-O0 -g -fcheck=all -ffpe-trap=zero,overflow -fbacktrace -Wall -Wextra -std=f2018"
```

## Jacobi elliptic validation (2026-09-24)

This batch adds the pure elemental real subroutine
`ellipj(u,m,sn,cn,dn,ph)`, returning Jacobi `sn`, `cn`, `dn`, and the
unwrapped amplitude for SciPy's real parameter domain `0 <= m <= 1`. The
implementation is original SciFort code based on the arithmetic-geometric-mean
inversion in NIST DLMF 22.20. No third-party implementation source or new
license was added.

The validation environment was:

```text
GNU Fortran (Debian 14.2.0-19) 14.2.0
Python 3.13.5
NumPy 2.3.5
SciPy 1.17.0
mpmath 1.3.0
```

`fpm` and `fprettify` were not installed. Direct GNU Fortran validation used:

```text
-std=f2018 -O0 -g -fcheck=all -Wall -Wextra
-Wimplicit-interface -Werror=implicit-interface
```

A deterministic ordinary-domain comparison used 1,755 points: a structured
set spanning `m=0` through `m=1`, `abs(u)` through 100, plus 1,500 seeded
random pairs in `0 < m < 1` and `abs(u) < 100`. Against finite SciPy 1.17.0
references, the largest absolute discrepancies were approximately
`4.8e-14` for `sn`, `4.2e-14` for `cn`, `1.0e-13` for `dn`, and `5.7e-14`
for the amplitude. Near zeros, relative errors can look larger while the
absolute errors remain near binary64 roundoff.

Independent mpmath checks at 100 decimal digits cover regimes where SciPy
1.17.0's Cephes-backed near-`m=1` approximation becomes nonphysical or returns
NaN. At `u=1000` and binary64 `m=1-1e-12`, SciFort agrees with mpmath to about
`3e-14` absolute or better in all four outputs, while SciPy 1.17.0 returns
NaNs. At `abs(u)=1e6`, high-precision checks for `m=0.5`, `m=0.9`, and
`m=1-1e-12` show absolute errors from about `2e-11` through `2e-10` for the
ordinary cases and about `1.6e-12` for the near-one `cn`/`dn` case. These
large-argument errors come primarily from binary64 period reduction. Since the
period itself depends on binary64 `K(m)`, phase accuracy necessarily degrades
as `abs(u)` spans more periods; arguments far beyond `1e6` can lose additional
phase digits and this limitation is documented in the public API.

The focused regression was rebuilt with the strict flags above and passed.
The final IEEE-trap run added:

```text
-ffpe-trap=invalid,zero,overflow -fbacktrace
```

and executed the eleven direct special-function programs: Airy, Airy
completion, classic Bessel, Bessel derivatives, real-order Bessel J/Y,
integer/modified Bessel order, Carlson elliptic, elementary, Jacobi elliptic,
Legendre elliptic, and spherical Bessel. Result: 11/11 passed.

The exact final executable tree compiled 193 library source files and 42 test
support/reference source files. All 72 Fortran test programs linked and passed
with the strict non-trap flags. The explicit returned ranges were tests 1--15,
16--30, 31--40, 41 alone, 42--43, 44 alone, 45--48, 49--54, 55 alone, and
56--72. Broader invocations of 42--55 and 42--48 exceeded the execution window
and were not counted; every program in those ranges was rerun to completion in
the smaller returned groups above. `test_merge_regressions`, aggregate
`test_scifort`, QMC/resampling, multivariate coverage, and every direct special
regression are included in the 72/72 total.

`python tools/check_ascii.py` passed. `python
 tools/generate_special_reference.py` completed successfully with Python
3.13.5, NumPy 2.3.5, SciPy 1.17.0, and mpmath 1.3.0, emitting ordinary SciPy
references plus 100-digit mpmath large-argument/near-one Jacobi references.
Git was 2.47.3. `fpm` and `fprettify` were unavailable in the validation
environment.

On a system with FPM and fprettify installed, run at minimum:

```text
fprettify -r src test
fpm test
fpm test --flag "-O0 -g -fcheck=all -ffpe-trap=zero,overflow -fbacktrace -Wall -Wextra -std=f2018"
```

## Exponential and trigonometric integral validation (2026-09-24)

This batch adds the pure elemental real interfaces `exp1(x)`, `expi(x)`,
`expn(n,x)`, `sici(x,si,ci)`, and `shichi(x,shi,chi)`.  The `exp1`, `expi`,
`expn`, and `shichi` algorithms are original SciFort code using convergent
series, a modified-Lentz continued fraction, asymptotic expansions, and the
standard Ei/E1 identities.  The `sici` Chebyshev kernels adapt the MATH77
`DSI`, `DCI`, and `DCPVAL` routines under the Caltech/JPL BSD-style license
already retained by SciFort for the real-order Bessel J/Y batch.  No new
license family or external dependency was introduced.

The validation environment was:

```text
GNU Fortran (Debian 14.2.0-19) 14.2.0
Python 3.13.5
NumPy 2.3.5
SciPy 1.17.0
mpmath 1.3.0
```

`fpm` and `fprettify` were not installed.  Direct GNU Fortran validation used:

```text
-std=f2018 -O0 -g -fcheck=all -Wall -Wextra
-Wimplicit-interface -Werror=implicit-interface
```

The focused regression was run with:

```text
python /mnt/data/manual_build_scifort.py \
  /mnt/data/scifort_next_work/scifort-main \
  /mnt/data/build_expint_focus \
  --run-tests --only test_special_exponential_integrals
```

It compiled 194 library source files and 42 test support/reference files and
passed `test_special_exponential_integrals` (1/1).  The direct special-function
IEEE-trap run added:

```text
-ffpe-trap=invalid,zero,overflow -fbacktrace
```

and executed all twelve `test_special_*` programs: Airy, Airy completion,
classic Bessel, Bessel derivatives, real-order Bessel J/Y, integer/modified
Bessel order, Carlson elliptic, elementary, exponential/trigonometric
integrals, Jacobi elliptic, Legendre elliptic, and spherical Bessel.  Result:
12/12 passed.

A deterministic probe of the exact compiled Fortran routines produced 2,020
comparisons against SciPy 1.17.0 over ordinary finite real domains and a broad
range of orders for `expn`.  The largest observed finite relative differences
were approximately `9.2e-15` for `exp1`, `5.8e-14` for `expi`, `3.5e-15` for
`expn`, `1.4e-15` for `Si`, `1.8e-15` for `Shi`, and `1.6e-15` for `Chi`.
The largest relative difference for `Ci` was about `1.1e-12` at a value close
to zero, corresponding to only about `1.3e-15` absolute error.  The regression
also checks `expi(715)` against an 80-digit mpmath reference; SciFort retains a
finite binary64 result there while SciPy 1.17.0 overflows earlier on that path.

The exact final source/test tree was validated with the strict non-trap flags.
An initial all-in-one invocation exceeded the execution window after compiling
all 194 library files, 42 support/reference files, and linking the first test
executables, so that partial invocation was not counted as a suite result.  The
same compiled objects were then used to link and execute the complete sorted
suite in returned chunks: tests 1--15, 16--19, 20--30, 31--40, 41 alone,
42--43, 44 alone, 45--50, 51--58, 59--70, and 71--73.  Every program returned
success: 73/73 passed.  This includes `test_merge_regressions`, aggregate
`test_scifort`, QMC/resampling, multivariate coverage, and all twelve direct
special-function regressions.

The repository checks:

```text
python tools/check_ascii.py
python tools/generate_special_reference.py
```

both completed successfully.  The reference generator uses SciPy 1.17.0 for
ordinary finite API comparisons and mpmath 1.3.0 at 80 digits for the large
positive Ei references.

The public implementation is intentionally real-only.  `exp1` and `expn`
return NaN for negative finite real arguments rather than complex
continuations.  `expi` implements the principal real branch.  The finite
negative `Ci` result uses the real even continuation `Ci(abs(x))`; for
`x=-infinity`, `sici` follows SciPy 1.17.0 and returns NaN for `Ci`.  At very
large finite arguments, `sici` ultimately depends on the platform libm sine
and cosine argument reduction.  The positive `expi` asymptotic path favors
range and remains finite closer to the true binary64 overflow threshold than
SciPy 1.17.0, with the ordinary large-x comparison error documented above.

On a system with FPM and fprettify installed, run at minimum:

```text
fprettify -r src test
fpm test
fpm test --flag "-O0 -g -fcheck=all -ffpe-trap=zero,overflow -fbacktrace -Wall -Wextra -std=f2018"
```

## Real gamma-family validation (2026-09-24)

This batch adds the pure elemental real interfaces `gamma(x)`, `rgamma(x)`,
`gammasgn(x)`, real-axis `loggamma(x)`, and `multigammaln(a,d)`.  It also
corrects the existing `digamma(-0.0)` signed-zero convention to match SciPy
1.17.0.  The implementation is original SciFort code built on the existing
`gammaln` reflection kernel and the real gamma reflection/sign structure in
NIST DLMF 5.5.3.  No third-party source or new license was added.

The validation environment was:

```text
GNU Fortran (Debian 14.2.0-19) 14.2.0
Python 3.13.5
NumPy 2.3.5
SciPy 1.17.0
mpmath 1.3.0
Git 2.47.3
```

`fpm` and `fprettify` were not installed.  Direct GNU Fortran validation used:

```text
-std=f2018 -O0 -g -fcheck=all -Wall -Wextra
-Wimplicit-interface -Werror=implicit-interface
```

The focused regression was built and run with:

```text
python3 /mnt/data/build_scifort_manual.py \
  /mnt/data/scifort_next_work/scifort-main \
  /mnt/data/build_gamma_focus \
  --test test_special_gamma_family
```

It compiled all 195 library source files and passed
`test_special_gamma_family` (1/1).  The focused IEEE-trap build added:

```text
-ffpe-trap=invalid,zero,overflow -fbacktrace
```

and passed the same regression.  The final direct special-function trap suite
executed all 13 `test_special_*` programs from the trap build and returned
13/13 passes.

A deterministic development sweep evaluated 1,400 real arguments: 700 points
across the finite negative non-pole region and 700 geometrically spaced
positive values from `1e-300` through 170.  Against finite SciPy 1.17.0
references, after excluding two independently confirmed SciPy range artifacts,
the largest relative discrepancies were approximately `1.9e-13` for both
`gamma` and `rgamma`, about `2.2e-16` for real `loggamma`, and zero for
`gammasgn`.  The exceptional extreme cases were checked with mpmath 1.3.0 at
80 decimal digits.  For example, at `x=-170.75` SciFort retains the finite
binary64 values for both Gamma and reciprocal Gamma while SciPy 1.17.0
prematurely underflows/overflows on that path.  The regression also covers true
subnormal Gamma/reciprocal-Gamma values, signed zeros, poles, infinities, NaNs,
multivariate log-gamma values, and elemental array invocation.

The exact final executable tree was validated with the strict non-trap flags.
The library build contained 195 source files and the test layer contained 42
support/reference source files and 74 program files.  The 74 executables were
linked successfully and run in explicit returned ranges 1--20, 21--40, 41
alone, 42--48, 49--55, and 56--74.  Every executable returned success:
74/74 passed.  This includes `test_merge_regressions`, aggregate
`test_scifort`, QMC/resampling, multivariate coverage, and all 13 direct
special-function regressions.

The repository checks:

```text
python3 tools/check_ascii.py
python3 tools/generate_special_reference.py
```

completed successfully.  The reference generator uses SciPy 1.17.0 for
ordinary finite API comparisons and mpmath 1.3.0 at 80 digits for extreme
representable Gamma/reciprocal-Gamma cases where SciPy's real path can lose
range.

The public implementation is intentionally real-only.  `loggamma(x)` returns
NaN for negative finite real arguments instead of a complex continuation.
`multigammaln(a,d)` returns NaN for invalid real-domain arguments rather than
raising Python's `ValueError`.  `gamma` returns NaN at negative integer poles,
while `rgamma` returns zero there, matching SciPy 1.17.0 real semantics.

On a system with FPM and fprettify installed, run at minimum:

```text
fprettify -r src test
fpm test
fpm test --flag "-O0 -g -fcheck=all -ffpe-trap=zero,overflow -fbacktrace -Wall -Wextra -std=f2018"
```

## Real gamma combinatorics validation (2026-09-24)

This batch adds the pure elemental real interfaces `beta(a,b)`, generalized
`binom(n,k)`, `poch(a,m)`, and integer-order `polygamma(n,x)`.  It also extends
`betaln(a,b)` from the earlier positive-shape statistical domain to the real
meromorphic continuation used by SciPy where a finite real result is defined.
The implementation is original SciFort code.  It uses the existing SciFort
Gamma/sign/log-Gamma and Hurwitz-zeta kernels together with the standard DLMF
Gamma/Beta/Pochhammer/binomial relations and DLMF 5.15 polygamma recurrence
and reflection identities.  No third-party source or new license was added.

The validation environment was:

```text
GNU Fortran (Debian 14.2.0-19) 14.2.0
Python 3.13.5
NumPy 2.3.5
SciPy 1.17.0
mpmath 1.3.0
Git 2.47.3
```

`fpm` and `fprettify` were not installed.  Direct GNU Fortran validation used:

```text
-std=f2018 -O0 -g -fcheck=all -Wall -Wextra
-Wimplicit-interface -Werror=implicit-interface
```

The focused regression was built and run with:

```text
python3 /mnt/data/build_scifort_manual.py \
  /mnt/data/scifort_gamma_next_work/scifort-main \
  /mnt/data/build_gamma_comb_focus \
  --test test_special_gamma_combinatorics
```

The focused regression passed.  A focused IEEE-trap build added:

```text
-ffpe-trap=invalid,zero,overflow -fbacktrace
```

and passed both `test_special_gamma_combinatorics` and the broadened
`test_special_elementary` regression.  The final direct special-function trap
suite linked the same trap-built library objects against all 14
`test_special_*` programs and returned 14/14 passes.

A deterministic compiled-Fortran development sweep evaluated 2,900 cases:
1,125 `polygamma`, 625 `beta`, 625 `betaln`, 525 `poch`, and 625 `binom`
evaluations.  Against SciPy 1.17.0 ordinary real references, the largest
observed finite relative differences were approximately `1.8e-14` for
`polygamma`, `2.9e-15` for `beta`, `1.1e-15` for `betaln`, `2.5e-15` for
`poch`, and `5.0e-15` for `binom`, with no material mismatches at the tested
tolerances.  Negative-axis polygamma reflection was also checked independently
against mpmath 1.3.0 at 80 decimal digits for derivative orders 1 through 8
and arguments including `-20.25`, `-20.5`, `-50.25`, `-100.25`, and
`-1000.125`.  The maximum observed relative discrepancy in that set was about
`9.0e-16`.  Exact reduced sine/cosine handling at integer and half-integer
arguments is regression-tested because it avoids cancellation in far-negative
reflection cases such as `polygamma(8,-20.5)`.

The exact final executable tree was validated with the strict non-trap flags.
It contained 196 library source files, 42 test support/reference files, and 75
test program files.  The one-shot compile/link step exceeded the execution
window after linking the first executables, without a compiler failure; the
same successfully compiled objects were then used to finish linking all 75
programs.  Executions were run in explicit returned ranges 1--15, 16--30,
31--40, 41 alone, 42--43, 44 alone, 45--55, and 56--75.  Every executable
returned success: 75/75 passed.  This includes `test_merge_regressions`,
aggregate `test_scifort`, QMC/resampling, multivariate coverage, and all 14
direct special-function regressions.

The repository checks:

```text
python3 tools/check_ascii.py
python3 tools/generate_special_reference.py
```

completed successfully before the full-suite run.  They are rerun on the
frozen source below as part of delivery verification.  The reference generator
uses SciPy 1.17.0 for ordinary finite API comparisons and mpmath 1.3.0 for
independent high-precision negative-axis polygamma references.

The public implementation is intentionally real-only.  `polygamma` takes an
integer derivative order; negative orders return NaN rather than raising a
Python exception.  At simultaneous numerator poles of the meromorphic Beta
continuation, a path-independent finite real value does not exist; SciFort
represents the singular case symmetrically as infinity rather than reproducing
argument-order-dependent signed-infinity artifacts.  Very high derivative
orders combined with very large negative arguments can ultimately be limited
by binary64 coefficient growth in the reflection polynomial, although the
ordinary and tested negative ranges are validated near double precision.

On a system with FPM and fprettify installed, run at minimum:

```text
fprettify -r src test
fpm test
fpm test --flag "-O0 -g -fcheck=all -ffpe-trap=zero,overflow -fbacktrace -Wall -Wextra -std=f2018"
```

## Floating factorial/combinatorial validation (2026-09-24)

This batch adds the floating SciPy-style interfaces `factorial`, `factorial2`,
`factorialk`, `comb`, and `perm`.  The implementation is original SciFort code.
Real factorial evaluation uses the Gamma identity but forms it by recurrence
from an anchor `Gamma(1+f)` with `f` in `[0,1)`, so large values do not depend on
exponentiating a nearly overflowing log-gamma result; generalized floating
combinations/permutations reuse the existing `binom` and `poch` kernels. Integer
factorial, multifactorial, and permutation actual arguments
use their defining finite products with explicit binary64 overflow guards so
small integer results do not inherit avoidable Gamma/Pochhammer rounding.  No
third-party implementation source or new license was added.

The public Fortran API corresponds to SciPy 1.17.0's floating `exact=False`,
`extend='zero'` behavior where it maps naturally to a statically typed scalar
interface.  `factorial` accepts either default-integer or binary64-real actual
arguments.  `factorial2` and `factorialk` use integer arguments for the default
zero-extension domain.  `comb` and `perm` accept either default integers or
binary64 reals; `comb` has the optional scalar logical `repetition`.  Python's
arbitrary-precision `exact=True` option and complex `extend='complex'`
continuations are not implemented.

The validation environment was:

```text
GNU Fortran (Debian 14.2.0-19) 14.2.0
Python 3.13.5
NumPy 2.3.5
SciPy 1.17.0
mpmath 1.3.0
```

`fpm` and `fprettify` were not installed.  Direct GNU Fortran validation used:

```text
-std=f2018 -O0 -g -fcheck=all -Wall -Wextra
-Wimplicit-interface -Werror=implicit-interface
```

The focused regression was built and run with the temporary out-of-tree build
helper:

```text
python3 /mnt/data/build_scifort_manual.py \
  /mnt/data/scifort_fact_work/scifort-main \
  /mnt/data/build_fact_focus test_special_factorials
```

It compiled all 197 library sources and passed
`test_special_factorials`.  Repeating the build with

```text
TRAPS=1 python3 /mnt/data/build_scifort_manual.py \
  /mnt/data/scifort_fact_work/scifort-main \
  /mnt/data/build_fact_trap test_special_factorials
```

added `-ffpe-trap=invalid,zero,overflow -fbacktrace` and also passed.  All 15
`test/test_special_*.f90` programs were subsequently linked against the
trap-built library and returned success: **15/15 passed**.

A deterministic compiled-Fortran development sweep evaluated 9,914 cases:
701 real `factorial`, 308 integer `factorial2`, 3,040 integer `factorialk`, and
1,955 each of ordinary `comb`, repeated `comb`, and `perm`.  Against SciPy
1.17.0 floating references, there were no material mismatches at a `3e-12`
relative criterion.  The largest observed finite relative differences were
approximately:

```text
factorial             1.5e-13
factorial2             1.2e-15
factorialk             3.9e-14
comb                   3.5e-14
comb(repetition=true)  5.8e-14
perm                   3.0e-14
```

The larger differences in that original sweep were ordinary Gamma/Pochhammer
last-bit accumulation, not domain disagreements. A 2026-09-25 Windows gfortran
report showed that the original real-factorial `exp(gammaln(n+1))` path missed
the committed `factorial(170.5)` reference by about `1.1e-13` relative because a
one-ulp large log-gamma error is amplified by the final exponential. The current
guarded recurrence returns the same value without that log-domain amplification;
the focused factorial regression passes both ordinary runtime checking and
`-ffpe-trap=invalid,zero,overflow` under GNU Fortran 14.2.0. Integer actual
arguments intentionally use direct products, so representable small integer
results such as `factorial(5)` and
`perm(10,3)` are exact in binary64 even where SciPy's approximate path can be
a few ulps away from the integer.

The exact source/test tree was then built with the strict non-trap flags.  It
contained 197 library source files, 42 test support/reference source files, and
76 test program files.  Executions were run in explicit returned ranges 1--20,
21--40, 41 alone (`test_merge_regressions`), 42--43, 44 alone
(`test_multivariate`), 45--50, 51--58, and 59--76.  The attempted combined
45--58 invocation exceeded the execution window before returning per-test
results and was not counted; the same executables were rerun in the smaller
returned ranges.  Every explicitly counted program returned success:
**76/76 passed**.  This includes merge regressions, aggregate `test_scifort`,
QMC/resampling, multivariate coverage, and all 15 direct special-function
regressions.

The repository/reference checks:

```text
python3 tools/check_ascii.py
python3 tools/generate_special_reference.py
```

completed successfully on the working tree after the API/reference updates.

Known compatibility limitations are intentional.  Arbitrary-precision
`exact=True` results cannot be represented by these binary64-returning Fortran
interfaces.  Complex factorial/multifactorial continuation is not implemented.
For `factorialk`, invalid `k <= 0` returns NaN rather than raising Python's
`ValueError`.  The default Fortran integer kind is used by the integer overloads;
callers needing integer widths beyond that range should use the real generalized
interfaces where applicable or a separate arbitrary-precision package.

On a system with FPM and fprettify installed, run at minimum:

```text
fprettify -r src test
fpm test
fpm test --flag "-O0 -g -fcheck=all -ffpe-trap=zero,overflow -fbacktrace -Wall -Wextra -std=f2018"
```

## Real zeta and Spence validation (2026-09-24)

This batch adds the real SciPy-style interfaces `spence(x)`, `zeta(s[,q])`,
and `zetac(s)`.  The implementation is original SciFort code.  `spence`
evaluates `Li_2(1-x)` from the DLMF dilogarithm series and connection formulas
25.12.1, 25.12.3, 25.12.4, and 25.12.6.  The existing positive-q Hurwitz-zeta
Euler-Maclaurin kernel remains the base for `zeta(s,q)` and `zetac`; real
Riemann analytic continuation uses DLMF 25.4.2.  No third-party implementation
source or new license was added.

The validation environment was:

```text
GNU Fortran (Debian 14.2.0-19) 14.2.0
Python 3.13.5
NumPy 2.3.5
SciPy 1.17.0
mpmath 1.3.0
Git 2.47.3
```

`fpm` and `fprettify` were not installed.  Direct GNU Fortran validation used:

```text
-std=f2018 -O0 -g -fcheck=all -Wall -Wextra
-Wimplicit-interface -Werror=implicit-interface
```

The focused regression was built from all 198 library sources with:

```text
python3 /mnt/data/build_scifort_manual.py \
  /mnt/data/scifort_next_work \
  /mnt/data/build_zeta_focus test_special_zeta_spence
```

and passed.  Rebuilding the same focused test with

```text
TRAPS=1 python3 /mnt/data/build_scifort_manual.py \
  /mnt/data/scifort_next_work \
  /mnt/data/build_zeta_trap test_special_zeta_spence
```

added `-ffpe-trap=invalid,zero,overflow -fbacktrace` and also passed.  All 16
`test/test_special_*.f90` programs were then linked against the trap-built
library and returned success: **16/16 passed**.

A deterministic compiled-Fortran development sweep evaluated 2,121 cases:
426 Spence values, 822 Riemann-zeta values, 804 `zetac` values, and 69
positive/negative-shift Hurwitz-zeta values.  Against SciPy 1.17.0 ordinary
real references, no material mismatch occurred.  Maximum observed finite
relative differences were approximately:

```text
spence       2.5e-15
zeta         7.6e-14
zetac        7.6e-14
zeta(s,q)    3.9e-14
```

The largest Riemann-zeta relative discrepancy was in the far-negative
reflection region near `s=-82.5`; ordinary positive and near-pole cases were
closer.  The generator independently checks representative negative-axis zeta
values and Spence values with mpmath at 80 decimal digits.  `zetac(1000)` is
checked with 400-digit mpmath because its correct value is about
`9.33263618503219e-302`; SciPy 1.17.0 returns zero there, while SciFort retains
the representable binary64 value by scaling the leading `2**(-s)` term.

The exact numerical source/test tree contained 198 library source files, 42
test-support/reference files, and 77 test programs.  All sources compiled and
all programs linked under the strict flags above.  The first one-shot build
command exceeded the execution window while linking tests after all library
and support sources had compiled successfully; linking resumed from those exact
objects.  Test execution was then completed in explicit returned ranges 1--20,
21--40, 41--44, 45--58, and 59--77.  Every executable returned success:
**77/77 passed**.

The repository/reference checks:

```text
python3 tools/check_ascii.py
python3 tools/generate_special_reference.py
```

completed successfully after the API/reference updates and are rerun on the
freshly extracted delivery archive below.

The public implementation is intentionally real-only.  `spence` returns NaN
on its negative-real branch cut and for infinite arguments.  `zeta(s,q)`
supports the SciPy real domain for positive `q`; for negative noninteger `q`,
real continuation is supplied when `s` is an integer by separating the finite
negative-base terms from positive-q Hurwitz tails.  Noninteger `s` with
negative `q` returns NaN, and nonpositive integer `q` is a pole.  Complex zeta
and dilogarithm branches are not yet implemented.  Far-negative Riemann-zeta
values are ultimately limited by binary64 reflection-phase and log-Gamma
accuracy, though the tested range through about `s=-100` remains near double
precision.

On a system with FPM and fprettify installed, run at minimum:

```text
fprettify -r src test
fpm test
fpm test --flag "-O0 -g -fcheck=all -ffpe-trap=zero,overflow -fbacktrace -Wall -Wextra -std=f2018"
```

## Magnitude-ordered Bessel J/J-prime zero validation (2026-09-24)

This batch extends the existing real integer-order Bessel-zero implementation
with `jnjnp_zeros(nt, zeros, orders, serials, kinds[, status])`. The routine is
original SciFort code over the existing `jn_zeros`/`jnp_zeros` root kernels. It
constructs a finite order/root candidate grid and expands it until the requested
last zero is below every first omitted root. The exact identity
`J_0'(x) = -J_1(x)` is coalesced before sorting, with deterministic J-before-J'
tie ordering. No third-party implementation source or new license was added.

The validation environment was:

```text
GNU Fortran (Debian 14.2.0-19) 14.2.0
Python 3.13.5
NumPy 2.3.5
SciPy 1.17.0
mpmath 1.3.0
Git 2.47.3
```

`fpm` and `fprettify` were not installed. The focused regression and the full
manual suite were built from all 199 library sources and 42 test-support sources
with:

```text
-std=f2018 -O0 -g -fcheck=all -Wall -Wextra
-Wimplicit-interface -Werror=implicit-interface -fbacktrace
```

The focused `test_special_bessel_zeros` regression passed. All 78 repository
test programs were then linked and executed in four ranges; **78/78 passed**.
A focused rebuild of `test_special_bessel_zeros` additionally used

```text
-ffpe-trap=invalid,zero,overflow
```

and passed. A whole-suite trap-mode probe is not a clean repository-wide check:
legacy tests intentionally exercise invalid IEEE/domain inputs, and among the
first 40 programs `test_beta`, `test_discrete`, `test_discrete_extended`, and
`test_gamma` raised SIGFPE in those pre-existing paths. The ordinary runtime-
checking suite above is therefore the complete whole-repository result for this
batch.

A compiled Fortran dump of all 1200 supported merged zeros was compared with
SciPy 1.17.0 `special.jnjnp_zeros(1200)`. The maximum absolute zero difference
was `5.684341886080802e-14`; the maximum difference scaled by
`max(1,abs(reference))` was `9.645314177720275e-16`. Order/serial/type metadata
matched at 1196 of 1200 array positions. The four positional differences are
two `J_1`/`J_0'` pairs at serial numbers 11 and 15. Those are exact mathematical
ties; SciPy's separately computed binary64 approximations reverse those two
pairs by last-bit roundoff, whereas SciFort coalesces each pair and applies its
documented deterministic tie order. The 1200th entry is
`68.56313856707132` with order 32, serial 9, derivative kind 1; the largest
order appearing in the table is 65.

The repository/reference checks

```text
python3 tools/check_ascii.py
python3 tools/generate_special_reference.py
```

also completed successfully. The generator records the first 20 merged entries,
the final 1200th entry, and the maximum order from SciPy 1.17.0.

On a system with FPM and fprettify installed, run at minimum:

```text
fprettify -r src test
fpm test
fpm test --flag "-O0 -g -fcheck=all -fbacktrace -Wall -Wextra -std=f2018"
```

## Real Y0/Y1 specialized zero-helper validation (2026-09-24)

This batch adds real-branch `y0_zeros`, `y1_zeros`, and `y1p_zeros` wrappers
using the existing SciFort Y/Y-prime zero solver and `yvp`/`yv` kernels. The
paired outputs are the documented mathematical quantities: `Y_0'` at `Y_0`
zeros, `Y_1'` at `Y_1` zeros, and `Y_1` at `Y_1'` zeros. NIST DLMF 10.6.3 gives
`Y_0'=-Y_1`; SciFort follows that identity. SciPy 1.17.0's real `y0_zeros`
second array has the opposite sign despite the same documented identity, so
that sign is an intentional compatibility difference.

The first four real roots and companion values are regression-tested against
SciPy 1.17.0's generic `yn_zeros`/`ynp_zeros` roots plus `yvp`/`yv`, and the
141st entries are checked as an extension point. In SciPy 1.17.0 the three
legacy specialized routines return finite entries through 140 and NaNs from
entry 141 onward; SciFort's real wrappers remain finite there because they use
its general real zero kernels. The optional complex branch of the SciPy helper
APIs remains unimplemented.

## Real Dawson/error/Fresnel special functions (2026-09-25)

This batch adds real `dawsn`, `erfi`, `erfcx`, and paired-output `fresnel`.
The implementation is original SciFort code based on the defining relations and
real asymptotic expansions in NIST DLMF Chapter 7. Dawson uses a local power
series, bounded transformed-integral quadrature, and its large-real asymptotic
series. `erfi` is formed from Dawson with a logarithmic overflow guard. `erfcx`
uses the Fortran intrinsic `erfc` in its safe central range, its standard real
asymptotic expansion for large positive inputs, and reflection for negative
inputs. Fresnel S/C use power series, composite 16-point Gauss-Legendre
quadrature, and DLMF auxiliary-function asymptotics.

The large-Fresnel branch reduces the quadratic phase from a compensated square.
For exact binary64 integers at and above `2**52`, parity is handled explicitly.
This is a mathematical-correctness choice: SciPy 1.17.0 can lose the low bits of
`x**2` for large noninteger arguments before trigonometric reduction. Regression
therefore uses SciPy only in the range where its public outputs retain the
phase, and uses 80-digit mpmath values for selected large binary64 arguments.
No third-party implementation source was imported.


Validation environment:

```text
GNU Fortran (Debian 14.2.0-19) 14.2.0
Python 3.13.5
NumPy 2.3.5
SciPy 1.17.0
mpmath 1.3.0
Git 2.47.3
```

`fpm` and `fprettify` were not installed. The baseline
`test_special_elementary` regression passed before editing. The final working
source was compiled from all 200 library sources with:

```text
-std=f2018 -O0 -g -fcheck=all -Wall -Wextra
-Wimplicit-interface -Werror=implicit-interface -fbacktrace
```

The focused `test_special_error_fresnel` regression passed. A separate focused
build adding

```text
-ffpe-trap=invalid,zero,overflow
```

also passed. The complete ordinary regression suite passed **79/79** test
programs. The repository checks

```text
python3 tools/check_ascii.py
python3 tools/generate_special_reference.py
```

also passed.

A compiled dense sweep over 4001 points on `[-10,10]` was compared with SciPy
1.17.0. Maximum errors scaled by `max(1,abs(reference))` were:

```text
dawsn   8.33e-16
erfi    7.39e-15
erfcx   5.25e-15
Fresnel S 2.11e-15
Fresnel C 2.45e-15
```

For the large Fresnel phase cases `1000000.125`, `100000000.25`,
`134262676.7054991`, `4503599627370497.0`, and `9007199254740992.0`, compiled
SciFort outputs agreed exactly after binary64 rounding with 80-digit mpmath
`fresnels`/`fresnelc` results evaluated at the exact binary64 input rationals.
The latter cases intentionally do not use SciPy 1.17.0 as the numerical oracle
because its direct quadratic phase loses low bits there.

On systems with FPM and fprettify installed, run at minimum:

```text
fprettify -r src test
fpm test
fpm test --flag "-O0 -g -fcheck=all -fbacktrace -Wall -Wextra -std=f2018"
```

## Complex Faddeeva and complex erfcx validation (2026-09-25)

This batch adds complex binary64 `wofz(z)` and overloads the top-level
`scifort_special` `erfcx` generic so existing real arguments retain the prior
real implementation while complex arguments evaluate `erfcx(z)=wofz(i*z)`.
The implementation is original SciFort code based on NIST DLMF Chapter 7:
the defining relation in 7.2.3, reflection identity 7.4.3, entire series 7.6.3,
upper-half-plane integral representation 7.7.2, and continued fraction 7.9.3.
The transition-region integral is evaluated with independently written
composite 16-point Gauss-Legendre quadrature. No third-party implementation
source was imported.

Validation environment:

```text
GNU Fortran (Debian 14.2.0-19) 14.2.0
Python 3.13.5
NumPy 2.3.5
SciPy 1.17.0
mpmath 1.3.0
Git 2.47.3
```

`fpm` and `fprettify` were not installed. Before editing, the untouched
Dawson/Fresnel baseline passed:

```text
python3 /mnt/data/scifort_manual_build.py \
  /mnt/data/scifort_complex_base/scifort-main \
  /mnt/data/build_complex_baseline_check test_special_error_fresnel
```

The final working source compiled all 201 library sources with:

```text
-std=f2018 -O0 -g -fcheck=all -Wall -Wextra
-Wimplicit-interface -Werror=implicit-interface -fbacktrace
```

The focused `test_special_faddeeva` regression passed. A separate focused build
adding

```text
-ffpe-trap=invalid,zero,overflow
```

also passed. The complete ordinary regression suite passed **80/80** test
programs. Repository checks

```text
python3 tools/check_ascii.py
python3 tools/generate_special_reference.py
```

both passed.

A compiled 97-by-97 grid covering complex `x+i*y` with `x,y` from -12 to 12 in
steps of 0.25 was compared with SciPy 1.17.0. Over all 9409 points, the maximum
error scaled by `max(1,abs(reference))` was:

```text
wofz          1.0362907589042387e-15
complex erfcx 1.0362907589042387e-15
```

The worst `wofz` grid point was `5.75-1.0i`; the worst complex-`erfcx` point
was `-1.0-5.75i`. The regression test also covers each numerical branch
(series, real-axis Taylor continuation, transition quadrature, large-argument
continued fraction, and lower-half-plane reflection), the real-axis Dawson
identity, Faddeeva reflection, elemental arrays, NaN handling, and the
positive-imaginary-infinity limit.

On systems with FPM and fprettify installed, run at minimum:

```text
fprettify -r src test
fpm test
fpm test --flag "-O0 -g -fcheck=all -fbacktrace -Wall -Wextra -std=f2018"
```


## Windows portability repairs: high-order J derivatives and near-overflow factorial

Validation date: 2026-09-25. A user-provided Windows `fpm test` log exposed two
platform-sensitive regressions that did not reproduce in the earlier Linux-only
suite: `jvp(0,1,1000)` returned NaN after internal overflow in the platform
integer-order `BESSEL_JN` path, and `factorial(170.5)` missed its committed
reference by about `1.1e-13` relative because the large log-gamma value lost one
ulp before exponentiation. The existing committed regressions already covered
both cases, so no tolerance was relaxed and no expected value was changed.

For exact nonnegative integer J order, SciFort now uses its existing defining
power series whenever

```text
x <= sqrt((n + 1) / 2)
```

so the initial term ratio `x**2/(4*(n+1))` is at most `1/8` and decreases
thereafter. This bypasses platform recurrence overflow in the small-argument,
large-order regime. The series underflow guard was also moved from the smallest
normal binary64 value to the binary64 subnormal floor, represented in logarithmic
form as `log(tiny)+log(epsilon)`, so representable subnormal J values are not
silently flushed to zero.

For nonnegative real factorial arguments, SciFort now evaluates the identity

```text
Gamma(n + 1) = Gamma(1 + f) * product(f + j, j=1,...,floor(n))
```

where `f = n-floor(n)`, with an overflow guard on each multiplication. Anchoring
Gamma in `[1,2)` avoids magnifying a one-ulp large log-gamma error through a
near-overflow exponential. In the repaired build `factorial(170.5)` is
`9.483367566824795e307`, versus the SciPy 1.17.0 reference
`9.483367566824803e307`.

The validation environment was:

```text
GNU Fortran (Debian 14.2.0-19) 14.2.0
Python 3.13.5
SciPy 1.17.0
```

`fpm` and `fprettify` were unavailable. All 201 library sources and 42 test
support/reference sources were compiled out of tree with GNU Fortran using
`-std=f2018 -O0 -g -fcheck=all -fbacktrace`. The focused
`test_special_bessel_derivatives` and `test_special_factorials` programs passed
both that ordinary build and a build linked with
`-ffpe-trap=invalid,zero,overflow`. The complete manually linked test sweep then
returned success for **81/81** test programs.

A compiled development sweep exercised 4,004 exact integer-order J evaluations
for orders 0 through 1000 at four arguments spanning the newly selected series
region. No SciFort result was NaN or infinite. The maximum absolute difference
from SciPy 1.17.0 was about `5.6e-17`; SciPy itself underflowed 3,176 of these
tiny values to zero, while spot checks of SciFort's nonzero results against
high-precision mpmath confirmed the retained tiny values. A second compiled
sweep covered 1,711 real factorial arguments from 0 through 170.9 (including an
explicit duplicate check at 170.5); the maximum finite relative difference from
SciPy 1.17.0 Gamma references was about `6.9e-14`, while the repaired 170.5 case
was within about `8.4e-16` relative.

The repository checks

```text
python3 tools/check_ascii.py
python3 tools/generate_special_reference.py
```

also completed successfully. On the Windows/FPM environment that exposed the
original failures, rerun at minimum:

```text
fpm test
```

## Complex Fresnel validation (2026-09-25)

This batch overloads the top-level `fresnel` subroutine for complex binary64
arguments while preserving the existing real implementation. Complex S/C values
are evaluated independently from NIST DLMF 7.5.8 through the two error-function
interrelations for `C(z)+i*S(z)` and `C(z)-i*S(z)`. A defining entire-series path
is used when `abs(z) <= 20` and `abs(abs(real(z))-abs(aimag(z))) <= 1.5`, where
that series is well conditioned; elsewhere the already validated complex
`erf`/Faddeeva foundation is used. Exact real and imaginary axes map to the real
Fresnel kernel, including `C(i*y)=i*C(y)` and `S(i*y)=-i*S(y)`.

Validation environment:

```text
GNU Fortran (Debian 14.2.0-19) 14.2.0
Python 3.13.5
NumPy 2.3.5
SciPy 1.17.0
mpmath 1.3.0
```

`fpm` and `fprettify` were not installed. All 201 library sources compiled with
GNU Fortran using:

```text
-std=f2018 -O0 -g -fcheck=all -fbacktrace -Wall -Wextra
```

The focused `test_special_complex_fresnel` regression passed with those runtime
checks and also with:

```text
-ffpe-trap=invalid,zero,overflow
```

The complete manually linked ordinary regression suite passed **82/82** test
programs. Repository checks

```text
python3 tools/check_ascii.py
python3 tools/generate_special_reference.py
```

also passed.

A compiled 97-by-97 grid over `[-12,12] x [-12,12]` in quarter-unit steps was
compared with SciPy 1.17.0. Away from the coordinate axes the maximum error
scaled by `max(1,abs(reference))` was approximately `1.71e-13`; the worst cases
have Fresnel magnitudes above `1e143`. SciPy 1.17.0 has an axis-only discrepancy
for larger pure-imaginary arguments (for example at `z=4*i`), so SciPy's complex
axis values were not treated as authoritative there.

An independent 1089-point subset of the same compiled grid, on 0.75-unit spacing
and including both axes, was compared with 70-digit mpmath values formed from the
DLMF 7.5.8 error-function identities. Maximum scaled errors were:

```text
Fresnel S  6.77e-14
Fresnel C  6.77e-14
```

At `z=4*i`, SciFort returns
`S=-0.42051575424692844*i` and `C=0.4984260330381776*i`, agreeing with direct
high-precision integration and the analytic-continuation identities above.
SciPy 1.17.0's complex path instead returns approximately
`S=-1.4161704677316096*i` and `C=1.6785488324958335*i` at that exact point,
even though nearby non-axis complex values converge to the DLMF continuation.

On systems with FPM and fprettify installed, run at minimum:

```text
fprettify -r src test
fpm test
fpm test --flag "-O0 -g -fcheck=all -fbacktrace -Wall -Wextra -std=f2018"
```
## Modified Fresnel validation (2026-09-25)

This batch adds real-argument `modfresnelp` and `modfresnelm`. The scaled positive
branch uses `K_+(x)=wofz(exp(i*pi/4)*x)/2` for nonnegative `x`; `F_+` is reconstructed
with a compensated binary64 square phase. Negative lower limits use the documented
integral symmetry rather than the inconsistent first component returned by SciPy 1.17.0.
`F_-` and `K_-` are the complex conjugates of the positive pair for real `x`.

Validation environment:

```text
GNU Fortran (Debian 14.2.0-19) 14.2.0
Python 3.13.5
NumPy 2.3.5
SciPy 1.17.0
mpmath 1.3.0
```

`fpm` and `fprettify` were unavailable in this environment. All 202 library sources
compiled out of tree with:

```text
-std=f2018 -O0 -g -fcheck=all -fbacktrace
```

The focused `test_special_modified_fresnel` program passed with those runtime checks
and also with `-ffpe-trap=invalid,zero,overflow`. The complete manually linked ordinary
regression suite passed **84/84** test programs.

A compiled 4,001-point sweep over `[-20,20]` in steps of `0.01` was compared with
SciPy 1.17.0 on the nonnegative branch. Maximum errors scaled by
`max(1,abs(reference))` were about `9.65e-15` for `F_+`/`F_-` and `5.47e-15` for
`K_+`/`K_-`. An independent 201-point subset spanning both signs was compared with
80-digit mpmath values formed from the error-function integral identity; maximum scaled
errors were about `3.22e-16` for `F` and `1.69e-16` for `K`. Additional exact-binary64
phase cases at `1e3`, `1000000.125`, `100000000.25`, and `10000000000.125` agreed
with the 80-digit references to scaled errors below `1.1e-19`.

The freshly extracted delivery tree also passed:

```text
python3 tools/check_ascii.py
python3 tools/generate_special_reference.py
```

That fresh extraction compiled all 202 library sources and passed the complete
**84/84** ordinary regression suite plus the focused strict-trap test.

On systems with FPM and fprettify installed, run at minimum:

```text
fprettify -r src test
fpm test
fpm test --flag "-O0 -g -fcheck=all -fbacktrace -Wall -Wextra -std=f2018"
```

## Complex gamma-family validation (2026-09-25)

This batch extends the existing pure elemental `gamma`, `rgamma`, and `loggamma`
interfaces to complex binary64 arguments without changing the real specifics.
Complex `loggamma` is original SciFort code based on NIST DLMF 5.5.3 recurrence
and reflection plus the 5.11.1 Stirling series. The implementation tracks the
principal log-gamma sheet explicitly so the only branch cut is the negative real
axis, as documented by SciPy. Complex `gamma` and `rgamma` exponentiate that
branch with explicit overflow handling; reciprocal Gamma returns zero at its
nonpositive-integer zeros. No third-party implementation source or new license
was added.

The validation environment was:

```text
GNU Fortran (Debian 14.2.0-19) 14.2.0
Python 3.13.5
NumPy 2.3.5
SciPy 1.17.0
```

`fpm` and `fprettify` were unavailable in this environment. Direct GNU Fortran
validation compiled all 202 library source files with:

```text
-std=f2018 -O0 -g -fcheck=all -Wall -Wextra
```

The focused complex-gamma regression was also run with:

```text
-ffpe-trap=invalid,zero,overflow
```

and passed. It covers ordinary values, upper/lower limits of the negative-real
branch cut (including signed zero), recurrence, conjugation, Gamma/reciprocal-Gamma
identities, nonpositive-integer poles, elemental arrays, near-pole behavior, and
overflowing complex magnitudes without triggering the trapped overflow exception.

A compiled 6,540-point grid on `[-20,20] x [-20,20]`, excluding exact poles,
was compared with SciPy 1.17.0 public outputs. Maximum scaled discrepancies were
approximately `5.13e-15` for `loggamma`, `1.44e-14` for `gamma`, and `4.55e-14`
for `rgamma`. Additional large points such as `-100+50i` and `100+200i` are
covered by the focused regression.

The complete manually linked repository suite contained 85 program tests and 42
test support/reference sources. All 85/85 program tests passed against the same
202-source library build. Repository ASCII and special-reference-generator checks
were also run after the final documentation update.

On a system with FPM and fprettify installed, run at minimum:

```text
fprettify -r src test
fpm test
fpm test --flag "-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Wall -Wextra -std=f2018"
```
## Complex digamma/psi validation (2026-09-25)

This batch extends the top-level pure elemental `digamma` and `psi` interfaces to
complex binary64 arguments without changing the existing real specifics. The
complex implementation is original SciFort code based on NIST DLMF 5.5.2
recurrence, 5.5.4 reflection, and 5.11.2 asymptotic expansion. The reflection
term uses a scaled complex `cot(pi*z)` evaluation so offsets as small as `1e-300`
from a negative-integer pole remain finite when mathematically representable and
do not trigger avoidable intermediate divide-by-zero/overflow exceptions.

Validation environment:

```text
GNU Fortran (Debian 14.2.0-19) 14.2.0
Python 3.13.5
NumPy 2.3.5
SciPy 1.17.0
```

`fpm` and `fprettify` were unavailable in this environment. All 203 library
sources compiled out of tree with GNU Fortran using:

```text
-std=f2018 -O0 -g -fcheck=all -Wall -Wextra
```

The focused `test_special_complex_digamma` regression passed normally and with:

```text
-ffpe-trap=invalid,zero,overflow
```

Coverage includes ordinary complex values, recurrence/reflection/conjugation
identities, complex pole conventions, exact real-axis dispatch, elemental
arrays, large complex arguments, and offsets of `1e-8` and `1e-300` from a
negative-integer pole.

A compiled 13,020-point finite grid on `[-20,20] x [-10,10]`, excluding exact
poles, was compared with SciPy 1.17.0. The maximum error scaled by
`max(1,abs(reference))` was approximately `2.22e-15`. Additional probes near the
positive and first two negative real digamma zeros remained within about
`9e-16` absolute error.

The complete manually linked ordinary regression suite contained 86 program
tests and 42 support/reference source files. All **86/86** program tests passed
against the same library build.

On systems with FPM and fprettify installed, run at minimum:

```text
fprettify -r src test
fpm test
fpm test --flag "-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Wall -Wextra -std=f2018"
```

## Complex exponential-integral validation (2026-09-25)

This batch extends the top-level pure elemental `exp1` and `expi` interfaces to
complex binary64 arguments while preserving the existing real specifics. The
complex implementation is original SciFort code based on NIST DLMF 6.2.1--6.2.7
and the E1 continued fraction in DLMF 6.9.1. Small and near-cut arguments use the
convergent entire Ein series; the complementary region uses a modified-Lentz
continued fraction with guarded complex exponential scaling. Exact real-axis
arguments dispatch to the existing real kernels, with explicit signed-zero logic
for SciPy-compatible branch-cut behavior.

Validation environment:

```text
GNU Fortran (Debian 14.2.0-19) 14.2.0
Python 3.13.5
NumPy 2.3.5
SciPy 1.17.0
```

`fpm` and `fprettify` were unavailable in this environment. All 204 library
sources compiled out of tree with GNU Fortran using:

```text
-std=f2018 -O0 -g -fcheck=all -Wall -Wextra
```

The focused `test_special_complex_exponential_integrals` regression passed
normally and with:

```text
-ffpe-trap=invalid,zero,overflow
```

Coverage includes ordinary complex values, exact upper/lower values on the E1
negative-real branch cut (including signed zero), SciPy's real-valued exact-cut
Ei convention, near-cut values, conjugation, the Ei/E1 interrelation, exact
real-axis dispatch, elemental arrays, NaN handling, and a large finite
left-half-plane E1 value.

A compiled 25,840-point grid on `[-40,40] x [-40,40]` with spacing 0.5,
excluding the E1 cut itself, was compared with SciPy 1.17.0. Maximum errors
scaled by `max(1,abs(reference))` were approximately:

```text
complex exp1  1.69e-14
complex expi  1.32e-14
```

The complete manually linked ordinary regression suite contained 87 program
tests and 42 support/reference source files. All **87/87** program tests passed
against the same library build. Repository checks also passed:

```text
python tools/check_ascii.py
python tools/generate_special_reference.py
```

The updated special-reference generator produced 551 lines.

On systems with FPM and fprettify installed, run at minimum:

```text
fprettify -r src test
fpm test
fpm test --flag "-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Wall -Wextra -std=f2018"
```

## Complex sine/hyperbolic-integral validation (2026-09-25)

This batch extends the top-level pure elemental `sici` and `shichi` subroutines
to complex binary64 arguments while preserving the existing real specifics.
The complex implementation is original SciFort code based on the defining
series and interrelations in NIST DLMF 6.2. Small arguments use the convergent
Si/Ci/Shi/Chi series directly; larger `sici` arguments reuse the validated
principal E1 implementation and larger `shichi` arguments reuse principal Ei.
Complex Si and Shi are entire. Complex Ci and Chi use the principal logarithm
and therefore have the standard negative-real branch cut.

Validation environment:

```text
GNU Fortran (Debian 14.2.0-19) 14.2.0
Python 3.13.5
NumPy 2.3.5
SciPy 1.17.0
```

`fpm` and `fprettify` were unavailable in this environment. All 205 library
sources compiled out of tree with GNU Fortran using:

```text
-std=f2018 -O0 -g -fcheck=all -Wall -Wextra
```

The focused `test_special_complex_sine_integrals` regression passed normally
and with:

```text
-ffpe-trap=invalid,zero,overflow
```

Coverage includes ordinary complex values, tiny series arguments, large complex
values, exact negative-real cut values (including both signed-zero inputs),
upper/lower near-cut continuation, exact real-axis dispatch, imaginary-axis
identities, E1/Ei interrelations, elemental arrays, zero, and NaN handling.

A compiled 6,560-point finite grid on `[-40,40] x [-40,40]` with unit spacing,
excluding only the origin, was compared with SciPy 1.17.0. Maximum errors scaled
by `max(1,abs(reference))` were approximately:

```text
complex Si   5.93e-15
complex Ci   8.61e-15
complex Shi  5.93e-15
complex Chi  8.47e-15
```

The dense sweep exposed and fixed a signed-zero canonicalization issue on the
negative imaginary axis before finalization. Negating an exact `-i*y` argument
can create a `-0` real component; the canonical right-half-plane E1 relation now
normalizes that component to `+0` so the inherited E1 cut is approached from the
correct side.

The complete manually linked ordinary regression suite contained 88 program
tests and 42 support/reference source files. All **88/88** program tests passed
against the final source tree.

Repository checks to execute for the final delivery are:

```text
python tools/check_ascii.py
python tools/generate_special_reference.py
```

On systems with FPM and fprettify installed, run at minimum:

```text
fprettify -r src test
fpm test
fpm test --flag "-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Wall -Wextra -std=f2018"
```

## Complex Spence validation (2026-09-25)

This batch extends the top-level pure elemental `spence` interface to complex
binary64 arguments while preserving the existing real specific. The complex
implementation is original SciFort code for the principal `Li_2(1-z)` branch,
based on NIST DLMF 25.12.1, 25.12.4, and the logarithmic expansion associated
with 25.12.12. Direct power series are used inside `abs(w) <= 0.85`, DLMF
inversion outside the reciprocal radius, and a 24-term logarithmic/Bernoulli
series in the narrow annulus around the unit circle. Exact negative-real
complex inputs preserve signed-zero side selection on the Spence branch cut.

Validation environment:

```text
GNU Fortran (Debian 14.2.0-19) 14.2.0
Python 3.13.5
NumPy 2.3.5
SciPy 1.17.0
mpmath 1.3.0
```

`fpm` and `fprettify` were unavailable in this environment. The two changed
library sources, `src/special/spence.f90` and `src/special.f90`, were rebuilt
against the previously validated 205-source complex-Si/Ci/Shi/Chi baseline
using:

```text
-std=f2018 -O0 -g -fcheck=all -fbacktrace -Wall -Wextra
```

The focused `test_special_complex_spence` regression passed normally and with:

```text
-ffpe-trap=invalid,zero,overflow
```

Coverage includes SciPy documentation examples, exact upper/lower cut values
including signed zero, near-cut continuation, positive-real dispatch to the
existing real kernel, the dilogarithm reflection identity, elemental arrays,
nonfinite inputs, and a `1e300 + 1e300*i` finite stress case.

A compiled 25,921-point grid on `[-20,20] x [-20,20]` with spacing 0.25 was
compared with SciPy 1.17.0. The maximum scaled error over the 25,760 points with
nonzero imaginary part was approximately `9.97e-16`. Including the exact real
axis, the maximum was approximately `9.99e-15`, from the unchanged real
`spence(4.75)` path. Independent 100-digit mpmath checks at near-cut, unit-circle,
`1e100`, and `1e300` complex arguments agreed at binary64 precision.

The complete manually linked ordinary regression suite contained 89 program
tests and 42 support/reference source files. All **89/89** program tests passed
against the updated library objects. Repository checks also passed:

```text
python tools/check_ascii.py
python tools/generate_special_reference.py
```

The updated special-reference generator produced 584 lines.

On systems with FPM and fprettify installed, run at minimum:

```text
fprettify -r src test
fpm test
fpm test --flag "-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Wall -Wextra -std=f2018"
```

## Complex Airy `airy` and `airye`

Validation date: 2026-09-26. GNU Fortran 14.2.0 on Linux with SciPy 1.17.0.
The real specifics are unchanged. The complex `airy` and `airye` specifics were
compiled with `-O0 -g -fcheck=all -fbacktrace -Wall -Wextra -std=f2018`; focused
real/complex Airy regressions passed, and the complex `airy`/`airye` tests also
passed when the driver was built with `-ffpe-trap=invalid,zero,overflow`.

For complex `airye`, direct scaling of the moderate Airy power series is used
unless Ai/Ai' are recessive enough for cancellation to be amplified. In that
sector the implementation starts from a scaled DLMF 9.7 asymptotic value at a
radial anchor and continues Ai/Ai' inward with a 24-term Taylor recurrence
derived from `y''=z*y`. Outside the moderate disk, scaled asymptotic Ai values
and DLMF 9.2 rotated-Ai connections are combined through exponent differences,
so the routine does not first form exponentially large unscaled intermediates.
The real-specific negative-axis convention remains unchanged, while a complex
negative-real argument follows the principal complex scaling branch.

A compiled 2401-point grid on `[-12,12] x [-12,12]` had maximum componentwise
scaled errors of about `2.8e-12`, `7.4e-12`, `3.0e-12`, and `7.9e-12` for
`eAi`, `eAi'`, `eBi`, and `eBi'`, respectively, against SciPy 1.17.0. The
largest Ai/Ai' discrepancy occurs near the same Stokes-sector transition around
`-4-5.5i` seen by the unscaled complex implementation. A `-200+300i` regression
also remains finite and agrees with SciPy to the test tolerance, exercising the
scaled outer-plane path where unscaled intermediate values would be impractical.

`fpm` and `fprettify` were unavailable in this validation environment. As a
substitute for `fpm test`, all 205 library source files and 42 test support files
were compiled manually with GNU Fortran 14.2.0, then all **91/91** test programs
were linked and run successfully. `python tools/check_ascii.py` also passed, and
`python tools/generate_special_reference.py` completed successfully with SciPy
1.17.0.


## Complex Carlson symmetric elliptic validation (2026-09-26)

This batch extends `elliprc`, `elliprd`, `elliprf`, `elliprg`, and `elliprj`
to complex binary64 arguments while preserving every existing real specific.
Complex RC/RF/RD use principal-root Carlson duplication formulas after
homogeneity scaling. Complex RG uses the RF/RD identity with separate principal
square-root factors so product simplification cannot change the branch. Complex
RJ uses the defining Carlson integral after `t=(u/(1-u))**2` and adaptive
16-point Gauss-Legendre integration; this path was chosen after a direct
complexification of the legacy real RJ recurrence showed incorrect branch
selection on unbalanced complex inputs.

Validation environment:

```text
GNU Fortran (Debian 14.2.0-19) 14.2.0
Python 3.13.5
NumPy 2.3.5
SciPy 1.17.0
```

`fpm` and `fprettify` were unavailable in this environment. All 205 library
sources were manually compiled with GNU Fortran using a checked build based on:

```text
-std=f2018 -O0 -g -fcheck=all -Wall -Wextra
```

The focused `test_special_complex_carlson_elliptic` regression passed both the
checked build and a Carlson rebuild with:

```text
-ffpe-trap=invalid,zero,overflow -fbacktrace
```

Coverage includes ordinary complex values, a conjugate-pair domain, exact and
near negative-real cuts, RC Cauchy principal values with a complex first
argument, a deliberately unbalanced RJ case, complex-pole RJ with real
symmetric arguments, real-axis RJ principal-value dispatch, homogeneity,
elemental arrays, NaN/infinity conventions, and insufficient-domain rejection.

A compiled-Fortran comparison used 400 deterministic random cases with
independent magnitudes from `1e-6` through `1e6` and arguments spanning the
open right half-plane. Against SciPy 1.17.0, the maximum errors scaled by
`max(1,abs(reference))` were approximately `4.8e-16` for RC, `5.8e-16` for
RF, `6.8e-16` for RD, `4.7e-15` for RG, and `7.6e-15` for RJ. All 400 SciPy
references were finite. Additional targeted SciPy references exercise the
negative-real RC second argument, the RJ conjugate-pair case, a complex RJ
pole with negative real part when the symmetric arguments are nonnegative
real, and the intentionally rejected generic RJ domain.

The final optimized manual regression build linked and ran all **92/92** test
programs successfully, including the C ABI and all pre-existing special,
distribution, statistics, and resampling tests. Repository checks also passed:

```text
python tools/check_ascii.py
python tools/generate_special_reference.py
```

The updated special-reference generator produced 598 lines. On systems with
FPM and fprettify installed, run at minimum:

```text
fprettify -r src test
fpm test
fpm test --flag "-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -Wall -Wextra -std=f2018"
```

## Generalized hypergeometric `hyp0f1` validation (2026-09-26)

This batch adds the pure elemental `hyp0f1(v,z)` interface for a real denominator
parameter and real or complex binary64 arguments. The origin disk uses the
defining 0F1 series with compensated accumulation. Larger real arguments use the
DLMF Bessel connection and SciFort's existing real J/I kernels. Larger complex
arguments start from the converged origin series and continue radially outward
with local Taylor expansions derived from the 0F1 differential equation.

Validation environment:

```text
GNU Fortran (Debian 14.2.0-19) 14.2.0
Python 3.13.5
NumPy 2.3.5
SciPy 1.17.0
```

`fpm` was unavailable in this environment. All 206 library source files compiled
successfully with a checked GNU Fortran build based on:

```text
-std=f2018 -O0 -g -fcheck=all -Wall -Wextra
```

The focused `test_special_hyp0f1` regression passed normally and after recompiling
the new numerical module with:

```text
-ffpe-trap=invalid,zero,overflow -fbacktrace
```

The focused cases cover the defining-series region, large positive real growth,
large negative-real oscillation, negative noninteger parameters, complex
continuation, small oscillatory complex results, elemental arrays, parameter
poles, NaNs, and infinities. A deterministic 180-case complex sweep used
`v = -10.5, -3.25, -0.5, 0.1, 0.5, 1, 2.5, 10, 25` and magnitudes from `1e-2`
through roughly `6.3e2`; the maximum error scaled by `max(1,abs(reference))`
against SciPy 1.17.0 was approximately `5.9e-13`. All **32/32** existing and new
special-function test programs were then compiled and run successfully against the
checked library build. `python tools/check_ascii.py` and `python tools/generate_special_reference.py` also passed.

The complex continuation has not yet been exhaustively characterized for very
large complex arguments combined with large negative denominator parameters;
that regime remains a future accuracy target rather than a claimed compatibility
range.



## Confluent hypergeometric `hyp1f1` validation (2026-09-26)

This batch adds pure elemental `hyp1f1(a,b,z)` overloads for real parameters and real or complex binary64 arguments. The implementation uses the defining Kummer series in the origin region. Moderate arguments use Taylor continuation with coefficients derived from `z*y'' + (b-z)*y' - a*y = 0`; left-half-plane continuation uses Kummer's transformation to reduce cancellation. Large positive/negative real arguments use separate log-domain DLMF 13.7 expansions, while large nonreal arguments use both asymptotic branches with the phase selected by the sign of the imaginary part. Negative-integer numerator parameters are evaluated as terminating polynomials.

The checked GNU Fortran build compiled all **206** library source files successfully with `-std=f2018 -O0 -fcheck=all -Wall -Wextra`. The focused `test_special_hyp1f1` regression passed normally and with `-ffpe-trap=invalid,zero,overflow`. All **33/33** special-function test programs then compiled and ran successfully.

A deterministic validation sweep compared 160 finite real cases and 160 finite complex cases with SciPy 1.17.0. Complex cases were retained only when the SciPy value also agreed with an 80/40-digit mpmath evaluation within `2e-10` scaled error, because SciPy 1.17.0's legacy complex 1F1 backend has isolated discrepancies. The SciFort maximum error scaled by `max(1,abs(reference))` was about **2.9e-13** for the real cases and **1.8e-10** for the complex cases. The sweep covered real parameters approximately from -8 to 12, real arguments from -150 to 150 plus a targeted value at 700, and complex components out to roughly +/-120.

Focused references additionally cover a removable real denominator pole, SciPy's complex denominator-pole convention, terminating polynomials, the `a=b` exponential identity, NaNs, positive-infinite `b`, elemental arrays, and large complex arguments on both sides of the imaginary axis. Extreme simultaneous large parameters and arguments, especially close to Stokes transitions or gamma poles, remain a future accuracy target.

## Gauss hypergeometric `hyp2f1` validation (2026-09-26)

This batch adds pure elemental `hyp2f1(a,b,c,z)` overloads for real parameters and real or complex binary64 arguments. The implementation evaluates the defining Gauss series in a small origin disk and then analytically continues the solution with local Taylor expansions. The coefficient recurrence is derived directly from the Gauss differential equation `z*(1-z)*y'' + (c-(a+b+1)*z)*y' - a*b*y = 0`; continuation steps are limited by the distance to the regular singular points at zero and one. Complex continuation uses an upper- or lower-half-plane detour around `z=1`, with signed zero selecting the branch-cut side. Exact `z=1`, terminating numerator polynomials, and exact `a=c`/`b=c` cancellations are handled separately.

GNU Fortran 14.2.0 compiled all **206** library source files successfully with `-std=f2018 -O0 -fcheck=all -Wall -Wextra`. The focused `test_special_hyp2f1` regression passed both normally and with `-ffpe-trap=invalid,zero,overflow`. All **34/34** `test_special*.f90` programs then compiled and ran successfully against the checked build, and `python tools/check_ascii.py` passed. FPM was not available in the execution environment, so `fpm test` could not be run here.

Regression references cover ordinary real values, continuation close to `z=1`, large negative real arguments, a finite real case beyond the branch point, Gauss' value at one, terminating polynomials, denominator-pole cancellation, upper/lower signed-zero branch-cut values, elemental arrays, NaNs, and a difficult complex case for which SciPy 1.17.0 differs from an 80-digit mpmath value by about `1.6e-7`; SciFort agrees with the high-precision value to about `4e-13` absolute error. Nonfinite nonpolynomial `z` values are currently returned as NaN rather than using a dedicated asymptotic limit formula.


## Tricomi confluent hypergeometric `hyperu` validation (2026-09-26)

This batch adds pure elemental real `hyperu(a,b,x)` for finite real parameters and the principal real branch `x >= 0`. The implementation anchors the recessive solution with the DLMF large-positive-argument expansion and continues it inward with the same local Taylor recurrence derived from `x*y'' + (b-x)*y' - a*y = 0` that underlies the 1F1 continuation code. Small positive arguments use the exact 1F1 connection formula in a log-gamma form; near integer `b`, symmetric evaluations at two parameter offsets are Richardson-extrapolated to avoid the removable gamma singularities of that formula. Negative-integer `a` is evaluated by the exact terminating generalized-Laguerre recurrence.

GNU Fortran 14.2.0 compiled all **206** library source files successfully with `-std=f2018 -O0 -fcheck=all -Wall -Wextra`. The focused `test_special_hyperu` regression passed both normally and with `-ffpe-trap=invalid,zero,overflow`. A deterministic 280-case sweep against 70-digit mpmath values covered `a` approximately from -10 to 12, `b` from -8 to 12, positive `x` from `1e-6` through about 200, and additional exact/near-integer-`b` small-argument cases. The maximum error scaled by `max(1,abs(reference))` was about **5.3e-10**; no case exceeded `5e-8`. In the same sweep, finite SciPy 1.17.0 values differed from mpmath by as much as about `2.8e-4`, so the high-precision mathematical values were used as the primary validation target.

Focused references cover ordinary values, small positive arguments with integer and noninteger `b`, near-integer regularization, large positive arguments, a difficult large-parameter continuation case, negative noninteger `a`, terminating negative-integer `a`, Kummer's transformation, elemental arrays, negative-axis NaNs, infinite arguments, and `x=0` limits. SciFort returns the mathematically signed one-sided divergence at `x=0`; SciPy 1.17.0 returns positive infinity for some negative-limit cases. Extremely large simultaneous parameters remain a future accuracy and overflow target.


## Parabolic-cylinder `pbdv` validation (2026-09-26)

This batch adds pure elemental real `pbdv(v,x,d,dprime)` with SciPy-compatible value/derivative output semantics. Positive arguments reuse `hyperu` through the DLMF relation between `D_v` and Kummer U. Negative arguments start from the exact gamma-ratio values `D_v(0)` and `D_v'(0)` and use local Taylor continuation of `y'' + (v + 1/2 - x**2/4)*y = 0`. Exact nonnegative integer orders on the negative axis use the Hermite three-term recurrence, because direct ODE propagation can otherwise amplify the mathematically absent dominant solution.

GNU Fortran 14.2.0 compiled all **207** library source files successfully with `-std=f2018 -O0 -fcheck=all -Wall -Wextra`. The focused `test_special_parabolic_cylinder` regression passed normally and with `-ffpe-trap=invalid,zero,overflow`. A deterministic 250-case sweep against 70-digit mpmath values covered orders approximately from -15 to 25, arguments from -15 to 20, and additional exact nonnegative-integer orders on the negative axis. Across both values and first derivatives, the maximum error scaled by `max(1,abs(reference))` was about **5.8e-13**.

Focused regression values include high-order negative-argument cases such as `v=10.2, x=-10` and `v=20.2, x=-10`, where SciPy 1.17.0's legacy backend differs substantially from high-precision mpmath; SciFort follows the high-precision mathematical values. Coverage also includes ordinary positive/negative arguments, the `hyperu` positive branch, exact Hermite orders, large positive arguments, origin values, elemental arrays, NaNs, and infinite arguments. After the new routines were integrated, all **36/36** `test_special*.f90` programs and all **97/97** program-bearing test files in `test/` compiled and ran successfully against the checked library build. `python tools/check_ascii.py` also passed. FPM was not available in the execution environment, so `fpm test` could not be run here. Extremely large order/positive-argument combinations that require a fully scaled Kummer-U intermediate remain a future target.

## Parabolic-cylinder `pbvv` and sequence validation (2026-09-26)

This batch extends the parabolic-cylinder module with real `pbvv(v,x,vv,vprime)`, real-order `pbdv_seq`/`pbvv_seq`, and analytic complex nonnegative-integer `pbdn_seq`. `pbvv` starts from the exact DLMF reciprocal-gamma values at zero, advances moderate positive arguments by the parabolic-cylinder ODE, switches to the DLMF growing-solution expansion at large positive arguments, and evaluates negative arguments through the exact D/V connection formula. Reduced `sin(pi*v)` and `cos(pi*v)` helpers preserve exact integer and half-integer cancellations. The sequence wrappers preserve SciPy's historical minimum-two-output convention.

A deterministic 63-point high-precision sweep covering seven noninteger orders from -7.3 through 10.2 and nine arguments from -10 through 15 compared both `pbvv` values and first derivatives with 50-digit mpmath values. The derivative references used the exact adjacent-order identity rather than numerical differencing. Across 126 value/derivative comparisons, the largest error scaled by `max(1,abs(reference))` was about **6.4e-12**, at `v=10.2, x=-10`. Focused tests also cover the exact origin formulas, the large-positive asymptotic path, integer-order negative arguments, NaNs, positive/negative sequence conventions, optional status values, and complex `pbdn_seq` values/derivatives against 80-digit mpmath references. In particular, `pbvv(1,-8)` is checked against the high-precision value about `1.1643327222791885e5`; SciPy 1.17.0's legacy `pbvv` backend returns a spurious result of order `1e261` for this case, so SciFort follows the high-precision mathematical value. Complex `pbdn_seq` likewise follows the analytic DLMF `D_n(z)` definition rather than reproducing SciPy 1.17.0's nonanalytic legacy complex behavior.


GNU Fortran 14.2.0 compiled all **207** library source files in dependency order with `-std=f2018 -O0 -fcheck=all`; all **36/36** `test_special*.f90` programs passed. The complete test tree contains **97** program-bearing test files and all **97/97** compiled and ran successfully against the checked library build. The terminal-oriented `test_multiple_comparisons` test was rerun with `TERM=xterm` because the parallel noninteractive harness did not set `TERM`; it then passed normally. `python tools/check_ascii.py` and the 132-character Fortran source-line audit also pass. FPM is not installed in this execution environment, so `fpm test` could not be run here.

## Scalar special helpers validation (2026-09-26)

This batch adds `agm`, `huber`, `pseudo_huber`, `ndtri_exp`, `owens_t`, and
`voigt_profile`, and moves skew-normal CDF evaluation onto the new public Owen-T
kernel. It also hardens `normal_logcdf` for finite negative arguments whose
logarithmic probability lies below binary64 range.

GNU Fortran 14.2.0 compiled all **208** library source files in dependency order
with `-std=f2018 -O0 -g -fcheck=all -ffpe-trap=zero,overflow -fbacktrace`. All
**37/37** `test_special*.f90` programs passed. The complete test tree now
contains **98** program-bearing files; all **98/98** compiled and ran against the
checked build. The terminal-oriented multiple-comparisons test was run with
`TERM=xterm`. FPM is not installed in this execution environment, so `fpm test`
could not be run here.

A deterministic 121-point `owens_t` grid used eleven `h` values from -20 through
20 and eleven slopes from -1000 through 1000. Against SciPy 1.17.0, the largest
relative difference for a nonzero reference was about **2.8e-14**. A 112-point
positive-parameter `voigt_profile` grid used seven `x` values from -10 through
10, four Gaussian widths from 0.1 through 2, and four Cauchy widths from zero
through 3. The largest relative difference from SciPy 1.17.0 was about
**1.1e-13**, occurring for a profile value of order `1e-195`.

For `ndtri_exp`, ordinary and subnormal-probability cases agree with SciPy to
roundoff. The deep-tail values were also solved independently with 80-digit
mpmath. At `y=-10000`, SciFort gives `-141.37983987312717`, agreeing with the
high-precision value `-141.3798398731271637...`; SciPy 1.17.0 differs by about
`2.6e-12`. At `y=-1e6`, SciFort gives `-1414.2077829910172`, agreeing with the
high-precision value `-1414.2077829910173269...`; SciPy 1.17.0 differs by about
`6.9e-10`. Focused tests also cover AGM sign/zero/infinity behavior, Huber and
pseudo-Huber branch semantics, Voigt Gaussian/Cauchy/degenerate limits, Owen-T
symmetries and infinite slope, elemental arrays, NaNs, and the extreme
`normal_logcdf(-1e200)` overflow-trap regression.


## Legacy distribution-special helper validation (2026-09-26)

This batch adds 21 SciPy-style probability helpers for binomial, negative-binomial, Poisson, chi-square, F, gamma, and Student-t laws. The implementations are thin numerical wrappers over SciFort's regularized incomplete beta/gamma and normal kernels; no new probability summation or integration algorithm is introduced. The F helpers form complementary beta coordinates in the log-ratio domain, and the Student-t helpers likewise avoid forming `t*t/df` when that intermediate can overflow.

GNU Fortran 14.2.0 compiled all **209** library source files in dependency order with `-std=f2018 -O0 -g -fcheck=all -Wall -Wextra -ffpe-trap=invalid,zero,overflow -fbacktrace`. All **38/38** `test_special*.f90` programs passed against that checked build. The complete test tree now contains **99** program-bearing files; all **99/99** compiled and ran successfully against the final checked library and test-support build (with `TERM=xterm` for the terminal-oriented multiple-comparisons test). A deterministic **426-case** comparison with SciPy 1.17.0 over ordinary and tail regimes had maximum error scaled by `max(1,abs(reference))` of about **1.65e-14**. The largest relative difference was about **1.33e-13**, for a gamma lower-tail probability near `5.21e-219`. Focused regressions also cover inverse round trips, support and parameter validation, simultaneous Poisson infinities, extreme scale products, asymmetric F tails, and deep Student-t tails.

SciFort intentionally differs from two SciPy 1.17.0 legacy endpoint results: `pdtri(k,1)` returns the exact inverse mean zero, and `stdtrit(df,0)` returns negative infinity. Both satisfy the documented inverse definitions. `python tools/check_ascii.py` and the 132-character Fortran line-length audit pass. FPM and `fprettify` are not installed in this execution environment, so `fpm test` and automatic formatting could not be run here.

## Legacy distribution parameter-inverse validation (2026-09-26)

This batch adds 13 SciPy-style inverse-parameter helpers: `bdtrik`, `bdtrin`,
`btdtria`, `btdtrib`, `chdtriv`, `fdtridfd`, `gdtrib`, `nbdtrik`, `nbdtrin`,
`pdtrik`, `stdtridf`, `nrdtrimn`, and `nrdtrisd`. The beta/gamma/Student-t
shape inversions use a safeguarded positive-parameter bracket expanded
multiplicatively and refined geometrically; the normal mean/standard-deviation
inverses use the exact `ndtri` transform. Count inverses use the continuous
incomplete-beta/gamma extension and return real values. `fdtridfd` scans from
the small-dfd branch and therefore returns the smaller positive root when the
F-CDF equation is not globally one-to-one.

GNU Fortran 14.2.0 compiled all **209** library source files from scratch in
dependency order with `-std=f2018 -O0 -g -fcheck=all
-ffpe-trap=invalid,zero,overflow -fbacktrace`. All **39/39**
`test_special*.f90` programs passed against that exact build with floating-point
traps enabled. The complete repository now contains **100** program-bearing test
files; all **100/100** compiled and passed with `-fcheck=all` and normal IEEE
exception handling (the general suite intentionally exercises NaN/invalid cases
that cannot be run with `-ffpe-trap=invalid`). A deterministic **423-case**
comparison against SciPy 1.17.0 public outputs over the 13 new interfaces had
maximum error scaled by `max(1,abs(reference))` of about **3.61e-11**. Focused
round-trip tests additionally verify each inverse against the underlying
incomplete-beta/gamma, F, Student-t, or normal CDF identity.

`python tools/check_ascii.py` and the 132-character Fortran line-length audit
pass. FPM and `fprettify` are not installed in this execution environment, so
`fpm test` and automatic formatting could not be run here.



## SciPy convenience-helper validation (2026-09-26)

This batch adds `cbrt`, `exp10`, `exp2`, `expm1`, `log1p`, `powm1`, `radian`,
`round`, normalized `sinc`, `sindg`, `cosdg`, `tandg`, `cotdg`, `softplus`, and
`diric`. The implementations are independent formulas using Fortran intrinsics
and existing SciFort stable elementary kernels; no SciPy implementation source
was copied or translated.

GNU Fortran 14.2.0 compiled all **209** library source files with
`-std=f2018 -O0 -g -fcheck=all`; the focused special-function build additionally
used `-ffpe-trap=invalid,zero,overflow`. All **40/40** `test_special*.f90`
programs passed with those traps enabled. The complete repository now contains
**101** program-bearing test files, and all **101/101** passed with `-fcheck=all`
and normal IEEE exception handling; `TERM=xterm` was supplied to the
terminal-oriented multiple-comparisons test.

A deterministic **1005-case** sweep compared the new routines with SciPy
1.17.0 public outputs. Excluding four deliberately different degree-trigonometric
values at the finite input `1.00000000000013e14`, where SciPy 1.17.0 returns
zero while SciFort performs modulo-360 reduction, the remaining **1001** direct
comparisons had maximum error scaled by `max(1,abs(reference))` of about
**8.98e-16**. The four large-angle SciFort values agree with direct reduction
of the represented binary64 degree input. Focused tests also cover overflow
guards, infinities and NaNs, exact quadrant values, ties-to-even rounding,
negative-base integer powers, cancellation near `x=1` in `powm1`, exact sinc
integer zeros, Dirichlet-kernel removable singularities, elemental arrays, and
tiny-angle `cotdg` reciprocal overflow.

`python tools/check_ascii.py` and the 132-character Fortran line-length audit
pass. FPM and `fprettify` are not installed in this execution environment, so
`fpm test` and automatic formatting could not be run here.

## Legacy noncentral-distribution special-helper validation (2026-09-26)

This batch adds the 13 legacy noncentral-distribution `scipy.special` interfaces
`chndtr`, `chndtridf`, `chndtrinc`, `chndtrix`, `ncfdtr`, `ncfdtri`,
`ncfdtridfd`, `ncfdtridfn`, `ncfdtrinc`, `nctdtr`, `nctdtridf`, `nctdtrinc`,
and `nctdtrit`. Forward probabilities and variate inverses reuse SciFort's
existing noncentral chi-square, noncentral F, and noncentral Student-t kernels;
parameter inverses use independently written safeguarded scalar searches.
Noncentral-F degree-of-freedom inverses scan upward from the small-positive
branch and therefore return the smallest positive solution when two roots are
present.

GNU Fortran 14.2.0 compiled all **210** library source files from scratch in
dependency order with `-std=f2018 -O0 -g -fcheck=all -fbacktrace`; the first
compile tranche also used `-Wall -Wextra`. All **41/41** `test_special*.f90`
programs passed against the exact rebuilt library with
`-ffpe-trap=invalid,zero,overflow`. The complete repository now contains
**102** program-bearing test files, and all **102/102** compiled and passed with
`-fcheck=all` and normal IEEE exception handling (`TERM=xterm` was supplied to
the terminal-oriented multiple-comparisons test).

A deterministic **575-case** comparison of the forward CDF and variate-inverse
interfaces with SciPy 1.17.0 public outputs covered noncentral chi-square and F
ordinary/tail regimes, noncentral-t CDFs including `df=0.5`, and noncentral-t
quantiles for `df >= 2` with probabilities from 0.01 through 0.99. The maximum
error scaled by `max(1,abs(reference))` was about **1.93e-12**, at an
`ncfdtri` value near 1024.44. A separate **84-case** parameter-inverse
round-trip sweep substituted each returned degree-of-freedom or noncentrality
parameter back into the corresponding SciFort CDF; the maximum absolute CDF
residual was about **6.66e-15**. Focused tests additionally exercise the
smallest-positive-root convention for both noncentral-F degree inverses.

The wrappers intentionally do not reproduce SciPy 1.17.0's observed
`nctdtrit(df,nc,0)=+infinity` legacy endpoint; SciFort returns mathematical
negative infinity. The wrapper also rejects infinite degrees of freedom to
match the documented `scipy.special` finite-parameter domain, even though
SciFort's higher-level noncentral-t distribution supports the shifted-normal
infinite-df limit. Very extreme low-df noncentral-t quantiles retain the
accuracy limits of the pre-existing `nct_ppf` kernel; this batch does not alter
that distribution algorithm.

`python tools/check_ascii.py` and the 132-character Fortran line-length audit
pass. FPM and `fprettify` are not installed in this execution environment, so
`fpm test` and automatic formatting could not be run here.

## Orthogonal-polynomial evaluator validation (2026-09-26)

This batch adds 15 nonnegative-integer-degree evaluator interfaces: `eval_legendre`,
`eval_jacobi`, `eval_gegenbauer`, `eval_laguerre`, `eval_genlaguerre`,
`eval_hermite`, `eval_hermitenorm`, `eval_chebyt`, `eval_chebyu`, `eval_chebyc`,
`eval_chebys`, `eval_sh_legendre`, `eval_sh_jacobi`, `eval_sh_chebyt`, and
`eval_sh_chebyu`. Real/complex polynomial evaluation uses independently written
three-term recurrences from DLMF 18.9; shifted families are standard affine
transformations. This checkpoint intentionally restricts Jacobi, Gegenbauer, and
generalized-Laguerre parameters to their classical orthogonality domains and does
not yet implement SciPy's noninteger-degree hypergeometric continuation.

GNU Fortran 14.2.0 compiled all **211** library sources from scratch with
`-std=f2018 -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow`. All **42/42**
`test_special*.f90` programs passed with those floating-point traps enabled. The
complete repository now contains **103** program-bearing tests; all **103/103**
compiled and passed with `-fcheck=all` and normal IEEE exception handling, with
`TERM=xterm` for the terminal-oriented multiple-comparisons test.

A deterministic **2,205-case** real sweep covered degrees 0 through 20, seven
arguments from -1.7 through 1.8, and all 15 new evaluator names. Against SciPy
1.17.0 public outputs, the maximum error scaled by `max(1,abs(reference))` was
about **1.72e-14**. A separate **845-case** complex sweep covered degrees 0
through 12 and five complex arguments for the 13 families with complex overloads;
against 70-digit mpmath values, the maximum scaled error was about **2.42e-15**.
SciPy 1.17.0's legacy complex Chebyshev path differs from high precision by up to
about `1.1e-8` absolute on that grid, while the direct SciFort recurrence remains
near binary64 roundoff. Focused tests also cover elemental arrays and invalid
degree/parameter domains.

FPM and `fprettify` are not installed in this execution environment, so `fpm test`
and automatic formatting could not be run here.

## Orthogonal-polynomial quadrature-root validation (2026-09-26)

This batch adds 15 SciPy-style Gaussian-quadrature interfaces: `roots_legendre`,
`roots_jacobi`, `roots_gegenbauer`, `roots_laguerre`, `roots_genlaguerre`,
`roots_hermite`, `roots_hermitenorm`, `roots_chebyt`, `roots_chebyu`,
`roots_chebyc`, `roots_chebys`, `roots_sh_legendre`, `roots_sh_jacobi`,
`roots_sh_chebyt`, and `roots_sh_chebyu`. General classical rules use the
Golub-Welsch Jacobi-matrix construction with the pre-existing SciFort symmetric
Jacobi eigensolver and exact zeroth moments; Chebyshev families use explicit
trigonometric nodes and weights. Symmetric rules average opposite node/weight
pairs so odd orders have an exact central zero.

GNU Fortran 14.2.0 compiled all **212** library sources from scratch with
`-std=f2018 -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow`. All **43/43**
`test_special*.f90` programs passed with those traps enabled. The complete
repository contains **104** program-bearing tests after this addition; all
**104/104** compiled and passed with `-fcheck=all` and normal IEEE exception
handling, with `TERM=xterm` for terminal-oriented tests.

A deterministic **5,757-value** comparison covered the 15 interfaces at orders
1, 2, 3, 5, 8, 12, 20, 32, and 64, including multiple Jacobi,
Gegenbauer, generalized-Laguerre, and shifted-Jacobi parameter sets. Against
SciPy 1.17.0 public outputs, the maximum error scaled by
`max(1,abs(reference))` was about **9.73e-12**, at an endpoint weight of the
64-point `roots_jacobi(64,-0.8,0.2)` rule. Restricting the sweep to order 32 or
below, the maximum scaled error was about **5.57e-13**. Focused tests cover all
15 APIs, exact/known zeroth moments, positive weights, ascending nodes, a
32-point asymmetric Jacobi rule, and invalid order/parameter domains.

The current general-family path uses a dense O(n^3) Jacobi eigensolver and is
therefore intended primarily for moderate quadrature orders. SciPy's specialized
large-order asymptotic Hermite algorithms are future work. FPM and `fprettify`
are not installed in this execution environment, so `fpm test` and automatic
formatting could not be run here.
## Historical statistical special-wrapper validation (2026-09-26)

This batch adds `smirnov`, `smirnovi`, `kolmogorov`, `kolmogi`, and
`tklmbda`. They are thin compatibility wrappers around the pre-existing
SciFort `ksone`, `kstwobign`, and Tukey-lambda distribution kernels, so no
independent probability approximation is duplicated in the special layer.

GNU Fortran 14.2.0 compiled all **213** library sources from scratch with
`-std=f2018 -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace`.
All **44/44** `test_special*.f90` programs passed against that exact build
with the floating-point traps enabled. The complete repository contains
**105** program-bearing tests after this addition; all **105/105** compiled
and passed with `-fcheck=all` and normal IEEE exception handling, with
`TERM=xterm` supplied to terminal-oriented tests.

A deterministic **170-case** comparison against SciPy 1.17.0 public outputs
covered finite-sample one-sided KS probabilities/inverses, limiting
two-sided KS probabilities/inverses, and Tukey-lambda CDF values across
negative, zero, and positive shape parameters. The maximum error scaled by
`max(1,abs(reference))` was about **2.40e-14**, at
`smirnov(100,0.03)`. Focused tests also cover support/probability endpoints,
invalid finite-sample arguments, inverse round trips, and Tukey-lambda
finite-support limits.

`python tools/check_ascii.py` and the 132-character Fortran line-length audit
pass. FPM and `fprettify` are not installed in this execution environment, so
`fpm test` and automatic formatting could not be run here.

## Legacy orthogonal compatibility validation (2026-09-26)

This batch adds SciPy's historical orthogonal compatibility names without a
second numerical implementation. `assoc_laguerre(x,n,k)` reverses the argument
order of `eval_genlaguerre(n,k,x)` and defaults `k` to zero. The 15 historical
quadrature names `p_roots`, `j_roots`, `cg_roots`, `l_roots`, `la_roots`,
`h_roots`, `he_roots`, `t_roots`, `u_roots`, `c_roots`, `s_roots`, `ps_roots`,
`js_roots`, `ts_roots`, and `us_roots` are generic aliases of the canonical
`roots_*` procedures. Focused tests verify bit-identical node, weight, zeroth-
moment, and status results for every alias.

A deterministic **468-case** real/complex `assoc_laguerre` sweep covered degrees
0 through 12, four generalized-Laguerre orders, six real arguments, and four
complex arguments. Against SciPy 1.17.0 public outputs, the maximum error scaled
by `max(1,abs(reference))` was about **3.64e-14**, at a degree-12 complex value.
The compatibility wrapper itself introduces no extra numerical error because it
directly calls the previously validated SciFort generalized-Laguerre recurrence.

GNU Fortran 14.2.0 compiled all **213/213** library sources with
`-std=f2008 -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow`. All **45/45**
`test_special*.f90` programs passed with those floating-point traps enabled.
The complete repository contains **106** program-bearing tests after this
addition; all **106/106** compiled and passed with `-fcheck=all` and normal IEEE
exception handling, with `TERM=xterm` supplied to the terminal-oriented
multiple-comparisons test.

`python tools/check_ascii.py` and the 132-character Fortran line-length audit
pass. FPM and `fprettify` are not installed in this execution environment, so
`fpm test` and automatic formatting could not be run here.

## Lambert W and Wright omega validation (2026-09-27)

This batch adds SciPy-style `lambertw` for real or complex arguments on arbitrary integer branches and real/complex `wrightomega`. The implementation is independent: initial values come from the DLMF square-root expansion around `-1/e` and the logarithmic large-argument expansion, followed by Halley refinement of the defining equations. Wright omega uses the branch-selected identity with Lambert W while `exp(z)` is numerically safe and switches to direct refinement of `w + log(w) = z` outside that range.

GNU Fortran 14.2.0 compiled all **214/214** library sources from scratch with `-std=f2018 -O0 -fcheck=all -Wall -Wextra -ffpe-trap=invalid,zero,overflow`. All **46/46** `test_special*.f90` programs passed with invalid/zero/overflow traps enabled. The complete repository contains **107** program-bearing tests after this addition; all **107/107** compiled and passed with `-fcheck=all` and normal IEEE exception handling, with `TERM=xterm` available to terminal-oriented tests.

A deterministic **2,108-case** Lambert-W sweep covered branches `k=-8,...,8`, magnitudes from about `1e-10` through `1e10`, complex quadrants, negative-real branch cuts, and points near `-1/e`. Against SciPy 1.17.0 public outputs, the maximum error scaled by `max(1,abs(reference))` was about **2.09e-16**. A separate **1,477-case** Wright-omega sweep covered real parts from `-1000` through `1000` and multiple unwinding strips; the maximum scaled difference from SciPy 1.17.0 was about **1.66e-12**, at a point only `1e-8` from the ill-conditioned `-1-i*pi` branch point. Away from that branch-point neighborhood the observed errors were near machine precision.

Focused tests additionally verify the defining identities `W*exp(W)=z` and `omega+log(omega)=z`, arbitrary positive and negative branch indices, real `k=-1`, exact `-1/e`, exact Wright-omega branch points `-1 +/- i*pi`, very large positive and noncentral negative Wright-omega arguments, nonprincipal zero behavior, and invalid convergence tolerances. SciFort intentionally returns the mathematical value `W(-1/e)=-1`, rather than reproducing the NaN observed from SciPy 1.17.0 at that exact point.

`python tools/check_ascii.py` and the 132-character Fortran line-length audit pass. FPM and `fprettify` are not installed in this execution environment, so `fpm test` and automatic formatting could not be run here.

## Associated Legendre and spherical-harmonic validation (2026-09-27)

The first associated-Legendre batch added real-degree/integer-order Ferrers `lpmv` together with normalized `sph_legendre_p` and complex `sph_harm_y`. The implementation is independent and follows DLMF hypergeometric formulas for noninteger degree, recurrence for integer degree, and directly normalized spherical recurrences. Positive noninteger orders use a differentiated hypergeometric representation to avoid large gamma-ratio cancellation near `x=-1`.

This extension adds integer-degree real/complex `assoc_legendre_p` with branch-cut conventions 2 and 3, optional spherical normalization, and Fortran derivative helpers for associated Legendre, spherical Legendre, and spherical harmonics. Branch 2 is constructed with `sqrt(1-z)*sqrt(1+z)` and branch 3 with `sqrt(z-1)*sqrt(z+1)`, retaining the principal-root continuation on each cut. Derivatives are generated by propagating truncated Taylor jets through the same recurrence rather than by finite differences.

GNU Fortran 14.2.0 compiled all **215/215** library sources from scratch with `-std=f2018 -O0 -fcheck=all -Wall -Wextra -ffpe-trap=invalid,zero,overflow`. All **49/49** `test_special*.f90` programs passed with invalid/zero/overflow traps enabled. The complete repository contains **110** program-bearing tests after the Legendre-Q addition; all **110/110** compiled and passed with `-fcheck=all` and normal IEEE exception handling, with `TERM=xterm` supplied where appropriate.

The original deterministic **616-case** `lpmv` sweep had maximum scaled error about **8.29e-13** against 70-digit mpmath. The earlier **1,105** `sph_legendre_p` and **3,315** `sph_harm_y` value comparisons against SciPy 1.17.0 had maximum scaled error about **3.46e-15**.

For the extension, **3,380** complex `assoc_legendre_p` comparisons covered degrees through 12, all orders `-n,...,n`, both branch conventions, normalized and unnormalized forms, and five nonreal arguments. The maximum error scaled by `max(1,abs(reference))` versus SciPy 1.17.0 was about **4.75e-15**. A separate **720-value** first/second derivative sweep had maximum scaled error about **2.66e-15**. Another **507** normalized spherical-Legendre derivative comparisons through degree 12 had maximum scaled error about **6.02e-15**, and mixed spherical-harmonic derivatives were likewise at a few `1e-15`. Exact third/fourth derivatives of the polynomial `P_4^2(z)` are included to exercise derivative orders beyond SciPy 1.17.0's currently implemented `diff_n <= 2` path.

The 2026-09-28 compatibility extension adds scalar real/complex `legendre_p` and the current first-kind table APIs `legendre_p_all`, `assoc_legendre_p_all`, `sph_legendre_p_all`, and `sph_harm_y_all`. A deterministic **5,205-entry** sweep compared values and first/second derivatives against SciPy 1.17.0 public outputs. It covered ordinary Legendre degrees through 12 at four real arguments; associated Legendre degrees through 8 and orders `-5:5` at four real and four complex arguments, including normalization; spherical Legendre degrees through 8 and orders `-5:5` at three angles; and all value/gradient/Hessian components of the corresponding spherical harmonics. With errors scaled by `max(1,abs(reference))`, the maxima were about **8.02e-15** for `legendre_p_all`, **3.79e-13** for real `assoc_legendre_p_all`, **2.48e-15** for normalized complex `assoc_legendre_p_all`, **5.65e-15** for `sph_legendre_p_all`, and **5.65e-15** for `sph_harm_y_all`. The permanent regression test also checks native Fortran bounds, branch convention 3, scalar/table consistency, exact higher derivatives of `P_4`, invalid shapes, NaN/infinity handling, and IEEE-trapping execution. A fresh GNU Fortran 14.2.0 build compiled all **224/224** library sources with `-Werror=implicit-interface`; all **125/125** repository tests passed, and all **64/64** `test_special*.f90` programs passed with invalid/zero/overflow traps enabled. FPM and `fprettify` were not installed in the execution environment.

The Legendre-Q extension adds real/complex `lqn` and `lqmn`. A deterministic **3,780-value/derivative** sweep through order 8 and degree 20 over ten real/complex points on both sides of the unit circle had maximum error scaled by `max(1,abs(reference))` of about **1.49e-12** versus SciPy 1.17.0. A harder **5,330-value/derivative** sweep through order 12 and degree 40 had maximum scaled difference about **2.09e-11**. The permanent tests cover real and complex interior/exterior values, derivatives, zero-based array bounds, invalid dimensions, singular endpoints, `lqn`/`lqmn` consistency, and the finite canonical unit-circle continuation `Q_0(i)=-i*pi/4`; SciPy 1.17.0's legacy backend returns infinities/NaNs at that nonreal unit-circle point.

Focused first-kind tests cover noninteger and negative degrees, positive and negative orders, exact integer-degree pole-cancellation cases, branch cuts 2/3, real-domain rejection when a selected branch is complex, normalized negative orders, derivatives through fourth order, spherical mixed derivatives, endpoint behavior for the value APIs, elemental arrays, invalid inputs, and complex spherical harmonics. Derivatives exactly at square-root branch points remain a documented limitation. `python tools/check_ascii.py` and the 132-character Fortran line-length audit pass. FPM and `fprettify` are not installed in this execution environment, so `fpm test` and automatic formatting could not be run here.


## Struve and Riccati-Bessel validation (2026-09-27)

This batch adds real-order `struve`/`modstruve`, the order-zero integral helpers
`itstruve0`, `it2struve0`, and `itmodstruve0`, plus zero-based
`riccati_jn`/`riccati_yn` sequences and derivatives. Struve H combines the
defining power series, DLMF Bessel-J expansion, and H-Y asymptotics; modified
Struve L uses the positive defining series and L-I asymptotics. Exact negative
half-integer Bessel identities are used before regional selection. The H0
integral helpers integrate the local series analytically and continue with
composite eight-point Gauss-Legendre panels.

A compiled deterministic sweep of **5,040** Struve/modified-Struve evaluations
against SciPy 1.17.0 had maximum error scaled by `max(1,abs(reference))` of
about **7.7e-13** for `struve(-12,x)` in the transition region. A separate
high-precision mpmath grid had maximum scaled errors about **7.0e-13** for H
and **6.3e-14** for L. In a cancellation-sensitive modified-Struve case near
`v=1.5, x=13.36`, the SciFort defining series agrees with high precision more
closely than the SciPy 1.17.0 public result. Permanent tests also use
high-precision references for the three integral helpers; in particular
`itstruve0(40)` is about `3.1484176274471858`, while the observed SciPy 1.17.0
legacy result is badly incorrect. Riccati-Bessel regression values and
derivatives are checked against the defining spherical-Bessel products at
positive, negative, and zero arguments, including exact x=0 limits and invalid
orders.

A separate **1,392-value/derivative** Riccati sweep covered both J and Y
sequences through order 25 at six positive arguments from `0.05` through `30`.
Against SciPy 1.17.0 public outputs, the maximum error scaled by
`max(1,abs(reference))` was about **1.22e-15**.

GNU Fortran 14.2.0 compiled all **216/216** library source files from scratch
with `-std=f2018 -O0 -g -Wall -Wextra -fcheck=all` and
`-ffpe-trap=invalid,zero,overflow`. All **51/51** `test_special*.f90` programs
passed with those traps enabled. The complete repository now contains **112**
program-bearing tests; all **112/112** compiled and passed with `-fcheck=all`
and normal IEEE exception handling, with `TERM=xterm` supplied to the
terminal-oriented test. `python tools/check_ascii.py` and the 132-character
Fortran line-length audit pass. FPM and `fprettify` are not installed in this
execution environment, so `fpm test` and automatic formatting could not be run.

## Order-zero Bessel integral validation (2026-09-27)

This batch adds the paired SciPy-style order-zero Bessel integral APIs
`itj0y0`, `it2j0y0`, `iti0k0`, and `it2i0k0`. The implementation is
independent and uses the DLMF defining/convergent series near the origin,
composite 16-point Gauss-Legendre continuation in the moderate region, and
Bessel-tail asymptotics for large oscillatory arguments. Modified-Bessel K
tails use direct finite-interval quadrature after exponential decay makes the
remaining tail negligible.

A 48-value high-precision reference set covered all eight returned integrals
at `x = 0.1, 1, 10, 20, 40, 100`, spanning the series, quadrature, and
asymptotic branches. References were obtained independently from high-precision
mpmath integration/series and DLMF formulas. The maximum error scaled by
`max(1,abs(reference))` was about **9.43e-15**, in the oscillatory Y0/t
integral at `x=10`; the large-argument x=40 asymptotic values were near
binary64 rounding after extending the recurrence through ten correction
terms. High-precision references are used instead of reproducing several
less-accurate SciPy 1.17.0 legacy outputs in the oscillatory moderate region.

Permanent tests also cover negative-argument parity/domain behavior, zero and
infinite limits, NaN propagation, divergent tail integrals as IEEE infinities,
and elemental-array calls. GNU Fortran compiled all **217/217** library source
files from scratch with `-std=f2018 -O0 -g -fcheck=all`. All **52/52**
`test_special*.f90` programs passed with invalid/zero/overflow traps enabled.
The wider repository regression run was also started and the first 44 program
tests passed before the execution environment's per-command time limit was
reached; no test failure was observed. `python tools/check_ascii.py` passes.
FPM and `fprettify` are not installed in this execution environment, so
`fpm test` and automatic formatting could not be run here.

## Bernoulli, Euler, and Stirling-number validation (2026-09-27)

This batch adds `bernoulli`, `euler`, and `stirling2` to `scifort_special`. The source is an
independent implementation from mathematical identities; SciPy source code was not used.
Bernoulli numbers use the even-index zeta identity, Euler numbers use the secant recurrence
through the representable range, and Stirling numbers use the positive dynamic recurrence.
The latter intentionally avoids SciPy 1.17/1.18's documented Temme approximation for `n>50`.

The validation environment was GNU Fortran 14.2.0 and SciPy 1.17.0. `fpm` was not installed.
Focused builds used

```text
-std=f2018 -O0 -g -fcheck=all -Wall -Wextra
-Wimplicit-interface -Werror=implicit-interface
```

and a second build added

```text
-ffpe-trap=invalid,zero,overflow -fbacktrace
```

`test/test_special_integer_sequences.f90` passed in both configurations. It checks Bernoulli
values through `B_100`, Euler values through `E_100`, the overflow sign at `E_188`, invalid
sequence sizes/status, exact odd-index zeros, scalar Stirling values including `(100,50)`,
invalid Stirling domains, and elemental array evaluation. The committed `(100,50)` Stirling
reference is the exact integer converted to binary64 (`4.3098323700936634e101`); SciPy 1.17.0
returns about `4.309847786526242e101` because it changes to the documented approximation.

The sequence routines return IEEE infinity deliberately when the exact Bernoulli/Euler value
exceeds binary64 range. The Stirling dynamic program guards multiply/add overflow so trapping
builds return positive infinity rather than raising a floating-point exception.

GNU Fortran 14.2.0 compiled all **218/218** library source files from scratch with the strict
flags above. All **53/53** `test_special*.f90` programs passed with
`-ffpe-trap=invalid,zero,overflow` enabled. The complete repository contains **114**
program-bearing tests after this batch; all **114/114** compiled and passed with `-fcheck=all`,
strict interface diagnostics, and normal IEEE exception handling. A separate attempt to apply
invalid-operation trapping to the entire repository stopped in the pre-existing incomplete-beta
NaN-validation path, where a chained Fortran logical expression can evaluate a comparison with
NaN; that path is outside this batch and was not modified. `python tools/check_ascii.py` passes.


## Wright generalized-Bessel validation (2026-09-27)

This batch adds `wright_bessel(a,b,x)` and `log_wright_bessel(a,b,x)` on the
nonnegative real parameter domain used by SciPy. The implementation is
independent and starts from the DLMF 10.46.1 positive defining series. It
locates the continuous maximum of the series term with the digamma derivative,
then sums outward from the neighboring largest integer term after subtracting
its logarithm. Compensated accumulation limits loss in the scaled positive sum.
For exceptionally large modal indices, where a direct tail sweep would be
unnecessarily expensive, the logarithmic API uses a local saddle/Laplace
expansion from polygamma derivatives.

Permanent regression values cover `(a,b,x) = (1,0,1)`, `(1,1,1)`,
`(1.5,1.1,2.5)`, `(0.5,0.25,10)`, and `(0.1,2,100)`, along with the exact
`a=0` reductions. High-precision references were generated independently from
the defining series with mpmath. The tests also check the standard Wright
recurrence, elemental-array evaluation, invalid negative/nonfinite parameters,
positive-infinite `x`, finite logarithms when the unlogged value overflows, and
the large-mode saddle branch. The large-log reference at `a=b=1`, `x=10^6`
was independently checked through `Phi(1,1;x) = I_0(2 sqrt(x))`; the same
identity supplies the `x=10^24` saddle reference.

GNU Fortran 14.2.0 compiled all **219/219** library source files from scratch
with `-std=f2018 -O0 -g -fcheck=all -Wall -Wextra -Wimplicit-interface
-Werror=implicit-interface`. All **54/54** `test_special*.f90` programs passed
with `-ffpe-trap=invalid,zero,overflow -fbacktrace` enabled. The complete
repository contains **115** program-bearing tests after this batch; all
**115/115** compiled and passed with the strict non-trapping flags above.
`python tools/check_ascii.py` and the 132-character audit of the new Fortran
sources pass. FPM and `fprettify` are not installed in this execution
environment, so `fpm test` and automatic formatting could not be run here.

## Lambda, weighted-Bessel-integral, and erf-zero validation (2026-09-27)

This batch adds `lmbda`, `besselpoly`, and `erf_zeros` to `scifort_special`.
`lmbda` evaluates the exact normalized-Bessel/hypergeometric identity, while
`besselpoly` evaluates its defining weighted Bessel-J integral with endpoint
smoothing and composite 16-point Gauss-Legendre quadrature. Near a convergent
singular endpoint the leading Bessel-series contribution is integrated
analytically and only the smoother residual is quadrature evaluated. The
`erf_zeros` routine starts from the NIST DLMF 7.13.1-7.13.2 asymptotic zero
expansion and applies Newton refinement through the existing complex error
function.

Positive-argument `lmbda` values were swept against SciPy 1.17.0 over orders
`0, 0.2, 0.7, 1, 2.3, 5, 10.2, 25.5` and arguments `0, 0.1, 2, 10, 30`; apart
from the documented exact-zero convention difference, the maximum scaled error
was about `7.69e-14`. Negative arguments and the `x=0` limit were independently
checked from `0F1(;v+1;-x^2/4)` and its analytic derivative.

A 95-case `besselpoly` sweep covering ordinary, endpoint-weighted,
noninteger-order, and oscillatory inputs was checked against independent
adaptive quadrature of the defining integral using SciPy's Bessel-J evaluator.
The maximum scaled difference was about `4.29e-13`. Permanent high-precision
regressions additionally include `a=100`, a noninteger-order large-scale case,
and integrable endpoint exponents as close as `-0.99` to the convergence
boundary after the Bessel leading power is included. These references were
checked with high-precision series/integration rather than with the legacy
`besselpoly` implementation itself.

For `erf_zeros`, the first 76 values agree with the finite SciPy 1.17.0 values
to a maximum scaled difference of `5.98e-15`. In that SciPy build, requesting
200 roots yields finite values only through root 76; the SciFort implementation
remains finite through the permanent 1000-root regression. Evaluating those
1000 roots with SciPy's complex `erf` gives a maximum residual of about
`1.08e-12`; the 1000th root is approximately
`56.02435373382717 + 56.06847001647009 i` with residual about `5.19e-13`.

GNU Fortran 14.2.0 compiled all **221/221** library source files from scratch
with `-std=f2018 -O0 -g -fcheck=all -Wall -Wextra -Wimplicit-interface
-Werror=implicit-interface -ffpe-trap=invalid,zero,overflow -fbacktrace`.
All **55/55** `test_special*.f90` programs passed against that trapped build.
A separate strict non-trapping rebuild also compiled all 221 sources; after
compiling the repository's 42 test-support/reference modules, all **116/116**
program-bearing repository tests passed. `python tools/check_ascii.py` and the
132-character Fortran line-length audit pass. FPM and `fprettify` are not
installed in this execution environment, so `fpm test` and automatic
formatting could not be run here.

## Kelvin-function validation (2026-09-27)

This batch adds the real Kelvin functions `ber`, `bei`, `ker`, and `kei`, their
first derivatives `berp`, `beip`, `kerp`, and `keip`, the combined `kelvin`
interface, and all eight zero-table routines exposed by SciPy. The implementation
is original SciFort code based on the NIST DLMF Kelvin/Bessel identities in
sections 10.61 and 10.67; SciPy implementation source was not copied.

The Be family is evaluated from the defining complex Bessel-J relation. A direct
power series is used at moderate argument and a scaled Hankel asymptotic is used
for large argument so exponentially growing components can reach the binary64
limit without intermediate overflow. The Ke family is evaluated from the
complex K0/K1 relation, using convergent I0/I1 harmonic-number series near the
origin and a complex K asymptotic at larger argument. Negative real arguments
use the exact even/odd parity of the Be family; the Ke family returns IEEE NaN
outside its positive-real branch. Exact zero-argument limits are handled
separately.

A 100-digit mpmath sweep checked the four functions and four derivatives at
arguments from `1e-8` through `900`, including both sides of the series/asymptotic
handoffs. The maximum scaled error, with scaling `max(1,abs(reference))`, was
about `2.59e-13`; it occurred for `keip(9.9)` immediately below the K-family
handoff. The first eight roots of each of the eight Kelvin functions agree with
SciPy 1.17.0 to displayed binary64 precision. The permanent tests also request
the 200th roots; for example `ber_zeros(200)` ends at approximately
`886.9106475809849`. The current zero-table API deliberately caps requests at
200 roots, which keeps root evaluation inside the tested binary64 dynamic
range, and reports a nonzero status for larger requests.

GNU Fortran 14.2.0 compiled all **222/222** library source files from scratch
with `-std=f2018 -O0 -g -fcheck=all -Wall -Wextra -Wimplicit-interface
-Werror=implicit-interface -ffpe-trap=invalid,zero,overflow -fbacktrace`.
All **56/56** `test_special*.f90` programs passed against that trapped build.
A separate strict non-trapping rebuild also compiled all 222 sources; after
compiling the repository's 42 shared test-support/reference modules, all
**117/117** program-bearing repository tests compiled and passed. The Kelvin
large-argument stress driver also completed with floating-point traps enabled at
arguments through `1e300`, returning finite values, IEEE infinities, or
underflowed zeros as appropriate without raising an invalid/zero/overflow trap.
FPM and `fprettify` are not installed in this execution environment, so
`fpm test` and automatic formatting could not be run here.

## Angular Mathieu validation (2026-09-27)

This batch adds integer-order angular Mathieu characteristic values
`mathieu_a`/`mathieu_b`, normalized Fourier coefficients
`mathieu_even_coef`/`mathieu_odd_coef`, and the angular functions
`mathieu_cem`/`mathieu_sem`. The implementation is original SciFort code based
on the four Fourier parity eigenproblems, normalizations, and parameter
reflection identities in NIST DLMF section 28.4. SciFort's existing symmetric
Jacobi eigensolver supplies the finite tridiagonal eigenpairs; no SciPy or
Zhang-Jin/specfun implementation source was copied or translated.

A direct comparison with SciPy 1.17.0 covered orders `m=0..10`, parameters
`q=-100,-10,-1,0,1,10,100`, and angular arguments
`x=-123,0,12.3,45,90` degrees. Across 1,617 characteristic values, function
values, and radian derivatives, the largest scaled difference
`abs(actual-reference)/max(1,abs(reference))` was about `4.92e-12`. It occurred
in a derivative with substantial Fourier cancellation; the characteristic
values themselves agree more closely. Permanent tests also check the DLMF
normalizations and Fourier recurrences, exact `q=0` trigonometric reductions,
negative-`q` reflection identities, invalid inputs, and elemental array calls.

One high-parameter regression deliberately does not copy a SciPy 1.17.0 output.
That version returns the identical value `580.0323280120538` for both
`mathieu_a(19,300)` and `mathieu_a(21,300)`. Direct diagonalization of the DLMF
odd-cosine Fourier matrix at dimensions 32, 40, 59, 80, 120, and 180 is stable
to displayed binary64 precision at approximately `512.9573263575256` and
`580.0323280120539`, respectively. SciFort therefore preserves the converged
same-parity eigenvalue ordering instead of reproducing the duplicated legacy
value.

GNU Fortran 14.2.0 compiled all **223/223** library source files from scratch
with `-std=f2018 -O0 -g -fcheck=all -Wall -Wextra -Wimplicit-interface
-Werror=implicit-interface -ffpe-trap=invalid,zero,overflow -fbacktrace`.
All **57/57** `test_special*.f90` programs passed against that trapped build.
A separate strict non-trapping rebuild also compiled all 223 sources; after
compiling the repository's 42 shared test-support/reference modules, all
**118/118** program-bearing repository tests compiled and passed. FPM and
`fprettify` are not installed in this execution environment, so `fpm test` and
automatic formatting could not be run here.

## Radial Mathieu validation (2026-09-27)

This batch adds `mathieu_modcem1`, `mathieu_modcem2`, `mathieu_modsem1`, and
`mathieu_modsem2`. The implementation is original SciFort code based on the DLMF
modified-Mathieu Bessel cross-product expansions and reuses the Fourier coefficients
from the angular implementation. No SciPy/specfun implementation source was copied.

A direct comparison with SciPy 1.17.0 covered orders `m=0..10`, parameters
`q=0.01,0.1,1,10,100`, and nonnegative radial arguments
`x=0,0.1,0.5,1,2,3` radians. Across 2,520 values and derivatives, the largest scaled
difference `abs(actual-reference)/max(1,abs(reference))` was about `3.0e-9`; it
occurred in a very small first-kind value. Adaptive selection of the Bessel
cross-product anchor removes the high-order/high-q cancellation seen when the largest
Fourier coefficient is used mechanically.

SciPy 1.17 numerical outputs agree with the DLMF radial variable in radians, despite
those documentation pages describing `x` as degrees. SciFort follows the numerical
behavior and documents the discrepancy. Negative first-kind arguments use exact
parity. For negative second-kind arguments, direct cross-product summation can be
catastrophically ill-conditioned; SciFort instead continues from the origin using the
same differential-equation symmetry and the first-kind companion solution.

GNU Fortran 14.2.0 compiled all **223/223** library source files from scratch
with `-std=f2018 -O0 -g -fcheck=all -Wall -Wextra -Wimplicit-interface
-Werror=implicit-interface`. All **58/58** `test_special*.f90` programs passed
with `-ffpe-trap=invalid,zero,overflow`. The remaining repository program tests
passed under the same strict runtime/interface checks without global IEEE traps;
combined repository coverage is **119/119** program tests. FPM and `fprettify`
are not installed in this execution environment.

## Spheroidal characteristic-value validation (2026-09-27)

This batch adds the prolate/oblate spheroidal characteristic-value functions
`pro_cv`/`obl_cv` and sequence interfaces `pro_cv_seq`/`obl_cv_seq`. The
implementation is original SciFort code based on the even/odd finite
tridiagonal eigenproblems in NIST DLMF sections 30.3 and 30.16. The DLMF
recurrence matrix is converted by a diagonal similarity transform to a real
symmetric tridiagonal matrix, and the requested ordered eigenvalue is isolated
with a Sturm-count bisection. No SciPy/specfun implementation source was copied.

A direct comparison with SciPy 1.17.0 covered orders `m=0..10`, degrees
`n=m..m+10`, and `c=0.1,1,3,10,25,50,100,200` for both prolate and oblate
families. Across 1,936 characteristic values, all results were finite and the
largest scaled difference `abs(actual-reference)/max(1,abs(reference))` was
about `8.12e-15`. The worst case was `pro_cv(2,2,100)`, where SciFort returned
approximately `103.26840553255833` versus SciPy's `103.2684055325575`.
Permanent tests also cover the exact `c=0` reduction, evenness in `c`, degree
ordering, elemental array evaluation, sequence interfaces, and invalid inputs.

GNU Fortran 14.2.0 compiled all **224/224** library source files from scratch
with `-std=f2018 -O0 -g -fcheck=all -Wall -Wextra -Wimplicit-interface
-Werror=implicit-interface -ffpe-trap=invalid,zero,overflow -fbacktrace`.
All **59/59** `test_special*.f90` programs passed against that trapped build.
A separate strict non-trapping rebuild compiled all 224 sources; after compiling
the repository's 42 shared test-support/reference modules, all **120/120**
program-bearing repository tests passed. FPM and `fprettify` are not installed
in this execution environment, so `fpm test` and automatic formatting could not
be run here.

## Spheroidal angular-function validation (2026-09-27)

This batch adds the first-kind prolate/oblate spheroidal angular functions
`pro_ang1`/`obl_ang1` and precomputed-characteristic interfaces
`pro_ang1_cv`/`obl_ang1_cv`. The implementation is original SciFort code based
on the Ferrers-function expansion and tridiagonal recurrence in NIST DLMF
sections 30.8 and 30.16, with the center normalization described in DLMF 30.1.
No SciPy/specfun implementation source was copied or translated.

A direct comparison with SciPy 1.17.0 covered `m=0..5`, `n=m..m+5`,
`c=0,0.1,0.5,1,2,5,10`, and
`x=-0.99,-0.9,-0.5,-0.1,0,0.1,0.5,0.9,0.99` for both prolate and oblate
families. Across 4,536 function calls, or 9,072 separately compared values and
derivatives, all results were finite. The largest scaled difference
`abs(actual-reference)/max(1,abs(reference))` was about `2.80e-11`, for
`obl_ang1(0,0,10,0)`: SciFort returns exactly `1` by construction while SciPy
1.17.0 returns approximately `1.0000000000279907`. Away from that legacy center
normalization drift, the largest observed discrepancies were substantially smaller;
for example, the largest absolute difference in the grid was about `2.06e-8` in a
derivative of magnitude about `2.71e5`, a scaled difference about `7.59e-14`.
Permanent tests also cover `_cv` agreement, exact `c=0` associated-Legendre
reduction, parity in `x`, center normalization at `c=20`, and invalid inputs.

GNU Fortran 14.2.0 compiled all **224/224** library source files from scratch
with `-std=f2018 -O0 -g -fcheck=all -Wall -Wextra -Wimplicit-interface
-Werror=implicit-interface -ffpe-trap=invalid,zero,overflow -fbacktrace`. All
**60/60** `test_special*.f90` programs passed against that trapped build. A
separate strict non-trapping rebuild compiled all 224 sources; after compiling
the repository's 42 shared test-support/reference modules, all **121/121**
program-bearing repository tests passed. FPM and `fprettify` are not installed
in this execution environment, so `fpm test` and automatic formatting could not
be run here.

## Spheroidal first-kind radial validation (2026-09-27)

This batch adds the first-kind prolate/oblate spheroidal radial functions
`pro_rad1`/`obl_rad1` and precomputed-characteristic interfaces
`pro_rad1_cv`/`obl_rad1_cv`. The implementation is original SciFort code based
on the spherical-Bessel radial expansion in NIST DLMF 30.11.3 and the
`a_n,k^-m` coefficient relation in DLMF 30.8.6. No SciPy/specfun
implementation source was copied, translated, or adapted.

The implementation reuses SciFort's tridiagonal spheroidal eigenvector
coefficients, converts them to the negative-order coefficients required by the
radial expansion with logarithmic scaling, and uses compensated summation for
both the Bessel series and the `A_n^-m` normalization. The oblate `x=0` endpoint
is evaluated analytically from the small-argument spherical-Bessel limit. The
prolate API follows the current SciPy domain `x > 1`; the oblate API accepts
`x >= 0`. This batch covers the first-kind functions only; second-kind radial
functions remain future work.

A direct comparison with SciPy 1.17.0 covered `m=0..5`, `n=m..m+5`,
`c=0.1,0.5,1,2,5,10`, prolate
`x=1.001,1.01,1.1,1.5,2,5,10`, and oblate
`x=0,0.01,0.1,0.5,1,2,5,10`. Across 3,240 function calls, or 6,480
separately compared values and derivatives, the largest scaled difference
`abs(actual-reference)/max(1,abs(reference))` was about `2.96e-12`, for
`obl_rad1(0,0,10,0)`: SciFort returned approximately
`0.1057735930149705` versus SciPy's `0.10577359301793428`. Permanent tests
also cover representative prolate/oblate values and derivatives, near-boundary
prolate evaluation, exact oblate-origin parity limits, `_cv` agreement, the
SciPy-compatible `c=0` limiting convention, and invalid inputs.

GNU Fortran 14.2.0 compiled all **224/224** library source files from scratch
with `-std=f2018 -O0 -g -fcheck=all -Wall -Wextra -Wimplicit-interface
-Werror=implicit-interface -ffpe-trap=invalid,zero,overflow -fbacktrace`. All
**61/61** `test_special*.f90` programs passed against that trapped build. A
separate non-trapping rebuild compiled all 224 sources; after compiling the
repository's 42 shared test-support/reference modules, all **122/122**
program-bearing repository tests passed. `python tools/check_ascii.py` and the
132-character source-line audit also pass. FPM and `fprettify` are not installed
in this execution environment, so `fpm test` and automatic formatting could not
be run here.

## Spheroidal second-kind radial validation (2026-09-27)

This batch adds the second-kind prolate/oblate spheroidal radial functions
`pro_rad2`/`obl_rad2` and precomputed-characteristic interfaces
`pro_rad2_cv`/`obl_rad2_cv`. The implementation is original SciFort code based
on the radial spherical-Bessel expansion in NIST DLMF 30.11.3, using the
second-kind spherical function specified by DLMF 30.11.1-30.11.2. No
SciPy/specfun implementation source was copied, translated, or adapted.

Direct spherical-Bessel-Y summation is well-conditioned only away from the
prolate singular boundary and the small-oblate-argument region. SciFort
therefore evaluates the DLMF series at an automatically selected stable anchor,
using the same independently reconstructed `a_n,k^-m` coefficients and
`A_n^-m` normalization as the first-kind implementation. When the requested
argument lies below that anchor, an independently written embedded
Dormand-Prince RK5(4) integrator continues the value and derivative with the
prolate or oblate spheroidal differential equation. The tableau is the published
Dormand-Prince 1980 formula (doi:10.1016/0771-050X(80)90013-3). This avoids importing the legacy
small-argument second-kind algorithms. Positive finite `c` is required; the
mathematically singular `c=0` second-kind normalization returns NaN rather than
a legacy finite sentinel.

A direct comparison with SciPy 1.17.0 in the stable spherical-Y expansion
region covered `m=0..5`, `n=m..m+5`, `c=0.2,0.5,1,2,5,10`, both prolate and
oblate functions, and three radial arguments beginning at
`max(2,(n+m+12)/c)`. Across 1,296 function calls, or 2,592 separately compared
values and derivatives, there were no nonfinite results and the largest scaled
difference `abs(actual-reference)/max(1,abs(reference))` was about `1.18e-12`,
for the derivative of `pro_rad2(0,0,10,2)`.

The difficult region was checked independently through the radial Wronskian.
For 4,104 cases spanning the same `m`, `n`, and `c` grid, prolate
`x=1.001,1.01,1.05,1.2,1.5,2,3,5,10`, and oblate
`x=0,0.01,0.1,0.5,1,1.5,2,3,5,10`, the largest error in the scaled identity
`c*(x*x-1)*(R1*R2'-R1'*R2)=1` for prolate functions or
`c*(1+x*x)*(R1*R2'-R1'*R2)=1` for oblate functions was about `8.48e-12`.
Permanent tests include stable SciPy regression values, near-boundary and
small-oblate Wronskian checks, `_cv` agreement, domain checks, and the `c=0`
NaN convention.

GNU Fortran 14.2.0 compiled all **224/224** library source files from scratch
with `-std=f2018 -O0 -g -fcheck=all -Wall -Wextra -Wimplicit-interface
-Werror=implicit-interface -ffpe-trap=invalid,zero,overflow -fbacktrace`. All
**62/62** `test_special*.f90` programs passed against that trapped build. A
separate non-trapping build compiled all 224 sources; after compiling the
repository's 42 shared test-support/reference modules, all **123/123**
program-bearing repository tests passed. The ASCII and 132-character Fortran
source-line audits also pass. FPM and `fprettify` are not installed in this
execution environment, so `fpm test` and automatic formatting could not be run
here.

## Parabolic-cylinder W validation (2026-09-28)

This batch adds `pbwa(a,x,w,wprime)`, the real DLMF parabolic-cylinder
`W(a,x)` function and its first derivative. The implementation is original
SciFort code based on NIST DLMF Chapter 12. Exact DLMF gamma-ratio values at
`x=0` are evaluated through complex `loggamma`, and a Taylor recurrence derived
directly from `y'' + (x*x/4-a)*y = 0` continues the solution to the requested
argument. No SciPy/specfun implementation source was copied, translated, or
adapted. SciPy's legacy real `pbwa` domain is retained: finite
`abs(a) <= 5` and `abs(x) <= 5`; other inputs return NaN.

Permanent regression values were generated independently at 80-digit precision
from the DLMF origin normalization and confluent-hypergeometric representation.
The tests cover the origin, positive and negative parameters and arguments, all
four domain edges, the exact Wronskian of `W(a,x)` and `W(a,-x)`, elemental
array evaluation, out-of-domain inputs, and NaN propagation. A debug run with
`-ffpe-trap=invalid,zero,overflow` also passes; this check caught and removed an
initial dependence on compiler short-circuit evaluation in the NaN/domain guard.

A 41-by-41 direct comparison grid against SciPy 1.17.0 covered 1,681 calls over
`a,x = -5,-4.75,...,5`. All SciFort results were finite. The largest scaled
value difference `abs(actual-reference)/max(1,abs(reference))` was about
`5.47e-11` at `a=0,x=-5`; the largest derivative difference was about
`2.33e-10` at `a=0,x=5`. The independent 80-digit references used by the
permanent endpoint tests agree with the SciFort values, so the edge differences
are not used to reproduce the lower-accuracy legacy SciPy backend.

GNU Fortran 14.2.0 compiled all **224/224** library source files from scratch
with Fortran 2018 mode, runtime checks, warnings, implicit-interface errors,
floating-point traps, and backtraces enabled. With the repository's **42**
shared test-support/reference module files, all **124/124** program-bearing
repository tests compiled and passed in the non-trapping runtime sweep. All
**63/63** `test_special*.f90` programs also passed when linked with
`-ffpe-trap=invalid,zero,overflow`. `python tools/check_ascii.py` and the
132-character Fortran source-line audit pass. FPM and `fprettify` are not
installed in this execution environment, so `fpm test` and automatic formatting
could not be run here.

## Ellipsoidal-harmonic validation (2026-09-28)

This batch adds the real ellipsoidal-harmonic family `ellip_harm`,
`ellip_harm_2`, and `ellip_normal`. The first-kind Lame-polynomial class
matrices, diagonal symmetrization, coefficient scaling, and class-prefactor
evaluation are adapted from SciPy 1.17.0 `scipy/special/_ellip_harm.pxd` under
its retained BSD-2-Clause notice. SciFort replaces the upstream LAPACK
tridiagonal eigensolver with its existing pure symmetric Jacobi eigensolver.
The second-kind and normalization routines are independently written from the
published defining integral identities. Their endpoint weights are removed by
trigonometric substitutions before adaptive 16-point Gauss-Legendre
quadrature, so no SciPy QUADPACK/Cython integration code is copied or
translated.

A deterministic direct comparison with SciPy 1.17.0 covered both focal-parameter
pairs `(h2,k2)=(1,4)` and `(5,8)`, every degree `n=0..8`, every mode
`p=1..2*n+1`, five first-kind coordinates `s=-3.2,0.7,2,3.2,7`, two exterior
second-kind coordinates for every mode, and every corresponding normalization
constant. All **1,296/1,296** finite comparisons passed a scaled error threshold
of `1e-10`. The largest scaled first-kind difference
`abs(actual-reference)/max(1,abs(reference))` was about `8.74e-13`, at
`ellip_harm(5,8,8,9,0.7)`. The largest second-kind difference was about
`5.55e-17`, and the largest normalization difference was about `5.21e-14`.
Permanent tests also cover all four Lame classes, the optional square-root sign
controls, negative exterior coordinates, elemental array evaluation, invalid
parameters, NaNs, infinities, and the second-kind exterior-domain restriction.

GNU Fortran 14.2.0 compiled all **225/225** library source files from scratch
with Fortran 2018 mode, runtime checks, warnings, implicit-interface errors,
and backtraces enabled. The repository's **42/42** shared test-support/reference
modules compiled, and all **126/126** program-bearing repository tests passed.
All **65/65** `test_special*.f90` programs also passed in executables with
`-ffpe-trap=invalid,zero,overflow`. `python tools/check_ascii.py` and the
132-character Fortran source-line audit pass. FPM and `fprettify` are not
installed in this execution environment, so `fpm test` and automatic formatting
could not be run here.

## Real-axis Hankel validation (2026-09-28)

This batch adds the real-argument subset of SciPy's complex Hankel family:
`hankel1`, `hankel2`, `hankel1e`, `hankel2e`, `h1vp`, and `h2vp`. Positive
arguments combine SciFort's existing independently validated real-order `jv`
and `yv` kernels. Negative real arguments use upper-edge principal-branch
continuation derived from NIST DLMF 10.11.3-10.11.4. Derivatives reuse the
existing DLMF 10.6.7 shifted-order `jvp`/`yvp` kernels. This checkpoint
established the real-axis baseline; the complex-argument continuation is
validated in the following section.

A deterministic comparison with SciPy 1.17.0 covered ten orders from `-3.7`
through `20.25`, nine positive and negative real arguments from `0.1` through
`100` in magnitude, all four base/scaled Hankel APIs, and derivative orders
`0`, `1`, `2`, `5`, and `20` for both derivative APIs. All **1,260/1,260**
finite comparisons passed a scaled-error threshold of `2e-10`. The largest
scaled error `abs(actual-reference)/max(1,abs(reference))` was about
`2.62e-14`, for `hankel2(0.2,-20)`. Permanent tests additionally cover the
positive-axis J/Y identities, negative-axis continuation, negative order,
scaled values, elemental arrays, zero/NaN/infinity behavior, the SciPy
`hankel2(0,0)` edge case, and a twentieth derivative.

GNU Fortran 14.2.0 compiled all **226/226** library source files from scratch
with Fortran 2018 mode, runtime checks, warnings, implicit-interface errors,
and backtraces enabled. The repository's **42/42** shared test-support modules
compiled, and all **127/127** program-bearing repository tests passed. All
**66/66** `test_special*.f90` programs also passed when the main programs were
compiled and linked with `-ffpe-trap=invalid,zero,overflow`. The repository
ASCII audit and the 132-character Fortran source-line audit pass. FPM and
`fprettify` are not installed in this execution environment, so `fpm test` and
automatic formatting could not be run here.

## Complex-argument Hankel validation (2026-09-28)

This batch extends `hankel1`, `hankel2`, `hankel1e`, `hankel2e`, `h1vp`, and
`h2vp` to complex arguments while preserving the established real-axis
overloads and principal-branch behavior. The implementation is original
SciFort code based on NIST DLMF connection and continuation formulas. In the
right half-plane, the exponentially favorable scaled Hankel solution is
obtained from the DLMF 10.27.8 relation to a scaled modified Bessel K solution.
That K solution is analytically continued from the positive real axis along a
constant-radius arc by integrating the scaled modified-Bessel differential
equation with an embedded Dormand-Prince 5(4) method. The complementary scaled
Hankel solution is continued directly from the positive real axis with the
scaled cylindrical-Bessel equation. Left-half-plane values use the DLMF 10.11
Hankel continuation formulas in scaled form, with exponentially decaying
factors selected to reduce avoidable overflow. Complex derivatives use the DLMF
10.6.7 shifted-order finite sum, accumulated in scaled form before restoring
the common exponential factor.

A deterministic comparison with SciPy 1.17.0 covered six real orders from
`-3.7` through `20.25` and eight complex arguments spanning both half-planes
and both sides of the negative-real branch cut, including `20+10i` and
`20-10i`. The sweep exercised all four base/scaled APIs plus first derivatives
of `hankel1` and second derivatives of `hankel2`, for **288/288** comparisons.
Every comparison passed a scaled-error threshold of `2e-9`; the largest scaled
error was about **3.58e-11**, for `hankel2(20.25,20+10i)`. Permanent regression
tests additionally cover upper- and lower-half-plane values, both sides of the
negative-axis branch cut, scaled values, negative order, first through third
derivatives, Schwarz reflection, the `n=0` derivative identity, exact
real-axis dispatch, elemental complex arrays, and zero/NaN/infinity/error
behavior.

GNU Fortran 14.2.0 compiled all **226/226** library source files from scratch
with Fortran 2018 mode, runtime checks, warnings, implicit-interface warnings,
and backtraces enabled. The repository's **42/42** shared test-support/reference modules
compiled, and all **128/128** program-bearing repository tests passed. All
**67/67** `test_special*.f90` programs also passed when their main programs
were compiled with `-ffpe-trap=invalid,zero,overflow`.
`python tools/check_ascii.py` passes, and the 132-character Fortran source-line audit
reports no overlength source/test lines. FPM and `fprettify` are not installed
in this execution environment, so `fpm test` and automatic formatting could
not be run here.


## Complex cylindrical and modified Bessel validation (2026-09-28)

This batch extends the public `scifort_special` overloads for `jv`, `jve`, `yv`,
`yve`, `jn`, `iv`, `ive`, `kv`, `kve`, `jvp`, `yvp`, `ivp`, and `kvp` to
complex arguments while preserving every established real-argument behavior.
J/Y values are reconstructed from the scaled Hankel pair except when a
nonnegative order (or equivalent negative integer order) is large compared with
`abs(z)`. In that cancellation-prone regime J is evaluated from the DLMF 10.2.2
power series directly. I uses the DLMF 10.27 quarter-turn connection to J, and K
uses the corresponding connection to the scaled Hankel functions. The
half-plane-dependent rotations preserve the principal branch; an exactly negative
real complex argument uses its upper edge. Complex derivatives reuse the DLMF
10.6.7/10.29.5 shifted-order finite sums with scaled accumulation.

A deterministic comparison with SciPy 1.17.0 covered seven real orders from
`-3.7` through `20.25` and nine complex arguments spanning both half-planes,
both sides of the negative-real branch region, the exact negative-real upper
edge, a small complex argument, and `20+/-10i`. The sweep exercised all eight
base/scaled J/Y/I/K APIs and first and second derivatives of all four families,
for **1008/1008** comparisons. Every comparison passed the scaled-error measure
`abs(actual-reference)/max(1,abs(reference)) < 2e-9`; the largest scaled error
was about **9.15e-10**, for `iv(20.25,20+10i)`. The permanent regression test
also pins the large-order/small-complex-argument case where naive Hankel
reconstruction catastrophically cancels, plus negative order, the principal
negative-axis branch, elemental arrays, `n=0`, zero/NaN/infinity handling, and
negative derivative orders.

GNU Fortran 14.2.0 compiled all **227/227** library source files from scratch
with Fortran 2018 mode, runtime checks, warnings, implicit-interface warnings,
and backtraces enabled. The repository's **42/42** shared test-support/reference
modules compiled, and all **129/129** program-bearing repository tests passed.
All **68/68** `test_special*.f90` programs passed when their main programs were
compiled with `-ffpe-trap=invalid,zero,overflow`. FPM and `fprettify` are not
installed in this execution environment, so `fpm test` and automatic formatting
could not be run here.

## Complex spherical Bessel validation (2026-09-28)

This batch extends `spherical_jn`, `spherical_yn`, `spherical_in`, and
`spherical_kn` to complex binary64 arguments, including the SciPy-style optional
first derivative. For genuinely nonreal arguments, values generally use the DLMF
10.47 half-integer definitions in terms of SciFort's complex cylindrical J/Y/I/K
kernels. High-order spherical J near its turning region uses normalized backward
Miller recurrence from the DLMF 10.51 identities. The shared complex-I kernel now
uses its direct defining power series when the order is large relative to
`abs(z)`, avoiding cancellation in the quarter-turn path. The scaled-Hankel arc
continuation tolerance was also tightened after stress testing exposed a
high-order spherical-Y accuracy limit. Nonzero complex arguments lying exactly
on the real axis continue to dispatch to the established real spherical kernels.

The broad deterministic comparison with SciPy 1.17.0 covered five integer orders
from 0 through 20 and six complex arguments spanning all four quadrants, a
small-magnitude complex point, a near-negative-axis point, and `10+5i`. All four
spherical families and both value/first-derivative modes passed, for **240/240**
comparisons. Every comparison satisfied
`abs(actual-reference)/max(1,abs(reference)) < 2e-9`; after the high-order
hardening, the largest scaled error in this sweep was about **4.02e-14**.

A second stress sweep covered orders 0, 1, 5, 20, 50, and 100 at six complex
arguments including `100+20i` and `-20+10i`. All **288/288** spherical
value/derivative comparisons passed the same `2e-9` criterion; the largest
scaled error was about **1.05e-9**, for `spherical_yn(100,100+20i)`. A separate
shared-kernel sweep exercised `iv`, `ive`, and first/second `ivp` derivatives at
orders through 100.5 over ten complex arguments. All **280/280** comparisons
passed, with worst scaled error about **3.21e-14**. Permanent regressions pin the
order-100 turning-region cases as well as large-order/small-complex values,
elemental arrays, exact real-axis dispatch, zero behavior, negative order, NaN,
and infinity handling.

GNU Fortran 14.2.0 compiled all **228/228** library source files from scratch
with Fortran 2018 mode, runtime checks, warnings, and implicit-interface errors
enabled. The repository's **42/42** shared test-support/reference module files
compiled, and all **130/130** program-bearing repository tests passed after the
high-order hardening. All **69/69** `test_special*.f90` programs also passed
when compiled and linked with `-ffpe-trap=invalid,zero,overflow`.
`python tools/check_ascii.py` passes, and the 132-character Fortran source/test
line audit reports no overlength lines. FPM and `fprettify` are not installed in
this execution environment, so `fpm test` and automatic formatting could not be
run here.



## Complex elementary and normal-CDF validation (2026-09-28)

This batch extends `ndtr`, `log_ndtr`, `expm1`, and `log1p` to complex binary64
arguments and adds all real/complex combinations for `xlogy` and `xlog1py`.
Complex `ndtr` uses SciFort's Faddeeva-based complex complementary error
function. `log_ndtr` evaluates the principal logarithm with three numerical
regimes: a scaled-`erfcx` identity in the left half-plane, direct evaluation in
the central region, and a reflected `log1p` form when the complementary normal
probability is small. Complex `expm1` and `log1p` use Cartesian rearrangements
that retain small increments, while `xlogy`/`xlog1py` preserve the SciPy-style
exact-zero multiplier convention before taking a logarithm. Established real
overloads remain unchanged.

A deterministic comparison with SciPy 1.17.0 used a 13-by-13 complex grid for
each of `ndtr`, `log_ndtr`, `expm1`, and `log1p`, spanning real and imaginary
parts from about -30 through +30, for **676** comparisons. A separate 8-by-8
complex-pair grid exercised both `xlogy` and `xlog1py`, including mixed-sign and
branch-sensitive values, for another **128** comparisons. All **804/804**
comparisons satisfied
`abs(actual-reference)/max(1,abs(reference)) < 2e-11`. The largest scaled errors
were about **1.79e-13** for `ndtr`, **3.72e-14** for `log_ndtr`, **0.00**
at displayed precision for `expm1`, **8.33e-17** for `log1p`, zero to displayed
precision for `xlogy`,
and **1.84e-16** for `xlog1py`. The permanent regression test also pins tiny
complex increments, both sides of the `log1p` negative-real branch cut, mixed
real/complex overload resolution, exact-zero multipliers, real-overload parity,
elemental arrays, and NaN propagation.

GNU Fortran 14.2.0 compiled all **229/229** library source files from scratch in
Fortran 2018 mode with implicit-interface errors and runtime checks enabled. The
repository's **42/42** shared test-support/reference modules compiled, and all
**131/131** program-bearing repository tests passed. All **70/70**
`test_special*.f90` programs passed when compiled and linked with
`-ffpe-trap=invalid,zero,overflow -fcheck=all -fbacktrace`.
`python tools/check_ascii.py` passes, and the 132-character Fortran source/test
line audit reports no overlength lines. FPM and `fprettify` are not installed in
this execution environment, so `fpm test` and automatic formatting could not be
run here.

## Orthogonal-polynomial constructor validation (2026-09-28)

This batch adds coefficient constructors for the 15 classical SciPy orthogonal-
polynomial families already covered by SciFort evaluators and Gaussian roots.
Constructors return descending-power binary64 coefficient arrays and support the
SciPy-style optional `monic` normalization. Coefficients are generated directly
from independent DLMF Chapter 18 recurrences; shifted Legendre/Chebyshev/Jacobi
families use recurrences in the shifted variable rather than an after-the-fact
monomial translation, avoiding avoidable cancellation. Monic Gegenbauer uses an
explicit coefficient ratio, including the analytic `alpha=0` limit.

A deterministic SciPy 1.17.0 comparison covered degrees 0 through 8, both monic
states, all 11 parameter-free families, five Jacobi parameter pairs, four ordinary
Gegenbauer parameters, five generalized-Laguerre parameters, and five shifted-
Jacobi parameter pairs. Across **2,700** coefficients, the largest error scaled by
`max(1,max(abs(reference coefficients)))` was about **7.99e-15**. The permanent
regression test additionally pins the exact `alpha=0` monic Gegenbauer limit,
constructor/evaluator consistency, degree-zero `chebyc`, monic normalization, and
invalid-domain empty-array behavior. High-degree monomial coefficient arrays are
inherently ill-conditioned for direct evaluation and SciPy's returned arrays can
contain small parity-forbidden roundoff terms, so validation does not treat those
artifacts as required compatibility behavior.

GNU Fortran 14.2.0 compiled all **230/230** library source files from scratch in
Fortran 2018 mode with warnings, implicit-interface errors, runtime checks, and
backtraces enabled. The repository's **42/42** shared test-support/reference
modules compiled, and all **132/132** program-bearing repository tests passed. All
**71/71** `test_special*.f90` programs passed when compiled and linked with
`-ffpe-trap=invalid,zero,overflow -fcheck=all -fbacktrace`. `python
tools/check_ascii.py` passes and the 132-character source/test audit reports no
overlength Fortran lines. FPM and `fprettify` are not installed in this execution
environment, so `fpm test` and automatic formatting could not be run here.

## Initial scipy.integrate quadrature validation (2026-09-28)

This batch introduces the first `scipy.integrate`-style API surface: sampled-data
`trapezoid`, `cumulative_trapezoid`, `simpson`, and `cumulative_simpson`, Romberg
integration with `romb`, fixed-order Gauss-Legendre integration with `fixed_quad`,
Newton-Cotes weights with `newton_cotes`, and scalar adaptive `quad`. The sampled-
data implementation follows SciPy 1.17.0 edge semantics where practical in a
one-dimensional Fortran API, including decreasing and duplicate abscissae for
`simpson`, strictly increasing explicit abscissae for `cumulative_simpson`, and a
zero-only explicit initial value for cumulative trapezoids.

A deterministic reference program compared the sampled-data and Newton-Cotes
routines with SciPy 1.17.0 across irregular and uniform grids, negative and zero
spacing, Romberg grids, built-in Newton-Cotes orders 1 through 14, and irregular
Newton-Cotes node sets. All **643/643** individual values and weights passed with
runtime checks and `-ffpe-trap=invalid,zero,overflow`. The permanent integration
regression additionally covers high-order Newton-Cotes moment/symmetry checks,
fixed Gauss-Legendre quadrature, finite/reversed adaptive bounds, semi-infinite
and whole-line dispatch, and invalid-argument status handling.

Adaptive `quad` is an original SciFort facade over the external BSD-3-Clause
`jacobwilliams/quadpack` package pinned at tag 2.1.0. Finite bounds dispatch to
`dqags`; semi-infinite and whole-line bounds dispatch to `dqagi`. The tagged
Quadpack2 source and public interfaces were inspected directly, but the external
dependency could not be downloaded or built in this execution environment:
`fpm` is not installed and shell downloads of the dependency are unavailable.
For local compilation and control-flow testing only, an API-compatible temporary
`quadpack` stub was used outside the repository. Therefore the local adaptive-
quadrature tests validate SciFort's facade, bound classification, sign handling,
defaults, and status propagation; they are not an independent numerical
validation of the real Quadpack2 implementation. The temporary stub is not
included in the project or release archive.

GNU Fortran 14.2.0 compiled all **233/233** SciFort library sources from scratch
in Fortran 2018 mode with warnings, implicit-interface errors, runtime checks,
and backtraces enabled. The repository's **42/42** shared test-support/reference
modules compiled and all **133/133** program-bearing repository tests passed.
All **71/71** `test_special*.f90` programs passed with
`-ffpe-trap=invalid,zero,overflow -fcheck=all -fbacktrace`, and the new
integration regression also passed with those floating-point traps. FPM and
`fprettify` are not installed in this execution environment, so `fpm test` and
automatic formatting could not be run here.

## Specialized QUADPACK modes and nested integration validation (2026-09-28)

This batch extends `quad` with finite-interval break points and SciPy-style
weighted integration. Break-point integration dispatches to Quadpack2 `dqagp`.
Finite cosine/sine weights dispatch to `dqawo`; one-sided infinite cosine/sine
weights dispatch to `dqawf`; algebraic/logarithmic endpoint weights dispatch to
`dqaws`; and Cauchy principal-value weights dispatch to `dqawc`. The same batch
adds original SciFort `dblquad`, `tplquad`, and `nquad`-style nesting helpers on
top of the public `quad` facade.

Fifteen representative numerical reference results were regenerated directly
with installed SciPy 1.17.0 before the final validation. They cover a supplied
break point, finite and one-sided-infinite cosine/sine weights (including a
left-infinite sine integral), all four algebraic/logarithmic endpoint weights, a
Cauchy principal value, variable-bound `dblquad` and `tplquad`, and the equivalent
three-dimensional `nquad` case. The permanent tests contain those numerical
references plus reversed-limit, duplicate/end-point break-point, invalid-option,
constant-bound, unit-box, and analytic polynomial checks.

As in the initial integration batch, the external Quadpack2 dependency itself
could not be downloaded or built in this execution environment because `fpm` is
not installed and shell network access is unavailable. The exact tagged 2.1.0
public interfaces for `dqagp`, `dqawo`, `dqawf`, `dqaws`, and `dqawc` were
inspected upstream. For local compilation and dispatch testing, an API-compatible
temporary `quadpack` harness was used outside the repository. It is not included
in the project or release archive. Consequently, the local weighted tests verify
SciFort's argument validation, workspace sizing, mode selection, reflection/sign
handling, nesting, and status propagation; they do not substitute for an
independent numerical validation of the real Quadpack2 implementation.

GNU Fortran 14.2.0 compiled all **234/234** SciFort library sources under
Fortran 2018 mode with implicit-interface errors enabled. All **42/42** shared
test-support/reference modules compiled, and all **135/135** program-bearing
repository tests passed. The **71/71** `test_special*.f90` programs passed with
`-ffpe-trap=invalid,zero,overflow -fcheck=all -fbacktrace`. All **3/3**
`test_integrate*.f90` programs also passed with those runtime checks; the
strict run exposed and led to the correction of a missing `recursive` attribute
on the internal callback used by recursive `nquad` levels. FPM and `fprettify`
are not installed in this execution environment, so `fpm test` and automatic
formatting could not be run here.

## Advanced integration validation (2026-09-28)

This batch adds the remaining one-dimensional/Monte-Carlo integration routines
before multidimensional cubature: `quad_vec`, `qmc_quad`, and `tanhsinh`.
`quad_vec` is intentionally a correctness-first componentwise facade over the
existing scalar QUADPACK wrapper, so its numerical values inherit the scalar
`quad` validation while its error estimate is the requested 2-norm or max-norm
of the component error estimates. `qmc_quad` uses independently scrambled
SciFort Halton engines initialized from an explicit `rng_state`; a two-dimensional
product integral over `[0,1] x [0,2]` with 8 estimates of 2048 points returned
`1.0001989963439561` for the exact value 1, with standard error
`1.7080241611021205e-4`. Constant-integrand regressions pin reversed-bound signs,
zero replicate variance, and the NaN standard error for one estimate.

The `tanhsinh` implementation was compared directly with SciPy 1.17.0 over 15
finite, endpoint-singular, semi-infinite, and whole-line cases. The set includes
powers `x**p` for `p=-0.5,-0.25,0,1,5,20`, `log(x)` on `[0,1]`, exponential
tails at three decay rates, Gaussian whole-line integrals at three rates, the
Cauchy integral, and an endpoint-singular beta integral. All 15 comparisons
agreed with scaled error below `1.13e-13`; the largest discrepancy was the
`exp(-5*x)` semi-infinite case, where SciPy 1.17.0 reported
`0.20000000000011237` and SciFort returned `0.19999999999999998`. A permanent
slow-tail regression integrates `x**(-1.01)` from 1 to infinity and confirms
that the transformed-tail guard returns `tanhsinh_status_maxlevel` rather than
claiming convergence when the fixed binary64 truncation boundary is still
material.

GNU Fortran 14.2.0 compiled all **237/237** library sources under Fortran 2018
mode with runtime checks, warnings, and implicit-interface errors enabled. The
**42/42** shared test-support/reference modules compiled, and all **136/136**
program-bearing repository tests passed. The four `test_integrate*.f90`
programs and all **71/71** `test_special*.f90` programs also passed when linked
with `-fcheck=all -ffpe-trap=invalid,zero,overflow`. The local build used an
API-compatible QUADPACK stand-in outside the repository because FPM and external
shell network access are unavailable in this execution environment; the stand-in
is not included in the project archive. Therefore these runs validate SciFort's
QUADPACK dispatch/interface layer but are not an independent numerical test of
the external Quadpack2 implementation. `fpm test` and `fprettify` could not be
run because those executables are not installed here.


## Multidimensional cubature validation (2026-09-28)

This batch adds adaptive `cubature` for scalar and vector-valued callbacks over
finite, one-sided-infinite, whole-line, and independently reversed
hyperrectangles. It provides tensor-product `gk15`/`gk21` rules and the
Genz-Malik degree-7/degree-5 embedded rule, user feature-point partitioning,
global componentwise error accumulation, and largest-error-region refinement.
The Genz-Malik formulas, Gauss-Kronrod tables, and high-level split strategy are
adapted from SciPy 1.17.0 under its BSD-3-Clause license; the Fortran storage,
callbacks, transformations, and status handling are SciFort code.

A deterministic comparison sweep evaluated 27 finite integrals against installed
SciPy 1.17.0: dimensions 2 through 4, three smooth integrand families, and each
of `genz-malik`, `gk15`, and `gk21`. All **27/27** estimates agreed; the worst
scaled estimate difference was `2.42e-14`, and all 27 subdivision counts matched
SciPy. Permanent regressions additionally cover scalar/vector output, reversed
bounds, a two-dimensional cusp split by `points=`, the whole-line Gaussian, a
left-infinite exponential integral, exact Genz-Malik/Gauss-Kronrod polynomial
cases, and max-subdivision termination. The whole-line and left-infinite tests
use analytic references so correctness of SciFort's transformed-domain mapping
does not depend on another implementation's infinite-limit transformation.

GNU Fortran 14.2.0 compiled all **238/238** library sources under Fortran 2018
mode with runtime checks, warnings, and implicit-interface errors enabled. The
**42/42** shared test-support/reference modules compiled, and all **137/137**
program-bearing repository tests passed. All **5/5** `test_integrate*.f90`
programs and all **71/71** `test_special*.f90` programs passed with
`-fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace`. As in prior
integration batches, an API-compatible QUADPACK stand-in outside the repository
was used only because FPM and external shell network access are unavailable in
this execution environment; it is not included in the archive. The new
`cubature` implementation itself does not depend on that stand-in. `fpm test`
and `fprettify` could not be run because those executables are not installed.


## Initial-value ODE validation (2026-09-28)

This batch starts the `solve_ivp` milestone with real-system RK23 and RK45.
The Bogacki-Shampine 3(2) and Dormand-Prince 5(4) tableaux are adapted from
Jacob Williams' BSD-3-Clause `rklib`; the adaptive driver, result type,
`t_eval`, backward integration, cubic-Hermite dense interpolation, and scalar
event handling are SciFort code. The permanent regression integrates
`y'=y`, a two-state harmonic oscillator, a backward exponential trajectory,
dense-output queries, a terminal event at `y=2`, and invalid method input.

A separate comparison with installed SciPy 1.17.0 used 11 RK45 exponential
`t_eval` points on `[0,2]`, 9 RK23 harmonic-oscillator points on `[0,pi/2]`,
and the terminal event for `y'=y`. The largest scaled RK45 sampled-output
difference was `3.46e-8`; this is dominated by SciFort's generic cubic-Hermite
interpolant whereas SciPy RK45 uses a quartic dense interpolant. The RK23
harmonic outputs differed by at most `1.0e-15`, and the located terminal-event
time differed from SciPy by about `3.3e-8`. Direct analytic endpoint tests are
tighter: RK45 gives `exp(2)` within `2e-9` scaled and RK23 gives the quarter-
period oscillator state within `2e-7` at the test tolerances.

GNU Fortran 14.2.0 compiled all **239/239** library sources under Fortran 2018
mode with runtime checks, warnings, and implicit-interface errors enabled. All
**42/42** shared test-support/reference modules compiled. The complete
program-bearing repository suite passed **138/138** tests in batches. All
**6/6** `test_integrate*.f90` programs and all **71/71**
`test_special*.f90` programs passed with
`-fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace`.

As in prior integration validation, the external Quadpack2 dependency could not
be downloaded because shell network access is unavailable, so an API-compatible
QUADPACK stand-in outside the repository was used only to satisfy/link the
existing quadrature wrapper during whole-project testing. It is not included in
the archive and is unrelated to the new ODE implementation. `fpm` and
`fprettify` are not installed in this environment, so `fpm test` and automatic
formatting could not be run.

## DOP853 validation (2026-09-28)

This batch adds `method='DOP853'` to the existing `solve_ivp` framework. The
DOP853 stage tableau, eighth-order weights, and embedded fifth-/third-order
error estimator are adapted from Jacob Williams' modern Fortran DOP853, which
retains Ernst Hairer's original permissive DOP853 notice. SciFort continues to
own the adaptive driver, initial-step selection, result/storage types, events,
and interpolation layer.

Against installed SciPy 1.17.0, a harmonic oscillator on `[0,pi]` with
`rtol=1e-10` and `atol=1e-12` produced the same 13 accepted nodes and 146 RHS
evaluations. The final state was `(-0.9999999999973785,
-2.265526655165218e-11)`, agreeing with SciPy at roughly roundoff in the first
component and about `1.3e-15` absolute in the second. On a mildly nonlinear Van
der Pol system over `[0,20]` with `rtol=1e-9`, `atol=1e-12`, both SciFort and
SciPy used 710 RHS evaluations and 54 stored nodes; SciFort's final state
`(0.6101970163495973, 0.3912935547253569)` agrees with the SciPy 1.17.0
reference `(0.6101970163495998, 0.3912935547253552)` to about `3e-15`.

The permanent DOP853 regressions at this checkpoint cover exponential growth,
the harmonic oscillator with `t_eval`, and the nonlinear Van der Pol problem.
This initial checkpoint still used the common cubic-Hermite interpolation layer;
the subsequent dense-output/event batch below replaces that DOP853 limitation
with the native seventh-order interpolant.

GNU Fortran 14.2.0 compiled all **240/240** library sources under strict
Fortran 2018 checking and all **42/42** shared test-support/reference
modules. A separate optimized clean build passed all **138/138**
program-bearing repository tests. The checked build passed all **6/6**
`test_integrate*.f90` programs and all **71/71** `test_special*.f90`
programs with `-fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace`.
The external Quadpack2 dependency was again represented only by an
API-compatible stand-in outside the repository for whole-project linking;
it is unrelated to DOP853 and is not packaged. `fpm` and `fprettify` remain
unavailable in this environment.

## DOP853 dense output and multi-event validation (2026-09-28)

This batch adds the native Hairer/Wanner DOP853 seventh-order dense interpolant.
For each accepted step that requires interpolation, SciFort computes the three
additional DOP853 dense stages and seven coefficient vectors; endpoint-only
DOP853 runs do not pay or store that extra work. The same interpolant is used
for `t_eval`, retained dense solutions, and DOP853 event root location. The
new `ivp_event_set` interface also supports several event functions through one
vector callback, with per-event terminal and direction controls and a flattened
`event_indices` result array.

Against installed SciPy 1.17.0, integrating `y'=y` on `[0,5]` with DOP853,
`rtol=1e-11`, `atol=1e-13`, and retained dense output gives the same **332** RHS
evaluations and **23** stored nodes. At query points 0.123, 0.777, 1.234, 2.345,
and 4.321, the maximum SciFort-versus-SciPy absolute difference is
`2.84e-14`, and the maximum error scaled by `max(1,abs(reference))` is
`4.08e-16`. A two-event constant-velocity test locates roots at approximately
1 and 2 in the correct order, records indices `[1,2]`, and stops at the second
terminal event. SciPy locates the same roots; SciFort reports one additional RHS
evaluation in that event-only case because it recomputes the derivative at the
terminal root before storing the truncated endpoint.

A clean GNU Fortran 14.2.0 build compiled all **240/240** library sources and
all **42/42** shared test-support/reference modules. The optimized build passed
all **138/138** program-bearing repository tests. The checked build passed all
**6/6** `test_integrate*.f90` programs and all **71/71** `test_special*.f90`
programs with `-fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace`. As in
the preceding integration checkpoints, an API-compatible QUADPACK stand-in was
used outside the repository solely because the pinned external Quadpack2
dependency and `fpm` are unavailable in this environment; it is unrelated to
the DOP853 implementation and is not packaged. `fprettify` is also unavailable.

## BDF initial-value solver validation (2026-09-29)

This batch adds real-system `method='BDF'` to `solve_ivp`. The variable-order
1--5 quasi-constant-step BDF/NDF recurrence, Newton convergence logic, dense
output, and dense finite-difference Jacobian strategy are adapted from SciPy
1.17.0 `scipy/integrate._ivp.bdf` and `common` under BSD-3-Clause. SciFort adds
a Fortran Jacobian callback, compact partial-pivoting dense LU kernels, and
integration with the existing result, `t_eval`, dense-solution, and event APIs.

Against installed SciPy 1.17.0, `y'=y` on `[0,1]` with `rtol=1e-8` and
`atol=1e-10` gives the same `nfev=94`, `njev=1`, and `nlu=12`. Values at
`[0, .25, .5, .75, 1]` and dense queries at `.123` and `.777` agree with the
SciPy references to about binary64 roundoff. A stiff scalar problem with exact
solution `cos(t)`, stiffness coefficient 1000, and an analytic Jacobian gives
the same `nfev=54`, `njev=2`, and `nlu=9` as SciPy and agrees at the sampled
points to roundoff. A two-component diagonal stiff system, using the numerical
Jacobian rather than a callback, likewise matches SciPy's `nfev=332`, `njev=1`,
and `nlu=34` and its sampled states to roundoff. Backward integration of `y'=y`
from 1 to 0 gives the same `nfev=92`, `njev=1`, and `nlu=12` as SciPy, with
matching sampled and dense
values to roundoff. A terminal BDF event for `y'=y`, `y=2` is located at
`0.6931471709305568`, matching SciPy's `0.6931471709305571`; the retained BDF
polynomial remains usable before the truncated terminal point.

GNU Fortran 14.2.0 compiled the BDF/IVP dependency set in Fortran 2018 mode.
The complete `test_integrate_ivp` program, including the new BDF regressions,
passed both an optimized build and a checked build with `-fcheck=all` and
`-ffpe-trap=invalid,zero,overflow`. `fpm` and `fprettify` are not installed in
this execution environment, so `fpm test` and automatic formatting could not
be run at this checkpoint.

## Radau initial-value solver validation (2026-09-29)

This batch adds real-system `method='Radau'` to `solve_ivp`. The fifth-order
three-stage Radau IIA collocation method, embedded third-order error estimator,
transformed Newton systems, step-size predictor, cubic collocation dense output,
and dense finite-difference Jacobian strategy are adapted from SciPy 1.17.0
`scipy/integrate._ivp.radau` and `common` under BSD-3-Clause. SciFort adds the
Fortran Jacobian callback, compact partial-pivoting real and complex dense LU
kernels, and integration with the existing result, `t_eval`, dense-solution,
backward-integration, and event APIs.

Against installed SciPy 1.17.0, `y'=y` on `[0,1]` with `rtol=1e-8` and
`atol=1e-10` gives the same `nfev=239`, `njev=2`, and `nlu=6`; values at
`[0, .25, .5, .75, 1]` and dense queries at `.123`, `.321`, and `.9` agree to
binary64 roundoff. The reverse integration from 1 to 0 gives the same work
counts and matching sampled/dense values. A stiff scalar problem with exact
solution `cos(t)`, stiffness coefficient 1000, and an analytic Jacobian gives
the same `nfev=84`, `njev=1`, and `nlu=16`. A two-component diagonal stiff
system using the numerical Jacobian matches SciPy's `nfev=770`, `njev=2`, and
`nlu=40` and its sampled states to roundoff. A terminal event for `y'=y`, `y=2`
is located at `0.6931471805888532` versus SciPy's `0.6931471805888538`, while
the retained Radau cubic remains usable before the truncated terminal point.

GNU Fortran 14.2.0 compiled the Radau/IVP dependency set in Fortran 2018 mode.
The complete `test_integrate_ivp` program, including the new Radau and existing
BDF regressions, passed both an optimized `-O2` build and a checked build with
`-fcheck=all -ffpe-trap=invalid,zero,overflow`. The unchanged public
`scifort_integrate` umbrella also compiled and passed the IVP test using
compile-time stand-ins outside the repository for unrelated integration modules.
`fpm` and `fprettify` are not installed in this execution environment, so
`fpm test` and automatic formatting could not be run.

## Initial LSODA compatibility validation (2026-09-29)

This batch adds a real-system `method='LSODA'` compatibility path without importing
ODEPACK source. The nonstiff branch uses adaptive RK45 startup to establish an
equal-step history, followed by an Adams-Bashforth 4 predictor and Adams-Moulton 5
corrector. Periodic dense-Jacobian norm estimates trigger a switch to SciFort's
existing variable-order BDF/NDF implementation; low stiffness in the BDF branch can
restart the Adams path. BDF intervals retain their native backward-difference dense
polynomial, while nonstiff intervals use endpoint Hermite interpolation.

Installed SciPy 1.17.0 was used to establish reference behavior, not implementation
code. For `y'=y` on `[0,1]` with `rtol=1e-8` and `atol=1e-10`, SciFort's sampled
values at `[0,.25,.5,.75,1]` agree with `exp(t)` within about `4e-9` relative error;
retained dense values at `.123` and `.777` agree within `4e-8`. The nonstiff case
performs no LU factorizations. For the stiff scalar problem
`y'=-1000*(y-cos(t))-sin(t)` with its analytic Jacobian, SciFort automatically
switches to BDF (`nlu > 0`) and samples the exact `cos(t)` solution within
`3e-7` absolute error over `[0,1]`. A terminal event at `y=2` for `y'=y` is located
at `0.6931471826737028`, within `3e-9` of `log(2)`, and retained dense output before
the terminal point remains valid. A two-component diagonal stiff system also passes
with a finite-difference Jacobian, exercising automatic switching without `jac=`.

GNU Fortran 14.2.0 compiled the IVP dependency set in Fortran 2018 mode. The complete
`test_integrate_ivp` program, including all existing RK23/RK45/DOP853/Radau/BDF and
new LSODA regressions, passed through a compile-time `scifort_integrate` shim that
reexports the unchanged IVP public API. The shim is outside the repository and is
not packaged; it is used only because the pinned external Quadpack2 dependency and
`fpm` are unavailable in this execution environment. `fpm test` and `fprettify`
therefore could not be run. At this initial checkpoint, order-12 Adams work, the
original switch heuristics/work counters, banded Jacobians, and `min_step` were
documented as follow-on work; later checkpoints below address the first, third, and
fourth items.

## Complex explicit `solve_ivp` validation (2026-09-29)

This batch overloads `solve_ivp` for complex-valued systems with `RK23`, `RK45`,
and `DOP853`. The complex steppers use the same real Runge-Kutta tableaux and
dense-output formulas as the established real implementations, with complex
stage/state vectors and RMS error norms based on complex magnitudes. Scalar and
vector event callbacks still return real values, while located event states are
stored as complex values. The real-system API is unchanged; the public
`scifort_integrate` generic resolves from the type of `y0`.

Reference values were generated with installed SciPy 1.17.0 for
`y'=(0.5+2i)y`, `y(0)=1+0.25i` on `[0,1]`. With the same tolerances, SciFort
matches SciPy's RHS counts exactly: `nfev=881` for RK23, `218` for RK45, and
`152` for DOP853. Values at `[0,.25,.5,.75,1]` and dense queries at
`[.123,.321,.9]` agree with the SciPy references within `5e-13` in complex
absolute error (typically much closer). Reverse-time DOP853 integration from
1 to 0 also matches SciPy's `nfev=152`, with sampled and dense values within
`5e-12`. A rotating solution `y'=i y` locates the falling real-part zero near
`pi/2` and stores the corresponding complex event state near `i`. A second
regression terminates on the second occurrence of the same real-valued event,
near `3*pi/2`, exercising complex event-state storage together with the
existing positive-integer occurrence-limit logic. A two-function vector event
set finds the expected roots near `pi/6` and `pi/3` with event indices 1 and 2
and the corresponding complex states. Complex `Radau` is rejected explicitly;
complex BDF is validated separately below.

GNU Fortran 14.2.0 compiled the real and complex IVP dependency set in Fortran
2018 mode with `-Werror=implicit-interface`. The new
`test_integrate_ivp_complex` passed both an optimized `-O2` build and a checked
build with `-fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace`. The
existing real `test_integrate_ivp` also passed through the new overloaded public
`scifort_integrate%solve_ivp` generic. Compile-time stand-ins outside the
repository were used only for unrelated integration modules because the pinned
external Quadpack2 dependency and `fpm` are unavailable in this environment.
`fprettify` is also unavailable, so `fpm test` and automatic formatting could
not be run.

## Complex BDF `solve_ivp` validation (2026-09-29)

This batch adds complex-system `method='BDF'` to the existing complex
`solve_ivp` overload. The variable-order 1--5 BDF/NDF recurrence, Newton
convergence logic, step/order adaptation, and native dense-output polynomial
are adapted from SciPy 1.17.0. The dense numerical Jacobian follows SciPy's
complex convention: finite-difference perturbations are made along the real
axis and the sign of the real part of the right-hand side selects the
perturbation direction. SciFort supplies a compact complex partial-pivoting LU
solver and also accepts an optional analytic complex `jac(t,y,J)` callback.

For `y'=(0.5+2i)y`, `y(0)=1+0.25i` on `[0,1]` with `rtol=1e-9` and
`atol=1e-12`, SciFort and installed SciPy 1.17.0 both report `nfev=174`,
`njev=1`, and `nlu=19`. Values at `[0,.25,.5,.75,1]` and retained-dense
queries at `[.123,.321,.9]` agree to binary64 roundoff. Supplying the exact
analytic Jacobian gives the same `174/1/19` work counts. A two-component
complex diagonal stiff system using the finite-difference Jacobian matches
SciPy's `nfev=562`, `njev=1`, and `nlu=51`, with sampled values agreeing within
`5e-12`. A second two-component system with a fully coupled complex matrix
exercises off-diagonal finite differences and complex LU solves; it matches
SciPy's `nfev=536`, `njev=1`, and `nlu=49`, with sampled and dense states within
`5e-12`. Reverse integration of the scalar problem matches SciPy's
`176/1/20` work counts and sampled/dense states within `5e-12`. A rotating
solution `y'=i y` also exercises complex BDF event interpolation and locates
the falling real-part zero at `pi/2` within `3e-10` absolute error.

GNU Fortran 14.2.0 compiled the complete real/complex IVP dependency set in
Fortran 2018 mode through the public `scifort_integrate` umbrella. Both
`test_integrate_ivp_complex` and the existing `test_integrate_ivp` passed an
optimized `-O2` build and a checked build with `-fcheck=all
-ffpe-trap=invalid,zero,overflow -fbacktrace -Wimplicit-interface
-Werror=implicit-interface`. Compile-time stand-ins outside the repository were
used only for unrelated integration modules. `fpm` and `fprettify` are not
installed in this execution environment, so `fpm test` and automatic formatting
could not be run.

## LSODA banded Jacobian and `min_step` validation (2026-09-29)

This batch extends the independent SciFort LSODA compatibility path with the public
SciPy `lband`, `uband`, and `min_step` controls. Supplying either bandwidth selects
banded Jacobian semantics; the omitted bandwidth defaults to zero. Analytic
Jacobians use packed storage with `uband+lband+1` rows, while finite-difference
Jacobians perturb columns separated by the total band width in one RHS call and
retain zeros outside the declared band. The BDF Newton solve still uses SciFort's
dense LU factorization after unpacking, so this change reduces Jacobian callback
work but does not yet add a dedicated banded linear solver.

New regressions integrate diagonal stiff systems with exact solution `cos(t)`. The
finite-difference regression compares the same eight-state system with dense and
zero-band Jacobian declarations and verifies that the banded path uses fewer actual
RHS callback invocations while still switching to BDF and meeting the solution
tolerance. A second six-state regression supplies an analytic packed diagonal
Jacobian and rejects an unexpected dense callback shape. A third regression uses a
tridiagonal stiff system with both finite-difference and analytic packed Jacobians,
verifying the `lband=1`, `uband=1` packing convention against a known `cos(t)`
solution. The `min_step` regression uses a deliberately oversized first trial with
tight tolerances and verifies that a required retry below the requested floor fails
cleanly; negative `min_step` and out-of-range bandwidths return
`ivp_status_invalid_argument`.

GNU Fortran 14.2.0 passes the complete `test_integrate_ivp` program in Fortran 2018
mode with `-fcheck=all -ffpe-trap=invalid,zero,overflow` and with implicit interfaces
promoted to errors. The same integration regression program also passes an optimized
`-O2` build. As in the preceding IVP checkpoints, a compile-time umbrella shim is
used only because `fpm` and the pinned external Quadpack2 dependency are unavailable
in this execution environment; the shim is not part of the repository.

## Higher-order Adams LSODA validation (2026-09-29)

The nonstiff LSODA branch now generalizes the original AB4/AM5 pair without importing
coefficient tables. For a history of `q` equally spaced derivatives, `q=4,...,11`,
SciFort constructs the Lagrange interpolating basis at normalized nodes
`0,-1,...,-(q-1)` and integrates each basis polynomial on `[0,1]` to obtain the
Adams-Bashforth predictor weights. It repeats the construction at nodes
`1,0,-1,...,-(q-1)` for the Adams-Moulton corrector, giving corrector orders 5
through 12. At that checkpoint the order was raised as equal-step history became
available and step-size changes rebuilt the history through RK45. The subsequent
variable-step checkpoint below removes that history-reset restriction while still not
claiming ODEPACK's Nordsieck representation or exact order-selection logic.

A direct internal-stepper regression holds `max_step=0.01` for `y'=y`, tracks the
public `lsoda_state%history_count`, and verifies that the history reaches 11 entries,
thereby exercising the order-12 Adams-Moulton path. It remains on the nonstiff branch
with no LU factorizations and agrees with `exp(t)` within `2e-8` relative error after
20 accepted steps. This regression also catches a subtle history-reset bug: when a
requested step increase is clipped by `max_step` to the existing step size, history
is now preserved instead of being needlessly discarded.


## Variable-step Adams LSODA validation (2026-09-29)

The variable-step LSODA checkpoint removes the equal-step restriction from the SciFort
Adams predictor/corrector without importing ODEPACK coefficient tables. Each
accepted derivative is retained together with its signed time offset from the current
endpoint. For every trial step, the Adams-Bashforth predictor and Adams-Moulton
corrector weights are regenerated by integrating the Lagrange basis on the actual
normalized nonuniform history nodes. Accepted step-size changes therefore preserve
usable history instead of forcing an RK45 rebuild, and `adams_order` is adapted
separately from the retained history depth.

At this checkpoint, the switching path was changed so ordinary nonstiff Adams work no
longer performed periodic Jacobian-norm probes. Accepted and rejected trials first used
a directional secant estimate of `h*||df||/||dy||`; an explicit dense or banded
Jacobian-norm estimate was still reserved as a fallback after repeated rejected trials.
The later roundoff/failure-recovery checkpoint below removes that rejected-step probe.
Once stiffness is detected,
the existing SciFort BDF/NDF branch and its analytic or finite-difference Jacobian
machinery are still used.

A new direct regression integrates `y'=y` from 0 to 1 starting with `h=0.001` and
`max_step=0.2`. It verifies that accepted step sizes change while the nonuniform
history is retained, that the history never falls below six entries after first
reaching six, that the Adams predictor order rises above the startup order, and that
no Jacobian or LU work occurs. Under the checked GNU Fortran 14.2.0 configuration
described above, the regression completed in 17 accepted calls with 59 RHS evaluations,
zero Jacobian evaluations, zero LU factorizations, and absolute final error about
`1.52e-8`.

The complete `test_integrate_ivp` program, including the existing stiff-switch,
banded-Jacobian, backward-integration, dense-output, event, and `min_step` LSODA
checks, also passes under the same checked GNU Fortran build through a temporary
compile-time `scifort_integrate` shim. `fpm test` could not be run in this environment
because the `fpm` executable is unavailable. This checkpoint did not yet claim
ODEPACK Nordsieck representation or ideal-step method-switch comparisons; the later
switch-comparison checkpoint below addresses the common finite-error comparison.

## ODEPACK-style LSODA stability and switch hysteresis validation (2026-09-29)

The next LSODA checkpoint adapts three method-selection details from the public-domain
ODEPACK `DSTODA` routine: the 12-entry Adams `SM1` stability envelope, the
20-successful-step delay between method-switch reconsiderations, and the rule that an
Adams-to-BDF switch is not considered above Adams method order 5. SciFort stores the
Adams-Bashforth predictor order, so `SM1` and the order-5 guard are indexed with the
corresponding Adams-Moulton method order, `predictor_order + 1`. SciFort continues to
use its nonuniform derivative/time history rather than ODEPACK's Nordsieck array.

For active Adams history, the already-computed predictor and corrected derivatives
provide a directional local Lipschitz estimate, `||f_new-f_predict|| /
||y_new-y_predict||`. This detects a hidden fast mode even when accepted-step secants
lie close to a slow manifold, without adding RHS or Jacobian evaluations. RK startup
and degenerate-correction cases retain the accepted-step secant estimate.

A new direct stiff-cosine regression initializes `h=1e-4`. It verifies zero Jacobian
and LU work through the first 20 accepted Adams steps, checks that the method-order-5
stability limit caps the active Adams proposal at approximately `0.35/1000`, observes
the BDF switch on the 21st accepted step, and sees the first BDF factorization on the
following call. The complete `test_integrate_ivp` program passes both the checked GNU
Fortran 14.2.0 build (`-fcheck=all`, floating-point invalid/zero/overflow traps, and
implicit interfaces promoted to errors) and an optimized `-O2` build. A smooth
`y'=y` benchmark remains at 96 RHS evaluations with zero Jacobian evaluations and zero
LU factorizations; the stiff cosine benchmark completes with 126 RHS evaluations, two
Jacobian evaluations, and eight LU factorizations. `fpm` and `fprettify` are not
available in this environment. The subsequent checkpoint below adds the common
finite-error ideal Adams/BDF step-size comparison and removes a redundant BDF endpoint
RHS call; Nordsieck step/order adaptation and exact counter parity remain follow-on
work.

## LSODA ideal-step switching and BDF endpoint derivative validation (2026-09-29)

This checkpoint extends the ODEPACK-derived method-selection logic without replacing
SciFort's independent variable-step Adams representation. For the common finite-error
switch path, SciFort now forms ideal next-step factors for the current and alternate
methods using the `CM1` and `CM2` error constants generated by ODEPACK `DCFODE` and the
ratio tests in `DSTODA`. Adams-to-BDF therefore requires the BDF ideal factor to be at
least five times the stability-limited Adams factor; BDF-to-Adams uses the corresponding
ODEPACK comparison. The implementation deliberately keeps SciFort's conservative
post-switch step initialization instead of directly rescaling a non-Nordsieck state by
ODEPACK's `RH` value.

The BDF compatibility path also removes one known source of RHS-counter inflation. The
converged BDF Newton equation already satisfies `d_corr = c*f_new - psi`, so
`bdf_advance` can return `f_new=(d_corr+psi)/c` to LSODA instead of reevaluating the RHS
at the accepted endpoint. A direct LSODA regression verifies that the first accepted
BDF step after the hysteresis-triggered switch uses two RHS evaluations rather than an
additional third endpoint call, and compares the implied derivative against a fresh
callback value to `2e-12` absolute tolerance.

With GNU Fortran 14.2.0, the complete `test_integrate_ivp` program passes the checked
Fortran 2018 build with `-fcheck=all -ffpe-trap=invalid,zero,overflow` and implicit
interfaces promoted to errors. The same suite also passes an optimized `-O2` build.
The smooth `y'=y` LSODA benchmark remains at 96 RHS evaluations, zero Jacobian
evaluations, and zero LU factorizations. The stiff cosine benchmark retains the same
accepted trajectory and eight LU factorizations while dropping from 126 to 104 counted
RHS evaluations because accepted BDF endpoints are no longer reevaluated. Exact
ODEPACK/SciPy counter parity is not claimed: the solver still has independent startup,
Nordsieck/order adaptation, Jacobian-refresh, and failure-recovery behavior. The
remaining switching gap is primarily ODEPACK's roundoff/failure branches together with
full Nordsieck step/order evolution.


## LSODA Adams Nordsieck predictor validation (2026-09-29)

This checkpoint introduces the first persistent Nordsieck representation in the
SciFort LSODA Adams path without replacing the existing order-4-and-higher
nonuniform-history Adams-Moulton corrector. After an accepted step has accumulated
usable history, SciFort interpolates the retained derivative samples at their actual
signed time offsets and converts that local derivative polynomial into scaled
solution derivatives. Column `j+1` therefore represents the local `j`-th solution
derivative scaled by the signed reference step to power `j` and divided by `j!`,
matching the `YH` convention documented by ODEPACK `DSTODA`.

When a later Adams trial uses a different step, the stored derivative columns are
rescaled in place by successive powers of the step ratio. The predictor is then
formed by the same Pascal-triangle shift used by `DSTODA`. SciFort still computes
the Adams-Moulton correction from its variable-node Lagrange integration weights and
reconstructs the Nordsieck state from accepted nonuniform derivative history; it does
not yet apply ODEPACK's `EL` correction directly to `YH`, nor does it use `IALTH` and
`TESCO` for exact ODEPACK order selection.

Normal history startup now uses Nordsieck Adams orders 1 through 3 rather than a
six-evaluation RK45 step. SciFort generates the low-order correction polynomial and
error-test constants with the public-domain ODEPACK `DCFODE` recurrence. The initial
Nordsieck-startup checkpoint used two fixed functional corrections; the later
functional-corrector checkpoint below replaces those with the `DSTODA` convergence
test while preserving an exact accepted endpoint derivative whenever genuine iteration
was required. If Nordsieck construction fails, the previous RK45 startup remains as a
guarded fallback. Rejected low-order startup trials do not use the accepted-step secant
as an immediate stiffness trigger, because explicit time dependence can make that
secant large even when the Jacobian is nonstiff. The later roundoff/failure-recovery
checkpoint below also removes the rejected-step Jacobian probe.

The direct constant-RHS regression now verifies that the first low-order startup step
is accepted in one attempt, extends derivative history from one to two entries, uses
one step RHS evaluation, and performs no Jacobian or LU work. A second direct regression
integrates the polynomial solution `y=t**2`. Once active higher-order Adams history is
available, it checks the Nordsieck value column, the scaled first derivative column
`h*y'`, and the scaled second-derivative column `h**2`. The test then forces the next
maximum step to one quarter of the stored reference step and rechecks the same
identities after the rescaled Nordsieck predictor is exercised. Because the derivative
is linear, these columns have exact polynomial reference values apart from binary64
rounding. A `solve_ivp` regression also integrates `y'=y` backward from 1 to 0 at
several `t_eval` points; it verifies the signed-step Nordsieck path against `exp(t)`
while requiring zero Jacobian evaluations and zero LU factorizations.

GNU Fortran 14.2.0 passes the complete `test_integrate_ivp` program in Fortran 2018
mode with `-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow` and
`-Werror=implicit-interface`. An optimized `-O2` build also passes. Relative to the
preceding checkpoint, the smooth `y'=y` benchmark drops from 96 to 78 RHS evaluations
with zero Jacobian evaluations and zero LU factorizations. The stiff cosine benchmark
drops from 104 to 95 RHS evaluations while retaining two Jacobian evaluations and
eight LU factorizations. Final errors remain within the established regression
bounds. `fpm` and `fprettify` are not installed in this execution environment, so
`fpm test` and the repository formatting command could not be run here.

## LSODA roundoff and repeated-failure recovery validation (2026-09-29)

This checkpoint adapts two additional control-flow details from the public-domain
ODEPACK `DSTODA` routine without changing SciFort's higher-order variable-node
Adams-Moulton corrector. Method-switch decisions remain an accepted-step operation:
rejected Adams error-test trials no longer probe a Jacobian or enter BDF. On accepted
Adams steps, SciFort retains the most recent nonzero correction-based Lipschitz
estimate as a `PDLAST` analogue for the stability envelope. When the local-error
information is too small to distinguish reliably from roundoff, the Adams-to-BDF
branch follows the `DSTODA` fallback and switches only if the Adams proposal is
already stability-limited.

`DSTODA` enters its special repeated-error-test recovery after `KFLAG <= -3`, reduces
the step by 0.1, resets to order 1, recomputes the first derivative, and retries.
SciFort uses the same restart-and-factor-of-ten pattern but conservatively starts it
on the fourth consecutive rejected trial. Applying it on the third rejection
regressed an established LSODA dense-output tolerance in SciFort's hybrid
variable-node controller. Because a rejected trial has not changed the retry origin,
SciFort already holds a valid derivative there and reuses it rather than issuing an
extra RHS callback.

A direct forced-failure regression first establishes normal Adams history for
`y'=y`, then requests an intentionally oversized step under `rtol=1e-14` and
`atol=1e-16`. The call exercises eight rejected Adams trials, recovers through the
low-order restart, finishes successfully, leaves two derivative-history entries, and
performs zero Jacobian evaluations and zero LU factorizations. The accepted value
matches `exp(t)` within the regression's `2e-8` relative tolerance.

A second direct regression manufactures the roundoff-polluted switch branch with the
constant system `y'=1`. It seeds a preceding Adams correction-rate estimate large
enough to make the proposed step stability-limited, forces the switch countdown to be
eligible, and verifies that the accepted Adams step selects BDF. The switch performs
one analytic Jacobian evaluation but no LU factorization until the following BDF
step; the retained Adams rate is cleared on the switch and both accepted states match
`y=t` to binary64 roundoff.

GNU Fortran 14.2.0 passes the complete `test_integrate_ivp` program in Fortran 2018
mode with `-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow` and
`-Werror=implicit-interface`. An optimized `-O2` build also passes. The normal smooth
benchmark remains at 78 RHS evaluations, zero Jacobian evaluations, and zero LU
factorizations. The stiff cosine benchmark remains at 95 RHS evaluations, two
Jacobian evaluations, and eight LU factorizations, so the added branches do not alter
the established normal trajectories. `fpm` and `fprettify` are not installed in this
execution environment, so `fpm test` and automatic formatting could not be run.

## LSODA low-order functional-corrector convergence validation (2026-09-29)

The low-order Nordsieck Adams startup no longer takes two functional corrections
unconditionally. The public-domain ODEPACK `DSTODA` corrector loop was inspected and
its convergence tests were adapted to SciFort's scalar tolerance interface. The first
weighted correction change may terminate at the `100*PNORM*UROUND` roundoff threshold.
Otherwise the startup forces another functional evaluation, estimates the contraction
rate, updates the `CRATE` analogue from its ODEPACK initial value 0.7, and tests the
scaled correction with `CONIT=0.5/(q+2)` and `TESCO(2,q)`. At most three functional
evaluations are attempted, and loss of contraction forces the outer Adams step to be
rejected with ODEPACK's factor-of-four retry reduction rather than allowing a
misleading local-error estimate to accept it.

SciFort still differs deliberately from pure `DSTODA` after a genuinely iterated
startup correction. ODEPACK can retain the converged derivative implicitly in its
persistent Nordsieck `YH` state, while SciFort's order-4-and-higher hybrid currently
retains explicit nonuniform endpoint derivatives. SciFort therefore reevaluates the
RHS at the accepted endpoint before extending that history. If the correction converges
at roundoff on the first evaluation, no refresh is needed. For the constant system
`y'=1`, the direct startup regression consequently uses one step RHS evaluation (two
cumulative evaluations including the initial derivative), rather than the preceding
checkpoint's three step evaluations.

A second regression uses `y'=100*y`, `first_step=0.1`, `min_step=0.09`, and deliberately
loose unit tolerances. With only a fixed number of corrections, the rapidly expanding
functional iteration can produce a tolerance-scaled accumulated correction that appears
acceptable even though the iterates are diverging. The new contraction test detects the
loss of convergence after two step RHS evaluations, rejects the trial, and then reports
the required retry as blocked by `min_step`; no Jacobian evaluation or LU factorization
is performed.

GNU Fortran 14.2.0 passes the complete `test_integrate_ivp` program in Fortran 2018
mode with `-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow` and
`-Werror=implicit-interface`. The smooth `y'=y` benchmark now uses 77 RHS evaluations,
zero Jacobian evaluations, and zero LU factorizations, compared with 78/0/0 in the
preceding checkpoint. The stiff cosine benchmark remains at 95 RHS evaluations, two
Jacobian evaluations, and eight LU factorizations. An optimized `-O2` build of the
complete IVP regression program also passes.


## LSODA persistent Nordsieck correction/order validation (2026-09-29)

This checkpoint replaces the remaining successful-step hybrid Adams path with a
persistent Nordsieck update through Adams method orders 1--12. The accepted scaled
solution-derivative array is now the active Adams state: trial-step changes power-rescale
its derivative columns, the Pascal shift predicts the new endpoint, functional iteration
forms the accumulated correction, and the accepted step applies the ODEPACK `EL`
correction directly to the predicted array. Nonuniform endpoint-derivative/time history
is still retained for state reconstruction and conservative recovery, but a successful
high-order step no longer rebuilds the active Nordsieck state from that history.

Successful-step order and step reconsideration now adapts the ODEPACK `IALTH`/`TESCO`
logic. `IALTH` suppresses unnecessary changes for the appropriate number of successful
steps; when reconsideration is allowed, SciFort compares the TESCO-scaled lower,
current, and higher order candidates, applies the `SM1` Adams stability limits, requires
a roughly 10 percent useful step change before changing order/step, and constructs the
additional scaled derivative when increasing order. A direct regression verifies that
successive order changes respect the `IALTH` delay and that the cooldown is reset after
an order change. A separate direct regression seeds the exact order-12 Nordsieck state
for `y'=y`, takes one genuine order-12 Adams step, and checks the accepted value,
derivative, and negative dense-order marker. This preserves explicit coverage of the
maximum Adams order without forcing the adaptive controller to make an artificial
order increase. The high-order regression also checks that accepted nonstiff LSODA
steps return a native Nordsieck dense polynomial instead of relying on the previous
Hermite fallback.

The functional corrector now applies at all Adams orders. On convergence, the derivative
from the last functional iteration is retained directly, as in the persistent ODEPACK
corrector, rather than issuing an extra RHS evaluation at the corrected endpoint. This
reduces the smooth `y'=y` benchmark from 77 to 72 RHS evaluations while retaining zero
Jacobian evaluations and zero LU factorizations. The stiff cosine benchmark changes from
95/2/8 to 94 RHS evaluations, two Jacobian evaluations, and eight LU factorizations.
The smooth dense-output errors at `t=[0.25,0.5,0.75,1]` remain about `1e-10` to `1e-8`.

The variable-step exponential regression now uses a `4e-7` relative endpoint tolerance.
For the same `rtol=1e-7`, `atol=1e-9`, `first_step=0.001`, and `max_step=0.2`, installed
SciPy 1.17.0 LSODA gives relative endpoint error about `3.07e-7`; SciFort gives about
`2.86e-7`. The preceding `2e-7` test bound was therefore tighter than the SciPy
compatibility reference even though both solvers satisfy the requested integration
tolerances.

GNU Fortran 14.2.0 passes the complete direct `test_integrate_ivp` program in Fortran
2018 mode with `-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow` and
`-Werror=implicit-interface`. The complete IVP regression program also passes an
optimized `-O2` build. `fpm` and `fprettify` are not installed in this execution
environment, so the repository-level `fpm test` and automatic formatting commands could
not be run. Subsequent LSODA batches completed the documented rejected-step
order-selection, BDF-to-Adams transfer/Jacobian-refresh, startup work-count,
and native stiff-controller milestones.

## Componentwise `solve_ivp` absolute-tolerance validation (2026-09-30)

Real and complex `solve_ivp` now accept a rank-one `atol` array with one value per
state component. A two-scale decay regression, `y1'=-y1`, `y2'=-10*y2`, on `[0,2]`
uses `rtol=1e-6` and compares `atol=[1e-9,1e-3]` with scalar `atol=1e-9`. For the
real methods, SciFort reports relaxed/tight RHS counts of **170/803** (RK23),
**74/302** (RK45), **74/194** (DOP853), **155/722** (Radau), **80/268** (BDF),
and **65/201** (LSODA). Installed SciPy 1.17.0 gives the same first five pairs and
**72/203** for LSODA. In every case the slowly decaying first component remains near
`exp(-2)` while the deliberately loose second-component tolerance materially reduces
work.

For the complex regression `y1'=(-1+i)y1`, `y2'=(-10+2i)y2`, SciFort and installed
SciPy 1.17.0 both report **86/308** RHS evaluations for RK45 and **92/272** for BDF
under the same relaxed/tight tolerance comparison. The regression suite also checks
that a vector with the wrong length, a rank-two tolerance array, and a vector with a
negative entry are rejected with `ivp_status_invalid_argument`; the complex driver
checks the wrong-length case as well. Existing scalar-tolerance IVP regressions remain
unchanged.

GNU Fortran 14.2.0 passes the complete real and complex IVP regression programs in
Fortran 2018 mode with `-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow` and
`-Werror=implicit-interface`. During this checked run, an existing LSODA stiff-path
expression `max_step/h_abs` was found to overflow when `max_step` has its default
`huge()` value. The same failure was reproduced from the untouched uploaded baseline.
The cap is now tested with an overflow-safe comparison before forming the quotient, and
both IVP regression programs pass with floating-point traps enabled. An independent
`-O2 -Werror=implicit-interface` build also passes the real/complex IVP programs and
the LSODA stiff-switch, transfer, and work-count regressions. `fpm` and `fprettify`
are not installed in this execution environment, so repository-level `fpm test` and
automatic formatting could not be run.

## Constant and sparsity-structured `solve_ivp` Jacobian validation (2026-09-30)

Real `BDF` and `Radau` now accept a constant dense `jac_matrix=` and a logical
`jac_sparsity=` finite-difference pattern; complex `BDF` accepts the same forms with
a complex constant matrix. The primary real regression is an eight-component diagonal
stiff problem with exact solution `y_i(t)=cos(t)` and diagonal rates from 500 through
1200 on `[0,0.5]`, using `rtol=2e-6` and `atol=1e-9`. With a constant matrix, SciFort
uses **38/0/7** (`nfev/njev/nlu`) for BDF and **44/0/12** for Radau. Installed SciPy
1.17.0 gives the same reported counts for the corresponding constant-array Jacobian
calls. The zero `njev` count verifies that a constant matrix is stored rather than
re-evaluated as a callback.

For numerical differentiation of the same eight-component diagonal problem, a dense
Jacobian requires 46 total RHS callback invocations in SciFort BDF and 52 in Radau.
Supplying an identity `jac_sparsity` pattern reduces those totals to **39** and **45**,
respectively: the eight diagonal columns are perturbed in one grouped RHS call instead
of eight independent calls. The solver-reported `nfev`, `njev`, and `nlu` are unchanged,
as expected because Jacobian finite-difference probes are not included in BDF/Radau
`nfev`. Installed SciPy 1.17.0 uses 47/40 total RHS calls for dense/sparse BDF and
52/45 for dense/sparse Radau on the equivalent problem; its reported solver counts are
38/1/7 and 44/1/12, respectively. The one-call dense-BDF difference is due to SciFort's
pre-existing numerical-Jacobian base-value reuse, while the sparsity grouping saves the
same seven perturbation calls.

For the existing two-component complex stiff BDF regression, the constant complex
matrix preserves **562/0/51** reported work counts and the finite-difference identity
pattern reduces total RHS callbacks from **564** to **563** while leaving the reported
562/1/51 counts unchanged. Tests also reject wrong matrix/pattern shapes, unsupported
explicit/LSODA use, and simultaneous `jac=` plus `jac_matrix=`.

GNU Fortran 14.2.0 passes the complete real and complex IVP regression programs in
Fortran 2018 mode with `-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow` and
`-Werror=implicit-interface`. The direct LSODA rejection, stiff-switch, transfer, and
work-count regression programs also pass with the checked build, confirming that the
new standalone BDF Jacobian paths leave LSODA behavior unchanged. An independent
`-O2 -Werror=implicit-interface` build passes the same real/complex IVP programs and
all four LSODA regressions. The IVP dependency slice was used because the execution
environment does not provide FPM. `fpm` and `fprettify` are not installed, so the
repository-level `fpm test` and automatic formatting commands could not be run.


## Vectorized `solve_ivp` RHS validation (2026-09-30)

Real `BDF` and `Radau` and complex `BDF` now accept an optional rank-two
`fun_vectorized(t,y,dydt)` callback for finite-difference Jacobian probes. The
ordinary rank-one RHS remains the callback used by accepted/trial integration
steps. On the existing eight-component diagonal stiff regression, both real
implicit methods invoke the vectorized callback with all eight perturbed columns
in one dense-Jacobian batch. With the identity `jac_sparsity` pattern, the greedy
coloring has one group and the maximum batch width is one. The vectorized and
scalar finite-difference runs have identical reported `nfev` and `njev`, while
the instrumented scalar callback count equals reported `nfev` in the vectorized
case; this verifies that hidden Jacobian probes are routed through the batch
callback and retain SciPy-style work-count accounting. Final states agree within
the existing stiff-problem tolerances.

The two-component complex BDF regression gives the corresponding dense batch
width of two and sparsity-group width of one, again preserving reported
`nfev`/`njev` and the final complex solution. Real and complex explicit RK45
calls that supply `fun_vectorized=` are rejected with
`ivp_status_invalid_argument` instead of silently ignoring an inapplicable
option. Constant `jac_matrix=` and analytic `jac=` paths remain independent of
the batched finite-difference callback.

GNU Fortran 14.2.0 passes the complete real and complex IVP regression programs
in Fortran 2018 mode with `-O0 -g -fcheck=all
-ffpe-trap=invalid,zero,overflow -Werror=implicit-interface`. The direct LSODA
rejection, stiff-switch, transfer, and work-count regression programs pass under
the same checked build, confirming that LSODA's separate Jacobian/controller
path is unchanged. An independent `-O2 -Werror=implicit-interface` build also
passes the real/complex IVP programs and all four LSODA regressions. The IVP
dependency slice was used because this execution environment does not provide
FPM. `fpm` and `fprettify` are not installed, so repository-level `fpm test` and
automatic formatting could not be run.

## Constant CSC `solve_ivp` Jacobian validation (2026-09-30)

SciFort now exposes minimal one-based compressed-sparse-column containers for
constant implicit-solver Jacobians. Real `BDF` and `Radau` accept
`ivp_sparse_matrix` through `jac_sparse_matrix=`, and complex `BDF` accepts the
corresponding `ivp_complex_sparse_matrix`. The `ivp_csc_matrix` constructor
stores `indptr`, `indices`, and `data` directly in CSC form; `solve_ivp` then
validates the structure and values. Row indices are one-based and strictly
increasing within each column. Solver
state does not allocate the ordinary dense Jacobian when this path is active.
The constant-CSC path now also avoids the ordinary dense LU workspace: shifted
Newton matrices are factored directly into sparse row L/U storage with partial
row pivoting and dynamically inserted fill-in.

On the existing eight-component diagonal stiff problem with exact solution
`y_i(t)=cos(t)`, the CSC constant Jacobian preserves **38/0/7**
(`nfev/njev/nlu`) for BDF and **44/0/12** for Radau, matching the already
validated dense constant-matrix path. Instrumented RHS counts equal reported
`nfev`, confirming that the constant sparse matrix causes no hidden Jacobian
probes. The two-component complex BDF regression likewise preserves
**562/0/51** and the same final complex state as the dense constant matrix.

The permanent regression suite additionally rejects malformed CSC column
pointers, explicit-method and LSODA use of `jac_sparse_matrix=`, and simultaneous
dense and CSC constant Jacobians. The solver also checks finite values and exact
`(size(y0),size(y0))` shape before initialization. LSODA retains its separate
dense/packed `jac=` plus `lband`/`uband` interface.

GNU Fortran 14.2.0 passes the complete real and complex IVP regression programs
in Fortran 2018 mode with `-O0 -g -fcheck=all
-ffpe-trap=invalid,zero,overflow -Werror=implicit-interface`. The direct LSODA
rejection, stiff-switch, transfer, and work-count programs also pass under the
same checked build. An independent `-O2 -Werror=implicit-interface` build passes
the real/complex IVP programs and all four LSODA regressions. The IVP dependency
slice was used because this execution environment does not provide FPM. `fpm`
and `fprettify` are not installed, so repository-level `fpm test` and automatic
formatting could not be run.

The sparse factorization itself is covered independently by
`test/test_integrate_sparse_lu.f90`. It checks a 40-by-40 tridiagonal real
system, a matrix that requires a first-column row interchange, Radau-style
complex diagonal shifts applied to a real CSC Jacobian, a fully complex CSC
factorization, and singular-pivot detection. The tridiagonal regression also
requires the stored L/U nonzero count to remain below one quarter of dense
storage, guarding against accidental dense fill or a dense-workspace fallback.

## Dynamic and numerical CSC `solve_ivp` Jacobian validation (2026-09-30)

Real `BDF` and `Radau` and complex `BDF` now accept state-dependent sparse
Jacobian callbacks through `jac_sparse=`. The callback returns the existing
one-based CSC container, is counted in `njev`, and is validated for square
shape, finite values, and valid CSC structure before each refreshed Jacobian is
used. Permanent regressions integrate the existing diagonal real stiff problem
with both real implicit methods and the two-component complex stiff problem with
complex BDF. In every case the callback invocation count equals reported
`njev`, reported `nfev` equals the instrumented scalar stepping-RHS count, the
native sparse LU performs nonzero factorization work, and final states retain
the established tolerances. Additional regressions reject explicit-method and
LSODA use, simultaneous dense and sparse callbacks, and malformed initial CSC
callback output.

The `jac_sparsity=` numerical path now retains sparse storage as well as sparse
probe grouping. BDF/Radau initialize a CSC matrix directly from the logical
pattern, write finite-difference values only into its declared structural
entries, and factor shifted Newton matrices through `scifort_integrate_sparse_lu`.
Direct state-level tests require `jacobian_sparse=.true.`, a valid CSC numerical
Jacobian with the expected nonzero count, and no allocated dense Jacobian array
for real BDF, real Radau, and complex BDF. The existing vectorized-RHS tests
continue to verify that coloring determines finite-difference batch width and
that hidden Jacobian probes are excluded from `nfev`.

While exercising rejected sparse Radau trials, this batch also closed a latent
fallback in the embedded-error estimator: its second solve now dispatches to the
same sparse real factorization instead of the zero-sized dense LU workspace.
Constant CSC behavior and LSODA's separate dense/packed Jacobian controller path
are unchanged.

GNU Fortran 14.2.0 passes `test_integrate_ivp`,
`test_integrate_ivp_complex`, `test_integrate_sparse_lu`, and the direct LSODA
rejection, stiff-switch, transfer, and work-count programs in Fortran 2018 mode
with `-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow` and
`-Werror=implicit-interface`. An independent `-O2 -Werror=implicit-interface`
build passes the same seven programs. The IVP dependency slice was used because
this execution environment does not provide FPM. `fpm` and `fprettify` are not
installed, so repository-level `fpm test` and automatic formatting could not be
run.
## Legacy `odeint` facade validation (2026-09-30)

The public `scifort_integrate` umbrella now exports `odeint`,
`odeint_tfirst`, and `odeint_result`. The facade delegates integration to the
validated real LSODA `solve_ivp` implementation while presenting the historical
`func(y,t)` callback order, time-major output, repeated monotone output times,
and legacy Jacobian/band and step-control names. A separate `odeint_tfirst`
entry point provides the `func(t,y)` convention without a runtime procedure-type
switch. The default scalar tolerances are `sqrt(epsilon(1.0_dp))`; both `rtol` and
`atol` also accept componentwise vectors.

`test/test_integrate_odeint.f90` checks selected SciPy 1.17.0 outputs for the
documented damped-pendulum example, backward integration with an exactly repeated
requested time, dense/packed LSODA Jacobian routing, `col_deriv`, the t-first
callback convention, `h0`/`hmax`/`hmin`/`mxstep` controls, and invalid grids or
control values. Single-point and all-equal repeated grids were also checked
directly and return the initial state without spurious evolution. The packed
Jacobian regression exercises a stiff problem and requires nonzero `njev` and
`nlu`.

The initial facade used a pooled `solve_ivp` attempt budget for `mxstep`; that
limitation has since been removed. `mxstep` now counts accepted internal LSODA
steps separately for each requested output call, with zero selecting the DLSODA
legacy default of 500. Printing controls and full-output histories remain roadmap
items.
The implementation and tests use SciPy 1.17.0 only for public semantics and
numerical references; the facade itself is original SciFort code.

GNU Fortran 14.2.0 passes `test_integrate_odeint`, `test_integrate_ivp`,
`test_integrate_ivp_complex`, `test_integrate_sparse_lu`, and the four direct
LSODA rejection/stiff-switch/transfer/work-count programs in Fortran 2018 mode
with `-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow` and
`-Werror=implicit-interface`. The same targeted programs are also exercised in
an independent `-O2 -Werror=implicit-interface` build. The public umbrella was
compiled against narrow temporary stubs for unrelated quadrature modules because
this execution environment does not provide the FPM-resolved QUADPACK dependency.
`fpm` and `fprettify` are not installed, so repository-level `fpm test` and
automatic formatting could not be run.

## Componentwise relative-tolerance validation (2026-09-30)

`solve_ivp` now accepts `rtol` as either a scalar or a rank-one vector of length
`size(y0)`, matching the documented componentwise tolerance form. The vector is
used in initial-step selection and every adaptive error weight for RK23, RK45,
DOP853, Radau, BDF, LSODA, and the complex RK/BDF variants. Radau/BDF retain a
scalar Newton convergence threshold; for vector `rtol` SciFort conservatively
uses `maxval(rtol)` for that threshold while leaving the local-error scale fully
componentwise. LSODA carries the vector through both Adams and BDF branches and
through ODEPACK-style finite-difference Jacobian weighting.

`test/test_integrate_ivp.f90` integrates a two-rate decay problem with
`rtol=[1e-3,1e-8]` and verifies all six real methods meet the tight fast-component
target and perform more RHS work than the corresponding scalar `rtol=1e-3` run.
It also rejects short, rank-two, and nonpositive `rtol` inputs.
`test/test_integrate_ivp_complex.f90` applies the same regression to complex
RK23/RK45/DOP853/BDF, including invalid shape checks. The `odeint` facade now
accepts the same scalar-or-vector `rtol`; `test/test_integrate_odeint.f90`
verifies that tightening only the second component improves the fast-mode result
and increases LSODA work, and rejects a short vector.

GNU Fortran 14.2.0 passes `test_integrate_ivp`, `test_integrate_ivp_complex`,
`test_integrate_odeint`, `test_integrate_sparse_lu`, and the four direct LSODA
rejection/stiff-switch/transfer/work-count programs in Fortran 2018 mode with
`-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow` and
`-Werror=implicit-interface`. An independent `-O2 -Werror=implicit-interface`
build passes the same eight targeted programs. FPM and `fprettify` are unavailable
in the execution environment, so those repository-level commands were not run.


## `odeint` critical-time and order-control validation (2026-09-30)

The legacy `odeint` facade now accepts `tcrit`, `mxordn`, and `mxords`. The
`tcrit` sequence is validated as monotone in the integration direction and is
forwarded to the real LSODA driver as an active hard step bound. The bound is
updated without reinitializing LSODA, so Adams/BDF Nordsieck history and method
state survive a critical point. The behavior follows the DLSODA ITASK=4 idea:
the current critical time prevents internal work from overshooting that bound
while an earlier requested output is being obtained. The nonstiff and stiff
order ceilings are stored in LSODA state; the Adams controller is capped at
`mxordn` (maximum 12) and every BDF initialization/transfer/order increase is
capped at `mxords` (maximum 5). Zero selects the legacy default maximum.

`test/test_integrate_odeint.f90` instruments a decay RHS around `tcrit=1` and
requires an exact visit to the critical time with no RHS evaluation beyond it
before that visit. It also compares ordinary decay with the default Adams cap
against `mxordn=1`, and a stiff cosine-tracking problem with the default BDF cap
against `mxords=1`; both capped runs must take more accepted steps while meeting
their numerical tolerances. Negative order limits and a nonmonotone `tcrit`
sequence are rejected with `ivp_status_invalid_argument`.

GNU Fortran 14.2.0 passes `test_integrate_odeint`, `test_integrate_ivp`,
`test_integrate_ivp_complex`, `test_integrate_sparse_lu`, and all four direct
LSODA rejection/stiff-switch/transfer/work-count programs in Fortran 2018 mode
with `-fcheck=all -ffpe-trap=invalid,zero,overflow` and
`-Werror=implicit-interface`. An independent `-O2 -Werror=implicit-interface`
build passes the same eight targeted programs. Existing compiler diagnostics are
the project's pre-existing real/complex equality and unused-dummy warnings.
FPM and `fprettify` are unavailable in this execution environment, so those
repository-level commands were not run.

## `odeint` per-output `mxstep` validation (2026-09-30)

`odeint` now implements DLSODA's one-call `MXSTEP` semantics directly in the
continuous SciFort LSODA driver. A successful internal step increments a call-local
counter. Whenever the current internal mesh point has reached or passed the next
requested `t` value, the counter resets exactly as a subsequent DLSODA call would,
without reinitializing the Adams/BDF state. Repeated requested times therefore cost
zero additional steps, and an internal step that overshoots more than one requested
output can satisfy those zero-step follow-up calls by dense interpolation. Rejected
trial steps do not consume `mxstep`, consistent with DLSODA's `NST-NSLAST` check.
`mxstep=0` selects 500.

`test/test_integrate_odeint.f90` uses `t=[0,0.1,1]`, `hmax=0.15`,
`rtol=1e-4`, and `atol=1e-8`. Local SciPy 1.17.0 reports excess work for
`mxstep=6` after cumulative `nst=10`, `nfe=21`, and succeeds for `mxstep=7` with
`nst=11`, `nfe=23`; SciFort now matches those counts exactly. The same regression
also verifies that `mxstep=0` and explicit `mxstep=500` produce identical work on a
simple decay problem.

GNU Fortran 14.2.0 passes `test_integrate_odeint`, `test_integrate_ivp`,
`test_integrate_ivp_complex`, `test_integrate_sparse_lu`, and all four direct LSODA
rejection/stiff-switch/transfer/work-count programs in Fortran 2018 mode with
`-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow` and
`-Werror=implicit-interface`. An independent `-O2 -Werror=implicit-interface`
build passes the same eight targeted programs. FPM and `fprettify` are unavailable
in the execution environment, so those repository-level commands were not run.

## `odeint` full-output diagnostic validation (2026-09-30)

The legacy `odeint` facade now accepts `full_output=.true.` and records LSODA
diagnostics at each requested output call after the initial time. The result
contains `hu`, `tcur`, `tolsf`, `tsw`, `nst`, `nfe`, `nje`, `nqu`, and `mused`
vectors plus scalar `imxer`, `lenrw`, and `leniw`. LSODA state retains the last
successful step size and order, the method actually used on that step, and the
latest Adams/BDF switch time; this avoids confusing a just-selected next method
with the method responsible for the returned solution. Cumulative work counters
are captured when the internal mesh first reaches each requested output.

For `y'=-y`, `y(0)=1`, and requested times `[0,0.1,0.2]`, a local SciPy 1.17.0
reference run reports `nst=[10,13]`, `nfe=[23,29]`, `nje=[0,0]`, `nqu=[4,4]`,
and `mused=[1,1]`; SciFort matches those integer histories exactly. Its forward
`hu` values are positive, while a backward repeated-time regression verifies
negative `hu` and deterministic initial diagnostics for an exactly repeated
initial output. A stiff cosine-tracking case verifies BDF `mused=2`, BDF orders
in the range 1--5, monotone cumulative work, and persistence of the recorded
method-switch time. For the one-equation dense case with default order ceilings,
`lenrw=36` and `leniw=21`, matching the DLSODA workspace formulas.

SciFort intentionally returns `tolsf=0` on normal successful integrations and
`imxer=-1`; component-resolved diagnostics for DLSODA error returns are not yet
modeled. Also, because SciFort's `solve_ivp` enforces the final `t_span` endpoint
as a hard bound, its final `tcur` can equal the requested endpoint instead of
overshooting it and interpolating back. When `full_output` is omitted or false,
the diagnostic vectors remain zero length and `lenrw=leniw=0`.

GNU Fortran 14.2.0 passes `test_integrate_odeint`, `test_integrate_ivp`,
`test_integrate_ivp_complex`, `test_integrate_sparse_lu`, and all four direct
LSODA rejection/stiff-switch/transfer/work-count programs in Fortran 2018 mode
with `-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow` and
`-Werror=implicit-interface`. An independent `-O2 -Werror=implicit-interface`
build passes the same eight targeted programs. Existing compiler diagnostics are
the project's pre-existing real/complex equality and unused-dummy warnings. FPM
and `fprettify` are unavailable in this execution environment, so those
repository-level commands were not run.



## `odeint` legacy reporting-control validation (2026-09-30)

The legacy facade now accepts the remaining documented reporting controls
`ixpr`, `mxhnil`, and `printmessg`. `ixpr=.true.` compares the LSODA method state
before and after every successful internal step and prints a message only when
an Adams/BDF switch actually occurs. `mxhnil=0` selects the DLSODA default
message limit of 10; a positive value supplies an explicit cap. SciFort retains
its existing policy that a step below its LSODA machine-progress floor is a
terminal failure, so unlike DLSODA's continuing `T+H=T` warning path, at most one
hnil warning can currently be observed in one integration. `printmessg=.true.`
prints the final status, using `Integration successful.` on successful runs.
These controls do not alter tolerances, step selection, work accounting, or the
returned result.

`test/test_integrate_odeint.f90` runs the stiff cosine-tracking problem once with
all reporting disabled and once with `ixpr=.true.`, `mxhnil=3`, and
`printmessg=.true.`. The verbose run emits one Adams-to-BDF switch line and the
final success message while matching the quiet run exactly in final state,
accepted-step count, RHS evaluations, and Jacobian evaluations. A negative
`mxhnil` is rejected with `ivp_status_invalid_argument`. A separate local
huge-time probe exercises the too-small-step reporting path and emits the hnil
warning before returning the pre-existing LSODA failure status.


## Rank-2 and weighted special reductions validation (2026-09-30)

This batch extends the original real one-dimensional `logsumexp`, `softmax`,
and `log_softmax` interfaces to real rank-2 arrays. `logsumexp` supports whole-
array reduction or a single Fortran `dim`, plus optional same-shape weights.
The new `logsumexp_signed` subroutine provides the log-absolute-value and sign
outputs corresponding to SciPy's `return_sign=True`. `softmax` and
`log_softmax` accept rank-2 arrays with an optional `dim`; omitting it
normalizes over the complete array. General NumPy-style broadcasting, tuples of
axes, and `keepdims` remain future work; complex reductions were added in the
subsequent checkpoint described below.

Focused references use the documented SciPy matrix example and local SciPy
1.17.0 weighted/sign outputs. Tests cover both dimensions, whole-array
normalization, negative weighted totals, exact zero totals, opposite-sign
infinite contributions, shape/dimension validation, and the rule that zero
weights are ignored before exponentiation. A deterministic 85-case audit used
logarithmic inputs spanning roughly `[-1000,1000]` and signed log-uniform
weights; every ordinary/infinite result agreed with the local SciPy reference.
One local SciPy 1.17.0 overflow corner case (`b=0` with a large finite positive
logarithm) produces NaN there; SciFort intentionally returns the mathematically
stable zero weighted sum because the disabled term is never exponentiated.

GNU Fortran 14.2.0 compiled the 67-source `scifort_special` dependency closure
from scratch with `-std=f2018 -O0 -g -fcheck=all -Wall -Wextra -fbacktrace`.
The focused `test_special_elementary` program also passed with
`-ffpe-trap=invalid,zero,overflow`, and all 71 `test_special*.f90` programs
passed in the checked configuration (the broad sweep traps zero/overflow but
allows intentional IEEE invalid cases). An independent `-O2` rebuild of the
same 67-source closure passed all 71 special tests again. A separate checked
compile of the complete 183-source `scifort_stats` dependency closure also
succeeded, verifying the umbrella dependency graph after the public reduction
export change. FPM and `fprettify` are unavailable in this execution
environment.


## Complex special reductions validation (2026-09-30)

This checkpoint extends `logsumexp`, `logsumexp_signed`, `softmax`, and
`log_softmax` to complex vectors and rank-2 arrays. Complex `logsumexp` accepts
same-shape real or complex weights. The signed form returns a real
log-magnitude and the complex unit phase, matching SciPy's documented
`return_sign=True` convention for complex inputs. The implementation shifts by
the largest logarithmic term magnitude and uses compensated real/imaginary
summation, so finite inputs with very large real offsets do not require forming
unscaled exponentials. Zero weights are skipped before exponentiation.

`test/test_special_reductions_complex.f90` checks vector, whole-matrix, and
single-dimension reductions; real and complex weights; complex phases;
softmax/log-softmax normalization; exact cancellation; zero-weight masking;
and invalid weight/dimension handling. Its ordinary reference values were
cross-checked against SciPy 1.17.0. An additional deterministic 80-case audit
used complex inputs with real offsets spanning approximately `[-700,700]` and
complex signed weights; SciFort agreed within `2e-13` with mpmath 90-digit
reference values. In a separate comparison, SciFort was more accurate than the
local SciPy result in a cancellation-sensitive case, so the permanent tests use
mathematical/high-precision references rather than forcing bit-level SciPy
agreement.

GNU Fortran 14.2.0 compiled the 67-source `scifort_special` dependency closure
from scratch with bounds/runtime checks and floating-point invalid/zero/overflow
traps. All 72 `test_special*.f90` programs passed. An independent `-O2` rebuild
of the same closure passed all 72 programs again. FPM and `fprettify` are not
available in this execution environment.


## Broadcast and keepdims special reductions validation (2026-10-01)

This checkpoint extends the real/complex vector and rank-2 reduction API with
fixed-rank NumPy-style weight broadcasting and multiple-axis behavior. Vector
weights may be scalar, length one, or full length. Rank-2 weights may be scalar,
a trailing vector, or a matrix with singleton/matching extents. Subsequent
checkpoints generalized this to joint broadcasting that can expand singleton
input dimensions. Complex weights are accepted for real as well as complex inputs.

`softmax` and `log_softmax` now accept an integer `dims(:)` list for rank-2
arrays. Empty lists implement NumPy/SciPy `axis=()` behavior, one-element lists
select a Fortran dimension, and `[1,2]` / `[2,1]` normalize jointly over both.
The explicit `logsumexp_keepdims` and `logsumexp_signed_keepdims` companions
retain reduced dimensions with extent one; separate names are used because
Fortran function result rank cannot depend on a runtime logical flag.

`test/test_special_reductions_broadcast.f90` checks scalar/vector/row/column
weight broadcasting, real-input complex weights, signed/phase outputs, empty and
two-axis reductions, keepdims shapes, duplicate-axis rejection, and real/complex
softmax normalization. Reference values were independently cross-checked against
local SciPy 1.17.0 public outputs. The implementation is an original SciFort
translation of the documented reduction/broadcast semantics; no SciPy source was
copied or adapted.

GNU Fortran 14.2.0 compiled the 67-source `scifort_special` dependency closure
from scratch with `-std=f2018 -O0 -g -fcheck=all -fbacktrace`; the full 73-program
`test_special*.f90` sweep passed with invalid/zero/overflow traps enabled. An
independent `-O2` rebuild passed all 73 programs again. A separate checked build
of the complete 183-source `scifort_stats` dependency closure also succeeded.
FPM and `fprettify` are unavailable in this execution environment.

## Rank-4 weighted special reductions validation (2026-10-01)

This checkpoint extends weighted `logsumexp` parity to rank-4 arrays. Real and
complex rank-4 inputs accept scalar through rank-4 weights using NumPy-style
trailing-dimension broadcasting. Broadcasting is joint: singleton dimensions of
`x` may expand to the corresponding weight extent. Real inputs may use real or
complex weights; complex inputs may use real or complex weights. The same
broadcasting applies to single-dimension reductions, `logsumexp_signed`,
`logsumexp_keepdims`, and `logsumexp_signed_keepdims`.

`test/test_special_reductions_rank4_weighted.f90` checks weight ranks zero through
four, expansion of `(2,1,3,1)` inputs to `(2,4,3,5)`, single-dimension reductions,
all 16 rank-4 keepdims axis subsets for real and complex data, negative real weights
and signs, real inputs with complex weights, complex inputs with real weights,
complex inputs with complex weights, empty broadcast extents, and
incompatible-shape rejection. Numerical references were generated from the
public SciPy 1.17.0 `logsumexp` interface; no SciPy implementation source was
copied or adapted.

GNU Fortran 14.2.0 compiled `scifort_kinds`, `scifort_math`, and
`scifort_special_reductions` with `-std=f2018 -O0 -Wall -Wextra
-Wimplicit-interface -fcheck=all`. The new weighted rank-4 regression and the six
pre-existing reduction regression programs all pass against that clean module
build. An independent `-O2` rebuild passes the same seven programs. FPM and
`fprettify` are unavailable in this execution environment.

## Rank-5 unweighted special reductions validation (2026-10-01)

This checkpoint extends unweighted real/complex `logsumexp`,
`logsumexp_signed`, `logsumexp_keepdims`, `logsumexp_signed_keepdims`, `softmax`,
and `log_softmax` to rank-5 arrays. A dimension bit mask validates tuple-style
Fortran dimension lists, and each selected rank-5 slice is reduced through the
existing stable one-dimensional kernels rather than duplicating 32 handwritten
axis combinations. Weighted broadcasting remains supported through rank 4.

`test/test_special_reductions_rank5.f90` checks whole-array and single-dimension
reductions, nonadjacent/reversed tuple dimensions, four-axis and empty-axis
keepdims behavior, signed empty-axis `-infinity`, real and complex softmax and
log-softmax normalization, duplicate-axis rejection, and complex phase results.
Reference values were generated from the public SciPy 1.17.0 interfaces. During
that comparison, a complex `log_softmax` case exposed a branch difference when
computing `x - logsumexp(x)` directly; SciFort now follows SciPy's max-shift
branch convention, with a dedicated rank-1 regression for the `2*pi*i` case.

GNU Fortran 14.2.0 compiled `scifort_kinds`, `scifort_math`, and
`scifort_special_reductions` with `-std=f2018 -O0 -g -fcheck=all -Wall -Wextra`.
The new rank-5 regression and all seven pre-existing reduction regression programs
pass against the clean module build using a local umbrella-module shim for tests
that normally import `scifort_special`. An independent `-std=f2018 -O2` rebuild
passes the same eight programs, including zero-extent rank-5 reductions. FPM and
`fprettify` are not installed in this environment.


## Rank-5 weighted special reductions validation (2026-10-01)

This checkpoint extends weighted `logsumexp` parity to rank-5 real and complex
arrays. Scalar through rank-5 weights use NumPy-style trailing-dimension joint
broadcasting, so singleton dimensions of either `x` or `b` may expand. Real
inputs accept real or complex weights, complex inputs accept real or complex
weights, and the same semantics apply to single-dimension, signed/phase, and
fixed-rank keepdims forms. Zero extents preserve the empty weighted-sum result.

`test/test_special_reductions_rank5_weighted.f90` checks weight ranks zero through
five, joint expansion from `(2,1,2,1,3)` to larger broadcast shapes, negative
weighted totals and sign output, single-dimension reduction, all 32 rank-5
keepdims axis subsets for real and complex data, real inputs with complex weights,
complex inputs with real weights, complex inputs with complex weights, unsupported
weight-rank rejection, and zero-extent broadcasting. Numerical references were
generated from the public SciPy 1.17.0 `logsumexp` interface; no SciPy
implementation source was copied or adapted for the rank-5 weighted code.

GNU Fortran 14.2.0 compiled `scifort_kinds`, `scifort_math`, and
`scifort_special_reductions` with `-std=f2018 -O0 -g -fcheck=all -Wall -Wextra
-Wimplicit-interface`. The new weighted rank-5 regression and all eight
pre-existing reduction regression programs pass against that checked module build,
using a local umbrella-module shim only for tests that normally import
`scifort_special`. A separate `-std=f2018 -O2 -fno-inline-functions` rebuild
passes the same nine programs. FPM and `fprettify` are not installed in this
environment.

## Rank-15 weighted special reductions validation (2026-10-01)

This checkpoint extends weighted `logsumexp` parity from rank 5 through rank 15,
which is the standard Fortran maximum array rank. Real inputs accept real or
complex weights, complex inputs accept real or complex weights, and scalar through
rank-15 weights use NumPy-style trailing-dimension joint broadcasting. Singleton
dimensions of either `x` or `b` may expand; zero extents and incompatible shapes
retain the established weighted reduction behavior. Signed/phase-returning,
single-dimension, and fixed-rank keepdims forms use the same broadcast engine.

The high-rank weighted implementation is split into a shared broadcast/index
support module plus rank-6--10 and rank-11--15 wrapper modules behind the unchanged
`scifort_special_reductions` facade. `test/test_special_reductions_rank15_weighted.f90`
checks SciPy 1.17.0 reference values at ranks 6 and 15, full weight ranks 6 through
15, singleton expansion at rank 15, negative weighted totals, all real/complex
input-weight combinations, a weight whose rank exceeds the input rank, incompatible
broadcasting, and zero-extent reductions. No SciPy implementation source was
copied or adapted for this high-rank weighted extension.

GNU Fortran 14.2.0 passes all 12 reduction regression programs with a checked
build using `-std=f2018 -O0 -g -fcheck=all -Wall -Wextra -Wimplicit-interface`.
An independent optimized build using `-std=f2018 -O2 -fno-inline-functions
-Wall -Wextra -Wimplicit-interface` passes the same 12 programs. Tests that
normally import the full `scifort_special` umbrella use only a local test shim to
import `scifort_kinds` and `scifort_special_reductions`; repository sources are
unchanged by that shim. FPM and `fprettify` are not installed in this environment.

## Scalar root finding

`test_optimize_root_scalar` compares the initial `scifort_optimize` scalar-root
surface against values generated with SciPy 1.17.0. For `x^3-1` and
`cos(x)-x`, `bisect`, `brentq`, `brenth`, `ridder`, and default-`k` `toms748`
are checked against SciPy root values and exact iteration/function-call counts.
Secant, derivative Newton, Halley, and the `root_scalar` numerical-derivative
Newton path receive the same accounting checks. Additional regressions cover
a one-iteration Brent nonconvergence result, non-bracketing endpoints, nonfinite
objective values, automatic method selection, and `toms748(k=2)` convergence.
The test is compiled with Fortran 2018, all warnings as errors, runtime checks,
and IEEE invalid/zero/overflow traps.

## Scalar minimization (2026-10-04)

`test_optimize_minimize_scalar` compares `minimize_scalar`, `brent`, `bounded`,
and `golden` with references generated from SciPy 1.17.0 for the polynomial
`(x-2)*x*(x+2)^2`. Two- and three-point Brent/golden brackets, default and
tight bounded tolerances, exact iteration/function-call counts, forced
nonconvergence, invalid bounds/methods, automatic method selection, and NaN
objectives are covered. GNU Fortran 14.2.0 passes the test with
`-std=f2018 -Wall -Wextra -Werror=implicit-interface -fcheck=all -fbacktrace`.
The scalar-minimization source and regression also compile cleanly with
`-Werror` in addition to the flags above.

## Multivariate minimization trust-krylov (2026-10-05)

`test_optimize_minimize` now includes SciPy 1.17.0 references for
`minimize(method='trust-krylov')` and direct `trust_krylov`. Dense-Hessian and
Hessian-product Rosenbrock paths, `inexact=.false.`, custom trust radii and
acceptance threshold, a one-iteration limit, a convex quadratic, and a
five-variable Rosenbrock Hessian-product case are covered. The native Fortran
solver retains its Krylov basis when a trial step is rejected so a repeated
subproblem at the smaller trust radius hot-starts like SciPy's TRLIB wrapper.
The deterministic cases match SciPy's final point/value, status, outer iteration
count, function/gradient evaluations, and Hessian/Hessian-product work counts.
The modified minimization module and regression compile with Fortran 2018 and
`-Werror=implicit-interface`; a checked manual build passes the complete
multivariate-minimization regression. Full `fpm test` remains unavailable in
this environment because FPM is not installed.

The `odeint` callback fix was independently validated by compiling the complete
local dependency chain through `ivp`, LSODA, BDF, Radau, and DOP853 with the
same implicit-interface error policy. The existing `test_integrate_odeint`
regression passes against those objects, including dense/banded Jacobian and
`col_deriv` cases. FPM itself is not installed in this environment, so the full
`fpm test` command could not be run here.

## Adaptive coordinate-scaled Krylov Jacobian differences (2026-10-10)

`test/test_optimize_root_krylov_adaptive_steps.f90` checks the new opt-in
`jacobian_step_mode='adaptive'` and `'componentwise'` interfaces in both public
Krylov entry points. With a variable of order `1e13` beside one of order
`0.1`, the default global Jacobian step can leave the residual domain, while
coordinate-scaled forward and centered probes recover the analytic
Jacobian-vector product. A second case uses residual offsets of order `1e12`
to expose cancellation from the global residual-normalized step. The test
checks sparse/zero directions, exact function-evaluation increments,
componentwise minimum representable steps, all six built-in Krylov methods,
and invalid option combinations. The existing centered-stencil and default
SciPy 1.17.0 work-count tests are retained unchanged.

A 12-program Krylov/nonlinear-root regression suite passed in both `-O0 -g`
and `-O2` builds with GNU Fortran 14.2, `-fcheck=all`,
`-ffpe-trap=invalid,zero,overflow`, and `-Werror=implicit-interface`.
These are manually linked, focused tests using a temporary MINPACK interface
shim; they are not a replacement for the dependency-backed full `fpm test`.

## Centered Jacobian-vector differences in Krylov root methods (2026-10-10)

`test/test_optimize_root_krylov_central_difference.f90` checks SciFort's new
opt-in `jacobian_diff_scheme='3-point'` against exact directional derivatives
of a smooth two-variable nonlinear residual. It checks that a custom inner
solver receives the same stencil as all six built-in methods
(GMRES, LGMRES, BiCGSTAB, CGS, MINRES, TFQMR), and that `newton_krylov` forwards the option.
The automatic centered-step scale uses cube-root machine epsilon, while the
unchanged default two-point scale uses square-root epsilon; explicit `rdiff`
remains supported. The test also checks exact `nfev` increments (one residual
call for forward products, two for centered products), the zero direction,
invalid scheme/method combinations, and rejection of nonfinite or
unrepresentable samples. A finite large-magnitude signed-sample example tests
that halving before subtracting avoids unnecessary intermediate overflow.

On GNU Fortran 14.2.0, the focused 11-program suite passes with `-std=f2018`
in both `-O0 -g` and `-O2` builds using `-fcheck=all`,
`-ffpe-trap=invalid,zero,overflow`, and `-Werror=implicit-interface`. The
manual public-root suite links a temporary MINPACK test shim; it does **not**
substitute for the project's full dependency-backed `fpm test`.

## Multidimensional root finding (2026-10-04)

`test_optimize_root_nonlin` compares the translated SciPy 1.17 inexact-Newton,
spectral-residual, and matrix-free Krylov root methods against generated SciPy
1.17.0 references. The suite covers `broyden1`, `broyden2`, `anderson`,
`diagbroyden`, `linearmixing`, `excitingmixing`, `df-sane`, and `krylov`, plus
the direct `broyden1` and `newton_krylov` entry points. It checks exact
iteration/function-evaluation counts for the default Armijo path, no line search,
strong-Wolfe search, generic tolerances, iteration/evaluation limits, DF-SANE
Cruz/Cheng searches, LGMRES and GMRES, disabled LGMRES recycling, custom
finite-difference `rdiff`, stricter residual tolerance, and Broyden restart,
simple, and SVD rank reduction. A deliberately difficult Broyden-2/simple-rank
case reproduces SciPy's nonconvergence status after 300 iterations and 1343
residual evaluations.

GNU Fortran 14.2.0 passes `test_optimize_root_nonlin`,
`test_optimize_root_scalar`, and `test_optimize_minimize_scalar` in a checked
manual build using `-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface
-fcheck=all`. The new `wolfe.f90`, `nonlin.f90`, `root.f90`, generated reference,
and regression program all compile under that warning-as-error policy. The
existing `src/linalg.f90` has a pre-existing real-equality warning when compiled
with blanket `-Werror`; it is compiled without blanket `-Werror` in this manual
harness. `test_optimize_root.f90`, which exercises the external MINPACK-backed
`hybr`/`lm` paths, compiles against an interface-accurate local MINPACK stub so
that callback/interface conformance is checked, but the stub is not used for a
numerical validation claim and is not packaged. FPM/dependency fetching is not
available in this environment, so the real MINPACK numerical test still needs a
normal local `fpm test` run.

The default LGMRES and GMRES Krylov reference cases reproduce SciPy's converged
iterates and work counts. A tiny two-variable nonlinear problem with SciPy 1.17's
default recycled LGMRES is pathological: both implementations reach the 300-step
limit, but the detailed recycled trajectory/evaluation count differs. Permanent
Krylov regressions therefore exercise recycling on a convergent linear system and
use `outer_k=0` for the nonlinear work-count cases rather than encoding that
pathological trajectory as a compatibility requirement.

## Custom inexact-Newton/Krylov termination norms (2026-10-10)

`test/test_optimize_root_tol_norm.f90` exercises `tol_norm=` for all six
inexact-Newton methods and matrix-free Krylov, as well as every direct public
wrapper. Fortran stopping points and iteration/evaluation counts are checked
against SciPy 1.17.0 affine-system reference runs with anisotropic residual and
step norms. The test also verifies default max-norm parity, zero-residual
initial states, invalid/NaN norm handling and unsupported-method diagnostics.

Local GNU Fortran 14.2 testing compiles `src/optimize/nonlin.f90` and
`src/optimize/root.f90` with Fortran 2018, bounds checks, and explicit
interfaces. The two numerical regressions, including existing
`test_optimize_root_nonlin`, pass. Because FPM and the external MINPACK
Fortran module are not installed here, the root-module tests link with a
strictly test-only MINPACK interface stub; no MINPACK numerical routine is
exercised or packaged. Full `fpm test` remains necessary on a provisioned
build system.

## Independent nonlinear constraint blocks for trust-constr (2026-10-07)

The new `nonlinear_constraint`/`constraint_objects=` adapter is covered by
`test/test_optimize_trust_constr_objects.f90`. Its reference file is generated
from SciPy 1.17.0 by
`tools/generate_optimize_trust_constr_objects_reference.py`. The separate
analytic equality/lower-inequality case agrees with SciPy in its final solution,
objective, 11 iterations, seven objective calls, and per-object
constraint work counts (7/7 function, 7/7 Jacobian, 12/12 Hessian).
With the second Jacobian numerical, the analogous SciPy counts are 7/21
function and 7/0 analytic Jacobian calls; the adapter reproduces these.
Additional regressions check generic dispatch, multiple constraint Hessians,
per-object multiplier routing, variable-bounds reporting, finite-difference
clipping at a retained bound, and invalid argument combinations.

GNU Fortran 14.2.0 compiled the affected module and both
`test_optimize_trust_constr` and `test_optimize_trust_constr_objects` with
`-std=f2018 -Werror=implicit-interface`. Both tests passed in a debug build
with `-O0 -g -fcheck=all`, and again in an optimized `-O2` build. These are
selected numerical regression tests, not the entire FPM suite. FPM and the
original external L-BFGS-B/SLSQP dependencies are unavailable in this
sandbox. The manual build used temporary, interface-compatible stubs for
these *unexecuted* algorithms; the stubs are not shipped. A normal `fpm test`
on a machine that can resolve the declared dependencies remains required.

## Unconstrained trust-constr (2026-10-07)

`test/test_optimize_trust_constr_unconstrained.f90` exercises direct
`trust_constr` and generic `minimize(method='trust-constr')` with no user
constraints and no variable bounds. References are generated against SciPy
1.17.0 by `tools/generate_optimize_trust_constr_unconstrained_reference.py`.
The selected analytic-Hessian quadratic, numerical-gradient quadratic,
initially stationary, one-iteration limit, and analytic Rosenbrock cases match
SciPy's termination status, solution, objective, and exact objective function,
gradient, Hessian, and CG iteration counts. The quadratic BFGS-Hessian case
converges to the same minimizer, but detailed quasi-Newton iteration counts
are not claimed to match. With `-O3 -march=native`, cancellation in the
finite-difference gradient close to the quadratic optimum can cause extra
rejected trial steps and successful `xtol` rather than `gtol` termination;
the regression permits this floating-point sensitivity while still checking
the final solution, derivative counts, and empty constraint diagnostics. A
zero-length `constraint_objects=` array is also
covered, both with and without finite variable bounds.

An unconstrained result uses an empty multiplier vector, a zero-row constraint
Jacobian, zero-length per-constraint evaluation arrays, and zero constraint
violation. The existing constrained tests cover the nonzero-row path. Checked
and optimized manual builds validate these focused tests using local interface
stubs for the unexecuted L-BFGS-B and SLSQP external dependencies. Those
stubs are not included in the package. The full `fpm test` suite still requires
FPM and the declared external packages.

## Three-point objective gradients in trust-constr (2026-10-07)

`test/test_optimize_trust_constr_objective_fd.f90` covers an explicit
`objective_fd_scheme='3-point'` option for direct `trust_constr`, generic
`minimize(method='trust-constr')`, and independent constraint blocks. Static
reference data are regenerated using installed SciPy 1.17.0 via
`tools/generate_optimize_trust_constr_objective_fd_reference.py`. The solver
matches SciPy's exact objective/gradient/Hessian evaluation counts and
iterations in unconstrained (3/15/3/3), linear-equality (3/15/3/3), and
nonlinear-equality with numeric constraint Jacobian (12/49/9/9) cases.
The nonlinear-constraint block separately retains zero explicit Jacobian
callback counts, despite its three-point constraint-function probes.

A retained lower bound is tested from an interior starting point, and a
separate one-iteration setup test starts exactly at both lower and upper
retained endpoints. The endpoint cases check that no objective callback
crosses the bound, the one-sided objective derivative remains accurate,
and the analytic-Hessian barrier initialization is finite. A finite initial
slack floor of `sqrt(epsilon(1.0_dp))` is used only for inequalities marked
`keep_feasible`, avoiding overflow from a tiny positive slack in the barrier
Hessian. The bounded optimizer's full iteration/evaluation counts can differ
from SciPy because the dense local barrier controller remains approximate.

The tests also cover an explicit `eps=` step, preservation of analytic
objective-gradient callbacks, a bogus scheme, and rejection of this
trust-constr-only option by other `minimize` methods. The four prior
`trust-constr` regression executables should continue to pass. Manual
compilation uses local L-BFGS-B/SLSQP interface stubs only for the unexecuted
external methods; the ZIP does not ship these stubs. The full FPM test suite
requires FPM and its external dependencies.


## Objective relative finite-difference steps in trust-constr (2026-10-08)

`test/test_optimize_trust_constr_objective_rel_step.f90` verifies the new
optional `objective_finite_diff_rel_step(:)` argument with both two- and
three-point numerical objective gradients. The independently generated
`test/optimize_trust_constr_objective_rel_step_reference.f90` stores SciPy
1.17.0 `_numdiff.approx_derivative` values for `[x1**3 + 2*x2**3]` at
`x=[2,-3]`, with unequal relative steps `[0.1,0.2]`. The reference generator
is `tools/generate_optimize_trust_constr_objective_rel_step_reference.py`.
Tests verify the signed two-point and centered three-point derivatives, a
second-order one-sided derivative at a retained lower bound, an exactly zero
coordinate, strict callback feasibility, objective evaluation counts, and
forwarding through generic `minimize` and `constraint_objects`. Additional
cases check exact-Jacobian precedence, invalid array shapes, nonfinite step
values, and rejection of this trust-constr-only keyword by BFGS. A full
unconstrained quadratic solve checks that the new vector persists through
later SQP iterations instead of affecting just initialization.

The public vector contains one entry per **original** variable, not per
barrier slack. The internal slack coordinates continue using the normal
finite-difference step, and bounded one-sided probes never cross retained
variable bounds. The absence of this keyword preserves earlier behavior.
The standalone targeted build uses local interface stubs for unexecuted
L-BFGS-B and SLSQP dependencies; these are not distributed with SciFort.
Full FPM integration tests still require the declared external libraries.

## Multiplier-weighted constraint Hessian products in trust-constr (2026-10-08)

`test/test_optimize_trust_constr_constraint_hessp.f90` verifies an optional
`nonlinear_constraint%hessp(x, v, p, hp)` callback, with `hp` equal to the
Hessian of `dot_product(v, c(x))` applied to direction `p`. SciPy 1.17.0
references are generated by
`tools/generate_optimize_trust_constr_constraint_hessp_reference.py`, using
`scipy.sparse.linalg.LinearOperator` from `NonlinearConstraint.hess(x, v)`.
The new tests exercise a nonconstant, off-diagonal Hessian and compare it to
an equivalent full `%hess` callback, full-Hessian precedence, generic
`minimize` dispatch, bound-interval inequalities (barrier path), simultaneous
objective and constraint products, mixed nonlinear/linear blocks, mixed exact
full/product constraint Hessians across different nonlinear blocks, two-output
multiplier vectors, and rejection of callbacks on constant-matrix blocks.

The dense SQP implementation assembles each nonlinear block Hessian using
`size(x0)` products and symmetrizes it before accumulation. `constr_nhev`
counts these assemblies rather than products. Full `%hess` takes precedence
over `%hessp` on the same block. Mixed exact and missing constraint Hessians
are supported by numerical differentiation of each missing block's Jacobian;
the legacy combined-constraint callback is unchanged.

The ten focused trust-constr programs
can be compiled against `scifort_optimize_minimize` independently of the
external FPM libraries. Such local builds use stubs for the unrelated
L-BFGS-B/SLSQP entry points; they do not replace the full `fpm test` suite.

## BFGS objective with exact nonlinear-constraint Hessians (2026-10-08)

`test/test_optimize_trust_constr_bfgs_exact_constraint.f90` exercises
the previously rejected combination of a missing objective `hess=` or
`hessp=` and an exact `constraint_hess=` callback or per-block `%hess`/
`%hessp`.
The equality-SQP kernel now updates objective curvature with BFGS using
the **objective gradient difference**, then adds the exact weighted
constraint Hessian at the current Lagrange multipliers. It must not update
BFGS with the full Lagrangian gradient difference in this case, because
that would double-count constraint curvature. This also covers a mixed
set of exact and numerically differentiated constraint Hessians.

The interior-point wrapper passes exact constraint Hessians to its barrier
equality-SQP subproblems even without an exact objective Hessian. The tests
cover direct equality and upper-inequality constraints, bound-constrained
nonlinear constraint objects, a numerically differentiated three-point
objective gradient, and two nonlinear equality blocks with different
Hessian sources, including a Hessian-vector product callback. The exact
callback's invocation count is checked against
`constr_nhev`; `nhev` remains zero for a BFGS objective. Reference
minimizers are generated by
`tools/generate_optimize_trust_constr_bfgs_exact_constraint_reference.py`
against SciPy 1.17.0, with no exact objective Hessian supplied to SciPy.

The ten targeted Fortran `trust-constr` test programs passed with GNU
Fortran 14.2 in the isolated build, where only test-only placeholders for
the unrelated external L-BFGS-B/SLSQP procedures were used. These stubs
are not part of the distribution. The complete FPM dependency build and
test suite were not available for this checkpoint.

## Root Krylov preconditioner and custom-inner-solver regression (2026-10-10)

New `test/test_optimize_root_krylov_hooks.f90` exercises matrix-free
right-preconditioned GMRES and LGMRES, adaptive inverse-Jacobian updates,
custom inner linear solves, the default identity action for a custom solver,
and rejection of invalid or nonfinite callbacks. Reproducible SciPy 1.17.0
reference values come from
`tools/generate_optimize_root_krylov_hooks_reference.py`:

| SciPy 1.17.0 case | `nit` | `nfev` | Preconditioner updates |
| --- | ---: | ---: | ---: |
| 4-dimensional diagonal GMRES, one inner iteration | 2 | 4 | -- |
| 4-dimensional diagonal LGMRES, one inner iteration | 2 | 3 | -- |
| 2-dimensional nonlinear adaptive LGMRES | 7 | 16 | 6 |

On GNU Fortran 14.2, the tests passed with `-std=f2018 -fimplicit-none
-Werror=implicit-interface -fcheck=all -ffpe-trap=invalid,zero,overflow`.
The focused root module tests used a test-only MINPACK interface shim because
FPM and the external MINPACK dependency were unavailable in the test
container; the shim is not part of the SciFort archive. The standalone
nonlinear-root and custom-norm regressions were rerun to check that default
solver behavior and its evaluation counts remain unchanged. A full `fpm test`
was not performed in this environment.


## Additional Krylov inner methods (2026-10-10)

`test/test_optimize_krylov_iterative.f90` validates new matrix-free
BiCGSTAB, CGS, MINRES, and TFQMR recurrences with finite-difference-compatible
operator callbacks, exact diagonal preconditioning, a symmetric-indefinite
MINRES system, incomplete iteration budgets, zero right-hand sides, and
failed/nonfinite operator outputs. It also checks rejection of an indefinite
MINRES preconditioner. The inner methods return a finite partial iterate if
only the iteration budget is exhausted, as `scipy.optimize.KrylovJacobian`
does; actual breakdown or invalid callback output is rejected.

`test/test_optimize_krylov_iterative_stress.f90` checks deterministic
16-variable nonsymmetric solves, positive preconditioning of a
symmetric-indefinite MINRES problem, and failures *after* successful earlier
matrix-vector evaluations. It uses identity preconditioning for TFQMR on the
larger nonsymmetric system: SciPy 1.17.0 and SciFort's equivalent
left-preconditioned TFQMR recurrence can terminate early based on the
preconditioned residual estimate with the chosen approximate inverse.
In a separate SciPy check, that case produced an original relative residual
near 1.72e-3 despite reporting successful inner convergence. The default
unpreconditioned TFQMR met the target residual tolerance.

`test/test_optimize_krylov_true_residual.f90` checks the new, opt-in
original-residual criterion. The default TFQMR still returns a finite iterate
for the deterministic 16-variable diagonally preconditioned case, but its true
relative residual is approximately `1.7183e-3` despite a requested `1e-10`
inner tolerance. The verified mode instead reports failure; with identity
preconditioning it succeeds and satisfies the original relative tolerance.
Further cases cover rejection of an inaccurate truncated iterate, failure
on the verification matrix-vector product, exact zero right-hand sides,
unsupported solver names, invalid output dimensions, and nonfinite inputs.
`tools/generate_optimize_krylov_true_residual_reference.py` reproduces the
SciPy 1.17.0 default behavior (`1.71831244e-3` for diagonal versus
`4.80472218e-12` for identity preconditioning). The guard is a SciFort
extension and is not claimed to reproduce SciPy's stopping behavior when
explicitly enabled. `test/test_optimize_krylov_large_system.f90` tests all
four inner solvers on 64-variable deterministic nonsymmetric and
symmetric-indefinite systems, with nonfinite preconditioner checks.

`test/test_optimize_krylov_augmentation_images.f90` checks the opt-in
current-Jacobian-image filter on orthogonal solution directions with nearly
parallel images, dependence on multiple images, FIFO retention, rechecking
an older span under a changed Jacobian, invalid and nonfinite inputs, and
huge/subnormal scales. `test/test_optimize_root_krylov_augmentation_images.f90`
checks forwarding through `root` and `newton_krylov`, nonlinear convergence
with analytic and finite-difference products, preconditioning and two-pass
reorthogonalization, exact two-point/three-point extra evaluation counts,
zero-tolerance parity, callback failures, and option validation. The image
filter deliberately requires extra J*v evaluations and is distinct from the
zero-additional-J*v solution-space filter. This snapshot passed 21 focused
programs under GNU Fortran 14.2 with `-O0 -g` and `-O2`, both with bounds
checks, explicit-interface enforcement, and floating-point traps. A test-only
MINPACK shim and public API facade were used; these were not packaged.

`test/test_optimize_krylov_augmentation.f90` exercises the optional LGMRES
history filter directly: combinations of two old vectors, near-duplicate
vectors, finite independent directions, FIFO eviction, preservation of the
stored vectors, zero/nonfinite direction rejection, and huge/subnormal
normalization. `test/test_optimize_root_krylov_augmentation.f90` checks
forwarding through `root` and `newton_krylov`, nonlinear convergence with
reorthogonalization, that `outer_prune_tol=0` exactly matches the default
solution and work counts, and preflight rejection of invalid values and
unsupported methods. The filter adds no Jacobian-vector evaluations, though
it can change subsequent Krylov work by changing the augmentation history.
In this snapshot, 19 focused programs (including the two new filter tests)
passed in both GNU Fortran 14.2 debug and `-O2` builds with bounds checking,
explicit-interface enforcement, and floating-point traps. These builds
used a test-only MINPACK interface shim and a test-only public API
re-export module; neither is part of the project archive. The full FPM
build and external MINPACK functionality are not thereby validated.

`test/test_optimize_root_krylov_reorthogonalization.f90` covers the optional
GMRES/LGMRES DGKS/second-pass orthogonalization modes. On a deterministic
48-variable diagonal system with eigenvalues from 1 down to `1e-12`, the
largest absolute off-diagonal Arnoldi Gram-matrix entry was approximately
`1.00e-7` (one pass), `1.53e-16` (selective), and `1.22e-16` (two passes)
with GNU Fortran 14.2 debug compilation. The three modes use identical numbers
of Jacobian-vector products; the default solution remains identical to an
explicit `'once'` call. This probes basis orthogonality, **not** a claimed
improvement in final nonlinear residual. Additional checks cover LGMRES with
an inverse-Jacobian preconditioner, both public entry points and rejection
of invalid/unavailable combinations. A test-only MINPACK interface shim is
used for the focused root-regression builds; it is not included in the ZIP.

`test/test_optimize_root_krylov_methods.f90` exercises both public entry
points with the four methods. The four-variable systems use an unsymmetric
linear matrix for BiCGSTAB/CGS/TFQMR, a symmetric-indefinite matrix for MINRES,
and a nonlinear diagonal-quadratic residual for all four. The following
counts reproduce `tools/generate_optimize_krylov_iterative_reference.py`
against SciPy 1.17.0 with `rdiff=1e-6`, `inner_maxiter=30`, no line search,
and the stated `fatol`:

| Method | Linear `nit` | Linear `nfev` | Quadratic `nit` | Quadratic `nfev` |
| --- | ---: | ---: | ---: | ---: |
| BiCGSTAB | 2 | 9 | 9 | 34 |
| CGS | 2 | 10 | 8 | 30 |
| MINRES | 2 | 6 | 9 | 27 |
| TFQMR | 2 | 9 | 9 | 37 |

The original six focused Krylov/nonlinear-root tests passed with GNU
Fortran 14.2 at `-O2`, `-std=f2018`, `-fimplicit-none`,
`-Werror=implicit-interface`, `-fcheck=all`, and
`-ffpe-trap=invalid,zero,overflow`. The verified-TFQMR extension adds two
standalone regression programs, extra cases in the public root-method and
custom-solver-hook tests, and a new pinned SciPy reference generator.
The focused root tests use an external **test-only** MINPACK interface shim;
it is not included in the release ZIP. A full dependency-backed `fpm test`
has not been run in this environment.

## LGMRES augmentation-policy comparison (2026-10-10)

`test/test_optimize_root_krylov_pruning_comparison.f90` exercises three
retained-vector policies on **changing** 16-variable nonlinear Jacobians.
For the exact matrix and nonlinear residual definition, see
`benchmarks/README.md`. The dense linear factor has explicitly controlled
singular values between `1` and `1/condition`, with two orthogonal
Householder factors. The nonlinear Jacobian changes along the Newton path;
its condition number is not asserted to be exactly that of the linear
factor. All three policies use an analytic J*v callback, Armijo line search,
`inner_maxiter=12`, `outer_k=4`, selective reorthogonalization,
`fatol=1e-9`, and at most 80 nonlinear iterations. The driver independently
re-evaluates the returned residual at the final iterate and compares
`nfev` with the actual function callback count.

With GNU Fortran 14.2, `-O2`, bounds checks and FP traps, the most difficult
case (nominal linear-factor condition `1e7`) gave:

| Policy | Success | Nonlinear iterations | `nfev` | Analytic J*v calls | Max actual residual |
| --- | --- | ---: | ---: | ---: | ---: |
| None | No | 80 | 159 | 238 | `1.39155e-3` |
| Solution space | Yes | 59 | 113 | 495 | `3.20157e-10` |
| Current Jacobian images | Yes | 59 | 113 | 777 | `5.12575e-10` |

With the same compiler in a checked debug build, the hard-case image-filtered
solve completed in 61 iterations with `2.56418e-18` maximum residual,
versus 59 iterations and `1.27264e-10` for solution-space filtering; the
unpruned case did not converge in 80 iterations. This difference between
build configurations is why the test asserts residual integrity and callback
accounting rather than hardcoded rankings.

At nominal condition `1e3` and `1e4`, all three policies converged in seven
nonlinear checks. For those two analytic-J*v problems, pruning in solution
space used 86 products (the same as no pruning); image-space pruning used
106, without improving convergence. For the `1e3` two-point finite-difference
case, all three policies converged in seven checks. Image-space pruning used
113 residual samples versus 93 for no pruning and solution-space pruning.
The full `-O2` output is in
`benchmarks/lgmres_pruning_gfortran14_o2.csv`.

In that snapshot, **22 of 22 focused Krylov/root Fortran test programs
passed** in both checked configurations. The full FPM build was unavailable.

## Automatic LGMRES Jacobian-image pruning (2026-10-10)

`outer_image_prune_auto=.true.` selects a relative image-span cutoff from
observed first-pass modified Gram-Schmidt projection loss, with an
`epsilon`- and dimension-scaled floor. This is an opt-in heuristic and is
not a replacement for conditioning/error analysis. Two new tests cover
near-dependent versus distinguishable images, dependence on combinations,
FIFO retention, huge and subnormal vectors, nonfinite images, public `root`
and `newton_krylov` forwarding, callback failures, conflicting or unsupported
options, and exact two- versus three-point residual-sample counts. An
explicit `.false.` retains identical solutions, iterations, and evaluation
counts to the default. The deterministic 16-variable comparison now exercises
four policies (none, solution-space, fixed-image, automatic-image). Its GNU
Fortran 14.2 `-O2` results show automatic and fixed-image modes identical on
these particular cases: at nominal condition `1e7`, both converged in 59
nonlinear checks, with 777 analytic J*v calls and a maximum residual of
`5.12575e-10`. This outcome is not asserted as a portable ranking.

The **24 focused Krylov/root programs passed** in both GNU Fortran 14.2
checked debug (`-O0 -g`) and optimized (`-O2`) configurations. Tests used
`-std=f2018`, `-fcheck=all`, `-ffpe-trap=invalid,zero,overflow`, and
`-Werror=implicit-interface`. The focused dependency stubs were external
and not included in the ZIP. The full `fpm test` with real dependencies
remains unrun.

## Reusing current-iterate LGMRES images (2026-10-10)

`outer_image_reuse=.true.` is an opt-in optimization for LGMRES with
Jacobian-image pruning. Previously the current-iterate pruning pass and the
following Arnoldi solve each evaluated the retained vectors' same `J(x)*v`.
The new inner-solve-only image cache eliminates the repeated calls. The
original retained solution vectors are passed unchanged to the cached
operator, and no image is retained across nonlinear iterates. The default
`outer_image_reuse=.false.` keeps prior work-count behavior.

The new `test_optimize_root_krylov_image_reuse` verifies analytic callback
savings with a nonlinear, changing Jacobian, fixed and automatic filtering,
preconditioning, reorthogonalization, one- and two-probe numerical products,
unchanged disabled-default behavior, and rejection of invalid combinations
before the first residual evaluation. The existing 16-variable stress test now
compares five policies, including automatic image pruning with reuse. It checks
that residual counts still agree with actual callback invocations; the test
does not require identical hard-case trajectories under the different policies.
The current build results and work counts are recorded in the accompanying
`benchmarks/lgmres_pruning_gfortran14_o2.csv`.

With GNU Fortran 14.2 `-O2`, the deterministic 16-variable comparison yielded:

| Nominal linear-factor condition | Automatic image pruning | With image reuse |
| --- | ---: | ---: |
| `1e3`, analytic J*v calls | 106 | 92 |
| `1e4`, analytic J*v calls | 106 | 92 |
| `1e7`, analytic J*v calls | 777 | 553 |
| `1e3`, two-point residual evaluations | 113 | 99 |

All these cases converged under both policies, although the most ill-conditioned
case followed a slightly different floating-point path; the comparison
regression deliberately checks the true final residual instead of identical
iterates. One-dimensional public API regressions saved two forward-difference
samples and four centered-difference samples during three nonlinear iterations.
**All 25 focused root/Krylov tests passed** under both debug (`-O0 -g`) and
optimized (`-O2`) builds, each with `-std=f2018`, `-fcheck=all`,
`-ffpe-trap=invalid,zero,overflow`, and `-Werror=implicit-interface`.
The MINPACK and top-level facade shims were build-only and are not shipped.
The complete dependency-backed `fpm test` was not run.

## LGMRES candidate-image reconstruction and scaling (2026-10-10)

Added opt-in `outer_candidate_image_reuse` for built-in LGMRES with fixed or
automatic Jacobian-image pruning and an analytic `jacobian_product` callback.
The raw `J*z` Arnoldi columns, saved before Gram-Schmidt, reconstruct the
candidate Newton step image by linear combination. The pruning code compares
*normalized* images, so no multiplication by the original RHS norm is needed.
The option does not reuse stale images across changing nonlinear Jacobians,
does not suppress explicit inner residual checks, and does not accept
finite-difference Jacobian products. Its default is disabled.

The new `test_optimize_root_krylov_candidate_image_reuse.f90` covers
changing nonlinear Jacobians, automatic and fixed image pruning, both public
entry points, independent and combined cache modes, unchanged disabled mode,
actual analytic callback savings, and invalid combinations rejected before
residual evaluation. The previous 16-variable ill-conditioning regression
now compares seven analytic policies plus five finite-difference policies.
For its `1e7` case, both caches together require 495 analytic J*v calls,
versus 777 without caching and 553 with retained-image caching alone; all
three image-filtered configurations converge in 59 nonlinear iterations.
This is a deterministic numerical example, not a universal speed guarantee.

The new `benchmark_lgmres_image_cache.f90` uses systems of 32, 64, and 128
variables and runs each caching policy 30 times. GNU Fortran 14.2 with `-O2`
and runtime checks produced the following analytic J*v call counts for
no caching versus both caches: 694 vs 534 (n=32), 1450 vs 950 (n=64),
and 967 vs 694 (n=128). All converged to maximum residual less than `1e-8`.
The last two problem sizes sometimes followed different nonlinear iteration
counts with retained-image caching. Measured timings are reference data only;
see `benchmarks/lgmres_image_cache_scaling_gfortran14_o2.csv`.

**All 26 focused root/Krylov regression programs passed** in GNU Fortran
14.2 checked debug (`-O0 -g`) and optimized (`-O2`) configurations with
`-std=f2018 -fcheck=all -ffpe-trap=invalid,zero,overflow`
and `-Werror=implicit-interface`. External MINPACK and API stubs used only
for focused linking are not in the project ZIP. A complete dependency-backed
`fpm test` has not been run.


## LGMRES image-cache compact workspace (2026-10-10)

The image-space pruning pass now compacts retained direction columns in place,
using a separate one-vector candidate to avoid Fortran dummy-argument aliasing.
When both `outer_image_reuse` and `outer_candidate_image_reuse` are enabled,
Arnoldi no longer duplicates the retained image-cache columns in its raw
operator-product table. It continues to store the actual products for the
remaining ordinary Krylov columns; no change was made to the default SciPy
path. A trial implementation based on rounded Arnoldi/Hessenberg coefficients
was **rejected** because it changed convergence on the ill-conditioned
`1e7` stress case. The retained-product sharing version preserves that
case's 59-iteration convergence and 495 analytic J*v calls with both caches.

`test_optimize_root_krylov_cache_compaction.f90` covers combined caching,
preconditioning, two-pass reorthogonalization, fixed and automatic image
pruning, one-column FIFO history, callback counts, and independent residuals.
All **27** focused root/Krylov test programs passed under GNU Fortran 14.2
with `-std=f2018 -fcheck=all -ffpe-trap=invalid,zero,overflow`
and `-Werror=implicit-interface`, in both debug (`-O0 -g`) and checked
optimized (`-O2 -g`) builds. The dependency-free focused builds used a
*temporary* MINPACK interface shim and a facade for unrelated dependencies;
those test-only shims are excluded from the archive.

The new O(n) matrix-free benchmark reports matching six nonlinear iterations,
six residual evaluations, and 46 Jacobian-vector calls for the previous and
updated versions. At 250,000 variables, measured process peak RSS decreased
from 184,992 KiB to 145,812 KiB. At 500,000 variables, peak RSS decreased
from 337,112 KiB to 294,192 KiB. These are single-run GNU Fortran 14.2
Linux references (checked `-O2`), not portable memory guarantees.
See `benchmarks/lgmres_cache_workspace_gfortran14_o2.csv` and the benchmark
README for reproduction. A full dependency-backed `fpm test` was not run
because FPM was unavailable.

## GMRES/LGMRES per-solve Arnoldi workspace reuse (2026-10-11)

`scifort_optimize_nonlin` now keeps its Arnoldi/Hessenberg/Givens and
small least-squares storage in a solver-local workspace across nonlinear
iterations instead of repeatedly allocating those arrays. The upcoming
unused Arnoldi basis column doubles as `J*v`/orthogonalization scratch.
For LGMRES, retained augmentation directions are read directly from the
existing history; only the `inner_maxiter` ordinary Krylov solution columns
are kept in the `z` array. The exact pre-orthogonalization Jacobian products
are retained when candidate-image reconstruction is requested. GMRES reuses
a finished basis column for its extra true linear-residual evaluation.

The new `test_optimize_root_krylov_workspace_reuse.f90` exercises 24- and
64-variable systems, GMRES, no-cache and combined-cache LGMRES, changing
augmentation-history capacities, multi-step convergence, independent true
residuals, and callback accounting. **All 28 focused root/Krylov programs
passed** in GNU Fortran 14.2 checked debug (`-O0 -g`) and optimized
(`-O2 -g`) builds using `-std=f2018 -fcheck=all`, explicit-interface
enforcement, and invalid/zero/overflow floating-point traps. The earlier
16-variable stress comparison matches the preceding snapshot *byte for
byte*, including the difficult nominal-condition-`1e7` problem.

The existing large, matrix-free memory benchmark was compiled separately
against the original compact-image-cache ZIP and this updated source. Single
GNU Fortran 14.2 `-O2` checked runs produced:

| Variables | Previous RSS KiB | New RSS KiB | New/old J*v calls |
| ---: | ---: | ---: | ---: |
| 250,000 | 145,820 | 141,876 | 46 / 46 |
| 500,000 | 294,168 | 282,528 | 46 / 46 |

Both versions had six nonlinear iterations, six residual evaluations and
the same independently reported maximum residual at each size. RSS/timing
are environment-specific whole-process measurements, not guaranteed
allocations or universal performance improvements. Raw observations are in
`benchmarks/lgmres_arnoldi_workspace_gfortran14_o2.csv`. The focused builds
used test-only MINPACK and top-level API facade shims **outside** the project
archive. Full external-dependency `fpm test` was not run.

## Lazy Krylov workspaces and large-history comparison (2026-10-11)

`test/test_optimize_root_krylov_lazy_workspace.f90` checks that the same
nonlinear problem with `outer_k=8` and `outer_k=96` produces the same
solution, nonlinear iteration count, residual evaluations, Jacobian-vector
callback count, and preconditioner callback count. It exercises both plain
LGMRES and automatic image pruning with retained/candidate-image reuse,
plus preconditioned GMRES. The pre-existing Krylov tests cover full and
one-vector histories, finite-difference evaluation counts, and near-singular
Jacobian systems.

All 29 focused Krylov/root test programs passed with GNU Fortran 14.2 in
`-O0 -g` and `-O2` builds, using `-std=f2018`, bounds checks, floating-point
traps, and explicit-interface enforcement. Temporary MINPACK interface
shims were used only for focused linking; they are not in the source archive.
These checks are not a full dependency-backed `fpm test`, which remains
outstanding.

The matrix-free memory benchmark has optional `outer_k` and preconditioner
arguments. See `benchmarks/lgmres_lazy_workspace_gfortran14_o2.csv` and
`benchmarks/README.md` for matched previous/updated peak-RSS and numerical
results. The table is descriptive, not a fixed performance requirement.


## Reusable LGMRES pruning workspaces (2026-10-11)

`test/test_optimize_krylov_augmentation_workspace.f90` compares the original
allocation-per-call filters against a reusable workspace and caller-owned
Arnoldi scratch basis. It exercises fixed/automatic Jacobian-image thresholds,
solution-space filtering, FIFO retention, rejected vectors, and changes in
problem dimension. Retained directions and insertion decisions agree exactly.

All 30 focused Krylov/root programs passed under GNU Fortran 14.2 with both
`-O0 -g` and `-O2`, `-std=f2018`, `-fcheck=all`,
`-Werror=implicit-interface`, and `-ffpe-trap=invalid,zero,overflow`.
Temporary MINPACK and public-module test-only shims (both outside the archive)
were used for the focused public-root tests. These checks do not replace a
full external-dependency `fpm test`, which remains outstanding.

The existing `benchmark_lgmres_cache_memory.f90` was compiled separately
against the prior `lgmres-lazy-workspace` snapshot and the revised sources
with identical GNU Fortran 14.2 `-O2` runtime-check flags. For `outer_k=64`:

| Variables | Previous peak RSS KiB | Revised peak RSS KiB | Nonlinear iterations | J*v calls |
| ---: | ---: | ---: | ---: | ---: |
| 150,000 | 151,220 | 145,392 | 6 / 6 | 46 / 46 |
| 250,000 | 251,264 | 241,568 | 6 / 6 | 46 / 46 |

The maximum residual and nonlinear residual-evaluation counts were unchanged.
The measured peak-RSS reductions were about 3.9% in these single Linux runs;
whole-process RSS and CPU times depend on allocation behavior and hardware.
See `benchmarks/lgmres_pruning_workspace_gfortran14_o2.csv` for the observations
and `benchmarks/benchmark_lgmres_pruning_workspace.f90` for a separate
standalone repeated-pruning timing experiment. Small-vector timing effects
are not consistently positive across all problem sizes.

The optimized 16-variable LGMRES pruning stress-comparison output also matched
the previous snapshot byte for byte, including all reported work counts.


## Incremental LGMRES Jacobian-image span screening (2026-10-11)

During the recheck of retained augmentation directions at a nonlinear point,
`scifort_optimize_krylov_augmentation` now extends a two-pass orthonormal
image basis as images arrive, rather than rebuilding that basis from all
previously accepted images on every candidate. The basis is reset when the
current nonlinear point changes. The original FIFO insertion routine and
SciPy-compatible unpruned path have not changed.

`test/test_optimize_krylov_augmentation_stream.f90` compares fixed and
automatic thresholds against the reconstruct-on-each-call helper, including
near-dependent and null images, very large and subnormal values, stored-vector
identity, scratch reuse, and resets with a changed Jacobian. Accepted
vectors and decisions matched exactly. All 31 focused root/Krylov programs
passed with GNU Fortran 14.2 under both `-O0 -g` and `-O2`, with
`-std=f2018`, `-Werror=implicit-interface`, bounds checking, and
floating-point traps. Temporary MINPACK and public-module test-only shims
were used for focused linking but are not packaged. The 16-variable
changing-Jacobian stress-comparison output of the optimized build matched
the previous snapshot byte for byte, including callback counts.

The standalone `benchmarks/benchmark_lgmres_incremental_pruning.f90`
measures only the repeated span-screening step for 96 proposed vectors of
length 192 (84 retained), repeating the pass eight times. It uses shared
scratch for both paths, so the difference is rebuilding the orthonormal
basis, not per-call allocation. One GNU Fortran 14.2 checked `-O2` run gave:

| Repeated screening method | Total CPU time (s) |
| --- | ---: |
| Reconstruct basis each candidate | 0.55291 |
| Incremental image basis | 0.01994 |

This is about 27.7 times faster *for that pruning microbenchmark*, not for
the complete root solver. Results are environment-dependent and are not
fixed regression thresholds; input data and flags are documented in
`benchmarks/README.md`.

## Incremental LGMRES solution-space span screening (2026-10-11)

Added `retain_lgmres_direction_stream` and
`reset_lgmres_direction_pruning_stream` to the internal augmentation module.
The new append-only routine keeps a two-pass orthonormal solution-space
basis across accepted candidates. It intentionally does not perform FIFO
eviction: when `outer_count` reaches `outer_k`, the root solver falls back to
the existing rank-reconstructing insertion routine. The cache is reset for
new solves. The streaming path is limited to eight MiB of extra basis storage
for each solve, preserving the previous behavior for large matrix-free cases.

`test/test_optimize_krylov_augmentation_solution_stream.f90` compares the
incremental and reconstruction implementations across multiple candidate
orders and tolerances, near dependence below the `32*epsilon` floor, exact
duplicates, zeros, large and subnormal magnitudes, reset behavior, and the
transition to full-history FIFO replacement. Retention decisions and stored
vectors were identical for all tested cases. All 32 focused root/Krylov
programs passed with GNU Fortran 14.2 in checked debug (`-O0 -g`) and
optimized (`-O2`) builds with `-std=f2018`, `-Werror=implicit-interface`,
`-fcheck=all` and floating-point traps for invalid, zero and overflow.
Focused public-root linking used test-only MINPACK and public-interface
shims; a full dependency-backed `fpm test` was not run.

The prior 16-variable changing-Jacobian pruning-comparison output matched
the previous snapshot byte for byte, including convergence flags, residuals,
nonlinear iterations and Jacobian-vector evaluations.

The standalone `benchmarks/benchmark_lgmres_solution_span.f90` times 96
proposed solution directions of length 192, with 70 retained, repeated 12
times. Both paths reuse ordinary scratch; this isolates repeated basis
reconstruction. A single GNU Fortran 14.2 `-O2`, checked Linux run gave:

| Solution-space pruning operation | Total CPU time (s) |
| --- | ---: |
| Reconstruct basis for each candidate | 0.505688 |
| Incremental basis, append-only | 0.035134 |

This is a 14.4x speedup for the isolated screening loop, not the nonlinear
solver. Timing is environment dependent and is not a CI pass/fail threshold.
The benchmark and instructions are in `benchmarks/`.

## LGMRES full-history FIFO QR downdating (2026-10-11)

`outer_fifo_prune_fast=.true.` enables an opt-in Givens QR-downdate path for
solution-space LGMRES augmentation pruning once `outer_k` is full. It accepts
only built-in LGMRES with `outer_prune_tol > 0` and `outer_k >= 2`. A factor
that loses numerical rank falls back to the previous rebuild-for-each-insert
path; an otherwise well-conditioned factor is refreshed after 32 accepted
FIFO replacements. The existing default is unchanged, and large problems
are constrained by the same eight-MiB basis limit as append-only screening.

`test/test_optimize_krylov_augmentation_fifo_qr.f90` compares 240 sequential
candidate decisions with the legacy filter, including many FIFO evictions,
exact duplicates, zeros, extreme magnitudes, reset/reordered histories, and
nearly dependent full histories that exercise the safe fallback. Retained
normalized directions agree exactly. The separate public regression
`test/test_optimize_root_krylov_fifo_qr.f90` exercises `root` and
`newton_krylov`, checks computed residuals independently, and verifies that
invalid keyword combinations are rejected before calling the residual.

All **34 focused root/Krylov tests** passed with GNU Fortran 14.2, both
`-O0 -g` and `-O2`, using `-std=f2018`, `-Werror=implicit-interface`,
`-fcheck=all`, and floating-point traps for invalid, zero, and overflow.
Focused linking used test-only MINPACK and public-interface shims; a full
FPM dependency-backed build was unavailable. On the existing 16-variable
changing-Jacobian comparison, enabling FIFO updates for solution-space
pruning produced identical reference rows (including the nominal `1e7`
condition-number case) in the optimized checked build.

The optional `benchmarks/benchmark_lgmres_fifo_pruning.f90` uses 192-value
directions, `outer_k=48`, 200 proposals, and 12 repetitions, with some
rejected duplicates. One GNU Fortran 14.2 `-O2` Linux checked run gave:

| Full-history augmentation pruning | Total CPU time (s) |
| --- | ---: |
| Rebuild complete QR span each proposal | 1.437279 |
| Incremental FIFO QR update with periodic rebuilding | 0.184281 |

This is a **7.8x speedup for the isolated pruning loop**, not the entire
root solver. CPU timings and convergence trajectories vary across systems;
these measurements are not CI acceptance thresholds.

## Reusing the current-Jacobian LGMRES image basis (2026-10-11)

`outer_image_basis_reuse=.true.` avoids reconstructing the image-space
orthonormal span for the candidate Newton direction. It is opt-in for
built-in LGMRES with fixed or automatic Jacobian-image pruning. The cache
is valid for only one insertion and is reset as the nonlinear Jacobian
changes. An 8 MiB per-solve basis cap selects the legacy code path for
large matrix-free systems; the default is unchanged.

`test/test_optimize_krylov_augmentation_candidate_basis.f90` compares
candidate acceptance, retained normalized vectors, full-history FIFO
replacement, near-dependence, extreme magnitudes, and one-shot invalidation
with the previous reconstruction. The public regression
`test/test_optimize_root_krylov_image_basis_reuse.f90` checks fixed and
automatic image pruning, retained-image caching, analytic and finite-
difference products, work counts, and option preflight. The expanded
`test/test_optimize_root_krylov_pruning_comparison.f90` verifies the same
convergence status, nonlinear iterations, residual calls, analytic J*v calls,
and solution accuracy on three changing-Jacobian 16-variable cases with
nominal conditioning from 1e3 to 1e7.

With GNU Fortran 14.2 on Linux, 36 focused root/Krylov tests passed in each
of the checked `-O0 -g` and `-O2` configurations, using `-std=f2018`,
`-Werror=implicit-interface`, `-fcheck=all`, and floating-point traps.
The focused public-root tests used test-only MINPACK and thin public-API
shims, which are not included in the project ZIP; a complete dependency-
backed `fpm test` was not run.

The standalone 192-variable, 48-history, 250-repetition benchmark measures
only image-span filtering and candidate insertion; it excludes evaluating
`J*v` and the Arnoldi solve. The single GNU Fortran 14.2 `-O2` checked run
measured 0.383 s versus 0.210 s for fixed tolerance and 0.367 s versus
0.199 s with the automatic cutoff, approximately 1.8x in both cases.
These timings are environment dependent and not CI acceptance criteria.


## In-place LGMRES image-basis workspace and end-to-end benchmark (2026-10-11)

`test/test_optimize_krylov_image_basis_storage.f90` compares the original
private QR basis with a caller-owned basis stored in the existing image
buffer. It covers fixed/automatic thresholds, dependent Jacobian images,
full FIFO candidate insertion, rejection, and resets. Both implementations
retain exactly the same normalized directions and decisions. The new
`basis_work` argument to the low-level cached-candidate routine is optional.

GNU Fortran 14.2 passed all 37 focused root/Krylov test programs in debug
and `-O2` checked builds. A test-only MINPACK module, root facade, and
reference fixture were compiled separately; none is shipped as a project
module. Full dependency-backed `fpm test` remains unrun.

`benchmarks/benchmark_lgmres_inplace_image_basis.f90` exercises the entire
nonlinear LGMRES solve with analytic Jacobian products, automatic image
pruning, and optional image-basis reuse, deliberately leaving
`outer_image_reuse` off. Compile the same benchmark against both the
previous image-basis snapshot and the updated source tree. The paired Linux
`/usr/bin/time` RSS measurements below used GNU Fortran 14.2 `-O2`,
`-fcheck=all`, and floating-point traps, with four repeated solves per run,
and two process runs per source variant. Results are in
`benchmarks/lgmres_inplace_image_basis_gfortran14_o2.csv`.

| Variables | Previous mean peak RSS (MiB) | In-place mean peak RSS (MiB) |
| ---: | ---: | ---: |
| 30,000 | 31.4 | 30.3 |
| 100,000 | 101.3 | 97.5 |
| 250,000 | 251.2 | 241.7 |

Every configuration converged in six nonlinear iterations per solve, with
six residual and 61 analytic Jacobian-vector evaluations per solve. All
repeated runs matched the previous source's checksum and maximum residual
for the same problem size. CPU times did not show a consistent gain; this
change is primarily a memory optimization. These measurements depend on the
compiler, allocator, OS, and workload and are not CI thresholds.


## Packed LGMRES retained-image buffers (2026-10-11)

`test/test_optimize_root_krylov_packed_image_buffers.f90` exercises the
zero-column initial image history, `outer_k=1` FIFO replacement, growing
histories, automatic and fixed image thresholds, combined retained/candidate
image reuse, and centered finite differences. It compares the cached analytic
solve with the same uncached solver for bitwise solutions, iteration counts,
and residual evaluations. All 38 focused root/Krylov programs passed under
GNU Fortran 14.2 with `-O0 -g` and `-O2`, `-std=f2018`, `-fcheck=all`,
`-Werror=implicit-interface`, and floating-point traps. The previous
16-variable changing-Jacobian stress output matched byte for byte. A full
FPM dependency-backed test remains outstanding.

The end-to-end matrix-free benchmark `benchmark_lgmres_cache_memory.f90` was
compiled from each source snapshot separately with identical checked `-O2`
flags. Analytic Jacobian products, automatic image pruning, retained-image
reuse, candidate-image reuse, `inner_maxiter=10`, and `outer_k=8` were enabled.
Two Linux `/usr/bin/time -f 'peak_kib=%M'` process runs per size produced:

| Unknowns | Previous peak RSS (MiB) | Updated peak RSS (MiB) | Reduction |
| ---: | ---: | ---: | ---: |
| 100,000 | 52.4 | 50.9 | 2.9% |
| 250,000 | 129.1 | 125.4 | 2.9% |
| 500,000 | 257.0 | 249.3 | 3.0% |

Each run converged in six nonlinear iterations with six residual and 46
analytic Jacobian-vector evaluations. The benchmark comparison is logged in
`benchmarks/lgmres_packed_images_gfortran14_o2.csv`. Peak RSS is sensitive
to the OS allocator and machine; no minimum saving or CPU speedup is asserted.


## Bounded LGMRES raw-product cache (2026-10-11)

`test/test_optimize_root_krylov_candidate_cache_budget.f90` covers
unlimited and sufficient budgets (identical Jv counts and solutions),
a deliberately small budget (the candidate image is recomputed through the
analytic Jv callback), the history crossing its budget between nonlinear
iterations, both retained-image caching policies, `newton_krylov` forwarding,
and rejection of unsupported/negative configurations before the residual
callback. All 39 focused root/Krylov regressions passed under GNU Fortran
14.2 with `-O0 -g` and `-O2`, `-std=f2018`, bounds checking and traps for
invalid floating-point operations, division by zero, and overflow. Focused
public-root linking used a test-only MINPACK placeholder and facade which
are **not** included in the project ZIP; full FPM/dependency tests remain
outstanding.

`benchmarks/benchmark_lgmres_candidate_cache_budget.f90` compares a
zero/unlimited cap with a 1 MiB cap on the same matrix-free nonlinear
problem. Runs used GNU Fortran 14.2 `-O2 -fcheck=all` and floating-point
traps, and `/usr/bin/time -f 'peak_kb=%M elapsed_s=%e'`. Results in
`benchmarks/lgmres_candidate_cache_budget_gfortran14_o2.csv`:

| Variables | Unlimited peak RSS (MiB) | 1 MiB cap peak RSS (MiB) | Jv calls (unlimited/bounded) |
| ---: | ---: | ---: | ---: |
| 100,000 | 70.7 | 55.5 | 57 / 62 |
| 250,000 | 174.9 | 136.6 | 56 / 61 |

Both policies took 6 nonlinear iterations and 6 residual evaluations for
each size, with matching reported maximum residuals. The bounded option
saves raw-image storage but adds one exact Jv callback per applicable
nonlinear step. These are single-run, platform-dependent peak-RSS values,
not guarantees or end-to-end runtime speedup claims.

### GMRES restarted-subspace extension (October 2026)

`test/test_optimize_root_krylov_gmres_restart.f90` checks the optional
`inner_restart` keyword through `root` and `newton_krylov`, right
preconditioning, different reorthogonalization modes, original-path work
counts, finite-difference accounting, and invalid-argument rejection.
Together with the existing focused root/Krylov tests, **39 of 39** selected
programs passed with GNU Fortran 14.2 in both checked debug and `-O2` builds.
The offline focused build used local MINPACK interface stubs that were
**not added to the project**; it is not a replacement for `fpm test`.

For the end-to-end benchmark, compile
`benchmarks/benchmark_krylov_restarted_gmres.f90` against the SciFort root
modules and dependencies, then run `benchmark_gmres_restart full 100000`
and `benchmark_gmres_restart restart 100000`. The same binary uses
`inner_maxiter=64` with no restarts or with `inner_restart=8`.
A checked Linux GNU Fortran 14.2 `-O2` run using `/usr/bin/time` gave:

| Unknowns | Mode | Peak KiB | CPU seconds | Jv calls | Iterations | Residuals |
|---:|---|---:|---:|---:|---:|---:|
| 50,000 | full | 56,472 | 0.0440 | 26 | 3 | 3 |
| 50,000 | restart | 13,476 | 0.0683 | 101 | 4 | 4 |
| 100,000 | full | 111,472 | 0.0800 | 26 | 3 | 3 |
| 100,000 | restart | 25,612 | 0.1458 | 101 | 4 | 4 |

Both versions converged to a true maximum residual below `1e-10`.
These are single-run process-memory and timing measurements, not universal
performance expectations. The iteration and callback differences reflect
restarted GMRES convergence behavior, not just the storage optimization.
The raw measurements are in
`benchmarks/krylov_gmres_restart_gfortran14_o2.csv`.


### LGMRES memory-bounded restart extension (October 2026)

`test/test_optimize_root_krylov_lgmres_restart.f90` covers public forwarding
through `root` and `newton_krylov`, unchanged results and work counts for
`inner_restart >= inner_maxiter`, matrix-free analytic and finite-difference
Jacobians, right preconditioning, changing nonlinear Jacobians, current-Jacobian
image pruning and image reuse, true inner-residual checks, rejection of
invalid arguments, and callback failures. The earlier GMRES restart test now
also verifies that built-in LGMRES accepts the restart argument.

The 40 selected Krylov/root tests passed in both GNU Fortran 14.2 `-O0` and
`-O2` checked builds with bounds checking, floating-point traps, and explicit
interface checking. This was an offline-focused build with temporary MINPACK
interface shims; it is **not** a complete `fpm test` and does not validate the
real external dependencies.

`benchmarks/benchmark_krylov_restarted_lgmres.f90` compares `inner_maxiter=64`
with full-width LGMRES and `inner_restart=8`, retaining `outer_k=4` vectors.
The paired GNU Fortran 14.2 `-O2` Linux measurements were:

| Unknowns | Variant | Peak KiB (two trials) | Jv calls | Nonlinear iterations |
|---:|---|---|---:|---:|
| 20,000 | Full | 24,204; 24,292 | 25 | 3 |
| 20,000 | Restarted | 6,968; 6,980 | 32 | 3 |
| 100,000 | Full | 116,252; 116,268 | 90 | 4 |
| 100,000 | Restarted | 30,188; 30,184 | 119 | 4 |
| 250,000 | Full | 288,484; 288,436 | 90 | 4 |
| 250,000 | Restarted | 73,484; 73,512 | 105 | 4 |

All runs converged; each mode's residual evaluation count is included in
`benchmarks/krylov_lgmres_restart_gfortran14_o2.csv` along with process CPU
times and the final true residual. Restarting is a memory-versus-work
tradeoff; elapsed times and nonlinear trajectories depend on the problem.

### Adaptive restarted-LGMRES extension (October 2026)

`test/test_optimize_root_krylov_adaptive_restart.f90` verifies adaptive
restart growth with true-residual-based stagnation detection, the default
fourfold ceiling, explicit capped-width equivalence to fixed restarting,
public `root` and `newton_krylov` forwarding, invalid option rejection
before residual evaluation, residual/Jacobian callback counts, and a
nonlinear system with preconditioning and Jacobian-image pruning. Together
with the prior focused root/Krylov suite, **41/41** test programs passed
under GNU Fortran 14.2 with `-O0` and `-O2`, bounds checking, floating-point
traps, and explicit-interface checking. The local MINPACK/public-interface
stubs used by the offline tests are not included in the distribution.
This does **not** replace the full dependency-backed `fpm test`, which was
not run in this environment.

`benchmarks/benchmark_krylov_adaptive_restart.f90` reproduces the following
optimized results; the raw records are in
`benchmarks/krylov_adaptive_restart_gfortran14_o2.csv`:

| System | Restart policy | Converged | True max residual | Analytic Jv calls |
|---|---|---|---:|---:|
| 2-variable skew | Fixed width 1 | No | 2.00 | 1 |
| 2-variable skew | Adaptive 1 to 2 | Yes | 4.44e-16 | 4 |
| 2-variable skew | Adaptive capped at 1 | No | 2.00 | 1 |
| 12-variable scaled diagonal | Fixed width 1 | No | 2.50e-2 | 360 |
| 12-variable scaled diagonal | Adaptive 1 to 4 | Yes | 1.53e-10 | 182 |
| 12-variable scaled diagonal | Adaptive capped at 1 | No | 2.50e-2 | 360 |

These are deterministic small tests of the adaptive algorithm, not
representative performance guarantees or a general comparison against SciPy.
In general a larger restart width may consume more memory, and some systems
can still stagnate or require more Jacobian-vector products.


### Large, ill-conditioned restarted-LGMRES comparison (2026-10-11)

`test/test_optimize_root_krylov_restart_scaling.f90` supplements the earlier
adaptive-restart regression with a family of changing-Jacobian problems whose
linear factor is an orthogonally rotated diagonal matrix. At the root the
2-norm Jacobian condition number is exactly the prescribed value. The test
verifies analytic Jv against centered differences, the actual returned
residual, a custom scaled convergence norm, preconditioner dispatch,
residual callback counts and solution accuracy at 64 and 256 variables.

With the new test, all **42 focused root/Krylov test programs passed** in
both GNU Fortran 14.2 debug `-O0` and optimized `-O2` configurations with
`-fcheck=all`, `-Werror=implicit-interface`, and
`-ffpe-trap=invalid,zero,overflow`. The offline test compilation used
MINPACK and public-interface facades for unrelated dependencies; none of
these test-only files is shipped in the archive. This is **not** a full
FPM dependency-backed build or `fpm test`.

The accompanying `benchmarks/benchmark_krylov_restart_scaling.f90` was run
in 78 configurations, with raw reference records saved in
`benchmarks/krylov_restart_scaling_gfortran14_o2.csv`. The benchmark checks
reported residuals against an independently recomputed nonlinear residual
and checks `nfev` against the callback counter. In 65,536 dimensions at
nominal condition `1e6`, full width/fixed/adaptive used 59,832/17,912/24,004
KiB peak RSS and 164/301/276 analytic Jv calls, respectively; all met
`fatol=1e-8`. At nominal condition `1e8`, all three returned success but
the forward errors were still approximately `0.0053`-`0.0095`, emphasizing
that a small absolute residual need not imply an accurate root for a badly
conditioned system. A scaled `tol_norm` may prevent such premature success,
although it may also require more work or iterations. Measured RSS and CPU
time are environment-specific, not correctness requirements.

### Automatic Krylov diagonal-preconditioner validation (2026-10-11)

`test/test_optimize_root_krylov_auto_diagonal.f90` adds analytic and
finite-difference diagonal-estimation tests to the existing root/Krylov
suite. It validates zero-diagonal identity fallback, a user-configurable
relative inverse floor, MINRES positive preconditioning for a symmetric
indefinite diagonal Jacobian, rejected conflicting options, failed analytic
Jv callbacks, and exact Jacobian evaluation accounting for refresh intervals
of 0, 1 and 2. On a coupled nonlinear system, an automatic initial-only
estimate agrees with a supplied inverse diagonal without additional
residual evaluations.

**44/44 focused test programs** passed with GNU Fortran 14.2 in both debug
`-O0` and optimized `-O2` configurations, using `-fcheck=all`,
`-Werror=implicit-interface` and `-ffpe-trap=invalid,zero,overflow`.
The focused offline build used temporary MINPACK and public-root compatibility
stubs, which were *not* included in the archive; a full external-dependency
`fpm test` remains outstanding.

`benchmarks/krylov_auto_diagonal_gfortran14_o2.csv` gives eight reference
cases (32 and 128 variables, four policies). On the 128-variable coupled
case, unpreconditioned LGMRES exhausted its iteration budget with 1,551
analytic Jv calls; a supplied inverse diagonal converged with 18, an
initial-only automatic estimate with 146, and per-iteration automatic
updates with 530. All three preconditioned cases had maximum forward error
below `1e-14`. This illustrates the tradeoff between ease of use, update
accuracy and diagonal-estimation overhead, not a general performance claim.

### Krylov seeded diagonal probing and exact fallback (2026-10-11)

`test/test_optimize_root_krylov_random_diagonal.f90` exercises seeded
Rademacher probes, unchanged deterministic Walsh behavior, reproducibility,
Fortran RNG isolation, direct diagonal inspection, batch-consistency
checks, spot-check tolerances, exact fallback, bad inputs and failed
Jacobian callbacks. Two-point and centered finite-difference evaluation
counts are checked separately. The benchmark
`benchmarks/krylov_random_diagonal_gfortran14_o2.csv` records deterministic
operator counts and measured errors for dense/irregularly coupled systems.

Focused builds use a temporary MINPACK module for unrelated hybrid and LM
methods, not included in the release. A full dependency-backed `fpm test`
has not been run in this offline environment.

The completed October 11 continuation passed **47/47 focused root/Krylov
programs** in both GNU Fortran 14.2 `-O0 -g` and `-O2` configurations with
`-std=f2018`, `-fcheck=all`, `-Werror=implicit-interface`, and
`-ffpe-trap=invalid,zero,overflow`. Temporary test-only MINPACK and umbrella
module facades were compiled outside the release tree. The earlier Walsh
snapshot and the revised default configuration produced byte-for-byte
identical output in the established ill-conditioned LGMRES pruning comparison.
The full dependency-backed `fpm test` has **not** been run.

### Krylov progressive verified diagonal probes (2026-10-12)

`test/test_optimize_root_krylov_progressive_diagonal.f90` verifies the new
`inner_m_diag_probe_max=` option through both public Krylov entry points.
A deterministic stride-eight Jacobian case first fails eight-probe
consistency verification, then accepts sixteen probes plus one coordinate
spot check (17 analytic products total) rather than exact fallback
(8+64 products). A more strongly coupled case fails after sixteen probes
and correctly uses 64 exact coordinate probes (80 total). A budget equal
to the initial probe count leaves the old estimator unchanged. Tests also
check two-point and centered finite-difference evaluation counts, callback
failure during expansion, error messages and invalid-option preflight.

The included dense rank-one nonlocal-coupling benchmark covers six sizes
(64 through 1,024), two coupling magnitudes and four probe policies, for
**48 recorded configurations** in
`benchmarks/krylov_progressive_diagonal_gfortran14_o2.csv`. At n=1,024 and
coupling strength 0.50, fixed verified eight-probe estimation fell back to
coordinate estimation and used 1,032 Jv calls; progressive probing used 36,
with mean/max inverse-diagonal relative errors 0.038/0.333. The latter is
approximate: verification only checks selected coordinates and probe batch
stability. Do not interpret this as a certified diagonal or a universal
solver-speed improvement. The benchmark measures estimation callback counts,
not completed nonlinear-solve convergence or CPU acceleration.

The root/Krylov-focused selection completed **33/33 tests** with GNU
Fortran 14.2 under both `-std=f2018 -O0 -g` and `-std=f2018 -O2`, plus
`-fcheck=all`, `-Werror=implicit-interface`, and floating-point traps.
MINPACK and the umbrella optimizer module were provided by temporary
**test-only facades**, omitted from the archive. This is a focused test
subset, not a complete dependency-backed `fpm test`.

### Krylov progressive spot-check reuse (2026-10-12)

`src/optimize/nonlin.f90` now caches exact diagonal coordinate checks during
progressive Rademacher verification. A checked coordinate is not recomputed
at later probe budgets or in the exact fallback. In contrast, omitting
`inner_m_diag_probe_max=` still follows the existing fixed-budget path and
retains the same evaluation counts.

The extended `test/test_optimize_root_krylov_progressive_diagonal.f90` covers
batch-consistent yet biased eight-probe estimates, growth to sixteen probes,
repeated spot checks, two distinct cached coordinates, exhausted-budget
exact fallback, and 2-point/3-point residual-evaluation counts. In the
single-coordinate example, adaptive acceptance takes 17 rather than 18
analytic Jacobian-vector products. The corresponding fallback takes 80
rather than 81 products. With three selected spots, accepting after growth
takes 19 rather than 21 products. Final estimated diagonal values agree
with the uncached algorithm.

A focused local build with **GNU Fortran 14.2.0** (`-std=f2018`,
`-fcheck=all`, `-Werror=implicit-interface`,
`-ffpe-trap=invalid,zero,overflow`) passed **45/45 root/Krylov tests** at
both `-O0 -g` and `-O2`. The non-Krylov MINPACK calls were compiled against
a disposable test-only interface facade, which is not included in this
archive and is not a substitute for the real dependency. The dependency-backed
`fpm test` was not run because `fpm` and its external dependencies are not
available in this environment. Testing with the actual FPM dependencies on
the target platform remains necessary.

### Krylov coverage-aware progressive verification (2026-10-12)

A rotated verification offset can be zero after integer rounding when the
number of exact-coordinate checks is greater than half the system size. In
that situation the earlier progressive implementation could check the same
set of sites again after probe expansion, missing a biased unexamined row.
The estimator now uses a path-compressed, per-estimate successor table to
choose the next unvisited coordinate (wrapping cyclically). At each stage
previously checked coordinates are still compared with the expanded probe
estimate without new Jacobian-vector evaluations. The fixed-probe path and
public `root` / `newton_krylov` interfaces are untouched.

`test/test_optimize_root_krylov_progressive_diagonal.f90` includes a
64-variable two-coupling example with seed 17 and 50 verification samples.
At eight probes, an exactly checked coupling in row 1 disagrees with the
estimate; at sixteen probes that row becomes acceptable. A second off-diagonal
coupling in row 5 produces a persistent 30% positive diagonal-estimate bias
at four, eight, and sixteen probes. The previous stage-2 verifier missed
row 5 because the stratum offset rounded to zero. The new coverage-aware
sampling checks row 5, rejects the biased estimate, and recovers the exact
diagonal in **16 + 64 = 80 analytic Jv calls**, including all cached checks.
The same test confirms **81 residual evaluations** for forward differences
and **161** for centered differences. Recompiling the original code against
the expanded regression reproduces both expected failures; the updated code
passes. This demonstrates detection for this constructed example, not a
statistical guarantee for every unexamined coordinate.

All 48 selected root/Krylov/linear-operator regression executables passed with
GNU Fortran 14.2.0 in checked debug (`-O0 -g -fcheck=all`, IEEE traps) and
optimized (`-O2`) builds. These focused builds used **temporary test-only
MINPACK and optimize facade modules**, not the real external dependencies;
`fpm test` could not run because FPM is absent in the environment. The
48-row dense-coupling progressive benchmark data are byte-for-byte identical
to the prior validated CSV, so the change has no effect on those cases.

### Krylov progressive rotated spot verification (2026-10-12)

`src/optimize/nonlin.f90` now selects previously unexamined exact diagonal
coordinates at successive verified Rademacher probe stages using deterministic
base-2 bit-reversed offsets within evenly spaced strata. Earlier exact checks
are revalidated against the newer estimate at zero additional Jacobian-vector
product cost. All coordinate spot products are cached for possible exact
fallback. The public API, default estimator, and fixed-budget behavior are
unchanged. For `s` stages with batch consistency, verification requires at
most `min(n,s*k)` distinct coordinate products; when exact fallback occurs,
the complete estimate still costs at most `q+n` products, including spots.

The extended `test/test_optimize_root_krylov_progressive_diagonal.f90` adds a
64-variable counterexample with one off-diagonal coupling in row 33. With
seed 17, the first eight probes fail batch consistency while the eight- and
sixteen-probe means at row 33 share a **30% bias**. The old method's repeated
spot at row 1 would miss this; the rotated spot at row 33 detects it at stage
2 and triggers exact fallback. The full 16-plus-64 analytic Jv budget and
81 forward-difference residual evaluations confirm that no product is
computed twice. Allowing 32 probes instead causes a third-stage acceptance
at **34 analytic Jv calls**, including a fresh row-17 spot check and cached
revalidation of row 33. This accepted result still has a 15% bias at row 33,
so rotated verification does not eliminate residual estimator error.

The tests additionally exercise: acceptance after a rejected original spot,
multiple cached coordinates, no repeated checks on exact fallback, full
coverage with `k=n`, callback failure at a shifted spot, and finite-difference
2-point and 3-point work counts. The prior 48-case dense rank-one benchmark
CSV was regenerated with GNU Fortran 14.2.0 at `-O2` and found to be
byte-for-byte unchanged. All **45/45** focused root/Krylov test programs
passed at both `-std=f2018 -O0 -g` and `-std=f2018 -O2`, using
`-ffree-line-length-132 -fcheck=all -Werror=implicit-interface` and
`-ffpe-trap=invalid,zero,overflow`. The temporary MINPACK interface-only
shim used to compile `root.f90` is outside the distribution archive; tests
of actual MINPACK implementations and full dependency-backed `fpm test`
remain outstanding because `fpm` and its Git dependencies are unavailable.

### Checked mean-field rank-one Krylov preconditioning (2026-10-12)

Added `inner_m_diag_rank_one=.true.` to `root(method='krylov')` and
`newton_krylov` as a SciFort-only, opt-in inverse-preconditioner extension.
The estimator requires the exact auto-diagonal path, computes `J*ones`, and
checks two sign-vector products before using a Sherman-Morrison inverse for
`J = diag(a) + u*ones^T`. The fit costs three Jv calls beyond the n exact
coordinate products; the accepted inverse action costs O(n) per application.
Inaccurate structure, a singular fitted diagonal, or a nearly singular
Sherman-Morrison denominator disables the correction and retains the original
safeguarded diagonal. A callback failure returns an unsuccessful solve.
This is a heuristic structural check and does not certify arbitrary dense
Jacobian couplings. The option is rejected for built-in MINRES and incompatible
probe/auto settings. Omitting it preserves all earlier behavior.

`test/test_optimize_root_krylov_rank_one.f90` checks both public entry points,
analytic fit/product counts, an exact mean-field inverse action, one-vector
GMRES convergence, default-path preservation, intentionally nonuniform
coupling, singular fit and update-denominator fallback, an ill-scaled matrix
whose ordinary inverse diagonal would be floored, forward/centered
finite-difference work counts, scalar dimensionality, validation before any
residual evaluation, and early termination on failed analytic Jv callbacks.

GNU Fortran 14.2.0 (`gfortran -std=f2018`) compiled all relevant modules.
The 47 `test_optimize_root_krylov*` / `test_optimize_krylov*` programs passed
with both configurations:

- Checked debug: `-O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow`
  and `-Werror=implicit-interface`.
- Optimized: `-O2 -Werror=implicit-interface`.

These focused builds used a **temporary MINPACK stub and optimize facade
module** for linking the independent Krylov tests, not the real FPM
external dependencies. The stubs and build artifacts are **not included** in
the ZIP. The full `fpm test` suite was not run because FPM and its dependency
checkout are not available here. Windows/FPM integration and numerical
performance benchmarks on more general dense Jacobians remain to be tested.

### Sketched rank-two Krylov preconditioning (2026-10-13)

Added the opt-in `inner_m_diag_rank_two_sketch` for known right factors
`V(n,2)`. Three bounded sign products estimate the true Jacobian diagonal
from a rowwise `3x3` solve; the existing rank-two fitting and two
independent structural probes then check the model. Accepted fits need
seven rather than `n+4` Jv calls per refresh; rejected fits use the
previous exact-diagonal and checked-rank-two path.

`test/test_optimize_root_krylov_rank_two_sketch.f90` verifies analytic
inverse action and Jv counts, public `root` and `newton_krylov` forwarding,
two-point and three-point finite differences with `nfev` accounting,
incompatible and dependent-column exact fallback, early Jacobian callback
failure, and invalid option handling. The adjacent benchmark
`benchmarks/benchmark_krylov_woodbury.f90` supplies 40 comparisons over
`n=32,128`, scaling exponents 2 and 6, and rank-one/rank-two systems.

GNU Fortran 14.2.0 (Debian, `-std=f2018`) checks used
`-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace`
and `-O2`, with `-Werror=implicit-interface`. The focused offline build
uses an explicitly labeled temporary MINPACK stub and `scifort_optimize`
forwarder for the otherwise independent Krylov programs; these test-only
files are not distributed. The full dependency-backed `fpm test` could
not be run because FPM and the dependency sources are unavailable here.

## Symmetric low-product rank-two Krylov discovery (2026-10-13)

The opt-in `inner_m_diag_rank_two_symmetric_base` extension was built on
Linux with GNU Fortran 14.2.0 and tested using two isolated configurations:

- checked debug: `-std=f2018 -O0 -g -fcheck=all
  -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface`;
- optimized: `-std=f2018 -O2 -Werror=implicit-interface`.

All **56 focused Krylov regression programs** passed in each configuration.
This selection includes the new small- and 128-variable symmetric rank-two
tests as well as the preceding root/Krylov regression programs. The new
analytic-Jacobian fit needs four Jv products on the accepted test models,
while rejected models invoke the existing n-coordinate exact-diagonal
fallback. Forward and centered finite differences, one-vector GMRES,
`newton_krylov` forwarding, six decades of base-diagonal scaling, invalid
options, nonfinite/zero base entries, rank-three contamination, and
callback failures are exercised.

The isolated `root` build used temporary MINPACK interfaces and a minimal
test-only umbrella-module facade; neither is included in this source
archive. Other non-Krylov MINPACK-dependent entry points were not
exercised. **The complete dependency-backed `fpm test` was not run:**
FPM and the full dependency set were unavailable in this environment.
The two fitted and two held-out directions remain heuristic structural
checks; they do not certify rank two for every arbitrary matrix.

## Adaptive-cost transpose checks for nonsymmetric rank two (2026-10-13)

Added `inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost=c` to
Krylov `root` and `newton_krylov`. For verified fits, the optional policy
uses one or two held-out transpose Jacobian products depending on the
modeled cost and fit quality, always retaining at least one check. The
existing default (no held-out adjoints) and fixed two-check option are
unchanged. The policy forces two checks on weak interpolation pivots
or borderline forward-model verification, regardless of declared cost.

`test/test_optimize_root_krylov_rank_two_transpose_adaptive.f90` covers
three- versus four-adjoint accepted setups, inverse action, a deliberately
borderline forward fit, inconsistent adjoint rejection, exact-diagonal
fallback, failing callbacks, `newton_krylov` forwarding and invalid
option combinations/costs. GNU Fortran 14.2: `-std=f2018 -O0 -g
-fcheck=all -ffpe-trap=invalid,zero,overflow -Werror=implicit-interface`
and `-std=f2018 -O2 -Werror=implicit-interface` each passed all **60/60**
Krylov-focused programs with temporary MINPACK and umbrella-module
stubs. A full dependency-backed `fpm test` was not available.

The updated `benchmarks/benchmark_krylov_rank_two_discovery.f90`
completed **56/56** cases with successful solves. In a well-conditioned
rank-two run using six decades of diagonal scaling, high-cost adaptive
verification (`c=16`) used 6 forward + 3 transpose products at n=32 and
6 forward + 4 transpose products at n=128 (counts include the linear
solve). For comparison, fixed two-check verification used 6+4 in both.
CPU times and cost-model columns are recorded in
`benchmarks/krylov_rank_two_discovery_gfortran14_o2.csv`; cost weights
are models, not timings. The held-out probes remain heuristic tests,
not exhaustive low-rank certificates.

## Optional transpose verification and discovery benchmark (2026-10-13)

Added `test/test_optimize_root_krylov_rank_two_transpose_verify.f90`
for the new opt-in transpose-consistency checks. It constructs a false
analytic transpose callback that agrees **exactly** on the two fitting
vectors yet differs on a held-out vector. The old fit accepts this
callback; optional verification rejects it and uses the original exact
Jacobian diagonal. Also checks successful four-Jv/four-transpose-product
setup, `newton_krylov` keyword forwarding, failed held-out callbacks,
and early invalid-option rejection.

GNU Fortran 14.2 checked debug (`-std=f2018 -O0 -g -fcheck=all
-ffpe-trap=invalid,zero,overflow -Werror=implicit-interface`) and
optimized (`-std=f2018 -O2 -Werror=implicit-interface`) runs each passed
**59 of 59 focused Krylov test programs** with temporary MINPACK and
umbrella-module stubs. The complete dependency-backed `fpm test` was
not run; the stubs are not distributed.

`benchmarks/benchmark_krylov_rank_two_discovery.f90` was compiled with
the optimized solver objects and run for 40 combinations of n=32/128,
rank-two/rank-three coupling, log10 diagonal scale=2/6, and five policies.
All 40 reported success, maximum scaled residual < 1.5e-13 and maximum
solution error < 1.5e-13 in this deterministic linear workload. The raw
CSV `benchmarks/krylov_rank_two_discovery_gfortran14_o2.csv` includes
CPU times and independent work counts. At n=128 with rank two, dense
and streamed methods each needed 130 forward products (including the
linear solve); basic low-product needed 6 Jv + 2 J^Tv, and transpose-
verified low-product needed 6 Jv + 4 J^Tv. If one J^Tv costs 16 Jv,
modeled costs are 130, 38 and 70 respectively. At n=32 the corresponding
costs were 34, 38 and 70, so the low-product path is not always cheaper.
For rank-three contamination all discovery policies reverted to diagonal
scaling in these test cases; low-product rejection incurred additional
setup work. The relative-cost columns are **models**, not measured times.

## Nonsymmetric low-product rank-two Krylov discovery (2026-10-13)

The opt-in `inner_m_diag_rank_two_nonsymmetric_base` plus analytic
`jacobian_transpose_product` implementation was tested on Linux with
GNU Fortran 14.2.0 using `-std=f2018 -Werror=implicit-interface` in
both checked debug (`-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow`)
and optimized (`-O2`) configurations.

All **58 focused Krylov test programs passed** in each configuration.
The selection includes the new 24-variable and
128-variable nonsymmetric rank-two programs. Accepted analytic fits use
four forward Jacobian products and two transpose products; failed
heuristic fits retain the existing n-coordinate exact-diagonal fallback.
The new tests check exact inverse action, single-basis GMRES convergence,
`newton_krylov` keyword forwarding, forward and centered finite differences,
residual evaluation counts, six decades of diagonal scaling, corrupted
base data, rank-three contamination, inconsistent transpose callbacks,
near-singular Woodbury fallback, option validation and failed callbacks.

The isolated build used temporary MINPACK and umbrella-module stubs;
these are not distributed. Full dependency-backed `fpm test` was not
run because FPM and its external dependencies were unavailable.

## Nonlinear nonsymmetric rank-two refresh (2026-10-13)

The new optional `inner_m_diag_rank_two_nonsymmetric_base_fun` interface
provides the uncoupled `D(x)` at each scheduled Krylov preconditioner refresh.
`test/test_optimize_root_krylov_rank_two_nonlinear_base.f90` constructs a
24-dimensional genuinely nonlinear model with a cubic diagonal term and two
varying nonsymmetric coupling coefficients. The exact target solution is
specified independently of the initial iterate. It tests Newton convergence,
residual and solution errors, multiple changing base-diagonal evaluations,
`inner_m_diag_update_interval=0` (one initial evaluation),
`newton_krylov` forwarding, an incorrect diagonal at a later step causing
exact-diagonal fallback, callback failure at a later refresh, nonfinite
callback output, and incompatible or missing option combinations.

With GNU Fortran 14.2.0 (`gfortran` on Linux) **61 of 61** focused
`test_optimize_root_krylov*.f90` and `test_optimize_krylov*.f90` programs
passed in both configurations:

```text
-std=f2018 -Werror=implicit-interface -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow
-std=f2018 -Werror=implicit-interface -O2
```

The isolated builds link test-only MINPACK and umbrella-module stubs from
the prior focused harness. These are not distributed. This is **not** a
complete dependency-backed `fpm test` run; FPM and the full declared
external dependencies were not available in this environment.


## Nonlinear Krylov refresh-interval benchmark (2026-10-13)

New original Fortran benchmark:
`benchmarks/benchmark_krylov_nonlinear_refresh.f90`. It constructs
nonlinear systems with known targets,
`F(x) = D0*(x + 0.5*r*x**3) + u1*sin(r*v1^T*x)/r +
u2*(v2^T*x+0.25*r*(v2^T*x)**3) - F(target)`, with exact
Jacobian and transpose-Jacobian products, and exact uncoupled diagonal
`D(x) = D0*(1+1.5*r*x**2)`. The benchmark independently checks the
forward derivative with centered differences and the transpose identity.

A full factorial of n=32/128, log10 diagonal scale=2/5, nonlinear
rate=0.15/1.5/6, and refresh interval=0/1/2/4/8 gives **60 cases**.
All 60 converged both in checked debug and optimized GNU Fortran 14.2
runs. All had scaled residual and solution error below `7.5e-11` in the
optimized run. The CSV at
`benchmarks/krylov_nonlinear_refresh_gfortran14_o2.csv` reports
iteration count, residual evaluations, diagonal callback count, forward
and transpose Jacobian products, CPU seconds, and modeled work at
transpose-product costs 1x/4x/16x. The work figures are counts weighted
by assumed relative cost; they are **not** measured timing equivalents.
CPU timings are short and noisy; compare product counts before relying
on the timings.

Averaged across the 12 workloads for each refresh interval:

| Interval | Mean Jv | Mean J^Tv | Work at 1x | Work at 4x | Work at 16x |
|---|---:|---:|---:|---:|---:|
| 0 | 53.58 | 3.50 | 57.08 | 67.58 | 109.58 |
| 1 | 48.42 | 17.00 | 65.42 | 116.42 | 320.42 |
| 2 | 37.67 | 9.08 | 46.75 | 74.00 | 183.00 |
| 4 | 41.08 | 5.58 | 46.67 | 63.42 | 130.42 |
| 8 | 53.58 | 3.50 | 57.08 | 67.58 | 109.58 |

At n=128, log10 scale=5 and rate=6, interval 4 needed 34 forward and
8 transpose products, versus 66+4 for a frozen preconditioner and
36+24 for interval 1. At 1x adjoint cost, modeled work was 42/70/60
respectively; at 16x, it was 162/130/420. All converged in 7 Newton
iterations. Thus a shorter interval sometimes reduces linear work but
increases total adjoint-weighted work. Since no benchmark run required
more than 7 Newton iterations, interval 8 and 0 coincided in this set.
These results are specific to the tested models and tolerances.

The existing nonlinear callback regression now checks intervals 2/3/4,
callback count `1+max(0,nit-2)/interval`, distinct accepted iterates,
and a deferred callback failure at interval 3. All **61/61** focused
Krylov programs passed in both configurations:

```text
gfortran -std=f2018 -Werror=implicit-interface -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow
gfortran -std=f2018 -Werror=implicit-interface -O2
```

The focused harness used compiled solver objects from the matching
preceding source snapshot and temporary MINPACK and umbrella-module
stubs, not a complete rebuilt dependency graph. The benchmark executable
was compiled and run under both sets of flags. `fpm test` was not run
because FPM and the declared external dependencies were unavailable.


### Work-aware adaptive Krylov refresh (October 2026)

Added `inner_m_diag_update_work_growth` and
`inner_m_diag_update_rebuild_cost` to `root` and `newton_krylov`. The
work signal counts preconditioner applications during each inner solve;
its baseline is measured immediately after a refresh. A later accepted
step requests an early refresh only when the increase in work exceeds
the selected growth factor plus the estimated rebuild cost. Strong
nonlinear progress does not suppress this independent work signal;
existing stall and maximum-age policies remain in force.

The new deterministic regression verifies early refresh on a single
costly inner solve, no refresh when rebuild cost is too high, recovery
after a refresh, zero baseline work, and early validation of invalid
parameters. All **63/63** files matching `test/*krylov*.f90` passed in
both GNU Fortran 14.2 configurations:

```text
checked:   -std=f2018 -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow
optimized: -std=f2018 -O2
```

The test harness compiled the solver and helper modules from the
current sources and used temporary MINPACK and umbrella forwarding
stubs. This is **not** a full dependency-backed `fpm test` (FPM and
external dependencies were unavailable). No real-world speedup is
claimed from this deterministic scheduling test.

### Adaptive nonlinear diagonal refresh (October 2026)

Added `inner_m_diag_update_adaptive` and optional finite threshold
`inner_m_diag_update_stall_ratio` to `root` and `newton_krylov`. The
preconditioner refreshes no later than its positive maximum interval,
and it may refresh earlier when the **accepted** residual's Euclidean norm
falls by less than the configured amount. The default scheduling is
unchanged. The test `test_optimize_root_krylov_adaptive_refresh.f90` uses
a deterministic damped inner solver to test known schedules, a single
weak step followed by strong steps, the maximum-age restart, forwarding,
and early rejection of invalid parameters.

A separate 84-case optimized benchmark CSV is saved as
`benchmarks/krylov_nonlinear_refresh_adaptive_gfortran14_o2.csv`.
It extends the preceding 60 fixed-interval cases with adaptive maximum
interval 8 at stall thresholds 0.8 and 0.95, retaining n=32/128,
log10 diagonal scales=2/5 and nonlinear rates=0.15/1.5/6.
All 84 cases converged under GNU Fortran 14.2 with scaled residuals below
`7.5e-11`. The new modes did not trigger earlier refreshes on these
well-progressing problems, so their modeled work was identical to a
fixed interval of 8 (which behaves like a frozen preconditioner for at
most 7 Newton iterations). **This benchmark does not demonstrate a speed
improvement**; the adaptive option is useful when progress actually
stalls, as confirmed by the deterministic schedule regression. Actual
preconditioner setup cost depends on the selected inverse estimator and
transpose-product expense; CPU timing for these small problems is noisy.

The 62 files matching `test/*krylov*.f90` compiled and passed in both:

```text
gfortran -std=f2018 -O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow
gfortran -std=f2018 -O2
```

The matching solver and helper modules were compiled from the updated
sources. Temporary MINPACK/umbrella forwarding stubs were used in the
isolated test harness; no stub files are included in the archive. Full
FPM dependency integration and `fpm test` were not available.

### Work-triggered refresh on changing nonlinear Jacobians (October 2026)

The new `benchmarks/benchmark_krylov_work_refresh.f90` exercises
`root(method='krylov')` with a *built-in* GMRES inner solver on genuinely
nonlinear nonsymmetric rank-two systems. Each Jacobian has the form
`J(x) = D(x) + U(x)*V^T`, with nonconstant `D(x)`. This deliberately
produces more challenging inner solves than the earlier fast-converging
84-case benchmark. The benchmark uses dimensions 32 and 128, diagonal
scale ranges of 2 and 5 decimal orders, and nonlinear rates 0.3, 0.9,
and 2.0. Each of the 12 systems is solved with nine policies (108 cases):
fixed intervals 0, 1, 2, 4, and 8; residual-stall adaptive refresh;
and work-triggered adaptive refresh using growth factor 1.2 with rebuild
costs 0 and 4, plus growth factor 2.0 with cost 0. The work-triggered
policies have maximum age 12. All use a known nonlinear base-diagonal
callback, analytic Jacobian and transpose products, and independent
residual/solution checks.

The reproducible optimized results are recorded in
`benchmarks/krylov_work_refresh_gfortran14_o2.csv`. Representative
*mean counts per 128-variable case* from GNU Fortran 14.2 `-O2` are:

| Policy | Base refreshes | Jv | J^Tv | Work at 1x | Work at 16x |
| --- | ---: | ---: | ---: | ---: | ---: |
| Initial only (`fixed0`) | 1.00 | 71.83 | 4.00 | 75.83 | 135.83 |
| Every iteration (`fixed1`) | 4.67 | 30.83 | 18.67 | 49.50 | 329.50 |
| Fixed interval 2 | 2.67 | 28.33 | 10.67 | 39.00 | 199.00 |
| Work-aware, growth 1.2, rebuild cost 0 | 2.67 | 28.33 | 10.67 | 39.00 | 199.00 |
| Work-aware, growth 1.2, rebuild cost 4 | 2.00 | 36.00 | 8.00 | 44.00 | 164.00 |

Here `Jv` and `J^Tv` are callback counts for the *complete solve*,
not only the preconditioner setup. Modeled work at 1x or 16x equals
`Jv + cost*J^Tv`; it excludes ordinary residual evaluations, setup
arithmetic, and other costs. Actual CPU times are included in the CSV,
but are noisy for these small cases. In all 108 optimized and all 108
checked-debug benchmark cases, the solver reported success, the scaled
residual remained below `1e-9`, and the recovered root was within
`1e-7` of the prescribed solution. All 64 Krylov-focused tests passed
in each configuration; the new nonlinear regression confirms early work
refresh, suppression by higher rebuild cost, and correct convergence.

The focused validation commands used freshly compiled SciFort solver
sources, GNU Fortran 14.2, and temporary test-only MINPACK and umbrella
module stubs:

```text
gfortran -std=f2018 -O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow
gfortran -std=f2018 -O2
```

The builds also used `-ffree-line-length-none` and
`-Werror=implicit-interface`. Full `fpm test` remains unverified because
FPM and the declared external dependencies were unavailable. No
workload-wide improvement or single best refresh interval is claimed.

## 2026-10-31 — tensor PCHIP and banded B-spline collocation

The original per-axis dense LU in `src/interpolate_regular_grid.f90` was
replaced by an original narrow-band factorization (half-bandwidth 1, 3, or 5).
SciFort now also implements real-valued 1–6D `RegularGridInterpolator` PCHIP.

Independently compared with locally installed SciPy 1.17.0, NumPy, and
GNU Fortran 14.2 in checked (`-O0 -g -fcheck=all -fbacktrace`) and optimized
(`-O3`) configurations:

- Tensor PCHIP: **232 SciPy cases, 3,528 scalar values**, on *each* of two
  independent random seeds (20261101 and 8675309), in both builds. Largest
  scale-normalized discrepancy across these cases: approximately `4e-14`.
- Existing 1–4D spline grid comparison: **189 cases, 2,997 comparisons**
  per seed and build; two seeds passed.
- Existing 1–6D linear/nearest comparison: **540 cases, 17,340 values**
  per build (default seed); both builds passed.
- Existing explicit 1D B-spline suite: **1,390 cases, 1,077,003 comparisons**
  per build; both passed.
- Existing 1D PCHIP/Akima/Makima/Hermite suite: **480 cases, 149,400 values**
  per build; both passed.
- All six native interpolation programs (`regular_grid`, `tensor_grid`,
  `pchip_grid`, `bspline`, `cubic`, and `cubic_spline`) passed in both builds.
- One run of the quintic 1D fit benchmark at 512 knots under GNU Fortran
  `-O3`: old dense solver `0.769019 s`, new banded solver `0.000104 s`;
  maximum nodal interpolation error `4.44e-16`. Treat these as machine- and
  build-specific measurements, not universal speedup claims.

The long randomized `compare_cubic_spline_scipy.py` run did not complete
within the execution limit, so it is **not** reported as passed here.
The native CubicSpline test and unaffected earlier PCHIP/Akima and B-spline
reference comparisons did pass. FPM was not installed locally; the
repository-wide `fpm test` suite was not executed.

Reproduction:

```sh
python tools/interpolate_reference/compare_tensor_pchip_scipy.py
python tools/interpolate_reference/compare_tensor_pchip_scipy.py --seed 8675309
python tools/interpolate_reference/compare_tensor_grid_scipy.py
python tools/interpolate_reference/compare_tensor_grid_scipy.py --seed 8675309
```
