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
- The C ABI does not yet expose the new functions.

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
- Only GNU Fortran was tested. The C ABI does not yet expose these functions.

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
- Only GNU Fortran was tested. The C ABI does not expose these functions.
