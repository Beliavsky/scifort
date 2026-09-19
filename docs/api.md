# Initial probability-distribution API

Every current distribution uses `loc` and `scale` parameterization. Optional
arguments default to `loc = 0` and `scale = 1`.

```fortran
use scifort_stats
```

## Normal

```text
normal_pdf(x [, loc, scale])
normal_logpdf(x [, loc, scale])
normal_cdf(x [, loc, scale])
normal_sf(x [, loc, scale])
normal_logcdf(x [, loc, scale])
normal_logsf(x [, loc, scale])
normal_ppf(p [, loc, scale])
normal_isf(p [, loc, scale])
```

## Uniform

The support is `[loc, loc + scale]`.

```text
uniform_pdf, uniform_logpdf, uniform_cdf, uniform_sf
uniform_logcdf, uniform_logsf, uniform_ppf, uniform_isf
```

## Exponential

The support is `[loc, infinity)` and `scale` is the inverse rate.

```text
exponential_pdf, exponential_logpdf, exponential_cdf, exponential_sf
exponential_logcdf, exponential_logsf, exponential_ppf, exponential_isf
```

## Laplace, logistic, and Cauchy

Each uses `loc` as its center and positive `scale` as its scale parameter.
They follow the same eight-function naming pattern.

## Gamma

Shape `a > 0` is a required argument after `x` or `p`. With
`z = (x - loc) / scale`, the density is
`z**(a - 1) exp(-z) / (Gamma(a) scale)` on `z >= 0`, matching
`scipy.stats.gamma`.

```text
gamma_pdf(x, a [, loc, scale])
gamma_logpdf, gamma_cdf, gamma_sf, gamma_logcdf, gamma_logsf
gamma_ppf(p, a [, loc, scale]), gamma_isf(p, a [, loc, scale])
```

At `x = loc` the density is `+infinity` for `a < 1`, `1 / scale` for
`a = 1`, and `0` for `a > 1`. A non-finite or nonpositive shape returns NaN.

## Chi-square

`chi2_*(x, df [, loc, scale])` is the gamma distribution with shape `df / 2`
and scale `2 * scale`, matching `scipy.stats.chi2`. Degrees of freedom must be
finite and positive; noninteger values are allowed.

## Beta

Shapes `a > 0` and `b > 0` follow `x` or `p`. With `z = (x - loc) / scale`,
the density is `z**(a - 1) (1 - z)**(b - 1) / (B(a, b) scale)` on
`0 <= z <= 1`, matching `scipy.stats.beta`.

```text
beta_pdf(x, a, b [, loc, scale])
beta_logpdf, beta_cdf, beta_sf, beta_logcdf, beta_logsf
beta_ppf(p, a, b [, loc, scale]), beta_isf(p, a, b [, loc, scale])
```

At `z = 0` the density is `+infinity` for `a < 1`, `b / scale` for `a = 1`,
and `0` for `a > 1`; the endpoint `z = 1` is symmetric in `b`.

## Student t

```text
t_pdf(x, df [, loc, scale])
t_logpdf, t_cdf, t_sf, t_logcdf, t_logsf
t_ppf(p, df [, loc, scale]), t_isf(p, df [, loc, scale])
```

Degrees of freedom must be positive; noninteger values are allowed and
`df = +infinity` gives the normal distribution, as in `scipy.stats.t`. The
tails use the regularized incomplete beta function with both `x` and `1 - x`
formed directly. For `|t| > 1e150 sqrt(df)`, and for quantiles that far out,
the exact leading term of the tail is used in scaled form, so tail
probabilities and quantiles remain accurate where SciPy 1.15.3 returns `0` for
`t.sf(1e200, 1)` and `-1e100` for `t.ppf(1e-300, 1)`.

## F

```text
f_pdf(x, dfn, dfd [, loc, scale])
f_logpdf, f_cdf, f_sf, f_logcdf, f_logsf
f_ppf(p, dfn, dfd [, loc, scale]), f_isf(p, dfn, dfd [, loc, scale])
```

Both degrees of freedom must be finite and positive. SciPy accepts
`dfd = +infinity` but SciPy 1.15.3 returns values inconsistent with the
limiting distribution (for example `f.cdf(1, 2, inf) = 0`), so SciFort returns
NaN for infinite degrees of freedom. The density at `x = loc` is `+infinity`
for `dfn < 2`, `1 / scale` for `dfn = 2`, and `0` for `dfn > 2`.

## Lognormal

Shape `s > 0` is the standard deviation of `log((X - loc) / scale)`. The
support is `x >= loc`. With `z = (x - loc) / scale > 0`, `log(z) / s` is
standard normal, matching `scipy.stats.lognorm`.

```text
lognormal_pdf(x, s [, loc, scale])
lognormal_logpdf, lognormal_cdf, lognormal_sf, lognormal_logcdf, lognormal_logsf
lognormal_ppf(p, s [, loc, scale]), lognormal_isf(p, s [, loc, scale])
```

## Weibull

Shape `c > 0` follows `x` or `p`. With `z = (x - loc) / scale >= 0`, the density is
`(c / scale) z**(c - 1) exp(-z**c)`, matching `scipy.stats.weibull_min`.

```text
weibull_pdf(x, c [, loc, scale])
weibull_logpdf, weibull_cdf, weibull_sf, weibull_logcdf, weibull_logsf
weibull_ppf(p, c [, loc, scale]), weibull_isf(p, c [, loc, scale])
```

At `x = loc`, the density is `+infinity` for `c < 1`, `1 / scale` for `c = 1`,
and `0` for `c > 1`.

## Pareto

Shape `b > 0` follows `x` or `p`. The support is `x >= loc + scale`, so
`z = (x - loc) / scale >= 1`. The density is `(b / scale) / z**(b + 1)`,
matching `scipy.stats.pareto`.

```text
pareto_pdf(x, b [, loc, scale])
pareto_logpdf, pareto_cdf, pareto_sf, pareto_logcdf, pareto_logsf
pareto_ppf(p, b [, loc, scale]), pareto_isf(p, b [, loc, scale])
```

## Rayleigh

The support is `x >= loc`. With `z = (x - loc) / scale >= 0`, the density is
`(z / scale) exp(-z**2 / 2)`, matching `scipy.stats.rayleigh`.

```text
rayleigh_pdf(x [, loc, scale])
rayleigh_logpdf, rayleigh_cdf, rayleigh_sf, rayleigh_logcdf, rayleigh_logsf
rayleigh_ppf(p [, loc, scale]), rayleigh_isf(p [, loc, scale])
```

## Discrete distributions

Discrete families provide `pmf` and `logpmf` in place of `pdf` and `logpdf`.
The count may be an `integer` or a `real(dp)`; the two forms are separate
specific procedures behind one generic name, so the count and the number of
trials must have the same type. Following `scipy.stats`:

- a non-integer count has probability zero, and `cdf`, `sf`, and their
  logarithms use its floor;
- `ppf(p)` is the smallest count whose `cdf` reaches `p`, and `isf(p)` the
  smallest count whose `sf` falls to `p` or below;
- `ppf(0)` returns `loc - 1`, the value just below the support, and `isf(1)`
  does the same. `ppf(1)` returns `+infinity` for the Poisson distribution
  and `loc + n` for the binomial one.

### Poisson

```text
poisson_pmf(k, mu [, loc])
poisson_logpmf, poisson_cdf, poisson_sf, poisson_logcdf, poisson_logsf
poisson_ppf(p, mu [, loc]), poisson_isf(p, mu [, loc])
```

The mean must satisfy `mu >= 0`; `mu = 0` places all mass at `loc` and
`mu = +infinity` pushes it beyond every finite count. The tails use
`P(X <= k) = Q(k + 1, mu)` and `P(X > k) = P(k + 1, mu)` (DLMF 8.4.8 with
8.4.11), so no sum over counts is formed.

### Binomial

```text
binomial_pmf(k, n, p [, loc])
binomial_logpmf, binomial_cdf, binomial_sf, binomial_logcdf, binomial_logsf
binomial_ppf(probability, n, p [, loc]), binomial_isf(probability, n, p [, loc])
```

The number of trials `n` must be a finite nonnegative integer value and the
success probability must lie in `[0, 1]`; other values give NaN, as in
`scipy.stats.binom`. The tails use `P(X > k) = I_p(k + 1, n - k)`
(DLMF 8.17.5).

## Special functions

```fortran
use scifort_special
```

```text
gammainc(a, x)       regularized lower incomplete gamma P(a, x)
gammaincc(a, x)      regularized upper incomplete gamma Q(a, x) = 1 - P(a, x)
gammaincinv(a, p)    x such that P(a, x) = p
gammainccinv(a, q)   x such that Q(a, x) = q
```

All four are elemental. `gammaincc` is computed directly rather than as
`1 - gammainc`. Special values follow `scipy.special`:

- `a < 0`, `x < 0`, or a NaN argument returns NaN;
- `a = 0` gives `P = 1` for `x > 0` and NaN for `x = 0`;
- `a = +infinity` gives `P = 0` for finite `x` and NaN for `x = +infinity`;
- the inverses return NaN for `a <= 0`, `a = +infinity`, or `p` outside
  `[0, 1]`, and map the probability endpoints to `0` and `+infinity`.

An inverse whose exact value is below the smallest positive normal number
returns `0`.

```text
betainc(a, b, x)       regularized incomplete beta I_x(a, b)
betaincc(a, b, x)      1 - I_x(a, b), computed directly
betaincinv(a, b, p)    x such that I_x(a, b) = p
betainccinv(a, b, q)   x such that 1 - I_x(a, b) = q
```

As in `scipy.special`, `a` and `b` must be finite and positive and `x`, `p`,
and `q` must lie in `[0, 1]`; otherwise the result is NaN. The endpoints map
to `0` and `1`. Roots above one half are computed through the complement
`1 - x`, so both `x` and `1 - x` keep their relative accuracy internally.

Accuracy: in validation against mpmath at 50 digits for `1e-8 <= a <= 1e6`,
the error of `P` and `Q` stayed below
`1e-15 * (1 + |log V| + x f(x) / V)`, where `V` is the computed tail and
`x f(x) / V` is the condition number with respect to `x`. The cost of an
evaluation near `x = a` grows like `sqrt(a)`. A uniform asymptotic expansion
for very large `a` is not yet implemented.

For the incomplete beta function the same bound, with `x f(x) / V` replaced by
`min(x, 1 - x) f(x) / V`, held in validation for parameters from `1e-7` to
`1e5`; see `docs/validation.md`. When one parameter is much larger than the
other, an expansion in incomplete gamma functions replaces the continued
fraction near the transition region, where the continued fraction would
lose accuracy in proportion to the large parameter.

The discrete distributions inherit these bounds, with one addition: a
probability is also sensitive to the rounding of `mu` or `p`, which moves it
by `mu pmf(k)` or `n p pmf(k; n - 1, p)`. In a comparison of 4000 random
cases with SciPy 1.15.3 the median relative difference was at most `4e-15`,
and in every one of the six largest differences SciPy was the less accurate
of the two (see `docs/validation.md`).

## Elemental operation

All current routines are scalar elemental functions. For example:

```fortran
real(dp) :: x(5)
real(dp) :: p(5)

p = normal_cdf(x, loc=1.0_dp, scale=2.0_dp)
```

Optional dummy arguments are scalar in an elemental invocation. Conformable
array arguments may be supplied for arguments that are present.
