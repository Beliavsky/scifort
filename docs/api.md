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

## Elemental operation

All current routines are scalar elemental functions. For example:

```fortran
real(dp) :: x(5)
real(dp) :: p(5)

p = normal_cdf(x, loc=1.0_dp, scale=2.0_dp)
```

Optional dummy arguments are scalar in an elemental invocation. Conformable
array arguments may be supplied for arguments that are present.
