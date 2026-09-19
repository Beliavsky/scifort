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

Accuracy: in validation against mpmath at 50 digits for `1e-8 <= a <= 1e6`,
the error of `P` and `Q` stayed below
`1e-15 * (1 + |log V| + x f(x) / V)`, where `V` is the computed tail and
`x f(x) / V` is the condition number with respect to `x`. The cost of an
evaluation near `x = a` grows like `sqrt(a)`. A uniform asymptotic expansion
for very large `a` is not yet implemented.

## Elemental operation

All current routines are scalar elemental functions. For example:

```fortran
real(dp) :: x(5)
real(dp) :: p(5)

p = normal_cdf(x, loc=1.0_dp, scale=2.0_dp)
```

Optional dummy arguments are scalar in an elemental invocation. Conformable
array arguments may be supplied for arguments that are present.
