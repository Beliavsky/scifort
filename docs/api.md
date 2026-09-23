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

## Arcsine

The standardized support is `[0, 1]` and the density is
`1 / (pi*sqrt(z*(1-z)))`, matching `scipy.stats.arcsine`. The density diverges
at both support endpoints while the CDF remains zero and one there.

```text
arcsine_pdf(x [, loc, scale])
arcsine_logpdf, arcsine_cdf, arcsine_sf, arcsine_logcdf, arcsine_logsf
arcsine_ppf(p [, loc, scale]), arcsine_isf(p [, loc, scale])
```

## Half-normal

The support is `x >= loc`. With `z = (x-loc)/scale`, the density is
`sqrt(2/pi)*exp(-z**2/2)/scale`, matching `scipy.stats.halfnorm`.

```text
halfnorm_pdf(x [, loc, scale])
halfnorm_logpdf, halfnorm_cdf, halfnorm_sf, halfnorm_logcdf, halfnorm_logsf
halfnorm_ppf(p [, loc, scale]), halfnorm_isf(p [, loc, scale])
```

## Half-Cauchy

The support is `x >= loc`. With `z = (x-loc)/scale`, the density is
`2/(pi*scale*(1+z**2))`, matching `scipy.stats.halfcauchy`.

```text
halfcauchy_pdf(x [, loc, scale])
halfcauchy_logpdf, halfcauchy_cdf, halfcauchy_sf
halfcauchy_logcdf, halfcauchy_logsf
halfcauchy_ppf(p [, loc, scale]), halfcauchy_isf(p [, loc, scale])
```

## Lomax (Pareto II)

The support is `x >= loc`. With shape `c > 0` and
`z = (x-loc)/scale`, the density is
`(c/scale)*(1+z)**(-c-1)`, matching `scipy.stats.lomax`.

```text
lomax_pdf(x, c [, loc, scale])
lomax_logpdf, lomax_cdf, lomax_sf, lomax_logcdf, lomax_logsf
lomax_ppf(p, c [, loc, scale]), lomax_isf(p, c [, loc, scale])
```

## Chi

Support is `x >= loc`. With `df > 0`, `z = (x-loc)/scale`, and
`a = df/2`, the standardized density is proportional to
`z**(df-1)*exp(-z**2/2)`, matching `scipy.stats.chi`.

```text
chi_pdf(x, df [, loc, scale])
chi_logpdf, chi_cdf, chi_sf, chi_logcdf, chi_logsf
chi_ppf(p, df [, loc, scale]), chi_isf(p, df [, loc, scale])
```

## Maxwell

The Maxwell family is the chi distribution with `df = 3` and supports
`x >= loc`, matching `scipy.stats.maxwell`.

```text
maxwell_pdf(x [, loc, scale])
maxwell_logpdf, maxwell_cdf, maxwell_sf, maxwell_logcdf, maxwell_logsf
maxwell_ppf(p [, loc, scale]), maxwell_isf(p [, loc, scale])
```

## Cosine

The cosine family has standardized support `[-pi, pi]`, shifted and scaled by
`loc` and positive `scale`, matching `scipy.stats.cosine`.

```text
cosine_pdf(x [, loc, scale])
cosine_logpdf, cosine_cdf, cosine_sf, cosine_logcdf, cosine_logsf
cosine_ppf(p [, loc, scale]), cosine_isf(p [, loc, scale])
```

## Semicircular

The semicircular family has standardized support `[-1, 1]`, shifted and scaled
by `loc` and positive `scale`, matching `scipy.stats.semicircular`.

```text
semicircular_pdf(x [, loc, scale])
semicircular_logpdf, semicircular_cdf, semicircular_sf
semicircular_logcdf, semicircular_logsf
semicircular_ppf(p [, loc, scale]), semicircular_isf(p [, loc, scale])
```

## Anglit

The Anglit family has standardized support `[-pi/4, pi/4]` and density
`cos(2*x)`, shifted and scaled by `loc` and positive `scale`, matching
`scipy.stats.anglit`.

```text
anglit_pdf(x [, loc, scale])
anglit_logpdf, anglit_cdf, anglit_sf, anglit_logcdf, anglit_logsf
anglit_ppf(p [, loc, scale]), anglit_isf(p [, loc, scale])
```

## Moyal

The Moyal family is supported on the real line with standardized density
`exp(-(x + exp(-x))/2) / sqrt(2*pi)`, shifted and scaled by `loc` and positive
`scale`, matching `scipy.stats.moyal`.

```text
moyal_pdf(x [, loc, scale])
moyal_logpdf, moyal_cdf, moyal_sf, moyal_logcdf, moyal_logsf
moyal_ppf(p [, loc, scale]), moyal_isf(p [, loc, scale])
```

## Hyperbolic secant

The hyperbolic-secant family is supported on the real line with standardized
density `sech(x)/pi`, shifted and scaled by `loc` and positive `scale`, matching
`scipy.stats.hypsecant`.

```text
hypsecant_pdf(x [, loc, scale])
hypsecant_logpdf, hypsecant_cdf, hypsecant_sf
hypsecant_logcdf, hypsecant_logsf
hypsecant_ppf(p [, loc, scale]), hypsecant_isf(p [, loc, scale])
```

## Half logistic

The half-logistic family has standardized support `[0, +infinity)` with
`loc` as its lower endpoint and positive `scale`, matching
`scipy.stats.halflogistic`.

```text
halflogistic_pdf(x [, loc, scale])
halflogistic_logpdf, halflogistic_cdf, halflogistic_sf
halflogistic_logcdf, halflogistic_logsf
halflogistic_ppf(p [, loc, scale]), halflogistic_isf(p [, loc, scale])
```

## Inverse gamma

The inverse-gamma family uses positive shape `a`, lower support endpoint `loc`,
and positive `scale`, matching `scipy.stats.invgamma`.

```text
invgamma_pdf(x, a [, loc, scale])
invgamma_logpdf, invgamma_cdf, invgamma_sf, invgamma_logcdf, invgamma_logsf
invgamma_ppf(p, a [, loc, scale]), invgamma_isf(p, a [, loc, scale])
```

## Inverse Gaussian

The inverse-Gaussian family uses SciPy's positive shape parameter `mu`, lower
support endpoint `loc`, and positive `scale`, matching `scipy.stats.invgauss`.

```text
invgauss_pdf(x, mu [, loc, scale])
invgauss_logpdf, invgauss_cdf, invgauss_sf, invgauss_logcdf, invgauss_logsf
invgauss_ppf(p, mu [, loc, scale]), invgauss_isf(p, mu [, loc, scale])
```

## Levy

The Levy family has lower support endpoint `loc` and positive `scale`, matching
`scipy.stats.levy`.

```text
levy_pdf(x [, loc, scale])
levy_logpdf, levy_cdf, levy_sf, levy_logcdf, levy_logsf
levy_ppf(p [, loc, scale]), levy_isf(p [, loc, scale])
```

## Log-Laplace

The log-Laplace family uses positive shape `c`, lower support endpoint `loc`,
and positive `scale`, matching `scipy.stats.loglaplace`.

```text
loglaplace_pdf(x, c [, loc, scale])
loglaplace_logpdf, loglaplace_cdf, loglaplace_sf
loglaplace_logcdf, loglaplace_logsf
loglaplace_ppf(p, c [, loc, scale]), loglaplace_isf(p, c [, loc, scale])
```


## Bradford

The Bradford family uses positive shape `c` on the finite standardized support
`0 <= z <= 1`, where `z = (x-loc)/scale`.

```text
bradford_pdf(x, c [, loc, scale])
bradford_logpdf, bradford_cdf, bradford_sf
bradford_logcdf, bradford_logsf
bradford_ppf(p, c [, loc, scale]), bradford_isf(p, c [, loc, scale])
```

## Truncated exponential

The truncated-exponential family uses positive shape `b` and standardized
support `0 <= z <= b`, matching `scipy.stats.truncexpon`.

```text
truncexpon_pdf(x, b [, loc, scale])
truncexpon_logpdf, truncexpon_cdf, truncexpon_sf
truncexpon_logcdf, truncexpon_logsf
truncexpon_ppf(p, b [, loc, scale]), truncexpon_isf(p, b [, loc, scale])
```

## Fisk / log-logistic

The Fisk family uses positive shape `c` and support `x > loc`, matching
`scipy.stats.fisk`.

```text
fisk_pdf(x, c [, loc, scale])
fisk_logpdf, fisk_cdf, fisk_sf, fisk_logcdf, fisk_logsf
fisk_ppf(p, c [, loc, scale]), fisk_isf(p, c [, loc, scale])
```

## Double Weibull

The double-Weibull family uses positive shape `c` on the real line, matching
`scipy.stats.dweibull`.

```text
dweibull_pdf(x, c [, loc, scale])
dweibull_logpdf, dweibull_cdf, dweibull_sf
dweibull_logcdf, dweibull_logsf
dweibull_ppf(p, c [, loc, scale]), dweibull_isf(p, c [, loc, scale])
```

## Alpha

The Alpha family uses positive shape `a` and support `x > loc`, matching
`scipy.stats.alpha`. With `z = (x-loc)/scale`, its CDF is
`Phi(a - 1/z) / Phi(a)`.

```text
alpha_pdf(x, a [, loc, scale])
alpha_logpdf, alpha_cdf, alpha_sf, alpha_logcdf, alpha_logsf
alpha_ppf(p, a [, loc, scale]), alpha_isf(p, a [, loc, scale])
```

## Birnbaum-Saunders / fatigue-life

The fatigue-life family uses positive shape `c` and support `x > loc`, matching
`scipy.stats.fatiguelife`.

```text
fatiguelife_pdf(x, c [, loc, scale])
fatiguelife_logpdf, fatiguelife_cdf, fatiguelife_sf
fatiguelife_logcdf, fatiguelife_logsf
fatiguelife_ppf(p, c [, loc, scale]), fatiguelife_isf(p, c [, loc, scale])
```

## Generalized logistic

The generalized-logistic family uses positive shape `c` on the real line,
matching `scipy.stats.genlogistic`.

```text
genlogistic_pdf(x, c [, loc, scale])
genlogistic_logpdf, genlogistic_cdf, genlogistic_sf
genlogistic_logcdf, genlogistic_logsf
genlogistic_ppf(p, c [, loc, scale]), genlogistic_isf(p, c [, loc, scale])
```

## Generalized normal

The generalized-normal family uses positive shape `beta` on the real line,
matching `scipy.stats.gennorm`; `beta = 2` is Gaussian-shaped and `beta = 1`
is Laplace-shaped after scale adjustment.

```text
gennorm_pdf(x, beta [, loc, scale])
gennorm_logpdf, gennorm_cdf, gennorm_sf, gennorm_logcdf, gennorm_logsf
gennorm_ppf(p, beta [, loc, scale]), gennorm_isf(p, beta [, loc, scale])
```

## Nakagami

The Nakagami family uses positive shape `nu` and support `x >= loc`, matching
`scipy.stats.nakagami`.

```text
nakagami_pdf(x, nu [, loc, scale])
nakagami_logpdf, nakagami_cdf, nakagami_sf
nakagami_logcdf, nakagami_logsf
nakagami_ppf(p, nu [, loc, scale]), nakagami_isf(p, nu [, loc, scale])
```

## Power normal

The power-normal family uses positive shape `c` on the real line, matching
`scipy.stats.powernorm`.

```text
powernorm_pdf(x, c [, loc, scale])
powernorm_logpdf, powernorm_cdf, powernorm_sf
powernorm_logcdf, powernorm_logsf
powernorm_ppf(p, c [, loc, scale]), powernorm_isf(p, c [, loc, scale])
```

## Log gamma

The log-gamma family uses positive shape `c` on the real line, matching
`scipy.stats.loggamma`.

```text
loggamma_pdf(x, c [, loc, scale])
loggamma_logpdf, loggamma_cdf, loggamma_sf
loggamma_logcdf, loggamma_logsf
loggamma_ppf(p, c [, loc, scale]), loggamma_isf(p, c [, loc, scale])
```

## Wald

The Wald family has support `x > loc` and is the inverse-Gaussian special case
with SciPy shape `mu = 1`, matching `scipy.stats.wald`.

```text
wald_pdf(x [, loc, scale])
wald_logpdf, wald_cdf, wald_sf, wald_logcdf, wald_logsf
wald_ppf(p [, loc, scale]), wald_isf(p [, loc, scale])
```

## Gumbel right and left

```text
gumbel_r_pdf(x [, loc, scale])
gumbel_r_logpdf(x [, loc, scale])
gumbel_r_cdf(x [, loc, scale])
gumbel_r_sf(x [, loc, scale])
gumbel_r_logcdf(x [, loc, scale])
gumbel_r_logsf(x [, loc, scale])
gumbel_r_ppf(p [, loc, scale])
gumbel_r_isf(p [, loc, scale])

gumbel_l_pdf(x [, loc, scale])
gumbel_l_logpdf(x [, loc, scale])
gumbel_l_cdf(x [, loc, scale])
gumbel_l_sf(x [, loc, scale])
gumbel_l_logcdf(x [, loc, scale])
gumbel_l_logsf(x [, loc, scale])
gumbel_l_ppf(p [, loc, scale])
gumbel_l_isf(p [, loc, scale])
```

For `z = (x - loc)/scale`, the right-skewed standard density is
`exp(-z - exp(-z))`; the left-skewed family is its reflection. Both have
unbounded real support and require finite `loc` and positive finite `scale`.

## Power-function

```text
powerlaw_pdf(x, a [, loc, scale])
powerlaw_logpdf(x, a [, loc, scale])
powerlaw_cdf(x, a [, loc, scale])
powerlaw_sf(x, a [, loc, scale])
powerlaw_logcdf(x, a [, loc, scale])
powerlaw_logsf(x, a [, loc, scale])
powerlaw_ppf(p, a [, loc, scale])
powerlaw_isf(p, a [, loc, scale])
```

The shape `a` must be finite and positive. Standardized support is `[0, 1]`
and the density is `a*z**(a-1)`. Following SciPy's endpoint convention, the
density at the lower endpoint is zero when `a < 1`, one for the standardized
`a = 1` case, and zero when `a > 1`; the limiting divergence for `a < 1` is
seen only from within the support.

## Triangular

```text
triang_pdf(x, c [, loc, scale])
triang_logpdf(x, c [, loc, scale])
triang_cdf(x, c [, loc, scale])
triang_sf(x, c [, loc, scale])
triang_logcdf(x, c [, loc, scale])
triang_logsf(x, c [, loc, scale])
triang_ppf(p, c [, loc, scale])
triang_isf(p, c [, loc, scale])
```

The standardized mode `c` lies in `[0, 1]`, with support `[loc, loc+scale]`.
The endpoint cases `c = 0` and `c = 1` are supported directly rather than by
perturbing the mode into the interior.

## Generalized Pareto

```text
genpareto_pdf(x, c [, loc, scale])
genpareto_logpdf(x, c [, loc, scale])
genpareto_cdf(x, c [, loc, scale])
genpareto_sf(x, c [, loc, scale])
genpareto_logcdf(x, c [, loc, scale])
genpareto_logsf(x, c [, loc, scale])
genpareto_ppf(p, c [, loc, scale])
genpareto_isf(p, c [, loc, scale])
```

For standardized `z >= 0`, the density is
`(1 + c*z)**(-1 - 1/c)`, with the `c = 0` limit evaluated as the exponential
distribution. For `c < 0`, support ends at `z = -1/c`; for `c >= 0`, support
is unbounded above. `c = -1` reduces to the standard uniform distribution.

## Gompertz

```text
gompertz_pdf(x, c [, loc, scale])
gompertz_logpdf(x, c [, loc, scale])
gompertz_cdf(x, c [, loc, scale])
gompertz_sf(x, c [, loc, scale])
gompertz_logcdf(x, c [, loc, scale])
gompertz_logsf(x, c [, loc, scale])
gompertz_ppf(p, c [, loc, scale])
gompertz_isf(p, c [, loc, scale])
```

The shape `c` must be finite and positive. Standardized support is `z >= 0`,
and the standard density is `c*exp(z)*exp(-c*(exp(z)-1))`. Tail calculations
use the cumulative hazard `c*expm1(z)` so the survival and logarithmic survival
remain stable without subtracting a CDF from one.

## Inverse Weibull

```text
invweibull_pdf(x, c [, loc, scale])
invweibull_logpdf(x, c [, loc, scale])
invweibull_cdf(x, c [, loc, scale])
invweibull_sf(x, c [, loc, scale])
invweibull_logcdf(x, c [, loc, scale])
invweibull_logsf(x, c [, loc, scale])
invweibull_ppf(p, c [, loc, scale])
invweibull_isf(p, c [, loc, scale])
```

For positive finite `c`, standardized support is `z > 0` with density
`c*z**(-c-1)*exp(-z**(-c))`. Quantiles and tails are evaluated in the log
domain where needed to avoid intermediate powers overflowing.

## Beta-prime

```text
betaprime_pdf(x, a, b [, loc, scale])
betaprime_logpdf(x, a, b [, loc, scale])
betaprime_cdf(x, a, b [, loc, scale])
betaprime_sf(x, a, b [, loc, scale])
betaprime_logcdf(x, a, b [, loc, scale])
betaprime_logsf(x, a, b [, loc, scale])
betaprime_ppf(p, a, b [, loc, scale])
betaprime_isf(p, a, b [, loc, scale])
```

Both shapes `a` and `b` must be finite and positive. For standardized `z > 0`,
the density is `z**(a-1)*(1+z)**(-a-b)/B(a,b)`. CDF/SF and inverse tails
reuse SciFort's regularized incomplete-beta and inverse-beta kernels with the
transformation `u=z/(1+z)`.

## Burr XII

```text
burr12_pdf(x, c, d [, loc, scale])
burr12_logpdf(x, c, d [, loc, scale])
burr12_cdf(x, c, d [, loc, scale])
burr12_sf(x, c, d [, loc, scale])
burr12_logcdf(x, c, d [, loc, scale])
burr12_logsf(x, c, d [, loc, scale])
burr12_ppf(p, c, d [, loc, scale])
burr12_isf(p, c, d [, loc, scale])
```

Both shapes `c` and `d` must be finite and positive. Standardized support is
`z > 0`, with density
`c*d*z**(c-1)*(1+z**c)**(-d-1)`. The implementation uses softplus/log-domain
identities for `log(1+z**c)` and explicit overflow guards for heavy-tail and
large-shape cases.

## Generalized half-logistic

```text
genhalflogistic_pdf(x, c [, loc, scale])
genhalflogistic_logpdf(x, c [, loc, scale])
genhalflogistic_cdf(x, c [, loc, scale])
genhalflogistic_sf(x, c [, loc, scale])
genhalflogistic_logcdf(x, c [, loc, scale])
genhalflogistic_logsf(x, c [, loc, scale])
genhalflogistic_ppf(p, c [, loc, scale])
genhalflogistic_isf(p, c [, loc, scale])
```

The shape `c` must be finite and positive. Standardized support is
`0 <= z <= 1/c`. The implementation evaluates the power term through
`log1p(-c*z)/c` and uses logarithmic probability transformations in the
quantiles to retain accuracy near both probability endpoints.

## Exponential-power

```text
exponpow_pdf(x, b [, loc, scale])
exponpow_logpdf(x, b [, loc, scale])
exponpow_cdf(x, b [, loc, scale])
exponpow_sf(x, b [, loc, scale])
exponpow_logcdf(x, b [, loc, scale])
exponpow_logsf(x, b [, loc, scale])
exponpow_ppf(p, b [, loc, scale])
exponpow_isf(p, b [, loc, scale])
```

For positive finite shape `b`, standardized support is `z >= 0` and the
density is `b*z**(b-1)*exp(1 + z**b - exp(z**b))`. CDF, survival, and inverse
tails use `expm1`/`log1p` rearrangements to avoid cancellation.

## Exponentiated Weibull

```text
exponweib_pdf(x, a, c [, loc, scale])
exponweib_logpdf(x, a, c [, loc, scale])
exponweib_cdf(x, a, c [, loc, scale])
exponweib_sf(x, a, c [, loc, scale])
exponweib_logcdf(x, a, c [, loc, scale])
exponweib_logsf(x, a, c [, loc, scale])
exponweib_ppf(p, a, c [, loc, scale])
exponweib_isf(p, a, c [, loc, scale])
```

Both shapes `a` and `c` must be finite and positive. On standardized support
`z >= 0`, the CDF is `(1-exp(-z**c))**a`. The implementation forms the base
Weibull CDF and its inverse in the logarithmic domain where needed.

## Power-lognormal

```text
powerlognorm_pdf(x, c, s [, loc, scale])
powerlognorm_logpdf(x, c, s [, loc, scale])
powerlognorm_cdf(x, c, s [, loc, scale])
powerlognorm_sf(x, c, s [, loc, scale])
powerlognorm_logcdf(x, c, s [, loc, scale])
powerlognorm_logsf(x, c, s [, loc, scale])
powerlognorm_ppf(p, c, s [, loc, scale])
powerlognorm_isf(p, c, s [, loc, scale])
```

Both `c` and `s` must be finite and positive. For `z > 0`, letting
`w=log(z)/s`, the survival function is `Phi(-w)**c`. Normal log-CDF and
quantile kernels are reused directly so extreme tails do not depend on
subtraction from one.


## Left-Levy

```text
levy_l_pdf(x [, loc, scale])
levy_l_logpdf(x [, loc, scale])
levy_l_cdf(x [, loc, scale])
levy_l_sf(x [, loc, scale])
levy_l_logcdf(x [, loc, scale])
levy_l_logsf(x [, loc, scale])
levy_l_ppf(p [, loc, scale])
levy_l_isf(p [, loc, scale])
```

`scale` must be finite and positive. The support is `x <= loc`. This family is
implemented as the reflection of the right-Levy law, so its lower and upper
tails delegate to the independently validated Levy kernels without subtracting
from one in the extreme tail.

## Weibull maximum

```text
weibull_max_pdf(x, c [, loc, scale])
weibull_max_logpdf(x, c [, loc, scale])
weibull_max_cdf(x, c [, loc, scale])
weibull_max_sf(x, c [, loc, scale])
weibull_max_logcdf(x, c [, loc, scale])
weibull_max_logsf(x, c [, loc, scale])
weibull_max_ppf(p, c [, loc, scale])
weibull_max_isf(p, c [, loc, scale])
```

The shape `c` and `scale` must be finite and positive. The support is
`x <= loc`. The standardized law is the reflection of Weibull-min, and the
implementation reuses the existing Weibull tail and inverse-tail kernels.

## R-distribution

```text
rdist_pdf(x, c [, loc, scale])
rdist_logpdf(x, c [, loc, scale])
rdist_cdf(x, c [, loc, scale])
rdist_sf(x, c [, loc, scale])
rdist_logcdf(x, c [, loc, scale])
rdist_logsf(x, c [, loc, scale])
rdist_ppf(p, c [, loc, scale])
rdist_isf(p, c [, loc, scale])
```

The shape `c` and `scale` must be finite and positive. Standardized support is
`-1 <= z <= 1`, so the physical support is `[loc-scale, loc+scale]`. Under the
transform `u=(z+1)/2`, `u` has a symmetric beta distribution with shapes
`c/2, c/2`; SciFort therefore reuses its incomplete-beta tails and inverses.

## Skew-Cauchy

```text
skewcauchy_pdf(x, a [, loc, scale])
skewcauchy_logpdf(x, a [, loc, scale])
skewcauchy_cdf(x, a [, loc, scale])
skewcauchy_sf(x, a [, loc, scale])
skewcauchy_logcdf(x, a [, loc, scale])
skewcauchy_logsf(x, a [, loc, scale])
skewcauchy_ppf(p, a [, loc, scale])
skewcauchy_isf(p, a [, loc, scale])
```

The skewness parameter must satisfy `-1 < a < 1`, and `scale` must be finite
and positive. The two half-axes use factors `1-a` and `1+a`. Direct two-argument
arctangent formulas are used for CDF/SF evaluation so large-tail probabilities
do not require subtracting nearly equal numbers.

## Double-gamma

```text
dgamma_pdf(x, a [, loc, scale])
dgamma_logpdf(x, a [, loc, scale])
dgamma_cdf(x, a [, loc, scale])
dgamma_sf(x, a [, loc, scale])
dgamma_logcdf(x, a [, loc, scale])
dgamma_logsf(x, a [, loc, scale])
dgamma_ppf(p, a [, loc, scale])
dgamma_isf(p, a [, loc, scale])
```

The shape `a` and `scale` must be finite and positive. In standardized
coordinates the density is proportional to `abs(z)**(a-1)*exp(-abs(z))`; the
implementation reuses the existing incomplete-gamma tails and inverse tails.

## Asymmetric Laplace

```text
laplace_asymmetric_pdf(x, kappa [, loc, scale])
laplace_asymmetric_logpdf(x, kappa [, loc, scale])
laplace_asymmetric_cdf(x, kappa [, loc, scale])
laplace_asymmetric_sf(x, kappa [, loc, scale])
laplace_asymmetric_logcdf(x, kappa [, loc, scale])
laplace_asymmetric_logsf(x, kappa [, loc, scale])
laplace_asymmetric_ppf(p, kappa [, loc, scale])
laplace_asymmetric_isf(p, kappa [, loc, scale])
```

`kappa` and `scale` must be finite and positive. The standardized density is
`exp(z/kappa)/(kappa+1/kappa)` for `z < 0` and
`exp(-kappa*z)/(kappa+1/kappa)` for `z >= 0`.

## Truncated normal

```text
truncnorm_pdf(x, a, b [, loc, scale])
truncnorm_logpdf(x, a, b [, loc, scale])
truncnorm_cdf(x, a, b [, loc, scale])
truncnorm_sf(x, a, b [, loc, scale])
truncnorm_logcdf(x, a, b [, loc, scale])
truncnorm_logsf(x, a, b [, loc, scale])
truncnorm_ppf(p, a, b [, loc, scale])
truncnorm_isf(p, a, b [, loc, scale])
```

This checkpoint supports finite standardized bounds `a < b`, applied to the
parent normal before `loc` and `scale`, matching SciPy's finite-bound
parameterization. Infinite standardized endpoints remain future work. Stable
normal log-tail differences are used for the truncation mass and bisection is
used for inverse probabilities.

## Log-uniform / reciprocal

```text
loguniform_pdf(x, a, b [, loc, scale])
loguniform_logpdf(x, a, b [, loc, scale])
loguniform_cdf(x, a, b [, loc, scale])
loguniform_sf(x, a, b [, loc, scale])
loguniform_logcdf(x, a, b [, loc, scale])
loguniform_logsf(x, a, b [, loc, scale])
loguniform_ppf(p, a, b [, loc, scale])
loguniform_isf(p, a, b [, loc, scale])
```

The standardized endpoints must satisfy `0 < a < b`; `scale` must be finite
and positive. The density is `1/(z*log(b/a))` on `[a,b]`. If all of `a`, `b`,
`loc`, and `scale` are fitted simultaneously the parameterization is
overidentified; fix `scale=1` to reproduce SciPy's default fit convention.

## Folded normal

```text
foldnorm_pdf(x, c [, loc, scale])
foldnorm_logpdf(x, c [, loc, scale])
foldnorm_cdf(x, c [, loc, scale])
foldnorm_sf(x, c [, loc, scale])
foldnorm_logcdf(x, c [, loc, scale])
foldnorm_logsf(x, c [, loc, scale])
foldnorm_ppf(p, c [, loc, scale])
foldnorm_isf(p, c [, loc, scale])
```

The folding offset `c` must be finite and nonnegative, and `scale` must be
finite and positive. Standardized support is `z >= 0`. Density and tail
probabilities are evaluated as stable mixtures of the two reflected normal
components; quantiles use safeguarded bisection.

## Folded Cauchy

```text
foldcauchy_pdf(x, c [, loc, scale])
foldcauchy_logpdf(x, c [, loc, scale])
foldcauchy_cdf(x, c [, loc, scale])
foldcauchy_sf(x, c [, loc, scale])
foldcauchy_logcdf(x, c [, loc, scale])
foldcauchy_logsf(x, c [, loc, scale])
foldcauchy_ppf(p, c [, loc, scale])
foldcauchy_isf(p, c [, loc, scale])
```

The folding offset `c` must be finite and nonnegative. Standardized support is
`z >= 0`. Direct two-argument arctangent formulas are used for the upper tail
to avoid cancellation, and quantiles use safeguarded bisection.

## Reciprocal inverse-Gaussian

```text
recipinvgauss_pdf(x, mu [, loc, scale])
recipinvgauss_logpdf(x, mu [, loc, scale])
recipinvgauss_cdf(x, mu [, loc, scale])
recipinvgauss_sf(x, mu [, loc, scale])
recipinvgauss_logcdf(x, mu [, loc, scale])
recipinvgauss_logsf(x, mu [, loc, scale])
recipinvgauss_ppf(p, mu [, loc, scale])
recipinvgauss_isf(p, mu [, loc, scale])
```

The shape `mu` and `scale` must be finite and positive. Standardized support is
`z > 0`. SciFort uses the reciprocal relationship to its existing
inverse-Gaussian distribution for the cumulative probabilities and quantiles.

## Truncated Pareto

```text
truncpareto_pdf(x, b, c [, loc, scale])
truncpareto_logpdf(x, b, c [, loc, scale])
truncpareto_cdf(x, b, c [, loc, scale])
truncpareto_sf(x, b, c [, loc, scale])
truncpareto_logcdf(x, b, c [, loc, scale])
truncpareto_logsf(x, b, c [, loc, scale])
truncpareto_ppf(p, b, c [, loc, scale])
truncpareto_isf(p, b, c [, loc, scale])
```

`b` must be finite and nonzero, `c` must be finite and greater than one, and
`scale` must be finite and positive. Standardized support is `1 <= z <= c`.
Both positive and negative `b` are supported; logarithmic normalizers and
quantile mixtures use stable `expm1`/log-domain rearrangements.

## Exponentially modified normal

```text
exponnorm_pdf(x, k [, loc, scale])
exponnorm_logpdf(x, k [, loc, scale])
exponnorm_cdf(x, k [, loc, scale])
exponnorm_sf(x, k [, loc, scale])
exponnorm_logcdf(x, k [, loc, scale])
exponnorm_logsf(x, k [, loc, scale])
exponnorm_ppf(p, k [, loc, scale])
exponnorm_isf(p, k [, loc, scale])
```

`k` must be finite and positive. The standardized law is the sum of a unit
normal variate and an independent exponential variate with mean `k`, matching
`scipy.stats.exponnorm`. Logarithmic tails are evaluated directly and the
quantiles use safeguarded bisection over the monotone CDF/SF kernels.

## Johnson SB

```text
johnsonsb_pdf(x, a, b [, loc, scale])
johnsonsb_logpdf(x, a, b [, loc, scale])
johnsonsb_cdf(x, a, b [, loc, scale])
johnsonsb_sf(x, a, b [, loc, scale])
johnsonsb_logcdf(x, a, b [, loc, scale])
johnsonsb_logsf(x, a, b [, loc, scale])
johnsonsb_ppf(p, a, b [, loc, scale])
johnsonsb_isf(p, a, b [, loc, scale])
```

`a` must be finite and `b` finite and positive. Standardized support is
`0 < z < 1`; with location and scale the support is `loc < x < loc+scale`.
The implementation applies the normal law to `a + b*logit(z)` and uses a
stable logistic inverse for quantiles.

## Johnson SU

```text
johnsonsu_pdf(x, a, b [, loc, scale])
johnsonsu_logpdf(x, a, b [, loc, scale])
johnsonsu_cdf(x, a, b [, loc, scale])
johnsonsu_sf(x, a, b [, loc, scale])
johnsonsu_logcdf(x, a, b [, loc, scale])
johnsonsu_logsf(x, a, b [, loc, scale])
johnsonsu_ppf(p, a, b [, loc, scale])
johnsonsu_isf(p, a, b [, loc, scale])
```

`a` must be finite and `b` finite and positive. Support is the whole real
line. The transformed normal coordinate is `a + b*asinh(z)`; inverse
transforms use guarded hyperbolic-sine evaluation for extreme probabilities.

## Trapezoid

```text
trapezoid_pdf(x, c, d [, loc, scale])
trapezoid_logpdf(x, c, d [, loc, scale])
trapezoid_cdf(x, c, d [, loc, scale])
trapezoid_sf(x, c, d [, loc, scale])
trapezoid_logcdf(x, c, d [, loc, scale])
trapezoid_logsf(x, c, d [, loc, scale])
trapezoid_ppf(p, c, d [, loc, scale])
trapezoid_isf(p, c, d [, loc, scale])
```

The standardized plateau edges satisfy `0 <= c <= d <= 1`; support is
`0 <= z <= 1`. The density rises linearly to the plateau, remains constant
from `c` through `d`, and falls linearly to zero. Exact piecewise inverses are
used. Likelihood scores are nonregular when observations lie exactly on moving
plateau or support boundaries and follow SciFort's boundary-NaN convention.

## Burr Type III

```text
burr_pdf(x, c, d [, loc, scale])
burr_logpdf(x, c, d [, loc, scale])
burr_cdf(x, c, d [, loc, scale])
burr_sf(x, c, d [, loc, scale])
burr_logcdf(x, c, d [, loc, scale])
burr_logsf(x, c, d [, loc, scale])
burr_ppf(p, c, d [, loc, scale])
burr_isf(p, c, d [, loc, scale])
```

`c` and `d` must be finite and positive. Standardized support is `z >= 0`.
The implementation uses direct log-domain CDF/survival formulas and closed-form
inverses, avoiding subtraction from one in the upper tail.

## Mielke beta-kappa / Dagum

```text
mielke_pdf(x, k, s [, loc, scale])
mielke_logpdf(x, k, s [, loc, scale])
mielke_cdf(x, k, s [, loc, scale])
mielke_sf(x, k, s [, loc, scale])
mielke_logcdf(x, k, s [, loc, scale])
mielke_logsf(x, k, s [, loc, scale])
mielke_ppf(p, k, s [, loc, scale])
mielke_isf(p, k, s [, loc, scale])
```

`k` and `s` must be finite and positive. Standardized support is `z >= 0`.
This parameterization is mathematically equivalent to Burr Type III with
`c=s` and `d=k/s`; SciFort evaluates it directly to preserve its public
parameterization and stable shape arithmetic.

## Gibrat

```text
gibrat_pdf(x [, loc, scale])
gibrat_logpdf(x [, loc, scale])
gibrat_cdf(x [, loc, scale])
gibrat_sf(x [, loc, scale])
gibrat_logcdf(x [, loc, scale])
gibrat_logsf(x [, loc, scale])
gibrat_ppf(p [, loc, scale])
gibrat_isf(p [, loc, scale])
```

Standardized support is `z > 0`. The Gibrat distribution is the lognormal law
with shape parameter one, so these routines reuse SciFort's lognormal kernels.
`scale` must be finite and positive.

## Wrapped Cauchy

```text
wrapcauchy_pdf(x, c [, loc, scale])
wrapcauchy_logpdf(x, c [, loc, scale])
wrapcauchy_cdf(x, c [, loc, scale])
wrapcauchy_sf(x, c [, loc, scale])
wrapcauchy_logcdf(x, c [, loc, scale])
wrapcauchy_logsf(x, c [, loc, scale])
wrapcauchy_ppf(p, c [, loc, scale])
wrapcauchy_isf(p, c [, loc, scale])
```

`c` must be finite with `0 < c < 1`. Standardized support is
`0 <= z <= 2*pi`; after location and scale it is
`loc <= x <= loc + 2*pi*scale`. The CDF and inverse use quadrant-preserving
`atan2` formulas, and the density denominator is evaluated in a cancellation-
resistant form near `c=1` and the angular origin.

## Generalized extreme-value

```text
genextreme_pdf(x, c [, loc, scale])
genextreme_logpdf(x, c [, loc, scale])
genextreme_cdf(x, c [, loc, scale])
genextreme_sf(x, c [, loc, scale])
genextreme_logcdf(x, c [, loc, scale])
genextreme_logsf(x, c [, loc, scale])
genextreme_ppf(p, c [, loc, scale])
genextreme_isf(p, c [, loc, scale])
```

`c` must be finite. SciFort follows SciPy's sign convention
`F(z)=exp(-(1-c*z)**(1/c))`, with the `c=0` limit equal to the right Gumbel
law. For `c>0` the standardized upper endpoint is `1/c`; for `c<0` the
standardized lower endpoint is `1/c`. Direct logarithmic tails and a
closed-form Gumbel transform are used.

## Kappa-3

```text
kappa3_pdf(x, a [, loc, scale])
kappa3_logpdf(x, a [, loc, scale])
kappa3_cdf(x, a [, loc, scale])
kappa3_sf(x, a [, loc, scale])
kappa3_logcdf(x, a [, loc, scale])
kappa3_logsf(x, a [, loc, scale])
kappa3_ppf(p, a [, loc, scale])
kappa3_isf(p, a [, loc, scale])
```

`a` must be finite and positive and standardized support is `z>0`. The CDF is
`(1+a*z**(-a))**(-1/a)`. Log-domain evaluation and the analytic inverse avoid
subtraction from one in the tails.

## Kappa-4

```text
kappa4_pdf(x, h, k [, loc, scale])
kappa4_logpdf(x, h, k [, loc, scale])
kappa4_cdf(x, h, k [, loc, scale])
kappa4_sf(x, h, k [, loc, scale])
kappa4_logcdf(x, h, k [, loc, scale])
kappa4_logsf(x, h, k [, loc, scale])
kappa4_ppf(p, h, k [, loc, scale])
kappa4_isf(p, h, k [, loc, scale])
```

`h` and `k` must be finite. The standardized law uses the nested transform
`T=(1-k*z)**(1/k)` and `F=(1-h*T)**(1/h)`, with exponential limits when either
shape is zero. Support depends on both shapes and is handled explicitly.
Useful special cases are `h=0` (generalized extreme-value with shape `k`),
`h=1` (generalized Pareto with shape `-k`), and `h=-1,k=0` (logistic).

## Doubly truncated Weibull minimum

```text
truncweibull_min_pdf(x, c, a, b [, loc, scale])
truncweibull_min_logpdf(x, c, a, b [, loc, scale])
truncweibull_min_cdf(x, c, a, b [, loc, scale])
truncweibull_min_sf(x, c, a, b [, loc, scale])
truncweibull_min_logcdf(x, c, a, b [, loc, scale])
truncweibull_min_logsf(x, c, a, b [, loc, scale])
truncweibull_min_ppf(p, c, a, b [, loc, scale])
truncweibull_min_isf(p, c, a, b [, loc, scale])
```

`c` must be finite and positive, with finite standardized truncation points
`0<=a<b`. Support is `a<=z<=b`. The implementation normalizes the Weibull
minimum law between the two truncation points, evaluates both tails directly
in the log domain, and uses a closed-form log-sum-exp inverse.

## Generalized gamma

```text
gengamma_pdf(x, a, c [, loc, scale])
gengamma_logpdf(x, a, c [, loc, scale])
gengamma_cdf(x, a, c [, loc, scale])
gengamma_sf(x, a, c [, loc, scale])
gengamma_logcdf(x, a, c [, loc, scale])
gengamma_logsf(x, a, c [, loc, scale])
gengamma_ppf(p, a, c [, loc, scale])
gengamma_isf(p, a, c [, loc, scale])
```

`a` must be finite and positive and `c` finite and nonzero. On standardized
support `z>0`, the density is
`abs(c)*z**(c*a-1)*exp(-z**c)/Gamma(a)`. Positive and negative `c` use the
appropriate lower or upper regularized incomplete-gamma tail directly, and the
quantiles use the matching inverse incomplete-gamma function. The lower support
density has its mathematical finite, zero, or singular limit when `c>0`; it is
zero for `c<0`.

## Half-generalized normal

```text
halfgennorm_pdf(x, beta [, loc, scale])
halfgennorm_logpdf(x, beta [, loc, scale])
halfgennorm_cdf(x, beta [, loc, scale])
halfgennorm_sf(x, beta [, loc, scale])
halfgennorm_logcdf(x, beta [, loc, scale])
halfgennorm_logsf(x, beta [, loc, scale])
halfgennorm_ppf(p, beta [, loc, scale])
halfgennorm_isf(p, beta [, loc, scale])
```

`beta` must be finite and positive and standardized support is `z>=0`. The
density is `beta*exp(-z**beta)/Gamma(1/beta)`. CDF, survival, and inverse
functions are evaluated through the regularized incomplete-gamma kernels.
`beta=1` is exponential, while `beta=2` is a half-normal with scale
`1/sqrt(2)` relative to the same standardized coordinate.

## ARGUS

```text
argus_pdf(x, chi [, loc, scale])
argus_logpdf(x, chi [, loc, scale])
argus_cdf(x, chi [, loc, scale])
argus_sf(x, chi [, loc, scale])
argus_logcdf(x, chi [, loc, scale])
argus_logsf(x, chi [, loc, scale])
argus_ppf(p, chi [, loc, scale])
argus_isf(p, chi [, loc, scale])
```

`chi` must be finite and positive. Standardized support is `0<=z<=1`. The
normalizer is evaluated as one half of the regularized lower incomplete gamma
`P(3/2, chi**2/2)`, avoiding cancellation in the equivalent normal-CDF
expression. The survival function is formed as a direct ratio of those
normalizers. PPF and ISF use safeguarded bisection on the finite support, with
the lower or upper tail selected according to the requested probability.

## Erlang

```text
erlang_pdf(x, a [, loc, scale])
erlang_logpdf(x, a [, loc, scale])
erlang_cdf(x, a [, loc, scale])
erlang_sf(x, a [, loc, scale])
erlang_logcdf(x, a [, loc, scale])
erlang_logsf(x, a [, loc, scale])
erlang_ppf(p, a [, loc, scale])
erlang_isf(p, a [, loc, scale])
```

The Erlang law delegates to the gamma kernels. Conventionally `a` is a positive
integer. For numerical compatibility with SciPy's computational behavior,
SciFort accepts any finite `a>0`; unlike SciPy's Python interface, this Fortran
numerical layer does not emit a warning when `a` is noninteger.

## Crystal Ball

```text
crystalball_pdf(x, beta, m [, loc, scale])
crystalball_logpdf(x, beta, m [, loc, scale])
crystalball_cdf(x, beta, m [, loc, scale])
crystalball_sf(x, beta, m [, loc, scale])
crystalball_logcdf(x, beta, m [, loc, scale])
crystalball_logsf(x, beta, m [, loc, scale])
crystalball_ppf(p, beta, m [, loc, scale])
crystalball_isf(p, beta, m [, loc, scale])
```

`beta` must be finite and positive and `m` finite and greater than one. The
standardized density is Gaussian for `z>-beta` and has the usual matched power
law below `-beta`. The normalization, lower power-law tail, Gaussian upper
tail, and both inverse branches are evaluated directly rather than by
subtracting nearly equal probabilities.

## Jones-Faddy skew-t

```text
jf_skew_t_pdf(x, a, b [, loc, scale])
jf_skew_t_logpdf(x, a, b [, loc, scale])
jf_skew_t_cdf(x, a, b [, loc, scale])
jf_skew_t_sf(x, a, b [, loc, scale])
jf_skew_t_logcdf(x, a, b [, loc, scale])
jf_skew_t_logsf(x, a, b [, loc, scale])
jf_skew_t_ppf(p, a, b [, loc, scale])
jf_skew_t_isf(p, a, b [, loc, scale])
```

`a` and `b` must be finite and positive. With
`r=z/sqrt(a+b+z**2)`, the lower-tail beta coordinate is `(1+r)/2`.
CDF, survival, and quantiles therefore use the regularized incomplete-beta
kernels and their direct inverses. When `a=b`, the law reduces to Student t
with `2*a` degrees of freedom.

## Pearson III

```text
pearson3_pdf(x, skew [, loc, scale])
pearson3_logpdf(x, skew [, loc, scale])
pearson3_cdf(x, skew [, loc, scale])
pearson3_sf(x, skew [, loc, scale])
pearson3_logcdf(x, skew [, loc, scale])
pearson3_logsf(x, skew [, loc, scale])
pearson3_ppf(p, skew [, loc, scale])
pearson3_isf(p, skew [, loc, scale])
```

`skew` must be finite. For nonzero skew the implementation maps the
standardized variate to a gamma law with `beta=2/skew` and
`alpha=beta**2`, reversing lower and upper tails when `skew<0`. A small-skew
normal branch matches SciPy's numerically stable limiting behavior, and
`skew=0` is the standard normal exactly at the API level.

## Relativistic Breit-Wigner

```text
rel_breitwigner_pdf(x, rho [, loc, scale])
rel_breitwigner_logpdf(x, rho [, loc, scale])
rel_breitwigner_cdf(x, rho [, loc, scale])
rel_breitwigner_sf(x, rho [, loc, scale])
rel_breitwigner_logcdf(x, rho [, loc, scale])
rel_breitwigner_logsf(x, rho [, loc, scale])
rel_breitwigner_ppf(p, rho [, loc, scale])
rel_breitwigner_isf(p, rho [, loc, scale])
```

`rho` must be finite and positive and standardized support is `z>=0`. The
density uses the normalized relativistic Breit-Wigner form. The CDF uses the
complex antiderivative employed by SciPy 1.17.0; its provenance and retained
BSD-3-Clause notice are recorded in `CODE_PROVENANCE.md` and
`THIRD_PARTY_LICENSES.md`. The survival function uses the reciprocal-argument
form of the same antiderivative to preserve upper-tail accuracy. PPF and ISF
use tail-aware safeguarded bisection.


## Generalized exponential

```text
genexpon_pdf(x, a, b, c [, loc, scale])
genexpon_logpdf(x, a, b, c [, loc, scale])
genexpon_cdf(x, a, b, c [, loc, scale])
genexpon_sf(x, a, b, c [, loc, scale])
genexpon_logcdf(x, a, b, c [, loc, scale])
genexpon_logsf(x, a, b, c [, loc, scale])
genexpon_ppf(p, a, b, c [, loc, scale])
genexpon_isf(p, a, b, c [, loc, scale])
```

`a`, `b`, and `c` must be finite and positive, and standardized support is
`z>=0`. The survival function is evaluated directly in the log domain. PPF and
ISF invert that survival law through an independently implemented safeguarded
principal-real Lambert-W solver on `[-1/e,0]`, avoiding a generic root search.

## Skew-normal

```text
skewnorm_pdf(x, a [, loc, scale])
skewnorm_logpdf(x, a [, loc, scale])
skewnorm_cdf(x, a [, loc, scale])
skewnorm_sf(x, a [, loc, scale])
skewnorm_logcdf(x, a [, loc, scale])
skewnorm_logsf(x, a [, loc, scale])
skewnorm_ppf(p, a [, loc, scale])
skewnorm_isf(p, a [, loc, scale])
```

`a` may be any finite shape value. The density is `2*phi(z)*Phi(a*z)` and
`a=0` is exactly the standard normal. CDF evaluation uses independent
Gauss-Legendre quadrature of Owen's T integral, with a direct transformed
lower-tail integral in the difficult `z<0, a>0` quadrant. Survival values use
the reflected lower tail. Quantiles use tail-aware safeguarded bisection.

## Tukey lambda

```text
tukeylambda_pdf(x, lam [, loc, scale])
tukeylambda_logpdf(x, lam [, loc, scale])
tukeylambda_cdf(x, lam [, loc, scale])
tukeylambda_sf(x, lam [, loc, scale])
tukeylambda_logcdf(x, lam [, loc, scale])
tukeylambda_logsf(x, lam [, loc, scale])
tukeylambda_ppf(p, lam [, loc, scale])
tukeylambda_isf(p, lam [, loc, scale])
```

`lam` may be any finite value. The standardized quantile is
`(p**lam-(1-p)**lam)/lam`, with its logarithmic limit used at `lam=0` and a
small-`lam` series to avoid cancellation. For positive `lam` the support is
`[-1/lam,1/lam]`; otherwise it is unbounded. CDF and density are obtained by
inverting/differentiating the monotone quantile. `lam=0` is logistic and
`lam=1` is uniform on `[-1,1]`.

## Rice

```text
rice_pdf(x, b [, loc, scale])
rice_logpdf(x, b [, loc, scale])
rice_cdf(x, b [, loc, scale])
rice_sf(x, b [, loc, scale])
rice_logcdf(x, b [, loc, scale])
rice_logsf(x, b [, loc, scale])
rice_ppf(p, b [, loc, scale])
rice_isf(p, b [, loc, scale])
```

`b` must be finite and nonnegative, and standardized support is `z>=0`. The
log-density uses a scaled modified-Bessel `I0` evaluation so `exp(b*z)` is not
formed separately. Lower and upper tails are computed directly through the
Poisson/incomplete-gamma representation of the equivalent noncentral
chi-square law. Quantiles use the requested tail directly in safeguarded
bisection. `b=0` reduces to the Rayleigh law.

## Double-Pareto lognormal

```text
dpareto_lognorm_pdf(x, u, s, a, b [, loc, scale])
dpareto_lognorm_logpdf(x, u, s, a, b [, loc, scale])
dpareto_lognorm_cdf(x, u, s, a, b [, loc, scale])
dpareto_lognorm_sf(x, u, s, a, b [, loc, scale])
dpareto_lognorm_logcdf(x, u, s, a, b [, loc, scale])
dpareto_lognorm_logsf(x, u, s, a, b [, loc, scale])
dpareto_lognorm_ppf(p, u, s, a, b [, loc, scale])
dpareto_lognorm_isf(p, u, s, a, b [, loc, scale])
```

`u` may be any finite log-location shape; `s`, `a`, and `b` must be finite
and positive. Standardized support is `z>0`. Density and tail expressions are
evaluated in log space with normal Mills ratios, and quantiles are inverted in
`log(z)` so the Pareto tails do not force premature overflow. When `u=0` and
`a=b`, the log distribution is symmetric and `CDF(x)=SF(1/x)`.

## Circular von Mises

```text
vonmises_pdf(x, kappa [, loc, scale])
vonmises_logpdf(x, kappa [, loc, scale])
vonmises_cdf(x, kappa [, loc, scale])
vonmises_sf(x, kappa [, loc, scale])
vonmises_logcdf(x, kappa [, loc, scale])
vonmises_logsf(x, kappa [, loc, scale])
vonmises_ppf(p, kappa [, loc, scale])
vonmises_isf(p, kappa [, loc, scale])
```

`kappa` must be finite and nonnegative. The density uses the scaled modified
Bessel function `i0e(kappa)`. Following SciPy's circular convention, the CDF is
periodically extended so `CDF(x + 2*pi*scale)=CDF(x)+1`; PPF and ISF return the
central-period representative. The central-period CDF is evaluated by adaptive
quadrature and quantiles use safeguarded bisection.

## Von Mises on a line

```text
vonmises_line_pdf(x, kappa [, loc, scale])
vonmises_line_logpdf(x, kappa [, loc, scale])
vonmises_line_cdf(x, kappa [, loc, scale])
vonmises_line_sf(x, kappa [, loc, scale])
vonmises_line_logcdf(x, kappa [, loc, scale])
vonmises_line_logsf(x, kappa [, loc, scale])
vonmises_line_ppf(p, kappa [, loc, scale])
vonmises_line_isf(p, kappa [, loc, scale])
```

`kappa` must be finite and nonnegative. Support is
`[loc-pi*scale, loc+pi*scale]`. Inside that interval the density and CDF agree
with the central period of `vonmises`; outside it the line distribution has
ordinary finite-support tails. At `kappa=0` it is uniform on the support.

## Asymptotic two-sided Kolmogorov

```text
kstwobign_pdf(x [, loc, scale])
kstwobign_logpdf(x [, loc, scale])
kstwobign_cdf(x [, loc, scale])
kstwobign_sf(x [, loc, scale])
kstwobign_logcdf(x [, loc, scale])
kstwobign_logsf(x [, loc, scale])
kstwobign_ppf(p [, loc, scale])
kstwobign_isf(p [, loc, scale])
```

Standardized support is `z>0`. This is SciPy's `kstwobign` asymptotic
Kolmogorov-Smirnov law. For small `z`, SciFort evaluates the Jacobi-theta
representation of the CDF; for larger `z`, it evaluates the alternating
Kolmogorov survival series. The switch avoids catastrophic tail cancellation.
PDF and analytic score derivatives use differentiated versions of the same
series, and PPF/ISF use the requested tail directly in safeguarded bisection.

## Irwin-Hall

```text
irwinhall_pdf(x, n [, loc, scale])
irwinhall_logpdf(x, n [, loc, scale])
irwinhall_cdf(x, n [, loc, scale])
irwinhall_sf(x, n [, loc, scale])
irwinhall_logcdf(x, n [, loc, scale])
irwinhall_logsf(x, n [, loc, scale])
irwinhall_ppf(p, n [, loc, scale])
irwinhall_isf(p, n [, loc, scale])
```

`n` must be a positive integer-valued `real(dp)`. Standardized support is
`0 <= z <= n`. The density is evaluated as a cardinal B-spline recurrence,
which avoids the large cancellation in the direct alternating-power formula;
the CDF is accumulated from the corresponding order-`n+1` spline identity.
Quantiles use safeguarded bisection. The continuous score contains only
location and scale because `n` is integral, and the bounded fitter enforces an
integer constraint on `n`. At `n=1` the standardized law is uniform on `[0,1]`.

## One-sided finite-sample Kolmogorov-Smirnov

```text
ksone_pdf(x, n [, loc, scale])
ksone_logpdf(x, n [, loc, scale])
ksone_cdf(x, n [, loc, scale])
ksone_sf(x, n [, loc, scale])
ksone_logcdf(x, n [, loc, scale])
ksone_logsf(x, n [, loc, scale])
ksone_ppf(p, n [, loc, scale])
ksone_isf(p, n [, loc, scale])
```

`n` must be a positive integer-valued `real(dp)` and standardized support is
`0 <= z <= 1`. SciFort evaluates the exact finite-sample Smirnov survival sum
with a shared logarithmic scale and differentiates the same sum for density
and location/scale score calculations. At the finite set of summation knots,
the density uses the symmetric limiting average while the derivative score is
reported as undefined. Quantiles use the requested tail directly in
safeguarded bisection. The fitter enforces integral `n`; the continuous score
contains location and scale only. For `n=1` the law is uniform on `[0,1]`.

## Two-sided finite-sample Kolmogorov-Smirnov

```text
kstwo_pdf(x, n [, loc, scale])
kstwo_logpdf(x, n [, loc, scale])
kstwo_cdf(x, n [, loc, scale])
kstwo_sf(x, n [, loc, scale])
kstwo_logcdf(x, n [, loc, scale])
kstwo_logsf(x, n [, loc, scale])
kstwo_ppf(p, n [, loc, scale])
kstwo_isf(p, n [, loc, scale])
```

`n` is a positive integer sample size. The standardized support is
`1/(2*n) <= x <= 1`. CDF/SF evaluation uses exact endpoint formulas, the
Durbin/Marsaglia-Tsang-Wang matrix method for the small finite-sample regime,
the exact `2*ksone` one-sided identity in the appropriate upper-tail regimes,
and the Pelz-Good expansion for large samples. CDF and SF select their own
regimes so small upper tails are not formed by subtracting a rounded CDF.
PPF/ISF bisect the requested direct tail. The density follows the derivative of
the finite-sample CDF; away from exact closed-form regions SciFort uses a
five-point local derivative, as does SciPy's finite-sample implementation in
spirit. Consequently the loc/scale score differentiates the numerical density
locally; `n` remains an integral fit parameter and is not differentiated.

## Levy-stable

```text
levy_stable_pdf(x, alpha, beta [, loc, scale])
levy_stable_logpdf(x, alpha, beta [, loc, scale])
levy_stable_cdf(x, alpha, beta [, loc, scale])
levy_stable_sf(x, alpha, beta [, loc, scale])
levy_stable_logcdf(x, alpha, beta [, loc, scale])
levy_stable_logsf(x, alpha, beta [, loc, scale])
levy_stable_ppf(p, alpha, beta [, loc, scale])
levy_stable_isf(p, alpha, beta [, loc, scale])
```

`0 < alpha <= 2` and `-1 <= beta <= 1`. SciFort implements SciPy's default
S1 parameterization. Evaluation converts to the Zolotarev/Nolan S0 form for the
piecewise density and CDF quadratures, including SciPy-compatible difficult
input rounding near `alpha=1` and the S0 shift point. The exact reductions
`alpha=2` to a normal with standardized standard deviation `sqrt(2)`,
`alpha=1, beta=0` to Cauchy, and `alpha=0.5, beta=1` to Levy are handled
directly. In the S1 parameterization the `alpha=1` loc/scale transform includes
the logarithmic scale shift `2*beta*scale*log(scale)/pi`. Quantiles bisect the
requested direct tail. Shape, location, and scale likelihood scores are
numerical finite differences.

## Studentized range

```text
studentized_range_pdf(x, k, df [, loc, scale])
studentized_range_logpdf(x, k, df [, loc, scale])
studentized_range_cdf(x, k, df [, loc, scale])
studentized_range_sf(x, k, df [, loc, scale])
studentized_range_logcdf(x, k, df [, loc, scale])
studentized_range_logsf(x, k, df [, loc, scale])
studentized_range_ppf(p, k, df [, loc, scale])
studentized_range_isf(p, k, df [, loc, scale])
```

`k > 1`, `df > 0`, and standardized support is `x >= 0`. For finite `df`,
SciFort conditions on `S=sqrt(chi-square(df)/df)` and integrates the normal
range CDF/PDF over that scale mixture. For `df >= 100000` it uses the
infinite-degrees-of-freedom normal-range integrals, matching SciPy's asymptotic
switch. Quantiles bisect the requested direct tail. The likelihood score with
respect to `k`, `df`, location, and scale is evaluated by finite differences.

## Noncentral chi-square

```text
ncx2_pdf(x, df, nc [, loc, scale])
ncx2_logpdf(x, df, nc [, loc, scale])
ncx2_cdf(x, df, nc [, loc, scale])
ncx2_sf(x, df, nc [, loc, scale])
ncx2_logcdf(x, df, nc [, loc, scale])
ncx2_logsf(x, df, nc [, loc, scale])
ncx2_ppf(p, df, nc [, loc, scale])
ncx2_isf(p, df, nc [, loc, scale])
```

`df` must be finite and positive and `nc` finite and nonnegative. Standardized
support is `z>=0`. Density and both tails use a Poisson mixture of central
chi-square terms centered near the modal Poisson index; lower and upper tails
use `gammainc` and `gammaincc` directly rather than subtraction. Analytic
scores differentiate the mixture with respect to `df`, `nc`, location, and
scale. Quantiles use direct-tail safeguarded bisection. `nc=0` reduces exactly
to the central chi-square implementation.

## Noncentral F

```text
ncf_pdf(x, dfn, dfd, nc [, loc, scale])
ncf_logpdf(x, dfn, dfd, nc [, loc, scale])
ncf_cdf(x, dfn, dfd, nc [, loc, scale])
ncf_sf(x, dfn, dfd, nc [, loc, scale])
ncf_logcdf(x, dfn, dfd, nc [, loc, scale])
ncf_logsf(x, dfn, dfd, nc [, loc, scale])
ncf_ppf(p, dfn, dfd, nc [, loc, scale])
ncf_isf(p, dfn, dfd, nc [, loc, scale])
```

`dfn` and `dfd` must be finite and positive and `nc` finite and nonnegative.
Standardized support is `z>=0`. The distribution is evaluated as a centered
Poisson mixture of beta-density and regularized incomplete-beta terms after
the usual F-to-beta transform. CDF and survival use `betainc` and `betaincc`
directly, analytic scores differentiate the mixture in all three shape
parameters plus location and scale, and PPF/ISF use direct-tail safeguarded
bisection. `nc=0` reduces exactly to the central F implementation.

## Noncentral Student t

```text
nct_pdf(x, df, nc [, loc, scale])
nct_logpdf(x, df, nc [, loc, scale])
nct_cdf(x, df, nc [, loc, scale])
nct_sf(x, df, nc [, loc, scale])
nct_logcdf(x, df, nc [, loc, scale])
nct_logsf(x, df, nc [, loc, scale])
nct_ppf(probability, df, nc [, loc, scale])
nct_isf(probability, df, nc [, loc, scale])
```

`df` must be positive and `nc` finite; `df=+infinity` reduces to a normal law
shifted by `nc`. Support is the real line. For finite `df`, SciFort writes the
noncentral-t law as an expectation over a `Gamma(df/2,1)` variate and integrates
in gamma-probability coordinates, avoiding direct quadrature of the singular
gamma density when `df<2`. A double-exponential map resolves both endpoints;
the quadrature spacing is refined for `df<1` and again for `df<0.3`. CDF and
SF are integrated directly in their requested tail, and PPF/ISF bisect on the
corresponding log tail. `nc=0` reduces to the central Student-t implementation.

## Gauss hypergeometric

```text
gausshyper_pdf(x, a, b, c, z [, loc, scale])
gausshyper_logpdf(x, a, b, c, z [, loc, scale])
gausshyper_cdf(x, a, b, c, z [, loc, scale])
gausshyper_sf(x, a, b, c, z [, loc, scale])
gausshyper_logcdf(x, a, b, c, z [, loc, scale])
gausshyper_logsf(x, a, b, c, z [, loc, scale])
gausshyper_ppf(probability, a, b, c, z [, loc, scale])
gausshyper_isf(probability, a, b, c, z [, loc, scale])
```

`a` and `b` must be positive, `c` finite, and `z>-1`; standardized support is
`0<x<1`. The density is proportional to
`x**(a-1)*(1-x)**(b-1)*(1+z*x)**(-c)`. Normalization and both tails are
evaluated after transforming the integral through the beta quantile, so the
possibly singular beta endpoint factors are absorbed analytically rather than
being sampled by the quadrature. PPF/ISF use direct-tail safeguarded bisection.
Either `c=0` or `z=0` reduces exactly to the beta family.

## Landau

```text
landau_pdf(x [, loc, scale])
landau_logpdf(x [, loc, scale])
landau_cdf(x [, loc, scale])
landau_sf(x [, loc, scale])
landau_logcdf(x [, loc, scale])
landau_logsf(x [, loc, scale])
landau_ppf(probability [, loc, scale])
landau_isf(probability [, loc, scale])
```

Support is the whole real line. The standard distribution matches
`scipy.stats.landau`; ordinary `loc`/`scale` standardization is used. SciFort
combines a left-tail asymptotic expansion, characteristic-function quadrature
in the central region, and a scaled nonoscillatory defining integral in the
heavy right tail. Lower and upper tails are evaluated directly, and PPF/ISF
bisect on log-CDF/log-SF for small requested probabilities.

## Generalized hyperbolic

```text
genhyperbolic_pdf(x, p, a, b [, loc, scale])
genhyperbolic_logpdf(x, p, a, b [, loc, scale])
genhyperbolic_cdf(x, p, a, b [, loc, scale])
genhyperbolic_sf(x, p, a, b [, loc, scale])
genhyperbolic_logcdf(x, p, a, b [, loc, scale])
genhyperbolic_logsf(x, p, a, b [, loc, scale])
genhyperbolic_ppf(probability, p, a, b [, loc, scale])
genhyperbolic_isf(probability, p, a, b [, loc, scale])
```

`a` must be positive. For `p >= 0`, `abs(b) < a`; for `p < 0`, the
closed boundary `abs(b) = a` is also valid. Support is the real line. The
density uses scaled-log modified-Bessel-K evaluation so the normalizer remains
stable when its Bessel argument is large. CDF/SF calculations use the
substitution `z=sinh(y)` and integrate the requested tail directly in log
space. `p=-0.5` reduces to the normal-inverse-Gaussian family.

## Generalized inverse Gaussian

```text
geninvgauss_pdf(x, p, b [, loc, scale])
geninvgauss_logpdf(x, p, b [, loc, scale])
geninvgauss_cdf(x, p, b [, loc, scale])
geninvgauss_sf(x, p, b [, loc, scale])
geninvgauss_logcdf(x, p, b [, loc, scale])
geninvgauss_logsf(x, p, b [, loc, scale])
geninvgauss_ppf(probability, p, b [, loc, scale])
geninvgauss_isf(probability, p, b [, loc, scale])
```

`p` may be any finite real value and `b` must be finite and positive. Standard
support is `z>0`. The density uses the real modified-Bessel-K normalization.
CDF/SF calculations integrate in the logarithmic coordinate `y=log(z)` and
accumulate the requested tail directly in log space. On a monotone deep tail,
the quadrature interval is localized from the endpoint log density before
integration; PPF/ISF compare logarithmic tails during safeguarded bisection, so
very small probabilities do not require forming a probability that underflows.

## Normal inverse Gaussian

```text
norminvgauss_pdf(x, a, b [, loc, scale])
norminvgauss_logpdf(x, a, b [, loc, scale])
norminvgauss_cdf(x, a, b [, loc, scale])
norminvgauss_sf(x, a, b [, loc, scale])
norminvgauss_logcdf(x, a, b [, loc, scale])
norminvgauss_logsf(x, a, b [, loc, scale])
norminvgauss_ppf(probability, a, b [, loc, scale])
norminvgauss_isf(probability, a, b [, loc, scale])
```

`a` must be finite and positive and `b` finite with `abs(b)<a`; support is the
real line. Density evaluation uses `K_1(a*sqrt(1+z**2))`. CDF/SF calculations
use the substitution `z=sinh(y)` and direct log-domain quadrature of the
requested tail. Tail intervals and panel widths adapt to the local transformed
density, and inverse functions bisect on log-CDF or log-SF for small requested
probabilities.

## Discrete distributions

Discrete families provide `pmf` and `logpmf` in place of `pdf` and `logpdf`.
The count may be an `integer` or a `real(dp)`; the two forms are separate
specific procedures behind one generic name. Following `scipy.stats`:

- a non-integer count has probability zero, and `cdf`, `sf`, and their
  logarithms use its floor;
- `ppf(p)` is the smallest count whose `cdf` reaches `p`, and `isf(p)` the
  smallest count whose `sf` falls to `p` or below;
- at a probability endpoint, `ppf(0)` and `isf(1)` return the integer just
  below the shifted support. `ppf(1)` returns the finite upper endpoint for
  Bernoulli and binomial and `+infinity` for the unbounded families.

### Bernoulli

```text
bernoulli_pmf(k, p [, loc])
bernoulli_logpmf, bernoulli_cdf, bernoulli_sf, bernoulli_logcdf, bernoulli_logsf
bernoulli_ppf(probability, p [, loc]), bernoulli_isf(probability, p [, loc])
```

The success probability must lie in `[0, 1]`. The distribution is exactly
the binomial distribution with one trial, and the implementation delegates to
that common kernel. Its shifted support is `{loc, loc + 1}`.

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

### Geometric

```text
geometric_pmf(k, p [, loc])
geometric_logpmf, geometric_cdf, geometric_sf, geometric_logcdf, geometric_logsf
geometric_ppf(probability, p [, loc]), geometric_isf(probability, p [, loc])
```

This is the number of Bernoulli trials through the first success, with
unshifted support `1, 2, ...` and success probability `0 < p <= 1`, matching
`scipy.stats.geom`. The survival function is evaluated directly as
`(1 - p)**floor(k - loc)`, using `log1p`- and `expm1`-based helpers for
stable logarithmic and complementary tails. The quantile has a closed-form
logarithmic starting value followed by an exact integer correction.

### Negative binomial

```text
negative_binomial_pmf(k, n, p [, loc])
negative_binomial_logpmf, negative_binomial_cdf, negative_binomial_sf
negative_binomial_logcdf, negative_binomial_logsf
negative_binomial_ppf(probability, n, p [, loc])
negative_binomial_isf(probability, n, p [, loc])
```

`k` counts failures before `n` successes. The shape `n` may be any finite
positive real value and `0 < p <= 1`, matching `scipy.stats.nbinom`. Its mass
is `Gamma(k + n) p**n (1 - p)**k / (Gamma(n) Gamma(k + 1))`. The tails use
`P(X <= k) = I_p(n, k + 1)` (DLMF 8.17.5), with CDF and survival probability
formed together by the incomplete-beta kernel rather than by subtraction.


### Discrete uniform (`randint`)

```text
randint_pmf(k, low, high [, loc])
randint_logpmf, randint_cdf, randint_sf, randint_logcdf, randint_logsf
randint_ppf(probability, low, high [, loc])
randint_isf(probability, low, high [, loc])
```

`low` and `high` must be finite integers with `low < high`. The unshifted
support is `low, ..., high - 1` and every support point has mass
`1/(high-low)`. A real-valued `loc` translates that lattice. CDF/SF routines
floor the shifted count, while PMF/log-PMF return zero/`-infinity` at
non-lattice points. As in SciPy, `ppf(0)` returns the point immediately below
the lower support and `ppf(1)` returns the finite upper support point.

### Planck

```text
planck_pmf(k, lambda [, loc])
planck_logpmf, planck_cdf, planck_sf, planck_logcdf, planck_logsf
planck_ppf(probability, lambda [, loc])
planck_isf(probability, lambda [, loc])
```

`lambda` must be finite and positive. On the unshifted support `k=0,1,...`,
`P(X=k)=(1-exp(-lambda))*exp(-lambda*k)`. CDF/SF and their logarithms are
formed directly with `expm1`/`log1p`-style helpers. This is exactly a geometric
distribution with success probability `1-exp(-lambda)`, shifted down by one.

### Discrete Laplace

```text
dlaplace_pmf(k, a [, loc])
dlaplace_logpmf, dlaplace_cdf, dlaplace_sf, dlaplace_logcdf, dlaplace_logsf
dlaplace_ppf(probability, a [, loc])
dlaplace_isf(probability, a [, loc])
```

`a` must be finite and positive. The unshifted support is every integer and
`P(X=k)=tanh(a/2)*exp(-a*abs(k))`. Lower and upper tails are evaluated from
closed forms on their respective half-lines rather than by complementing one
another. The distribution is symmetric around `loc`.

### Logarithmic series

```text
logser_pmf(k, p [, loc])
logser_logpmf, logser_cdf, logser_sf, logser_logcdf, logser_logsf
logser_ppf(probability, p [, loc])
logser_isf(probability, p [, loc])
```

`p` must satisfy `0 < p < 1`; the unshifted support is `1,2,...` and
`P(X=k)=-p**k/(k*log(1-p))`. The implementation sums the shorter lower tail
when appropriate and otherwise evaluates the upper tail directly by a
positive-term recurrence. PPF/ISF use exponential bracketing followed by an
integer binary search, so the returned value obeys the discrete SciPy
quantile convention without converting a continuous approximation to an
integer.


### Beta-negative-binomial

```text
betanbinom_pmf(k, n, a, b [, loc])
betanbinom_logpmf, betanbinom_cdf, betanbinom_sf
betanbinom_logcdf, betanbinom_logsf
betanbinom_ppf(probability, n, a, b [, loc])
betanbinom_isf(probability, n, a, b [, loc])
```

`n` is a positive integer success count and `a,b > 0` are beta shapes. The
unshifted support is `0,1,...`. The mass is evaluated with log-gamma/log-beta
identities and the lower tail is accumulated by a stable probability-ratio
recurrence. Quantiles use exponential bracketing followed by integer binary
search.

### Yule-Simon

```text
yulesimon_pmf(k, alpha [, loc])
yulesimon_logpmf, yulesimon_cdf, yulesimon_sf
yulesimon_logcdf, yulesimon_logsf
yulesimon_ppf(probability, alpha [, loc])
yulesimon_isf(probability, alpha [, loc])
```

`alpha` must be finite and positive, with unshifted support `1,2,...`.
`P(X=k)=alpha*B(k,alpha+1)` and the survival function is evaluated directly as
`k*B(k,alpha+1)`, avoiding subtraction in the upper tail.

### Zipf / zeta

```text
zipf_pmf(k, a [, loc])
zipf_logpmf, zipf_cdf, zipf_sf, zipf_logcdf, zipf_logsf
zipf_ppf(probability, a [, loc]), zipf_isf(probability, a [, loc])
```

`a > 1` and the unshifted support is `1,2,...`. The normalization is the
Riemann zeta value `zeta(a)=zeta(a,1)`. Upper tails are evaluated directly as
`zeta(a,k+1)/zeta(a)` using SciFort's real Euler-Maclaurin Hurwitz-zeta kernel;
CDF and log-CDF are formed from that direct tail with cancellation-resistant
complements.

### Finite Zipfian

```text
zipfian_pmf(k, a, n [, loc])
zipfian_logpmf, zipfian_cdf, zipfian_sf
zipfian_logcdf, zipfian_logsf
zipfian_ppf(probability, a, n [, loc])
zipfian_isf(probability, a, n [, loc])
```

`a >= 0`, `n` is a positive integer, and the unshifted support is `1,...,n`.
Masses are proportional to `k**(-a)` and normalized by the generalized
harmonic number. Lower and upper finite tails are accumulated independently;
`a=0` reduces exactly to the discrete uniform law on `1,...,n`.

### Fisher noncentral hypergeometric

```text
nchypergeom_fisher_pmf(k, M, n, N, odds [, loc])
nchypergeom_fisher_logpmf, nchypergeom_fisher_cdf, nchypergeom_fisher_sf
nchypergeom_fisher_logcdf, nchypergeom_fisher_logsf
nchypergeom_fisher_ppf(probability, M, n, N, odds [, loc])
nchypergeom_fisher_isf(probability, M, n, N, odds [, loc])
```

`M`, `n`, and `N` are integer-valued population/sample parameters with the
usual hypergeometric constraints, and `odds > 0`. The PMF is proportional to
`choose(n,k)*choose(M-n,N-k)*odds**k` over the finite hypergeometric support.
Normalization and both tails are accumulated with log-sum-exp arithmetic.
`odds=1` reduces exactly to the ordinary hypergeometric distribution. Only
`odds` has a continuous analytic score; the population parameters and `loc`
remain integral fit parameters.

### Wallenius noncentral hypergeometric

```text
nchypergeom_wallenius_pmf(k, M, n, N, odds [, loc])
nchypergeom_wallenius_logpmf, nchypergeom_wallenius_cdf
nchypergeom_wallenius_sf, nchypergeom_wallenius_logcdf
nchypergeom_wallenius_logsf
nchypergeom_wallenius_ppf(probability, M, n, N, odds [, loc])
nchypergeom_wallenius_isf(probability, M, n, N, odds [, loc])
```

`M`, `n`, and `N` are integer-valued with the ordinary hypergeometric
constraints, and `odds > 0`. The unshifted support is
`max(0,N-(M-n)),...,min(N,n)`. SciFort evaluates the law by exact finite-state
dynamic programming over the defining sequential biased-urn process and sums
each requested tail directly. `odds=1` reduces to the ordinary hypergeometric
distribution. The analytic score contains the single continuous `odds`
component; `M`, `n`, `N`, and `loc` are integral fit parameters. RNG follows
the same sequential biased-urn experiment exactly.

### Poisson-binomial

```text
poisson_binom_pmf(k, p(:) [, loc])
poisson_binom_logpmf(k, p(:) [, loc])
poisson_binom_cdf(k, p(:) [, loc])
poisson_binom_sf(k, p(:) [, loc])
poisson_binom_logcdf(k, p(:) [, loc])
poisson_binom_logsf(k, p(:) [, loc])
poisson_binom_ppf(probability, p(:) [, loc])
poisson_binom_isf(probability, p(:) [, loc])
```

Each `p(i)` must lie in `[0,1]`; support is `0,...,size(p)` before `loc` is
applied. Because the shape is a rank-one vector, these evaluators are pure
scalar procedures rather than elemental functions. PMF and tails use dynamic
programming for the convolution of independent Bernoulli variables. The score
is a vector of length `size(p)` computed from leave-one-out convolutions. The
fit parameter vector is `[p(:), loc]`; permutations of `p` are statistically
equivalent, so an unrestricted fit is nonidentifiable up to permutation. RNG
returns the exact sum of independent Bernoulli draws. The C ABI accepts a
`const double *p` and `size_t n` and returns scalar PMF/CDF/PPF values.

### Skellam

```text
skellam_pmf(k, mu1, mu2 [, loc])
skellam_logpmf(k, mu1, mu2 [, loc])
skellam_cdf(k, mu1, mu2 [, loc])
skellam_sf(k, mu1, mu2 [, loc])
skellam_logcdf(k, mu1, mu2 [, loc])
skellam_logsf(k, mu1, mu2 [, loc])
skellam_ppf(probability, mu1, mu2 [, loc])
skellam_isf(probability, mu1, mu2 [, loc])
```

`mu1` and `mu2` must be finite and strictly positive. The distribution is the
difference of two independent Poisson counts. The mass uses a scaled integer-
order modified-Bessel-I recurrence; lower and upper tails use equivalent
noncentral-chi-square identities directly, avoiding subtraction. Real-valued
counts follow the same shifted-lattice floor convention as the other discrete
families.

## Random generation

Random-number state lives in a separate module so stochastic routines never
hide process-global mutable state:

```fortran
use scifort_random, only : rng_state, rng_seed, rng_uniform, &
    rng_fill_uniform, rng_get_state, rng_set_state
use scifort_stats, only : normal_rvs, normal_rvs_array, gamma_rvs, &
    poisson_rvs, binomial_rvs
```

`rng_state` defaults to the six-component MRG32k3a seed vector with every
component equal to `12345`. `rng_seed(state, seed)` deterministically maps any
default-integer or `int64` scalar seed onto a valid repeated state component;
positive seeds in the valid component range are preserved exactly. Use
`rng_get_state` and `rng_set_state` to checkpoint and restore the six `int64`
recurrence components. `rng_set_state` reports `rng_status_ok` or
`rng_status_invalid_state` and leaves the existing state unchanged on error.

`rng_uniform(state)` returns an open-interval binary64 value. Each 26-bit
chunk is formed by rejection sampling a range whose size is an exact multiple
of `2**26`, partitioning that range into equal consecutive buckets. Two chunks
are combined into a 52-bit integer and shifted by half a binary64 unit on that
grid, so neither zero nor one is returned. `rng_fill_uniform` applies the same
operation to a caller-provided array.

Every current distribution has:

```text
<family>_rvs(state, parameters... [, loc] [, scale])
<family>_rvs_array(state, samples, parameters... [, loc] [, scale])
```

Continuous samplers accept the same shape, `loc`, and `scale` parameters as
their PPF. Discrete samplers accept the same shape/probability parameters and
optional `loc`. The implementation is currently inverse-transform sampling:
each valid draw evaluates the family's PPF at `rng_uniform(state)`. Invalid
parameters return NaN (or fill the array with NaN) without advancing the RNG
state. This makes sampling conventions exactly track the public quantile API,
but it is a correctness-first implementation and can be substantially slower
than specialized rejection or transformation algorithms.

## Likelihood evaluation

For independent observations, `scifort_likelihood` provides one summed
log-likelihood and one negative log-likelihood for every current family. The
same names are re-exported by `scifort_stats`:

```fortran
use scifort_stats, only : normal_loglikelihood, normal_nnlf, &
    gamma_loglikelihood, poisson_loglikelihood
```

The naming pattern is:

```text
<family>_loglikelihood(data, parameters... [, loc] [, scale])
<family>_nnlf(data, parameters... [, loc] [, scale])
```

`*_loglikelihood` is the sum of the family's `logpdf` or `logpmf` over the
input array. `*_nnlf` is its negative, matching the negative-log-likelihood
objective convention used by SciPy's distribution fitting infrastructure.
The routines are `pure`, allocate no temporary state, and return zero for an
empty observation array. Invalid distribution parameters propagate NaN;
observations outside a distribution's support contribute negative infinity.

For Bernoulli, Poisson, geometric, binomial, and negative-binomial families,
the public likelihood names are generic over integer and `real(dp)` data
arrays. The binomial integer-data overload takes an integer trial count, just
as the underlying integer-count PMF interface does.

These functions assume independent, uncensored, equally weighted
observations. Censoring and observation weights are not yet implemented.

## Score functions

`scifort_score` provides analytic gradients of the summed log-likelihood for
all current families, and `scifort_stats` re-exports the same names:

```fortran
use scifort_stats, only : normal_score, gamma_score, beta_score, &
    poisson_score, negative_binomial_score
```

The naming pattern is `<family>_score(data, parameters... [, loc] [, scale])`.
Continuous-family result arrays follow the same canonical order as fitting:
shape parameters first, then `loc`, then `scale`. Families without shape
parameters return `[loc, scale]`. The shape-family score orders are
`[a, loc, scale]` for gamma, `[df, loc, scale]` for chi, chi-square, and Student t,
`[s, loc, scale]` for lognormal, `[c, loc, scale]` for Weibull, `[b, loc,
scale]` for Pareto, `[a, b, loc, scale]` for beta, `[dfn, dfd, loc, scale]`
for F, `[a, loc, scale]` for power-function, `[c, loc, scale]` for triangular
and generalized Pareto, Lomax, log-Laplace, Bradford, Fisk, double-Weibull,
fatigue-life, generalized logistic, power-normal, log-gamma, Gompertz,
inverse-Weibull, generalized half-logistic, exponential-power, Weibull maximum, or R-distribution,
`[beta, loc, scale]` for generalized normal, `[nu, loc, scale]`
for Nakagami, `[b, loc, scale]` for truncated exponential, `[a, loc, scale]`
for Alpha, inverse gamma, or double-gamma, `[kappa, loc, scale]` for asymmetric Laplace,
`[a, b, loc, scale]` for beta-prime, truncated normal, or log-uniform, `[b, c, loc, scale]` for
truncated Pareto, `[a, c, loc, scale]` for exponentiated-Weibull, `[c, s, loc, scale]` for power-lognormal,
`[c, d, loc, scale]` for Burr XII, `[a, loc, scale]` for skew-Cauchy, `[c, loc, scale]` for folded normal or folded Cauchy,
`[mu, loc, scale]` for inverse Gaussian or reciprocal inverse-Gaussian,
and `[loc, scale]` for both Gumbel,
arcsine, half-normal, half-Cauchy, Maxwell, cosine, semicircular, Anglit, Moyal, Landau,
hyperbolic-secant, half-logistic, Levy, left-Levy, and Wald families.

Discrete scores include differentiable distribution parameters only. Bernoulli,
Poisson, geometric, and binomial return `[p]`, `[mu]`, `[p]`, and `[p]`,
respectively. Negative binomial returns `[n, p]` for the continuous positive
extension of `n`. Wallenius returns `[odds]`; Poisson-binomial returns one score component for
each entry of `p(:)`. Integer support shifts and integral population/trial
parameters are not differentiated. Integer and `real(dp)` observation arrays
are both accepted.

For distributions whose support moves with a fitted parameter, an observation
exactly on a support boundary can make the ordinary two-sided likelihood
derivative nonregular. Affected score components return NaN at those boundary
cases rather than returning a one-sided derivative.

## Bounded maximum-likelihood fitting

`scifort_fit` provides one `<family>_fit` routine for every current family and
re-exports them through `scifort_stats`. Each routine minimizes the existing
`*_nnlf` objective over explicit finite bounds:

```fortran
use scifort_stats, only : fit_result, normal_fit, gamma_fit, binomial_fit

type(fit_result) :: result
real(dp) :: lower(2), upper(2)

lower = [-10.0_dp, 0.01_dp]
upper = [ 10.0_dp, 10.0_dp]
call normal_fit(data, lower, upper, result)
```

All current fitting bounds must be finite. Setting `lower(i) == upper(i)` fixes
parameter `i`. An optional `guess` array supplies the starting point; otherwise
the midpoint of each effective bound is used. Optional `max_iter` and
`tolerance` arguments control the deterministic projected coordinate/pattern
search. The optimizer is local rather than globally guaranteed, so bounds and
a scientifically reasonable starting guess remain important for multimodal or
weakly identified problems.

The canonical fit parameter vectors are:

| Family | Parameter vector |
| --- | --- |
| normal, uniform, exponential, Laplace, logistic, Cauchy, Rayleigh, Gumbel right/left, arcsine, half-normal, half-Cauchy, Maxwell, cosine, semicircular, Anglit, Moyal, Landau, hyperbolic-secant, half-logistic, Levy, left-Levy, Wald | `[loc, scale]` |
| gamma, inverse gamma, Alpha, double-gamma | `[a, loc, scale]` |
| inverse Gaussian | `[mu, loc, scale]` |
| asymmetric Laplace | `[kappa, loc, scale]` |
| chi, chi-square, Student t | `[df, loc, scale]` |
| lognormal | `[s, loc, scale]` |
| Weibull | `[c, loc, scale]` |
| Pareto | `[b, loc, scale]` |
| power-function | `[a, loc, scale]` |
| triangular, generalized Pareto, Lomax, log-Laplace, Bradford, Fisk, double-Weibull, fatigue-life, generalized logistic, Gompertz, inverse-Weibull, generalized half-logistic, Weibull maximum, R-distribution | `[c, loc, scale]` |
| exponential-power | `[b, loc, scale]` |
| truncated exponential | `[b, loc, scale]` |
| generalized normal | `[beta, loc, scale]` |
| Nakagami | `[nu, loc, scale]` |
| power-normal, log-gamma | `[c, loc, scale]` |
| generalized extreme-value, wrapped Cauchy | `[c, loc, scale]` |
| kappa-3 | `[a, loc, scale]` |
| Burr Type III | `[c, d, loc, scale]` |
| Mielke beta-kappa / Dagum | `[k, s, loc, scale]` |
| kappa-4 | `[h, k, loc, scale]` |
| doubly truncated Weibull minimum | `[c, a, b, loc, scale]` |
| Gibrat | `[loc, scale]` |
| skew-Cauchy | `[a, loc, scale]` |
| beta, beta-prime, truncated normal, log-uniform | `[a, b, loc, scale]` |
| exponentiated-Weibull | `[a, c, loc, scale]` |
| power-lognormal | `[c, s, loc, scale]` |
| Burr XII | `[c, d, loc, scale]` |
| F | `[dfn, dfd, loc, scale]` |
| Bernoulli | `[p, loc]` |
| Poisson | `[mu, loc]` |
| geometric | `[p, loc]` |
| binomial | `[n, p, loc]` |
| negative binomial | `[n, p, loc]` |
| double-Pareto lognormal | `[u, s, a, b, loc, scale]` |
| circular von Mises, von Mises on a line | `[kappa, loc, scale]` |
| asymptotic two-sided Kolmogorov | `[loc, scale]` |
| finite-sample two-sided Kolmogorov-Smirnov | `[n, loc, scale]` |
| noncentral Student t | `[df, nc, loc, scale]` |
| Gauss hypergeometric | `[a, b, c, z, loc, scale]` |
| generalized inverse Gaussian | `[p, b, loc, scale]` |
| normal inverse Gaussian | `[a, b, loc, scale]` |
| generalized hyperbolic | `[p, a, b, loc, scale]` |
| Fisher noncentral hypergeometric | `[M, n, N, odds, loc]` |
| Wallenius noncentral hypergeometric | `[M, n, N, odds, loc]` |
| Poisson-binomial | `[p(:), loc]` |
| Skellam | `[mu1, mu2, loc]` |

The fitting layer constrains Bernoulli/Poisson/geometric/Skellam `loc`,
binomial `n` and `loc`, negative-binomial `n` and `loc`, and Wallenius
`M`, `n`, `N`, and `loc`, plus finite-sample two-sided KS `n`, to integer
values. Poisson-binomial constrains only
`loc` to an integer while each probability remains continuous in `[0,1]`. This does not
change the lower-level negative-binomial density API, which continues to allow
a positive real `n` for evaluation and scoring.

`fit_result%params` contains the final parameter vector, `nllf` the final
negative log-likelihood, `iterations` and `evaluations` optimizer diagnostics,
`status` one of `fit_status_success`, `fit_status_max_iter`,
`fit_status_invalid_input`, or `fit_status_no_finite_objective`, and `success`
a convenience logical flag. Score and fitting interfaces are currently native
Fortran only.

## Descriptive statistics

`scifort_descriptive` is re-exported by `scifort_stats` and currently provides:

```text
mean(x)
variance(x [, ddof])
standard_deviation(x [, ddof])
central_moment(x, order)
quantile(x, q)
median(x)
covariance(x, y [, ddof])
pearson_correlation(x, y)
rankdata(x [, method])
```

`variance`, `standard_deviation`, and `covariance` default to `ddof = 1`;
`ddof = 0` gives the population divisor. Mean and second-order reductions use
scaled/compensated calculations so large finite observations are not summed or
squared naively. Variance/covariance-style reductions require finite data and
return NaN when their divisor is nonpositive. Empty reductions return NaN.

`quantile` implements linear interpolation at position `(n - 1) q`, equivalent
to the common Hyndman-Fan type 7 / NumPy `method="linear"` convention. `q`
may be a scalar or an array. NaN observations or probabilities outside
`[0, 1]` return NaN. `median(x)` is `quantile(x, 0.5)`.

`rankdata` returns `real(dp)` ranks beginning at one. The default method is
`average`; `min`, `max`, `dense`, and stable `ordinal` are also supported.
NaN input propagates to every returned rank. Infinities are ordered normally.

## Classical hypothesis tests

`scifort_hypothesis` and `scifort_stats` expose:

```text
ttest_1samp(x, popmean [, alternative])
ttest_ind(x, y [, equal_var] [, alternative])
ttest_rel(x, y [, alternative])
pearsonr(x, y [, alternative])
spearmanr(x, y [, alternative])
mannwhitneyu(x, y [, alternative] [, use_continuity])
```

`alternative` is `two-sided` by default and may also be `less` or `greater`.
The t-test routines return `ttest_result` with `statistic`, `pvalue`, and `df`.
`ttest_ind` uses the pooled-variance test by default; `equal_var=.false.` selects
Welch's test and the Welch-Satterthwaite degrees of freedom.

`pearsonr` returns the ordinary sample correlation and the exact beta-law
p-value under the independent bivariate-normal null, with the `n = 2` p-value
special case handled explicitly. `spearmanr` uses average ranks and the usual
Student-t large-sample transformation. Both return `correlation_test_result`.

`mannwhitneyu` returns `mannwhitneyu_result`. The reported statistic is the U
statistic for the first sample. Its current p-value implementation is the
normal asymptotic approximation with tie correction and, by default, a 0.5
continuity correction. Exact small-sample Mann-Whitney p-values are not yet
implemented.

All current hypothesis routines propagate invalid/non-finite reductions as
NaN rather than issuing warnings or stopping the program.

### Extended hypothesis and contingency tests

`scifort_hypothesis_extended` and `scifort_stats` additionally expose:

```text
make_sample_group(values)
ranksums(x, y [, alternative])
kruskal(groups)
friedmanchisquare(samples)
wilcoxon(x [, y, zero_method, correction, alternative, method])
power_divergence(f_obs [, f_exp, ddof, lambda_])
power_divergence_named(f_obs [, f_exp, ddof, lambda_name])
chisquare(f_obs [, f_exp, ddof])
fisher_exact(table [, alternative])
contingency_expected_freq(observed)
chi2_contingency(observed [, correction, lambda_])
f_oneway(groups [, equal_var])
bartlett(groups)
levene(groups [, center, proportiontocut])
fligner(groups [, center, proportiontocut])
```

`sample_group` stores one independent sample, allowing unequal group sizes in
Kruskal-Wallis, one-way ANOVA, Bartlett, Levene, and Fligner-Killeen calls.
`friedmanchisquare` instead accepts a rectangular array with samples/conditions
by row and matched blocks/subjects by column. `ranksums` is the large-sample
Wilcoxon rank-sum test without tie correction; `wilcoxon` is the paired or
one-sample signed-rank test and supports `wilcox`, `pratt`, and `zsplit` zero
methods plus `auto`, `exact`, `asymptotic`, and small-sample exhaustive
`permutation` methods.

`power_divergence_named` accepts `pearson`, `log-likelihood`, `freeman-tukey`,
`mod-log-likelihood`, `neyman`, and `cressie-read`. `chisquare` is the Pearson
special case. `fisher_exact` currently implements the exact 2-by-2 fixed-margin
hypergeometric test for all three alternatives. `chi2_contingency` operates on
two-dimensional tables and applies Yates' correction by default when the table
has one degree of freedom.

`f_oneway(..., equal_var=.false.)` selects Welch's heteroscedastic one-way
ANOVA. `levene` and `fligner` accept `median`, `mean`, or `trimmed` centers; the
default is `median`. As elsewhere in the Fortran API, invalid argument shapes or
non-finite inputs are represented by NaN/status-like result values rather than
Python exceptions or warnings.

## Exact contingency tests, effect measures, and p-value utilities

`scifort_contingency_meta` and `scifort_stats` expose:

```text
binomtest(k, n [, p, alternative])
result%proportion_ci([confidence_level, method])
barnard_exact(table [, alternative, pooled, n])
boschloo_exact(table [, alternative, n])
odds_ratio(table [, kind])
result%confidence_interval([confidence_level, alternative])
relative_risk(exposed_cases, exposed_total, control_cases, control_total)
result%confidence_interval([confidence_level])
association(observed [, method, correction])
margins(observed, row_sums, column_sums)
combine_pvalues(pvalues [, method, weights])
false_discovery_control(ps [, method])
```

`binomtest` implements exact binomial p-values for all three alternatives. Its
result object provides `exact`, `wilson`, and `wilsoncc` proportion confidence
intervals. `barnard_exact` supports pooled and unpooled Wald statistics;
`boschloo_exact` uses Fisher's one-sided p-value as its statistic. Their nuisance
success probability is maximized deterministically in one dimension with a
dense global grid followed by local golden-section refinement. This replaces
SciPy's SHGO/Sobol optimizer while preserving the same test statistic and
rejection-region definitions.

`odds_ratio(..., kind='conditional')` returns the conditional maximum-likelihood
estimate under Fisher's noncentral hypergeometric model and exact conditional
confidence intervals. `kind='sample'` returns the ordinary cross-product ratio.
`relative_risk` returns the sample risk ratio and a Katz log-scale confidence
interval.

`association` supports `cramer`, `tschuprow`, and `pearson`. `margins` currently
accepts two-dimensional integer or real tables and returns row and column sums;
SciPy's arbitrary-rank margin helper is not yet mirrored. `combine_pvalues`
currently operates on one-dimensional inputs and supports `fisher`, `pearson`,
`mudholkar_george`, `tippett`, and `stouffer`; Stouffer weights are optional.
`false_discovery_control` is currently one-dimensional and supports `bh` and
`by`. Axis/nan-policy wrappers remain future API work.

## Multiple comparisons and heterogeneous-mean tests

`scifort_multiple_comparisons` and `scifort_stats` expose:

```text
alexandergovern(groups)
tukey_hsd(groups [, equal_var])
result%confidence_interval([confidence_level])
dunnett(groups, control [, alternative, maxpts])
result%confidence_interval([confidence_level, maxpts])
poisson_means_test(k1, n1, k2, n2 [, diff, alternative])
```

`alexandergovern` implements the Alexander-Govern chi-square approximation for
comparing two or more means without assuming equal variances. `tukey_hsd` uses
Tukey HSD for equal-sized homoscedastic groups, Tukey-Kramer for unequal group
sizes when `equal_var=.true.`, and Games-Howell when `equal_var=.false.`. The
returned pairwise statistic matrix is the matrix of mean differences; adjusted
p-values and simultaneous confidence intervals reuse SciFort's translated
studentized-range distribution.

`dunnett` compares each treatment sample with a common control using the pooled
within-group variance and Dunnett correlation matrix. Adjusted p-values and
simultaneous confidence intervals are evaluated with SciFort's multivariate
Student-t rectangular-probability integrator. `maxpts` controls that numerical
integration budget; unlike SciPy, SciFort's default integration stream is
deterministic rather than NumPy-RNG driven, so last digits can differ within
the integration tolerance.

`poisson_means_test` implements the Krishnamoorthy-Thomson E-test for two
Poisson rates, including nonzero nonnegative null differences and all three
alternatives. The infinite double sum is truncated at the same extreme Poisson
quantiles used by SciPy 1.17.0. Invalid scalar inputs return NaN-like result
values instead of raising Python exceptions.

## Additional nonparametric and multi-sample tests

`scifort_nonparametric_extended` and `scifort_stats` expose:

```text
ansari(x, y [, alternative])
mood(x, y [, alternative])
epps_singleton_2samp(x, y [, t])
anderson_ksamp(groups [, variant])
median_test(groups [, ties, correction, lambda_])
median_test_named(groups [, ties, correction, lambda_name])
```

`ansari` implements the Ansari-Bradley scale test. Untied samples with both
sample sizes below 55 use an exact finite distribution computed by dynamic
programming; otherwise the tie-corrected normal approximation is used. `mood`
implements Mood's two-sample scale test with the Mielke tie correction. Both
support `two-sided`, `less`, and `greater` alternatives.

`epps_singleton_2samp` implements the characteristic-function two-sample test.
The optional positive, distinct `t` vector replaces SciPy's default two-point
frequency grid. The covariance pseudoinverse is formed with SciFort's symmetric
eigensolver, so tiny last-bit differences from SciPy's LAPACK-backed path are
possible for ill-conditioned or larger custom grids.

`anderson_ksamp` accepts two or more `sample_group` values and supports
`midrank`, `right`, and `continuous` statistics. It returns the normalized
statistic, the seven standard critical values, and SciPy-style asymptotic
p-value interpolation with 0.25/0.001 cap/floor behavior. A permutation-method
argument is not currently provided.

`median_test` returns the power-divergence statistic, p-value, grand median,
and the 2-by-k contingency table (above/below). Median ties can be placed
`below`, `above`, or `ignore`d. `median_test_named` accepts `pearson`,
`log-likelihood`, `freeman-tukey`, `mod-log-likelihood`, `neyman`, and
`cressie-read`; the numeric form accepts an explicit Cressie-Read lambda.
Invalid inputs return NaN/status-like results rather than Python exceptions.

## Association, robust regression, and ordered trend tests

`scifort_association_extended` and `scifort_stats` expose:

```text
pointbiserialr(x, y)
kendalltau(x, y [, method, variant, alternative])
linregress(x, y [, alternative])
theilslopes(y, x [, alpha, method])
siegelslopes(y, x [, method])
brunnermunzel(x, y [, alternative, distribution])
page_trend_test(data [, ranked, predicted_ranks, method])
```

`kendalltau` supports tau-b and tau-c, all three alternatives, exact untied
finite-sample probabilities, and the tie-corrected asymptotic approximation.
`linregress` reports slope, intercept, Pearson r, p-value, slope standard error,
and intercept standard error. `theilslopes` supports `separate` and `joint`
intercepts and returns a confidence interval for the median slope;
`siegelslopes` supports `hierarchical` and `separate` intercepts.

`brunnermunzel` supports Student-t and normal reference distributions and all
three alternatives. `page_trend_test` ranks within rows unless `ranked=.true.`,
accepts an explicit permutation of predicted ranks, and provides exact or
asymptotic one-sided p-values. The exact Page implementation is limited to at
most ten conditions; `auto` selects exact only within the smaller SciPy-style
region, so ordinary automatic calls do not encounter that implementation cap.
Invalid inputs return NaN/status-like results rather than Python exceptions.

## Goodness-of-fit and normality tests

`scifort_goodness_of_fit` and `scifort_stats` expose:

```text
skewtest(x [, alternative])
kurtosistest(x [, alternative])
normaltest(x)
jarque_bera(x)
ks_1samp(x, cdf [, alternative, method])
ks_2samp(x, y [, alternative, method])
kstest(...)
cramervonmises(x, cdf)
cramervonmises_2samp(x, y [, method])
shapiro(x)
anderson(x [, dist, interpolate_pvalue])
```

The KS result includes `statistic_location` and `statistic_sign`, matching the
corresponding SciPy result fields. One-sample KS supports `auto`, `exact`,
`approx`, and `asymp`; two-sample KS supports `auto`, `exact`, and `asymp`.
For bounded Fortran work/memory, exact two-sample KS is limited to sample-size
products at most 25,000,000; larger explicit-exact requests use the asymptotic
path. Two-sample Cramer-von Mises uses exact enumeration while the number of
label combinations is at most 2,000,000 and otherwise uses the asymptotic
formula.

`shapiro_result` contains the W statistic, Royston p-value approximation, and a
status code. `anderson` currently supports `norm`, `expon`, `logistic`,
`gumbel_l` (also `gumbel`/`extreme1`), and `gumbel_r`; SciPy's `weibull_min`
Anderson-Darling fitting path is not yet implemented. Legacy critical values and
fitted parameters are returned; optional `interpolate_pvalue=.true.` provides a
table-interpolated p-value clipped to the tabulated range.

## Resampling

`scifort_resampling` provides explicit-state scalar-statistic interfaces:

```fortran
result = bootstrap(x, statistic, state [, n_resamples] [, confidence_level] &
    [, alternative] [, method])
result = bootstrap(x, y, statistic, state [, n_resamples] [, confidence_level] &
    [, alternative] [, method] [, paired])
result = bootstrap_percentile(x, statistic, state [, n_resamples] [, confidence_level])
result = permutation_test(x, y, statistic, state [, n_resamples] [, alternative] &
    [, permutation_type])
result = monte_carlo_test(x, rvs, statistic, state [, n_resamples] [, alternative])
result = monte_carlo_test(x, y, rvs_x, rvs_y, statistic, state &
    [, n_resamples] [, alternative])
result = power(rvs, test, n_observations, state [, significance] [, n_resamples])
result = power(rvs_x, rvs_y, test, n_observations_x, n_observations_y, state &
    [, significance] [, n_resamples])
```

`bootstrap` defaults to 9999 draws, 95% confidence, a two-sided interval, and
`method='bca'`. Methods `percentile`, `basic`, and `bca` are supported.
One-sided intervals follow SciPy's convention: the finite endpoint of a
`confidence_level` one-sided interval equals the corresponding endpoint of a
two-sided interval whose confidence level is twice as far from one. For two
samples, `paired=.true.` resamples common indices; otherwise each sample is
resampled independently. `bootstrap_percentile` is retained as a compatibility
wrapper for the original one-sample percentile API.

`bootstrap_result` contains the observed statistic, sample standard deviation
of the bootstrap distribution, confidence endpoints, number of draws, and the
full bootstrap distribution. BCa acceleration is computed by leave-one-out
jackknife influence values; degenerate BCa problems return NaN endpoints.

`permutation_test` accepts `permutation_type='independent'`, `'samples'`, or
`'pairings'`. If the requested number of resamples is at least the number of
distinct permutations, SciFort enumerates the exact null distribution and does
not advance the RNG state. Otherwise it samples random permutations and applies
the conservative plus-one correction `(b + 1)/(B + 1)`. The result records the
observed statistic, p-value, null distribution, draw count, and whether the test
was exact. Floating-point comparisons use the same `100*epsilon` neighborhood
around the observed statistic used by SciPy 1.17.

`monte_carlo_test` accepts one or two null-sample generator callbacks. It stores
the simulated null distribution and uses the plus-one p-value correction for
`less`, `greater`, and twice-the-smaller-tail `two-sided` alternatives.

`power` accepts one or two alternative-sample generator callbacks and a test
callback that returns a p-value. It stores all simulated p-values and reports
the fraction strictly below `significance` (default 0.01), matching SciPy's
comparison convention.

All resampling APIs use caller-owned `rng_state` objects. Invalid arguments are
validated before random sampling and therefore do not advance the state.
Bootstrap and random-permutation indices are selected by rejection from the
52-bit uniform grid, avoiding modulo bias when the number of choices does not
divide `2**52`.

## C ABI distribution slice

The interoperable declarations are in `include/scifort.h`. Every currently
implemented distribution has scalar C entry points for its density or mass,
CDF, and PPF. The naming pattern is:

```text
scifort_<family>_pdf_f64(...)   continuous density
scifort_<family>_pmf_f64(...)   discrete mass
scifort_<family>_cdf_f64(...)
scifort_<family>_ppf_f64(...)
```

All scalar inputs and outputs use C `double`. Continuous entry points take
explicit `loc` and `scale` arguments after any required shape parameters.
Discrete entry points take explicit `loc`; count-valued arguments are passed
as `double` so the ABI can preserve the native real-count behavior. Families
that mathematically require an integer count or number of trials still reject
noninteger values according to the native Fortran rules.

Invalid scalar distribution parameters or probabilities return NaN exactly as
in the Fortran API. The normal distribution additionally provides bulk
`*_vec_f64` PDF, CDF, and PPF routines with explicit lengths, caller-owned
output arrays, and an integer status code for invalid `loc` or `scale`.

The C ABI does not yet expose survival, logarithmic, inverse-survival, or
special-function routines.

## Multivariate distributions

The first multivariate slice uses vector and matrix arguments and is therefore
not elemental. Samples are stored by column in the Fortran array APIs.

### Multivariate normal

```text
multivariate_normal_pdf(x [, mean, cov, allow_singular])
multivariate_normal_logpdf(x [, mean, cov, allow_singular])
multivariate_normal_cdf(x [, mean, cov, allow_singular, maxpts, abseps, releps, lower_limit, state])
multivariate_normal_logcdf(x [, mean, cov, allow_singular, maxpts, abseps, releps, lower_limit, state])
multivariate_normal_marginal(dimensions, mean, cov, marginal_mean, marginal_cov, status)
multivariate_normal_entropy([mean, cov])
multivariate_normal_rvs(state, sample [, mean, cov, status])
multivariate_normal_rvs_array(state, samples [, mean, cov, status])
multivariate_normal_fit(data, mean, cov, status [, fixed_mean, fixed_cov])
```

Absent `mean` and `cov` mean a zero vector and identity covariance. Only the
lower triangle of `cov` is authoritative, matching SciPy's ordinary dense
covariance path. Numerical rank uses the binary64 cutoff
`1e6*epsilon*max(abs(eigenvalues))`; singular covariance matrices require
`allow_singular=.true.`. Density is zero away from the affine support of a
singular covariance.

`cdf` evaluates rectangular probabilities with a Genz-style conditional
transform and randomized tent-transformed Halton integration. For scalar or
exactly diagonal covariance, `logcdf` instead sums exact one-dimensional normal
log-interval probabilities, using the smaller tail and `expm1` for interval
differences. This preserves finite log probabilities after ordinary probabilities
underflow, including shifted/scaled coordinates; this exact path does not advance
an optional RNG state. For correlated covariance, `logcdf` retains the existing
probability-domain integration and takes its logarithm: extremely small box
probabilities can still underflow to zero and yield negative infinity. General
correlated extreme-tail log integration remains a limitation. The
default lower limit is negative infinity componentwise. If a lower endpoint is
greater than its upper endpoint, that pair is swapped and the overall signed
box probability is adjusted as in SciPy. The default `maxpts` is approximately
`1000000*dimension` and the default absolute target is `1e-5`; `releps` is
accepted for SciPy-compatible call structure but is informational in this
implementation. Supplying `state` advances an explicit RNG state for randomized
shifts; omitting it uses deterministic local seeds, which makes regression runs
reproducible. A negative signed box probability has a complex logarithm in
SciPy; because SciFort's `logcdf` result is real, that case returns NaN.

`multivariate_normal_marginal` takes one-based, unique variable indices and
returns the corresponding mean vector and covariance submatrix in the requested
order. `fit` expects observations by column,
`data(dimension,n_observations)`, and uses the maximum-likelihood covariance
divisor `n_observations`.

### Multivariate Student t

```text
multivariate_t_pdf(x [, loc, shape, df, allow_singular])
multivariate_t_logpdf(x [, loc, shape, df])
multivariate_t_cdf(x [, loc, shape, df, allow_singular, maxpts, lower_limit, state])
multivariate_t_entropy([loc, shape, df])
multivariate_t_rvs(state, sample [, loc, shape, df, status])
multivariate_t_rvs_array(state, samples [, loc, shape, df, status])
multivariate_t_marginal(dimensions, loc, shape, marginal_loc, marginal_shape, status)
```

Absent `loc` and `shape` mean a zero vector and identity shape matrix; absent
`df` means one degree of freedom. The `shape` matrix is not generally the
covariance: for `df > 2`, covariance is `df/(df-2)*shape`. Only the lower
triangle of `shape` is authoritative. `pdf` rejects a singular shape by default
and accepts it with `allow_singular=.true.`; `logpdf` follows SciPy's public
semantics and evaluates with the pseudo-determinant/pseudo-inverse path. Infinite
`df` reduces to the multivariate normal.

`cdf` accepts a finite `lower_limit`, an optional explicit RNG `state`, and a
default point budget of approximately `1000*dimension`. It uses the same
conditional randomized low-discrepancy integration foundation as the normal
CDF, with the Student-t radial transformation. Marginal indices are one-based
and must be unique. Random draws use a multivariate normal vector divided by an
independent chi-square scale.

### Dirichlet

```text
dirichlet_pdf(x, alpha)
dirichlet_logpdf(x, alpha)
dirichlet_mean(alpha)
dirichlet_var(alpha)
dirichlet_cov(alpha)
dirichlet_entropy(alpha)
dirichlet_rvs(state, sample, alpha [, status])
dirichlet_rvs_array(state, samples, alpha [, status])
```

`alpha` must be positive. As in SciPy, `x` may contain either all `K` simplex
coordinates or the first `K-1`, in which case the final coordinate is inferred.
The simplex sum tolerance is `1e-9`.

### Multinomial

```text
multinomial_pmf(x, n, p)
multinomial_logpmf(x, n, p)
multinomial_mean(n, p)
multinomial_cov(n, p)
multinomial_entropy(n, p)
multinomial_rvs(state, sample, n, p [, status])
multinomial_rvs_array(state, samples, n, p [, status])
```

`n` is a nonnegative integer. For SciPy 1.17 compatibility, if `p` does not sum
to one within `10*epsilon`, the final probability is replaced by the residual
`1-sum(p(1:K-1))` before validation. This behavior is version-specific: SciPy
1.17 warns that it is scheduled to change in SciPy 1.18.

### Dirichlet-multinomial

```text
dirichlet_multinomial_pmf(x, alpha, n)
dirichlet_multinomial_logpmf(x, alpha, n)
dirichlet_multinomial_mean(alpha, n)
dirichlet_multinomial_var(alpha, n)
dirichlet_multinomial_cov(alpha, n)
```

The concentration vector must be positive and counts must be nonnegative
integers. A valid count vector whose total differs from `n` has zero probability.

### Multivariate hypergeometric

```text
multivariate_hypergeom_pmf(x, m, n)
multivariate_hypergeom_logpmf(x, m, n)
multivariate_hypergeom_mean(m, n)
multivariate_hypergeom_var(m, n)
multivariate_hypergeom_cov(m, n)
multivariate_hypergeom_rvs(state, m, n, sample [, status])
multivariate_hypergeom_rvs_array(state, m, n, samples [, status])
```

`m` contains nonnegative category population counts and `n` is the total sample
size without replacement. A valid `x` must contain nonnegative counts bounded by
`m` and sum to `n`. Valid zero- and one-object populations have zero finite-sample
variance/covariance, matching SciPy 1.17 behavior. Array random variates are stored
by column.

### Normal-inverse-gamma

```text
normal_inverse_gamma_pdf(x, s2 [, mu, lmbda, a, b])
normal_inverse_gamma_logpdf(x, s2 [, mu, lmbda, a, b])
normal_inverse_gamma_mean([mu, lmbda, a, b], mean_x, mean_s2)
normal_inverse_gamma_var([mu, lmbda, a, b], var_x, var_s2)
normal_inverse_gamma_rvs(state, x, s2 [, mu, lmbda, a, b])
normal_inverse_gamma_rvs_array(state, x, s2 [, mu, lmbda, a, b])
```

Defaults are `mu=0`, `lmbda=1`, `a=1`, and `b=1`; `lmbda`, `a`, and `b` must
be positive. The inverse-gamma component `s2` must be positive. Its mean exists
for `a > 1`, and its variance exists for `a > 2`; the marginal variance of `x`
exists for `a > 1`.

### Matrix normal

```text
matrix_normal_pdf(x [, mean, rowcov, colcov])
matrix_normal_logpdf(x [, mean, rowcov, colcov])
matrix_normal_entropy([rowcov, colcov])
matrix_normal_rvs(state, sample [, mean, rowcov, colcov, status])
matrix_normal_rvs_array(state, samples [, mean, rowcov, colcov, status])
```

Absent covariance factors are identities and an absent mean is zero. `rowcov`
and `colcov` must be SPD; only their lower triangles are authoritative for the
factorization, as in SciPy's dense matrix-distribution path. Array samples are
indexed by the third dimension.

### Wishart and inverse-Wishart

```text
wishart_pdf(x, df, scale)
wishart_logpdf(x, df, scale)
wishart_mean(df, scale)
wishart_mode(df, scale)
wishart_var(df, scale)
wishart_entropy(df, scale)
wishart_rvs(state, df, scale, sample [, status])
wishart_rvs_array(state, df, scale, samples [, status])

invwishart_pdf(x, df, scale)
invwishart_logpdf(x, df, scale)
invwishart_mean(df, scale)
invwishart_mode(df, scale)
invwishart_var(df, scale)
invwishart_entropy(df, scale)
invwishart_rvs(state, df, scale, sample [, status])
invwishart_rvs_array(state, df, scale, samples [, status])
```

For dimension `p`, both families require `df > p-1` and an SPD `scale`; only
the lower triangle of the scale matrix is used for factorization. Wishart draws
use the Bartlett decomposition. Inverse-Wishart draws currently use the exact
distributional identity obtained by inverting a Wishart draw with inverse scale;
therefore they are distribution-compatible but are not intended to reproduce
SciPy's random stream bit-for-bit. The inverse-Wishart mean requires `df > p+1`
and its elementwise variance requires `df > p+3`.

### Matrix t

```text
matrix_t_pdf(x [, mean, row_spread, col_spread, df])
matrix_t_logpdf(x [, mean, row_spread, col_spread, df])
matrix_t_rvs(state, sample [, mean, row_spread, col_spread, df, status])
matrix_t_rvs_array(state, samples [, mean, row_spread, col_spread, df, status])
```

The SciPy 1.17 matrix-t parameterization uses positive `df`, a mean matrix, and
SPD row/column spread matrices; defaults are zero mean, identity spreads, and
`df=1`. Density evaluation uses the two-sided determinant formula. Random draws
use the equivalent inverse-Wishart/matrix-normal mixture representation.

### Uniform direction and von Mises-Fisher

```text
uniform_direction_rvs(state, sample [, status])
uniform_direction_rvs_array(state, samples [, status])

vonmises_fisher_pdf(x, mu, kappa)
vonmises_fisher_logpdf(x, mu, kappa)
vonmises_fisher_entropy(mu, kappa)
vonmises_fisher_fit(data, mu, kappa, status)
vonmises_fisher_rvs(state, mu, kappa, sample [, status])
vonmises_fisher_rvs_array(state, mu, kappa, samples [, status])
```

`uniform_direction` infers the sphere dimension from the output vector and normalizes
an iid Gaussian vector. Von Mises-Fisher requires unit vectors of dimension at least
two and `kappa > 0`; use `uniform_direction` for the zero-concentration case. Array
samples and fit observations are stored by column. The sampler uses a circular
construction in two dimensions, a stable direct sphere method in three, and Wood's
rejection method in dimension four and above.

### Haar matrix groups and random correlations

```text
ortho_group_rvs(state, sample [, status])
ortho_group_rvs_array(state, samples [, status])
special_ortho_group_rvs(state, sample [, status])
special_ortho_group_rvs_array(state, samples [, status])
unitary_group_rvs(state, sample [, status])
unitary_group_rvs_array(state, samples [, status])
random_correlation_rvs(state, eigs, sample [, status, tol, diag_tol])
```

The group generators infer dimension from square output matrices and generate Haar
orthogonal, determinant-`+1` special-orthogonal, or complex unitary draws from
Gaussian QR constructions. `random_correlation_rvs` requires more than one
nonnegative eigenvalue whose sum equals the dimension within `tol` (default
`1e-13`); it uses the Davies-Higham orthogonal-similarity/Givens construction and
checks the final unit diagonal within `diag_tol` (default `1e-7`).

### Random contingency tables

```text
random_table_pmf(x, row, col)
random_table_logpmf(x, row, col)
random_table_mean(row, col)
random_table_rvs(state, row, col, sample [, status])
random_table_rvs_array(state, row, col, samples [, status])
```

`row` and `col` are nonnegative integer marginals with equal totals. A table with
valid dimensions but mismatched marginals has zero mass. Sampling is exact via
sequential row-conditional multivariate-hypergeometric draws. Matching SciPy 1.17,
the mean of a zero-total table is NaN. Array draws are indexed by the third
dimension.

## Gaussian kernel density estimation

```text
gaussian_kde_init(kde, dataset, status [, weights, bw_method, bw_factor])
kde%pdf(points [, status])
kde%logpdf(points [, status])
kde%scotts_factor()
kde%silverman_factor()
kde%set_bandwidth(method_or_factor, status)
kde%integrate_gaussian(mean, cov, status)
kde%integrate_box_1d(low, high, status)
kde%integrate_box(low, high [, maxpts, state, status])
kde%integrate_kde(other, status)
kde%resample(state, samples, status)
kde%marginal(dimensions, status)
kde%inv_cov(result, status)
```

`type(gaussian_kde)` stores observations by column, `dataset(d,n)`, matching
SciPy's mathematical data orientation. A rank-one initializer is also available
for univariate data. Optional weights are normalized and the effective sample size
is `1/sum(weights**2)`; weighted covariance uses the unbiased frequency/analytic-
weight normalization used by NumPy/SciPy. `bw_method` accepts `scott` or
`silverman`, while `bw_factor` supplies a positive constant covariance factor.
Python callable bandwidth rules are intentionally not part of the Fortran API.

Multidimensional box integration reuses SciFort's multivariate-normal rectangular
probability kernel. Marginal indices are one-based and unique. Resampling uses the
project's explicit `rng_state`; streams are deterministic within SciFort but are not
intended to match NumPy random streams bit-for-bit.

## Quasi-Monte Carlo

The QMC layer is exported through `scifort_stats` and `scifort_qmc`:

```text
qmc_scale(sample, lower, upper, scaled, status [, reverse])
qmc_to_integers(sample, lower, upper, values, status [, endpoint])
qmc_discrepancy(sample [, method, iterative, status])
qmc_update_discrepancy(x_new, sample, initial_disc [, status])
qmc_geometric_discrepancy(sample [, method, metric, status])
qmc_van_der_corput(n, base, sequence, status [, start_index, scramble, state])
qmc_halton_init(engine, d, status [, scramble, state])
engine%random(n, sample, status)
engine%fast_forward(n [, status])
engine%reset()
qmc_sobol_init(engine, d, status [, scramble, bits, state])
engine%random(n, sample, status)
engine%random_base2(m, sample, status)
engine%fast_forward(n [, status])
engine%reset()
qmc_latin_hypercube_init(engine, d, status [, scramble, strength])
engine%random(state, n, sample, status)
engine%reset()
qmc_multinomial_init(dist, pvals, n_trials [, state, status, scramble, bits])
dist%random(n, sample, status)
dist%reset()
qmc_multivariate_normal_init(dist, mean, status [, state, cov, cov_root, inv_transform, scramble, bits])
dist%random(n, sample, status)
dist%reset()
qmc_poisson_disk_init(engine, d, status [, radius, hypersphere, ncandidates, lower, upper])
engine%random(state, n, sample, status [, n_drawn])
engine%fill_space(state, sample, status)
engine%reset()
```

Unlike the multivariate distribution samplers, QMC samples use SciPy's row layout
`sample(n,d)`. `qmc_discrepancy` accepts `CD`, `WD`, `MD`, or `L2-star` and follows
SciPy's convention that `iterative=.true.` evaluates the existing sums with `n+1`
denominators so a candidate can be added with `qmc_update_discrepancy`.
`qmc_geometric_discrepancy` provides `mindist` and `mst` with Euclidean or
`cityblock` distance.

Unscrambled Halton uses successive prime bases and is sequence-compatible with
SciPy, including zero-based sequence indexing, reset, and fast-forward. Scrambling
uses the same finite sequence of independent base-digit permutations implied by
binary64 precision, but those permutations are drawn from SciFort's explicit RNG,
so a SciFort seed is reproducible within SciFort rather than NumPy-stream identical.
A scrambled initializer therefore requires an explicit `rng_state`.

Latin hypercubes require an explicit RNG state at draw time. Strength one places
one point in every one-dimensional stratum. Strength two requires `n=p**2` for a
prime `p` and `d<=p+1`, constructs a randomized orthogonal array over `p` symbols,
and refines every symbol block with a one-dimensional Latin hypercube.

Sobol uses the Joe-Kuo search-criterion-6 direction-number table distributed by
SciPy 1.17.0. Unscrambled sequences are numerically identical to SciPy through the
full supported dimension 21,201; `bits` may be 1 through 64 and limits the sequence
to `2**bits` distinct points. `random_base2` enforces the cumulative power-of-two
balance rule. LMS+digital-shift scrambling uses SciFort's explicit RNG state and is
therefore deterministic within SciFort but not NumPy-bitstream compatible.

`qmc_multinomial` consumes a one-dimensional Sobol stream to categorize repeated
uniform draws. `qmc_multivariate_normal` uses Sobol in `d` dimensions for inverse-
normal transformation, or the next even dimension for Box-Muller, and accepts a
symmetric PSD covariance (including singular matrices) or an explicit square root.

`qmc_poisson_disk` implements Bridson-style active-pool sampling in arbitrary finite
rectangular bounds. `hypersphere` is `volume` or `surface`; `random` may return fewer
than the requested number when the active pool is exhausted, with the actual count
reported by `n_drawn`, and `fill_space` continues until exhaustion. Neighbor checks
are exact against the accepted set rather than grid-accelerated, favoring a compact
dependency-free implementation over SciPy's spatial acceleration. The optional
SciPy `random-cd` and Lloyd post-sampling optimizers are not yet implemented; Lloyd
is deferred until SciFort has a Voronoi/Qhull-capable geometry layer.

## Special functions

```fortran
use scifort_special
```

Milestone 6 adds the following elemental scalar interfaces:

```text
erf(x), erfc(x)
erfinv(y), erfcinv(y)
gammaln(x)
betaln(a, b)
digamma(x), psi(x)
ndtr(x), log_ndtr(x), ndtri(p)
expit(x), logit(p), log_expit(x)
xlogy(x, y), xlog1py(x, y)
entr(x), rel_entr(x, y), kl_div(x, y)
exprel(x), cosm1(x)
boxcox(x, lambda), boxcox1p(x, lambda)
inv_boxcox(y, lambda), inv_boxcox1p(y, lambda)
i0e(x), i1e(x)
log_besseli_scaled(nu, x), besseli_ratio(nu, x)
log_besselk(nu, x)
log_besselk_order_derivative(nu, x), besselk_log_derivative_x(nu, x)
hurwitz_zeta(s, q), hurwitz_zeta_derivative(s, q)
```

`erf` and `erfc` wrap the Fortran intrinsic implementations. `erfinv` and
`erfcinv` are formed through the standard-normal inverse identity and therefore
inherit the accuracy and extreme-tail range of the current `ndtri`/normal-PPF
kernel. `gammaln` returns `log(abs(Gamma(x)))` on the real axis and uses the
reflection formula for negative nonintegers; nonpositive integer poles map to
`+infinity`. `digamma` uses recurrence/asymptotics for positive arguments and
the reflection formula for negative nonintegers; negative integer poles return
NaN. `psi` is an alias of `digamma`. The current public `betaln` contract is
limited to finite positive `a` and `b`, which is the parameter domain needed by
SciFort's statistical distributions.

The logistic and entropy helpers follow SciPy's real scalar conventions at
zero and at ordinary finite domain boundaries. `exprel` and `cosm1` use
cancellation-resistant formulas near zero. The Box-Cox routines use
`expm1`/`log1p`-style kernels so very small `lambda` does not collapse to a
naive subtractive formula.

`i0e(x)` and `i1e(x)` evaluate `exp(-abs(x))*I0(x)` and
`exp(-abs(x))*I1(x)`, respectively. The shared dependency-free implementation
uses convergent power series at moderate arguments and scaled asymptotics for
large magnitudes; `i0e` is even and `i1e` is odd. These kernels are currently
used by the Rice and von Mises families. `log_besselk` evaluates
`log(K_nu(x))` for real finite order and positive real argument using a scaled
integral representation with adaptive Gauss-Legendre panels and a large-`x`
asymptotic branch. `log_besselk_scaled(nu,x)` returns `log(exp(x)*K_nu(x))`
directly, avoiding cancellation and overflow in Bessel-normalized heavy-tail
families. The two derivative routines return derivatives of
`log(K_nu(x))` with respect to order and argument. `hurwitz_zeta` evaluates the real
Hurwitz zeta for finite `s > 1`, `q > 0` using Euler-Maclaurin summation;
`hurwitz_zeta_derivative` returns its derivative with respect to `s`.

Three one-dimensional pure reductions are also available:

```text
logsumexp(x)
softmax(x)
log_softmax(x)
```

They use max-shifting to avoid avoidable overflow. The current API is
unweighted and one-dimensional: SciPy's `axis`, `b`, `keepdims`, and
`return_sign` options are not yet implemented. `logsumexp` of an empty array is
`-infinity`; `softmax` and `log_softmax` return an empty allocatable array for
an empty input.

The regularized incomplete-gamma API remains:

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

Most current distribution routines are scalar elemental functions.
`poisson_binom_*` is the exception because its distribution parameter is the
rank-one probability vector `p(:)`. For example:

```fortran
real(dp) :: x(5)
real(dp) :: p(5)

p = normal_cdf(x, loc=1.0_dp, scale=2.0_dp)
```

Optional dummy arguments are scalar in an elemental invocation. Conformable
array arguments may be supplied for arguments that are present.
