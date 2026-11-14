## Complex Schur decomposition and reordering (2026-11-08)

`scifort_linalg_schur.schur(a,t,z,status,sort,sdim,selector)` accepts real
or complex matrices and returns a **complex** upper-triangular Schur form
with a unitary basis. The `sort` values are `lhp`, `rhp`, `iuc`, and `ouc`;
`sdim` is the selected invariant-subspace dimension. A user-defined
`selector(w)` can replace `sort`. This is SciPy's **complex-output** Schur
subset; SciPy's real quasi-triangular form is not yet implemented.
[Full API, status codes, validation and limitations](linalg_schur.md).

## Packed-band linear systems (`scipy.linalg`)

`scifort_linalg_banded` provides real and complex `solve_banded`,
`solveh_banded`, `cholesky_banded`, and `cho_solve_banded` for packed general
and Hermitian/symmetric positive-definite banded matrices. Both vectors and
multiple right-hand sides are supported; the Cholesky factor is reusable.
See [storage conventions, examples and limits](linalg_banded.md) and
[validation results](linalg_banded_validation_2026-11-07.md).

## Dense matrix equations (`scipy.linalg`)

Use `scifort_linalg_matrix_equations` for real or complex
`solve_sylvester`, `solve_continuous_lyapunov`, and
`solve_discrete_lyapunov`. These routines return a SciFort integer status
and are designed for small and moderate matrix equations. See
[API and limitations](linalg_matrix_equations.md). These now offer `method="auto"`, `"direct"`, or `"schur"` for
Schur-reduced solutions up to matrix order 40 (see the linked page).

## Complex sparse matrices, direct and iterative linear solvers

See [Complex CSR/CSC and Krylov solvers](sparse_complex.md) for complex128
compressed sparse storage, conversions, conjugate transpose, and
`cg_complex`, `bicgstab_complex`, and `gmres_complex`. See also the
[complex sparse direct LU API](sparse_complex_direct.md).

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
observations. Separate normal, exponential, Weibull, lognormal, gamma, beta, logistic, and both Gumbel censored likelihoods
and censored-data fitting overloads are described in
[censored likelihoods and fitting](censored_fitting.md). Observation weights
and censored fitting for other distribution families are not yet implemented.

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
gmean(x [, weights, nan_policy])
hmean(x [, weights, nan_policy])
pmean(x, p [, weights, nan_policy])
sem(x [, ddof, nan_policy])
variation(x [, ddof, nan_policy])
iqr(x [, rng, scale, nan_policy])
median_abs_deviation(x [, scale, nan_policy])
moment(x [, order, center, nan_policy])
skew(x [, bias, nan_policy])
kurtosis(x [, fisher, bias, nan_policy])
quantile(x, q)
median(x)
covariance(x, y [, ddof])
pearson_correlation(x, y)
rankdata(x [, method, axis, nan_policy])
tiecorrect(rankvals)
call mode(a, mode_value, count_value [, axis, nan_policy])
call mode_all(a, mode_value, count_value [, nan_policy])
call mode_keepdims(a, mode_value, count_value [, axis, nan_policy])
call describe(a, nobs, min, max, mean, variance, skewness, kurtosis
              [, axis, ddof, bias, nan_policy])
call describe_all(a, nobs, min, max, mean, variance, skewness, kurtosis
                  [, ddof, bias, nan_policy])
percentileofscore(a, score [, kind, nan_policy])
scoreatpercentile(a, per [, limit, interpolation_method, axis])
scoreatpercentile_all(a, per [, limit, interpolation_method])
cumfreq(a [, numbins, defaultreallimits, weights])
relfreq(a [, numbins, defaultreallimits, weights])
zscore(x [, axis, ddof, nan_policy])
gzscore(x [, axis, ddof, nan_policy])
gstd(x [, axis, ddof, nan_policy])
gstd_keepdims(x [, axis, ddof, nan_policy])
```

`variance`, `standard_deviation`, and `covariance` default to `ddof = 1`;
`ddof = 0` gives the population divisor. Mean and second-order reductions use
scaled/compensated calculations so large finite observations are not summed or
squared naively. Variance/covariance-style reductions require finite data and
return NaN when their divisor is nonpositive. Empty reductions return NaN.

`gmean`, `hmean`, and `pmean` follow SciPy's scalar/vector definitions. Optional
weights must have the same length as the data, be finite, and have nonzero total
weight; an invalid weight vector returns NaN. `gmean` and `hmean` require
nonnegative observations. `pmean` likewise requires nonnegative observations,
with `p=0` reducing to `gmean`. The power-mean implementation detects power
overflow/underflow before exponentiation so floating-point trapping builds return
the corresponding infinity/zero rather than raising an exception.

`sem` defaults to `ddof=1`; `variation` defaults to `ddof=0` and preserves the
sign of the ordinary arithmetic mean. `iqr` accepts percentile endpoints in
`rng` (default `[25,75]`), sorts reversed endpoints, and divides the percentile
range by `scale`. `median_abs_deviation` uses the sample median as its center and
also divides by a numeric `scale`.

These ten SciPy-style helpers accept vectors and native Fortran arrays through
rank 15. Array calls use a one-based Fortran `axis`; axis 1 is the analogue of
SciPy/NumPy axis 0. Separate `*_keepdims` names retain the reduced dimension
because a Fortran function result cannot change rank at runtime. Ranks 5 through
15 use shared assumed-rank slice engines and thin fixed-result-rank wrappers,
while the statistical formulas remain in the vector kernels. Rank-2 `gmean`,
`hmean`, and `pmean` accept rank-1/rank-2 trailing-broadcast weights. Rank-3
mean-family wrappers accept rank-1/rank-2/rank-3 trailing broadcast weights.
Rank-4 through rank-15 mean-family wrappers additionally accept scalar or
trailing-rank weights up to the input rank. These wrappers can expand singleton
data dimensions when the joint NumPy broadcast shape is larger than `x`. `iqr`
without an axis reduces all elements for rank-2 through rank-15 input, matching
SciPy's `axis=None` default.

For these helpers, `nan_policy='propagate'` preserves the existing NaN behavior,
`nan_policy='omit'` removes NaN observations before each vector/slice reduction
through rank 15, and `nan_policy='raise'` performs Fortran 2018 `error stop` when
a NaN is encountered. Weighted means apply the same policy when either an
observation or its corresponding broadcast weight is NaN. Infinities are not
omitted. If omission leaves too few observations for the requested
statistic/divisor, the result is NaN. `iqr` and `median_abs_deviation` accept a numeric scale or
the case-insensitive character scale `scale='normal'` through rank 15, including
axis, keepdims, and NaN-policy forms. The normal conversion divides IQR by
`2 * Phi^-1(0.75)` and MAD by `Phi^-1(0.75)`, matching SciPy. Unsupported named
scales return NaN in the native pure-function API. Unsupported `nan_policy`
names also return NaN. Dtype selection, custom MAD center callbacks,
and alternate percentile interpolation modes are not yet exposed.

`moment` defaults to the first moment about the sample mean and accepts an
optional finite scalar `center`. `skew` returns the Fisher-Pearson coefficient;
`bias=.false.` applies SciPy's finite-sample correction when `n>2`. `kurtosis`
uses Fisher's definition by default, `fisher=.false.` returns Pearson kurtosis,
and `bias=.false.` applies the unbiased correction when `n>3`. Constant samples
return NaN for skewness and kurtosis, matching SciPy.

`quantile` implements linear interpolation at position `(n - 1) q`, equivalent
to the common Hyndman-Fan type 7 / NumPy `method="linear"` convention. `q`
may be a scalar or an array. NaN observations or probabilities outside
`[0, 1]` return NaN. `median(x)` is `quantile(x, 0.5)`.

`rankdata` returns `real(dp)` ranks beginning at one for inputs through rank 15. The
default method is `average`; `min`, `max`, `dense`, and stable `ordinal` are also
supported. Without `axis`, multidimensional input is flattened in Fortran storage order
and a rank-1 result is returned, the natural Fortran analogue of SciPy `axis=None`. An
explicit one-based `axis` ranks each slice independently and preserves the input shape.
`nan_policy='propagate'` returns all NaNs for an affected slice, `'omit'` ignores NaNs
while retaining NaNs at their original positions, and `'raise'` terminates with
`error stop`. The historical invalid-method behavior is retained: an unrecognized tie
method returns NaNs rather than changing existing callers to an exception. Infinities
are ordered normally. `tiecorrect(rankvals)` returns the standard correction factor used
by Mann-Whitney U and Kruskal-Wallis H calculations; NaNs count toward the sample size
but are not treated as tied with one another, matching SciPy behavior.

`mode` and `describe` use subroutines because their SciPy result objects contain multiple
fields whose result rank depends on the selected axis. The output arguments are allocatable
and must be declared with rank one less than the input for `mode`/`describe`;
`mode_keepdims` outputs have the same rank as the input. `mode_all` and `describe_all` use
scalar outputs. `mode` chooses the most frequent value and the smallest value on a count
tie; with `nan_policy='propagate'`, all NaNs in a slice form one countable category, while
`'omit'` removes them and `'raise'` terminates. Empty mode samples return NaN with count zero.
`describe` returns nobs, min/max, arithmetic mean, variance, skewness, and Fisher kurtosis.
Its variance follows SciPy's unrestricted integer `ddof` convention, so `ddof == n` can
produce positive infinity and `ddof > n` can produce a negative value. The `bias` argument
controls finite-sample corrections for skewness and kurtosis. For SciPy 1.17 compatibility,
when `nan_policy='omit'` is activated by a NaN, the masked-statistics path excludes all
non-finite values in the input before describing each slice.

`percentileofscore(a, score)` compares scalar or rank-1 through rank-15 scores
against a one-dimensional reference sample. The default `kind='rank'` averages the
percentage rankings of matching observations; `weak`, `strict`, and `mean` implement
the other SciPy conventions. `nan_policy='propagate'`, `'omit'`, and `'raise'` are
supported. The result has the same shape as an array-valued `score`.

`scoreatpercentile(a, per)` reduces rank-1 through rank-15 data over the whole array
by default. An explicit one-based `axis` performs slice-wise reduction. Scalar
percentiles return the reduced shape, while a percentile vector becomes the leading
result dimension, matching NumPy/SciPy layout. Optional `limit=[lo,hi]` includes
only observations in the closed interval, and `interpolation_method` accepts
`fraction`, `lower`, or `higher`. `scoreatpercentile_all` is an explicit whole-array
alias.

`cumfreq` and `relfreq` flatten rank-1 through rank-15 input into an equal-width
histogram. They return `cumfreq_result` and `relfreq_result`, respectively; each
record contains an allocatable `cumcount` or `frequency` vector plus `lowerlimit`,
`binsize`, and `extrapoints`. Optional `numbins`, `defaultreallimits=[lo,hi]`, and
shape-matched weights are supported. `relfreq` intentionally follows SciPy's current
normalization convention for multidimensional input rather than renormalizing the
flattened histogram to one.

`zscore` and `gzscore` return arrays with the same rank and shape as the input.
They standardize along one one-based Fortran `axis` (default 1) with `ddof=0` by
default. `gzscore` standardizes natural logarithms, matching SciPy's geometric
z-score definition. `zmap(scores,compare,axis,...)` uses the comparison array for
its location and scale; its reduced-axis extent may differ from `scores`, and
singleton non-reduced dimensions broadcast jointly. `gstd` returns
`exp(std(log(x)))`, defaults to `ddof=1`, and uses `gstd_keepdims` when the
reduced dimension must be retained. All four support ranks 1 through 15 and
`nan_policy='propagate'` / `'omit'` / `'raise'`.

The explicit `zscore_all`, `gzscore_all`, `zmap_all`, and `gstd_all` companions
are the fixed-signature Fortran analogue of SciPy `axis=None`: statistics are
computed over every element, while the score functions retain the score-array
shape and `gstd_all` returns a scalar. `zmap_all` currently requires `scores` and
`compare` to have the same rank, although their extents may differ because the
comparison statistics are scalar. `gstd_keepdims(x,dims)` accepts a nonempty
list of unique one-based dimensions and reduces them jointly while retaining
each selected dimension with extent one. This provides tuple-axis `gstd`
without requiring a Fortran function result whose rank depends on the runtime
length of `dims`.

Under `omit`, a NaN input remains NaN in the `zscore`/`gzscore` output while the
other values in that slice are standardized from the non-NaN observations.
`gstd` follows SciPy's documented finite, strictly positive input domain; zero
or infinite observations return NaN. Negative observations become omitted
logarithmic NaNs when `nan_policy='omit'`; with `nan_policy='raise'`, they
terminate because the logarithmic transform produces NaN. For `zmap`, the
policy applies to the comparison data, matching SciPy: a NaN in `scores`
propagates only to the corresponding output element.

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
bws_test(x, y, state [, alternative, n_resamples])
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

`bws_test` implements the two-sample Baumgartner-Weiss-Schindler rank test and
returns a `permutation_test_result`. It supports all three alternatives.
`n_resamples` defaults to 9999; when the number of distinct independent-sample
splits does not exceed that limit, the null distribution is enumerated exactly.
Otherwise SciFort performs randomized permutations with the caller-supplied
`rng_state`, avoiding hidden mutable random state.

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
qmc_random_cd(sample, state, status [, n_iters, n_nochange, n_accepted, n_performed])
qmc_lloyd(sample, status [, tol, maxiter, n_iterations])
qmc_van_der_corput(n, base, sequence, status [, start_index, scramble, state])
qmc_halton_init(engine, d, status [, scramble, state, optimization])
engine%random(n, sample, status [, state])
engine%fast_forward(n [, status])
engine%reset()
qmc_sobol_init(engine, d, status [, scramble, bits, state, optimization])
engine%random(n, sample, status [, state])
engine%random_base2(m, sample, status [, state])
engine%fast_forward(n [, status])
engine%reset()
qmc_latin_hypercube_init(engine, d, status [, scramble, strength, optimization])
engine%random(state, n, sample, status)
engine%reset()
qmc_multinomial_init(dist, pvals, n_trials [, state, status, scramble, bits])
dist%random(n, sample, status)
dist%reset()
qmc_multivariate_normal_init(dist, mean, status [, state, cov, cov_root, inv_transform, scramble, bits])
dist%random(n, sample, status)
dist%reset()
qmc_poisson_disk_init(engine, d, status [, radius, hypersphere, ncandidates, lower, upper, optimization])
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

`qmc_random_cd` performs random coordinate descent for centered discrepancy by
proposing exchanges of two rows within one coordinate column. Accepted exchanges
are verified against the full centered discrepancy, so the returned sample never
increases CD apart from roundoff. The defaults are at most 10000 proposals and a
stop after 100 consecutive non-improving proposals. `n_accepted` and
`n_performed` are optional SciFort diagnostics. Pairwise coordinate exchange
preserves every column value exactly.

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

Halton, Sobol, LatinHypercube, and PoissonDisk accept
`optimization='random-cd'`, `'lloyd'`, or `'none'`. Because SciFort does not hide
mutable RNG state inside an engine, random-CD Halton/Sobol draws with more than one
point and dimension require the optional draw-time `state`; Lloyd itself is
deterministic and does not add an RNG requirement. Coordinate-exchange
post-processing preserves ordinary Latin strata, while either optimizer can
invalidate stronger construction properties such as strength-two OA balance or
Poisson-disk minimum separation. PoissonDisk optimization is restricted to the
unit cube; its internal pool remains based on the unoptimized construction.

`qmc_lloyd` follows SciPy's centroidal-Voronoi relaxation convention: finite
Voronoi vertices are averaged for each site, the displacement is multiplied by an
exponentially decaying factor, proposals outside the unit cube are discarded, and
iteration stops when the change in minimum city-block spacing is below `tol`
(default `1e-5`) or after `maxiter` iterations (default 10). SciFort reconstructs
finite Voronoi vertices with its own dependency-free N-dimensional Bowyer-Watson
Delaunay triangulation rather than Qhull. For generic nondegenerate inputs this is
validated against unjoggled Qhull geometry. SciPy requests Qhull's `QJ` joggle;
SciFort deliberately does not, so collinear, cospherical, or otherwise degenerate
designs can follow a different valid triangulation or return an invalid-parameter
status. As with SciPy's finite-vertex algorithm, Voronoi cells are not clipped to
the unit-cube boundary before centroids are formed.

`qmc_multinomial` consumes a one-dimensional Sobol stream to categorize repeated
uniform draws. `qmc_multivariate_normal` uses Sobol in `d` dimensions for inverse-
normal transformation, or the next even dimension for Box-Muller, and accepts a
symmetric PSD covariance (including singular matrices) or an explicit square root.

`qmc_poisson_disk` implements Bridson-style active-pool sampling in arbitrary finite
rectangular bounds. `hypersphere` is `volume` or `surface`; `random` may return fewer
than the requested number when the active pool is exhausted, with the actual count
reported by `n_drawn`, and `fill_space` continues until exhaustion. Neighbor checks
are exact against the accepted set rather than grid-accelerated, favoring a compact
dependency-free implementation over SciPy's spatial acceleration. SciPy-style
`random-cd` and `lloyd` post-sampling optimization are implemented. Lloyd uses
SciFort's native dependency-free Delaunay/Voronoi reconstruction rather than a
Qhull dependency.

## Numerical differentiation

The `scipy.differentiate` compatibility layer is exposed by:

```fortran
use scifort_differentiate
```

The adaptive real finite-difference API is:

```text
derivative(f, x [, atol, rtol, maxiter, order, initial_step, step_factor, step_direction])
jacobian(f, x [, atol, rtol, maxiter, order, initial_step, step_factor, step_direction])
hessian(f, x [, atol, rtol, maxiter, order, initial_step, step_factor])
```

`derivative` accepts scalar or rank-1 arrays of abscissae. Scalar results contain
`df`, `error`, `success`, `status`, `nit`, `nfev`, and `x`; array results expose
the same quantities elementwise. The default controller follows SciPy 1.17: order
8, initial step 0.5, step factor 2, at most 10 iterations, central differences by
default, and one-sided differences for negative/positive `step_direction`. Odd
orders are promoted to the next even stencil order.

`jacobian` accepts either a scalar-valued or allocatable-vector-valued callback of
a real vector argument and returns output-by-input matrices of estimates and
convergence metadata. `hessian` accepts a scalar-valued vector callback and uses
nested adaptive differentiation with a tighter inner relative tolerance, following
SciPy's construction. Because Fortran callbacks are scalar calls rather than
NumPy-vectorized calls, exact `nfev` values can differ in difficult entries even
when estimates/statuses agree; the reported count is the number of original scalar
function calls made by SciFort. Arbitrary Python `args`, `preserve_shape`, and
iteration callbacks are represented naturally by Fortran closures or are not
currently exposed as separate keyword features.

SciFort also provides analytic-function extensions:

```text
complex_step_derivative
complex_step_gradient
complex_step_jacobian
```

These use an imaginary perturbation (default `1e-30`) and are not part of SciPy's
public `differentiate` method selection. The scalar formula is adapted from
ForDiff 0.1.1; the vector forms are SciFort implementations.

## Optimization

The initial `scipy.optimize` compatibility layer is exposed by:

```fortran
use scifort_optimize
```

### Quadratic assignment (FAQ, 2-opt)

`quadratic_assignment(a, b, method='faq', maximize=.false., ...)`
approximately solves the permutation objective
`sum(a(i,j) * b(p(i),p(j)))`. Two independent methods are implemented:
FAQ (Frank-Wolfe relaxation and Hungarian rounding) and best-improvement
2-opt. Fixed seeded matches are accepted via `partial_match(:,2)`, with
one-based indices; optional `partial_guess(:,2)` initializes 2-opt and
optional doubly stochastic `p0` initializes FAQ. Both methods return
`qap_result` containing `col_ind`, `fun`, `nit`, `status`, and `success`.
Use `quadratic_assignment_score(a,b,p)` for a separate objective evaluator.
The algorithms are approximate, and their deterministic initializations
and iteration counters differ from SciPy. See
[the detailed QAP documentation](optimize_quadratic_assignment.md).

### Linear sum assignment and isotonic regression

Two additional SciPy-style optimization algorithms are available without
external dependencies:

```fortran
use scifort_optimize, only : dp, linear_sum_assignment, isotonic_regression, &
    isotonic_result, assignment_status_success
integer, allocatable :: rows(:), cols(:)
integer :: status
type(isotonic_result) :: fit
real(dp) :: costs(2, 3)

costs(1, :) = [4.0_dp, 1.0_dp, 6.0_dp]
costs(2, :) = [2.0_dp, 0.0_dp, 5.0_dp]
call linear_sum_assignment(costs, rows, cols, status)
! rows = [1, 2], cols = [2, 1], objective = 3.
fit = isotonic_regression([3.0_dp, 1.0_dp, 2.0_dp, 5.0_dp, 4.0_dp])
! fit%x = [2, 2, 2, 4.5, 4.5]
! fit%weights = [3, 2], fit%blocks = [1, 4, 6].
```

`linear_sum_assignment(cost_matrix, row_ind, col_ind, status, maximize)`
accepts a rank-two `real(dp)` cost matrix. When `maximize` is omitted it finds
minimum-cost pairs; use `maximize=.true.` for maximum weight matching. The
number of pairs is `min(size(cost_matrix,1),size(cost_matrix,2))`; row indices
are sorted and unique, and column indices are unique. A `+infinity` cost
forbids an edge in minimization; `-infinity` forbids an edge in maximization.
Other nonfinite costs are rejected. An empty dimension produces empty index
arrays and a successful status. Status constants are
`assignment_status_success=0`, `assignment_status_invalid_input=1`, and
`assignment_status_infeasible=2`. Unlike SciPy's Jonker-Volgenant kernel,
SciFort uses a rectangular primal-dual Hungarian shortest-augmenting-path
implementation with O(min(n,m)^2 max(n,m)) time, O(n+m) auxiliary memory,
and scaled costs to avoid overflow for very large finite entries. With tied
optima, the pairing need not match SciPy's choice.

`fit = isotonic_regression(y [, weights, increasing])` applies the weighted
pool-adjacent-violators algorithm in O(n) time and auxiliary memory.
Weights must be strictly positive and finite and default to one.
`increasing=.false.` fits a nonincreasing sequence. The returned
`isotonic_result` contains `x` (length n), `weights` (pooled sums), and
`blocks` (the 1-based start of each pool followed by n+1). For example
`fit%blocks=[1,4,6]` denotes pools 1:3 and 4:5. A length-zero input
succeeds with `fit%blocks=[1]`. `fit%status` is zero on success, 1 on invalid
input (including nonfinite observations or weights), or 2 on numerical
failure (such as pooled weight overflow). Invalid inputs return empty
allocated arrays rather than raising Python exceptions.

Both algorithms have regression cases checked against SciPy 1.17.0 and
independent randomized comparisons. No external source code was imported;
see `CODE_PROVENANCE.md`.

Multivariate `root` also supports SciPy-style `tol_norm=` for
`method='broyden1'`, `'broyden2'`, `'anderson'`, `'linearmixing'`,
`'diagbroyden'`, `'excitingmixing'`, and `'krylov'`. The same keyword is
supported by the direct functions `broyden1`, `broyden2`, `anderson`,
`linearmixing`, `diagbroyden`, `excitingmixing`, and `newton_krylov`.
A norm callback accepts a one-dimensional `real(dp)` vector and returns a
nonnegative finite `real(dp)` scalar. It controls absolute and relative
residual, iterate, and accepted-step convergence tests; it does not change
the Jacobian update, line search, or Krylov inner linear norm. Without a
callback, the maximum absolute element norm is used. A nonfinite/negative
norm results in unsuccessful status 0, not an exception. The initial point
can only satisfy convergence if its residual is exactly zero.

Example:

```fortran
function weighted_norm(v) result(value)
    real(dp), intent(in) :: v(:) !! vector to measure
    real(dp) :: value !! infinity norm after component scaling
    value = maxval(abs(v) * [0.1_dp, 1.0_dp])
end function weighted_norm
```

```fortran
r = root(f, [0.0_dp, 0.0_dp], method='broyden1', &
    fatol=0.08_dp, tol_norm=weighted_norm)
```

`tol_norm=` is rejected for methods such as `hybr`, `lm`, and `df-sane`.


### Matrix-free Krylov preconditioning and custom inner solvers

`root(method='krylov', ...)` and `newton_krylov(...)` accept three optional
Fortran procedure callbacks modeled on SciPy 1.17's `inner_M` and callable
inner methods:

- `inner_m(vec, product, ok)` applies an approximate inverse Jacobian to
  a vector. Its `real(dp)` input/output are rank-one arrays and its logical
  output indicates whether the application succeeded. SciFort uses **right
  preconditioning** for its matrix-free GMRES and LGMRES Arnoldi steps;
  LGMRES augmentation vectors, already in solution coordinates, are not
  preconditioned a second time. Without `inner_m=`, the identity is used.
- `inner_m_update(x, fvec)` updates user-held preconditioner state after **each
  accepted nonlinear step**, with the accepted iterate and its residual.
  The user initializes the preconditioner for the initial iterate; it is not
  updated before the first linear solve. This callback does not perform a
  residual evaluation and does not increment `nfev`.
- `inner_solver(matvec, preconditioner, rhs, rtol, maxiter, solution, ok)`
  replaces the built-in inner solver. The first two arguments are procedure
  callbacks with the same `(vec, product, ok)` interface as `inner_m`.
  `matvec` applies the finite-difference Jacobian at the current nonlinear
  iterate and **counts each residual evaluation in `nfev`**. `preconditioner`
  applies `inner_m` or the identity if omitted. The requested `rtol` is the
  nonlinear forcing-term linear tolerance; `maxiter` is `inner_maxiter`.
  The solver returns an approximate solution to `J * solution = rhs` and
  sets `ok` to `.false.` on failure.

The public abstract callback interfaces are `vector_root_linear_operator`,
`vector_root_preconditioner_update` and `vector_root_inner_solver`. Failed
operators or inner solvers, and nonfinite outputs, yield an unsuccessful root
result rather than an uninitialized numerical iterate. These keywords are
rejected on non-Krylov `root` methods.

```fortran
program example_krylov_preconditioner
    use scifort_kinds, only : dp
    use scifort_optimize, only : root, vector_root_result
    implicit none

    type(vector_root_result) :: r

    r = root(residual, [0.0_dp, 0.0_dp], method='krylov', &
        inner_method='gmres', inner_maxiter=1, inner_m=inverse_diagonal)
    if (.not. r%success) error stop 'Krylov solve failed'
    print *, r%x

contains

    subroutine residual(x, f)
        real(dp), intent(in) :: x(:) !! unknown vector
        real(dp), intent(out) :: f(:) !! residual

        f = [2.0_dp, 5.0_dp] * x - [4.0_dp, 10.0_dp]
    end subroutine residual

    subroutine inverse_diagonal(vec, product, ok)
        real(dp), intent(in) :: vec(:) !! input direction
        real(dp), intent(out) :: product(:) !! preconditioned direction
        logical, intent(out) :: ok !! true when application succeeds

        product = vec / [2.0_dp, 5.0_dp]
        ok = .true.
    end subroutine inverse_diagonal
end program example_krylov_preconditioner
```

Six built-in `inner_method=` strings are supported: `lgmres` (default),
`gmres`, `bicgstab`, `cgs`, `minres`, and `tfqmr`. The latter four are
implemented by `scifort_optimize_krylov_iterative`, using matrix-free
finite-difference products (all counted in `nfev`) and the optional `inner_m`
action. For BiCGSTAB and CGS, `inner_m` acts on each search vector before the
Jacobian product, corresponding to right preconditioning. TFQMR retains
SciPy's `M*J` preconditioned recurrence by default, with an opt-in `J*M`
variant described below.
MINRES requires a symmetric Jacobian and a positive-definite symmetric
preconditioner; an indefinite preconditioner causes an unsuccessful result.
These properties are the caller's responsibility and are not checked
exhaustively at runtime. `outer_k=` applies only to LGMRES.

TFQMR uses an internal residual *bound*, not a freshly evaluated original
linear residual for its inner stopping criterion. As with SciPy 1.17.0,
general approximate left preconditioners can make this estimate optimistic
about the unpreconditioned linear residual; a nonlinear caller should check
`res%fun`/`res%success` or use BiCGSTAB/GMRES if tight linear residual control
is important. The outer nonlinear convergence test still uses the actual
residual vector evaluated at its accepted iterate.

**SciFort extension -- verified TFQMR:** pass `inner_verify_residual=.true.`
with built-in `inner_method='tfqmr'` to `root(method='krylov', ...)` or direct
`newton_krylov(...)`. The corresponding low-level switch is
`krylov_iterative_solve(..., verify_residual=.true.)`. In this mode the solver
computes the original, *unpreconditioned* `norm2(rhs - J*step)` by making an
additional matrix-free operator call when TFQMR's estimator first indicates
convergence. If the actual residual exceeds the requested inner tolerance, the
recurrence continues. On exhaustion of `inner_maxiter`, it also checks the
original residual before reporting success; a failed check or failed operator
call reports an unsuccessful inner solve instead of returning an unchecked
finite partial iterate. Verification is enabled only for the built-in TFQMR
method; requesting it with another method, a custom `inner_solver=`, or a
non-Krylov `root` method is an invalid input. Additional Jacobian-vector
calls are included in `nfev`, so work counts intentionally differ from SciPy
when this option is enabled.

```fortran
res = root(fun, x0, method='krylov', inner_method='tfqmr', &
    inner_verify_residual=.true., inner_maxiter=80)
```

The verification option is an **accuracy guard**, not by itself a different
preconditioner algorithm. With the default TFQMR recurrence, some approximate
preconditioners still fail to attain the requested original residual regardless
of the iteration limit. The default `inner_verify_residual=.false.` preserves
SciPy 1.17 work counts and acceptance of finite incomplete iterates.

**SciFort extension -- consistent right-preconditioned TFQMR:** the optional
`inner_tfqmr_right_preconditioned=.true.` selects the composite matrix-free
operator `J*M` in the TFQMR recurrence, where `M` is the inverse-Jacobian
action passed as `inner_m=`. Since the solution update applies `M` to each
Krylov correction, this means the solver consistently solves `(J*M)*y=b` and
returns `step=M*y`. In contrast, SciPy 1.17.0's default TFQMR uses the
`M*J` recurrence while still applying `M` to the iterate updates; for a
nonidentity approximate preconditioner, that can yield a small estimated
residual yet a large original residual. SciFort **retains that SciPy behavior
by default**, so existing work-count regressions and users are unaffected.

```fortran
res = root(fun, x0, method='krylov', inner_method='tfqmr', &
    inner_m=my_inverse_jacobian, inner_maxiter=80, &
    inner_tfqmr_right_preconditioned=.true., &
    inner_verify_residual=.true.)
```

`inner_tfqmr_right_preconditioned` can also be passed to `newton_krylov`, or as
`right_preconditioned=.true.` to the low-level
`krylov_iterative_solve(...)`. It is valid only for built-in TFQMR; passing a
true value with a different method, a custom `inner_solver=`, or a non-Krylov
`root` method returns an unsuccessful result. Passing `.false.` has no effect.
It is valid without `inner_m=` (identity preconditioning). Use together with
`inner_verify_residual=.true.` to certify the original linear residual; the
verification requires additional Jacobian-vector products and may reject an
iterate if finite-difference Jacobian products lose accuracy near convergence.
Increasing `rdiff=` moderately can help with cancellation in those products.

**SciFort extension -- centered matrix-free Jacobian products:**
`root(method='krylov', ...)` and `newton_krylov(...)` accept
`jacobian_diff_scheme='3-point'`. The default `'2-point'` keeps SciPy 1.17.0's
forward directional difference, its automatic `rdiff=sqrt(epsilon)` scale, and
its residual-evaluation counts. The opt-in three-point scheme approximates
`J*v` by `(f(x+h*v) - f(x-h*v))/(2*h)`, evaluates the function twice per
nonzero product, and records both calls in `nfev`. For an omitted `rdiff`, it
uses `epsilon(1.0_dp)**(1.0_dp/3.0_dp)` rather than `sqrt(epsilon)` to balance
centered-difference truncation and floating-point cancellation. The step also
retains the existing SciPy-style scaling by the maximum norms of `x` and
`f(x)`, with normalization by the Euclidean norm of `v`. An explicitly supplied
positive, finite `rdiff=` replaces the automatic scale in either scheme.
The same stencil is used by all six built-in Krylov methods, including
GMRES/LGMRES extra residual checks, and by the matrix-free action passed to a
custom `inner_solver=`. It can be combined with preconditioning and verified
TFQMR:

```fortran
res = root(fun, x0, method='krylov', inner_method='tfqmr', &
    inner_m=approximate_inverse, inner_tfqmr_right_preconditioned=.true., &
    inner_verify_residual=.true., jacobian_diff_scheme='3-point')
```

The three-point mode rejects unrepresentable perturbations rather than
returning a spurious zero Jacobian-vector product. It also rejects nonfinite
samples; unlike an analytic Jacobian, the centered stencil requires *both*
perturbed points to lie inside the function's domain. Since the stencil costs
two residual evaluations, enabling it intentionally changes work counts from
SciPy's default. Specifying the scheme for a non-Krylov `root` method or
specifying anything other than `'2-point'`/`'3-point'` is invalid.
A centered stencil improves smooth derivative accuracy but does not guarantee
convergence near machine precision for ill-conditioned nonlinear systems.

**SciFort extension -- Jacobian step scaling:**
`jacobian_step_mode='global'` (default), `'adaptive'`, and `'componentwise'`
are accepted by `root(method='krylov')` and `newton_krylov`. The default
preserves the SciPy 1.17.0 matrix-free product and evaluation accounting.
`'componentwise'` forms the directional derivative as the sum of the
independently estimated Jacobian columns weighted by nonzero direction
components, using step `h_j = max(rdiff,8*epsilon)*max(1,abs(x_j))`.
It avoids the numerical problems caused by scaling every direction by the
largest, unrelated component of `x`. `'adaptive'` uses the original single
directional probe for well-scaled variables, and uses the coordinatewise
procedure when the active coordinate scales span more than five decimal
orders relative to the largest `x` component, when a proposed perturbation
rounds to the original coordinate, or when its relative magnitude falls below
`1e-3*rdiff`. This decision is made for each Jacobian-vector product.

Coordinatewise mode requires `nnz(v)` extra function calls per product for
`'2-point'`, or `2*nnz(v)` for `'3-point'`; a zero direction uses none.
Choosing it for dense high-dimensional problems may be substantially slower.
It is not an adaptive error estimator or a substitute for a well-scaled
nonlinear problem, and need not improve all ill-conditioned problems. A
positive finite `rdiff` is required if supplied. Extra calls are included in
`nfev`, and the mode is applied uniformly to all six built-in methods and
custom `inner_solver` callbacks. Modes are case-insensitive. Unsupported
mode names, or use with non-Krylov `root` methods, fail before evaluation.

```fortran
res = newton_krylov(fun, x0, jacobian_step_mode='adaptive', &
    jacobian_diff_scheme='3-point')
```

**SciFort extension -- analytic matrix-free Jacobian-vector products:**
`root(method='krylov', jacobian_product=...)` and
`newton_krylov(..., jacobian_product=...)` accept a direct action `J(x)*v`.
This optional Fortran-specific callback avoids numerical differentiation when
an exact product is available, including cases where the residual contains
large offsets or different scales that make subtraction inaccurate. Its
explicit interface is exported as `vector_root_jacobian_product` by
`scifort_optimize_root` (and through `scifort_optimize`). The corresponding
lower-level interface is `nonlin_jacobian_product` in
`scifort_optimize_nonlin`:

```fortran
subroutine my_jacobian_product(x, vec, product, ok)
    real(dp), intent(in) :: x(:) !! current nonlinear iterate
    real(dp), intent(in) :: vec(:) !! multiplication direction
    real(dp), intent(out) :: product(:) !! Jacobian at x times vec
    logical, intent(out) :: ok !! callback succeeded

    product = 2.0_dp * x * vec ! Example: F(x) = x*x - target*target
    ok = .true.
end subroutine my_jacobian_product
```

The callback must return the product for the **current** `x`, not necessarily
the initial point. It is called by the built-in LGMRES, GMRES, BiCGSTAB, CGS,
MINRES, and TFQMR inner methods, including GMRES's final residual check and
TFQMR's optional original-residual verification; custom `inner_solver=` hooks
also receive it as their matrix-free operator. A zero direction produces a
zero product without calling the callback. A failed callback (`ok=.false.`)
or nonfinite product terminates the inner solve without accepting its step.
The callback does not evaluate the residual `fun`, so its calls do **not**
increment `nfev`. Count analytic products in the callback if needed.

`jacobian_product=` cannot be combined with explicit finite-difference
options `rdiff=`, `jacobian_diff_scheme=`, or `jacobian_step_mode=`. It is valid
only for the Krylov `root` method, not the MINPACK `hybr`/`lm` Jacobian `jac=`
or Broyden-like methods. It is an opt-in **SciFort extension**, not a SciPy
1.17.0 argument. The existing finite-difference defaults and reference
function-evaluation counts remain unchanged.
Test coverage: `test/test_optimize_root_krylov_analytic_product.f90`.

**SciFort extension -- verify the final Krylov inner residual:** pass
`inner_check_residual=.true.` to either `root(method='krylov', ...)` or
`newton_krylov(...)`. It evaluates the final *unpreconditioned* linear
residual using the same matrix-free `J(x)*p` operator as the inner solver:

```fortran
res = newton_krylov(fun, x0, jacobian_product=my_jacobian_product, &
    inner_method='lgmres', inner_check_residual=.true.)
```

This rejects an inner solution before the nonlinear line search or accepted-
step preconditioner update if `||F(x) - J(x)*p||_2` exceeds the requested
inner relative tolerance times `||F(x)||_2`. Comparisons are scaled before
subtraction to avoid intermediate overflow on extreme finite residuals; a
`32*epsilon(1.0_dp)` relative floor accounts for roundoff when the requested
forcing tolerance is below working precision. Failed or nonfinite J*v
callbacks are rejected explicitly. The check applies to all six built-in
methods **and custom `inner_solver=` callbacks**; it can also be combined
with the TFQMR-specific `inner_verify_residual=.true.` (which instead checks
candidate TFQMR stopping points *during* its iteration). Unlike that TFQMR
option, `inner_check_residual` verifies only the final returned inner step
and does not automatically refine an inaccurate inner solution.

This is a stricter SciFort-specific option, **disabled by default** to retain
SciPy's accepted incomplete inner iterates and evaluation counts. It can
reject steps when the inner iteration budget is too small or when numerical
J*v approximations are inaccurate. A finite-difference check adds one residual
evaluation for the global two-point stencil, two for a centered three-point
stencil, or more for coordinatewise differences; the extra calls count toward
`nfev`. An analytic `jacobian_product=` adds no residual evaluations.

The outer Krylov controller also avoids squaring extremely large or tiny
residual norms when updating its forcing term, uses a scaled Euclidean norm
near the floating-point exponent limits, and detects a zero Newton step by
checking its components rather than squaring a subnormal vector. These are
numerical safeguards, not changes to the default finite-difference scheme.
Regression: `test/test_optimize_root_krylov_strict_inner.f90`.

**SciFort extension -- diagonal inverse-Jacobian preconditioning:**
`root(method='krylov', ...)` and `newton_krylov(...)` accept
`inner_m_diag=`, a rank-one real array containing the **diagonal of an
approximation to J(x) inverse**, not the diagonal of J itself. Its entries
act as a right-preconditioner `M*v = inner_m_diag*v` (for the appropriate
built-in Krylov method) and are static during one nonlinear solve:

```fortran
! For a diagonally dominated Jacobian with diagonal d(:):
! Construct approximate inverse entries in user code, excluding zeros.
res = root(fun, x0, method='krylov', inner_method='lgmres', &
    inner_maxiter=40, inner_restart=4, inner_m_diag=invdiag, &
    jacobian_product=my_jacobian_product)
```

This removes the need to write an `inner_m=` callback for a static diagonal
preconditioner. All six built-in Krylov methods, restarted GMRES/LGMRES,
analytic or numerical J*v, and a custom `inner_solver=` are supported.
`inner_m_diag` must have `size(x0)` finite, nonzero entries. For built-in
MINRES, all entries must be **strictly positive**, since its preconditioner
must be positive definite; other methods can accept signed entries.
It cannot be combined with `inner_m=` or `inner_m_update=`. Invalid
configurations are rejected before residual evaluation. Multiplication
checks for unrepresentably large products before they can overflow.
The outer convergence norm remains the norm of the **original residual**;
row scaling and solution scaling are not implied by this option. A diagonal
approximation is useful when it captures the main scaling or coupling of
the Jacobian, but does not generally improve a poorly conditioned dense
system. If the approximate inverse must change after a nonlinear step,
continue using `inner_m=` with `inner_m_update=` instead.

This is an optional **SciFort extension**, not a SciPy 1.17.0 keyword.
Omitting it leaves solver defaults and evaluation counts unchanged.
Regression: `test/test_optimize_root_krylov_diagonal_preconditioner.f90`;
benchmark: `benchmarks/benchmark_krylov_diagonal_preconditioner.f90`.

**SciFort extension -- automatically updated diagonal preconditioning:**
`root(method='krylov')` and `newton_krylov()` also accept
`inner_m_diag_auto=.true.`. Instead of supplying `inner_m_diag=`, SciFort
estimates `diag(J(x))` by applying the matrix-free Jacobian to each unit
coordinate vector, and then uses a safeguarded inverse diagonal as a
preconditioner. With `jacobian_product=` the estimation uses **n analytic Jv
calls per refresh** and no extra residual calls. With finite differences it
requires at least **n additional residual evaluations** for two-point
sampling or **2n** for centered three-point sampling (n = `size(x0)`),
plus the ordinary Krylov solve costs. This is an **O(n)** number of Jv
probes per refresh; for large matrix-free systems it can be prohibitively
expensive. No full Jacobian is assembled.

```fortran
res = root(fun, x0, method='krylov', inner_method='lgmres', &
    inner_restart=4, inner_maxiter=40, &
    inner_m_diag_auto=.true., inner_m_diag_update_interval=2, &
    jacobian_product=my_jacobian_product)
```

`inner_m_diag_update_interval=` defaults to `1`: estimate at the first
nonlinear step and again at steps `2, 3, ...`. An interval `k>1` refreshes
at steps `1, 1+k, 1+2*k, ...`; `0` estimates at the initial point only.
`inner_m_diag_floor=` defaults to `sqrt(epsilon(1.0_dp))` and must be finite
and in `(0,1]`. The absolute floor for nonzero estimated diagonal entries
is the larger of this fraction times the largest absolute diagonal entry
and the smallest positive normal `real(dp)` value. Entries that are
**exactly zero** are mapped to a preconditioner factor of `1` rather than
inverted. Very small **nonzero** entries are clipped to the floor while
preserving their signs; built-in MINRES takes their absolute values to keep
its diagonal preconditioner positive definite. Positive definiteness of
the *preconditioner* does not establish that the user-provided Jacobian is
symmetric; MINRES still requires a symmetric linear operator.

The automatic mode cannot be combined with `inner_m_diag=`, `inner_m=`, or
`inner_m_update=`. Invalid controls are rejected before the first residual
call. A failed Jv probe or a nonfinite inverse-diagonal action causes the
solve to report failure rather than using a partially estimated diagonal.
All auto-specific controls require `inner_m_diag_auto=.true.`. The mode is
optional and absent-option SciPy-compatible defaults remain unchanged.

**SciFort extension -- checked rank-one coupling correction:**
`root(method='krylov')` and `newton_krylov()` accept
`inner_m_diag_rank_one=.true.` together with `inner_m_diag_auto=.true.` and
**exact** diagonal estimation (`inner_m_diag_probes` omitted or zero).
By default, the opt-in model is

`J = diag(a) + u * transpose(ones)`.

For a known, fixed right coupling vector `v`, also pass
`inner_m_diag_rank_one_right=v`, a rank-one `real(dp)` array of length
`size(x0)`. This extends the model to

`J = diag(a) + u * transpose(v)`.

The right vector may have mixed signs or zero components, but must have
at least one nonzero entry, finite values, and the correct length. Supplying
it without `inner_m_diag_rank_one=.true.` is an input error. The solver
normalizes `v` by its largest absolute entry, which makes its scale
irrelevant and avoids overflow for very large finite weights. The default
right factor is `ones`, preserving the original mean-field interface.

The estimator first obtains the unclipped, exact Jacobian diagonal and
one sign-weighted directional product `J*sign(v)`. For the normalized
right vector and `s=sum(abs(v))`, it fits each component by

`u(i) = ( (J*sign(v))(i) - J(i,i)*sign(v(i)) ) / (s - abs(v(i)))`,

followed by `a(i) = J(i,i) - u(i)*v(i)`. Two deterministic sign-vector
products then check the candidate model with a normalized discrepancy
threshold of `1.0e-5`. For three or more variables, duplicate or sign-reversed
check directions are replaced deterministically; this does not make the
checks exhaustive. If the checks pass, SciFort applies

`P(z) = z/a - (u/a) * dot_product(v, z/a) / (1 + dot_product(v, u/a))`

without assembling the Jacobian. Application costs O(n), and fitting
costs **three additional Jv calls per diagonal refresh**, with O(n)
additional storage. For global forward/centered finite differences,
these calls cost three/six additional residual evaluations.

```fortran
res = root(fun, x0, method='krylov', inner_method='gmres', &
    inner_m_diag_auto=.true., inner_m_diag_rank_one=.true., &
    inner_m_diag_rank_one_right=v, jacobian_product=my_jvp, &
    inner_maxiter=20)
```

If the fitted uncoupled diagonal is singular or too small, there is too
little off-diagonal information to identify `u` (e.g. a right factor with
only one nonzero entry), the Sherman-Morrison denominator is near
singular, arithmetic may be unsafe, or the verification disagrees,
SciFort **falls back to the existing safeguarded diagonal preconditioner**.
A failed or nonfinite Jv call instead causes solver failure. The model is
refreshed with `inner_m_diag_update_interval`, and positive
`inner_m_diag_probes` remain unsupported. The built-in MINRES solver rejects
the option because this inverse can be nonsymmetric; in one dimension the
existing diagonal preconditioner is used.

Verification of the supplied-right model is **heuristic**, not proof that
the Jacobian has this structure. Untested directions can still have
significant errors, so use that approach when the right coupling pattern
is known from the model. It is a SciFort extension, not a SciPy keyword.
With the options omitted, preconditioner behavior and work counts are
unchanged. Regression tests:
`test/test_optimize_root_krylov_rank_one.f90` and
`test/test_optimize_root_krylov_weighted_rank_one.f90`.

**SciFort extension -- discover an unknown rank-one right factor:**
For a potentially dense Jacobian whose coupling is approximately rank one,
set `inner_m_diag_auto=.true., inner_m_diag_rank_one=.true.,
inner_m_diag_rank_one_discover=.true.` without supplying
`inner_m_diag_rank_one_right`. For example:

```fortran
res = root(fun, x0, method='krylov', inner_m_diag_auto=.true., &
    inner_m_diag_rank_one=.true., inner_m_diag_rank_one_discover=.true., &
    jacobian_product=my_jvp)
```

This **opt-in, full-Jacobian diagnostic** applies the Jacobian to all `n`
coordinate directions, stores its `n` columns, and reuses those products
for exact automatic diagonal scaling. It identifies two off-diagonal
pivot positions `p /= q` and reconstructs `u` and `v` from
`J(i,j) = u(i)*v(j)` for `i /= j`. Independent rows/columns determine the
missing pivot entries. Unlike the known-right method, it checks **every
off-diagonal matrix element**, within a normalized tolerance of `1e-5`,
before building the guarded Sherman-Morrison inverse. Singular fitted
bases or denominators, non-identifiable sparse patterns, incompatible
couplings, and unsafe arithmetic retain the existing diagonal scaling;
nonfinite or failed Jv calls remain fatal.

Setup uses **n Jv products, O(n^2) memory, and O(n^2) fitting work** per
refresh, versus n+3 Jv products and O(n) memory for the known-right
model. Forward/centered finite differences cost n/2n additional residual
perturbations, respectively, before ordinary inner-solver work. It is
**not** a cheap general matrix-free preconditioner: this option can be
valuable when a small but strongly coupled system amortizes the full
Jacobian measurement across subsequent Krylov solves, but has high
memory and setup costs for large systems. In dimensions smaller than
four or when a rank-one fit cannot be identified safely, the algorithm
falls back to diagonal scaling. The option is mutually exclusive with
`inner_m_diag_rank_one_right`, requires exact diagonal estimation, and
shares the rank-one option's incompatibility with built-in MINRES. It
is disabled by default and is not a SciPy API option. Regression tests:
`test/test_optimize_root_krylov_rank_one_discover.f90`.

**SciFort extension -- checked rank-two Woodbury correction:**
For a known two-column right factor `V`, `root(method='krylov')` and
`newton_krylov()` also accept
`inner_m_diag_rank_two=.true., inner_m_diag_rank_two_right=V` with
`V` a real(dp) `(size(x0),2)` array. This targets

`J = diag(a) + U * transpose(V)`

with an unknown `(n,2)` left factor `U`. For example:

```fortran
res = root(fun, x0, method='krylov', inner_method='gmres', &
    inner_m_diag_auto=.true., inner_m_diag_rank_two=.true., &
    inner_m_diag_rank_two_right=right_factors, &
    jacobian_product=my_jvp)
```

Both columns must be finite and nonzero; their scales are normalized
independently. `V` must be supplied explicitly, and the rank-one and rank-two
options cannot be used together. As with rank one, this requires exact
`inner_m_diag_auto=.true.` and rejects randomized/Walsh diagonal probing and
built-in MINRES. The right factors should represent the *known model
structure*, not an arbitrary attempted fit to a general dense Jacobian.

Two bounded sign-vector Jacobian products identify the two unknown
coefficients in each row of `U`, using the exact `J(i,i)` and equations

`(J*p)(i) - J(i,i)*p(i) = sum_k U(i,k) * (dot_product(V(:,k),p) - V(i,k)*p(i))`.

The estimator selects two fitting directions by maximizing the minimum
normalized per-row determinant over deterministic candidate sign vectors.
Two additional independent directions verify the fitted matrix to a
normalized discrepancy threshold of `1.0e-5`. If verified, preconditioning
uses the Woodbury identity:

`P(z) = D^(-1)*z - D^(-1)*U * (I + V^T*D^(-1)*U)^(-1) * V^T*D^(-1)*z`,

where `D = diag(a)`. This uses **four additional Jv products** per refresh
beyond `n` exact diagonal-coordinate products, or eight additional residual
evaluations for centered differences. Storage and application cost O(n).

If `n < 4`, right columns are dependent or not identifiable row-by-row,
`D` is singular, the two-by-two Woodbury factor is nearly singular, any
floating-point safety bound fails, or the structural checks disagree, the
original safeguarded diagonal preconditioner is retained. A failed Jv
callback instead stops the solve. Verification is heuristic, not a guarantee
of accurate inversion in unseen directions. The new options are disabled
by default and **are not SciPy keywords**. Test coverage:
`test/test_optimize_root_krylov_rank_two.f90`.

**SciFort extension -- sketched rank-two setup (opt-in):**
Add `inner_m_diag_rank_two_sketch=.true.` to the same `root` or
`newton_krylov` call to try a faster setup when the known `V` is reliable.
This option requires `inner_m_diag_rank_two=.true.` and inherits the
rank-two restrictions above. It is disabled unless explicitly requested.

With at least five variables, the sketch first selects three bounded sign
vectors whose moment coefficients give a nonsingular, well-conditioned
3-by-3 linear system *for every row*. In the structured model, three
Jacobian-vector products determine `D(i,i)`, `U(i,1)`, and `U(i,2)`:

`(J*p)(i) = D(i,i)*p(i) + U(i,1)*dot_product(V(:,1),p) + U(i,2)*dot_product(V(:,2),p)`.

The sketch reconstructs `J(i,i)`, then invokes the existing checked
rank-two estimator. Its two fitting and two independent verification
products give **seven Jv calls per accepted refresh**, independent of
problem dimension, instead of `n+4`. For finite differences this is seven
(two-point) or fourteen (three-point) residual evaluations, on top of
ordinary nonlinear-solver work. The added sign-vector search is O(n) in
memory and O(n) in work with a fixed 12-vector candidate family.

If row geometry is degenerate, the approximation fails structural checks,
or its guarded Woodbury inverse is unsafe, the solver returns to the
original exact n-coordinate diagonal estimation and rank-two verification.
This retry may cost *more* than the original path on incompatible systems.
A failed Jacobian callback aborts the solve, as in the existing exact
path. The structural test is heuristic, not a certificate of low rank;
use the sketch only for a model with known right factors. The fast path
does not change the default or the behavior of the existing
`inner_m_diag_rank_two` option. Regression:
`test/test_optimize_root_krylov_rank_two_sketch.f90`; reproducible
comparison: `benchmarks/benchmark_krylov_woodbury.f90`.

**SciFort extension -- discover two unknown rank-two factors (opt-in):**
For nonlinear Jacobians with unknown but identifiable two-column coupling,
set `inner_m_diag_auto=.true., inner_m_diag_rank_two=.true.,
inner_m_diag_rank_two_discover=.true.` and **omit**
`inner_m_diag_rank_two_right`. This is a SciFort-specific option, not part of
SciPy's `root` or `newton_krylov` interface:

```fortran
res = root(fun, x0, method='krylov', inner_method='gmres', &
    inner_m_diag_auto=.true., inner_m_diag_rank_two=.true., &
    inner_m_diag_rank_two_discover=.true., jacobian_product=my_jvp)
```

The algorithm probes each coordinate once to assemble the Jacobian and
reuse its exact diagonal estimate. It searches for an independent
*off-diagonal* 2-by-2 anchor, estimates both factors `U`, `V` without
using diagonal entries to fit them, and recovers the masked anchor rows
and columns using independent off-diagonal observations. It then checks
**all off-diagonal entries** against `U*transpose(V)` to a normalized
relative tolerance of `1e-5`, and obtains the uncoupled diagonal `D` from
the difference between the measured diagonal and the fitted coupling.
For safe fits it uses the guarded inverse

`P(z) = D^(-1)*z - D^(-1)*U*(I + V^T*D^(-1)*U)^(-1)*V^T*D^(-1)*z`.

The method requires **at least six variables** and independent row/column
observations. Rank-deficient, sparsely identifiable, contaminated, or
poorly conditioned coupling, singular fitted `D` or the 2-by-2 Woodbury
factor, and unsafe arithmetic all **fall back to the existing exact
Jacobian-diagonal inverse**. Failed Jacobian-vector callbacks instead
stop the solve. A full off-diagonal check guards against incompatible
models, but agreement within tolerance is not a general certification
of Jacobian structure or floating-point accuracy.

Setup uses **n Jacobian-vector products**, O(n^2) memory for the full
Jacobian, and O(n^2) additional fitting/checking work. With directional
forward or centered differences, it normally costs n or 2n residual
perturbations respectively (in addition to the usual nonlinear work).
This option may be useful for small, strongly coupled systems where
preconditioner construction is amortized; for large matrix-free systems
it often defeats the purpose of Krylov methods. It is disabled by default,
mutually exclusive with explicitly supplied rank-two right factors and
`inner_m_diag_rank_two_sketch`, and still requires exact automatic
diagonal estimation (`inner_m_diag_probes=0`). The built-in MINRES solver
is not supported. Regression:
`test/test_optimize_root_krylov_rank_two_discover.f90`.

**SciFort extension -- streamed unknown rank-two discovery (opt-in):**
To use the same exhaustive off-diagonal model check without allocating
an `n`-by-`n` Jacobian, add `inner_m_diag_rank_two_stream=.true.` to the
unknown-factor discovery options above:

```fortran
res = root(fun, x0, method='krylov', inner_m_diag_auto=.true., &
    inner_m_diag_rank_two=.true., &
    inner_m_diag_rank_two_discover=.true., &
    inner_m_diag_rank_two_stream=.true., jacobian_product=my_jvp)
```

At most eight evenly spaced coordinate columns are retained as pilot
vectors to identify a well-conditioned, disjoint off-diagonal 2-by-2
pivot. Each remaining coordinate direction is evaluated once; two
uncertain rows and columns are cached in O(n) memory until the factors
can be completed. Every off-diagonal entry is compared against the
rank-two model, including entries outside the pilot set and the missing
anchor rows and columns. If a fit is rejected, the same n products
provide the existing exact diagonal preconditioner. Analytic callbacks
need exactly n Jv calls per refresh (n or 2n finite-difference residual
perturbations), without a second checking pass. Computational setup
work remains O(n^2) because all n-by-n off-diagonal entries are checked,
but peak *discovery storage* is O(n), rather than O(n^2). The Krylov
iteration has its own separate storage requirements.

This is an optional **memory reduction, not a reduction in Jacobian
products**. Unlike the full-Jacobian anchor search, streamed discovery
can fail to find two independent right columns when the pilot set is
unrepresentative. That case safely retains diagonal scaling; the
streaming option is not a guaranteed low-rank detector. Other guards,
including singular Woodbury rejection, IEEE finite checks and callback
failure handling, remain in force. The option requires
`inner_m_diag_rank_two_discover=.true.` and is invalid for non-Krylov
methods, but does not change the default discovery algorithm or its
API. The public keyword is a SciFort extension, absent from SciPy.
Regressions: `test/test_optimize_root_krylov_rank_two_stream.f90` and
`test/test_optimize_root_krylov_rank_two_stream_large.f90`.

**SciFort extension -- four-product symmetric rank-two discovery (opt-in):**
When the uncoupled diagonal `a` is **already known** and the Jacobian
remainder is symmetric and rank two, pass
`inner_m_diag_rank_two_symmetric_base=a` with automatic rank-two discovery:

```fortran
res = root(fun, x0, method='krylov', inner_m_diag_auto=.true., &
    inner_m_diag_rank_two=.true., &
    inner_m_diag_rank_two_discover=.true., &
    inner_m_diag_rank_two_symmetric_base=a, &
    jacobian_product=my_jvp)
```

Here `a` is the diagonal of the *uncoupled* matrix `D` in
`J = D + L`, **not** `diag(J)` or the diagonal of an inverse
preconditioner. This method assumes `L = transpose(L)` and
`rank(L) = 2`. It does not require the coupling factors or a transpose
Jacobian callback. Two deterministic orthonormal directions `S` produce
`Y = (J-D)S`. When the 2-by-2 matrix `B = transpose(S)*Y` is
nonsingular, the symmetric rank-two identity
`L = Y*inverse(B)*transpose(Y)` reconstructs the factors using O(n)
workspace. Two additional held-out, orthogonal directions test the
reconstruction. The existing safeguarded Woodbury inverse is used only
if the checks pass. On a failed structural check, the implementation
reverts to **n exact coordinate Jacobian products** to obtain the
ordinary diagonal preconditioner. A failed Jacobian evaluation is
propagated as a solver error, not treated as a rejected model.

For accepted analytic-Jacobian fits, discovery uses **four Jv products**
per refresh, independent of dimension. With standard two-point numerical
products it requires four residual perturbations; centered three-point
products require eight. It uses O(n) discovery storage and O(n)
reconstruction/verification arithmetic. The accepted inverse action
is O(n). Very small systems (n < 6), singular or nonfinite base
diagonals, poorly conditioned interpolation matrices, and unsafe
Woodbury matrices cannot use the new fast fit. A base diagonal supplied
for a nonlinear Jacobian must remain valid at each preconditioner
refresh; `inner_m_diag_update_interval=0` can freeze the initial fit.

**Limit:** The two verification directions are deterministic heuristic
checks, not an exhaustive certificate that an arbitrary Jacobian is
symmetric or rank two. In particular, an untested perturbation can
escape detection. Use this method only with independently justified
symmetry and diagonal structure, not as a general-purpose automatic
rank-two detector. The new option conflicts with streamed discovery
and with supplied rank-two right factors; the ordinary discovery mode
remains the default when this keyword is absent. The existing MINRES
restriction on rank-two preconditioning still applies. The option is a
SciFort extension, not part of SciPy. Tests:
`test/test_optimize_root_krylov_rank_two_symmetric.f90` and
`test/test_optimize_root_krylov_rank_two_symmetric_large.f90`.

**SciFort extension -- six-product nonsymmetric rank-two discovery (opt-in):**
When the uncoupled diagonal `a` of `D` is known, but the remainder
`L=J-D` is **nonsymmetric** with rank two, callers can provide analytic
transpose-Jacobian-vector products in addition to normal products:

```fortran
res = root(fun, x0, method='krylov', inner_m_diag_auto=.true., &
    inner_m_diag_rank_two=.true., &
    inner_m_diag_rank_two_discover=.true., &
    inner_m_diag_rank_two_nonsymmetric_base=a, &
    jacobian_product=my_jvp, &
    jacobian_transpose_product=my_jtv)
```

`my_jtv(x, vec, product, ok)` has the same signature as `my_jvp`, but
computes `transpose(J(x))*vec`. The **Jacobian**, evaluated at the
current nonlinear iterate, must match the matrix used for the ordinary
Jacobian-vector products. The uncoupled diagonal `a` is **not**
`diag(J)`. The transpose callback is required; `jacobian_product` itself
is optional (if absent, forward products use the requested finite-
difference stencil). The provided base diagonal must remain accurate
at every refresh.

**Evolving uncoupled diagonal (SciFort extension):** For nonlinear systems
where `D` varies with `x`, use the alternate optional callback rather than
passing a constant array:

```fortran
res = root(fun, x0, method='krylov', inner_m_diag_auto=.true., &
    inner_m_diag_rank_two=.true., inner_m_diag_rank_two_discover=.true., &
    inner_m_diag_rank_two_nonsymmetric_base_fun=base_diagonal, &
    inner_m_diag_update_interval=1, &
    jacobian_product=my_jvp, jacobian_transpose_product=my_jtv)
```

The user-defined `base_diagonal` routine has the explicit interface

```fortran
subroutine base_diagonal(x, diagonal, ok)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: diagonal(:)
    logical, intent(out) :: ok
end subroutine base_diagonal
```

The public interface `vector_root_base_diagonal` can be imported from
`scifort_optimize_root`. The output is the *uncoupled* `D(x)`, not
`diag(J(x))`; the callback must return finite, nonzero entries and `ok=.true.`.
SciFort calls it at the current accepted iterate **once per scheduled
preconditioner refresh**, before evaluating fitting Jacobian products. By
default `inner_m_diag_update_interval=1` refreshes it every nonlinear step;
`0` evaluates it only at the initial step and deliberately holds that
preconditioner fixed even if `D(x)` changes. For positive intervals `k`,
the initial preconditioner is refreshed again every `k` Newton iterations,
not on a final convergence-only check.

**Adaptive diagonal refresh (SciFort extension, opt-in):**
`inner_m_diag_update_adaptive=.true.` retains the update interval as a
**maximum age**, but refreshes earlier when the previous accepted step
reduces the Euclidean residual norm insufficiently. The optional
`inner_m_diag_update_stall_ratio=` accepts a finite value in `(0,1]`
(default `0.8`); an accepted step with `||F(x_new)||_2 >= ratio*||F(x_old)||_2`
requests a refresh at the next Newton iteration. After every refresh, the
maximum-age clock restarts. Unsuccessful trial line-search points and final
convergence-only checks do not trigger a refresh. Without an explicit
interval, adaptive mode defaults to a maximum age of four iterations;
fixed mode still defaults to an interval of one. Adaptive mode requires
`inner_m_diag_auto=.true.` and an interval greater than zero; a custom
stall ratio requires adaptive mode. This is a *progress-based heuristic*,
not an estimate of Jacobian change or a guarantee of reduced work.
It is available for all automatically estimated diagonal preconditioners,
including rank-one and rank-two extensions, and through `newton_krylov`.
For example:

```fortran
res = root(fun, x0, method='krylov', inner_m_diag_auto=.true., &
    inner_m_diag_update_interval=4, &
    inner_m_diag_update_adaptive=.true., &
    inner_m_diag_update_stall_ratio=0.85_dp, &
    jacobian_product=jvp)
```

**Work-aware adaptive refresh (SciFort extension, opt-in):** With adaptive
refresh enabled, supply `inner_m_diag_update_work_growth=g` (finite `g > 1`)
to watch the number of preconditioner applications in each inner solve.
The count at the last diagonal refresh is the baseline. After a later
accepted Newton step, request an early refresh if

```text
(inner applications - estimated rebuild cost) / baseline >= g
```

The optional nonnegative integer `inner_m_diag_update_rebuild_cost=c`
(default `0`) expresses the estimated cost of rebuilding in equivalent
inner preconditioner applications. It requires the work-growth option;
that option in turn requires `inner_m_diag_update_adaptive=.true.`.
Both options are appended to the existing Fortran argument lists, so
previous positional calls remain compatible. The existing residual-stall
and maximum-age refresh conditions continue to apply independently.

```fortran
res = root(fun, x0, method='krylov', inner_m_diag_auto=.true., &
    inner_m_diag_update_adaptive=.true., &
    inner_m_diag_update_interval=6, &
    inner_m_diag_update_work_growth=2.0_dp, &
    inner_m_diag_update_rebuild_cost=4, &
    jacobian_product=jvp)
```

**Measured automatic-preconditioner rebuild cost (SciFort extension, opt-in):**
For any automatic diagonal or low-rank Krylov preconditioner, add
`inner_m_diag_update_measured_cost=.true.` to use the **observed number**
of attempted Jacobian-vector products performed by the current preconditioner
setup, including probes, fitting, structural verification, retries, and
exact-diagonal fallback. Covered methods include the ordinary exact or
randomized diagonal, known-right rank one and rank two, rank-two sketch,
dense or streamed unknown-factor discovery, and symmetric/nonsymmetric
low-product discovery. The optional finite positive
`inner_m_diag_update_transpose_cost=c` (default `1.0_dp`) weights each
measured transpose product relative to a forward product. The cost is

```text
cost = forward_fit_and_fallback_Jv + c * transpose_fit_and_verify_Jv
```

This substitutes for `inner_m_diag_update_rebuild_cost=`; specifying
both is an error. Measured mode requires `inner_m_diag_auto=.true.`,
`inner_m_diag_update_adaptive=.true.`, and the work-growth threshold.
An analytic transpose callback is needed **only** by the nonsymmetric
low-product discovery method, not by other measured-cost modes.

```fortran
res = root(fun, x0, method='krylov', inner_m_diag_auto=.true., &
    inner_m_diag_rank_two=.true., &
    inner_m_diag_rank_two_discover=.true., &
    inner_m_diag_rank_two_nonsymmetric_base_fun=base_diagonal, &
    inner_m_diag_update_interval=8, &
    inner_m_diag_update_adaptive=.true., &
    inner_m_diag_update_work_growth=1.2_dp, &
    inner_m_diag_update_measured_cost=.true., &
    inner_m_diag_update_transpose_cost=4.0_dp, &
    jacobian_product=jvp, jacobian_transpose_product=jtv)
```

Measurements are of **attempted products**, not elapsed time: a fallback
from an unsuccessful rank-two fit includes exact-diagonal re-estimation
(`n` coordinate Jv calls unless exact values were already cached).
A rejected sketch can additionally retry rank-two fitting with that exact
diagonal; these fitting products are counted too. The fit may make fewer
products if verification terminates early. Finite-difference Jacobian-vector
evaluations count as one directional product apiece, not per residual
evaluation: two-point uses one residual evaluation and three-point two.
Base-diagonal callbacks are not costed, and setup work is compared
to *preconditioner applications* during the inner solver, not actual
Jacobian-vector products in that solver. The comparison therefore remains
an approximate work-equivalence, particularly for non-GMRES or custom
solvers. The option adds no Jacobian products and keeps default behavior
unchanged. See `test/test_optimize_root_krylov_measured_rebuild.f90`,
`test/test_optimize_root_krylov_measured_setup_modes.f90`, and
`benchmarks/benchmark_krylov_work_refresh.f90`.
**CPU-timed Krylov preconditioner refresh (SciFort extension, opt-in):**
The `inner_m_diag_update_timed_cost=.true.` option compares the CPU time
used by the current inner linear solve with the CPU time measured for the
first inner solve immediately after the most recent preconditioner rebuild.
The refresh is requested early when

```text
(current_inner_cpu_seconds - last_rebuild_cpu_seconds) /
    post_rebuild_inner_cpu_seconds >= inner_m_diag_update_work_growth
```

The policy uses the standard Fortran `cpu_time` intrinsic. It measures
*process CPU time*, not wall time; concurrent threads, scheduling, and
clock resolution can influence the result. The default
`inner_m_diag_update_timed_min_seconds=0.001_dp` requires a measurable
post-refresh baseline and a measurable current solve. The threshold can
be tuned to the workload and must be finite and nonnegative. Unavailable,
nonmonotone, or sub-resolution measurements never trigger an early refresh;
residual-stall checks and the maximum refresh interval remain active.

```fortran
res = root(fun, x0, method='krylov', inner_m_diag_auto=.true., &
    inner_m_diag_update_adaptive=.true., &
    inner_m_diag_update_interval=6, &
    inner_m_diag_update_work_growth=2.0_dp, &
    inner_m_diag_update_timed_cost=.true., &
    inner_m_diag_update_timed_min_seconds=0.001_dp, &
    jacobian_product=jvp)
```

The timed option requires adaptive work-growth scheduling and an automatic
preconditioner. It is mutually exclusive with manual rebuild cost,
`inner_m_diag_update_measured_cost=.true.`, and
`inner_m_diag_update_transpose_cost`. It also works with `newton_krylov`.
It does not add Jacobian products, and it does not alter default behavior.
Timing measurements can be noisy, particularly for very fast linear solves;
this is a heuristic scheduling option, not a guarantee of speedup. See
`test/test_optimize_root_krylov_timed_refresh.f90` for the forwarding,
invalid-input, and early-refresh regressions.

**Wall-clock Krylov refresh (SciFort extension, opt-in):**
Set `inner_m_diag_update_wall_clock=.true.` together with
`inner_m_diag_update_timed_cost=.true.` to use elapsed seconds from the
standard Fortran `system_clock` intrinsic instead of process CPU seconds.
The same setup and inner-solve timings, work-growth decision, and
`inner_m_diag_update_timed_min_seconds` minimum apply. The default remains
CPU time. The wall-clock switch without the timed mode is invalid; all
existing manual and product-count cost exclusions still apply. A non-Krylov
method also rejects the option.

```fortran
res = root(fun, x0, method='krylov', inner_m_diag_auto=.true., &
    inner_m_diag_update_adaptive=.true., &
    inner_m_diag_update_interval=6, &
    inner_m_diag_update_work_growth=2.0_dp, &
    inner_m_diag_update_timed_cost=.true., &
    inner_m_diag_update_wall_clock=.true., &
    inner_m_diag_update_timed_min_seconds=0.001_dp, &
    jacobian_product=jvp)
```

Wall-clock intervals are computed from an `integer(int64)` counter and
positive counter rate. Unavailable readings, a decreasing counter (including
rollover), nonfinite durations, and intervals too short to meet the minimum
do not trigger early refresh. Fixed maximum-age and residual-stall refreshes
remain available. Wall time reflects thread scheduling and machine load; it
is generally more directly relevant to elapsed-runtime optimization with
OpenMP callbacks, but remains noisy for short inner solves. It is not a
speedup guarantee, and the Fortran standard does not promise monotonic
`system_clock` behavior on all implementations. See
`test/test_optimize_root_krylov_wall_clock_refresh.f90` and
`benchmarks/benchmark_krylov_quadrature_refresh.f90`.

The expanded 132-case optimized benchmark compares manual costs with measured
weights of `1` and `16`. For the six 128-variable nonlinear cases, measured
1x averages 46.2 forward and 6.7 transpose products per solve, versus
71.8 and 4.0 with no refresh. Measured 16x avoids early rebuilds on these
cases (71.8 and 4.0). Neither setting is globally optimal; the weights
are cost-model inputs, not empirical runtime measurements.

For built-in GMRES, one preconditioner application generally corresponds
to an Arnoldi iteration. Other inner solvers can apply the preconditioner
multiple times per iteration; a custom solver may make no such calls, in
which case no work signal is available and only the existing conditions
control refresh. The count excludes preconditioner setup products and
independent strict residual checks. The rebuild cost is a *user estimate*;
this heuristic does not learn timing costs or guarantee a performance gain.
It uses constant additional memory, adds no Jacobian-vector products
solely for scheduling, and leaves defaults unchanged. See
`test/test_optimize_root_krylov_work_adaptive.f90`.

Interval choice trades fitting and
verification products against inner-solver work and depends on the rate of
Jacobian change and the relative cost of transpose actions. See
`benchmarks/benchmark_krylov_nonlinear_refresh.f90` and its checked
60-configuration CSV for example workloads; no interval is optimal in general.
A failed callback or invalid
entry returns unsuccessful solver status with a descriptive message; an
inaccurate but finite diagonal that fails the structural checks triggers
the existing exact-diagonal fallback. This callback works with both `root`
and `newton_krylov`, with fixed or adaptive held-out transpose verification.
Do not supply the callback and `inner_m_diag_rank_two_nonsymmetric_base`
together. All existing fixed-array calls and defaults remain unchanged.
It uses O(n) persistent storage and adds no Jacobian-vector products beyond
those already required for the selected verification policy. Regression:
`test/test_optimize_root_krylov_rank_two_nonlinear_base.f90`.

For `L=J-D`, the algorithm uses two orthonormal right probes `Q` to
form `C=L*Q` and two orthonormal left probes `P` to form
`R=transpose(L)*P`. If `B=transpose(P)*C` is numerically nonsingular,
it reconstructs `L=C*inverse(B)*transpose(R)`. A bilinear check compares
`transpose(P)*C` and `transpose(Q)*R`, and two further forward probes
orthogonal to `Q` check the fitted coupling. Only after verification
and finite/conditioning checks does the existing guarded rank-two
Woodbury inverse apply. Rejected fits use the **original n-coordinate
exact-diagonal** estimator, not the supplied uncoupled diagonal; a
failed forward or transpose callback instead terminates the solve.

An accepted refresh uses **four forward Jv plus two transpose Jv
products**, independent of n, and O(n) additional storage. For additional
adjoint consistency checking, set
`inner_m_diag_rank_two_nonsymmetric_verify_transpose=.true.` alongside
the nonsymmetric base. This tests two independent held-out `J^T v`
directions, orthogonal to the transpose fitting probes, against the
reconstructed model. An accepted setup then uses **four forward and four
transpose products**, while the default remains unchanged. Rejected fits
use the exact Jacobian diagonal; failed held-out adjoint callbacks abort
the solve. The fixed-check Boolean option is also forwarded by `newton_krylov`.

**Adaptive adjoint verification (SciFort extension):** As an alternative
specify the positive, finite ratio of one transpose-product cost to
one forward-product cost:

```fortran
inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost=16.0_dp
```

This selects at least **one** held-out transpose direction. If two
transpose products have a modeled cost no greater than one quarter of
the `n` forward products needed for exact diagonal discovery (`2*c <= n/4`),
it selects **two**; otherwise it selects one. It always selects two if
the normalized interpolation determinant is at most `1.0e-4` or if
forward verification discrepancy is at least 10% of its acceptance
limit. These thresholds are heuristics, not accuracy guarantees.
Accepted setups thus cost four forward and either **three or four**
transpose products, depending on fit and cost. Incorrect adjoints
detected by the selected checks cause exact-diagonal fallback; failed
adjoint callbacks terminate the solve. The adaptive option is available
in both `root` and `newton_krylov`, conflicts with specifying the
fixed-check Boolean (even `.false.`), and rejects zero, negative,
NaN, and infinite ratios before evaluating the residual. Omitting both
options retains the original two fitting transpose products without
additional held-out checks. Each forward
numerical product costs one residual perturbation for two-point or two
for centered three-point differences. Analytic transpose calls do not
count as residual function evaluations (`nfev`). This option is mutually
exclusive with `inner_m_diag_rank_two_symmetric_base`, the streamed
unknown-factor mode, rank-two sketching and caller-supplied rank-two
right factors. It is disabled by default, does not support built-in
MINRES, and is a **SciFort extension**, not part of SciPy.

**Limit:** The two held-out directions are heuristic checks, not a
proof of low rank. Rank-three perturbations can be missed if they are
invisible to the chosen probes. Poorly conditioned projection matrices,
incorrect diagonal assumptions and unsafe Woodbury inverses cause
fallback. Even two extra adjoint checks do not exhaustively certify
rank two or callback consistency on all directions. Use this path only
when the rank-two model is supported by problem structure. Regression
programs: `test/test_optimize_root_krylov_rank_two_nonsymmetric.f90`,
`test/test_optimize_root_krylov_rank_two_nonsymmetric_large.f90`,
`test/test_optimize_root_krylov_rank_two_transpose_verify.f90`, and
`test/test_optimize_root_krylov_rank_two_transpose_adaptive.f90`.

The reproducible `benchmarks/benchmark_krylov_rank_two_discovery.f90`
compares exact-diagonal, dense, streamed, low-product, fixed-adjoint-checked
low-product and adaptive-cost policies (`c=1`, `c=16`) on 32/128-variable
rank-two and rank-three systems
with 2/6 decades of diagonal scaling. Its CSV includes forward and
transpose counts, measured CPU time, solution error, and modeled work
`njvp + c*njtv` for `c=1,4,16` (relative adjoint cost, **not** seconds).
For a single accepted rank-two setup, low-product work is `4 Jv + 2 J^Tv`
versus `4 Jv + 4 J^Tv` with fixed enhanced verification and either
`4 Jv + 3 J^Tv` or `4 Jv + 4 J^Tv` adaptively; the complete-solve
counts include Krylov iterations and refreshes. A rejected rank-three fit
can cost *more* than ordinary diagonal discovery because it performs
structural checks before exact-diagonal fallback.

**SciFort extension -- lower-cost structured diagonal probes:**
`inner_m_diag_probes=p` optionally estimates `diag(J)` from `p` Walsh-Hadamard
sign-vector products, rather than the default `n` coordinate products. The
probe vectors have entries `+1` or `-1`, and the estimated diagonal is
`sum(z * (J*z)) / p`. The estimator is **deterministic** (it neither uses nor
modifies the Fortran random-number generator) and stores O(n) working data.

```fortran
res = root(fun, x0, method='krylov', inner_method='lgmres', &
    inner_m_diag_auto=.true., inner_m_diag_probes=8, &
    inner_m_diag_update_interval=0, jacobian_product=my_jvp)
```

Omitting the argument, or setting it to zero, retains the previous exact
coordinate-by-coordinate estimator and evaluation counts. A positive count
`p >= n` also uses the exact method. For `1 <= p < n`, estimation requires
`p` analytic Jacobian-vector calls, or normally `p` extra two-point / `2p`
extra centered three-point residual evaluations. With
`jacobian_step_mode='componentwise'`, each dense probe itself costs up to
`n` coordinate evaluations, so analytic or global-directional Jv actions are
preferable for large systems.

**Accuracy is structure-dependent.** In particular, for a power-of-two
probe count, Walsh signs exactly cancel off-diagonal contributions between
coordinates with distinct indices modulo that count. A sufficiently narrow
banded Jacobian can therefore yield an exact diagonal with few probes, but
more distant or dense couplings can alias into a biased estimate. Non-power-
of-two counts are accepted but generally lose the exact cancellation
property. The existing small/zero-diagonal safety rules still apply to the
*estimated* entries and cannot correct a biased estimate. Start with 8 or
16 probes on locally coupled systems, and compare against the exact mode
before relying on the approximation.

This is an original SciFort feature, not a SciPy 1.17.0 option. See
`test/test_optimize_root_krylov_diagonal_probes.f90` and
`benchmarks/benchmark_krylov_diagonal_probes.f90` for regression and
performance comparisons.


**SciFort extension -- seeded and checked diagonal probes:**
When `inner_m_diag_probes=p` is positive and smaller than `size(x0)`, set
`inner_m_diag_probe_scheme='rademacher'` to use locally generated seeded
pseudorandom `+1`/`-1` vectors instead of Walsh columns. An optional integer
`inner_m_diag_probe_seed=` (default `1`) reproduces the sign sequence across
runs without changing the state of Fortran's `random_number`. Different
seeds generate different estimates. For ideal independent Rademacher signs,
`average(z*(J*z))` is an unbiased diagonal estimator; SciFort uses a
reproducible integer pseudorandom sequence, **not** independent truly random
samples. Sampling variance can be large on dense or ill-conditioned systems.
The preexisting `walsh` scheme is the default.

To guard against obvious probe errors, `inner_m_diag_verify_samples=k`
(default `0`, disabled) spot-checks up to `k` evenly spaced Jacobian diagonal
entries through exact coordinate Jv calls. `inner_m_diag_verify_rtol=t`
(default `0.25`, must be finite in `(0,1]`) controls the maximum normalized
spot-check discrepancy. In the `rademacher` scheme with at least two probes,
the guarded path also compares the first-half probe estimate against the
full estimate at every coordinate; this costs **no additional Jv products**.
If either check disagrees beyond `t`, the entire diagonal is recomputed by
the exact `n`-product method before it can affect a Krylov iteration. A
successful spot-check does **not** certify unexamined coordinates: the
reported full-solve error still depends on conditioning.

```fortran
res = root(fun, x0, method='krylov', inner_method='lgmres', &
    inner_m_diag_auto=.true., inner_m_diag_probes=16, &
    inner_m_diag_probe_scheme='rademacher', &
    inner_m_diag_probe_seed=17, &
    inner_m_diag_verify_samples=8, &
    inner_m_diag_verify_rtol=0.15_dp, &
    jacobian_product=my_jvp)
```

The extra verification costs between zero and `k` coordinate Jv calls
before any fallback. A failed spot-check or probe-batch check adds `n`
exact Jv calls, so a guarded approximation can be *more expensive* than
choosing exact estimation initially. With finite differences, each Jv
may require one/two additional residual evaluations for the global
2-point/3-point stencil; coordinatewise differencing may be much more
expensive. Use the exact method when reliable diagonals matter more than
operator cost. A zero probe count or `p >= n` always uses the original exact
method. New options require `inner_m_diag_auto=.true.`; a seed requires
`rademacher`, and a custom verification tolerance requires positive `k`.
All inputs are checked before the first residual evaluation.

In the recorded 2,048-variable irregular-coupling benchmark, 8 Walsh probes
had mean relative preconditioner-diagonal error about `0.309` and 8 seeded
Rademacher probes about `0.140`, but the latter's **maximum** relative
error was about `1.35` (versus `0.310` for Walsh). The verified variants
fell back to exact estimation, using more Jv calls than the exact mode.
This is a warning against treating lower *average* error as a guarantee.
See `test/test_optimize_root_krylov_random_diagonal.f90` and
`benchmarks/benchmark_krylov_random_diagonal.f90`. The algorithm and API
are original SciFort extensions, not SciPy 1.17.0 keywords.

**SciFort extension -- progressive verified Rademacher probes:**
`inner_m_diag_probe_max=q` optionally permits a checked estimate that fails
at the initial `inner_m_diag_probes=p` budget to use additional seeded
Rademacher vectors, growing the total budget roughly as `p, 2*p, 4*p, ...`
without repeating previously evaluated Jacobian-vector products. The budget
never exceeds `q`. After each stage, the algorithm checks disagreement with
an earlier probe group; if this passes, it checks the selected exact diagonal
coordinates. An estimate that still fails at `q` falls back to the original
coordinatewise exact estimator. This can avoid an expensive exact fallback
on dense or irregular couplings when a modest increase in probes suffices.

```fortran
res = root(fun, x0, method='krylov', inner_method='lgmres', &
    inner_m_diag_auto=.true., inner_m_diag_probes=8, &
    inner_m_diag_probe_scheme='rademacher', &
    inner_m_diag_probe_seed=17, inner_m_diag_verify_samples=4, &
    inner_m_diag_verify_rtol=0.25_dp, inner_m_diag_probe_max=32, &
    jacobian_product=my_jvp)
```

The cap requires `inner_m_diag_auto=.true.`, the `rademacher` scheme,
positive verification samples and `1 <= p < n`, and must satisfy
`p <= q <= n`. Invalid combinations are rejected **before evaluating the
residual**. Setting `q=p` retains the preexisting fixed-budget path; omitting
`q` is the default and does not change the old estimator. For progressive
stages, an online average is used to avoid excessive intermediate sums; it
can give a slightly different rounded estimate from a single fixed-budget
run. The estimator has O(n) workspace. Progressive coordinate verification
starts with the evenly spaced `k` samples and **rotates the sample positions**
within those strata as the probe budget grows, using reproducible base-2
bit-reversed offsets (0, 1/2, 1/4, 3/4, ... of the stratum width). If a rotated
position was checked earlier, progressive verification instead selects the
next unchecked coordinate, with cyclic wraparound. A local path-compressed
successor table keeps this search efficient even for `k` near `n`; if all
coordinates have been checked, the stage needs no additional spot Jv calls.
Previously checked coordinates are revalidated against each expanded diagonal
estimate using cached exact Jacobian-vector products. This avoids a blind spot
when an offset rounds to zero for narrow strata and guarantees new spatial
coverage while unchecked coordinates remain. Verification is still heuristic:
accepted estimates can have significant errors at unexamined coordinates.

With `s` stages whose split-batch check passes, the algorithm makes at most
`min(n, s*k)` **unique** exact coordinate spot calls, instead of at most `k`
in the previous implementation. If it must fall back to an exact diagonal,
all cached products are reused, so the worst-case total remains `q + n`
Jacobian-vector products: previous spot evaluations count toward the `n`
coordinate products rather than adding to them. Progressive probing can still
cost more than exact estimation, and some stage expansions will require more
spot calls than before. Omitting `q`, or setting `q=p`, retains the earlier
fixed-budget evaluation behavior and public APIs.

A dense, nonlocally coupled 1,024-variable test used 1,032 analytic Jv
calls with fixed eight-probe verification and 36 with an eight-to-32
progressive budget. However, the progressive diagonal had **maximum
coordinatewise relative error about 0.333**: successful spot and batch
verification cannot guarantee accuracy at unexamined coordinates. Use the
exact estimator or a validated preconditioner when stronger guarantees are
needed. See `test/test_optimize_root_krylov_progressive_diagonal.f90`,
`benchmarks/benchmark_krylov_progressive_diagonal.f90`, and
`benchmarks/krylov_progressive_diagonal_gfortran14_o2.csv`.

On a constructed coupled nonlinear system with 128 variables and six orders
of diagonal scaling, automatic initial-only preconditioning converged with
146 Jv calls, refresh-every-step with 530, and a supplied static diagonal
with 18; unpreconditioned LGMRES stalled at 1,551 calls. These are benchmark
cases, not performance guarantees. The automatic interface, floor policy,
and update schedule are original SciFort extensions rather than SciPy API
keywords. See `test/test_optimize_root_krylov_auto_diagonal.f90` and
`benchmarks/benchmark_krylov_auto_diagonal.f90`.

**SciFort extension -- memory-bounded restarted GMRES and LGMRES:** Built-in
`root(method='krylov', inner_method='gmres', ...)`,
`root(method='krylov', inner_method='lgmres', ...)`, and `newton_krylov(...)`
accept `inner_restart=` as an optional *positive integer*:

```fortran
res = root(fun, x0, method='krylov', inner_method='gmres', &
    inner_maxiter=80, inner_restart=8, &
    jacobian_product=my_jacobian_product)
```

When `inner_restart < inner_maxiter`, GMRES and LGMRES use at most
`inner_restart` **new** Krylov directions per Arnoldi cycle. `inner_maxiter` becomes the
total number of Arnoldi steps allowed **across all cycles**, not the number
of cycles. The unpreconditioned linear residual `rhs - J(x)*solution` is
recomputed with a matrix-free Jacobian-vector product after each cycle,
and a fresh Arnoldi basis is built from that residual. This adds one
Jacobian-vector product per cycle (and the corresponding residual function
calls with numerical finite differences). It avoids storing an
`O(n*inner_maxiter)` Arnoldi basis: the reusable basis and preconditioned
search-direction workspace scale as `O(n*inner_restart)` instead.
Right-preconditioners and `inner_orthogonalization=` are supported.

When the keyword is omitted, or `inner_restart >= inner_maxiter`, the
**original one-cycle GMRES/LGMRES implementation is used unchanged**. Invalid
nonpositive sizes, methods other than built-in GMRES/LGMRES, custom
`inner_solver=`, and
non-Krylov outer methods reject this option before residual evaluation.
Like the original implementation, exhausting the linear iteration budget
returns the best finite correction obtained; set
`inner_check_residual=.true.` to require the nonlinear controller to reject
steps that fail its true linear-residual requirement. Restarted GMRES/LGMRES
may converge more slowly or stagnate if the restart width is too small.

For **LGMRES**, the retained `outer_k` augmentation vectors remain available
throughout all restart cycles of one linear solve. They are screened against
the current Jacobian before the first cycle and the final accumulated Newton
step is considered for the retained history only after the last cycle.
Each restart checks the actual unpreconditioned residual `rhs - J(x)*step`;
its last Jacobian-vector product also supplies the candidate image for
Jacobian-image pruning without a duplicate evaluation. Consequently, with
restarting enabled `outer_candidate_image_reuse` does not allocate a raw
Arnoldi-product cache; the final true-residual check makes it unnecessary.
Image pruning, `outer_image_reuse`, preconditioning, and all existing
reorthogonalization modes continue to work. The memory bound concerns the
Arnoldi workspace and does not include retained `outer_k` vectors.

This is an *optional SciFort extension*, not a SciPy 1.17.0 keyword.
Regressions: `test/test_optimize_root_krylov_gmres_restart.f90` and
`test/test_optimize_root_krylov_lgmres_restart.f90`.

**SciFort extension -- adaptive restarted LGMRES:** For built-in
`inner_method='lgmres'` with `1 <= inner_restart < inner_maxiter`, the
optional `inner_restart_adaptive=.true.` lets the next Arnoldi cycle grow
when the **true unpreconditioned linear residual** falls by less than 20%
in the preceding cycle. The width grows by at most a factor of two per weak
cycle, up to `inner_restart_max=`. If no maximum is supplied, the ceiling is
`min(inner_maxiter, 4*inner_restart)`. The accumulated Newton correction and
retained augmentation vectors are not discarded when the width grows.
The total `inner_maxiter` budget is unchanged, and the solver retains the
existing early exit on exact stagnation when no further growth is possible.

```fortran
res = root(fun, x0, method='krylov', inner_method='lgmres', &
    inner_maxiter=48, inner_restart=4, &
    inner_restart_adaptive=.true., inner_restart_max=16, &
    jacobian_product=my_jacobian_product)
```

This is an opt-in heuristic, not a guarantee of faster convergence. The
steady-state Arnoldi storage is bounded by the selected maximum restart
width, rather than the initial width; workspace resizing can transiently
allocate additional storage. Enabling adaptation without a genuinely
restarted built-in LGMRES solver is rejected before evaluating the residual.
`inner_restart_max=` also requires adaptation, must be at least the initial
restart width, and cannot exceed `inner_maxiter`. The fixed-width default,
non-restarted solves, and built-in GMRES are unchanged. This adaptive policy
is **not** a SciPy keyword. Regression:
`test/test_optimize_root_krylov_adaptive_restart.f90`.

Benchmark: `benchmarks/benchmark_krylov_restarted_gmres.f90`.

**SciFort extension -- Arnoldi reorthogonalization:** Both
`root(method='krylov', ...)` and `newton_krylov(...)` accept
`inner_orthogonalization=` when `inner_method='gmres'` or `'lgmres'`:

```fortran
res = newton_krylov(fun, x0, inner_method='gmres', &
    inner_orthogonalization='selective', jacobian_product=my_jacobian_product)
```

- `'once'` (default) uses one pass of modified Gram-Schmidt (MGS), retaining
  SciFort's previous SciPy-compatible Krylov solver and work counts.
- `'selective'` applies a second MGS pass when the norm of the projected
  vector is at most half its original norm (DGKS orthogonality criterion).
- `'twice'` unconditionally applies two MGS passes.

The second-pass inner products are **added** to the existing Arnoldi
Hessenberg coefficients, and the projected-vector norm is recomputed before
normalization and happy-breakdown detection. No extra `J*v` or `fun` calls
occur, although the extra dot products and vector updates cost additional
floating-point operations. Reorthogonalization may improve the basis for
ill-conditioned or nearly linearly dependent Krylov vectors, but does not
ensure overall nonlinear convergence or fix an inadequate inner iteration
budget. Case-insensitive names are accepted; unknown names, non-Krylov
methods, the four non-Arnoldi inner methods and custom `inner_solver=` callbacks
reject this option before evaluating `fun`. This is a SciFort extension,
**not** a SciPy 1.17.0 keyword. The original default remains unchanged.

**SciFort extension -- prune redundant LGMRES augmentation:**
`root(method='krylov', ...)` and `newton_krylov(...)` accept the optional
real scalar `outer_prune_tol=` with built-in `inner_method='lgmres'` and
`outer_k > 0`:

```fortran
res = newton_krylov(fun, x0, inner_method='lgmres', outer_k=10, &
    outer_prune_tol=1.0e-8_dp, inner_orthogonalization='selective')
```

By default `outer_prune_tol=0` disables filtering and preserves the original
SciPy-compatible FIFO augmentation history, iterates and work counts. With a
positive tolerance less than 1, the normalized new Newton step is projected
onto a *temporary* orthonormal basis of previously retained solution vectors
using two modified Gram-Schmidt passes. If the remaining component has norm
at most `outer_prune_tol`, the new direction is considered redundant and is
not retained. This detects dependence on a **combination** of old vectors,
not just nearly parallel pairs. New independent directions use the original
FIFO replacement policy when `outer_k` is full; no existing direction is
orthogonalized or changed in place. Filtering requires no additional
Jacobian-vector or residual function calls. The solution-space span is built
incrementally during the append-only part of a solve, reducing the screening
work from roughly O(n*k**2) to O(n*k) per proposed direction with k retained
vectors; after history capacity is reached, the original reconstructed-span
path is used to handle FIFO eviction. To avoid large additional allocations,
the incremental basis is enabled only when n*outer_k does not exceed
1,048,576 `real(dp)` values (8 MiB in double precision). Larger problems
use the original implementation. The existing `outer_prune_tol=` interface
and default behavior are unchanged.

**SciFort extension -- fast full-history FIFO pruning:** for solution-space
pruning, the optional `outer_fifo_prune_fast=.true.` enables QR-factor
updating after the history is full. The solver removes the oldest retained
column from a small upper-triangular factor, applies Givens rotations to
restore triangularity, and updates the matching orthonormal basis. It then
projects the proposed new direction, using two-pass modified Gram-Schmidt.
The original normalized augmentation vectors remain unchanged except for
ordinary FIFO replacement. No additional residual or Jacobian-vector calls
are needed, and the option is disabled by default.

```fortran
res = root(fun, x0, method='krylov', inner_method='lgmres', &
    outer_k=12, outer_prune_tol=1.0e-8_dp, &
    outer_fifo_prune_fast=.true.)
```

This opt-in path requires built-in LGMRES, `outer_k >= 2`, and positive
`outer_prune_tol`; invalid combinations are rejected before the first
residual evaluation. Its extra orthonormal basis is limited to the existing
8 MiB workspace threshold; larger problems continue to use the previous
reconstruction. Near-singular histories also fall back to the previous
method rather than attempting an unstable downdate. The QR factor is
rebuilt periodically to limit roundoff accumulation. Floating-point
trajectories can differ slightly from legacy pruning, so performance and
convergence gains are problem dependent. The option is a SciFort extension,
not a SciPy keyword.

A positive tolerance requires built-in LGMRES with `outer_k > 0`; values
outside `[0,1)` or nonfinite values are rejected before evaluating `fun`.
It is invalid with other inner solvers (including user `inner_solver=`) or
non-Krylov root methods. The filter measures dependence in the *solution*
space, not in the Jacobian image space; very small but important directions
can be discarded if the threshold is set too aggressively. It does not
replace reorthogonalization or the rank-revealing least-squares fallback.
This is a SciFort extension, not a SciPy 1.17.0 keyword. Regressions:
`test/test_optimize_krylov_augmentation.f90` and
`test/test_optimize_root_krylov_augmentation.f90`.

**SciFort extension -- prune in current Jacobian-image space:**

For built-in LGMRES, set `outer_image_prune_tol=` to a positive real in
`(0, 1)` to detect a different kind of numerical redundancy:

```fortran
res = root(fun, x0, method='krylov', inner_method='lgmres', &
    outer_k=5, outer_image_prune_tol=1.0e-8_dp, &
    jacobian_product=my_jvp, inner_orthogonalization='selective')
```

At each nonlinear iterate, SciFort evaluates the *current* Jacobian images
`J(x)*v` of the stored augmentation vectors. A two-pass modified
Gram-Schmidt independence test removes vectors whose normalized images
have distance at most `outer_image_prune_tol` from the span of earlier
images. It also tests the new Newton direction before retaining it. The
stored vectors remain in **solution coordinates**; neither the Jacobian
images nor the orthonormalized basis replaces the original LGMRES directions.
No image is reused at a new nonlinear point, where the Jacobian may differ.
The remaining independent directions retain their original order and FIFO
replacement behavior.

For an established retained history, the current-Jacobian filtering pass
incrementally extends a temporary orthonormal image basis rather than
rebuilding it for every candidate. That scratch basis is reset before each
nonlinear iterate, and the new Newton step uses the ordinary rank-checking/FIFO insertion
path by default; optional `outer_image_basis_reuse` uses the retained basis. The separate low-level
`retain_lgmres_image_direction` API is unchanged. The streaming low-level
helper is append-only (no eviction within a pass), requires an explicit
`reset_lgmres_image_pruning_stream` before each new Jacobian, and is intended
for library implementations screening an existing history. No new public
`root` keyword or change in function-evaluation counts is involved.

This filtering can reject vectors that are independent in solution space
but nearly dependent under the Jacobian. It costs up to `outer_count + 1`
additional Jacobian-vector products per nonlinear step (`outer_count` for
reexamining stored history and one for the new candidate). With two-point
finite differences that entails one extra residual evaluation per product;
three-point finite differences take two, and componentwise differencing
can take more. The `nfev` counter includes all extra residual samples.
With analytic `jacobian_product=`, no extra *residual* evaluations are
needed, though the analytic callback is called more often.

The feature is disabled by default (`outer_image_prune_tol=0`), preserving
SciPy-compatible work counts, and requires built-in LGMRES with `outer_k>0`.
It cannot be used together with positive `outer_prune_tol=` (which
filters in solution space), a custom inner solver, or a non-Krylov method.
Nonfinite values and thresholds outside `[0, 1)` are rejected before the
initial residual is evaluated; invalid image callbacks fail the solve.
The threshold is a *relative directional* cutoff, not an absolute bound on
Jacobian singular values. Too-aggressive pruning and inaccurate numerical
Jacobian products can discard useful augmentation directions. Prefer an
analytic Jacobian-vector callback or conservative thresholds on ill-scaled
problems. Regressions: `test/test_optimize_krylov_augmentation_images.f90`
and `test/test_optimize_root_krylov_augmentation_images.f90`.

**SciFort extension -- reuse the current-Jacobian image basis:**

With positive `outer_image_prune_tol=` or `outer_image_prune_auto=.true.`,
`outer_image_basis_reuse=.true.` avoids rebuilding the image-space orthonormal
span for the candidate Newton direction after Arnoldi:

```fortran
res = root(fun, x0, method='krylov', inner_method='lgmres', &
    outer_k=8, outer_image_prune_auto=.true., &
    outer_image_basis_reuse=.true., jacobian_product=my_jvp)
```

The original filtering pass constructs a two-pass modified Gram-Schmidt
basis at the *current* nonlinear iterate. The opt-in path preserves that
basis and tests the new direction against it. When the history is full, the
candidate is tested **before** removing the oldest retained direction,
maintaining FIFO semantics. The cached basis is invalidated after one
candidate and reset whenever `x` changes. Reuse does not skip any Jacobian
products or residual evaluations; it saves the orthogonalization work of
one candidate insertion per nonlinear iteration.

The default is `.false.`; incompatible option combinations are rejected
before any residual evaluation. Reuse requires built-in LGMRES and fixed or
automatic current-Jacobian image pruning. When `outer_image_reuse=.false.`,
the solver uses the already allocated raw-image array as the orthonormal
basis workspace, since Arnoldi does not need the raw images in this case.
No second `n*outer_count` basis allocation or additional 8 MiB cap is needed;
if incremental cache state is unexpectedly invalid, raw images are rebuilt
before the original rank-revealing insertion fallback is used.

When `outer_image_reuse=.true.`, raw images must remain intact for Arnoldi;
therefore a separate cached basis remains capped at 1,048,576 `real(dp)`
entries (8 MiB). Larger cases use the original reconstruction method.
Both paths preserve numerical defaults and residual/Jacobian evaluation
counts in tested cases. Whole-solver timings may be dominated by residual
and Jacobian work, so no end-to-end speedup is guaranteed. The low-level
streaming routines and original reconstruction procedures remain available.
Tests: `test/test_optimize_krylov_augmentation_candidate_basis.f90`,
`test/test_optimize_krylov_image_basis_storage.f90`, and
`test/test_optimize_root_krylov_image_basis_reuse.f90`.
This is a SciFort-specific extension, not a SciPy keyword.

**SciFort extension -- automatic Jacobian-image pruning:**

For built-in LGMRES, `outer_image_prune_auto=.true.` enables an opt-in
*projection-dependent* tolerance without requiring the user to select
`outer_image_prune_tol=`:

```fortran
res = root(fun, x0, method='krylov', inner_method='lgmres', &
    outer_k=5, outer_image_prune_auto=.true., &
    jacobian_product=my_jvp, inner_orthogonalization='selective')
```

Each candidate Jacobian image is normalized before two-pass modified
Gram-Schmidt. Let `r` be the remaining norm after the first pass and let
`n` be the number of components. The relative rejection threshold is

`min(0.25, max(32*epsilon(1.0_dp)*sqrt(n), sqrt(epsilon(1.0_dp))*sqrt(max(0,1-r**2))))`.

Thus, nearly parallel images (small `r`) receive a larger cutoff, while
nearly orthogonal images receive only a floating-point roundoff floor. The
cutoff is recomputed for each old image and the new direction, using the
*current* Jacobian at every nonlinear iteration; it is not a global fixed
numerical rank. It is a conservative heuristic rather than a bound on
Jacobian conditioning or finite-difference truncation error. In particular,
noisy Jacobian-vector callbacks can require an explicit fixed threshold
instead. Unlike solution-space pruning, automatic image-space pruning
requires the same extra Jacobian-vector products as fixed image filtering.

The default `.false.` leaves SciPy-compatible behavior and evaluation
counts unchanged. Automatic pruning requires built-in LGMRES with
`outer_k > 0`; it cannot be combined with positive `outer_image_prune_tol`
or `outer_prune_tol`, or with a custom inner solver. Explicitly passing
`outer_image_prune_tol=0` is permitted. Invalid options are checked before
initial residual evaluation. The direct and public regressions are
`test/test_optimize_krylov_augmentation_auto.f90` and
`test/test_optimize_root_krylov_augmentation_auto.f90`. The comparison
regression also includes automatic filtering as a fourth policy.

**SciFort extension -- reuse current-iterate pruning images in Arnoldi:**

Set `outer_image_reuse=.true.` with either `outer_image_prune_auto=.true.`
or positive `outer_image_prune_tol=`:

```fortran
res = root(fun, x0, method='krylov', inner_method='lgmres', &
    outer_k=5, outer_image_prune_auto=.true., &
    outer_image_reuse=.true., jacobian_product=my_jvp)
```

Image filtering evaluates each retained `J(x)*v` at the current nonlinear
iterate before Arnoldi. With reuse enabled, the same already-validated images
are passed to Arnoldi directly, avoiding its redundant operator calls for
retained directions. The corresponding solution-space directions are preserved
*exactly*, rather than renormalized again, to ensure cached images correspond
to the very same input vectors. The candidate step's separate image check is
not removed. Cache storage is temporary and is discarded on return from each
nonlinear inner solve; all retained images are recomputed on the next iterate,
when `J(x)` may have changed.

This can save up to one Jacobian-vector product per retained direction per
nonlinear step (one or two finite-difference residual probes for 2-point or
3-point directional differences, respectively). Analytic Jacobian callbacks
do not affect `nfev`, but are called fewer times. For deterministic functions,
reuse avoids mathematically redundant linear-operator actions, though tiny
changes in floating-point normalization can change nonlinear trajectories on
very ill-conditioned problems. Callbacks with side effects or state-dependent
results should leave reuse disabled. With no stored directions (the first
step), there is nothing to reuse.

The default `.false.` keeps previous residual/Jacobian evaluation counts and
SciPy-compatible behavior unchanged. `.true.` requires built-in LGMRES with
active fixed or automatic Jacobian-image pruning; invalid combinations are
rejected before evaluating the residual. Tests:
`test/test_optimize_root_krylov_image_reuse.f90` and
`test/test_optimize_root_krylov_pruning_comparison.f90`.

**SciFort extension -- reuse Arnoldi products for the proposed LGMRES step:**

`outer_candidate_image_reuse=.true.` avoids a separate Jacobian-vector call
when evaluating the new augmentation candidate in Jacobian-image space:

```fortran
res = root(fun, x0, method='krylov', inner_method='lgmres', &
    outer_k=5, outer_image_prune_auto=.true., &
    jacobian_product=my_jvp, outer_image_reuse=.true., &
    outer_candidate_image_reuse=.true.)
```

At a fixed nonlinear iterate, an analytic Jacobian product is linear in its
second argument. If Arnoldi constructs a step `p = sum_i y_i*z_i`, then its
image `J*p = sum_i y_i*(J*z_i)` can be reconstructed from the products already
calculated before Gram-Schmidt. The raw-product workspace is allocated **only** with this option. The
reconstructed image is used for the *candidate's pruning test*, not to skip a requested independent
true-residual check. As the filter normalizes its image directions, it need
not multiply the reconstructed normalized-step image by the possibly extreme
right-hand-side norm. No operator images survive across Newton iterations.

**Memory use:** With both `outer_image_reuse` and
`outer_candidate_image_reuse` enabled, the candidate reconstruction reads
retained columns directly from the already validated current-iterate image
cache. Only the *new* Krylov columns need a separate raw-product workspace,
with shape `(n, inner_maxiter)` instead of `(n, inner_maxiter + outer_count)`.
The image filter also compacts retained directions in place using one
`n`-element candidate rather than a second `(n, outer_k)` array. For a
saturated history the two changes avoid about
`8*n*(outer_k + outer_count - 1)` bytes of transient real64 storage, before
allocator overhead. Both changes are internal and require no new arguments.
For robustness, reconstruction continues to use the actual products saved
*before* Arnoldi orthogonalization: reconstructing them instead from rounded
Hessenberg coefficients can change the convergence path of ill-conditioned
systems. See `benchmarks/benchmark_lgmres_cache_memory.f90`.

**SciFort extension -- bound candidate-image cache memory:**

For built-in LGMRES with analytic `jacobian_product`, Jacobian-image
pruning, and `outer_candidate_image_reuse=.true.`, the optional integer
`outer_candidate_cache_max_mib=` specifies a maximum capacity, in mebibytes,
for the separate raw Arnoldi-product array. The default `0` means unlimited,
matching the original candidate-image reuse behavior. Positive values are
validated before the first residual evaluation and require candidate-image
reuse; negative values and use with other root methods are rejected.

```fortran
res = root(fun, x0, method='krylov', inner_method='lgmres', &
    outer_image_prune_auto=.true., jacobian_product=my_jvp, &
    outer_image_reuse=.true., outer_candidate_image_reuse=.true., &
    outer_candidate_cache_max_mib=32)
```

Before each inner LGMRES solve, SciFort checks the worst-case raw-product
capacity, `8*n*inner_maxiter` bytes when retained images are shared, or
`8*n*(inner_maxiter+outer_count)` bytes otherwise. If the projected capacity
exceeds the positive limit, it does **not** allocate the raw-product cache.
Instead it evaluates `J(x)*p` directly for the candidate Newton direction
using one additional analytic callback. This is a memory-for-computation
tradeoff, not a different numerical approximation. Prior raw storage is
released if a growing history causes the budget to be exceeded later.

The limit applies **only to the separate raw-image cache**, not to total
process RSS, Arnoldi bases, preconditioner storage, or retained vectors.
The final candidate image may differ by rounding from a sum of cached raw
products, but it is evaluated directly and is not reconstructed from the
rounded Hessenberg factors. It also works with `newton_krylov()` and does
not change residual evaluation counts for analytic Jv callbacks. For
reproducible memory benchmarks, see
`benchmarks/benchmark_lgmres_candidate_cache_budget.f90` and the matching
CSV. New regression: `test/test_optimize_root_krylov_candidate_cache_budget.f90`.

**Arnoldi workspace reuse (internal optimization):** GMRES/LGMRES now keep
basis, Hessenberg and small linear-solve arrays across Newton iterations in
one solver-local workspace. The Jacobian-product and reorthogonalization
scratch uses the next (unused) Arnoldi basis column. LGMRES retains
augmentation directions directly in its existing `outer_vec` history and
stores only ordinary Krylov solution directions in its separate `z` basis;
this does not change the selected `inner_method`, callback API, optional
image-cache logic, or arithmetic/evaluation-count semantics. GMRES also reuses
a basis column for its independent final linear-residual check. Workspace
storage is local to each nonlinear solve, not shared module state, so two
separate solves use independent arrays. Benchmark reproduction details are in
`benchmarks/README.md`.

**Pruning workspace reuse (internal optimization):** Image-space and
solution-space augmentation filters can reuse temporary modified
Gram-Schmidt vectors. The built-in LGMRES path keeps its Jacobian-image arrays for the
lifetime of one nonlinear solve and borrows the dormant Arnoldi `z`
workspace for its prefilter candidate direction. Image storage is exactly
sized to the retained history, without an additional candidate column.
The old Arnoldi basis columns can also be reused for orthogonalization. This avoids a second dense retained-basis buffer.
The low-level `retain_lgmres_direction` and
`retain_lgmres_image_direction` routines retain their allocation-per-call
behavior when the optional `workspace=` argument is absent. The latter can
also be given an optional caller-owned `basis_work=` scratch matrix. Neither
`root` nor `newton_krylov` needs a new argument; independent root calls use
independent storage. See `test/test_optimize_krylov_augmentation_workspace.f90`.

The option is deliberately restricted to a user-supplied **analytic**
`jacobian_product=` callback and active fixed or automatic image pruning with
built-in LGMRES (`outer_k > 0`). With finite-difference products, the action
is not exactly linear in the direction owing to truncation and cancellation;
the separate candidate call is retained. A callback with side effects, noisy
or direction-dependent approximations may also break the assumed linearity;
leave candidate-image reuse disabled in that case. It can be enabled
independently of the earlier `outer_image_reuse=` option, or together with it
for maximum savings. Both defaults are `.false.`, preserving previous solver
outputs and evaluation counts when omitted. Invalid combinations are rejected
before evaluating the nonlinear residual.

Tests: `test/test_optimize_root_krylov_candidate_image_reuse.f90` and the
seven-policy analytic comparison in
`test/test_optimize_root_krylov_pruning_comparison.f90`. The
`benchmarks/benchmark_lgmres_image_cache.f90` scaling benchmark covers
32-, 64-, and 128-variable systems with each option separately and combined.

A reproducible five-policy stress comparison is available in
`test/test_optimize_root_krylov_pruning_comparison.f90`, with a GNU Fortran
optimized reference CSV at `benchmarks/lgmres_pruning_gfortran14_o2.csv`.
Its analytic-product runs report *both* `nfev` (residual evaluations) and the
number of `jacobian_product` callbacks. The fifth policy measures automatic
image pruning with reuse, which removes duplicate current-iterate J*v work. The same test includes
finite-difference runs, in which those extra image products increase `nfev`.
Neither filter is guaranteed to improve convergence for every nonlinear
system. See `docs/validation.md` for the chosen stress parameters and results.

**Rank-deficient LGMRES least squares:** If the small, rotated Arnoldi
triangular system has a tiny diagonal (as can occur with nearly dependent
augmentation vectors), the built-in GMRES/LGMRES solvers now use a direct
one-sided Jacobi SVD to find a truncated minimum-norm least-squares solution.
This replaces the earlier normal-equations (`H**T*H`) eigenvalue fallback.
Unlike normal equations, the direct SVD preserves singular directions down to
approximately `max(m,n)*epsilon` relative to the largest singular value;
smaller singular values are treated as numerical rank deficiency. The usual
full-rank triangular backsolve is unchanged. This internal improvement adds
no residual evaluations, does not require BLAS or LAPACK, and introduces no
new public keywords. It can change computed steps for ill-conditioned
systems, so SciPy parity of evaluation counts does not imply bitwise parity
of the iterates. Regressions: `test/test_optimize_krylov_rank_ls.f90` and
`test/test_optimize_root_krylov_rank_ls.f90`.

`test/test_optimize_root_krylov_reorthogonalization.f90` compares all three
modes on a 48-variable diagonal system spanning 12 orders of magnitude and
checks their Arnoldi basis overlaps, work counts, preconditioned LGMRES,
nonlinear reuse of outer augmentation vectors, public wrapper forwarding,
and invalid-option handling.

Right-mode reference comparisons are generated by
`tools/generate_optimize_krylov_right_preconditioner_reference.py`: SciPy
1.17.0 solves the equivalent *unpreconditioned* transformed system `(J*M)*y=b`
and maps back to `M*y`. On deterministic 16- and 64-variable test systems,
this yields 10-11 matrix-vector products and true relative residuals below
`1e-10`; the default SciPy `M=` path can leave a residual around `1e-3` for
the same example. Regression tests:
`test/test_optimize_krylov_right_preconditioner.f90` and
`test/test_optimize_root_krylov_right_preconditioner.f90`.

The inner tolerance is driven by the existing inexact-Newton forcing term.
If the configured inner-iteration limit is reached without a breakdown, the
finite incomplete iterate is returned to the nonlinear controller, as in
SciPy's `KrylovJacobian.solve`. A failed operator application, nonfinite
output, or numerical breakdown yields a failed root result. The name in
`inner_method=` is ignored when a custom `inner_solver=` callback is provided.
The Fortran procedure interface is not a literal Python `LinearOperator` API.

Reference tests: `test/test_optimize_krylov_iterative.f90` and
`test/test_optimize_root_krylov_methods.f90`, with reference generator
`tools/generate_optimize_krylov_iterative_reference.py` (SciPy 1.17.0).

References and reproducible parity checks: `test/test_optimize_root_krylov_hooks.f90`
and `tools/generate_optimize_root_krylov_hooks_reference.py` (SciPy 1.17.0).

The scalar root-finding surface is:

```text
root_scalar(f [, method, bracket, x0, x1, fprime, fprime2, xtol, rtol, maxiter, k])
bisect(f, a, b [, xtol, rtol, maxiter])
brentq(f, a, b [, xtol, rtol, maxiter])
brenth(f, a, b [, xtol, rtol, maxiter])
ridder(f, a, b [, xtol, rtol, maxiter])
toms748(f, a, b [, xtol, rtol, maxiter, k])
newton(f, x0 [, fprime, x1, fprime2, tol, rtol, maxiter])
```

All routines return `root_result`, containing `root`, `iterations`,
`function_calls`, integer `status`, logical `converged`, a human-readable `flag`,
and the method actually used. The public status constants are
`root_status_success`, `root_status_sign_error`,
`root_status_convergence_error`, `root_status_value_error`, and
`root_status_invalid_input`. Unlike Python, invalid brackets and nonfinite
objective values are represented in the result rather than by exceptions.

With no explicit method, `root_scalar` follows SciPy's selection rule: a bracket
selects `brentq`; an initial point plus first and second derivatives selects
Halley; an initial point plus a first derivative selects Newton; two initial
points select secant; and a lone initial point selects Newton with a numerical
derivative from `scifort_differentiate`. Bracketing methods use SciPy's default
`xtol=2e-12`, `rtol=4*epsilon`, and `maxiter=100`; direct `newton` uses SciPy's
`tol=1.48e-8` and `maxiter=50` defaults.

The bisection, Brent, and Ridders numerical kernels are adapted from
roots-fortran 1.0.0 by Jacob Williams under BSD-3-Clause. SciPy 1.17 supplies
the compatible dispatch/result semantics, open methods, and TOMS748
organization under BSD-3-Clause. SciFort keeps the implementation self-contained
rather than adding roots-fortran as a runtime FPM dependency.

`bracket_root(f, xl0 [, xr0, xmin, xmax, factor, maxiter])` supplies a
scalar counterpart to SciPy's `elementwise.bracket_root`. It expands the
starting interval in both directions (within optional bounds) and returns
`bracket_root_result` with `xl`, `xr`, `fl`, `fr`, `nit`, `nfev`, `status`,
and `success`. The default right endpoint is `xl0+1`, or
`xl0 + min((xmax-xl0)/8, 1)` for finite `xmax`; the default factor is 2.
Status 0 is success, -1 an exhausted interval (also available as
`bracket_root_status_limit`), -2 the iteration cap, -3 a nonfinite value,
and -5 invalid input. The default `maxiter` is 1000. Rank-one and rank-two
overloads also bracket multiple initial values, returning arrays of results. This routine
finds brackets; call `root_scalar` with
`bracket=[br%xl, br%xr]` to refine one.

The scalar minimization surface is:

```text
minimize_scalar(f [, method, bracket, bounds, tol, xtol, xatol, maxiter])
bracket_minimum(f, xm0 [, xl0, xr0, xmin, xmax, factor, maxiter])
find_minimum(f, bracket [, xatol, xrtol, fatol, frtol, maxiter, callback])
brent(f [, bracket, xtol, maxiter])
bounded(f, lower, upper [, xatol, maxiter])
golden(f [, bracket, xtol, maxiter])
```

`bracket_minimum` returns `bracket_minimum_result` with `xl`, `xm`, `xr`,
`fl`, `fm`, `fr`, `nit`, `nfev`, `status`, and `success`. It expands from
the initial middle point to obtain an ordered triple with a lower middle
objective value. The default starting endpoints are `xm0 +/- 0.5`, reduced
to 1/16 of the distance to a finite limit, with a growth factor of 2 and
at most 1000 iterations. With bounds, steps approach the finite limit
geometrically. Status 0 means a valid bracket; -1 means a search limit was
reached (`bracket_minimum_status_limit`); -2 is the iteration cap; -3 means
nonfinite values; and -5 means invalid input (checked before objective
evaluation). This is a scalar counterpart to SciPy 1.17
`elementwise.bracket_minimum`. Rank-one and rank-two overloads bracket independent
starting points; optional bounds and endpoints support broadcasting.
Pass
`[br%xl, br%xm, br%xr]` to `brent` to refine the minimizer.

`find_minimum` refines a valid three-point bracket using Chandrupatla's
quadratic-fit/golden-section minimization scheme, adapted from SciPy 1.17
`scipy.optimize.elementwise.find_minimum`. `bracket` must contain exactly
three finite real numbers, whose sorted points have a middle objective
no greater than either endpoint. It may be passed directly from
`bracket_minimum` as `[br%xl, br%xm, br%xr]`. The result type
`find_minimum_result` contains `x`, `f_x`, `xl`, `xm`, `xr`, `fl`, `fm`,
`fr`, `nit`, `nfev`, `status`, `success`.

Convergence occurs when the wider half-bracket width is at most twice
`abs(xm)*xrtol + xatol`, or the three-point function curvature is at
most twice `abs(fm)*frtol + fatol`. Defaults are
`xrtol=sqrt(epsilon(1.0_dp))`, the other tolerances `tiny(1.0_dp)`,
and `maxiter=100`. The optional callback interface is
`subroutine callback(state, stop)` with `type(find_minimum_result),
intent(in) :: state` and `logical, intent(out) :: stop`; it is invoked
at the initial state and after each iteration. Status 0 indicates
convergence, -1 an invalid bracket, -2 iteration exhaustion, -3
nonfinite evaluations, -4 callback termination
(`find_minimum_status_callback`), and -5 malformed input. Invalid
settings are reported in a result rather than raising an exception.
A rank-one overload accepts `(3,n)` starting triples and returns an
`n`-element array of `find_minimum_result`.


#### Rank-one and rank-two elementwise overloads

`bracket_root`, `find_root`, `bracket_minimum`, and `find_minimum` have generic
scalar, rank-one, and rank-two independent-problem interfaces. The rank-one bracketing
input `xl0(:)` or `xm0(:)` selects a result array of size `n`; its optional
`xr0`, `xl0`, `xmin`, and `xmax` accept scalars or length-`n` arrays. `factor`
and `maxiter` remain shared scalars.

`find_root(f, bracket, ...)` with `bracket(2,n)` returns
`find_root_result(n)`; `find_minimum(f, bracket, ...)` with `bracket(3,n)`
returns `find_minimum_result(n)`. Optional `xatol`, `xrtol`, `fatol`, and
`frtol` accept either scalars or length-`n` arrays. Scalar `f` and optional
scalar-state callback procedures are invoked independently for each problem.
The implementation loops over the existing scalar kernels: it does not
vectorize `f`, reuse evaluations between problems, or broadcast arbitrary
higher-rank shapes. Invalid element inputs have individual statuses; invalid
batch argument shapes yield all-invalid results without objective calls.

The rank-one regression and SciPy 1.17 reference generator are
`test/test_optimize_elementwise_arrays.f90` and
`tools/generate_optimize_elementwise_arrays_reference.py`.

The rank-two bracketing overloads take `xl0(m,n)` or `xm0(m,n)` and return
result objects with shape `(m,n)`. The refinement overloads accept
`bracket(2,m,n)` for `find_root` and `bracket(3,m,n)` for `find_minimum`;
results have shape `(m,n)`. Optional real bounds, starting endpoints and
refinement tolerances follow NumPy-style trailing-axis broadcasting to this
grid: scalar `()`, vector `(n)` or `(1)`, and matrix `(m,n)`, `(m,1)`,
`(1,n)`, or `(1,1)`. Optional `factor` and `maxiter` are shared scalars.
A shape mismatch returns invalid-input status for every grid element
without invoking `f`. Each scalar callback runs independently for each
grid point, including at the initial state; work counts are per point.
Temporary expression arguments are supported, including those with
non-unit lower bounds in compiler-generated descriptors. Zero-length axes
produce correctly shaped empty results. Rank-three and higher scalar-callback problem grids are supported;
packed vector-callback **function overloads** are available through rank four;
assumed-rank `*_vectorized_into` subroutines extend them to higher ranks.

`test/test_optimize_elementwise_rank2.f90` checks SciPy 1.17 rank-two
statuses, iterations, values and evaluation counts alongside scalar
consistency, all broadcast forms, callbacks, invalid arguments and empty
grids. Its reference file is generated by
`tools/generate_optimize_elementwise_rank2_reference.py`.


#### Packed vector-callback elementwise solvers

The four functions `find_root_vectorized`, `bracket_root_vectorized`,
`find_minimum_vectorized`, and `bracket_minimum_vectorized` accept
`subroutine f(x, fx)` with rank-one `real(dp)` arrays of equal length.
They evaluate only active candidates in each callback, so the callback must
apply one common elementwise function and cannot rely on fixed lane indices.
The function overloads support **rank-one through rank-four** grids;
the companion assumed-rank subroutines support higher-rank grids. All
return independent status, iteration, and evaluation counts.

For rank-one inputs, bracketing takes `xl0(n)` or `xm0(n)` and refinement
uses brackets `(2,n)` or `(3,n)`. For rank-two inputs, bracketing takes
`xl0(m,n)` or `xm0(m,n)` and refinement uses `(2,m,n)` or `(3,m,n)`.
Results retain the `(m,n)` grid shape. Optional real parameters broadcast
on trailing dimensions from a scalar, `(n)`, `(m,n)`, `(m,1)`, `(1,n)`, or
singleton equivalents. For bracketing these include `factor` as well as
bounds and initial endpoints; `maxiter` remains a common scalar. Refinement
has four broadcastable tolerances and an optional **scalar-state** callback,
which may stop individual lanes independently.

These wrappers flatten the problem grid in Fortran element order and
forward it to the tested rank-one active-lane controller. Rank-three grids
`(m,n,k)` and rank-four grids `(m,n,k,p)` are accepted by all four solvers;
refinement brackets prepend an axis of size 2 for roots or 3 for minima.
Optional real parameters broadcast from scalars or lower-rank arrays by
aligning trailing axes, expanding singleton dimensions. Incompatible shapes
reject the whole grid without evaluating `f`, while invalid per-lane inputs
or nonfinite objective outputs affect only the corresponding lanes. Empty
grids preserve every original dimension without calling `f`.

For example, `find_root_vectorized(f, bracket(2,2,3,4),
xatol=tolerances(4))` returns `find_root_result(2,3,4)`.
The callback continues to receive a **packed rank-one vector** whose size
may shrink during iteration, not a whole n-dimensional array. The
scalar-state stop callback also remains per lane. For ranks five and above, four additional public **subroutines** use
Fortran 2018 assumed-rank inputs and outputs:

| Routine | Input | Caller-provided result |
| --- | --- | --- |
| `find_root_vectorized_into(f, bracket, res, ...)` | `(2, grid...)` | `find_root_result(grid...)` |
| `find_minimum_vectorized_into(f, bracket, res, ...)` | `(3, grid...)` | `find_minimum_result(grid...)` |
| `bracket_root_vectorized_into(f, xl0, res, ...)` | `xl0(grid...)` | `bracket_root_result(grid...)` |
| `bracket_minimum_vectorized_into(f, xm0, res, ...)` | `xm0(grid...)` | `bracket_minimum_result(grid...)` |

`grid...` can be rank 0 through 14 for refinement and rank 0 through 15 for
bracketing. The bracketing starting-point grid must have the same shape as `res`;
refinement inputs have one additional leading bracket axis of extent 2 or 3.
For scalar outputs, rank-one bracket inputs have lengths 2 or 3.
The result argument has `intent(inout)` to work around an internal compiler
error in GNU Fortran 14.2 for default-initialized derived-type assumed-rank
`intent(out)` arguments. It is explicitly initialized to invalid-input
status at the beginning of each call. Invalid dimensions/broadcast shapes
cause no objective calls. Empty grids retain their shape and cost no work.
Noncontiguous inputs/results are supported. The scalar `maxiter`, packed
active-point callbacks, per-lane stop callbacks and trailing-axis optional
broadcast rules match the existing function overloads.

Example with a five-dimensional problem grid:

```fortran
use scifort_optimize, only : dp, find_root_result, find_root_vectorized_into
real(dp) :: intervals(2,2,1,2,1,2), atol(2)
type(find_root_result) :: answers(2,1,2,1,2)
! Set intervals(:,i,j,k,l,m) for each independent root problem.
atol = [1.0e-12_dp, 1.0e-10_dp]
call find_root_vectorized_into(root_batch, intervals, answers, xatol=atol)
```

The assumed-rank API avoids generating eleven additional families of
rank-specific wrappers. Only the `select rank` handling of flattened
inputs/output assignment is rank-specific. Rank-one through rank-four
function calls remain unchanged.

The companion assumed-rank **argument-aware** entry points are:

| Routine | Additional required arguments | Input/result grid |
| --- | --- | --- |
| `find_root_vectorized_args_into(f, bracket, args, res, ...)` | stacked real `args(..)` | `(2, grid...)` / `find_root_result(grid...)` |
| `find_minimum_vectorized_args_into(f, bracket, args, res, ...)` | stacked real `args(..)` | `(3, grid...)` / `find_minimum_result(grid...)` |
| `bracket_root_vectorized_args_into(f, xl0, args, res, ...)` | stacked real `args(..)` | `xl0(grid...)` / `bracket_root_result(grid...)` |
| `bracket_minimum_vectorized_args_into(f, xm0, args, res, ...)` | stacked real `args(..)` | `xm0(grid...)` / `bracket_minimum_result(grid...)` |

The vectorized callback interface is `subroutine f(x, args, fx)`, with
`real(dp), intent(in) :: x(:), args(:,:)` and
`real(dp), intent(out) :: fx(:)`. Each packed callback receives
`size(args,1) = nargs` and `size(args,2) = size(x)`; the parameter
columns are selected **by active lane**, not position in the previous
callback batch. Thus independently converged or invalid lanes do not
misalign other lanes' parameters. The stored parameter stack can have any
rank from 1 to 15: `args(nargs, [broadcast_grid_axes...])`. The axes
after `nargs` right-align to the full result-grid dimensions; each may be
one or equal to the matching grid dimension. For example, `args(2,1,3)`
broadcasts over a `(2,3)` result grid, while `args(2)` broadcasts both
parameters to every problem. Extra `args` do **not** need to be repeated
for the two or three initial abscissae of a refinement problem. A rank-zero
or zero-parameter stack is invalid. Invalid broadcast dimensions or
mismatched input/output grid shapes produce invalid-input results without
calling `f`. As in the original `*_into` APIs, output is `intent(inout)`
with explicit resetting; scalar and noncontiguous result sections work.
The parameter stack consumes an array dimension, so a rank-15 problem
grid cannot also have a full-rank, lane-varying argument stack.

The same four `*_vectorized_args_into` names are now **generic** over
the `args` argument. The original real-valued stacked array continues to
work; alternatively, pass a rank-one array of `type(elementwise_argument)`.
Construct each element with `elementwise_argument_from(value)`; `value` may
be scalar or a real array of rank up to 15. The argument-list constructor
copies data and its shape, preserving the Fortran first-index-fastest order.
Every element right-aligns its *own* dimensions to the result grid, with
singleton-axis broadcasting. Unlike the stack form, no extra parameter
axis is needed in the original argument, so rank-15 problem grids can have
independently varying rank-15 parameters:

```fortran
use scifort_optimize, only : dp, elementwise_argument, &
    elementwise_argument_from, find_root_vectorized_args_into, &
    find_root_result
real(dp) :: bracket(2,2,3), centers(3), offsets(2,1), scales(1,3)
type(elementwise_argument) :: args(3)
type(find_root_result) :: result(2,3)
args = [elementwise_argument_from(centers), &
        elementwise_argument_from(offsets), &
        elementwise_argument_from(scales)]
! Set bracket(:,i,j) to endpoints enclosing the desired root.
call find_root_vectorized_args_into(f, bracket, args, result)
! f(x, packed_args, fx):
! fx = packed_args(3,:)*(x - packed_args(1,:) - packed_args(2,:))
```

This overload applies equally to `find_minimum_vectorized_args_into`,
`bracket_root_vectorized_args_into`, and
`bracket_minimum_vectorized_args_into`. Each callback still receives a
`real(dp)` parameter matrix of shape `(nargs,nactive)`. An empty list,
an uninitialized argument object, or any incompatible dimension is an
invalid input and produces no objective evaluations. A direct call to a
low-level `*_vectorized_into` dispatcher may supply either `args` (stack)
or `args_list`, with `f_args`, but not both. Results are explicitly reset
even for invalid input, and noncontiguous output slices are supported.

A third `*_vectorized_args_into` overload accepts a rank-one list of
`type(elementwise_typed_argument)` rather than a real stack or real
argument list. Build each argument with
`elementwise_typed_argument_from(value)`, where `value` may be any Fortran
intrinsic or derived type (including a scalar, an assumed-shape section,
or an array of rank 15). The original dynamic type and shape are retained;
arguments broadcast independently on the trailing axes as above.

The packed callback has the *different* explicit interface:

```fortran
subroutine f(x, args, fx)
    real(dp), intent(in) :: x(:)
    type(elementwise_typed_argument), intent(in) :: args(:)
    real(dp), intent(out) :: fx(:)
    integer :: k
    do k = 1, size(x)
        select type(v => args(1)%values(k)%value)
        type is (integer)
            fx(k) = x(k) - real(v, dp)
        class default
            error stop 'Unexpected objective argument type'
        end select
    end do
end subroutine f
```

Here `args(i)%values(k)%value` is an allocatable unlimited-polymorphic
scalar preserving the value and dynamic type of objective argument `i` at
active point `k`. Each callback receives exactly `size(x)` scalar values
per argument, including duplicate copies for the two or three initial
refinement endpoints of the same lane. The numerical solver continues to
use the same packed Chandrupatla/bracketing controllers, so invalid or
converged lanes are removed without shifting parameters between problems.
All four overloads are exported by the public `scifort_optimize` module.

Example construction with independently shaped mixed arguments:

```fortran
type(elementwise_typed_argument) :: args(4)
args = [elementwise_typed_argument_from(center_vector), &
        elementwise_typed_argument_from(integer_offset_grid), &
        elementwise_typed_argument_from(.true.), &
        elementwise_typed_argument_from(user_defined_context)]
call find_root_vectorized_args_into(my_typed_callback, brackets, args, results)
```

The callback needs `select type` to interpret heterogeneous values; the
numeric `(nargs,nactive)` callback is still preferred when all inputs are
`real(dp)`. This original typed callback boxes and copies payloads on every
objective call; the indexed alternative below avoids these copies. Invalid
shapes, empty lists, uninitialized payloads and malformed objects are rejected
before objective evaluation; results are explicitly initialized with invalid-
input statuses. Callers may instead pass `typed_args=` and `f_typed=` directly
to the lower-level `*_vectorized_into` dispatchers, not together with numeric
`args` or `args_list`. The original `f(x,fx)` and numeric-stack interfaces
remain source-compatible. See
`test/test_optimize_elementwise_typed_args.f90`,
`test/test_optimize_elementwise_vectorized_arg_list.f90`,
`test/test_optimize_elementwise_vectorized_args.f90`,
`test/test_optimize_elementwise_vectorized_assumed_rank.f90`,
`test/test_optimize_elementwise_vectorized_rank2.f90`, and
`test/test_optimize_elementwise_vectorized_rank34.f90`.

#### Indexed, zero-copy heterogeneous objective arguments

For a performance-sensitive objective with mixed-type arguments, use
`find_root_vectorized_indexed_args_into`, `find_minimum_vectorized_indexed_args_into`,
`bracket_root_vectorized_indexed_args_into`, or
`bracket_minimum_vectorized_indexed_args_into` instead of the copied typed
variants. These four explicitly named procedures accept the same
`elementwise_typed_argument` lists and the same bounds, tolerance, iteration
limit and callback options as their corresponding typed variants. They preserve
all original typed values throughout the solve, without allocating fresh
polymorphic scalars at every objective call.

Their callback has the interface `f(x, args, indices, fx)`:

```fortran
subroutine f(x, args, indices, fx)
    real(dp), intent(in) :: x(:)
    type(elementwise_typed_argument), intent(in) :: args(:)
    integer, intent(in) :: indices(:,:)
    real(dp), intent(out) :: fx(:)
    integer :: k
    do k = 1, size(x)
        select type(v => args(1)%values(indices(1,k))%value)
        type is (integer)
            fx(k) = x(k) - real(v, dp)
        class default
            fx(k) = 0.0_dp
        end select
    end do
end subroutine f
```

`indices` has shape `(size(args),size(x))`. Entry `indices(i,k)` is the
one-based offset into the **original** `args(i)%values` array for active
point `k`. The same source entry may appear repeatedly when broadcasting or
when evaluating two/three initial bracket endpoints. Index columns track
active-lane compaction, so callers must not infer lane identity from `k`.
The original arguments remain read-only and are not reshaped or copied inside
the indexed callback. The integer map is gathered for each evaluation, so
"zero-copy" describes the typed payloads, not all solver temporaries.

For example, with `args=[elementwise_typed_argument_from([1,2,3])]`,
`find_root_vectorized_indexed_args_into(f, brackets, args, results)` can
solve three independent roots. Each callback reads the integer at
`args(1)%values(indices(1,k))%value`, never a boxed copy. All four indexed
solvers are also exported from `scifort_optimize`.

The lower-level `*_vectorized_into` dispatchers accept `typed_args=` together
with `f_indexed=` as an alternative to `f_typed=`. Passing both callback
forms, or combining an indexed typed callback with the numeric callback or
numeric arguments, returns invalid-input statuses without calling a callback.
Regression coverage is in `test/test_optimize_elementwise_indexed_args.f90`.
Internally, the indexed controllers reuse a preallocated integer map for the
current packed lanes. This avoids repeated vector-subscripted integer-map
array expressions but does not change the public callback interface. The map
remains `(nargs,nactive)`; each callback must read only its supplied columns.

A reproducible performance comparison of the copied and indexed typed APIs is
provided in `benchmarks/benchmark_elementwise_typed_args.f90`, with instructions
in `benchmarks/README.md`. Because copied callbacks allocate and duplicate
polymorphic payloads for every active point, the indexed versions are especially
useful for derived types containing substantial allocatable components. The
benchmark deliberately excludes one-time construction of the typed argument
list from both timings; it measures repeated solver calls only. There is no
performance-threshold assertion, as absolute timing depends on the platform.

The minimization routines `minimize_scalar`, `brent`, `bounded`, and `golden`
return `minimize_result`, containing `x`, `fun`, `iterations`,
`function_calls`, integer `status`, logical `success`, a human-readable `flag`,
and the method used. With no explicit method, `minimize_scalar` selects
`bounded` when bounds are present and `brent` otherwise. Brent and golden accept
two starting points, which trigger downhill bracketing, or a valid three-point
bracket. The defaults follow SciPy 1.17: Brent `xtol=1.48e-8`, bounded
`xatol=1e-5`, golden `xtol=sqrt(epsilon)`, with iteration/function-call
accounting matching SciPy on the regression cases. SciFort reports invalid
inputs and nonfinite objectives through status values rather than Python
exceptions.

The bounded Brent kernel is adapted from Jacob Williams' fmin 1.1.2 under
BSD-3-Clause. SciPy 1.17 supplies downhill bracketing, unbounded Brent,
golden-section search, dispatch, defaults, and result/count semantics under
BSD-3-Clause.

## Integration

The pure `lebedev_rule(n,x,w,status)` subroutine supplies all 32 supported
unit-sphere quadrature rules, orders 3-131, with allocatable Cartesian
`x(3,m)` and surface weights `w(m)` summing to `4*pi`. Orders 13, 25, and 27
contain negative weights. See the [Lebedev guide](lebedev.md) for supported
orders, status codes, examples, and numerical limits.

The native real/complex `solve_bvp` subroutine solves two-point boundary-value
ODEs with adaptive fourth-order collocation, unknown parameters, optional
singular endpoint terms and analytic Jacobians, and C1 cubic dense output.
See the [boundary-value solver guide](bvp.md) for callbacks, shapes, statuses,
numerical contracts, and SciPy compatibility limitations.

```fortran
use scifort_integrate
```

The `scipy.integrate` compatibility layer currently provides rank-one
sampled-data quadrature, scalar and vector adaptive integration, nested scalar
multidimensional integration, explicit-state QMC integration, and independent
double-exponential quadrature. The public surface is:

```text
trapezoid(y[, x, dx, status])
cumulative_trapezoid(y[, x, dx, initial, status])
simpson(y[, x, dx, status])
cumulative_simpson(y[, x, dx, initial, status])
romb(y[, dx, status])
fixed_quad(f, a, b[, n, status])
newton_cotes(n, weights, error_coeff[, status])
newton_cotes(rn, weights, error_coeff[, equal, status])
quad(f, a, b, result, abserr[, epsabs, epsrel, limit, neval, status, &
    points, weight, wvar, wvar2, limlst, maxp1])
dblquad(f, a, b, gfun, hfun, result, abserr[, epsabs, epsrel, limit, status])
tplquad(f, a, b, gfun, hfun, qfun, rfun, result, abserr[, epsabs, epsrel, limit, status])
nquad(f, ndim, ranges, result, abserr[, epsabs, epsrel, limit, status])
quad_vec(f, a, b, result, abserr[, epsabs, epsrel, limit, points, norm, neval, status])
qmc_quad(f, a, b, result, standard_error, state[, n_estimates, n_points, neval, status])
tanhsinh(f, a, b, result, abserr[, maxlevel, minlevel, atol, rtol, nfev, status])
cubature(f, a, b, result, abserr[, rule, rtol, atol, max_subdivisions, &
    points, neval, subdivisions, status])
res = solve_ivp(fun, t_span, y0[, method, t_eval, dense_output, events, &
    event_terminal, event_direction, rtol, atol, first_step, max_step, max_steps, &
    event_set, jac, event_max_occurrences, min_step, lband, uband, &
    jac_matrix, jac_sparsity, fun_vectorized, jac_sparse_matrix, jac_sparse]) ! real y0
cres = solve_ivp(cfun, t_span, cy0[, method, t_eval, dense_output, events, &
    event_terminal, event_direction, rtol, atol, first_step, max_step, max_steps, &
    event_set, event_max_occurrences, jac, jac_matrix, jac_sparsity, &
    fun_vectorized, jac_sparse_matrix, jac_sparse])       ! complex cy0
ores = odeint(func, y0, t[, rtol, atol, dfun, col_deriv, ml, mu, h0, &
    hmax, hmin, mxstep, tcrit, mxordn, mxords, full_output, ixpr, mxhnil, printmessg])
ores = odeint_tfirst(func, y0, t[, rtol, atol, dfun, col_deriv, ml, mu, h0, &
    hmax, hmin, mxstep, tcrit, mxordn, mxords, full_output, ixpr, mxhnil, printmessg])
```

Unlike SciPy's NumPy interfaces, the sampled-data routines currently operate on
one-dimensional arrays rather than accepting an arbitrary `axis`. `x`, when
present, must have the same length as `y`. `cumulative_trapezoid` follows SciPy
1.17.0 in accepting only zero as an explicit `initial`;
`cumulative_simpson` permits an arbitrary scalar initial value. Its explicit `x`
coordinates must be strictly increasing. `romb` requires `2**k + 1` samples.

`fixed_quad` uses SciFort's existing Gauss-Legendre roots. `newton_cotes` returns
SciPy-compatible weights and the leading error coefficient; orders above the
tabulated rules use a centered `[-1,1]` moment system with residual refinement to
avoid the conditioning loss of a raw-power Vandermonde solve.

`quad` is a Fortran-native facade over the double-precision routines in Jacob
Williams' modern QUADPACK 2.1.0 FPM package. Defaults are `epsabs=1.49e-8`,
`epsrel=1.49e-8`, `limit=50`, `limlst=50`, and `maxp1=50`. Finite ordinary
integrals use `dqags`, semi-infinite and whole-line ordinary integrals use
`dqagi`, and finite integrals with `points=` use `dqagp`; break points may be
unordered or repeated and endpoints are ignored after validation.

Weighted `quad` accepts the following SciPy-style names:

- `weight='cos'` or `'sin'`, with `wvar=omega`: `dqawo` on finite intervals or
  `dqawf` on one-sided infinite intervals. Whole-line oscillatory weighting is
  not exposed by QUADPACK and returns `quad_status_invalid_argument`.
- `weight='alg'`, `'alg-loga'`, `'alg-logb'`, or `'alg-log'`, with
  `wvar=alpha` and `wvar2=beta`: finite-interval `dqaws`, requiring
  `alpha > -1` and `beta > -1`.
- `weight='cauchy'`, with `wvar=c`: finite-interval principal-value `dqawc`.

`points=` and `weight=` are mutually exclusive in the Fortran API; providing
both returns `quad_status_invalid_argument` rather than issuing a Python-style
warning and ignoring one option. Precomputed Chebyshev moments (`wopts`), `full_output`, and complex-valued
callbacks are not yet exposed by `quad`. Vector-valued integration is available
through the separate `quad_vec` interface. `status` retains the QUADPACK 0--6
termination codes and adds
`quad_status_cycle_failure=7` for a failed Fourier cycle.


`quad_vec` accepts a callback `f(x, values)` and integrates each element of the
caller-sized `result(:)` independently through `quad`. This correctness-first
implementation supports finite or infinite bounds and optional break points; it
does not yet share callback evaluations across components as SciPy's dedicated
adaptive vector algorithm does. `norm='2'` (default) or `norm='max'` controls
how component error estimates are combined. `neval` is the total number of
scalar callback evaluations across components.

`qmc_quad` integrates a scalar callback `f(x)` over a finite hyperrectangle.
Bounds `a(:)` and `b(:)` may be reversed independently, and zero-width axes give
an exact zero integral. The routine uses `n_estimates=8` independently scrambled
Halton estimates of `n_points=1024` points each by default. Unlike SciPy's Python
API, SciFort requires an explicit mutable `rng_state`; no hidden or process-global
random state is created. `standard_error` is the sample standard error of the
replicate estimates and is NaN when only one estimate is requested. The current
Fortran-native interface evaluates one multidimensional point per callback and
does not yet expose SciPy's custom-QMC-engine or logarithmic-integrand options.

`tanhsinh` implements the Takahasi-Mori double-exponential rule independently
for real scalar callbacks. Finite endpoint singularities are never sampled
exactly; finite, one-sided infinite, whole-line, and reversed intervals are
supported. Defaults are `maxlevel=10`, `minlevel=2`, `atol=0`, and
`rtol=epsilon(1.0_dp)**0.75`. Convergence uses both successive-level change and
the transformed-integrand magnitude at the truncation boundary. A slowly
convergent integral that exhausts the refinement budget returns
`tanhsinh_status_maxlevel=4` with the last estimate and error indicator. SciPy's
array broadcasting, `log=True`, callback, and `preserve_shape` options are not
yet exposed.

`cubature` adaptively integrates over an `ndim`-dimensional hyperrectangle.
Scalar callbacks have the form `f(x) result(value)`; vector callbacks use
`f(x, values)` with caller-sized `result(:)` and `abserr(:)`. Bounds may be
finite, one-sided infinite, whole-line, or independently reversed. The default
`rule='gk21'` is the tensor-product 21-point Gauss-Kronrod rule; `gk15` and
`gauss-kronrod` (an alias of `gk21`) are also available. `rule='genz-malik'`
uses the embedded degree-7/degree-5 Genz-Malik rule and requires at least two
dimensions. Product Gauss-Kronrod cost grows exponentially with dimension, so
Genz-Malik is normally preferable beyond low dimension.

The adaptive driver maintains a global componentwise error estimate, refines the
active region with largest maximum component error, and bisects that region
along every coordinate. `points(dim,npoint)` supplies physical-space feature or
singularity locations; each point is used to partition any region containing it
before numerical evaluation. Infinite intervals are mapped to finite transformed
intervals, and a doubly-infinite coordinate is split at the transformation
branch point. Defaults are `rtol=1e-8`, `atol=0`, and
`max_subdivisions=10000`. Exhausting the refinement budget returns
`cubature_status_max_subdivisions=4` with the current estimate and error.
SciPy's `workers` and Python `args` options are not needed in the current
Fortran callback model; host association or derived-type state can carry
parameters. Callback points are evaluated one at a time rather than in Python-
style batches.

`solve_ivp` supports real first-order systems with `method='RK45'` (default),
`method='RK23'`, `method='DOP853'`, `method='Radau'`, `method='BDF'`, or
`method='LSODA'`. LSODA is currently real-system only. It is
also overloaded for complex-valued systems with `RK23`, `RK45`, `DOP853`, and
`BDF`, matching the SciPy methods in this interface that support complex-domain
integration. The generic is selected by the type of `y0`. Real runs return
`solve_ivp_result`; complex runs return `solve_ivp_complex_result`. Both expose
`t`, `y`, `status`, `success`, `message`, `nfev`, `njev`,
`nlu`, accepted/rejected step counts, optional `t_events`/`y_events`,
`event_indices`, per-event `event_counts`, `n_event_functions`, and `sol` when
`dense_output=.true.`. Complex results store complex `y`, `y_events`, and dense
states while event times and event-function values remain real.
`event_times_for(i)` and `event_states_for(i)` return the
roots and states for one 1-based event index without requiring callers to unpack
the flattened storage. Defaults follow
SciPy for the explicit methods: `rtol=1e-3`, `atol=1e-6`, adaptive first-step
selection, and an unlimited `max_step`. Both `rtol` and `atol` may be scalars,
which are broadcast to every state component, or rank-one arrays with `size(y0)`
entries. All supported real and complex methods use the componentwise error scale
`atol + rtol*abs(y)`. BDF and Radau reduce a vector `rtol` to `maxval(rtol)` only
for their scalar Newton stopping heuristic; adaptive error weights remain
componentwise. Implicit-method finite-difference Jacobians use the same
componentwise absolute tolerances when choosing perturbation sizes. `t_eval` must
be monotone in the integration direction and lie inside `t_span`; backward
integration is supported.

The RK23 Bogacki-Shampine 3(2) and RK45 Dormand-Prince 5(4) tableaux are
adapted from Jacob Williams' BSD-3-Clause `rklib`. RK45's quartic dense-output
coefficient matrix and stage-combination formula are adapted from SciPy 1.17.0
under BSD-3-Clause. DOP853's eighth-order tableau, embedded 5th/3rd-order error
estimator, and native seventh-order dense interpolant are adapted from Jacob
Williams' permissively licensed modern Fortran edition of Hairer/Wanner DOP853.
SciFort implements its own explicit-method adaptive driver, result storage, event
orchestration, and root bookkeeping. `Radau` adapts SciPy 1.17.0's fifth-order
three-stage Radau IIA collocation algorithm, embedded third-order error control,
Newton stopping logic, and cubic collocation dense output under BSD-3-Clause.
`BDF` adapts SciPy 1.17.0's variable-order 1--5 quasi-constant-step BDF/NDF
algorithm, Newton stopping logic, dense-output polynomial, and dense
finite-difference Jacobian strategy under the same license for both real and
complex systems. The LSODA compatibility path combines SciFort orchestration and BDF/NDF stepping with
an Adams Nordsieck implementation adapted from public-domain ODEPACK `DSTODA` and
`DCFODE`. Nonstiff regions start at Adams method order 1 and carry the accepted
scaled-derivative Nordsieck state persistently through method order 12. Trial-step
changes power-rescale the derivative columns, prediction applies the Pascal-triangle
shift, functional correction uses ODEPACK's `EL` and `TESCO` constants, and accepted
steps update the persistent state with the accumulated correction. Successful-step
order/step reconsideration follows the `IALTH`/`TESCO` lower/current/upper-order
comparison with the 10 percent change threshold and the order-dependent `SM1`
stability limits. SciFort still retains nonuniform derivative/time history for state
reconstruction, diagnostics, and conservative recovery rather than as the active
high-order Adams representation. The correction gives a directional local Lipschitz
estimate without extra RHS evaluations; degenerate paths use an accepted-step secant
estimate. Method-switch tests are delayed for 20 successful steps after initialization
or a switch, and Adams-to-BDF switching is skipped above Adams method order 5. Rejected
Adams error-test trials do not perform Jacobian probes or method switches; repeated
failures conservatively restart low-order history with a factor-of-ten step reduction.
For Radau and real/complex BDF, `jac=` associates a
callback `jac(t,y,J)` that fills an `n` by `n` matrix of the matching scalar type.
LSODA uses the same dense callback convention when `lband` and `uband` are omitted.
If either bandwidth is supplied, the omitted bandwidth defaults to zero and an
LSODA analytic callback instead receives packed storage with shape
`(lband+uband+1,n)`, where `Jpacked(uband+i-j+1,j)=J(i,j)`. With no LSODA `jac=`,
finite differences honor the declared band and perturb nonoverlapping columns in
groups, reducing callback work for narrow Jacobians. `lband` and `uband` must each
lie between zero and `size(y0)-1`.
Real Radau/BDF and complex BDF also accept `jac_matrix=` with shape
`(size(y0),size(y0))` as SciFort's Fortran spelling of SciPy's constant
array-like Jacobian form. The matrix is copied into solver state once and does
not increment `njev`. A constant sparse matrix may instead be supplied through
`jac_sparse_matrix=`. Construct real `ivp_sparse_matrix` or complex
`ivp_complex_sparse_matrix` values with
`ivp_csc_matrix(nrows,ncols,indptr,indices,data)`. The CSC representation is
one-based: `indptr` has length `ncols+1`, begins at one, ends at `nnz+1`, and
row indices must be strictly increasing within each column. Values must be
finite. The sparse matrix remains in CSC form in solver state and, like a dense
constant matrix, reports zero Jacobian evaluations. Shifted Newton matrices are
factored directly from CSC storage by SciFort's dependency-free sparse LU path,
which uses partial row pivoting and stores fill-in sparsely. Constant-CSC runs do
not allocate the ordinary dense LU workspace.

State-dependent sparse Jacobians may be supplied through `jac_sparse=`. The
callback has the form `jac_sparse(t,y,J)` with `J` an `intent(out)` real
`ivp_sparse_matrix` or complex `ivp_complex_sparse_matrix`; assigning
`J=ivp_csc_matrix(...)` is the normal construction pattern. Real BDF/Radau and
complex BDF validate every returned CSC matrix before use, so callback outputs
must remain finite, structurally valid, and exactly `size(y0)` by `size(y0)`.
Each callback invocation increments `njev`; unlike `jac_sparse_matrix=`, this
form is treated as state dependent and may be refreshed after Newton difficulty.

When none of `jac=`, `jac_matrix=`, `jac_sparse_matrix=`, or `jac_sparse=` is
supplied, `jac_sparsity=` may provide a logical matrix of the same shape
describing the finite-difference nonzero pattern. Columns whose declared
nonzero rows do not overlap are greedily perturbed together, reducing RHS probes.
The numerical Jacobian values are written directly into a CSC matrix having the
declared structure and its shifted Newton systems use the native sparse LU path;
the ordinary dense Jacobian and LU workspaces are not allocated. A supplied
sparsity pattern is ignored when an analytic callback or constant matrix is
present. LSODA does not accept `jac_matrix=`, `jac_sparse_matrix=`, `jac_sparse=`,
or `jac_sparsity=`; its dense/packed `jac=` and `lband`/`uband` controls remain
the supported Jacobian interface. Complex BDF uses real-axis perturbations and
the sign of the real part of the right-hand side when selecting finite-difference steps,
following SciPy's complex numerical Jacobian convention. Real Radau/BDF and
complex BDF additionally accept an optional `fun_vectorized(t,y,dydt)` callback
for numerical-Jacobian probes, where
`y` and `dydt` have shape `(size(y0),k)` and each column is one independent RHS
state. The ordinary scalar `fun` remains required for stepping. Dense finite
differences use one `k=size(y0)` perturbation batch; with `jac_sparsity=`, `k` is
the number of active column groups. Hidden calls through `fun_vectorized` follow
SciPy's Jacobian accounting convention and are not included in `nfev`. This
separate rank-two callback is SciFort's Fortran-native equivalent of SciPy's
`vectorized=True` contract; explicit Runge-Kutta methods and LSODA reject it
rather than silently ignoring it. Analytic `jac=`, sparse callback `jac_sparse=`,
constant `jac_matrix=`, or constant `jac_sparse_matrix=` paths do not require it.

The legacy `events=` argument accepts one scalar event together with
`event_terminal` and `event_direction`. For multiple real-state events, construct an `ivp_event_set`; for a complex
system use `ivp_complex_event_set`. Set `n_events`, associate the corresponding
vector callback `fun(t,y,values)`,
and optionally provide `terminal(:)` and `direction(:)` controls. SciPy-style
occurrence limits are available through scalar `event_max_occurrences` or vector
`event_set%max_occurrences(:)`: zero means unlimited and a positive integer `N`
terminates after the Nth detected occurrence. An explicitly supplied occurrence
limit overrides the corresponding logical terminal control; `terminal=.true.`
remains equivalent to a limit of one. Negative limits are invalid. Located roots
are stored in flattened `t_events`/`y_events`; `event_indices` identifies the
1-based event function for each root. `event_counts` and the type-bound grouping
helpers provide the SciPy-style per-event view. If several roots occur in one
accepted step, they are ordered in the integration direction and roots after the
first event whose occurrence limit has been reached are discarded, matching
SciPy 1.17.0 event handling. The scalar and vector event interfaces are mutually
exclusive in one call. Event crossings use SciPy's inclusive endpoint/direction
tests, so a zero at the initial point is filtered according to the direction of
the first accepted step.

DOP853 uses its native seventh-order interpolant for `t_eval`, retained dense
solutions, and event-root states. The three extra DOP853 dense stages are only
evaluated for accepted steps that actually need interpolation, so ordinary
endpoint-only runs keep the base stepper evaluation count. RK45 now uses SciPy's
quartic continuous extension for `t_eval`, retained dense solutions, and event
root states without additional right-hand-side evaluations; RK23 retains its
cubic-Hermite interpolant. The complex RK23/RK45/DOP853 overloads use the same
dense formulas with complex state and coefficient arrays. Radau uses its native
cubic collocation polynomial, and real and complex BDF use the native
backward-difference polynomial, for `t_eval`, retained dense solutions, and event
roots. LSODA uses the BDF native polynomial on stiff intervals and its persistent
Adams Nordsieck polynomial on nonstiff intervals. LSODA `min_step` defaults to zero
and is used as the lower bound when an error or convergence failure requires retrying
an Adams or BDF step; it does not force otherwise acceptable startup or endpoint
steps upward. The Adams path starts at method order 1, power-rescales the persistent
Nordsieck array when the trial step changes, Pascal-shifts it to the new endpoint, and
uses the `DSTODA` functional-corrector convergence test. Roundoff-size corrections may
finish after one RHS evaluation; otherwise up to three functional evaluations are used
and a noncontracting corrector forces a retry at one quarter of the failed step. On
success, the accepted state is updated with the ODEPACK `EL` correction and
`IALTH`/`TESCO` controls when step/order changes are reconsidered. The last functional
iterate derivative is retained without an additional endpoint RHS refresh, matching
the persistent-state corrector semantics. The Adams stability estimate retains the
most recent nonzero correction-based Lipschitz estimate, analogous to ODEPACK's
`PDLAST`. LSODA uses ODEPACK's common finite-error ideal-step comparison when choosing
between Adams and BDF; its BDF branch also recovers the accepted endpoint derivative
from the converged BDF equation without an extra RHS call. If the accepted Adams error
estimate is roundoff-polluted, the switch branch selects BDF only when the Adams
proposal is already stability-limited. Rejected Adams error-test trials do not perform
Jacobian probes or method switches; after four consecutive failures SciFort
conservatively restarts low-order Adams history and cuts the retry step by a factor of
ten. Constant and callback CSC Jacobians use the native sparse LU path;
`jac_sparsity=` numerical Jacobians also retain CSC storage through factorization.
Radau and LSODA remain real-only.

`odeint` is a legacy compatibility facade over SciFort's real LSODA path. Its
required callback uses the historical `func(y,t,dydt)` order;
`odeint_tfirst` provides a separate type-safe `func(t,y,dydt)` interface rather
than a runtime `tfirst` flag. The requested `t(:)` grid must be monotone
increasing or decreasing, with repeated values allowed. `odeint_result%y` is
time-major with shape `(size(t),size(y0))`, matching SciPy's output orientation,
and the result also exposes final `nfev`, `njev`, `nlu`, accepted/rejected-step
counts, status, success, and message fields. Fortran host association or explicit
application state replaces Python's `args=` tuple.

The default scalar `rtol` and `atol` are `sqrt(epsilon(1.0_dp))`, the binary64
value corresponding to the historical approximately `1.49012e-8` defaults.
Both tolerances may be scalar or length `size(y0)`, with componentwise LSODA
error weights. `dfun`
may return a dense Jacobian, or with `ml`/`mu` a packed banded Jacobian using the
same storage convention as LSODA. `col_deriv=.true.` accepts the transposed
Jacobian convention. `h0=0` selects automatic startup, `hmax=0` means no
explicit upper bound, and `hmin` is forwarded to LSODA's retry floor. `mxstep`
limits the number of accepted internal LSODA steps allowed before each requested
output is satisfied; zero selects the legacy default of 500. The allowance resets
at each output call while the same Adams/BDF history is retained. `tcrit` accepts a
monotone sequence in the integration direction and imposes the corresponding LSODA
hard step bound when
a critical point is active, preventing internal RHS/Jacobian work from stepping
past that bound while preserving the accumulated Adams/BDF history. `mxordn=0`
and `mxords=0` select the legacy defaults 12 and 5; positive values cap the
Adams and BDF orders respectively, with values above those method maxima clamped.
The underlying `solve_ivp` entry point exposes these only as LSODA-specific
compatibility keywords (`lsoda_tcrit`, `lsoda_max_order_adams`,
`lsoda_max_order_bdf`, and `lsoda_output_max_steps`).

With `full_output=.true.`, `odeint_result` also exposes SciPy-style LSODA
diagnostics for each requested output after the initial time. The allocatable
vectors `hu`, `tcur`, `tolsf`, `tsw`, `nst`, `nfe`, `nje`, `nqu`, and `mused`
therefore have length `max(0,size(t)-1)`. `hu` is the signed size of the last
successful internal step, `tcur` is the internal mesh time that satisfied the
output call, `tsw` is the most recent Adams/BDF switch time, `nst`/`nfe`/`nje`
are cumulative work counters, `nqu` is the order used on the last successful
step, and `mused` is 1 for Adams or 2 for BDF. `imxer`, `lenrw`, and `leniw` are
provided as scalar result fields. On ordinary successful runs `tolsf` is zero
and `imxer=-1`; SciFort does not yet reproduce DLSODA's component-resolved
error-return diagnostics. Repeated initial output times receive deterministic
initial diagnostics (`hu=0`, zero steps, one initial RHS evaluation) rather than
exposing wrapper-dependent uninitialized values. Because `solve_ivp` treats the
final `t_span` endpoint as a hard bound, the final `tcur` may equal the requested
endpoint rather than overshooting it before interpolation. Diagnostic collection
is disabled by default, in which case these vectors are zero length and the
legacy work-length fields remain zero.

`ixpr=.true.` prints an informational line whenever LSODA changes between
Adams and BDF. `mxhnil` controls the maximum number of warnings associated with
a step too small to advance the independent variable reliably; zero selects the
legacy DLSODA default cap of 10. SciFort currently treats that condition as a
terminal LSODA failure, so at most one such warning is observable in one
integration even when a larger cap is requested. `printmessg=.true.` prints the
final integration status; successful runs use the legacy text `Integration
successful.`. These controls only affect reporting and do not alter the
accepted-step sequence or work counters. The underlying `solve_ivp` interface
exposes the first two reporting hooks as the LSODA-specific compatibility
keywords `lsoda_print_switches` and `lsoda_max_hnil_messages`.

New applications should normally use `solve_ivp`; `odeint` is intended for
source/API migration where its legacy conventions are useful.

`dblquad` expects `f(y,x)` with inner bounds `gfun(x)` and `hfun(x)`. `tplquad`
expects `f(z,y,x)`, middle bounds `gfun(x)`/`hfun(x)`, and inner bounds
`qfun(x,y)`/`rfun(x,y)`. The Fortran-native `nquad` accepts `f(x)`, where
`x(1)` is the innermost coordinate and `x(ndim)` the outermost; its
`ranges(axis,x,lower,upper)` callback supplies each axis bounds and may inspect
already assigned outer coordinates. A single tolerance/limit set is currently
applied to every nesting level. Nested `abserr` is the largest QUADPACK error
estimate observed across the levels and `status` is the first nonzero nested
status, or the outer status when that fails.

## Special functions

```fortran
use scifort_special
```

Milestone 6 adds the following elemental scalar interfaces:

```text
erf(x), erfc(x), dawsn(x), erfi(x), erfcx(x)
erf(z), erfc(z), dawsn(z), erfi(z), erfcx(z), wofz(z)  ! complex z
erfinv(y), erfcinv(y)
fresnel(x, s, c), fresnel(z, s, c)  ! complex z
modfresnelp(x, fp, kp), modfresnelm(x, fm, km)
erf_zeros(nt, zeros[, status])
gamma(x), rgamma(x), gammasgn(x)
gamma(z), rgamma(z), loggamma(z)  ! complex z
gammaln(x), loggamma(x)
multigammaln(a, d)
beta(a, b), betaln(a, b)
binom(n, k), poch(a, m)
polygamma(n, x)
factorial(n), factorial2(n), factorialk(n, k)
comb(n, k[, repetition]), perm(n, k)
bernoulli(n, values[, status]), euler(n, values[, status]), stirling2(n, k)
wright_bessel(a, b, x), log_wright_bessel(a, b, x)
spence(x), spence(z), zeta(s[, q]), zetac(s)
lambertw(z[, k, tol]), wrightomega(x), wrightomega(z)  ! complex Lambert W / omega
lpmv(m, v, x), legendre_p(n, z), assoc_legendre_p(n, m, z[, branch_cut, norm])
legendre_p_all(n, z, values[, status, diff_n])
assoc_legendre_p_all(n, m, z, values[, status, branch_cut, norm, diff_n])
lqn(n, z, qn, qnd[, status]), lqmn(m, n, z, qmn, qmnd[, status])
sph_legendre_p(n, m, theta), sph_harm_y(n, m, theta, phi)
sph_legendre_p_all(n, m, theta, values[, status, diff_n])
sph_harm_y_all(n, m, theta, phi, values[, status, diff_n])
legendre_p_derivatives(...), assoc_legendre_p_derivatives(...)
sph_legendre_p_derivatives(...), sph_harm_y_derivatives(...)
hyp0f1(v, x), hyp0f1(v, z)  ! complex z
hyp1f1(a, b, x), hyp1f1(a, b, z)  ! complex z
hyp2f1(a, b, c, x), hyp2f1(a, b, c, z)  ! complex z
hyperu(a, b, x)
pbdv(v, x, d, dprime)
pbvv(v, x, vv, vprime)
pbwa(a, x, w, wprime)
ellip_harm(h2, k2, n, p, s[, signm, signn])
ellip_harm_2(h2, k2, n, p, s)
ellip_normal(h2, k2, n, p)
pbdv_seq(v, x, dv, dpdx[, status])
pbvv_seq(v, x, vv, vp[, status])
pbdn_seq(n, z, dn, dpdz[, status])  ! complex z, n >= 0
digamma(x), psi(x)
digamma(z), psi(z)  ! complex z
ndtr(x), log_ndtr(x)  ! real or complex x
ndtri(p), ndtri_exp(y)  ! real inverse-normal helpers
agm(a, b), owens_t(h, a), voigt_profile(x, sigma, gamma)
huber(delta, r), pseudo_huber(delta, r)
bdtr(k, n, p), bdtrc(k, n, p), bdtri(k, n, y), bdtrik(y, n, p), bdtrin(k, y, p)
nbdtr(k, n, p), nbdtrc(k, n, p), nbdtri(k, n, y), nbdtrik(y, n, p), nbdtrin(k, y, p)
pdtr(k, m), pdtrc(k, m), pdtri(k, y), pdtrik(p, m)
chdtr(v, x), chdtrc(v, x), chdtri(v, p), chdtriv(p, x)
fdtr(dfn, dfd, x), fdtrc(dfn, dfd, x), fdtri(dfn, dfd, p), fdtridfd(dfn, p, x)
gdtr(a, b, x), gdtrc(a, b, x), gdtria(p, b, x), gdtrib(a, p, x), gdtrix(a, b, p)
stdtr(df, t), stdtrit(df, p), stdtridf(p, t)
btdtria(p, b, x), btdtrib(a, p, x)
nrdtrimn(p, std, x), nrdtrisd(mn, p, x)
chndtr(x, df, nc), chndtrix(p, df, nc), chndtridf(x, p, nc), chndtrinc(x, df, p)
ncfdtr(dfn, dfd, nc, f), ncfdtri(dfn, dfd, nc, p)
ncfdtridfd(dfn, p, nc, f), ncfdtridfn(p, dfd, nc, f), ncfdtrinc(dfn, dfd, p, f)
nctdtr(df, nc, t), nctdtrit(df, nc, p), nctdtridf(p, nc, t), nctdtrinc(df, p, t)
smirnov(n, d), smirnovi(n, p), kolmogorov(y), kolmogi(p), tklmbda(x, lmbda)
legendre(n[, monic]), jacobi(n, alpha, beta[, monic]), gegenbauer(n, alpha[, monic])
laguerre(n[, monic]), genlaguerre(n, alpha[, monic])
hermite(n[, monic]), hermitenorm(n[, monic])
chebyt(n[, monic]), chebyu(n[, monic]), chebyc(n[, monic]), chebys(n[, monic])
sh_legendre(n[, monic]), sh_jacobi(n, p, q[, monic])
sh_chebyt(n[, monic]), sh_chebyu(n[, monic])
eval_legendre(n, x), eval_jacobi(n, alpha, beta, x), eval_gegenbauer(n, alpha, x)
eval_laguerre(n, x), eval_genlaguerre(n, alpha, x), assoc_laguerre(x, n[, k])
eval_hermite(n, x), eval_hermitenorm(n, x)
eval_chebyt(n, x), eval_chebyu(n, x), eval_chebyc(n, x), eval_chebys(n, x)
eval_sh_legendre(n, x), eval_sh_jacobi(n, p, q, x)
eval_sh_chebyt(n, x), eval_sh_chebyu(n, x)
roots_legendre(n, nodes, weights[, mu, status])
roots_jacobi(n, alpha, beta, nodes, weights[, mu, status])
roots_gegenbauer(n, alpha, nodes, weights[, mu, status])
roots_laguerre(n, nodes, weights[, mu, status])
roots_genlaguerre(n, alpha, nodes, weights[, mu, status])
roots_hermite(n, nodes, weights[, mu, status])
roots_hermitenorm(n, nodes, weights[, mu, status])
roots_chebyt(n, nodes, weights[, mu, status])
roots_chebyu(n, nodes, weights[, mu, status])
roots_chebyc(n, nodes, weights[, mu, status])
roots_chebys(n, nodes, weights[, mu, status])
roots_sh_legendre(n, nodes, weights[, mu, status])
roots_sh_jacobi(n, p, q, nodes, weights[, mu, status])
roots_sh_chebyt(n, nodes, weights[, mu, status])
roots_sh_chebyu(n, nodes, weights[, mu, status])
p_roots/j_roots/cg_roots/l_roots/la_roots/h_roots/he_roots(...)  ! legacy aliases
t_roots/u_roots/c_roots/s_roots(...)  ! legacy Chebyshev aliases
ps_roots/js_roots/ts_roots/us_roots(...)  ! legacy shifted aliases
expit(x), logit(p), log_expit(x)
xlogy(x, y), xlog1py(x, y)  ! real/complex combinations
entr(x), rel_entr(x, y), kl_div(x, y)
exprel(x), cosm1(x)
cbrt(x), exp10(x), exp2(x), expm1(x), log1p(x), powm1(x, y)  ! expm1/log1p also complex
radian(d, m, s), round(x), sinc(x)
sindg(x), cosdg(x), tandg(x), cotdg(x), softplus(x), diric(x, n)
boxcox(x, lambda), boxcox1p(x, lambda)
inv_boxcox(y, lambda), inv_boxcox1p(y, lambda)
airy(x, ai, aip, bi, bip), airy(z, ai, aip, bi, bip)  ! complex z
airye(x, eai, eaip, ebi, ebip), airye(z, eai, eaip, ebi, ebip)  ! complex z
itairy(x, apt, bpt, ant, bnt)
elliprc(x, y), elliprd(x, y, z), elliprf(x, y, z)  ! real or complex
elliprg(x, y, z), elliprj(x, y, z, p)  ! real or complex
exp1(x), expi(x), expn(n, x)
exp1(z), expi(z)  ! complex z
sici(x, si, ci), shichi(x, shi, chi)
sici(z, si, ci), shichi(z, shi, chi)  ! complex z
i0(x), i1(x)
i0e(x), i1e(x)
j0(x), j1(x), y0(x), y1(x)
itj0y0(x, ij0, iy0), it2j0y0(x, ij0, iy0)
iti0k0(x, ii0, ik0), it2i0k0(x, ii0, ik0)
ber(x), bei(x), ker(x), kei(x), berp(x), beip(x), kerp(x), keip(x)
kelvin(x, be, ke, bep, kep)
ber_zeros/bei_zeros/ker_zeros/kei_zeros(nt, zeros[, status])
berp_zeros/beip_zeros/kerp_zeros/keip_zeros(nt, zeros[, status])
kelvin_zeros(nt, berz, beiz, kerz, keiz, berpz, beipz, kerpz, keipz[, status])
mathieu_a(m, q), mathieu_b(m, q)
mathieu_even_coef(m, q, coeff[, status]), mathieu_odd_coef(m, q, coeff[, status])
mathieu_cem(m, q, x, value, derivative), mathieu_sem(m, q, x, value, derivative)
mathieu_modcem1/modcem2(m, q, x, value, derivative)
mathieu_modsem1/modsem2(m, q, x, value, derivative)
pro_cv(m, n, c), obl_cv(m, n, c)
pro_cv_seq(m, n, c, values[, status]), obl_cv_seq(m, n, c, values[, status])
lmbda(v, x, values, derivatives[, status]), besselpoly(a, lmb, nu)
jv(v, z), jve(v, z), yv(v, z), yve(v, z), jn(v, z)  ! real/complex z
jvp(v, z[, n]), yvp(v, z[, n]), ivp(v, z[, n]), kvp(v, z[, n])  ! real/complex z
hankel1(v, z), hankel2(v, z), hankel1e(v, z), hankel2e(v, z)  ! real/complex z
h1vp(v, z[, n]), h2vp(v, z[, n])  ! real/complex z
k0(x), k0e(x), k1(x), k1e(x)
yn(n, x), kn(n, x)
iv(v, z), ive(v, z), kv(v, z), kve(v, z)  ! real/complex z
spherical_jn(n, z[, derivative]), spherical_yn(n, z[, derivative])  ! real/complex z
spherical_in(n, z[, derivative]), spherical_kn(n, z[, derivative])  ! real/complex z
riccati_jn(n, x, values, derivatives[, status])
riccati_yn(n, x, values, derivatives[, status])
struve(v, x), modstruve(v, x)
itstruve0(x), it2struve0(x), itmodstruve0(x)
log_besseli_scaled(nu, x), besseli_ratio(nu, x)
log_besselk(nu, x)
log_besselk_order_derivative(nu, x), besselk_log_derivative_x(nu, x)
hurwitz_zeta(s, q), hurwitz_zeta_derivative(s, q)
```

Complex `ndtr(z)` evaluates `0.5*erfc(-z/sqrt(2))` through SciFort's Faddeeva
foundation. `log_ndtr(z)` returns the principal logarithm and switches between
scaled-`erfcx` evaluation in the left half-plane, direct evaluation in the central
region, and a reflected `log1p` form when the complementary normal CDF is small.
This avoids real-tail underflow and loss of the tiny logarithm when `ndtr(z)` is
close to one. Complex `expm1(z)` and `log1p(z)` use cancellation-safe Cartesian
forms near the origin while preserving the principal logarithm branch. `xlogy`
and `xlog1py` accept real/complex arguments in either position; an exact zero
multiplier returns zero for non-NaN logarithm arguments, matching SciPy's special
zero convention. Existing real overloads retain their prior kernels and behavior.

`lmbda(v,x,values,derivatives[,status])` returns the Jahnke-Emden Lambda sequence for
orders `v-int(v), ..., v`, together with derivatives in `x`. Finite `0 <= v <= 100000`
and finite real `x` are required; the implementation evaluates the exact hypergeometric identity
`Lambda_v(x) = 0F1(;v+1;-x**2/4)`. Thus all sequence values are exactly one at `x=0`,
and derivative parity follows the mathematical definition for negative `x`. Invalid input
returns empty arrays and `bessel_misc_status_invalid_argument`.

`besselpoly(a,lmb,nu)` evaluates the defining weighted integral
`integral_0^1 x**lmb J_nu(2*a*x) dx` for finite real parameters with `a >= 0` whenever
the endpoint integral converges. It is elemental. Endpoint smoothing and composite
16-point Gauss-Legendre quadrature handle ordinary and oscillatory cases; near an
integrable singular endpoint, the leading Bessel-series term is integrated analytically
and quadrature is applied only to the smoother residual. Divergent endpoints return a
signed IEEE infinity and invalid/nonfinite or negative-scale inputs return NaN. Work grows
roughly linearly with `a` because increasingly oscillatory integrands use more panels.

`ber(x)` and `bei(x)` are the real and imaginary parts of
`J_0(x*exp(3*pi*i/4))`; `ker(x)` and `kei(x)` are the real and imaginary parts of
`K_0(x*exp(pi*i/4))`. `berp`/`beip`/`kerp`/`keip` are derivatives with respect to
`x`, and `kelvin(x,be,ke,bep,kep)` returns the four corresponding complex pairs in
one elemental call. The B-family accepts all finite real `x` with even values and
odd derivatives. The K-family follows SciPy's real domain: negative `x` returns NaN,
`ker(0)=+infinity`, `kei(0)=-pi/4`, `kerp(0)=-infinity`, and `keip(0)=0`.

SciFort evaluates the convergent complex Bessel series at small/moderate arguments
and switches to DLMF large-argument expansions on the fixed Kelvin rays. The
large-B branch factors the exponentially growing magnitude from the oscillatory
complex phase so finite arguments near binary64 overflow do not raise an avoidable
floating-point exception. The large-K branch decays to zero once the binary64 result
is below the representable range.

The eight `*_zeros` procedures return strictly positive roots and omit the
mathematical zero of `bei`, `berp`, `beip`, or `keip` at the origin, matching SciPy's
root-table convention. `kelvin_zeros` returns all eight tables together. The current
independent sign-bracketing implementation supports `1 <= nt <= kelvin_max_zeros`,
where `kelvin_max_zeros=200`; invalid counts return empty arrays and
`kelvin_status_invalid_argument`. `kelvin_status_success` and
`kelvin_status_numerical_failure` complete the status set.

`mathieu_a(m,q)` and `mathieu_b(m,q)` return the integer-order Mathieu
characteristic values for finite real `q`; `m >= 0` is required for `a` and
`m >= 1` for `b`. `mathieu_even_coef` and `mathieu_odd_coef` return normalized
Fourier coefficients for the even (`ce_m`) and odd (`se_m`) angular families.
The coefficient APIs follow SciPy in requiring finite `q >= 0`; invalid input
returns an empty array and `mathieu_status_invalid_argument`. The returned
coefficient array contains the full internal Fourier truncation rather than
SciPy's implementation-specific shorter heuristic length.

`mathieu_cem` and `mathieu_sem` preserve SciPy's historical angle convention:
`x` is supplied in degrees, while `derivative` is the derivative with respect
to the corresponding radian argument. They accept finite positive or negative
`q`; negative parameters use the exact DLMF parity/reflection identities. At
`q=0`, for example, `ce_0=1/sqrt(2)`, `ce_m(x,0)=cos(m*x)`, and
`se_m(x,0)=sin(m*x)` for positive `m`, under the standard normalized Mathieu
convention.

The implementation forms one of the four symmetric tridiagonal Fourier
eigenproblems from DLMF section 28.4 and reuses SciFort's pure symmetric Jacobi
eigensolver. The truncation grows with order and approximately with `sqrt(abs(q))`;
the current binary64 implementation caps the matrix dimension at 256. Requests
that exceed this implementation ceiling return NaN (or an empty coefficient
array) with `mathieu_status_numerical_failure`. `mathieu_status_success` is zero.
`mathieu_modcem1`/`mathieu_modcem2` and `mathieu_modsem1`/`mathieu_modsem2`
return the first- and second-kind even/odd radial Mathieu functions. They require
`q > 0`, `m >= 0` for the even pair, and `m >= 1` for the odd pair. The radial
argument `x` is interpreted in radians. This matches SciPy 1.17 numerical behavior
and the DLMF radial variable, even though the corresponding SciPy 1.17 documentation
pages describe the input as degrees. Derivatives are with respect to `x`.

The implementation uses the DLMF Bessel cross-product expansions with an adaptive
Fourier anchor chosen among well-scaled coefficients. Negative first-kind arguments
use exact even/odd parity. Negative second-kind arguments are continued from the
origin by the modified-Mathieu differential equation symmetry, avoiding the severe
cancellation of direct cross-product evaluation there. Invalid parameters return NaN.

`pro_cv(m,n,c)` and `obl_cv(m,n,c)` return the characteristic values of prolate and
oblate spheroidal wave functions for nonnegative integer `m`, integer `n >= m`, and
finite real `c`. The functions are even in `c`, and reduce exactly to `n*(n+1)` at
`c=0`. `pro_cv_seq` and `obl_cv_seq` return all degrees `m..n`; `n < m` produces an
empty successful sequence, matching SciPy's sequence behavior. Invalid scalar inputs
return NaN. The sequence routines report `spheroidal_status_success`,
`spheroidal_status_invalid_argument`, or `spheroidal_status_numerical_failure`.

SciFort forms the even or odd DLMF 30.16 tridiagonal eigenproblem with signed
`gamma**2 = c**2` for prolate and `-c**2` for oblate cases. The nonsymmetric DLMF
recurrence matrix is converted to its diagonally similar symmetric tridiagonal form,
and the requested ordered eigenvalue is isolated by a Sturm-count bisection. SciPy's
characteristic-value convention equals the DLMF eigenvalue plus the signed
`gamma**2` shift, which is incorporated directly in the tridiagonal diagonal to avoid
large-c cancellation. Matrix dimensions grow adaptively and are capped at 256; cases
that do not converge before that ceiling return NaN or numerical-failure status.

The first-kind angular interfaces are pure elemental subroutines:

```text
call pro_ang1(m, n, c, x, value, derivative)
call pro_ang1_cv(m, n, c, cv, x, value, derivative)
call obl_ang1(m, n, c, x, value, derivative)
call obl_ang1_cv(m, n, c, cv, x, value, derivative)
```

They require `m >= 0`, `n >= m`, finite `c`, and `abs(x) < 1`; `derivative` is
with respect to `x`. The `_cv` forms reuse a precomputed characteristic value.
As for the characteristic-value routines, the mathematical dependence is through
`c**2`, so SciFort treats positive and negative `c` identically. Invalid input or
failure to converge within the current 256-term ceiling returns NaN for both
outputs.

SciFort evaluates the DLMF Ferrers-function expansion rather than a center Taylor
series. Relative tridiagonal-eigenvector coefficients are reconstructed by stable
continued-fraction recurrences from both sides of the target degree, and the series
is accumulated with compensated summation. The overall scale follows the
Flammer/SciPy convention: at `x=0`, the value for even `n-m`, or the derivative for
odd `n-m`, equals the corresponding phase-free associated Legendre normalization.
SciFort enforces that normalization directly, including in high-parameter cases
where legacy SciPy outputs can lose center normalization through cancellation.

The first-kind radial spheroidal interfaces are also pure elemental subroutines:

```text
call pro_rad1(m, n, c, x, value, derivative)
call pro_rad1_cv(m, n, c, cv, x, value, derivative)
call obl_rad1(m, n, c, x, value, derivative)
call obl_rad1_cv(m, n, c, cv, x, value, derivative)
```

They require `m >= 0`, `n >= m`, finite nonnegative `c`, and finite radial
arguments. The prolate forms require `x > 1`; the oblate forms accept `x >= 0`.
The `_cv` forms reuse a precomputed characteristic value. At `c=0`, SciFort
returns `(value,derivative)=(1,0)` for `m=n=0` and `(0,0)` otherwise, matching
SciPy's limiting convention. Invalid inputs or failure of the adaptive finite
expansion return NaN for both outputs.

SciFort evaluates DLMF 30.11.3 with the existing tridiagonal eigenvector
coefficients. The angular coefficients `a_n,k^m` are transformed to the required
`a_n,k^-m` using DLMF 30.8.6 in logarithmically scaled form, and the resulting
spherical-Bessel series and normalization sum are accumulated with compensated
summation. The oblate origin is evaluated from the analytic small-argument limit,
avoiding the cancelling product of a singular coordinate prefactor and a vanishing
spherical Bessel function. The current batch implements the first-kind radial
family.

The second-kind radial interfaces follow the same value/derivative convention:

```text
call pro_rad2(m, n, c, x, value, derivative)
call pro_rad2_cv(m, n, c, cv, x, value, derivative)
call obl_rad2(m, n, c, x, value, derivative)
call obl_rad2_cv(m, n, c, cv, x, value, derivative)
```

They require `m >= 0`, `n >= m`, finite positive `c`, and finite radial
arguments. The prolate forms require `x > 1`; the oblate forms accept `x >= 0`.
The `_cv` variants reuse a precomputed characteristic value. Because the
second-kind normalization is singular as `c -> 0`, `c=0` returns NaN rather
than the very large finite sentinel used by some legacy implementations.

For numerically safe large arguments SciFort evaluates the DLMF 30.11.3
spherical-Bessel-Y series with the same `a_n,k^-m` coefficients and normalization
used by the first-kind family. Near the prolate singular boundary and for small
oblate arguments, where direct summation of that Y series is badly conditioned,
SciFort evaluates the series at an automatically selected stable anchor and
continues both the value and derivative to the requested argument with an
adaptive embedded Dormand-Prince method applied to the spheroidal differential
equation. This construction preserves the first/second-kind radial Wronskian
without importing the legacy SciPy/specfun small-argument implementation.

The integer-order Bessel-zero tables are array-producing pure subroutines:

```text
call jn_zeros(n, nt, zeros[, status])
call jnp_zeros(n, nt, zeros[, status])
call yn_zeros(n, nt, zeros[, status])
call ynp_zeros(n, nt, zeros[, status])
call jnyn_zeros(n, nt, jzeros, jpzeros, yzeros, ypzeros[, status])
call jnjnp_zeros(nt, zeros, orders, serials, kinds[, status])
call y0_zeros(nt, zeros, derivatives[, status])
call y1_zeros(nt, zeros, derivatives[, status])
call y1p_zeros(nt, zeros, values[, status])
```

`bessel_zeros_status_success`, `bessel_zeros_status_invalid_argument`, and
`bessel_zeros_status_numerical_failure` provide the corresponding status codes.
The fixed-order routines require `nt >= 1`. `jnjnp_zeros` follows SciPy's
1200-entry compatibility ceiling and requires `1 <= nt <= 1200`; `orders`
contains the integer Bessel order, `serials` contains the one-based root number,
and `kinds` is zero for a J root and one for a J-prime root. The leading entry is
the special `J_0'(0)` root with `serials(1)=0`. Because `J_0'(x)=-J_1(x)`, those
positive roots are exact mathematical ties; SciFort coalesces each tied pair to
one binary64 value and deterministically places the J_1 entry before the J_0'
entry. SciPy's compiled routine can reverse a tied pair when its independently
computed approximations differ by last-bit roundoff, so tie ordering is not a
numerical compatibility requirement.

`y0_zeros`, `y1_zeros`, and `y1p_zeros` have real and complex overloads,
selected by the types of both allocatable output arrays. Real outputs retain
the positive-real-root behavior. Complex outputs return the upper-half-plane
root sequence corresponding to SciPy's `complex=True`, requiring `1 <= nt <= 1000`.
The associated outputs are, respectively, `Y_0'`, `Y_1'`, and `Y_1` at the roots.
SciFort follows `Y_0'=-Y_1` for both types; SciPy's observed numerical second
array has the opposite sign despite its documented derivative identity.
The first complex `Y_1'` root has positive real part, as in SciPy's numerical
output; the rest have negative real part. Complex numerical failure returns
status 2 and two empty arrays. See [complex Bessel zeros](complex_bessel_zeros.md)
for ordering, convergence limits, and reference discrepancies. The existing real
routines retain their behavior without an artificial 140-root ceiling.

`dawsn`, `erf`, `erfc`, `erfi`, and `erfcx` are overloaded at the top-level
`scifort_special` API for real and complex binary64 arguments. Real arguments retain
the existing real kernels. Complex `erfcx` uses `erfcx(z)=wofz(i*z)`; complex
`erfc`, `erf`, `erfi`, and Dawson reuse the Faddeeva foundation together with
standard error-function identities, local entire series near the origin, and symmetry
rearrangements selected to reduce cancellation. `wofz(z)` returns
`exp(-z**2)*erfc(-i*z)` for complex binary64 `z`. The complex implementation is
independent SciFort code using DLMF series, integral, continued-fraction, and
reflection identities. No complex C ABI symbols are provided.

`erfcx` evaluates the real scaled complementary error function without forming
`exp(x**2)` in the positive large-argument tail. `fresnel` is a pure elemental
subroutine overloaded for real and complex binary64 arguments and returns
SciPy-ordered outputs `(S,C)` through its `s` and `c` arguments. The real
implementation combines convergent series, Gauss-Legendre quadrature, and DLMF
auxiliary-function asymptotics. Large finite real arguments use compensated
reduction of `x**2 modulo 4`; this preserves the mathematically correct binary64
phase in cases where SciPy 1.17.0 loses low quadratic-phase bits. Complex Fresnel
evaluation uses DLMF 7.5.8 error-function interrelations over the existing
complex `erf`/Faddeeva foundation, plus the defining entire series in a
well-conditioned diagonal sector. Exact real and imaginary axes reuse the real
Fresnel kernel through the analytic continuation identities
`C(i*y)=i*C(y)` and `S(i*y)=-i*S(y)`. This also avoids the discontinuous
large-imaginary-axis values returned by SciPy 1.17.0's complex path.

The modified Fresnel interfaces accept a real lower limit `x` and return complex pairs:

```text
call modfresnelp(x, fp, kp)
call modfresnelm(x, fm, km)
```

They implement the documented integrals `F_+(x)=integral_x^inf exp(i*t**2) dt` and
`F_-(x)=integral_x^inf exp(-i*t**2) dt`, together with
`K_+(x)=exp(-i*(x**2+pi/4))*F_+(x)/sqrt(pi)` and the conjugate negative form.
For nonnegative arguments the scaled result is evaluated through the existing `wofz` kernel,
then an error-free split square reconstructs the oscillatory phase without first rounding
`x**2` to a single binary64 value. For negative real arguments SciFort follows the integral
definition and symmetry exactly; this intentionally differs from SciPy 1.17.0's first `F`
component, while its `K` component agrees. Nonfinite inputs, and finite magnitudes whose
square exceeds binary64 range, return complex NaNs.

The Fresnel zero-table interfaces are:

```fortran
call fresnel_zeros(nt, zeros_sine, zeros_cosine[, status])
call fresnels_zeros(nt, zeros[, status])
call fresnelc_zeros(nt, zeros[, status])
```

They return the first `nt` zeros in the first quadrant, ordered by increasing magnitude, matching SciPy's convention. `fresnel_zeros_status_success`, `fresnel_zeros_status_invalid_argument`, and `fresnel_zeros_status_numerical_failure` provide status codes. `nt` must be at least one. The roots are seeded from NIST DLMF 7.13.5-7.13.7 and refined with Newton iterations using the existing complex Fresnel kernel and the exact derivatives `S'(z)=sin(pi*z**2/2)` and `C'(z)=cos(pi*z**2/2)`. SciFort does not impose SciPy 1.17.0's observed high-index NaN limit.

The error-function zero-table interface is:

```fortran
call erf_zeros(nt, zeros[, status])
```

It returns the first `nt` nonzero zeros of `erf(z)` in the first quadrant, ordered by
increasing magnitude. `erf_zeros_status_success`, `erf_zeros_status_invalid_argument`,
and `erf_zeros_status_numerical_failure` provide status codes, and `nt` must be at least
one. Initial values use the NIST DLMF 7.13.1-7.13.2 large-index expansion and are refined
by Newton iterations with SciFort's existing complex `erf` kernel and the exact derivative
`2*exp(-z**2)/sqrt(pi)`. The implementation remains finite in regression tests through at
least the 1000th zero.

`erf` and `erfc` wrap the Fortran intrinsic implementations. `erfinv` and
`erfcinv` are formed through the standard-normal inverse identity and therefore
inherit the accuracy and extreme-tail range of the current `ndtri`/normal-PPF
kernel. The real `gamma` and `rgamma` specifics use the existing real `gammaln` magnitude plus
real-axis sign intervals, with explicit overflow guards and signed zero handling.
Negative integer poles of real `gamma` return NaN, matching SciPy 1.17.0, while
real `rgamma` returns zero there. The complex overloads use principal-branch
`loggamma(z)`, evaluated by recurrence/reflection plus a DLMF Stirling series;
complex `rgamma` is entire and returns zero at nonpositive integer poles, while
complex `gamma` and `loggamma` return complex NaNs at those poles. On the
negative real axis, complex `loggamma(x +/- 0i)` takes the corresponding upper/lower
limit of the single documented branch cut. `gammasgn` remains real-only and
preserves the sign of zero. `gammaln` remains real-only and returns
`log(abs(Gamma(x)))`, with nonpositive integer poles mapped to `+infinity`. The
real `loggamma` specific is unchanged: positive arguments agree with `gammaln`,
zero maps to `+infinity`, and negative real arguments return NaN. No complex C
ABI is exposed.
`multigammaln(a,d)` implements the standard multivariate log-gamma sum for
integer `d >= 1` and `a > (d-1)/2`; invalid Fortran arguments return NaN rather
than raising SciPy's Python `ValueError`. `beta` and `betaln` implement the real
meromorphic continuation, including finite cancellations when a negative-integer
Gamma pole is canceled by a positive-integer argument and signed zero when only
the denominator Gamma factor has a pole. Simultaneous numerator poles are
singular and are represented as infinity. `poch(a,m)` is the generalized rising
factorial and handles integer pole cancellations explicitly; `binom(n,k)` is the
generalized real binomial coefficient and follows SciPy's NaN convention when
`n` is a negative integer. `polygamma(n,x)` accepts a nonnegative integer
derivative order: `n=0` delegates to the real `digamma` specific, positive
arguments use the Hurwitz-zeta identity, and negative nonintegers use
recurrence/reflection with exact integer/half-integer phase reduction. Top-level
`digamma` and `psi` are overloaded for real and complex binary64 arguments. The
real specific uses recurrence/asymptotics for positive arguments and the
reflection formula for negative nonintegers; negative integer poles return NaN
and signed zero follows SciPy (`digamma(+0)=-infinity`,
`digamma(-0)=+infinity`). The complex specific is meromorphic with no branch cut,
uses DLMF 5.5.2 recurrence and 5.5.4 reflection plus the 5.11.2 asymptotic
series, and evaluates `cot(pi*z)` in scaled form near negative-real poles to
avoid intermediate underflow/overflow. Complex nonpositive-integer poles return
complex NaNs. Exact real-axis complex nonpoles reuse the real kernel. `psi` is
an alias of `digamma` for both real and complex arguments. No complex C ABI is
exposed.

`factorial(n)` provides SciPy's floating `exact=False`, `extend='zero'` real
semantics: negative real arguments return zero and nonnegative real arguments
represent `Gamma(n+1)`. The real implementation evaluates that identity by a
guarded recurrence from `Gamma(1+fractional_part(n))`, avoiding accuracy loss
from exponentiating a near-overflow log-gamma value. Integer actual arguments use
a guarded finite product so small integer factorials are exact when representable
in binary64. `factorial2(n)` and
`factorialk(n,k)` take integer Fortran arguments, matching the default real-zero
extension's integer domain; negative `n` returns zero, `0!(k)=1`, and invalid
`k <= 0` returns NaN instead of raising Python's `ValueError`. `comb(n,k)` and
`perm(n,k)` expose SciPy's floating combinatorial paths and accept either default
integers or binary64 reals. Their ordinary domain returns zero when `n < 0`,
`k < 0`, or `k > n`; `comb(...,repetition=.true.)` follows SciPy's repeated-choice
transformation through generalized `binom`. Python's arbitrary-precision
`exact=True` option and the complex factorial extensions are not implemented.

`bernoulli(n,values[,status])` and `euler(n,values[,status])` are pure array-producing
subroutines with zero-based output bounds `0:n`, so the mathematical index is also the Fortran
index. They require `n >= 0`; invalid input returns an empty array and
`integer_sequences_status_invalid_argument`. `bernoulli` uses the zeta/factorial identity for
even indices with the `B_1=-1/2` convention and exact zeros at odd indices above one. `euler`
uses the secant-number binomial recurrence while the result is representable, then a guarded
logarithmic asymptotic identity once binary64 overflow is imminent. `stirling2(n,k)` is a pure
elemental integer-input function returning binary64 values. Negative indices and `k>n` return
zero, and the implementation uses the direct positive recurrence in binary64 rather than SciPy's
large-`n` Temme approximation; therefore some `n>50` results can differ from SciPy's approximate
values. Arbitrary-precision `exact=True` Stirling values are not provided.

`wright_bessel(a,b,x)` evaluates Wright's generalized Bessel function
`sum(x**k/(k!*Gamma(a*k+b)), k=0..infinity)` on the same nonnegative real parameter
domain documented by SciPy (`a >= 0`, `b >= 0`, `x >= 0`). `log_wright_bessel` returns
the natural logarithm without first forming the potentially overflowing function value. The
implementation locates the largest positive series term from the log-derivative, sums outward
in scaled log form with compensated accumulation, and switches to a Laplace/saddle expansion
only when the modal index is exceptionally large. Invalid negative or nonfinite parameters
return NaN; a positive-infinite `x` returns positive infinity.

`spence` is overloaded for real and complex binary64 arguments and follows SciPy's
`Li_2(1-z)` convention. The existing real specific is unchanged: it is defined for
finite `x >= 0`, returns `pi**2/6` at zero, and returns NaN on the negative-axis
branch cut and at infinite arguments. The complex specific uses the principal DLMF
dilogarithm branch, combining the power series, inversion formula, and a logarithmic
series around `abs(1-z)=1`. Complex negative-real arguments return the sided cut
limit; signed zero in the imaginary component selects the upper/lower side, matching
SciPy's complex convention. Nonfinite complex arguments return complex NaNs. No
complex C ABI is exposed. `zeta(s)` provides real and complex
Riemann zeta continuation, including exact real-axis trivial zeros and the pole
at `s=1`; `zetac(s)` evaluates `zeta(s)-1` directly for `s>1` so tiny positive
values are not lost to subtraction. The optional `q` in `zeta(s,q)` exposes
the real Hurwitz form for `s >= 1`: positive `q` uses the existing
Euler-Maclaurin kernel, while negative noninteger `q` is supported when `s`
is an integer by separating the finite negative-base terms from positive-q
Hurwitz tails. Noninteger `s` with negative `q` returns NaN, and nonpositive
integer `q` is a pole. The pure elemental complex overload supports the
one-argument Riemann form only. For nonreal arguments it requires
`abs(real(s)) <= 10000` and `abs(aimag(s)) <= 100000`, and returns complex NaNs
outside these limits or when guarded near-pole/reflection arithmetic would
overflow. Exactly real complex inputs retain the existing real behavior.
See [complex zeta](complex_zeta.md) for accuracy, range limits, and provenance.
For large positive real `s`, SciFort
`zetac` can retain representable values such as `zetac(1000)` where SciPy
1.17.0 returns zero.

`ndtri_exp(y)` inverts `log_ndtr(x)=y` for `y <= 0`. Ordinary representable
probabilities reuse `ndtri`; upper-tail values form `1-exp(y)` with `expm1`,
and very negative log probabilities use Newton correction of the normal-tail
asymptotic without exponentiating `y`. Thus values such as `y=-1e6` remain
finite even though `exp(y)` underflows. `y=0` maps to positive infinity,
negative infinity maps to negative infinity, and positive or NaN inputs return
NaN.

`owens_t(h,a)` evaluates Owen's T function. The implementation is even in `h`
and odd in `a`; it integrates the tangent-substituted defining integral with a
32-point Gauss-Legendre rule for `abs(a) <= 1` and uses the reciprocal-slope
identity for larger magnitudes. Infinite `a` uses the exact half-normal-tail
limit. `voigt_profile(x,sigma,gamma)` evaluates the Voigt profile through the
existing Faddeeva `wofz` kernel when `sigma /= 0`; the exact Cauchy, Gaussian,
and doubly-degenerate limits are handled explicitly. As in SciPy's ufunc,
negative `sigma` or `gamma` are evaluated algebraically rather than rejected as
distribution parameters.

`agm(a,b)` is the real arithmetic-geometric mean. Equal-sign negative inputs
return the negative of the positive AGM, mixed strict signs return NaN, and a
finite zero input returns zero. `huber(delta,r)` and `pseudo_huber(delta,r)`
implement SciPy's real loss functions, including positive infinity for negative
`delta`; their finite paths are rearranged to avoid avoidable overflow and
cancellation.

The legacy distribution-special helpers are implemented directly from SciFort's regularized incomplete beta/gamma kernels and standard-normal transforms rather than by calling the higher-level distribution modules. `bdtr` and `bdtrc` use `floor(k)` and require integer `n`; `bdtri` inverts the Bernoulli probability. `nbdtri` likewise inverts the negative-binomial success probability, while `pdtri` inverts the Poisson mean. `chdtri` follows SciPy's historical convention of inverting the *survival* probability. `gdtr` uses the SciPy rate/shape ordering `(a,b)` and `gdtria`/`gdtrix` invert the rate and variate respectively. `stdtrit` is the ordinary lower-tail Student-t quantile.

The companion parameter-inverse helpers solve the same continuous incomplete-beta/gamma identities without rounding their inverted count or shape parameter. `bdtrik`, `bdtrin`, `nbdtrik`, `nbdtrin`, and `pdtrik` therefore return real-valued continuous inverses even though the corresponding probability laws are discrete. `btdtria`/`btdtrib`, `chdtriv`, `gdtrib`, and `stdtridf` invert shape or degrees-of-freedom parameters with safeguarded geometric bracketing. `nrdtrimn` and `nrdtrisd` use the exact normal-CDF transform. `fdtridfd` is not globally one-to-one for every `(dfn,x)` pair; when two positive denominator-df solutions exist, SciFort returns the smaller root and documents that branch explicitly instead of using a large failure sentinel. Invalid or nonunique finite-parameter cases return NaN under SciFort's numerical API policy. Two SciPy 1.17 legacy endpoint quirks are intentionally not copied: `pdtri(k,1)` returns the mathematically exact zero, and `stdtrit(df,0)` returns negative infinity rather than positive infinity.

The noncentral legacy helpers reuse the standard SciFort `ncx2`, `ncf`, and `nct` distribution kernels. `chndtrix`, `ncfdtri`, and `nctdtrit` are their lower-tail quantiles; the remaining `*idf`, `*idfn`, and `*inc` routines invert a degrees-of-freedom or noncentrality parameter with safeguarded searches. The noncentral-F CDF is not globally one-to-one in either degrees of freedom, so `ncfdtridfd` and `ncfdtridfn` deliberately return the smallest positive root when more than one exists. `nctdtridf` uses the same smallest-positive-root convention if its equation is nonunique. `nctdtrit(df,nc,0)` returns negative infinity, following the mathematical lower-tail inverse rather than SciPy 1.17.0's observed positive-infinity legacy result. The wrappers preserve SciPy-special validation of finite positive degrees of freedom even though the higher-level SciFort noncentral-t distribution also supports an infinite-df shifted-normal limit.

`smirnov(n,d)` and `smirnovi(n,p)` are the historical one-sided
Kolmogorov-Smirnov survival function and inverse survival function; they
delegate to `ksone_sf` and `ksone_isf` with unit location/scale and require
a positive integer sample size. `kolmogorov(y)` and `kolmogi(p)` likewise
delegate to the asymptotic two-sided `kstwobign` survival/inverse-survival
kernels. `tklmbda(x,lmbda)` is the standardized Tukey-lambda CDF. These
special-function names are compatibility wrappers only; SciFort's higher-level
distribution APIs remain the preferred interfaces when PDF/log-CDF or fitted
location/scale behavior is needed.

The orthogonal-polynomial coefficient constructors implement the same classical
nonnegative-integer-degree families as SciPy's `orthopoly1d` constructors, but use an
idiomatic Fortran return value: an allocatable real coefficient array in descending
power order, matching SciPy's `.coeffs` field. The optional logical `monic` defaults
to false and, when true, scales the leading coefficient to one. Invalid degrees or
family parameters return an empty array. The `gegenbauer` constructor handles the
`alpha=0` monic limit analytically, where the unnormalized polynomial is identically
zero for positive degree. Evaluation and quadrature metadata are deliberately kept in
the existing `eval_*` and `roots_*` APIs instead of emulating Python's polynomial
object. The constructors are independent DLMF recurrence implementations.

The orthogonal-polynomial evaluators currently implement the nonnegative-integer-degree
polynomial branch. Legendre, Jacobi, Gegenbauer, Laguerre/generalized Laguerre,
Chebyshev, and shifted-family interfaces accept real or complex binary64 `x`;
`eval_hermite` and `eval_hermitenorm` follow SciPy's real-argument interface.
Jacobi parameters are currently restricted to the classical orthogonality domain
`alpha > -1`, `beta > -1`; Gegenbauer requires `alpha > -1/2`; generalized
Laguerre requires `alpha > -1`; shifted Jacobi requires `q > 0` and `p-q > -1`.
Negative degree, nonfinite arguments, or parameters outside those implemented
domains return NaN. The recurrences are independently implemented from NIST DLMF
18.9; `eval_chebyc` uses `C_n(x)=2*T_n(x/2)`, `eval_chebys` uses
`S_n(x)=U_n(x/2)`, and the shifted families use the standard affine mappings.
SciPy's noninteger-degree hypergeometric continuation remains future work.
`assoc_laguerre(x,n,k)` is the SciPy compatibility wrapper around
`eval_genlaguerre(n,k,x)` and defaults `k` to zero; both real and complex
arguments are supported.

The matching `roots_*` procedures return allocatable ascending node and weight
arrays for Gaussian quadrature of order `n >= 1`; optional `mu` is the exact
zeroth moment (sum of the weights), and optional `status` uses
`orthogonal_roots_status_success`, `orthogonal_roots_status_invalid_argument`, or
`orthogonal_roots_status_numerical_failure`. Legendre/Jacobi/Gegenbauer,
Laguerre/generalized Laguerre, and Hermite/Hermite-normalized rules use the
Golub-Welsch Jacobi-matrix identity with SciFort's existing dense symmetric
Jacobi eigensolver. Chebyshev first/second-kind rules and their `C`/`S` and
shifted variants use closed-form nodes and weights. The same classical parameter
domains as the evaluator interfaces are enforced. Invalid input returns empty
node/weight arrays, NaN `mu`, and invalid-argument status. This initial
implementation favors dependency-free correctness for moderate quadrature order;
its dense eigensolver has O(n^3) work and O(n^2) storage, so SciPy's specialized
large-order Hermite asymptotics are not yet reproduced. The historical
SciPy names `p_roots`, `j_roots`, `cg_roots`, `l_roots`, `la_roots`, `h_roots`,
`he_roots`, `t_roots`, `u_roots`, `c_roots`, `s_roots`, `ps_roots`, `js_roots`,
`ts_roots`, and `us_roots` are exact generic aliases of their corresponding
canonical `roots_*` procedures, so nodes, weights, optional `mu`, and status
behavior are identical.

The logistic and entropy helpers follow SciPy's real scalar conventions at
zero and at ordinary finite domain boundaries. `exprel` and `cosm1` use
cancellation-resistant formulas near zero. The Box-Cox routines use
`expm1`/`log1p`-style kernels so very small `lambda` does not collapse to a
naive subtractive formula.

The convenience helpers are pure elemental real interfaces except for the
integer order of `diric`. `cbrt` uses a logarithmic seed plus Newton correction;
`exp10` and `exp2` use guarded real powers; public `expm1`/`log1p` expose the
existing cancellation-resistant SciFort kernels; and `powm1` switches between a
log-domain `expm1` path near one and direct real powers away from cancellation.
`round` implements ties-to-even and returns a real result. `sinc` is normalized
as `sin(pi*x)/(pi*x)` and returns exact zero at nonzero binary64 integers.
Degree trigonometric functions reduce every finite input modulo 360 degrees and
return exact quadrant values. This intentionally differs from SciPy 1.17.0 at
some extremely large finite angles where its legacy backend returns zero rather
than reducing the represented angle. Nonfinite degree arguments return NaN.
`softplus` uses the stable `log(1+exp(x))` branch already present in SciFort, and
`diric(x,n)` implements the periodic Dirichlet-kernel ratio for integer `n > 0`.


`airy` and `airye` are pure elemental subroutines because the SciPy operation returns
four values. `airy(x, ai, aip, bi, bip)` returns Ai(x), Ai'(x), Bi(x), and Bi'(x).
The overloaded `airy(z, ai, aip, bi, bip)` accepts a complex binary64 argument.
Exact complex values on the real axis reuse the existing real kernel. Elsewhere SciFort
uses the Airy-ODE power series for moderate magnitude and DLMF 9.7 asymptotics plus
DLMF 9.2 rotated-Ai connection formulas in the outer plane; Bi and Bi' are reconstructed
from stable rotated Ai/Ai' evaluations. `airye` is overloaded for real and complex
binary64 arguments. For complex `z`, it follows SciPy's scaling
`eAi=Ai*exp(zeta)`, `eAip=Ai'*exp(zeta)`, and
`eBi/eBip=(Bi/Bi')*exp(-abs(real(zeta)))`, where `zeta=2*z*sqrt(z)/3` on the
principal square-root branch. Moderate complex Ai/Ai' values in the recessive
sector use Airy-ODE Taylor continuation from a scaled asymptotic anchor to avoid
cancellation in the unscaled power series; outer values apply scaling before the
rotated-Ai connection formulas. Exact nonnegative complex real-axis inputs reuse
the real scaled kernel. Negative complex real-axis inputs retain the complex branch
rather than the NaN convention of the real specific. No complex Airy C ABI is
currently exposed.
For `abs(x) <= 1`, the real kernel sums the two independent power-series solutions of
`y'' = x*y`; outside that interval it uses the DLMF real-axis relations to SciFort's
real-order J/I/K Bessel kernels. `airye` returns the exponentially scaled forms used
by SciPy: for positive real `x`, Ai/Ai' are multiplied by `exp(2*x**1.5/3)` and
Bi/Bi' by the reciprocal factor. Its positive branch is evaluated directly with
`kve`/`ive`, so large finite arguments do not require forming the exponential.
For negative real `x`, SciPy 1.17.0's real `airye` returns NaN for the scaled Ai/Ai'
components; SciFort follows that behavior, while the Bi/Bi' components equal their
unscaled real values. Nonfinite real arguments return NaN. For very large positive
finite arguments where `2*x**1.5/3` itself would overflow, `airye` uses the leading
DLMF asymptotic forms; at those magnitudes the omitted corrections are below binary64
resolution. The unscaled `airy` returns zero for underflowing Ai/Ai' and positive
infinity for overflowing Bi/Bi', which is a deliberate mathematical range extension
over SciPy/AMOS cases that can return NaN on overflow.

`itairy(x, apt, bpt, ant, bnt)` is a pure elemental subroutine matching SciPy's
four real integrals: `apt = integral_0^x Ai(t) dt`, `bpt = integral_0^x Bi(t) dt`,
`ant = integral_0^x Ai(-t) dt`, and `bnt = integral_0^x Bi(-t) dt`. Near the
origin SciFort integrates the Airy-ODE power series term by term; a composite
16-point Gauss-Legendre bridge covers the moderate oscillatory region; and DLMF
9.10 asymptotic antiderivative expansions cover the tails. For infinite real
limits SciFort returns the mathematical DLMF limits (`1/3`, `2/3`, and `0`) and
the expected signed divergence of the positive-Bi integral, rather than
propagating implementation NaNs.

`ai_zeros(nt, a, ap, ai, aip[, status])` and
`bi_zeros(nt, b, bp, bi, bip[, status])` allocate and return the first `nt`
negative real zeros of the function and its derivative together with the
cross-values returned by SciPy. DLMF 9.9 leading asymptotic zero locations seed
Newton refinement through the existing `airy` kernel. `nt` must be positive;
invalid input returns zero-length arrays and status `airy_status_invalid_argument`
when `status` is present. The status constants exported by `scifort_special` are
`airy_status_success`, `airy_status_invalid_argument`, and
`airy_status_numerical_failure`. Complex `airy` and `airye` arguments are supported;
complex Bi zeros remain outside the current API.


`exp1(x)` evaluates the real exponential integral E1 on `x >= 0`; negative
real arguments return NaN. `expi(x)` evaluates the real principal exponential
integral Ei(x), using `Ei(x) = -E1(-x)` on the negative axis. The overloaded
`exp1(z)` and `expi(z)` accept complex binary64 arguments. Complex E1 uses the
principal branch with its cut on the negative real axis; exact `-x+0i` and
`-x-0i` values select the lower and upper limiting values respectively, matching
SciPy's signed-zero convention. Complex Ei has the same analytic branch cut away
from the axis, while exact complex values on the negative real axis return the
real principal-value branch, matching SciPy's documented convention.

`expn(n,x)` provides the integer-order generalized exponential integral for
nonnegative integer `n` and real `x >= 0`; negative order or argument returns
NaN. The real SciFort kernels use a convergent small-argument series and a
modified-Lentz continued fraction for E_n, plus a positive-Ei power
series/asymptotic expansion with log-domain overflow guarding. The complex E1/Ei
kernel is an independent implementation of the DLMF Ein series and E1 continued
fraction, with guarded complex exponential scaling in the large-argument region.

`sici` and `shichi` are pure elemental subroutines overloaded for real and complex
binary64 arguments and preserve SciPy's two-result ordering. The real `sici`
specific adapts the MATH77 DSI/DCI Chebyshev approximations, retaining the real
even continuation of Ci for finite negative arguments; at negative infinity it
returns NaN for Ci to match SciPy 1.17.0. Real `shichi` uses the standard
identities in terms of Ei(abs(x)) and E1(abs(x)); Shi is odd and Chi is even on
finite real arguments. The complex specifics use the DLMF defining power series
for small arguments and reuse the validated principal E1/Ei kernels outside that
disk. Complex Si and Shi are entire. Complex Ci and Chi use the principal logarithm
and therefore carry a negative-real-axis branch cut; exact complex negative-real
arguments return the upper-cut value with `+i*pi`, matching SciPy's convention.
No complex C ABI is exposed.

`lambertw(z,k,tol)` solves `w*exp(w)=z` for integer branch `k` (default zero) and returns a complex binary64 result for either real or complex `z`. The implementation follows the DLMF branch structure, uses the square-root expansion around `z=-1/e` on the principal and real `k=-1` branches, uses the standard logarithmic large-argument expansion elsewhere, and applies Halley refinement to the defining equation. A positive optional `tol` may be supplied, but SciFort never relaxes convergence below a machine-scale floor. Exact `z=-1/e` returns the mathematical branch-point value `-1`; this intentionally differs from the NaN observed from SciPy 1.17.0 at that exact binary64 input. At `z=0`, the principal branch is zero and nonprincipal branches return negative infinity in the real component.

`wrightomega(z)` solves `w + log(w) = z` using the principal logarithm. Real and complex overloads are available. In the safe exponential range SciFort evaluates the equivalent branch-selected Lambert-W relation; for very large positive or noncentral negative real parts it refines the defining equation directly to avoid overflow or destructive underflow. The exact branch points `-1 +/- i*pi` return `-1`. Values extremely close to those branch points are ill-conditioned in binary64 and may lose several digits.

`lpmv(m,v,x)` evaluates the real Ferrers associated Legendre function for an integer order `m`, real binary64 degree `v`, and `-1 <= x <= 1`, including the Condon-Shortley phase. Exact integer degrees use the classical associated-Legendre recurrence, including `P_{-n-1}^m=P_n^m`, so gamma-pole cancellations are avoided. For noninteger degrees, nonpositive orders use the regularized DLMF Gauss-hypergeometric representation and positive orders use the differentiated hypergeometric representation; this is substantially more accurate near `x=-1` than forming a large positive/negative-order gamma ratio. Arguments outside the real Ferrers interval return NaN.

`legendre_p(n,z)` evaluates the ordinary Legendre polynomial `P_n(z)` for a nonnegative integer degree and real or complex binary64 `z`. `legendre_p_derivatives(n,z,diff_n,values)` returns derivatives with respect to `z` in `values(0:diff_n)`. The derivative implementation propagates a truncated Taylor series through the three-term Legendre recurrence rather than using finite differences. For polynomial degree `n`, derivatives above order `n` are exactly zero; invalid negative degrees return NaN.

The four current SciPy-style table APIs are also available. `legendre_p_all(n,z,values[,status,diff_n])` allocates `values(0:diff_n,0:n)`. `assoc_legendre_p_all(n,m,z,values[,status,branch_cut,norm,diff_n])` and `sph_legendre_p_all(n,m,theta,values[,status,diff_n])` allocate `values(0:diff_n,0:n,-m:m)`, so the final Fortran subscript is the mathematical signed order directly rather than NumPy's wrapped storage order. `sph_harm_y_all(n,m,theta,phi,values[,status,diff_n])` allocates `values(0:diff_n,0:diff_n,0:n,-m:m)`; `values(i,j,l,k)` is the mixed derivative `d**(i+j)Y_l^k/(dtheta**i*dphi**j)` when `i+j <= diff_n`, with the unused upper triangle set to zero. This uniform tensor corresponds to SciPy's value/gradient/Hessian tuple for derivative orders zero through two. SciFort additionally permits higher derivative orders. `diff_n` defaults to zero. Invalid table dimensions, negative derivative order, or an invalid associated-Legendre branch choice produce an empty allocatable result and `associated_legendre_status_invalid_argument` when `status` is present.

`assoc_legendre_p(n,m,z[,branch_cut,norm])` evaluates the integer-degree associated Legendre polynomial for real or complex binary64 `z`. `branch_cut=2` (default) uses the square-root factors `sqrt(1-z)*sqrt(1+z)`, placing the cut on the real axis outside `[-1,1]`; `branch_cut=3` uses `sqrt(z-1)*sqrt(z+1)`, placing the cut on `(-1,1)`. This factorized form preserves the upper/lower-cut continuation and the signed-zero behavior of the principal square roots. `norm=.true.` multiplies by `sqrt((2*n+1)*(n-m)!/(2*(n+m)!))`; negative orders use the corresponding branch-dependent connection relation. Real overloads return NaN when the selected branch requires a complex value.

`sph_legendre_p(n,m,theta)` returns the fully normalized integer-degree spherical associated Legendre value used by SciPy's current spherical-harmonic API. `sph_harm_y(n,m,theta,phi)` multiplies that amplitude by `exp(i*m*phi)`. Both preserve the Condon-Shortley phase, return zero for `n < 0` or `abs(m) > n`, and are elemental.

Fortran-friendly derivative procedures are provided because SciPy's tuple/tensor `diff_n` return convention has no direct scalar-function analogue in Fortran. `assoc_legendre_p_derivatives(...,diff_n,values,...)` stores the value and derivatives with respect to `z` in `values(0:diff_n)`. `sph_legendre_p_derivatives(...,diff_n,values)` stores derivatives with respect to `theta`. `sph_harm_y_derivatives(...,diff_n,values)` stores mixed derivatives as `values(i,j) = d**(i+j)Y/(dtheta**i*dphi**j)` for `i+j <= diff_n`; other entries are zero. The implementation propagates truncated Taylor coefficients through the associated-Legendre recurrence, so ordinary non-branch-point inputs support derivative orders beyond two without finite differencing. Derivatives exactly at square-root branch points can be singular or branch-sided and are not a portability target in this checkpoint.

`lqn(n,z,qn,qnd[,status])` returns `Q_k(z)` and `dQ_k(z)/dz` for degrees `k=0,...,n`; `lqmn(m,n,z,qmn,qmnd[,status])` returns associated second-kind values and derivatives for orders `0,...,m` and degrees `0,...,n`. Output arrays use zero lower bounds so the Fortran indices are the mathematical order/degree. Real and complex overloads are available. Inside the unit disk the implementation uses the Ferrers seed `Q_0(z)=atanh(z)` and DLMF degree/order recurrences. Outside the unit disk each order-zero `Q_n` is anchored independently with the DLMF 14.3.7 hypergeometric representation before the order recurrence is applied, avoiding the unstable forward degree recurrence for the recessive exterior solution. Derivatives use the exact DLMF 14.10 relation rather than finite differences. The singular points `z=+/-1` return NaNs. On nonreal points exactly on the unit circle SciFort evaluates the canonical exterior complex-Q continuation; this intentionally avoids the infinities/NaNs observed from SciPy 1.17.0's legacy `lqn`/`lqmn` backend there. No C ABI is currently exposed for these routines.

`hyp0f1(v,z)` provides the generalized hypergeometric function 0F1 with a real denominator parameter and either real or complex binary64 argument. Nonpositive integer `v` values are poles and return NaN. The implementation evaluates the defining series in the origin disk. Larger real arguments use the DLMF 10.16.9 Bessel connection together with SciFort's existing real J/I kernels, combining Gamma and power factors in the log domain. Larger complex arguments continue the regular solution outward with local Taylor series whose coefficient recurrence is derived from the DLMF 16.8 generalized-hypergeometric differential equation. This avoids the cancellation of summing the origin series directly at large oscillatory complex arguments. The current complex validation emphasizes finite `v` from about -10.5 through 25 and `abs(z)` through about 630; more extreme negative parameters with very large complex arguments remain a future accuracy target.

`hyp1f1(a,b,z)` provides Kummer's confluent hypergeometric function 1F1 for real binary64 parameters `a` and `b` and either real or complex binary64 arguments. The denominator parameter is singular at nonpositive integers. For the real overload, SciPy-compatible terminating cases with a negative-integer numerator that end before the denominator pole are evaluated as finite polynomials; the complex overload follows SciPy's convention and reports the denominator pole directly. The origin region uses the defining series. Moderate arguments use local Taylor continuation of the Kummer differential equation, with Kummer's transformation used in the left half-plane to reduce cancellation. Large real arguments use log-domain DLMF asymptotics, including a terminating transformed polynomial when the algebraic large-negative branch vanishes. Large nonreal arguments use the two-branch DLMF 13.7 expansion with upper/lower-half-plane phase selection. The current validation covers real parameters roughly from -8 through 12 and real arguments through magnitude 700, plus complex arguments through magnitude about 170; extreme parameter/argument combinations near Stokes transitions remain a future accuracy target.

`hyp2f1(a,b,c,z)` provides Gauss' hypergeometric function 2F1 for real binary64 parameters and either real or complex binary64 arguments. The origin disk uses the defining Gauss series. Outside that disk SciFort continues the regular solution with local Taylor expansions whose coefficient recurrence is derived directly from `z*(1-z)*y'' + (c-(a+b+1)*z)*y' - a*b*y = 0`. Each step stays inside the disk bounded by the regular singular points at zero and one. Complex targets on or beyond the positive-real branch point are reached through an upper- or lower-half-plane detour; signed zero in the imaginary component selects the corresponding side of the cut. Exact `z=1` uses Gauss' gamma-ratio value when `c-a-b > 0` and the signed divergent limit otherwise. Negative-integer numerator parameters use terminating polynomials, while exact `a=c` or `b=c` denominator-pole cancellations use the identity `(1-z)**(-b)` or `(1-z)**(-a)`. The real overload returns SciPy-style positive infinity beyond `z=1` when the principal continuation is genuinely complex. Nonfinite nonpolynomial arguments currently return NaN rather than attempting an asymptotic infinity limit.

`hyperu(a,b,x)` provides Tricomi's confluent hypergeometric function U for real binary64 parameters and the principal real branch `x >= 0`. Moderate and large positive arguments start from the DLMF large-`x` asymptotic solution `U(a,b,x) ~ x**(-a)` and continue inward with local Taylor expansions of `x*y'' + (b-x)*y' - a*y = 0`. The continuation step is one percent of the current positive abscissa, keeping every Taylor disk away from the regular singular point at zero. Small positive arguments use the DLMF connection formula in terms of `hyp1f1`; values close to integer `b` use symmetric Richardson extrapolation in `b` so the removable singularities of the connection formula are not evaluated directly. Negative-integer `a` uses the exact terminating generalized-Laguerre recurrence. Negative real `x` returns NaN because SciPy exposes only the real principal branch across the negative-axis cut. At `x=0`, SciFort returns the mathematical one-sided limit, including its sign when divergent; this intentionally differs from SciPy 1.17.0 cases that return positive infinity even when the limit is negative. No complex `hyperu` overload or C ABI is currently exposed.

`pbdv(v,x,d,dprime)` evaluates the real parabolic-cylinder function `D_v(x)` and its first derivative. For positive `x`, SciFort uses the DLMF confluent-hypergeometric relation `D_v(x) = 2**(v/2)*exp(-x*x/4)*U(-v/2,1/2,x*x/2)`, reusing the `hyperu` kernel; the derivative follows the adjacent-order identity `D_v'(x) = x*D_v(x)/2 - D_(v+1)(x)`. Negative arguments are continued from the exact DLMF gamma-ratio values at zero with Taylor coefficients derived directly from the parabolic-cylinder differential equation. Nonnegative integer orders on the negative axis instead use the Hermite three-term recurrence, preventing roundoff from exciting the dominant non-Hermite solution. The interface is a pure elemental subroutine, so scalar and conformable array arguments are supported. Extremely large simultaneous order and positive argument can still encounter an intermediate overflow in the Kummer-U representation even when a scaled final value would underflow, and remain a future scaled-kernel target.

`pbvv(v,x,vv,vprime)` evaluates SciPy's real parabolic-cylinder `V_v(x)` normalization and its first derivative. Exact reciprocal-gamma values at zero seed Taylor continuation on the positive axis; sufficiently large positive arguments use the DLMF asymptotic expansion, and negative arguments use the exact D/V connection formula with reduced `sin(pi*v)` and `cos(pi*v)` evaluation so integer and half-integer cancellations are exact in binary64. `pbdv_seq` and `pbvv_seq` return allocatable value/derivative arrays following SciPy's sequence-order convention, including the historical minimum length of two. `pbdn_seq` returns the analytic complex integer-order sequence `D_0(z),...,D_n(z)` and derivatives using the Hermite recurrence; `n=0` likewise returns orders zero and one for sequence-shape compatibility. Negative `n` is rejected with `parabolic_cylinder_status_invalid_argument`. The complex `pbdn_seq` values intentionally follow the analytic DLMF definition rather than a nonanalytic behavior observed in SciPy 1.17.0's legacy complex backend.

`pbwa(a,x,w,wprime)` evaluates the real DLMF parabolic-cylinder function `W(a,x)` and its first derivative on SciPy's legacy domain `abs(a) <= 5`, `abs(x) <= 5`. Exact DLMF gamma-ratio values at the origin seed Taylor continuation whose recurrence is derived directly from `y'' + (x*x/4-a)*y = 0`. Inputs outside the legacy domain, NaNs, and infinities return NaN for both outputs. The interface is pure elemental. Regression values are generated independently from the DLMF confluent-hypergeometric representation rather than treating the legacy SciPy backend as the mathematical reference; this matters near the edge of the documented accuracy range.


`ellip_harm(h2,k2,n,p,s[,signm,signn])` evaluates the real first-kind
ellipsoidal harmonic (Lame function) for `h2 > 0`, `k2 > h2`, degree
`n >= 0`, and mode `1 <= p <= 2*n+1`. Optional `signm` and `signn` are
each `+1` or `-1` and select the square-root signs used by the L/M/N
classes, matching SciPy's convention. The finite class eigenproblem and
coefficient scaling are adapted from SciPy 1.17.0's BSD-2-Clause
`_ellip_harm.pxd` implementation, but SciFort solves the symmetric
eigenproblem with its existing pure Jacobi eigensolver rather than LAPACK.
The routine is pure elemental. Invalid parameters return NaN.

`ellip_harm_2(h2,k2,n,p,s)` evaluates the corresponding real second-kind
harmonic for finite exterior arguments `abs(s) >= sqrt(k2)`. SciFort uses
the defining integral after the substitution `u = sin(theta)/abs(s)`,
which removes the finite endpoint and is integrated by an independently
written adaptive 16-point Gauss-Legendre rule. `ellip_normal(h2,k2,n,p)`
returns the ellipsoidal-harmonic normalization constant. Its two
endpoint-weighted defining integrals are likewise transformed to smooth
angle integrals before quadrature. These two quadrature implementations
do not copy SciPy's QUADPACK/Cython integration code. Inputs outside their
real domains return NaN; all three interfaces use binary64 arithmetic and
are pure elemental.


`elliprc`, `elliprd`, `elliprf`, `elliprg`, and `elliprj` provide binary64
Carlson symmetric elliptic forms RC, RD, RF, RG, and RJ as pure elemental
functions for both real and complex arguments. The real specifics retain their
existing domains: `elliprc(x,y)` accepts `x >= 0`, `y /= 0` and evaluates
negative `y` as a Cauchy principal value; RF/RD/RG use the nonnegative real
Carlson domains; and negative real `p` in RJ uses Carlson's real
Cauchy-principal-value continuation.

For complex arguments, RC/RF/RD/RG follow the principal square-root branches
on the Carlson cut plane. RF/RD/RG reject symmetric arguments lying exactly on
the negative-real cut. RC permits a negative-real second argument through the
Cauchy-principal-value transformation, including when its first argument is
complex. RJ accepts the sufficient complex domains used by SciPy, including
nonnegative-real symmetric arguments, a conjugate pair plus a nonnegative-real
third argument, and the generic right-half-plane case with `Re(p) > 0`;
unsupported branch configurations return complex NaN rather than selecting an
unverified continuation. Mathematically divergent boundaries return positive
infinity. Infinite-magnitude limits follow SciPy's complex conventions: RC,
RF, RD, and RJ tend to zero, while RG tends to positive infinity when the
infinite direction is well defined. No complex C ABI is exposed.

The real RC/RD/RF/RJ kernels adapt the fifth-order duplication algorithms in
Jacob Williams' BSD-3-Clause `carlson-elliptic-integrals` release 1.0.0,
itself a modern refactoring of the public-domain SLATEC routines. SciFort
normalizes inputs by homogeneity before evaluation. Complex RC/RF/RD use
independently written Carlson duplication formulas with principal roots, and
complex RG uses the RF/RD identity while keeping principal square roots as
separate factors to preserve the branch. Complex RJ uses the defining Carlson
integral after homogeneity scaling and the substitution `t=(u/(1-u))**2`,
with adaptive 16-point Gauss-Legendre integration; this avoids the additional
branch bookkeeping required by a direct complex extension of the legacy RJ
recurrence. Extremely unbalanced finite arguments can require more adaptive RJ
subintervals than ordinary inputs.

`ellipk(m)`, `ellipkm1(p)`, and `ellipe(m)` provide the real complete
Legendre elliptic integrals using SciPy's parameter convention `m = k**2`;
`ellipkm1(p)` accepts the complementary parameter `p = 1 - m` and is the
preferred interface when `m` is extremely close to one. `ellipkinc(phi,m)`
and `ellipeinc(phi,m)` provide the corresponding real incomplete integrals
for arbitrary real amplitudes in radians. All five interfaces are pure
elemental functions. The real branch requires `m <= 1`, while `ellipkm1`
requires `p >= 0`; invalid real-domain inputs return NaN.

The implementation uses the DLMF 19.25 Carlson connection formulas and
SciFort's existing RF/RD/RG kernels. Complete E switches to the convergent
DLMF 19.12.2 complementary-parameter series near `m = 1`, preserving
binary64 accuracy for tiny `1-m` and after the negative-parameter
transformation. Incomplete amplitudes are reduced by `atan2(sin(phi),
cos(phi))` plus exact pi-period bookkeeping, which remains stable for very
large finite amplitudes. At `m = 1`, crossings of the logarithmic K
singularity return infinity with the orientation sign; this differs from
SciPy 1.17.0's positive-infinity result for negative crossings but follows
the signed real integral. Complex continuation for `m > 1` is not provided.

`ellipj(u,m,sn,cn,dn,ph)` is a pure elemental subroutine returning the real
Jacobi elliptic functions `sn(u|m)`, `cn(u|m)`, `dn(u|m)`, and the Jacobi
amplitude `ph = am(u|m)`. The parameter follows SciPy's convention and must
satisfy `0 <= m <= 1`; invalid parameters, NaN inputs, and infinite `u` return
NaN outputs. The implementation uses the DLMF 22.20 AGM inversion after
reducing finite `u` modulo the real period. Exact `m=0` and `m=1` branches use
trigonometric and hyperbolic identities, respectively. `dn` is reconstructed
as `sqrt(cn**2 + (1-m)*sn**2)` to avoid cancellation near `m=1`. The amplitude
is unwrapped so that `ellipkinc(ph,m)` recovers `u` modulo ordinary floating-
point roundoff. Because the real period depends on the binary64 value of
`ellipk(m)`, phase accuracy gradually degrades when `abs(u)` contains many
periods; 100-digit reference checks at `abs(u)=1e6` show absolute errors of
order `1e-10`, while much larger arguments can lose correspondingly more
phase digits. Complex arguments and parameters outside `[0,1]` are not yet
supported.

`i0e(x)` and `i1e(x)` evaluate `exp(-abs(x))*I0(x)` and
`exp(-abs(x))*I1(x)`, respectively. The shared dependency-free implementation
uses convergent power series at moderate arguments and scaled asymptotics for
large magnitudes; `i0e` is even and `i1e` is odd. These kernels are currently
used by the Rice and von Mises families. `i0` and `i1` reconstruct the
unscaled functions from those kernels in the log domain, which delays overflow
until the binary64 result itself is too large. Consequently, SciFort can return
finite `i0`/`i1` values for some large finite arguments where SciPy 1.17.0's
Cephes-backed fast functions already overflow.

`j0`, `j1`, `y0`, and `y1` wrap the corresponding Fortran standard Bessel
intrinsics for finite real arguments. The Y functions accept only nonnegative
arguments, return `-infinity` at zero, and return NaN for negative arguments.
The fast `k0`, `k1`, `k0e`, and `k1e` interfaces reuse SciFort's existing
real-order log-K kernels; the `e` variants return `exp(x)*K_n(x)`. Infinite J,
Y, and unscaled-I inputs return NaN to match the SciPy fast-function behavior,
while K functions tend to zero at positive infinity.

`itj0y0(x,ij0,iy0)` returns `integral_0^x J0(t) dt` and
`integral_0^x Y0(t) dt`. `it2j0y0(x,ij0,iy0)` returns
`integral_0^x (1-J0(t))/t dt` and `integral_x^infinity Y0(t)/t dt`.
The modified-Bessel counterparts `iti0k0` and `it2i0k0` return, respectively,
`integral_0^x I0(t) dt` with `integral_0^x K0(t) dt`, and
`integral_0^x (I0(t)-1)/t dt` with `integral_x^infinity K0(t)/t dt`.
These are pure elemental subroutines because each SciPy API produces a pair of
values. The J/I first components are extended by their real parity where
applicable; Y/K components return NaN for negative arguments. Divergent x=0
tail integrals return signed IEEE infinity. Convergent infinite-argument limits
use their mathematical values.

`jv(v,z)` and `yv(v,z)` accept finite real order with real or complex binary64
arguments through `scifort_special`. Real arguments preserve the established MATH77-based
kernels and their real-domain restrictions. Complex arguments use the principal branch,
with the negative real axis interpreted from above. Away from the large-order/small-argument
regime, J and Y are reconstructed from the validated scaled Hankel pair so that `jve` and
`yve` apply SciPy's `exp(-abs(aimag(z)))` scaling without first forming exponentially large
unscaled Hankel values. When a nonnegative order (or an equivalent negative integer order)
is large relative to `abs(z)`, J instead uses its defining power series directly; this avoids
catastrophic cancellation between huge Hankel terms when the true J value is tiny. `jn(v,z)`
remains a real-order alias of `jv(v,z)`, matching SciPy 1.17.0; unlike SciFort `yn`, its order
dummy is real, not integer. For real noninteger orders the adapted MATH77 phase-reduction
limit on very large real `x` remains unchanged.

`jvp`, `yvp`, `ivp`, and `kvp` accept the same real or complex argument forms. Their optional
integer `n` is the derivative order (default 1); `n=0` returns the corresponding base function.
Both real and complex implementations use the DLMF 10.6.7 and 10.29.5 shifted-order finite
sums with normalized binomial weights and scaled accumulation. Because Fortran has no scalar
exception analogue to SciPy's `ValueError`, a negative `n` returns NaN (complex NaN for a
complex argument).

`hankel1(v,z)` and `hankel2(v,z)` accept finite real order and finite real or
complex nonzero arguments and return complex binary64 values. Real arguments reuse
SciFort's independently validated `jv`/`yv` kernels, with the negative real axis
interpreted as the upper edge of the principal branch. For nonreal arguments,
SciFort evaluates the scaled Hankel functions first. The recessive half-plane
solution is connected to a scaled complex `K_v` through DLMF 10.27.8; the scaled
K and complementary scaled Hankel solution are continued from positive-real
anchors by the corresponding Bessel differential equations along constant-radius
arcs. Left-half-plane values use DLMF 10.11.3-10.11.5 continuation from `-z`,
which avoids integrating around the branch point. `hankel1e` and `hankel2e`
apply the SciPy scale factors `exp(-i*z)` and `exp(+i*z)` directly and are the
preferred interfaces when an unscaled value would overflow or underflow.
`h1vp` and `h2vp` accept an optional nonnegative derivative order `n` (default
1); complex derivatives use the DLMF 10.6.7 shifted-order finite sum with
normalized binomial weights. A negative derivative order returns a complex NaN.
`jn_zeros`, `jnp_zeros`, `yn_zeros`, and `ynp_zeros` return the first `nt`
positive real zeros of integer-order J, J', Y, and Y', respectively;
`jnyn_zeros` returns all four tables together. Negative integer orders use the
same zeros as their absolute orders. The zero at the origin is excluded,
matching SciPy's positive-interval convention. Successive roots are bracketed
with a fixed step smaller than the real-root spacing and refined by safeguarded
Newton iteration using the existing Bessel derivative kernels. `nt < 1` returns
an empty array (or arrays) with `bessel_zeros_status_invalid_argument` rather
than raising Python's `ValueError`. The implementation has no artificial
1200-root ceiling; compatibility has been regression-tested through SciPy's
documented 1200-root combined-table scale.
`y0_zeros`, `y1_zeros`, and `y1p_zeros` are convenience interfaces over these
same real zero kernels and additionally evaluate the companion derivative or
function values at each root.

`yn(n,x)` and `kn(n,x)` add
integer-order Y and K with an integer Fortran order argument; negative integer orders
follow the standard parity/even-order identities. The `scifort_special` overloads
`iv(v,z)`, `ive(v,z)`, `kv(v,z)`, and `kve(v,z)` accept finite real order and either real or
complex binary64 arguments. Real arguments preserve the existing real-domain behavior.
Complex I uses the DLMF J/I quarter-turn connection, selecting the rotation by half-plane so
the principal branch is preserved; `ive` applies SciPy's `exp(-abs(real(z)))` scale through
the correspondingly scaled J evaluation. Complex K uses the DLMF K/Hankel quarter-turn
connection and the already scaled Hankel functions, so `kve` computes `exp(z)*K_v(z)`
directly. The negative real axis uses the upper edge of the principal branch, and `K_-v=K_v`
is inherited from the connection formulas. The spherical interfaces take a nonnegative integer
order `n`, a real or complex binary64 argument, and an optional logical `derivative`
(default false). They are pure elemental functions. Real arguments preserve the
existing specialized kernels: `spherical_jn` and `spherical_in` use small-argument
power series and recurrence/logarithmic kernels chosen to preserve tiny high-order
values; `spherical_yn` uses the standard J/Y-style recurrence; and `spherical_kn`
evaluates its terminating half-integer form in logarithmic space. For genuinely
nonreal arguments, the public overloads generally use the DLMF half-integer
definitions in terms of the complex cylindrical `jv`/`yv`/`iv`/`kv` kernels.
High-order spherical J near its turning region uses a normalized backward Miller
recurrence, and the shared complex-I kernel switches to its direct defining series
where quarter-turn reconstruction would lose accuracy. First derivatives use the
standard adjacent-order spherical identities. A nonzero complex argument on the
real axis delegates to the real kernel, preserving its parity and negative-axis
behavior exactly. Negative orders and nonfinite complex arguments return NaN. At
complex zero, J/I use their regular limits while Y/K return complex NaN, matching the
SciPy complex-argument convention; the real overloads retain their established signed
infinities at singular Y/K endpoints. `riccati_jn` and `riccati_yn` return
allocatable zero-based arrays indexed `0:n` containing `x*j_k(x)` or
`x*y_k(x)` and their first derivatives. Invalid order or nonfinite input returns
empty arrays and `spherical_bessel_status_invalid_argument` when `status` is
present. At `x=0`, `riccati_jn` uses its exact regular limits; `riccati_yn(0)`
is `-1` with zero derivative and higher orders use the mathematical signed
IEEE infinities rather than finite compatibility sentinels.

`struve(v,x)` and `modstruve(v,x)` implement real-order Struve H and modified
Struve L for real binary64 arguments. Negative real `x` is supported for
integer order through the parity identities; a negative argument with
noninteger order returns NaN because the result is generally complex. Exact
negative half-integer orders use their Bessel J/I identities. The implementation
uses convergent defining series near the origin, a Bessel-J expansion for the H
transition region, and DLMF large-argument difference asymptotics.
`itstruve0(x)` is `integral_0^x H_0(t) dt`, `it2struve0(x)` is
`integral_x^infinity H_0(t)/t dt`, and `itmodstruve0(x)` is
`integral_0^x L_0(t) dt`; the H integrals combine exactly integrated local
series with composite eight-point Gauss-Legendre continuation. The routines are
pure elemental scalar functions.

`log_besselk` evaluates
`log(K_nu(x))` for real finite order and positive real argument using a scaled
integral representation with adaptive Gauss-Legendre panels and a large-`x`
asymptotic branch. `log_besselk_scaled(nu,x)` returns `log(exp(x)*K_nu(x))`
directly, avoiding cancellation and overflow in Bessel-normalized heavy-tail
families. The two derivative routines return derivatives of
`log(K_nu(x))` with respect to order and argument. `hurwitz_zeta` evaluates the real
positive-q Hurwitz zeta for finite `s > 1`, `q > 0` using Euler-Maclaurin
summation; `hurwitz_zeta_derivative` returns its derivative with respect to
`s`. The public `zeta(s[,q])` wrapper supplies Riemann analytic continuation
and the additional real Hurwitz conventions described above.

Stable real/complex reductions are available for vectors through rank-15 arrays;
weighted `logsumexp` covers the same rank range:

```text
logsumexp(x(:) [, b])                                -> scalar
logsumexp(x(rank 2-15) [, b])                        -> scalar
logsumexp(x(rank 2-15), dim [, b])                   -> rank(x)-1 array
logsumexp_keepdims(x(rank 2-15), dims(:) [, b])      -> same rank as x
call logsumexp_signed(x, logabs, sign_or_phase [, b])
call logsumexp_signed(x(rank 2-15), dim, logabs, sign_or_phase [, b])
call logsumexp_signed_keepdims(x(rank 2-15), dims, logabs, sign_or_phase [, b])
softmax(x(:))                                        -> rank-1 array
softmax(x(rank 2-15) [, dim])                        -> same shape as x
softmax(x(rank 2-15), dims(:))                       -> same shape as x
log_softmax(x(:))                                    -> rank-1 array
log_softmax(x(rank 2-15) [, dim])                    -> same shape as x
log_softmax(x(rank 2-15), dims(:))                   -> same shape as x
```

The array forms cover fixed ranks 1 through 15 where the result rank is
representable by a standard Fortran generic interface; explicit `dim`, `dims`, and
keepdims forms start at rank 2. Weighted `logsumexp` forms cover the same ranks.
`softmax` and `log_softmax` are unweighted. All forms use magnitude/real-part
shifting before exponentiation. Complex `log_softmax` follows SciPy's max-shift
branch convention, so its imaginary part is not forced to the principal branch of
`x - logsumexp(x)`.
Weighted `logsumexp` evaluates `log(sum(b*exp(x)))`. For real inputs and real
weights, negative totals return NaN, matching SciPy's default
`return_sign=False` behavior. `logsumexp_signed` is the idiomatic Fortran
equivalent of `return_sign=True`: real inputs with real weights return a sign in
`{-1,0,+1}`, while complex inputs or complex weights return a unit complex phase
and a real log-magnitude. Complex weights are accepted with either real or
complex `x`; a real `x` with complex weights therefore returns a complex
`logsumexp` value.

Weight broadcasting follows NumPy trailing-dimension rules through rank 15 and is
joint between `x` and `b`. A scalar weight broadcasts everywhere. A lower-rank
weight is aligned with the trailing dimensions of `x`; each aligned extent must
match or one of the two extents must be one. Singleton dimensions of `x` may
therefore expand. For example, rank-5 `x(2,1,2,1,3)` and rank-4 `b(3,1,4,3)` jointly
broadcast to `(2,3,2,4,3)`. Zero weights are ignored before the corresponding
exponential is evaluated, so zero weights remain suitable for masking.

A scalar `dim` uses Fortran's one-based dimension numbering: `dim=1` corresponds
to NumPy `axis=0`. The `dims(:)` forms are the Fortran analogue of an axis tuple.
An empty list means no axes, and any unique subset/order of the valid dimensions
is accepted. `softmax` and `log_softmax` keep the input shape naturally; empty
`dims` therefore returns all ones or all zeros, respectively. Because a Fortran
function result rank cannot depend on a runtime `keepdims` flag, SciFort exposes
`logsumexp_keepdims` and `logsumexp_signed_keepdims` as explicit fixed-rank
companions. Reduced dimensions have extent one. Invalid or duplicate dimension
lists return zero-size keepdims results, while invalid `softmax`/`log_softmax`
dimension lists return NaNs.

`logsumexp` of an empty vector is `-infinity`; `softmax` and `log_softmax` retain
the existing empty-array behavior. Rank 15 is the standard Fortran maximum array
rank, so the native fixed-rank reduction API reaches the language ceiling.

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


## Signal processing

The `medfilt` and `wiener` generics support ranks 1-3, scalar or
per-axis window sizes and zero padding. `medfilt2d` is the rank-2 alias.
`wiener` supports real, integer or complex data, with estimated or supplied
real/complex noise pseudo-variance. See
[Signal processing](signal.md) for examples and limitations.

The new `scifort_signal` module implements 1D `find_peaks` (including
position-dependent lower/upper interval bounds for all five constrained
properties), `peak_prominences`, `peak_widths`, `detrend`, `convolve`, `correlate`,
and `correlation_lags`, together with `savgol_coeffs`, `savgol_filter`, `savgol_filter_nd`,
`lfilter`, `lfilter_nd`, and `lfilter_zi`, and with `fftconvolve`, `oaconvolve`,
`periodogram`, `welch`, `csd`, `coherence`, `stft`, and `spectrogram`
for real and complex 1D data. See [Signal processing](signal.md) for signatures,
examples, one-based indexing rules, and SciPy compatibility limitations.


### Polyphase, Fourier, and IIR/FIR resampling

`scifort_signal` additionally exports `upfirdn`, `resample_poly`, and
`decimate`, and Fourier `resample`. The functions accept one-dimensional real(dp) or
complex(dp) signals; both `upfirdn` and `resample_poly` support real/complex FIR coefficients.
Their optional `status=` returns `resample_status_success` (0) or
`resample_status_invalid_input` (1); invalid requests produce a zero-length
result. `resample_poly` reduces the rational ratio by its greatest common
divisor and aligns the first output sample to the first input sample. The
default Kaiser-window prototype uses beta=5, as in SciPy.

```fortran
use scifort_signal, only : dp, upfirdn, resample_poly, decimate
real(dp), allocatable :: y(:)
y = upfirdn([0.25_dp, 0.5_dp, 0.25_dp], x, up=3, down=2)
y = resample_poly(x, 3, 2, window='kaiser')
y = decimate(x, 4, ftype='fir', zero_phase=.true.)
```

For `resample_poly`, complex `taps=` produce complex output even with real
input. Boundary modes include `'antisymmetric'` and `'antireflect'`;
`padtype` additionally supports `'mean'`, `'median'`, `'minimum'`, and
`'maximum'`. The rank-2/3 overloads allow an optional one-based or negative
`axis` argument; see the signal guide for examples and reference tests.

For higher-dimensional data, `scifort_signal` also exports four Fortran
2018 assumed-rank, caller-allocated subroutine generics:
`resample_into`, `resample_poly_into`, `upfirdn_into`, and `decimate_into`.
They support input/output ranks 1-15 with `axis=` and check the output
shape before processing. A nonzero `status=` leaves the output unchanged.
For real input with complex FIR taps, declare a complex output array.

```fortran
real(dp) :: x(2, 1, 9, 3), y(2, 1, 14, 3)
call resample_poly_into(x, y, 3, 2, axis=3, status=status)
```

For real and complex rank-1–15 spectral analysis, `periodogram_into`,
`welch_into`, `csd_into` and `coherence_into` use caller-allocated results
and selected positive-one-based or negative axes, with median Welch/CSD
averaging and mixed-type cross-spectra. See
[`docs/signal_spectral_nd.md`](signal_spectral_nd.md) for examples and testing.

`stft_into` and `spectrogram_into` support rank-1 through rank-14 real/complex
inputs and caller-allocated rank-plus-one time-frequency output arrays.
The legacy `stft` and `stft_into` functions accept zero, even-reflected,
odd-reflected, constant-endpoint, or no boundary extension via `boundary=`.
Both support `axis=`, return frequency/time coordinates, and provide
transactional failure behavior. Spectrogram modes include `psd`, `complex`,
`magnitude`, `angle` and `phase`; see
[`docs/signal_timefreq_nd.md`](signal_timefreq_nd.md).

See `docs/signal.md#assumed-rank-multirate-transforms-ranks-1-15` for the
full signatures, defaults, memory considerations, and regression tests.

The `decimate` implementation now defaults to `ftype='iir'` (eighth-order
Chebyshev-I with ripple 0.05 dB and cutoff 0.8/q), using existing SOS
forward or forward/backward filtering; `ftype='fir'` remains available.
The FFT-based `resample` assumes periodic input and handles even-length
Nyquist bins explicitly. Data-padding, multidimensional axes, and window
restrictions are detailed in `docs/signal.md`. Standalone checked/optimized GNU Fortran
tests and a 1,120-case randomized SciPy 1.17.0 comparison validate the
supported operations. No extra FPM dependencies are introduced.

### Digital filter design and frequency response

The public `scifort_signal` module now also exports `firwin`, `butter`,
`freqz`, `sosfreqz`, `cheby1`, `cheby2`, `iirfilter`, `iirdesign`,
`buttord`, `cheb1ord`, `cheb2ord`, `iirnotch`, `iirpeak`, and the result types `iir_order_result`, `butter_result`, and
`frequency_response`. The routines support windowed-sinc FIR coefficients,
low/high/bandpass/bandstop Butterworth, Chebyshev-I, and Chebyshev-II filters,
notch/peak biquads, scalar lowpass/highpass minimum-order design, and real/complex
frequency-response calculations. See `docs/signal.md` for exact signatures,
status codes, numerical limitations, and examples. SOS output is recommended
for high-order Butterworth filters. No new FPM dependencies are required.

## Hierarchical clustering (`scifort_cluster_hierarchy`)

The `linkage` generic implements the seven classical hierarchical
agglomerative methods from `scipy.cluster.hierarchy`; both raw observations
and condensed pairwise distances are accepted. Companion functions provide
flat clustering by distance or maximum count, the unified `fcluster`
dispatcher, direct observation-to-label `fclusterdata`, cophenetic
distances/correlation, and leaf ordering. See [cluster_hierarchy.md](cluster_hierarchy.md) for API
signatures, one-based index conventions, examples and limitations.


## Compressed-sparse graph algorithms (`scifort_sparse_csgraph`)

`min_weight_full_bipartite_matching(graph,row_ind,col_ind,status,maximize,total_weight)`
solves sparse rectangular weighted assignment, minimizing costs by default or
maximizing them with `maximize=.true.`. Both final arguments are optional;
the total-cost output is a SciFort extension. The existing
`minimum_weight_full_bipartite_matching(graph,row_ind,col_ind,total_weight,status,maximize)`
remains supported. Both use `csr_weighted_bipartite_graph`, return sorted
one-based rows, and return `csgraph_no_full_matching` when infeasible.

The one-based CSR graph module supports `connected_components`,
`breadth_first_order`, `depth_first_order`, `breadth_first_tree`,
`depth_first_tree`, single-source `dijkstra` and `bellman_ford`, all-pairs
`floyd_warshall` and `johnson`, the `shortest_path` dispatcher, Kruskal
`minimum_spanning_tree`, dense `laplacian`, `reconstruct_path` (CSR tree),
and `path_from_predecessors` (SciFort vertex-sequence helper). Johnson and
the dispatcher accept optional one-based source `indices`, and both all-pairs
algorithms return predecessor matrices with zero as the unavailable sentinel.
See [sparse_csgraph.md](sparse_csgraph.md) for examples, exact signatures,
negative-cycle behavior, directedness, normalization, and test commands.

## Mixed-integer linear programming (`scifort_optimize_milp`)

`milp(c, integrality, a, b_l, b_u, lower, upper, node_limit, mip_rel_gap,
tol, lp_maxiter, time_limit)` minimizes `dot_product(c,x)`, with dense
two-sided linear constraints and variable types 0 continuous, 1 integer,
2 semi-continuous, 3 semi-integer. The module and result/status constants are
re-exported by `scifort_optimize`. See the
[detailed MILP API, limitations, and independent tests](optimize_milp.md).

## Spherical splines (`scipy.interpolate`)

`smooth_sphere_bivariate_spline(theta,phi,r,w,s,eps,ntest,npest,status)`
provides native weighted scattered-data smoothing corresponding to SciPy's
`SmoothSphereBivariateSpline`. It selects knots automatically for the weighted
squared-residual target `s`, which defaults to zero. Optional knot-capacity
controls are SciFort extensions. The shared result type supports periodic
queries and coordinate derivatives. Deficient fits remain valid with status 2;
the initial polynomial has rank -1 (not estimated). See
[scattered spherical smoothing](interpolate_smooth_sphere.md) for signatures,
limits, rank semantics, and compatibility details.

`rect_sphere_bivariate_spline(u,v,r,s,pole_continuity,pole_values,pole_exact,pole_flat,status)`
provides native rectangular-grid interpolation (`s=0`) and adaptive smoothing
(`s>0`) corresponding to SciPy's `RectSphereBivariateSpline`.
Each coordinate axis requires at least
four sites; pole options are north/south two-element arrays. It returns the
same `sphere_bivariate_spline_type` described below, with a squared grid
residual and rank -1 (not estimated). Failed convergence returns status 3
and an invalid object. See [rectangular spherical splines](interpolate_rect_sphere.md)
for pole estimation, grid bounds, diagnostics, and compatibility differences.

`lsq_sphere_bivariate_spline(theta,phi,r,tt,tp,w,eps,status)` returns a native
`sphere_bivariate_spline_type`, corresponding to SciPy's
`LSQSphereBivariateSpline`. The function is exported by `scifort_interpolate`
and `scifort_interpolate_sphere`. It fits real weighted scattered data with
bicubic periodic longitude and FITPACK pole coefficient constraints.
Angles are in radians; `tt` and `tp` contain interior knots. At least three
longitude knots are required. Weights and `eps` are optional.

Methods provide `%value`, paired `%values`, Cartesian-product `%grid`,
`%get_knots`, `%coefficients`, and `%residual`. Coordinate derivative orders
`dtheta,dphi` independently range from 0 through 2. `%rank` and
`%degrees_of_freedom` are SciFort extension diagnostics. A rank-deficient
minimum-norm fit remains valid and returns status 2. Dense workspace and
256-coefficient limits apply; query longitude wrapping is a SciFort extension.
See [spherical spline fitting](interpolate_sphere.md) for exact signatures,
status codes, solver differences, limits, and provenance.

## AAA rational approximation (`scipy.interpolate`)

`scifort_interpolate` exports `aaa`, `aaa_type`, and `complex_aaa_type` for
real/complex scalar sampled data, adaptive support selection, scalar/vector
evaluation, poles, zeros, simple-pole residues, and Froissart cleanup.
See the [AAA guide](interpolate_aaa.md) for controls, statuses, convergence,
SciPy compatibility differences, and the 512-support-point spectral limit.

### Native scalar adaptive smoothing

`splrep(x,y,w=...,k=...,s=positive)` now selects knots and fits a nonperiodic
scalar smoothing spline without the optional FITPACK package. Degrees 1-5,
positive residual weights, and enclosing `bbox` bounds are supported. The
existing default `s=0`, tuple layout, evaluation/calculus routines, and status
names are preserved. See the [TCK guide](interpolate_fitpack_tck_backend.md)
for residual tolerance, numerical guards, and remaining parametric/periodic gaps.

### Native parametric adaptive smoothing

`splprep(points,u=...,w=...,k=...,s=positive)` now fits a nonperiodic curve
with shared automatic knots and one total weighted residual target across
channels. Parameters may be supplied or generated from chord lengths. Degrees
1-5 and enclosing `ub`/`ue` bounds are supported; unpadded coefficient storage,
default `s=0`, and existing status names are preserved. More than ten coordinate
channels remain supported as a SciFort extension. See the
[TCK guide](interpolate_fitpack_tck_backend.md) for scaling, conditioning,
iteration-path differences from SciPy, and remaining periodic/modern-interface gaps.
