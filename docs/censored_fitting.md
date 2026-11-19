# Censored likelihoods and distribution fitting

SciFort supports exact, left-censored, right-censored, and interval-censored
likelihoods and bounded maximum-likelihood fits for normal, exponential,
minimum/maximum Weibull, Rayleigh, lognormal, log-Laplace, Gibrat, gamma, Erlang, beta, Cauchy, symmetric/asymmetric Laplace,
Pareto, Lomax, Fisk, Burr III/XII, uniform, power-function, arcsine, Bradford, triangular, trapezoidal, double-gamma, double-Weibull, log-gamma, inverse-gamma, beta-prime, Maxwell, power-normal, power-lognormal, skew-Cauchy, exponential-power, exponentiated Weibull, Mielke, fatigue-life, Johnson SU, generalized exponential, cosine, anglit, semicircular, truncated exponential, Johnson SB, kappa3, log-uniform, truncated Pareto, generalized half-logistic, folded normal, folded Cauchy, generalized Pareto, generalized normal, half generalized normal, R distribution, chi-squared, chi, Nakagami, Student t, F, Jones-Faddy skew-t, inverse Gaussian, Wald, reciprocal inverse Gaussian, generalized gamma, truncated normal, truncated Weibull, alpha, generalized extreme value, Crystal Ball, ARGUS, Pearson III, kappa4, Gauss hypergeometric, Irwin-Hall, real-line von Mises, exponentially modified normal, skew-normal, double Pareto-lognormal, Rice, noncentral chi-squared, noncentral F, wrapped Cauchy, asymptotic Kolmogorov, relativistic Breit-Wigner, generalized inverse Gaussian, normal-inverse Gaussian, generalized hyperbolic, noncentral Student t, Landau, Levy-stable, studentized range, Tukey lambda, one-sided and two-sided finite-sample Kolmogorov-Smirnov, fixed-branch von Mises (a SciFort extension), generalized logistic, Gompertz, inverse Weibull, Moyal, Levy and reflected Levy, logistic, half-Cauchy, half-logistic, half-normal, hyperbolic-secant, and both Gumbel distributions.
Lognormal, Gibrat, gamma, and Erlang censored fits currently require fixed location.
Beta censored fits require fixed location and scale. The definitions follow SciPy's
[`CensoredData`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.CensoredData.html)
and [censored fitting likelihood](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.rv_continuous.fit.html).
Existing array-based likelihood and fitting calls are unchanged.

## Data representation

```fortran
use scifort_stats, only : dp, censored_data, make_censored_data
type(censored_data) :: data
integer :: status

data = make_censored_data(uncensored=[13.0_dp, 22.0_dp, 17.0_dp, 15.0_dp], &
    right=[20.0_dp, 18.0_dp], status=status)
```

The existing `censored_data` type gains two optional allocatable components,
`left(:)` and `interval(:,:)`, after its original `uncensored(:)` and `right(:)`
components. Existing keyword construction and the `n_observations()` method
remain available. Unallocated or zero-length groups represent no observations.
`make_censored_data` copies supplied arrays. Its optional status is 0 for valid
data and 1 for invalid data; invalid input is retained for inspection.
`censored_data_valid(data)` checks validity, and `censored_size(data)` counts rows
across all groups, as does `data%n_observations()`.

- `uncensored`: finite exact measurements.
- `left`: finite upper limits, representing `X <= limit`.
- `right`: finite lower limits, representing `X > limit`.
- `interval`: shape `(n,2)`, lower endpoints in column 1 and upper endpoints
  in column 2. Endpoints must be ordered and not NaN. An infinite lower or
  upper endpoint represents left or right censoring, respectively. At least
  one endpoint must be finite. Equal finite endpoints represent an exact
  observation, rather than a zero-probability interval.

ECDF and logrank still support only exact/right-censored observations. Passing
nonempty left or interval groups to these procedures returns NaN results;
these groups are never silently discarded. This does not add a general
interval-censored empirical-distribution estimator.

## Likelihoods

The procedures are available from `scifort_censored` and `scifort_stats`:

```text
levy_stable_censored_loglikelihood(data, alpha, beta [, loc, scale])
studentized_range_censored_loglikelihood(data, k, df [, loc, scale])
kstwo_censored_loglikelihood(data, n [, loc, scale])
nct_censored_loglikelihood(data, df, nc [, loc, scale])
landau_censored_loglikelihood(data [, loc, scale])
tukeylambda_censored_loglikelihood(data, lam [, loc, scale])
ksone_censored_loglikelihood(data, n [, loc, scale])
vonmises_censored_loglikelihood(data, kappa [, loc, scale])
geninvgauss_censored_loglikelihood(data, p, b [, loc, scale])
norminvgauss_censored_loglikelihood(data, a, b [, loc, scale])
genhyperbolic_censored_loglikelihood(data, p, a, b [, loc, scale])
kstwobign_censored_loglikelihood(data [, loc, scale])
wrapcauchy_censored_loglikelihood(data, c [, loc, scale])
rel_breitwigner_censored_loglikelihood(data, rho [, loc, scale])
rice_censored_loglikelihood(data, b [, loc, scale])
ncx2_censored_loglikelihood(data, df, nc [, loc, scale])
ncf_censored_loglikelihood(data, dfn, dfd, nc [, loc, scale])
exponnorm_censored_loglikelihood(data, k [, loc, scale])
skewnorm_censored_loglikelihood(data, a [, loc, scale])
dpareto_lognorm_censored_loglikelihood(data, u, s, a, b [, loc, scale])
gausshyper_censored_loglikelihood(data, a, b, c, z [, loc, scale])
irwinhall_censored_loglikelihood(data, n [, loc, scale])
vonmises_line_censored_loglikelihood(data, kappa [, loc, scale])
argus_censored_loglikelihood(data, chi [, loc, scale])
pearson3_censored_loglikelihood(data, skew [, loc, scale])
kappa4_censored_loglikelihood(data, h, k [, loc, scale])
alpha_censored_loglikelihood(data, a [, loc, scale])
genextreme_censored_loglikelihood(data, c [, loc, scale])
crystalball_censored_loglikelihood(data, beta, m [, loc, scale])
gengamma_censored_loglikelihood(data, a, c [, loc, scale])
truncnorm_censored_loglikelihood(data, a, b [, loc, scale])
truncweibull_min_censored_loglikelihood(data, c, a, b [, loc, scale])
invgauss_censored_loglikelihood(data, mu_shape [, loc, scale])
wald_censored_loglikelihood(data [, loc, scale])
recipinvgauss_censored_loglikelihood(data, mu_shape [, loc, scale])
t_censored_loglikelihood(data, df [, loc, scale])
f_censored_loglikelihood(data, dfn, dfd [, loc, scale])
jf_skew_t_censored_loglikelihood(data, a, b [, loc, scale])
chi2_censored_loglikelihood(data, df [, loc, scale])
chi_censored_loglikelihood(data, df [, loc, scale])
nakagami_censored_loglikelihood(data, nu [, loc, scale])
gennorm_censored_loglikelihood(data, beta [, loc, scale])
halfgennorm_censored_loglikelihood(data, beta [, loc, scale])
rdist_censored_loglikelihood(data, c [, loc, scale])
foldnorm_censored_loglikelihood(data, c [, loc, scale])
foldcauchy_censored_loglikelihood(data, c [, loc, scale])
genpareto_censored_loglikelihood(data, c [, loc, scale])
loguniform_censored_loglikelihood(data, a, b [, loc, scale])
truncpareto_censored_loglikelihood(data, b, c [, loc, scale])
genhalflogistic_censored_loglikelihood(data, c [, loc, scale])
truncexpon_censored_loglikelihood(data, b [, loc, scale])
johnsonsb_censored_loglikelihood(data, a, b [, loc, scale])
kappa3_censored_loglikelihood(data, a [, loc, scale])
cosine_censored_loglikelihood(data [, loc, scale])
anglit_censored_loglikelihood(data [, loc, scale])
semicircular_censored_loglikelihood(data [, loc, scale])
fatiguelife_censored_loglikelihood(data, c [, loc, scale])
johnsonsu_censored_loglikelihood(data, a, b [, loc, scale])
genexpon_censored_loglikelihood(data, a, b, c [, loc, scale])
exponpow_censored_loglikelihood(data, b [, loc, scale])
exponweib_censored_loglikelihood(data, a, c [, loc, scale])
mielke_censored_loglikelihood(data, k, s [, loc, scale])
powernorm_censored_loglikelihood(data, c [, loc, scale])
powerlognorm_censored_loglikelihood(data, c, s [, loc, scale])
skewcauchy_censored_loglikelihood(data, a [, loc, scale])
normal_censored_loglikelihood(data [, loc, scale])
moyal_censored_loglikelihood(data [, loc, scale])
levy_censored_loglikelihood(data [, loc, scale])
levy_l_censored_loglikelihood(data [, loc, scale])
rayleigh_censored_loglikelihood(data [, loc, scale])
cauchy_censored_loglikelihood(data [, loc, scale])
halfcauchy_censored_loglikelihood(data [, loc, scale])
halfnorm_censored_loglikelihood(data [, loc, scale])
hypsecant_censored_loglikelihood(data [, loc, scale])
halflogistic_censored_loglikelihood(data [, loc, scale])
genlogistic_censored_loglikelihood(data, c [, loc, scale])
gompertz_censored_loglikelihood(data, c [, loc, scale])
invweibull_censored_loglikelihood(data, c [, loc, scale])
laplace_censored_loglikelihood(data [, loc, scale])
laplace_asymmetric_censored_loglikelihood(data, kappa [, loc, scale])
logistic_censored_loglikelihood(data [, loc, scale])
gumbel_r_censored_loglikelihood(data [, loc, scale])
gumbel_l_censored_loglikelihood(data [, loc, scale])
exponential_censored_loglikelihood(data [, loc, scale])
weibull_censored_loglikelihood(data, c [, loc, scale])
weibull_max_censored_loglikelihood(data, c [, loc, scale])
uniform_censored_loglikelihood(data [, loc, scale])
powerlaw_censored_loglikelihood(data, a [, loc, scale])
arcsine_censored_loglikelihood(data [, loc, scale])
invgamma_censored_loglikelihood(data, a [, loc, scale])
betaprime_censored_loglikelihood(data, a, b [, loc, scale])
maxwell_censored_loglikelihood(data [, loc, scale])
dgamma_censored_loglikelihood(data, a [, loc, scale])
dweibull_censored_loglikelihood(data, c [, loc, scale])
loggamma_censored_loglikelihood(data, c [, loc, scale])
bradford_censored_loglikelihood(data, c [, loc, scale])
triang_censored_loglikelihood(data, c [, loc, scale])
trapezoid_censored_loglikelihood(data, c, d [, loc, scale])
fisk_censored_loglikelihood(data, c [, loc, scale])
burr_censored_loglikelihood(data, c, d [, loc, scale])
burr12_censored_loglikelihood(data, c, d [, loc, scale])
lomax_censored_loglikelihood(data, c [, loc, scale])
pareto_censored_loglikelihood(data, b [, loc, scale])
loglaplace_censored_loglikelihood(data, c [, loc, scale])
lognormal_censored_loglikelihood(data, s [, loc, scale])
gibrat_censored_loglikelihood(data [, loc, scale])
gamma_censored_loglikelihood(data, a [, loc, scale])
erlang_censored_loglikelihood(data, a [, loc, scale])
beta_censored_loglikelihood(data, a, b [, loc, scale])
```

The corresponding `*_censored_nnlf` procedures return the negative
log-likelihood. Location defaults to 0 and scale to 1. Weibull/log-Laplace/Lomax/Fisk/generalized-logistic/Gompertz/inverse-Weibull shape `c`,
lognormal shape `s`, gamma/Erlang shape `a`, asymmetric Laplace shape `kappa`,
power-function shape `a`, Burr III/XII shapes `c` and `d`, Pareto shape `b`, and beta shapes `a` and `b` must be supplied.
Bradford shape `c`, triangular mode `c`, trapezoidal edges `c,d`,
double-gamma/inverse-gamma shape `a`, beta-prime shapes `a,b`, and
double-Weibull/log-gamma shape `c` are also required.
Scale and shapes with positive domains must be positive and finite; location
must be finite. Triangular/trapezoidal shapes use the bounded domains below.
Invalid data or parameters return quiet NaN, including for empty data.
An empty valid sample has log-likelihood zero. An impossible observation
contributes negative infinity. No finite penalty is substituted for impossible
observations; this differs from SciPy's fitting penalty policy.

For independent observations, contributions are `logpdf(x)` for exact values,
`logcdf(u)` for left censoring, `logsf(l)` for right censoring, and
`log(F(u)-F(l))` for intervals. No observation weights are supported.

Interval masses use log differences in the better-conditioned tail.
Exponential and Weibull intervals additionally use cumulative-hazard
increments in the log domain, avoiding underflow in tiny interval masses.
For narrow normal intervals, a three-point Gauss-Legendre density integral
avoids subtraction of nearly equal tail logs. It is used when the standardized
width times the larger of 1 and the absolute standardized lower endpoint is
at most 0.01; the sixth-order quadrature error is below binary64 roundoff there.
Logistic interval masses use the exact identity
`F(b)-F(a) = F(b)*S(a)*(1-exp(-(b-a)/scale))`, evaluated with log tails,
`expm1`, and the log interval width. This retains finite log masses when the
standardized width underflows. Logistic standardization avoids intermediate
subtraction or division overflow; log magnitudes beyond binary64 range return
negative infinity. No location or scale fixing is required for logistic fits.
The distribution kernels' accuracy and representability limits still apply.

`censored_loglikelihood(data, params, logpdf, logcdf, logsf [, loginterval])`
is a **SciFort extension** for user-supplied continuous distributions. The
three required callbacks are pure scalar functions of `(x, params(:))`.
They must validate their parameter domains and provide stable log tails.
The optional pure callback `(lower, upper, params(:))` supplies a specialized
log interval probability. Without it, intervals too narrow to distinguish
the callback's tail logs can return negative infinity. Empty samples do not
invoke callbacks, so callback-specific parameter validity is the caller's
responsibility for the generic procedure.

## Bounded fits

The existing `exponnorm_fit`, `skewnorm_fit`, `dpareto_lognorm_fit`, `gausshyper_fit`, `irwinhall_fit`, `vonmises_line_fit`, `argus_fit`, `pearson3_fit`, `kappa4_fit`, `alpha_fit`, `genextreme_fit`, `crystalball_fit`, `gengamma_fit`, `truncnorm_fit`, `truncweibull_min_fit`, `invgauss_fit`, `wald_fit`, `recipinvgauss_fit`, `t_fit`, `f_fit`, `jf_skew_t_fit`, `chi2_fit`, `chi_fit`, `nakagami_fit`, `gennorm_fit`, `halfgennorm_fit`, `rdist_fit`, `foldnorm_fit`, `foldcauchy_fit`, `genpareto_fit`, `loguniform_fit`, `truncpareto_fit`, `genhalflogistic_fit`, `truncexpon_fit`, `johnsonsb_fit`, `kappa3_fit`, `cosine_fit`, `anglit_fit`, `semicircular_fit`, `fatiguelife_fit`, `johnsonsu_fit`, `genexpon_fit`, `exponpow_fit`, `exponweib_fit`, `mielke_fit`, `powernorm_fit`, `powerlognorm_fit`, `skewcauchy_fit`, `invgamma_fit`, `betaprime_fit`, `maxwell_fit`, `dgamma_fit`, `dweibull_fit`, `loggamma_fit`, `bradford_fit`, `triang_fit`, `trapezoid_fit`, `moyal_fit`, `levy_fit`, `levy_l_fit`, `normal_fit`, `exponential_fit`, `weibull_fit`, `weibull_max_fit`,
`genlogistic_fit`, `gompertz_fit`, `invweibull_fit`, `rayleigh_fit`, `uniform_fit`, `powerlaw_fit`, `arcsine_fit`, `fisk_fit`, `burr_fit`, `burr12_fit`, `lomax_fit`, `pareto_fit`, `loglaplace_fit`, `lognormal_fit`, `gibrat_fit`, `gamma_fit`, `erlang_fit`,
`beta_fit`, `cauchy_fit`, `laplace_fit`, `laplace_asymmetric_fit`, `logistic_fit`, `halfcauchy_fit`,
`halflogistic_fit`, `halfnorm_fit`, `hypsecant_fit`, `gumbel_r_fit`, and `gumbel_l_fit` generics accept
`censored_data` in place of the observation array:

```fortran
use scifort_stats, only : dp, censored_data, make_censored_data, &
    fit_result, exponential_fit
type(censored_data) :: data
type(fit_result) :: fit

data = make_censored_data(uncensored=[13.0_dp, 22.0_dp, 17.0_dp, 15.0_dp], &
    right=[20.0_dp, 18.0_dp])
call exponential_fit(data, [0.0_dp, 0.1_dp], [0.0_dp, 100.0_dp], fit)
! fit%params is approximately [0, 26.25].
```

Bounds and parameter order follow the existing API: `[loc, scale]` for normal,
uniform, arcsine, Gibrat, exponential, Rayleigh, Cauchy, symmetric Laplace, logistic, half-Cauchy, half-logistic, half-normal, hyperbolic-secant, and both Gumbel families, `[c, loc, scale]` for both Weibull families, log-Laplace, Lomax, and Fisk, `[a,loc,scale]` for the power function, `[c, d, loc, scale]` for Burr III/XII, `[b, loc, scale]` for Pareto, and `[s, loc, scale]` for
lognormal, `[a, loc, scale]` for gamma and Erlang, `[kappa,loc,scale]` for asymmetric Laplace, and `[a, b, loc, scale]` for beta. Equal bounds fix a parameter.
The optional `guess`, `max_iter`, and `tolerance` arguments retain their
existing meanings. Empty samples are rejected. Results use `fit_result` and
the existing optimizer status codes.

These are bounded local fits using SciFort's coordinate/pattern search.
Success indicates convergence within the specified box, not uniqueness,
identifiability, global optimality, or existence of an unconstrained MLE.
All-censored samples can converge at a box boundary. Fully free location
fits can also be ill-posed; use scientifically justified fixed parameters and
bounds, and inspect results. SciPy and SciFort need not use the same optimizer
or return identical iteration counts.

## Bounded uniform, power-function, and arcsine details

`uniform_censored_loglikelihood(data [, loc, scale])`,
`powerlaw_censored_loglikelihood(data, a [, loc, scale])`, and
`arcsine_censored_loglikelihood(data [, loc, scale])` have `_nnlf` companions.
Location is the finite lower support endpoint (default zero); scale is the
positive finite support width (default one). The rounded upper edge
`loc+scale` must also be finite; an overflowing edge returns quiet NaN, even
for empty data. This matches the existing beta likelihood's support contract.
The two-term upper-edge calculation retains its rounding correction, so a
small width or endpoint offset is not discarded when it rounds away in the
sum. The density support includes both mathematical endpoints.

For `z=(x-loc)/scale` on `[0,1]`, uniform density is `1/scale` and the
power-function density is `a*z**(a-1)/scale`, with finite `a>0`.
Uniform is the `a=1` special case. Arcsine density is
`1/(pi*scale*sqrt(z*(1-z)))`, the beta `(1/2,1/2)` special case.
The power-function CDF is `z**a`; arcsine CDF is `2*asin(sqrt(z))/pi`.
All three have zero density outside the support, zero CDF at/below the lower
edge, and zero survival at/above the upper edge.

At the power-function lower endpoint, density is infinite for `a<1`,
`1/scale` for `a=1`, and zero for `a>1`; upper-endpoint density is `a/scale`.
The scalar `powerlaw_pdf` and `powerlaw_logpdf` now return the correct infinite
limit at the lower endpoint for `a<1`, correcting their previous zero result.
Arcsine density is infinite at either endpoint. These infinities are preserved
in likelihoods, and the existing optimizer rejects nonfinite objectives.
Equal finite interval endpoints still mean exact observations.

Support clipping and normalized distances reuse the beta helpers. The power
function has direct guarded log-CDF and complementary-tail calculations;
subnormal or underflowed normalized upper distances retain their logarithms
so a large shape can compensate without losing the tail. Its
intervals retain `a*log1p((b-a)/(a-loc))` in the log domain instead of
subtracting rounded CDFs. Uniform interior intervals use log width divided by
support width. Nearly full-support intervals complement both excluded tails.
Arcsine reuses the beta density, tails, and narrow-interval quadrature, with a
separate complement branch that retains both small excluded endpoint masses.
These guarantees apply to the new censored kernels, with the documented
finite-upper-edge restriction and binary64 representability limits.

The existing fit generics accept `censored_data`: bounds/results are
`[loc,scale]` for uniform/arcsine and `[a,loc,scale]` for power-function fits.
Every parameter may be free or fixed by equal bounds; the fixed-support
restriction of the general beta fit does not apply to these overloads.
Fits use the existing bounded local optimizer. In particular, a uniform fit
with free location and width can stop on a flat location direction before
reaching the global minimum-width support. Fixing location, supplying an
informed initial guess, or comparing multiple bounded starts can help; success
alone does not prove global optimality. Arcsine and small-shape power-function
fits can approach singular endpoint densities when support parameters are free.

```fortran
use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
    uniform_fit, powerlaw_fit, arcsine_fit
type(censored_data) :: data
type(fit_result) :: fit

data = make_censored_data(uncensored=[0.25_dp,0.5_dp,0.75_dp])
call uniform_fit(data, [0.0_dp,0.75_dp], [0.0_dp,2.0_dp], fit)
call powerlaw_fit(data, [0.1_dp,0.0_dp,1.0_dp], [10.0_dp,0.0_dp,1.0_dp], fit)
call arcsine_fit(data, [0.0_dp,1.0_dp], [0.0_dp,1.0_dp], fit)
```

See the [validation report](../validation/bounded_family_censored_2026-11-16.md).

## Burr family details

`fisk_censored_loglikelihood(data, c [, loc, scale])`,
`burr_censored_loglikelihood(data, c, d [, loc, scale])`, and
`burr12_censored_loglikelihood(data, c, d [, loc, scale])` have `_nnlf`
companions. All shapes and scale must be positive and finite; `loc` is the
finite lower support endpoint. Location defaults to zero and scale to one.
The names follow SciPy: `burr` is Type III and `burr12` is Type XII.

For `z=(x-loc)/scale > 0` and `t=c*log(z)`, Burr XII has survival
`exp(-d*log1p(exp(t)))`; Burr III has CDF `exp(-d*log1p(exp(-t)))`.
Fisk is the `d=1` special case of both. At `c=1`, Burr XII reduces to Lomax;
reciprocation of standardized Burr III gives Burr XII. The likelihoods use
these mathematical identities in log coordinates without forming powers or
reciprocal data. Exact densities include the logarithmic-transform Jacobian.

Below `loc`, density and CDF are zero and survival is one. At `loc`, Fisk and
Burr XII have infinite density for `c<1`, finite density `d/scale` for `c=1`
(with `d=1` for Fisk), and zero density for `c>1`. Burr III uses `c*d` for
this classification, with finite endpoint density `c*d/scale` when `c*d=1`.
The endpoint CDF is zero and survival is one in all three families. Equal
finite interval endpoints remain exact observations. Nonfinite likelihoods
are preserved and rejected as fit objectives by the existing optimizer.

Interval probabilities use increments of `d*softplus(t)` (Type XII) or
`d*softplus(-t)` (Type III). Log widths and stable softplus differences retain
adjacent intervals even when endpoint tail logs coincide. Nearly full-support
intervals complement both excluded tails. Guarded products and log magnitudes
allow compensating shapes even when `c*log(z)` exceeds binary64 range.
Log likelihoods outside binary64 range return negative infinity. These are
new censored-kernel guarantees; scalar routines retain their existing limits.
The scalar Fisk/Burr overflow guards now also handle shapes below one with
floating-point traps enabled, without relying on short-circuit evaluation.

The existing fit generics accept `censored_data`, with `[c,loc,scale]` bounds
and results for Fisk and `[c,d,loc,scale]` for Burr III/XII. Every parameter
may be free or fixed by equal bounds. Fits are bounded local searches; shape,
location and scale can be weakly identified. A successful fit does not prove
an unconstrained MLE exists or is unique. Singular endpoint densities can make
free-location fits ill-posed. Use scientifically justified bounds and inspect
optimizer status; harder four-parameter fits may need more than 1000 sweeps.

```fortran
use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
    fisk_fit, burr_fit, burr12_fit
type(censored_data) :: data
type(fit_result) :: fit

data = make_censored_data(uncensored=[0.5_dp,1.0_dp,2.0_dp])
call fisk_fit(data, [0.1_dp,0.0_dp,1.0_dp], [10.0_dp,0.0_dp,1.0_dp], fit)
call burr_fit(data, [1.0_dp,0.1_dp,0.0_dp,1.0_dp], [1.0_dp,10.0_dp,0.0_dp,1.0_dp], fit)
call burr12_fit(data, [1.0_dp,0.1_dp,0.0_dp,1.0_dp], [1.0_dp,10.0_dp,0.0_dp,1.0_dp], fit)
! With c=1 and fixed location/scale, d=3/sum(log(1+x)) for the reciprocal sample above.
```

See the [validation report](../validation/burr_family_censored_2026-11-16.md).

## Lomax details

`lomax_censored_loglikelihood(data, c [, loc, scale])` and its `_nnlf`
companion use finite shape `c > 0`, finite lower endpoint `loc` (default 0),
and finite `scale > 0` (default 1), matching the
[SciPy Lomax parameterization](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.lomax.html).
For `z=(x-loc)/scale >= 0`, density is `c*(1+z)**(-c-1)/scale`, survival is
`(1+z)**(-c)`, and `log1p(z)` is exponential with scale `1/c`.
The equivalent Pareto location is `loc-scale`; the implementation does not
form that potentially unrepresentable shifted location.

Below `loc`, density and CDF are zero and survival is one. At `loc`, density
is `c/scale`, CDF is zero, and survival is one. Intervals crossing `loc` are
clipped to the support. Cumulative hazards retain their logarithms when
standardized displacements or shape products underflow. Subnormal
displacements also use log distances so a huge shape does not amplify their
relative rounding error. Narrow intervals use
the exact hazard increment `c*log1p((b-a)/(scale+a-loc))`, evaluated through
log distances and positive log sums. Large hazards use direct guarded
multiplication to retain accuracy in nearly unit probabilities. Log magnitudes
beyond binary64 range return negative infinity. These guarantees apply to the
new censored kernels; the existing scalar routines are unchanged.

`lomax_fit` accepts bounds, guesses, and returns parameters in `[c, loc, scale]`
order. Any parameter can be free or fixed by equal bounds. Fits are bounded
local searches; an unconstrained shape/scale fit can approach an exponential
limit rather than have a finite optimum. Fix scientifically justified
parameters or choose meaningful bounds. For example:

```fortran
use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, lomax_fit
type(censored_data) :: data
type(fit_result) :: fit

data = make_censored_data(uncensored=[0.0_dp,1.0_dp,3.0_dp])
call lomax_fit(data, [0.1_dp,0.0_dp,1.0_dp], [10.0_dp,0.0_dp,1.0_dp], fit)
! Fixed location 0 and scale 1: fitted shape is 1/log(2).
```

See the [validation report](../validation/lomax_censored_2026-11-15.md)
for independent references and numerical checks.

## Pareto details

`pareto_censored_loglikelihood(data, b [, loc, scale])` and its `_nnlf`
companion use finite shape `b > 0`, finite location shift `loc` (default 0),
and finite `scale > 0` (default 1), matching the
[SciPy Pareto parameterization](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.pareto.html).
The support is `x >= loc+scale`; `loc` alone is not the lower support endpoint.
With `z=(x-loc)/scale`, density is `b*z**(-b-1)/scale`, survival is `z**(-b)`,
and `log(z)` follows an exponential distribution with scale `1/b`.
Exact log densities include the transformation Jacobian `-log(x-loc)`.

Below the support, density and CDF are zero and survival is one. At the lower
endpoint the density is `b/scale`, CDF is zero, and survival is one.
Intervals crossing the support boundary are clipped to the support. The
implementation avoids forming `loc+scale`, overflowing raw displacements,
or overflowing standardized ratios. Within-support interval masses use
`S(a)*(1-exp(-b*log1p((upper-a)/(a-loc))))`, with the increment retained in
the log domain. This preserves finite logs for narrow intervals and
underflowed probabilities. Log magnitudes outside binary64 range return
negative infinity. Existing scalar kernels are unchanged.

`pareto_fit` accepts bounds, guesses, and returns parameters in `[b, loc, scale]`
order. Any parameter may be free or fixed by equal bounds; no fixed-location
restriction is imposed. As with the other censored fits, impossible data
contribute negative infinity, rather than SciPy's finite fitting penalty.
The fit is a bounded local search; support-changing parameters require
scientifically justified bounds. For example:

```fortran
use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, pareto_fit
type(censored_data) :: data
type(fit_result) :: fit

data = make_censored_data(uncensored=[1.0_dp,2.0_dp,4.0_dp])
call pareto_fit(data, [0.1_dp,0.0_dp,1.0_dp], [10.0_dp,0.0_dp,1.0_dp], fit)
! Fixed location 0 and scale 1: fitted shape is 1/log(2).
```

See the [validation report](../validation/pareto_censored_2026-11-15.md)
for independent references and numerical checks.

## Log-Laplace details

`loglaplace_censored_loglikelihood(data, c [, loc, scale])` and its `_nnlf`
companion use finite `c > 0`, finite lower support endpoint `loc` (default 0),
and finite `scale > 0` (default 1). For `z=(x-loc)/scale > 0`,
`log(z)` follows Laplace with center zero and scale `1/c`. The exact density
includes the Jacobian `1/(x-loc)`; censoring contributions use the transformed
probabilities without a density Jacobian.

Below `loc`, density and CDF are zero and survival is one. At `loc`, density
is infinite for `c < 1`, `1/(2*scale)` for `c = 1`, and zero for `c > 1`.
The likelihood preserves these endpoint limits. Fits whose objective is
nonfinite at every candidate return `fit_status_no_finite_objective`.
Free-location fits can approach an unbounded endpoint density when `c < 1`;
choose bounds and fixed parameters with care.

The implementation reuses the stable Weibull log-distance/scale calculation,
with guarded shape multiplication. Same-side interval masses retain
`c*log1p((b-a)/(a-loc))` in the log domain. Intervals crossing `z=1` add positive
pieces or complement both omitted tails, preserving tiny masses and tiny
near-unit deficits. Log likelihoods beyond binary64 range return negative
infinity. These guarantees apply to the new censored kernels, without changing
the existing scalar kernels.

`loglaplace_fit` accepts bounds, guesses, and returns parameters in
`[c, loc, scale]` order. Any parameter can be free or fixed by equal bounds;
no fixed-location restriction is imposed. For example:

```fortran
use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, loglaplace_fit
type(censored_data) :: data
type(fit_result) :: fit

data = make_censored_data(uncensored=[1.0_dp,2.0_dp,4.0_dp])
call loglaplace_fit(data, [0.1_dp,0.0_dp,2.0_dp], [10.0_dp,0.0_dp,2.0_dp], fit)
! Fixed location 0 and scale 2: fitted shape is 3/(2*log(2)).
```

See the [validation report](../validation/loglaplace_censored_2026-11-15.md)
for independent references and numerical checks.

## Asymmetric Laplace details

`laplace_asymmetric_censored_loglikelihood(data, kappa [, loc, scale])`
and its `_nnlf` companion require finite positive asymmetry `kappa`, finite
cusp location `loc` (default zero), and positive finite `scale` (default one).
Support is the entire real line. With `z=(x-loc)/scale`, density is
`exp(-kappa*z)/(scale*(kappa+1/kappa))` for `z>=0`, and
`exp(z/kappa)/(scale*(kappa+1/kappa))` for `z<0`. The probability to the left
of the cusp is `kappa**2/(1+kappa**2)`. Shape one recovers symmetric Laplace;
reflection exchanges `kappa` with `1/kappa` without changing scale.

The censored kernels evaluate shape weights, normalization, decay rates, and
interval widths in the log domain, avoiding explicit reciprocals, squared
shapes, and overflowing raw displacements. When a complementary probability
would underflow, its log is obtained by adding the positive masses on each
side rather than subtracting a near-unit tail. Intervals entirely on one side
use exponential increments; crossing intervals split at the cusp or complement
both omitted tails, preserving narrow and nearly full-support probabilities.
These robustness measures belong to the new censored API; existing scalar
functions retain their current numerical behavior.

`laplace_asymmetric_fit` accepts bounds and returns parameters in
`[kappa,loc,scale]` order. Any parameter may be free or fixed using equal bounds.
Fits are bounded local searches, without uniqueness or global-optimality
guarantees; use scientifically justified bounds.

```fortran
use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
    laplace_asymmetric_fit
type(censored_data) :: data
type(fit_result) :: fit

data = make_censored_data(uncensored=[-2.0_dp,-1.0_dp,1.0_dp,3.0_dp])
call laplace_asymmetric_fit(data, [2.0_dp,0.0_dp,0.05_dp], [2.0_dp,0.0_dp,10.0_dp], fit)
! Fixed shape 2 and location zero; fitted scale is approximately 2.375.
```

## Hyperbolic-secant details

`hypsecant_censored_loglikelihood(data [, loc, scale])` and its `_nnlf`
companion use finite location `loc` (default zero) and positive finite scale
`scale` (default one). Support is the entire real line. For standardized
`z=(x-loc)/scale`, density is `1/(pi*cosh(z)*scale)` and
`F(z)=2*atan(exp(z))/pi`. The survival function follows by reflection.

The new callbacks use bounded exponentials, log-domain arctangents, and guarded
location/scale standardization. Interval probabilities use the arctangent
subtraction identity with positive numerator and denominator. For narrow
intervals the numerator retains the log of `1-exp(-(b-a)/scale)`, including
when the standardized width underflows. Wide intervals crossing the center
use the complement of both omitted tails, preserving nearly unit probabilities.
These measures apply to the new censored routines; scalar routines retain
their existing numerical behavior.

`hypsecant_fit` accepts bounds and returns parameters in `[loc,scale]` order.
Either parameter can be free or fixed using equal bounds. Fits remain bounded
local searches, without uniqueness or global-optimality guarantees.

```fortran
use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, hypsecant_fit
type(censored_data) :: data
type(fit_result) :: fit

data = make_censored_data(uncensored=[-1.0_dp,0.0_dp,1.0_dp])
call hypsecant_fit(data, [-2.0_dp,1.0_dp], [2.0_dp,1.0_dp], fit)
! Fixed scale one; symmetric data has fitted location approximately zero.
```

## Half-normal details

`halfnorm_censored_loglikelihood(data [, loc, scale])
hypsecant_censored_loglikelihood(data [, loc, scale])` and its `_nnlf`
companion use lower support endpoint `loc` (default zero) and positive normal
scale `scale` (default one). For `z=(x-loc)/scale >= 0`, density is
`sqrt(2/pi)*exp(-z*z/2)/scale`, CDF is `erf(z/sqrt(2))`, and SF is
`erfc(z/sqrt(2))`. At the support endpoint the density is finite and positive,
CDF is zero, and SF is one. Below support the density and CDF are zero.

The likelihood uses direct log tails. Near support it evaluates the normal
interval from `loc` to `x`, retaining its log width even when standardization
underflows. The complementary tail uses `log1p` so tiny excluded probabilities
survive near unit mass. Narrow intervals reuse the normal density quadrature;
wide intervals use stable survival differences. Raw displacement overflow is
avoided by the existing guarded location/scale standardization helper.
These robustness measures belong to the new censored API; the existing scalar
functions retain their current numerical behavior.

`halfnorm_fit` accepts bounds and returns parameters in `[loc,scale]` order.
Either parameter may be free or fixed using equal bounds. Exact observations
at the support endpoint are allowed. Fits remain bounded local searches,
without a guarantee of a global optimum or an unconstrained MLE.

```fortran
use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, halfnorm_fit
type(censored_data) :: data
type(fit_result) :: fit

data = make_censored_data(uncensored=[1.0_dp,2.0_dp,3.0_dp])
call halfnorm_fit(data, [0.0_dp,0.05_dp], [0.0_dp,10.0_dp], fit)
! Fixed location zero; fitted scale is approximately sqrt(14/3).
```

## Gibrat details

`gibrat_censored_loglikelihood(data [, loc, scale])` and its `_nnlf`
companion reuse lognormal with log-shape fixed at 1. Location is the lower
support endpoint (default zero), and scale must be positive (default one).
The existing stable lognormal log tails and interval calculation apply,
including when standardized displacements underflow or raw displacements
overflow. At or below the support endpoint, density and CDF are zero and SF
is one. Invalid inputs return NaN; valid empty data has log-likelihood zero.

`gibrat_fit` accepts bounds in `[loc,scale]` order. The censored overload
requires fixed location, as does the current lognormal overload. Scale may be
free or fixed; shape is always 1 and is absent from the bounds and results.
An exact observation at or below the support endpoint has zero likelihood,
so a fully fixed fit containing one returns `fit_status_no_finite_objective`.
Existing array-based calls retain their original behavior, including support
for free location.

```fortran
use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, gibrat_fit
type(censored_data) :: data
type(fit_result) :: fit

data = make_censored_data(uncensored=[1.0_dp,2.0_dp,4.0_dp])
call gibrat_fit(data, [0.0_dp,0.05_dp], [0.0_dp,10.0_dp], fit)
! Fixed location zero; fitted scale is approximately the geometric mean, 2.
```

## Erlang details

`erlang_censored_loglikelihood(data, a [, loc, scale])` and its `_nnlf`
companion reuse the gamma censored kernels, including stable interval masses
and extreme log tails. Shape `a` must be finite and positive; location defaults
to zero and scale to one. As in the existing Erlang scalar/array API, positive
noninteger shapes are accepted without warnings. The mathematical Erlang
special case has integer shape; fitting searches shape continuously, as SciPy
does. Fix shape with equal bounds when an integer shape is required.

`erlang_fit` accepts bounds in `[a,loc,scale]` order and requires fixed location.
It excludes exact observations at the support endpoint, including equal-endpoint
intervals, matching the current gamma fit restriction. Likelihood evaluation
itself retains the gamma endpoint density contract: positive infinity for
`a<1`, `1/scale` for `a=1`, and zero for `a>1` at `x=loc`.

```fortran
use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, erlang_fit
type(censored_data) :: data
type(fit_result) :: fit

data = make_censored_data(uncensored=[1.0_dp,2.0_dp,3.0_dp],right=[4.0_dp])
call erlang_fit(data, [2.0_dp,0.0_dp,0.05_dp], [2.0_dp,0.0_dp,10.0_dp], fit)
! Fixed integer shape 2, fixed location 0, free scale.
```

## Rayleigh and Weibull maximum details

`rayleigh_censored_loglikelihood(data [, loc, scale])` uses finite lower
support endpoint `loc` and positive finite `scale`, defaulting to zero and
one. `weibull_max_censored_loglikelihood(data, c [, loc, scale])` uses positive
finite shape `c`, finite upper support endpoint `loc`, and positive finite
`scale`, with the same defaults. Corresponding `*_censored_nnlf` procedures
return negative log-likelihoods. Invalid data or parameters return NaN,
including for empty data; a valid empty sample returns zero.

Rayleigh has cumulative hazard `H(x)=((x-loc)/scale)**2/2` above support.
Its SF is `exp(-H)` and its CDF is `1-exp(-H)`. The callbacks form log
hazards from log distances and log scales, avoiding overflowing subtraction,
squaring, or conversion to Weibull scale `sqrt(2)*scale`. Its density is zero
at and below `loc`. The CDF is zero and SF is one at/below that endpoint.
For finite intervals, the implementation clips to support and evaluates
`-H(lower) + log(1-exp(-delta_H))`, with the positive increment obtained from

```text
delta_H = (upper-lower)*((upper-loc)+(lower-loc))/(2*scale**2)
```

using log factors rather than subtracting squared endpoints. Tiny increments
use `expm1`; large increments use a complementary `log1p` expression, which
retains tiny excluded probabilities for nearly full-support intervals.

For Weibull maximum, `H(x)=((loc-x)/scale)**c` below support and
`logcdf=-H`. Its SF uses the stable log hazard mass above. Intervals are
reflections of minimum-Weibull intervals: endpoints reverse, and left and
right censoring interchange. Narrow intervals reuse the existing stable
minimum-Weibull hazard-increment kernel. Wide intervals use a direct log CDF
difference to retain near-one log masses. Density evaluation also uses log
distance and log hazard, extending the range beyond explicitly forming the
standardized displacement. At `loc`, density is positive infinity for `c<1`,
`1/scale` for `c=1`, and zero for `c>1`, matching existing Weibull behavior.
The CDF is one and SF zero at/above `loc`. Equal finite interval endpoints
still mean exact observations, including endpoint density behavior.

```fortran
use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
    rayleigh_fit, weibull_max_fit
type(censored_data) :: data
type(fit_result) :: fit

data = make_censored_data(uncensored=[0.5_dp,1.2_dp,2.5_dp], right=[3.5_dp])
call rayleigh_fit(data, [0.0_dp,0.05_dp], [0.0_dp,10.0_dp], fit)
data = make_censored_data(uncensored=[-0.5_dp,-1.2_dp,-2.5_dp], left=[-3.5_dp])
call weibull_max_fit(data, [0.1_dp,0.0_dp,0.05_dp], [6.0_dp,0.0_dp,10.0_dp], fit)
! These examples fix loc=0. Bounds may free either location.
```

Bounds follow `[loc,scale]` for Rayleigh and `[c,loc,scale]` for Weibull
maximum. All parameters can be free or fixed; `guess`, `max_iter`, and
`tolerance` retain the existing meanings. The bounded local optimizer is
unchanged, as are the scalar and array APIs. Infinite endpoint likelihoods
are not replaced with a finite likelihood penalty. The optimizer rejects
nonfinite objective candidates; fully fixed nonfinite objectives return
`fit_status_no_finite_objective`. A finite local fit does not establish a
global finite MLE, especially when free support locations and shapes can
produce singular endpoint likelihoods. Sample sums must fit in binary64;
log magnitudes beyond its range return negative infinity where appropriate.
Existing minimum-Weibull interval accuracy limits apply to the reflected
narrow-interval kernel.

## Half-Cauchy and half-logistic details

`halfcauchy_censored_loglikelihood` and `halflogistic_censored_loglikelihood`
accept `(data [, loc, scale])`; the corresponding `*_censored_nnlf` procedures
return the negative log-likelihood. Location is the finite lower support
endpoint, default zero, and scale must be positive and finite, default one.
The support includes `loc`. Their `halfcauchy_fit` and `halflogistic_fit`
censored overloads use `[loc, scale]` bounds with the existing `guess`,
`max_iter`, and `tolerance` options. Either parameter can be free or fixed;
no fixed-location restriction is imposed.

For `z=(x-loc)/scale >= 0`, these are normalized positive halves of Cauchy
and logistic distributions. Their densities and upper tails are twice the
parent density and upper tail. The lower-tail identities are

```text
half-Cauchy:   F(z) = 2*atan(z)/pi
half-logistic: F(z) = tanh(z/2) = (1-exp(-z))/(1+exp(-z))
```

Below support, the density and CDF are zero and the SF is one. At `loc`, the
CDF is zero, the SF is one, and densities are respectively `2/(pi*scale)`
and `1/(2*scale)`. Exact endpoint observations are valid. As elsewhere in
this API, equal finite interval endpoints mean an exact observation.
Left censoring at or below `loc` and intervals wholly below support have
zero probability and contribute negative infinity, rather than an invalid
input status. Intervals crossing `loc` are clipped to support. Infinite
endpoints follow the usual one-sided censoring rules.

Narrow intervals wholly above support reuse the stable parent interval kernels
with log normalization `log(2)`. Wide intervals use
`logsf(lower) + log1p(-exp(logsf(upper)-logsf(lower)))` when the upper
survival probability is at most half the lower survival probability. This
retains tiny excluded probabilities for intervals covering almost the whole
support. Clipped intervals use the half-family CDF.
Near the boundary, direct arctangent or exponential-increment CDF formulas
retain finite log probabilities when the standardized distance underflows.
Near-one SF values use `log1p` of the small CDF. In the opposite tail, the
implementation directly evaluates the log SF and uses `log1p` for the CDF.
This avoids subtracting a parent CDF near one half or losing tiny logarithms
of probabilities near one. The half-logistic endpoint log density is evaluated
directly to avoid cancellation in the normalization.

The half-Cauchy callbacks retain finite density and SF logs beyond the
representable standardized range. Half-logistic log densities and SFs whose
magnitudes exceed binary64 return negative infinity. Existing scalar and
array APIs retain their behavior. Sample sums must still fit in binary64.
The bounded local optimizer makes no global-optimality or identifiability
guarantee. Endpoint data can favor the lower scale bound; all-censored data
can have flat regions or boundary solutions. Inspect diagnostics and bounds.

```fortran
use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
    halfcauchy_fit, halflogistic_fit
type(censored_data) :: data
type(fit_result) :: fit

data = make_censored_data(uncensored=[0.5_dp, 1.2_dp, 2.5_dp], right=[3.5_dp])
call halfcauchy_fit(data, [0.0_dp, 0.05_dp], [0.0_dp, 10.0_dp], fit)
call halflogistic_fit(data, [0.0_dp, 0.05_dp], [0.0_dp, 10.0_dp], fit)
! Both calls fix loc=0; fit%params = [loc, scale].
```

## Cauchy details

`cauchy_censored_loglikelihood(data [, loc, scale])` uses finite median `loc`
and positive finite half width at half maximum `scale`. `cauchy_fit` accepts
censored data with bounds ordered `[loc, scale]`; either parameter may be
free or fixed with equal bounds. Optional `guess`, `max_iter`, and `tolerance`
retain their existing meanings. The support is the whole real line.

The density is `1/(pi*scale*(1+z*z))`, where `z=(x-loc)/scale`.
Censored callbacks compute the log absolute displacement before scaling.
This retains finite log densities and small-tail log probabilities even
when forming `z` or its square would overflow, or when the ordinary tail
probability would underflow. The small tail uses the reciprocal-arctangent
identity. The near-one log tail uses `log1p` directly.

For finite intervals on one side of the median, the positive angle between
the endpoint arctangents satisfies

```text
tan(angle) = ((upper-lower)/scale) / (1 + z_lower*z_upper)
log_mass = log(angle) - log(pi)
```

The calculation uses log widths and log products, without forming endpoint
products or subtracting CDFs. Narrow intervals crossing the median add two
positive arctangent contributions. Crossing intervals whose endpoints are
at least one scale from the median instead use `log1p` of the negative sum
of the two missing tails, preserving near-one log probabilities.

`log(atan(exp(t)))` is evaluated without overflowing or underflowing its
exponential argument. When `t < log(epsilon)/2`, it returns `t`: the relative
error of replacing `atan(u)` by `u` is at most `u*u/3`, below binary64
roundoff there. For positive `t`, it uses `atan(u)=pi/2-atan(1/u)`.

Invalid data or parameters return NaN; an empty valid sample returns zero.
Equal interval endpoints represent exact observations; one infinite endpoint
represents one-sided censoring. A sum must still fit in binary64. Existing
array likelihood/fit calls and scalar Cauchy routines are unchanged; their
original standardization limits still apply. The censored overload extends
the numerical range of its private callbacks, without changing those APIs.

Fits use the existing bounded local search. Cauchy likelihoods can have
multiple stationary points; convergence does not establish a global optimum,
identifiability, or a finite unconstrained MLE. Check scientifically justified
bounds and consider multiple starting guesses.

```fortran
use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, cauchy_fit
type(censored_data) :: data
type(fit_result) :: fit

data = make_censored_data(uncensored=[-1.0_dp, 0.0_dp, 2.0_dp], right=[3.0_dp])
call cauchy_fit(data, [-5.0_dp, 0.05_dp], [8.0_dp, 10.0_dp], fit)
! fit%params = [loc, scale]; inspect fit%success and fit%nllf.
```

## Laplace details

`laplace_censored_loglikelihood(data [, loc, scale])
laplace_asymmetric_censored_loglikelihood(data, kappa [, loc, scale])` uses center `loc` and
positive finite mean absolute deviation `scale`, matching SciPy's Laplace
parameterization. `laplace_fit` accepts censored data with bounds ordered
`[loc, scale]`. Neither parameter must be fixed. Equal bounds fix a parameter;
`guess`, `max_iter`, and `tolerance` retain their existing meanings.

For finite intervals wholly on one side of the center, the mass is the
nearer endpoint's exponential tail multiplied by
`1-exp(-(upper-lower)/scale)`. Both terms are evaluated in the log domain.
Intervals crossing the center use the average of the positive exponential
masses on the two sides, evaluated with a stable log average. This avoids subtracting nearly equal
CDFs, preserves adjacent and subnormal-width intervals, and retains finite
log probabilities when the standardized width underflows. The shared
standardization and log-distance helpers avoid intermediate displacement
and width overflow.

Invalid data or parameters return NaN; an empty valid sample returns zero.
Unrepresentable log densities or log tails return negative infinity. The
existing binary64 likelihood accumulation limits still apply. Bounded fits
use the existing local optimizer: convergence does not establish a global
optimum or a finite unconstrained MLE. In particular, exact observations at
the center can drive scale to its lower bound, and one-sided censored samples
can drive location to a box boundary. Existing array-based Laplace likelihood
and fit procedures retain their behavior.

```fortran
use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, laplace_fit
type(censored_data) :: data
type(fit_result) :: fit

data = make_censored_data(uncensored=[-1.0_dp, 0.0_dp, 2.0_dp], right=[3.0_dp])
call laplace_fit(data, [-5.0_dp, 0.05_dp], [8.0_dp, 10.0_dp], fit)
! fit%params = [loc, scale]; inspect fit%success and fit%nllf.
```

## Gumbel details

`gumbel_r` and `gumbel_l` follow SciPy's respective right- and left-skewed
parameterizations, with finite location and positive finite scale. Both
parameters may be free in censored fits; equal bounds fix either parameter.
For the right-skewed family, `F(x)=exp(-exp(-z))`, where
`z=(x-loc)/scale`. For a finite interval with standardized endpoints `za < zb`,

```text
F(b)-F(a) = F(b) * (1-exp(-delta))
delta = exp(-za) * (1-exp(-(b-a)/scale))
```

The implementation evaluates the increment in the log domain using `expm1`
and the log interval width. This preserves narrow and extreme-tail masses
without subtracting nearly equal tail probabilities. Left-skewed likelihoods
share the same kernels through reflection, reversing interval endpoints and
interchanging the two tails. The standardization helper shared with logistic
censoring checks binary exponents to avoid intermediate overflow.
Unrepresentable log densities or log tails return negative infinity; finite
sample sums must still fit in binary64. Fits use the same bounded local
optimizer and have the convergence limitations described above.

## Lognormal details

Lognormal uses SciPy's parameterization: `log((x-loc)/scale)` is normal with
mean zero and standard deviation `s`. Exact observations contribute the
normal log density minus `log(x-loc)` (the Jacobian); censored observations
contribute probabilities without a Jacobian. Support is strictly `x > loc`.

The censored `lognormal_fit` overload currently requires fixed location:
set `lower(2) == upper(2)`. Unequal location bounds return
`fit_status_invalid_input`; the existing uncensored array overload is unchanged.
For example:

```fortran
call lognormal_fit(data, [0.1_dp, 0.0_dp, 0.1_dp], &
    [5.0_dp, 0.0_dp, 100.0_dp], fit)
! fit%params = [s, loc, scale], with loc fixed at zero.
```

Likelihoods accept any finite location. They compute log distances before
normalization, avoiding overflow/underflow in `(x-loc)/scale`. Narrow intervals
retain their log-transformed width using `log1p((b-a)/(a-loc))`, then reuse
the normal quadrature above. Tail subtraction handles wider intervals.
This avoids losing adjacent endpoints when their separately computed
logarithms round to the same value. The normal kernels' extreme standardized
tail limits still apply; bounded convergence does not establish identifiability.

## Gamma details

`gamma_censored_loglikelihood(data, a [, loc, scale])
erlang_censored_loglikelihood(data, a [, loc, scale])` uses positive shape
`a` and inverse-rate scale, matching the existing gamma API. At `a=1`, it
calls the exponential censored likelihood, preserving that identity exactly.
An exact observation at `loc` contributes positive infinity for `a<1`,
`-log(scale)` for `a=1`, and negative infinity for `a>1`. A sum containing
both positive and negative infinite contributions is undefined and returns
quiet NaN, without an invalid floating-point operation.

```fortran
call gamma_fit(data, [0.1_dp, 0.0_dp, 0.1_dp], &
    [5.0_dp, 0.0_dp, 100.0_dp], fit)
! fit%params = [a, loc, scale], with loc fixed at zero.
```

The censored `gamma_fit` overload requires equal location bounds and excludes
exact observations at the fixed support endpoint, including equal-endpoint
interval rows. Violations return `fit_status_invalid_input`. Endpoint densities
can make the free-shape likelihood infinite; the likelihood API retains these
mathematical values, while this initial fitting interface excludes such samples.
The existing uncensored array fitting overload is unchanged.

Intervals touching or crossing `loc` use the incomplete-gamma lower tail.
Positive narrow intervals integrate the density with three-point
Gauss-Legendre quadrature after factoring out the lower-endpoint density.
With `q=(upper-lower)/(lower-loc)` and `d=(upper-lower)/scale`, the quadrature
is used when `q*max(1,abs(a-1)) <= 0.001` and `d <= 0.001`. These bounds
limit relative distance and log-density variation. The log-density ratio at
fraction `t` of the interval is `(a-1)*log1p(t*q)-t*d`. Other intervals
subtract log tails in the better-conditioned tail.

If the standardized argument falls below the smallest normal binary64 value,
logarithmic distance and the leading incomplete-gamma series term preserve
finite log probabilities. The neglected relative correction is bounded by
the standardized argument. Ordinary cases reuse the existing gamma kernels;
their numerical and very-large-shape limitations still apply.

## Beta details

Beta uses positive shapes `a` and `b`, location `loc`, and positive support
width `scale`. The support runs from `loc` to the mathematical sum
`loc+scale`. The rounded upper edge must be finite for this initial interface;
parameters outside that arithmetic range return NaN (invalid input for fitting).
The small rounding residual of the upper-edge sum is retained with FastTwoSum.
Lower and upper distances are normalized independently, avoiding a rounded
`1-z` that loses an upper-tail probability. Consequently, extreme-tail results
can be more accurate than the existing scalar beta API, which forms `1-z`.
Ordinary results agree within numerical tolerances.

```fortran
call beta_fit(data, [0.1_dp, 0.1_dp, 0.0_dp, 1.0_dp], &
    [5.0_dp, 5.0_dp, 0.0_dp, 1.0_dp], fit)
! fit%params = [a, b, loc, scale], with support fixed to [0,1].
```

The censored overload requires equal lower/upper bounds for both location and
scale. Exact support-endpoint observations, including equal-endpoint interval
rows, are excluded from fitting; violations return `fit_status_invalid_input`.
The existing uncensored array overload is unchanged. This initial fit is a
bounded local optimization of the two shapes, with no global-optimum or
identifiability guarantee.

Likelihoods retain mathematical endpoint densities. At the lower edge the
log-density is positive infinity for `a<1`, `log(b)-log(scale)` for `a=1`, and
negative infinity for `a>1`; the upper edge exchanges `a` and `b`. Opposing
infinite contributions return NaN. With `a=b=1`, the likelihood agrees with
the uniform distribution, including censoring probabilities.

Intervals crossing a support edge use a single log tail. For a narrow interval
strictly inside support, let `ql=width/lower_distance` and
`qr=width/upper_distance`, both distances evaluated at its lower endpoint.
Three-point Gauss-Legendre integration of the factored density is used when
`ql*max(1,abs(a-1)) <= 0.001` and `qr*max(1,abs(b-1)) <= 0.001`. Its
log-density ratio at fraction `t` is
`(a-1)*log1p(t*ql)+(b-1)*log1p(-t*qr)`. Other intervals use the
better-conditioned difference of log tails. Edge-touching intervals avoid
quadrature over singular densities.

When a normalized endpoint distance rounds to zero, the leading incomplete-beta
term is evaluated from its finite logarithm. It is used only when the far
shape times that coordinate is at most `epsilon*exp(-2)`. Otherwise this
underflow regime returns NaN rather than an uncontrolled approximation.
Gamma-difference and scaled gamma-ratio helpers retain the normalization and
small complementary probability when one shape is tiny. Other shape/range
limitations of the existing incomplete-beta and log-beta kernels still apply.

## Analytic right-censored exponential fit

```fortran
call exponential_fit_right_censored(data, loc=0.0_dp, result=fit)
```

This **SciFort extension** computes the fixed-location exponential MLE without
optimization. With `d` exact observations, its scale is the total exposure
divided by `d`: each exact value contributes `x-loc`, and each right-censoring
limit contributes `max(limit-loc,0)`. The routine accepts exact/right groups
and equivalent exact or right-unbounded interval rows. Left censoring and
finite-width intervals are invalid for this analytic routine.

The result has `params=[loc,scale]`, zero optimizer iterations, and one objective
evaluation on success. Exposure is normalized before summation to avoid
overflow or premature underflow when the MLE is representable.

| Status | Meaning |
| --- | --- |
| `fit_status_success` (0) | Positive finite scale and finite objective. |
| `fit_status_invalid_input` (2) | Invalid data/location, empty sample, or unsupported censoring. |
| `fit_status_no_finite_objective` (3) | An exact observation is below support, or the fitted objective is nonfinite. |
| `fit_status_no_finite_mle` (4) | No exact observations, zero exposure, or scale outside binary64's positive finite range. |

Status 4 leaves parameters and objective as NaN. For example, all right-censored
data with positive exposure have their likelihood supremum at infinite scale;
exact observations all at the fixed location instead give an unbounded
likelihood as scale tends to zero. Neither is reported as a successful fit.

## Verification

`test/test_censored.f90` and `test/test_censored_gumbel.f90` cover the new likelihoods and fits, compatibility with
uncensored calls, invalid inputs, infinite bounds, support endpoints, narrow
and extreme-tail intervals, analytic MLEs, and unsupported ECDF/logrank groups.
Regime-specific tolerances are used. Numerical constants can be reproduced by:

```text
python tools/generate_censored_reference.py
fpm test test_censored test_censored_gumbel test_fit test_likelihood test_logrank test_empirical_distribution
```

The reference generator uses mpmath for independent high-precision defining
formulas and SciPy for fit comparisons. The initial constants use mpmath 1.3.0
at 100 to 700 decimal digits and SciPy 1.15.3. See
[provenance](../CODE_PROVENANCE.md) for implementation sources and attribution.

After an FPM build, randomized comparisons against 100-digit formulas can
also be run with `validation/check_censored.py`. Supply `--module-dir` for
the directory containing `scifort_kinds.mod` and `--library` for the matching
`scifort/libscifort.a`. The script builds its input-driven reference driver
outside `test/`, in checked and optimized configurations. `--seed` and
`--cases` control the reproducible sampling; `--cases` is the number per
distribution, followed by one large-shape Weibull, three lognormal, nine gamma,
nine beta, seven logistic, and sixteen Gumbel stress cases (eight per direction).
Reference tail selection is based on conditioning,
and the checker rejects nonfinite expected values.

See the [dated validation record](../validation/censored_fitting_2026-11-12.md)
for measured results, compiler flags, and reproducible commands.

See also the [lognormal validation record](../validation/lognormal_censored_2026-11-12.md)
for the extension checks and commands.

See the [gamma validation record](../validation/gamma_censored_2026-11-12.md)
for the gamma extension checks and commands.

See the [beta validation record](../validation/beta_censored_2026-11-12.md)
for the beta extension checks and commands.

See the [logistic validation record](../validation/logistic_censored_2026-11-12.md)
for the logistic extension checks and commands.

See the [Gumbel validation record](../validation/gumbel_censored_2026-11-12.md)
for the two Gumbel families, reference checks, and reproducible commands.

## Generalized logistic, Gompertz, and inverse Weibull details

All three families use finite shape `c>0`, finite `loc` (default zero), and
finite `scale>0` (default one). Generalized logistic has support on the whole
real line; Gompertz and inverse Weibull have lower support endpoint `loc`.
For standardized `z=(x-loc)/scale`, their direct exponential tails are
`F_genlogistic=exp(-c*log1p(exp(-z)))`,
`S_gompertz=exp(-c*expm1(z))`, and `F_invweibull=exp(-z**(-c))`.
Their complements and interval increments retain log hazards and original
interval widths, avoiding subtraction of rounded tail probabilities.
Intervals outside the positive-support families have zero probability;
crossing intervals are clipped at `loc`. At that endpoint Gompertz density
is `c/scale`, while inverse Weibull density is zero. Equal interval endpoints
still represent exact measurements.

Each `*_fit` accepts `censored_data` and bounds in `(c,loc,scale)` order;
equal bounds fix any parameter. All parameters can be free, with no additional
location restriction. Existing array-specific fit procedures remain available.
The existing local optimizer may encounter boundary optima or ill-posed MLEs;
success does not establish global optimality. Nonfinite objectives are rejected.

```fortran
use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
    genlogistic_fit, gompertz_fit, invweibull_fit
type(censored_data) :: data
type(fit_result) :: fit

data = make_censored_data(uncensored=[0.5_dp,1.2_dp,2.5_dp], &
    left=[0.3_dp],right=[3.5_dp])
call genlogistic_fit(data,[0.1_dp,0.0_dp,1.3_dp],[10.0_dp,0.0_dp,1.3_dp],fit)
call gompertz_fit(data,[0.1_dp,0.0_dp,1.3_dp],[10.0_dp,0.0_dp,1.3_dp],fit)
call invweibull_fit(data,[0.1_dp,0.0_dp,1.3_dp],[10.0_dp,0.0_dp,1.3_dp],fit)
```

Generalized logistic shape one agrees with logistic; reciprocal standardized
inverse Weibull observations follow minimum Weibull, with the density Jacobian
included for exact data. The new likelihoods have additional extreme-range
guards; this does not change other scalar kernels' representability limits.
The scalar inverse Weibull small-shape multiplication guard is corrected.
See the [validation report](../validation/hazard_family_censored_2026-11-16.md).

## Moyal, Levy, and reflected Levy details

All three families use finite `loc` (default zero) and finite `scale>0`
(default one). Moyal has support on the whole real line. Levy has lower
support endpoint `loc`; reflected Levy (`levy_l`) has upper endpoint `loc`.
Both Levy densities are zero at their endpoint. Intervals are clipped to
support, and equal finite endpoints continue to represent exact data.

For standardized `z`, the transformations `t=exp(-z/2)` for Moyal and
`t=1/sqrt(z)` for Levy give a half-normal variable. Their CDF is
`erfc(t/sqrt(2))` and survival probability is `erf(t/sqrt(2))`. Reflected
Levy swaps these tails after reflection. Log transformed coordinates retain
tiny survival probabilities, and guarded normal exponents retain finite
log densities/CDFs when an unscaled exponential would overflow. Moyal also
retains a finite half coordinate when the full standardized coordinate
overflows. True log magnitudes outside binary64 range return negative infinity.

Intervals retain their transformed width in logarithmic form before using
the existing narrow normal-interval quadrature. Wider intervals use the
better-conditioned tail difference; nearly full intervals complement both
omitted tails. The existing normal-kernel accuracy limits still apply.

Each `*_fit` overload accepts `censored_data` and bounds in `(loc,scale)`
order. Either or both parameters can be free; equal bounds fix a parameter.
Existing array-specific procedures remain available. Fits use the existing
bounded local optimizer, so success does not establish a global optimum.

```fortran
use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
    moyal_fit, levy_fit, levy_l_fit
type(censored_data) :: data
type(fit_result) :: fit

data = make_censored_data(uncensored=[0.5_dp,1.2_dp,2.5_dp],right=[3.5_dp])
call moyal_fit(data,[0.0_dp,0.1_dp],[0.0_dp,10.0_dp],fit)
call levy_fit(data,[0.0_dp,0.1_dp],[0.0_dp,10.0_dp],fit)
data = make_censored_data(uncensored=[-0.5_dp,-1.2_dp,-2.5_dp],left=[-3.5_dp])
call levy_l_fit(data,[0.0_dp,0.1_dp],[0.0_dp,10.0_dp],fit)
```

The new likelihoods have additional extreme-range guards; other scalar
kernels' behavior is unchanged. See the
[validation report](../validation/normal_transform_censored_2026-11-16.md).

## Bradford, triangular, and trapezoidal details

All three families have support `[loc,loc+scale]`, with finite `loc` and
finite `scale>0` (defaults zero and one). These interfaces require a finite
rounded upper support edge, like the existing beta/uniform interfaces;
invalid parameters return quiet NaN even for empty data. Compensated endpoint
distances retain small offsets in `loc+scale` and log normalized distances
when ordinary division underflows.

Bradford requires finite `c>0`. Its standardized density is
`c/((1+c*z)*log1p(c))`, including finite endpoint densities. Interval mass
is `log1p(c*(z_b-z_a)/(1+c*z_a))/log1p(c)`, evaluated in logarithmic form.
For `c<epsilon(1.0_dp)*exp(-2)`, the uniform limit is used: the relative
correction is below binary64 roundoff. The public parameter domain still
excludes `c=0`.

Triangular mode `c` may range from zero to one, including endpoint modes.
Trapezoidal plateau edges satisfy `0<=c<=d<=1`. With
`h=2/(1+d-c)`, its standardized density rises as `h*z/c`, is constant `h`
on the plateau, and falls as `h*(1-z)/(1-d)`. A zero-length ramp is omitted.
The endpoint density is zero unless that endpoint belongs to the plateau.
Triangular is exactly trapezoidal with `d=c`; `(c,d)=(0,1)` is uniform.
Triangular endpoint modes are beta `(1,2)` and `(2,1)`.

Intervals within one ramp integrate its linear density exactly. Intervals
crossing corners integrate the constant plateau height minus triangular
deficits, preserving the original interval width. Separate lower/upper log
coordinates retain subnormal breakpoints and upper distances. Nearly full
intervals complement both omitted tails. No quadrature or observation weights
are added. Equal finite interval endpoints retain the exact-data convention.

The new `*_fit` overloads accept `censored_data`. Parameter bounds have order
`(c,loc,scale)` for Bradford/triangular and `(c,d,loc,scale)` for trapezoidal.
All parameters may be free; equal bounds fix any parameter. The optimizer
rejects candidates violating shape ordering or support. Existing array-specific
procedures remain available. Fits retain the existing local optimizer's
limits, including boundary/flat optima; success does not prove a global MLE.

```fortran
use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
    bradford_fit, triang_fit, trapezoid_fit
type(censored_data) :: data
type(fit_result) :: fit

data = make_censored_data(uncensored=[0.1_dp,0.3_dp,0.6_dp],left=[0.05_dp],right=[0.8_dp])
call bradford_fit(data,[0.01_dp,0.0_dp,1.0_dp],[10.0_dp,0.0_dp,1.0_dp],fit)
call triang_fit(data,[0.0_dp,0.0_dp,1.0_dp],[1.0_dp,0.0_dp,1.0_dp],fit)
call trapezoid_fit(data,[0.0_dp,0.7_dp,0.0_dp,1.0_dp], &
    [0.7_dp,0.7_dp,0.0_dp,1.0_dp],fit)
```

Other scalar kernels' extreme-range behavior is unchanged. See the
[validation report](../validation/piecewise_censored_2026-11-16.md).

## Double-gamma, double-Weibull, and log-gamma details

All three distributions have real-line support, finite positive shape and
scale, and finite location. Shapes follow SciPy: `a` for `dgamma`, and `c`
for `dweibull` and `loggamma`. With `z=(x-loc)/scale`, double-gamma and
double-Weibull put half of a gamma or Weibull density on each side of zero.
Their CDF at location is exactly one half. Exact observations at location
have infinite density for shape below one, zero density above one, and
`1/(2*scale)` at shape one. Both shape-one families equal Laplace.

Same-side intervals reuse the underlying gamma/Weibull interval kernels.
Intervals crossing location sum two positive radial CDF contributions;
nearly full intervals instead retain both omitted tails. Reflection uses
negated observations/location, avoiding an overflowing `2*loc-x`.
Log distances retain normalized radii below the smallest positive real.

Log-gamma means `exp(Z)` is unit-scale gamma with shape `c`; it differs
from lognormal and from the logarithm of the gamma function. Its density is
`exp(c*z-exp(z))/(Gamma(c)*scale)`. Tails use incomplete gamma functions
and log-domain small-argument expansions. Very negative standardized
coordinates retain `c*z` when small shape rescues an overflowing coordinate.
The large-shape density uses a Stirling saddle form around `z=log(c)`.
Narrow intervals integrate relative densities with three-point Gauss-Legendre
quadrature when both `c*delta` and `exp(z)*delta` are small; broader intervals
use the better tail difference. Beyond representable `exp(z)`, an endpoint
expansion retains representable log survival probabilities.

`dgamma_fit`, `dweibull_fit`, and `loggamma_fit` accept `censored_data` with
bounds in `(shape,loc,scale)` order. Any parameter may be free or fixed by
equal bounds. Exact center singularities are retained in likelihoods;
a fully fixed singular or zero-density objective returns
`fit_status_no_finite_objective`. Free-location fits for shape below one can
have unbounded likelihoods as location approaches an exact observation.
A finite bounded optimizer result is not a guarantee of a global MLE.

```fortran
use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
    dgamma_fit, dweibull_fit, loggamma_fit
type(censored_data) :: data
type(fit_result) :: fitted

data = make_censored_data(uncensored=[-1.2_dp,0.5_dp,2.5_dp], &
    left=[-2.0_dp],right=[3.5_dp],interval=reshape([-0.2_dp,0.3_dp],[1,2]))
call dgamma_fit(data,[1.1_dp,0.0_dp,0.2_dp],[5.0_dp,0.0_dp,5.0_dp],fitted)
call dweibull_fit(data,[1.1_dp,0.0_dp,0.2_dp],[5.0_dp,0.0_dp,5.0_dp],fitted)
call loggamma_fit(data,[0.1_dp,0.0_dp,0.2_dp],[5.0_dp,0.0_dp,5.0_dp],fitted)
```

The incomplete-gamma kernel's binary64 accuracy limits carry through to
double-gamma and log-gamma. At extremely large shapes, rounding in the
coordinate and `log(shape)` can limit accuracy near the sharply concentrated
mode; subnormal results retain binary64 resolution limits. Existing scalar
APIs and original array fit procedures keep their signatures. C ABI and
wrapper additions remain deferred. See the
[validation report](../validation/reflected_censored_2026-11-16.md).

## Inverse-gamma, beta-prime, and Maxwell details

These distributions have standardized support `z=(x-loc)/scale >= 0`.
Location is a finite lower endpoint and scale is finite and positive.
Inverse-gamma requires finite positive shape `a`; beta-prime requires finite
positive shapes `a,b`. The current beta-prime censored interface additionally
requires a finite rounded `a+b`, matching the incomplete-beta kernel's
parameter range; otherwise it returns quiet NaN, including for empty data.
Maxwell has no shape parameter and equals chi with three degrees of freedom.

Inverse-gamma uses gamma shape `a` at `y=1/z`, with CDF `Q(a,y)` and survival
`P(a,y)`. Maxwell uses gamma shape `3/2` at `y=z*z/2`. Transformed arguments
are retained as logarithms so a reciprocal, square, or normalized distance
can underflow or overflow without losing a representable log probability.
Exact densities include the analytic Jacobians. At the lower endpoint both
inverse-gamma and Maxwell have zero density. Inverse-gamma with shape one
equals inverse-Weibull with shape one.

Beta-prime maps to beta through `u=z/(1+z)`, `v=1/(1+z)`. Both complementary
coordinates and their logarithms are retained independently. Endpoint-tail
expansions retain probabilities when a transformed coordinate rounds to
zero. At location, density is infinite for `a<1`, `b/scale` for `a=1`, and
zero for `a>1`. At `a=1` beta-prime equals Lomax with shape `b`.

Narrow intervals integrate relative densities with three-point Gauss-Legendre
quadrature after bounding log-density variation. Interval widths come from
the original observations, preserving adjacent and subnormal widths. Broader
intervals use the better tail difference; nearly full intervals complement
both omitted tails. Intervals below support have zero mass, and intervals
crossing location are clipped to support. Equal finite endpoints represent
exact observations under the shared `censored_data` convention.

The `*_fit` overloads take `censored_data` and bounds in `(a,loc,scale)`,
`(a,b,loc,scale)`, or `(loc,scale)` order for inverse-gamma, beta-prime, or
Maxwell respectively. Any parameter can be free or fixed by equal bounds.
Beta-prime free-location fits with `a<1` can have unbounded likelihood as
location approaches an exact observation. A fully fixed singular or
zero-density objective returns `fit_status_no_finite_objective`; a finite
bounded optimizer result does not guarantee a global MLE.

```fortran
use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
    invgamma_fit, betaprime_fit, maxwell_fit
type(censored_data) :: data
type(fit_result) :: fitted

data = make_censored_data(uncensored=[0.5_dp,1.2_dp,2.5_dp], &
    left=[0.3_dp],right=[3.5_dp],interval=reshape([1.1_dp,1.4_dp],[1,2]))
call invgamma_fit(data,[0.1_dp,0.0_dp,0.2_dp],[5.0_dp,0.0_dp,5.0_dp],fitted)
call betaprime_fit(data,[1.1_dp,0.5_dp,0.0_dp,0.2_dp], &
    [5.0_dp,5.0_dp,0.0_dp,5.0_dp],fitted)
call maxwell_fit(data,[0.0_dp,0.2_dp],[0.0_dp,5.0_dp],fitted)
```

Binary64 and the existing incomplete-gamma/beta kernels' numerical limits
apply. For very large shapes, coordinate rounding can limit near-mode
accuracy in sharply concentrated distributions. Subnormal values retain
binary64 resolution limits. Existing scalar and array-fit APIs are preserved;
C ABI and wrapper additions remain deferred. See the
[validation report](../validation/positive_transform_censored_2026-11-16.md).

## Power-normal, power-lognormal, and skew-Cauchy details

`powernorm` uses finite `c > 0`, location `loc`, and positive scale.
Writing `z=(x-loc)/scale`, its survival is `Phi(-z)**c` and its density is
`c*phi(z)*Phi(-z)**(c-1)/scale`. At `c=1` it is the normal distribution.
`powerlognorm` has shapes `(c,s)` in that order, both positive and finite.
It applies the same power-normal law at `w=log((x-loc)/scale)/s`, with
Jacobian `1/(s*(x-loc))`. Its support is `x > loc`; density and CDF at or
below `loc` are zero and survival is one. At `c=1` it equals lognormal with
shape `s`. All parameters must be finite, and invalid parameters return NaN
even for empty samples. Optional location and scale default to zero and one.

The power kernels retain the logarithm of cumulative hazard
`H=-c*log(Phi(-w))`. Normal upper-tail asymptotics avoid overflowing `w*w`
before multiplying by a tiny power and avoid subtracting two large log
probabilities in the density. For narrow intervals the hazard increment
is integrated in normal coordinates before evaluating
`-H(a)+log(1-exp(-(H(b)-H(a))))`. Original interval widths are retained
when transformed endpoint coordinates round together. Broader intervals
use the better tail difference and complement both omitted tails for
nearly full intervals. These routines do not change the existing scalar
distribution implementations, whose extreme-parameter accuracy can differ.

`skewcauchy` uses finite `-1 < a < 1`. Its left and right scale factors are
`1-a` and `1+a`, respectively, and its density is
`1/(pi*scale*(1+(z/h)**2))` on each side with the appropriate `h`.
The probability to the left of location is `(1-a)/2`; to the right it is
`(1+a)/2`. At `a=0` it equals Cauchy. Both sides have unbounded support.
Logarithmic side scales avoid overflowing their products with `scale`;
positive arctangent differences retain adjacent and crossing interval masses.
Center probabilities and small inner masses are added directly when
complementing a nearly unit tail would lose accuracy near `a=+/-1`.

Each family has matching `_censored_nnlf` and a censored-data `_fit`
overload, exported by `scifort_stats`. Bounds follow `(c,loc,scale)`,
`(c,s,loc,scale)`, and `(a,loc,scale)`, respectively. Any parameter may be
free or fixed with equal bounds. The bounded local optimizer can require
more than its default 1,000 sweeps when power, location, and scale are
jointly free; check `success` and `status`, and adjust `max_iter` or bounds
when needed. Success does not certify a global maximum or identifiability.
Existing array-based fits and specific procedure names remain available.

```fortran
use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
    powernorm_fit, powerlognorm_fit, skewcauchy_fit
type(censored_data) :: sample
type(fit_result) :: fitted
sample = make_censored_data(uncensored=[0.5_dp,1.2_dp,2.5_dp], &
    left=[0.3_dp],right=[3.5_dp],interval=reshape([1.1_dp,1.4_dp],[1,2]))
call powernorm_fit(sample,[0.05_dp,0.0_dp,1.3_dp], &
    [3.0_dp,0.0_dp,1.3_dp],fitted)
call powerlognorm_fit(sample,[0.05_dp,0.8_dp,0.0_dp,1.3_dp], &
    [3.0_dp,0.8_dp,0.0_dp,1.3_dp],fitted)
call skewcauchy_fit(sample,[-0.95_dp,0.0_dp,1.3_dp], &
    [0.95_dp,0.0_dp,1.3_dp],fitted)
```

The formulas follow the mathematical specifications for
[power-normal](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.powernorm.html),
[power-lognormal](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.powerlognorm.html),
and [skew-Cauchy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.skewcauchy.html).
See [validation and commands](../validation/power_skew_censored_2026-11-16.md).

## Exponential-power, exponentiated Weibull, and Mielke details

All three families use finite `loc` (default zero), positive finite `scale`
(default one), and positive finite shapes. Their support is `x >= loc`;
CDF is zero at and below location, survival is one, and intervals are
clipped to the support. A nonempty interval entirely below the support
has log probability negative infinity. Invalid parameters return NaN even
for empty samples. Shape order matches SciPy and existing SciFort APIs.

- `exponpow`, shape `b`: write `z=(x-loc)/scale`, `t=z**b`. Survival is
  `exp(-expm1(t))`, with log density
  `log(b)-log(x-loc)+log(t)+t-expm1(t)` above location. At location the
  density is infinite for `b<1`, `1/scale` for `b=1`, and zero for `b>1`.
  At `b=1` it equals Gompertz with shape one. This is SciPy's
  exponential-power family, distinct from generalized normal (`gennorm`).
- `exponweib`, shapes `(a,c)`: `t=z**c` and CDF is
  `(1-exp(-t))**a`. At location, the binary64 product `a*c` determines
  infinite, finite `1/scale`, or zero density as it is below, equal to,
  or above one. An overflowing positive product gives zero endpoint
  density. At `a=1` the family equals minimum Weibull with shape `c`.
- `mielke`, shapes `(k,s)`: CDF is `(1+z**(-s))**(-k/s)` above location.
  At location, density is infinite for `k<1`, `1/scale` for `k=1`, and
  zero for `k>1`. It equals Burr III with `c=s,d=k/s` when the ratio is
  representable, and Fisk with shape `s` when `k=s`. The censored kernels
  also handle shape ratios outside the binary64 range in log form.

The new kernels keep powers, cumulative hazards, and shape ratios in log
form until safe to evaluate. Exponentiated Weibull combines opposing
shapes before an overflowing transformed log can discard a finite result.
Its density uses `a*c/(x-loc)` times `t/expm1(t)` times its CDF, which
avoids cancellation between large shape-dependent terms. Mielke uses
the equivalent Burr III identities directly in `(k,s)` coordinates.

Narrow intervals retain the original endpoint gap when transformed
endpoints round together. If `delta=c*log((b-loc)/(a-loc))`, the increment
in `t` is `t(a)*expm1(delta)`. Exponential-power then increments its
hazard by `exp(t(a))*expm1(t(b)-t(a))`. Exponentiated Weibull and Mielke
use positive differences of reversed cumulative hazards. Broader and
nearly full intervals preserve the smaller tail or both omitted tails.

Matching `_censored_nnlf` and `_fit` overloads are exported from
`scifort_stats`. Bounds follow `(b,loc,scale)`, `(a,c,loc,scale)`, and
`(k,s,loc,scale)`, respectively. Any parameter can be free or fixed with
equal bounds; existing array-based fits remain available. The optimizer
is bounded and local. Jointly free shapes/location/scale can need more
than the default 1,000 sweeps, so inspect `success` and `status`.
Unrestricted location fits may have unbounded likelihood at exact
observations when the endpoint density is singular; a successful bounded
local fit does not establish a finite global MLE.

The scalar exponentiated-Weibull log-CDF now uses a complementary `log1p`
identity in the upper tail, fixing the observed loss in scalar log-SF.
Other scalar extreme-parameter limitations are unchanged; the new censored
kernels are independently validated for their documented ranges.

```fortran
use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
    exponpow_fit, exponweib_fit, mielke_fit
type(censored_data) :: sample
type(fit_result) :: fitted
sample = make_censored_data(uncensored=[0.5_dp,1.2_dp,2.5_dp], &
    left=[0.3_dp],right=[3.5_dp],interval=reshape([1.1_dp,1.4_dp],[1,2]))
call exponpow_fit(sample,[0.05_dp,0.0_dp,1.3_dp], &
    [3.0_dp,0.0_dp,1.3_dp],fitted)
call exponweib_fit(sample,[0.05_dp,1.2_dp,0.0_dp,1.3_dp], &
    [3.0_dp,1.2_dp,0.0_dp,1.3_dp],fitted)
call mielke_fit(sample,[0.05_dp,1.2_dp,0.0_dp,1.3_dp], &
    [3.0_dp,1.2_dp,0.0_dp,1.3_dp],fitted)
```

The formulas follow the mathematical specifications for
[exponential-power](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.exponpow.html),
[exponentiated Weibull](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.exponweib.html),
and [Mielke](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mielke.html).
See [validation and commands](../validation/lifetime_censored_2026-11-16.md).

## Fatigue-life, Johnson SU, and generalized exponential details

`fatiguelife` uses positive finite shape `c`, finite `loc`, and positive
finite `scale`. With `z=(x-loc)/scale`, its normal coordinate is
`w=(sqrt(z)-1/sqrt(z))/c`. Support is `x > loc`; density and CDF at or
below location are zero and survival is one. The log density is the
standard normal log density at `w` plus
`log((sqrt(z)+1/sqrt(z))/(2*c*(x-loc)))`. Log-distance and `expm1`
identities retain near-scale displacements and normalized distances
outside the binary64 range. At zero location/unit scale, reciprocal
coordinates exchange CDF and survival.

`johnsonsu` uses finite `a` of either sign and positive finite `b`, with
finite location and positive finite scale. Its support is the real line.
The normal coordinate is `q=a+b*asinh(z)`, with log Jacobian
`log(b)-log(scale)-log(hypot(1,z))`. Log coordinates retain very large
or underflowing normalized displacements; guarded arithmetic prevents
overflowing intermediate products and sums when the result is finite.
Negating `a`, the observation, and location reflects the distribution.
Very large opposing shape terms make `q` intrinsically ill-conditioned:
rounding in `asinh(z)` can be greatly amplified. Binary64 accuracy near
such cancellation is limited and is not established by ordinary cases.

`genexpon` uses three positive finite shapes `(a,b,c)`. It has lower
endpoint `loc`, with density `a/scale` there. CDF is zero and survival
is one at and below the endpoint. With `v=c*z`, cumulative hazard is
`H=a*z+(b/c)*(v+expm1(-v))`, and density is
`(a+b*(1-exp(-v)))*exp(-H)/scale`. The kernels keep the two nonnegative
hazard terms in log form, using a series for the small remainder
`v+expm1(-v)`. This permits compensating extreme rates without forming
an overflowing `a+b`, an underflowing ratio, or subtracting large hazards.
As `b` tends to zero, it approaches the exponential family with scale
`scale/a`.

All optional locations/scales default to zero/one. Invalid parameters
return NaN even for empty data. The matching `_censored_nnlf` and
censored-data `_fit` overloads are available through `scifort_stats`.
Bounds follow `(c,loc,scale)`, `(a,b,loc,scale)`, and `(a,b,c,loc,scale)`,
respectively. Any parameter can be free or fixed using equal bounds;
original array-based fit procedures remain available.

For narrow fatigue-life intervals, the normal-coordinate width is
`expm1(log(zb/za)/2)*(sqrt(za)+1/sqrt(zb))/c`. Johnson SU uses
`asinh(zb)-asinh(za)=2*atanh((zb-za)/(hypot(1,za)+hypot(1,zb)))`.
The original endpoint gap is retained independently of rounded endpoint
coordinates. Both widths feed the existing narrow normal integral;
broader intervals use the better tail difference and preserve both
omitted tails for nearly full intervals.

Generalized exponential integrates its hazard with a sum of nonnegative
terms. For lower standardized point `z`, width `d`, and
`f(v)=v+expm1(-v)`, the increment is
`a*d+b*d*(1-exp(-c*z))+(b/c)*exp(-c*z)*f(c*d)`.
This keeps adjacent interval probabilities even when both endpoint
cumulative hazards round to the same value.

```fortran
use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
    fatiguelife_fit, johnsonsu_fit, genexpon_fit
type(censored_data) :: sample
type(fit_result) :: fitted
sample = make_censored_data(uncensored=[0.5_dp,1.2_dp,2.5_dp], &
    left=[0.3_dp],right=[3.5_dp],interval=reshape([1.1_dp,1.4_dp],[1,2]))
call fatiguelife_fit(sample,[0.05_dp,0.0_dp,1.3_dp], &
    [3.0_dp,0.0_dp,1.3_dp],fitted)
call johnsonsu_fit(sample,[-3.0_dp,1.2_dp,0.0_dp,1.3_dp], &
    [3.0_dp,1.2_dp,0.0_dp,1.3_dp],fitted)
call genexpon_fit(sample,[0.05_dp,1.2_dp,0.8_dp,0.0_dp,1.3_dp], &
    [3.0_dp,1.2_dp,0.8_dp,0.0_dp,1.3_dp],fitted)
```

Fits are bounded local optimizations. Jointly free shapes, location, and
scale can be weakly identifiable and may need more than 1,000 sweeps;
check `success` and `status`. Success does not certify a global MLE.
The new censored kernels are independently validated; existing scalar
distribution implementations and their extreme-parameter limits are unchanged.

The formulas follow the specifications for
[fatigue-life](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.fatiguelife.html),
[Johnson SU](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.johnsonsu.html),
and [generalized exponential](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.genexpon.html).
See [validation and commands](../validation/transform_hazard_censored_2026-11-16.md).


## Cosine, anglit, and semicircular details

The three likelihood/NNLF pairs accept optional finite `loc` and positive
finite `scale`. Their fit generics accept bounds in `[loc, scale]` order;
equal bounds fix either parameter. Original scalar and array-fit calls
remain available. Invalid parameters return NaN even for empty data.
Empty valid likelihoods are zero. Exact observations outside the support,
and exact observations at a mathematical support endpoint, contribute
negative infinity because the endpoint density is zero.

These censored likelihoods use a **rounded physical half-width** contract.
Compute `r=h*scale` in binary64, where `h=pi`, `pi/4`, or `1` for cosine,
anglit, or semicircular, respectively. The probability law uses center
`loc` and this `r`, with density normalization based on `r`. This explicitly
retains a normalized law even when multiplication rounds a subnormal scale.
The rounded half-width must be finite and positive; both rounded support
edges `loc-r` and `loc+r` must be finite and strictly straddle `loc`.
Otherwise likelihood/NNLF returns NaN, including for empty samples.
This restricts representable supports compared with the scalar APIs.

Endpoint gaps include the error-free residual of the center/half-width
sum. Consequently a point at a *rounded* edge may lie just inside or
outside the mathematical endpoint `loc +/- r`; it is not automatically
an endpoint observation. For ordinary scales, the rounded-half-width law
agrees with the scalar parameterization to rounding accuracy. At subnormal
scales, multiplication can have a large relative rounding error, so scalar
and censored likelihood values need not agree closely. This numerical
contract is intentional and limited to these censored likelihoods/fits.

Let `q=(r-abs(x-loc))/r` be the nearer endpoint gap, formed with compensated
distances before division. The smaller tail probabilities are
`(q-sin(pi*q)/pi)/2`, `sin(pi*q/4)**2`, and
`(acos(1-q)-(1-q)*sqrt(q*(2-q)))/pi`, respectively. The implementation uses
integrated sine and binomial series where these expressions cancel, and
retains logarithmic gaps when a standardized gap underflows. Densities use
positive sine or square-root identities instead of cancellation-prone
trigonometric sums. The other tail uses `log1p` of the smaller probability.

For narrow intervals, three-point Gauss-Legendre integration operates in
endpoint-gap coordinates, preserving widths between adjacent floating-point
observations. Broad intervals select a tail difference; nearly full
intervals complement both omitted tails. Fits are bounded local searches,
with the same iteration and identifiability limits as other censored fits;
they do not guarantee a global maximum.

For example, fit anglit location and scale to a mixed sample:

```fortran
use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
    anglit_fit, anglit_censored_nnlf
type(censored_data) :: sample
type(fit_result) :: fit

sample = make_censored_data(uncensored=[-0.2_dp,0.1_dp,0.35_dp],left=[-0.4_dp], &
    right=[0.45_dp],interval=reshape([-0.3_dp,0.2_dp,-0.1_dp,0.25_dp],[2,2]))
call anglit_fit(sample,[-0.2_dp,0.7_dp],[0.2_dp,2.0_dp],fit,max_iter=10000)
if (fit%success) print *, fit%params, &
    anglit_censored_nnlf(sample,fit%params(1),fit%params(2))
```

See the [validation record](../validation/compact_censored_2026-11-16.md)
for checked/optimized results, independent reference checks, and commands.


## Truncated exponential, Johnson SB, and kappa3 details

These families support likelihoods and matching `_censored_nnlf` procedures
from `scifort_censored` and `scifort_stats`. Their fit generics accept
`censored_data` with bounds in `[b, loc, scale]`, `[a, b, loc, scale]`, or
`[a, loc, scale]` order, respectively. Equal bounds fix parameters. The
original scalar and array-fitting APIs remain available. Invalid parameters
return NaN even for empty samples; valid empty likelihoods are zero.

### Truncated exponential

Shape `b` is positive and finite. The censored implementation uses a
**rounded physical width**: `r=b*scale` is evaluated in binary64, and the
probability density is proportional to `exp(-(x-loc)/scale)` on
`[loc,loc+r]`. Its exact mathematical normalizer uses `r/scale`, rather
than replacing that ratio by `b`. Thus the law stays normalized when the
width multiplication rounds. This intentional numerical difference from
the scalar API is most visible for subnormal widths.

The rounded width must be positive and finite; the rounded upper support
edge must be finite and strictly greater than `loc`. Underflowed,
overflowed, and collapsed supports return NaN, including for empty samples.
Endpoint gaps retain the error-free residual of `loc+r`, so a point at
the rounded edge need not be the exact mathematical endpoint.

The endpoint density is positive, with log density
`-(x-loc)/scale-log(scale)-log(1-exp(-r/scale))`. Lower/upper support
endpoints have CDF 0/1 and SF 1/0. Exact observations outside support
contribute negative infinity. Tails and interval probabilities use positive
exponential mass increments. Nearly full intervals complement both omitted
tails. A centered hyperbolic-sine identity retains the log density in the
nearly uniform limit, including representable corrections as small as
`2.5e-301` that would be lost by subtracting large logarithms. Subnormal
widths avoid forming an incorrectly rounded half-width midpoint.

### Johnson SB

The first shape `a` is any finite real; `b` is positive and finite.
Support is `[loc,loc+scale]`, with positive finite width `scale`.
The rounded upper edge must be finite and strictly above `loc`.
Otherwise the likelihood/NNLF returns NaN. Compensated distances follow
the same upper-edge convention as the beta likelihood. Endpoint and
outside-support densities are zero. Endpoint CDF/SF behavior follows the
usual bounded-support conventions.

For `z=(x-loc)/scale`, the normal coordinate is
`q=a+b*log(z/(1-z))`. Independent endpoint gaps preserve both sides of
the logit. Near the midpoint, a compensated center displacement preserves
changes that a large `b` can amplify, even when `z` rounds to one half.
Odd subnormal widths use the difference between the full endpoint gaps
instead of a rounded midpoint. Products and sums are guarded against
overflow, including an opposing `a` that rescues an overflowing product.
Severe cancellation between large opposing shape terms remains
ill-conditioned; uniform relative accuracy through that cancellation is
not promised.

For an interior interval `(x1,x2)`, the exact logit width is
`log1p((x2-x1)/(x1-loc))+log1p((x2-x1)/(loc+scale-x2))`.
The implementation retains that width logarithmically and integrates in
normal coordinates before falling back to stable tail differences.

### Kappa3

Shape `a` is positive and finite, with support `[loc,infinity)`.
There is no finite upper-edge or representable `loc+scale` requirement.
The endpoint density is `a**(-1/a)/scale`; its log can be positive infinity
for a sufficiently small valid shape. Endpoint CDF/SF values are 0/1.
Below-support densities and probabilities follow the usual support rules.

The implementation computes the reversed cumulative hazard
`H=log1p(a*z**(-a))/a`, where `z=(x-loc)/scale`. Logarithmic softplus
identities retain tiny shapes and shape-times-log-distance products that
overflow individually. Compensated displacement at the scale point also
preserves large-shape effects when the ordinary standardized distance
rounds to one. For `a=1`, kappa3 agrees with Lomax shape one. As `a`
tends to zero, finite positive standardized observations have CDF tending
to `exp(-1)`; the code retains this limit instead of underflowing a
reciprocal shape or normalizing constant.

Interval probabilities use positive reversed-hazard increments. The
large-positive softplus branch has a relative correction below `5e-18`.
Allowing free location with exact endpoint observations can produce
unbounded likelihoods as parameters approach singular limits. Bound the
parameters appropriately and inspect fit diagnostics; bounded local
searches do not establish a global maximum or a unique identifiable fit.

For example, fit a Johnson SB shape with the other parameters fixed:

```fortran
use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
    johnsonsb_fit, johnsonsb_censored_nnlf
type(censored_data) :: sample
type(fit_result) :: fit

sample = make_censored_data(uncensored=[0.2_dp,0.4_dp,0.6_dp],left=[0.1_dp], &
    right=[0.7_dp],interval=reshape([0.25_dp,0.5_dp,0.35_dp,0.55_dp],[2,2]))
call johnsonsb_fit(sample,[-3.0_dp,1.2_dp,0.0_dp,1.3_dp], &
    [3.0_dp,1.2_dp,0.0_dp,1.3_dp],fit)
if (fit%success) print *, fit%params, &
    johnsonsb_censored_nnlf(sample,fit%params(1),fit%params(2), &
        fit%params(3),fit%params(4))
```

See the [validation record](../validation/bounded_shape_censored_2026-11-17.md)
for checked/optimized results, independent references, and commands.


## Log-uniform, truncated Pareto, and generalized half-logistic details

The likelihoods and matching `_censored_nnlf` functions are available from
`scifort_censored` and `scifort_stats`. The existing fit generics also accept
`censored_data`, with bounds ordered `[a,b,loc,scale]` for log-uniform,
`[b,c,loc,scale]` for truncated Pareto, and `[c,loc,scale]` for generalized
half-logistic. Equal bounds fix parameters. Existing scalar and array-fit
APIs remain available. Invalid parameters return NaN even for empty samples;
valid empty likelihoods are zero.

### Log-uniform and truncated Pareto

Log-uniform requires finite shapes `0<a<b`. Truncated Pareto requires finite
`b!=0` (either sign) and finite `c>1`. Both use a normalized law on rounded
**physical cutoffs** `L,U` for `t=x-loc`:

- Log-uniform: `L=fl(a*scale)`, `U=fl(b*scale)`, density `1/(t*log(U/L))`.
- Truncated Pareto: `L=scale`, `U=fl(c*scale)`, density proportional to
  `t**(-b-1)` on `[L,U]`.

Here `fl` denotes binary64 rounding. Normalization uses the rounded physical
cutoffs, so rounding the support does not destroy normalization. This is an
intentional numerical difference from the original scalar APIs, especially
for subnormal cutoffs. Require `0<L<U`, both finite, with finite rounded
translated edges satisfying `loc<fl(loc+L)<fl(loc+U)`. Underflowed,
overflowed, or collapsed supports return NaN. Compensated endpoint gaps
retain addition residuals, so a rounded translated edge need not equal the
exact mathematical endpoint.

Density is positive at both exact support endpoints. CDF values there are
0 and 1; SF values are 1 and 0. Exact observations outside support contribute
negative infinity. The Pareto law tends to log-uniform as `b` approaches
zero from either side, but `b=0` itself remains invalid for compatibility.
Negative shapes concentrate toward the upper cutoff and are fully supported
by the native censored likelihood.

The implementation integrates the power law in logarithmic coordinates.
Positive endpoint and interval increments use `log1p`, `expm1`, and
logarithmic hazards. It selects the smaller tail and obtains its complement
stably; nearly full intervals complement both omitted tails. A centered
hyperbolic-sine normalization avoids cancellation for nearly zero powers.
Adjacent support cutoffs and adjacent floating-point interval endpoints do
not require subtracting nearly equal CDF values.

For fitting truncated Pareto with unknown shape sign, search negative-only
and positive-only bounds separately and compare finite objectives. A local
search spanning the excluded point `b=0` can stop on the poorer side.
Log-uniform's shape/scale parameterization is redundant when all three are
free; fixing scale or another appropriate parameter improves identifiability.

### Generalized half-logistic

Shape `c` is positive and finite. Define the rounded physical width
`r=fl(scale/c)`, the support `[loc,loc+r]`, and `u=(x-loc)/r`.
The normalized law is
`F=(1-(1-u)**(1/c))/(1+(1-u)**(1/c))`.
Width `r` must be positive and finite, and the rounded upper edge must be
finite and strictly greater than `loc`; otherwise return NaN. The original
scalar API remains unchanged. Width rounding can produce numerical
differences, especially for subnormal scales. Endpoint gaps retain addition
residuals, as in the other bounded censored kernels.

At the lower endpoint the density is `1/(2*c*r)`. At the upper endpoint it
is zero for `c<1`, `2/r` for `c=1`, and positive infinity for `c>1`.
Endpoint CDF/SF values follow the bounded-support conventions. For `c=1`,
the interior density is `2/(r*(2-u)**2)` and the CDF is `u/(2-u)`.

The transformed hazard `H=-log1p(-u)/c` is retained logarithmically, avoiding
an overflowing reciprocal shape or an underflowed standardized observation.
Replacing `-log1p(-u)` by `u` below `epsilon/exp(2)` has relative correction
below `epsilon/(2*exp(2))`. Interior interval probabilities use the positive
increment `DeltaH=log1p((x2-x1)/(loc+r-x2))/c`, rather than a CDF subtraction.
This retains tiny shapes, extreme tails, and subnormal interval widths.

Endpoint singularities or movable support edges can make unrestricted
likelihood maximization unbounded. Fits are bounded local searches; success
does not establish a global maximum or unique identifiable parameters.

For example, fit a positive truncated Pareto shape with fixed support:

```fortran
use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
    truncpareto_fit, truncpareto_censored_nnlf
type(censored_data) :: sample
type(fit_result) :: fit

sample = make_censored_data(uncensored=[1.2_dp,1.4_dp,1.6_dp],left=[1.1_dp], &
    right=[1.7_dp],interval=reshape([1.25_dp,1.5_dp,1.35_dp,1.55_dp],[2,2]))
call truncpareto_fit(sample,[0.01_dp,3.0_dp,0.0_dp,1.0_dp], &
    [5.0_dp,3.0_dp,0.0_dp,1.0_dp],fit)
if (fit%success) print *, fit%params, &
    truncpareto_censored_nnlf(sample,fit%params(1),fit%params(2), &
        fit%params(3),fit%params(4))
```

See the [validation record](../validation/bounded_power_censored_2026-11-17.md)
for checked/optimized tests, independent references, and exact commands.


## Folded normal, folded Cauchy, and generalized Pareto details

These families provide `_censored_loglikelihood` and `_censored_nnlf` from
`scifort_censored` and `scifort_stats`. Their fit generics accept
`censored_data` with bounds in `[c,loc,scale]` order. Equal bounds fix
parameters. The original scalar and array-fit APIs remain available.
All parameters must be finite, scale positive. Invalid parameters return
NaN even for empty samples; valid empty likelihoods are zero.

### Folded normal and folded Cauchy

The folding offset `c` is nonnegative, including zero. In standardized
coordinates `z=(x-loc)/scale`, support is `[0,infinity)`. The density is
the sum of the base normal or Cauchy densities at `z-c` and `z+c`, divided
by scale. At the lower endpoint it is twice the base density at `c`,
divided by scale. CDF/SF values there are 0/1. Below support the density
and CDF are zero, and SF is one. At zero folding offset these laws agree
with half-normal and half-Cauchy.

The lower tail is a base-law interval `[-z-c,z-c]`; the upper tail is the
sum of the two base survival probabilities. Censoring intervals add the
positive masses in the two folded branches. Normal intervals use the
existing narrow-interval quadrature or stable normal-tail differences.
Cauchy intervals use positive arctangent identities and logarithmic
coordinates. Both retain the physical interval width independently of
rounded standardized endpoints, including subnormal widths. Nearly full
intervals complement both omitted tails.

Guarded standardized coordinates and logarithmic Cauchy distances avoid
overflow in squared distances, sums, and ratios. A normal log density or
log tail becomes negative infinity when its magnitude exceeds the binary64
magnitude range. Cauchy log densities and tails remain finite even when
the standardized distance overflows. Halving can recover a representable
normal displacement from an overflowing standardized distance. However,
subtraction of a large folding offset from a nearly equal standardized
distance is ill-conditioned; uniform relative accuracy through such
cancellations is not promised.

### Generalized Pareto

Shape `c` is any finite real. For `c>=0`, support is `[loc,infinity)`.
For positive shape the cumulative hazard is
`H=log1p(c*(x-loc)/scale)/c`; the zero-shape limit is exponential with
`H=(x-loc)/scale`. The implementation retains H logarithmically, avoiding
an overflowing reciprocal shape or an underflowed shape-distance product.
There is no representable location-plus-scale requirement in these cases.

For `c<0`, define a **rounded physical width** `r=fl(scale/(-c))`, with
binary64 rounding `fl`. Support is `[loc,loc+r]` and the normalized SF is
`(1-(x-loc)/r)**(-1/c)`. Width r must be positive and finite, and the
rounded upper edge must be finite and strictly greater than loc.
Underflowed, overflowed, or collapsed supports return NaN. Compensated
endpoint gaps retain the addition residual of loc+r, so a rounded edge
need not equal the exact mathematical endpoint. Normalization uses r;
this intentional numerical difference from the scalar API is especially
visible for subnormal widths.

The bounded lower endpoint density is `1/((-c)*r)`. At the upper endpoint
it is zero for `-1<c<0`, `1/r` for `c=-1`, and positive infinity for
`c<-1`. The unbounded lower endpoint density is `1/scale`. Endpoint tails
follow the support conventions. At c=-1 the law is uniform on the bounded
support; at c=1 it agrees with Lomax shape one.

An interior interval has probability `exp(-H1)*(1-exp(-DeltaH))`.
For positive c,
`DeltaH=log1p(c*(x2-x1)/(scale+c*(x1-loc)))/c`;
for negative c,
`DeltaH=log1p((x2-x1)/(loc+r-x2))/(-c)`.
Logarithmic positive increments retain adjacent intervals without
subtracting CDF values. Tiny bounded coordinates reuse the independently
validated logarithmic hazard transformation from generalized half-logistic.

Free support endpoints and singular limits can make likelihood maximization
unbounded. Fits are bounded local searches and do not establish global
optimality or unique identifiable parameters.

For example, fit a generalized Pareto shape with location and scale fixed:

```fortran
use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
    genpareto_fit, genpareto_censored_nnlf
type(censored_data) :: sample
type(fit_result) :: fit

sample = make_censored_data(uncensored=[0.2_dp,0.4_dp,0.6_dp],left=[0.1_dp], &
    right=[0.7_dp],interval=reshape([0.25_dp,0.5_dp,0.35_dp,0.55_dp],[2,2]))
call genpareto_fit(sample,[-1.0_dp,0.0_dp,1.3_dp],[3.0_dp,0.0_dp,1.3_dp],fit)
if (fit%success) print *, fit%params, &
    genpareto_censored_nnlf(sample,fit%params(1),fit%params(2),fit%params(3))
```

See the [validation record](../validation/fold_pareto_censored_2026-11-17.md)
for checked/optimized results, independent references, and exact commands.


## Generalized normal, half generalized normal, and R details

The likelihoods and matching `_censored_nnlf` functions are available from
`scifort_censored` and `scifort_stats`. Their fit generics accept
`censored_data` with bounds ordered `[beta,loc,scale]` for generalized normal
and its half law, or `[c,loc,scale]` for R. Equal bounds fix parameters.
Original scalar and array-fit APIs remain available. Invalid parameters
return NaN even for empty data; valid empty likelihoods are zero.

### Generalized normal and half generalized normal

Power beta is positive and finite, and its binary64 reciprocal must be
finite. Location is finite and scale positive and finite. A reciprocal
that overflows returns NaN; extending that range remains a compatibility
gap. An overflowing gamma log normalizer for a valid reciprocal produces
the appropriate infinite log likelihood, while valid empty data still
returns zero.

With `z=(x-loc)/scale`, the half-law density on `[loc,infinity)` is
`exp(-z**beta)/(scale*Gamma(1+1/beta))`. The symmetric law has half this
density, with z replaced by abs(z), on the whole real line. Its CDF/SF at
loc are both one half. The half-law CDF/SF at loc are 0/1, with positive
endpoint density `1/(scale*Gamma(1+1/beta))`; below loc its density/CDF are
zero and SF one.

The powered distance `Y=abs(z)**beta` follows a unit-scale gamma law with
shape `1/beta`. Independent gamma lower/upper logarithmic tails retain
small probabilities. Log-gamma-one-plus and a grouped Stirling expression
avoid cancellation and intermediate overflow in the normalization. The
implementation retains the tiny positive endpoint log density in the
nearly uniform half-law limit. Compensated scale-point displacements
preserve effects amplified by large beta even if the ordinary standardized
coordinate rounds to one.

For very small Y, the identity `Y**(1/beta)=abs(z)` is evaluated using
the original log distance. This preserves finite CDFs even when the
intermediate product beta*log(abs(z)) overflows negatively. Subnormal
log CDFs at the largest finite power remain representable.

Intervals on one side use gamma tail differences or three-point quadrature
in the original coordinate. Quadrature requires small relative width,
powered-coordinate change, and density change. Independently stored log
widths retain adjacent floating-point intervals. If the powered coordinate
is uniformly below epsilon/exp(2), a constant-density integral has relative
error below that threshold. Symmetric intervals reflect negative ranges,
split at loc when needed, and complement omitted tails for nearly full
support. At beta=1 these laws agree with Laplace and exponential; at
beta=2 they agree with normal and half-normal with scale divided by sqrt(2).

### R distribution

The current censored implementation supports finite `0<c<=1e6` and requires
the rounded half-shape `a=fl(c/2)>0`. Larger c and an underflowed half-shape
return NaN, including for empty data. This explicit upper limit bounds the
validated beta continued-fraction regime; extending it is a remaining
numerical compatibility gap.

The implementation uses a normalized symmetric beta law with shapes a,a
on physical support `[L,L+w]`, where `L=fl(loc-scale)` and `w=2*scale`.
Require positive finite scale, finite positive w, finite L and rounded upper
edge, and `L<loc<fl(L+w)`. Underflowed/overflowed/collapsed transforms return
NaN. The rounded lower edge defines the numerical law; subnormal or shifted
supports can differ from the scalar API. Compensated distances retain the
residual of L+w. R scale denotes a half-width; its loc denotes the requested
center before lower-edge rounding.

Endpoint density is infinite for c<2, `1/w` for c=2, and zero for c>2.
Endpoint CDF/SF values are 0/1 and 1/0. R with c=1,2,3 agrees with the
corresponding arcsine, uniform, and semicircular laws when their physical
supports coincide.

For half-shape at most epsilon/4096, positive density and logit-interval
formulas avoid intrinsic gamma overflow at a subnormal shape. Interior
tails round to one half with absolute error below epsilon. The leading
density/interval formulas have relative errors below binary64 roundoff
over representable endpoint coordinates. The remaining range reuses the
existing beta tails and interval kernels. At an exactly symmetric stored
coordinate, both tails are one half. All parameters can vary in bounded
fits; movable endpoints and shape/scale limits can make unrestricted
likelihoods unbounded. Fit success denotes a local bounded search.

For example, fit a half generalized-normal power with location and scale fixed:

```fortran
use scifort_stats, only : dp, censored_data, make_censored_data, fit_result, &
    halfgennorm_fit, halfgennorm_censored_nnlf
type(censored_data) :: sample
type(fit_result) :: fit

sample = make_censored_data(uncensored=[0.2_dp,0.4_dp,0.6_dp],left=[0.1_dp], &
    right=[0.7_dp],interval=reshape([0.25_dp,0.5_dp,0.35_dp,0.55_dp],[2,2]))
call halfgennorm_fit(sample,[0.3_dp,0.0_dp,1.3_dp],[5.0_dp,0.0_dp,1.3_dp],fit)
if (fit%success) print *, fit%params, &
    halfgennorm_censored_nnlf(sample,fit%params(1),fit%params(2),fit%params(3))
```

See the [validation record](../validation/generalized_normal_censored_2026-11-17.md)
for checked/optimized results, references, and exact commands.


## Chi-squared, chi, and Nakagami details

The likelihoods and matching `*_censored_nnlf` functions are exported by
`scifort_censored` and `scifort_stats`. The existing fit generics accept
`censored_data` with three bounds in `[df,loc,scale]` order for `chi2`/`chi`,
or `[nu,loc,scale]` for `nakagami`. Equal bounds fix any parameter; all three
may vary. Location defaults to zero and scale to one. Shape and scale must
be positive and finite, location finite. Chi and chi-squared additionally
require the rounded binary64 half-shape `fl(df/2)>0`; an underflowed half-shape
returns NaN, including for empty data. Valid empty likelihoods are zero.

The current numerical contract limits gamma shape a to at most 1e6:
`df<=2e6` for chi/chi-squared and `nu<=1e6` for Nakagami. Larger shapes
return NaN even for empty data. Extending this range needs a validated
large-shape treatment; it remains a compatibility gap.

All three laws have support `[loc,infinity)`. In standardized coordinates
`z=(x-loc)/scale`, the transformed variable `Y=rate*z**power` is gamma:

| Distribution | Gamma shape a | Rate | Power |
| --- | --- | --- | --- |
| Chi-squared | fl(df/2) | 1/2 | 1 |
| Chi | fl(df/2) | 1/2 | 2 |
| Nakagami | nu | nu | 2 |

The positive half-shape defines the numerical law at subnormal degrees of
freedom, where rounding can differ from the ideal real-number half-shape.
The density is `power*rate**a*z**(power*a-1)*exp(-Y)/(Gamma(a)*scale)`.
Its endpoint density is infinite, finite, or zero according as `power*a`
is below, equal to, or above one. CDF/SF at loc are 0/1; below loc the
density/CDF are zero and SF one. The finite endpoint density is
`power*rate**a/(Gamma(a)*scale)`.

Compensated log distances avoid forming overflowing physical differences
or squared standardized coordinates. Log-gamma density and independent
gamma tails handle extreme transformed distances. The implementation retains
tiny-shape survival logarithms through `Q(a,Y)=a*E1(Y)*(1+O(a*(1+abs(log(Y)))))`
when a<=epsilon/4096 and Y<=1.5. The omitted relative correction is below
roundoff for finite binary64 observation/scale log distances. For underflowed
Y, `E1(Y)=-EulerGamma-log(Y)+O(Y)` uses the original log distance. This avoids
rounding a subnormal survival probability before taking its logarithm.
Narrow intervals retain
their original physical log width and use three-point Gauss-Legendre
integration only when relative width and both log-density variation terms
are small. The analytic density ratios retain subnormal increments. Other
intervals use the better-conditioned gamma tail, or complement both omitted
tails for nearly full support. Existing special-function accuracy and
representability limits apply; this adds no new large-shape asymptotic kernel.

Chi with df=1/2/3 agrees with half-normal/Rayleigh/Maxwell respectively.
Chi-squared is gamma with shape df/2 and scale twice its supplied scale;
the implementation need not form that doubled scale. Nakagami with nu=1/2
is half-normal, and nu=1 is Rayleigh with scale divided by sqrt(2).

Fits use the existing local bounded optimizer. Check success, status, and
objective; an iteration limit can leave a useful finite candidate without
declaring convergence. Moving endpoints and singular limits may make a fit
unbounded. Bounds and a successful local fit do not establish global
optimality or identifiability. See the
[validation record](../validation/radial_gamma_censored_2026-11-17.md).


## Student t, F, and Jones-Faddy skew-t details

The likelihoods and corresponding `*_censored_nnlf` functions are exported
by `scifort_censored` and `scifort_stats`. Their existing fit generics accept
`censored_data`, with bounds in `[df,loc,scale]`, `[dfn,dfd,loc,scale]`, or
`[a,b,loc,scale]` order respectively. Equal bounds fix any parameter; all
parameters may vary. Original array-based fitting calls remain available.
Invalid parameters return NaN even for empty data; valid empty likelihoods
are zero. Location defaults to zero and scale to one.

The current numerical contract requires positive normal binary64 beta
parameters, each at least `tiny(1.0_dp)`, and their sum at most 1e6:

| Distribution | Beta parameters | Native shape limits |
| --- | --- | --- |
| Student t | df/2, df/2 | 2*tiny <= df <= 1e6 |
| F | dfn/2, dfd/2 | dfn,dfd >= 2*tiny; (dfn+dfd)/2 <= 1e6 |
| Jones-Faddy | a, b | a,b >= tiny; a+b <= 1e6 |

Student t additionally
accepts `df=+infinity` in its likelihood/NNLF as the normal limit. The
bounded fitter requires finite bounds, as elsewhere. Subnormal beta shapes
and larger finite shape sums remain explicit compatibility gaps. Existing
beta-kernel accuracy and representability limits apply.

Student t and Jones-Faddy have support on the whole real line. With
`z=(x-loc)/scale`, `s=a+b`, and `h=sqrt(s+z*z)`, Jones-Faddy maps monotonically
to a beta variable `u=(1+z/h)/2`, `v=1-u`. Student t is the symmetric case
`a=b=df/2`. Independent log distances and the small-coordinate identity
`s/(2*h*(h+abs(z)))` retain both tails when z or z*z would overflow.
Endpoint-rounded beta coordinates use logarithmic beta-tail limits.

For F, support is `[loc,infinity)`, and `y=(dfn/dfd)*z` is beta-prime with
shapes dfn/2 and dfd/2. The implementation retains log(dfn/dfd), avoiding
overflow in an effective scale or ratio. CDF/SF at loc are 0/1. Endpoint
density is infinite for dfn<2, `1/scale` for dfn=2, and zero for dfn>2.
Below loc the density/CDF are zero and SF one.

Narrow Student/Jones-Faddy intervals use a three-point physical-coordinate
density integral when log width, scale, h, and shape sum bound log-density
variation by 0.001. The bound is `abs(d log(f)/dz)<=(2*s+1)/h`, using the
smallest h on the interval. F reuses the validated beta-prime interval
formula. Other intervals use the better-conditioned logarithmic tail, or
complement both omitted tails for nearly full support.

When both beta shapes are at most epsilon/16384, finite-point CDF/SF use
the endpoint mixture weights b/(a+b) and a/(a+b) with stable logarithms.
The neglected relative terms are below binary64 roundoff for the finite
input log-coordinate range. Finite interval masses are still retained:
F integrates its leading log-ratio density, and Student/Jones-Faddy use
`2*a*b/(a+b)` times the increment of `asinh(z/sqrt(a+b))`. Stable increment
identities preserve widths even when rounded CDFs are equal.

Student t with df=1 agrees with Cauchy; symmetric Jones-Faddy agrees with
Student t; F agrees with the corresponding scaled beta-prime law. Fits
use the existing local bounded optimizer. Check status, success, and the
objective; iteration limits can leave finite improving candidates without
convergence. Success does not prove global optimality or identifiability.
See the [validation record](../validation/beta_t_censored_2026-11-17.md).


## Inverse Gaussian, Wald, and reciprocal inverse Gaussian details

`invgauss_censored_loglikelihood`, `wald_censored_loglikelihood`, and
`recipinvgauss_censored_loglikelihood`, and their `*_censored_nnlf`
counterparts, are exported by `scifort_censored` and `scifort_stats`.
Each accepts mixed exact, left-, right-, and interval-censored data.
The existing fit generics accept `censored_data`, with bounds in
`[mu_shape,loc,scale]` order for inverse Gaussian and its reciprocal,
and `[loc,scale]` for Wald. Equal bounds fix any parameter; all may vary.
The original array-fit calls remain available. SciPy names this shape `mu`;
SciFort preserves its established `mu_shape` dummy argument name.

Shape and scale must be positive and finite; loc must be finite. Defaults
are loc=0 and scale=1. Invalid parameters return NaN even for empty data;
valid empty likelihoods are zero, and empty fits are invalid. With
`z=(x-loc)/scale`, the support is x>loc: density and CDF are zero at/below
loc, and SF is one. Equal finite interval endpoints represent an exact
observation. Infinite one-sided interval endpoints follow the shared
conventions above; infinite exact/left/right observations are invalid.

For inverse Gaussian, the log density is
`-log(sqrt(2*pi))-1.5*log(z)-0.5*(z/mu_shape-1)**2/z-log(scale)`.
Wald is exactly the mu_shape=1 case. Reciprocal inverse Gaussian is
`1/Y` for standardized `Y~invgauss(mu_shape)`, then shifted and scaled;
its CDF/SF swap inverse-Gaussian SF/CDF at 1/z, and its standardized log
density has the Jacobian correction -2*log(z). These transformations
are applied through log distances, avoiding a rounded/overflowed 1/z.

Let `A=sqrt(z)/mu_shape-1/sqrt(z)`, `B=sqrt(z)/mu_shape+1/sqrt(z)`,
and `R(t)=Phi(-t)/phi(t)` be the normal Mills ratio. The second CDF term
`exp(2/mu_shape)*Phi(-B)` equals `phi(A)*R(B)`. This identity avoids
overflow and exponent cancellation. For A>=12 the SF uses the
inverse-power expansion of `R(A)-R(B)`, factoring out `(B-A)/B` before
summing so a nearly equal A and B do not erase the tail. The expansion
stops at binary64 roundoff or before terms grow. For nearly coincident
moderate coordinates, a three-point integral of `1-t*R(t)` computes the
difference. Other regions use the better normal/log-tail identity;
small complementary log probabilities use log1p. Near the center,
scaled two-component sums and a binary-split product retain the signed
normal coordinate when log(z)-log(mu_shape) rounds to zero.

Narrow intervals use three-point Gauss-Legendre integration in physical
coordinates when a bound on log-density variation is below 0.001.
For inverse Gaussian the bound is based on
`abs(d log(f)/d log(z)) <= 1.5+(z/mu_shape**2+1/z)/2`;
for the reciprocal it is `0.5+(z+1/(mu_shape**2*z))/2`.
Other intervals subtract the better-conditioned logarithmic tail;
nearly full intervals complement both omitted tails.

For extremely steep monotone intervals narrower than 1e-6 times the lower
support distance, an endpoint Laplace integral avoids coincident rounded
log tails. It is used only when one signed density-slope term dominates
the other by more than 1e16, the absolute log-coordinate slope exceeds
1e16, and slope times interval width exceeds 1e8. The neglected relative
curvature correction is then of order the reciprocal log-coordinate slope.
The integral retains the endpoint log density, logarithmic physical slope,
and exponential interval-mass factor without forming their probabilities.

These likelihood callbacks accept every positive finite binary64 shape,
including subnormals, without widening the scalar distribution API's
existing numerical guarantees. Finite mathematical log probabilities
beyond binary64 range return negative infinity; exponentiated masses
may underflow while their logarithms remain finite. Extremely concentrated
affine models can be ill-conditioned relative to input rounding, despite
the compensated center calculation. Bounds and fitted parameter
identifiability remain the caller's responsibility; optimizer success is
local convergence and does not establish a global optimum.
See the [validation record](../validation/inverse_gaussian_censored_2026-11-17.md).


## Generalized gamma and truncated normal/Weibull details

The new likelihoods and corresponding `*_censored_nnlf` functions are
exported by `scifort_censored` and `scifort_stats`. Existing fit generics
accept `censored_data` with bounds in these orders:

| Family | Parameter order | Current shape/bound limits |
| --- | --- | --- |
| `gengamma` | a, c, loc, scale | tiny<=a<=1e6; 1e-6<=abs(c)<=1e6 |
| `truncnorm` | a, b, loc, scale | finite -100<=a<b<=100 |
| `truncweibull_min` | c, a, b, loc, scale | 1e-6<=c<=1e6; finite 0<=a<b |

Here tiny means `tiny(1.0_dp)`. All parameters must be finite, scale must
be positive, and omitted loc/scale default to 0/1. Equal fit bounds fix
any parameter; all parameters can vary. Original array-based fitting
calls remain available. Invalid parameters return NaN even for empty
data; valid empty likelihoods are zero, and empty fits are invalid.
These shape and truncation limits are explicit compatibility gaps,
separate from named-family coverage.

For generalized gamma, `z=(x-loc)/scale` and `Y=z**c` has a unit-scale
Gamma(a) law. Positive c preserves CDF/SF order; negative c reverses it.
The density is evaluated through the existing log-gamma density and the
Jacobian `abs(c)/(x-loc)`, retaining log(Y) when Y underflows or overflows.
Tiny shapes reuse the existing E1 gamma-tail limit and retain small
complementary log probabilities with log1p. With c>0, endpoint density
at loc is infinite, finite, or zero as a*c is less than, equal to, or
greater than one. Negative c has zero endpoint density. Density/CDF
below loc are zero, and SF is one.

Adjacent generalized-gamma intervals use a three-point relative-density
integral when both signed power and gamma/hazard slope terms vary by less
than 0.001. Other intervals use logarithmic tail differences. Extremely
steep monotone narrow intervals use the existing endpoint-expansion
approach: Y must dominate a by 1e16, abs(c)*Y must exceed 1e16,
rate*width must exceed 1e8, and relative width must be below 1e-6.
The omitted relative curvature correction is of order 1/Y or 1/(abs(c)*Y).
The c=1, c=-1, and a=1 cases agree with gamma, inverse gamma, and Weibull.

Truncation bounds a,b are standardized parent coordinates, before loc
and scale. Both physical endpoints `loc+scale*a` and `loc+scale*b` must
be finite, distinct, and representable without a nonzero scaled bound
underflowing to zero. The current implementation conservatively requires
abs(loc) and abs(scale*a/b) to be at most huge/4. Endpoints are included
for density evaluation; CDF is zero/one at the lower/upper endpoint and
SF is one/zero. Finite intervals are intersected with the support.

If `b-a < 1e-4*max(1,abs(a),abs(b))`, the current contract additionally
requires loc=0 and an exact power-of-two scale. All truncation configurations
require nonzero scaled bounds and nonzero physical endpoints to be normal
binary64 values, and `abs(loc)/(scale*(b-a))<=1e4`. These restrictions keep
subnormal and strongly rounded affine endpoints from changing the support.
The additional narrow-window restriction keeps very narrow
standardized windows from being changed by rounded affine endpoints;
general translated/nonbinary-scaled narrow windows remain a documented gap.
Ordinary affine windows can still be ill-conditioned near their endpoints;
finite noncollapsed support alone does not establish parameter identifiability.

Truncated normal divides parent normal densities/interval masses by
`Phi(b)-Phi(a)`, using the existing better-tail and small-interval normal
integral rather than subtracting rounded probabilities. Truncated Weibull
uses the shifted hazard `H(z)-H(a)`, H(z)=z**c. Its density is
`c*z**(c-1)*exp(-(H(z)-H(a)))/(scale*(1-exp(-(H(b)-H(a)))))`.
Logarithmic power increments preserve adjacent boundaries and avoid
subtracting large rounded powers. For very small normalizers, factored
power ratios cancel the common logarithmic hazard before subtraction,
retaining conditional densities and interval masses even when H(b)
underflows. Complementing the smaller tail retains tiny log-CDF/log-SF.

Likelihoods accept the shared finite observation and infinite one-sided
interval conventions. Equal finite interval endpoints represent exact
observations. Log probabilities beyond binary64 range return negative
infinity; exponentiated probabilities may underflow while logarithms remain
finite. Fit success denotes local bounded convergence, not global optimality.
See the [validation record](../validation/power_trunc_censored_2026-11-17.md).


## Alpha, generalized extreme value, and Crystal Ball details

The likelihoods and their `*_censored_nnlf` partners are exported by
`scifort_censored` and `scifort_stats`. The corresponding fit generics
accept `censored_data`; existing real-array calls remain available.

| Family | Fit bounds order | Current numerical shape limits |
| --- | --- | --- |
| `alpha` | a, loc, scale | 1e-6<=a<=100 |
| `genextreme` | c, loc, scale | abs(c)<=100; c=0 or abs(c)>=1e-6 |
| `crystalball` | beta, m, loc, scale | 1e-6<=beta<=100; 1<m<=1e6 |

All parameters must be finite and scale positive. Omitted loc and scale
default to zero and one. Equal fit bounds fix any parameter; all parameters
may vary within these limits. Valid empty likelihoods are zero, invalid
parameters return NaN even for empty data, and empty fits are invalid.
These finite shape ranges are compatibility limitations of the new censored
methods, separate from the scalar distribution APIs and named-family count.
In particular, very small nonzero GEV shapes remain unsupported; c=0 is
evaluated by its exact Gumbel limit. Bounds spanning zero can contain this
unsupported region; likelihoods return NaN there and the optimizer rejects
such candidates.

Alpha has support x>loc, with zero endpoint density/CDF and unit SF.
Writing z=(x-loc)/scale, its CDF is Phi(a-1/z)/Phi(a). The SF and interval
probabilities are normal interval integrals, retaining reciprocal interval
widths in logarithmic physical coordinates. Small SFs therefore remain
available even when a-1/z rounds to a. Extremely steep negative normal
intervals use the upper CDF when the lower contribution is negligible:
the normal log-density variation must exceed 1e8.

GEV follows SciPy's shape sign convention: T=(1-c*z)**(1/c),
F=exp(-T), with T=exp(-z) at c=0. Its support satisfies 1-c*z>=0.
For positive c the finite upper endpoint density is zero for c<1,
1/scale for c=1, and infinite for c>1. The negative-c lower endpoint
density is zero. Tail masses and intervals use log(T) and physical-width
increments before exponential differencing, retaining underflowed survival
probabilities in logarithmic form. As in the scalar API, endpoint detection
uses binary64 standardized support arithmetic; evaluation extremely close
to a shape-dependent endpoint can be ill-conditioned relative to input
rounding. This does not provide additional precision for affine endpoints.

Crystal Ball joins a power-law left tail at z=-beta to a Gaussian core.
Its log normalizer is the log sum of the two integrated pieces. The smaller
normalization weight and its complement are evaluated separately to retain
tiny log probabilities. Tail intervals factor a positive power-law ratio;
Gaussian intervals reuse the existing stable normal integral. Intervals
crossing the transition add the two positive pieces in log space. Density
is continuous at the transition, and the support is the whole real line.
For steep Gaussian intervals above standardized 12 with normal log-density
variation exceeding 1e8, the upper contribution is negligible and the lower
SF is retained directly, including near binary64's finite log-probability limit.
For crossing intervals with standardized width times max(1,beta)<=1e-6,
the shared density and slope permit a local Gaussian extension over the
entire retained physical width. Its relative truncation error is bounded
by order (1+beta**2/m)*width**2, below 1e-12 in this regime. This avoids
changing an adjacent interval's total width by rounding its two pieces.

All three retain the shared exact/left/right/interval conventions, including
infinite one-sided interval endpoints and equal finite endpoints as exact
observations. Mathematical log probabilities beyond binary64 range return
negative infinity. A finite logarithm does not imply its exponentiated
probability is representable. Very concentrated affine models and finite
support endpoints remain subject to input rounding; fit success denotes
local bounded convergence, not global optimality or identifiability.

Definitions: [alpha](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.alpha.html),
[GEV](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.genextreme.html),
and [Crystal Ball](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.crystalball.html).
See the [validation record](../validation/shape_censored_2026-11-17.md).


## ARGUS, Pearson III, and kappa4 details

The likelihoods and corresponding `*_censored_nnlf` functions are exported
by `scifort_censored` and `scifort_stats`. The existing fit generics accept
`censored_data`, preserving their original real-array overloads.

| Family | Fit bounds order | Current numerical shape limits |
| --- | --- | --- |
| `argus` | chi, loc, scale | 1e-6<=chi<=100 |
| `pearson3` | skew, loc, scale | abs(skew)<1.6e-5, or 0.002<=abs(skew)<=2000 |
| `kappa4` | h, k, loc, scale | abs(h),abs(k)<=10; each zero or magnitude>=1e-6; if h>0, k*log(h)<=log(10) |

All parameters must be finite and scale positive. Omitted loc and scale
default to zero and one. Equal bounds fix any fit parameter; all parameters
may vary within the supported ranges. Valid empty likelihoods are zero,
invalid parameters return NaN even for empty data, and empty fits are invalid.
These shape restrictions are numerical compatibility gaps of the censored
methods; they do not narrow the scalar distribution APIs. Unsupported
candidate parameters in fitting are rejected by the optimizer. Bounds
crossing the documented gaps do not establish convergence across them.

ARGUS has standardized support 0<=z<=1, z=(x-loc)/scale, with zero
density at both endpoints. Set A=chi**2/2 and Y=A*(1-z**2).
Y has a Gamma(3/2) law conditioned on 0<Y<A, with reversed ordering.
The SF is P(3/2,Y)/P(3/2,A); its CDF and finite intervals use positive
gamma interval integrals, not subtraction of probabilities near one.
Logarithmic increments retain A*(zb+za)*(b-a)/scale even for adjacent
physical endpoints. Complementing the smaller tail retains small log-SF.

For Pearson III, beta=2/skew and alpha=beta**2. The positive coordinate
Y=beta*(z+beta) follows Gamma(alpha); negative skew reverses CDF/SF.
The support endpoint is z=-beta, below the support for positive skew
and above it for negative skew. Its density is infinite for alpha<1,
abs(beta)/scale for alpha=1, and zero for alpha>1. The gamma coordinate
and Jacobian are retained logarithmically, including when Y overflows.
Both the normal approximation threshold abs(skew)<1.6e-5 and its use
for small nonzero skew follow the existing scalar API. The gap between
that threshold and abs(skew)=0.002 remains unsupported in the censored
method; the nonnormal branch currently has 1e-6<=alpha<=1e6.

Gamma interval masses use a three-point relative-density integral when
the relative width and both gamma log-density terms vary by at most 0.001.
Otherwise they use the better logarithmic tail difference. When the lower
gamma coordinate exceeds max(1,alpha) by 1e16 and the width exceeds 1e8,
the upper contribution is negligible and the lower gamma SF is retained.
The Pearson normal branch similarly retains a dominant normal tail when
the log-density variation exceeds 1e8 beyond standardized magnitude 12.

Kappa4 uses U=1-k*z, T=U**(1/k), V=1-h*T and F=V**(1/h), with
continuous formulas for h=0 and k=0. Its standardized lower endpoint is
(1-h**(-k))/k when h>0 (log(h) for k=0), 1/k when h<=0 and k<0,
and negative infinity otherwise. The upper endpoint is 1/k for k>0
and positive infinity otherwise. Density at a finite upper endpoint is
zero, 1/scale, or infinite as k is below, equal to, or above one. At the h>0
lower endpoint the same classification uses h. For h<0,k<0, it uses
h*k: equality gives (-h)**(1/h-1)/scale.

Logarithmic reverse hazards and increments retain finite interval masses
without subtracting rounded CDFs. For negative h, the density factors
T/V before combining it with the hazard, avoiding overflow in cancelling
intermediate terms. The additional h>0 restriction is equivalent to
h**(-k)>=0.1; more compressed outer endpoints remain unsupported.
The h=0, (h=-1,k=0), (h=1,k=0), and h=k=1 cases agree with GEV,
logistic, exponential, and uniform, respectively.

Finite support detection uses binary64 standardized coordinates, as do
the existing scalar APIs. Evaluation extremely close to a shape-dependent
or affine endpoint can be ill-conditioned relative to input rounding;
these routines do not provide additional precision for endpoints or make
unrepresentable quantiles representable. Infinite one-sided intervals and
equal finite endpoints follow the shared data conventions. Log probabilities
beyond binary64 range return negative infinity. Fit success denotes local
bounded convergence, not a global optimum or parameter identifiability.

Definitions: [ARGUS](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.argus.html),
[Pearson III](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.pearson3.html),
and [kappa4](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.kappa4.html).
See the [validation record](../validation/gamma_kappa_censored_2026-11-17.md).


## Gauss hypergeometric, Irwin-Hall, and line von Mises details

The likelihoods and `*_censored_nnlf` functions are exported from
`scifort_censored` and `scifort_stats`. The existing fit generics accept
`censored_data`, retaining their original real-array overloads.

| Family | Fit bounds order | Current censored shape limits | Standardized support |
| --- | --- | --- | --- |
| `gausshyper` | a, b, c, z, loc, scale | 0.1<=a,b<=20; abs(c)<=10; -0.5<=z<=10 | [0,1] |
| `irwinhall` | n, loc, scale | integer 1<=n<=32 | [0,n] |
| `vonmises_line` | kappa, loc, scale | 0<=kappa<=100 | [-pi,pi] |

Here pi denotes the library's binary64 constant. These are numerical limits
of the censored methods, not new restrictions on scalar distribution APIs.
All parameters must be finite, scale must be at least `tiny(1.0_dp)`,
abs(loc)/scale must not exceed 1e4, and scale*extent must not exceed
huge/2, where extent is 1, n, or pi. Also require
abs(loc)+scale*extent<=huge, ensuring finite affine support endpoints. These bounds avoid unrepresentable supports and severe affine
cancellation. Omitted loc/scale default to zero/one. Invalid parameters
return NaN even for empty data; a valid empty likelihood is zero and an
empty fit is invalid. Equal bounds fix any parameter. Irwin-Hall fits
enforce integrality of n, including when real bounds contain several integers.
Fit success denotes local bounded convergence, not a global optimum.

Gauss hypergeometric uses the beta-weighted tilted density
u**(a-1)*(1-u)**(b-1)*(1+z*u)**(-c), normalized on [0,1]. Direct
logarithmic interval integrals reuse the existing double-exponential
endpoint transformation and preserve singular endpoint densities. At zero,
density is infinite/finite/zero according as a is below/equal to/above one;
the upper endpoint has the analogous classification using b.
The finite lower endpoint log density for a=1 is -lognorm-log(scale);
the finite upper endpoint for b=1 also includes -c*log1p(z).

`scifort_gausshyper` additionally exports
`gausshyper_logkernel_integral(lo,hi,a,b,c,z)`, a **SciFort extension**.
This is the unnormalized log integral in standard coordinates with
0<=lo<=hi<=1, a,b>0, finite c, and finite z>-1; equal endpoints return
negative infinity and invalid arguments return NaN. It reuses the scalar
API's existing quadrature and does not promise accuracy throughout the
unrestricted shape domain. The censored contract above states the newly
validated range. This helper is available directly from its family module.

Irwin-Hall evaluates its positive cardinal B-spline recurrence in the log
domain. The first and last polynomial pieces use direct powers in log
form, with analytic log integrals on those outer pieces, including
subnormal positive endpoints. Other intervals split at integer knots and
use a 16-point Gauss-Legendre rule, exact for polynomials through degree 31 in exact arithmetic. Thus
tiny probabilities can retain finite logarithms after ordinary probability
underflow. For n=1 both endpoint densities are 1/scale; for n>1 they
are zero. Its n=1 case is uniform and n=2 is triangular.

Line von Mises integrates exp(kappa*(cos(u)-1)), with the existing scaled
Bessel normalization. A composite 16-point Gauss-Legendre rule limits
each panel to width at most 0.25. Its support is a fixed real interval:
there is no wrapping or circular arc interpretation. Both endpoint
densities are finite. Kappa=0 gives uniform on [-pi,pi]. The log kernel
uses -2*kappa*sin(u/2)**2 to avoid cancellation at the center.

Tails integrate the requested support segment directly and complement the
smaller tail when appropriate. Narrow interior intervals retain
log(b-a)-log(scale), using physical endpoints before standardization can
collapse their width. They use GL16 when width is less than 1e-5 times
the distance to either endpoint. Wider intervals use the family-specific
integrator. Interval clipping and support classification use binary64
standard coordinates. If standardization underflows or rounds onto a
support endpoint, that endpoint behavior applies. Thresholds extremely
close to a rounded or affine endpoint remain ill-conditioned relative to
input rounding; extra endpoint precision is not provided. One-sided infinite intervals and equal finite
endpoints follow the common data conventions.

Definitions: [Gauss hypergeometric](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gausshyper.html),
[Irwin-Hall](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.irwinhall.html),
and [line von Mises](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.vonmises_line.html).
See the [validation record](../validation/bounded_final_censored_2026-11-18.md).


## Normal-exponential, skew-normal, and double Pareto-lognormal details

The likelihoods and corresponding `*_censored_nnlf` functions are exported
by `scifort_censored` and `scifort_stats`. Existing fit generics accept
`censored_data` and retain their real-array overloads.

| Family | Fit bounds order | Current censored shape limits | Standard support |
| --- | --- | --- | --- |
| `exponnorm` | k, loc, scale | 0.01<=k<=100 | real line |
| `skewnorm` | a, loc, scale | abs(a)<=10 | real line |
| `dpareto_lognorm` | u, s, a, b, loc, scale | abs(u)<=100; 0.01<=s<=10; 0.01<=a,b<=100 | positive half-line |

All parameters must be finite and scale positive, with omitted loc/scale
defaulting to zero/one. These numerical limits apply to the new censored
methods, not the existing scalar APIs. Invalid parameters return NaN even
for empty data; a valid empty likelihood is zero and an empty fit invalid.
Equal bounds fix any parameter. Fit success denotes bounded local
convergence, not a global optimum. Double Pareto-lognormal u and outer
scale are confounded through u+log(scale); fixing one is advisable when
interpreting parameters. All-free candidates remain supported.

For the exponentially modified normal, r=1/k and
H(z,r)=exp(r*r/2-r*z)*Phi(z-r), z=(x-loc)/scale. The standard PDF is
r*H, CDF is Phi(z)-H, and SF is Phi(-z)+H. In the Gaussian tail,
H=phi(z)*R(r-z), where R is the normal Mills ratio. When z<=-12,
the CDF uses phi(z)*(R(-z)-R(r-z)), retaining r independently in the
existing factored Mills-ratio difference. This avoids cancellation of
large Gaussian exponents and nearly equal Mills ratios. The exponential
right tail is evaluated directly, with guarded products for log overflow.
CDF/SF values near one complement the smaller logarithmic tail, retaining
representable tiny log probabilities rather than rounding them to zero.

Skew-normal uses 2*phi(z)*Phi(a*z). Its a=0 case is normal and
F(0)=1/2-atan(a)/pi. For z<=0 and a>=0, direct relative-density
integration uses t=z-v/rate with rate=1+a+(1+a*a)*abs(z). A composite
GL16 rule uses 24 panels on 0<=v<=48. The omitted integral relative to
the endpoint density is bounded by exp(-47); the normalized omitted
fraction is below 4*exp(-47). Negative shapes use
F(z,a)=2*Phi(z)-F(z,-a), and upper tails use reflection. Normal Mills
ratios retain relative density changes when separate Gaussian log densities
would round equal. Logarithms remain finite after ordinary tail underflow.

For double Pareto-lognormal, let y=(x-loc)/scale and
z=(log(y)-u)/s. Its log coordinate is a positive mixture of an EMG with
rate a*s and weight b/(a+b), and a reflected EMG with rate b*s and
weight a/(a+b). CDF, SF, density, and interval masses use these positive
mixtures rather than signed probability sums. When the requested tail
exceeds one half, its log uses the complementary smaller mixture. This
retains tiny negative log-CDF/log-SF values. The physical distance
log(x-loc) is retained before scaling, including when x-loc or y would
overflow or y would underflow. As in the existing scalar API, exact
observations at x=loc have zero density; CDF/SF there are zero/one.

Narrow EMG/skew-normal intervals use three-point Gauss-Legendre integration
of log densities with the physical width retained before standardization.
The small-window criterion bounds width times a conservative density-slope
scale by 0.001. Other intervals use the better logarithmic tail difference.
EMG intervals beginning above r+40 use the exponential hazard increment
directly; omitted Gaussian corrections are below ordinary binary64 precision.
Double Pareto-lognormal intervals retain log(log(1+(b-a)/(a-loc)))-log(s)
before separately integrating the two normal-exponential components.

These functions do not provide additional precision for input coordinates.
For EMG/skew-normal, a standard coordinate that overflows is treated as
infinite; a finite log tail can be lost if scaling overflows before a small
exponential rate is applied. Gaussian logs outside binary64 range return
negative infinity. Very steep intervals can be ill-conditioned relative
to input rounding; scaled log-error tests are not uniform relative
probability guarantees. Wider shape ranges remain numerical coverage gaps.
Equal finite and infinite one-sided intervals follow the shared conventions.

Definitions: [EMG](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.exponnorm.html),
[skew-normal](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.skewnorm.html),
and [double Pareto-lognormal](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.dpareto_lognorm.html).
See the [validation record](../validation/normal_mix_censored_2026-11-18.md).

## Rice and noncentral chi-squared/F censoring

```fortran
value = rice_censored_loglikelihood(data,b,loc,scale)
value = rice_censored_nnlf(data,b,loc,scale)
value = ncx2_censored_loglikelihood(data,df,nc,loc,scale)
value = ncx2_censored_nnlf(data,df,nc,loc,scale)
value = ncf_censored_loglikelihood(data,dfn,dfd,nc,loc,scale)
value = ncf_censored_nnlf(data,dfn,dfd,nc,loc,scale)
call rice_fit(data,lower,upper,result,guess,max_iter,tolerance)
call ncx2_fit(data,lower,upper,result,guess,max_iter,tolerance)
call ncf_fit(data,lower,upper,result,guess,max_iter,tolerance)
```

All likelihoods use the existing exact/left/right/interval observation
conventions, with support x>=loc. Optional loc defaults to zero and scale
to one; both must be finite, with scale positive. Parameter order for
fit bounds and guesses is [b,loc,scale], [df,nc,loc,scale], or
[dfn,dfd,nc,loc,scale]. Equal bounds fix parameters. Existing real-array
fit overloads remain available; no C ABI additions are included.

Validated parameter contracts are 0<=b<=10 for Rice,
0.1<=df<=2000 and 0<=nc<=100 for noncentral chi-squared, and
0.1<=dfn,dfd<=2000 and 0<=nc<=100 for noncentral F.
Noncentral F additionally requires scale*dfd/dfn to be normal and finite.
Invalid parameters return NaN even with empty data; valid empty likelihoods
return zero. Empty fits are invalid. Noncentrality zero is evaluated by
the corresponding central gamma/beta-prime law, without a Poisson sum.

Positive mixtures retain log densities, both log tails, and interval
masses without subtracting rounded ordinary probabilities. Rice maps
z=(x-loc)/scale to a noncentral chi-squared variable z*z with df=2 and
nc=b*b; its density includes the transformation Jacobian. Transformations
use logarithms of physical distances and widths, retaining positive
subnormal inputs and adjacent intervals even when direct scaling or
squaring would underflow or overflow. A tail larger than one half is
formed by complementing the directly evaluated smaller mixture tail.
Large interval probabilities similarly complement the two omitted tails.

The Poisson sum stops when a decreasing bound on successive component
ratios bounds the remaining mass by 2e-16 times the accumulated mass.
At most 10,001 components (indices 0 through 10000) are evaluated. Failure
to meet this bound returns NaN rather than a partial sum. Very remote
noncentral chi-squared/Rice upper tails and densities can reach this
limit, even when their log probability is representable; callers and fit
bounds must account for such numerical failures. Centrality zero does
not have this series restriction. The remainder bound controls truncation,
not all binary64 rounding or pre-existing special-function error.

At x=loc, Rice density is zero. Noncentral chi-squared density is infinite
for df<2, exp(-nc/2)/(2*scale) for df=2, and zero for df>2.
Noncentral F has the analogous cases according to dfn, with density
exp(-nc/2)/scale for dfn=2. Left and right probabilities at the lower
endpoint are zero and one. Log results outside binary64's range become
negative infinity. No claim of uniform relative probability accuracy is
made for extreme logs: the validation measures scaled log errors and
separately checks tiny nonzero complementary logs with relative error.

See the [validation record](../validation/noncentral_censored_2026-11-18.md)
for checked/optimized native tests, independent high-precision references,
fit benchmarks, exact commands, and remaining numerical limitations.

## Kolmogorov, wrapped Cauchy, and relativistic Breit-Wigner censoring

```fortran
value = kstwobign_censored_loglikelihood(data,loc,scale)
value = kstwobign_censored_nnlf(data,loc,scale)
value = wrapcauchy_censored_loglikelihood(data,c,loc,scale)
value = wrapcauchy_censored_nnlf(data,c,loc,scale)
value = rel_breitwigner_censored_loglikelihood(data,rho,loc,scale)
value = rel_breitwigner_censored_nnlf(data,rho,loc,scale)
call kstwobign_fit(data,lower,upper,result,guess,max_iter,tolerance)
call wrapcauchy_fit(data,lower,upper,result,guess,max_iter,tolerance)
call rel_breitwigner_fit(data,lower,upper,result,guess,max_iter,tolerance)
```

Each family supports exact, left/right, and interval censoring through
`censored_data`. The location and scale are optional, defaulting to zero
and one. Fit bounds and guesses use [loc,scale], [c,loc,scale], or
[rho,loc,scale]. Equal bounds fix parameters. Existing real-array fits
remain available. Valid empty likelihoods are zero; invalid parameters
return NaN even for empty data, and empty fits are invalid.
No C ABI or wrapper additions are included.

Kolmogorov (`kstwobign`) has no shape parameter and support x>=loc, with
zero density at loc. Location must be finite and scale positive/finite.
Its two theta-series representations are evaluated after factoring their
leading exponentials, including the differentiated density series.
This retains log tails even when ordinary probabilities underflow.
Small-coordinate and large-coordinate interval formulas retain physical
endpoint increments before exponentiation; adjacent observations do not
require distinguishable standardized coordinates. Log probabilities
beyond the finite binary64 range become negative infinity.

Wrapped Cauchy uses SciPy's fixed standardized interval [0,2*pi], with
physical upper=loc+(2*pi)*scale computed in binary64. This method requires
0<c<=0.99, normal positive scale, abs(loc)/scale<=10,
scale<=huge/(4*pi), and a finite noncollapsed upper endpoint.
Observations outside this interval have zero density; endpoint densities
are finite. CDF/SF endpoint probabilities are zero/one and one/zero.
Finite censoring intervals are clipped to the support. An interval
crossing the angular boundary must not be supplied with its endpoints
reversed: that remains invalid input. A future crossing-arc method would
be a **SciFort extension** with a distinct observation convention.
There is no periodic reduction of exact data or censoring bounds.

Wrapped-Cauchy lower and upper tails use half-angle arctangent formulas,
with independently retained physical distances to both support endpoints.
Near the upper endpoint, the method reflects about the computed finite
support, rather than subtracting rounded standardized angles. This
binary64 endpoint convention, also used for endpoint classification,
can differ from an ideal real-arithmetic affine endpoint by rounding.
The c and location/scale limits bound sensitivity; broader concentration
and affine ranges remain numerical coverage gaps. Narrow interval masses
use physical widths and three-point Gaussian quadrature.

Relativistic Breit-Wigner requires 0.1<=rho<=100, finite loc, and positive
finite scale. Its support is x>=loc, with finite density at loc. The
ordinary tails call the existing family module's antiderivative, retaining
its documented provenance. The far upper tail instead integrates a
convergent reciprocal-coordinate expansion; tiny lower tails use the
constant endpoint density. Densities use log-domain quartic denominators,
and far-tail intervals integrate each expansion term with its retained
log-relative physical width. This avoids cancellation in the complex
antiderivative at extreme coordinates. For physical resonance modelling,
fix loc=0 as recommended in SciPy's documentation; the mathematical fit
API permits other locations within the stated finite-parameter contract.

Ordinary narrow intervals use three-point Gaussian quadrature when a
conservative density-variation criterion is small. Quadrature coordinates
are formed from log distances and physical width increments, preserving
narrow interval masses even when adding a node to a large location would
round it onto an endpoint. Wider intervals use the better log-tail difference, complementing omitted tails for nearly
whole-support probabilities. Scaled log errors and tiny complementary-log
relative errors are validated separately. Very steep tails remain
sensitive to rounding of physical observations; no uniform relative
probability accuracy is claimed after exponentiation.
See the [validation record](../validation/three_tail_censored_2026-11-18.md).


## Generalized inverse Gaussian and hyperbolic censoring

`geninvgauss`, `norminvgauss`, and `genhyperbolic` support mixed exact,
left/right/interval censored likelihoods, matching `_censored_nnlf` functions,
and bounded `_fit(censored_data, lower, upper, result, ...)` overloads.
All are exported by `scifort_stats`; existing real-array fit calls remain.
Parameter/bound order is respectively `(p,b,loc,scale)`, `(a,b,loc,scale)`,
and `(p,a,b,loc,scale)`. Equality of lower/upper bounds fixes a parameter.

The current censored methods require `-5<=p<=5`, `0.5<=a<=10`, and
`abs(b)<=0.8*a` for generalized hyperbolic; normal-inverse Gaussian uses
that same tail/skew range and the identity `p=-0.5`. Generalized inverse
Gaussian requires `-5<=p<=5` and `0.5<=b<=10`. Location is finite and scale
is positive finite. Values outside these shape ranges return NaN, including
for empty data. Extending these ranges and the generalized-hyperbolic
`abs(b)=a, p<0` boundary remains future work. These restrictions apply to
the censored methods, not the existing scalar distribution APIs.

The GIG standardized support is `(0,infinity)`; observations at/below the
location have zero density/CDF and unit survival. Hyperbolic families have
real-line support. Valid empty likelihoods equal zero; empty fits are invalid.
Equal finite interval endpoints denote exact observations. Infinite interval
endpoints use the corresponding one-sided tail; two infinite endpoints and
reversed intervals are invalid, consistent with the common data contract.

Positive log-density quadrature uses `y=log(z)` for GIG and `y=asinh(z)`
for the hyperbolic families. The smaller probability determines its
complement through stable log identities. Narrow intervals integrate the
density while retaining physical widths separately from location/scale.
The normal-inverse Gaussian implementation shares the generalized-hyperbolic
kernel at `p=-0.5` rather than computing two different approximations.

For `z>1e8` and `0<z<1e-8`, GIG uses leading endpoint Laplace tails;
for `abs(z)>1e8`, hyperbolic tails use the leading exponential endpoint
formula. These preserve finite logarithms after probability underflow.
Their error is assessed relative to the large magnitude of the log tail;
they do not promise correctly rounded probabilities throughout those remote
ranges. When the true negative log probability exceeds representable
binary64 magnitude, the result is negative infinity. Tiny complementary
log probabilities can round to zero when the probability is below the
smallest representable positive number. For the hyperbolic families, a
standardized coordinate that itself overflows uses the corresponding
infinite-coordinate limit. Retaining finite log tails when `(x-loc)/scale`
is unrepresentable remains an additional numerical gap.

See [Bessel-family numerical validation](../validation/bessel_censored_2026-11-18.md)
for the independent references, commands, and measured accuracy.


## Tukey lambda censoring

`tukeylambda_censored_loglikelihood`, its `_censored_nnlf`, and the
`tukeylambda_fit(censored_data, lower, upper, result, ...)` overload use
`(lam,loc,scale)` parameter order. They are exported by `scifort_stats`.
The current censored methods require `abs(lam)<=5` and either `lam=0` or
`abs(lam)>=1e-6`. Location is finite and scale is positive finite.
For positive lambda, scale must be normal, `abs(loc)/scale<=10`, and
`extent=scale/lam<=huge/4`; both finite endpoints `loc-/+extent` must be
representable and distinct from the location. Unsupported parameters return
NaN, even for empty data. These bounds apply only to the censored methods.

Negative/zero lambda has real-line support. Positive lambda has the finite
computed affine support `[loc-scale/lam,loc+scale/lam]`. Tail evaluation near
these endpoints retains the physical distance to the computed endpoint.
Central inversion instead uses the quantile itself, avoiding subtraction of
a very large support radius when positive lambda is small. Inversion solves
for log probability, retaining extremely small probabilities and densities.
The lambda-zero case uses the logistic law directly. Some SciPy versions
approximate sufficiently small nonzero lambda as zero; the new methods use
the nonzero-lambda quantile formula over the documented envelope.

At finite endpoints the censored density uses its mathematical one-sided
limit: zero for `0<lam<1`, `1/(2*scale)` for `lam=1`, and `1/scale` for
`lam>1`. This intentionally differs from the current scalar SciFort
`tukeylambda_logpdf`, which returns negative infinity when its computed CDF
is an endpoint probability. Existing scalar APIs are unchanged. In
particular, lambda 1 is uniform on the affine `[-1,1]` branch and lambda 2
is uniform on the affine `[-0.5,0.5]` branch.

## One-sided finite-sample KS censoring

`ksone_censored_loglikelihood`, `_censored_nnlf`, and the censored `ksone_fit`
overload use `(n,loc,scale)` order. The current methods require integer
`1<=n<=100`, normal positive `scale<=huge/4`, `abs(loc)/scale<=10`, and a
finite noncollapsed upper endpoint `loc+scale`. The fit treats `n` as an
integer parameter; equality of bounds can fix it. Existing real-array fit
calls remain available. Extending the sample-size range remains future work.

The finite computed support is `[loc,loc+scale]`. Positive Smirnov survival
sums provide ordinary tails; direct endpoint polynomials and log-distance
formulas retain tiny complementary probabilities and intervals after
probability underflow. Narrow intervals use positive density quadrature,
split at the first density knot `loc+scale/n`. At that exact computed knot,
exact observations use the midpoint of the one-sided densities, consistent
with SciFort's existing knot convention. Intervals use the appropriate
one-sided density, rather than the midpoint. Some SciPy values at knots use
a different representative density; this does not change interval probability.

For `n=1` this is uniform on the affine `[0,1]` support, including finite
endpoint density. For `n>1`, exact observations at either support endpoint
retain SciFort's existing zero-density convention. Values outside support
have zero density; CDF/SF endpoint limits are zero/one in the usual order.

## Fixed-branch von Mises censoring: SciFort extension

The **censored method and censored fit overload are a SciFort extension**:
`vonmises_censored_loglikelihood`, `_censored_nnlf`, and
`vonmises_fit(censored_data, lower, upper, result, ...)` interpret all
observations on the fixed affine branch
`[loc-pi*scale,loc+pi*scale]`. Parameter order is `(kappa,loc,scale)`.
They delegate to the existing `vonmises_line` censored implementation and
inherit its numerical limits: `0<=kappa<=100`, normal positive scale,
`abs(loc)/scale<=1e4`, and representable noncollapsed finite support with
`pi*scale<=huge/2`. Existing real-array circular fits are unchanged.

This method does not evaluate SciPy's periodic, unwrapped circular CDF as
an ordinary probability. Angles are not reduced modulo `2*pi`; reversed
intervals are invalid, including intervals intended to cross the branch
boundary. Censoring of crossing arcs needs a separate representation and
remains future work. `scale` retains SciFort's affine-scale convention;
SciPy documents different scale handling for the circular distribution.

All three new methods retain the common data contract: valid empty
likelihoods equal zero; empty fits are invalid; equal finite interval
endpoints are exact observations; one infinite interval endpoint uses a
one-sided probability; two infinite endpoints and NaNs are invalid.
See [numerical validation](../validation/quantile_ks_censored_2026-11-18.md).

## Noncentral Student t and Landau censoring

`nct_censored_loglikelihood`, `nct_censored_nnlf`, and the censored
`nct_fit` overload use `[df,nc,loc,scale]` bounds order. The censored method
currently accepts **0.5<=df<=100** and **abs(nc)<=10**, with finite location
and positive finite scale. Infinite df and shapes outside these limits
return NaN, including for empty likelihood data; fit bounds must remain
within the supported domain. These restrictions are specific to the
censored method, and do not change the existing scalar distribution API.

The noncentral-t likelihood integrates the positive normal/chi mixture on
`log(sqrt(V/df))`, using a composite 16-point Gauss-Legendre rule with
maximum panel width 0.25. The lower integration bound is shifted with
`log(abs(z))`; a lower cutoff of `(100+nc**2)/df` and an upper normal-numerator
cutoff of abs(nc)+40 keep discarded contributions below the tested log
error target in this shape domain. Complementary tails are evaluated from
the smaller tail. This avoids cancellation and a fixed gamma-probability
floor in remote tails. Log absolute standardized distances retain finite
power-law log densities and tails even when direct division overflows.

`landau_censored_loglikelihood`, `landau_censored_nnlf`, and censored
`landau_fit` use `[loc,scale]` bounds order, with finite location and
positive finite scale. The parameterization is the same as SciFort's
existing scalar Landau API and SciPy's location-scale interface; there is
no extra scale-dependent location shift.

For standardized `z<=-1`, the likelihood independently inverts the Landau
Laplace transform on a vertical contour through its real saddle. Scaling
out the saddle exponent avoids both tiny-probability underflow and the
existing scalar left-tail asymptotic switch error near -2.7. At saddle
magnitude greater than 1e12, it uses the leading saddle terms; the omitted
relative expansion corrections are O(1/saddle). Central and right tails
reuse the existing Landau kernels. At positive standardized magnitudes
at least 1e100, the leading `2/(pi*z)` survival law and `2/(pi*z**2)` density
are evaluated through log absolute distances. A log probability below
binary64's finite range returns negative infinity; its complementary log
probability is zero after underflow.

For both families, intervals use well-separated log-tail differences or a
positive 16-point density integral for nearby endpoints. The latter uses
standardized offsets and the physical interval's log width, so adjacent
and subnormal widths survive physical-location rounding. At standardized
magnitudes at least 1e100, same-side intervals instead use the leading
power-law tail with physical distance ratios: exponent df for noncentral t,
and exponent 1 for Landau's right tail. This retains finite log masses and
relative interval widths even when one or both standardized endpoints
overflow. Corrections to those limiting laws are negligible at that cutoff
within the supported shape domain. All four observation types may be mixed
in one dataset. Valid empty likelihoods
are zero; empty fits and invalid data are rejected by the existing shared
contract. Existing real-array fit overloads remain available.

The tested target is `3e-11*max(1,abs(reference_log_value))`, with separate
relative checks for tiny complementary logs. For very large negative log
values, this does not promise small relative error in the underlying
probability or last-bit monotonicity. The new censored noncentral-t
integration and Landau contour are deliberately independent of the
existing scalar quadratures, so their last bits and difficult tails may
differ. See [validation](../validation/specialist_censored_2026-11-18.md)
for independent references, commands, and results. No C ABI or wrapper
entry points are added in this batch.


## Finite-sample two-sided Kolmogorov-Smirnov censoring

`kstwo_censored_loglikelihood(data,n,loc,scale)` and
`kstwo_censored_nnlf(data,n,loc,scale)` accept mixed exact, left-, right-, and
interval-censored observations. The `kstwo_fit` censored-data overload uses
parameter order `[n,loc,scale]`, with `n` constrained to integers. Existing
real-array fit calls remain available. These procedures are exported from
`scifort_censored` / `scifort_fit` and `scifort_stats`.

The current censored method requires integer `1<=n<=32`, a normal positive
scale no larger than `huge/2`, `abs(loc)/scale<=100`, and a finite noncollapsed
affine support. Unsupported parameters return NaN, including for empty data;
a valid empty likelihood is zero and an empty fit is invalid. Wider sample
sizes remain a numerical coverage gap even though the family is now covered.

The computed support is `[loc+scale/(2*n),loc+scale]`. Exact observations at
both endpoints have zero density. At the computed first knot `loc+scale/n`,
point observations use the left density, consistent with the existing scalar
small-statistic branch. For n>1 the right density is smaller by `1-1/n`.
Intervals crossing that knot integrate both sides; they do not use the point
value as a constant approximation.

An independent positive ordered-uniform rectangle-volume recurrence evaluates
the CDF and sums disjoint first-failure volumes for the SF. A positive boundary
flux evaluates the analytic density, avoiding finite differences. Endpoint
polynomials are evaluated in logs using distances from the computed physical
endpoints. Narrow intervals use endpoint power differences or positive GL16
quadrature split at every polynomial knot; each density piece has degree at
most 31. This retains interval logs even for adjacent binary64 endpoints.

Rounding an affine endpoint or knot can shift its mathematical position by a
few ulps. The computed endpoint conventions above deliberately take precedence;
last-bit monotonicity is not promised across those switches. The tested target
is `3e-11*max(1,abs(reference_log_value))`; it is a log-value target rather than
a guarantee of relative accuracy after exponentiation of an extreme log.
See the [validation record](../validation/kstwo_censored_2026-11-18.md).
C ABI and wrapper additions remain deferred.


## Levy-stable and studentized-range censoring

The two remaining named continuous families now expose censored log likelihoods,
NNLFs, and `censored_data` bounded-fit overloads. They accept all four observation
kinds and are exported from `scifort_censored` / `scifort_fit` and `scifort_stats`.
Fit parameter orders are `[alpha,beta,loc,scale]` and `[k,df,loc,scale]`.
Existing real-array fit overloads remain available. Valid empty likelihoods
are zero; empty fits and invalid data are rejected by the shared contract.
C ABI and wrapper work remains deferred.

Both new methods require normal positive scale<=huge/100 and
abs(loc)/scale<=10000. Unsupported parameters or observation coordinates return
NaN; a fit must stay within the supported domain. The named-family count is
now 109/109, but this is not complete numerical or SciPy compatibility.

### Levy-stable numerical domain

The stable method follows the default S1 parameterization. It supports:

- alpha=2, any beta in [-1,1]: the existing normal censored method with
  standard deviation sqrt(2)*scale;
- alpha=1, beta=0: the existing Cauchy censored method;
- alpha=0.5, beta=1 or -1: the existing Levy or reflected-Levy censored method;
- 1.25<=alpha<=1.75 and abs(beta)<=0.5, including remote tails and finite
  observations whose standardized distance overflows binary64.

The first three cases inherit the delegated family's tail and interval
contracts. The fourth case independently inverts the characteristic function
in the central region. The substitution t=v**4 smooths the origin, and a
composite GL16 rule integrates through t**alpha=48. Density, CDF, and reflected
CDF for the SF are evaluated separately through abs(z)<=6. Beyond that join,
a rotated characteristic-function contour computes the small tail directly.
Analytic subtraction and factoring of the leading power retain density and
tail logarithms when ordinary probabilities underflow. Reflection supplies
the left tail; stable logarithmic complementation supplies the large tail.
Narrow intervals integrate the density with physical widths, retaining logs
for adjacent binary64 endpoints. See the
[tail validation and derivation](../validation/stable_censored_tails_2026-11-19.md).

General skew alpha=1 and other stability/skewness ranges remain gaps. In
particular, no alpha=1 skew location-scale correction is
approximated: that branch is unsupported unless beta=0. The new independent
kernel does not inherit the existing scalar routine's near-one/near-mode input
rounding. Small differences from that scalar routine can therefore be expected.
The tested target is 3e-10*max(1,abs(reference_log_value)).

### Studentized-range numerical domain

The studentized-range method supports real 2<=k<=100 and either 1<=df<=100 or
df>=100000. Values at and above 100000, including positive infinity, use the
infinite-df normal-range approximation, matching the existing asymptotic switch.
The gap 100<df<100000 and 1<k<2 remain unsupported. Finite df=100
still uses the chi-scale mixture; it does not switch to the normal-range limit.
Fit bounds must be finite;
fix df=100000 to fit the infinite-df approximation.

For finite df, positive log-domain normal-range/chi-scale mixture integrals
retain small CDFs, survival probabilities, and densities. The SF integrates
first-failure probabilities of the normal range directly. A range-coordinate
change retains finite power-law log tails even when the standardized physical
coordinate would overflow. For the infinite-df approximation, positive
coordinates are limited to (x-loc)/scale<=40; larger coordinates return NaN.
No such coordinate cap applies to the supported finite-df branch.
For k>10 the inner normal-range integration uses half-width panels to resolve
its narrower peak. The chi-scale upper limit includes the (k-1)th moment
that controls the lower-end power law. Survival integration uses log1p when
the normal-tail ratio is small, avoiding cancellation in the hazard increment.

Support is x>=loc. At loc the CDF is zero and SF is one. Density is positive
for k=2 and zero for k>2; below loc it is zero. The k=2 law is exactly
sqrt(2)*abs(Student_t(df)), or sqrt(2)*abs(Normal(0,1)) in the infinite-df
limit. Narrow positive intervals use density integration; intervals beginning
at loc use the direct CDF to preserve the lower-end power behavior.
The tested log target is 3e-10*max(1,abs(reference_log_value)) for independent
references. Ordinary SciPy comparisons allow 3e-8 because its absolute-tolerance
integration can lose relative accuracy in small probabilities. This target does
not promise relative probability accuracy after exponentiating extremely large
negative logs, or last-bit monotonicity across numerical switches.

See the [validation record](../validation/final_censored_2026-11-18.md) for
independent Fourier, folded-t, endpoint-power, adaptive integral, and fit checks.

Validation of the extension through df=100 is recorded in the
[finite-df validation record](../validation/studentized_range_df100_2026-11-18.md).

Validation of the extension through k=100, including noninteger k, is recorded
in the [larger-k validation record](../validation/studentized_range_k100_2026-11-19.md).
