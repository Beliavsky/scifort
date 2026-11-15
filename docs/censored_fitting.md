# Censored likelihoods and distribution fitting

SciFort supports exact, left-censored, right-censored, and interval-censored
likelihoods and bounded maximum-likelihood fits for normal, exponential,
minimum/maximum Weibull, Rayleigh, lognormal, Gibrat, gamma, Erlang, beta, Cauchy, symmetric/asymmetric Laplace,
logistic, half-Cauchy, half-logistic, half-normal, hyperbolic-secant, and both Gumbel distributions.
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
normal_censored_loglikelihood(data [, loc, scale])
rayleigh_censored_loglikelihood(data [, loc, scale])
cauchy_censored_loglikelihood(data [, loc, scale])
halfcauchy_censored_loglikelihood(data [, loc, scale])
halfnorm_censored_loglikelihood(data [, loc, scale])
hypsecant_censored_loglikelihood(data [, loc, scale])
halflogistic_censored_loglikelihood(data [, loc, scale])
laplace_censored_loglikelihood(data [, loc, scale])
laplace_asymmetric_censored_loglikelihood(data, kappa [, loc, scale])
logistic_censored_loglikelihood(data [, loc, scale])
gumbel_r_censored_loglikelihood(data [, loc, scale])
gumbel_l_censored_loglikelihood(data [, loc, scale])
exponential_censored_loglikelihood(data [, loc, scale])
weibull_censored_loglikelihood(data, c [, loc, scale])
weibull_max_censored_loglikelihood(data, c [, loc, scale])
lognormal_censored_loglikelihood(data, s [, loc, scale])
gibrat_censored_loglikelihood(data [, loc, scale])
gamma_censored_loglikelihood(data, a [, loc, scale])
erlang_censored_loglikelihood(data, a [, loc, scale])
beta_censored_loglikelihood(data, a, b [, loc, scale])
```

The corresponding `*_censored_nnlf` procedures return the negative
log-likelihood. Location defaults to 0 and scale to 1. Weibull shape `c`,
lognormal shape `s`, gamma/Erlang shape `a`, asymmetric Laplace shape `kappa`,
and beta shapes `a` and `b` must be supplied. Scale and shape must be positive and finite; location must be
finite. Invalid data or parameters return quiet NaN, including for empty data.
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

The existing `normal_fit`, `exponential_fit`, `weibull_fit`, `weibull_max_fit`,
`rayleigh_fit`, `lognormal_fit`, `gibrat_fit`, `gamma_fit`, `erlang_fit`,
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
Gibrat, exponential, Rayleigh, Cauchy, symmetric Laplace, logistic, half-Cauchy, half-logistic, half-normal, hyperbolic-secant, and both Gumbel families, `[c, loc, scale]` for both Weibull families, and `[s, loc, scale]` for
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
