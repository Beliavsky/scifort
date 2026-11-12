# Censored likelihoods and distribution fitting

SciFort supports exact, left-censored, right-censored, and interval-censored
likelihoods and bounded maximum-likelihood fits for normal, exponential,
minimum-Weibull, lognormal, gamma, beta, logistic, and both Gumbel distributions.
Lognormal and gamma censored fits currently require fixed location. Beta censored fits require fixed location
and scale. The definitions follow SciPy's
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
logistic_censored_loglikelihood(data [, loc, scale])
gumbel_r_censored_loglikelihood(data [, loc, scale])
gumbel_l_censored_loglikelihood(data [, loc, scale])
exponential_censored_loglikelihood(data [, loc, scale])
weibull_censored_loglikelihood(data, c [, loc, scale])
lognormal_censored_loglikelihood(data, s [, loc, scale])
gamma_censored_loglikelihood(data, a [, loc, scale])
beta_censored_loglikelihood(data, a, b [, loc, scale])
```

The corresponding `*_censored_nnlf` procedures return the negative
log-likelihood. Location defaults to 0 and scale to 1. Weibull shape `c`,
lognormal shape `s`, gamma shape `a`, or beta shapes `a` and `b` must be supplied. Scale and shape must be positive and finite; location must be
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

The existing `normal_fit`, `exponential_fit`, `weibull_fit`, `lognormal_fit`, `gamma_fit`,
`beta_fit`, `logistic_fit`, `gumbel_r_fit`, and `gumbel_l_fit` generics accept
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
exponential, logistic, and both Gumbel families, `[c, loc, scale]` for Weibull, and `[s, loc, scale]` for
lognormal, `[a, loc, scale]` for gamma, and `[a, b, loc, scale]` for beta. Equal bounds fix a parameter.
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

`gamma_censored_loglikelihood(data, a [, loc, scale])` uses positive shape
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
