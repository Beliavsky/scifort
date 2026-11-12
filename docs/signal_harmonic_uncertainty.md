# Harmonic-fit coefficient covariance and mean-curve uncertainty

This is a **SciFort extension**, not an existing `scipy.signal` API. It builds on
[`lombscargle_harmonic_fit`](signal_harmonic_fit.md), which fits weighted Fourier
series to irregular observations. Specifying `estimate_covariance=.true.`
computes an estimated covariance matrix for every successful trial frequency,
standard errors of the fitted coefficients, and the standard error of the
**fitted mean curve** at new times.

## Fortran example

```fortran
use scifort_signal, only : dp, harmonic_fit_result, lombscargle_harmonic_fit, &
    lombscargle_harmonic_predict, lombscargle_harmonic_predict_std, &
    irregular_status_success

real(dp) :: t(100), y(100), w(100), omega(2), future(20)
real(dp), allocatable :: forecast(:), forecast_se(:)
type(harmonic_fit_result) :: fit
integer :: status

! Initialize t, y, w, omega (angular frequencies), future.
fit = lombscargle_harmonic_fit(t, y, omega, 3, &
    weights=w, floating_mean=.true., estimate_covariance=.true.)

if (fit%valid(1)) then
    forecast = lombscargle_harmonic_predict(fit, future, 1, status)
    forecast_se = lombscargle_harmonic_predict_std(fit, future, 1, status)
    if (status == irregular_status_success) then
        ! forecast is the fitted mean; forecast_se is its standard error.
    end if
end if
```

The array-valued uncertainty fields are allocated **only** with
`estimate_covariance=.true.`. The scalar fields `degrees_of_freedom` and
`intercept_fitted` are always present:

| Result field | Shape | Meaning |
|---|---|---|
| `covariance` | `(2*H+1,2*H+1,nfreq)` | Covariance of fitted parameters |
| `standard_error` | `(2*H+1,nfreq)` | Square roots of marginal parameter variances |
| `residual_variance` | `(nfreq)` | Fitted weighted noise scale |
| `degrees_of_freedom` | scalar | `count(weights>0) - 2*H - (1 if floating_mean)` |
| `intercept_fitted` | logical | Whether the model estimated an intercept |

The first covariance row/column and standard-error entry refer to the
intercept; remaining entries are `cos(1), sin(1), cos(2), sin(2), ...`.
When `floating_mean=.false.`, the intercept is **not** a fitted parameter:
its covariance row/column and standard error are NaN. The harmonic block
remains available and `predict_std` omits the fixed intercept. `precenter`
adds a shift to the reported intercept but does not turn it into a fitted
parameter.

An invalid/rank-deficient frequency has `valid(k)=.false.` and NaN covariance
and uncertainties; other valid frequencies are retained. The new predictor
returns an all-NaN array and `irregular_status_invalid_input` for invalid
frequencies, missing covariance, malformed result arrays, nonfinite times,
or floating-point overflow. Empty time arrays are permitted for valid fits.

## Weighted least-squares model

At a fixed angular frequency, with `H` harmonics, the model is

\[
 y_i = c + \sum_{h=1}^H[a_h\cos(h\omega t_i)+b_h\sin(h\omega t_i)] + \epsilon_i.
\]

The covariance formulas assume **independent errors with variance inversely
proportional to the supplied positive observation weights**. Weights are
normalized to sum to one, just as for the original fit; rescaling all
weights leaves the estimates, standard errors, and covariances unchanged.
Observations with zero weight are excluded from the degrees of freedom.

Let \(A\) be the weighted-mean-centered Fourier design (or the uncentered
design when there is no floating intercept), \(W=\operatorname{diag}(w_i)\)
and \(A_w=W^{1/2}A=QR\). The estimated noise scale and harmonic
coefficient covariance are

\[
 s^2 = \frac{\sum_i w_i(y_i-\hat y_i)^2}{n_+-p},\qquad
 \operatorname{Cov}(\hat\beta)=s^2 R^{-1}R^{-T},
\]

where \(n_+\) is the number of positive weights and \(p=2H+1\) for a
floating intercept or \(p=2H\) otherwise. The code obtains \(R^{-1}\)
with triangular solves; it does **not** form/invert normal equations.
For a floating mean, with column weighted means \(m\), the intercept's
covariance blocks are

\[
 \operatorname{Var}(\hat c)=s^2+m^T\operatorname{Cov}(\hat\beta)m,
 \quad \operatorname{Cov}(\hat c,\hat\beta)=-m^T\operatorname{Cov}(\hat\beta).
\]

The standard error of the fitted **mean** at a new time is
\(\sqrt{v(t)^T\operatorname{Cov}(\hat\theta)v(t)}\), where \(v(t)\)
contains the intercept and Fourier basis values. It **does not** include
future observation noise; it is not a prediction interval.

## Statistical limitations

- This is conditional on the supplied **fixed trial frequency**, the chosen
  number of harmonics, and the inverse-variance-weight assumption. It does
  not account for choosing the frequency based on the same observations.
- When `precenter=.true.` and `floating_mean=.false.`, the precentering
  shift is treated as fixed. Its sampling variability is not estimated;
  choose `floating_mean=.true.` when inference on the intercept matters.
- The residual scale is estimated separately at each valid frequency; no
  heteroskedasticity-robust/sandwich covariance, correlated-noise model,
  multiple-testing correction, or confidence-interval coverage guarantee
  is provided. Inference assumes the model and noise structure are adequate.
- Severe harmonic collinearity is rejected by the existing QR rank test.
  `n_+ > p` is required; arbitrary high harmonic counts are not allocated
  when the original parameter validation rejects them.

## Independent validation

```bash
python tools/run_signal_harmonic_uncertainty.py
```

The runner performs checked and `-O3` GNU Fortran builds, native regression
checks, and 256 seeded weighted NumPy least-squares comparisons, validating
covariance matrices, coefficient standard errors, residual noise estimates,
and fitted-mean standard errors. It checks all combinations of floating and
fixed intercepts, precentering, omitted or present weights, zero-weight
observations, and weight rescaling. The original fit and multiharmonic
regressions are rerun. There are 59,796 finite reference values per build;
maximum observed absolute error was below `1.6e-15`. A separate public API
test checks exports via `scifort_signal`.
