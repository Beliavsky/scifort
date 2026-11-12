# Weighted harmonic coefficient recovery and prediction (SciFort extension)

`lombscargle_harmonic_fit` complements the multiharmonic Lomb–Scargle
periodogram by returning the **fitted Fourier-series coefficients** at each
trial angular frequency. It also reports intercepts, residual RMS errors,
unnormalized and normalized periodogram powers, and a validity flag for every
frequency. `lombscargle_harmonic_predict` evaluates a fitted series at new
times without refitting.

These interfaces are original **SciFort extensions** to the `scipy.signal`
API, not functions in SciPy 1.17. The underlying model matches the weighted
multiharmonic least-squares formulation described in
[Multiharmonic Lomb–Scargle](signal_multiharmonic.md).

## Example

```fortran
use scifort_signal, only : dp, harmonic_fit_result, &
    lombscargle_harmonic_fit, lombscargle_harmonic_predict, &
    irregular_status_success

real(dp) :: time(100), signal(100), weights(100), freqs(2), new_times(20)
real(dp), allocatable :: forecast(:)
type(harmonic_fit_result) :: fit
integer :: status

! Initialize time, signal, weights, freqs (angular frequency), new_times.
fit = lombscargle_harmonic_fit(time, signal, freqs, 3, &
    weights=weights, floating_mean=.true.)

if (fit%valid(2)) then
    forecast = lombscargle_harmonic_predict(fit, new_times, 2, status)
    if (status == irregular_status_success) then
        ! forecast is the 3-harmonic reconstructed signal at new_times.
    end if
end if
```

The fitted model for frequency index `k` is

\[
\hat y_k(t)=c_k+\sum_{h=1}^{H} a_{hk}\cos(h\omega_k t)
                         + b_{hk}\sin(h\omega_k t).
\]

The result contains:

| Field | Shape | Meaning |
|---|---|---|
| `frequencies` | `(nfreq)` | Input angular frequencies (radians/time) |
| `nterms` | scalar | Number of fitted harmonics |
| `offset` | `(nfreq)` | `c_k` on the **original** data scale |
| `cosine`, `sine` | `(nterms,nfreq)` | `a_hk` and `b_hk`, indexed by harmonic and frequency |
| `power` | `(nfreq)` | `N/2` times weighted reduction in squared residuals |
| `normalized_power` | `(nfreq)` | Fractional explained weighted variance |
| `rms` | `(nfreq)` | Square root of weighted mean squared reconstruction error |
| `valid` | `(nfreq)` | Whether the harmonic design is full rank at that frequency |
| `covariance`, `standard_error`, `residual_variance` | optional | Available with `estimate_covariance=.true.`; see [uncertainty](signal_harmonic_uncertainty.md) |
| `status` | scalar | Zero for full success, or invalid/singular status |

All returned indices are Fortran **one-based**. `predict` accepts any finite
rank-one time vector, including nonmonotonic times and empty arrays. Output
NaNs and `irregular_status_invalid_input` indicate an invalid frequency index,
invalid fit, or invalid times.

## Weighted least squares and intercepts

Weights must be finite and nonnegative, with at least one positive weight.
The algorithm normalizes them to sum to one; multiplying every weight by a
positive constant does not affect the fit. A zero weight excludes the
observation without changing array shapes.

`precenter=.true.` subtracts the **unweighted** sample mean before fitting.
The reported intercept adds this shift back, so `predict` always reproduces
the **original** observation scale. `floating_mean=.true.` adds an independently
estimated intercept at every frequency; the harmonic basis columns are
weighted-mean-centered before solution. When no floating mean is used,
`offset` is zero (or the precentering shift, if requested).

For each frequency the implementation builds a weighted Fourier design,
performs **twice-reorthogonalized modified Gram–Schmidt QR**, collects the
upper-triangular `R` coefficients, and solves `R*c=Q^T*y` by back
substitution. This avoids directly forming or inverting normal equations.
Coefficients and predictions are computed only for full-rank frequencies.
There is no high-order polynomial expansion or other external dependency.

## Error behavior and limitations

The input must contain more **positive-weight samples** than fitted harmonic
parameters (`2*nterms`, plus one when `floating_mean=.true.`). Nonfinite
inputs, mismatched dimensions, negative weights, or impossible harmonic
counts return `irregular_status_invalid_input`. No huge coefficient array is
allocated when the requested count exceeds the number of observations.

At a frequency with collinear harmonics (for example zero frequency), the
coefficients, intercept, diagnostics, and prediction are NaN and `valid(k)`
is false. Other frequency bins remain usable, while the aggregate status is
`irregular_status_singular`. Constant signals with a floating intercept have
zero reference variance and return the singular status. The normalized power
is defined relative to the weighted floating mean when that option is used,
or zero after optional precentering otherwise.

This routine fits **one real-valued time series** across an arbitrary
rank-one frequency grid. Independent channels sharing one irregular sample
and weight grid can be fitted efficiently with cached factors; see
[Prepared multiharmonic fitting](signal_harmonic_plan.md). Joint multichannel
parameter fitting, regularization, and nonlinear frequency optimization are not
provided. Optional covariance and fitted-mean standard errors are described in
[Harmonic uncertainty](signal_harmonic_uncertainty.md). Numerically close trial frequencies can give ill-conditioned
harmonic designs; the QR rank check reports invalid bins rather than
silently returning unstable coefficients.

## Verification

```bash
python tools/run_signal_harmonic_fit.py
```

The test runner builds checked and optimized GNU Fortran binaries, runs native
coefficient and forecasting regressions, reruns the earlier multiharmonic
regression, and compares **256 seeded fits** with independent NumPy weighted
least squares. It checks the returned intercept, every sine/cosine
coefficient, both power normalizations, residual RMS, and predictions at
previously unseen times. The one-harmonic cases are also checked against
SciPy 1.17 `signal.lombscargle`. There are 20,061 numeric comparisons per
build. A separate public-module test exercises the exported result type and
prediction routine through `scifort_signal`.
