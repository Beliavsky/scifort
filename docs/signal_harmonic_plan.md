# Prepared multiharmonic fitting for many irregularly sampled signals

**SciFort extension.** `harmonic_fit_plan` caches the weighted, twice-reorthogonalized
QR factorization of a Fourier-series design for a fixed irregular time grid,
frequency grid, harmonic count, and optional observation weights. Many signal
channels can then be fitted without rebuilding the time-dependent harmonic
basis or refactorizing it. There is no corresponding public SciPy 1.17 API.

The existing unprepared `lombscargle_harmonic_fit` and periodogram functions
are unchanged. Prepared fits return the **same** `harmonic_fit_result` type,
so the existing `lombscargle_harmonic_predict` and
`lombscargle_harmonic_predict_std` functions apply unchanged.

## Example

```fortran
use scifort_signal, only : dp, harmonic_fit_plan, harmonic_fit_result, &
    lombscargle_harmonic_prepare, lombscargle_harmonic_fit_batch, &
    lombscargle_harmonic_predict, lombscargle_harmonic_predict_std

type(harmonic_fit_plan) :: plan
type(harmonic_fit_result), allocatable :: fits(:)
real(dp) :: t(100), y(100,5), weights(100), freqs(3), future(10)
real(dp), allocatable :: forecast(:), mean_se(:)

! Supply actual times, signals, weights, frequencies, and future times.
plan = lombscargle_harmonic_prepare(t, freqs, 3, &
    weights=weights, floating_mean=.true.)
fits = lombscargle_harmonic_fit_batch(plan, y, &
    precenter=.true., estimate_covariance=.true.)
if (fits(2)%valid(1)) then
    forecast = lombscargle_harmonic_predict(fits(2), future, 1)
    mean_se = lombscargle_harmonic_predict_std(fits(2), future, 1)
end if
```

`y` has shape `(nsamples,nchannels)` (time/sample index first).
`fits(j)` is a complete independent result for channel `j`. The same plan
may be reused in subsequent calls using new signals sampled at exactly the
same times. To fit one channel with the same plan:

```fortran
single_fit = lombscargle_harmonic_fit_prepared(plan, y(:,2), &
    precenter=.true., estimate_covariance=.true.)
```

The plan construction and fitting API are:

- `lombscargle_harmonic_prepare(x, freqs, nterms, weights, floating_mean)`
  returns `type(harmonic_fit_plan)`.
- `lombscargle_harmonic_fit_prepared(plan, y, precenter, estimate_covariance)`
  returns `type(harmonic_fit_result)` for one channel.
- `lombscargle_harmonic_fit_batch(plan, y, precenter, estimate_covariance)`
  returns an allocatable array of `harmonic_fit_result`, one per channel.

All `x`, `y`, `freqs`, and `weights` values are double-precision real; frequencies
are angular frequencies in radians per unit time. The plan accepts optional
nonnegative weights shared by **all** channels. A different weight vector
requires building a different plan. Multiplying the shared weights by a
positive constant has no effect because weights are normalized to sum one.
The `floating_mean` choice is fixed at plan construction; `precenter` is
independent for each channel and can be changed on every fit call. Computing
coefficient covariance and fitted-mean standard errors is opt-in. Covariance
assumptions and formulas are described in
[Harmonic uncertainty](signal_harmonic_uncertainty.md).

## Results and failure behavior

`plan%status`, `plan%valid(:)`, `plan%frequencies(:)`,
`plan%nterms`, and `plan%floating_mean` are public diagnostic fields.
The underlying QR factors and normalized weights cannot be modified through
the public plan API. A zero-frequency bin or another rank-deficient design
marks `plan%valid(k)=.false.` and produces the aggregate
`irregular_status_singular`, but good bins remain available. An invalid plan
has `irregular_status_invalid_input`. Very large or impossible harmonic counts
are rejected before allocating high-dimensional QR matrices.

Each channel has its own `valid(:)` and `status`. An invalid channel (for
example, with NaN observations) returns NaN-valued outputs and invalid status
without corrupting other channels. A constant channel with a floating
intercept has a singular reference variance. Bad sample dimensions return
invalid results, not truncated fits. A zero-channel matrix produces an empty
array of fits. Strided time inputs and noncontiguous channel arrays work.

## Complexity and memory

Write `N` for observations, `M` for frequency bins, `H` for harmonics,
`P=2H` for the number of sine/cosine coefficients, and `C` for channels.
Preparing a plan performs approximately `O(M*N*P**2)` QR work once; each
channel then needs approximately `O(M*N*P + M*P**2)` operations for weighted
projections and residuals plus triangular coefficient solves. Optional
covariance adds a small factorization of the cached upper-triangular matrices
per frequency/channel. The principal cached allocation is the weighted `Q`
array with shape `(N,P,M)`, plus `R(P,P,M)` and weighted harmonic means.
Consequently the plan trades extra memory for lower costs when many series
share the same sample geometry. For just one signal or a very large frequency
grid, the unprepared fitting routine may be preferable.

The implementation uses real weighted QR, not normal-equations inversion.
The root/sine/cosine factors follow the original `lombscargle_harmonic_fit`
ordering and floating-intercept conventions. Individual channel offsets,
coefficients, normalized/unnormalized powers, RMS, and uncertainty estimates
are computed separately.

## Verification

```bash
python tools/run_signal_harmonic_plan.py
```

The standalone runner performs checked and optimized GNU Fortran builds,
executes new and previous native regressions, and checks randomized
multichannel results against independent NumPy weighted least squares and
SciPy's one-harmonic Lomb–Scargle output. It verifies zero-weight observations,
precentering, both intercept choices, optional covariance, invalid frequency
bins, and real-valued coefficients. A public-module regression is
`test/test_signal_harmonic_plan_public.f90`.

An optional CPU-time benchmark is supplied in
`tools/benchmark_signal_harmonic_plan.f90`. It compares 20 channels on a
160-sample/64-frequency grid with three harmonics and ten repetitions,
including the one-time preparation cost each repetition. On the GNU Fortran
14.2 `-O3` test environment, an example run took about 0.93 s for repeated
ordinary fits and 0.049 s for prepared batches (roughly 19 times faster).
This is an illustrative, hardware- and problem-dependent result, **not** a
promised speedup. The extra `Q(N,2H,M)` allocation is part of that tradeoff.
