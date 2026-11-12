# Multiharmonic Lomb–Scargle periodograms (SciFort extension)

SciFort provides `lombscargle_multiharmonic` and an optional trailing
`nterms=` keyword on the existing rank-one `lombscargle` function. They
fit a **joint Fourier-series model** to unevenly spaced observations,
with more than one harmonic at every trial angular frequency:

\[
y(t_i) \simeq c + \sum_{h=1}^{H}
  \left[a_h\cos(h\omega t_i)+b_h\sin(h\omega t_i)\right].
\]

`H` is the positive integer `nterms`. The constant `c` is fitted only when
`floating_mean=.true.`. This is a **SciFort extension**, not a keyword of
SciPy 1.17's public `scipy.signal.lombscargle`, which fits one harmonic.

## API

```fortran
use scifort_signal, only : dp, lombscargle, lombscargle_multiharmonic

real(dp) :: t(100), y(100), omega(256), w(100)
real(dp), allocatable :: p(:), pn(:)
integer :: status

! Initialize t, y, omega, w first.
p = lombscargle(t, y, omega, nterms=3, weights=w, &
                floating_mean=.true., status=status)

pn = lombscargle_multiharmonic(t, y, omega, 3, weights=w, &
                 floating_mean=.true., normalize='normalize', status=status)
```

The dedicated function takes `nterms` as its fourth required argument;
`lombscargle` accepts it as the **last optional keyword**. Omitting `nterms`
or setting it to one follows the original single-harmonic implementation,
so previous positional calls and their outputs remain unchanged.

Both functions accept `precenter`, `weights`, `floating_mean`, `normalize`
and `status`. The `freqs` argument contains **angular** frequencies in
radians per unit of `t`, not cycles per unit time. `precenter` subtracts
the unweighted sample mean *before* fitting; for a nonzero unknown offset,
`floating_mean` is generally more appropriate.

## Normalizations and numerical method

Let \(w_i\) be nonnegative weights normalized to sum to one, \(\hat y_i\)
be the fitted Fourier series, and \(y_i^*\) the observations after optional
precentering. For a floating mean, the reference model is the weighted mean
\(\bar y_w\); otherwise its reference value is zero. The model improvement is

\[
\Delta = \sum_i w_i(y_i^*-b)^2
       - \sum_i w_i(y_i^*-\hat y_i)^2,
\]

where \(b=\bar y_w\) with a floating mean and zero otherwise. The results are:

* `normalize='power'` (default): \(N\Delta/2\), consistent with SciFort's
  existing single-harmonic convention.
* `normalize='normalize'`: \(\Delta/\sum_i w_i(y_i^*-b)^2\), the fraction
  of weighted variance explained by all harmonics jointly, in [0, 1].

The implementation centers harmonic columns when fitting an intercept,
then uses twice-reorthogonalized modified Gram–Schmidt QR to evaluate the
least-squares projection. It does **not** form or invert normal equations.
Its per-frequency cost is \(O(NH^2)\), with \(O(NH)\) workspace.

## Error handling and limitations

* Nonfinite observations, times, or frequencies; negative/nonfinite weights;
  mismatched arrays; unsupported normalization; and `nterms < 1` return
  `irregular_status_invalid_input` and NaNs.
* There must be more positive-weight observations than fitted parameters
  (two per harmonic, plus the optional intercept).
* A frequency with an ill-conditioned or rank-deficient Fourier design
  returns NaN at that frequency and sets `irregular_status_singular`. Other
  valid frequencies are still computed.
* A constant signal with `floating_mean=.true.` has no reference variance;
  normalized power is undefined (NaN, singular status).
* This first implementation operates on **one real-valued time series**.
  The existing planned/rank-15 Lomb–Scargle APIs remain single-harmonic.
  The complex `lombscargle_amplitude` entry point also remains
  single-harmonic; it does not return an array of Fourier coefficients.
* This is a direct frequency-sweep algorithm rather than an FFT-based fast
  Lomb–Scargle implementation.

## Verification

Run from the repository root:

```bash
python tools/run_signal_multiharmonic.py
```

The runner builds checked and `-O3 -march=native` executables with
GNU Fortran, runs the new deterministic Fortran tests, and checks 384
randomized cases (2,493 frequencies) against **independent NumPy weighted
least squares**, including weights, DC offsets, and several harmonic counts.
Cases with `nterms=1` are additionally compared to SciPy 1.17. The separate
public-module test is `test/test_signal_multiharmonic_public.f90` and is
included in normal project testing. No NumPy or SciPy dependency is needed
to use the Fortran routine.
