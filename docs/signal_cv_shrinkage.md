# Leave-one-taper-out cross-validated spectral shrinkage (SciFort extension)

The static and moving-window multichannel multitaper estimators now provide
**opt-in, frequency-local cross-validation** of diagonal-target shrinkage.
Unlike SciFort's existing Gaussian-Wishart `auto_shrinkage=.true.` mode,
`cv_shrinkage=.true.` works with **Thomson adaptive taper weighting**, as well
as equal and DPSS eigenvalue weighting. It fits a *predictive* coefficient using
individual Fourier spectra, **not** a plug-in effective taper count. This is
an original SciFort extension, not an API from `scipy.signal`.

## Signal-level use

```fortran
use scifort_signal, only : dp, partial_coherence_result, &
    partial_coherence_tf_result, multitaper_partial_coherence, &
    multitaper_tf_partial_coherence

real(dp) :: x(512,4)  ! populate: (sample, channel)
type(partial_coherence_result) :: static
type(partial_coherence_tf_result) :: moving

static = multitaper_partial_coherence(x,3.5_dp,kmax=6, &
    weighting='adaptive',cv_shrinkage=.true.,ridge=0.02_dp, &
    return_coherency=.true.)
print *, static%shrinkage_used  ! length: nfrequency

moving = multitaper_tf_partial_coherence(x,3.5_dp,nperseg=128, &
    noverlap=64,kmax=6,weighting='adaptive', &
    cv_shrinkage=.true.,ridge=0.02_dp)
print *, shape(moving%shrinkage_used)  ! (nfrequency,nwindow)
```

Both calls work for real and complex channel arrays; `axis=2` selects
channels on rows instead of columns. Each frequency/time cell gets a separate
coefficient between 0 (no shrinkage) and 1 (diagonal target). Shrinkage
preserves each channel's power spectrum and is applied **before** separately
requested trace-relative ridge loading. Result precision, partial coherence,
complex partial coherency, `valid` masks, and numerical status use the
existing Hermitian-Cholesky implementation.

To reuse a precomputed spectrum *with its cross-validation metadata*:

```fortran
use scifort_signal, only : multitaper_spectral_matrix, &
    partial_coherence_matrix, multitaper_matrix_result, partial_coherence_result

type(multitaper_matrix_result) :: spectral
type(partial_coherence_result) :: conditioned

spectral = multitaper_spectral_matrix(x,3.5_dp,kmax=6, &
    weighting='adaptive',cv_shrinkage=.true.)
conditioned = partial_coherence_matrix(spectral,cv_shrinkage=.true.)
! spectral%cv_shrinkage(:) matches conditioned%shrinkage_used(:)
```

The analogous `multitaper_spectral_matrix_spectrogram` result has
`%cv_shrinkage(f,t)`; `partial_coherence_tf_matrix` accepts
`cv_shrinkage=.true.` and uses these values. Without explicit
`cv_shrinkage=.true.` at spectral estimation time, the metadata array is
**unallocated**, and requesting CV at the inversion stage returns
`partial_coherence_invalid_input` rather than silently using the Wishart
estimator. Spectral power is unchanged (up to roundoff) whether CV is
requested or not.

## Objective and computation

Let `J_k(f)` be the vector of tapered, normalized complex Fourier coefficients
for taper `k` and frequency `f`, with the same conjugation convention as
SciPy CSD: `P_ij = conjg(J_i)*J_j`. For every held-out taper `k`, recompute
the spectral weights using **only the other tapers' periodograms** (for
adaptive weighting, re-run the Thomson frequency-local iterations). Denote
the resulting multichannel training cross-spectrum by `S^(-k)(f)` and the
held-out rank-one spectrum by `H^(k)(f)`. For channel-specific weights
`q_{l,i}`, the training cross-spectrum is

\[
S_{ij}^{(-k)}(f) = \sum_{\ell \ne k} \sqrt{q_{\ell i}^{(-k)}(f)
q_{\ell j}^{(-k)}(f)}\,\overline{J_{\ell i}(f)}J_{\ell j}(f),
\qquad H_{ij}^{(k)}(f)=\overline{J_{ki}(f)}J_{kj}(f).
\]

Only the off-diagonals contribute to the held-out predictive loss:

\[
L_f(\beta)=\sum_k\sum_{i<j}
\left|\beta S_{ij}^{(-k)}(f)-H_{ij}^{(k)}(f)\right|^2,
\qquad \beta=1-\alpha.
\]

The unconstrained real slope is obtained by scalar least squares and clipped:

\[
\hat\beta_f=\operatorname{clip}_{[0,1]}\frac{
\sum_{k,i<j}\operatorname{Re}\{\overline{S_{ij}^{(-k)}}H_{ij}^{(k)}\}}{
\sum_{k,i<j}|S_{ij}^{(-k)}|^2},
\qquad \hat\alpha_f=1-\hat\beta_f.
\]

If the denominator vanishes, choose full diagonal shrinkage. For an entirely
zero spectral matrix the normal matrix inversion still marks that cell
singular, with a NaN `shrinkage_used` value. The risk calculation normalizes
spectral vectors frequency-by-frequency to prevent the physical power units
from overflowing the numerator/denominator. One-sided real FFT doubling
cancels from the cross-validated slope and is applied to the final CSD as
usual.

## Important statistical and API limits

- This criterion selects shrinkage to **predict a held-out taper**, not to
  minimize an exactly unbiased estimate of full-sample covariance error.
  It is not an oracle, a Ledoit-Wolf estimator, or a calibrated confidence
  interval/p-value. When only a few tapers are available, the coefficient can
  have considerable variability.
- Thomson's segment noise/variance estimate still uses the full signal
  segment, even though frequency-local adaptive weights are refitted after
  each taper's removal. Thus, adaptive CV is not perfectly independent
  data splitting, and its reported coefficient has no exact unbiased-risk
  guarantee.
- DPSS tapers are mathematically orthogonal but **not always statistically
  independent**, particularly for colored/nonstationary signals, line
  components or non-Gaussian observations. Hold-out spectra with low DPSS
  concentration can have different bias/variance than training spectra.
  DC/Nyquist bins for real data need additional caution.
- CV requires at least **three tapers** and **two channels**. Fewer tapers or
  channels produce `spectral_status_invalid_input` at spectral estimation.
  `cv_shrinkage=.true.` cannot be combined with `shrinkage=` or with
  `auto_shrinkage=.true.`; matrix-level `effective_tapers=` also conflicts.
  Invalid coefficient metadata (wrong shape, NaN, or values outside [0,1])
  is rejected.
- CV is more expensive than the Gaussian-Wishart plug-in: each omitted taper
  requires re-estimation of weights and accumulation of cross-spectra. Its
  additional workspace grows as `O(frequency * taper * channels)` per window;
  no full recording of FFTs across time windows is retained.
- The original `auto_shrinkage=.true.` implementation is unchanged and still
  **rejects adaptive taper weighting**, since the fixed-weight Wishart moment
  identity does not hold for adaptive weights.

## Regressions

```bash
python tools/run_signal_partial_shrinkage_cv.py
python tools/run_signal_partial_shrinkage_auto.py
python tools/check_signal_partial_shrinkage_auto_public.py
```

The CV suite builds native tests with checked and optimized GNU Fortran,
including the previous shrinkage, coherency and multitaper tests. An
independent Python/NumPy/SciPy reference separately computes DPSS windows,
FFT coefficients, equal/eigen/adaptive weights, leave-one-out refits,
cross-spectra and shrinkage for **12 static or time-frequency real/complex
configurations (1,677 frequency/time cells)**. This comparison is strictly
numerical validation of the implementation, **not** evidence of statistical
calibration or guaranteed prediction improvement.
