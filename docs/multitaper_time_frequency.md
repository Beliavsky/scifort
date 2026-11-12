# Time-frequency multitaper estimation and harmonic-line detection

The procedures described here are SciFort **extensions** built on the existing
`scifort_signal_dpss` and `scifort_signal_spectral` kernels. They are not
one-for-one translations of public SciPy functions.

## Whole-record multitaper detrending

The following multitaper APIs accept an optional trailing
`detrend_mode='none'|'constant'|'linear'`, which overrides the older logical
`detrend` argument:

- `multitaper_psd`, including real/complex and multichannel overloads;
- `multitaper_spectral_matrix`, `multitaper_csd`, `multitaper_coherence`;
- `multitaper_psd_ci` (the jackknife intervals use the same detrended input);
- `multitaper_harmonic_ftest` and `multitaper_harmonic_scan`.

The default is `'constant'` for PSD/cross-spectrum/CI APIs and `'none'` for
the harmonic-line estimators, matching their previous logical defaults.
`'linear'` fits and subtracts one least-squares line per channel over the
**entire record** *before* tapering. For example:

```fortran
psd = multitaper_psd(x,3.5_dp,kmax=5,detrend_mode='linear')
mat = multitaper_spectral_matrix(channels,3.5_dp,kmax=5, &
    weighting='adaptive',axis=1,detrend_mode='linear')
```

The streaming `multitaper_stream%initialize` accepts the same keyword but
fits a separate line to each complete window, not the entire recording.
Both modes use the standard `scipy.signal.detrend(type='linear')` least-squares
convention. SciFort's multitaper estimators are extensions; SciPy DPSS,
linear detrending, and NumPy FFTs provide independent numerical checks.
Invalid `detrend_mode` values return the existing invalid-input status.

## `multitaper_spectrogram`

`multitaper_spectrogram(x, nw, nperseg, noverlap, kmax, fs, nfft, detrend, weighting, axis, detrend_mode)`
returns a `multitaper_spectrogram_result` for rank-1 real/complex signals,
or a `multitaper_spectrogram_channels_result` for rank-2 real/complex data.
The `axis` argument is accepted on rank-2 arrays only. Its default is 1
(samples are in the first dimension), and `axis=2` means that samples are
in the second dimension. Results contain `frequencies`, `times`, `power`,
and `status`.

- `nw`: positive time-half-bandwidth product; must be smaller than `nperseg/2`.
- `nperseg`: segment length in samples, at least 2 and at most `size(x,axis)`.
- `noverlap`: number of overlapping samples, from 0 to `nperseg-1`; default
  `nperseg/2` (integer division).
- `kmax`: DPSS count; default `max(1,int(2*nw)-1)`, capped at `nperseg`.
- `fs`: sampling frequency (default 1), with `power` in squared signal units/Hz.
- `nfft`: FFT length; must be at least `nperseg`, default `nperseg`.
- `detrend`: legacy logical per-segment arithmetic-mean removal, default true.
- `detrend_mode`: optional character selector `'none'`, `'constant'`, or
  `'linear'`. Each segment and channel is detrended independently **before**
  the DPSS taper is applied. `linear` removes its least-squares mean and
  slope, including independent real/imaginary slopes for complex data.
  When specified it overrides `detrend`; otherwise `detrend=.false.`
  corresponds to `none` and the default `detrend=.true.` to `constant`.
  Unknown values return invalid-input status and empty outputs.
  The selector is the last argument for both rank-1 and rank-2 overloads,
  after `axis` when the latter is applicable.
- `weighting`: `equal` (default), `eigen`, or frequency-dependent `adaptive`.

Time bins use the centres `(window_start_zero_based + nperseg/2)/fs`, matching
the legacy `scipy.signal.spectrogram` time convention. Only complete windows
are used; no end-padding or boundary extension is applied. Real inputs produce
one-sided PSDs (positive non-DC/non-Nyquist bins doubled), while complex
inputs produce unshifted signed two-sided FFT frequencies. Rank-1 results
have `power(frequency,time)` and rank-2 results have
`power(frequency,time,channel)`.

The DPSS eigensystem is calculated **once per call**. FFTs are then computed
per taper, segment, and channel. This removes redundant DPSS calculation
but is not a streaming implementation and does not cache transform plans.


## Time-frequency cross-spectra and coherence

Three more SciFort-specific APIs extend the cached-DPSS spectrogram engine:

* `multitaper_cross_spectrogram(x,y,nw,nperseg,noverlap,kmax,fs,nfft,detrend,weighting,detrend_mode)`
  takes two equally long rank-1 signals of real and/or complex type. It returns
  `multitaper_cross_spectrogram_result`: `frequencies`, `times`, the complex
  cross-density `power(frequency,time)`, the real diagonal `auto_x` and
  `auto_y` PSDs, and `status`.
* `multitaper_tf_coherence(x,y,nw,nperseg,noverlap,kmax,fs,nfft,detrend,weighting,detrend_mode)`
  takes the same inputs and returns `multitaper_tf_coherence_result` with
  the magnitude-squared time-frequency `coherence(frequency,time)` and
  corresponding frequency/time axes.
* `multitaper_spectral_matrix_spectrogram(x,nw,nperseg,noverlap,kmax,fs,nfft,detrend,weighting,axis,detrend_mode)`
  takes rank-2 real or complex channel data. Its result is
  `multitaper_spectral_matrix_spectrogram_result` with complex
  `power(frequency,time,x_channel,y_channel)` and `frequencies`, `times`,
  and `status`. Samples are in dimension 1 by default or dimension 2 with
  `axis=2` (one-based Fortran axis numbering).

These use the same settings, default values and input validation as
`multitaper_spectrogram`. A real-real pair or an all-real multichannel matrix
uses a one-sided spectrum and doubles positive non-DC bins except an even-
length Nyquist bin; any complex input uses unshifted two-sided FFT bins.

Cross-density follows SciPy's `csd` convention:

`P_xy = sum_k sqrt(q_x(k,f)*q_y(k,f))*conjg(X_k(f))*Y_k(f)`.

All four windowed time-frequency estimators now accept the trailing
`detrend_mode` string. The streaming and whole-record multitaper APIs also
accept that keyword, while retaining their previous logical `detrend`
argument for backward compatibility.
For example, a drifting pair can be analyzed using:

```fortran
xy = multitaper_cross_spectrogram(x,y,3.5_dp,64, &
    noverlap=32,kmax=5,detrend_mode='linear')
```

The `'linear'` option agrees with `scipy.signal.detrend(...,type='linear')`
on each complete window. The implementation uses centered sample indices
to remove a segment's least-squares slope and intercept. Tests include
`test/test_signal_multitaper_linear_detrend.f90` and the independent
`tools/signal_reference/compare_signal_multitaper_linear_detrend.py`
reference (64 seeded examples, 36,787 coefficients, maximum absolute
error 8.92e-15). The spectral-matrix API is a SciFort extension, so the
reference assembles the cross-periodograms using SciPy DPSS and NumPy FFTs.

The `q_x` and `q_y` are normalized nonnegative taper weights fitted **separately**
for each channel. This agrees with the existing `multitaper_spectrogram` PSDs
on the diagonal and enforces a Cauchy-Schwarz-consistent coherence estimate:

`C_xy = |P_xy|**2/(P_xx*P_yy)`.

Coherence is clipped to `[0,1]` to eliminate tiny roundoff excursions;
zero-power pairs return zero rather than a NaN. The full complex
cross-spectral matrix is Hermitian at each frequency and time bin.

**Important interpretation:** Coherence estimated from only one taper
is trivially one wherever both powers are nonzero. Use several tapers for
meaningful coherence estimates. Time-frequency coherence is a local,
smoothed estimate; it does not provide significance probabilities or a
confidence interval. The three APIs are **SciFort extensions**, not names
from the public SciPy API.

Example:

```fortran
use scifort_signal, only : dp, multitaper_cross_spectrogram, &
    multitaper_tf_coherence, multitaper_spectral_matrix_spectrogram, &
    multitaper_cross_spectrogram_result, multitaper_tf_coherence_result, &
    multitaper_spectral_matrix_spectrogram_result
real(dp) :: x(256), y(256), multichannel(256,2)
type(multitaper_cross_spectrogram_result) :: xy
type(multitaper_tf_coherence_result) :: coh
type(multitaper_spectral_matrix_spectrogram_result) :: matrix
! ... initialize x and y with data ...
multichannel(:,1) = x
multichannel(:,2) = y
xy = multitaper_cross_spectrogram(x,y,3.5_dp,64,noverlap=32, &
    kmax=5,fs=100.0_dp,weighting='adaptive')
coh = multitaper_tf_coherence(x,y,3.5_dp,64,noverlap=32,kmax=5, &
    fs=100.0_dp,weighting='adaptive')
matrix = multitaper_spectral_matrix_spectrogram(multichannel,3.5_dp,64, &
    noverlap=32,kmax=5,fs=100.0_dp,weighting='adaptive')
```

Verification: `test/test_signal_multitaper_tf_cross*.f90` tests mixed types,
real and complex modes, auto-spectrogram agreement, Hermitian symmetry,
axis selection, perfect/zero-power coherence and invalid inputs. The new
independent test `python test/compare_signal_multitaper_tf_cross.py` uses
SciPy 1.17 DPSS windows and NumPy FFTs in 120 randomized cases, with
maximum absolute spectral and coherence differences of 1.30e-11 and
2.90e-11, respectively. The complete standalone signal suite passes **89/89**
at both runtime-checked `-O0` and optimized `-O3` using GNU Fortran 14.2.0;
the targeted new regression also passes `-O3` with bounds checks and
floating-point traps. Project-wide `fpm test` was not run (FPM unavailable).
No additional Fortran dependencies or upstream numerical code were imported.

## `multitaper_noise_floor`: local colored-noise background

`multitaper_noise_floor(x,nw,median_bins,kmax,fs,nfft,detrend,weighting,detrend_mode)`
is a **SciFort extension** accepting rank-one real or complex signals. It calls
`multitaper_psd` and estimates a local spectral noise background with an
odd-width running median:

`background(f) = median(power(f-h), ..., power(f+h))`, with
`h=(median_bins-1)/2`. Real-input, one-sided spectra use half-sample
**reflect** boundaries (the `scipy.ndimage.median_filter(mode='reflect')`
convention). Complex-input, unshifted two-sided spectra use cyclic **wrap**
boundaries, consistent with their periodic Fourier-bin indexing.

The result type `multitaper_noise_floor_result` contains `frequencies`,
`power`, `background`, `whitened=power/background`, and `status`. All array
outputs are empty on failure. `median_bins` must be odd, at least 3, and at
most the number of returned spectral bins. Invalid PSD arguments are reported
through the same `spectral_status_invalid_input` convention. A zero power
sample gives `whitened=0` even if its background vanishes; a positive power
with zero background saturates at `huge(1.0_dp)` instead of dividing by zero.

```fortran
use scifort_signal, only : dp, multitaper_noise_floor, &
    multitaper_noise_floor_result, spectral_status_success
type(multitaper_noise_floor_result) :: local
real(dp) :: x(256)
! Initialize x before calling the estimator.
local = multitaper_noise_floor(x,3.5_dp,9,kmax=6, &
    weighting='adaptive',detrend_mode='linear')
if (local%status == spectral_status_success) then
    print *, maxval(local%whitened), 'largest local power ratio'
end if
```

This is a **descriptive** robust background estimator, not a statistical
hypothesis test or spectral-whitening filter applied to the time signal.
The window should span more bins than a narrow spectral line while remaining
short enough to track broadband spectral curvature. Strong broad peaks or
multiple nearby lines can contaminate the median. Zero-padding changes the
number of grid bins per physical bandwidth, so `median_bins` must be scaled
accordingly. No uniform false-alarm probability is assigned to the ratios.
Use `multitaper_harmonic_scan` for approximate line p-values and multiplicity
correction, with its stated null-model limitations.

`test/compare_signal_multitaper_noise_floor.py` validates 144 seeded real and
complex colored-AR(1) / sinusoid cases against SciPy 1.17 DPSS tapers,
NumPy FFT, SciPy linear detrending, and SciPy's independent median filter.
The 10,889 compared bins agree in checked and optimized GNU Fortran builds;
maximum absolute PSD, background, and normalized errors are approximately
7.85e-13, 4.94e-13, and 1.90e-12, respectively. No new Fortran dependency
is introduced.

## Windowed and streaming colored-noise backgrounds (SciFort extension)

`multitaper_noise_floor_spectrogram(x,nw,median_bins,nperseg, &
noverlap,kmax,fs,nfft,detrend,weighting,axis,detrend_mode)` adds local-median
noise backgrounds to overlapping multitaper PSD windows. It follows the
parameters and **complete-window only** convention of
`multitaper_spectrogram`. `axis` is available only for rank-2 inputs and is
one-based (`axis=1` by default, or `axis=2` for channels-by-samples).
`median_bins` is required and must be odd, from 3 through the number of
frequency bins. A zero-width time-frequency output is possible for a valid
streaming chunk that completes no new windows.

For rank-1 real or complex input, the
`multitaper_noise_floor_spectrogram_result` has `power`, `background`, and
`whitened` arrays of shape `(frequency,time)`. Rank-2 inputs return the
`multitaper_noise_floor_spectrogram_channels_result` with those three arrays
of shape `(frequency,time,channel)`. Both also contain `frequencies`, `times`,
and `status`. For every frame and channel the local median runs **along the
frequency axis only**, not across time or other channels. As in the
whole-record estimator, one-sided (real) spectra reflect their frequency
boundaries, while unshifted two-sided (complex) spectra wrap periodically.
The zero/zero ratio is zero and positive power with zero background is
saturated to `huge(1.0_dp)`.

```fortran
use scifort_signal, only : dp, multitaper_noise_floor_spectrogram, &
    multitaper_noise_floor_spectrogram_channels_result, spectral_status_success
real(dp) :: channels(300,2)
type(multitaper_noise_floor_spectrogram_channels_result) :: floor
! Populate channels with samples for two synchronized channels.
floor = multitaper_noise_floor_spectrogram(channels,3.5_dp,7,64, &
    noverlap=32,kmax=5,nfft=128,weighting='adaptive',detrend_mode='linear')
if (floor%status == spectral_status_success) then
    print *, shape(floor%whitened) ! frequency, complete windows, channels
end if
```

For ongoing signals, configure a `multitaper_stream` in the usual way, then
call `stream%push_noise_floor(chunk,median_bins,axis)` with real/complex rank-1
or rank-2 chunks. Each call returns a
`multitaper_noise_floor_spectrogram_channels_result` **even for one
channel**. This includes only the *newly completed* time windows; no past
noise-floor frames are stored in the stream. The `axis` argument applies to
rank-2 chunks only. An invalid width returns invalid-input status without
consuming the chunk. Other invalid inputs retain the existing stream's
transactional behavior. Calls to `push()` and `push_noise_floor()` consume
samples from the same underlying stream state; do not feed the same chunk to
both.

`stream%mean_noise_floor(median_bins)` returns the separate
`multitaper_noise_floor_channels_result`: `power(frequency,channel)`,
`background(frequency,channel)`, `whitened(frequency,channel)`, `frequencies`
and `status`. It applies the median **after averaging complete-window power
spectra**, not to the average of framewise normalized ratios. Before the
first complete window, it reports invalid-input status and empty arrays.
`reset()` clears the running mean, preserving tapers and detrend settings.

These ratios are descriptive, not formal spectral-line p-values. Strong
wide-band peaks, multiple nearby lines, and zero-padding may affect a
frequency-domain local median. The same numerical caveats as for
`multitaper_noise_floor` apply. Streaming state storage remains bounded by
the FFT/segment/channel dimensions, although each returned result grows
with the number of new complete windows in its input chunk.

`test/test_signal_multitaper_noise_floor_window_stream.f90` checks batch vs.
chunked output, mean-spectrum semantics, zero spectra, invalid arguments,
nonfinite chunks, transposed channels, all three weighting methods, and
real/complex input. `test/compare_signal_multitaper_noise_floor_window.py`
compares 120 seeded cases (47,098 time-frequency PSD bins) against independent
SciPy 1.17 periodic DPSS, `signal.detrend`, NumPy FFT, and
`ndimage.median_filter`. Checked and optimized builds passed, with maximum
absolute power, background, and normalized differences of about 2.14e-11,
2.14e-11, and 1.11e-10, respectively.

## `multitaper_harmonic_ftest`

`multitaper_harmonic_ftest(x, nw, kmax, fs, nfft, detrend, detrend_mode)` accepts rank-1
real or complex inputs and returns `multitaper_harmonic_result` with
`frequencies`, `f_statistic`, `pvalue`, `amplitude`, and `status`.
`kmax` must be at least 2. The default is
`min(size(x),max(2,int(2*nw)-1))`; `detrend` defaults to false.

Each DPSS taper produces an eigencoefficient `Y_k(f)`. For its zero-frequency
taper response `h_k = sum_n v_k(n)`, fit the line coefficient

`a(f) = sum_k h_k * Y_k(f) / sum_k h_k**2`.

At each FFT frequency, compute

`F(f) = (K-1) * |a(f)|**2 * sum_k h_k**2 / sum_k |Y_k(f)-a(f)*h_k|**2`.

Under Gaussian, locally smooth-spectrum noise, the test is approximately
`F(2, 2K-2)` distributed at **interior** Fourier frequencies. The returned
`pvalue` is its survival probability. The test assumes individual frequency
hypotheses; it **does not** correct for inspecting many bins. Real-signal
DC and Nyquist endpoints do not generally have this same null distribution.
For an isolated sinusoid at a positive interior frequency, the complex
coefficient `a` is about half its physical real-sinusoid amplitude; for a
complex exponential, `a` is the exponential coefficient.

If no line is present in a zero signal, the implementation returns an
F-statistic of 0 and p-value of 1. An exact line with zero residual produces
a saturated, very large F statistic and p-value zero.

## `multitaper_harmonic_scan`: multiplicity-corrected line detection

`multitaper_harmonic_scan(x,nw,kmax,fs,nfft,detrend,method,alpha,include_edges,test_bins,detrend_mode)`
is a **SciFort extension** wrapping the existing real/complex Thomson
`multitaper_harmonic_ftest`. Its result, `multitaper_harmonic_scan_result`,
contains the original `frequencies`, `f_statistic`, `pvalue`, `amplitude`,
plus `adjusted_pvalue`, `significant`, `num_tests`, and `status`.

- `method='bh'` (default) applies Benjamini–Hochberg false-discovery-rate
  correction to the selected bins. `method='by'` applies the more conservative
  Benjamini–Yekutieli correction, appropriate under arbitrary dependence;
  `method='holm'` applies Holm's step-down family-wise error correction under
  arbitrary dependence. `method='bonferroni'` provides another family-wise
  error correction; `method='none'` leaves probabilities unchanged.
- `alpha=0.05_dp` by default is the threshold for the adjusted probabilities.
  Only values strictly between zero and one are accepted.
- By default, DC is excluded, as is Nyquist for a real input with even FFT
  length; these endpoints do not obey the usual interior F null model.
  `include_edges=.true.` includes those bins at the caller's discretion.
- `test_bins` overrides the default selection and accepts unique **one-based**
  frequency indices, for example `[4,7,9]`, for a prespecified family of
  comparisons. Unselected entries have adjusted p-value 1 and cannot be
  marked significant. Empty selections are valid.

```fortran
type(multitaper_harmonic_scan_result) :: lines
real(dp) :: signal(256)
! Signal must be initialized before use.
lines = multitaper_harmonic_scan(signal,3.5_dp,kmax=6,method='by',alpha=0.01_dp)
if (lines%status == spectral_status_success) then
    print *, count(lines%significant), 'significant Fourier bins'
end if
```

The multiple-testing algorithms are independently implemented and use
`O(n log n)` stable sorting. Correction does not fix inaccuracies in the
underlying approximate Thomson F null model. In particular, BH controls FDR
under independence or suitable positive dependence, which must **not** be
assumed automatically for nearby/overlapping spectral bins; BY is more
conservative under dependent valid p-values. Searching adaptively for bins
and then correcting *only* those bins is not a valid global scan without
accounting for the selection process.

The independent `test/compare_signal_multitaper_harmonic_scan.py` checks 225
seeded cases against SciPy DPSS/NumPy FFT and SciPy's F distribution,
SciPy's `false_discovery_control` for BH/BY, and Statsmodels `multipletests`
for Holm. All passed; maximum corrected-p-value error was about 5.03e-13.
Statsmodels is a **test-only** dependency. No new Fortran dependencies were
introduced.

## Verification

The complete standalone signal suite passed **87/87** tests with `-O0` and
runtime checks, and **87/87** at `-O3`, using GNU Fortran 14.2.0. New targeted
optimized tests passed with bounds checks and floating-point traps.

`python test/compare_signal_multitaper_time_frequency.py` independently
checks 90 real/complex, single/multichannel spectrogram cases and 100
harmonic-test cases using SciPy 1.17 DPSS tapers, NumPy FFT and SciPy F
survival probabilities. All 190 cases passed. Maximum observed absolute
PSD discrepancy: 8.83e-12; maximum p-value discrepancy: 2.48e-13.
Project-wide `fpm test` was not run (FPM not installed).

The algorithms are independently written under SciFort's MIT license;
no third-party Fortran routines or new runtime dependencies were added.
