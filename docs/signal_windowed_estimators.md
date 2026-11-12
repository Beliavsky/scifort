# Parameterized and caller-supplied windows in legacy signal estimators

The one-dimensional `scifort_signal` methods `periodogram`, `welch`, `csd`,
`coherence`, `stft`, and `spectrogram` and their rank-independent caller-output
counterparts support SciPy-compatible window families that require numerical
parameters. Both real and complex input overloads are supported. The new
optional keywords are **appended to all existing argument lists** so earlier
positional and keyword calls are source-compatible.

| Keyword | Meaning | Examples |
| --- | --- | --- |
| `window` | Window name, as before | `'hann'`, `'kaiser'`, `'tukey'` |
| `window_parameter` | First named-window shape parameter | Kaiser beta; Gaussian standard deviation; Tukey alpha; general Gaussian power |
| `window_parameter2` | Second parameter where supported | General Gaussian width; Taylor sidelobe attenuation |
| `window_coefficients` | Real one-dimensional cosine weights | `window='general_cosine'`, `window_coefficients=[0.5_dp,0.4_dp,0.1_dp]` |
| `window_values` | Real one-dimensional array of explicit window samples | `window_values=weights` |

These correspond to SciPy's `window=('kaiser', beta)`,
`window=('general_gaussian', p, sig)`, `window=('general_cosine', a)`, and
`window=np.asarray(weights)` conventions respectively. Direct Fortran window
arrays use the new `window_values=` keyword; the existing `window=` argument
remains a character string. `window_parameter2` is supported for
`general_gaussian` and `taylor`. Other families accepting `window_parameter`
include Tukey, Gaussian, exponential, Kaiser, Chebyshev, Taylor,
Kaiser-Bessel-derived, general Hamming, and DPSS.

## Examples

```fortran
use scifort_signal, only : dp, welch, stft_into, spectral_result, spectral_status_success
real(dp) :: x(128), channels(128,3)
real(dp) :: taper(32)
complex(dp) :: z(17,3,7)
real(dp), allocatable :: freq(:), time(:)
type(spectral_result) :: estimate
integer :: st

! A parameterized Gaussian window, analogous to SciPy ('gaussian', 5.0).
estimate = welch(x, fs=100.0_dp, window='gaussian', &
    window_parameter=5.0_dp, nperseg=32, noverlap=16)
if (estimate%status /= spectral_status_success) error stop 'Welch failed'

! An explicit 32-point array determines nperseg when omitted.
taper = [(0.5_dp-0.5_dp*cos(2.0_dp*acos(-1.0_dp)*real(st-1,dp)/32.0_dp), st=1,32)]
call stft_into(channels, z, freq, time, fs=100.0_dp, &
    window_values=taper, noverlap=16, nfft=32, &
    boundary='none', padded=.false., axis=1, status=st)
if (st /= spectral_status_success) error stop 'STFT failed'
```

A numerical `window_values` array may **not** be combined with `window`,
`window_parameter`, `window_parameter2`, or `window_coefficients`.
When `window_values` is present, its length determines `nperseg`; an explicitly
supplied `nperseg` must equal that length, and the window must fit the available
samples. Window weights must be finite. Incompatible parameter families and
zero-energy normalization return `spectral_status_invalid_input`.

Spectrogram's **omitted-window default** now uses the SciPy legacy Tukey window
with `alpha=0.25`. Explicit `window='tukey'` without a parameter retains the
ordinary `get_window` default of `alpha=0.5`. This distinction also applies
to `spectrogram_into(mode='magnitude'|'complex'|'angle'|'phase')`.

`stft_into` and `spectrogram_into` accept input ranks 1 through 14; the
frequency axis replaces the requested sample axis and a time axis is
appended. `periodogram_into`, `welch_into`, `csd_into`, and `coherence_into`
accept ranks 1 through 15 and preserve the rank. Their transaction semantics
remain unchanged: when a channel fails validation, no caller-provided output
samples are changed and coordinate arrays are returned empty.

## Standalone validation

The script below builds only seven SciFort modules, so it does not need FPM
or external Fortran dependencies. It runs checked and optimized GNU Fortran
builds and compares 96 window configurations against installed SciPy 1.17.

```sh
python tools/run_signal_windowed_estimators.py
```

The numerical cases cover six estimators, eight window configurations,
real/complex signals, 3,264 Fourier/spectral coefficients, frequency grids,
and time grids. `test/test_signal_windowed_estimators.f90` additionally checks
array-window length inference, conflicting options, nonfinite samples,
rank-14 and rank-15 inputs, noncontiguous sections, and unchanged outputs on
failure. The previous STFT boundary test is rerun in the same harness.

## Legacy spectrogram modes (2026-10-28)

The one-dimensional real and complex `spectrogram` now also offer all five
SciPy output modes, including complex STFT coefficients and unwrapped phase.
The rank-independent `spectrogram_into` uses the same implementation. See
[Legacy spectrogram output modes](signal_spectrogram_modes.md).
