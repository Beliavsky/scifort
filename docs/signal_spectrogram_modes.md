# Legacy spectrogram output modes

The real and complex rank-one `scifort_signal::spectrogram` functions implement
all five output modes of `scipy.signal.spectrogram` (tested against SciPy 1.17):

| `mode` | Mathematical result | SciFort result component |
| --- | --- | --- |
| `'psd'` (default) | Power spectral density / spectrum | `result%power` (real) |
| `'complex'` | Windowed, detrended **complex STFT** coefficients | `result%values` (complex) |
| `'magnitude'` | Absolute value of those complex coefficients | `result%power` (real) |
| `'angle'` | Principal complex argument, in radians | `result%power` (real) |
| `'phase'` | Complex argument unwrapped **along frequency bins** | `result%power` (real) |

`mode` is an optional **trailing** character argument. Previously compiled
call sites and existing positional source calls are not altered. The default
`%power` output remains unchanged. Because Fortran cannot change the return
type based on a runtime string, `%values` was added as an allocatable component
to the existing `spectrogram_result` type. It is allocated only for
`mode='complex'`; in that case `%power` is allocated with shape `(0,0)`.
For all other valid modes, `%values` stays unallocated. Invalid options
return `spectral_status_invalid_input` and empty coordinate/power arrays.

```fortran
use scifort_signal, only : dp, spectrogram, spectrogram_result, &
    spectral_status_success
real(dp) :: signal(256)
type(spectrogram_result) :: s, p

! Initialize signal first.
s = spectrogram(signal, fs=4000.0_dp, nperseg=32, noverlap=4, &
                nfft=64, window='hann', mode='complex')
if (s%status /= spectral_status_success) error stop 'complex spectrogram'
! s%values has shape (33, n_frames) for the default one-sided real FFT.

p = spectrogram(signal, fs=4000.0_dp, nperseg=32, noverlap=4, &
                nfft=64, window='hann', mode='phase')
if (p%status /= spectral_status_success) error stop 'phase spectrogram'
! p%power has shape (33, n_frames), with phase values in radians.
```

The default window when omitted is **periodic Tukey, alpha=0.25**, not the
standard `window='tukey'` default alpha=0.5. When omitted, the overlap is
`nperseg/8` (integer division) and `nperseg` is inferred as before. The
`scaling='density'` and `'spectrum'` modes apply to PSD and to the amplitude
normalization used by the four other modes. For one-sided real FFTs,
only PSD bins are doubled; complex STFT coefficients are **not** doubled.
Full two-sided spectra retain NumPy FFT bin order. For complex signals, the
spectrum is always two-sided.

For real input at negative-valued DC and even-length Nyquist bins, the
`angle`/`phase` branch matches SciPy's signed-zero convention: a complete
complex FFT selects `-pi`, while a one-sided real FFT selects `+pi`. This is
important when the phase is subsequently unwrapped. Unwrapping traverses
**frequency**, not successive time frames.

## Multidimensional signals

`scifort_signal::spectrogram_into` already exposes all five modes for input
ranks 1–14. It now shares the same kernel as rank-one `spectrogram`, including
window defaults and phase handling. Use a real output array for `psd`,
`magnitude`, `angle`, or `phase`, and a complex array for `complex`. The
selected axis is one-based (or negative, with `-1` denoting the last input
axis). Frequency replaces the selected input axis and time is appended last.
The output must be correctly sized; the routine does not partially modify a
caller-owned output on errors. Both real and complex input types are accepted.

## Validation

```sh
python tools/run_signal_spectrogram_modes.py
```

The script compiles the isolated spectral modules in GNU Fortran checked and
`-O3` modes, runs the native one-dimensional and rank-independent regression
tests, and compares **440 independent SciPy 1.17 cases** across all five
modes. The comparison checks 58,254 spectral values and both coordinate grids;
maximum observed absolute difference was `8.31e-14`. Tested choices include
real/complex data, default/Hann/parameterized Gaussian/explicit windows,
odd/even FFT lengths, one-sided/two-sided output, three detrending settings,
both scaling settings, and omitted sizing arguments. SciPy is needed only to
run this optional comparison, not to use the Fortran library.

The complete 91-module checked GNU Fortran public-module build also passed
`test_signal_timefreq_nd.f90`. The older independent multidimensional
STFT/spectrogram reference suite passed 396 SciPy cases and 727,408
coefficients with maximum absolute error `2.13e-12`. The full `fpm test`
suite has not been run in this environment.
