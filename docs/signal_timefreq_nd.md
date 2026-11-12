# Rank-independent STFT and spectrogram interfaces

The `scifort_signal` module exports Fortran 2018 caller-output generics
`stft_into` and `spectrogram_into`. They wrap SciFort's existing legacy
`stft` and `spectrogram` kernels to process real and complex arrays of
**ranks 1–14** along any selectable signal axis. The frequency dimension
**replaces** the selected input dimension, and a new **time dimension is
appended last**. The result rank is input rank plus one, with a maximum of
Fortran's standard rank 15.

## Usage

```fortran
use scifort_signal, only : dp, stft_into, spectrogram_into, &
    spectral_status_success
real(dp) :: data(12,4), power(9,4,2)
complex(dp) :: coefficients(9,4,2)
real(dp), allocatable :: freq(:), time(:)
integer :: st

! Four real channels, with 12 samples per channel.
call stft_into(data, coefficients, freq, time, fs=100.0_dp, &
    nperseg=6, noverlap=2, nfft=16, axis=1, boundary='none', &
    padded=.false., scaling='psd', status=st)
if (st /= spectral_status_success) error stop 'STFT failed'

call spectrogram_into(data, power, freq, time, fs=100.0_dp, &
    nperseg=6, noverlap=2, nfft=16, axis=1, mode='psd', status=st)
if (st /= spectral_status_success) error stop 'spectrogram failed'
```

### Interfaces

- `stft_into(x, result, frequencies, times, fs, window, nperseg,
  noverlap, nfft, detrend, boundary, padded, scaling, axis, status, ...)` always
  writes a **complex** result; it accepts real or complex input. The
  `scaling` choices are `'spectrum'` (default) and `'psd'`. It supports
  `boundary='zeros'` (default), `'even'` (symmetric reflection), `'odd'`
  (reflection around endpoint values), `'constant'` (repeat endpoints), or
  `'none'` (no extension), and `padded=.true.` (default) or `.false.`.
  The one-dimensional real and complex kernels use the same boundary rules.
  These modes also work with arbitrary `axis=` choices and rank-14 inputs.
- `spectrogram_into(x, result, frequencies, times, fs, window,
  nperseg, noverlap, nfft, detrend, scaling, mode, axis, status, ...)` supports
  `'psd'` (default with real-valued output), `'magnitude'`, `'angle'`,
  and `'phase'` with real-valued output; for `'complex'`, declare a complex
  output array (this is the default for a complex-valued output). `'phase'`
  unwraps along the frequency dimension. The one-dimensional and
  multidimensional modes now share the same spectral kernel (including
  SciPy's real FFT DC/Nyquist phase-branch convention); see
  [Legacy spectrogram output modes](signal_spectrogram_modes.md).
  The default window is Tukey (alpha 0.25), the
  overlap is one eighth of the segment length, and no time-boundary padding
  is used. Scaling is `'density'` (default) or `'spectrum'`.
- A real input produces a one-sided result of length `nfft/2+1` by default;
  pass `return_onesided=.false.` to retain the full two-sided result of length
  `nfft`, with signed FFT-ordered frequencies (negative Nyquist for even
  `nfft`). Complex input always produces a two-sided result regardless of
  the `return_onesided` setting. The `frequencies` and `times` arrays are
  allocated and returned on success.
- The `axis` argument defaults to the **last input dimension**. Positive
  indices are **one-based** (`axis=1` refers to the first axis), unlike
  SciPy's zero-based positive axes. Negative axes count from the end.
- The caller supplies an output of the correct rank and extents. On an
  invalid axis, shape mismatch, failed individual channel, or other error,
  the `result` array remains completely unchanged, while `frequencies` and
  `times` are allocated as zero-length arrays. `status` is optional;
  `spectral_status_success` is zero.
- Noncontiguous input and output array sections are supported. The
  implementation currently stages flattened complex input/output copies to
  ensure all-or-nothing output writes; this uses extra memory for large data.
  Zero-size input extents are currently rejected. Input rank 15 is not
  supported because the legacy STFT representation adds a time axis.

For input with shape `(2, 12, 3)` and `axis=2`, the output has shape
`(2, number_of_frequencies, 3, number_of_time_slices)`. All other input
extents are unchanged. Existing rank-one `stft` and `spectrogram` functions
and the separate centered `short_time_fft` derived type are unchanged.

## Verification

`test/test_signal_stft_boundaries.f90` checks all five STFT boundary
conventions, real/complex inputs, rank-14 data, strided arrays, and preservation
of caller outputs after invalid options. The independent SciPy comparison is
`test/compare_signal_stft_boundaries_scipy.py`, using the separately compiled
`tools/signal_reference/stft_boundaries_driver.f90` executable.

`test/test_signal_timefreq_nd.f90` exercises the public interface for
ranks 2, 3 and 14, real/complex inputs, multiple axis placements,
noncontiguous arrays, STFT and spectrogram outputs, and all-or-nothing
handling of invalid options. The comparison driver
`validation/reference_drivers/reference_signal_timefreq_nd.f90` is reproducibly generated by
`scripts/gen_signal_timefreq_nd_reference.py`.

For an executable linked to SciFort's public module and the driver:

```sh
python test/compare_signal_timefreq_nd_scipy.py /path/to/reference_signal_timefreq_nd
```

This compares 396 cases covering ranks 2, 3, 6 and 14 and 727,408 complex
coefficients with SciPy 1.17, including the five spectrogram modes, window
families, real and complex signals, STFT padding/boundary choices, two
normalizations, and all supported signal axes in the fixture. To avoid
higher-rank padding bugs in SciPy 1.17, the reference applies SciPy's public
rank-one transform to each channel and then checks the expected output
layout. Phase comparisons omit bins with numerically vanishing STFT
magnitude, where phase is undefined, and compare unwrapped phase modulo
integer turns of `2*pi` to accommodate differences in the branch chosen
across such bins. No SciPy code is included at runtime; there are no new
Fortran library dependencies.

## Inverse STFT for rank-2 through rank-15 coefficients

The `scifort_signal` module also exports `istft_into`, the counterpart to
`stft_into`. It takes a complex rank-2 through rank-15 time-frequency
array, with **time and frequency on independently selectable axes**, and
writes a caller-allocated rank-1 through rank-14 real or complex signal. The selected frequency
dimension is replaced by samples and the time dimension disappears.
All other axes retain their positions. By default, time is still the final
input axis, so existing `istft_into` calls are unchanged.

```fortran
use scifort_signal, only : dp, stft_into, istft_into, &
    istft_status_success
real(dp) :: x(16,3), y(16,3)
complex(dp) :: z(5,3,5)
real(dp), allocatable :: frequencies(:), frames(:), sample_times(:)
integer :: st

! Initialize x with signals on the first dimension.
call stft_into(x, z, frequencies, frames, nperseg=8, noverlap=4, &
    nfft=8, axis=1, status=st)
call istft_into(z, y, sample_times, nperseg=8, noverlap=4, &
    nfft=8, axis=1, status=st)
if (st /= istft_status_success) error stop 'inverse STFT failed'
```

The transform axes can also be placed in any order. For example, the
following moves the STFT frame axis from last to first without changing the
reconstructed signal shape:

```fortran
complex(dp) :: z_time_first(5,5,3)
integer :: frame, bin, channel

do channel=1,3
    do bin=1,5
        do frame=1,5
            z_time_first(frame,bin,channel)=z(bin,channel,frame)
        end do
    end do
end do
call istft_into(z_time_first,y,sample_times,nperseg=8,noverlap=4, &
    nfft=8,time_axis=1,freq_axis=2,status=st)
if (st /= istft_status_success) error stop 'inverse STFT failed'
```

The generic supports the following arguments:

`istft_into(zxx, result, times, fs, window, nperseg, noverlap, nfft,
boundary, scaling, axis, status, window_values, input_onesided, time_axis,
freq_axis, window_parameter, window_parameter2, window_coefficients)`

- `zxx` is always complex. `result` is real or complex; its rank is exactly
  `rank(zxx)-1`. All extents except the frequency/sample extent must match.
  `times` is allocated on success and contains `(i-1)/fs` for each sample.
- `time_axis` selects the frame dimension of **`zxx`**; it defaults to the
  final input dimension. `freq_axis` selects the frequency dimension of
  **`zxx`**; it defaults to the dimension corresponding to `axis` after
  deleting `time_axis`. These dimensions must be distinct.
- `axis` selects the sample dimension of **`result`**, defaulting to its
  last dimension unless `freq_axis` is supplied. When `freq_axis` is
  specified and `axis` is absent, the sample dimension is inferred from
  the input frequency position after removing `time_axis`. When both are
  present, they must identify the same surviving dimension.
- Positive axes are **one-based**, unlike Python's zero-based axes; negative
  axes count from the end of the corresponding array (`axis` for the output,
  `time_axis` and `freq_axis` for the input). SciPy's legacy `istft` may put
  the sample dimension at its input **time** location; SciFort instead puts
  the sample dimension at the input **frequency** location. This exactly
  reverses the layout used by SciFort's `stft_into`.
- For real output, `input_onesided` defaults to true (and cannot be false).
  For complex output it defaults to false. Complex output can also
  reconstruct one-sided spectra by setting `input_onesided=.true.`.
  The inverse accepts named windows with optional `window_parameter`,
  `window_parameter2`, and `window_coefficients`, matching forward `stft_into`,
  or an explicit `window_values` array (not both), `scaling='spectrum'` or `'psd'`, and a
  logical `boundary` flag (default true, removing center extensions).
- On bad dimensions/parameters, noninvertible overlap-add normalization,
  invalid axes, or a channel failure, `result` is **unchanged** and `times`
  is allocated to length zero. `status` is `istft_status_success` only if all
  channels reconstruct; other codes are the underlying inverse STFT codes.
  The transaction uses a flattened input copy and private output staging,
  so memory overhead is proportional to the input and output sizes.

`test/test_signal_istft_nd.f90` tests real, complex, one-sided and two-sided
reconstruction, all relative positions of the independent frequency and
frame axes, rank-15 input, strided array sections, and transactional failure
semantics. The independent comparison driver is
`tools/reference_signal_istft_nd.f90`; it is used with
`test/compare_signal_istft_nd_scipy.py` to check SciFort against SciPy 1.17
on 240 cases and 10,008 output samples, including odd windows, FFT padding,
two scalings, and both output types. The maximum absolute error was
`1.76e-14` in the checked GNU Fortran 14.2 build.

For the independently permuted-axes extension, the driver
`tools/reference_signal_istft_axes.f90` and the script
`test/compare_signal_istft_axes_scipy.py` verify **1,296** additional SciPy
1.17 cases across all six orderings of frequency, time and channel axes,
three invocation forms (`freq_axis`/`time_axis`, `axis`/`time_axis`, or all
three), Hann/Hamming/boxcar windows, real/complex inversions, odd and padded
FFTs, one- and two-sided frequency spectra, boundary trimming and PSD or
spectrum scaling. The **32,832** compared samples differ by at most
**1.83e-14**. SciPy 1.17 occasionally reports a *time-vector length* that
does not match the reconstructed array when both axis locations are moved;
in that case the comparison tests SciFort's time vector against the analytic
`(0:nout-1)/fs` grid rather than the inconsistent SciPy coordinate array.
Run with:

```sh
python test/compare_signal_istft_axes_scipy.py /path/to/reference_signal_istft_axes
```

To compile and run just the native inverse-STFT tests using GNU Fortran
(without FPM), the standalone runner accepts a substring filter:

```sh
python tools/run_signal_standalone.py checked ../build-signal --match istft
python tools/run_signal_standalone.py optimized ../build-signal-O3 --match istft
```

The standalone runner still compiles SciFort's transitive Fortran module
dependencies, but it only links and runs matching tests. The full FPM test
suite depends on the FPM package manager and its external dependencies.

For recordings too large for a complete STFT array, see the independent
[streaming inverse STFT API](signal_istft_stream.md), which reconstructs
one frame at a time while retaining only the overlap-add state.

## Parameterized and explicit windows (2026-10-27)

The legacy `*_into` methods now accept four optional trailing keywords:
`window_parameter`, `window_parameter2`, `window_coefficients`, and
`window_values`. These provide named-window shape values and full user-defined
window samples, including `nperseg` inference for the latter. See
[Parameterized windows in legacy spectral estimators](signal_windowed_estimators.md)
for examples, validation rules, and a standalone test command.

## Full two-sided real-input STFT and spectrogram (2026-10-27)

Both public rank-one kernels and the rank-independent `stft_into` and
`spectrogram_into` interfaces accept a final optional keyword
`return_onesided`. For real input it defaults to `.true.`; set
`return_onesided=.false.` for the full `nfft` frequency bins with negative
frequencies and no one-sided PSD bin doubling. The argument is supported for
all spectrogram output modes (`'psd'`, `'complex'`, `'magnitude'`, `'angle'`,
`'phase'`), for any selected axis, and for real or complex arrays. Complex
input is always two-sided, even when `.true.` is explicitly requested.

For example, if `signal` is a `real(dp)` array of shape `(2,16)` with the
sample dimension last, a complex output `z(2,8,3)` can receive a two-sided
8-point STFT along that dimension:

```fortran
call stft_into(signal, z, frequencies, times, nperseg=8, nfft=8, &
    noverlap=4, axis=2, boundary='none', padded=.false., &
    return_onesided=.false., status=st)
```

The native `test/test_signal_twosided.f90` also checks stride handling,
frequency parity, rank 14 input/15 output, and no partial writes on invalid
shapes. The independent SciPy 1.17 suite runs alongside four earlier
regression executables with both checked and optimized GNU Fortran:

```sh
python tools/run_signal_twosided.py
```
