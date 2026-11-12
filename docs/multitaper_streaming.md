# Streaming multitaper spectral estimation

The SciFort `multitaper_stream` type incrementally estimates DPSS multitaper
cross-spectral matrices from incoming real or complex data. Unlike the batch
`multitaper_spectral_matrix_spectrogram`, the estimator **does not retain the
whole time series**. It caches the periodic DPSS tapers and retains only an
incomplete segment, FFT work arrays, channel weights, and a running mean
spectral matrix. This is a SciFort extension, not a SciPy function.

The stream also retains a reusable Fourier transform plan, so FFT chirps,
roots and scratch storage are not re-created for every DPSS taper or window.
`reset()` preserves both tapers and the plan. `initialize()` clears the old
plan before preparing one for the new FFT length. This changes neither FFT
sign/normalization nor the output memory bound; a stream state must not be
shared concurrently between threads.

## Example

```fortran
use scifort_signal, only : dp, multitaper_stream, &
    multitaper_spectral_matrix_spectrogram_result, multitaper_matrix_result, &
    spectral_status_success
implicit none
type(multitaper_stream) :: stream
type(multitaper_spectral_matrix_spectrogram_result) :: new_windows
type(multitaper_matrix_result) :: running_mean
integer :: status
real(dp) :: chunk1(150,2), chunk2(125,2)

! Fill chunk1 and chunk2 with the next samples from two channels.
call stream%initialize(3.5_dp,64,2,status,noverlap=32,kmax=5, &
    fs=100.0_dp,nfft=128,weighting='adaptive',detrend_mode='linear')
if (status /= spectral_status_success) error stop 'invalid stream settings'
new_windows = stream%push(chunk1)
! Process new_windows%power(:, :, channel_x, channel_y), if desired.
new_windows = stream%push(chunk2)
running_mean = stream%mean_matrix()
! To restart while reusing the cached tapers:
call stream%reset()
```

## API

`call stream%initialize(nw,nperseg,nchannels,status,noverlap,kmax,fs,nfft, &
 detrend,weighting,complex_input,detrend_mode)` has the same window conventions and defaults
as `multitaper_spectral_matrix_spectrogram`:

- `nw`: time-half-bandwidth product (`0 < nw < nperseg/2`).
- `nperseg`: integer window size, at least two; `nchannels`: positive channel count.
- `noverlap`: zero through `nperseg-1`, default `nperseg/2`; `kmax` defaults to
  `min(nperseg,max(1,int(2*nw)-1))`.
- `fs`: samples per second (default 1); `nfft`: transform length, at least
  `nperseg` (default `nperseg`).
- `detrend`: remove each window's per-channel mean (default `.true.`).
- `detrend_mode`: optional `'none'`, `'constant'` (default), or `'linear'`;
  overrides logical `detrend` when present. Linear removes each complete
  window's fitted intercept and slope **independently for every channel**
  before applying DPSS tapers, including in adaptive weighting.
- `weighting`: `equal` (default), `eigen`, or `adaptive`.
- `complex_input=.true.` selects complex two-sided spectra; the default uses
  real one-sided spectra. The input type at `push` must match this choice.

`result = stream%push(chunk,axis)` accepts real/complex rank-one chunks when
`nchannels=1` and real/complex rank-two chunks for all channel counts.
Rank-two chunks normally have shape `(new_samples,nchannels)`. Set the
one-based `axis=2` to supply `(nchannels,new_samples)`. Empty chunks are valid.
Each returned `multitaper_spectral_matrix_spectrogram_result` contains:

- `power(frequency,new_window,channel_x,channel_y)`, complex DPSS cross-spectra
  satisfying `P_xy = conjg(FFT(x))*FFT(y)`;
- `times(new_window)`, the absolute window-centre times since the last
  initialization/reset, and `frequencies(frequency)`;
- `status = spectral_status_success` for valid input. If the current chunk
  finishes no complete windows, the time dimension of `power` is zero.

The number of new windows depends only on cumulative input length, never on
chunk boundaries. Only complete windows are emitted; incomplete final samples
are discarded if the stream is abandoned. No padding is introduced at the
beginning or end. Cross-channel adaptive weights are computed independently
per channel and combined geometrically, just like the batch estimator.
The linear-detrending mode uses centered least-squares coordinates, and
**does not** fit one trend over the full input stream. Use the whole-record
multitaper functions when one global trend estimate is wanted. The stream's
chosen mode survives `reset()`, but can be changed by `initialize()`. Unknown
modes produce `spectral_status_invalid_input` at initialization.

`running_mean = stream%mean_matrix()` returns the average of the cross-spectral
matrices of **all complete windows** processed since initialization/reset.
It is invalid before the first complete window, and its `power` shape is
`(frequency,channel_x,channel_y)`. The mean is not an average of the input
samples or a full-record PSD estimate.

`stream%window_count()` and `stream%sample_count()` return `integer(int64)`
counters. `call stream%reset()` clears both counters and accumulated spectra
but keeps the configuration, tapers, and work buffers. Reinitialization
replaces those settings. Invalid input (including mismatched channel counts,
wrong real/complex type, nonfinite data, or invalid axis) returns an invalid
status without consuming any samples. The methods are stateful and should not
be called concurrently on the same stream object.

## Streaming local noise-floor estimates

The same state also supplies `stream%push_noise_floor(chunk,median_bins,axis)`.
The output type is `multitaper_noise_floor_spectrogram_channels_result` with
arrays `power`, `background`, and `whitened`, each shaped
`(frequency,new_complete_window,channel)`, plus time and frequency axes and a
status. The median acts only over frequency, separately for each window and
channel, using reflected one-sided or periodic two-sided boundary handling.
Odd `median_bins` values from 3 through the number of frequency bins are
valid. For rank-1 chunks, omit `axis` and initialize a one-channel stream.

```fortran
use scifort_signal, only : multitaper_noise_floor_spectrogram_channels_result, &
    multitaper_noise_floor_channels_result
type(multitaper_noise_floor_spectrogram_channels_result) :: next_floor
type(multitaper_noise_floor_channels_result) :: mean_floor
! Initialize stream as above, then feed a real rank-2 sample chunk:
next_floor = stream%push_noise_floor(chunk1,7)
mean_floor = stream%mean_noise_floor(7)
```

`mean_noise_floor()` applies the local median to the *mean* of completed
window PSDs. This is not the average of time-frequency `whitened` values;
its result has `power(frequency,channel)`, `background(frequency,channel)`,
`whitened(frequency,channel)` and no artificial time coordinate. Calling it
before any complete window gives invalid-input status with empty outputs.
No additional historical samples are retained. Both `push` and
`push_noise_floor` advance **the same stream**, so submit each chunk only
once. Invalid median widths do not consume data.

## Complexity and verification

Internal memory is bounded by the segment size, channel count, taper count and
FFT length, **not the total number of processed samples**. Returning a batch
of all newly completed windows necessarily allocates a result proportional to
that batch's size; callers can choose small chunks when they need tight memory
bounds or use `mean_matrix()` for history-free summaries.

`test/test_signal_multitaper_stream*.f90` tests chunk-boundary independence,
real/complex spectra, selected axes, all weighting methods, odd and even FFT
sizes, one- and multichannel arrays, empty chunks, running means, reset, and
invalid inputs against the existing batch implementation. An additional
`python test/compare_signal_multitaper_stream.py` checks 80 randomized
streamed cross-spectra against an independent SciPy DPSS/NumPy FFT reference.

## Independent linear-detrending reference

`test/test_signal_multitaper_whole_stream_detrend.f90` checks chunking,
real/complex signals, mean spectra, all three taper-weighting methods,
invalid arguments, mode precedence, and preserved configurations after reset.
`tools/signal_reference/compare_signal_multitaper_whole_stream_detrend.py`
independently assembles reference spectral matrices using SciPy's
`signal.detrend(type='linear')`, periodic DPSS tapers and NumPy FFTs. This
checks whole-record matrices as well as *independently detrended complete
streaming windows* for equal and eigen weighting. Adaptive weighting is
covered by the native batch-versus-stream regression.
