# Streaming forward STFT and spectrogram

The public `scifort_signal` module exposes `stft_stream` with incremental
forward analysis, complementing the existing streaming inverse transform:

- `stft_stream_init`, `stft_stream_push`, `stft_stream_finish`
- `spectrogram_stream_init`, `spectrogram_stream_push`,
  `spectrogram_stream_finish`
- `stft_stream_frequencies`

A caller passes **successive time-domain chunks**. Each `push` returns the
complete frequency-domain frames available at that point (possibly zero
frames), and `finish` emits the final frames after applying the requested
right-boundary extension and final padding. Unlike `stft`/`spectrogram`, no
allocation for an entire recording is required. The internal ring buffer,
window and FFT workspace have persistent size `O(nperseg + nfft)`; each
returned allocation grows only with the frames emitted by that call.

## Example: compute and store STFT frames

```fortran
use scifort_signal, only : dp, stft_stream, stft_stream_init, &
    stft_stream_push, stft_stream_finish, stft_stream_frequencies, &
    spectral_status_success
implicit none
type(stft_stream) :: encoder
real(dp) :: samples(1000)
real(dp), allocatable :: frequencies(:), times(:)
complex(dp), allocatable :: frames(:,:)
integer :: i, st

samples = [(sin(0.01_dp*real(i,dp)),i=1,1000)]
call stft_stream_init(encoder,nperseg=64,status=st,fs=1000.0_dp, &
    window='hann',noverlap=48,nfft=128,boundary='zeros', &
    padded=.true.,scaling='spectrum')
if (st/=spectral_status_success) error stop 'STFT init'
call stft_stream_frequencies(encoder,frequencies,st)
if (st/=spectral_status_success) error stop 'STFT bins'

do i=1,size(samples),113
    call stft_stream_push(encoder,samples(i:min(i+112,size(samples))), &
        frames,st,times)
    if (st/=spectral_status_success) error stop 'STFT push'
    ! Save/process frames(:,j) at times(j).
end do
call stft_stream_finish(encoder,frames,st,times)
if (st/=spectral_status_success) error stop 'STFT finish'
! Process final frames; `finish` can be called only once.
```

Call `spectrogram_stream_init` in place of `stft_stream_init` and
`spectrogram_stream_push` / `_finish` for **real-valued** power-spectrum
matrices. The corresponding power-normalization arguments are
`scaling='density'` (the default) or `'spectrum'`, with constant detrending,
no boundary extension and no final padding by default. Real signals return
one-sided power with the appropriate doubled interior frequency bins;
complex signals return two-sided power. For complex input call initialization
with `input_complex=.true.` and pass `complex(dp)` chunks. All frame matrices
use `(frequency,frame)` order and each returned `times(:)` has the same number
of entries as the second matrix dimension.

## Interface and conventions

- **Window length:** `nperseg` is required at initialization, rather than
  inferred from a recording length that is not yet known. Default overlap is
  `nperseg/2` for STFT and `nperseg/8` for spectrogram. `nfft` defaults to
  `nperseg`; zero padding the FFT is supported.
- **Analysis:** default periodic window `'hann'` for STFT and periodic Tukey
  with alpha 0.25 for spectrogram, matching the batch spectrogram and SciPy's
  legacy spectrogram default. Explicit `window='tukey'` uses alpha 0.5,
  matching `get_window('tukey',nperseg)`. Supply `window_values(:)` instead
  of a name for custom weights or another Tukey alpha, for example
  `get_window('tukey',nperseg,parameter=0.25_dp)`.
- **Detrending:** `'none'`, `'false'`, `'constant'`, or `'linear'` are supported.
  The stream uses SciFort's existing detrend-then-window-then-FFT convention.
- **STFT edges:** `boundary='zeros'` (default) adds `nperseg/2` zeros on each
  side, and `boundary='none'` does not extend either edge. `padded=.true.`
  (default) completes the last hop. The spectrogram API follows legacy
  SciPy's no-boundary, no-padding convention. These routines do not yet
  support `'even'`, `'odd'`, or `'constant'` boundary modes.
- **Scaling:** STFT uses `'spectrum'` or `'psd'` (with amplitude units);
  spectrogram uses `'density'` or `'spectrum'` (with power units). Custom
  windows and zero-padded FFT lengths are supported. Real inputs are
  one-sided unless `return_onesided=.false.`; complex inputs are always
  two-sided, consistent with SciPy's handling of complex time series.
- **Output:** `stft_stream_push` / `_finish` return allocatable `complex(dp)`
  matrices, and `spectrogram_stream_push` / `_finish` return allocatable
  `real(dp)` matrices. `stft_stream_frequencies` returns bin centers once
  at initialization, using negative-frequency ordering for two-sided FFTs.
- **Error handling:** all routines return `spectral_status_success` (0) or
  `spectral_status_invalid_input` (1). Malformed initialization preserves
  the previous state; nonfinite chunks, wrong input type, or calls after
  `finish` return empty outputs without advancing the stream. Zero-length
  pushes are allowed. `finish` with no samples produces no frames.
- **Scope:** this `stft_stream` interface analyzes one signal per state.
  For synchronized channels, use the new `stft_stream_multi` and
  `spectrogram_stream_multi_*` interfaces in
  [multichannel streaming](signal_stream_multi.md). The rank-independent
  batch `stft_into` and `spectrogram_into` interfaces are unchanged. Chunk
  boundaries may be arbitrary, but supply chunks in sample order.

For streams shorter than `nperseg`, the requested window size stays fixed;
SciPy's whole-array convenience functions may instead shorten an explicitly
requested `nperseg` to fit the available input length. For direct equivalence
in that case, choose matching window lengths explicitly.

## Validation

`test/test_signal_stft_stream.f90` runs **1,128** deterministic comparisons
against SciFort's existing `stft` and `spectrogram` batch implementations.
Additional tests cover rejected inputs, initialization preservation,
custom windows, real two-sided output, constant-sized rolling state across
20,000 samples, and combined forward/inverse STFT reconstruction of real
and complex signals.

The added window regressions cover the omitted spectrogram window, explicit
Tukey names, and custom Tukey(0.25) weights with real and complex input,
odd/even segment and FFT lengths, both power scalings, and irregular chunks.

`tools/reference_signal_stft_stream.f90` and
`test/compare_signal_stft_stream_scipy.py` perform **768 independent
comparisons**, including omitted windows and explicit `boxcar`, `hann`, and
`tukey` names. Checked (`-O0 -fcheck=all`) and optimized (`-O3`) GNU Fortran
17.0.0 experimental builds pass against SciPy 1.15.3 over **99,396 spectral
bins**, with maximum absolute error `4.77e-15` in both builds. The cases exercise
real and complex data, chunk boundaries, STFT and spectrogram mode, FFT lengths,
detrending, edges, and scaling.

With FPM installed, use `fpm test`. Without FPM:

```text
python tools/run_signal_standalone.py checked ../signal-build --match stft_stream
```

Compile `tools/reference_signal_stft_stream.f90` with the signal module
objects and pass its executable path to
`test/compare_signal_stft_stream_scipy.py` to rerun the SciPy comparisons.
There are no new library dependencies.
