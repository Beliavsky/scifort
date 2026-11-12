# Multichannel streaming STFT, spectrogram and inverse STFT

`scifort_signal` provides synchronized multichannel variants of the incremental
short-time Fourier transform routines. They accept one state for a fixed number
of channels and return results for every channel on the same time grid. The
single-channel `stft_stream` and `istft_stream` interfaces remain unchanged.

| Operation | State | Routines |
|---|---|---|
| Forward STFT | `stft_stream_multi` | `stft_stream_multi_init`, `_push`, `_finish` |
| Power spectrogram | `stft_stream_multi` | `spectrogram_stream_multi_init`, `_push`, `_finish` |
| Inverse STFT | `istft_stream_multi` | `istft_stream_multi_init`, `_push`, `_finish` |
| Frequency bins | `stft_stream_multi` | `stft_stream_multi_frequencies` |

The full routine names have prefixes as shown above, for example
`spectrogram_stream_multi_push`.

## Array shapes and axis conventions

- Forward input: `(sample, channel)` (`real(dp)` or `complex(dp)`). The
  first dimension can contain an arbitrary number of consecutive samples,
  including zero. The **second dimension must equal `nchannels`** supplied
  at initialization, including when the input chunk is empty.
- STFT output: `complex(dp)` array `(frequency, channel, frame)`.
- Power spectrogram output: `real(dp)` array `(frequency, channel, frame)`.
- Inverse input: `complex(dp)` matrix `(frequency, channel)` containing
  **one frame** for each channel.
- Inverse output: `real(dp)` or `complex(dp)` matrix `(sample, channel)`.
- Every call can return shared frame/sample coordinates through the optional
  allocatable `times(:)` argument. `stft_stream_multi_frequencies` returns
  a common frequency grid. An empty output is valid when no samples or
  frames are finalized by a call.

Fortran 2018 assumed-shape arguments allow array sections with non-unit
strides, including noncontiguous sample and channel selections. No channels
are added or dropped after initialization. These wrappers intentionally fix
sample, frequency and channel positions; use `stft_into`, `spectrogram_into`,
or `istft_into` when batch arrays require selectable axes or higher rank.

## Example: real multichannel STFT into streaming inverse STFT

```fortran
use scifort_signal, only : dp, stft_stream_multi, istft_stream_multi, &
    stft_stream_multi_init, stft_stream_multi_push, stft_stream_multi_finish, &
    istft_stream_multi_init, istft_stream_multi_push, istft_stream_multi_finish, &
    spectral_status_success, istft_status_success
implicit none
type(stft_stream_multi) :: analyzer
type(istft_stream_multi) :: synthesizer
real(dp) :: x(32,3)
real(dp), allocatable :: samples(:,:), time(:)
complex(dp), allocatable :: frames(:,:,:)
integer :: status, k

! Fill x(:,1:3) from three synchronized real signals.
x = 0.0_dp
call stft_stream_multi_init(analyzer,nchannels=3,nperseg=64,status=status, &
    window='hann',noverlap=48,nfft=64,fs=48000.0_dp)
if (status/=spectral_status_success) error stop 'forward init'
call istft_stream_multi_init(synthesizer,nchannels=3,nfreq=33,status=status, &
    nperseg=64,noverlap=48,nfft=64,fs=48000.0_dp)
if (status/=istft_status_success) error stop 'inverse init'

call stft_stream_multi_push(analyzer,x,frames,status,time)
if (status/=spectral_status_success) error stop 'forward push'
do k=1,size(frames,3)
    call istft_stream_multi_push(synthesizer,frames(:,:,k),samples,status)
    if (status/=istft_status_success) error stop 'inverse push'
    ! Consume all completed samples in samples(:,1:3).
end do
call stft_stream_multi_finish(analyzer,frames,status,time)
if (status/=spectral_status_success) error stop 'forward finish'
do k=1,size(frames,3)
    call istft_stream_multi_push(synthesizer,frames(:,:,k),samples,status)
    if (status/=istft_status_success) error stop 'inverse push'
end do
call istft_stream_multi_finish(synthesizer,samples,status)
if (status/=istft_status_success) error stop 'inverse finish'
```

In a real application, call `stft_stream_multi_push` repeatedly with
successive input chunks and pass each emitted spectral frame into the
inverse state. An inverse stream may emit no samples until enough frames
have arrived. The example initializes an input shorter than `nperseg`;
finalization emits the padded forward frames, matching the single-channel
streaming conventions.

## Parameters and memory

`stft_stream_multi_init` accepts the same options as `stft_stream_init`,
with a required positive `nchannels` inserted before `nperseg`. The
spectrogram and inverse initializers follow their respective single-channel
signatures in the same way. `input_complex=.true.` is needed when forward
chunks are complex; complex input always produces two-sided spectra.
The real-valued inverse output overload requires `input_onesided=.true.`;
two-sided inverse spectra require the complex-valued overload.

The public API has only one multichannel state, internally managing one
small rolling transform/overlap state per channel. Persistent memory is
`O(nchannels * (nperseg + nfft))`, independent of the total samples and
frames processed. Each call allocates its own output, proportional to the
number of newly completed frames or samples. The implementation does not
allocate the complete input recording, complete STFT cube, or complete
reconstructed signal. Cross-channel spectra, beamforming, channel-varying
windows and independently sampled channels are outside this interface;
`multitaper_stream` provides separate cross-spectral functionality.

Window options, padding, detrending and normalization follow the existing
single-channel stream implementations. Forward STFT offers `boundary='zeros'`
or `'none'` and optional final `padded`; spectrogram mode has no boundary
extension or final padding and defaults to constant detrending. Inverse
STFT `boundary` is logical and controls trimming of boundary padding.

Invalid initialization leaves an existing stream unchanged. Wrong channel
counts, nonfinite chunks/frames, wrong real/complex output choices and
calls after finalization are rejected before advancing any channel. Error
outputs are allocated with zero elements. The inverse returns
`istft_status_noninvertible` when **any channel** contains finalized samples
whose overlap-add normalization is insufficient; the diagnostic remains
sticky as in the scalar inverse stream. On valid input, all channels have
identical frame/sample counts and times.

## Validation

`test/test_signal_stream_multi.f90` compares the multichannel routines
against independent single-channel SciFort streams over 74 deterministic
cases, including varied overlap/FFT lengths, real and complex signals,
spectrogram PSDs, irregular chunk sizes, noncontiguous sections, incorrect
inputs, an inverse NOLA warning, round-trip reconstruction and a 10,000-sample
four-channel long-stream stress test.

`tools/reference_signal_stream_multi.f90` and
`test/compare_signal_stream_multi_scipy.py` provide independent numerical
comparisons against SciPy 1.17. They cover 128 forward STFT/spectrogram
configurations (12,304 spectral values; max absolute difference `1.33e-15`)
and 128 inverse STFT configurations (4,352 output samples; max difference
`4.78e-14`). The inverse comparison checks each channel against SciPy's
one-dimensional inverse routine, avoiding time-coordinate ambiguities in
SciPy's legacy higher-rank inverse convenience interface.

When FPM is installed run `fpm test`. The provided standalone runner works
with GNU Fortran:

```text
python tools/run_signal_standalone.py checked build_signal --match stream_multi
```

Compile the reference driver against the SciFort signal objects and supply
its executable path to `test/compare_signal_stream_multi_scipy.py`. These
multichannel wrappers are MIT-licensed and add no dependencies.
