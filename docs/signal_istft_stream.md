# Streaming inverse STFT

`scifort_signal` now exports `istft_stream`, `istft_stream_init`,
`istft_stream_push`, and `istft_stream_finish`. This API reconstructs a
signal **one STFT frame at a time**, without allocating an array for the
complete time-frequency input or reconstructed recording. The legacy
`istft` and rank-independent `istft_into` remain unchanged. For synchronized
multichannel reconstruction from `(frequency,channel)` frames, use
[`istft_stream_multi`](signal_stream_multi.md).

The stream reuses the same inverse FFT, overlap-add normalization, and
options as `istft`. Its persistent memory is proportional to `nfft +
nperseg` rather than the number of frames, and each pushed frame returns
only the finalized samples (which can be an empty array). `finish` emits
the remaining samples and removes the final boundary extension.

## Example: process frames incrementally

```fortran
use scifort_signal, only : dp, istft_stream, istft_stream_init, &
    istft_stream_push, istft_stream_finish, istft_status_success
implicit none
type(istft_stream) :: decoder
complex(dp) :: frame(5)
real(dp), allocatable :: chunk(:), sample_times(:)
integer :: status, j

call istft_stream_init(decoder, nfreq=5, status=status, fs=8000.0_dp, &
    window='hann', nperseg=8, noverlap=4, nfft=8)
if (status /= istft_status_success) error stop 'initialization failed'

do j = 1, 1000
    ! Read the next frequency-domain frame from a file or live source.
    ! This is a placeholder; substitute actual STFT coefficients.
    frame = cmplx(0.0_dp,0.0_dp,dp)
    call istft_stream_push(decoder, frame, chunk, status, times=sample_times)
    if (status /= istft_status_success) error stop 'incomplete reconstruction'
    ! Write/play chunk (some initial chunks may be empty).
end do
call istft_stream_finish(decoder, chunk, status, times=sample_times)
if (status /= istft_status_success) error stop 'incomplete reconstruction'
! Write/play the final chunk.
```

Use a `complex(dp), allocatable` output for two-sided complex STFTs and pass
`input_onesided=.false.` to initialization. Complex output is also supported
for one-sided frames and has zero imaginary component, matching `istft`.
The **real** output overload accepts only one-sided frames. Explicit
`window_values` may be passed instead of `window`.

## API and streaming contract

- `istft_stream_init(state, nfreq, status, fs, window, nperseg, noverlap,
  nfft, boundary, scaling, window_values, input_onesided)` checks the
  parameters and constructs a fresh state. The default is `input_onesided
  =.true.`, and defaults for `fs`, window, overlap, FFT length, scaling, and
  boundary follow the legacy inverse STFT. Supply `nfft` explicitly for
  odd-length one-sided FFTs when the bin count is ambiguous.
- `istft_stream_push(state, frame, samples, status, times)` consumes a
  complex rank-1 STFT frame with `size(frame)==nfreq`. `samples` is an
  allocatable rank-1 **real or complex** array, selected by its declared
  type. The optional `times` array contains the corresponding `(i-1)/fs`
  sample coordinates in the complete reconstructed signal.
- `istft_stream_finish(state, samples, status, times)` finalizes the stream,
  flushes its tail and removes the right boundary extension when requested.
  At least one valid `push` is required. Calling `push`/`finish` after
  `finish` is invalid; reinitialize the state to start another stream.
- A frame with the wrong number of bins or any nonfinite value is **rejected
  without advancing the state**. Likewise an invalid `init` preserves any
  previously initialized stream. Rejected operations return empty outputs.
- Status constants match `istft`: `istft_status_success=0`,
  `istft_status_invalid_input=1`, `istft_status_noninvertible=2`.
  With noninvertible overlap normalization, the affected samples are still
  produced (as in the batch `istft` routine); status 2 remains sticky
  through subsequent `push` and `finish` calls. The caller decides whether
  to use the samples.
- With `boundary=.true.`, the implementation deliberately holds back
  enough completed samples to remove **both** edge extensions, even if the
  final frame is not yet known. Results from successive calls concatenate
  to the corresponding batch `istft` output. The output stream is
  **single-channel**; users may keep an independent state for each channel.

## Validation

- `test/test_signal_istft_stream.f90`: 960 batch-versus-stream cases across
  real/complex input choices, odd and even FFT lengths, overlaps, windows,
  normalization modes, boundary options and different frame counts. Also
  exercises user-provided windows, invalid/rejected frames, reinitialization,
  finalization and a 20,000-frame fixed-memory scenario.
- `tools/reference_signal_istft_stream.f90` and
  `test/compare_signal_istft_stream_scipy.py`: 1,350 independent cases
  against SciPy 1.17, covering 17,100 output samples (maximum absolute
  difference about `4.41e-13`). The SciPy 1.17 `nperseg=1` / `boundary=True`
  case is excluded because its boundary removal produces an empty result;
  SciFort's batch and streaming inverses retain the one sample.

Run with FPM when available:

```text
fpm test
```

Without FPM, the repository's standalone script can compile the full signal
module dependency chain and run filtered tests:

```text
python tools/run_signal_standalone.py checked ../signal-build --match istft_stream
```

No new library dependencies are required.
