## Canonical spectral coherence modes (SciFort extension, 2026-10-29)

The `canonical_coherence_matrix` and `canonical_coherence_tf_matrix` APIs
return the sorted squared canonical correlations between channel groups,
optionally conditioning on another group, together with the existing
Wilks-determinant aggregate. Frequency-band and direct real/complex
multitaper wrappers are also available. See
[Canonical spectral coherence modes](signal_canonical_coherence.md) for
mathematics, regularization, examples, and independent NumPy SVD tests.

## Multiple coherence (SciFort extension, 2026-10-29)

The new `multiple_coherence_matrix`, `multiple_coherence_tf_matrix`,
`multitaper_multiple_coherence`, and `multitaper_tf_multiple_coherence`
functions compute joint multichannel spectral multiple coherence, conditional
power, and explained power. They reuse the precision-matrix solver and its
singularity, ridge, and shrinkage rules. See
[Multiple coherence](signal_multiple_coherence.md) for mathematics, examples,
limitations and NumPy Schur-complement regression coverage.

## 2026-10-28: Mixed-type one-dimensional convolution

The rank-one `convolve`, `correlate`, `fftconvolve` and `oaconvolve`
interfaces now support real/complex and complex/real operand pairs without
manual casts. `choose_conv_method` accepts the same mixed input types.
All modes, both operand orders, and conjugation of the second correlation
input match SciPy. See `docs/signal_mixed_convolution.md` for examples,
interfaces and reproducible numerical comparisons.

## Legacy spectrogram output modes (2026-10-28)

Real and complex one-dimensional `spectrogram` now supports all SciPy legacy
`mode=` options (`psd`, `complex`, `magnitude`, `angle`, `phase`), consistent with
existing rank-independent `spectrogram_into` modes. See
[Legacy spectrogram output modes](signal_spectrogram_modes.md) for the
`%power` versus `%values` result convention and validation.

## FFT performance for arbitrary lengths

SciFort's shared `scifort_signal_spectral::forward_fft` uses the forward
DFT convention `exp(-2*pi*i*j*k/n)`, without normalization. Power-of-two
lengths use the original iterative radix-2 FFT. Non-power-of-two lengths
at most 16 use the direct DFT to avoid setup overhead. Larger composite
lengths with no prime factors greater than seven use a **mixed-radix
Cooley-Tukey FFT** (radices 2, 3, 5, and 7). All other larger lengths,
including prime lengths, use **Bluestein's chirp-convolution FFT**, with
an internal power-of-two transform of length at least `2*n-1`.
Both accelerated paths are `O(n log n)`, replacing the previous quadratic
DFT fallback; their temporary storage is `O(n)`. Inverse transforms
continue to use conjugation and division by `n` at their call sites,
so API behavior is unchanged.

This implementation is dependency-free and shares the spectral FFT used
by STFT, periodogram/Welch, CWT/CZT, convolution and Fourier resampling.
The original `forward_fft(values)` interface remains unchanged.

### Reusable FFT plans for repeated transforms

For multiple transforms of one length, `fft_plan` caches either Bluestein's
chirp and transformed convolution kernel or the mixed-radix root table, as
well as reusable scratch space. A plan is set up once rather than
recomputing the transform-specific data for every call:

```fortran
use scifort_kinds, only : dp
use scifort_signal, only : fft_plan, forward_fft_planned, &
    spectral_status_success

type(fft_plan) :: plan
integer :: code
complex(dp) :: samples(1009)

call plan%initialize(size(samples), code)
if (code /= spectral_status_success) error stop 'FFT plan initialization failed'
! Fill samples before executing, and refill between repeated transforms.
call forward_fft_planned(samples, plan, code)
if (code /= spectral_status_success) error stop 'FFT plan execution failed'
call plan%inverse(samples, code)  ! inverse is normalized by 1/n
call plan%clear()                 ! optionally release cached storage
```

The same methods work for small direct DFTs and radix-2 sizes. Mixed-radix
plans allocate two `n`-element scratch arrays plus an `n`-element root table
and do not allocate during `execute`. Their radix-3/5/7 stages now
compute symmetric-pair sum/difference butterflies rather than summing all
`r` terms afresh at every output bin. This preserves the same forward FFT
sign, natural-order output and inverse `1/n` convention, and is transparent
to callers. Bluestein plans retain their padded
convolution buffers. Reinitializing a plan for another length replaces the
previous algorithm's cached arrays. `plan%length()`
returns zero for an uninitialized or cleared plan. Invalid lengths or failed
plan allocations leave an existing plan intact; an invalid or mismatched
transform leaves input values unchanged and returns an error status.
`forward_fft_planned` is an **opt-in SciFort extension**, not a public
`scipy.signal` function. Plans are not global; each has mutable scratch
storage and therefore **must not be used concurrently by multiple threads**.
Use one plan per concurrent worker. Plans support strided complex arrays.

Internal reusable plans are automatically prepared in real-valued batch
periodogram/Welch/CSD/STFT/spectrogram and forward/inverse streaming STFT,
and now also the multidimensional FFT convolution and overlap-add kernels.
`fftconvolve_nd` and `oaconvolve_nd` prepare one plan per selected axis and
reuse it across every column-major line, operand, inverse transform, and
(overlap-add) tile. These kernels use a cost-aware choice between power-of-two and selected
2/3/5/7-smooth padding; plan reuse avoids setup, but speedups depend on
transform shape and radix cost. The rank-one and rank-2/3 convolution interfaces also reuse plans.

Every multitaper whole-record/windowed/harmonic estimator now prepares one
plan for the chosen `nfft` and reuses it across tapers, channels, and windows.
`multitaper_stream` owns its FFT plan; `reset()` retains it, while a new
`initialize()` replaces it. The persistent streaming memory bound remains
independent of recording duration. There are no global caches and no API
changes; do not share a single streaming state across simultaneous workers.

Additional regression `test/test_signal_fft_plan_reuse_consumers.f90` covers
rank-15 convolution and stream reset/reinitialization for prime-length and
mixed-radix transforms. Independent SciPy checks include the existing 1,008
convolution/correlation comparisons, 232 tiled overlap-add comparisons, and
120 noise-floor, 80 streamed-spectra, 90 spectrogram plus 100 harmonic-test
multitaper comparisons. In the local optimized streaming benchmark at
`nfft=263`, reusing a plan reduced execution time by roughly 30% (~1.4x
speedup); for `nfft=315`, timings were nearly identical. The benchmark is
`benchmarks/bench_signal_multitaper_fft_plan_reuse.f90`.

Validation: `test/test_signal_fft_bluestein.f90` checks the original FFT.
`test/test_signal_fft_plans.f90` checks repeated planned/unplanned agreement,
inverse transforms, strided arrays, invalid initialization, and rejection
without modifying inputs. `tools/signal_reference/compare_signal_fft_plans.py`
compares 198 cases and 297,189 coefficients with SciPy 1.17 (maximum
magnitude-scaled difference `8.96e-13`). The optimized benchmark
`benchmarks/bench_signal_fft_plans.f90` measures plan reuse relative to
repeated unplanned calls; speedups depend on compiler and hardware.
`benchmarks/bench_signal_fft_non_power_two.f90` retains the original direct
DFT comparison. `test/test_signal_fft_mixed_radix.f90` covers 38 lengths,
repeated forward/inverse transforms, strided input, and invalid operations.
`tools/signal_reference/compare_signal_fft_mixed_radix.py` checks 246 SciPy
FFT comparisons (321,702 complex coefficients; largest magnitude-scaled
difference `3.77e-14`). `benchmarks/bench_signal_fft_mixed_radix.f90`
measures planned and unplanned smooth-length transforms. In the local GNU
Fortran 14.2 `-O3` benchmark, reusable mixed-radix plans were roughly
1.4-2.5 times faster than the previous planned Bluestein implementation
for selected composite sizes between 225 and 5040. A later butterfly
optimization improves planned transform time over that mixed-radix baseline
by approximately 1.36-2.60x on ten measured 2/3/5/7-smooth lengths.
`test/test_signal_fft_butterflies.f90` checks 66 additional native cases;
`tools/signal_reference/compare_signal_fft_butterflies.py` checks 138
independent SciPy 1.17 cases (206,016 complex values) with maximum
magnitude-scaled error `1.78e-15` in the checked build. The previous
246-case FFT reference remains green. Clean public builds
passed 16 FFT-related and 16 neighboring regression programs in both checked
and optimized configurations. Performance varies with transform length,
compiler, and CPU.

# SciPy-style signal algorithms

```fortran
use scifort_signal
```

The first `scifort_signal` release adds seven independent, MIT-licensed
numerical routines corresponding to functions in `scipy.signal`. These
algorithms use real(dp) 1D input unless complex support is noted below.
There is no new FPM dependency. Unlike NumPy/SciPy indices, all sample
indices and peak bases in this Fortran API are **one-based**.

## ShortTimeFFT: invertible centered time-frequency analysis

The public derived type `short_time_fft` covers the central numerical
functionality of `scipy.signal.ShortTimeFFT`: centered window slices, real or
complex data, inversion using a canonical dual window, `stft_detrend`,
spectrograms, cross-spectrograms, and dual-window construction. It operates on
rank-one through rank-fourteen sampled signals. A rank-one signal returns a
frequency-by-time rank-two array, while higher ranks insert the frequency
axis at the selected signal axis and append a time-slice axis.
Unlike a static windowed FFT, the time axis also contains slices with centers
outside the data interval; `p_min` and `p_max(n)` give the zero-based slice
indices matching SciPy's convention. Fortran array indexing remains one-based.

```fortran
use scifort_signal, only : dp, short_time_fft, short_time_fft_result, &
    short_time_spectrogram_result, short_fft_success

type(short_time_fft) :: sft
type(short_time_fft_result) :: spectrum
type(short_time_spectrogram_result) :: power
real(dp) :: x(256)
real(dp), allocatable :: reconstruction(:)
integer :: status

! Supply the data in x before these calls.
call sft%initialize_window('hann',64,16,100.0_dp,mfft=128,scale_to='psd')
if (sft%status /= short_fft_success) error stop 'invalid STFT plan'
call sft%stft(x,spectrum,padding='zeros')
if (spectrum%status /= short_fft_success) error stop 'invalid data'
call sft%istft(spectrum%values,reconstruction,status,k1=size(x))
if (status /= short_fft_success) error stop 'inverse STFT failed'
call sft%spectrogram(x,power)
! spectrum%frequencies, spectrum%times, spectrum%values
! power%power equals abs(spectrum%values)**2
```

The type supports the following type-bound procedures:

- `initialize(real_or_complex_window,hop,fs,...)` and
  `initialize_window(name,nperseg,hop,fs,...)`, with optional FFT length
  `mfft`, `'onesided'`, `'onesided2X'`, `'twosided'`, or `'centered'` FFT
  layout, `'none'`, `'magnitude'`, or `'psd'` window scaling, phase shift,
  and an explicit dual window. `onesided2X` requires magnitude/PSD scaling.
- `initialize_from_dual(dual_window,hop,fs,...)` constructs the canonical
  analysis window for a supplied real or complex synthesis window.
  `initialize_equal_dual(desired_window,hop,fs,...)` builds a tight window
  whose unscaled analysis and synthesis windows are equal. The latter also
  supports `scale_to='unitary'`, where the dual equals `mfft*win`.
- `stft(x,spectrum,p0,p1,k_offset,padding,axis)` with optional zero-based
  time-slice bounds, sample offset, 1-based `axis` for rank-2--14 inputs, and `'zeros'`, `'edge'`, `'even'`, or `'odd'`
  boundary padding; `stft_detrend(x,spectrum,detrend,...)` supports `'none'`,
  `'constant'`, and `'linear'` detrending before window multiplication.
  For rank-1/2/3 real or complex inputs, `stft_detrend` also accepts an
  in-place user procedure as `detrend=`. The callback receives one
  unwindowed, boundary-extended segment before FFT window multiplication.
  It is called independently for every time slice and every batch channel.
  The input array is not changed. If the callback produces a NaN or infinity,
  the result has status `short_fft_invalid_input` and no allocated coefficient
  array. Callback and string forms of `detrend=` share the same type-bound
  generic name; they differ by the type of the argument.
- `istft(spectrum,samples,status,k0,k1,p0,f_axis,t_axis)` reconstructs
  rank-1 through rank-14 sample arrays. For higher-rank spectra, `f_axis` and `t_axis`
  identify frequency and time in the input coefficients; by default the
  last two axes are frequency and time. The time axis is removed on output
  and the frequency axis becomes the sample axis. `k0` defaults to zero;
  specify `k1` as the original signal length along the transform axis.
  One-sided spectra return real samples; two-sided/centered spectra return
  complex samples, even for real input.
- `spectrogram(x,result,...)` returns squared magnitudes and the coordinate
  axes, while `cross_spectrogram(x,y,cross,status,...)` returns `Sx*conjg(Sy)`.

### Custom ShortTimeFFT detrending

For custom nonlinear or polynomial detrending, supply an in-place Fortran
subroutine that accepts a rank-one `real(dp)` or `complex(dp)` slice. A
callback for real inputs **must** accept real samples; complex inputs require
a complex callback. To remove the per-slice mean without applying a window:

```fortran
use scifort_signal, only : dp, short_time_fft, short_time_fft_result

type(short_time_fft) :: sft
type(short_time_fft_result) :: spec
real(dp) :: x(128)

! Initialize sft and fill x before calling.
call sft%stft_detrend(x, spec, detrend=remove_mean, padding='even')

contains
    subroutine remove_mean(samples)
        real(dp), intent(inout) :: samples(:) !! one real STFT time slice
        samples = samples - sum(samples)/real(size(samples), dp)
    end subroutine remove_mean
```

This follows `scipy.signal.ShortTimeFFT.stft_detrend` callable semantics for
rank-1 inputs, with rank-2 and rank-3 Fortran extensions using an optional
1-based/negative `axis`. The existing string-based `detrend` remains available
through rank 14. The current callable interface is limited to rank 1--3;
that does not limit ordinary `stft`, `istft`, or named detrending.

Run the focused SciPy verification with
`python tools/run_signal_stft_callable_detrend.py` (GNU Fortran and
SciPy/NumPy required; FPM is not required).

The routines use SciFort's independent radix-two FFT with a direct DFT fallback
for other lengths; power-of-two `mfft` is faster. The higher-rank `stft`,
`stft_detrend`, `spectrogram`, and `cross_spectrogram` methods accept real or
complex rank-2 through rank-14 inputs. Rank-one `cross_spectrogram` additionally
accepts mixed real/complex inputs on a two-sided or centered FFT plan; see
[Mixed cross-spectral inputs](signal_mixed_cross.md). `axis` defaults to the last input axis;
positive axis indices use Fortran's 1-based convention, while negative axis
indices follow Python's convention (`-1` is the last dimension).
The rank-3 through rank-15 inverse accepts arbitrary distinct `f_axis` and `t_axis`,
including negative indices. `short_time_fft_result_3d` through `short_time_fft_result_15d` contain
rank-3 through rank-15 complex coefficients, and the corresponding
`short_time_spectrogram_result_3d` through `_15d` contain power. Since the
STFT adds one axis and the Fortran maximum rank is 15, a conventional
rank-15 input would need an unavailable rank-16 output. The input limit is
therefore **rank 14**; rank-15 input may be handled externally by flattening
batch dimensions and reshaping the transform result. Higher-rank interfaces
are generated by `scripts/gen_short_time_fft_high_rank.py` from the shared
rank-one analysis and synthesis kernels. The generator is retained for
reproducibility, rather than maintaining nine independent numerical kernels.
Results use error status codes, not Python exceptions.

For a batch of three real channels sampled at 256 time points:

```fortran
use scifort_signal, only : dp, short_time_fft, short_time_fft_result_3d
real(dp) :: channels(3,256)
real(dp), allocatable :: recovered(:,:)
type(short_time_fft) :: sft
type(short_time_fft_result_3d) :: coeff
integer :: ierr

call sft%initialize_window('hann',32,8,100.0_dp,mfft=64)
call sft%stft(channels,coeff,axis=2)
! coeff%values has shape [3,33,sft%p_num(256)]
call sft%istft(coeff%values,recovered,ierr,k1=256,f_axis=2)
```

### ShortTimeFFT coordinate and slice-range utilities

The class now includes the following query methods (matching SciPy's sampling
and slice-index conventions, but returning Fortran status codes instead of
raising Python exceptions):

- `p_range(n,first,last,status,p0,p1)` validates a slice interval; bounds are
  zero-based and the last slice is exclusive. `p_num(n)` counts time slices.
- `k_min()` and `k_max(n)` return the covered signal-index bounds;
  `lower_border_end()` and `upper_border_begin(n)` each return a two-integer
  vector `[sample_index, slice_index]` marking the padding transitions.
- `nearest_k_p(k,left)` finds the previous or next hop-aligned sample,
  including correct flooring for negative indices.
- `frequency_count()`, `frequency_points()`, `delta_f()` and `delta_t()`
  report frequency bins and frequency/time spacings.
- `time_points(n,status,p0,p1,k_offset)` returns slice centers in seconds.
  Time metadata from `stft(...,k_offset=...)` now incorporates that offset.
- `extent(n,bounds,status,axes_seq,center_bins)` returns plot bounds in
  time-frequency (`'tf'`) or frequency-time (`'ft'`) order. As in SciPy,
  unshifted `'twosided'` spectra have no single monotone frequency extent.
- `invertible()` tests whether a valid synthesis window is available.

For example:

```fortran
integer :: first, last, ierr
real(dp), allocatable :: f(:), t(:)
real(dp) :: bounds(4)

call sft%p_range(256,first,last,ierr)
f = sft%frequency_points()
t = sft%time_points(256,ierr,k_offset=5)
call sft%extent(256,bounds,ierr,axes_seq='tf',center_bins=.true.)
```

Even-length two-sided transforms now correctly label the Nyquist bin at
`-fs/2`, not `+fs/2`, and STFT result timestamps include `k_offset`.

Test programs are `test_signal_short_time_fft_high_rank.f90`,
`test_signal_short_time_fft_high_rank_public.f90` and
`test_signal_short_time_fft_helpers.f90`. The optional SciPy tests are
`compare_signal_short_time_fft_high_rank.py` (with companion Fortran reference
producer) and `compare_signal_short_time_fft_helpers.py`. The high-rank
comparisons check 216 configurations and 72,920 complex values (maximum
absolute difference 7.10e-15); coordinate comparisons check 216 cases
(maximum difference 3.55e-15). SciPy 1.17 has high-dimensional padding
limitations, so the high-rank reference applies SciPy's one-dimensional STFT
independently to each signal channel.

`scifort_signal` also exports `freqz_sos` as a SciPy-style alias of the
previously available `sosfreqz` routine.

Native tests are `test_signal_short_time_fft.f90` and
`test_signal_short_time_fft_public.f90`. The optional SciPy 1.17 comparison is
`test/compare_signal_short_time_fft.py` with the companion Fortran reference
producer `validation/reference_drivers/reference_signal_short_time_fft.f90`. It checks 312 diverse
rank-1 cases (44,058 complex coefficients) including data type, odd/even FFT
size, scaling, padding, phase, detrending, slice bounds, and sample offsets.
The rank-2/3 regression `test/compare_signal_short_time_fft_nd.py` with
`validation/reference_drivers/reference_signal_short_time_fft_nd.f90` compares another 92 cases
(62,091 coefficients). SciPy 1.17 fails some multidimensional zero-padding
cases internally; where necessary, the reference applies SciPy's rank-1 STFT
independently to each channel, preserving the documented axis layout.
No new library dependencies are introduced.

Rank-4 and rank-5 signals are supported by the same transform, detrending,
inverse, auto-spectrogram, and cross-spectrogram overloads. The new
`short_time_fft_result_5d` and `_6d` / `short_time_spectrogram_result_5d`
and `_6d` types are exported from `scifort_signal`. Frequency is inserted at
the input signal axis and time remains the final output axis. Inverse transforms
can use any two distinct frequency and time axes; removal of the time axis
preserves the relative order of other dimensions. A native regression is
`test/test_signal_short_time_fft_rank45.f90` and SciPy comparison is
`test/compare_signal_short_time_fft_rank45.py`, using
`validation/reference_drivers/reference_signal_short_time_fft_rank45.f90`. The reference suite covers
108 rank-4/5 configurations and 208,000 complex values; its largest absolute
error against SciPy 1.17 per-channel transforms was 8.15e-15. SciPy 1.17 has
an internal high-rank padding issue, so these references independently apply
SciPy's rank-one transform to each channel. The rank-4/5 Fortran implementation
also checks forward/inverse reconstruction and arbitrary inverse axis placement.

## Newer signal APIs: smoothing, normalization and convolution selection

The following additions are also public members of `scifort_signal`:

- `whittaker_henderson(signal, lamb, order, weights)` returns a
  `whittaker_result` with `x`, `lamb` and `status`. The optional fixed penalty
  `lamb` must be nonnegative; when absent, an automatic restricted-maximum-
  likelihood (REML) optimization selects it, following the profiled Gaussian
  criterion in SciPy 1.18. `order` defaults to 2 (supported
  1--12), and `weights` defaults to ones. Setting weights to zero marks
  missing values, including NaN observations. A positive-definite banded
  Cholesky solver solves the weighted least-squares difference-penalty
  problem in O(n * order**2) time and O(n * order) memory per objective
  evaluation. Automatic selection scans log-lambda over [-24,24] and
  refines a bracketed minimum with golden-section search, requiring dozens
  of banded factorizations. For nearly polynomial/degenerate data the optimum
  may lie outside this finite interval; this is an explicit numerical limit.
  With insufficient observed data the matrix can be singular
  (`whittaker_status_singular`).
- `normalize(b, a)` returns `normalized_tf_result` containing complex
  `b(:,:)`, complex `a(:)`, `trimmed` and `status`. All real/complex and mixed
  variants are supported, with a rank-1 or rank-2 numerator and rank-1
  denominator. Rows of `b(:,:)` are output channels; a rank-1 input still
  returns one row. Leading *exact* zero denominator coefficients and leading
  numerator columns with magnitude at most `1e-14` after normalization are
  removed, leaving at least one numerator column.
- `abcd_normalize(a, b, c, d)` returns `state_space_result` with
  `a(:,:)`, `b(:,:)`, `c(:,:)`, `d(:,:)` and `status`. The optional arguments
  are rank-2 intrinsic real, complex or integer Fortran arrays; a scalar D
  must be passed as `reshape([d],[1,1])`. Absent or zero-size matrices are
  supplied as zeros once state, input and output dimensions are inferred.
- `choose_conv_method(x, h, mode, measure, times, status)` returns a
  6-character `'direct'` or `'fft'` selector for real/complex rank-1 arrays, including mixed real-complex inputs.
  `measure=.true.` benchmarks both methods using `cpu_time` and optionally
  returns timings in `times(1:2)=[direct,fft]`. The default heuristic is tuned
  to SciFort's own radix-2 FFT and need not produce the same choice as
  SciPy's platform-dependent selector. `convolve` and `correlate` now accept
  `'direct'`, `'fft'`, or `'auto'` (default), with `'auto'` using the selector.
  Inputs with NaNs/Infs always select direct convolution under `'auto'`.

```fortran
use scifort_signal, only : dp, whittaker_henderson, whittaker_result, &
    normalize, normalized_tf_result, abcd_normalize, state_space_result, &
    convolve, choose_conv_method
real(dp) :: y(6), t(6)
type(whittaker_result) :: smooth
type(normalized_tf_result) :: tf
type(state_space_result) :: matrices
character(6) :: method

y = [1.0_dp, 2.0_dp, 5.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
smooth = whittaker_henderson(y, 10.0_dp, order=2)  ! fixed penalty
smooth = whittaker_henderson(y, order=2)             ! REML penalty
tf = normalize([0.0_dp, 2.0_dp, 4.0_dp], [2.0_dp, 1.0_dp])
matrices = abcd_normalize(b=reshape([1.0_dp, 2.0_dp], [2,1]), &
    c=reshape([3.0_dp, 4.0_dp], [1,2]))
method = choose_conv_method(y, y)
t = convolve(y, [1.0_dp, 1.0_dp], mode='same', method='fft')
```

The regression driver `validation/reference_drivers/test_signal_missing_extra_reference.f90` and
`test/compare_signal_missing_extra_scipy.py` compare numerical results with
SciPy 1.17 and a sparse linear-system reference for Whittaker-Henderson.
Automatic REML tuning is implemented and checked in
`test/compare_signal_whittaker_reml.py`, which compares 27 weighted/unweighted
cases (972 fitted values) against an independent dense-Cholesky likelihood
optimizer, following the SciPy 1.18 criterion. The rank-independent
`convolve(method='auto')` remains future work. The `short_time_fft` type
now supports input ranks 1-3 (see above).

## Newly covered signal algorithms (October 24, 2026)

The following seven SciPy-compatible routines are exported by the public
`scifort_signal` module. They use binary64 real/complex values and return an
optional integer `status` (`0` indicates success). Invalid inputs produce
empty allocated arrays. No new FPM dependencies are required.

**Short-time Fourier transforms:**
`closest_stft_dual_window(win, hop, dual, alpha, desired_dual, scaled, status)`
is a subroutine returning the dual window closest to `desired_dual` (ones by
default), optionally optimizing `alpha` (`scaled=.true.` by default). Fortran
has separate real and complex overloads; the 1-based hop length must be
between 1 and `size(win)`. The implementation groups window samples by
residue class modulo `hop`, using a diagonal overlap energy rather than dense
least squares. A noninvertible overlap returns `stft_dual_status_noninvertible`.

**Two-dimensional FIR design:**
`firwin_2d(hsize, window, fc, window2, circular, fs, pass_zero, scale, status,
window_parameter, window_parameter2, fc2, fs2, pass_zero2, scale2,
ellipse_axes, ellipse_angle)` returns a real two-dimensional array.
`hsize=[nrow,ncol]`, and `fc` is either a scalar or a vector of band edges.
`window` and optional `window2` are names recognized by `firwin`/`get_window`;
parameterized windows (for example Kaiser) accept the corresponding real
shape parameters. With `circular=.false.`, the result is a separable outer
product, optionally with highpass/bandpass/bandstop response. With
`circular=.true.`, the implementation interpolates an oversampled radial
profile. For non-square `hsize`, its output shape is `[hsize(2),hsize(1)]`,
matching SciPy 1.17's `meshgrid` convention; this differs from the dimensions
advertised in SciPy's documentation. Unlike SciPy 1.17's implementation,
SciFort actually forwards `pass_zero` to the underlying `firwin` filters.
Large designs have a bounded allocation guard.

**SciFort extension: independent-axis and elliptical FIR designs.** For
separable designs, `fc2` supplies a **rank-one vector** of second-axis
cutoff edges, even when there is only one cutoff (write `fc2=[0.35_dp]`).
`fs2`, `pass_zero2`, and `scale2` optionally select the second-axis sample
rate, passband type, and normalization independently. By default `fs2` and
`scale2` inherit their first-axis values. When `pass_zero` is explicitly
supplied, that is also the default for `pass_zero2`; otherwise each axis
infers its default passband from its own number of cutoffs, as `firwin` does.
The first-axis window and
cutoffs are unchanged, and the result remains the outer product of two
independently designed one-dimensional FIR filters. Unlike SciPy 1.17's
`firwin_2d`, this can design distinct lowpass, highpass, bandpass, or
bandstop responses along each axis. Each axis must independently meet the
usual `firwin` restrictions, including odd length for a passband reaching
Nyquist.

For radial designs (`circular=.true.`), `ellipse_axes=[a,b]` and optional
`ellipse_angle` (radians, positive counterclockwise in normalized row/column
coordinates) define the interpolated elliptical radius. The normalized
coordinates run from -1 to +1 along each requested dimension; `a` and `b`
are strictly positive and default to 1. This extends SciPy's circular
construction to rotated ellipses, but is a *coefficient-domain radial
interpolation*, not an exact elliptical frequency-response specification.
The existing transposed output convention for non-square radial designs is
preserved. Independent-axis options cannot be combined with radial designs;
ellipse options cannot be combined with separable designs. Conflicting,
nonfinite, or otherwise invalid arguments return `design_status_invalid_input`
and an allocated empty `(0,0)` array.

```fortran
use scifort_signal, only : dp, firwin_2d
real(dp), allocatable :: fir(:,:)
integer :: status

! Different cutoff frequencies and passband types on each axis.
fir = firwin_2d([9,13], 'hamming', 0.18_dp, &
    fc2=[0.16_dp,0.48_dp], pass_zero2='bandpass', status=status)

! Rotated elliptical radial profile (output is [13,9] here).
fir = firwin_2d([9,13], 'hamming', 0.2_dp, circular=.true., &
    ellipse_axes=[0.65_dp,1.4_dp], ellipse_angle=0.4_dp, status=status)
```

The independent SciPy numerical harness is
`tools/signal_reference/compare_signal_fir2d_anisotropic_scipy.py`.
It compares each separable design with an outer product of two SciPy `firwin`
results and each elliptical result with independently rotated NumPy grids and
`numpy.interp` applied to a SciPy `firwin` radial profile. The original
`firwin_2d` comparison is retained separately to verify defaults.

**Recursive smoothing with mirror-symmetric boundaries:**
`symiirorder1(signal, c0, z1, precision, status)` and
`symiirorder2(signal, r, omega, precision, status)` use the existing spline
recurrences. They accept real and complex rank-1 arrays, as well as rank-2
batches filtered along the second dimension. First-order `z1` is a real pole,
with transfer function `c0/((1-z1/z)*(1-z1*z))`. The second-order filter has
stable conjugate poles with radius `r` and angle `omega`. Poles must satisfy
`abs(z1)<1` or `abs(r)<1`; the second-order recurrence requires nonzero
`sin(omega)` except when `r=0`. The optional `precision` controls boundary
series truncation; the default is `1e-11`. If a finite record is too short
for the recursive tail to decay, the result has status
`symiir_status_no_convergence`. Complex input to `symiirorder2` is an extension
to SciPy 1.17, which only accepts real input for that routine.

**Analog grids, repeated poles, and chirp-Z points:**
`findfreqs(num, den, n, kind, status)` constructs the logarithmic analog
frequency grid for transfer-function coefficients (`kind='ba'`, default)
or directly supplied roots (`kind='zp'`). Both real and complex inputs are
supported. `unique_roots(p, roots, multiplicities, tol, rtype, status)` is a
subroutine grouping approximate roots greedily, in input order, with absolute
Euclidean tolerance (`tol=0.001` by default). `rtype` selects `'min'`,
`'max'`, or `'avg'` (including SciPy aliases). Its straightforward quadratic
grouping is intended for modest filter orders. `czt_points(m, w, a, status)`
returns the complex locations `a*w**(-k)`, with defaults forming the unit
circle; it complements the existing `czt` and `zoom_fft` functions.

```fortran
use scifort_signal, only : dp, closest_stft_dual_window, firwin_2d, &
    findfreqs, unique_roots, czt_points
real(dp), allocatable :: synthesis(:), mask(:,:), freq(:), roots(:)
complex(dp), allocatable :: spiral(:)
integer, allocatable :: mult(:)
real(dp) :: alpha
integer :: status

call closest_stft_dual_window([1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], &
    2, synthesis, alpha, status=status)
mask = firwin_2d([9,9], 'hamming', 0.15_dp, circular=.true., status=status)
freq = findfreqs([1.0_dp,0.0_dp],[1.0_dp,8.0_dp,25.0_dp],20,status=status)
call unique_roots([1.0_dp,1.0005_dp,3.0_dp],roots,mult,status=status)
spiral = czt_points(16,status=status)
```

Native regressions live in `test/test_signal_missing_algorithms.f90` and
`test/test_signal_missing_algorithms_public.f90`. Run
`python tools/signal_reference/compare_signal_missing_algorithms_scipy.py`
for 499 SciPy 1.17 reference comparisons covering 11,610 numeric results
(maximum observed absolute difference `2.49e-14`), or add `--optimized`
for an optimized test build.

## Auditory and periodic comb filter design

`gammatone(freq, ftype, order, numtaps, fs)` produces the coefficients of
an auditory gammatone filter. `ftype='fir'` (default order 4) uses the
Glasberg/Moore equivalent rectangular bandwidth fit and Slaney's impulse
response; `ftype='iir'` produces a fourth-order gammatone approximation as
an **eighth-order** real-coefficient transfer function. The IIR form ignores
`order` and `numtaps`, as in SciPy. Frequencies must lie strictly between
zero and `fs/2` (default `fs=2`). FIR orders range from 1 to 24 and FIR
lengths must be positive. The FIR default length is
`max(int(0.015*fs),15)`.

`iircomb(w0, qfactor, ftype, fs, pass_zero)` designs a notch or peak comb.
`ftype` defaults to `'notch'`, `fs` defaults to 2, and `pass_zero` defaults
to `.false.`. The spacing `w0` must exactly divide the sample rate within
SciPy's `1e-14` relative tolerance, and the Q factor must specify a stable
positive tangent-bandwidth response. To prevent inadvertent enormous
allocations, the returned comb order is capped at one million samples.
Both filter types are case-insensitive ASCII strings.

Each returns a `butter_result`: `%b(:)` and `%a(:)` are real coefficients,
`%status` is `design_status_success` or `design_status_invalid_input`, and
`%sos` is unallocated. Filter coefficients use the standard descending
powers of z^-1 employed by `lfilter` and `freqz`. In narrowband IIR
applications, direct-form eighth-order filtering may be poorly conditioned:
consider conversion to second-order sections where appropriate.

```fortran
use scifort_signal, only : dp, butter_result, gammatone, iircomb, &
    design_status_success
type(butter_result) :: auditory, harmonic

auditory = gammatone(440.0_dp, 'fir', numtaps=128, fs=16000.0_dp)
harmonic = iircomb(50.0_dp, 30.0_dp, ftype='notch', fs=1000.0_dp)
if (auditory%status /= design_status_success) error stop 'gammatone'
if (harmonic%status /= design_status_success) error stop 'comb'
```

**Reference:** SciPy 1.17.0 `scipy.signal.gammatone` and `iircomb`
(whose formulas reference M. Slaney, 1993, and S. Orfanidis, 1996).
See `CODE_PROVENANCE.md` for attribution. Reproduce 1,350 SciPy
coefficient comparisons with
`python tools/signal_reference/compare_signal_auditory_comb.py`.

## Wavelet peak detection

`find_peaks_cwt(x, widths, max_distances, gap_thresh, min_length, min_snr,
noise_perc, window_size, status, wavelet)` implements SciPy-style continuous
wavelet ridge detection for a finite, real-valued one-dimensional signal.
The optional procedure `wavelet=` accepts two `real(dp)` scalar arguments,
`points` and `width`, and returns an allocatable one-dimensional `real(dp)`
array. The first argument is `min(10*width, size(x))` and **may be fractional**;
the kernel must have `ceiling(points)` finite samples. Kernels should normally
be symmetric and appropriately normalized for meaningful ridge detection.
This interface currently supports **real-valued wavelet callbacks only**.

The default is the normalized Ricker (Mexican-hat) wavelet. The public helper
`cwt_ricker(x,widths,status)` is retained for compatibility, and the new
`cwt_wavelet(x,widths,wavelet,status)` returns the transform for any valid real
callback. Both transform results have shape `(size(widths),size(x))`, with
wavelet scales in rows. `cwt_wavelet_callback` is the published abstract
interface used for callback type checking. `ricker_wavelet(points,width)`
returns the default wavelet samples. As in SciPy, the kernel is reversed for
convolution; this matters for asymmetric user kernels.

Detected indices are **one-based** and sorted ascending. `max_distances`, when
supplied, must cover every width. Invalid widths, inputs, callback lengths,
or nonfinite callback outputs return empty arrays and status
`cwt_peak_invalid_input`. A failure at any width discards the entire CWT.
On nearly constant signals, strict extrema comparisons can be sensitive to
summation roundoff and can differ by one sample from SciPy.

```fortran
real(dp) :: signal(100), widths(8)
real(dp), allocatable :: coefficients(:,:)
integer, allocatable :: peaks(:)
integer :: status
widths = [1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp,7.0_dp,8.0_dp]
peaks = find_peaks_cwt(signal,widths,min_snr=1.0_dp,status=status)
coefficients = cwt_wavelet(signal,widths,scaled_ricker,status=status)
peaks = find_peaks_cwt(signal,widths,wavelet=scaled_ricker,status=status)

contains
    function scaled_ricker(points,width) result(w)
        real(dp), intent(in) :: points !! possibly fractional number of samples
        real(dp), intent(in) :: width !! positive wavelet scale
        real(dp), allocatable :: w(:) !! sampled custom kernel
        w = 0.5_dp*ricker_wavelet(points,width)
    end function scaled_ricker
```

Reproduce SciPy 1.17 comparison of default Ricker behavior with
`python tools/signal_reference/verify_cwt_place_poles_scipy.py`, and custom
wavelet behavior with
`python tools/signal_reference/compare_signal_cwt_custom_scipy.py`.

## Complex continuous-wavelet transforms

`cwt_complex_wavelet(x,widths,wavelet,status)` extends the existing real
continuous-wavelet calculation to **complex-valued wavelets**, with real or
complex rank-one input `x`. `cwt_complex_wavelet_callback` describes the
`wavelet(points,width)` Fortran procedure, which must return an allocated,
finite rank-one `complex(dp)` wavelet of length `ceiling(points)`.
The unrounded `points=min(10*width,size(x))` value is passed to the callback,
even for fractional scales. Wavelet rows are ordered by the caller's widths;
columns correspond to input sample positions.

`cwt_morlet2(x,widths,omega0,status)` provides the normalized Morlet wavelet
with default angular center frequency `omega0=5.0_dp`; it accepts real or
complex input, and optionally adjusts `omega0` for different time-frequency
localization. `morlet2_wavelet(points,width)` provides the default kernel
as a reusable complex callback. Both output interfaces return allocatable
`complex(dp)` arrays with shape `(size(widths),size(x))` on success, or an
empty `(0,0)` array and `cwt_peak_invalid_input` on errors.

The sampled Morlet family is

\[
\psi_s(k) = \pi^{-1/4}s^{-1/2}\exp(-t_k^2/2)\exp(i\omega_0t_k),
\qquad t_k = \frac{k-(L-1)/2}{s}, \quad
L=\min(10s,N).
\]

The CWT uses **conjugate and time-reversed** complex wavelets in the
`same` zero-padded convolution, not the ordinary unconjugated convolution.
Invalid scales, nonfinite input samples, nonfinite Morlet frequencies,
malformed/nonfinite callbacks, and nonfinite computed outputs are rejected.
The existing real `cwt_wavelet`, `find_peaks_cwt`, and Ricker APIs retain
signed-coefficient ridge maxima. Complex coefficients are handled by the
**distinct magnitude-based** peak finder described below.

```fortran
use scifort_signal, only : dp, cwt_morlet2, cwt_complex_wavelet, morlet2_wavelet
real(dp) :: x(128), widths(3)
complex(dp), allocatable :: transform(:,:), default_transform(:,:)
widths = [1.5_dp, 3.0_dp, 6.0_dp]
transform = cwt_morlet2(x, widths, omega0=4.0_dp)
default_transform = cwt_complex_wavelet(x, widths, morlet2_wavelet)
```

SciPy removed its public `signal.cwt` and `signal.morlet2` functions before
1.17. The maintained CWT semantics are checked with the public SciPy 1.17
`signal.convolve(data, conj(wavelet[::-1]), mode='same')` reference,
using independently generated Morlet and asymmetric complex wavelets.
`python tools/run_signal_cwt_complex.py` compiles a checked and optimized
standalone suite, including 240 SciPy comparisons and native invalid-input
regressions, without requiring FPM.

## Magnitude-based complex wavelet peak detection

`find_peaks_cwt_complex(x,widths,max_distances,gap_thresh,min_length,min_snr,`
`noise_perc,window_size,status,wavelet,omega0)` accepts rank-one real or
complex input and returns a sorted, one-based vector of peak positions.
By default it computes `cwt_morlet2(x,widths)` with `omega0=5.0_dp`.
Pass `omega0=` for a different Morlet angular center frequency **or**
`wavelet=` for a caller-supplied `cwt_complex_wavelet_callback`, but not both.
The callback receives `(points,width)`, where `points` may be fractional,
and must return an allocatable `complex(dp)` vector of `ceiling(points)`
finite coefficients.

The ridge criterion is defined explicitly as **local maxima of
`abs(CWT(scale,position))`**, not maxima of the real part. Ridge lines are
connected from coarser to finer scales using `max_distances` (default
`widths/4`) and `gap_thresh` (default `widths(1)`). Filtering requires
`min_length` connected maxima (default `ceiling(nwidths/4)`) and a
magnitude signal-to-noise ratio (default threshold `min_snr=1`). Noise at
position `j` is the `noise_perc` percentile (default 10) of
`abs(CWT(1,:))` in a local `window_size` neighborhood (default
`ceiling(nsamples/20)`). Scales should normally be ordered from narrow to
wide for meaningful ridge linking. A scale row with no strict maxima
cannot by itself contribute a peak.

`find_peaks_cwt_magnitude(coefficients,widths,max_distances,gap_thresh,`
`min_length,min_snr,noise_perc,window_size,status)` performs **the same
magnitude ridge selection** on an already computed complex CWT matrix.
Here `coefficients` has shape `(size(widths),nsamples)`. This avoids an
extra CWT when testing multiple detection thresholds or using a different
complex wavelet transform. Invalid or nonfinite input, inconsistent dimensions,
invalid ridge parameters, conflicting `omega0`/`wavelet`, or malformed
callback outputs result in an empty peak array with
`cwt_peak_invalid_input`. Successful processing (including no peaks)
sets `cwt_peak_success`.

```fortran
use scifort_signal, only : dp, cwt_morlet2, find_peaks_cwt_complex, &
    find_peaks_cwt_magnitude
real(dp) :: widths(5), x(256)
complex(dp), allocatable :: cwt(:,:)
integer, allocatable :: peaks(:), matched(:)
integer :: status
widths = [2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp]
peaks = find_peaks_cwt_complex(x,widths,omega0=4.0_dp,status=status)
cwt = cwt_morlet2(x,widths,omega0=4.0_dp)
matched = find_peaks_cwt_magnitude(cwt,widths,status=status)
```

This is a SciFort extension, **not** a claim that SciPy exposes a public
complex `find_peaks_cwt` equivalent. The independent test constructs the
complex transform with SciPy's public `signal.convolve` and applies SciPy
1.17's private real-valued ridge functions to its magnitude. Those private
functions are **test-only** and can change in future SciPy releases.
Run `python tools/run_signal_cwt_complex_peaks.py` for checked and optimized
native tests, 280 independent magnitude-ridge cases, 300 legacy real-Ricker
comparisons, and the preceding 240 complex-CWT coefficient comparisons.

## Real state-feedback pole placement (SISO and MIMO)

`place_poles(a,b,poles)` accepts real state and input matrices `a(n,n)`,
`b(n,m)` and either real or complex desired eigenvalues `poles(n)`. It returns
`pole_placement_result` containing `gain_matrix(m,n)`, `requested_poles(:)`,
`status`, `nb_iter`, `rtol`, and `det_transfer`. The intended closed-loop system is `a - matmul(b,gain_matrix)`.
Complex target poles must occur in conjugate pairs. The following cases are handled:

- `m=1`, with a nonsingular controllability matrix and distinct target poles;
  this uses the analytic Ackermann formula.
- `m>=n`, when `b` has full row rank; a real minimum-norm input gain
  realizes a block-diagonal real target with the requested poles.
- Partially actuated `1<m<n` systems, and `m>=n` with rank-deficient `b`,
  using a pivoted-nullspace eigenvector assignment. It selects feasible
  eigenvectors from `(a-lambda*I)*x=b*u` and solves `K*X=U` with a real
  transfer matrix. Complex-conjugate pairs and repeated poles (multiplicity
  no greater than `rank(b)`) are supported **when independent eigenvectors
  can be found**. Real KNV0 rank-one projections and Yang--Tits real-pair
  singular-vector/complex-pair Hermitian updates improve the transfer matrix.

Optional parameters `method='YT'` (default) or `method='KNV0'`,
`rtol=1.0e-3_dp` (default) and `maxiter=30` (default) control determinant
optimization. `KNV0` accepts **real** poles only. Negative `rtol` forces
all `maxiter` iterations, as in SciPy. The result's `nb_iter` counts
iterations, `rtol` gives their achieved relative determinant change and
`det_transfer` is the absolute determinant of the chosen **real** transfer
matrix. The diagnostics are zero for SISO/fully actuated direct solutions.

The YT updates use SciPy's recommended pair scheduling but the SciFort
nullspace basis, deterministic initialization, singular-vector tie handling,
and real transfer basis are not identical to SciPy's. Consequently the gain,
iteration count and raw transfer determinant need not agree with SciPy even
when poles are correctly assigned. In particular, `det_transfer` is not
numerically identical to SciPy's complex-basis `abs(det(X))` when complex
conjugate pairs are present.

**Not implemented:** exact parity with all SciPy YT/KNV0 numerical tie cases
or generalized Jordan chains for defective closed-loop assignments. The nullspace method may fail for highly conditioned
systems or repeated complex poles that require a defective eigenstructure;
this is not a promise of success for every mathematically controllable pair.
The Fortran gain need not equal SciPy's (MIMO solutions are nonunique).

Invalid dimensions, non-conjugate roots or excessive pole multiplicity
return `pole_status_invalid_input`. Uncontrollable input matrices return
`pole_status_uncontrollable` when detected; failure to find an independent
transfer matrix or unstable floating-point results return
`pole_status_numerical_failure`. The result contains *requested*, not
independently computed, poles. Users should check the eigenvalues and
conditioning of `a-matmul(b,gain_matrix)` before relying on the design.

```fortran
real(dp) :: a(2,2), b(2,1)
type(pole_placement_result) :: placed
a(1,:) = [0.0_dp,1.0_dp]
a(2,:) = [-2.0_dp,-3.0_dp]
b(:,1) = [0.0_dp,1.0_dp]
placed = place_poles(a,b,[-4.0_dp,-5.0_dp])
! placed%gain_matrix(1,:) is approximately [18.0_dp, 6.0_dp].
! For a multi-input system you can select iterative optimization:
! placed = place_poles(a,b,poles,method='YT',rtol=1.0e-4_dp,maxiter=50)
```

References: SciPy 1.17 `signal.find_peaks_cwt` and `signal.place_poles`.
Reproduce the additional 470 MIMO comparisons with
`python tools/signal_reference/verify_mimo_place_poles_scipy.py`.
Reproduce the 666 YT/KNV0 method-and-iteration comparisons with
`python tools/signal_reference/verify_iterative_pole_placement_scipy.py`.
The source provenance and applicable SciPy BSD-3-Clause attribution are in
`CODE_PROVENANCE.md` and `THIRD_PARTY_LICENSES.md`. No new FPM dependencies.

## Order statistics and separable spatial FIR filtering

`order_filter(x, domain, rank, status)` evaluates the sorted order statistic
selected by a **zero-based** `rank` over the nonzero/true elements of the
centered `domain` footprint. Both inputs must have the same rank (1, 2, or 3).
`x` may be `real(dp)` or default `integer`; `domain` may be logical, default
integer, or `real(dp)` with **positive odd extents**. The output has the same
shape and element type as `x`; samples outside the input are zero. `rank=0`
selects the minimum; `rank=count(domain /= 0)-1` selects the maximum.
Noncontiguous inputs are supported. The implementation uses a three-way
quickselect rather than sorting all neighbor values, and does not specify
ordering for NaNs.

`sepfir2d(x, hrow, hcol, status)` applies an odd-length FIR filter `hrow`
along the **second** array dimension, followed by `hcol` along the **first**
dimension (matching SciPy's row/column naming). Coefficients are applied as
convolution, not correlation; the boundary extension repeats edge samples
when reflected, and also handles tap lengths longer than image dimensions.
Real and complex `real(dp)`/`complex(dp)` images and coefficients may be mixed;
default-integer images are promoted to real or complex results. Where SciPy
rejects complex coefficients with a real or integer image, SciFort deliberately
supports that combination by promotion.

```fortran
use scifort_signal, only : dp, order_filter, sepfir2d
real(dp) :: image(80,120), hrow(3), hcol(5)
real(dp), allocatable :: ranked(:,:), filtered(:,:)
integer :: status

ranked = order_filter(image, reshape([1,0,1, 0,1,0, 1,0,1],[3,3]), &
    rank=2, status=status)
filtered = sepfir2d(image, hrow, hcol, status=status)
```

`spatial_status_success=0` and `spatial_status_invalid_input=1` identify
successful operations and invalid window lengths, masks, or ranks; invalid
calls return a zero-filled result of the input's shape. There are no new
runtime dependencies. A standalone reference script,
`tools/signal_reference/compare_signal_spatial_scipy.py`, checks SciPy and
independent `ndimage.generic_filter` baselines. SciPy 1.17's optimized
`sparse-footprint` median path returned different values on three of the
504 randomized masked-order cases; SciFort follows the documented footprint
semantics (as `generic_filter` does) instead. For safety, `sepfir2d` SciPy
comparisons restrict filter lengths to the image dimensions because some
larger-filter/small-image combinations caused a crash in the installed SciPy
version. SciFort handles such filters through repeated reflection; an extra
Fortran regression tests this extension.
The standalone GNU Fortran 14.2 signal suites pass 72/72 tests with runtime
checks and 72/72 optimized; project-wide FPM testing was unavailable.

To repeat the randomized comparison after compiling the Fortran driver:

```bash
gfortran -O3 -c src/kinds.f90 src/signal_spatial.f90
gfortran -O3 kinds.o signal_spatial.o \
    tools/signal_reference/spatial_reference_driver.f90 -o spatial_driver
python tools/signal_reference/compare_signal_spatial_scipy.py ./spatial_driver
```

## Median and Wiener noise reduction

The generic `medfilt`, `medfilt2d`, and `wiener` procedures accept **rank-1,
rank-2, and rank-3** arrays. `medfilt2d` is restricted to rank 2; it shares the
same general-purpose filtering engine as `medfilt`, without SciPy's optimized
2D C kernel. Median filters accept `real(dp)` and default `integer` inputs,
returning the same type; Wiener filters accept both types and return `real(dp)`.
`wiener` also accepts `complex(dp)` arrays (and returns complex arrays), with
real or complex explicitly supplied noise. Other integer kinds, complex
median filtering, and ranks 4-15 are not yet implemented.

```fortran
use scifort_signal, only : dp, medfilt, medfilt2d, wiener
real(dp) :: image(40,60)
real(dp), allocatable :: clean(:,:), smooth(:,:)
integer :: status

clean = medfilt2d(image, kernel_shape=[3,5], status=status)
smooth = wiener(image, kernel_size=5, noise=0.02_dp, status=status)
```

Each function supports optional `kernel_size=` (one scalar applied to every
axis) or `kernel_shape=` (one integer size per dimension); omit both to use
3 on every axis. Supplying both is invalid. Median sizes must be positive
and odd, but can exceed the input extents; `wiener` accepts positive odd or
even sizes. Outside the array, the input is treated as zero. For even Wiener
windows, the center matches SciPy's `correlate(..., mode='same')` convention.

`wiener` estimates the noise power as the arithmetic mean of local variances
unless the optional `noise=` is provided. Its variances use the population
second moment, `mean(x**2) - mean(x)**2`. For complex signals this is a
*complex pseudo-variance*, not `mean(abs(x)**2)-abs(mean(x))**2`, matching
SciPy. The complex noise estimate and the comparison of complex values
follow NumPy conventions; see [Complex Wiener filtering](signal_wiener_complex.md)
for the numerical details and examples. When the local variance is below
the noise power, the result is the local mean. A zero-variance region with
zero noise power yields a quiet NaN, matching SciPy's indeterminate gain.

Optional `status=` is `denoise_status_success` (0) or
`denoise_status_invalid_input` (1). Invalid parameter combinations produce
a shape-preserving, zero-filled result for legacy real/integer inputs, or
a complex array of quiet NaNs for complex Wiener inputs. Input/output may be noncontiguous
since the routines copy into local column-major working arrays. Median
selection uses a three-way quickselect; the implementation is intended for
finite inputs, and does not specify the ordering of NaNs.

Standalone Fortran regressions are in `test/test_signal_denoise.f90` and
`test/test_signal_denoise_public.f90`; 1,050 randomized comparisons with
SciPy 1.17 are reproduced by
`tools/signal_reference/compare_signal_denoise_scipy.py`.
The standalone 70-test signal suite passes with GNU Fortran 14.2.0 in both
runtime-checked and `-O3` configurations; the new direct regression and
1,050 SciPy comparisons also pass with floating-point traps enabled. FPM was
unavailable in the validation environment. To repeat the reference comparison
from the project root:

```bash
gfortran -O3 -fcheck=all -c src/kinds.f90 src/signal_denoise.f90
gfortran -O3 -fcheck=all kinds.o signal_denoise.o \
    tools/signal_reference/denoise_reference_driver.f90 -o denoise_driver
python tools/signal_reference/compare_signal_denoise_scipy.py ./denoise_driver
```


## Multidimensional FFT convolution and axis-selectable filtering

`fftconvolve` now accepts **rank-1, rank-2, and rank-3** arrays of `real(dp)`
or `complex(dp)` values. All three ranks accept mixed real/complex
operands. In addition to `mode='full'` (the default), `same` and
`valid` are supported. The optional `axes=[...]` argument selects which axes
are convolved; if omitted, all axes are convolved. **Axes are one-based**, as
in all other SciFort indexing APIs: `[1,2]` means the first two dimensions.
Nonconvolved dimensions must match or be broadcastable from an extent of one.
SciPy's `same` centering across broadcast dimensions is also reproduced.

```fortran
real(dp) :: frames(64,64,8), kernel(5,5,1)
real(dp), allocatable :: blurred(:,:,:)
integer :: status

blurred = fftconvolve(frames, kernel, mode='same', &
    axes=[1,2], status=status)
! shape(blurred) = [64,64,8]; each frame is filtered independently
```

For ranks **1 through 15**, use the Fortran 2018 caller-output subroutine
`fftconvolve_nd(in1, in2, y, mode, axes, status)`. The caller allocates `y`
with the requested convolution shape. Its real/complex type is real only when
both operands are real; mixed types produce complex output. `axes` uses the
same one-based convention as `fftconvolve` and defaults to all axes.

```fortran
real(dp) :: frames(32,24,2,3,10), kernel(3,3,1,1,1)
real(dp) :: blurred(32,24,2,3,10)
integer :: status

call fftconvolve_nd(frames, kernel, blurred, mode='same', &
    axes=[1,2], status=status)
```

Each selected axis has output length `n+m-1` in `full` mode, `n` in `same`
mode, and `abs(n-m)+1` in `valid` mode. Unselected axes broadcast singleton
extents and return their maximum extent in `full`/`valid` mode; `same` mode
returns the first input's extent (with center cropping). Inputs must have the
same rank, and the caller's output must match the expected shape exactly.
On invalid rank, shape, axes, empty input, nonfinite input, or unsupported
mode, `status` is `convolution_status_invalid_input` and provided output
cells are filled with quiet NaNs. Input/output actual arguments must not
alias each other. Strided sections are supported through Fortran copy-in/out.
The internal column-major storage remapping requires `dp` to interoperate
with C double on the compiler, as for the other SciFort assumed-rank APIs.

The multidimensional implementation pads each transformed dimension to a
power of two, applies separable complex FFTs, multiplies the spectra, and
performs inverse transforms. It uses SciFort's existing FFT and **adds no
external dependency**. The one-dimensional array-returning `oaconvolve` API now also accepts mixed
real/complex operands (both orders), without changing existing calls.
For ranks 1-15, the analogous `oaconvolve_nd(in1,in2,y,mode,axes,status)`
subroutine supports multidimensional overlap-add (details below). Existing
rank-1/2/3 `fftconvolve` functions return empty arrays on invalid input, while
`fftconvolve_nd` uses the status code and NaN-filled caller output.
In `valid` mode, one operand must dominate in all convolved dimensions.

### Rank-independent direct convolution and cross-correlation

`convolve_nd` and `correlate_nd` provide the direct and FFT variants of
SciPy's N-dimensional convolution and cross-correlation through **rank 15**.
Unlike `fftconvolve_nd`, they accept a `method=` argument with `'direct'`,
`'fft'`, or `'auto'` (default). The direct kernel computes the explicit
multidimensional convolution sum; the FFT kernel reuses SciFort's existing
separable-transform engine. `correlate_nd` conjugates the second operand
and reverses it on each selected convolution axis, like SciPy's
`signal.correlate`. `convolve_nd` does not conjugate either operand.

```fortran
use scifort_signal, only : dp, correlate_nd, convolution_shape_nd, &
    convolve_nd_success
implicit none
real(dp) :: image(8,6), pattern(3,2)
real(dp), allocatable :: response(:,:)
integer, allocatable :: dims(:)
integer :: status

! Set image and pattern from your data.
dims = convolution_shape_nd(shape(image),shape(pattern),mode='same')
allocate(response(dims(1),dims(2)))
call correlate_nd(image,pattern,response,mode='same',method='auto',status=status)
if (status /= convolve_nd_success) error stop 'correlation failed'
```

With a third channel dimension, use `axes=[1,2]` and a kernel of shape
`(height,width,1)` to process all channels independently. Rank-independent
inputs may be real or complex, and mixed operand types promote the result to
complex. The output rank equals the input rank. Arrays can be noncontiguous
Fortran sections, although compilers may create contiguous temporaries.
As with `fftconvolve_nd`, the implementation requires a C-interoperable `dp`
for its packed column-major views, and input/output arguments must not alias.

`convolution_shape_nd(shape(in1),shape(in2),mode,axes,status)` returns the
required extents, or an empty integer array and
`convolve_nd_invalid_input` for invalid shapes, modes or axes. Each
selected axis has output size `n+m-1` (`full`), `n` (`same`) or
`abs(n-m)+1` (`valid`), and valid mode requires one operand to dominate
all selected axes. Unselected axes permit equal extents or broadcasting
from singleton extents; in same mode they are center-cropped when needed.
Axes are **one-based**, and an omitted axes argument means all axes.
Duplicate, empty or out-of-range axis lists are rejected.

These two new routines use *transactional* errors: invalid shape, rank,
method, mode, axes, nonfinite operand values or empty input returns
`convolve_nd_invalid_input` **without modifying the caller's output**.
This differs deliberately from the NaN-filled error output of the older
`fftconvolve_nd` interface. The current `auto` heuristic chooses direct
when the approximate number of complex multiply-adds does not exceed
15,000, otherwise FFT. The heuristic is not a measured performance model.
For very large inputs, explicitly selecting FFT is generally preferable;
for short kernels, direct avoids FFT setup. Neither routine implements
`convolve2d`/`correlate2d`'s fill, wrap or symmetric boundary modes.

`test/test_signal_convolve_nd.f90` checks public interfaces, strided arrays,
rank 1/2/3/15, correlation conjugation, FFT/direct parity, mixed types,
invalid inputs and copy-out behavior. Run
`python test/compare_signal_convolve_nd_scipy.py --checked` for the 1008
independent, seeded SciPy comparisons using actual rank-15 operands
(including singleton dimensions representing lower-dimensional signals),
or omit `--checked` to test with `-O3`. The matching reference driver is
`tools/signal_reference/convolve_nd_driver.f90`.

### Multidimensional overlap-add convolution

For highly unequal signal and filter lengths, `oaconvolve_nd` uses blocked
FFT convolution rather than transforming the entire long signal at once.
It uses the **same call signature, output layout, argument promotion, axis
selection, and error codes** as `fftconvolve_nd`:

```fortran
use scifort_signal, only : dp, oaconvolve_nd, convolution_status_success
real(dp) :: samples(4096,3), kernel(17,1), smoothed(4096,3)
integer :: status

call oaconvolve_nd(samples,kernel,smoothed,mode='same', &
    axes=[1],status=status)
if (status /= convolution_status_success) error stop 'oaconvolve_nd'
```

On sufficiently asymmetric selected axes, the implementation splits the
longer operand into blocks, **including multiple axes in the same call**.
A common zero-padded FFT layout and the unchanged operand's transformed
spectrum are reused for every block. Block results are accumulated directly
into the requested `full`/`same`/`valid` crop, avoiding a full-size temporary
convolution. Broadcast dimensions are handled without transforming them.
If different operands are longer on different axes, the routine chooses one
side to split and leaves other axes intact, preserving spectrum reuse. When
inputs have comparable lengths, or blocking would not save work, it uses the
ordinary full multidimensional FFT instead.

The 2D tiled implementation was checked in 232 seeded cases against SciPy 1.17,
plus the original 400 comparisons spanning ranks 1/2/3/5/15. A reproducible
before/after benchmark lives in `benchmarks/bench_signal_oaconvolve_tiling.f90`.
On one GNU Fortran 14.2.0 environment with `-O3`, two evaluations of a
513x257 input with a 5x5 kernel took 0.079 s before and 0.025 s after the
change; timings depend on input sizes and compiler.

The routine inherits `fftconvolve_nd`'s current C-interoperable `dp`
storage-view requirement. Caller arrays may be noncontiguous sections
(compiler copies in and out); inputs must not alias the output.

The native regressions are `test_signal_oaconvolve_assumed_rank.f90` and
`test_signal_oaconvolve_assumed_rank_public.f90`. The optional SciPy 1.17
comparison is in `tools/signal_reference/verify_oaconvolve_assumed_rank.py`
with the Fortran driver `oaconvolve_assumed_rank_driver.f90`. The current
400 seeded cases cover ranks 1, 2, 3, 5, and 15 and all four operand types,
with maximum absolute discrepancy 4.55e-14.

`lfilter` now supports rank-2 and rank-3 real and complex signals. Specify
`axis=1`, `2`, or `3` to filter along a particular dimension; by default it
filters the final dimension. Each orthogonal lane is filtered independently.
The optional `zi=` and allocatable `zf=` have the same shape as the input,
except that the selected axis is replaced by
`max(size(a),size(b))-1` (the number of delay states). Other `lfilter`
semantics and the rank-1 interfaces are unchanged.

```fortran
real(dp) :: frames(16,32,128), b(2), a(1)
real(dp), allocatable :: smoothed(:,:,:), end_states(:,:,:)

b = [0.5_dp,0.5_dp]
a = [1.0_dp]
smoothed = lfilter(b,a,frames,axis=3,zf=end_states)
! shape(end_states) = [16,32,1]
```

Verification: `test_signal_fftconvolve_nd.f90`, `test_signal_lfilter_axes.f90`,
and `test_signal_lfilter_rank3.f90` exercise these APIs. The comparison scripts
`verify_fftconvolve_nd.py`, `verify_lfilter_axes.py`, and
`verify_lfilter_rank3.py` respectively compare 1,000, 400, and 300 seeded
cases with SciPy 1.17, including selected axes, singleton broadcasting,
mixed real/complex inputs, filter histories, and final states. Array ranks
higher than three use the separate caller-output `lfilter_nd` and
`fftconvolve_nd` and `oaconvolve_nd` interfaces (described above and below). `test_signal_fftconvolve_assumed_rank.f90`
and `verify_fftconvolve_assumed_rank.py` additionally validate `fftconvolve_nd`
against direct convolution and 576 seeded rank-15 SciPy 1.17 cases.

### Assumed-rank IIR/FIR filtering (ranks 1-15)

`lfilter_nd(b, a, x, y, zi, zf, status, axis)` is a Fortran 2018
**subroutine** exported from `scifort_signal` and `scifort_signal_lfilter`.
It accepts real(dp), complex(dp), and all six mixed real/complex combinations
of numerator coefficients, denominator coefficients, and input data. The
caller supplies the output `y` with exactly the same rank and extents as
`x`; ranks **1 through 15** are supported. The optional `axis` is **one-based**
and defaults to the last dimension. The arrays `zi` and `zf`, when present,
have the input's shape with the filtered axis replaced by
`max(size(a),size(b))-1`. For any mixed-type call, the output and states
must be `complex(dp)`; otherwise they match the real/complex input type.

```fortran
real(dp) :: x(4,3,2,5,100), y(4,3,2,5,100)
real(dp) :: zi(4,3,2,5,2), zf(4,3,2,5,2)
integer :: status

zi = 0.0_dp
call lfilter_nd([0.25_dp,0.5_dp,0.25_dp], [1.0_dp,-0.4_dp], &
    x, y, zi=zi, zf=zf, axis=5, status=status)
```

The implementation normalizes coefficient vectors once per call and applies
one transposed direct-form-II kernel per independent channel. The full
multidimensional input does **not** need to be copied into a temporary
array for filtering; mixed calls with real signals and complex coefficients
use one reusable complex line buffer. Streaming is supported by passing a
previous call's `zf` to a subsequent call's `zi`, with appropriately shaped
output/state arrays. Input arrays are never modified.

Invalid coefficients, axis, rank or state shapes return
`lfilter_status_invalid_input`, with the available output and final state
values set to quiet NaNs. Zero-length filtered axes preserve their starting
states. Noncontiguous input/output array sections work via Fortran's
copy-in/copy-out of `contiguous` dummy arguments, so input/output arguments
must not overlap. Standard `iso_c_binding` pointer remapping requires
`dp == c_double` on the target compiler. All calculations are otherwise
Fortran-only; no new library dependency is introduced.

Verification: `test_signal_lfilter_nd.f90` and
`tools/signal_reference/verify_lfilter_nd.py`, and
`test/compare_signal_lfilter_nd_scipy.py` exercise ranks 5/15, all eight
real/complex combinations, streaming, noncontiguous sections, and SciPy 1.17
numerical agreement. The reference driver passes 240 checked and optimized
comparisons; `test/compare_signal_lfilter_nd_scipy.py` checks 300 additional
cases with varying coefficient lengths. Existing array-returning rank-1/2/3 `lfilter` APIs
are unaffected.

## Peak location and measurement

```fortran
real(dp) :: x(9)
type(peak_result) :: res
real(dp), allocatable :: p(:), width(:), heights(:), left(:), right(:)
integer, allocatable :: lb(:), rb(:)
integer :: status

x = [0.0_dp, 1.0_dp, 0.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, &
     0.0_dp, 2.0_dp, 0.0_dp]
res = find_peaks(x, prominence_min=1.0_dp, width_min=0.5_dp, &
    plateau_size_min=1)
! res%peaks = [2, 5, 8]
! res%prominences = [1, 2, 2]
! res%widths = [1, 3, 1]

call peak_prominences(x, res%peaks, p, lb, rb, status)
call peak_widths(x, res%peaks, width, heights, left, right, status)
```

`find_peaks` accepts the following optional **scalar** inclusive bounds:

- `height_min`, `height_max` for amplitude;
- `threshold_min`, `threshold_max` for the drop to both immediate neighbors;
- `distance` for minimum center spacing (integer >= 1); stronger peaks win;
- `prominence_min`, `prominence_max` for contour prominence;
- `width_min`, `width_max` for width in sample intervals;
- `plateau_size_min`, `plateau_size_max` for level-top length;
- `wlen` to restrict prominence searches (integer >= 2, even values rounded
  upward to an odd centered window, consistent with SciPy);
- `rel_height` for the width level, default 0.5 (>= 0).

### Position-dependent peak constraints

Each of the five bounded properties also accepts an **array** of `real(dp)`
values, with exactly `size(x)` elements. The ten optional keywords are
`height_min_samples`, `height_max_samples`, `threshold_min_samples`,
`threshold_max_samples`, `prominence_min_samples`, `prominence_max_samples`,
`width_min_samples`, `width_max_samples`, `plateau_size_min_samples`, and
`plateau_size_max_samples`. Bounds are looked up at the **one-based center
index** of each candidate peak, including flat-topped peaks. A scalar bound
and a sample-dependent bound on the same endpoint are intersected; both must
be satisfied. The arrays can be noncontiguous sections and plateau bounds
can be fractional real values. Distance remains a single scalar.

For example, to detect peaks that exceed a varying baseline by two units:

```fortran
real(dp), allocatable :: baseline(:), x(:)
type(peak_result) :: res
! ... allocate and populate x and baseline, both of the same length ...
res = find_peaks(x, height_min_samples=baseline + 2.0_dp, &
    prominence_min=1.0_dp, plateau_size_max=3)
```

The dynamic bounds use SciPy's inclusive `find_peaks` interval semantics.
A supplied array requests the associated properties, even if its values
exclude no peaks. Every array must be finite and have length `size(x)`;
otherwise `peak_status_invalid_input` is returned, with empty result arrays.
Sample-dependent plateau bounds, like scalar plateau bounds, must be
nonnegative. This extension preserves all existing positional scalar calls.

### Open-interval requests (property-only measurements)

In SciPy, passing `height=(None, None)` requests peak heights without excluding
any candidates. SciFort represents this behavior with optional logical flags:
`height_open`, `threshold_open`, `prominence_open`, `width_open`, and
`plateau_size_open`. Set a flag to `.true.` to request the corresponding
properties, while retaining all peak locations unless **other** constraints
filter them. An absent or `.false.` flag preserves the old behavior. Existing
scalar and per-sample bounds can be supplied at the same time; they still
filter candidates and the open flag guarantees the property is returned.

```fortran
type(peak_result) :: res
res = find_peaks(x, height_open=.true., prominence_open=.true., &
    width_open=.true., plateau_size_open=.true.)
! res%peaks, res%peak_heights, res%prominences, res%widths,
! res%plateau_sizes and related geometric properties are allocated.
```

`width_open=.true.` also returns `prominences`, `left_bases`, and `right_bases`,
matching SciPy `width=(None, None)`. These are logical keyword arguments,
not Python-style intervals; peak locations and bases remain **one-based**.
The public test `test_signal_peaks_open_intervals.f90` covers all 32 flag
combinations, mixed scalar/array bounds, empty and noncontiguous arrays, and
fractional widths. The driver `peak_open_intervals_driver.f90` and
`compare_signal_peaks_open_intervals_scipy.py` cover 2,400 independently
seeded SciPy 1.17 cases (157,239 values, maximum error `2.84e-14`) in both
checked and optimized peak-kernel builds.

Properties are held in `peak_result`. Arrays remain allocated with length zero
when the relevant bound was not requested. Prominences and bases are also
returned when width is requested, consistent with SciPy. Plateau edge indices
are inclusive, one-based. `left_ips`/`right_ips` are fractional **one-based**
positions; widths measure differences in those positions.

`peak_prominences` reports the prominence and integer bases for specified
one-based peak indices. `peak_widths` calculates interpolated contour
intersections for those peaks; it can optionally accept precomputed arrays
`prominences=`, `left_bases=`, `right_bases=` together. When none are supplied,
it computes them using `peak_prominences`.

The implementation uses explicit status codes (`peak_status_success=0`,
`peak_status_invalid_input=1`), rather than exceptions/warnings. Nonfinite
signals, nonfinite bounds, invalid peak indices, negative relative heights,
and bad windows are rejected. Empty and monotone signals return zero peaks.

**Current differences from SciPy:** The Fortran API uses explicit logical
open-interval flags instead of Python's `(None, None)` interval syntax.
`distance` uses a deterministic rightmost-on-tie
priority; SciPy's internal sort may select different but equally high peaks
when distances conflict. The distance filter currently uses a straightforward
O(k^2) priority selection for k candidate peaks; other passes are typically
linear in signal length, apart from prominence searches.

The geometric definitions follow the public contour-prominence and fractional
width specifications in the SciPy manual, not its source code. Comparison with
SciPy 1.17.0 included 700 original randomized signals and 2,200 additional
seeded cases with all ten position-dependent bounds (55,439 peak indices and
property values; maximum absolute numerical difference `2.84e-14`), covering
distance and threshold constraints, local windows, plateaus, prominences,
and widths. Exact tie behavior is the
known exception. Dedicated self-contained Fortran reference tests also run
without SciPy installed.

## Detrending

```fortran
real(dp), allocatable :: residual(:)
integer :: status

residual = detrend([1.0_dp, 3.0_dp, 4.0_dp, 8.0_dp], status=status)
! residual = [0.3, 0.1, -1.1, 0.7]
residual = detrend(x, kind='constant')
residual = detrend(x, kind='linear', breaks=[5])
```

`detrend` removes the mean (`kind='constant'`) or fits and subtracts the
least-squares straight line (`kind='linear'`, default). Linear detrending
supports an unsorted array of one-based segment-start breakpoints, including
1 and size(x)+1. Duplicate breakpoint values are ignored. Breakpoints do
not change constant detrending. Invalid arguments return a quiet NaN array
with `detrend_status_invalid_input=1`; valid empty input gives an empty
array. In addition to the original real 1-D API, overloads now accept real
and complex arrays of ranks 1–3. For rank 2/3, the optional `axis` uses
**one-based** dimension indices and defaults to the last dimension (like
SciPy's default `axis=-1`). `breaks` are also **one-based**, so SciPy's
zero-based `bp` values must be increased by 1. The constant mode ignores
breakpoints. Noncontiguous sections are accepted and outputs preserve shape.
All 720 seeded real/complex rank-1/2/3 comparisons with SciPy 1.17 passed;
maximum absolute difference was approximately 5.35e-15.

```fortran
real(dp) :: channels(100,4)
complex(dp) :: complex_cube(8,40,6)
real(dp), allocatable :: centered(:,:)
complex(dp), allocatable :: residuals(:,:,:)
! Input arrays must be filled before calling these functions.
centered = detrend(channels,kind='constant',axis=1)
residuals = detrend(complex_cube,kind='linear',axis=2,breaks=[21])
```

## Direct convolution and cross-correlation

```fortran
real(dp), allocatable :: c(:), r(:)
integer, allocatable :: lags(:)

c = convolve([1.0_dp, 2.0_dp, 3.0_dp], [4.0_dp, 5.0_dp], mode='same')
! c = [4, 13, 22]
r = correlate([1.0_dp, 2.0_dp, 3.0_dp], [4.0_dp, 5.0_dp])
! r = [5, 14, 23, 12]
lags = correlation_lags(3, 2, mode='full')
! lags = [-1, 0, 1, 2]
```

`convolve` and `correlate` provide generic real(dp) and complex(dp) 1D
implementations, in modes `'full'` (default), `'same'` (size of first signal),
and `'valid'` (full overlap, permitting either input to be longer). For complex
correlation, the *second* argument is conjugated, following SciPy. The
optional `method` accepts `'auto'` and `'direct'`, both performing a direct
O(n*m) calculation; `method='fft'` is not supported by these direct routines and
returns a nonzero status; use `fftconvolve` for FFT-based convolution. Invalid modes or zero-length inputs also return
`convolution_status_invalid_input` and a zero-length result. All status
arguments are optional and are `intent(out)`.

`correlation_lags` returns integer shifts; these are signed lag offsets,
not array indices. It follows SciPy's `same` lag-centering convention, which
can differ by one from the direct correlation's array trimming for odd first
lengths and even second lengths, just as in SciPy 1.17.

## Two-dimensional convolution and correlation

`convolve2d` and `correlate2d` accept rank-two arrays of `real(dp)` or
`complex(dp)` (including mixed real/complex input pairs). They compute
**direct** two-dimensional sums and share SciPy's `mode='full'` (default),
`'same'` (first input's shape), and `'valid'` (complete overlap) options.
`'valid'` requires one operand to be at least as large as the other **on both
axes**; crossed dimensions are rejected. The `boundary=` argument accepts
`'fill'` (default), `'wrap'` (periodic), and `'symm'` (symmetric extension).
`fillvalue=` defaults to zero and applies only with `boundary='fill'`.

```fortran
real(dp) :: image(3,3), kernel(2,2)
real(dp), allocatable :: filtered(:,:)
integer :: status

image = 1.0_dp
kernel = 0.25_dp
filtered = convolve2d(image, kernel, mode='same', boundary='symm', &
    status=status)
```

Complex correlation conjugates the **second** operand. For compatibility
with SciPy 1.17's operand-swapping convention, `correlate2d` with
`mode='valid'` and a larger **second** complex input also conjugates the
result relative to the corresponding `full` slice. SciPy centers `same`
correlation one element to the right/down for even-sized kernels compared
with `same` convolution; SciFort follows this 2D convention. Inputs with
zero dimensions, unsupported modes/boundaries, or invalid `valid` shapes
return an empty `(0,0)` result and
`convolution_status_invalid_input` (`1`) when `status=` is requested.

The routines have O(N1*N2*M1*M2) runtime, use no FFT, and do not support
higher-rank convolution, batching, or sparse kernels. Other SciFort 1D
convolution functions are unchanged. With mixed real/complex operands,
`fillvalue=`, if supplied, must be a `complex(dp)` scalar.

To validate against SciPy 1.17, compile the standalone test driver with
`kinds.f90`, `signal_convolve.f90`, and `signal_convolve2d.f90` and run:

```text
python tools/signal_reference/compare_signal_convolve2d_scipy.py --exe path/to/convolve2d_driver
```

The seeded test set contains 1,600 real-and-complex comparisons across all
three modes and boundary rules, both `valid` operand orderings, even and
odd kernels, and nonzero complex fill values. A separate Fortran regression
is in `test/test_signal_convolve2d.f90`.

## Testing and limitations

Self-contained FPM test programs:

```text
test/test_signal_peaks.f90
test/test_signal_detrend.f90
test/test_signal_convolve.f90
```

For an isolated dependency-free build, compile `src/kinds.f90` followed by
`src/signal_convolve.f90`, `src/signal_detrend.f90`,
`src/signal_peaks.f90`, and `src/signal.f90`; link the three test programs
against those objects. The overall project still requires its normal FPM
dependencies to run the full test suite.

An optional cross-language regression can be run when Python, NumPy, SciPy,
and GNU Fortran are available:

```text
python tools/signal_reference/verify.py
```

The script compiles three test-only drivers in a temporary directory and
checks 700 peak/property cases, 450 detrending cases, and 900 complex/real
convolution/correlation/lag cases against the installed SciPy version. No
SciPy runtime dependency is introduced for SciFort itself.

## Savitzky-Golay polynomial smoothing and derivatives

```fortran
real(dp), allocatable :: coeff(:), smooth(:), rate(:)
integer :: status

coeff = savgol_coeffs(5, 2, deriv=1, use='conv', status=status)
! coeff = [0.2, 0.1, 0, -0.1, -0.2]
smooth = savgol_filter(x, 7, 3)
rate = savgol_filter(x, 7, 3, deriv=1, delta=0.1_dp, mode='mirror')
```

`savgol_coeffs(window_length, polyorder, deriv, delta, pos, use, status)`
returns least-squares polynomial evaluation/derivative weights. `deriv`
defaults to zero, `delta` to one, and `use` to `'conv'` (reversed
coefficients); `'dot'` is also supported. Unlike SciPy's zero-based
`pos`, **Fortran `pos` is a one-based integer in `[1, window_length]`**;
omitting it selects the middle, including half-sample centering for even
windows. The derivative order may exceed the polynomial degree, in which
case the weights are zero.

`savgol_filter(x, window_length, polyorder, deriv, delta, mode, cval, status)`
accepts real(dp) rank-1 through rank-4 samples. **Ranks 2 through 4** also
accept the optional **one-based `axis=`**, defaulting to the last dimension.
Rank-1 through rank-4 **complex(dp)** samples are also supported as a
**SciFort extension**, with real polynomial weights applied to complex values
without discarding the imaginary component; the optional `cval` is complex
for complex inputs. (SciPy 1.17 instead casts complex input to real.)
These axis overloads preserve the input shape and process each unfiltered
channel independently. Center weights and the independent left/right
`'interp'` edge weights are precomputed **once per array call** and reused
across channels (also for complex input), rather than fitted again for every
line. The weights are not cached between separate calls.
The boundary mode defaults to `'interp'`,
which fits the first and last windows separately rather than padding them.
Other modes are `'mirror'` (edge sample excluded), `'nearest'`, `'wrap'`,
and `'constant'` with configurable `cval` (default zero). Odd and even
window lengths are supported. For `'interp'`, the window must not exceed
the signal length; other modes support shorter signals.

Weights are computed from a scaled polynomial Vandermonde matrix using
reorthogonalized QR and a triangular solve, not normal equations. This
avoids squaring the least-squares condition number, though very high-order
polynomial fits can remain ill-conditioned. The API returns a quiet NaN
array for invalid or numerically singular parameters, with
`savgol_status_invalid_input=1` or `savgol_status_singular=2` (success is 0).
No array-API broadcasting is provided. For ranks 5 through 15, use the
assumed-rank `savgol_filter_nd` interface described below.

```fortran
real(dp) :: images(32,64,8)
real(dp), allocatable :: smoothed(:,:,:)
complex(dp) :: complex_data(20,40)
complex(dp), allocatable :: complex_smoothed(:,:)
real(dp) :: stack(10,20,30,40)
real(dp), allocatable :: smoothed_stack(:,:,:,:)

! Smooth across dimension 2; preserve other dimensions.
smoothed = savgol_filter(images,7,3,axis=2)
complex_smoothed = savgol_filter(complex_data,5,2, &
    axis=1,mode='constant',cval=(0.0_dp,0.0_dp))
! Rank-4 smoothing, along dimension 3.
smoothed_stack = savgol_filter(stack,7,2,axis=3)
```

### Irregularly sampled Savitzky-Golay filtering (SciFort extension)

`savgol_irregular` is provided as `savgol_irregular` from
both `scifort_signal` and `scifort_signal_savgol`. Unlike SciPy's uniform-grid
`savgol_filter`, it fits a **weighted local polynomial directly to the actual
observation times**, and optionally evaluates its time derivatives at
independent query times. It is particularly useful for jittered sampling,
dropped observations, and heteroscedastic measurements.

```fortran
real(dp) :: sample_time(100), measurements(100), precision(100)
real(dp) :: query_time(150)
real(dp), allocatable :: smoothed(:), derivative(:)
integer :: status

! Fill time, measurements, weights, and query locations first.
smoothed = savgol_irregular(sample_time, measurements, 9, 3, &
    weights=precision, query=query_time, status=status)
derivative = savgol_irregular(sample_time, measurements, 9, 3, &
    deriv=1, weights=precision, query=query_time, status=status)
```

The arguments are `savgol_irregular(x, y, window_length, polyorder, deriv,
weights, query, status)`; both real(dp) and complex(dp) observed series `y`
are accepted. The strictly increasing and finite `x` locations need not be
uniformly spaced. `query` defaults to `x`, and all query locations must be
finite. Each output selects the `window_length` *nearest consecutive samples*,
with ties going to the earlier time. The endpoints are handled by shifting
this fitting window inside the available observations. The degree must be
`0 <= polyorder < window_length` and the window length must be at least two.
The `deriv` argument defaults to zero and refers to derivatives with respect
to the physical units of `x`, not sample number. For orders above `polyorder`
the derivative is zero if the local fit is otherwise valid.

Optional nonnegative `weights` are **relative inverse-variance/precision
weights**, not standard deviations; only their ratios matter. Zero weights
exclude samples locally, and insufficient positive-weight points can cause a
singular fit. Every window is solved using twice-reorthogonalized QR of a
centered and scaled weighted Vandermonde matrix, avoiding normal equations.
An invalid request returns NaNs and `savgol_status_invalid_input`; a
rank-deficient or numerically singular local window returns an entirely NaN
result with `savgol_status_singular`. The routine evaluates all requested
locations transactionally, rather than returning a partly valid array.

This **does not add irregular-grid support to** SciPy's `savgol_filter`:
`savgol_irregular` is an original SciFort extension. The
special case of regular spacing and an odd window matches SciPy's
`mode='interp'` at the original sample times. There is no SciPy runtime
dependency. The NumPy weighted-least-squares comparison and SciPy uniform-grid
check run with `python tools/run_signal_savgol_irregular.py`.

### Assumed-rank Savitzky-Golay filtering (ranks 1-15)

`savgol_filter_nd(x, y, window_length, polyorder, deriv, delta, mode,
cval, status, axis)` is a Fortran 2018 **subroutine**, exported from both
`scifort_signal` and `scifort_signal_savgol`. The input and caller-supplied
output can have any rank from 1 to 15, real(dp) or complex(dp); they must
have **identical rank, extents and type**. The selected `axis` is **one-based**
and defaults to the last axis. Existing array-returning `savgol_filter`
functions at ranks 1-4 remain available and now use this same internal
rank-independent path. The polynomial weights are computed once per call
and shared by all channels, including the interpolated boundary fits.

```fortran
real(dp) :: x(4,3,2,5,10), y(4,3,2,5,10)
integer :: status

call savgol_filter_nd(x, y, 5, 2, axis=5, status=status)
! y has the same rank and shape as x
```

The output must already have its final shape. On invalid arguments, the
routine returns `savgol_status_invalid_input`, initializing the available
output elements to quiet NaNs; a singular polynomial fit returns
`savgol_status_singular`. The input and output must not overlap. Both
contiguous arrays and noncontiguous Fortran array sections are accepted;
Fortran may create a copy-in/copy-out temporary for noncontiguous sections.
Internally, contiguous column-major storage is remapped using standard
`iso_c_binding` intrinsics, which require a C-interoperable double-precision
kind (`dp == c_double` on supported compilers); **no C runtime or external
library is needed**. Extents and array data are not broadcast. The maximum
standard Fortran rank is 15; large arrays are limited by available memory.

`lfilter_nd`, `sosfilt_nd`, `filtfilt_nd`, `sosfiltfilt_nd`, and
`fftconvolve_nd` support ranks 1-15. The array-returning `sosfilt`,
`filtfilt`, `sosfiltfilt`, and `fftconvolve` APIs still stop at rank 3.

Validation includes `test_signal_savgol_nd.f90`, 300 seeded rank-5/rank-15
SciPy comparisons in `test/compare_signal_savgol_nd_scipy.py`, and the
existing rank-1 through rank-4 tests and reference comparisons.

## Direct-form digital filters

```fortran
real(dp), allocatable :: y(:), zi(:), zf(:)
integer :: status

zi = lfilter_zi([0.1_dp, 0.2_dp], [1.0_dp, -0.5_dp])
y = lfilter([0.1_dp, 0.2_dp], [1.0_dp, -0.5_dp], x, &
    zi=zi*x(1), zf=zf, status=status)
```

The generic `lfilter(b, a, x, zi, zf, status)` implements 1D FIR/IIR filtering
using the transposed direct-form-II difference equation and automatically
normalizes coefficients by `a(1)`. It accepts real(dp), complex(dp), and
**mixed-type coefficient/input arrays** for ranks 1 through 3. If any of
`b`, `a`, or `x` is complex, the result and optional `zi`/`zf` states are
complex(dp); otherwise they remain real(dp). Mixed types use SciFort's
existing complex filtering implementation, without discarding imaginary
parts. Inputs `b` and `a` are coefficients in ascending **lag** order. Optional `zi` supplies
`max(size(a), size(b))-1` initial delay values (zero by default); optional
allocatable output `zf` receives final delay values, even without `zi`.
Passing `zf` into the next `lfilter` call as `zi` supports blockwise filtering.

`lfiltic(b, a, y, x, status)` computes the same direct-form-II initial state
from **past outputs and inputs**, rather than from a forward filtering run.
Both `y` and optional `x` are in *most-recent-first* order, corresponding to
`[y(-1),y(-2),...]` and `[x(-1),x(-2),...]`. The optional `x` defaults to
zero history; missing samples in either history are zero-padded and excess
samples are ignored. The returned state has length `max(size(a),size(b))-1`.
The coefficients are normalized by `a(1)`, including for complex filters.
For example:

```fortran
real(dp), allocatable :: history_state(:), remaining(:)
! Most recent output was 0.4 and input was 0.3:
history_state = lfiltic([0.2_dp,0.1_dp], [1.0_dp,-0.5_dp], &
    [0.4_dp], [0.3_dp])
remaining = lfilter([0.2_dp,0.1_dp], [1.0_dp,-0.5_dp], &
    [1.0_dp,2.0_dp], zi=history_state)
```

Both real(dp) and complex(dp) coefficient/history forms are supported.
Invalid or nonfinite parameters return `lfilter_status_invalid_input=1` and
NaN state values, just as the existing `lfilter_zi` API does.

The generic `lfilter_zi(b, a, status)` computes delay values for a steady
unit-step input. It returns `lfilter_status_no_steady_state=2` when the
transfer function has a zero/near-zero denominator at unit frequency;
otherwise status is 0. Empty coefficient arrays, zero `a(1)`, nonfinite
filter coefficients, or invalid initial-state lengths return
`lfilter_status_invalid_input=1` and a quiet NaN output or empty final
state. Rank-2/3 `lfilter` accepts optional one-based `axis=` (default last
axis), with the filtered dimension of `zi`/`zf` replaced by the number of
delay states. `lfiltic` and `lfilter_zi` still use 1D coefficient/history
interfaces. As in SciPy, second-order sections
are recommended for numerically sensitive high-order IIR filters.

For ranks 4-15, `lfilter_nd` provides the caller-output assumed-rank
interface described in **Assumed-rank IIR/FIR filtering (ranks 1-15)**
above, with caller-supplied initial and final states. The original rank-1/2/3
array-valued `lfilter` functions are unchanged.

## New filter regression tests

```text
test/test_signal_savgol.f90
test/test_signal_savgol_axes.f90
test/test_signal_savgol_rank4.f90
test/test_signal_lfilter.f90
test/test_signal_lfilter_axes.f90
test/test_signal_lfilter_rank3.f90
test/test_signal_lfilter_mixed.f90
test/test_signal_lfilter_nd.f90
test/test_signal_lfiltic.f90
```

For an isolated dependency-free build, also compile
`src/signal_savgol.f90` and `src/signal_lfilter.f90` before `src/signal.f90`.
An optional randomized comparison is available using a local SciPy install:

```text
python tools/signal_reference/verify_filters.py
python test/compare_signal_savgol_axes_scipy.py
python test/compare_signal_savgol_rank4_scipy.py
python test/compare_signal_savgol_nd_scipy.py
python test/compare_signal_lfilter_mixed_scipy.py
python test/compare_signal_lfilter_nd_scipy.py
```

The additional drivers compare 600 seeded rank-1/2/3 Savitzky-Golay,
400 rank-4 Savitzky-Golay, 300 rank-5/rank-15 Savitzky-Golay, and 540 mixed-type IIR cases against SciPy 1.17,
including selected axes, complex data, and IIR final states. The Savitzky-Golay
complex comparisons use separate SciPy real/imaginary references.
Each compiles its own temporary
Fortran driver and requires only gfortran, NumPy and SciPy. The first script
compiles a test-only Fortran driver in a temporary directory and runs
2,380 deterministic comparisons (440 coefficient, 600 polynomial-filter,
520 real and 520 complex `lfilter`, 150 real and 150 complex `lfilter_zi`
cases) against the installed `scipy.signal` public functions. This reference
harness requires Python, NumPy, SciPy and gfortran, but the numerical Fortran
modules have no new runtime dependency.

A CPU benchmark comparing independent fits for every channel to shared weights
is provided in `benchmarks/benchmark_savgol_multichannel.f90`. With 768 channels
of length 160, 21-point windows, cubic interpolation, and 30 repetitions,
GNU Fortran 14.2 with `-O3` measured about 0.53 s before and 0.11 s after
this change (approximately 4.8x faster). These times are machine-dependent.

## Rank-independent periodogram, Welch, CSD and coherence

SciFort also provides `periodogram_into`, `welch_into`, `csd_into`, and
`coherence_into` for caller-allocated real or complex rank-1–15 arrays with
selectable (one-based or negative) axes. These interfaces retain SciFort's
existing window, FFT, scaling, and mean/median spectral-estimation kernels,
and add mixed real/complex cross-spectra. See the dedicated
[rank-independent spectral guide](signal_spectral_nd.md) for examples,
expected shapes, status handling, and SciPy comparisons.

## Rank-independent STFT and spectrogram (ranks 1–14)

`stft_into` and `spectrogram_into` expose caller-provided time-frequency
arrays with a configurable one-based or negative processing axis and one
extra, trailing time axis. They support one-sided real and two-sided complex
spectra, window selection, detrending, and all five legacy spectrogram modes.
Legacy `istft_into` additionally accepts `time_axis` and `freq_axis` for
arbitrarily permuted input transform dimensions, retaining a selectable
sample axis in the output.

See the [rank-independent time-frequency guide](signal_timefreq_nd.md)
for dimensions, signatures, status handling, and the SciPy reference suite.

## Spectral analysis (SciPy-style real and complex 1D subset)

This increment provides `get_window`, `periodogram`, `welch`, `csd`,
`coherence`, `stft`, and `spectrogram` via `scifort_signal`. The segment
accumulation and spectral-estimator layout were informed by the MIT-licensed
[`spectrum`](https://github.com/jchristopherson/spectrum) library by Jason
Christopherson, with notice recorded in `CODE_PROVENANCE.md` and
`THIRD_PARTY_LICENSES.md`. These are separately implemented SciFort interfaces
and do not require the external `spectrum` package or BLAS/LAPACK.

```fortran
use scifort_signal, only : dp, welch, periodogram, csd, coherence, &
    stft, spectrogram, get_window, spectral_result, cross_spectral_result, &
    stft_result, spectrogram_result

type(spectral_result) :: psd
type(cross_spectral_result) :: cross
type(stft_result) :: z
type(spectrogram_result) :: spect

psd = welch(x, fs=200.0_dp, window='hann', nperseg=128, noverlap=64)
! psd%frequencies: nonnegative Hz; psd%power: one-sided PSD
cross = csd(x, y, fs=200.0_dp, nperseg=128)
z = stft(x, fs=200.0_dp, nperseg=128)
! z%values(frequency_index, time_index) is complex
spect = spectrogram(x, fs=200.0_dp, nperseg=128)
! spect%power(frequency_index, time_index) is a real PSD estimate
```

The public Fortran functions return derived types with allocatable arrays and
`status` (`spectral_status_success=0` or
`spectral_status_invalid_input=1`) instead of Python tuples/exceptions.
`periodogram`, `welch`, and `coherence` return `spectral_result` with
`frequencies(:)` and `power(:)`. `csd` returns a `cross_spectral_result` with
complex `power(:)`; cross-spectral sign follows `conjg(FFT(x))*FFT(y)`.
`stft` and `spectrogram` also return `times(:)` and two-dimensional
`values(:,:)` or `power(:,:)`, respectively. A first FFT bin corresponds to
zero frequency and Fortran indices are one-based.

- `get_window(name,n,fftbins,status,parameter,parameter2)` supports 23 named
  window families, including `chebwin`, `taylor`, and `kaiser_bessel_derived`.
  Default windows are periodic; set `fftbins=.false.` for symmetric output.
  Kaiser, Gaussian and other parameterized windows use `parameter` and,
  where needed, `parameter2`. KBD requires an even length and symmetric output.
  See `docs/signal_window_families.md` for details and direct functions.
- `periodogram` defaults to a boxcar window, constant detrending and PSD
  density scaling, and allows a larger zero-padded or shorter/truncated DFT.
- `welch` and `csd` default to a periodic Hann window,
  `nperseg=min(256,n)`, half-window overlap, constant detrending and
  one-sided density scaling. Optional `average='mean'` (default) or
  `'median'` uses the SciPy median-bias correction; density/spectrum scaling
  is supported. Shorter of two CSD signals is zero-padded.
- `coherence` is `abs(Pxy)**2/(Pxx*Pyy)` using Welch spectra. Bins with
  zero/underflowed auto power yield quiet NaN instead of division by zero.
- `stft` defaults to a periodic Hann window, 50% overlap, zero-valued boundary
  extension and trailing zero-padding, with amplitude (`'spectrum'`) or PSD
  (`'psd'`) scaling. The `boundary` option supports `'zeros'`, `'even'`,
  `'odd'`, `'constant'`, and `'none'` (the Fortran equivalent of Python `None`);
  each extension prepends and appends `nperseg/2` samples, rounded down.
  `padded=.true.` adds zeros after the extension to complete the final frame.
  `padded` is a logical flag.
- `spectrogram` defaults to the SciPy legacy periodic Tukey(0.25) window,
  overlap `nperseg/8`, constant detrending and density scaling. It returns
  PSD for the rank-one derived-type interface; the `spectrogram_into`
  overloads additionally support complex, magnitude, angle, and phase modes.
  Explicit `window='tukey'` defaults to Tukey(0.5) unless a parameter is given.

The real-valued variants accept **real(dp) 1D signals** and produce
**one-sided frequency bins** by default. Set `return_onesided=.false.`
for full two-sided `nfft`-bin spectra on real input, with FFT ordering from
zero through positive and negative frequencies and a negative Nyquist bin
for even `nfft`. One-sided power estimates double interior bins; two-sided
estimates do not. The complex-valued overloads accept **complex(dp)** arrays
and always return two-sided FFT-order bins. Supported detrending choices are
`'constant'`, `'linear'`, and `'none'`; SciPy's Python `detrend=False` maps
to `'none'` in Fortran. Named parameterized and explicit sample windows are
supported by all six estimators; the caller-output assumed-rank adapters
support selectable axes. `coherence` has no `return_onesided` option in its
SciPy legacy API. Mixed real/complex rank-one CSD and coherence are now supported; see
[Mixed cross-spectral inputs](signal_mixed_cross.md). Exponential averaging
is not implemented. For invertible dual-window analysis use the
`short_time_fft` derived type documented above.
The FFT is O(n log n) for power-of-two lengths and uses an O(n**2) direct DFT
for other sizes. For large non-power-of-two transforms a mixed-radix FFT would
be significantly faster. Window and normalization conventions were tested
against SciPy 1.17.0; minor rounding differences at mathematically zero
periodogram/coherence bins can be magnified by normalization.

Testing:

```text
python tools/signal_reference/verify_spectral.py
python tools/signal_reference/verify_spectral_complex.py
python tools/signal_reference/verify_fftconvolve.py
```

The test-only Python scripts build independent Fortran comparison drivers,
require SciPy/NumPy and GNU Fortran, and check 1,640 real spectral cases,
960 complex spectral cases, and 1,000 FFT/overlap-add convolution cases
against SciPy.
The project itself adds **no new FPM dependencies**. Full `fpm test` still
requires the project's normal external dependencies.


## Inverse short-time Fourier transforms and overlap-add checks

`istft` complements the existing `stft` algorithm. It reconstructs a sampled
signal from the complex STFT matrix using an inverse DFT followed by
window-weighted overlap-add. Normalization divides by the sum of overlapping
**squared windows**. The routine supports ordinary real one-sided spectra and
complex two-sided spectra, odd-length and zero-padded FFTs, periodic named, parameterized, or
explicit window samples, `spectrum`/`psd` scaling, sampling frequency and
optional removal of zero-boundary extensions.

```fortran
use scifort_signal, only : dp, stft, stft_result, istft, &
    istft_status_success, check_cola, check_nola
real(dp) :: x(32)
real(dp), allocatable :: y(:), time(:)
type(stft_result) :: z
integer :: j, status

x = [(sin(real(j, dp)), j=1,32)]
z = stft(x, nperseg=8, noverlap=4)
call istft(z%values, y, status, nperseg=8, noverlap=4, times=time)
if (status /= istft_status_success) error stop 'istft'
! size(y) == 32, reconstructed to rounding precision
```

A named window can use the same three trailing shape arguments as the
forward STFT: `window_parameter` (e.g. Gaussian width or Kaiser beta),
`window_parameter2` (e.g. general Gaussian scale or Taylor sidelobe level),
and `window_coefficients` (general cosine series). All existing positional
arguments retain their order. For example, to reconstruct an STFT calculated
with a Gaussian window of standard deviation 2.5 samples:

```fortran
z = stft(x, window='gaussian', window_parameter=2.5_dp, &
    nperseg=8, noverlap=6)
call istft(z%values, y, status, window='gaussian', &
    window_parameter=2.5_dp, nperseg=8, noverlap=6)
if (status /= istft_status_success) error stop 'Gaussian ISTFT'
```

The inverse requires the same window weights, hop length, FFT size and scaling
as the forward transform to reconstruct the original signal. For inverse
STFT, explicit window arrays still use `window_values`, and cannot be mixed
with a named or parameterized window. Invalid window-family parameters or
nonfinite weights return `istft_status_invalid_input`.

The generic `istft(zxx,samples,status,...)` distinguishes **real** and
**complex** allocatable output arrays. Real output assumes a one-sided
Hermitian spectrum; complex output assumes a two-sided spectrum. To return
complex samples from a one-sided STFT, pass `input_onesided=.true.`. For
real output, `input_onesided=.false.` is invalid; use a complex output
instead. The frequency axis is dimension 1 and frame axis is dimension 2.
`nperseg`, `noverlap`, and `nfft` follow SciPy's inference rules, including
odd segment lengths when explicitly supplied. The `window='hann'` default
uses the periodic window from `get_window`; `window_values(:)` accepts
explicit weights as an alternative. Named-window `check_cola` and `check_nola`
accept the same `window_parameter`, `window_parameter2`, and
`window_coefficients` options, for example
`check_nola('kaiser',8,6,window_parameter=4.0_dp)`. `boundary=.true.` (the default) removes
`nperseg/2` extension samples from **both** ends. `times(:)` is an optional,
allocatable vector of reconstructed sample times. As in SciPy, an STFT padded
at the end can produce extra reconstructed zeros when the original sample
length is not a multiple of the hop.

Status values are `istft_status_success=0`,
`istft_status_invalid_input=1` and
`istft_status_noninvertible=2`. A noninvertible status means some sample
positions have overlapping window-square weight at or below `1e-10`:
SciFort still returns the overlap-add output but cannot guarantee those
samples were reconstructible. SciPy uses a warning for the same condition.
Inversion assumes the forward transform was not detrended, or the detrended
trend need not be recovered.

`check_cola` and `check_nola` are separate generic logical functions, accepting
either a window name or a real window vector. Each requires `nperseg` and
`noverlap`, supports optional `tol`, and returns optional input-validation
`status`. `check_nola` is the condition required for weighted overlap-add
inversion; `check_cola` tests the more restrictive constant-overlap sum.
For example, `check_cola('hann',8,4)` and `check_nola('hann',8,4)` are both
true, while `check_nola('hann',8,0)` is false.

The new `test/test_signal_istft_parameter_windows.f90` regression exercises
real and complex round trips, rank-14 input, noncontiguous array sections,
and unchanged caller outputs after invalid-window errors. Run independent
SciPy 1.17 comparisons of 192 coefficient-driven real and complex cases
(5,613 samples, maximum absolute error `6.57e-15`) in checked and optimized
GNU Fortran builds without FPM using:

```sh
python tools/run_signal_istft_parameter_windows.py
```

The independent regression programs are `test_signal_istft.f90` and
`test_signal_istft_public.f90`. The optional reference scripts
`test/compare_signal_istft_scipy.py` and
`test/compare_signal_overlap_checks_scipy.py` compile Fortran binary-stream
drivers and compare their results to the public SciPy APIs. These three
algorithms are implemented without new Fortran dependencies.

## Complex two-sided spectral analysis

The same generic names accept `complex(dp)` input and produce full two-sided
spectra in standard FFT order: frequencies start at zero, reach positive
frequencies, then continue with negative frequencies. For even `nfft`, the
Nyquist bin is negative, matching `numpy.fft.fftfreq`. All six spectral
functions `periodogram`, `welch`, `csd`, `coherence`, `stft`, and `spectrogram`
support complex input; `csd` and `coherence` currently require **both** input
arrays to be complex or both to be real.

```fortran
use scifort_signal, only : dp, periodogram, csd, stft, spectral_result, &
    cross_spectral_result, stft_result
complex(dp) :: x(64), y(64)
type(spectral_result) :: p
type(cross_spectral_result) :: xy
type(stft_result) :: z

p = periodogram(x, fs=64.0_dp)
xy = csd(x, y, nperseg=32, average='median')
z = stft(x, nperseg=16, boundary='zeros')
! size(p%frequencies) == 64; two-sided complex-input result
```

The optional logical `return_onesided` argument is accepted on complex
`periodogram`, `welch`, `csd`, `stft`, and `spectrogram` calls but cannot
force one-sided output: complex inputs always produce full spectra, just as
SciPy does.
A SciPy warning is not emitted by the Fortran API. The routines support
constant, linear, or no detrending of complex values, mean/median Welch
averaging, PSD/spectrum scaling, unequal-length CSD zero padding, and the
existing two STFT boundary modes. Complex `spectrogram` returns a real PSD.

Independent `test/test_signal_spectral_complex.f90` and
`tools/signal_reference/verify_spectral_complex.py` cover both FFT parities,
complex phase, padding, derivative-free detrending, scaling and spectral
estimation. The SciPy 1.17 comparison script passed 960 randomized cases.
Near-zero spectral bins can have unstable coherence ratios because machine
roundoff dominates numerator and denominator; the comparison excludes these
mathematically undefined bins.

## FFT and overlap-add convolution

`fftconvolve` performs zero-padded FFT linear convolution. `oaconvolve`
uses overlap-add with a reused kernel FFT for long/short filter combinations;
for similar input lengths it uses a single FFT. Both are generic for real and
complex 1D arrays and accept `mode='full'`, `'same'`, or `'valid'`.

```fortran
use scifort_signal, only : dp, fftconvolve, oaconvolve
real(dp), allocatable :: a(:), b(:), filtered(:)
filtered = fftconvolve(a, b, mode='same')
filtered = oaconvolve(a, b, mode='same')
```

Empty arrays, unsupported modes, and nonfinite data return an allocated
zero-length result and `convolution_status_invalid_input`, if requested.
This intentionally differs from SciPy, which may propagate nonfinite values
through its FFT with a warning. Both implementations use a radix-2 FFT, zero
padding to a power of two, and inverse transforms normalized by the transform
length. Large-input performance will differ from SciPy's pocketfft backend.
The array-returning `oaconvolve` supports 1D inputs; `oaconvolve_nd`
accepts multidimensional inputs through rank 15, just like `fftconvolve_nd`.
The legacy `fftconvolve` overloads cover ranks 1-3.
Neither implementation offers FFT worker counts.

The standalone `test/test_signal_fftconvolve.f90` validates against direct
convolution, including long-signal overlap-add. An optional 1,000-case
SciPy 1.17 comparison is in
`tools/signal_reference/verify_fftconvolve.py`.

## Cascaded second-order filtering and zero-phase filtering

The `scifort_signal` module provides the following generic `real(dp)` and
`complex(dp)` filtering functions, in addition to `lfilter`. They accept
rank-one, rank-two, or rank-three inputs:

```fortran
real(dp) :: sos(2,6), x(100)
real(dp), allocatable :: zi(:,:), zf(:,:), y(:), zero_phase(:)
integer :: status

! Each row is [b0,b1,b2,a0,a1,a2] with a0 exactly one.
sos(1,:) = [0.5_dp,0.5_dp,0.0_dp,1.0_dp,-0.3_dp,0.0_dp]
sos(2,:) = [0.5_dp,0.5_dp,0.0_dp,1.0_dp,-0.2_dp,0.0_dp]

zi = sosfilt_zi(sos, status=status)
y = sosfilt(sos, x, zi=zi*x(1), zf=zf, status=status)
zero_phase = sosfiltfilt(sos, x, padtype='odd', status=status)
```

- `sosfilt(sos, x, zi, zf, status, axis)` runs cascaded transposed
  direct-form-II sections. `sos` has shape `(n_sections, 6)`. For rank-one
  input, `zi`/`zf` have shape `(n_sections, 2)`; at higher rank they have
  shape `(n_sections, *shape(x))` with the filtered signal axis replaced by
  length 2. `zi` is optional (default zero state); `zf` is an optional
  allocatable output for continued streaming. Section index is **always**
  the first dimension of `zi` and `zf`.
- `sosfilt_zi(sos, status)` returns section states corresponding to the
  steady-state response to a unit step, accounting for the gain of each
  preceding section. Filters with a pole at one do not have a finite unit-step
  steady state and return `sos_status_no_steady_state`.
- `filtfilt(b, a, x, padtype, padlen, status, axis)` runs `lfilter`
  forwards and backwards, initializing each pass from scaled `lfilter_zi`.
- `sosfiltfilt(sos, x, padtype, padlen, status, axis)` performs the same
  operation using cascaded second-order sections; this is generally preferable
  for higher-order IIR filtering.

For rank-two and rank-three signals, `axis` is **one-based** and defaults to
`rank(x)`; all other dimensions are independent channels. For example:

```fortran
real(dp) :: sos(2,6), records(8,100,3)
real(dp), allocatable :: filtered(:,:,:), zf(:,:,:,:)

sos(1,:) = [0.5_dp,0.5_dp,0.0_dp,1.0_dp,-0.3_dp,0.0_dp]
sos(2,:) = [0.5_dp,0.5_dp,0.0_dp,1.0_dp,-0.2_dp,0.0_dp]

! Filter the 100 samples along the second axis.
filtered = sosfilt(sos,records,axis=2,zf=zf)
! shape(zf) == [2,8,2,3]
filtered = sosfiltfilt(sos,records,axis=2,padtype='even')
```

The rank-one overloads retain their existing argument order. In ranks two and
three, the new `axis` argument follows `status`, so previous keyword-based
calls are unaffected. Unlike `sosfilt`, the forward/backward routines do not
expose filter-state arguments, as in SciPy.

### Assumed-rank SOS filtering (ranks 1-15)

`call sosfilt_nd(sos, x, y, zi=..., zf=..., axis=..., status=...)` is a
Fortran 2018 subroutine exported by `scifort_signal` and `scifort_signal_sos`.
It filters real or complex inputs of ranks **1 through 15** along an optional
one-based axis (default: last). `sos` always has shape `(nsections,6)` with
`a0` exactly one for every section. Real and complex SOS coefficients may be
combined with either real or complex input data; any complex operand requires
**complex output and state arrays**. `y` must have exactly `shape(x)`.

Initial `zi` and final `zf` are optional caller-owned arrays. The normal
SOS state layout is `(nsections, d1, ..., d_axis=2, ..., d_rank)`, where
`d1 ... d_rank` are the data extents, with the selected axis replaced by 2.
That layout works for input ranks **1 through 14** (states have one extra
leading dimension). There is also a packed one-dimensional representation:
flatten the conceptual state layout in Fortran column-major order. It can be
used **at any rank**, and is necessary for rank-15 inputs because Fortran
2018 does not permit rank-16 arrays. Each packed buffer contains
`2*nsections*product(shape(x))/size(x,axis)` elements when the filtered
axis is nonempty (or, more generally,
`2*nsections*product(shape(x) with axis removed)`). `zi` and `zf` may
independently use packed or normal layouts as long as their shapes are valid.

```fortran
real(dp) :: sos(2,6), samples(4,3,8,2,5), filtered(4,3,8,2,5)
real(dp) :: zi(2,4,3,2,2,5), zf(2,4,3,2,2,5)
integer :: status

! SOS coefficients have already been initialized with a0 = 1.
zi = 0.0_dp
call sosfilt_nd(sos,samples,filtered,zi=zi,zf=zf,axis=3,status=status)
! shape(zf) == [2,4,3,2,2,5]
! In the next block, use the previous zf as the new zi.
```

On invalid coefficient, rank, state shape, or axis, the subroutine returns
`sos_status_invalid_input` with NaN output and final states. Nonfinite
initial states are invalid; input samples themselves may be nonfinite, as
in ordinary SciPy filtering. An empty filtering axis leaves each state
unchanged. Strided actual array sections are supported through compiler
copy-in/out for `contiguous` dummy arguments. The algorithm reuses one
`(nsections,2)` work array across all lanes; mixed-type filtering does not
allocate a complex copy of the entire input array.

The existing rank-1/2/3 array-returning `sosfilt` generic is unchanged.
Verification: `test/test_signal_sos_nd.f90` and
`test/compare_signal_sos_nd_scipy.py` (320 seeded rank-2/5/15 numerical
comparisons for all four real/complex combinations, including final states).

The rank-1/2/3 `filtfilt` generic also supports `method='gust'` and an
optional `irlen=` for Gustafsson forward/backward initial states. For this
method padding arguments are ignored, including on very short signals;
rank-2/3 calls retain their existing one-based `axis=` argument. The
standalone `filtfilt_gust` function additionally accepts mixed
real/complex coefficient/data combinations. See
[**Gustafsson zero-phase filtering**](signal_filtfilt_gust.md) for API,
examples, numerical behavior and validation details.

For zero-phase second-order-section filtering, set `method='gust'`
on `sosfiltfilt(sos,x,...)` (ranks 1-3) or use the rank-1 generic
`sosfiltfilt_gust(sos,x,irlen,status)`, including mixed real and complex
coefficients/data. This SOS-native alternative avoids the numerical hazards
of multiplying high-order polynomial coefficients. It is implemented by
SciFort, not a parameter supported by SciPy's `sosfiltfilt`. See
[the Gustafsson method documentation](signal_filtfilt_gust.md).

`filtfilt` defaults to `padtype='odd'` and a padding length of
`3*max(size(a),size(b))` samples on **each** side. For `sosfiltfilt`, the
SciPy-style default is
`3*(2*n_sections+1-min(count(sos(:,3)==0),count(sos(:,6)==0)))`.
The optional `padtype` supports `'odd'`, `'even'`, `'constant'`, and `'none'`;
`'none'` represents SciPy's `padtype=None`. An explicit `padlen=0` also
performs unpadded filtering. `padlen` must be strictly shorter than the length
of the **selected filtering axis**, not the overall array size.

### Assumed-rank forward/backward filtering (ranks 1-15)

Use `call filtfilt_nd(b,a,x,y,axis=...,padtype=...,padlen=...,status=...)`
or `call sosfiltfilt_nd(sos,x,y,...)` from `scifort_signal`
or `scifort_signal_filtfilt_nd`.
These Fortran 2018 subroutines extend the existing zero-phase
forward/backward digital filtering functions to real or complex arrays
of ranks **1 through 15**. They preserve the original rank-1/2/3 array-returning
calls, whose optional trailing arguments now also accept Gustafsson states.

```fortran
real(dp) :: records(4,3,24,2,2), zero_phase(4,3,24,2,2)
real(dp) :: sos(2,6)
integer :: code

! Provide a normalized, valid SOS matrix with sos(:,4) == 1.0_dp.
call sosfiltfilt_nd(sos,records,zero_phase,axis=3, &
    padtype='even',padlen=4,status=code)
call filtfilt_nd([0.25_dp,0.5_dp,0.25_dp], &
    [1.0_dp,-0.3_dp],records,zero_phase,axis=3, &
    padtype='odd',padlen=4,status=code)
```

The caller owns `y`, whose **rank and shape must match** `x`. The
optional `axis` is one-based and defaults to `rank(x)`; remaining axes
identify independent channels. Both variants support `'odd'` (default),
`'even'`, `'constant'`, and `'none'` padding and a caller-set `padlen`;
`'none'` corresponds to SciPy's `padtype=None`. With the default
padtype, `padlen` is `3*max(size(a),size(b))` for BA filters and
`3*(2*nsections+1-min(count(sos(:,3)==0),count(sos(:,6)==0)))`
for SOS. A valid filtered axis must contain more than `padlen` samples.
An explicit `padlen=0` or `padtype='none'` removes padding, but the
steady-state initialization still must exist for padding. Negative
Python-style axes are not supported. Both `filtfilt_nd` and
`sosfiltfilt_nd` support `method='gust'` and positive `irlen`; the SOS
variant uses direct-form-II cascade states rather than expanding SOS rows
into a high-order transfer polynomial. Extremely short, ill-conditioned
SOS state systems can return `filtfilt_status_numerical_failure` with a
quiet-NaN result. This SOS Gustafsson method is a SciFort extension (SciPy's
`sosfiltfilt` currently supports only padding). For example:

```fortran
call filtfilt_nd(b, a, records, zero_phase, axis=3, &
    method='gust', irlen=16, status=code)
call sosfiltfilt_nd(sos, records, zero_phase, axis=3, &
    method='gust', irlen=16, status=code)
```

Real and complex coefficients and data may be mixed: `filtfilt_nd`
accepts all eight real/complex `(b,a,x)` permutations;
`sosfiltfilt_nd` accepts all four `(sos,x)` combinations. If any operand
is complex, the caller's output must be `complex(dp)`. There are no
initial/final-state outputs because filtering is forward/backward, not
streaming. Input/output need not be explicitly contiguous: compiler
copy-in/out allows sliced array sections. Internally, memory is remapped
in column-major order and one 1D lane is processed at a time; thus no
whole-array complex promotion is necessary. The implementation reuses
the existing SciFort 1D routines to preserve SciPy-style pad and
steady-state semantics. As with those routines, the result is undefined
if caller input and output overlap.

Return codes are `filtfilt_status_success`,
`filtfilt_status_invalid_input`, `filtfilt_status_no_steady_state`
(padding only), and `filtfilt_status_numerical_failure` (Gustafsson
linear-system solve). Failed validation and per-channel failures initialize
the entire selected output to quiet NaNs. Empty outer dimensions are permitted when the
selected filtering axis and filter are valid. On all-scalar arguments,
the subroutine returns the invalid-input code and NaN output.

Verification: `test/test_signal_filtfilt_nd.f90`,
`test/test_signal_filtfilt_nd_public.f90`, the checked standalone
signal suite, and `python test/compare_signal_filtfilt_nd_scipy.py`,
which exercises **684** SciPy 1.17 reference cases (ranks 2/5/15,
multiple axes, padding modes, all type combinations, and both filter
representations). The maximum observed absolute error was `3.33e-16`
for these padded cases using GNU Fortran 14.2. In addition,
`test/test_signal_filtfilt_nd_gust.f90` and
`tools/signal_reference/compare_signal_filtfilt_nd_gust_scipy.py` cover
960 Gustafsson cases (48,384 output samples, including 240 direct SciPy
comparisons and 720 independent complex-coefficient references), passing
in checked and optimized builds with maximum error <= `2.24e-16`.
No external Fortran libraries are required.

The original array-returning filtering interfaces support ranks one
through three with matching coefficient/input types. Mixed real/complex
coefficient promotion is available through `sosfilt_nd`,
`filtfilt_nd` and `sosfiltfilt_nd` but not through the
array-returning `sosfilt`, `filtfilt`, or `sosfiltfilt` overloads. Negative Python-style axis indices are not implemented. Inputs with no SOS rows are rejected by `sosfilt` and
`sosfiltfilt`. `sosfilt_zi` requires at least one row, and, unlike SciPy,
rejects an empty matrix. Each SOS section requires `a0` **exactly equal to
one** (as SciPy's SOS routines do); use a normalized coefficient matrix.

All routines return an allocated output, with optional integer status:
`sos_status_success`, `sos_status_invalid_input`,
`sos_status_no_steady_state` and the analogous `filtfilt_status_*` codes.
Invalid filters or padding return quiet NaN arrays rather than raising a
Python exception. The operation is deterministic and needs no external
numerical dependency.

An independent implementation of the standard cascaded transposed DF-II
recurrence and endpoint-reflected forward/backward filtering was verified
against the documented SciPy 1.17 behavior. The included optional reference
harness checks 2,100 randomized cases covering real/complex filter states,
step initialization, and all four padding types:

```text
python tools/signal_reference/verify_sos.py
```

The fixed regression programs are `test_signal_sos.f90`,
`test_signal_filtfilt.f90`, `test_signal_sos_axes.f90`, and
`test_signal_filtfilt_axes.f90`. The additional SciPy 1.17 reference command:

```text
python tools/signal_reference/verify_filter_axes.py
```

checks 480 seeded rank-two/three SOS and forward/backward cases, including
real/complex data, all axes, padding modes, input states, and final states.
The largest observed absolute discrepancy was about `5.55e-16`.

## Windowed FIR and Butterworth IIR filter design

These routines are available from `use scifort_signal`. The algorithms are
independent SciFort implementations, with no new dependencies:

```fortran
use scifort_signal, only : dp, firwin, butter, butter_result, &
    freqz, sosfreqz, sos2tf, frequency_response
real(dp), allocatable :: taps(:)
type(butter_result) :: lowpass
type(frequency_response) :: h

taps = firwin(31, [0.20_dp], window='hamming')
lowpass = butter(4, 0.30_dp, btype='lowpass')
h = sosfreqz(lowpass%sos, worN=256)
lowpass = sos2tf(lowpass%sos)
! lowpass%b, lowpass%a: transfer-function polynomial coefficients
! h%frequencies: radians/sample; h%response: complex transfer function
```

`firwin(numtaps, cutoff, pass_zero, window, fs, scale, status)` returns
windowed-sinc coefficients, with the cutoff frequencies strictly between
zero and `fs/2` (`fs=2` by default). `cutoff` may be a `real(dp)` scalar or
an array of strictly increasing cutoff frequencies. `pass_zero` accepts
logical `.true.`/`.false.` or the strings `'lowpass'`, `'highpass'`,
`'bandpass'`, and `'bandstop'`. If omitted, a single cutoff gives lowpass and
multiple cutoffs give bandstop. Multiband designs accept **two or more**
cutoff frequencies, including odd counts (e.g. three or five edges); bands
alternate between pass and stop at each edge. Whether Nyquist is in a passband
depends on the number of edges and `pass_zero`; an even tap count is rejected
when Nyquist is in a passband. The
symmetric window name defaults to `'hamming'`; supported names include
`'boxcar'`, `'hann'`, `'hamming'`, `'blackman'`, and `'bartlett'` (plus the
other names accepted by `get_window`). `scale=.true.` normalizes the first
passband's reference gain to one. SciFort reports invalid input if that
normalization is singular, rather than returning SciPy's occasional NaN.

### Advanced FIR filter design and minimum-phase conversion

The `scifort_signal` module also exports `firwin2`, `firls`, and
`minimum_phase`. All three return allocatable real(dp) coefficient arrays.
Supply the optional `status` argument to inspect input validation; on errors,
the returned array is empty. The `fir_advanced_status_success`,
`fir_advanced_status_invalid_input`, and `fir_advanced_status_singular` constants
are public. No additional FPM dependency is needed.

```fortran
use scifort_signal, only : dp, firwin2, firls, minimum_phase
real(dp), allocatable :: h(:), h_ls(:), h_min(:)

! Frequency and gain breakpoints, with fs=2 (Nyquist=1).
h = firwin2(31, [0.0_dp,0.3_dp,0.4_dp,1.0_dp], &
    [1.0_dp,1.0_dp,0.0_dp,0.0_dp], window='hamming')

! Least-squares linear-phase FIR, with relative pass/stop weights.
h_ls = firls(31, [0.0_dp,0.3_dp,0.4_dp,1.0_dp], &
    [1.0_dp,1.0_dp,0.0_dp,0.0_dp], weight=[1.0_dp,3.0_dp])

! Homomorphic FFT-cepstrum conversion: ceil(31/2) coefficients.
h_min = minimum_phase(h, method='homomorphic')
```

- `firwin2(numtaps, freq, gain, nfreqs, window, antisymmetric, fs, status)`
  uses linear interpolation of gains and an inverse real frequency transform.
  Frequency points must run from 0 to `fs/2` in nondecreasing order;
  an interior point may be repeated exactly twice to model a gain jump.
  The interpolation mesh defaults to one more than the next power of two
  at least as large as `numtaps`; `nfreqs` must exceed `numtaps`.
  The default symmetric window is `'hamming'`; `window='none'` bypasses
  windowing. With `antisymmetric=.false.` (default), odd/even taps yield
  type I/II filters, respectively. `antisymmetric=.true.` yields
  type III/IV filters. The usual zeros at DC and/or Nyquist are enforced.
- `firls(numtaps, bands, desired, weight, fs, status)` uses analytically
  integrated weighted cosine normal equations for **odd** tap counts.
  `bands` contains nonoverlapping `(low, high)` pairs and `desired` contains
  one nonnegative desired gain per band endpoint. `weight`, if specified,
  is one nonnegative weight per band. The cosine integrals use
  cancellation-resistant midpoint formulas for narrow bands. A Cholesky solve is tried first; if the normal-equation
  matrix loses numerical positive definiteness, a dependency-free symmetric
  Jacobi eigensolver computes a rank-revealing pseudoinverse solution. Band
  weights are scaled by their maximum to avoid intermediate overflow, so
  equivalent weight scalings produce the same design. An all-zero weight
  vector still returns `fir_advanced_status_singular`. Unlike exact SciPy
  coefficient comparisons, nearly rank-deficient filters may have different
  equally accurate tap vectors because the singular-value cutoff differs;
  compare the weighted in-band response instead. The Gram formulation can
  still lose accuracy for very high orders or extremely narrow bands.
- `minimum_phase(h, method, n_fft, half, status)` provides the FFT-cepstrum
  `'homomorphic'` default and modified discrete-Hilbert `'hilbert'` method.
  `half=.true.` (default) returns `(size(h)+1)/2` coefficients with
  approximately the square-root magnitude response; `half=.false.`
  returns `size(h)` coefficients with approximately the original magnitude,
  and is only supported for `'homomorphic'`. The FFT length defaults to
  the next power of two at least `200*(size(h)-1)`; explicit `n_fft` must
  be at least `size(h)`. Non-power-of-two lengths are supported through the
  existing slower DFT fallback. The Hilbert method can be numerically
  sensitive to cancellation at near-zero spectral samples, so comparisons
  with SciPy may differ by several `1e-3` for some filter designs even
  though the standard transform construction is the same. Prefer
  `'homomorphic'` for general FIR filters.

The independently written numerical tests are in
`test/test_signal_fir_advanced.f90`, `test/test_signal_firls_rank.f90`,
`tools/signal_reference/verify_fir_advanced.py`, and
`tools/signal_reference/verify_firls_conditioned.py`. The Python drivers
compare frequency-sampling filters of all four types, weighted least-squares
filters with both well-conditioned and nearly rank-deficient systems, and
homomorphic/Hilbert minimum-phase transforms against installed SciPy. No SciPy source code
is included in the Fortran implementation.

### Equiripple FIR design (`remez`)

`scifort_signal` now exports `remez`, independently implemented using a dense
frequency grid and barycentric polynomial Remez exchange. It produces
least-maximum-error FIR filters with the usual type-I/II symmetric or type-III/IV
antisymmetric linear phase. Supported `filter_type` values are `'bandpass'`
(default), `'hilbert'`, and `'differentiator'`. Unlike `firwin2`/`firls`, there is
**one desired gain per band**, not one per frequency edge. The `bands` argument
holds disjoint, increasing low/high frequency pairs; `weight`, if present, holds
one strictly positive weight per band.

```fortran
use scifort_signal, only : dp, remez, remez_status_success
integer :: status
real(dp), allocatable :: h(:)

! 31-tap equiripple lowpass; frequencies measured with fs=2.
h = remez(31, [0.0_dp,0.3_dp,0.4_dp,1.0_dp], [1.0_dp,0.0_dp], &
    weight=[1.0_dp,2.0_dp], fs=2.0_dp, status=status)
if (status /= remez_status_success) error stop 'Remez did not converge'
```

Signature: `remez(numtaps, bands, desired, weight, filter_type, maxiter,
grid_density, fs, status, iterations, ripple)`. The optional `fs` defaults to
**1** (matching SciPy `remez`, unlike some other SciFort signal routines),
`maxiter` to 25, and `grid_density` to 16. In differentiator mode the desired
response increases linearly with **cycles/sample** (`frequency/fs`), with a
relative-error weight within nonzero-gain bands. The `iterations` and `ripple`
optional outputs describe the last successful exchange and the maximum sampled
weighted error. Use a larger `grid_density` to reduce frequency-grid error.

The returned allocatable coefficients are empty on failure. Status constants
are `remez_status_success`, `remez_status_invalid_input`, and
`remez_status_no_convergence`. Numerical exchange can fail on high-order,
very narrow, highly disparate, or otherwise ill-conditioned specifications;
this is reported as nonconvergence rather than returned as an unreliable
filter. Tap vectors need not exactly equal SciPy's due to differences in
extremal-frequency grids and numerical stopping criteria. Prefer checking
weighted passband/stopband ripple when comparing implementations.

The exchange, barycentric interpolation, and Chebyshev reconstruction use
**at least 30-digit intermediate arithmetic** if the compiler offers a suitable
real kind (`selected_real_kind(30,400)`); otherwise they fall back to `dp`.
The public interface and returned FIR taps remain double precision. Higher
intermediate precision can be slower but considerably reduces cancellation
in high-order Hilbert and differentiator designs. The alternating-extrema
selection uses parity-indexed prefix maxima rather than quadratic searches.
Eight previously failing high-order cases now return finite coefficients with
weighted ripple close to SciPy 1.17; one of the cases does not converge in
SciPy itself. Very ill-conditioned designs can still return
`remez_status_no_convergence`, especially on compilers without an extended
real kind. These changes introduce no external dependencies.

`test/test_signal_remez.f90` contains fixed SciPy-based test cases, and
`tools/signal_reference/verify_remez.py` runs an independent randomized
SciPy comparison. No SciPy code is bundled or required for FPM builds.

`butter(order, cutoff, btype, fs)` supports integer prototype order 1..20,
with scalar or length-one cutoff for lowpass/highpass and a length-two array
for bandpass/bandstop. Frequencies are measured in the units of `fs`
(default `fs=2`). The returned `butter_result` provides `b`, `a`, `sos`,
and `status`. Its SOS array has shape `(nsections,6)`, with columns
`[b0,b1,b2,a0,a1,a2]`, as expected by `sosfilt`/`sosfiltfilt` and
`sosfreqz`. For high-order and narrowband filters, **prefer SOS**: transfer-
function polynomial coefficients may be badly conditioned even when the
section representation remains accurate. Section ordering and scaling can
be different from SciPy's; the transfer function is the intended invariant.

`freqz(b, a, worN, whole, fs, include_nyquist, w)` calculates the complex
transfer-function response. `a` is optional, default `[1]`. The routine
is generic over matching real(dp) or complex(dp) arrays. `worN` defaults
512; `whole=.true.` generates the full circle. `include_nyquist=.true.`
includes the Nyquist endpoint for the half-circle grid. The optional `w`
array supplies arbitrary physical frequencies (including negative values)
and overrides grid generation. Frequencies default to radians/sample with
`fs=2*pi`; all arrays of response values use a `frequency_response` result
with `frequencies`, `response`, and `status` components.

### Analog and zero-pole-gain frequency responses

`freqs(b, a, worN, w)` evaluates the **analog** transfer function at
`s = i*w`. Coefficients are in **descending** powers of `s`, unlike digital
`freqz` coefficients, which are in ascending powers of unit delay. `a` is
optional and defaults to `[1]`. Arrays `b` and `a` may both be `real(dp)` or
both be `complex(dp)`. The optional `w` overrides automatic grid selection
and contains angular frequencies in **rad/s**, with no `fs` argument. If
`w` is absent, `worN` (default 200) requests a logarithmic grid around
characteristic filter frequencies, selected from numerator and denominator
roots.

`freqs_zpk(z, p, k, worN, w)` evaluates the same analog response using
zeros, poles, and scalar gain `k`, i.e. `k * product(i*w-z)/product(i*w-p)`.
`freqz_zpk(z, p, k, worN, whole, fs, w)` instead evaluates a **digital**
transfer function with `z = exp(i*2*pi*w/fs)`. Digital `worN` defaults to 512,
`fs` defaults to `2*pi`, and `whole=.true.` requests a full frequency circle.
All zero/pole inputs may be real or complex independently, with real or
complex gain. Empty zero/pole arrays are allowed. All three functions
return the existing `frequency_response` type, with the `frequencies`,
`response`, and `status` components and `design_status_*` codes.

```fortran
type(frequency_response) :: analog, digital
analog = freqs([1.0_dp], [1.0_dp, 2.0_dp, 2.0_dp], worn=200)
analog = freqs_zpk([real(dp) ::], &
    [cmplx(-1.0_dp,1.0_dp,dp),cmplx(-1.0_dp,-1.0_dp,dp)], &
    1.0_dp,w=[0.0_dp,0.5_dp,1.0_dp,2.0_dp])
digital = freqz_zpk([1.0_dp,-1.0_dp],[0.5_dp],0.25_dp,fs=1000.0_dp)
```

By default, **analog** frequency grids are logarithmic and positive; the
**digital** grids are evenly spaced. Explicit `w` allows negative values and
an empty array. Invalid coefficients/frequencies return
`design_status_invalid_input`; poles exactly on an evaluated imaginary-axis
or unit-circle point return `design_status_singular` with NaN at that point.
For analog transfer-function automatic grids, an independent iterative
polynomial root estimator selects characteristic frequency *ranges*; like
other root-finding methods, decade boundaries can differ from SciPy for
clustered or repeated high-order roots. Supplying explicit `w` bypasses this
step and avoids that ambiguity. These functions do not implement SciPy's
optional plotting callbacks or multidimensional coefficient broadcasting.

The deterministic test is `test/test_signal_analog.f90`. An optional
1,000-case SciPy 1.17 comparison (16,519 evaluated frequencies) is provided by
`python tools/signal_reference/verify_analog.py`; the largest observed
absolute response difference was about `3.8e-9` in the randomized sample.
No new FPM dependency was added.

## Normalized analog lowpass prototypes

The public `scifort_signal` module also exports five SciPy-compatible analog
lowpass prototype generators; they can be imported directly from
`scifort_signal_analog_prototypes` (or `scifort_signal_ellipap` for the elliptic prototype) when a small dependency set is preferred.
The routines return a common `analog_prototype_result` type with
`complex(dp)` allocatable `%zeros(:)` and `%poles(:)`, `real(dp)` `%gain`,
and integer `%status`:

| Routine | Arguments | Prototype |
|---|---|---|
| `buttap(n)` | nonnegative order | Butterworth, unit cutoff |
| `cheb1ap(n, rp)` | passband ripple in dB | Chebyshev I, cutoff at `-rp` dB |
| `cheb2ap(n, rs)` | stopband attenuation in dB | Chebyshev II, cutoff at `-rs` dB |
| `besselap(n, norm)` | `phase`, `delay`, or `mag` | Bessel/Thomson, default `phase` |
| `ellipap(n, rp, rs)` | passband ripple and stopband attenuation in dB | Elliptic/Cauer |

`%zeros` lists **finite** transmission zeros (the empty array for
Butterworth, Chebyshev I, and Bessel). The analog transfer function is
`gain * product(s-zeros) / product(s-poles)`, where `s` is measured in
radians/second. Pole and zero order is not semantically significant.
Order zero returns a constant prototype; Chebyshev I returns
`10**(-rp/20)` gain, as SciPy does for that special case.

```fortran
use scifort_signal, only : dp, cheb2ap, analog_prototype_result, prototype_success
type(analog_prototype_result) :: prototype

prototype = cheb2ap(5, 40.0_dp)
if (prototype%status /= prototype_success) error stop 'invalid analog prototype'
! prototype%zeros, prototype%poles, and prototype%gain can feed ZPK transforms.
```

The Bessel normalization options have distinct interpretations:
`phase` makes the lowpass asymptotic gain coefficient 1 (the SciPy/MATLAB
default); `delay` gives unit low-frequency group delay; and `mag` places
the -3 dB amplitude point at 1 rad/s. The Bessel poles are calculated
using Campos-Calderon starting points and simultaneous Aberth-Ehrlich
refinement of the scaled complex modified-Bessel-K function. This avoids
high-order polynomial evaluation cancellation. For validated **double precision**
accuracy, `besselap` currently accepts orders 0 through **20**;
higher orders return `prototype_too_expensive`. The trigonometric prototypes
accept orders 0 through **64**. The elliptic prototype accepts orders 0 through **20**, with `0 < rp < rs <= 100` dB; it reuses SciFort's existing AGM/Jacobi elliptic-filter routines. Inputs with invalid order, normalization, or
attenuation return `prototype_invalid_input`, and Bessel failures return
`prototype_no_convergence`. The Fortran APIs do not raise exceptions for
ordinary invalid parameters. No new FPM dependencies are introduced.

The native regression is `test/test_signal_analog_prototypes.f90`; the
independent comparison with SciPy is run with
`python tools/signal_reference/compare_signal_analog_prototypes_scipy.py`.
The original 339 deterministic reference cases cover poles, finite zeros, gains, and
continuous-time frequency responses. The extended regression is
`test/test_signal_bilinear_ellipap.f90`; run
`python tools/signal_reference/compare_signal_bilinear_ellipap.py` to compare elliptic
prototypes and Bessel orders 11--20 with SciPy.

## Analog lowpass frequency transformations

SciFort now implements eight analog-frequency transformation functions in
`scifort_signal` (and `scifort_signal_frequency_transforms`):

| Transfer-function form | Zero-pole-gain form | Frequency substitution |
|---|---|---|
| `lp2lp(b,a,wo)` | `lp2lp_zpk(z,p,k,wo)` | `s / wo` |
| `lp2hp(b,a,wo)` | `lp2hp_zpk(z,p,k,wo)` | `wo / s` |
| `lp2bp(b,a,wo,bw)` | `lp2bp_zpk(z,p,k,wo,bw)` | `(s**2+wo**2)/(bw*s)` |
| `lp2bs(b,a,wo,bw)` | `lp2bs_zpk(z,p,k,wo,bw)` | `(bw*s)/(s**2+wo**2)` |

`wo` and `bw` are optional, positive real angular frequencies, both
in **rad/s**, defaulting to 1.0. Transfer-function coefficients are ordered
by **descending powers of s**, not descending powers of `z**(-1)` as in
SciFort's digital filtering functions. The BA functions accept both real and
complex coefficient arrays and return `analog_tf_transform_result` with
complex allocatable `%b`, `%a`, and integer `%status`. Input leading zeros are
ignored, and the transformed denominator is normalized to have unit leading
coefficient. The implementation composes rational polynomials directly,
without factoring high-degree polynomials into roots. The computational
order limit is 32; exceeding it yields `frequency_transform_too_expensive`.

The ZPK functions accept complex arrays of zeros and poles and a real or
complex scalar gain. They return `analog_zpk_transform_result` with `%zeros`,
`%poles`, `%gain`, and `%status`. As in SciPy, poles must be at least as
numerous as finite zeros. Transformations requiring reciprocals reject zeros
or poles at the origin. SciFort preserves complex gain ratios even for
non-conjugate root sets, a deliberate mathematical extension beyond SciPy's
real-prototype assumptions. Conjugate-paired real-filter prototypes match
SciPy's results.

```fortran
use scifort_signal, only : dp, lp2bp, analog_tf_transform_result, &
    frequency_transform_success
type(analog_tf_transform_result) :: transformed

transformed = lp2bp([1.0_dp], [1.0_dp,2.0_dp,1.0_dp], &
    wo=20.0_dp, bw=4.0_dp)
if (transformed%status /= frequency_transform_success) error stop 'lp2bp'
! transformed%b and transformed%a hold complex(dp) analog coefficients.
```

Input validation uses `frequency_transform_invalid_input`, and no exceptions
are thrown by the Fortran API for ordinary invalid parameters. The native
regression is `test/test_signal_frequency_transforms.f90`; the SciPy 1.17
numerical comparison uses
`python tools/signal_reference/compare_signal_frequency_transforms_scipy.py`.

## Analog-to-digital bilinear transformations

`bilinear(b, a, fs)` and `bilinear_zpk(z, p, k, fs)` apply
`s = 2*fs*(z-1)/(z+1)` to analog transfer functions. Both are
exported from `scifort_signal`, or directly from `scifort_signal_bilinear`.
`fs` is a **positive sampling frequency in Hz**, not an analog cutoff in
rad/s. **No frequency prewarping** is performed. `bilinear` defaults to
`fs=1`; `bilinear_zpk` requires `fs`. The former accepts real or complex
coefficient arrays in descending powers of `s`, and returns
`analog_tf_transform_result` with complex `%b(:)` and `%a(:)` in descending
powers of `z` (equivalently digital `z**(-1)` coefficient order), normalized
so that `%a(1) = 1`. The returned coefficient arrays have equal lengths,
`max(size(b),size(a))`. Analog coefficients of degree above 32 return
`frequency_transform_too_expensive`.

`bilinear_zpk` takes complex zero and pole arrays and a real or complex gain,
and maps finite poles/zeros using `(2*fs+s)/(2*fs-s)`. Missing zeros at
infinity become zeros at `z=-1`. Improper systems and finite singularities
at `s=2*fs` return `frequency_transform_invalid_input`. The ZPK routines
preserve the **complete complex gain ratio** even for non-conjugate poles and
zeros; SciPy 1.17 intentionally takes the real part of that ratio. For
conjugate-paired real-coefficient systems the outputs agree with SciPy.
Neither SciFort variant performs implicit frequency prewarping.

```fortran
use scifort_signal, only : dp, bilinear, analog_tf_transform_result
real(dp) :: sampling_rate
integer :: status
type(analog_tf_transform_result) :: digital

sampling_rate = 100.0_dp
digital = bilinear([1.0_dp], [1.0_dp, 10.0_dp], fs=sampling_rate)
status = digital%status
if (status /= 0) error stop 'bilinear transform failed'
! digital%b and digital%a can be passed to digital filtering routines.
```

The native test is `test/test_signal_bilinear_ellipap.f90`. The comparison
`python tools/signal_reference/compare_signal_bilinear_ellipap.py` tests
complex polynomial transfer functions, ZPK root mappings and gain identities,
elliptic prototypes, and Bessel normalization against SciPy 1.17. No new
runtime or FPM dependencies are introduced.


## Bessel filter design: analog and digital

The `bessel(order, cutoff, btype, analog, norm, fs)` generic covers
lowpass, highpass, bandpass, and bandstop Bessel/Thomson filters through
prototype order 20. Single-edge designs accept a scalar cutoff; band
filters take two ascending cutoff frequencies. `norm='phase'` is the default;
`'delay'` selects the maximally flat low-frequency group-delay normalization,
and `'mag'` selects a -3 dB analog magnitude cutoff. The underlying
`besselap` poles are frequency-transformed in ZPK form. Digital designs
prewarp the edges and apply the bilinear transform, as in SciPy; the
transform is not exactly group-delay-preserving at high digital frequencies.

```fortran
use scifort_signal, only : dp, bessel, iirfilter, butter_result, design_status_success

type(butter_result) :: filt

! 14th-order digital bandpass at sampling rate 200 Hz.
filt = bessel(14, [15.0_dp, 35.0_dp], btype='bandpass', fs=200.0_dp)
if (filt%status /= design_status_success) error stop 'invalid Bessel design'
! filt%sos(nsections,6) is recommended over direct-form b/a at high orders.

! Analog filter with 1 s low-frequency group delay at a 1 rad/s reference.
filt = bessel(8, 1.0_dp, norm='delay', analog=.true.)

! Default phase-normalized Bessel filter through the IIR dispatcher.
filt = iirfilter(6, [0.2_dp, 0.4_dp], ftype='bessel', btype='bandstop')
```

Digital designs return real `%b(:)`, `%a(:)`, and `%sos(:,:)` arrays using
the existing `butter_result`. For analog designs, `%sos` remains unallocated.
Near high-order bandstop notches, analog BA polynomial evaluation can become
severely ill-conditioned, even in SciPy. For robust analog designs,
`iirfilter_zpk(order, cutoff, ftype='bessel', btype=...)` and
`iirfilter_sos(order, cutoff, ftype='bessel', btype=...)` now accept
Bessel/Thomson filters through order 20 with phase normalization; use these
representations to avoid expanded high-order polynomials. Numerical comparisons are in
`tools/signal_reference/compare_signal_bessel_design_scipy.py` and the
native test is `test/test_signal_bessel_design.f90`.

## Analog Butterworth and Chebyshev filter design

The `butter`, `cheby1`, `cheby2`, `ellip`, `bessel`, and `iirfilter` generics
support the final optional argument `analog=.true.`. In this mode the cutoff
is an **analog angular frequency in rad/s**, rather than a digital frequency
relative to Nyquist. The implemented families are Butterworth, Chebyshev I
(passband ripple `rp` in dB), Chebyshev II (stopband attenuation `rs`
in dB), elliptic (both `rp` and `rs`), and Bessel/Thomson for lowpass, highpass,
bandpass, and bandstop transfer functions.

```fortran
use scifort_signal, only : dp, butter, cheby1, cheby2, iirfilter, &
    butter_result, freqs, frequency_response

type(butter_result) :: filt
type(frequency_response) :: response

! 3rd-order Butterworth analog lowpass with cutoff at 100 rad/s.
filt = butter(3, 100.0_dp, analog=.true.)
response = freqs(filt%b, filt%a, w=[50.0_dp,100.0_dp,200.0_dp])

! 4th-order Chebyshev I analog bandpass with 1 dB ripple.
filt = cheby1(4, 1.0_dp, [20.0_dp,80.0_dp], &
    btype='bandpass', analog=.true.)

! 5th-order inverse Chebyshev analog highpass with 50 dB rejection.
filt = cheby2(5, 50.0_dp, 15.0_dp, &
    btype='highpass', analog=.true.)

! Equivalent selection using the shared IIR dispatcher.
filt = iirfilter(5, 15.0_dp, ftype='cheby2', rs=50.0_dp, &
    btype='highpass', analog=.true.)
```

Results use the existing `butter_result`: on success `b(:)` and `a(:)`
contain **descending powers of the Laplace variable** `s`, as required by
`freqs`; unlike digital IIR designs, `sos` is **unallocated** in analog
mode. Only real-valued polynomial coefficients are returned. This analog
mode implements the transfer-function (`ba`) representation, not the analog
`zpk` or `sos` output formats supported by SciPy. Digital calls that omit
`analog` behave exactly as before and still return SOS as well as BA.

Orders 1 through 20 are accepted; bandpass and bandstop transforms produce
order `2*N`. Invalid cutoff arrays, zero/negative/nonfinite frequencies,
invalid ripple/attenuation, and `fs` supplied with `analog=.true.` return
`design_status_invalid_input` with unallocated output arrays. As in SciPy,
`fs` is meaningful only for **digital** designs. High-order transfer-function
polynomials can be ill-conditioned even when their roots are stable, so
compare coefficients at moderate orders and use a root-based representation
elsewhere for numerically delicate calculations. Analog minimum-order
selection through `buttord`/`cheb1ord`/`cheb2ord` and `iirdesign` is available,
as described below.

Implementation: `src/signal_analog_design.f90`. Fixed SciPy-reference
regressions: `test/test_signal_analog_design.f90`. Run an optional 1,200-case
SciPy 1.17 randomized comparison (all three families, four modes, and
orders 1 to 8):

```text
python tools/signal_reference/verify_analog_iir.py
```

That comparison passed all 1,200 cases (maximum scaled coefficient
error `4.55e-13`, maximum absolute frequency-response discrepancy
`4.69e-9`). There are no new FPM dependencies.

`group_delay(b, a, worN, whole, fs, w)` computes **digital filter group
delay** in samples, defined as the negative derivative of transfer-function
phase with respect to radians/sample. `b` and optional `a` (default `[1]`)
are real(dp) or complex(dp) coefficient arrays. It returns a
`group_delay_result` containing `frequencies`, real(dp) `delay`, and `status`.
The frequency-grid behavior matches `freqz`, except `include_nyquist` is not
a separate option: use explicit `w` to request the Nyquist frequency. The
algorithm evaluates the analytic derivative of the numerator and denominator
polynomials, avoiding finite differences of wrapped phase, and scales each
polynomial for coefficient-magnitude invariance. Near a numerator or
denominator zero, the corresponding delay is NaN and the overall status is
`design_status_singular=2`, unlike SciPy's convention of substituting zero
in singular bins. This distinction is intentional. Example:

```fortran
type(group_delay_result) :: delay_result
delay_result = group_delay([0.0_dp, 1.0_dp], worN=128)
! delay_result%delay = 1 sample at every frequency
```

`group_delay_zpk(z, p, k, worn, whole, fs, w)` adds stable zero-pole-gain
analysis for mixed real/complex root arrays, with no transfer polynomial
expansion. It also works as `group_delay(z, p, k, ...)`, including unequal
root counts and per-bin singularity diagnostics. See
[ZPK group delay](signal_group_delay_zpk.md) for the formula, conventions,
status codes, and reproducible numerical tests.

`group_delay_sos(sos, worN, whole, fs, w)` is a **SciFort extension** for
stable group delay of real(dp) and complex(dp) SOS cascades. It sums each
second-order section's analytic delay without forming a high-order
transfer-function polynomial. The `group_delay` generic also accepts a
rank-two SOS matrix directly. It uses the existing `group_delay_result`
type, frequency units, status codes, and singular-bin NaN convention.
See [SOS group delay](signal_group_delay_sos.md) for examples, input
validation and SciPy regression instructions.

`sos2tf(sos)` converts a real(dp) SOS matrix of shape `(nsections,6)`
into numerator and denominator polynomials using coefficient convolution.
It returns a `butter_result` with the `b`, `a`, and `sos` components populated
and `status=design_status_success`. An empty `(0,6)` cascade produces
`b=[1]`, `a=[1]`. As in SciPy, the polynomials have length `2*nsections+1`,
so an odd-order filter may include a trailing zero arising from a padded
first-order section. Section coefficients are not normalized by this routine;
nonunit `a0` is preserved. Malformed or nonfinite sections are rejected.
Only real coefficients are currently supported, and forming high-order
transfer-function polynomials can substantially degrade numerical stability.

`sosfreqz(sos, worN, whole, fs, w)` evaluates a real(dp) or complex(dp)
SOS cascade using the same conventions. The leading denominator of each
section must equal one. Evaluation uses polynomial Horner recurrences and
multiplies individual section responses. These routines do not include
SciPy's `plot` option, multidimensional coefficient broadcasting, arbitrary
filter-design families, or analog second-order-section evaluation. Invalid arguments
return `design_status_invalid_input=1` and unallocated result arrays;
frequency points at denominator zeros report `design_status_singular=2`.

To validate with SciPy 1.17 and GNU Fortran:

```text
python tools/signal_reference/verify_design.py
```

For SOS-to-polynomial conversion, `test/test_signal_sos2tf.f90` checks
reference coefficients, empty cascades, padded first-order sections, and
coefficient agreement with the other IIR designs. The optional
`python tools/signal_reference/verify_sos2tf.py` script compares 545 cases
against SciPy 1.17.0, including random sections and designed filters.

The self-contained `test/test_signal_design.f90` also verifies reference
coefficients, filter responses, argument validation, and real/complex
frequency evaluation without requiring Python. The optional Python script
covers 1,620 randomized SciPy reference cases, including wide bandstop
and narrow bandpass Butterworth designs.

## Polyphase multirate filtering and downsampling

```fortran
use scifort_signal, only : dp, upfirdn, resample_poly, decimate, &
    resample_status_success
real(dp), allocatable :: y(:)
integer :: status

! Upsample by three, filter, retain every second high-rate sample.
y = upfirdn([0.25_dp, 0.5_dp, 0.25_dp], x, up=3, down=2, status=status)

! Automatic zero-phase Kaiser(5) anti-alias filter, rational 3:2 rate change.
y = resample_poly(x, 3, 2, status=status)

! Explicit FIR downsampling (IIR is now the default).
y = decimate(x, 4, ftype='fir', zero_phase=.true., status=status)
```

- `upfirdn(h, x [, up, down, mode, cval, status])`: sparse polyphase
  computation (no explicit zero-stuffed array), supporting real/complex
  signals and real/complex FIR kernels. `up` and `down` default to 1.
  Boundary `mode` may be `'constant'` (default, zero unless `cval` is
  given), `'edge'`, `'wrap'`, `'reflect'`, `'symmetric'`, `'line'`,
  `'smooth'`, `'antisymmetric'`, or `'antireflect'`. The boundary modes extend original input samples, not the
  zero-stuffed signal.
- `resample_poly(x, up, down [, taps, window, beta, padtype, cval, status])`:
  rational resampling with a zero-phase low-pass FIR, one-based arrays,
  and `ceiling(n*up/down)` output samples. `taps=` provides the FIR
  impulse response *before* up-scaling. If omitted, a symmetric
  windowed-sinc FIR with `20*max(up,down)+1` taps is designed **after
  reducing the ratio by gcd**. The default `window='kaiser'` and
  `beta=5.0_dp`; non-Kaiser windows use the symmetric `firwin` interface
  (`hann`, `hamming`, `blackman`, `boxcar`, etc.). `padtype` accepts the
  `upfirdn` boundary modes and `'mean'`, `'median'`, `'minimum'`, and
  `'maximum'` (background subtraction and restoration). For complex samples,
  median/minimum/maximum use the same lexicographic complex ordering as
  NumPy. `cval` applies only to `'constant'` padding.
- `decimate(x, q [, n, ftype, zero_phase, status])`: FIR or IIR
  downsampling. SciPy's default `ftype='iir'` now uses an eighth-order
  Chebyshev-I lowpass SOS filter. For `ftype='fir'`, the default order is
  `n=20*q` with symmetric Hamming taps; `zero_phase=.true.` compensates
  for group delay through `resample_poly`, while false performs causal FIR
  decimation. The IIR design is described below.

Input arrays must be nonempty and factors positive; invalid arguments
return allocated zero-length arrays with `resample_status_invalid_input`.
The first two operations support complex(dp) input; both `upfirdn` and
`resample_poly` support complex FIR taps. `resample_poly` with real input
and complex `taps=` returns a complex(dp) array. With real taps and real
input it continues to return real(dp). FIR `decimate` uses real taps.

Example: multi-channel polyphase processing along the second dimension:

```fortran
real(dp) :: channels(4, 24)
real(dp), allocatable :: filtered(:,:)
! The 24-sample dimension is axis 2 in Fortran, or -1 from the end.
filtered = resample_poly(channels, 3, 2, axis=2, status=status)
! shape(filtered) == [4, 36]
```

Complex custom FIR taps preserve their imaginary parts:

```fortran
complex(dp) :: complex_taps(3)
complex(dp), allocatable :: complex_filtered(:,:)
complex_taps = cmplx([0.2_dp, 0.6_dp, 0.2_dp], &
                     [0.1_dp, -0.2_dp, 0.1_dp], dp)
complex_filtered = resample_poly(channels, 3, 2, taps=complex_taps, &
                                 padtype='median', axis=-1, status=status)
! shape(complex_filtered) == [4, 36]
```

The rank-2/3 overloads call the established one-dimensional kernels for each
independent channel; no new FPM dependencies are needed. The 53-case SciPy
reference driver and generator are at
`tools/signal_reference/{multirate_nd_driver.f90,check_multirate_nd.py}` and
`tools/generate_signal_multirate_nd.py`. The regression test is
`test/test_signal_multirate_nd.f90`. All 53 deterministic cases agree with
SciPy 1.17 to approximately `4.3e-15` (6,105 complex-valued samples).

The complex-tap and extra-padding regression driver is
`tools/signal_reference/{multirate_complex_driver.f90,check_multirate_complex.py}`;
`test/test_signal_multirate_complex_taps.f90` checks the Fortran interface
and invalid inputs. These tests cover all ranks 1-3, real/complex inputs,
and both anti-extension modes. SciFort uses a constant extension for a
single-sample signal in `antireflect` mode, avoiding an undefined slope.

**Current SciPy differences:** These APIs accept rank-one through rank-three
arrays. The rank-two and rank-three overloads add `axis` before `status`;
positive axes are **one-based** (Fortran convention) and negative axes count
from the last dimension (`axis=-1` selects the last). In the multidimensional
overloads, `upfirdn` and `decimate` default to the last axis, whereas
`resample_poly` defaults to the first axis, matching SciPy's respective
default positions. The allocatable-result function overloads stop at rank three. For ranks
1-15, use the separate caller-allocated assumed-rank `*_into` subroutines
below. `decimate` now supports caller-provided **real FIR taps** (`taps=`) and
**real IIR second-order sections** (`sos=`) instead of the default designed
filters. Direct Python DLTI objects, general complex IIR coefficients and
transfer-function `(b,a)` inputs remain unsupported. See
[Custom FIR and SOS decimation](signal_decimate_custom_filters.md).
The Kaiser parameter is provided
through `beta=` instead of Python's `('kaiser', beta)` tuple. Operations
use a direct polyphase multiply-accumulate kernel, avoiding a zero-stuffed
array while retaining linear work in the necessary filter taps and samples;
performance may differ from SciPy's native polyphase implementation.

### Assumed-rank multirate transforms (ranks 1-15)

Fortran 2018 callers can use four new subroutine generics through
`scifort_signal`:

- `call resample_into(x, y, num [, window, t, output_times, domain, axis, status,
  frequency_window, complex_window])` (`domain=` and `complex_window=` are
  for complex input only).
- `call resample_poly_into(x, y, up, down [, taps, window, beta, padtype, cval, axis, status])`.
- `call upfirdn_into(h, x, y [, up, down, mode, cval, axis, status])`.
- `call decimate_into(x, y, q [, n, ftype, zero_phase, axis, status, taps, sos])`.
  Complex `sos=` matrices select a separate complex-output overload.

Here `x` and `y` are *assumed-rank* arguments; the caller must allocate or
otherwise provide `y` with the same rank and unchanged non-transform extents.
For example:

```fortran
use scifort_signal, only : dp, resample_poly_into, resample_status_success
real(dp) :: channels(2, 1, 9, 3)
real(dp) :: filtered(2, 1, 14, 3)
integer :: status
call resample_poly_into(channels, filtered, 3, 2, axis=3, status=status)
if (status /= resample_status_success) error stop 'polyphase resampling failed'
```

For rank-15 arrays this caller-allocated subroutine approach does not add a
new array dimension, so the rank-15 Fortran limit is supported. Positive
axis indices are one-based, and negative axes count from the last dimension.
`resample_into` and `resample_poly_into` default to axis 1; `upfirdn_into` and
`decimate_into` default to the last axis. Real and complex signals are
supported, including all four tap/signal combinations for `upfirdn_into` and
real-input/complex-tap `resample_poly_into` yielding complex output.

These subroutines return a nonzero `status` without changing **any** elements
of `y` if the axis or output shape is invalid, or if a lane fails. They use
Fortran-order temporary buffers and reuse the existing one-dimensional
kernels, so memory use temporarily increases by approximately two complex
arrays of input/output size. Strided input and output sections are supported.
The existing rank-1/2/3 functions with allocatable results remain unchanged.

`src/signal_multirate_nd_into.f90` contains the shared lane traversal, and
`tools/generate_multirate_assumed_rank.py` regenerates only the portable
`select rank` pack/unpack cases. The native public-module regression is
`test/test_signal_multirate_assumed_rank.f90`. Randomized SciPy comparisons
for ranks 4, 5, 6, 8, 14, and 15 use
`tools/signal_reference/{multirate_assumed_rank_driver.f90,verify_multirate_assumed_rank.py}`.
No additional FPM dependencies are required.

Independent reference verification is available in
`tools/signal_reference/verify_resample.py`. This compares 1,120 randomized
real/complex 1D cases against SciPy 1.17.0, including rational-ratio
reduction, extended input boundaries, explicit FIR taps, automatic Kaiser
and standard windows, zero-phase and causal FIR decimation. The maximal
absolute error in that set was under 8e-15. This reference utility uses
SciPy only for validation, not as a SciFort runtime dependency.


## Hilbert analytic signals (1D and 2D)

`hilbert(x, n, status)` computes the analytic signal for a **real(dp)** 1D
input. The optional positive `n` selects the FFT size, truncating or
zero-padding the samples; by default, the transform size equals `size(x)`.
Its real part recovers the original signal (or the truncated/zero-padded
signal), and its imaginary part is the periodic Hilbert transform. The
result is a complex(dp) allocatable vector.

`hilbert2(x, n, status)` computes SciPy's **single-orthant** 2D analytic
signal for a real(dp) rank-two input, returning a complex(dp) matrix. The
optional `n` accepts either `[nrows, ncols]` or a scalar for a square FFT.
The 2D operation doubles positive-frequency bins along **both** axes;
it is not a transform of each original row or column independently, and
its real part need not equal the original matrix.

```fortran
use scifort_signal, only : dp, hilbert, hilbert2
real(dp) :: x(4) = [0.2_dp, 1.0_dp, -0.5_dp, 2.0_dp]
real(dp) :: img(2,3)
complex(dp), allocatable :: analytic(:), analytic_image(:,:)

analytic = hilbert(x)              ! four complex values
analytic = hilbert(x, n=8)         ! pad to eight samples
analytic_image = hilbert2(img, n=[4,6])
```

Both routines use `hilbert_status_success=0` and
`hilbert_status_invalid_input=1`. Invalid transform sizes or nonfinite
samples return allocated zero-length results. Explicit positive FFT lengths
also permit zero-length input arrays and produce zero-padded results.
Only rank-one and rank-two inputs are exposed respectively; there is no
multi-axis `hilbert` overload for higher-rank arrays. As of SciPy 1.17,
`hilbert2` excludes the even-length Nyquist frequency in each axis; SciFort
matches this convention, while the 1D `hilbert` retains its Nyquist bin.

The implementation uses the already available SciFort radix-2 FFT and
arbitrary-length direct DFT fallback. Thus non-power-of-two lengths are
correct but can be slower (O(n**2) per transform) than SciPy's optimized
arbitrary-length FFT. Run `test/test_signal_hilbert.f90` with `fpm test`; an
optional 780-case randomized comparison against SciPy 1.17 is available via
`python tools/signal_reference/verify_hilbert.py`. No new runtime dependency
was added.

## Chebyshev-I IIR design, IIR decimation, and Fourier resampling

The public generic `cheby1(order, rp, cutoff, btype, fs)` returns the
existing `butter_result` type, with `b`, `a`, `sos`, and `status` fields.
`rp` is strictly positive ripple in dB; digital critical frequencies must
be strictly within `(0, fs/2)`. Scalar cutoff is accepted for lowpass or
highpass design, and a two-element array for bandpass/bandstop. Prototype
orders 1-20 are supported. SOS output is preferable for high orders because
transfer-function coefficients can be ill-conditioned. Even-order filters
retain Chebyshev-I ripple attenuation at the passband reference frequency.

`decimate(x, q, n, ftype, zero_phase, status [, taps, sos])` now defaults to `ftype='iir'`,
matching SciPy's eighth-order Chebyshev-I antialias filter (`rp=0.05 dB`,
critical normalized frequency `0.8/q`). Set `zero_phase=.false.` for
causal SOS filtering or leave it true to use the existing padded
`sosfiltfilt` implementation. Set `ftype='fir'` to retain the previously
available FIR method. Real and complex signals are supported, and invalid
or too-short zero-phase inputs produce a status code and empty output.
For user-supplied filters, pass `taps=` (real FIR coefficients) or `sos=`
(real, normalized second-order sections) instead of `ftype=`/`n=`. These
trailing keywords work for vector, axis-selected rank-2/3 and caller-allocated
rank-1-through-15 decimators. Complex `sos=` arrays are now supported via
complex-output overloads, including for real input. See
[Custom FIR and SOS decimation](signal_decimate_custom_filters.md) and
[Complex SOS decimation](signal_decimate_complex_sos.md).

`resample` and `resample_into` now also accept `frequency_window(:)` for
user-supplied **real FFT-order spectral weights**, or `complex_window(:)` for
complex weights with complex input. These arguments come after the original
`status` argument and do not change positional calls. Named windows and array
windows are mutually exclusive. The array length must equal the input length
along the resampled axis, and all weights must be finite. See
[Fourier resampling frequency-bin windows](signal_resample_frequency_window.md).

`resample(x, num, window, t, output_times, status)` implements ideal periodic
Fourier resampling by copying low positive and negative FFT bins and handling
the Nyquist coefficient explicitly on even-length up/down sampling. Real and
complex 1D signals are supported; complex inputs additionally accept
`domain='time'` (default) or `domain='freq'` (already transformed). `window`
can be any named periodic window implemented by `get_window`, such as
`'hann'`, `'hamming'`, `'blackman'`, or `'boxcar'`. If `t` has at least two
samples, the optional `output_times` array is filled at the new uniform
sampling interval. The real overload deliberately omits `domain='freq'`
because its Fourier coefficients would generally be complex.

```fortran
use scifort_signal, only : dp, cheby1, butter_result, decimate, resample
real(dp) :: x(100)
real(dp), allocatable :: y(:), tnew(:)
type(butter_result) :: filter
filter = cheby1(4, 0.5_dp, 0.3_dp)
y = decimate(x, 3)                  ! IIR, zero-phase by default
y = decimate(x, 3, ftype='fir')     ! explicit FIR decimator
y = resample(x, 160, window='hann', t=[0.0_dp, 0.01_dp], &
    output_times=tnew)
```

The Fourier implementation shares SciFort's radix-2 FFT and uses a direct
DFT for other lengths, with worst-case O(n**2) time. This can be significantly
slower than SciPy's arbitrary-length optimized FFT backend. Real and complex
**rank-two and rank-three** inputs are now supported via an optional
`axis=` argument. The default is the first axis; positive axis indices are
one-based and negative axes count from the end. The multidimensional complex
overload supports `domain='time'` and `domain='freq'`; the shared
`output_times=` vector describes the selected axis only. Rank-four and higher input require `resample_into`; caller-defined array
spectral windows are now supported for all these ranks. Python-style
frequency-window callable objects and tuple returns are not supported. The Chebyshev and resampling routines
are independent MIT-licensed mathematical implementations, not translated
copies of SciPy or Christopherson's `spectrum` code.

Optional randomized cross-language validation:

```text
python test/compare_signal_cheby_resample_scipy.py
```


## Chebyshev-II and elliptic IIR designs, minimum orders, and resonant biquads

The public `scifort_signal` module additionally exports `cheby2`, `ellip`,
`iirfilter`, `iirdesign`, `buttord`, `cheb1ord`, `cheb2ord`, `ellipord`, `iirnotch`,
`iirpeak`, and `iir_order_result`. They use the same `butter_result` type as
`butter` and `cheby1`, holding `b`, `a`, `sos`, and `status`.

```fortran
use scifort_signal, only : dp, cheby2, ellip, iirfilter, iirdesign, iirnotch, &
    buttord, ellipord, butter_result, iir_order_result

type(butter_result) :: filt, notch
type(iir_order_result) :: ord

filt = cheby2(4, 40.0_dp, [0.2_dp, 0.65_dp], btype='bandpass')
filt = iirfilter(4, 0.3_dp, ftype='cheby2', rs=40.0_dp)
ord = buttord(0.2_dp, 0.34_dp, 1.0_dp, 42.0_dp)
filt = iirdesign(0.2_dp, 0.34_dp, 1.0_dp, 42.0_dp, ftype='cheby2')
filt = ellip(4, 1.0_dp, 45.0_dp, [0.2_dp, 0.6_dp], btype='bandpass')
ord = ellipord(0.2_dp, 0.3_dp, 1.0_dp, 40.0_dp)
notch = iirnotch(60.0_dp, 30.0_dp, fs=200.0_dp)
```

`cheby2(order, rs, cutoff, btype, fs, analog)` accepts prototype orders 1 through
20, positive stopband attenuation up to 100 dB, scalar digital cutoff for
lowpass/highpass, or a two-element increasing cutoff array for
bandpass/bandstop. Frequencies must lie strictly between zero and `fs/2`.
`btype` defaults to `'lowpass'`; `fs` defaults to 2. Its gain normalization
is unity at the selected passband reference. Output SOS section ordering is
not intended to match SciPy's exact pairing, and SOS is preferred over BA
for higher orders. Set `analog=.true.` to use rad/s cutoffs and return
analog transfer-function BA coefficients instead (with `sos` unallocated).

`ellip(order, rp, rs, cutoff, btype, fs, analog)` designs elliptic (Cauer)
filters with both passband ripple `rp > 0` and stopband attenuation
`rs > rp` in dB (currently `rs <= 100`). Orders 1 to 20, low/high/bandpass/
bandstop modes, analog rad/s cutoffs, and digital `fs` are supported. Even
prototype orders have gain `-rp` dB at the passband reference (DC for
lowpass); odd orders have 0 dB. Elliptic transmission zeros are generated
from Jacobi elliptic functions and poles from the elliptic prototype
relations. Digital designs return BA coefficients and second-order sections
(`sos`), using zero/pole bilinear mapping and reference-passband scaling;
analog designs return BA coefficients only, with `sos` unallocated.

`iirfilter(order, cutoff, btype, ftype, rp, rs, fs, analog)` selects
`ftype='butter'` (default), `'cheby1'`, `'cheby2'`, or `'ellip'` (also
`'elliptic'` and `'cauer'`). `rp` is required for Chebyshev I, `rs` for
Chebyshev II, and **both** for elliptic designs.

`buttord(wp, ws, gpass, gstop, fs, analog)`, `cheb1ord(...)`,
`cheb2ord(...)`, and `ellipord(...)` find the minimum prototype order and
critical frequencies for Butterworth, Chebyshev I, Chebyshev II, and elliptic
filters. Arguments `wp` and
`ws` can be **scalars** (lowpass if `wp < ws`, otherwise highpass) or
**rank-one two-element arrays**. For bandpass use
`ws(1) < wp(1) < wp(2) < ws(2)`; for bandstop use
`wp(1) < ws(1) < ws(2) < wp(2)`. These inequalities are enforced.

The optional `analog=.true.` accepts all four shapes with rad/s frequencies,
without a sampling rate. Otherwise `fs` defaults to 2, all frequencies are
in the units of `fs` and must be strictly below Nyquist, `fs/2`.
Analog mode rejects any `fs` argument. Positive `gpass`, `gstop > gpass`,
stopband attenuation at most 100 dB, and selected order at most 20 are
required by the current filter-design implementation.

`iir_order_result` contains `order`, `status`, and two cutoff representations:
`cutoff` is the scalar field preserved for older callers; the new allocatable
`cutoffs(:)` contains **one or two ascending cutoff frequencies**. For a
scalar specification, `cutoff == cutoffs(1)`. For bandpass/bandstop designs,
use the two-element `cutoffs` array. A nonzero status indicates an invalid
or unsupported specification; on failure `cutoffs` remains unallocated.

```fortran
use scifort_signal, only : dp, iir_order_result, butter_result, &
    buttord, cheb2ord, butter, iirdesign

type(iir_order_result) :: order_spec
type(butter_result) :: filter

order_spec = buttord([20.0_dp,50.0_dp], [14.0_dp,60.0_dp], &
    3.0_dp,40.0_dp,analog=.true.)
! order_spec%order == 13; order_spec%cutoffs has two rad/s frequencies
filter = butter(order_spec%order,order_spec%cutoffs, &
    btype='bandpass',analog=.true.)

order_spec = cheb2ord([0.1_dp,0.6_dp], [0.2_dp,0.5_dp], &
    3.0_dp,60.0_dp)
filter = iirdesign([0.1_dp,0.6_dp],[0.2_dp,0.5_dp], &
    3.0_dp,60.0_dp,ftype='cheby2')
```

`iirdesign(wp,ws,gpass,gstop,ftype,fs,analog)` selects the minimum-order
filter and calls the matching `iirfilter` routine, including bandpass,
bandstop, and analog designs. SciFort retains `ftype='butter'` as the
default for backward compatibility; specify `ftype='ellip'` to obtain the
elliptic design (SciPy's `iirdesign` defaults to elliptic).
`iirdesign` returns `butter_result` BA and digital SOS data;
for analog designs only BA is available.

Bandstop order selection minimizes the continuous prototype order by bounded
one-dimensional searches on the two passband boundaries. SciPy uses a
different stopping tolerance for those searches, so bandstop critical
frequencies can differ slightly while meeting the same specifications.
High-order BA coefficients, particularly for analog and narrowband designs,
are ill-conditioned even if the underlying filter is stable. Prefer digital
SOS when available.

Fixed regressions: `test/test_signal_iir_order_bands.f90` and
`test/test_signal_elliptic.f90`. For independent randomized verification
against SciPy 1.17 of the older families and the new elliptic family, run:

```text
python tools/signal_reference/verify_iir_order_bands.py
python tools/signal_reference/verify_elliptic.py
```

The elliptic verifier covers 1,200 analog/digital filter designs, 800
minimum-order selections, and 231 `iirdesign` cases across four filter
shapes. SciFort's bandstop edge optimization has a different stopping
tolerance from SciPy, so individual optimized bandstop cutoffs and BA
coefficients can differ slightly. SOS coefficients are not guaranteed to
have the same section order or per-section gain distribution as SciPy;
compare the total transfer response. For high orders, SOS avoids the
polynomial-conditioning problems of BA output. Analog ZPK and SOS output are available through the separate APIs below.

### Analog zero-pole-gain and second-order-section output

The independent analog-design entry points `iirfilter_zpk` and `iirfilter_sos`
return `analog_zpk_result` (`zeros(:)`, `poles(:)`, `gain`, `status`) and
`analog_sos_result` (`sos(:,:)`, `status`). They support `butter`, `cheby1`,
`cheby2`, and `ellip` (`elliptic`/`cauer`), lowpass, highpass, bandpass, and
bandstop. Scalar `cutoff` and one-/two-edge array overloads are available.
Analog cutoff frequencies are in rad/s, with no `fs` argument.

```fortran
use scifort_signal, only : dp, analog_zpk_result, analog_sos_result, &
    frequency_response, iirfilter_zpk, iirfilter_sos, sosfreqs, freqs_zpk

type(analog_zpk_result) :: zp
type(analog_sos_result) :: sec
type(frequency_response) :: response

zp = iirfilter_zpk(5, [2.0_dp, 5.0_dp], btype='bandpass', &
    ftype='ellip', rp=1.0_dp, rs=45.0_dp)
sec = iirfilter_sos(5, [2.0_dp, 5.0_dp], btype='bandpass', &
    ftype='ellip', rp=1.0_dp, rs=45.0_dp)
if (sec%status == 0) then
    response = sosfreqs(sec%sos, [1.0_dp, 2.0_dp, 3.0_dp, 5.0_dp])
end if
```

`zpk2sos_analog(zeros,poles,gain)` converts arbitrary real-coefficient,
proper analog transfer functions from their finite roots, including odd-order
systems and constant gains. The roots must occur in real or complex-conjugate
pairs. Section numerator/denominator coefficients are **descending powers of
`s`** in `[b2,b1,b0,a2,a1,a0]` layout. Missing highest-degree terms are
explicitly zero. Unlike digital SOS arrays these may have a denominator's
first coefficient zero for first-order sections; **do not use** `sosfilt` or
`sosfreqz` to process/evaluate analog SOS. Use `sosfreqs(sos,w)` with explicit
angular frequencies in rad/s, or `freqs_zpk(z,p,k,w=w)` for ZPK responses.

The original `butter`, `cheby1`, `cheby2`, `ellip`, and `iirfilter` APIs
remain backward compatible and still return analog **BA coefficients only**
when `analog=.true.`. A separate design API avoids changing the meaning of
`butter_result%sos`, whose existing values are digital sections.
For high analog orders, prefer ZPK or SOS representation over polynomial BA.
For Butterworth and Chebyshev designs, these APIs avoid forming BA
polynomials internally, so they can return finite root/section data for
large frequencies and orders even when expanded BA coefficients overflow.
Output section order and per-section scaling need not match SciPy's exactly;
the full frequency response and transfer-function roots are validated instead.

Deterministic tests are in `test/test_signal_analog_sos.f90`. A randomized
SciPy 1.17 cross-language verifier with 1,280 designs and 12,800 response
points is available via `python test/compare_signal_analog_sos_scipy.py`
(NumPy, SciPy, and gfortran required; FPM not required).

`iirnotch(w0, qfactor, fs)` and `iirpeak(w0, qfactor, fs)` implement real
second-order notch/rejection and resonant/peak filters. Both return BA and a
single SOS section; `w0` is the center frequency and `qfactor` is center
frequency divided by the -3 dB bandwidth. For numerical safety this API
requires `0 < w0 < fs/2`, positive Q, and normalized bandwidth below pi.

Invalid inputs leave the result status at `design_status_invalid_input`
rather than raising an exception. The original MIT source contains no new
FPM dependencies. An optional cross-language comparison script runs from the
project root with `python tools/signal_reference/verify_iir_extra.py` and
checks 1,350 randomized cases against SciPy 1.17.0; deterministic cases are
in `test/test_signal_iir_extra.f90`.

## Filter histories and group delay reference validation

The new routines `lfiltic` and `group_delay` have standalone regression tests:
`test/test_signal_lfiltic.f90` and `test/test_signal_group_delay.f90`.
They check nonunit `a(1)`, real/complex filters, true continuation from
prior samples, sampled frequency grids, sample-rate units, singularity
status, and nonfinite input validation. An optional 600-case randomized
comparison against SciPy 1.17 is available via
`python tools/signal_reference/verify_filter_diagnostics.py`.

## Frequency sweeps and Gaussian-modulated pulses

`scifort_signal` now exports the SciPy-style real-valued waveform generators
`chirp(t, f0, t1, f1, method, phi, vertex_zero, status)`,
`sweep_poly(t, poly, phi, status)`, and
`gausspulse(t, fc, bw, bwr, quadrature, envelope, status)`.
All take **real(dp) rank-one arrays of sample times** and return allocated
real(dp) arrays; the optional `status` is an integer `intent(out)` argument.

```fortran
use scifort_signal, only : dp, chirp, chirp_complex, sweep_poly, &
    gausspulse, gausspulse_cutoff
real(dp), allocatable :: t(:), y(:), q(:), env(:)
complex(dp), allocatable :: analytic(:)
real(dp) :: cutoff
integer :: i

t = [(0.01_dp*real(i,dp),i=0,200)]
y = chirp(t, f0=3.0_dp, t1=2.0_dp, f1=12.0_dp)
y = chirp(t, 3.0_dp, 2.0_dp, 12.0_dp, method='quadratic', &
    phi=90.0_dp, vertex_zero=.false.)
analytic = chirp_complex(t, 3.0_dp, 2.0_dp, 12.0_dp, method='logarithmic')
y = sweep_poly(t, [0.5_dp,-1.0_dp,3.0_dp])  ! f(t) = 0.5*t**2 - t + 3
y = gausspulse(t-1.0_dp, fc=8.0_dp, quadrature=q, envelope=env)
cutoff = gausspulse_cutoff(fc=8.0_dp, tpr=-60.0_dp)
```

`chirp` follows SciPy's four instantaneous-frequency laws:
`'linear'` (`'lin'`, `'li'`), `'quadratic'` (`'quad'`, `'q'`),
`'logarithmic'` (`'log'`, `'lo'`), and `'hyperbolic'` (`'hyp'`).
The default method is linear. `phi` is in degrees, whereas `f0` and `f1`
are in **cycles per unit of time**, not radians per unit. A quadratic chirp
uses a vertex at time zero unless `vertex_zero=.false.`, which puts the vertex
at `t1`. Logarithmic chirps require two nonzero frequencies of the same sign;
hyperbolic chirps require both frequencies to be nonzero and cannot evaluate at
their singular time. This API requires `t1` to be nonzero for every method.

SciPy's `chirp(..., complex=True)` has a separate Fortran entry point,
`chirp_complex`, because a Fortran function cannot choose its result type from
a runtime logical flag. It returns the analytic complex exponential with
unit magnitude, not the cosine output of `chirp`. Apart from the result type,
its arguments have the same names and meanings.

`sweep_poly` accepts polynomial frequency coefficients in **descending powers**
(as `numpy.polyval` does), integrates the polynomial with zero constant of
integration, then evaluates its cosine. Empty `poly` means identically zero
frequency; negative frequencies and arbitrary polynomial degree are allowed.

`gausspulse` computes a Gaussian-envelope cosine. `fc` defaults to 1000,
`bw` to 0.5, and `bwr` to -6 dB. Instead of SciPy's tuple-based `retquad`
and `retenv` flags, request optional `quadrature=` and/or `envelope=` output
arrays. `gausspulse_cutoff(fc, bw, bwr, tpr, status)` corresponds to
`scipy.signal.gausspulse('cutoff', ...)`; `tpr` defaults to -60 dB. All
specified frequencies and bandwidths must be finite, with `fc >= 0`,
`bw > 0`, and `bwr < 0`; cutoff also requires `tpr < 0` and `fc > 0`.
Zero center frequency is allowed for the sampled pulse but gives no finite
cutoff, so `gausspulse_cutoff` reports invalid input rather than returning
infinity.

The status constants are `waveform_status_success=0` and
`waveform_status_invalid_input=1`. Invalid inputs return quiet NaN arrays
(or a NaN scalar for cutoff) and an invalid status; valid empty time arrays
produce valid empty results. This initial API is limited to rank-one arrays
(no multi-axis broadcasting or scalar time overloads), and it does not support
SciPy's polynomial object argument to `sweep_poly` directly; pass its
coefficients as an array.

`test/test_signal_waveforms.f90` checks all four chirp methods, complex chirps,
quadratic vertex placement, polynomial sweeps, pulse in-phase/quadrature and
envelope components, pulse cutoff, and invalid inputs. The additional script

```text
python tools/signal_reference/verify_waveforms.py
```

compiles its own comparison driver with GNU Fortran and checks **1,210**
randomized cases against installed SciPy (tested with SciPy 1.17.0).
The largest observed absolute difference was approximately `1.8e-12`.
All routines are independent MIT-licensed Fortran implementations with no
additional FPM dependencies.


## Periodic waves, impulses, and maximum-length sequences

The independent waveform generators `sawtooth`, `square`, `unit_impulse`,
and `max_len_seq` are exported through `scifort_signal`.

```fortran
use scifort_signal, only : dp, sawtooth, square, unit_impulse, &
    max_len_seq, max_len_seq_result
real(dp), allocatable :: t(:), y(:), u(:), image(:,:)
type(max_len_seq_result) :: bits, continuation
integer :: i

t = [(0.1_dp*real(i,dp), i=0,199)]
y = sawtooth(t, width=0.5_dp)     ! triangular wave, period 2*pi
u = square(t, duty=0.25_dp)      ! 25-percent positive phase
u = square(t, duty=(sin(t)+1.0_dp)/2.0_dp) ! variable-width pulses
u = unit_impulse(9, idx=3)       ! one at Fortran index 3
image = unit_impulse([3,5], centered=.true.) ! one at (2,3)
bits = max_len_seq(8, length=100)
continuation = max_len_seq(8, state=bits%state, length=100)
```

### Periodic `sawtooth` and `square`

`*wave*(t, width_or_duty, status)` accepts `real(dp)` rank-one phase arrays
in **radians** (the period is `2*pi`). The fraction can be a scalar or a
rank-one array **of the same length as `t`**. For `sawtooth`, the default
`width=1` is a rising ramp from -1 to +1; `width=0` is a falling ramp,
and `width=0.5` is triangular. `square` defaults to `duty=0.5`, returning
+1 before the switching point and -1 afterwards. Both accept fractions
in the closed interval `[0,1]`, including the constant-wave square limits.
These functions are not band-limited and can produce aliasing, as in SciPy.

Bad fractions, nonfinite input phases, and array-length mismatches produce an
all-NaN return and `periodic_status_invalid_input=1` (success is 0). Unlike
SciPy's elementwise NaNs for a partly invalid fraction array, **SciFort
rejects the entire array**. There is currently no NumPy-style multidimensional
broadcasting. Empty arrays are valid.

### `unit_impulse`: one-based Fortran indices

`unit_impulse(n, idx, centered, status)` returns a one-dimensional `real(dp)`
array. `unit_impulse(shape, idx, centered, status)` returns a rank-two array
when `shape` is an integer pair `[rows, columns]`; its `idx` is also an
integer pair. Indices are **one-based** in the Fortran API, unlike SciPy's
zero-based `idx`. When `idx` is absent, the impulse is at index 1 in each
dimension, unless `centered=.true.` selects `floor(shape/2)+1`. Explicit
`idx` takes priority over `centered`. The zero-sized output is allowed when a
dimension is zero and no explicit index is supplied; negative dimensions,
out-of-range indices, and indices in an empty array are rejected. Invalid
nonempty positions return quiet NaNs plus nonzero `periodic_status`.

Only `real(dp)` rank-one and rank-two impulses are implemented. SciPy accepts
arbitrary dimensions and its `dtype` parameter is not currently supported.

### `max_len_seq`: linear-feedback shift registers

`max_len_seq(nbits, state, length, taps, skip)` returns `max_len_seq_result`,
with integer arrays `seq` (0/1 output bits) and `state` (final register)
and an integer `status`. `nbits` must be in `[2,32]`. The default initial
register contains all ones; a supplied integer `state` has length `nbits`,
with nonzero entries treated as ones. An all-zero initial register is
invalid. Pass `bits%state` into the next invocation to generate the sequence
in blocks without losing samples.

Default primitive-polynomial feedback taps are provided for all 2 to 32 bits,
following the table linked from the SciPy `max_len_seq` documentation. Custom
`taps` are the same **polynomial tap positions** used by SciPy: positive
integer positions in `[1,nbits-1]` (not ordinary one-based array indices).
Each custom tap location is counted once if duplicated. These taps need not
be primitive, so a custom choice might not produce a maximum-length cycle.
Unlike SciPy's ability to accept tap endpoint values `0` or `nbits`, this
Fortran interface rejects them to avoid degenerate conventions.

An optional `skip=` advances or rewinds the Fibonacci register **before**
producing any output bits. The signed argument has kind `integer(int64)`;
positive values move forward, negative values move backward, and zero leaves
the initial state alone. For example:

```fortran
use, intrinsic :: iso_fortran_env, only : int64
use scifort_signal, only : max_len_seq, max_len_seq_result
type(max_len_seq_result) :: segment, rewound

segment = max_len_seq(16, length=128, skip=1000000000000_int64)
rewound = max_len_seq(16, state=segment%state, length=0, &
                       skip=-1000000000128_int64)
! rewound%state is the original all-ones register
```

This **SciFort extension** is not a keyword of SciPy's `signal.max_len_seq`.
Repeated squaring of the binary state-transition matrix takes
`O(nbits**2 * log(abs(skip)))` bit operations and `O(nbits)` packed-word
workspace; it does **not** produce or allocate skipped samples. The feedback
recurrence has an invertible one-clock transition because the first register
bit always participates in feedback; rewind therefore works even if a custom
tap set is not primitive. All signed 64-bit distances, including the most
negative `int64`, are handled without overflowing `abs(skip)`. A call with
`length=0` returns just the advanced/rewound state. The final register state
always corresponds to `skip + length` clock steps relative to the supplied
initial state. The sequence output and saved-state semantics are otherwise
identical to the earlier API.

The default `length` is `2**nbits-1` for `nbits<=20`. For `nbits>20`,
**explicit `length=` is required** to avoid implicitly allocating enormous
full-cycle arrays. `length=0` succeeds without advancing the state. Negative
lengths, invalid state dimensions, invalid taps, and invalid bit widths return
`sequence_status_invalid_input=1`; allocation errors return
`sequence_status_allocation_failed=2`. A zero `status` denotes success.

The binary recurrence implements an ordinary Fibonacci linear-feedback shift
register. The primitive-tap table follows the mathematical reference linked
from [SciPy max_len_seq](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.max_len_seq.html).
No SciPy source code is imported, and no additional FPM dependency is needed.

The deterministic tests are `test/test_signal_periodic.f90`,
`test/test_signal_sequences_skip.f90`, and `test/test_signal_sequences.f90` (including full-period balance,
two-level circular autocorrelation, resumability, default tap tables,
and error statuses). A reproducible reference check is available with:

```text
python tools/signal_reference/verify_periodic_sequences.py
```

Using GNU Fortran 14.2.0 and SciPy 1.17.0, all 1,300 randomized SciPy
comparisons passed. Real waveform values agreed to within `1.62e-14`
absolute error; MLS samples and final register states agreed exactly.
The full 22-program standalone signal suite passed with bounds checking.
For skip-ahead, `python tools/run_signal_sequence_skip.py` runs the previous
native tests, new jump/rewind and maximal-`int64` tests, and 698 independent
SciPy 1.17 comparisons in both checked and optimized GNU Fortran builds.
The reference constructs the one-step transition from SciPy output, inverts
it by bitwise Gauss-Jordan elimination, and exponentiates it independently.
All emitted bits and final states agree exactly.


## Spectra and phase synchrony from irregularly spaced samples

`scifort_signal` now exports independent `lombscargle`,
`lombscargle_amplitude`, and `vectorstrength` routines.

```fortran
use scifort_signal, only : dp, lombscargle, lombscargle_amplitude, &
    vectorstrength, vector_strength_result
real(dp) :: times(5), observations(5), omega(3)
real(dp), allocatable :: power(:), fractional_power(:)
complex(dp), allocatable :: amplitude(:)
type(vector_strength_result) :: synchrony

times = [0.0_dp,0.3_dp,0.9_dp,1.5_dp,2.4_dp]
observations = [0.0_dp,0.5_dp,-0.3_dp,0.6_dp,-0.1_dp]
omega = [1.0_dp,2.0_dp,3.0_dp]  ! radians per time unit
power = lombscargle(times,observations,omega,floating_mean=.true.)
fractional_power = lombscargle(times,observations,omega, &
    normalize='normalize',floating_mean=.true.)
amplitude = lombscargle_amplitude(times,observations,omega, &
    floating_mean=.true.)
synchrony = vectorstrength(times, [1.0_dp,2.0_dp])
```

`lombscargle(x, y, freqs, precenter, normalize, weights, floating_mean,
status)` accepts **rank-one real(dp)** observation times and measurements;
`freqs` are angular frequencies in **radians per unit time** (not cycles).
It fits independent weighted sinusoidal least-squares models
`a*cos(omega*x) + b*sin(omega*x) + c`, with `c=0` by default or estimated
per frequency if `floating_mean=.true.`. Weights are nonnegative and
normalized to sum one, so rescaling all weights does not change the result.
The default `normalize='power'` matches SciPy's convention of
`N/4 * A**2` for a well-sampled sinusoid of amplitude `A`; use
`normalize='normalize'` for the fraction of weighted reference-model
variance explained by the fitted sinusoid. The optional legacy
`precenter=.true.` subtracts the **unweighted** arithmetic mean of `y`
first. For an unknown baseline, prefer `floating_mean` to `precenter`.

`lombscargle_amplitude` shares `x`, `y`, `freqs`, `precenter`, `weights`,
`floating_mean` and `status`, but returns `complex(dp)` best-fit amplitude
coefficients. This is a separate entry point because Fortran cannot
select a function's real or complex return type using a runtime string.
A fitted sinusoid can be reconstructed as
`real(amplitude)*cos(omega*x) + aimag(amplitude)*sin(omega*x)` (plus the
fitted intercept if `floating_mean` is enabled). SciPy's
`normalize='amplitude'` corresponds to this separate function.

`vectorstrength(events, periods)` returns a `vector_strength_result`
with `strength(:)`, `phase(:)` and integer `status`. It computes the
magnitude and argument of the mean complex event phase at each test
period. `periods` may be a rank-one array or a scalar; in both cases
**the result fields are rank-one** (length one for a scalar input).
Input periods must be finite and strictly positive.

Invalid lengths, nonfinite observations or frequencies, negative weights,
nonpositive weight totals, and unsupported normalization modes are
rejected with `irregular_status_invalid_input=1` and NaN results. Empty
sample arrays or frequency/period arrays are invalid, as in SciPy.
For a constant reference signal, normalized Lomb power is undefined:
`irregular_status_singular=2` is returned and the normalized result is
NaN; ordinary power/amplitude still compute. Success is
`irregular_status_success=0`. As with SciPy, amplitude estimation near
zero frequency or for a severely degenerate sampling grid can be
ill-conditioned. This is a direct O(n*m) frequency-sweep implementation,
not a fast grid/FFT Lomb-Scargle algorithm.

### Weighted multiharmonic Lomb–Scargle (SciFort extension)

`lombscargle(x, y, freqs, nterms=H)` and
`lombscargle_multiharmonic(x, y, freqs, H)` jointly fit `H` sinusoidal
harmonics. Weighted floating means, `precenter`, power/normalized outputs,
and singular-frequency diagnostics are supported. The fit uses
reorthogonalized QR rather than the normal equations. With `nterms` absent or
one, `lombscargle` retains its original SciPy-compatible implementation.
This `nterms` keyword is specific to SciFort; SciPy 1.17 does not implement
it. See [weighted multiharmonic periodograms](signal_multiharmonic.md) for
mathematical conventions and independent numerical verification.
`lombscargle_harmonic_fit` returns the fitted sine/cosine coefficients and
intercepts, and optionally (`estimate_covariance=.true.`) their estimated
covariance and standard errors. The companion
`lombscargle_harmonic_predict_std` calculates fitted-mean standard errors,
conditional on fixed trial frequencies. See
[harmonic covariance and mean-curve uncertainty](signal_harmonic_uncertainty.md).

### Multichannel Lomb-Scargle through rank 15 (SciFort extension)

`lombscargle_into` and `lombscargle_amplitude_into` accept real observations
of rank 1 through 15 and one chosen sample axis. Every other axis indexes
independent channels. The times `x(:)` and angular frequencies `freqs(:)` are
common to all channels; this multichannel API is a **SciFort extension**
(SciPy's `lombscargle` accepts only rank-one observations).

```fortran
use scifort_signal, only : dp, lombscargle_into, lombscargle_amplitude_into
real(dp) :: times(120), observations(4,120,3), omega(81)
real(dp) :: shared_weights(120), channel_weights(4,120,3)
real(dp) :: power(4,81,3)
complex(dp) :: amplitude(4,81,3)
integer :: status

! Common weights: reuse the trigonometric frequency basis for all channels.
call lombscargle_into(times, observations, omega, power, axis=2, &
    weights=shared_weights, floating_mean=.true., &
    normalize='normalize', status=status)

! Different sampling confidence weights for each channel.
call lombscargle_amplitude_into(times, observations, omega, amplitude, &
    axis=-2, weights_nd=channel_weights, floating_mean=.true., &
    status=status)

! Omit each channel's nonfinite observations and zero-weight samples.
call lombscargle_into(times, observations, omega, power, axis=2, &
    weights_nd=channel_weights, nan_policy='omit', &
    normalize='normalize', floating_mean=.true., status=status)
```

**Weights and missing observations:** `weights(:)` is an optional shared
nonnegative weight vector of length `size(x)`; `weights_nd` is an optional
nonnegative, finite real array with the **same rank and shape as the input**,
including its sample axis. Use either, not both. With the default
`nan_policy='raise'`, nonfinite observations are invalid, including where
weights are zero, preserving the original SciFort behavior. With
`nan_policy='omit'`, every channel independently excludes observations whose
`y` is NaN or infinite and observations with zero weight. Nonfinite times,
frequencies, or weights always remain invalid. All-positive-weight samples
are retained. An empty channel or a later singular weighted fit rejects the
**entire** call. Omitting samples changes the sample count used by SciPy's
unnormalized Lomb-Scargle power; this is the behavior obtained by masking each
channel and then calling SciPy's rank-one function. `precenter` subtracts the
unweighted mean of the **retained** observations when omission is enabled.
This is a SciFort convenience option rather than SciPy's own `nan_policy`.

**Shapes and error handling:** Replace the sample-axis extent of the input
with `size(freqs)` to obtain the exact output shape, keeping every other
extent and the rank unchanged. `axis` uses one-based indices or negative
indices from the end; its default is the last dimension. Caller output must
be allocated. Noncontiguous array sections are supported for input, weights,
and output. All channel results are staged before assignment: invalid
parameters, missing data when `nan_policy='raise'`, or a singular channel
leave **every element of the provided output unchanged**. Errors are reported
via `status=`, consistent with the rank-one kernel.

**Numerical method:** With shared weights (or no weights) and the default
`nan_policy='raise'`, the implementation evaluates each sine and cosine at
most once per sample/frequency combination, reusing weighted Gram terms
across channels. The weighted least-squares projection is evaluated separately
for every channel, preserving the original normalization, precentering,
floating-mean, and amplitude conventions. With `weights_nd` or sample
omission, the existing rank-one weighted solver is applied to each channel,
including its independent list of retained samples. Computational complexity
remains `O(n_samples*n_freqs*n_channels)`; shared trigonometry reduces the
constant factor, not the asymptotic order. This is **not** an FFT-based
Lomb-Scargle transform. Additional cached workspace is
`O(n_samples*n_channels)` in the shared-weight path.

**Reusable plans (SciFort extension):** For repeated spectra on the same
observation times, angular-frequency grid, and common sample weights, a
`lombscargle_plan` stores the weighted phase-rotated sine/cosine basis and
least-squares denominators. Plan construction costs
`O(n_samples*n_freqs)` time and storage (approximately 16 bytes per
sample/frequency pair for the two double-precision basis arrays). Each
subsequent evaluation still costs `O(n_samples*n_freqs*n_channels)` but
avoids reevaluating trigonometric functions and Gram sums; speedups depend
on the number of channels, number of evaluations, and array layout.
The `floating_mean` setting is fixed at initialization, whereas `precenter`
and `normalize` are specified for each evaluation.

```fortran
use scifort_signal, only : dp, lombscargle_plan, &
    lombscargle_plan_into, lombscargle_plan_amplitude_into
real(dp) :: times(120), omega(81), shared_weights(120)
real(dp) :: first_record(4,120,3), next_record(4,120,3)
real(dp) :: power(4,81,3)
complex(dp) :: amplitude(4,81,3)
type(lombscargle_plan) :: plan
integer :: status

! Populate the input arrays, times, weights and frequencies first.
call plan%initialize(times, omega, weights=shared_weights, &
    floating_mean=.true., status=status)
call lombscargle_plan_into(plan, first_record, power, axis=2, &
    normalize='normalize', status=status)
call lombscargle_plan_amplitude_into(plan, next_record, amplitude, &
    axis=-2, status=status)
call plan%clear()
```

`initialize` validates all inputs and preserves an already initialized
plan if a later initialization fails. `is_initialized()` reports its
state. The prepared plan is read-only during evaluation, so concurrent
calls may share it, provided each call has its own input and output
arrays. Outputs are staged until every channel succeeds; errors preserve
all caller data. A plan is valid only for its initialized time/frequency
and shared-weight grids. Per-channel `weights_nd` and
`nan_policy='omit'` are **not supported by these planned interfaces**:
they may change the effective least-squares design for each channel.
Use the existing unplanned adapters for those cases. Precomputation also
requires more memory than the on-the-fly shared-weight path; avoid plans
when the product `n_samples*n_freqs` is too large or when only one
spectrum is needed.

Independent numerical comparisons and regression tests:

```text
test/test_signal_irregular_plan.f90
python test/compare_signal_irregular_plan_scipy.py checked
python test/compare_signal_irregular_plan_scipy.py optimized
benchmarks/bench_signal_lombscargle_plan.f90
test/test_signal_irregular_nd.f90
test/test_signal_irregular_nd_weights.f90
python test/compare_signal_irregular_nd_scipy.py checked
python test/compare_signal_irregular_nd_scipy.py optimized
python test/compare_signal_irregular_nd_weights_scipy.py checked
python test/compare_signal_irregular_nd_weights_scipy.py optimized
benchmarks/bench_signal_lombscargle_shared_trig.f90
```

The standalone Fortran regression is `test/test_signal_irregular.f90`.
An additional optional SciPy reference check is:

```text
python tools/signal_reference/verify_irregular.py
```

With GNU Fortran 14.2.0 and SciPy 1.17.0, **900 randomized cases**
(11,201 frequency bins) passed, covering weights, precentering,
floating means, normalized power, complex amplitudes, and vector
strength/phase. Maximum observed absolute error was `5.55e-13`.
No new FPM dependencies or third-party implementation code were added.


## Chirp Z-transform and zoom FFT

```fortran
use scifort_signal, only : dp, czt, zoom_fft
complex(dp), allocatable :: z(:), focused(:)
real(dp) :: x(4)
x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
z = czt(x) ! four output values identical to the ordinary DFT
focused = zoom_fft(x, [0.2_dp, 0.4_dp], m=32, fs=2.0_dp)
```

`czt(x, m, w, a, status)` computes the one-dimensional chirp Z-transform:

`y(k+1) = sum(x(n+1)*a**(-n)*w**(n*k), n=0..size(x)-1)`

for zero-based frequency index `k` from `0` to `m-1`. Inputs can be
`real(dp)` or `complex(dp)`; the result is always `complex(dp)`. `m` defaults
to `size(x)`, `w` defaults to `exp(-2*pi*i/m)`, and `a` defaults to `1`.
Provide `w=` and `a=` as complex scalars when choosing another logarithmic
spiral. The one-shot API uses ordinary Fortran rank-one arrays and explicit
status integers. Reusable derived-type plans are documented below.

`zoom_fft(x, fn, m, fs, endpoint, status)` returns samples of the
Z-transform at uniformly spaced frequencies in cycles per unit time.
`fn` may be a scalar (from zero to `fn`) or a two-element real array
`[f_start, f_stop]`; `m` defaults to the signal length, `fs` defaults to
`2`, and `endpoint=.false.` excludes `f_stop`, matching SciPy's conventions.
With `endpoint=.true.`, `f_stop` is included when `m > 1`.

For unit-circle contours these use Bluestein convolution with SciFort's
existing radix-2 FFT, costing roughly O((n+m) log(n+m)) operations. For
contours with `abs(w)` or `abs(a)` appreciably different from one, the
implementation uses direct Horner evaluation (O(n*m)), preventing the
large cancelling intermediate chirp factors that can arise in the FFT
formula. Extremely long or very rapidly growing/decaying nonunit spirals
may still yield nonfinite outputs. Status values are
`czt_status_success=0`, `czt_status_invalid_input=1`, and
`czt_status_nonfinite_result=2`; invalid calls return an empty array.

The Fortran test is `test/test_signal_czt.f90`. To compare against SciPy:

```text
python tools/signal_reference/verify_czt.py
```

On GNU Fortran 14.2 and SciPy 1.17, 480 deterministic randomized cases
passed for standard and custom CZTs, off-unit spirals, and scalar/pair
zoom frequency intervals, with maximum scaled error `2.055e-13`.

### Reusable CZT and zoom-FFT plans

For repeated transforms with the same input length, output length, and spiral,
use `czt_plan` or `zoom_fft_plan`. Initialization caches the Bluestein chirps
and Fourier-transformed convolution kernel; later calls reuse them.
Initialization and transforms report the same status codes as the original
functions. A plan can be reinitialized; invalid initialization makes it
unusable until a successful reinitialization.

```fortran
use scifort_signal, only : dp, czt_plan, zoom_fft_plan
implicit none
type(czt_plan) :: chirp
type(zoom_fft_plan) :: zoom
complex(dp), allocatable :: spectrum(:), channels_out(:,:)
real(dp) :: observations(128), channels(128,4)
integer :: status

! Populate observations and channels before calling transform.
call chirp%initialize(128, m=256, status=status)
if (status == 0) spectrum = chirp%transform(observations, status=status)
if (status == 0) channels_out = chirp%transform(channels, axis=1, status=status)

call zoom%initialize(128, [0.1_dp, 0.4_dp], m=256, fs=2.0_dp, status=status)
if (status == 0) spectrum = zoom%transform(observations, status=status)
```

`czt_plan%initialize(n, m, w, a, status)` takes the same spiral options as
`czt`, except it receives the input length `n` instead of the data. The
`zoom_fft_plan%initialize(n, fn, m, fs, endpoint, status)` generic accepts a
scalar or a two-element frequency interval, just like `zoom_fft`.
Both types expose `transform` for real and complex vectors and rank-two arrays,
with an optional **one-based** `axis=1` or `axis=2` for matrices (default 1).
The transformed axis keeps its position, with its length replaced by `m`.
`czt_plan%points()` returns the complex evaluation contour, and
`zoom_fft_plan%frequencies()` returns its frequency coordinates.
`is_valid()`, `input_length()`, and `output_length()` report plan metadata.
Input-length mismatches, invalid axes, and nonfinite samples return empty
output and a nonzero status.

The FFT-based precomputation is used for near-unit-magnitude `w` and `a`,
which cover ordinary DFT and zoom-FFT transforms. Other contours use cached
spiral ratios with stable Horner evaluation; these remain O(n*m) per call,
so off-unit contours do not gain the FFT complexity reduction.

The new native test is `test/test_signal_czt_plans.f90`, and comparison with
SciPy's reusable `CZT` and `ZoomFFT` classes can be run using:

```text
python tools/signal_reference/verify_czt_plans.py
```

This script exercises repeated transforms, complex/real inputs, unit and
nonunit spirals, zoom intervals, and both matrix axes.

### Rank-independent CZT and zoom FFT (ranks 1–15)

`czt_into(x, y, m, w, a, axis, status)` and
`zoom_fft_into(x, y, fn, m, fs, endpoint, axis, status)` expose the
SciPy-style arbitrary-axis behavior for caller-allocated arrays of rank
1 through 15. Real or complex input yields **complex output**. The
transform axis is **one-based** (negative axes count from the end) and
defaults to the **last** axis, matching SciPy's default axis in Fortran
numbering. The caller must allocate `y` with the same rank and all the
same non-transform dimensions; the transform-axis length must equal `m`
(or the input-axis length if `m` is omitted).

```fortran
use scifort_signal, only : dp, czt_into, zoom_fft_into, &
    czt_plan, czt_plan_into, zoom_fft_plan, zoom_fft_plan_into
implicit none
real(dp) :: samples(4,128,3)
complex(dp) :: spiral(4,256,3), focused(4,64,3)
type(czt_plan) :: chirp
type(zoom_fft_plan) :: zoom
integer :: st

! Initialize samples first. The second axis is the signal-time axis.
call czt_into(samples,spiral,m=256,axis=2,status=st)
call zoom_fft_into(samples,focused,[0.15_dp,0.35_dp], &
    m=64,fs=2.0_dp,axis=2,status=st)

! Reuse precomputed chirp and FFT kernels over successive datasets:
call chirp%initialize(128,m=256,status=st)
if (st == 0) call czt_plan_into(chirp,samples,spiral,axis=2,status=st)
call zoom%initialize(128,[0.15_dp,0.35_dp],m=64,fs=2.0_dp,status=st)
if (st == 0) call zoom_fft_plan_into(zoom,samples,focused,axis=2,status=st)
```

`zoom_fft_into` accepts scalar and pair `fn` forms. Both one-shot
adapters construct **one** transform plan and reuse it for all signal
lanes. The plan-based adapters accept previously initialized plans,
so repeated batches need not recompute chirps and transform kernels.
The output is staged in an internal flat array and written only when all
lanes succeed. Invalid shapes, invalid axes, invalid plans, and
nonfinite input/outputs return a nonzero status **without changing `y`**.
Noncontiguous input/output array sections and zero-sized *channel*
dimensions are supported; the transform-axis length must be positive.
The whole-batch staging design requires O(size(x)+size(y)) additional
memory. Avoid passing overlapping input/output sections as separate
Fortran actual arguments.

The new standalone regression is `test/test_signal_czt_nd.f90`.
Independent comparisons with SciPy 1.17 are available through:

```text
python test/compare_signal_czt_nd_scipy.py checked
python test/compare_signal_czt_nd_scipy.py optimized
```

Each mode passes 466 SciPy CZT/zoom comparisons (36,060 coefficients)
with maximum magnitude-scaled difference 2.19e-12. The SciPy 1.17
`zoom_fft` implementation raises a `ZeroDivisionError` when
`m=1` and `endpoint=True`; those cases are omitted from the independent
SciPy comparison rather than incorrectly reported as matched.

### Remez portability

`remez` prefers a >=30-digit real kind for internal exchange arithmetic,
falling back to a >=18-digit kind (if provided by the compiler), then to
`dp`. Its coefficients and status API remain `real(dp)` and integer.
The 80-bit path reduces the number of nonconverged cases in the supplied
240-design stress test from 21 to four compared with pure double precision;
full binary128 passed all 240. Compilers exposing only `dp` may return
`remez_status_no_convergence` for difficult high-order designs.


## Partial-fraction decomposition of analog and digital filters

`residue`, `residuez`, `invres`, and `invresz` follow SciPy's
partial-fraction convention with explicit Fortran result types:

```fortran
use scifort_signal, only : dp, residue, residuez, invres, invresz, &
    partial_fraction_result, partial_fraction_tf_result

type(partial_fraction_result) :: pf
type(partial_fraction_tf_result) :: tf

! Analog: coefficients of descending powers of s.
pf = residue([1.0_dp, 2.0_dp], [1.0_dp, -2.0_dp, 1.0_dp])
! pf%residues = [1, 3], pf%poles = [1, 1]
tf = invres(pf%residues, pf%poles, pf%direct)

! Digital: coefficients of ascending powers of z**(-1).
pf = residuez([1.0_dp, 2.0_dp], [1.0_dp, -0.5_dp])
! pf%residues = [5], pf%poles = [0.5], pf%direct = [-4]
tf = invresz(pf%residues, pf%poles, pf%direct)
```

`residue(b, a, tol, rtype)` and `residuez(b, a, tol, rtype)` accept
rank-one `real(dp)` or `complex(dp)` coefficients. Their
`partial_fraction_result` contains `residues(:)`, `poles(:)`,
`direct(:)` (all `complex(dp)`), and an integer `status`. The direct
polynomial has the same coefficient ordering as the input. Poles are
sorted by increasing modulus and near-equal poles are grouped; within
a group the residues correspond to successive inverse powers from 1
to multiplicity. `tol` is an absolute pole-clustering distance (default
`1e-3`), and `rtype` may be `avg` (the default), `min`, or `max` (with
their SciPy long-form aliases).

`invres(r, p, k, tol, rtype)` and `invresz(...)` return
`partial_fraction_tf_result` with `b(:)`, `a(:)` (complex, denominator
monic), and `status`. The direct coefficient array `k` is optional and
may be real or complex; omitted `k` is an empty polynomial. Repeated
poles must be **adjacent** on input to the inverse functions, as in SciPy.
A zero numerator may be represented as a single zero coefficient.
Leading/trailing zero trimming can produce slightly different output
array lengths from NumPy's `polyadd`, while preserving the same transfer
function.

The implementation reuses SciFort's polynomial-root solver, clusters
numerically indistinguishable poles, obtains the polynomial quotient by
long division, and solves a square complex linear system for the residues.
The inverse functions multiply the polynomial factors explicitly. No
new external dependencies are introduced. Very high-order denominators,
closely spaced distinct poles, or large repeated-pole multiplicities can
be ill-conditioned: this version targets moderate-order rational filters,
not arbitrary-precision symbolic partial fractions. A failed root solve,
singular residue system, or nonfinite arithmetic returns
`partial_fraction_singular=2`. Invalid inputs return
`partial_fraction_invalid_input=1`; normal completion returns
`partial_fraction_success=0`. These are status codes rather than SciPy
exceptions.

Regression and cross-library validation:

```text
python tools/signal_reference/run_signal_suite.py
python tools/signal_reference/compare_signal_partial_fraction_scipy.py
```

The second script compiles a Fortran probe and compares both decompositions
and inverse transfer-function evaluations across 540 deterministic
randomized cases (analog/digital, real/complex coefficients, repeated
poles, proper and improper rational functions) with SciPy 1.17. The
reconstruction comparison uses a `2e-5` relative/absolute tolerance to
accommodate cancellation in poorly conditioned near-repeated designs;
individual residues are compared with looser relative tolerances than
well-separated simple poles. The regression test is
`test/test_signal_partial_fraction.f90`.

## State-space and transfer-function conversions

`scifort_signal` also provides SciPy-style `tf2ss`, `ss2tf`, `zpk2ss`,
and `ss2zpk`. The new result types store complex(dp) arrays, including
when all inputs are real, so the same interface handles real and complex
systems without discarding imaginary components.

```fortran
use scifort_signal, only : dp, tf2ss, ss2tf, zpk2ss, ss2zpk, &
    state_space_result, state_space_tf_result, state_space_success

type(state_space_result) :: model
type(state_space_tf_result) :: transfer

model = tf2ss([1.0_dp,3.0_dp,3.0_dp],[1.0_dp,2.0_dp,1.0_dp])
if (model%status /= state_space_success) error stop 'tf2ss failed'
! model%a = [[-2,-1],[1,0]] and model%c = [[1,2]].
transfer = ss2tf(model%a,model%b,model%c,model%d)
! transfer%b(1,:) = [1,3,3], transfer%a = [1,2,1].
```

* `tf2ss(b,a)` accepts 1D real(dp) or complex(dp) transfer-function
  coefficients in **descending powers of s or z** and returns a single-input,
  single-output controller-canonical realization. A rank-2 numerator with
  shape `(number_of_outputs, number_of_coefficients)` constructs a
  single-input/multiple-output realization with a common denominator.
  Leading zero coefficients are trimmed and `a(1)` is normalized. An
  improper transfer function or a zero denominator is rejected. A constant
  transfer function uses SciPy's one-dummy-zero-state convention.
* `ss2tf(A,B,C,D,input)` accepts either all-real or all-complex matrices,
  sized `(n,n)`, `(n,nin)`, `(nout,n)`, and `(nout,nin)`. Its optional input
  is **one-based** (default 1). It returns `state_space_tf_result%b(nout,n+1)`
  and `%a(n+1)`. The coefficient algorithm uses the Faddeev-LeVerrier
  characteristic polynomial and is intended for moderate state dimensions:
  poorly scaled high-order matrices can suffer numerical cancellation.
* `zpk2ss(z,p,k)` accepts real or complex roots and a real or complex gain;
  it expands polynomial factors with the existing `zpk2tf` routine and then
  creates the canonical state-space realization. Improper transfer
  functions are rejected.
* `ss2zpk(A,B,C,D,input,output)` returns the existing
  `filter_zpk_result`, including `zeros(:)`, `poles(:)`, `gain` and `status`.
  `input` and `output` are optional **one-based** indices (both default 1).
  Selecting one output is a Fortran convenience for multiple-output models;
  SciPy's `ss2zpk` is principally a single-output conversion.

`state_space_result` contains allocatable complex matrices `%a`, `%b`, `%c`,
`%d` and an integer `%status`. `state_space_tf_result` contains complex
`%b(:,:)`, `%a(:)` and `%status`. Status constants are
`state_space_success=0`, `state_space_invalid_input=1`, and
`state_space_numerical_failure=2`. Invalid dimensions, invalid channel
indices and nonfinite inputs are rejected without raising exceptions.
The conversion routines apply to both continuous and discrete LTI systems;
there is no sample-time parameter.

Validation commands (GNU Fortran 14.2, SciPy 1.17):

```text
python tools/signal_reference/run_signal_suite.py
python tools/signal_reference/run_signal_suite.py --optimized
python tools/signal_reference/compare_signal_state_space_scipy.py
```

The seeded reference comparison checks 280 `tf2ss`, 270 `ss2tf`,
138 `ss2zpk`, and 144 `zpk2ss` cases, including complex coefficients,
zero-state and multi-channel systems. On these cases the maximum scaled
matrix/coefficient discrepancy from SciPy 1.17 is approximately `1.08e-13`.

## Continuous-to-discrete LTI conversion

`cont2discrete` converts continuous-time state-space matrices or a
single-input transfer function (SISO or SIMO numerator rows) to a discrete system. The
corresponding ZPK entry point is `cont2discrete_zpk`: in Fortran, a ZPK
call with real-valued gain conflicts with the transfer-function generic
because their scalar and array argument ranks/types overlap. Both entry
points accept real or complex coefficients; result arrays are complex(dp)
so no imaginary information is lost.

```fortran
use scifort_signal, only : dp, cont2discrete, cont2discrete_zpk, &
    discrete_state_space_result, discrete_tf_result, discrete_zpk_result

type(discrete_state_space_result) :: discrete_model
type(discrete_tf_result) :: discrete_filter
type(discrete_zpk_result) :: roots

! State-space matrices: A is (n,n), B is (n,nin), C is (nout,n),
! D is (nout,nin); arbitrary positive real dt (sampling period).
discrete_filter = cont2discrete([1.0_dp], [1.0_dp, 2.0_dp, 1.0_dp], &
    0.1_dp, method='bilinear')
roots = cont2discrete_zpk([-0.5_dp], [-1.0_dp, -2.0_dp], 2.0_dp, &
    0.1_dp, method='zoh')
```

Supported methods:

- `zoh` (default): exact zero-order-hold state-space discretization;
- `foh`: exact first-order-hold discretization, including the additional
  input and feedthrough terms;
- `gbt` with required `alpha` between 0 and 1: generalized bilinear transform;
- `bilinear`/`tustin`: generalized bilinear with alpha=0.5;
- `euler`/`forward_diff`: generalized bilinear with alpha=0;
- `backward_diff`: generalized bilinear with alpha=1;
- `impulse`: impulse-invariant mapping for **strictly proper** systems
  (the continuous-time feedthrough matrix must be exactly zero).

Each result includes the sample interval `%dt` and `%status`. State-space
results provide `%a`, `%b`, `%c`, `%d`; TF results provide `%b(nout,n+1)` and
`%a(n+1)` in descending powers of z; ZPK results provide `%zeros`, `%poles`,
and `%gain`. Status constants are `cont2discrete_success=0`,
`cont2discrete_invalid_input=1`, `cont2discrete_numerical_failure=2`.
Nonpositive/nonfinite sample intervals, incompatible dimensions, unsupported
methods, invalid `gbt` weights and singular linear systems are rejected.

The ZOH and FOH algorithms use an independent scaled Padé(13) matrix
exponential; GBT uses a pivoted complex linear solver. GBT ZPK results use
an **algebraic root mapping** for poles, zeros, and gain, avoiding numerical
root-splitting of artificial repeated zeros. ZOH/FOH/impulse ZPK conversion
still uses the project polynomial-root solver, and can be less accurate on
high-order systems with clustered zeros. The matrix exponential and
Faddeev-LeVerrier TF reconstruction have the usual conditioning limitations
for large or highly nonnormal systems. No BLAS, LAPACK or other new
external dependency is required.

Validation (GNU Fortran 14.2.0, SciPy 1.17.0):

```text
python tools/signal_reference/run_signal_suite.py
python tools/signal_reference/run_signal_suite.py --optimized
python tools/signal_reference/compare_signal_cont2discrete_scipy.py
```

The standalone regression is `test/test_signal_cont2discrete.f90`. It covers
all methods, complex-valued systems, TF and ZPK overloads, zero-state models,
and invalid inputs. The 390-case batched seeded comparison checks 120 state-space, 130 SISO TF,
60 SIMO TF and 80 ZPK conversions against SciPy 1.17. Its largest scaled
response/matrix difference is about `8.51e-11`. Full FPM integration still requires `fpm test`
on a machine with dependencies installed.


## Discrete-time LTI system objects

Use the `dlti_system` derived type to retain a discrete-time LTI model in
transfer-function (TF), zeros/poles/gain (ZPK), or state-space (SS) form. The
separate `dlti_tf`, `dlti_zpk`, and `dlti_ss` constructors accept both real and
complex data, with optional positive sampling interval `dt` (default `1`, a
unit sampling interval). Polynomials use **descending positive powers of z**,
not delay-filter coefficients in ascending powers of z^-1. Fortran's generic
resolution cannot always distinguish `dlti(num,den,dt)` from
`dlti(zeros,poles,gain)` when `dt` and `gain` are both real, so the `dlti`
shorthand supports only TF and SS constructors; use `dlti_zpk` explicitly for
ZPK systems.

```fortran
use scifort_signal, only : dp, dlti_tf, dlti_zpk, dlti_ss, &
    dlti_system, discrete_simulation_result, frequency_response, dlti_success
implicit none

type(dlti_system) :: filter, realization, stable_filter
type(discrete_simulation_result) :: samples
type(frequency_response) :: spectrum
real(dp) :: input(128), frequencies(3)

filter = dlti_tf([0.2_dp,0.1_dp], [1.0_dp,-0.6_dp], dt=0.01_dp)
if (filter%status /= dlti_success) error stop 'invalid system'
realization = filter%to_ss()
stable_filter = filter%to_zpk()
! Fill input with 128 observations first.
samples = filter%output(input)
frequencies = [0.1_dp,0.5_dp,1.0_dp]
spectrum = filter%freqresp(w=frequencies)
! spectrum%frequencies is in radians per sample.
! The bode method instead returns frequencies in radians per second.
```

`to_tf()`, `to_zpk()`, and `to_ss()` return new `dlti_system` instances. The
original representation is retained until explicitly converted. Transfer
functions and ZPK represent **SISO** systems; the state-space constructor
accepts SISO, MISO, SIMO, and MIMO arrays `(A,B,C,D)`. Converting MIMO state
space to TF or ZPK is deliberately rejected, since selecting one input/output
channel would silently discard the other channels. Use the procedural
`ss2tf` or `ss2zpk` APIs when a particular channel is wanted.

Available type-bound methods:

- `output(u, t, x0)`: real/complex input vectors for SISO or rank-two
  `(samples,inputs)` matrices for multiple inputs; optional timestamps and
  complex initial state. Returns a `discrete_simulation_result` with
  `t`, `y(samples,outputs)`, `x(samples,states)` and `status`.
- `step(n, t, x0)` and `impulse(n, t, x0)`: response to each input channel
  separately, with `y(samples,outputs,inputs)`.
- `freqresp(w,n,whole)` and `bode(w,n)`: SISO frequency responses.
  `w` is in radians per sample; `bode%frequencies` is scaled by `1/dt`.
- `freqresp_mimo(w,n,whole,solver,inputs,outputs)` and
  `bode_mimo(w,n,solver,inputs,outputs)`: MIMO state-space frequency analysis,
  optional channel selections (one-based) and direct/Hessenberg/Schur solvers.
- `ninputs()`, `noutputs()`, `nstates()`, and `stable(tol)`. Stability checks
  require all poles, including internal state-space modes, to lie strictly
  inside the unit circle. `tol` optionally reserves a margin inside it.
- `series(upstream)`: cascade `self(upstream(u))` in that order, meaning
  `H_total = H_self * H_upstream`. `parallel(other)` sums the output of
  systems driven by identical inputs. Both use state-space connections and
  retain nonzero `D` feedthrough, with channel-shape and `dt` validation.
- `append(other)`: construct an independent-channel block-diagonal
  MIMO state-space realization. `select_io(inputs,outputs)` selects/reorders
  channels with one-based index arrays, including repeated channels.
- `feedback(controller,sign)`: close a MIMO loop with `u = r + sign*K*y`,
  `y = P*u` where `self=P`, `controller=K` and `sign` defaults to `-1`
  (negative feedback). The returned system maps the external input `r`
  to the plant output `y`. Handles nonzero plant/controller feedthrough
  using `(I-sign*Dk*Dp)^-1`; singular loops produce numerical-failure status.
  Both systems must use the same positive `dt`.
- `controllability()` and `observability()` return the complex Kalman
  matrices `[B,AB,...,A^(n-1)B]` and `[C;CA;...;CA^(n-1)]`, respectively.
- `minreal(tol)` removes numerically unreachable and unobservable modes with
  two orthogonal projections. Returns a new state-space `dlti_system`;
  default relative rank cutoff is `1e-10`. This is Kalman realization
  reduction, **not** approximate balanced model-order reduction. As with any
  numerical rank decision, poorly separated modes require care.
- `gramians(tol,max_iter)` computes the infinite-horizon discrete
  controllability Gramian `Wc`, observability Gramian `Wo`, and H2 impulse
  energy norm `sqrt(trace(C Wc C^H + D D^H))`. The returned
  `dlti_gramian_result` has allocatable `.wc`, `.wo`, scalar `.h2_norm`,
  `.iterations`, and `.status` members. It uses a doubling (Smith)
  iteration with a relative tolerance of `1e-12` by default and fails
  rather than returning a Gramian for an unstable/nonconvergent system.

### Balanced state-space truncation (SciFort extension)

`dlti_system%balanced_truncate(order, tol)` reduces a stable discrete-time
SISO or MIMO model to the requested state count using square-root balanced
truncation. `order` is required, with `0 <= order <= nstates()`. The optional
relative numerical cutoff `tol` defaults to `1e-10`. The result has fields:

- `model`: the reduced `dlti_system`, preserving input/output dimensions,
  direct feedthrough, and sampling interval;
- `hsv`: all Hankel singular values in descending order, padded with zeros
  for states numerically unreachable/unobservable;
- `error_bound`: `2*sum(hsv(order+1:))`, the standard infinity-norm error
  bound for **stable minimal** systems;
- `status`: `dlti_success`, `dlti_invalid_input`, or
  `dlti_numerical_failure` (including unstable and unreliable-rank systems);
- `transform(n,order)` and `inverse(order,n)`: biorthogonal state
  projections, with `inverse*transform` approximately equal to identity.
  For the reduced realization, `Ar=inverse*A*transform`,
  `Br=inverse*B`, and `Cr=C*transform`. These projections are also
  used by the continuous-time balanced-truncation implementation.

```fortran
use scifort_signal, only : dp, dlti_system, dlti_ss, &
    dlti_balanced_result, dlti_success
implicit none
type(dlti_system) :: original
type(dlti_balanced_result) :: reduced
real(dp) :: a(2,2),b(2,2),c(2,2),d(2,2)
a = reshape([0.4_dp,0.0_dp,0.0_dp,0.2_dp],[2,2])
b = reshape([1.0_dp,0.0_dp,0.0_dp,0.1_dp],[2,2])
c = b
d = 0.0_dp
original = dlti_ss(a,b,c,d,dt=0.01_dp)
reduced = original%balanced_truncate(order=1)
if (reduced%status == dlti_success) then
    print *, 'Hankel singular values:', reduced%hsv
    print *, 'Error bound:', reduced%error_bound
end if
```

The algorithm diagonalizes complex Hermitian Gramians and uses an eigenvalue
problem for the squared singular values of their factor product. It does
not require an external LAPACK library, but very small Hankel singular
values (roughly below the square root of machine precision relative to
those retained) may be inaccurate. The implementation refuses requested
orders beyond the estimated numerical rank. For high-order or strongly
non-normal systems, use a LAPACK-based SVD in a future backend or check the
reduced frequency response against the original. The full state-space
realization must be asymptotically stable, including hidden modes.

For example:

```fortran
use scifort_signal, only : dlti_system, dlti_gramian_result, dlti_success
! Assume plant and controller are previously constructed compatible systems.
type(dlti_system) :: plant, controller, closed_loop, minimal
type(dlti_gramian_result) :: measures
closed_loop = plant%feedback(controller) ! negative feedback
if (closed_loop%status /= dlti_success) error stop 'invalid closed loop'
minimal = closed_loop%minreal()
measures = minimal%gramians()
if (measures%status == dlti_success) print *, 'H2 norm = ', measures%h2_norm
```

All methods return invalid status or `false` on invalid data; they do not
silently infer a compatible sampling interval. Access to the stored arrays
is supported, but treat them as read-only to preserve constructor invariants.
ZPK and TF conversion and internal eigenvalue-based stability detection can
be poorly conditioned for high orders; state-space MIMO frequency analysis
avoids polynomial expansion.

Regression commands (GNU Fortran and SciPy/NumPy required):

```shell
python tools/signal_reference/verify_dlti_objects.py
python tools/signal_reference/verify_dlti_objects.py --optimized
python tools/signal_reference/verify_dlti_compositions.py
python tools/signal_reference/verify_dlti_compositions.py --optimized
# Compile dlti_interconnections_driver.f90 and dlti_minreal_driver.f90 first:
python tools/signal_reference/verify_dlti_interconnections.py path/to/dlti_interconnections_driver
python tools/signal_reference/verify_dlti_minreal.py path/to/dlti_minreal_driver
python tools/signal_reference/run_signal_suite.py --jobs 4
```

The reference suites cover 168 deterministic real/complex SISO and MIMO
systems plus 110 randomized MIMO series/parallel compositions. It uses SciPy's state-space realizations to independently evaluate
state trajectories, and SciPy's public frequency-response methods where
usable. SciPy 1.17 contains complex-casting limitations in selected ZPK
frequency and discrete simulation paths; these cases are independently
checked using the equivalent complex state-space recurrences and SciPy's
`zpk2tf` coefficients.

## SOS-native discrete step and impulse responses (SciFort extension)

Digital SOS matrices can now be passed directly to the `dstep`/`dimpulse`
Fortran generics, or to `dstep_sos`/`dimpulse_sos`. Both real and complex
coefficients are supported. They use the existing stable section-wise
`sosfilt` kernel; there is no high-order polynomial expansion. The optional
`zi(nsections,2)` stores direct-form-II section delays (not canonical `x0`).
See [SOS response details and examples](signal_sos_responses.md).

## Discrete-time LTI simulation: `dlsim`, `dstep`, and `dimpulse`

The `scifort_signal_discrete_simulation` module exports the SciPy-style
`dlsim`, `dstep`, and `dimpulse` routines through `scifort_signal`. All three
accept either real(dp) or complex(dp) discrete **state-space** matrices
`(A,B,C,D,dt)` or SISO **transfer-function** coefficients `(b,a,dt)`;
the corresponding zero-pole-gain entry points are `dlsim_zpk`,
`dstep_zpk`, and `dimpulse_zpk`. The ZPK functions have their own names
to avoid Fortran generic ambiguities. Polynomial coefficients use descending
powers of **z**, as in SciPy, and must describe a proper transfer function.

```fortran
use scifort_signal, only : dp, dlsim, dstep, dimpulse, &
    discrete_simulation_result, discrete_response_result, discrete_sim_success

type(discrete_simulation_result) :: simulation
type(discrete_response_result) :: response

! Single-input state-space: A=(1,1), B=(1,1), C=(1,1), D=(1,1),
! U=(number of input samples, number of inputs).
simulation = dlsim(reshape([0.5_dp],[1,1]), &
    reshape([1.0_dp],[1,1]),reshape([1.0_dp],[1,1]), &
    reshape([0.0_dp],[1,1]),reshape([1.0_dp,2.0_dp,0.0_dp],[3,1]),0.2_dp)

! Transfer function H(z) = 1/(z-0.5): step and impulse responses.
response = dstep([0.0_dp,1.0_dp],[1.0_dp,-0.5_dp],0.2_dp,n=25)
response = dimpulse([0.0_dp,1.0_dp],[1.0_dp,-0.5_dp],0.2_dp,n=25)
```

The **SciFort extension** `dlsim_sos(sos,u,dt,t,zi)` (also available as
`dlsim(sos,u,dt,t,zi)`) simulates arbitrary SISO input through real or
complex SOS sections without transfer-function expansion. The returned
`%zf(nsections,2)` contains final DF-II states for continuation. See
[`signal_dlsim_sos.md`](signal_dlsim_sos.md) for examples and state layout.

The **SciFort extension** `dlsim_sos_bank` (and new `dlsim_sos`/`dlsim`
rank-two-input overloads) simulates multiple independent channels with
shared SOS coefficients `(nsections,6)` or separate channel coefficients
`(nsections,6,nchannels)`. It provides channel-wise filter-state histories
and final states for blockwise continuation. See
[`signal_dlsim_sos_bank.md`](signal_dlsim_sos_bank.md).

`discrete_simulation_result` has real `%t(nt)` and complex `%y(nt,nout)`
and `%x(nt,nstates)` arrays, where **the state at sample k is recorded
before the k-th input update**. `discrete_response_result` provides real
`%t(nt)` and complex `%y(nt,nout,nin)` arrays: the third index selects the
input channel excited by a unit step or unit impulse. This representation
provides all MIMO experiments in one call. For SISO `dlsim` with transfer
functions or ZPK, the input is a rank-1 array; for state-space `dlsim`,
the input is a rank-2 array with one column per input.

Optional arguments: `dlsim(...,t=t_input,x0=x_initial)` accepts a strictly
increasing real input-time array and an initial state vector. With input
times, SciPy-style piecewise-linear interpolation evaluates inputs at
`linspace(0,t_input(end),floor(t_input(end)/dt)+1)` and returns the same
output grid. Without `t`, samples are spaced exactly `dt` apart. `dstep`
and `dimpulse` accept `n=` (default 100), `t=` (overrides `n`), and `x0=`.
When neither `t` nor `x0` is given, the initial state is zero. State-space
calls support MIMO and zero-dimensional states; TF and ZPK entry points
handle SISO systems. Real and complex inputs are both supported, and the
results retain complex parts even when the state matrices are real.

Status constants are `discrete_sim_success=0`,
`discrete_sim_invalid_input=1`, and `discrete_sim_numerical_failure=2`.
Nonfinite coefficients, invalid dimensions, nonpositive `dt`, wrong state
lengths, and invalid/duplicate input times return status 1. Nonfinite
state/output trajectories return status 2. Transfer-function and ZPK initial
states use SciFort's controller-canonical realization; changing realization
changes the interpretation of `x0` but not zero-state input/output responses.

SciPy 1.17 has two minor edge cases: `dlsim` casts **complex** trajectories
to real, losing imaginary parts; SciFort preserves them. `dstep`/`dimpulse`
can return one fewer time sample than the requested `n` for some fractional
sampling periods because of `floor()` after floating-point endpoint rounding;
SciFort instead returns the requested `n` when no explicit `t` is supplied.

Validation:

```text
python tools/signal_reference/run_signal_suite.py
python tools/signal_reference/run_signal_suite.py --optimized
python tools/signal_reference/compare_signal_discrete_simulation_scipy.py
```

The fixed Fortran regression is `test/test_signal_discrete_simulation.f90`.
The comparison driver checks **498 seeded cases** against SciPy 1.17 for
real systems and the same recurrence applied to SciPy's state-space matrices
for complex systems. The maximum absolute discrepancy in that comparison was
approximately `1.14e-14`. The routines require no new FPM dependencies.

## Continuous-time LTI simulation: `lsim`, `step`, and `impulse`

The `scifort_signal_continuous_simulation` module adds continuous-time
`lsim`, `step`, and `impulse`. Like the discrete simulation routines,
`lsim` accepts either state-space matrices `(A,B,C,D,U,T)` (with
`U(nt,nin)`) or SISO transfer-function polynomials `(b,a,U,T)` (with
`U(nt)`). Zeros/poles/gain are accepted through the distinct
`lsim_zpk(z,p,k,U,T)` entry point; likewise `step_zpk` and `impulse_zpk`.
Transfer-function polynomials use descending positive powers of **s**.
All functions support real and complex state-space/TF/ZPK inputs.

```fortran
use scifort_signal, only : dp, lsim, step, impulse, &
    continuous_simulation_result, continuous_response_result, continuous_sim_success

type(continuous_simulation_result) :: trajectory
type(continuous_response_result) :: response
real(dp) :: t(101), u(101)
integer :: i

do i = 1, size(t)
    t(i) = real(i-1,dp)*0.01_dp
end do
u = t  ! linear input ramp

! H(s) = 1/(s+1). The default input interpolation is linear.
trajectory = lsim([1.0_dp],[1.0_dp,1.0_dp],u,t)
if (trajectory%status /= continuous_sim_success) error stop 'lsim failed'

response = step([1.0_dp],[1.0_dp,1.0_dp],t)
response = impulse([1.0_dp],[1.0_dp,1.0_dp],t)

! Optional automatic time grid, with n=100 points by default.
response = step([1.0_dp],[1.0_dp,1.0_dp])
response = impulse([1.0_dp],[1.0_dp,1.0_dp],n=200)
```

`continuous_simulation_result` contains real `%t(nt)` and complex
`%y(nt,nout)` and `%x(nt,nstates)` arrays. `continuous_response_result`
contains real `%t(nt)` and complex `%y(nt,nout,nin)` arrays: the third
index selects the input channel receiving a unit step or Dirac impulse.
The response routines evaluate all MIMO input channels independently,
and accept irregularly spaced response times. Both types expose `%status`:
`continuous_sim_success=0`, `continuous_sim_invalid_input=1`, and
`continuous_sim_numerical_failure=2`.

`lsim` requires a nonempty, nonnegative, strictly increasing and (for two
or more intervals) **uniformly spaced** time grid, as in SciPy. The
optional `interp=` argument is `.true.` by default, giving piecewise-linear
(FOH) inputs. Set `interp=.false.` for zero-order-held input samples.
Transition matrices use the same Padé matrix exponential as `cont2discrete`,
including an augmented matrix for exact FOH propagation. Optional `x0=`
specifies the continuous state at time **zero**; when the first requested
sample is at `t(1)>0`, it is propagated with zero input beforehand.
TF/ZPK initial states refer to SciFort's controller-canonical realization.
`step`, `impulse`, `step_zpk` and `impulse_zpk` allow `t` to be omitted.
With no `t`, the optional Fortran integer `n=` specifies the
number of equally spaced output samples (default: 100, must be positive).
The default grid starts at zero and ends at `7/r`, where `r` is the
smallest absolute real part of a state-matrix eigenvalue. If a pole lies
on the imaginary axis or an eigenvalue has zero real part, `r=1` is used,
as in SciPy. One requested sample gives `[0]`; when `t` is provided,
`n=` is ignored and the exact requested (possibly irregular) grid is
retained. A zero-dimensional state-space model uses a seven-unit range.
The implementation uses exact diagonal roots for triangular matrices
and a characteristic-polynomial root solver otherwise; high-order or
nearly defective state matrices may yield less accurate automatic grid
endpoints. Supply explicit `t` for applications requiring precise horizons.

For ZPK calls, zeros and poles must have the same Fortran element type,
but the gain may independently be real or complex; `lsim_zpk` also accepts
real or complex input arrays independently of root/gain types. Complex
initial states can be passed with `x0=` to the mixed-type variants.

A unit-step response includes the direct-feedthrough matrix `D`. `impulse`
returns the *regular* impulse response `C exp(A t) B` (plus any response to
`x0`): a nonzero `D` represents an additional Dirac delta at `t=0`, which
is **not** represented by finite sampled values, matching SciPy's impulse
routine convention. Empty state matrices (static gains) are supported.
Wrong matrix dimensions, improper transfer functions, nonfinite values,
nonmonotone/nonuniform time grids, and invalid initial-state lengths
return invalid status without an exception.

Validation (SciPy 1.17 reference, GNU Fortran):

```text
python tools/signal_reference/run_signal_suite.py
python tools/signal_reference/run_signal_suite.py --optimized
python tools/signal_reference/compare_signal_continuous_simulation_scipy.py
python tools/signal_reference/compare_signal_continuous_auto_time_scipy.py
```

The fixed regression `test/test_signal_continuous_simulation.f90`
covers analytic real and complex responses, zero/first-order holds,
nonzero initial times, MIMO channels, static systems and error handling.
Additional tests `test_signal_continuous_auto_time.f90` and
`test_signal_continuous_mixed_zpk.f90` check automatic grid selection,
zero/imaginary poles, invalid sample counts, and mixed real/complex
ZPK gain/input overloads. The 184-case automatic-grid comparison uses
SciPy 1.17 and independent matrix-exponential reference responses.
The earlier randomized driver checks 522 cases against SciPy's LTI routines or
SciPy's matrix exponential (for complex input/MIMO variants not reliably
handled by the SciPy LTI convenience functions). No new FPM dependencies
are required.

## Transfer-function, ZPK, and digital SOS conversions

`scifort_signal` exports five SciPy-style conversion routines with typed
outputs and explicit status codes:

| Fortran routine | Result type | Purpose |
| --- | --- | --- |
| `tf2sos(b,a,pairing,analog)` | `filter_sos_result` | Transfer-function coefficients to digital or analog SOS |
| `tf2zpk(b,a)` | `filter_zpk_result` | Transfer-function polynomial coefficients to zeros, poles, gain |
| `zpk2tf(z,p,k)` | `filter_tf_result` | Zeros, poles, gain to monic-denominator polynomial coefficients |
| `zpk2sos(z,p,k,pairing)` | `filter_sos_result` | Digital zero/pole/gain representation to biquadratic sections |
| `sos2zpk(sos)` | `filter_zpk_result` | Cascade of digital SOS sections to zeros, poles, gain |

```fortran
use scifort_signal, only : dp, tf2sos, tf2zpk, zpk2tf, zpk2sos, sos2zpk, &
    filter_zpk_result, filter_tf_result, filter_sos_result, zpk_conversion_success

type(filter_zpk_result) :: zp
type(filter_tf_result) :: tf
type(filter_sos_result) :: sections

zp = tf2zpk([1.0_dp, 0.4_dp, 0.1_dp], [1.0_dp, -0.7_dp, 0.2_dp])
if (zp%status /= zpk_conversion_success) error stop 'tf2zpk failed'
tf = zpk2tf(zp%zeros, zp%poles, zp%gain)
sections = zpk2sos(zp%zeros, zp%poles, zp%gain)
if (sections%status /= zpk_conversion_success) error stop 'zpk2sos failed'
zp = sos2zpk(sections%sos)
if (zp%status /= zpk_conversion_success) error stop 'sos2zpk failed'
sections = tf2sos([1.0_dp], [1.0_dp,-0.8_dp,0.2_dp])
if (sections%status /= zpk_conversion_success) error stop 'tf2sos failed'
! Analog SOS sections use descending powers of s instead of z**(-1).
sections = tf2sos([2.0_dp], [1.0_dp,3.0_dp], analog=.true.)
```

`tf2sos` accepts real(dp) or complex(dp) coefficient vectors. By default it
returns digital sections with `pairing='nearest'`. `pairing='keep_odd'` is
also supported, with a nearest-pair fallback when unequal polynomial orders
prevent a one-to-one match of first-order factors. `pairing='minimal'`
returns polynomial SOS factors that must not be used by digital `sosfilt`.
For analog transfer functions specify `analog=.true.`; this selects minimal
pairing and rejects `nearest`/`keep_odd` as SciPy does. Analog sections are
polynomials in descending positive powers of `s`; evaluate them using
Horner's rule at `s=iw`, not `sosfreqz`. Constant analog systems have the row
`[0,0,k,0,0,1]`. Like SciPy, `tf2sos` removes leading zero numerator
coefficients during TF normalization. For numerically sensitive high-order
filters, design directly in ZPK/SOS form rather than converting from an
expanded BA polynomial when possible.

`filter_zpk_result` contains complex `zeros(:)`, `poles(:)` and `gain`;
`filter_tf_result` contains complex `b(:)` and `a(:)`; `filter_sos_result`
contains complex `sos(nsections,6)`. To use an SOS cascade with *real*
inputs, check that its imaginary coefficients are negligible and pass
`real(sections%sos,dp)` to the existing real `sosfilt` or `sosfreqz`.
All routines accept real or complex roots and coefficients through specific
Fortran interfaces; an explicit complex gain is supported.

**Polynomial convention:** `tf2zpk` and `zpk2tf` use polynomial coefficients
in *descending positive powers* of `s` or `z`. Digital SOS rows represent
*increasing negative powers* of `z`. As in SciPy's default digital
`zpk2sos(pairing='nearest')`, unequal numerator and denominator root counts
are padded at `z=0`; the resulting SOS expresses the usual digital
negative-power convention, which need not match the original positive-power
rational expression if its degrees differ. When the original degrees match,
`sosfreqz` reproduces `freqz_zpk` to numerical precision.

Digital `zpk2sos` accepts `pairing='nearest'` (default), `keep_odd`, or
`minimal` with conjugate-aware group pairing. `minimal` does **not** add roots
at zero and requires `size(p)>=size(z)` (otherwise invalid input). For odd or
unequal orders, its SOS rows may have leading zero numerator/denominator
coefficients, including `sos(:,4)=0`. Such rows are polynomial factors in
*descending positive powers of z*, not normal causal digital sections; they
**must not** be passed to `sosfilt` or `sosfreqz`. To analyze them, use
`[b0,b1,b2]` and `[a0,a1,a2]` as numerator and denominator polynomials in
`z`, or convert back using `sos2zpk` and remove bookkeeping origin roots.
`nearest`/`keep_odd` remain appropriate for digital filtering. Section order
and individual coefficients can differ from SciPy even when responses agree.

`sos2zpk` accepts real or complex SOS matrices shaped `(nsections,6)` and
returns two complex zeros and poles per section, including **zeros at the
origin to pad absent polynomial roots**, exactly as SciPy does. Each section's
leading denominator coefficient is normalized; the section gains are
multiplied. SOS matrices with zero denominators, nonfinite values, or the
wrong shape return invalid input. Gain overflow returns numerical failure.
For direct analog ZPK-to-SOS conversion, use `zpk2sos_analog` instead of
`zpk2sos`; `tf2sos(analog=.true.)` is the TF-to-SOS convenience route. High polynomial orders and clustered roots may be ill-conditioned;
the current root solver is intended for moderate orders.

Statuses: `zpk_conversion_success=0`, `zpk_conversion_invalid_input=1`,
`zpk_conversion_numerical_failure=2`. Failed calls return the status rather
than raising a Python-style exception.

Numerical verification (SciPy 1.17 and GNU Fortran):

```text
python tools/signal_reference/run_signal_suite.py
python tools/signal_reference/run_signal_suite.py --optimized
python tools/signal_reference/compare_signal_zpk_conversion_scipy.py
python tools/signal_reference/compare_signal_sos2zpk_minimal_scipy.py
python tools/signal_reference/compare_signal_tf2sos_scipy.py
```

The seeded randomized comparison covers 840 distinct real/complex designs,
including unequal degrees, odd orders, coefficient/root reconstruction, and
SOS frequency-response equivalence. The extended comparison tests 792
minimal-pairing designs and 759 SOS inverses (including section scaling and
leading zeros) against SciPy 1.17. The TF-to-SOS comparison covers 992 valid
analog/digital coefficient sets plus 64 deliberately improper minimal-pairing
requests. SciPy produces reference SOS matrices for 903 valid sets; for 89
cases where SciPy's own pole pairing raises an error, the Fortran response is
compared directly against the input polynomial ratio. Across all valid sets,
the largest scaled frequency-response difference was about `5.04e-7`,
primarily in the higher-order cases.

## Continuous-time LTI system objects

`clti_system` is the continuous-time counterpart of `dlti_system`. Use
`clti_tf`, `clti_zpk` or `clti_ss` to construct it from analog transfer
coefficients (descending powers of **s**), Laplace-plane zeros/poles/gain,
or complex/real MIMO state-space matrices `(A,B,C,D)`. The TF/SS-only
shorthand `clti` avoids ambiguous Fortran generic overloads. The type
returns numerical status from `%status()`, rather than a public component.

```fortran
use scifort_signal, only : dp, clti_system, clti_tf, dlti_system, &
    clti_success, bode_result
implicit none
type(clti_system) :: analog
type(dlti_system) :: sampled
type(bode_result) :: curve
analog = clti_tf([2.0_dp],[1.0_dp,3.0_dp])
if (analog%status() /= clti_success) error stop 'invalid analog model'
curve = analog%bode(w=[0.1_dp,1.0_dp,10.0_dp])
sampled = analog%to_discrete(dt=0.02_dp,method='zoh')
```

Type-bound methods include `to_ss`, `to_tf`, `to_zpk`, `freqresp`, `bode`,
`freqresp_mimo`, `bode_mimo`, `output`, `step`, `impulse`, `ninputs`,
`noutputs`, `nstates`, `stable`, `series`, `parallel`, `feedback`, `append`,
`select_io`, `minreal`, `controllability`, `observability`, `gramians`,
`balanced_truncate`, and `to_discrete`. `output(u,t,x0,interp)` accepts
real/complex SISO vectors or `(samples,inputs)` matrices and a mandatory
nonnegative, evenly spaced time grid; `interp` defaults to linear input
interpolation, while `interp=.false.` requests zero-order hold. SISO and
MIMO analog frequencies are in **radians per second**, not per sample.

Continuous stability checks require strictly negative real parts of the
poles (including hidden state-space modes), with optional margin `tol`.
Connections use the same state-space algebra as discrete systems but no
sample interval is involved. The `to_discrete` method accepts the same
`method` and `alpha` options as the procedural `cont2discrete` routine.
MIMO TF/ZPK conversion is not provided; select a channel or use `ss2tf`.
`clti_system` does not expose *discrete* Gramians or discrete
Bode methods.

### Continuous-time Gramians and balanced truncation (SciFort extensions)

`clti_system%gramians(tol,max_iter)` solves the two continuous-time
Lyapunov equations for a **Hurwitz-stable** realization (including hidden
modes):

```text
A*Wc + Wc*A^H + B*B^H = 0
A^H*Wo + Wo*A + C^H*C = 0
```

It returns a `clti_gramian_result` with complex `wc`, `wo`, real
`h2_norm`, the doubling `iterations`, and `status`. For strictly proper
models (`D=0`), `h2_norm=sqrt(trace(C*Wc*C^H))`; when any direct
feedthrough element is nonzero, the **continuous-time** H2 norm is infinite
because the impulse response contains a Dirac distribution. A zero-state
system with `D=0` has zero H2 norm. The continuous-time convention differs
from the discrete-time H2 norm, which includes the finite impulse at step 0.

`clti_system%balanced_truncate(order,tol)` returns a `clti_balanced_result`
containing a reduced continuous-time `model`, descending Hankel singular
values `hsv`, `error_bound=2*sum(hsv(order+1:))`, and `status`. The bound is
valid for asymptotically stable minimal models. `order` may range from zero
to the original state count. The direct feedthrough matrix and input/output
channel counts are preserved.

```fortran
use scifort_signal, only : dp, clti_system, clti_ss, &
    clti_gramian_result, clti_balanced_result, clti_success
implicit none
type(clti_system) :: system
type(clti_gramian_result) :: gr
type(clti_balanced_result) :: red
real(dp) :: a(2,2),b(2,2),c(2,2),d(2,2)
a = reshape([-1.0_dp,0.0_dp,0.0_dp,-4.0_dp],[2,2])
b = reshape([1.0_dp,0.0_dp,0.0_dp,0.3_dp],[2,2])
c = b
d = 0.0_dp
system = clti_ss(a,b,c,d)
gr = system%gramians()
if (gr%status /= clti_success) error stop 'Gramian solve failed'
red = system%balanced_truncate(order=1)
if (red%status /= clti_success) error stop 'balanced reduction failed'
print *, 'Hankel singular values:', red%hsv
print *, 'continuous H2 norm:', gr%h2_norm
```

The solver maps the continuous Lyapunov equations *exactly* to equivalent
discrete Lyapunov equations using a scaled bilinear/Cayley transform. It then
reuses SciFort's Smith doubling solver. The balancing projections are applied
to the **original analog** `A`, `B`, and `C` matrices, rather than applying
an inverse Tustin transform to the reduced discrete realization. This is
continuous-time square-root balanced truncation, not discretization followed
by model reduction. No new external dependencies are required.

As with the discrete implementation, a Hermitian eigendecomposition of a
normal matrix replaces a direct SVD, so tiny or poorly separated Hankel
singular values are less accurate. Extremely stiff or marginally stable
systems can fail the doubling convergence test. Invalid, unstable, and
numerically nonconvergent cases return an error status rather than silently
providing incorrect Gramians or reduced models.

Verification of these extensions: `test/test_signal_clti_gramians.f90` and
`tools/signal_reference/verify_clti_gramians.py`, compared with SciPy's
continuous Lyapunov solver and an independent square-root SVD balance.

Verification: `test_signal_clti.f90`, `test_signal_dlti_balanced.f90`,
`tools/signal_lti_objects_reference_driver.f90`, and
`tools/compare_signal_lti_objects.py` (100 random complex MIMO cases of each
kind, with independent SciPy/NumPy references).

## Continuous/discrete LTI frequency responses and Bode plots

`scifort_signal` now exports `freqresp`, `bode`, `dfreqresp`, and `dbode`.
Each generic accepts `(b,a)` transfer-function polynomials or `(A,B,C,D)`
state-space matrices, with either all-real or all-complex coefficients. A
system must have **one input and one output**. For root-based models, use
`freqresp_zpk`, `bode_zpk`, `dfreqresp_zpk`, or `dbode_zpk` with `(z,p,k)`.
Real or complex zeros/poles and gains are supported (the root arrays must
currently have the same Fortran kind/type). Transfer-function and ZPK wrappers reuse
SciFort's existing analog/digital response evaluators. State-space wrappers
evaluate the resolvent directly, solving `(i*w*I-A)*x=B` (analog) or
`(exp(i*w)*I-A)*x=B` (digital), and return `C*x+D` without expanding
the characteristic polynomial. No additional dependencies are required.

```fortran
use scifort_signal, only : dp, freqresp, bode, dfreqresp, dbode, &
    frequency_response, bode_result, design_status_success

type(frequency_response) :: analog, digital
type(bode_result) :: response

analog = freqresp([1.0_dp],[1.0_dp,2.0_dp],n=128)
response = bode([1.0_dp],[1.0_dp,2.0_dp])
digital = dfreqresp([1.0_dp],[1.0_dp,-0.5_dp],n=256)
response = dbode([1.0_dp],[1.0_dp,-0.5_dp],dt=0.01_dp)
if (response%status /= design_status_success) error stop 'dbode failed'
```

**Coefficient and frequency conventions:** Analog TF coefficients are in
**descending powers of s** and their frequencies are rad/s. Digital TF
coefficients are in **descending powers of z** (such as `z - 0.5` encoded
by `[1,-0.5]`). They are padded to a common polynomial degree and converted
to delay coefficients before reusing `freqz`. `dfreqresp%frequencies` are in
**radians per sample** regardless of sampling interval, whereas
`dbode%frequencies` are scaled to **radians per second**, dividing by
optional positive `dt` (default 1). For the discrete response, optional
`whole=.true.` evaluates the full unit circle instead of just `[0,pi)`.

Each frequency-response result has `%frequencies(:)`, `%response(:)`
(complex), and `%status`. `bode_result` has `%frequencies(:)`,
`%magnitude(:)` (dB), `%phase(:)` (degrees, unwrapped), and `%status`.
At exactly zero amplitude the dB magnitude is negative infinity. Invalid
inputs or singular frequency evaluations are signalled through nonzero
status codes; inspect `%status` before accessing outputs. The existing
`design_status_success`, `design_status_invalid_input`, and
`design_status_singular` constants apply.

With no explicit `w`, `freqresp` uses an automatically chosen logarithmic
analog grid with `n=10000` points, while `bode` uses `n=100` points.
`dfreqresp` and `dbode` likewise default to 10000 and 100 equally spaced
frequency points respectively. Explicit `w` takes precedence over `n`.
For state-space inputs without an explicit `w`, the analog logarithmic
frequency grid is estimated from nonzero state-matrix entry scales, without
forming the characteristic polynomial. This grid need not match the
transfer-function/ZPK automatic grid point-for-point; supply `w` for exact
comparison. The direct solve uses complex double precision and partial
pivoting; nearly singular resolvents remain ill-conditioned, and dense
frequency sweeps cost O(nfreq*nstate**3) arithmetic. The original generic
response still represents a single input and single output; separate MIMO
interfaces are documented below. If an
evaluation matrix is singular or produces nonfinite results, the response
status is `design_status_singular` and the affected response entries are
NaN. In comparison, high-order polynomial conversion can overflow even
when the original state-space frequency response is well behaved.

Verification (GNU Fortran and SciPy 1.17):

```text
python tools/signal_reference/run_signal_suite.py
python tools/signal_reference/run_signal_suite.py --optimized
python tools/signal_reference/compare_signal_lti_frequency_scipy.py
python tools/signal_reference/compare_signal_lti_resolvent_scipy.py
```

The new standalone Fortran test checks state-space, TF and ZPK responses,
real and complex systems, phase unwrapping, automatic frequencies, full
unit-circle grids, and invalid inputs. The seeded reference comparison
covers 720 independent analog/digital TF/ZPK and Bode evaluations.
The additional direct-resolvent test checks high-order state matrices,
complex transfer paths, singular frequencies, and 560 randomized real/complex
state-space systems (5,600 frequency points) against independent SciPy
linear solves.

## MIMO state-space frequency responses and Bode plots

SciFort extends the SISO `scipy.signal.freqresp`/`bode` interface with four
**separate** state-space entry points: `freqresp_mimo`, `dfreqresp_mimo`,
`bode_mimo`, and `dbode_mimo`. Both real `real(dp)` and complex `complex(dp)`
matrices are accepted. Using separate names preserves all existing SISO calls
and return types; SciPy's corresponding convenience frequency-response
functions do not directly support general MIMO state-space systems.

A state-space model has `A(nstate,nstate)`, `B(nstate,ninput)`,
`C(noutput,nstate)`, and `D(noutput,ninput)`. The response type is
`mimo_frequency_response` with `%frequencies(:)`,
`%response(noutput,ninput,nfreq)`, and `%status`. The Bode result type is
`mimo_bode_result` with `%magnitude(noutput,ninput,nfreq)` in decibels,
`%phase(noutput,ninput,nfreq)` in independently unwrapped degrees,
`%frequencies(:)` and `%status`. Indices are **one-based**, so the frequency
response from input 2 to output 1 at frequency index `k` is
`result%response(1,2,k)`.

```fortran
use scifort_signal, only : dp, freqresp_mimo, dbode_mimo, &
    mimo_frequency_response, mimo_bode_result, design_status_success

real(dp) :: a(2,2), b(2,2), c(2,2), d(2,2)
type(mimo_frequency_response) :: fr
type(mimo_bode_result) :: br

a = reshape([-1.0_dp,0.0_dp,0.0_dp,-2.0_dp],[2,2])
b = reshape([1.0_dp,0.0_dp,0.0_dp,1.0_dp],[2,2])
c = b
d = 0.0_dp

fr = freqresp_mimo(a,b,c,d,w=[0.0_dp,1.0_dp,2.0_dp])
if (fr%status /= design_status_success) error stop 'frequency response'
print *, fr%response(1,1,2)

br = dbode_mimo(a,b,c,d,n=128,dt=0.01_dp)
if (br%status /= design_status_success) error stop 'Bode response'
print *, br%magnitude(1,1,1)
```

Analog frequencies are angular frequencies in **rad/s**. Digital
`dfreqresp_mimo` frequencies use **radians per sample**, independently of
`dt`. Digital `dbode_mimo` returns **rad/s**, dividing by its positive `dt`
(default 1). The optional `whole=.true.` is available on `dfreqresp_mimo` and
covers the full unit circle. Explicit `w` takes precedence over `n`;
frequency-response defaults use `n=10000`, Bode defaults use `n=100`.
Analog automatic grids are estimated from state-matrix entry scales.

Each MIMO entry point now accepts an optional `solver=` keyword:

- `solver='auto'` (default): use reusable Hessenberg reduction for
  `nstate >= 8` and at least eight frequencies, else the direct solver.
- `solver='direct'`: use the original pivoted dense shifted solve.
- `solver='hessenberg'`: apply complex unitary Householder similarities
  once and solve the resulting upper Hessenberg shifted systems using
  adjacent-row pivoting at each frequency. This method also supports
  purely real systems without a separate API.
- `solver='schur'`: additionally reduce Hessenberg form to **complex Schur
  triangular form** by shifted, unitary QR iterations, then use triangular
  back substitutions for each frequency. The preparation cost is greater,
  but repeated or long sweeps may be faster. If QR iteration does not converge
  within its bounded iteration budget, the implementation falls back to the
  valid Hessenberg solver rather than returning inaccurate Schur results.

For example, `fr = freqresp_mimo(a,b,c,d,n=600,solver='schur')`.
`auto` deliberately retains its earlier Hessenberg/direct selection, avoiding
Schur setup cost for typical one-off requests. Unrecognized solver names return
`design_status_invalid_input`.

### Select channels and reuse a prepared reduction

The four MIMO routines accept optional `inputs(:)` and `outputs(:)` integer
arrays. Indices are **one-based**. Channel selection is applied before solving
and producing results; the result's first two dimensions follow the selected
output/input order. Indices may be repeated. Both arrays must be nonempty
and in range when specified. Invalid selectors return
`design_status_invalid_input` without allocating response data.

To avoid another Householder reduction for each frequency sweep, use the
independent `prepare_mimo_frequency(A,B,C,D)` API, which accepts real or
complex state-space matrices and returns a `mimo_frequency_plan`. By default
it stores Hessenberg form; specify `solver='schur'` to request a reusable Schur
form instead. Inspect `plan%method` to see whether Schur converged or the
safe Hessenberg fallback was used. The same plan can be passed as the
**first argument** to `freqresp_mimo`,
`dfreqresp_mimo`, `bode_mimo`, or `dbode_mimo` to evaluate any number of
frequency grids and channel subsets, without repeating the reduction.
A prepared plan is intended to be treated as immutable. To change system
coefficients, call `prepare_mimo_frequency` again.

```fortran
use scifort_signal, only : dp, prepare_mimo_frequency, mimo_frequency_plan, &
    freqresp_mimo, dbode_mimo, mimo_frequency_response, mimo_bode_result, &
    design_status_success

real(dp) :: a(2,2), b(2,2), c(2,2), d(2,2)
type(mimo_frequency_plan) :: plan
type(mimo_frequency_response) :: response
type(mimo_bode_result) :: plot

a = reshape([-1.0_dp,0.0_dp,0.0_dp,-2.0_dp],[2,2])
b = reshape([1.0_dp,0.0_dp,0.0_dp,1.0_dp],[2,2])
c = b
d = 0.0_dp
plan = prepare_mimo_frequency(a,b,c,d)
if (plan%status /= design_status_success) error stop 'invalid system'

! Output 2, input 1 only; size(response%response) is [1,1,50].
response = freqresp_mimo(plan,n=50,inputs=[1],outputs=[2])
! The same plan supports discrete-time Bode responses.
plot = dbode_mimo(plan,n=80,dt=0.1_dp,outputs=[1,2])
```

These plan calls use the same status codes, singular-frequency handling,
frequency units, and default grids as ordinary MIMO calls. Analog automatic
grid bounds are retained from the original state matrix, even though the
internally transformed matrix is different. A plan retains either a unitary
Hessenberg or a unitary complex-Schur reduction and the transformed input/output
matrices; no frequency-dependent factorization is cached. Per-frequency solves
are still performed as needed. The optional Schur path requires a larger
one-time QR setup but triangular solves need no pivoting. Unlike an eigenvector
basis, both choices use unitary similarity transformations.
The plan is particularly useful when a model is queried many times with
short frequency grids, different channel subsets, or both. For a single
sweep, `solver='auto'` remains the recommended default. In one GNU Fortran
`-O3` benchmark (`benchmarks/benchmark_signal_mimo_plan.f90`), an 80-state
model with two selected inputs, two selected outputs, and 25 sweeps of 20
frequencies took 0.0467 CPU seconds with separate reductions versus 0.0130
seconds using one prepared plan (3.6x faster, including preparation). This
is machine- and compiler-dependent, not a guaranteed speedup.

The underlying frequency-domain equation remains `(i*w*I-A)*X=B`
(analog) or `(exp(i*w)*I-A)*X=B` (digital). Both paths avoid
characteristic-polynomial expansion. The direct path requires roughly
`O(nstate**3)` work per frequency; Hessenberg reduction requires
`O(nstate**3)` preparation once and `O(nstate**2*(1+ninput))` work per
frequency, plus output multiplication. The similarity transform is unitary
and does not require computing potentially ill-conditioned eigenvectors.
With `solver='schur'`, the prepared reduction is a complex Schur form;
with the default solver it is Hessenberg. For `nstate=0`, the response is
exactly `D`. All methods store
`noutput*ninput*nfreq` complex response values. Large sweeps should use
`n` or explicit `w` to control that storage.

On one machine with GNU Fortran 14.2 at `-O3`, the benchmark
`benchmarks/benchmark_signal_mimo_hessenberg.f90` measured 0.174 CPU
seconds for dense evaluation versus 0.015 seconds for Hessenberg on an
80-state, 600-frequency, two-output/three-input case (about 11.6x faster).
Performance depends on matrix structure, compiler, and frequency count.
Use `solver='direct'` to compare results or avoid setup cost for short grids.

In a separate 48-state, 1,500-frequency, four-sweep GNU Fortran `-O3`
benchmark (`benchmarks/benchmark_signal_mimo_schur.f90`), Schur preparation
cost 0.00234 CPU seconds versus 0.00037 for Hessenberg; afterward, Schur
sweeps took 0.0376 seconds versus 0.0578 for Hessenberg (about 1.54x faster
in evaluation). These timings are not general guarantees and should be measured
on the intended workload. An independent seeded comparison against SciPy 1.17
covered 150 real/complex MIMO systems (4,140 channel values) with maximum
scaled discrepancy `5.69e-14`, and all 150 Schur plans converged.

Invalid dimensions, nonfinite coefficients/frequencies, invalid sample
intervals, and nonpositive automatic `n` return
`design_status_invalid_input` (and do not allocate result arrays).
A singular/nonfinite frequency evaluation yields `design_status_singular`
with NaNs **for that frequency slice only**, preserving nonsingular slices
for inspection. The Bode result likewise preserves valid frequencies and
sets invalid slices to NaN; phase unwrapping restarts after a singularity.
Zero response channels have magnitude `-Inf` dB. Singular shifted matrices
are reported even if their transfer-function singularities would cancel.
Near-singular systems remain ill-conditioned in double precision.

Validation (GNU Fortran and SciPy 1.17):

```text
python tools/signal_reference/run_signal_suite.py
python tools/signal_reference/compare_signal_lti_mimo_scipy.py
python tools/signal_reference/compare_signal_lti_mimo_channels_scipy.py
# Optional: compile benchmarks/benchmark_signal_mimo_plan.f90 with module objects.
# Optional: build and run the Hessenberg benchmark with module objects.
```

The deterministic Hessenberg regression verifies both numerical paths,
including singular points, automatic solver selection, real/complex MIMO
models, and Bode magnitudes/phases. The baseline regression includes a
60-state system for which characteristic-polynomial coefficients overflow,
complex and static systems,
multiple inputs/outputs, time-unit conversion, phase unwrapping, invalid
inputs, and singular-frequency recovery. The independent seeded comparison
covers 480 MIMO systems across 3,840 frequency evaluations and 22,648
input-output channel responses, using `scipy.linalg.solve` as the numerical
reference rather than SciPy's SISO `freqresp` convenience function.
The separate channel-selection comparison tests 280 additional real/complex
systems, analog/digital plan overloads, reordered and repeated channel
indices, and 18,753 selected channel-frequency values. The deterministic
plan regression also covers both Bode variants, invalid selectors, static
feedthrough, singular frequencies, and automatic-grid parity.


## Polynomial deconvolution and FFT-based signal envelopes

The new `deconvolve(signal,divisor,quotient,remainder,status)` is a subroutine
implementing polynomial long division of sequences in convolution order.
It accepts rank-one `real(dp)` or `complex(dp)` signals and divisors in all
four combinations; mixed input types yield complex outputs. Quotient and
remainder are allocatable output arrays. The quotient has length
`max(0, size(signal)-size(divisor)+1)` and the remainder has the original
signal length. When the divisor is longer than the signal, the quotient is
empty and the remainder equals the signal. Empty/nonfinite inputs and a zero
leading divisor coefficient return `deconvolution_status_invalid_input` with
empty quotient/remainder arrays; valid calls return
`deconvolution_status_success`.

`envelope(z,bp_low,bp_high,n_out,squared,residual,axis,remainder,status)` is an
array-returning generic for real/complex input arrays of ranks **1 through 3**.
It computes an FFT-domain band-selected analytic/baseband envelope. The result
is always real and has the shape of the input, except that the selected axis
may be resized by `n_out`. `bp_low` and `bp_high` are **signed integer FFT-bin
indices**, with inclusive lower and exclusive upper bounds and the constraint
`-(n/2) <= bp_low < bp_high <= (n+1)/2`, where `n` is the signal length.
Defaults are 1 and `(n+1)/2`; unlike SciPy's `bp_in=(low,high)` tuple,
Fortran uses separate scalar arguments. To include negative frequencies, pass
an explicit negative `bp_low`.

The optional `axis=` is **one-based** and defaults to the last dimension;
`n_out=` Fourier-resamples the selected axis and `squared=.true.` returns the
squared envelope. The optional `remainder=` output is allocatable and real or
complex to match the input. `residual='lowpass'` (default) returns the rejected
frequencies below the selected band, `residual='all'` returns all rejected
frequencies, and `residual='none'` returns no residual (an empty allocated
array if `remainder` is nevertheless supplied). Invalid parameters or
nonfinite input produce `envelope_status_invalid_input`; valid input produces
`envelope_status_success`. SciPy stacks the envelope and residual on an extra
axis, but SciFort uses the separate optional argument to retain ordinary
Fortran-shaped arrays.

```fortran
use scifort_signal, only : dp, envelope, deconvolve
real(dp) :: samples(128)
real(dp), allocatable :: amplitude(:), drift(:), quotient(:), remainder(:)
integer :: status

amplitude = envelope(samples, bp_low=2, residual='lowpass', &
    remainder=drift, status=status)
call deconvolve([1.0_dp,3.0_dp,3.0_dp,1.0_dp], &
    [1.0_dp,1.0_dp], quotient, remainder, status)
```

The implementation uses the existing dependency-free SciFort FFT, including
its mixed-radix, radix-2, and Bluestein paths for larger transforms.

`envelope_into(z, y, bp_low, bp_high, n_out, squared, residual, axis,
remainder, status)` adds a caller-output form for **real or complex rank 1-15**
arrays. `y` is always real, `remainder` (when provided) has the same kind as
`z`, and both have the same rank as `z`. All non-sample axes retain their
input extents; the selected axis has length `n_out` (or the original length).
The sample axis defaults to the last and also accepts negative indexing
(`axis=-1` selects the last dimension). The bandpass and residual semantics
match those of the rank-1/2/3 `envelope` generic above.

The caller must allocate the output and any residual array at the required
shape. A successful `residual='none'` call leaves an explicitly supplied
`remainder` **unchanged**; when the optional argument is absent, residual
FFT work is skipped. If any input sample is nonfinite or any validation fails,
**neither output is modified**, even if earlier channels succeeded.
Noncontiguous input/output sections are supported. This separate-output
caller-owned form is a **SciFort extension** to SciPy's stacked return convention.

```fortran
use scifort_kinds, only : dp
use scifort_signal, only : envelope_into
real(dp) :: x(4,128,3), amplitude(4,64,3), drift(4,64,3)
integer :: status

call envelope_into(x, amplitude, axis=2, n_out=64, bp_low=2, &
    residual='lowpass', remainder=drift, status=status)
```

`test/test_signal_envelope_rank15_into.f90` checks both scalar types, rank 15,
negative axes, strided sections, empty batches, and transactional error
handling. `reference/compare_signal_envelope_nd_scipy.py` compares rank-four
outputs with independent SciPy 1.17 calculations (360 seeded cases). The
`reference/run_checks.py` script builds the standalone comparison drivers and
reproduces 1,900 seeded SciPy 1.17 cases, while `test/test_signal_envelope_deconvolve*.f90`
provide native regression tests and public API smoke tests.

The complete GNU Fortran 14.2 standalone signal suite passes 74/74 tests with
runtime checking and 74/74 at `-O3`. The targeted regression passes with
floating-point traps as well. The project-wide FPM test remains unverified
because FPM was not installed in this environment.

## DPSS/Slepian windows and multitaper power spectra

`dpss(m,nw,kmax,sym,norm,solver,ratio_method)` generates Discrete Prolate
Spheroidal Sequences using self-contained symmetric-tridiagonal eigenproblems. `nw` is the
positive time-half-bandwidth product (strictly smaller than `m/2`).
It returns a `dpss_result` with `status`, `tapers(m,k)` and `ratios(k)`.
**Fortran layout:** the first index is the sample and the second is the
window, whereas SciPy's explicit-`Kmax` result is `(Kmax,M)`. The real-valued
windows follow SciPy's sign conventions (symmetric windows have positive sum;
antisymmetric windows have an initially positive main lobe).

- `kmax` defaults to one; if omitted, normalization defaults to
  `'approximate'` (the SciPy single-window default). If explicitly supplied,
  normalization defaults to `'l2'`, matching SciPy's multitaper default.
- `sym=.true.` (default) requests symmetric filter-design windows;
  `.false.` requests periodic spectral-analysis windows.
- `norm='l2'`, `'approximate'`, and `'subsample'` are supported.
  The `subsample` mode uses the exact periodic Fourier-interpolation peak at
  the half-sample center of an even-length first taper, instead of the
  empirical even-length factor used by `approximate`. The correction is
  applied to every requested taper with the same scale factor, as in SciPy.
  Odd internal lengths use simple maximum normalization. It is evaluated
  using an O(m) Dirichlet-kernel identity with no additional FFT dependency.
- Concentration ratios are measured on the eigenvectors *before* any
  periodic trimming and peak renormalization, as in SciPy. Optional
  `ratio_method='auto'` (default), `'direct'` or `'fft'` is a SciFort
  extension. `auto` uses the FFT implementation when the internal
  window length is at least 256, otherwise the original direct sum.
- Optional `solver='auto'` uses the selected tridiagonal eigensolver for
  internal lengths at least 256 when at most 32 tapers are requested;
  `solver='ql'` forces the full-eigenvector method and `solver='partial'`
  forces selected eigenpairs. Eigenvector memory is O(m*k) for the partial
  method versus O(m^2) for QL. Both solver modes can use FFT ratios.
- FFT-based concentration ratios evaluate the *same* sinc Toeplitz
  Rayleigh quotient as the direct method, using zero-padded autocorrelation.
  Their work is O(k*m*log(m)) instead of O(k*m^2), with O(m) complex
  scratch space. Both ratio methods use O(m) storage per taper.

```fortran
use scifort_signal, only : dp, dpss, dpss_result, &
    multitaper_psd, spectral_result
real(dp) :: x(256)
type(dpss_result) :: win
type(spectral_result) :: psd

! Supply x from your sampled signal.
win = dpss(256, 3.5_dp, kmax=5, sym=.false.)
! Force the scalable eigensolver and FFT-based concentration calculation.
win = dpss(8192, 3.5_dp, kmax=5, solver='partial', ratio_method='fft')
! A peak-normalized symmetric design taper with half-sample correction.
win = dpss(64, 3.5_dp, norm='subsample')
psd = multitaper_psd(x, 3.5_dp, kmax=5, fs=100.0_dp, &
    nfft=512, weighting='eigen')
```

`multitaper_psd(x,nw,kmax,fs,nfft,detrend,weighting)` is a **SciFort extension**,
not a direct `scipy.signal` function. It returns `spectral_result`
(`frequencies`, `power`, and `status`) for finite, one-dimensional real input.
It generates periodic unit-energy DPSS windows, averages one-sided
`density` periodograms, and doubles all interior positive-frequency bins.
`detrend=.true.` by default subtracts the raw arithmetic mean before
windowing. `weighting` is `'equal'` (default), `'eigen'` (weighted by DPSS
concentration ratios), or `'adaptive'` (Thomson's frequency-dependent
weights). The default `kmax` is `max(1,int(2*nw)-1)`, capped
at the signal length. `nfft` defaults to signal length, and `fs` to 1.
The original real 1-D behavior remains unchanged. Complex one-dimensional
signals use the same name and return two-sided spectra in unshifted FFT
frequency order (nonnegative frequencies followed by negative frequencies).

For rank-2 real or complex arrays, `multitaper_psd` returns
`multitaper_channels_result` with `power(nfrequency,nchannels)`; by default
samples are along the first axis, or `axis=2` may be specified. The separate
`multitaper_spectral_matrix` function returns `multitaper_matrix_result` with
`power(nfrequency,nchannels,nchannels)` and the convention
`Pxy = conjg(FFT(x))*FFT(y)`. Every frequency slice is Hermitian and positive
semidefinite to roundoff, and its diagonal equals the channelwise PSD.
Real inputs have one-sided, interior-bin-doubled spectra; complex inputs
have two-sided spectra.

`multitaper_csd(x,y,nw,...)` estimates the cross-spectral density of two
equal-length signals, supporting real, complex or mixed real/complex inputs;
when either input is complex its frequency grid is two-sided.
`multitaper_coherence(x,y,nw,...)` estimates magnitude-squared coherence
for two real or two complex signals. Its values are between zero and one;
zero-power frequency bins return zero instead of an undefined ratio.
All spectral estimators support mean removal (`detrend`, default true),
`fs`, `nfft` (at least the input length), `kmax`, and `'equal'`,
`'eigen'`, or `'adaptive'` weighting. Invalid shapes/parameters return nonzero `status` rather than
stopping. All of these multitaper estimators are **SciFort extensions**, not
standard `scipy.signal` entry points.

```fortran
real(dp) :: samples(256,3)
complex(dp) :: quadrature(256,3)
type(multitaper_channels_result) :: psd_by_channel
type(multitaper_matrix_result) :: cross_matrix
type(cross_spectral_result) :: pxy

! Fill the sampled arrays first.
psd_by_channel = multitaper_psd(samples,3.5_dp,kmax=5,fs=100.0_dp)
cross_matrix = multitaper_spectral_matrix(quadrature,3.5_dp, &
    kmax=5,nfft=512,weighting='eigen')
pxy = multitaper_csd(samples(:,1),samples(:,2),3.5_dp)
! cross_matrix%power(:,1,2) stores channel 1 vs. channel 2.
```

`weighting='adaptive'` uses the DPSS concentration ratios and an iterative
Thomson-style estimate of the frequency-specific broadband-leakage bias. The
noise floor is estimated as the centered input mean-square divided by `fs`.
Adaptive weights are computed separately for each channel. Cross-spectra use
geometric-mean amplitude weights, preserving positive semidefiniteness and
the Hermitian property. With fewer than three tapers, adaptive weighting
falls back to concentration weighting; `equal` and `eigen` preserve their
previous numerical behavior. This is a **SciFort extension** motivated by
published multitaper methodology, not a claim that SciPy has the same API.

`multitaper_psd_ci(x,nw,kmax,fs,nfft,detrend,weighting,confidence)` computes
a real or complex **one-dimensional** PSD and pointwise log-jackknife
confidence intervals. The result type `multitaper_confidence_result` contains
`frequencies`, `power`, `lower`, `upper`, `stderr_log`, `equivalent_dof`, and
`status`. For each omitted taper, adaptive weighting is re-estimated, and the
sample variance of the delete-one **log** spectra is scaled by `(K-1)/K`.
The interval is `PSD * exp(+-t_critical * stderr_log)`, where the critical
value is the two-sided Student-t quantile with `K-1` degrees of freedom.
`confidence` defaults to `0.95` and must be strictly between zero and one.
The returned approximate spectral equivalent degrees of freedom at a
frequency is `2/sum(q_k**2)`, with normalized squared taper weights `q_k`.
At least **two** tapers are required for jackknife intervals. For a signal
with exactly zero power at a frequency, its PSD and interval are zero.

```fortran
use scifort_signal, only : dp, multitaper_psd_ci, multitaper_confidence_result
type(multitaper_confidence_result) :: estimate
real(dp) :: x(256)
! Fill x with observations first.
estimate = multitaper_psd_ci(x,3.5_dp,kmax=5,fs=100.0_dp, &
    weighting='adaptive',confidence=0.95_dp)
! estimate%lower and estimate%upper are in PSD units.
```

These confidence limits are **pointwise**, not simultaneous bands. Jackknife
tapers are not independent observations, so their coverage is approximate;
K=2 gives especially wide and unstable limits. The estimate does not account
for uncertainty from estimating `nw` or choosing frequencies. Exact
zero-power bins have zero-width intervals by convention.

The spectral-matrix calculation stores the complete `nchannels` squared
array for each frequency. For many channels, this needs proportionately
more memory than the channelwise `multitaper_psd` implementation, which
computes only its diagonal spectra. The QL DPSS option stores a dense eigenvector matrix;
the selected-eigenpair solver (`solver='partial'` or suitable `auto` calls)
does not.

For numerical validation, with SciPy and NumPy installed:

```shell
python test/compare_signal_dpss_scipy.py
python test/compare_signal_multitaper.py
python test/compare_signal_multitaper_channels.py
python test/compare_signal_multitaper_adaptive.py
```

The DPSS test checks 344 symmetric/periodic, normalized, and randomized cases
against `scipy.signal.windows.dpss`; the multitaper test uses SciPy tapers
and independent NumPy FFT power accumulation in 240 randomized cases.
The new multichannel comparison includes 160 randomized real/complex cases,
non-power-of-two FFT sizes, selected time axes, and Hermitian/positive-
semidefinite spectral matrix checks. The adaptive comparison covers 140
randomized real/complex matrices and one-dimensional confidence intervals,
using SciPy DPSS tapers, NumPy FFTs and Student-t quantiles as an
independent oracle.

## Legacy two-sided spectra for real input (2026-10-27)

`periodogram`, `welch`, `csd`, `stft`, and `spectrogram` now accept the
trailing optional logical `return_onesided` in both the rank-one and
caller-output assumed-rank interfaces. The default remains `.true.` for
real input; `.false.` exposes the entire FFT frequency grid, preserving
negative-frequency bins and avoiding the one-sided power factor of two.
For complex inputs the spectra remain two-sided regardless of the keyword.
The interfaces preserve all pre-existing positional arguments.

```fortran
use scifort_signal, only : dp, welch, spectral_result
real(dp) :: samples(128)
type(spectral_result) :: spectrum

spectrum = welch(samples, nperseg=32, noverlap=16, &
    window='hann', return_onesided=.false.)
! size(spectrum%frequencies) == 32, with signed FFT-order frequencies.
```

The reproducible checked/optimized GNU Fortran test suite is:

```sh
python tools/run_signal_twosided.py
```

It includes comparisons against 192 independent SciPy 1.17 reference
configurations and 7,248 spectral coefficients, plus native tests for full
frequency parity, PSD symmetry, strided views, high-rank inputs, and
all-or-nothing error handling. SciPy/NumPy are optional test-only
requirements; no extra Fortran library is required.


## Conditional coherence from multichannel cross-spectral matrices

The SciFort extension `multitaper_partial_coherence` computes squared
partial coherence across pairs of real or complex channels, conditioning on
all remaining channels. `partial_coherence_matrix` accepts a previously
computed Hermitian cross-spectral-density matrix. The implementation uses
complex Cholesky solves with optional scale-invariant diagonal loading,
returns spectral precision matrices, and marks singular FFT frequencies.
See [Partial coherence](signal_partial_coherence.md) for the full API and
limitations.

The additional `multitaper_tf_partial_coherence` and
`partial_coherence_tf_matrix` interfaces compute conditional coherence at
each time-frequency location from moving-window multitaper CSD matrices.
They return four-dimensional conditional-coherence/precision arrays and
independent `(frequency,time)` validity masks. See
[Time-frequency partial coherence](signal_tf_partial_coherence.md).
Both partial-coherence result types optionally include complex conditional
coherency when `return_coherency=.true.`; diagonal-target shrinkage is
available via `shrinkage=` and is distinct from trace-ridge loading.
The opt-in `auto_shrinkage=.true.` estimates a Gaussian-Wishart diagonal
shrinkage coefficient with a fixed effective number of tapers (equal/eigen
weighting only). Separately, `cv_shrinkage=.true.` recomputes taper weights
in leave-one-out Fourier-spectrum fits, including Thomson adaptive weights,
and selects per-frequency (or per-time-frequency) predictive shrinkage. The
latter requires at least three tapers and should not be interpreted as a
calibrated confidence procedure. See
[cross-validated diagonal shrinkage](signal_cv_shrinkage.md).


## Conditional coherence between channel groups

`group_coherence_matrix`, `group_coherence_tf_matrix` and the direct
`multitaper_group_coherence` variants estimate the determinant coherence of
two groups of spectral channels, optionally conditioning on a third group.
Band-integrated versions compute the statistic after integrating the CSD.
This is a SciFort extension, not a `scipy.signal` function. See
[Group spectral coherence](signal_group_coherence.md) for definitions,
examples, limitations, and native/NumPy tests.
