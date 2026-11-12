# Time-frequency partial coherence (SciFort extension)

`multitaper_tf_partial_coherence` estimates **squared partial coherence**
between multichannel signals, conditional on all other measured channels, in
overlapping time windows. This is a SciFort extension, **not a public
`scipy.signal` function**. It builds on the existing multitaper cross-spectral
matrix spectrogram and partial-coherence numerical kernel.

## Example

```fortran
use scifort_signal, only : dp, partial_coherence_tf_result, &
    multitaper_tf_partial_coherence, partial_coherence_success

real(dp) :: signals(512,4)
type(partial_coherence_tf_result) :: result

! Populate signals(sample, channel).
result = multitaper_tf_partial_coherence(signals, 3.5_dp, 128, &
    noverlap=96, kmax=6, fs=200.0_dp, ridge=0.05_dp)
if (result%status == partial_coherence_success) then
    ! Conditional association of channels 1 and 3, removing channels 2 and 4:
    print *, result%coherence(:,:,1,3)
end if
```

The public generic accepts either `real(dp) x(:,:)` or `complex(dp) x(:,:)`.
Its required inputs are `x, nw, nperseg`, followed by the optional keywords
`noverlap, kmax, fs, nfft, detrend, weighting, axis, detrend_mode, ridge,
shrinkage, return_coherency, auto_shrinkage`.
These arguments have the same meanings and validation rules as
`multitaper_spectral_matrix_spectrogram`, except for the added `ridge`.
The default sample axis is 1; use `axis=2` for `(channel, sample)` data.
The overlap defaults to `nperseg/2`. Real input uses one-sided frequencies;
complex input uses the full FFT frequency ordering. No padding or partial
final windows are used.

## Direct spectral-matrix interface

If the cross-spectral-density matrices have already been estimated, reuse them:

```fortran
use scifort_signal, only : dp, multitaper_spectral_matrix_spectrogram, &
    multitaper_spectral_matrix_spectrogram_result, &
    partial_coherence_tf_matrix, partial_coherence_tf_result

type(multitaper_spectral_matrix_spectrogram_result) :: spectral
type(partial_coherence_tf_result) :: result

spectral = multitaper_spectral_matrix_spectrogram(signals, 3.5_dp, 128, &
    noverlap=96, kmax=6)
result = partial_coherence_tf_matrix(spectral, ridge=0.02_dp)
```

`partial_coherence_tf_matrix` accepts the complex Hermitian array
`power(nfrequency, ntime, nchannel, nchannel)` and checks the supplied
frequency/time coordinates and status. If supplied spectral matrices are
non-Hermitian beyond rounding tolerance or contain nonfinite values, the
result has `partial_coherence_invalid_input` status.

## Returned arrays

`partial_coherence_tf_result` contains:

| Component | Dimensions | Meaning |
|---|---|---|
| `frequencies` | `(nf)` | Fourier frequencies in Hz when `fs` is supplied |
| `times` | `(nt)` | Window-center timestamps in seconds |
| `coherence` | `(nf,nt,nc,nc)` | Squared partial coherence, from 0 to 1 where valid |
| `precision` | `(nf,nt,nc,nc)` | Inverse **regularized** cross-spectral matrix |
| `coherency` | `(nf,nt,nc,nc)` when requested; otherwise empty | Complex partial coherency, including conditional phase |
| `valid` | `(nf,nt)` | Invertibility mask, independently set for every cell |
| `shrinkage_used` | `(nf,nt)` | Coefficient applied at each time-frequency cell |
| `status` | scalar | `partial_coherence_success` (0), `partial_coherence_invalid_input` (1), or `partial_coherence_singular` (2) |

For a Hermitian spectral matrix `S(f,t)`, the algorithm forms

\[
S_\lambda(f,t) = S(f,t) + \lambda\frac{\operatorname{tr}S(f,t)}{n_c}I,
\quad P(f,t) = S_\lambda(f,t)^{-1},
\]

then computes

\[
C_{ij\mid\mathrm{others}}(f,t) =
  \frac{|P_{ij}(f,t)|^2}{P_{ii}(f,t)P_{jj}(f,t)}.
\]

To retain conditional **phase**, request `return_coherency=.true.`. This
computes `coherency(f,t,i,j) = -P(i,j)/sqrt(P(i,i)*P(j,j))` for `i /= j`,
with unit diagonal, so `abs(coherency)**2` matches squared coherence.
The imaginary part is signed partial imaginary coherency, and `atan2` of its
imaginary and real parts yields the phase (undefined for zero magnitude).
Reversing the channel order reverses phase. Unrequested coherency arrays
remain zero-sized to avoid additional memory use. This is a **SciFort
extension**, not a public SciPy function.

The optional `shrinkage=alpha` with `0 <= alpha <= 1` replaces each spectral
matrix by `(1-alpha)*S + alpha*diag(S)` before trace-based `ridge` loading.
This separately preserves the spectral power of each channel, suppresses
cross terms, and can make low-rank matrices invertible when all channel
powers are positive. Defaults are `shrinkage=0`, `ridge=0`. Shrinkage is
biased and is **not** an automatic significance procedure. An optional
`auto_shrinkage=.true.` mode estimates a separate per-cell coefficient under
a complex Gaussian-Wishart approximation; it cannot be combined with
`shrinkage=`. The direct matrix interface requires `effective_tapers=`, while
the signal-level interface supports both equal and eigenvalue-weighted
DPSS tapers using their fixed-weight effective taper count. Adaptive
weights remain excluded from **Gaussian-Wishart automatic shrinkage**.
A separate `cv_shrinkage=.true.` method supports adaptive weighting by
refitting taper weights when each tapered spectrum is held out. See
[automatic diagonal shrinkage](signal_auto_shrinkage.md) and
[cross-validated shrinkage](signal_cv_shrinkage.md) for the different
objectives and limitations.

```fortran
result = multitaper_tf_partial_coherence(signals,3.5_dp,128, &
    noverlap=96,kmax=6,shrinkage=0.10_dp, &
    return_coherency=.true.)
! Complex conditional association across all time-frequency bins:
print *, aimag(result%coherency(:,:,1,3))
```

A Hermitian Cholesky factorization with triangular solves is used separately
for each time-frequency cell. The diagonal of valid coherence matrices is 1.

**Numerical and statistical considerations:**

- `ridge` is nonnegative, dimensionless, and defaults to zero. Positive values
  can stabilize rank-deficient *nonzero* matrices. A zero-power matrix remains
  singular even with ridge loading. Small ridge values can produce unstable
  statistical estimates and should not be treated as automatic tuning.
- If any cell is singular, `status` is `partial_coherence_singular`; only that
  cell's coherence, requested coherency and precision are NaN, while other valid cells are retained.
  Check `valid(f,t)` when processing results in this case.
- The number of independent tapers per window should comfortably exceed the
  number of channels when possible. Local time resolution and spectral
  resolution trade off with each other. All conditional associations are
  descriptive, not causal measures or significance tests.
- Result storage is **O(nf * nt * nc²)**; requesting complex coherency
  additionally allocates an `nf*nt*nc*nc` complex array. It does not keep a
  duplicate complete multitaper spectrogram inside the returned result.

## Numerical verification

```sh
python tools/run_signal_tf_partial_coherence.py
```

The self-contained runner builds checked and targeted optimized GNU Fortran
executables, exercises existing multitaper regressions, and compares
120 independent real/complex and matrix-only cases with NumPy inversions and
SciPy 1.17 DPSS/FFT references. SciPy and NumPy are **test-only** dependencies;
no new runtime Fortran package or system BLAS/LAPACK dependency is required.
The complete FPM suite must be run separately when FPM is available.

The independent complex-coherency/shrinkage regression suite is
`python tools/run_signal_partial_coherency.py` (add `--legacy-references`
to rerun the earlier independent reference suites).
