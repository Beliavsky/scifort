# Exact-band fixed-weight multitaper cross-spectral matrices

**SciFort extension, not a `scipy.signal` function.**
`fixed_taper_spectral_matrix_spectrogram` constructs time-resolved **complex
cross-spectral matrices** from a multichannel real or complex time series
using a caller-supplied deterministic bank of tapers. It evaluates arbitrary
point frequencies by direct finite Fourier sums and integrates entire
frequency bands **analytically**, rather than by summing FFT bins or applying
trapezoidal quadrature. The resulting spectrum is compatible with the
existing `canonical_moving_block_bootstrap` API.

## Native Fortran example

```fortran
use scifort_signal, only : dp, dpss, dpss_result, &
    fixed_taper_spectral_matrix_spectrogram, &
    multitaper_spectral_matrix_spectrogram_result, &
    canonical_moving_block_bootstrap, &
    canonical_moving_block_bootstrap_result, spectral_status_success

type(dpss_result) :: bank
type(multitaper_spectral_matrix_spectrogram_result) :: spectra
type(canonical_moving_block_bootstrap_result) :: inference
real(dp) :: x(640,5), bands(2,2)

! Fill x: samples first, channels second. Here 40 complete windows of 16.
bank = dpss(16,2.5_dp,kmax=3,sym=.false.)
if (bank%status /= 0) error stop 'DPSS generation failed'
bands(:,1) = [0.0_dp,0.125_dp]
bands(:,2) = [0.125_dp,0.5_dp]

spectra = fixed_taper_spectral_matrix_spectrogram( &
    x,bank%tapers,[0.1375_dp],bands=bands, &
    taper_weights=bank%ratios,normalized=.true.,average=.true.)
if (spectra%status /= spectral_status_success) error stop 'bad CSD input'

! The resulting three bins are one point frequency and two EXACT band means.
inference = canonical_moving_block_bootstrap(spectra,[1,2],[3,4], &
    conditioned_on=[5],block_length=4,n_resamples=399, &
    circular=.true.,seed=321,interval_method='percentile')
```

The public generic accepts both `real(dp) x(:,:)` and
`complex(dp) x(:,:)`. The returned existing
`multitaper_spectral_matrix_spectrogram_result` contains
`power(nfrequency+n_bands,n_windows,n_channels,n_channels)`, frequencies
(point frequencies followed by band midpoints), window-center times **in
sample units**, `effective_tapers`, and `status`. No C ABI or external
numerical-library dependencies were added.

## Inputs and conventions

- The first dimension of `x` is time; the second is channel. The window
  length `L` is determined by `size(tapers,1)` and the taper count `K` by
  `size(tapers,2)`; both must be positive. `size(x,1) >= L`.
- `frequencies(:)` is an arbitrary finite set of point frequencies **in
  cycles/sample**; negative and out-of-Nyquist values are accepted (Fourier
  aliases). `bands(2,nband)` lists finite, positive-width bands, in the same
  units. Either array may have zero entries, but not both. Band order and
  overlap are unrestricted.
- `taper_weights(:)` is fixed independently of time and frequency;
  defaults to all ones. Relative nonnegative weights are normalized to sum
  to one, allowing exact zeros, but not an all-zero array. All selected
  tapers must have finite, positive squared energy; an unused zero-energy
  taper is allowed. All data, tapers, frequencies and weights must be finite.
- `normalized=.true.` (default) divides each taper's outer product by
  **that taper's own squared energy**, including for nonorthogonal/unequal-
  norm banks. `normalized=.false.` returns the weighted average of unscaled
  outer products.
- `average=.true.` (default) computes each band-**mean** spectral matrix.
  `average=.false.` computes the band **integral**, retaining a factor of
  bandwidth. Neither option produces a one-sided doubling or applies an
  `fs` frequency-unit conversion. Detrending is the caller's responsibility.
- `noverlap=0` (default) creates complete nonoverlapping windows. Optional
  `0 <= noverlap < L` shifts windows by `L-noverlap` samples. Trailing
  incomplete windows are dropped. **Use `noverlap=0` when passing this
  output to the moving-block bootstrap**; overlapping original windows do
  not meet the bootstrap's spectral-window assumption.
- The complex cross-spectrum convention is `F_i(f)*conjg(F_j(f))`, matching
  the SciFort exact finite-window VAR(p) population target. This is the
  conjugate of the `conj(F_i)*F_j` convention of `scipy.signal.csd`; users
  comparing phases must account for this difference.

For weights `eta_k` and tapers `w(t,k)`, define

\[
Z_{wk,i}(f)=\sum_{t=0}^{L-1} w(t,k)x_{w,t,i}e^{-2\pi ift},\quad
S_w(f)=\sum_{k=1}^{K}\frac{\eta_k}{\sum_j\eta_j}
\frac{Z_{wk}(f)Z_{wk}(f)^H}{d_k},
\]

where `d_k=sum_t w(t,k)^2` if `normalized` and one otherwise. For a band
`[a,b]`, the exact **mean** kernel is

\[
K_{ts}(a,b)=\operatorname{sinc}\big((b-a)(t-s)\big)
 e^{-\pi i(a+b)(t-s)},
\quad \operatorname{sinc}(u)=\frac{\sin(\pi u)}{\pi u}.
\]

Multiplying this by `(b-a)` gives the integral kernel. The Fourier sum and
kernel retain off-diagonal complex phase and Hermitian positive
semidefiniteness (up to rounding). Reported `effective_tapers` equals
`1/sum_k(eta_k/sum_j eta_j)^2`, a fixed-weight equivalent taper count; it
is **not** a claim that finite-window DPSS periodograms are independent.
The kernel costs `O(K * n_windows * n_bands * L**2 * n_channels)` for its
matrix product plus channel operations, and is more expensive than FFT-bin
integration for very long windows. It is useful when exact finite-window
frequency-band targets are essential.

## Verification

```bash
python tools/run_signal_varp_multitaper_coverage.py
python tools/check_signal_varp_population_public.py
```

The first runner builds checked and optimized GNU Fortran, tests invalid
inputs/real-complex overloads/window overlap, and compares 60 random
scenarios (788 complex CSD matrices per build) with independent NumPy
finite Fourier sums and **96-node Gauss-Legendre quadrature** using SciPy
DPSS or arbitrary random tapers. It also runs the VAR(p) conditional
coverage smoke tests. All 60 cases pass in each build with observed maximum
complex entry error about `3.27e-13`. The public closure compiles 106
modules and passes 25 signal API regressions, including the new public
exact-band-to-bootstrap test.

This function **does not implement adaptive/data-dependent taper weights**;
its expectation is directly comparable to the deterministic bank
`varp_expected_multitaper_spectrum` and
`varp_expected_band_multitaper_spectrum` formulas only under fixed weights.
It does not by itself supply uncertainty intervals or calibrated inference.
