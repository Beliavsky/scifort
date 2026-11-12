# Window-jackknife uncertainty for canonical spectral coherence

**SciFort extension**; this is not a public `scipy.signal` function.

For temporally dependent windows, see the complementary
[equal-block jackknife](signal_canonical_block_jackknife.md). The window
entry points also accept an optional `block_length=` (default `1`); named
block entry points make this choice explicit. Blocks group contiguous,
nonoverlapping windows and require at least four equal-size deletion groups.


`canonical_window_jackknife` estimates **pointwise** variability in the
squared canonical spectral coherence modes between two disjoint channel groups
(optionally conditioning on another disjoint group). Instead of pretending
that neighboring Fourier bins or adaptively weighted DPSS tapers are independent,
it deletes **entire spectral-estimation windows**. The method is suitable for
reasonably stationary signal segments whose **nonoverlapping windows can be
treated as approximately independent sampling units**. It does not establish
that assumption from the data.

The input `multitaper_spectral_matrix_spectrogram_result` stores a complex
Hermitian cross-spectral matrix at `(frequency, window, channel, channel)`. Its
`times` array must be finite and strictly increasing. Provide **at least four
windows**; more independent windows are preferable. The signal-level generic
`multitaper_canonical_window_jackknife` always uses `noverlap=0` and returns
one spectral-window jackknife result. Its real and complex overloads support
the multitaper `equal`, `eigen`, and `adaptive` weighting modes.

## Fortran examples

```fortran
use scifort_signal, only : dp, canonical_window_jackknife_result, &
    multitaper_canonical_window_jackknife

type(canonical_window_jackknife_result) :: result
real(dp) :: signals(4096,5)

! Populate signals(:,1:5) before calling.
result = multitaper_canonical_window_jackknife(signals,3.5_dp,256, &
    [1,2],[3,4],conditioned_on=[5],kmax=6, &
    weighting='adaptive',ridge=0.05_dp,shrinkage=0.1_dp)

if (result%status == 0) then
    ! result%squared_modes(f,mode): full-data estimate
    ! result%stderr(f,mode): delete-one-window jackknife SE
    ! result%lower95(f,mode), result%upper95(f,mode): approximate pointwise limits
    ! result%valid(f,mode): whether uncertainty for an isolated mode is meaningful
end if
```

For a precomputed collection of nonoverlapping spectral windows (including
results from the public `multitaper_spectral_matrix_spectrogram` function):

```fortran
result = canonical_window_jackknife(window_spectra, [1,2], [3,4], &
    conditioned_on=[5],ridge=0.05_dp,shrinkage=0.1_dp)
```

For exact piecewise-linear **band integration before the canonical analysis**:

```fortran
result = canonical_window_jackknife(window_spectra, [1,2], [3,4], &
    conditioned_on=[5],ridge=0.05_dp, &
    bands=reshape([8.0_dp,13.0_dp,18.0_dp,30.0_dp],[2,2]))
```

The array `bands(2,nbands)` contains frequency intervals in the units of the
supplied frequency grid. Passing `average=.true.` divides the integrated
matrix by its band width; fixed relative ridge and shrinkage make the
canonical squared coherences invariant to this positive scalar division.

## Estimator and returned quantities

For each frequency or band, let the independent window matrices be
\(S_1,\ldots,S_W\) and let \(g(S)\) compute the ordered canonical squared
coherences of the **regularized** matrix. The point estimate is

\[
\hat\theta=g\!\left(\frac{1}{W}\sum_{w=1}^{W}S_w\right).
\]

For each held-out window \(w\), refit the same estimator to the other \(W-1\)
matrices, obtaining \(\hat\theta_{(-w)}\). The ordinary delete-one jackknife
standard error is

\[
\widehat{\mathrm{SE}}(\hat\theta)=
\sqrt{\frac{W-1}{W}\sum_{w=1}^W
    (\hat\theta_{(-w)}-\overline{\theta}_{(-\cdot)})^2}.
\]

The estimated bias is
\((W-1)(\overline{\theta}_{(-\cdot)}-\hat\theta)\), also returned. There
is no automatic correction of the point estimate for this bias.

An optional interpretation-oriented approximate **pointwise** 95% interval
uses \(z=\operatorname{atanh}(\sqrt{\theta})\):

\[
\theta_{\rm lower}=
\tanh^2\bigl(\max(0,z-1.959964\,\widehat{\mathrm{SE}}_z)\bigr),
\qquad
\theta_{\rm upper}=
\tanh^2\bigl(z+1.959964\,\widehat{\mathrm{SE}}_z\bigr).
\]

Here \(\widehat{\mathrm{SE}}_z\) is the **same delete-one jackknife formula**
applied to each refitted Fisher-transformed root. The transformation keeps
the interval in \([0,1]\). It is a **normal approximation, not an exact or
calibrated 95% confidence procedure**; small numbers of windows, weak
coherence, strong spectral dependence, and root multiplicity can all spoil
coverage. For values extremely close to one the transform clips the square
to \(1-10^{-12}\) to keep its arithmetic finite.

Returned arrays (first dimension always frequency or band):

| Field | Meaning |
| --- | --- |
| `frequencies` | Spectral frequencies or integrated band centers |
| `squared_modes(f,k)` | Descending squared canonical coherences from the full window mean |
| `stderr(f,k)` | Ordinary jackknife SE on the squared-root scale |
| `stderr_fisher(f,k)` | Jackknife SE on the Fisher scale |
| `lower95(f,k)`, `upper95(f,k)` | Approximate normal/Fisher pointwise 95% interval |
| `bias(f,k)` | Jackknife bias estimate on squared-root scale |
| `valid(f,k)` | Whether that mode's inference met all validity checks |
| `coherence(f)` | Aggregate group coherence of the window mean |
| `coherence_stderr(f)` | Ordinary jackknife SE of that aggregate |
| `coherence_valid(f)` | Every leave-one fit was valid for that frequency/band |
| `n_windows` | Number of supplied nonoverlapping resampling units |
| `status` | 0 success, 1 invalid input, 2 singular/unusable spectral cell |

`gap_tolerance` (default `1e-6`) controls rejection of individual modes with
squared strengths within that distance of either 0 or 1, or too close to
an adjacent canonical root. **All full and delete-one fits must pass** this
check for a mode to have `valid=.true.`. Rejected mode standard errors,
biases, and interval limits are quiet NaNs, not spuriously zero. The point
estimates and aggregate coherence remain available wherever the underlying
canonical calculation is valid. An aggregate SE can still be valid for tied
individual canonical roots.

## Assumptions and limitations

- **Window independence is an assumption.** `noverlap=0` avoids literal shared
  samples, but does not eliminate temporal dependence in a colored or
  nonstationary process. Strong dependence calls for larger independent
  blocks, explicit block resampling, or a different statistical model.
- **Do not jackknife neighboring frequency bins.** Each Fourier band is
  aggregated *inside each window first*, preserving within-window spectral
  dependence. Reported intervals are pointwise and do not account for
  searching over many frequencies or bands (no simultaneous coverage).
- **Adaptive tapers are allowed within windows.** Each window's complete
  adaptively weighted CSD is a sampling unit. The method does not treat
  individual adaptive taper weights as independent or recompute them when
  a whole window is removed. Dependence between windows remains possible.
- **Only fixed `shrinkage` and `ridge` are supported.** Both are held constant
  across refits. Automatic Gaussian-Wishart and cross-validated shrinkage
  coefficients are not transported across windows or reused as though they
  were known fixed parameters. If you tune regularization using the same
  windows, the reported jackknife omits that selection uncertainty.
- **Canonical modes must be distinguishable.** Close roots, zero and perfect
  roots make ordinary rank-by-rank jackknife inference nonregular. No inference
  for canonical coefficient directions or tracked subspaces is provided.
- **A window is the unit of resampling, not a taper.** Increasing `kmax` may
  improve each window's CSD estimate, but does *not* increase `n_windows`.
  At least four windows are required by the API; considerably more are
  advisable for normal-approximation intervals.

## Reproducibility

```bash
python tools/run_signal_canonical_jackknife.py
python tools/check_signal_canonical_jackknife_public.py
```

The native driver exercises invalid inputs, tied/zero modes, real/complex
signals, adaptive weighting, and band-integrated matrices. The independent
NumPy/SciPy oracle uses complex Cholesky and SVD, separate leave-one window
means, and linear interpolation/trapezoidal band integration. The scripts
require GNU Fortran and Python with NumPy and SciPy; the Fortran library
requires no new dependencies.
