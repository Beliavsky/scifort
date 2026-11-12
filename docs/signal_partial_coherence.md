# Multichannel partial spectral coherence (SciFort extension)

SciFort exposes squared **partial coherence** between every pair of channels,
conditioned on the remaining channels. This estimates the frequency-domain
association between two signals after the linear contributions of other
observed signals have been accounted for. It is distinct from ordinary
pairwise coherence and is **not a public `scipy.signal` function**.

## Example

```fortran
use scifort_signal, only : dp, multitaper_partial_coherence, &
    partial_coherence_result, partial_coherence_success, partial_coherence_singular

real(dp) :: signals(256,4)
type(partial_coherence_result) :: result

! Fill signals with four time series in (sample,channel) order.
result = multitaper_partial_coherence(signals,3.5_dp, &
    kmax=6,fs=100.0_dp,ridge=0.02_dp)
if (result%status == partial_coherence_success) then
    print *, result%frequencies(2)
    print *, result%coherence(2,1,3) ! Squared coherence of channels 1 and 3,
                                    ! conditioned on channels 2 and 4.
end if
```

The real and complex input overloads accept the same optional keywords as
`multitaper_spectral_matrix`: `kmax`, `fs`, `nfft`, `detrend`, `weighting`,
`axis`, and `detrend_mode`, plus the optional `ridge`. Samples are rows
(`axis=1`) by default. Use `axis=2` when channels are rows. All channel pairs
are estimated on the same Fourier frequency grid and returned in
`result%coherence(nfrequencies,nchannels,nchannels)`; real-valued inputs
produce one-sided spectra, and complex inputs produce two-sided FFT spectra.
`result%precision` has the same shape with complex values. Add
`return_coherency=.true.` to receive the optional complex-valued
`result%coherency` array; omitting it leaves that array empty (size zero).

## Algorithm and interpretation

Let `S(f)` be a Hermitian cross-spectral-density matrix with convention
`S_ij(f) = conj(FFT_i(f))*FFT_j(f)`. For each frequency, compute its inverse
`P(f) = S(f)^(-1)`, also called the spectral precision matrix. For channels
`i` and `j`, the squared partial coherence is

\[
C_{ij\mid \mathrm{others}}(f) =
\frac{|P_{ij}(f)|^2}{P_{ii}(f)P_{jj}(f)}.
\]

The diagonal is one and off-diagonal values lie between zero and one.

Like ordinary squared coherence, these values do not convey the direction
of causality. Reliable conditional estimates require enough independent
spectral information relative to the number of channels.

## Complex partial coherency, phase, and imaginary coupling

The additional **SciFort extension** `return_coherency=.true.` returns the
complex, phase-preserving partial coherency

\[
\gamma_{ij\mid\mathrm{others}}(f) =
  -\frac{P_{ij}(f)}{\sqrt{P_{ii}(f)P_{jj}(f)}} \quad(i\ne j),
\qquad \gamma_{ii}=1.
\]

Its squared magnitude equals the existing partial coherence (up to numerical
roundoff): `abs(result%coherency)**2 == result%coherence`. The minus sign is
part of the inverse-spectral conditional-correlation convention. The complex
coherency matrix is Hermitian, with unit real diagonal. For example:

```fortran
result = multitaper_partial_coherence(signals,3.5_dp, &
    kmax=6,shrinkage=0.1_dp,return_coherency=.true.)
if (result%valid(2)) then
    ! Imaginary (signed) partial coherency for the two channels:
    print *, aimag(result%coherency(2,1,3))
    ! Phase in radians; meaningful only when coherency is nonzero:
    print *, atan2(aimag(result%coherency(2,1,3)), &
        real(result%coherency(2,1,3),dp))
end if
```

The Fourier cross-spectrum convention `S_ij = conj(FFT_i)*FFT_j`
determines the **sign of the reported phase**, and reversing channel order
reverses it. Phase is meaningless at zero partial coherency and can be
unstable near zero. Neither complex nor squared partial coherence establishes
causality. The optional complex array costs an additional
`O(nfrequency*nchannels**2)` storage; ordinary calls do not allocate it.

## Diagonal-target shrinkage

The new optional `shrinkage=alpha` parameter, with `0 <= alpha <= 1`, forms

\[
S_\alpha(f)=(1-\alpha)S(f)+\alpha\operatorname{diag}(S(f)).
\]

It preserves **individual channel powers** while shrinking every
cross-spectrum by `1-alpha`. Unlike trace-based ridge loading, this rule is
invariant to independent *positive real* rescalings of channels, when no
additional ridge is applied. `alpha=1` makes every nonzero-diagonal channel
spectrally independent in the regularized estimate, giving zero off-diagonal
coherency. By default `shrinkage=0.0_dp`. For a separate, opt-in Gaussian-Wishart
plug-in estimate, use `auto_shrinkage=.true.`; see
[automatic diagonal shrinkage](signal_auto_shrinkage.md). Alternatively,
`cv_shrinkage=.true.` selects leave-one-taper-out predictive shrinkage for
equal, eigen **or adaptive** taper weighting, without assuming fixed weights.
The spectral-matrix call must also have `cv_shrinkage=.true.` to populate
its coefficients before calling `partial_coherence_matrix`. See
[cross-validated shrinkage](signal_cv_shrinkage.md). If `ridge` is also supplied, trace-based
ridge loading is applied **after** shrinkage, so the two regularizers can be
combined. Shrinkage reduces variance/conditioning problems at the cost of
biasing conditional relationships toward zero. It does **not** automatically
establish statistically reliable inference.

Directly measured cross-spectral matrices may have zero or nearly zero
eigenvalues, for example when two channels are identical or when fewer
tapers than channels contribute to the estimate. The optional `ridge`
argument applies *relative diagonal loading* before inversion:

\[
S_{\lambda}(f) = S(f) + \lambda\,
\frac{\mathrm{trace}(S(f))}{n_\mathrm{channels}}I.
\]

`ridge=0.0_dp` (default) performs no regularization; positive `ridge` values
make rank-deficient but nonzero spectral matrices invertible and bias the
result toward smaller conditional associations. This convention is invariant
to multiplying all channels by the same nonzero real factor.
There is **no automatic tuning of the ridge**. In particular, small
regularization does not imply statistically reliable estimates.

For numerical stability, SciFort symmetrizes small roundoff differences in
Hermitian matrices and uses a complex Hermitian Cholesky factorization followed
by triangular solves (no determinant expansions and no normal equations).

## Direct calculation from an estimated spectral matrix

If the spectral matrix has already been computed, reuse it without repeating
the DPSS and FFT calculations:

```fortran
use scifort_signal, only : dp, multitaper_spectral_matrix, &
    multitaper_matrix_result, partial_coherence_matrix, partial_coherence_result

type(multitaper_matrix_result) :: spectral
type(partial_coherence_result) :: result

spectral = multitaper_spectral_matrix(signals,3.5_dp,kmax=6)
result = partial_coherence_matrix(spectral,ridge=0.01_dp)
```

`partial_coherence_matrix` takes a `multitaper_matrix_result`, including
its frequency array and `power(nfrequencies,nchannels,nchannels)`. Any
Hermitian cross-spectral matrix may be supplied in that type; it does not
have to originate from multitaper estimation.

## Status, singular bins, and validation

The result contains

- `frequencies(:)` – frequency values in the input's units;
- `coherence(:,:,:)` – real squared partial coherence;
- `coherency(:,:,:)` – optional complex coherency, empty unless explicitly requested;
- `precision(:,:,:)` – complex inverse of each regularized matrix;
- `valid(:)` – whether each frequency could be inverted;
- `shrinkage_used(:)` – actual per-frequency shrinkage (NaN for zero spectrum);
- `status` – `partial_coherence_success` (0),
  `partial_coherence_invalid_input` (1), or
  `partial_coherence_singular` (2).

An invalid input (wrong shape, nonfinite spectrum, non-Hermitian matrix,
negative or nonfinite ridge, shrinkage outside `[0,1]` or nonfinite,
missing/invalid effective taper count, or combining automatic and manual shrinkage) returns status 1; malformed inputs initially
produce empty result arrays. Structurally valid spectral data with one or more
singular frequency bins returns status 2, preserves valid neighboring bins,
sets the invalid bin's coherence, requested coherency and precision entries to NaN, and identifies
these bins with `valid(f)=.false.`. A zero spectral matrix cannot be rescued
by relative diagonal loading because its trace is zero. The Cholesky test
also rejects matrices that are not positive definite to working precision.

The public routines support real and complex `real64` data and noncontiguous
rank-2 actual arguments. They do not automatically select the sample axis
from the input shape. There is no added Fortran package dependency.

## Validation

Run from the project root:

```bash
python tools/run_signal_partial_coherence.py
```

The standalone runner performs the independent NumPy matrix-inversion and
SciPy 1.17 DPSS/FFT reference comparisons, native singularity/validation
tests, and existing multitaper regressions, with a checked compiler build and
a *targeted* optimized build of the new module. An independent public-module
smoke test is also provided in
`test/test_signal_partial_coherence_public.f90`.

The new partial-coherence formulas and Cholesky implementation are
MIT-licensed SciFort contributions. SciPy DPSS and NumPy linear algebra
are used solely in the reference Python script, not in the Fortran library.

Complex coherency and shrinkage are also tested with independent NumPy
inversions and SciPy DPSS/FFT matrix estimates by running
`python tools/run_signal_partial_coherency.py` (add `--legacy-references`
to repeat both older independent reference suites).
