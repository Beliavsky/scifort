# Fixed-weight multitaper VAR(p) population spectral targets

**SciFort extension, not a `scipy.signal` API.** These routines extend the
[exact finite-taper VAR(p) population targets](signal_varp_population_targets.md)
to a bank of deterministic tapers. They calculate the expected *complex
spectral matrix*, not the expected nonlinear coherence statistic. Band
integration is analytic and is performed **before** conditional canonical
coherence estimation.

## Public interface and example

```fortran
use scifort_signal, only : dp, dpss, dpss_result, &
    varp_taper_spectrum_result, canonical_coherence_result, &
    varp_expected_multitaper_spectrum, &
    varp_expected_band_multitaper_spectrum, &
    varp_expected_multitaper_canonical_coherence, &
    varp_expected_band_multitaper_canonical_coherence, &
    varp_spectrum_success

type(dpss_result) :: tapers
type(varp_taper_spectrum_result) :: point, integrated
type(canonical_coherence_result) :: modes, band_modes
real(dp) :: a(5,5,2), sigma(5,5), bands(2,2)

! Initialize a stable real VAR(2) model in a and symmetric positive-
! definite innovation covariance in sigma before calling these routines.
tapers = dpss(24, 3.5_dp, kmax=5)
if (tapers%status /= 0) error stop 'DPSS tapers failed'
bands(:,1) = [-0.20_dp,0.0_dp]
bands(:,2) = [0.0_dp,0.25_dp]

point = varp_expected_multitaper_spectrum(a,sigma,tapers%tapers, &
    [0.0_dp,0.125_dp,0.25_dp],taper_weights=tapers%ratios, &
    normalized=.true.)
integrated = varp_expected_band_multitaper_spectrum(a,sigma, &
    tapers%tapers,bands,taper_weights=tapers%ratios, &
    average=.true.,normalized=.true.)

modes = varp_expected_multitaper_canonical_coherence(a,sigma, &
    tapers%tapers,[0.0_dp,0.125_dp,0.25_dp],[1,2],[3,4], &
    conditioned_on=[5],taper_weights=tapers%ratios,normalized=.true.)
band_modes = varp_expected_band_multitaper_canonical_coherence(a,sigma, &
    tapers%tapers,bands,[1,2],[3,4],conditioned_on=[5], &
    taper_weights=tapers%ratios,normalized=.true.,average=.true.)

if (point%status == varp_spectrum_success) then
    print *, real(point%power(1,1,1),dp)
end if
```

The four routines are exported from both
`scifort_signal_varp_multitaper_target` and the public umbrella module
`scifort_signal`. They return the existing `varp_taper_spectrum_result`
or `canonical_coherence_result`, without introducing a separate result type,
external dependency, or C ABI change.

## Input and normalization contract

- `a(n,n,p)` is a real, stable VAR(p) with `p>=1`; the finite, symmetric
  `innovation_covariance(n,n)` must be strictly positive definite. The
  model has a zero mean and is initialized in stationarity.
- `tapers(length,K)` has one taper per **column**. For SciFort's `dpss`
  routine this is exactly the shape of `dpss_result%tapers`. Each selected
  taper must be finite and have positive, finite squared energy.
- Optional `taper_weights(K)` supplies **fixed, nonnegative** weights.
  The default is equal weighting. The relative weights are normalized to
  sum to one; an all-zero vector, negative/nonfinite weights, or a length
  mismatch is invalid. A zero-weight taper is ignored, including an all-zero
  unused column. Taper values themselves must be finite, even if unused.
- With `normalized=.false.` (the default), each window contributes its raw
  expected tapered Fourier outer product. With `normalized=.true.`, **each
  taper is individually divided by its own squared energy** before
  averaging. This is important if the tapers have different norms.
- Frequencies and band endpoints are in **cycles/sample** and have the same
  rules as the single-taper VAR(p) API: finite, possibly negative/unsorted
  point frequencies, and finite positive-width, possibly overlapping bands
  in shape `(2,nband)`. `average=.true.` divides the exact band integral by
  band width. Neither option computes a one-sided PSD nor applies an `fs`
  unit conversion.
- Results return `frequencies`, optional `bands`, `lag0_covariance`,
  complex `power(freq,channel,channel)` and the existing status values.
  The population covariance is independent of the taper bank. Invalid
  calls do not return allocated power arrays.
- Canonical wrappers retain the single-taper group's `group_a`, `group_b`,
  optional `conditioned_on`, `ridge`, `shrinkage`, and `return_vectors`
  semantics. Group IDs are one-based and mutually disjoint. Reported
  validity and mode-stability flags are those of the existing canonical
  solver.

## Mathematical definition

For the stationary real VAR(p) process described in the linked document,
let `w(t,k)` be taper `k` with squared energy `e_k>0`, and let
`S_k(f) = E[Z_k(f) Z_k(f)^H]` be its exact finite-taper population
cross-spectral matrix. For fixed weights `eta_k>=0` with positive sum,

\[
 S_{\mathrm{multi}}(f)
 =\frac{\sum_{k=1}^K \eta_k S_k(f)/d_k}{\sum_{k=1}^K \eta_k},
 \qquad
 d_k=\begin{cases}e_k & \text{if normalized},\\1&\text{otherwise.}\end{cases}
\]

Every `S_k(f)` uses the exact stationary lag covariances of the VAR(p)
companion process and the taper autocorrelations. This formula can also be
interpreted as a weighted sum of taper autocorrelations; the current
implementation reuses the already validated single-taper solver for each
selected column. The computational cost consequently scales with `K` as
well as model dimension and window length.

A frequency-band result replaces `S_k(f)` with its **analytically integrated**
frequency-band matrix, and optionally divides by band width. For fixed
weights, integration and taper averaging commute. The order of **nonlinear**
canonical coherence computation does not commute with either operation.
The conditional canonical coherence is computed from the **combined**
spectral matrix, not by averaging separate taper coherences.

The ratio weights from `dpss_result%ratios` give one useful deterministic
choice. Data-adaptive Thomson/iterative weights, including frequencies where
weights depend on the observed periodograms, are **not represented**: an
expectation cannot generally pass through a data-dependent reweighting rule.
Similarly, these are population plug-in targets, not the finite-sample
expectation of estimated canonical coherence, a confidence interval, or an
assertion of calibrated bootstrap coverage.

## Validation

```bash
python tools/run_signal_varp_multitaper.py
python tools/check_signal_varp_population_public.py
```

The first command builds the minimal source closure with checked and
optimized GNU Fortran, executes native regression and invalid-input tests,
and compares independent SciPy/NumPy finite-window targets for 40 randomized
VAR(1)/VAR(2)/VAR(3) models. The reference uses SciPy's
`solve_discrete_lyapunov`, SciPy DPSS windows, an independent
`O(K*length**2)` pairwise Fourier sum and exact band kernels, then complex
Cholesky/SVD canonical modes **after** summing all taper matrices.

The second command checks the full `scifort_signal` dependency closure,
including the public multitaper entry points and earlier canonical/coherence
tests. NumPy and SciPy are verification-only dependencies, not Fortran
runtime dependencies. The Fortran Package Manager's full project suite
remains the recommended test where FPM and its project dependencies are
available.
