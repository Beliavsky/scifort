# Conditional two-mode VAR(p) canonical coherence: coverage experiments

**SciFort research extension, not a `scipy.signal` API.** This is a
reproducible *simulation and reference-test framework*, not a change to the
published Fortran library interfaces. It evaluates the existing five
moving-block bootstrap confidence-interval methods against appropriately
matched finite-window population spectral functionals for **VAR(2)** and
**VAR(3)**. It studies two canonical modes between channels `[1,2]` and
`[3,4]`, **conditioning on channel 5**.

## Reproduce

```bash
# Small checked/optimized smoke tests; NumPy/SciPy and gfortran required.
python tools/run_signal_varp_conditional_coverage.py

# Longer, deterministic study (the included CSV was generated this way).
python tools/run_signal_varp_conditional_coverage.py --checked-only \
    --study --trials 64 --resamples 99 \
    --csv docs/data/signal_varp_conditional_coverage_64.csv

# Alternative dependence lengths / taper lengths can be studied directly.
python tools/run_signal_varp_conditional_coverage.py --checked-only \
    --study --trials 64 --resamples 99 --block-length 8 \
    --samples-per-window 24 --seed 20261031 --csv /tmp/varp_alternative.csv
```

The default test is deliberately short (two independent datasets per model,
39 draws). The full script prints per-scalar population targets, valid-case
coverage, Monte Carlo 95% Wilson score intervals, and mean bootstrap interval
widths. The CSV also records the mean **point-estimator bias** and RMSE
relative to the matched target, the population gap between the two canonical
roots, and numbers of valid intervals, for all 180 combinations of four
models, three frequency selections, three scalar statistics, and five interval
methods.

Scripts and tests:

- `tools/signal_reference/study_canonical_varp_conditional.py` simulates the
  original time series and constructs finite-window spectral matrices.
- `tools/signal_reference/canonical_varp_coverage_driver.f90` calls the
  existing native Fortran canonical moving-block-bootstrap implementation.
- `tools/signal_reference/test_study_canonical_varp_conditional.py` checks the
  stationary model covariance, Fourier-band kernels, input/output protocol,
  Hermitian positive semidefiniteness, and canonical-mode algebra.
- `tools/run_signal_varp_conditional_coverage.py` builds checked and optimized
  GNU Fortran executables and runs all of the above tests.

There is no new Fortran runtime dependency or change in FPM/C-ABI exports.

## Matching the population estimand

For the stationary real Gaussian VAR(p),

\[
X_t=\sum_{j=1}^p A_j X_{t-j}+\varepsilon_t,\qquad
E[\varepsilon_t\varepsilon_t^T]=\Sigma,
\]

initialize its companion state with the **exact stationary Gaussian state
covariance**, computed by SciPy's discrete Lyapunov solve. This avoids the
finite-burn-in approximation that could contaminate a coverage experiment.
Each dataset comprises 40 **nonoverlapping 16-sample Hann-tapered windows**.
They remain temporally dependent because the underlying VAR process is
stationary and correlated across boundaries.

For the tapered Fourier vector
\(Z(f)=\sum_{t=0}^{L-1} w_t X_t e^{-2\pi i f t}\), a single-frequency
window matrix is \(Z(f)Z(f)^*/\sum w_t^2\). The two frequency-band window
matrices are the exact integrals of this outer product, **divided by band
width**, using the pairwise kernel

\[
K_{ts}(a,b)=\operatorname{sinc}[(b-a)(t-s)]
 e^{-\pi i (a+b)(t-s)}.
\]

No interpolated FFT grid enters the calculation. Independently, the
population-target calculation uses SciPy `solve_discrete_lyapunov` and the
direct time-sample-pair autocovariance sum from
`verify_varp_taper_spectrum.reference`, followed by an independent SciPy
Cholesky/SVD implementation of conditional canonical coherence. It computes
\(C(E[ZZ^*])\), or \(C(E[\overline S_{[a,b]}])\) for the bands, **not**
\(E[C(ZZ^*)]\), not coherence at frequency zero, and not the expectation of
the finite-sample nonlinear canonical statistic.

The three frequency selections are point frequency `0.1375` cycles/sample,
mean spectrum over band `[0,0.125]`, and mean spectrum over `[0.125,0.5]`.
The two squared canonical roots are sorted descending; aggregate group
coherence is \(1-(1-r_1)(1-r_2)\). The bootstrap resamples contiguous blocks
of **four complete spectral windows** with circular wrapping and 99 draws
per independent dataset. There is no ridge or shrinkage for this study,
so its population and sample spectra use identical regularization.
The underlying study uses one fixed Hann taper rather than a DPSS multitaper
or adaptive spectral weighting; results must not be generalized to those
estimators without a separate simulation.

## Illustrative results (64 datasets per model)

The following table gives empirical 95%-interval coverage of the **leading
squared canonical mode at f=0.1375**. All listed intervals were numerically
valid, but that does not establish their inferential reliability.

| Model | Population mode 1 | Mode 1-2 gap | Percentile | Basic | BCa | Studentized |
|---|---:|---:|---:|---:|---:|---:|
| VAR(2), separated modes | 0.4555 | 0.2478 | 65.6% | 82.8% | 60.9% | 90.6% |
| VAR(3), cross-lag coupling | 0.4980 | 0.2342 | 62.5% | 73.4% | 70.3% | 89.1% |
| VAR(2), weak second mode | 0.4809 | 0.4501 | 81.2% | 85.9% | 81.2% | 95.3% |
| VAR(3), nearly tied modes | 0.4606 | 0.0005 | 28.1% | 73.4% | 40.6% | 87.5% |

Mean leading-mode point-estimator biases at that frequency were approximately
`+0.0501`, `+0.0528`, `+0.0249`, and `+0.1151`, respectively. The final
model's second sorted root had mean bias around `-0.0837`; this split is a
consequence of **ordering nearly tied noisy modes**. Mode 1 and mode 2 are
not reliably separately identified at the near tie. Its numeric `valid`
flags should *not* be interpreted as proving population-mode separation.
Aggregate group coherence is invariant to exchanging the two roots and is
often a more interpretable target in this situation.

Across all 36 model/frequency/statistic combinations, studentized intervals
were valid in at least 63 of 64 datasets per combination; other methods were
valid in all 64 for this particular experiment. Empirical average coverage
across those combinations was around 90.5% for studentized and 77.2% for
percentile. These averages weight heterogeneous statistics equally; they are
not formal performance rankings and no method demonstrated guaranteed 95%
coverage. Finite-sample bias, near-tied roots, block selection and mode
ranking all affect results.

## Limits and next experiments

- **64 datasets give substantial Monte Carlo uncertainty.** Use Wilson score
  ranges and multiple independent seeds before inferring differences between
  methods. Only 99 bootstrap replicates were drawn per dataset.
- Pointwise intervals are not simultaneous frequency/mode confidence regions.
  Repeated or nearly tied canonical roots require inference on their joint
  *subspace* rather than rank-labeled eigenvectors/roots.
- The population functional `C(E[S_window])` is the population plug-in
  target of the mean-window estimator, but the finite-sample estimator's
  nonlinear bias need not vanish at the sample sizes used here.
- Although windows do not overlap, the underlying process couples them.
  Circular bootstrap blocks do not automatically recover nominal coverage.
  The selected block length is fixed, not optimized for each model.
- This study deliberately holds taper, regularization, bandwidth and sample
  sizes constant. Vary them in follow-up experiments, along with VAR
  coefficients, block length, seed, and independent model replications.

No calibration claim follows from passing the exact numerical references.
