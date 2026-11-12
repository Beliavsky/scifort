# VAR(p) conditional canonical coherence: fixed DPSS multitaper coverage study

**SciFort exploratory simulation, not a `scipy.signal` function or a
calibrated statistical method.** This study extends the earlier
[one-Hann-taper VAR(p) coverage experiment](signal_varp_conditional_coverage.md)
using the native [exact-band fixed-taper spectral estimator](signal_fixed_taper_spectrogram.md)
and the [fixed-weight VAR(p) population targets](signal_varp_multitaper_targets.md).

## Reproduction

```bash
# Checked and optimized builds, randomized spectral reference and quick study
python tools/run_signal_varp_multitaper_coverage.py

# Reproduce the included 1,620-row exploratory diagnostic CSV:
python tools/run_signal_varp_multitaper_coverage.py --checked-only --study \
    --trials 8 --resamples 59 \
    --settings 2.5:1:equal,2.5:3:eigen,3.5:5:eigen \
    --block-lengths 2,4,8 \
    --csv docs/data/signal_varp_multitaper_coverage_8.csv

# Larger recommended starting point for a scientific coverage investigation:
python tools/run_signal_varp_multitaper_coverage.py --checked-only --study \
    --trials 256 --resamples 999 \
    --csv docs/data/multitaper_coverage_large.csv
```

The longer setting is computationally expensive and **has not been run**.
The 8-dataset experiment above is a reproducibility example, not enough
Monte Carlo replication to rank bootstrap methods. Dependencies for the
validation/study scripts are NumPy, SciPy and `gfortran`; the Fortran
library itself requires neither NumPy/SciPy nor external BLAS/LAPACK.

The study tests 4 stable five-channel Gaussian VAR(2)/VAR(3) models:
separated modes, cross-lag dependence, a weak second canonical mode, and
near-tied canonical modes. Each dataset starts from the **exact stationary
companion covariance**, rather than from a finite burn-in, and has 40
**nonoverlapping 16-sample windows**. Canonical groups are channels `[1,2]`
and `[3,4]`, conditional on channel `[5]`.

For each model, the **same simulated sample paths** are reused for every
DPSS configuration and bootstrap block length. Tapers are SciPy periodic
DPSS (`sym=False`, L2 norm) with fixed either equal or concentration-ratio
weights. They are passed as raw numeric arrays to native Fortran, which
computes the complete complex cross-spectral matrix separately for every
window. No FFT grid or Python-computed spectral matrices enter the
Fortran estimator. The existing Fortran bootstrap resamples contiguous
blocks of whole spectral windows with circular wraparound and computes five
methods: percentile, basic, BC, BCa, studentized.

## Population target and independent checks

The target for a fixed bank is the conditional canonical functional of
**the expected multitaper spectral matrix**:

\[
\theta(f)=C\left(\sum_k \eta_k E[S_{w,k}(f)]/\sum_k\eta_k\right).
\]

`C` denotes the squared conditional canonical roots (sorted) or the
aggregate determinant-based coherence, computed **after** taper averaging
and after exact band integration. It is not `E[C(S)]` or the nonlinear
finite-sample estimator's expectation. The independent Python target uses
SciPy `solve_discrete_lyapunov`, explicit time-pair autocovariance sums,
then Schur complements, Cholesky and singular-value decomposition.
The separate independent sample-path check recomputes every spectral
matrix from Fourier sums and the analytic pairwise band kernel; each band
is also verified against 96-point direct Fourier quadrature. The full
Fortran canonical point estimate agrees with independent Python calculations
within approximately `1e-15` in the smoke suite.

Three spectral statistics are computed for each frequency selection:
leading squared canonical mode, second squared mode, and aggregate group
coherence. Frequency selections are the point `0.1375` cycles/sample and
exact means over `[0,0.125]` and `[0.125,0.5]`. Each scalar/method row in
`docs/data/signal_varp_multitaper_coverage_8.csv` reports population target,
leading/second gap, number of valid intervals, observed hits, coverage among
**valid** intervals, Wilson Monte Carlo bounds, mean interval width and
full-estimator bias/RMSE.

## Illustrative 8-dataset results

For the leading root at `f=0.1375`, pooled equally over the four models:

| DPSS bank | Block | Percentile coverage | Studentized coverage |
|---|---:|---:|---:|
| NW=2.5, K=1, equal | 2 | 62.5% | 100.0% |
| NW=2.5, K=1, equal | 4 | 71.9% | 100.0% |
| NW=2.5, K=1, equal | 8 | 59.4% | 100.0% |
| NW=2.5, K=3, eigen | 2 | 81.2% | 90.6% |
| NW=2.5, K=3, eigen | 4 | 68.8% | 84.4% |
| NW=2.5, K=3, eigen | 8 | 75.0% | 90.6% |
| NW=3.5, K=5, eigen | 2 | 84.4% | 93.8% |
| NW=3.5, K=5, eigen | 4 | 75.0% | 90.6% |
| NW=3.5, K=5, eigen | 8 | 75.0% | 96.9% |

Each row pools at most **32 very small-sample coverage indicators**, with
only 59 bootstrap resamples per underlying dataset. The rows also reuse
sample paths across configurations and are not independent comparisons.
The near-tied model is especially problematic: numeric convergence and
`valid` flags cannot by themselves establish reliable rank-specific mode
identification. These percentages must not be interpreted as a ranking or
evidence of nominal 95% calibration.

## Limitations and next steps

- Higher simulation replication and bootstrap resample counts, different
  independent seeds, longer time-series and window lengths are essential
  before inferring calibration or comparing block-selection methods.
- Taper count, time-bandwidth product and weights affect both the finite
  estimator **and its matched population target**. Different targets should
  not be compared as if they were a common parameter.
- Mode crossings and near-degeneracies require subspace-level confidence
  regions; sorted mode-specific scalar intervals alone are inadequate.
- The study omits data-adaptive weights, smoothing/detrending, overlapping
  windows, ridge/shrinkage tuning and automatic block-length estimation.
- Wilson ranges describe Monte Carlo uncertainty *conditional on interval
  validity*, which must also be assessed separately. A high validity rate
  does not imply coverage.
