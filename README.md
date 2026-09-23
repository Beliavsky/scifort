# SciFort

SciFort is a proposed modern Fortran statistics library with an idiomatic
Fortran API, an FPM build, and a stable C ABI for future Python, R, MATLAB, and
Octave bindings.

This repository is an independent project. It is not affiliated with or
endorsed by SciPy, NumFOCUS, or the Python Software Foundation.

## Current scope

Cloning or downloading the current `main` branch provides all the features
described below. For release history, see [CHANGELOG.md](CHANGELOG.md).

SciFort provides a broad collection of continuous, discrete, and multivariate
probability distributions, covering most distribution families available in
`scipy.stats`. Examples include normal, gamma, beta, binomial, Poisson, and
multivariate normal. Supported operations include density and mass functions,
cumulative and survival probabilities, quantiles, and random sampling. See the
[API reference](docs/api.md) for the complete list, supported operations,
parameterizations, and limitations.

The dependency-free special-function layer provides regularized incomplete
gamma and beta functions (`gammainc`, `gammaincc`, `betainc`, `betaincc`) with
their inverses. `scifort_special` also provides error/inverse-error
functions, `gammaln`, positive-shape `betaln`, real `digamma`/`psi`, normal-CDF
transforms, logistic transforms, entropy helpers, Box-Cox transforms, and
one-dimensional `logsumexp`, `softmax`, and `log_softmax`. The
special-function layer also exposes scaled modified-Bessel `i0e` and `i1e` kernels shared by the Rice and von Mises
families, plus a real Hurwitz-zeta kernel and its exponent derivative used by
the Zipf family, real modified-Bessel-K log/derivative kernels, including a scaled-log form, used by the
generalized hyperbolic, generalized inverse-Gaussian, and normal-inverse-Gaussian families, and
real-order scaled modified-Bessel-I log and adjacent-order ratio kernels used by von Mises-Fisher.

Multivariate normal exposes PDF/log-PDF, CDF/log-CDF with optional finite lower limits,
marginal extraction, entropy, explicit-state RNG, and maximum-likelihood fitting,
including SciPy-compatible positive-semidefinite covariance handling and lower-triangle
semantics. Multivariate Student t exposes PDF/log-PDF, CDF with optional finite lower
limits, marginal extraction, entropy, and explicit-state RNG. Their multidimensional
CDFs share a Genz-style conditional transformation with randomized low-discrepancy
integration. Dirichlet, multinomial, multivariate hypergeometric, and
normal-inverse-gamma provide their SciPy-style density/mass, moment, and random-variate
surfaces where SciPy exposes them. The matrix-valued layer provides matrix normal,
Wishart, inverse-Wishart, and matrix-t densities together with moments/entropy where
applicable and explicit-state RNG. The shared `scifort_linalg` and matrix helper
modules provide small dense Cholesky, triangular-solve, symmetric-eigensystem, PSD,
and SPD-solve kernels without an external BLAS/LAPACK link requirement. The directional/random-matrix
layer adds unit-sphere and von Mises-Fisher sampling/densities, Haar orthogonal/special-orthogonal/unitary
generators, prescribed-spectrum random correlation matrices, and fixed-margin random contingency tables.


The statistics layer also includes `gaussian_kde`, with weighted or
unweighted data, Scott/Silverman/constant bandwidth factors, PDF/log-PDF evaluation,
Gaussian and box integration, KDE-product integration, marginals, inverse covariance,
and explicit-state resampling. The `scipy.stats.qmc` layer provides sample scaling
and integer mapping, centered/wrap-around/mixture/L2-star discrepancies, geometric
minimum-distance/MST discrepancy, iterative discrepancy updates, Van der Corput,
Halton, and Joe-Kuo Sobol sequences, ordinary or strength-two orthogonal-array Latin
hypercubes, Poisson-disk sampling, and QMC multinomial and multivariate-normal
samplers. Sobol supports up to 21,201 dimensions, 64 direction-number bits, reset,
fast-forward, base-two blocks, and SciFort-state LMS+digital-shift scrambling. QMC
samples follow SciPy's `(n,d)` layout, unlike the column-oriented random-variate arrays
used by multivariate distributions.

Each scalar univariate distribution currently provides:

- `pdf`
- `logpdf`
- `cdf`
- `sf`
- `logcdf`
- `logsf`
- `ppf`
- `isf`

Most scalar distribution functions are `pure elemental`, so the same routines operate on
scalars and conformable arrays. The Poisson-binomial family takes a rank-one
probability vector `p(:)` and is therefore pure scalar rather than elemental.
Invalid location or scale parameters return a quiet NaN. Probability arguments
outside `[0, 1]` also return a quiet NaN.

Continuous distributions expose PDF, log-PDF, CDF, survival function,
log-CDF, log-survival, PPF, and ISF functions. Discrete distributions expose
PMF, log-PMF, CDF, survival function, log-CDF, log-survival, PPF, and ISF
functions. See `docs/api.md` for parameterizations, support, endpoint
behavior, and the complete public API.

SciFort also has an explicit-state random-number interface in
`scifort_random` and scalar or array `*_rvs` sampling routines for every
current distribution. Sampling is deterministic for a given `rng_state` and
never relies on Fortran's process-global `random_number` state. Most family samplers use inverse transforms through the existing PPF kernels,
prioritizing distributional consistency and reproducibility over sampling
speed. A few finite/discrete families use their defining exact samplers instead:
Wallenius draws sequentially from the biased urn, and Poisson-binomial sums
independent Bernoulli draws.

Array log-likelihood and negative-log-likelihood (`*_nnlf`) evaluators are
available for all current families through `scifort_likelihood` and are also
re-exported by `scifort_stats`. Discrete likelihoods accept either integer or
`real(dp)` observation arrays. `scifort_score` adds likelihood scores for all 130 families, analytic where a
closed-form derivative is practical. The finite-sample two-sided KS density and
its loc/scale score use local high-order differentiation, matching the numerical
nature of the finite-sample density itself. `scifort_fit` adds deterministic bounded maximum-
likelihood fits with optional starting guesses, fixed parameters represented by
equal lower and upper bounds, integer constraints for discrete fit parameters,
and explicit convergence diagnostics. The first fitter is a local projected
coordinate/pattern search and requires finite bounds for every parameter.

The normal CDF uses the standard `erfc` identity. The initial normal quantile
implementation uses safeguarded bisection against the CDF or log-CDF. It is a
correctness-first implementation intended to be replaced or supplemented by
a faster, carefully licensed approximation after accuracy and provenance
reviews.

SciFort includes native descriptive statistics and inference. `scifort_descriptive`
provides mean, sample/population variance through configurable `ddof`, standard
deviation, central moments, linear quantiles, median, covariance, Pearson
correlation, and SciPy-style `rankdata` tie methods. `scifort_hypothesis`
provides t-tests, Pearson/Spearman correlation tests, and asymptotic
Mann-Whitney U. `scifort_hypothesis_extended` adds Wilcoxon rank-sum and
signed-rank tests, Kruskal-Wallis and Friedman tests, Cressie-Read power
divergence and chi-square tests, 2-by-2 Fisher exact tests, contingency-table
chi-square, classical and Welch one-way ANOVA, and Bartlett, Levene, and
Fligner-Killeen variance tests. `scifort_nonparametric_extended` adds the
Ansari-Bradley and Mood scale tests, Epps-Singleton two-sample test,
k-sample Anderson-Darling test, and Mood median test. `scifort_association_extended`
adds Kendall tau-b/tau-c, point-biserial correlation, ordinary linear regression,
Theil-Sen and Siegel robust slopes, the Brunner-Munzel test, and Page's ordered-
alternative trend test. `scifort_contingency_meta` adds exact binomial tests,
Barnard and Boschloo unconditional 2-by-2 tests, conditional/sample odds ratios,
relative risk, contingency association/margins, p-value combination, and
Benjamini-Hochberg/Benjamini-Yekutieli false-discovery control.
`scifort_multiple_comparisons` adds Alexander-Govern heterogeneous-variance
mean testing, Tukey HSD/Tukey-Kramer and Games-Howell all-pairs comparisons
with simultaneous confidence intervals, Dunnett treatment-versus-control
comparisons with simultaneous intervals from the multivariate Student-t law,
and the Krishnamoorthy-Thomson Poisson means E-test.
`scifort_goodness_of_fit` adds one- and two-sample
Kolmogorov-Smirnov and Cramer-von Mises tests,
Shapiro-Wilk, D'Agostino skew/kurtosis/omnibus normality tests, Jarque-Bera,
and Anderson-Darling tests for normal, exponential, logistic, and Gumbel
families. `scifort_resampling` adds explicit-state percentile/basic/BCa
bootstrap intervals, exact-or-randomized permutation tests (`independent`,
`samples`, and `pairings`), Monte Carlo hypothesis tests, and simulated power
estimation through native Fortran callbacks. These interfaces are re-exported
by `scifort_stats` and remain native-Fortran-only.

The experimental C ABI now exposes scalar density or mass, CDF, and PPF
routines for every current distribution, including the Milestone 7
continuous-family additions. Normal PDF, CDF, and PPF also have
bulk array entry points. The Fortran API remains the primary API; survival,
logarithmic, inverse-survival, and special-function entry points are not yet in
the C ABI.

## Requirements

- A Fortran 2018 compiler
- FPM 0.12.0 or newer

The default FPM build produces a static library. FPM 0.13.0 or newer is needed
if the manifest is changed to build static and shared library targets together.

## Build and test

```text
fpm build
fpm test
fpm run
fpm run --example distributions_demo
```

For a debug build with current FPM:

```text
fpm test --profile debug
```

## Fortran use

```fortran
program example
    use scifort_stats, only : dp, normal_cdf, normal_ppf
    implicit none

    real(dp) :: p
    real(dp) :: x(3)

    x = [-1.0_dp, 0.0_dp, 1.0_dp]
    p = normal_cdf(1.96_dp)

    print *, p
    print *, normal_cdf(x)
    print *, normal_ppf(0.975_dp)
end program example
```

To use SciFort from another FPM project:

```toml
[dependencies]
scifort = { git = "https://github.com/Beliavsky/scifort.git" }
```

This tracks the default branch. After a release is tagged, add
`tag = "vX.Y.Z"` to pin that version.

## Repository rules

Coding agents and contributors must read `AGENTS.md`. Claude Code must also
read `CLAUDE.md`. The most important rules are:

- preserve exact code provenance and original license notices
- do not import GPL or LGPL implementation code into the permissive core
- never claim tests were run when they were not
- prioritize numerical correctness, especially in tails and at boundaries
- keep the native Fortran API separate from the stable C ABI
- do not break public APIs or C symbols without an explicit versioning decision

See `CONTRIBUTING.md`, `CODE_PROVENANCE.md`, and the files under `docs/`.

## Validation

The repository is tested with GNU Fortran and FPM on Linux, macOS, and
Windows in GitHub Actions. The current validation record also includes local
debug, release, runtime-checking, and numerical reference tests for the
special functions and distributions. See `docs/validation.md` for exact
commands, reference comparisons, and known limitations.

## C ABI

The experimental header is `include/scifort.h`. The initial ABI uses only
standard C scalar types, explicit array lengths, caller-owned output storage,
and integer status codes. It does not expose compiler-specific Fortran module
symbols or descriptors.

The scalar ABI covers PDF or PMF, CDF, and PPF for all current distributions.
Discrete counts, trial counts, and support shifts are passed as `double`; the
native validation rules still require integer-valued arguments where the
Fortran API does. Only the normal distribution currently has bulk C entry
points. Random-number state, random-variate, likelihood, score, and fitting
routines are currently native-Fortran-only. Python, R, MATLAB, and Octave wrappers should call this
common ABI rather than independently binding to compiler-specific Fortran
interfaces.

## License

Original SciFort code is licensed under the MIT License. Third-party code may
be included only when its license is compatible with distribution in this
project, and every imported or translated component must retain its applicable
notices. See `THIRD_PARTY_LICENSES.md` and `CODE_PROVENANCE.md`.
