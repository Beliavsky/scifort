# Changelog

All notable changes will be documented here.

## Unreleased

### Fixed

- Added merge-regression tests before fixing Zipf NaN/infinity handling,
  default-integer overflow in hypergeometric/discrete-uniform observations and
  discrete-uniform quantiles, and large-offset log-softmax normalization.
- Preserve scalar and diagonal-covariance multivariate-normal log probabilities
  with exact log-domain interval products. Correlated integration retains its
  documented probability-domain limits.
- Reconciled SciPy adaptation notices and provenance, retained the full CC BY
  4.0 research license, and excluded verbatim third-party license texts from
  project ASCII checks.

- Made the Fisher lower-tail regression tolerance portable across log-gamma
  implementations, with an independently verified exact reference of 427/429.

### Added

- Expanded `scifort_resampling` to SciPy-style percentile/basic/BCa bootstrap
  intervals (including paired two-sample resampling), exact-or-randomized
  permutation tests for `independent`, `samples`, and `pairings`, one-/two-sample
  Monte Carlo hypothesis tests, and one-/two-sample simulated power estimation.
  Added SciPy 1.17.0 deterministic references and
  `test/test_resampling_inference.f90`; the legacy permutation regression now
  follows SciPy's exact-switching rule.
- Added `scifort_multiple_comparisons` with Alexander-Govern, Tukey HSD/
  Tukey-Kramer, Games-Howell, Dunnett, and the Poisson means E-test. Tukey/
  Games-Howell results provide simultaneous studentized-range confidence
  intervals; Dunnett p-values and confidence intervals reuse SciFort's
  multivariate Student-t box integrator. Added deterministic SciPy 1.17.0
  reference generation and `test/test_multiple_comparisons.f90`.
- Added `scifort_contingency_meta` with exact binomial tests and confidence
  intervals, Barnard and Boschloo unconditional 2-by-2 exact tests, sample and
  conditional odds ratios with confidence intervals, relative risk, Cramer/
  Tschuprow/Pearson contingency association, 2D margins, Fisher/Pearson/
  Mudholkar-George/Tippett/Stouffer p-value combination, and Benjamini-
  Hochberg/Benjamini-Yekutieli false-discovery control. Added generated
  SciPy 1.17.0 references and `test/test_contingency_meta.f90`.
- Added `scifort_association_extended` with Kendall tau-b/tau-c, point-biserial
  correlation, SciPy-style least-squares `linregress`, Theil-Sen and Siegel
  robust slope estimators, Brunner-Munzel, and Page's ordered-alternative trend
  test. Kendall includes exact untied inversion-count probabilities; Page exact
  p-values use convolution of the one-row permutation distribution. Added
  generated SciPy 1.17.0 references and `test/test_association_extended.f90`.
- Added `scifort_nonparametric_extended` with Ansari-Bradley and Mood scale
tests, the Epps-Singleton two-sample test, k-sample Anderson-Darling with
midrank/right/continuous variants, and the Mood median test with below/above/
ignore tie conventions. Small untied Ansari-Bradley problems use an exact
subset dynamic program. Added generated SciPy 1.17.0 references and
`test/test_nonparametric_extended.f90` covering exact, tied, one-sided,
custom-grid, capped/interpolated, and contingency-table cases.
- Added `scifort_goodness_of_fit`: one- and two-sample Kolmogorov-Smirnov,
  one- and two-sample Cramer-von Mises, Shapiro-Wilk, D'Agostino skew/kurtosis
  and omnibus normality tests, Jarque-Bera, and Anderson-Darling for normal,
  exponential, logistic, and left/right Gumbel families. Added SciPy 1.17.0
  generated references and `test/test_goodness_of_fit.f90`.
- Expanded SciPy-style hypothesis/statistical inference with Wilcoxon rank-sum
  and signed-rank tests, Kruskal-Wallis and Friedman tests, Cressie-Read power
  divergence/chi-square, 2-by-2 Fisher exact and contingency chi-square tests,
  classical and Welch one-way ANOVA, and Bartlett, Levene/Brown-Forsythe, and
  Fligner-Killeen scale tests. Added unequal-size `sample_group` support, exact
  signed-rank subset-sum probabilities, small tied/zero sign enumeration, and
  `test/test_hypothesis_extended.f90` with SciPy 1.17.0 reference values.
- Added `gaussian_kde` with weighted/unweighted initialization, Scott/Silverman or
  constant bandwidths, PDF/log-PDF evaluation, Gaussian/box/KDE-product integration,
  inverse covariance, marginal extraction, and deterministic explicit-state resampling.
  Added SciPy 1.17.0 generated reference data and `test/test_gaussian_kde.f90`.
- Added the first `scipy.stats.qmc` slice: unit-cube scaling and integer mapping,
  CD/WD/MD/L2-star discrepancy, centered-discrepancy updating, geometric minimum-
  distance and MST criteria, Van der Corput, stateful Halton with optional digit
  scrambling, and ordinary/strength-two OA Latin-hypercube sampling. Unscrambled
  sequence and discrepancy values are validated directly against SciPy 1.17.0;
  randomized paths use SciFort's explicit RNG and are tested by replay/invariants.

- Added the advanced QMC slice: Joe-Kuo Sobol sequences through 21,201 dimensions
  and 64 bits with base-two blocks, reset/fast-forward, and LMS+digital-shift
  scrambling; `MultinomialQMC`; inverse-transform or Box-Muller
  `MultivariateNormalQMC`; and bounded Poisson-disk sampling with volume/surface
  proposals and fill-space exhaustion. Unscrambled Sobol, multinomial-QMC, and
  multivariate-normal-QMC outputs are validated directly against SciPy 1.17.0,
  including deep skipped sequences, singular covariance, and bit exhaustion.

- Added the directional/random-matrix SciPy block: `uniform_direction`,
  `vonmises_fisher`, `ortho_group`, `special_ortho_group`, `unitary_group`,
  `random_correlation`, and `random_table`. The block includes explicit-state RNG,
  von Mises-Fisher PDF/log-PDF, entropy and fitting, Haar orthogonal/unitary draws,
  Davies-Higham prescribed-spectrum correlation matrices, and fixed-margin table
  PMF/log-PMF, mean, and exact conditional-hypergeometric sampling.
- Added real-order scaled modified-Bessel-I log and adjacent-order ratio kernels for
  von Mises-Fisher normalization/entropy, with positive-series and large-argument
  asymptotic regimes. Added generated SciPy 1.17.0 references and
  `test/test_directional_random_matrix.f90` covering densities, fitting, all sampler
  branches, structural matrix invariants, edge cases, and deterministic RNG replay.
- Extended the SciPy-style multivariate layer with `multivariate_hypergeom` and
  `normal_inverse_gamma`, including PMF/PDF, first and second moments, SciPy-compatible
  support/degenerate-case handling, and explicit-state random variates.
- Added the matrix-valued distribution core: `matrix_normal`, `wishart`, `invwishart`,
  and SciPy 1.17's `matrix_t`, with PDF/log-PDF, moments and entropy where applicable,
  lower-triangle SPD semantics, and deterministic explicit-state RNG. Shared matrix
  helpers provide SPD factorization/solve, multivariate-gamma, and multivariate-digamma
  building blocks without an external BLAS/LAPACK dependency.
- `test/test_matrix_distributions.f90`, generated SciPy 1.17.0 references, and
  `tools/generate_matrix_distributions_reference.py`, covering densities/masses,
  moments, entropy, lower-triangle behavior, degenerate hypergeometric populations,
  matrix-t documentation values, support/error cases, and RNG replay.

- Expanded the SciPy-style multivariate distribution layer with multivariate-normal
  CDF/log-CDF, finite lower integration limits, marginal extraction, and a new
  `multivariate_t` family providing PDF/log-PDF, CDF, finite lower limits, entropy,
  marginal extraction, and explicit-state RNG. The normal and Student-t CDFs share
  a Genz-style conditional transform with randomized tent-transformed Halton
  integration, including singular positive-semidefinite matrices and deterministic
  replay when no RNG state is supplied.
- Initial SciPy-style multivariate distribution layer: `multivariate_normal`,
  `dirichlet`, `multinomial`, and `dirichlet_multinomial`, re-exported through
  `scifort_stats`. Multivariate normal provides PDF/log-PDF, entropy, singular-PSD
  support, explicit-state RNG, and MLE fitting; Dirichlet and multinomial provide
  densities/masses, moments, entropy, and RNG; Dirichlet-multinomial provides
  PMF/log-PMF and moments.
- `scifort_linalg`, a compact dependency-free dense linear-algebra foundation with
  lower Cholesky factorization, lower-triangular solve, cyclic-Jacobi symmetric
  eigendecomposition, and SciPy-compatible covariance PSD/rank decomposition.
- `test/test_multivariate.f90`, SciPy 1.17.0 numerical references, and
  `tools/generate_multivariate_reference.py`, covering covariance lower-triangle
  semantics, singular support, fitting, moments/entropy, deterministic RNG,
  two/three/four-dimensional normal and Student-t CDFs, finite boxes, marginals,
  and higher-dimensional pivot paths.

- Levy-stable (`levy_stable`) and studentized-range (`studentized_range`)
  continuous distributions with complete eight-function APIs, likelihood/NNLF,
  finite-difference likelihood scores, bounded/fixed fitting, deterministic
  inverse-transform RNG, `scifort_stats`, and scalar C-ABI integration.
- `levy_stable` matches SciPy's default S1 parameterization and adapts the
  SciPy 1.17.0 Nolan/Zolotarev piecewise S1-to-S0 formulas and difficult-input
  handling, while using SciFort-owned Gauss-Legendre/adaptive quadrature.
  `studentized_range` independently evaluates the normal-range/chi-square
  scale-mixture integrals and switches to the infinite-df range formula at
  `df >= 100000`.
- `test/test_levy_stable_studentized_range.f90`, generated SciPy 1.17.0
  references, and `tools/generate_levy_stable_studentized_range_reference.py`,
  covering special-case reductions, loc/scale behavior, direct-tail inverses,
  numerical scores, deterministic RNG, fixed fitting, and all six C ABI calls.

- Finite-sample two-sided Kolmogorov-Smirnov (`kstwo`) distribution with
  complete eight-function API, likelihood/NNLF, loc/scale score, bounded/fixed
  fitting with integral sample-size constraints, deterministic inverse-transform
  RNG, `scifort_stats`, and scalar C-ABI integration.
- The `kstwo` CDF/SF adapts SciPy 1.17.0's BSD-licensed regime selection,
  Durbin/Marsaglia-Tsang-Wang matrix algorithm, and Pelz-Good expansion; exact
  endpoint formulas and SciFort's one-sided KS implementation handle the outer
  regimes. The PDF and loc/scale score use high-order local differentiation.
- `test/test_kstwo.f90`, generated SciPy 1.17.0 references, and
  `tools/generate_kstwo_reference.py`, including exact `n=1`, one-sided-tail,
  deep direct-tail inversion, score, RNG, fixed-fit, and C-ABI checks.

- Wallenius noncentral-hypergeometric (`nchypergeom_wallenius`) and
  Poisson-binomial (`poisson_binom`) discrete families with complete
  eight-function PMF APIs, likelihood/NNLF, analytic continuous-parameter
  scores, bounded/fixed fitting, exact explicit-state RNG, `scifort_stats`, and
  C-ABI integration.
- Wallenius probabilities and the odds score are propagated by an exact
  finite-state dynamic program for the defining sequential biased-urn process;
  no BiasedUrn implementation source is used. Poisson-binomial probabilities use
  Bernoulli convolution with direct log tails, while its vector score uses
  leave-one-out convolutions.
- `test/test_wallenius_poisson_binom.f90`, generated SciPy 1.17.0 references,
  and `tools/generate_wallenius_poisson_binom_reference.py`, including
  Wallenius/hypergeometric and Poisson-binomial/binomial reduction identities,
  score finite differences, exact RNG replay, fixed-fit, shifted-lattice, and
  array-aware C-ABI coverage.

- Landau (`landau`) distribution with complete eight-function API,
  likelihood/NNLF, analytic loc/scale score, bounded/fixed fitting, deterministic
  inverse-transform RNG, `scifort_stats`, and scalar C-ABI integration.
- Hybrid Landau numerics using a Koelbig-Schorr/Yoshimura left-tail asymptotic
  expansion, exponentially damped characteristic-function quadrature centrally,
  and a scaled nonoscillatory defining integral for the heavy right tail.
- `test/test_landau.f90`, generated SciPy 1.17.0 references, and
  `tools/generate_landau_reference.py`, including deep direct-tail round trips.

- Noncentral Student t (`nct`) and Gauss hypergeometric (`gausshyper`) families
  with complete eight-function APIs, likelihood/NNLF, analytic scores,
  bounded/fixed fitting, deterministic inverse-transform RNG, `scifort_stats`,
  and scalar C-ABI integration.
- Adaptive double-exponential gamma-probability quadrature for noncentral t,
  with finer grids for degrees of freedom below one, and beta-quantile-transform
  quadrature for Gauss hypergeometric normalization and direct lower/upper tails.
- `test/test_nct_gausshyper.f90`, generated SciPy 1.17.0 references, and
  `tools/generate_nct_gausshyper_reference.py`, including independent
  high-precision stress constants for low-df noncentral t and concentrated/deep
  Gauss-hypergeometric tails where SciPy 1.17.0 loses accuracy.

- Generalized hyperbolic (`genhyperbolic`) and Fisher noncentral-hypergeometric
  (`nchypergeom_fisher`) families with complete eight-function APIs,
  likelihood/NNLF, analytic scores for continuous parameters, bounded/fixed
  fitting, deterministic inverse-transform RNG, `scifort_stats`, and scalar
  C-ABI integration.
- A scaled-log modified-Bessel-K evaluator in `scifort_special`, used to keep
  generalized-hyperbolic normalization and far-tail density calculations stable
  without forming `exp(x)*K_nu(x)` explicitly.
- `test/test_genhyperbolic_fisher.f90`, generated SciPy 1.17.0 references, and
  `tools/generate_genhyperbolic_fisher_reference.py`, including GH/NIG and
  Fisher/hypergeometric identity checks plus skewed-tail stress cases.

- Generalized inverse-Gaussian (`geninvgauss`), normal-inverse-Gaussian
  (`norminvgauss`), and Skellam (`skellam`) families with eight-function APIs,
  likelihood/NNLF, analytic scores, bounded/fixed fitting, deterministic
  inverse-transform RNG, `scifort_stats`, and scalar C-ABI integration.
- A dependency-free real modified-Bessel-K logarithm and log-derivatives with
  respect to order and argument in `scifort_special`. Adaptive quadrature panel
  sizing preserves accuracy for large order and argument.
- Direct log-domain, locally truncated tail quadrature for generalized
  inverse-Gaussian and normal-inverse-Gaussian CDF/SF evaluation, including
  log-tail inversion for extreme quantiles.
- `test/test_gig_nig_skellam.f90`, generated SciPy 1.17.0 references, and
  `tools/generate_gig_nig_skellam_reference.py`, with additional independent
  high-precision stress references for deep tails and wide-order Bessel-K.

- Beta-negative-binomial (`betanbinom`), Yule-Simon (`yulesimon`), Zipf/zeta
  (`zipf`), and finite Zipfian (`zipfian`) discrete families with eight-function
  PMF APIs, likelihood/NNLF, analytic continuous-shape scores, bounded/fixed
  fitting, deterministic inverse-transform RNG, `scifort_stats`, and scalar
  C-ABI integration.
- A dependency-free real Hurwitz-zeta evaluator and derivative in
  `scifort_special`, implemented with an Euler-Maclaurin expansion and used for
  stable direct Zipf upper tails and the Zipf score.
- `test/test_betanbinom_yulesimon_zipf_zipfian.f90`, generated SciPy 1.17.0
  references, and `tools/generate_betanbinom_yulesimon_zipf_zipfian_reference.py`
  for reproducible eight-method regression coverage.

- Beta-binomial (`betabinom`), hypergeometric (`hypergeom`), negative
  hypergeometric (`nhypergeom`), and Boltzmann truncated-geometric (`boltzmann`)
  discrete families with eight-function PMF APIs plus likelihood/NNLF,
  analytic shape scores where differentiable, bounded/fixed fitting,
  deterministic inverse-transform RNG, `scifort_stats`, and scalar C-ABI
  integration. Hypergeometric and negative-hypergeometric shapes are integral
  and therefore expose zero-length continuous score vectors.
- `test/test_betabinom_hypergeom_nhypergeom_boltzmann.f90`, generated SciPy
  1.17.0 numerical references, and
  `tools/generate_betabinom_hypergeom_nhypergeom_boltzmann_reference.py` for
  reproducible full-method validation. Two stress-tail references are computed
  directly at high precision because SciPy 1.17.0 loses several digits there.
- Direct independently summed lower and upper log tails for the three finite
  combinatorial laws, and `expm1`-stable closed-form Boltzmann tails for very
  small rates.

- Discrete uniform (`randint`), Planck (`planck`), discrete Laplace (`dlaplace`),
  and logarithmic-series (`logser`) families with eight-function PMF APIs plus
  likelihood/NNLF, analytic shape scores where differentiable, bounded/fixed
  fitting, deterministic inverse-transform RNG, `scifort_stats`, and scalar
  C-ABI integration. `randint` has only integral shape bounds and therefore no
  continuous shape-score component.
- `test/test_randint_planck_dlaplace_logser.f90`, generated SciPy 1.17.0
  numerical references, and
  `tools/generate_randint_planck_dlaplace_logser_reference.py` for reproducible
  full-method validation of this discrete four-family slice.
- Stable direct tails for Planck and discrete Laplace and direct upper-tail
  recurrence for logarithmic-series probabilities, avoiding `1-cdf`
  cancellation in their unbounded discrete tails.

- Irwin-Hall (`irwinhall`), one-sided finite-sample Kolmogorov-Smirnov
  (`ksone`), noncentral chi-square (`ncx2`), and noncentral F (`ncf`)
  continuous families with eight-function APIs plus likelihood, analytic-score,
  bounded/fixed fitting, deterministic RNG, `scifort_stats`, and scalar C-ABI
  integration. Integer-valued shape parameters for Irwin-Hall and `ksone` are
  preserved as integral fit parameters and omitted from continuous scores.
- `test/test_irwinhall_ksone_ncx2_ncf.f90`, generated SciPy 1.17.0 numerical
  references, and `tools/generate_irwinhall_ksone_ncx2_ncf_reference.py` for
  reproducible full-method validation of this four-family slice.
- Stable cardinal-B-spline evaluation for Irwin-Hall, scaled exact Smirnov
  sums for one-sided KS tails/densities, and centered Poisson mixtures of
  regularized incomplete-gamma/beta kernels for noncentral chi-square/F.

- Double-Pareto lognormal (`dpareto_lognorm`), circular von Mises
  (`vonmises`), von Mises-on-a-line (`vonmises_line`), and asymptotic
  two-sided Kolmogorov (`kstwobign`) continuous families with eight-function
  APIs plus likelihood, analytic-score, bounded/fixed fitting, deterministic
  RNG, `scifort_stats`, and scalar C-ABI integration.
- `test/test_dpareto_vonmises_kstwobign.f90`, generated SciPy 1.17.0 numerical
  references, and `tools/generate_dpareto_vonmises_kstwobign_reference.py` for
  reproducible full-method validation of this four-family slice.
- Shared scaled modified-Bessel `i0e` and `i1e` kernels in `scifort_special`;
  the Rice implementation now reuses these kernels with the new von Mises
  family rather than carrying private copies.

- Generalized exponential (`genexpon`), skew-normal (`skewnorm`), Tukey lambda
  (`tukeylambda`), and Rice (`rice`) continuous families with eight-function
  APIs plus likelihood, analytic-score, bounded/fixed fitting, deterministic
  RNG, `scifort_stats`, and scalar C-ABI integration.
- `test/test_genexpon_skewnorm_tukeylambda_rice.f90`, generated SciPy 1.17.0
  numerical references, and
  `tools/generate_genexpon_skewnorm_tukeylambda_rice_reference.py` for
  reproducible full-method validation of this four-family slice.
- Dependency-free internal numerical kernels for the principal real Lambert-W
  branch used by generalized-exponential quantiles, Owen-T quadrature used by
  skew-normal tails, and scaled modified-Bessel ratios plus a Poisson/gamma
  mixture used by the Rice family.

- Crystal Ball (`crystalball`), Jones-Faddy skew-t (`jf_skew_t`), Pearson III
  (`pearson3`), and relativistic Breit-Wigner (`rel_breitwigner`) continuous
  families with eight-function APIs plus likelihood, analytic-score,
  bounded/fixed fitting, deterministic RNG, `scifort_stats`, and scalar C-ABI
  integration.
- `test/test_crystalball_jf_pearson_breitwigner.f90`, generated SciPy 1.17.0
  numerical references, and
  `tools/generate_crystalball_jf_pearson_breitwigner_reference.py` for
  reproducible full-method validation of this four-family slice.

- Generalized-gamma (`gengamma`), half-generalized-normal (`halfgennorm`),
  ARGUS (`argus`), and Erlang (`erlang`) continuous families with eight-function
  APIs plus likelihood, analytic-score, bounded/fixed fitting, deterministic RNG,
  `scifort_stats`, and scalar C-ABI integration.
- `test/test_gengamma_halfgennorm_argus_erlang.f90`, generated SciPy 1.17.0
  numerical references, and
  `tools/generate_gengamma_halfgennorm_argus_erlang_reference.py` for
  reproducible full-method validation of this four-family slice.

- Generalized extreme-value (`genextreme`), kappa-3 (`kappa3`), kappa-4
  (`kappa4`), and doubly truncated Weibull-minimum (`truncweibull_min`)
  continuous families with eight-function APIs plus likelihood, analytic-score,
  bounded/fixed fitting, deterministic RNG, `scifort_stats`, and scalar C-ABI
  integration.
- `test/test_genextreme_kappa_truncweibull.f90`, generated SciPy 1.17.0
  numerical references, and
  `tools/generate_genextreme_kappa_truncweibull_reference.py` for reproducible
  full-method validation of the four-family slice.

- Burr Type III (`burr`), Mielke beta-kappa/Dagum (`mielke`), Gibrat (`gibrat`),
  and wrapped Cauchy (`wrapcauchy`) continuous families with eight-function APIs
  plus likelihood, analytic-score, bounded/fixed fitting, deterministic RNG,
  `scifort_stats`, and scalar C-ABI integration.
- `test/test_burr_mielke_gibrat_wrapcauchy.f90`, generated SciPy 1.17.0 numerical
  references, and `tools/generate_burr_mielke_gibrat_wrapcauchy_reference.py` for
  reproducible full-method validation of the four-family slice.
- Exponentially modified normal (`exponnorm`), Johnson SB (`johnsonsb`),
  Johnson SU (`johnsonsu`), and trapezoid (`trapezoid`) continuous families
  with eight-function APIs plus likelihood, analytic-score, bounded/fixed
  fitting, deterministic RNG, `scifort_stats`, and scalar C-ABI integration.
- `test/test_normal_transform_four.f90`, generated SciPy 1.17.0 numerical
  references, and `tools/generate_normal_transform_reference.py` for
  reproducible validation of the four-family normal/transform slice.
- Folded normal (`foldnorm`), folded Cauchy (`foldcauchy`), reciprocal
  inverse-Gaussian (`recipinvgauss`), and truncated Pareto (`truncpareto`)
  continuous families with eight-function APIs plus likelihood, analytic-score,
  bounded/fixed fitting, deterministic RNG, `scifort_stats`, and scalar C-ABI
  integration.
- `test/test_folded_reciprocal_four.f90`, generated SciPy 1.17.0 numerical
  references, and `tools/generate_folded_reciprocal_reference.py` for the
  four-family folded/reciprocal regression slice.
- Double-gamma (`dgamma`), asymmetric Laplace (`laplace_asymmetric`), finite-bound
  truncated normal (`truncnorm`), and log-uniform/reciprocal (`loguniform`)
  continuous families with eight-function APIs plus likelihood, analytic-score,
  bounded/fixed fitting, deterministic RNG, `scifort_stats`, and scalar C-ABI
  integration.
- `test/test_scipy_transform_four.f90`, generated SciPy 1.17.0 numerical
  references, and `tools/generate_scipy_transform_reference.py` for the
  four-family regression slice.
- Left-Levy (`levy_l`), Weibull maximum (`weibull_max`), R-distribution
  (`rdist`), and skew-Cauchy (`skewcauchy`) continuous families with
  eight-function APIs plus likelihood, analytic-score, bounded/fixed fitting,
  deterministic RNG, `scifort_stats`, and scalar C-ABI integration.
- `test/test_reflected_shape_four.f90` with SciPy 1.17.0 PDF/CDF/PPF and
  logarithmic-tail references, score finite differences, deterministic RNG,
  fixed-fit, and C-ABI checks for the reflected/shape slice.
- Generalized half-logistic (`genhalflogistic`), exponential-power
  (`exponpow`), exponentiated-Weibull (`exponweib`), and power-lognormal
  (`powerlognorm`) continuous families with eight-function APIs plus
  likelihood, analytic-score, bounded/fixed fitting, deterministic RNG,
  `scifort_stats`, and scalar C-ABI integration.
- `test/test_transform_shape_four.f90` with SciPy 1.17.0 PDF/CDF/PPF
  references, score finite differences, deterministic RNG, fixed-fit, and
  C-ABI checks for the four-family transform/shape slice.
- Gompertz (`gompertz`), inverse-Weibull (`invweibull`), beta-prime
  (`betaprime`), and Burr XII (`burr12`) continuous families with eight-function
  APIs plus likelihood, analytic-score, bounded/fixed fitting, deterministic
  RNG, `scifort_stats`, and scalar C-ABI integration.
- `test/test_positive_shape_continuous.f90`, SciPy 1.17.0 numerical references,
  and `tools/generate_positive_shape_reference.py` for reproducible validation
  of the four-family positive-support slice.
- Nakagami (`nakagami`), power-normal (`powernorm`), log-gamma (`loggamma`),
  and Wald (`wald`) continuous families with eight-function APIs plus
  likelihood, analytic-score, bounded/fixed fitting, deterministic RNG,
  `scifort_stats`, and scalar C-ABI integration.
- `test/test_transform_family_continuous.f90`, SciPy 1.17.0 reference values,
  and `tools/generate_transform_family_reference.py` for reproducible
  validation of the four-family transform/shape slice.
- Alpha (`alpha`), Birnbaum-Saunders/fatigue-life (`fatiguelife`), generalized
  logistic (`genlogistic`), and generalized normal (`gennorm`) continuous
  families with eight-function APIs plus likelihood, analytic-score,
  bounded/fixed fitting, deterministic RNG, `scifort_stats`, and scalar C-ABI
  integration.
- `test/test_shape_family_continuous.f90`, SciPy 1.17.0 numerical references,
  and `tools/generate_shape_family_reference.py` for reproducible validation of
  the four-family one-shape slice.
- Bradford (`bradford`), truncated exponential (`truncexpon`), Fisk/log-logistic
  (`fisk`), and double-Weibull (`dweibull`) continuous families with eight-function
  distribution APIs plus likelihood, analytic-score, bounded/fixed fitting,
  deterministic RNG, `scifort_stats`, and scalar C-ABI integration.
- `test/test_closed_form_continuous.f90`, SciPy 1.17.0 numerical references, and
  `tools/generate_closed_form_reference.py` for reproducible parity checks of the
  four-family closed-form slice.
- Inverse-gamma (`invgamma`), inverse-Gaussian (`invgauss`), Levy (`levy`),
  and log-Laplace (`loglaplace`) continuous families with eight-function
  distribution APIs plus likelihood, analytic-score, bounded/fixed fitting,
  deterministic RNG, `scifort_stats`, and scalar C-ABI integration.
- `test/test_inverse_positive.f90` with SciPy 1.17.0 PDF/CDF/PPF references,
  tail round trips, score finite differences, RNG, fit, and C-ABI checks.
- Hyperbolic-secant (`hypsecant`) and half-logistic (`halflogistic`) continuous
  families, including eight-function distribution APIs plus likelihood,
  analytic-score, fitting, RNG, `scifort_stats`, and scalar C-ABI integration.
- Anglit (`anglit`) and Moyal (`moyal`) continuous families, including the
  eight-function distribution API plus likelihood, analytic-score, fitting,
  RNG, `scifort_stats`, and scalar C-ABI integration.
- Cosine (`cosine`) and semicircular (`semicircular`) bounded continuous
  families, including eight-function distribution APIs plus likelihood,
  analytic-score, fitting, RNG, `scifort_stats`, and scalar C-ABI integration.
- Milestone 7 continuous distributions: right Gumbel (`gumbel_r`), left Gumbel
  (`gumbel_l`), power-function (`powerlaw`), triangular (`triang`), and
  generalized Pareto (`genpareto`), each with PDF/log-PDF, CDF/SF, logarithmic
  tails, PPF, and ISF.
- Full likelihood/NNLF, analytic-score, bounded/fixed fitting, explicit-state
  random-variate, and scalar C-ABI integration for all five new families.
- `test/test_extended_continuous.f90`, SciPy 1.17.0 reference tables, and
  `tools/generate_extended_continuous_reference.py` for reproducible Milestone 7
  parity checks.
- Arcsine (`arcsine`) and half-normal (`halfnorm`) continuous distributions,
  including likelihood/NNLF, analytic scores, bounded/fixed fitting,
  explicit-state random variates, scalar C ABI, and `test_more_continuous`.
- Half-Cauchy (`halfcauchy`) and Lomax/Pareto-II (`lomax`) continuous
  distributions with the same eight-function distribution API and full
  likelihood/score/fit/RNG/scalar-C-ABI integration.
- Chi (`chi`, positive real degrees of freedom) and Maxwell (`maxwell`)
  continuous distributions with the same eight-function API and full
  likelihood/score/fit/RNG/scalar-C-ABI integration.
- Regularized incomplete gamma functions `gammainc` and `gammaincc` and
  their inverses `gammaincinv` and `gammainccinv` in the new
  `scifort_special` module.
- Gamma distribution (`gamma_*`, shape `a` with `loc` and `scale`) and
  chi-square distribution (`chi2_*`, `df` with `loc` and `scale`), each with
  `pdf`, `logpdf`, `cdf`, `sf`, `logcdf`, `logsf`, `ppf`, and `isf`.
- High-precision reference data generated by
  `tools/generate_gamma_reference.py` with conditioning-based tolerances,
  and the `test/test_gamma.f90` test program.
- Regularized incomplete beta functions `betainc` and `betaincc` and their
  inverses `betaincinv` and `betainccinv` in `scifort_special`.
- Beta (`beta_*`), Student t (`t_*`, including `df = +infinity`), and F
  (`f_*`) distributions with the same eight functions.
- `tools/generate_beta_reference.py` and the `test/test_beta.f90` test
  program.
- Poisson (`poisson_*`) and binomial (`binomial_*`) distributions with
  `pmf`, `logpmf`, `cdf`, `sf`, `logcdf`, `logsf`, `ppf`, and `isf`. Counts
  may be given as `integer` or `real(dp)` through generic interfaces.
- `tools/generate_discrete_reference.py` and the `test/test_discrete.f90`
  test program.
- Bernoulli (`bernoulli_*`), geometric (`geometric_*`), and negative-binomial
  (`negative_binomial_*`) distributions with the same eight-function discrete
  API. Negative-binomial shape `n` may be any finite positive real value.
- High-precision reference data generated by
  `tools/generate_discrete_extended_reference.py` and the
  `test/test_discrete_extended.f90` test program.
- Lognormal (`lognormal_*`), Weibull (`weibull_*`), Pareto (`pareto_*`), and
  Rayleigh (`rayleigh_*`) distributions, each with `pdf`, `logpdf`, `cdf`,
  `sf`, `logcdf`, `logsf`, `ppf`, and `isf`.
- High-precision reference data generated by
  `tools/generate_continuous_reference.py` and the `test/test_continuous.f90`
  test program.
- Scalar C ABI entry points for PDF or PMF, CDF, and PPF across every
  currently implemented distribution, with matching declarations in
  `include/scifort.h` and regression coverage in `test/test_c_api.f90`.
- Explicit-state `scifort_random` generator API with deterministic seeding,
  state serialization/restoration, scalar open-interval uniforms, and bulk
  uniform filling. The recurrence is MRG32k3a; binary64 uniforms combine two
  unbiased 26-bit buckets into a 52-bit open grid.
- Scalar `*_rvs` functions and `*_rvs_array` subroutines for every current
  distribution. They use the existing PPF kernels so sampling follows the same
  parameter validation, support, and quantile conventions as the distribution
  API.
- `test/test_random.f90` covering generator reproducibility, state management,
  uniform smoke statistics, all distribution samplers, and invalid-parameter
  state preservation.
- `scifort_likelihood` with array `*_loglikelihood` and `*_nnlf`
  evaluators for all current distributions. Discrete families provide generic
  integer- and real-observation interfaces matching their PMF APIs.
- `test/test_likelihood.f90` covering all likelihood families, discrete
  overloads, empty samples, invalid parameters, and the `nnlf = -loglike`
  convention.
- `scifort_score` with analytic summed-log-likelihood scores for all 20
  distributions, plus an internal positive-argument digamma implementation for
  shape derivatives. Discrete scores provide integer and real-observation
  overloads and omit nondifferentiable integer support/trial parameters.
- `test/test_score.f90` checking digamma reference values, analytic scores
  against independent central differences, discrete overload equivalence, and
  nonregular support-boundary behavior.
- `scifort_fit` with deterministic bounded maximum-likelihood interfaces for
  all 20 distributions. Equal bounds fix parameters, discrete integral
  parameters are projected to integers, and `fit_result` reports objective,
  convergence status, iteration count, and evaluation count.
- `test/test_fit.f90` covering closed-form MLE checks, fixed-parameter fits,
  integral discrete parameters, all family wrappers, overloads, and failure
  diagnostics.
- `scifort_descriptive` with numerically scaled mean, configurable-`ddof`
  variance and standard deviation, central moments, scalar/array linear
  quantiles, median, covariance, Pearson correlation, and `rankdata` methods
  `average`, `min`, `max`, `dense`, and `ordinal`.
- `scifort_hypothesis` with one-sample, pooled/Welch independent, and paired
  t-tests; exact-null Pearson correlation p-values, Spearman rank-correlation
  p-values, and asymptotic Mann-Whitney U tests with tie and continuity
  corrections. All tests accept `two-sided`, `less`, and `greater` alternatives.
- `scifort_resampling` with generic percentile bootstrap confidence intervals
  and randomized two-sample permutation tests over pure statistic callbacks.
  Both use explicit `rng_state`; finite-grid index draws are unbiased by
  rejection and invalid calls do not advance the generator.
- `test/test_descriptive.f90`, `test/test_hypothesis.f90`, and
  `test/test_resampling.f90`, plus `tools/generate_statistics_reference.py`
  for reproducible NumPy/SciPy comparison values.
- Milestone 6 `scifort_special` scalar functions: `erf`, `erfc`, `erfinv`,
  `erfcinv`, `gammaln`, positive-shape `betaln`, `digamma`/`psi`, `ndtr`,
  `log_ndtr`, `ndtri`, `expit`, `logit`, `log_expit`, `xlogy`, `xlog1py`,
  `entr`, `rel_entr`, `kl_div`, `exprel`, `cosm1`, and Box-Cox forward and
  inverse transforms.
- One-dimensional `logsumexp`, `softmax`, and `log_softmax`, plus
  `test/test_special_elementary.f90` and `tools/generate_special_reference.py`.

### Changed

- Every dummy argument in `src/` and the test programs now carries a
  trailing FORD documentation comment (`!!`) describing its meaning, range,
  and default. `AGENTS.md` requires these for new code.

### Fixed

- `source_order_generated.txt` now includes the already-present `kstwo`,
  Wallenius noncentral-hypergeometric, and Poisson-binomial modules as well as
  the new Levy-stable and studentized-range modules, so the listed order covers
  the complete source graph.

- `log_beta` now evaluates `log(Gamma(small))` through the existing
  `log_gamma_one_plus(small) - log(small)` path when a beta shape is below
  `0.5` and the other shape is in the Stirling region, avoiding reliance on
  compiler behavior arbitrarily close to the gamma pole at zero.
- `log1p_safe` lost up to about `5e-13` relative accuracy for arguments
  near `1e-4`, and `expm1_safe` up to about `3e-12` near `1e-5`, because
  they formed `log(1 + x)` and `exp(x) - 1` directly outside a very small
  series region. Both now use series forms that keep the relative error
  below `6e-16` without relying on `(1 + x) - 1` being evaluated exactly.
  This improves the exponential, logistic, Laplace, and normal routines that
  use them.

## 0.1.0 - 2026-08-20

### Added

- Bradford (`bradford`), truncated exponential (`truncexpon`), Fisk/log-logistic
  (`fisk`), and double-Weibull (`dweibull`) continuous families with eight-function
  distribution APIs plus likelihood, analytic-score, bounded/fixed fitting,
  deterministic RNG, `scifort_stats`, and scalar C-ABI integration.
- `test/test_closed_form_continuous.f90`, SciPy 1.17.0 numerical references, and
  `tools/generate_closed_form_reference.py` for reproducible parity checks of the
  four-family closed-form slice.
- Cosine (`cosine`) and semicircular (`semicircular`) bounded continuous
  families, including eight-function distribution APIs plus likelihood,
  analytic-score, fitting, RNG, `scifort_stats`, and scalar C-ABI integration.
- FPM static and shared library build.
- MIT license and provenance policy.
- Repository instructions for Codex and Claude Code.
- Normal, uniform, exponential, Laplace, logistic, and Cauchy distributions.
- Scalar and elemental-array Fortran API.
- Experimental scalar and vector C ABI for normal PDF, CDF, and PPF.
- Tests, example program, C header, and continuous integration workflow.
- Direct logarithmic Cauchy tails that avoid cancellation near one.
- A documented validation record for the starter archive.
