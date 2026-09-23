# Code provenance

This file records the origin and licensing of numerical implementation code.
Every imported, translated, or substantially adapted routine must be entered
before it is merged.

## Current implementation inventory

Original SciFort implementation code is MIT licensed. The unreleased additions
also contain the explicitly identified SciPy BSD-3-Clause adaptations and Sobol
data below. Mixed-source files retain both MIT and BSD-3-Clause notices;
third-party notices and complete license texts are recorded in
`THIRD_PARTY_LICENSES.md` and `LICENSES/`. The Landau research coefficient-table
reference is separately attributed under CC BY 4.0.

- `src/math.f90`: stable elementary helpers implemented from standard
  power-series and algebraic identities. `floor_real` implements the mathematical
  floor using real truncation and a negative-fraction correction, without an
  intermediate integer conversion; nonfinite values pass through unchanged.
- `src/stats/normal.f90`: standard normal `erfc` identity, a
  Mills-ratio asymptotic expansion for log-CDF tails, and safeguarded
  bisection for quantiles.
- `src/stats/uniform.f90`: direct implementation from the mathematical
  definition.
- `src/stats/exponential.f90`: direct implementation with stable
  `expm1` and `log1p` helpers.
- `src/stats/laplace.f90`: direct piecewise implementation from the
  mathematical definition.
- `src/stats/logistic.f90`: direct implementation using stable
  logistic and softplus identities.
- `src/stats/cauchy.f90`: direct implementation using reciprocal-angle
  identities and cancellation-resistant logarithmic tails.
- `src/c_api/api.f90`: interoperable wrapper layer over the native
  Fortran API.
- `src/special/incomplete_gamma.f90`: regularized incomplete gamma
  functions and inverses. Independently implemented from these
  specifications, not from any library source:
  - NIST Digital Library of Mathematical Functions (DLMF),
    https://dlmf.nist.gov/, equations 8.7.1 (series for `P`), 8.7.3 (series
    used for `Q` when `a < 1` and `x <= 1.5`), 8.9.2 (continued fraction for
    `Q`), 5.7.3 (series for `log(Gamma(1 + a))`), 5.11.1 (Stirling series),
    and 4.6.1 (series for `log(1 + t) - t`);
  - the modified Lentz algorithm for continued fractions: W. J. Lentz,
    Applied Optics 15(3) (1976) 668, doi:10.1364/AO.15.000668;
    I. J. Thompson and A. R. Barnett, Journal of Computational Physics 64
    (1986) 490-509, doi:10.1016/0021-9991(86)90046-X;
  - the Wilson-Hilferty approximation for inverse starting values:
    E. B. Wilson and M. M. Hilferty, "The Distribution of Chi-Square",
    Proceedings of the National Academy of Sciences 17 (1931) 684-688,
    doi:10.1073/pnas.17.12.684.

  The region selection, safeguarded logarithmic Newton iteration, and
  tolerances are original. The 29 constants `zeta(k) - 1`, `k = 2..30`, were
  generated with mpmath 1.3.0 at 40 digits and checked against the closed
  forms of `zeta(2)` and `zeta(4)`. The Stirling coefficients are the exact
  rationals `B(2k) / (2k (2k - 1))` for `k = 1..7`, checked with mpmath.
- `src/special/log_gamma.f90`: log-gamma and log-beta helpers shared by the
  special functions, moved from `incomplete_gamma.f90` and extended with
  `log_beta`, `log_gamma_difference`, and `log_gamma_ratio_scaled`, all
  independently derived from DLMF 5.7.3 and 5.11.1 (Stirling series with the
  large terms cancelled analytically).
- `src/special/incomplete_beta.f90`: regularized incomplete beta function and
  inverses. Independently implemented from DLMF 8.17.4 (symmetry), 8.17.7
  (series used for the small-parameter complement), 8.17.22 and 8.17.23
  (continued fraction, evaluated by the modified Lentz method cited above),
  and 5.11.1. The expansion in incomplete gamma functions for one large and
  one small parameter (`gamma_expansion`) was derived for this project from
  the substitution `t = exp(-v)` in the integral for `B_x(a, b)`, the
  factorization `1 - exp(-v) = v exp(-v/2) sinh(v/2) / (v/2)`, the
  Bernoulli generating function of DLMF 24.2.1, and the recurrence of DLMF
  8.8.2. No implementation of a similar expansion was consulted. The 30
  coefficients `B(2n) / (2n (2n)!)` were generated with mpmath 1.3.0 at 40
  digits and checked against `log(sinh(v/2) / (v/2))` to `1e-40` at `v = 1`.
- `src/stats/beta.f90`, `src/stats/student_t.f90`,
  `src/stats/f_distribution.f90`: direct implementations from the density
  definitions over the incomplete beta function.
- `src/stats/poisson.f90`: direct implementation using DLMF 8.4.8 with
  8.4.11, which express the finite Poisson sum through the incomplete gamma
  functions, and the kernel of `src/special/incomplete_gamma.f90`.
- `src/stats/binomial.f90`: direct implementation using DLMF 8.17.5, which
  expresses the binomial tail as a regularized incomplete beta function, and
  the kernel of `src/special/incomplete_beta.f90`. The quantile search
  (a normal starting approximation with a continuity correction, a doubling
  bracket, and bisection over integers) is original in both modules.
- `src/stats/bernoulli.f90`: exact specialization of the binomial distribution
  to one trial, implemented by delegating to the SciFort binomial kernel.
- `src/stats/geometric.f90`: direct implementation from the mathematical
  definition. Stable CDF and survival functions use the project-local
  `log1p_safe` and `expm1_safe` helpers. The analytic logarithmic quantile and
  integer correction are original SciFort code.
- `src/stats/negative_binomial.f90`: direct implementation from the mass
  definition and DLMF 8.17.5, which gives the CDF as a regularized incomplete
  beta function. The log-mass rearrangement through the existing beta kernel,
  and the normal-start/doubling/bisection quantile search, are original
  SciFort code.
- `src/special.f90`: re-export module for special functions.
- `src/stats/gamma.f90`, `src/stats/chi2.f90`: direct implementations from
  the density definitions over the incomplete gamma functions.
- `src/stats/lognormal.f90`: direct implementation using standard normal
  distribution functions over standardized log-transformed coordinates.
- `src/stats/weibull.f90`: direct implementation from the mathematical definition
  of the Weibull minimum distribution with stable logarithmic tails and `expm1_safe`.
- `src/stats/pareto.f90`: direct implementation from the Pareto Type I definition
  with stable logarithmic tails and `expm1_safe`.
- `src/stats/rayleigh.f90`: direct implementation from the Rayleigh distribution
  definition with stable logarithmic tails and `expm1_safe`.
- `src/stats/gumbel_r.f90` and `src/stats/gumbel_l.f90`: original SciFort
  direct implementations from the standard Gumbel density/CDF definitions.
  Stable logarithmic and inverse-tail rearrangements were derived for this
  project. SciPy documentation was consulted only for public parameterization
  and numerical parity; no SciPy source was copied or adapted.
- `src/stats/powerlaw.f90`: original SciFort direct implementation of the
  power-function density `a*x**(a-1)` on standardized `[0,1]`. The SciPy
  `scipy.stats.powerlaw` documentation was used to confirm parameterization and
  endpoint conventions; no implementation source was copied or translated.
- `src/stats/triang.f90`: original SciFort piecewise implementation of the
  triangular density/CDF on standardized `[0,1]`, including the `c=0` and
  `c=1` endpoint modes and cancellation-resistant inverse formulas. SciPy's
  `scipy.stats.triang` documentation was used only as an API/parity reference.
- `src/stats/genpareto.f90`: original SciFort implementation of the generalized
  Pareto density `(1+c*x)**(-1-1/c)` and its exponential `c=0` limit. Direct
  logarithmic tails and quantiles use project-local `log1p_safe`/`expm1_safe`.
  SciPy's `scipy.stats.genpareto` documentation was used only to confirm public
  parameterization/support and for numerical validation; no SciPy source code
  was copied, translated, or adapted.
- `src/stats/arcsine.f90`: original SciFort implementation from the standard
  arcsine density/CDF identities on `[0,1]`, with complementary-tail formulas
  chosen directly from the same identities. SciPy 1.17.0 was used only for
  numerical/API parity checks; no SciPy source was copied or translated.
- `src/stats/halfnorm.f90`: original SciFort implementation from the standard
  half-normal density and error-function CDF, reusing SciFort's normal tail and
  quantile kernels. SciPy 1.17.0 was used only for numerical/API parity checks;
  no SciPy source was copied or translated.
- `src/stats/halfcauchy.f90`: original SciFort implementation from the
  standard half-Cauchy density and arctangent CDF identities, with direct
  complementary-tail formulas for numerical stability. SciPy 1.17.0 was used
  only for numerical/API parity checks; no SciPy source was copied or translated.
- `src/stats/lomax.f90`: original SciFort implementation of the Pareto-II/Lomax
  density and survival identities, using project-local `log1p_safe` and
  `expm1_safe` helpers for stable tails and inverse transforms. SciPy 1.17.0 was
  used only for numerical/API parity checks; no SciPy source was copied or translated.
- `src/stats/chi.f90`: original SciFort implementation using the standard chi
  density and the regularized incomplete-gamma identities for the CDF,
  survival function, and quantiles. SciPy 1.17.0 was used as an API and
  numerical-validation reference; no SciPy implementation source was copied.
- `src/stats/maxwell.f90`: original SciFort specialization of the chi
  distribution at three degrees of freedom, matching the standard Maxwell
  parameterization used by `scipy.stats.maxwell`.
- `src/stats/cosine.f90`: original SciFort implementation of the bounded
  cosine density and distribution formulas used by `scipy.stats.cosine`;
  quantiles use deterministic bisection on the monotone closed-form CDF.
- `src/stats/semicircular.f90`: original SciFort implementation of the
  Wigner semicircle density/CDF used by `scipy.stats.semicircular`; quantiles
  use deterministic bisection on the monotone closed-form CDF.
- `src/stats/anglit.f90`: original SciFort implementation of the standard
  Anglit density, trigonometric CDF, and inverse-CDF identities. SciPy 1.17.0
  was used for numerical and API parity checks; no SciPy implementation code
  was copied or translated.
- `src/stats/moyal.f90`: original SciFort implementation of the standard Moyal
  density and error-function tail identities, with explicit log-tail handling
  and normal-quantile inversion for stable probabilities. SciPy 1.17.0 was used
  for numerical and API parity checks; no SciPy implementation code was copied
  or translated.
- `src/stats/hypsecant.f90`: original SciFort implementation from the standard
  hyperbolic-secant density `sech(x)/pi`, its arctangent CDF, and analytic
  inverse-CDF identity. Stable log tails and `log(cosh(x))` are independently
  rearranged to avoid overflow. SciPy 1.17.0 was used only for public
  parameterization and numerical parity checks; no SciPy source was copied.
- `src/stats/invgamma.f90`: original SciFort implementation from the standard
  inverse-gamma density and regularized incomplete-gamma identities. SciPy
  1.17.0 was used only for public parameterization and numerical parity checks;
  no SciPy implementation source was copied or translated.
- `src/stats/invgauss.f90`: original SciFort implementation of the
  inverse-Gaussian density and normal-CDF representation, with deterministic
  tail-aware bisection for quantiles. SciPy 1.17.0 was used only for numerical
  and API parity checks; no SciPy source was copied or translated.
- `src/stats/levy.f90`: original SciFort implementation from the standard Levy
  density, error-function CDF, and inverse-normal/error-function quantile
  identities. No third-party implementation code was copied or translated.
- `src/stats/loglaplace.f90`: original SciFort implementation from the
  piecewise power-law form of the log-Laplace density/CDF and its analytic
  inverse. SciPy 1.17.0 was used only for parameterization and parity checks.
- `src/stats/halflogistic.f90`: original SciFort implementation from the
  standard half-logistic density and hyperbolic-tangent CDF identities. Direct
  survival/log-survival and inverse-survival formulas avoid cancellation in the
  upper tail. SciPy 1.17.0 was used only for public parameterization and
  numerical parity checks; no SciPy source was copied.

- `src/random.f90`: clean implementation of the MRG32k3a combined multiple
  recursive generator from P. L'Ecuyer, "Good Parameters and Implementations
  for Combined Multiple Recursive Random Number Generators", Operations
  Research 47(1) (1999) 159-164, doi:10.1287/opre.47.1.159. The published
  moduli and recurrence coefficients were transcribed from the paper; no
  software implementation was copied or translated. The scalar-seed mapping,
  state API, rejection into equal 26-bit buckets, and two-chunk 52-bit
  binary64 conversion are original SciFort code.
- `src/random_variates.f90`: original SciFort inverse-transform sampling
  wrappers. They combine `src/random.f90` with the existing distribution PPF
  routines and do not incorporate code from SciPy or another RNG library.
- `src/stats/bradford.f90`, `src/stats/truncexpon.f90`, `src/stats/fisk.f90`,
  and `src/stats/dweibull.f90`: original SciFort implementations from the standard
  mathematical definitions of the Bradford, finite truncated-exponential,
  log-logistic/Fisk, and symmetric double-Weibull distributions. Stable tail and
  inverse formulas were independently rearranged using the project-local
  `log1p_safe`, `expm1_safe`, and `log1pexp` helpers. SciPy 1.17.0 was used only
  as an independent numerical validation reference; no SciPy source code was
  copied, translated, or adapted.
- `src/stats/alpha.f90`, `src/stats/fatiguelife.f90`,
  `src/stats/genlogistic.f90`, and `src/stats/gennorm.f90`: original SciFort
  implementations from the published mathematical definitions of the Alpha,
  Birnbaum-Saunders, generalized-logistic, and generalized-normal families.
  The generalized-normal CDF and inverse reuse SciFort's independently
  implemented regularized incomplete-gamma kernel. Stable tail and inverse
  rearrangements were derived independently. SciPy 1.17.0 was used only for
  public parameterization and numerical validation; no SciPy source code was
  copied, translated, or adapted.
- `src/likelihood.f90`: original SciFort reductions over the existing
  distribution `logpdf`/`logpmf` routines. The family-specific
  `*_loglikelihood` and `*_nnlf` interfaces contain no imported numerical
  algorithm or translated third-party implementation.
- `src/special/digamma.f90`: original SciFort positive-argument digamma
  helper, independently implemented from DLMF 5.5.2 (recurrence) and DLMF
  5.11.2 (asymptotic expansion). The Bernoulli-number coefficients are exact
  rational constants. No third-party implementation was copied or translated.
- `src/score.f90`: original SciFort analytic derivatives of the existing
  distribution log-density and log-mass formulas. Shape derivatives use the
  project-local positive-argument digamma helper. No SciPy or third-party
  implementation was copied, translated, or adapted.
- `src/fit.f90`: original SciFort deterministic projected coordinate/pattern
  search and family-specific fitting wrappers over `scifort_likelihood`.
  Equal-bound fixed parameters and integer-constrained fit coordinates follow
  the documented public fitting conventions used by SciPy, but the optimizer
  itself is independently implemented and does not copy or translate SciPy's
  default differential-evolution implementation or another optimizer.
- `src/math.f90` (`log1p_safe`, `expm1_safe`): series from DLMF 4.6.4 and
  4.2.19.

- `src/descriptive.f90`: original SciFort implementations of standard
  descriptive-statistics definitions. The quantile interpolation is the
  Hyndman-Fan type 7 convention from R. J. Hyndman and Y. Fan, "Sample
  Quantiles in Statistical Packages", The American Statistician 50(4)
  (1996) 361-365, doi:10.1080/00031305.1996.10473566. Scaled compensated
  summation, scaled centered powers, merge sorting, and rank grouping were
  implemented independently for this project; no NumPy, SciPy, or Fortran
  stdlib implementation was copied or translated.
- `src/hypothesis.f90`: original SciFort implementations from the standard
  definitions of Student/Welch t statistics, the exact beta null law of the
  Pearson correlation coefficient, the rank-transform approximation for
  Spearman correlation, and the tie-corrected normal approximation to the
  Mann-Whitney U statistic. The Mann-Whitney test traces to H. B. Mann and
  D. R. Whitney, "On a Test of Whether one of Two Random Variables is
  Stochastically Larger than the Other", Annals of Mathematical Statistics
  18(1) (1947) 50-60. SciPy 1.17.0 was used only as an independent numerical
  validation reference; no SciPy source code was copied or adapted.
- `src/hypothesis_extended.f90`: original Fortran implementations of standard
  rank, contingency, analysis-of-variance, and scale-test formulas. SciPy
  1.17.0 source and documentation were consulted for public argument/default
  semantics, Wilcoxon method-selection behavior, named Cressie-Read lambda
  conventions, and edge-case validation. No SciPy implementation code or
  tables were copied or translated. Exact signed-rank probabilities use an
  independently implemented subset-sum dynamic program; small tied/zero
  signed-rank cases use exhaustive sign enumeration.
- `src/contingency_meta.f90`: original SciFort Fortran implementations of exact
  binomial testing, Barnard/Boschloo unconditional 2-by-2 tests, odds ratio,
  relative risk, contingency association, p-value combination, and FDR control.
  SciPy 1.17.0 `_binomtest.py`, `_hypotests.py`, `_odds_ratio.py`,
  `_relative_risk.py`, `contingency.py`, `_stats_py.py`, and `_morestats.py`
  were consulted for public defaults, formulas, numerical tolerances, and edge-
  case semantics. No SciPy implementation source was copied or translated.
  Barnard/Boschloo nuisance maximization is independently implemented as a
  deterministic grid search plus golden-section refinement rather than SHGO.
- `src/multiple_comparisons.f90`: original SciFort Fortran implementations of
  the Alexander-Govern approximation, Tukey HSD/Tukey-Kramer, Games-Howell,
  Dunnett's single-step treatment-control procedure, and the
  Krishnamoorthy-Thomson Poisson means E-test. SciPy 1.17.0 `_stats_py.py`,
  `_hypotests.py`, and `_multicomp.py` were consulted for public defaults,
  formulas, numerical conventions, and edge semantics. The implementation
  reuses SciFort's independently implemented studentized-range, multivariate-t,
  Poisson, and chi-square kernels; no SciPy implementation source was copied
  verbatim. Dunnett confidence limits use an independently implemented
  monotone bisection against the SciFort multivariate-t box probability rather
  than SciPy's noisy Brent minimization.
- `test/test_multiple_comparisons.f90` and
  `test/multiple_comparisons_reference.f90`: SciFort-owned regression logic and
  public SciPy 1.17.0 numerical outputs regenerated by
  `tools/generate_multiple_comparisons_reference.py`.
- `test/test_contingency_meta.f90` and
  `test/contingency_meta_reference.f90`: SciFort-owned regression logic and
  public SciPy 1.17.0 numerical outputs regenerated by
  `tools/generate_contingency_meta_reference.py`.
- `src/nonparametric_extended.f90`: original SciFort implementations of the
  published Ansari-Bradley, Mood, Epps-Singleton, k-sample Anderson-Darling,
  and median-test formulas. SciPy 1.17.0 `_morestats.py` and `_hypotests.py`
  were consulted for public defaults, tie conventions, finite/asymptotic method
  selection, correction factors, critical values, and edge-case semantics; no
  SciPy implementation code was copied or translated. The exact small-sample
  Ansari-Bradley distribution uses an independently implemented subset-count
  dynamic program. Epps-Singleton uses SciFort's existing symmetric eigensolver
  for the covariance pseudoinverse.
- `test/test_nonparametric_extended.f90` and
  `test/nonparametric_extended_reference.f90`: SciFort-owned regression logic
  and public SciPy 1.17.0 numerical outputs; the reference module is regenerated
  by `tools/generate_nonparametric_extended_reference.py`.
- `src/association_extended.f90`: SciFort Fortran implementations of Kendall
  tau, point-biserial correlation, ordinary and robust line fitting, the
  Brunner-Munzel test, and Page's ordered-alternative test. Public defaults,
  tie conventions, method selection, and finite-sample formulas were checked
  against SciPy 1.17.0. The exact Kendall inversion-count recurrence is adapted
  from SciPy's `_kendall_p_exact` implementation; the applicable SciPy
  BSD-3-Clause license is reproduced in `THIRD_PARTY_LICENSES.md`. Page exact
  probabilities are independently implemented by convolution of the single-row
  permutation distribution.
- `test/test_association_extended.f90` and
  `test/association_extended_reference.f90`: SciFort-owned regression logic and
  public SciPy 1.17.0 numerical outputs regenerated by
  `tools/generate_association_extended_reference.py`.
- `src/resampling.f90`: original SciFort explicit-state resampling code.
  Percentile/basic bootstrap intervals follow the ordinary nonparametric
  bootstrap framework, and BCa uses the standard bias-correction plus
  leave-one-out jackknife acceleration formulas described by Efron and
  Tibshirani. Exact permutation enumeration, randomized Fisher-Yates shuffles,
  Monte Carlo p-values, and simulated power were independently implemented.
  SciPy 1.17.0 documentation/source behavior was used to verify exact-switching,
  the `independent`/`samples`/`pairings` permutation definitions, plus-one
  randomized p-values, the `100*epsilon` comparison tolerance, one-sided
  bootstrap conventions, and the strict `p < significance` power definition.
  Selection from the finite 52-bit RNG grid uses an independently implemented
  rejection step to remove modulo bias. No third-party resampling source was
  copied or translated.

- `test/resampling_inference_reference.f90` and
  `tools/generate_resampling_inference_reference.py`: deterministic public
  SciPy 1.17.0 outputs for exact permutation, Monte Carlo, and power regression
  cases; these files contain numerical references, not SciPy implementation
  source.

- `src/special/elementary.f90`: original SciFort implementation. `erf` and
  `erfc` are thin wrappers around the Fortran standard intrinsics. The inverse
  error functions use the identities with the standard-normal quantile;
  `gammaln` and negative-argument `digamma` use the classical gamma/digamma
  reflection formulas (DLMF 5.5), while positive digamma delegates to the
  existing independently implemented recurrence/asymptotic kernel. Logistic,
  entropy, `exprel`, and `cosm1` routines follow their mathematical
  definitions with independently chosen stable rearrangements. Box-Cox and
  inverse Box-Cox functions follow G. E. P. Box and D. R. Cox, "An Analysis
  of Transformations", Journal of the Royal Statistical Society Series B
  26(2) (1964) 211-252, using the project-local `log1p_safe` and
  `expm1_safe` kernels. SciPy 1.17.0 was used only for numerical validation;
  no SciPy implementation code was copied, translated, or adapted.
- `src/special/reductions.f90`: original SciFort max-shift implementations of
  one-dimensional log-sum-exp, softmax, and log-softmax. No third-party
  implementation was copied or translated.


- `src/stats/nakagami.f90`, `src/stats/powernorm.f90`,
  `src/stats/loggamma.f90`, and `src/stats/wald.f90`: original SciFort
  implementations from the documented probability definitions and standard
  relationships to the gamma, normal, and inverse-Gaussian distributions.
  SciPy 1.17.0 is used only as an independent numerical validation reference;
  no SciPy implementation source was copied, translated, or adapted.
- `src/stats/gompertz.f90`, `src/stats/invweibull.f90`,
  `src/stats/betaprime.f90`, and `src/stats/burr12.f90`: original SciFort
  implementations from the published mathematical definitions of these
  distributions. Beta-prime reuses SciFort's independently implemented
  incomplete-beta kernel; the other three use project-local stable elementary
  helpers. SciPy 1.17.0 is used only as an independent numerical validation
  reference; no SciPy implementation source was copied, translated, or adapted.

- `src/stats/genhalflogistic.f90`, `src/stats/exponpow.f90`,
  `src/stats/exponweib.f90`, and `src/stats/powerlognorm.f90`: original
  SciFort implementations from the documented probability definitions and
  standard transform relationships. The implementations use project-local
  stable elementary and normal-distribution kernels. SciPy 1.17.0 is used
  only as an independent numerical validation reference; no SciPy
  implementation source was copied, translated, or adapted.


- `src/stats/levy_l.f90` and `src/stats/weibull_max.f90`: original SciFort
  implementations using reflection identities over SciFort's independently
  implemented Levy and Weibull-min kernels. SciPy 1.17.0 is used only as an
  independent numerical validation reference; no SciPy implementation source
  was copied, translated, or adapted.
- `src/stats/rdist.f90`: original SciFort implementation using the standard
  affine relationship between the R-distribution and a symmetric beta law;
  it reuses SciFort's independently implemented incomplete-beta and inverse
  incomplete-beta kernels. No third-party implementation source was copied or
  translated.
- `src/stats/skewcauchy.f90`: original SciFort implementation from the public
  skew-Cauchy probability definition, using independently selected stable
  arctangent and logarithmic rearrangements. SciPy 1.17.0 is used only for
  numerical validation; no SciPy implementation source was copied, translated,
  or adapted.
- `test/test_reflected_shape_four.f90`: contains only independently generated
  SciPy 1.17.0 numerical outputs for left-Levy, Weibull maximum,
  R-distribution, and skew-Cauchy regression cases; it contains no SciPy
  implementation source.

- `src/stats/dgamma.f90`: original SciFort implementation of the symmetric
  reflected-gamma definition, reusing SciFort's independently implemented
  incomplete-gamma kernels. No SciPy implementation source was copied or adapted.
- `src/stats/laplace_asymmetric.f90`: original SciFort implementation from the
  public piecewise asymmetric-Laplace density and its integrated probabilities.
  No SciPy implementation source was copied or adapted.
- `src/stats/truncnorm.f90`: original SciFort finite-bound implementation of a
  normal law truncated between standardized limits `a < b`, using SciFort's
  normal log-tail kernels and an independently written bisection inverse.
- `src/stats/loguniform.f90`: original SciFort implementation of the reciprocal
  density `1/(x*log(b/a))` on standardized support `[a,b]`.
- `test/test_scipy_transform_four.f90` and
  `test/scipy_transform_reference.f90`: contain only independently generated
  SciPy 1.17.0 numerical outputs for regression validation; no SciPy
  implementation source is included. The reference module is regenerated by
  `tools/generate_scipy_transform_reference.py`.

- `src/stats/foldnorm.f90` and `src/stats/foldcauchy.f90`: original SciFort
  implementations of folded normal and folded Cauchy laws from their public
  reflection identities. Stable mixture/arctangent rearrangements and
  bisection inverses were independently selected; no SciPy source was copied.
- `src/stats/recipinvgauss.f90`: original SciFort implementation using the
  reciprocal relationship to SciFort's independently implemented
  inverse-Gaussian law; no third-party implementation source was copied.
- `src/stats/truncpareto.f90`: original SciFort implementation from the
  normalized truncated-Pareto density on standardized support `[1,c]`, with
  independently selected log-domain normalizer and inverse formulas.
- `test/test_folded_reciprocal_four.f90` and
  `test/folded_reciprocal_reference.f90`: contain only independently generated
  SciPy 1.17.0 numerical outputs for regression validation. The reference module
  is regenerated by `tools/generate_folded_reciprocal_reference.py`; no SciPy
  implementation source is included.

- `src/stats/exponnorm.f90`: original SciFort implementation of the
  exponentially modified normal law from the convolution/normal-tail
  identities, with independently derived direct log-tail rearrangements and
  safeguarded bisection inverses. SciPy 1.17.0 was used only for public
  parameterization and numerical validation; no implementation source was
  copied or adapted.
- `src/stats/johnsonsb.f90` and `src/stats/johnsonsu.f90`: original SciFort
  implementations from the standard Johnson normalizing transformations
  `a+b*logit(z)` and `a+b*asinh(z)`, reusing SciFort's normal and logistic
  kernels. No third-party implementation source was copied.
- `src/stats/trapezoid.f90`: original SciFort piecewise implementation of the
  trapezoidal density, cumulative law, and exact inverse on standardized
  `[0,1]`. SciPy documentation was used only to confirm public shape
  conventions and numerical parity.
- `test/test_normal_transform_four.f90` and
  `test/normal_transform_reference.f90`: contain only independently generated
  SciPy 1.17.0 numerical outputs for regression validation. The reference
  module is regenerated by `tools/generate_normal_transform_reference.py`; no
  SciPy implementation source is included.

- `src/stats/burr.f90`: original SciFort implementation of the Burr Type III
  law from its public density and cumulative definitions, with independently
  derived log-domain tails and inverse formulas. SciPy 1.17.0 is used only to
  confirm public parameterization and numerical parity; no SciPy implementation
  source was copied, translated, or adapted.
- `src/stats/mielke.f90`: original SciFort implementation of the Mielke
  beta-kappa/Dagum law from its public density, with independently derived stable
  probability and quantile formulas. Its documented equivalence to Burr Type III
  is used as an independent regression identity. No third-party source was copied.
- `src/stats/gibrat.f90`: original SciFort implementation using the mathematical
  identity between the Gibrat law and SciFort's independently implemented
  lognormal distribution with shape one. No third-party implementation source was
  copied or adapted.
- `src/stats/wrapcauchy.f90`: original SciFort implementation of the wrapped
  Cauchy law on `[0,2*pi]` from the public density and circular CDF identities,
  with independently selected `atan2` and stable denominator rearrangements.
  SciPy 1.17.0 is used only for numerical validation.
- `test/test_burr_mielke_gibrat_wrapcauchy.f90` and
  `test/burr_mielke_gibrat_wrapcauchy_reference.f90`: contain only independently
  generated SciPy 1.17.0 numerical outputs and SciFort-owned regression logic. The
  reference module is regenerated by
  `tools/generate_burr_mielke_gibrat_wrapcauchy_reference.py`; no SciPy
  implementation source is included.

- `src/stats/genextreme.f90`: original SciFort implementation of the
  generalized extreme-value law from its public CDF and density equations,
  following SciPy's documented sign convention for the shape parameter. Stable
  logarithmic tails and the quantile transform were derived independently.
- `src/stats/kappa3.f90`: original SciFort implementation of the kappa-3 law
  from its public mathematical definition, with independently derived log-domain
  probability and inverse formulas. SciPy 1.17.0 is used only for parameter and
  numerical parity checks.
- `src/stats/kappa4.f90`: original SciFort implementation of the kappa-4 nested
  transform, including its shape-dependent support and exponential limiting
  cases. The implementation was derived from the public mathematical definition;
  no third-party implementation source was copied or adapted.
- `src/stats/truncweibull_min.f90`: original SciFort implementation of a Weibull
  minimum distribution conditioned to finite standardized bounds `a` and `b`,
  with direct log-domain normalization, tails, and inverse formulas.
- `test/test_genextreme_kappa_truncweibull.f90` and
  `test/genextreme_kappa_truncweibull_reference.f90`: contain only independently
  generated SciPy 1.17.0 numerical outputs and SciFort-owned regression logic.
  The reference module is regenerated by
  `tools/generate_genextreme_kappa_truncweibull_reference.py`; no SciPy
  implementation source is included.

- `src/stats/gengamma.f90`: original SciFort implementation of the generalized
  gamma law from its public density and regularized incomplete-gamma
  representation. SciPy 1.17.0 source and documentation were inspected to
  confirm parameterization, sign conventions, and tail selection; no SciPy
  implementation source was copied, translated, or adapted.
- `src/stats/halfgennorm.f90`: original SciFort implementation of the positive
  half generalized-normal law from its public density and incomplete-gamma
  CDF. Quantile and upper-tail formulas were derived from the same mathematical
  representation. No third-party implementation source was copied.
- `src/stats/argus.f90`: original SciFort implementation of the ARGUS law from
  its public density and normalization identity. The incomplete-gamma
  normalization and finite-support inverse were implemented independently;
  SciPy 1.17.0 was used to confirm numerical parity and public conventions.
- `src/stats/erlang.f90`: original SciFort specialization of the existing
  SciFort gamma kernels. It mirrors SciPy's numerical acceptance of positive
  noninteger shape values but does not reproduce the Python warning layer.
- `test/test_gengamma_halfgennorm_argus_erlang.f90` and
  `test/gengamma_halfgennorm_argus_erlang_reference.f90`: contain SciFort-owned
  regression logic and independently generated SciPy 1.17.0 numerical outputs.
  The reference module is regenerated by
  `tools/generate_gengamma_halfgennorm_argus_erlang_reference.py`; no SciPy
  implementation source is included.

- `src/stats/crystalball.f90`: original SciFort implementation of the Crystal
  Ball law from its public piecewise density. Normalization, direct logarithmic
  tails, and closed-form inverse branches were independently derived. SciPy
  1.17.0 was used to confirm parameterization and numerical behavior.
- `src/stats/jf_skew_t.f90`: original SciFort implementation of the Jones-Faddy
  skew-t law using its public beta-transform representation. Incomplete-beta
  tails and inverse-beta quantiles reuse SciFort kernels; no third-party
  implementation source was copied or adapted.
- `src/stats/pearson3.f90`: original SciFort implementation of Pearson III by
  mapping the public skewness parameterization to SciFort gamma and normal
  kernels. SciPy 1.17.0 was inspected to confirm its small-skew normal
  transition and public conventions; the implementation is independently
  written.
- `src/stats/rel_breitwigner.f90`: the density, stable survival form, inverse
  solver, likelihood, score, fit, and RNG integration are SciFort-owned. The
  complex CDF antiderivative is adapted from `rel_breitwigner_gen._cdf` in
  SciPy 1.17.0, `scipy/stats/_continuous_distns.py`,
  https://github.com/scipy/scipy/blob/v1.17.0/scipy/stats/_continuous_distns.py.
  Upstream copyright is Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy
  Developers, under BSD-3-Clause. The notice is retained in the source file and
  the complete applicable license text is in `THIRD_PARTY_LICENSES.md`. The
  adaptation replaces NumPy complex operations with Fortran intrinsic complex
  arithmetic, adds explicit loc/scale and endpoint handling, and derives a
  reciprocal-argument survival form for upper-tail stability.
- `test/test_crystalball_jf_pearson_breitwigner.f90` and
  `test/crystalball_jf_pearson_breitwigner_reference.f90`: contain SciFort-owned
  regression logic and generated SciPy 1.17.0 numerical outputs. The reference
  module is regenerated by
  `tools/generate_crystalball_jf_pearson_breitwigner_reference.py`; no SciPy
  implementation source is included in the generated reference data.


- `src/stats/genexpon.f90`: original SciFort implementation of the generalized
  exponential law from its public density and survival representation. The
  principal-real Lambert-W inverse is independently implemented from the
  defining equation with safeguarded Halley/Newton-style refinement; no
  third-party Lambert-W source was copied or adapted. SciPy 1.17.0 was used to
  confirm public parameterization and numerical outputs.
- `src/stats/skewnorm.f90`: original SciFort skew-normal implementation from
  `2*phi(z)*Phi(a*z)` and Owen's T integral. The fixed Gauss-Legendre Owen-T
  quadrature and transformed extreme-lower-tail quadrature are independently
  written from the integral definitions; no SciPy implementation source was
  copied or adapted.
- `src/stats/tukeylambda.f90`: original SciFort implementation from the public
  Tukey-lambda quantile definition. CDF inversion, the small-shape expansion,
  support handling, density derivative, and direct symmetric ISF were derived
  independently.
- `src/stats/rice.f90`: original SciFort implementation of the Rice law. The
  scaled modified-Bessel series/asymptotics and the Poisson/incomplete-gamma
  representation of the noncentral-chi-square tails are independently written
  from standard mathematical identities. No third-party Bessel or noncentral
  chi-square implementation source was copied or adapted.
- `test/test_genexpon_skewnorm_tukeylambda_rice.f90` and
  `test/genexpon_skewnorm_tukeylambda_rice_reference.f90`: contain SciFort-owned
  regression logic and generated SciPy 1.17.0 numerical outputs. The reference
  module is regenerated by
  `tools/generate_genexpon_skewnorm_tukeylambda_rice_reference.py`; no SciPy
  implementation source is included in the generated data.

- `src/stats/dpareto_lognorm.f90`: original SciFort implementation of the
  double-Pareto lognormal law from its published density/CDF identities. Tail
  expressions are evaluated in log space and quantiles are inverted on
  `log(x)` to preserve the heavy tails. SciPy 1.17.0 was used only to verify
  public parameterization and numerical outputs; no SciPy source was copied.
- `src/stats/vonmises.f90` and `src/stats/vonmises_line.f90`: original SciFort
  implementations from the public von Mises density. The line CDF is computed
  by independently written adaptive quadrature and the circular CDF extends it
  periodically. Quantiles use safeguarded bisection. No third-party CDF or
  quantile source was copied or adapted.
- `src/stats/kstwobign.f90`: original SciFort implementation of the asymptotic
  two-sided Kolmogorov law using complementary Jacobi-theta and alternating
  Kolmogorov series derived from the standard mathematical identities. The
  density and score derivative are differentiated forms of those series.
- `src/special/elementary.f90` (`i0e`, `i1e`): the independently written scaled
  modified-Bessel series/asymptotic kernels previously private to Rice were
  generalized to real arguments and promoted to the shared special-function
  layer. No external Bessel implementation was copied or adapted.
- `test/test_dpareto_vonmises_kstwobign.f90` and
  `test/dpareto_vonmises_kstwobign_reference.f90`: contain SciFort-owned
  regression logic and generated SciPy 1.17.0 numerical outputs. The reference
  module is regenerated by
  `tools/generate_dpareto_vonmises_kstwobign_reference.py`; no SciPy
  implementation source is included in the generated data.

- `src/stats/irwinhall.f90`: original SciFort implementation of the Irwin-Hall
  law. The density uses the standard cardinal-B-spline recurrence and the CDF
  the corresponding accumulated B-spline identity, chosen independently to
  avoid cancellation in the alternating polynomial representation. SciPy
  1.17.0 was used only to confirm public parameterization and numerical output.
- `src/stats/ksone.f90`: original SciFort implementation of the one-sided
  finite-sample Kolmogorov-Smirnov law from the classical
  Birnbaum-Tingey/Smirnov finite-sample sum. The logarithmically scaled survival
  and differentiated density/score sums were independently written; no SciPy
  implementation source was copied or adapted.
- `src/stats/ncx2.f90`: original SciFort noncentral chi-square implementation
  from its Poisson mixture of central chi-square laws. Mixture centering,
  direct incomplete-gamma tails, analytic score accumulation, and safeguarded
  quantile inversion were independently implemented.
- `src/stats/ncf.f90`: original SciFort noncentral F implementation from the
  Poisson mixture of beta-transformed central F components. Direct
  incomplete-beta tails, analytic mixture scores, and quantile inversion were
  independently implemented.
- `test/test_irwinhall_ksone_ncx2_ncf.f90` and
  `test/irwinhall_ksone_ncx2_ncf_reference.f90`: contain SciFort-owned
  regression logic and generated SciPy 1.17.0 numerical outputs. The reference
  module is regenerated by `tools/generate_irwinhall_ksone_ncx2_ncf_reference.py`;
  no SciPy implementation source is included in the generated data.

- `src/stats/randint.f90`: original SciFort implementation of the discrete
  uniform law on `low, ..., high-1`, including direct discrete tails and exact
  integer quantile correction. No SciPy implementation source was copied or
  translated.
- `src/stats/planck.f90`: original SciFort implementation of the Planck
  discrete exponential law from its geometric mass formula, using stable
  `expm1`/`log1p` identities for tails and normalization.
- `src/stats/dlaplace.f90`: original SciFort implementation of the symmetric
  discrete Laplace law from its closed-form mass and half-line tail sums.
- `src/stats/logser.f90`: original SciFort implementation of the logarithmic
  series law from its defining mass function. Lower and upper tail summations,
  recurrence termination, and integer quantile search were independently
  implemented.
- `test/test_randint_planck_dlaplace_logser.f90` and
  `test/randint_planck_dlaplace_logser_reference.f90`: contain SciFort-owned
  regression logic and generated SciPy 1.17.0 numerical outputs. The reference
  module is regenerated by
  `tools/generate_randint_planck_dlaplace_logser_reference.py`; no SciPy
  implementation source is included in the generated data.

- `src/stats/betabinom.f90`: original SciFort implementation of the
  beta-binomial law from the beta-binomial mass formula. Finite lower and upper
  tails are summed independently in log space and integer quantiles use binary
  search; no SciPy implementation source was copied or translated.
- `src/stats/hypergeom.f90`: original SciFort implementation of the
  hypergeometric law from products of binomial coefficients, evaluated through
  log-gamma identities with independently accumulated finite tails.
- `src/stats/nhypergeom.f90`: original SciFort implementation of the negative
  hypergeometric law from its combinatorial mass formula, with independent
  finite-tail log summation and integer quantile search.
- `src/stats/boltzmann.f90`: original SciFort implementation of the finite
  truncated-geometric/Boltzmann law. Stable `expm1`/`log1p` identities are used
  for normalization and both tails, including the small-rate regime.
- `test/test_betabinom_hypergeom_nhypergeom_boltzmann.f90` and
  `test/betabinom_hypergeom_nhypergeom_boltzmann_reference.f90`: contain
  SciFort-owned regression logic and generated public SciPy 1.17.0 numerical
  outputs. The generator additionally computes two stress tails directly from
  exact/high-precision defining formulas where SciPy 1.17.0 is less accurate;
  no SciPy implementation source is included.


- `src/special/zeta.f90`: original SciFort implementation of the real Hurwitz
  zeta function for `s > 1, q > 0` and its exponent derivative, using the
  Euler-Maclaurin representation in DLMF 25.11.43. No SciPy source was copied
  or adapted.
- `src/stats/betanbinom.f90`: original SciFort implementation of the
  beta-negative-binomial mass from binomial/beta identities, with a
  probability-ratio lower-tail recurrence and discrete quantile search.
- `src/stats/yulesimon.f90`: original SciFort implementation from
  `P(X=k)=alpha B(k,alpha+1)` and the exact beta-function survival identity.
- `src/stats/zipf.f90`: original SciFort implementation from the zeta-law mass
  and Hurwitz-zeta tail identity.
- `src/stats/zipfian.f90`: original SciFort implementation of the finite
  generalized-harmonic power law with independently accumulated finite tails.
- `test/test_betanbinom_yulesimon_zipf_zipfian.f90` and
  `test/betanbinom_yulesimon_zipf_zipfian_reference.f90`: SciFort-owned test
  logic and generated public SciPy 1.17.0 numerical outputs. The reference
  module is regenerated by
  `tools/generate_betanbinom_yulesimon_zipf_zipfian_reference.py`; no SciPy
  implementation source is included.

- `src/special/bessel_k.f90`: original SciFort implementation of real
  modified-Bessel `K_nu(x)` in logarithmic form from its standard integral
  representation, with composite Gauss-Legendre quadrature after logarithmic
  scaling and a standard large-argument asymptotic expansion. The mathematical
  specifications were taken from NIST DLMF sections 10.32 (integral
  representations) and 10.40 (asymptotic expansions). Adaptive panel sizing,
  truncation rules, and the order/argument log-derivative implementations are
  original SciFort code; no SciPy or other numerical-library source was copied
  or adapted.
- `src/stats/geninvgauss.f90`: original SciFort implementation of the
  generalized inverse-Gaussian density from its defining Bessel-K formula.
  Direct lower/upper tails use a logarithmic-coordinate integral with localized
  monotone-tail truncation and log-sum-exp quadrature. SciPy 1.17.0 public
  behavior was used for API and ordinary-value parity only.
- `src/stats/norminvgauss.f90`: original SciFort implementation of the normal
  inverse-Gaussian density from its standard Bessel-K formula. Direct tails use
  the substitution `z=sinh(y)`, localized log-domain quadrature, and direct
  log-tail inversion. SciPy 1.17.0 public behavior was used for API and
  ordinary-value parity only.
- `src/stats/skellam.f90`: original SciFort implementation of the difference of
  independent Poisson counts. The PMF uses the standard modified-Bessel-I
  identity and the tails use the standard noncentral-chi-square representation
  over SciFort's `ncx2` kernel. No SciPy source was copied or adapted.
- `test/test_gig_nig_skellam.f90` and
  `test/gig_nig_skellam_reference.f90`: SciFort-owned regression logic and
  generated public SciPy 1.17.0 numerical outputs. Additional deep-tail and
  wide-order Bessel-K constants were evaluated independently at 50-80 digit
  precision from the defining integral formulas. The SciPy reference module is
  regenerated by `tools/generate_gig_nig_skellam_reference.py`; no SciPy
  implementation source is included.

- `src/stats/genhyperbolic.f90`: original SciFort implementation of the
  generalized-hyperbolic density from its standard modified-Bessel-K formula.
  Direct CDF/SF evaluation uses the `z=sinh(y)` transformation and log-domain
  Gauss-Legendre tail integration. The exact `p<0, abs(b)=a` limiting
  normalizer is handled analytically. SciPy 1.17.0 public outputs were used only
  for API/numerical parity tests; no SciPy implementation source was copied or
  adapted.
- `src/stats/nchypergeom_fisher.f90`: original SciFort implementation of
  Fisher's noncentral hypergeometric law from its defining finite-support PMF,
  with log-sum-exp normalization and direct lower/upper sums. No code from
  SciPy's or BiasedUrn's implementations was copied, translated, or adapted.
- `test/test_genhyperbolic_fisher.f90` and
  `test/genhyperbolic_fisher_reference.f90`: SciFort-owned regression logic and
  generated public SciPy 1.17.0 numerical outputs. The reference module is
  regenerated by `tools/generate_genhyperbolic_fisher_reference.py`; no SciPy
  implementation source is included.

- `src/stats/nct.f90`: original SciFort implementation of the noncentral
  Student-t distribution from the representation `(Z+nc)/sqrt(2U/df)` with
  `U~Gamma(df/2,1)`. CDF, SF, and density are expectations in gamma-probability
  coordinates evaluated with a double-exponential endpoint map and adaptive
  spacing for small `df`. No SciPy/Boost implementation source was copied or
  adapted; SciPy 1.17.0 public outputs were used for ordinary parity only.
- `src/stats/gausshyper.f90`: original SciFort implementation of the Gauss
  hypergeometric density from its defining beta-weighted tilted density.
  Normalization and direct tails are integrated after a beta-quantile transform,
  avoiding direct endpoint-singular beta-density quadrature. No SciPy source was
  copied or adapted.
- `test/test_nct_gausshyper.f90` and `test/nct_gausshyper_reference.f90`:
  SciFort-owned regression logic and generated public SciPy 1.17.0 numerical
  outputs. Additional difficult-tail constants are independent defining-integral
  evaluations. The SciPy reference module is regenerated by
  `tools/generate_nct_gausshyper_reference.py`; no SciPy implementation source
  is included.

- `src/stats/landau.f90`: original SciFort implementation of the Landau
  distribution. The central density/CDF use the classical exponentially damped
  characteristic-function representation; the right tail uses a scaled form of
  Landau's nonoscillatory defining integral. The left-tail PDF/CDF asymptotic
  coefficients are the classical Koelbig-Schorr coefficients as tabulated in
  Takuma Yoshimura, "Numerical Evaluation and High Precision Approximation
  Formula for Landau Distribution" (2024), CC BY 4.0. No SciPy, Boost, or
  Yoshimura implementation source code was copied or adapted. SciPy 1.17.0
  public numerical outputs are used only for generated regression references.
- `test/test_landau.f90` and `test/landau_reference.f90`: SciFort-owned
  regression logic and generated public SciPy 1.17.0 numerical outputs. The
  reference module is regenerated by `tools/generate_landau_reference.py`.

- `src/stats/nchypergeom_wallenius.f90`: original SciFort implementation of
  Wallenius' noncentral hypergeometric distribution from the defining sequential
  biased-urn process. PMF/tails are obtained by exact finite-state Markov
  propagation, and the odds derivative is propagated through the same recursion.
  No SciPy or BiasedUrn source code was copied, translated, or adapted.
- `src/stats/poisson_binom.f90`: original SciFort implementation of the
  Poisson-binomial distribution from independent Bernoulli convolution. Direct
  probability/log-tail dynamic programs and leave-one-out convolutions provide
  the distribution and analytic vector score. No SciPy implementation source
  was copied or adapted.
- `test/test_wallenius_poisson_binom.f90` and
  `test/wallenius_poisson_binom_reference.f90`: SciFort-owned regression logic
  and generated public SciPy 1.17.0 numerical outputs. The reference module is
  regenerated by `tools/generate_wallenius_poisson_binom_reference.py`; only
  numerical reference values enter SciFort.

- `src/stats/kstwo.f90`: the finite-sample two-sided Kolmogorov-Smirnov
  regime selection, Durbin/Marsaglia-Tsang-Wang matrix CDF algorithm, and
  Pelz-Good large-sample expansion are translated/adapted from
  `scipy/stats/_ksstats.py` in SciPy 1.17.0, copyright Enthought, Inc. and the
  SciPy Developers, under BSD-3-Clause. SciFort omits SciPy's Pomeranz branch
  and instead extends the exact Durbin/MTW branch through the small-`n` region
  where both are accurate, while reusing SciFort's independent `ksone` kernel
  for the exact one-sided-tail regimes. The retained SciPy license is in
  `THIRD_PARTY_LICENSES.md`.
- `test/test_kstwo.f90` and `test/kstwo_reference.f90`: SciFort-owned
  regression logic and generated public SciPy 1.17.0 numerical outputs. The
  reference module is regenerated by `tools/generate_kstwo_reference.py`.

- `src/stats/levy_stable.f90`: the default S1 parameterization, S1-to-S0
  conversion, difficult-input rounding, analytic special cases, and
  Nolan/Zolotarev piecewise density/CDF formulas are translated/adapted from
  SciPy 1.17.0 `scipy/stats/_levy_stable/__init__.py` and
  `scipy/stats/_levy_stable/c_src/levyst.c`, copyright Enthought, Inc. and the
  SciPy Developers, under BSD-3-Clause. SciFort supplies its own fixed
  Gauss-Legendre and adaptive quadrature, tail inversion, likelihood, fitting,
  score, RNG, and ABI integration. The retained SciPy license is in
  `THIRD_PARTY_LICENSES.md`.
- `src/stats/studentized_range.f90`: original SciFort implementation from the
  defining normal-range/chi-square scale-mixture representation. The finite-df
  integrals correspond to the formulas described by Batista et al. (2017), and
  the infinite-df range integrals to Lund & Lund (1983). SciPy 1.17.0 was
  inspected for API semantics, the `df >= 100000` asymptotic switch, and
  numerical regression values; no SciPy implementation source is copied or
  adapted in this module.
- `test/test_levy_stable_studentized_range.f90` and
  `test/levy_stable_studentized_range_reference.f90`: SciFort-owned regression
  logic and generated public SciPy 1.17.0 numerical outputs. The reference
  module is regenerated by
  `tools/generate_levy_stable_studentized_range_reference.py`.

- `src/linalg.f90`: original SciFort compact dense-linear-algebra implementation.
  The user-supplied `Beliavsky/fortran-lapack` archive at commit
  `9b43beb7644479af46a56e743f61dbc8ae63ac7a` (BSD-3-Clause, copyright
  Federico Perini) was inspected for Cholesky interface/failure conventions and
  as a numerical-design reference. No fortran-lapack implementation source was
  copied, translated, or adapted; SciFort uses its own unblocked Cholesky,
  triangular solve, cyclic Jacobi eigensolver, and PSD wrapper.
- `src/stats/multivariate_integration.f90`: SciFort integration plumbing and
  randomized tent-transformed Halton rule around a conditional multivariate-box
  probability transform. The conditional normal/t transformation and the
  pivoted/permuted Cholesky logic are adapted from SciPy 1.17.0
  `scipy/stats/_qmvnt.py` and `scipy/stats/_qmvnt_cy.pyx` under SciPy's
  BSD-3-Clause license. SciFort does not copy SciPy's FFT/CBC lattice generator;
  the Halton sequence, randomized shifts, batching, stopping rules, singular-row
  handling, and explicit-state integration interface are SciFort-owned.
- `src/stats/multivariate_normal.f90`: original SciFort implementation of the
  multivariate-normal density, entropy, sampling, MLE, marginal-extraction, and
  public CDF/log-CDF wiring. The scalar/diagonal-covariance log-CDF path uses
  independence and logarithmic differences of the normal CDF/SF, derived from
  the error-function definition in DLMF 7.2.1 (https://dlmf.nist.gov/7.2.E1).
  This path and its identity-based regression tests are original SciFort code.
  SciPy 1.17.0 `scipy/stats/_multivariate.py` and
  `scipy/stats/_covariance.py` were inspected for dense-covariance API semantics,
  lower-triangle behavior, `_PSD` numerical rank cutoff, singular-support
  tolerance, fitting conventions, lower-limit semantics, and CDF defaults. The
  multidimensional box-probability kernel is the separately documented adapted
  code in `src/stats/multivariate_integration.f90`.
- `src/stats/multivariate_t.f90`: original SciFort implementation of the
  multivariate Student-t density/log-density, entropy, explicit-state sampling,
  marginal extraction, and public CDF wiring from the standard defining
  formulas. SciPy 1.17.0 was inspected for shape-matrix semantics, singular-matrix
  behavior, entropy asymptotics, default degrees of freedom, CDF point limits,
  and public method conventions. Its multidimensional CDF uses the separately
  documented adapted integration kernel.
- `src/stats/dirichlet.f90`, `src/stats/multinomial.f90`, and
  `src/stats/dirichlet_multinomial.f90`: original SciFort implementations from
  the standard defining probability formulas. SciPy 1.17.0 was inspected for
  parameter validation, the Dirichlet `K-1` coordinate convention, multinomial
  final-probability adjustment, entropy definition, and public method semantics.
- `test/test_multivariate.f90` and `test/multivariate_reference.f90`: SciFort-owned
  regression logic and generated public SciPy 1.17.0 numerical outputs. The
  reference module is regenerated by `tools/generate_multivariate_reference.py`.
- `src/stats/multivariate_hypergeom.f90`: original Fortran implementation of the
  standard multivariate-hypergeometric mass and moments. SciPy 1.17.0
  `scipy/stats/_multivariate.py` was inspected for parameter/support semantics and
  its sequential univariate-hypergeometric sampling decomposition; the translated
  sampling control flow is treated as a SciPy BSD-3-Clause adaptation. The
  retained SciPy license is in `THIRD_PARTY_LICENSES.md`.
- `src/stats/normal_inverse_gamma.f90`: original SciFort implementation from the
  standard normal/inverse-gamma joint density and moments. SciPy 1.17.0
  `scipy/stats/_multivariate.py` was inspected for defaults, existence conditions,
  support semantics, and the inverse-gamma/conditional-normal sampling definition.
- `src/stats/matrix_distribution_helpers.f90`, `src/stats/matrix_normal.f90`,
  `src/stats/wishart.f90`, and `src/stats/invwishart.f90`: original SciFort
  implementations of standard SPD, multivariate-gamma/digamma, matrix-normal,
  Wishart, and inverse-Wishart formulas. SciPy 1.17.0 `_multivariate.py` was
  inspected for lower-triangle semantics, parameter validation, moment/entropy
  conventions, and public API behavior. Wishart RNG uses the standard Bartlett
  construction; inverse-Wishart RNG uses a Wishart-inversion identity rather than
  SciPy's direct Axen implementation.
- `src/stats/matrix_t.f90`: SciFort density implementation from the standard
  two-sided matrix-t determinant formula. SciPy 1.17.0 `_multivariate.py` was
  inspected for the new public API, defaults, and the equivalent inverse-Wishart/
  matrix-normal mixture; the sampling decomposition is treated as a SciPy
  BSD-3-Clause adaptation. The retained SciPy license is in
  `THIRD_PARTY_LICENSES.md`.
- `test/test_matrix_distributions.f90` and
  `test/test_matrix_distributions_reference.f90`: SciFort-owned regression logic
  and generated public SciPy 1.17.0 numerical outputs. The reference module is
  regenerated by `tools/generate_matrix_distributions_reference.py`.
- `src/special/bessel_i.f90`: original SciFort implementation of the positive-term
  series for real-order modified Bessel I and its large-argument asymptotic
  expansion, based on DLMF 10.25.2 and 10.40.1. SciPy `special.ive` is used only
  for numerical regression values.
- `src/stats/uniform_direction.f90`: original normalized-Gaussian unit-sphere
  sampler. SciPy 1.17.0 `_multivariate.py` was inspected for public shape and
  parameter semantics; no SciPy implementation source is copied.
- `src/stats/random_matrix_helpers.f90`, `src/stats/ortho_group.f90`,
  `src/stats/special_ortho_group.f90`, and `src/stats/unitary_group.f90`: original
  SciFort Gaussian-QR/Haar implementations using re-orthogonalized modified
  Gram-Schmidt, based on the standard construction described by F. Mezzadri,
  Notices AMS 54 (2007), 592-604. SciPy 1.17.0 was used as an API/behavioral
  reference only.
- `src/stats/vonmises_fisher.f90`: original SciFort density, entropy, fit, and
  explicit-state sampling implementation. The dimension-specific sampling strategy
  follows J. Wenzel, "Numerically stable sampling of the von Mises Fisher
  distribution on S2", in dimension three and A. T. A. Wood, Communications in Statistics - Simulation and Computation 23 (1994),
  157-164, in dimensions four and above. SciPy 1.17.0 `_multivariate.py` was
  inspected for parameter validation and public behavior; no SciPy source code is
  copied.
- `src/stats/random_correlation.f90`: original Fortran implementation of the
  prescribed-eigenvalue correlation-matrix construction of P. I. Davies and
  N. J. Higham, BIT 40 (2000), 640-651. SciPy 1.17.0 `_multivariate.py` was
  inspected for tolerance defaults and validation semantics; no SciPy source code
  is copied.
- `src/stats/random_table.f90`: original fixed-margin contingency-table mass, mean,
  and exact row-conditional multivariate-hypergeometric sampler. SciPy 1.17.0 was
  inspected for public support semantics and numerical references; the sampler is
  intentionally independent of SciPy's compiled Boyett/Patefield implementations.
- `test/test_directional_random_matrix.f90` and
  `test/test_directional_random_matrix_reference.f90`: SciFort-owned regression
  logic and generated public SciPy 1.17.0 numerical outputs. The reference module
  is regenerated by `tools/generate_directional_random_matrix_reference.py`.

- `src/stats/gaussian_kde.f90`: original SciFort implementation of the
  standard Gaussian-mixture KDE formulas. SciPy 1.17.0 `scipy/stats/_kde.py` was
  inspected for public API behavior, data orientation, effective-sample-size and
  weighted-covariance conventions, bandwidth factors, integration semantics, and
  marginal behavior. No SciPy implementation source was copied or translated.
- `test/test_gaussian_kde.f90` and `test/gaussian_kde_reference.f90`: SciFort-owned
  regression logic and generated public SciPy 1.17.0 numerical outputs. The
  reference module is regenerated by `tools/generate_gaussian_kde_reference.py`.
- `src/stats/qmc.f90`: original Fortran implementation of standard QMC formulas and
  constructions. SciPy 1.17.0 `scipy/stats/_qmc.py`, `_qmc_cy.pyx`, and `_sobol.pyx`
  were inspected for public semantics, discrepancy formulas, Halton/Sobol indexing,
  LMS+digital-shift structure, strength-two OA-LHS constraints, QMC distribution
  transforms, and Poisson-disk behavior. The discrepancy equations, sequence engines,
  scrambling code, multinomial/normal transforms, Bridson sampler, Fisher-Yates
  permutations, Prim MST, and orthogonal-array construction were independently
  written in Fortran; no SciPy implementation source was copied or translated.
  Randomized paths intentionally use SciFort's explicit MRG32k3a state instead of
  NumPy RNGs.
- `src/stats/qmc_sobol_data.f90`: generated Fortran representation of the Joe-Kuo
  search-criterion-6 Sobol direction-number initialization table distributed in
  SciPy 1.17.0 as `_sobol_direction_numbers.npz`. This data file is carried under
  SciPy's BSD-3-Clause terms; the original direction-number source is the UNSW
  Joe-Kuo table cited in SciPy's `_sobol.pyx`.
- `test/test_qmc.f90` and `test/qmc_reference.f90`: SciFort-owned regression logic
  and generated public SciPy 1.17.0 numerical outputs, regenerated by
  `tools/generate_qmc_reference.py`.

- `src/goodness_of_fit.f90`: SciFort implementation of goodness-of-fit and
  normality tests. The Shapiro-Wilk W/p-value kernel is adapted from SciPy 1.17.0
  `scipy/stats/_ansari_swilk_statistics.pyx` (Royston AS R94 / AS 181 lineage);
  the finite-sample Cramer-von Mises correction and Anderson-Darling critical
  value tables/public fitting conventions were adapted from SciPy 1.17.0
  `scipy/stats/_hypotests.py` and `_morestats.py`. The KS lattice probability,
  moment tests, optimization code, and remaining Fortran implementation were
  independently written. SciPy's BSD-3-Clause notice is retained in the source
  header and `THIRD_PARTY_LICENSES.md`.
- `test/test_goodness_of_fit.f90` and `test/goodness_of_fit_reference.f90`:
  SciFort-owned regression logic and public SciPy 1.17.0 numerical outputs; the
  reference module is regenerated by `tools/generate_goodness_of_fit_reference.py`.

Test reference data in `test/gamma_reference.f90`,
`test/beta_reference.f90`, `test/discrete_reference.f90`,
`test/discrete_extended_reference.f90`, and `test/continuous_reference.f90`
is generated by the corresponding scripts under `tools/` with mpmath, which
is BSD licensed; only computed numerical values, not mpmath source, enter the
repository. `test/extended_continuous_reference.f90`, `test/closed_form_reference.f90`,
`test/shape_family_reference.f90`, `test/transform_family_reference.f90`,
`test/positive_shape_reference.f90`, `test/scipy_transform_reference.f90`, and
`test/folded_reciprocal_reference.f90`, `test/normal_transform_reference.f90`, and
`test/burr_mielke_gibrat_wrapcauchy_reference.f90`, and
`test/genextreme_kappa_truncweibull_reference.f90`, and
`test/gengamma_halfgennorm_argus_erlang_reference.f90`, and
`test/crystalball_jf_pearson_breitwigner_reference.f90`, and
`test/genexpon_skewnorm_tukeylambda_rice_reference.f90`, and
`test/dpareto_vonmises_kstwobign_reference.f90`, and
`test/irwinhall_ksone_ncx2_ncf_reference.f90`, and
`test/randint_planck_dlaplace_logser_reference.f90`, and
`test/betabinom_hypergeom_nhypergeom_boltzmann_reference.f90`, and
`test/betanbinom_yulesimon_zipf_zipfian_reference.f90`, and
`test/gig_nig_skellam_reference.f90`, and
`test/genhyperbolic_fisher_reference.f90`,
`test/nct_gausshyper_reference.f90`, `test/landau_reference.f90`,
`test/wallenius_poisson_binom_reference.f90`, `test/kstwo_reference.f90`,
`test/levy_stable_studentized_range_reference.f90`,
`test/multivariate_reference.f90`, `test/test_matrix_distributions_reference.f90`,
`test/test_directional_random_matrix_reference.f90`, `test/gaussian_kde_reference.f90`,
and `test/qmc_reference.f90` are generated from public SciPy 1.17.0
numerical outputs by their corresponding
generator scripts under `tools/`;
only numerical reference values enter SciFort, and no SciPy implementation
source is copied or adapted. `test/test_transform_shape_four.f90` likewise
contains only independently generated SciPy 1.17.0 numerical outputs for the
generalized half-logistic, exponential-power, exponentiated-Weibull, and
power-lognormal regression cases; it contains no SciPy implementation source.

The documented SciPy source adaptations are relativistic Breit-Wigner,
finite-sample `kstwo`, Levy-stable, multivariate integration, multivariate-
hypergeometric sampling, matrix-t sampling, the Kendall exact recurrence, and
the identified goodness-of-fit kernels. The Sobol initialization data are also
redistributed from SciPy. These entries must not be described as wholly original
MIT-only implementations. The remaining entries state their individual origins.

## Merge-regression references

`test/test_merge_regressions.f90` contains original SciFort regression checks
based on support endpoints, NaN propagation, translation invariance, equal-weight
softmax normalization, uniform-distribution quantiles, and independence of normal
coordinates. It reuses the existing scalar normal log tails as an identity check;
no third-party implementation or numerical tables were copied for these tests.
The verbatim `LICENSES/CC-BY-4.0.txt` was retrieved on 2026-09-23 from the
research repository and the canonical Creative Commons legal-code URL recorded
in `THIRD_PARTY_LICENSES.md`. Its SHA-256 is
`9ba9550ad48438d0836ddab3da480b3b69ffa0aac7b7878b5a0039e7ab429411`.

## Required entry format for future imports

Add a row or section containing:

- SciFort file and procedure names;
- upstream project, file path, and URL;
- exact version, tag, or commit hash;
- original authors and copyright holders;
- original license and location of retained license text;
- whether the code was copied, translated, adapted, or independently
  reimplemented from a publication;
- meaningful modernization changes;
- tests used to establish parity.
