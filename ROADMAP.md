# Roadmap

## Milestone 1: distribution foundation

- stabilize the current six distributions;
- add high-precision generated reference tables;
- test GNU Fortran, Intel `ifx`, and LLVM Flang;
- establish documented accuracy targets;
- complete one production-quality Python wrapper through the C ABI.

## Milestone 2: special functions

- incomplete gamma and inverse;
- incomplete beta and inverse;
- inverse error functions;
- stable gamma ratios and log-beta;
- reusable root-solving and continued-fraction infrastructure.

Every imported algorithm requires a license and provenance review.

## Milestone 3: core SciPy-style distributions

- [x] beta, gamma, chi-square, Student t, and F;
- [x] Bernoulli, binomial, geometric, negative binomial, and Poisson;
- [x] direct stable survival and logarithmic functions;
- [x] documented C ABI coverage beyond the initial normal-distribution slice.

## Milestone 4: random variates and fitting

- [x] explicit-state RNG interface;
- [x] distribution sampling;
- [x] array log-likelihood and negative-log-likelihood evaluation;
- [x] score functions;
- [x] bounded and fixed-parameter fitting interfaces.

## Milestone 5: descriptive statistics and tests

- [x] moments, quantiles, covariance, correlation, ranks;
- [x] t-tests, correlation tests, and Mann-Whitney U;
- [x] rank-sum/signed-rank, Kruskal-Wallis, and Friedman tests;
- [x] power-divergence/chi-square/Fisher contingency tests;
- [x] classical/Welch one-way ANOVA and Bartlett/Levene/Fligner variance tests;
- [x] one-/two-sample KS and Cramer-von Mises goodness-of-fit tests;
- [x] Shapiro-Wilk, D'Agostino/Pearson, Jarque-Bera, and Anderson-Darling normality/fit tests;
- [x] Ansari-Bradley, Mood, Epps-Singleton, k-sample Anderson-Darling, and median tests;
- [x] Kendall/point-biserial association, linear/robust slopes, Brunner-Munzel, and Page trend tests;
- [x] binomial, Barnard/Boschloo exact, odds-ratio/relative-risk, contingency association, p-value combination, and FDR control;
- [x] Alexander-Govern, Tukey/Tukey-Kramer/Games-Howell, Dunnett, and Poisson means E-tests;
- [x] resampling infrastructure: percentile/basic/BCa bootstrap, exact/randomized permutation tests, Monte Carlo tests, and power simulation.

Breadth should expand only after the numerical foundations and provenance
process are demonstrably reliable.

## Milestone 6: core SciPy-special compatibility

- [x] error functions and inverse error functions;
- [x] log-gamma, positive-shape log-beta, and real digamma/psi;
- [x] normal-CDF and logistic transforms;
- [x] entropy, relative-entropy, and stable elementary helpers;
- [x] Box-Cox transforms and inverses;
- [x] one-dimensional log-sum-exp, softmax, and log-softmax reductions.

The first reduction API is intentionally unweighted and one-dimensional.
Broadcasting, axis selection, weighted/sign-returning log-sum-exp, and complex
special functions remain future work.

## Milestone 7: additional continuous distributions

- [x] right- and left-skewed Gumbel distributions;
- [x] power-function distribution;
- [x] triangular distribution;
- [x] generalized Pareto distribution;
- [x] likelihood, analytic-score, bounded-fit, random-variate, and scalar C-ABI
  integration for every new family;
- [x] arcsine and half-normal distributions with the same eight-function and
  cross-layer interfaces.
- [x] half-Cauchy and Lomax distributions with the same likelihood, score,
  fitting, RNG, and scalar C-ABI integration.
- [x] chi and Maxwell distributions with the same eight-function and
  cross-layer interfaces.
- [x] cosine and semicircular bounded distributions with full likelihood,
  score, fitting, RNG, and scalar C-ABI integration.
- [x] Anglit and Moyal loc/scale distributions with the same eight-function
  distribution API and full likelihood, score, fitting, RNG, and scalar C-ABI
  integration.
- [x] Hyperbolic-secant and half-logistic distributions with the same
  eight-function API and full likelihood, score, fitting, RNG, and scalar C-ABI
  integration.
- [x] Inverse-gamma, inverse-Gaussian, Levy, and log-Laplace positive-support
  families with full likelihood, analytic-score, fitting, RNG, and scalar
  C-ABI integration.
- [x] Bradford, truncated exponential, Fisk/log-logistic, and double-Weibull
  closed-form families with full eight-function, likelihood, analytic-score,
  fitting, RNG, and scalar C-ABI integration.
- [x] Alpha, Birnbaum-Saunders/fatigue-life, generalized logistic, and
  generalized normal one-shape families with the same full-stack integration.
- [x] Nakagami, power-normal, log-gamma, and Wald transform/shape families with
  full eight-function, likelihood, analytic-score, bounded-fit, RNG, and scalar
  C-ABI integration.
- [x] Gompertz, inverse-Weibull, beta-prime, and Burr XII positive-support
  families with full eight-function, likelihood, analytic-score, bounded-fit,
  RNG, and scalar C-ABI integration.
- [x] Generalized half-logistic, exponential-power, exponentiated-Weibull, and
  power-lognormal transform/shape families with full eight-function,
  likelihood, analytic-score, bounded-fit, RNG, and scalar C-ABI integration.
- [x] Left-Levy, Weibull maximum, R-distribution, and skew-Cauchy reflected or
  transform families with full eight-function, likelihood, analytic-score,
  bounded-fit, RNG, and scalar C-ABI integration.
- [x] Double-gamma, asymmetric Laplace, finite-bound truncated normal, and
  log-uniform/reciprocal families with full eight-function, likelihood,
  analytic-score, bounded-fit, RNG, and scalar C-ABI integration.
- [x] Folded normal, folded Cauchy, reciprocal inverse-Gaussian, and truncated
  Pareto families with full eight-function, likelihood, analytic-score,
  bounded-fit, RNG, and scalar C-ABI integration.
- [x] Exponentially modified normal, Johnson SB, Johnson SU, and trapezoid
  families with full eight-function, likelihood, analytic-score, bounded-fit,
  RNG, and scalar C-ABI integration.
- [x] Burr Type III, Mielke beta-kappa/Dagum, Gibrat, and wrapped Cauchy
  families with full eight-function, likelihood, analytic-score, bounded-fit,
  RNG, and scalar C-ABI integration.
- [x] Generalized extreme-value, kappa-3, kappa-4, and doubly truncated
  Weibull-minimum families with full eight-function, likelihood, analytic-score,
  bounded-fit, RNG, and scalar C-ABI integration.
- [x] Generalized-gamma, half-generalized-normal, ARGUS, and Erlang families with
  full eight-function, likelihood, analytic-score, bounded-fit, RNG, and scalar
  C-ABI integration.
- [x] Crystal Ball, Jones-Faddy skew-t, Pearson III, and relativistic
  Breit-Wigner families with full eight-function, likelihood, analytic-score,
  bounded-fit, RNG, and scalar C-ABI integration.
- [x] Generalized exponential, skew-normal, Tukey lambda, and Rice families
  with full eight-function, likelihood, analytic-score, bounded-fit, RNG, and
  scalar C-ABI integration.
- [x] Double-Pareto lognormal, circular von Mises, von Mises-on-a-line, and
  asymptotic two-sided Kolmogorov families with full eight-function,
  likelihood, analytic-score, bounded-fit, RNG, and scalar C-ABI integration.
- [x] Irwin-Hall, one-sided finite-sample Kolmogorov-Smirnov (`ksone`),
  noncentral chi-square, and noncentral F families with full eight-function,
  likelihood, analytic-score, bounded-fit, RNG, and scalar C-ABI integration.
- [x] Discrete uniform (`randint`), Planck, discrete Laplace, and logarithmic-series
  families with full eight-function, likelihood, analytic-shape-score where
  differentiable, bounded-fit, RNG, and scalar C-ABI integration.
- [x] Beta-binomial, hypergeometric, negative-hypergeometric, and Boltzmann
  discrete families with full eight-function, likelihood, analytic-shape-score
  where differentiable, bounded-fit, RNG, and scalar C-ABI integration.
- [x] Beta-negative-binomial, Yule-Simon, Zipf/zeta, and finite Zipfian
  discrete families with full eight-function, likelihood, analytic-shape-score,
  bounded-fit, RNG, and scalar C-ABI integration, backed by an independent
  Euler-Maclaurin real Hurwitz-zeta kernel.
- [x] Generalized inverse-Gaussian, normal-inverse-Gaussian, and Skellam
  families with full eight-function, likelihood, analytic-score, bounded-fit,
  RNG, and scalar C-ABI integration, backed by an independent real
  modified-Bessel-K log/derivative kernel.
- [x] Generalized hyperbolic and Fisher noncentral-hypergeometric families with
  full eight-function, likelihood, analytic-score, bounded-fit, RNG, and scalar
  C-ABI integration. Generalized hyperbolic reuses a scaled-log Bessel-K kernel;
  Fisher uses exact finite-support log-sum-exp normalization.
- [x] Noncentral Student t and Gauss hypergeometric families with full
  eight-function, likelihood, analytic-score, bounded-fit, RNG, and scalar C-ABI
  integration. Noncentral t uses gamma-probability quadrature with an adaptive
  low-df grid; Gauss hypergeometric integrates over a beta-quantile transform so
  endpoint-singular beta densities are never integrated directly.
- [x] Landau distribution with full eight-function, likelihood, analytic-score,
  bounded-fit, RNG, and scalar C-ABI integration. The implementation combines
  left-tail asymptotics, central characteristic-function quadrature, and a
  scaled nonoscillatory right-tail integral.
- [x] Wallenius noncentral-hypergeometric and Poisson-binomial discrete families
  with full eight-function, likelihood, analytic-score, bounded-fit, RNG, and
  C-ABI integration. Wallenius uses an exact finite-state sequential biased-urn
  dynamic program, including propagated odds derivatives and exact sequential
  sampling. Poisson-binomial uses log-domain Bernoulli convolution, a vector
  leave-one-out score, exact Bernoulli-sum sampling, and an array-aware C ABI.
- [x] Finite-sample two-sided Kolmogorov-Smirnov (`kstwo`) with full
  eight-function, likelihood/score, bounded-fit, RNG, and scalar C-ABI
  integration. The CDF/SF implementation adapts SciPy 1.17.0's BSD-licensed
  Simard-L'Ecuyer regime selection, Durbin/MTW matrix method, and Pelz-Good
  expansion, while reusing SciFort's independent exact one-sided KS kernel.
- [x] Levy-stable (`levy_stable`) and studentized-range (`studentized_range`)
  continuous families with full eight-function, likelihood, numerical-score,
  bounded-fit, deterministic inverse-transform RNG, `scifort_stats`, and scalar
  C-ABI integration. Levy-stable uses the SciPy 1.17.0 Nolan/Zolotarev S1/S0
  piecewise formulation with SciFort-owned quadrature. Studentized range uses
  the defining normal-range/chi-square scale-mixture integrals and the
  infinite-degrees-of-freedom range formula for `df >= 100000`.

These additions bring the scalar univariate distribution inventory to 130 families.
The multivariate milestone below now adds eighteen further named distribution families.
Further breadth should continue in coherent groups and must include the same
cross-layer validation rather than adding density functions in isolation.

## Milestone 8: multivariate distributions

- [x] Add a dependency-free small dense linear-algebra layer with lower Cholesky,
  lower-triangular solve, symmetric eigendecomposition, and PSD decomposition.
- [x] Add multivariate normal PDF/log-PDF, entropy, positive-semidefinite and
  singular-support handling, explicit-state RNG, and maximum-likelihood fitting.
- [x] Add Dirichlet PDF/log-PDF, moments, entropy, and explicit-state RNG.
- [x] Add multinomial PMF/log-PMF, moments, entropy, and explicit-state RNG.
- [x] Add Dirichlet-multinomial PMF/log-PMF and moments.
- [x] Add multivariate-normal CDF/log-CDF and marginal extraction.
- [x] Add multivariate Student t.
- [x] Add multivariate hypergeometric and normal-inverse-gamma families.
- [x] Add matrix normal, Wishart, inverse-Wishart, and matrix-t families.
- [x] Add directional and random-matrix distributions: von Mises-Fisher, uniform
  direction, random correlation, orthogonal/special-orthogonal/unitary groups, and
  fixed-margin random tables.



## Milestone 9: density estimation and quasi-Monte Carlo

- [x] Add weighted/unweighted Gaussian KDE with Scott, Silverman, and constant
  bandwidth factors; PDF/log-PDF evaluation; Gaussian, one-dimensional and
  multidimensional box, and KDE-product integrals; marginal extraction; inverse
  covariance; and explicit-state resampling.
- [x] Add QMC unit-cube scaling/reverse scaling and integer mapping.
- [x] Add centered, wrap-around, mixture, and L2-star discrepancies, iterative
  centered-discrepancy updates, and geometric minimum-distance/MST criteria.
- [x] Add Van der Corput and Halton sequences, including SciFort explicit-state
  Owen-style digit scrambling, reset, and fast-forward behavior.
- [x] Add ordinary Latin-hypercube sampling and prime-square strength-two
  orthogonal-array Latin hypercubes.
- [x] Add Sobol sequences with bundled Joe-Kuo direction numbers and LMS+shift scrambling.
- [x] Add Poisson-disk sampling with volume/surface proposals, arbitrary finite bounds, reset, and fill-space behavior.
- [ ] Add optional post-sampling `random-cd` and Lloyd optimization; Lloyd requires a future Voronoi/Qhull-capable geometry layer.
- [x] Add `MultinomialQMC` and `MultivariateNormalQMC` using Sobol as the default engine, including singular PSD covariance support and inverse-transform/Box-Muller normal generation.
