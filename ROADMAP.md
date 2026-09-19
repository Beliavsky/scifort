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

- beta, gamma, chi-square, Student t, and F;
- Bernoulli, binomial, geometric, negative binomial, and Poisson;
- direct stable survival and logarithmic functions;
- documented C ABI coverage.

## Milestone 4: random variates and fitting

- explicit-state RNG interface;
- distribution sampling;
- likelihood evaluation and score functions;
- bounded and fixed-parameter fitting interfaces.

## Milestone 5: descriptive statistics and tests

- moments, quantiles, covariance, correlation, ranks;
- classical parametric and nonparametric tests;
- resampling infrastructure.

Breadth should expand only after the numerical foundations and provenance
process are demonstrably reliable.
