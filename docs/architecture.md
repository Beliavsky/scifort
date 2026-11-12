# Architecture

## Layers

SciFort separates four concerns:

1. **Numerical core**: pure modern Fortran algorithms.
2. **Native API**: idiomatic module procedures for Fortran users.
3. **C ABI**: a small interoperable and versioned procedural interface.
4. **Language wrappers**: Python, R, MATLAB, and Octave conventions built over
   the common C ABI.

The numerical core must not depend on Python, R, MATLAB, Octave, Cython, MEX,
or compiler-specific Fortran descriptors.

## Source layout

```text
src/
    kinds.f90
    constants.f90
    math.f90
    descriptive.f90
    descriptive_lowrank.f90
    descriptive_highrank_means.f90
    descriptive_weighted_engine.f90
    descriptive_highrank_weights_4_9.f90
    descriptive_highrank_weights_10_15.f90
    descriptive_highrank_stats_engine.f90
    descriptive_highrank_stats_5_9.f90
    descriptive_highrank_stats_10_15.f90
    descriptive_named_scale_engine.f90
    descriptive_named_scale_lowrank.f90
    descriptive_named_scale_highrank.f90
    descriptive_standardize_engine.f90
    descriptive_standardize_1_9.f90
    descriptive_standardize_10_15.f90
    hypothesis.f90
    hypothesis_extended.f90
    contingency_meta.f90
    multiple_comparisons.f90
    nonparametric_extended.f90
    goodness_of_fit.f90
    random.f90
    resampling.f90
    special.f90
    special/
        <function family>.f90
    stats.f90
    stats/
        <distribution>.f90
    c_api/
        api.f90
```

## Native API

Distribution procedure names use a family prefix, for example:

```text
normal_pdf
normal_logpdf
normal_cdf
normal_sf
normal_logcdf
normal_logsf
normal_ppf
normal_isf
```

Scalar arguments are elemental. Array reductions, fitting routines, random
sampling, and algorithms requiring workspaces should use separate procedural
interfaces rather than overloading scalar distribution evaluation. The same
rule applies to special functions: scalar transforms in `scifort_special` are
elemental, while `logsumexp`, `softmax`, and `log_softmax` are explicit reduction
interfaces. They support real or complex vectors through rank-15 arrays; array calls
may reduce/normalize one Fortran dimension, a tuple-style list of dimensions,
or, where meaningful, the whole array. Weighted log-sum-exp accepts the
fixed-rank broadcast forms documented in the native API through rank 15, and
explicit keepdims companions preserve rank when dimensions are reduced.
For descriptive `gmean`, `hmean`, and `pmean`, ranks 5 through 15 share an
assumed-rank dispatcher that flattens the source and indexes reduction slices.
Weighted rank-4 through rank-15 calls use a second assumed-rank engine that
computes trailing joint broadcast indices without materializing full broadcast
arrays. `sem`, `variation`, `iqr`, `median_abs_deviation`, `moment`, `skew`, and
`kurtosis` ranks 5 through 15 use a third shared assumed-rank slice engine with
the same thin-wrapper pattern. Thin rank-specific functions exist only because
Fortran function result rank is static. The numerical formulas remain in the
vector kernels.
`zscore`, `gzscore`, and `gstd` use a separate assumed-rank standardization
engine through rank 15. It gathers one-dimensional slices in column-major order,
computes scaled mean/standard-deviation quantities without materializing rank-
specific work arrays, and scatters standardized values back to their original
positions. Rank-specific wrappers only fix the static result rank; `gstd` uses
fixed-rank reduction wrappers plus `gstd_keepdims` companions.

## C ABI

The C ABI uses:

- `double` for initial public real values;
- `size_t` for array lengths;
- caller-owned arrays;
- `int` status outputs;
- stable lower-case symbol names beginning with `scifort_`.

Future opaque state, workspace, fitted-model, or RNG objects must be represented
by `c_ptr` handles with explicit create and destroy functions.

## Wrapper policy

Language wrappers own language-specific behavior:

- Python: NumPy conversion, broadcasting, exceptions, wheel packaging.
- R: vector recycling decisions, `NA` policy, registration, package metadata.
- MATLAB: MEX arrays, column-major dimension handling, error translation.
- Octave: oct-file or compatible MEX interface and package metadata.

Wrappers must not independently reimplement numerical formulas.
