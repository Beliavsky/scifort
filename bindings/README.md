# Language bindings

Language bindings are intentionally not implemented as separate numerical
libraries. They will wrap the stable C ABI declared in `include/scifort.h`.

Binding directories (Python implemented; others planned):

```text
bindings/python/
bindings/R/
bindings/matlab/
bindings/octave/
```

The current development version provides scalar C ABI entry points for PDF
or PMF, CDF, and PPF across every implemented distribution. Normal PDF, CDF,
and PPF additionally provide bulk array entry points with explicit lengths and
status codes. Generated normal-distribution scalar and bulk entry points now
also cover SF, log-PDF, log-CDF, log-SF, and ISF; see the
[generator configuration and coverage](cabi_config/README.md). Existing PDF,
CDF, and PPF symbols retain their original implementations and behavior.
Outside the eleven generated families listed below, survival, logarithmic and
inverse-survival methods remain native-Fortran-only, as do special-function,
random-state, random-variate and likelihood entry points.

The generator also provides an explicit-shape C adapter for the native
assumed-shape lower-triangular solver, `scifort_solve_lower_f64`, plus Cholesky,
symmetric eigenpair and PSD decomposition adapters. It uses ordinary
Fortran array declarations and documented dimensions; see the
[linear-algebra configuration](cabi_config/linalg.json) and
[array interface contract](cabi_config/README.md).

Wrapper implementations must preserve the numerical results of the Fortran
core and must document language-specific choices involving shapes, missing
values, exceptions, and memory ownership.

The generated ABI also covers all eight PDF/log-PDF/CDF/SF/log-CDF/log-SF/PPF/ISF
methods for exponential, logistic, Laplace, Cauchy, Weibull, lognormal, gamma,
beta, Student t and chi-square. All eleven families also provide elementwise
C entry points for array-valued parameters; Python applies NumPy broadcasting.
Together with normal, these provide bulk entry points for all 56 distribution
array workloads in the SciPy timing comparison. See the
[family manifests and regeneration commands](cabi_config/README.md#benchmark-distribution-families).

A generated [NumPy/ctypes Python interface](python/README.md) now exposes these
eleven families, with scalar calls, same-shape arrays, explicit defaults and optional
caller-owned output memory. The focused development build and practical
Python-to-SciPy timing commands are documented there. R, MATLAB and Octave wrappers
remain planned.

The [discovery workflow](../tools/cabi_generator/README.md#api-inventory-and-draft-manifests)
reports public interfaces, existing C bindings and configured C/Python coverage.
The latest additions are available as `api.gamma`, `api.beta` and `api.t`, with
required shapes `a`, `a, b` and degrees of freedom `df`, respectively.

The [unified workflow](../tools/cabi_generator/README.md#unified-reviewed-interface-workflow)
now generates the dependency-ordered development build plan, aggregate C header
and [reviewed interface coverage](cabi_config/interface_coverage.md). Adding the
reviewed chi-square manifest did not require editing the development builder.

The generated C ABI also exposes four [matrix-distribution helper adapters](cabi_config/matrix_helpers.json):
SPD factorization with log determinant, lower-triangular inversion, and triangular
or SPD solves with multiple right-hand sides. Their native procedures pilot
[structured shape contracts](../docs/matrix_helper_shapes.md), checked against
reviewed wrapper dimensions and rendered as Markdown tables.

The Python owner now exposes all eight reviewed matrix adapters under `api.linalg`.
Their [Python interface contract](python/README.md) documents exact shapes,
Fortran storage, named outputs, native statuses and caller-provided buffers.

Multivariate normal PDF/log-PDF are available through generated explicit-shape C
adapters and `api.multivariate_normal.pdf/logpdf` in Python. See the
[density input/default/error contract](python/README.md#multivariate-normal-densities)
and [native shape tables](../docs/multivariate_normal_shapes.md).

Multivariate Student t PDF/log-PDF also use generated C adapters and
`api.multivariate_t.pdf/logpdf`. See their
[shape/default/singularity contract](python/README.md#multivariate-student-t-densities)
and [native shape tables](../docs/multivariate_t_shapes.md).

Dirichlet PDF, log-PDF, entropy and moments are available through generated C
adapters and `api.dirichlet`. Mean/variance/covariance return arrays with optional
caller-owned buffers. See the [Python contract](python/README.md#dirichlet-distributions)
and [native shape tables](../docs/dirichlet_shapes.md).

Generated multinomial and multivariate hypergeometric interfaces provide masses
and moments under `api.multinomial` and `api.multivariate_hypergeom`. Their
[integer/buffer contract](python/README.md#multinomial-and-multivariate-hypergeometric-distributions)
includes explicit overflow guards; see the
[multinomial](../docs/multinomial_shapes.md) and
[hypergeometric](../docs/multivariate_hypergeom_shapes.md) shape tables.

Dirichlet-multinomial PMF, log-PMF and moments are available through five generated
C adapters and `api.dirichlet_multinomial`. See the
[Python/error contract](python/README.md#dirichlet-multinomial-distribution) and
[shape tables](../docs/dirichlet_multinomial_shapes.md).

Wishart and inverse-Wishart each provide six generated C/Python methods: PDF,
log-PDF, entropy, mean, mode and entrywise variance. Read the
[matrix and moment contract](python/README.md#wishart-and-inverse-wishart-distributions)
and the [Wishart](../docs/wishart_shapes.md) and
[inverse-Wishart](../docs/invwishart_shapes.md) shape tables. Sampling continues
to use the native explicit RNG-state API.

Matrix-normal PDF/log-PDF/entropy and matrix-t PDF/log-PDF are generated for
rectangular evaluation matrices. See the [Python shape/default contract](python/README.md#matrix-normal-and-matrix-t-distributions),
[matrix-normal shapes](../docs/matrix_normal_shapes.md), and
[matrix-t shapes](../docs/matrix_t_shapes.md). The Python entropy interface uses
independent 1-by-1 identity defaults for absent row/column covariance matrices.

Poisson-binomial mass, tail, quantile and analytic log-mass score methods are available under
`api.poisson_binom`. New C symbols use `_explicit_f64` and preserve the three
existing scalar-returning symbols. See the
[probability-vector and endpoint contract](python/README.md#poisson-binomial-distribution)
and [shape tables](../docs/poisson_binom_shapes.md). Their manifest also exercises
reviewed selection of specific procedures behind public Fortran generics.

Multivariate normal and Student-t also expose entropy and marginal parameter
selection in C/Python. Python entropy infers dimensions from optional arrays;
Python marginal indices are zero-based while C/Fortran remain one-based. See the
[shapes, outputs and status contract](python/README.md#multivariate-entropy-and-marginal-parameters).

Multivariate normal maximum-likelihood fitting is available through C/Python,
with optional fixed mean/covariance and caller-owned outputs. See the
[column-oriented fitting and error contract](python/README.md#multivariate-normal-fitting).

Explicit RNG streams now support native bulk uniforms and generated Dirichlet
sampling through six caller-owned int64 state components. Python exposes
`api.rng(seed)` and `api.dirichlet.rvs_array(stream, alpha, n_samples=...)`.
Read the [state, shape and failure contract](python/README.md#explicit-rng-streams-and-dirichlet-sampling).

Generated integer-output batch samplers now cover multinomial and multivariate
hypergeometric distributions using the same caller-owned RNG stream as Dirichlet.
See the [sampling and integer buffer contract](python/README.md#integer-count-batch-sampling).

Multivariate normal and Student t batch samplers are also available through the
state bridge, with optional location/matrix parameters and explicit output
extents. See [vector sampling](python/README.md#multivariate-vector-batch-sampling).
