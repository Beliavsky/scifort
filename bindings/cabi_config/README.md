# SciFort C ABI configuration

The reusable implementation is in [tools/cabi_generator](../../tools/cabi_generator/README.md).
This directory contains SciFort-specific interface and compatibility decisions.

`normal.json` preserves the five generated scalar/uniform-bulk pairs and the
three PDF/CDF/PPF scalar/bulk pairs maintained in `src/c_api/api.f90`. It adds
eight new `scifort_normal_*_elementwise_f64` symbols. `gamma.json` adds the same
elementwise pattern for its eight methods. All arguments to these new wrappers
are explicit-shape arrays sharing one length; invalid parameter elements retain
native results with status 0. They do not apply the uniform bulk parameter checks.

- [Normal wrappers](../../src/c_api/generated/normal_generated.f90).
- [Normal header](../../include/scifort_generated_normal.h).
- [Normal coverage](normal_coverage.json).
- [Gamma header](../../include/scifort_generated_gamma.h).

Regenerate after changing manifests or native declarations, review the diff, and
run the unified pipeline with `--check` before committing. Optional arguments
are required at the C boundary. Existing uniform bulk calls still validate
finite location and finite positive scale before touching output. Buffer,
length and status contracts appear in the generated headers and generator guide.
All eleven configured distribution families now provide eight elementwise
symbols each (88 total). Python methods support NumPy broadcasting of shapes,
location and scale. Student t retains its positive-infinite-df normal limit;
chi-square and other shapes retain their finite-positive requirements.

`linalg.json` uses schema version 2 to wrap the existing assumed-shape
`solve_lower(l(:,:),b(:),x(:),status)` procedure. Its generated adapter declares
`l(n,n)`, `b(n)`, and `x(n)` and calls the native routine directly. It adds
`scifort_solve_lower_f64` without modifying the numerical implementation.

- [Explicit-shape Fortran adapter](../../src/c_api/generated/linalg_generated.f90).
- [C header](../../include/scifort_generated_linalg.h), included by `scifort.h`.
- [Coverage record](linalg_coverage.json), with explicit exclusions for other
  routines and native status constants.

The solver accepts a by-value C `int` order and contiguous column-major buffers;
all pointers must be valid and non-NULL. Orders below one return status 1 and
leave output unchanged. Otherwise, native status codes apply: 0 success and
5 singular/nonfinite diagonal, with the native routine's partial-output behavior.
No allocation, buffer copying, or pointer conversion occurs in this adapter.

The native `solve_lower` declaration documents shapes `(n,n)`, `(n)` and `(n)`,
with `n = size(b), n >= 1` in its procedure description. These match the reviewed
manifest. The [coverage record](linalg_coverage.json) has an empty
`native_shape_documentation_review` list. Missing or mismatched native shape
comments produce an advisory during generation/checking; native code is never
rewritten automatically. See the [documentation convention](../../tools/cabi_generator/README.md).


Schema-version-1 elemental manifests retain their pointer-based default. New
version-2 elemental manifests default to explicit-shape bulk wrappers. The
normal manifest pins `array_style: "pointer"` to retain its existing NULL-buffer
contract; do not remove that policy when updating its schema version.

The generated [Python interface](../python/README.md) supplies explicit defaults
and array conventions over this ABI. [python.json](python.json) configures it;
its generation and read-only check are separate from the C manifest batch.

## Benchmark distribution families

[manifests.json](manifests.json) lists all eleven configured modules. Regenerate
or check them together from the repository root:

```text
python tools/cabi_generator/generate.py bindings/cabi_config/manifests.json --batch
python tools/cabi_generator/generate.py bindings/cabi_config/manifests.json --batch --check
```

The [exponential](exponential.json), [logistic](logistic.json), [Laplace](laplace.json),
[Cauchy](cauchy.json), [Weibull](weibull.json), and [lognormal](lognormal.json)
manifests each account for all eight public probability methods: PDF, log-PDF,
CDF, SF, log-CDF, log-SF, PPF and ISF. Each emits five new scalar symbols and
eight bulk symbols, for 78 new symbols across these six families. PDF/CDF/PPF
scalar symbols remain in `src/c_api/api.f90`; `existing_scalar_symbol` records
them without emitting duplicate implementations or declarations.

Together with normal, all 56 distribution-array benchmark workloads now have
bulk entry points. New bulk adapters use explicit-shape Fortran arrays and require
valid non-NULL, nonoverlapping buffers even for zero-length calls. All use
`size_t` counts and return status 0 on success or 1 for invalid lengths/uniform
parameters, preserving output on call-level errors. Shapes (`c` for Weibull,
`s` for lognormal) and scale must be finite and positive; location must be finite.
Scalar calls retain native validation and NaN/endpoint behavior. The header for
each family documents the complete call contract and is included by `scifort.h`.

The focused family test compiles verbatim production scalar wrapper bodies
separately to avoid building unrelated distribution modules. It compares generated
bulk results to those scalar wrappers and direct native calls in checked and
optimized builds. This bounds validation to the changed ABI integration; it does
not replace full-project validation before release or merge.

## Gamma, beta and Student t

The [gamma](gamma.json), [beta](beta.json) and [Student t](student_t.json)
manifests were reviewed from discovery drafts. Each adds five scalar and eight
bulk symbols (39 new C symbols total); existing PDF/CDF/PPF scalar definitions
remain unchanged in `src/c_api/api.f90`. All three appear in the C batch and
[Python policy](python.json), bringing coverage to 80 scalar/bulk method pairs
and Python methods across ten families. This is configured family coverage,
not a claim of full-project parity. The current Python timing benchmark still
covers its original 56 workloads across seven families.

Gamma shape `a` and beta shapes `a, b` must be finite and positive. Student t
`df` must be positive but **may be positive infinity**, which delegates to the
normal distribution. Bulk validation intentionally omits `df` from the finite
list. All three require finite location and finite positive scale. Python
requires the shapes and defaults only location/scale to 0/1; it exposes Student
t as `api.t`, preserving the native `t_*` and existing `scifort_t_*` names.

A report for these three families can be regenerated without writing wrappers:

```text
python tools/cabi_generator/discover.py src/stats/gamma.f90 src/stats/beta.f90 src/stats/student_t.f90 --abi-source src/c_api --manifest-batch bindings/cabi_config/manifests.json --python-policy bindings/cabi_config/python.json --output build/cabi-discovery/after.json
```

It writes JSON and a sibling Markdown report. The reviewed scope has 24 public
methods, 24 observed scalar C wrappers, 24 observed procedural bulk wrappers,
24 configured scalar symbols, 24 configured bulk symbols and 24 Python methods.
Observed bindings come from syntactic imported references, not name matching;
configured coverage does not prove that generated artifacts or binaries are current.
Use the separate generation drift checks and focused integration tests.

A wider `src/stats` scan is available locally in
`build/cabi-discovery/distributions.json` and `.md`. It attempted 162 files,
parsed 156, recorded five parser failures and explicitly skipped
`qmc_sobol_data.f90` under a 1,000,000-byte limit. One parsed default-public module
also remains unresolved. It therefore marks itself incomplete and exits nonzero.
The resolved scope contains 1,273 public exports, with 909 supported scalar-function signatures
and 906 supported elemental signatures; these include methods beyond the eight
probability methods and are not counts of distributions or a SciPy parity score.
The report records 80 configured scalar/bulk/Python method triples. Regenerate:

```text
python tools/cabi_generator/discover.py src/stats --progress --max-source-bytes 1000000 --abi-source src/c_api --manifest-batch bindings/cabi_config/manifests.json --python-policy bindings/cabi_config/python.json --output build/cabi-discovery/distributions.json
```

The five parser failures are `foldcauchy.f90` and `skewcauchy.f90` (fparser's
intrinsic ATAN argument-count handling), plus `power_transform_engine.f90`,
`probability_plot.f90` and `sample_distance.f90` (Fortran 2018 `error stop ...,
quiet=` syntax outside the configured parser grammar). `scifort_power_transform`
has default-public access, which the current public-export resolver deliberately
rejects. These are discovery limitations; numerical source has not been rewritten
to satisfy the parser. The three reviewed families have no such inventory errors.

## Unified regeneration and chi-square

Use the reviewed [pipeline policy](pipeline.json) to regenerate/check all related
outputs together:

```text
python tools/cabi_generator/pipeline.py bindings/cabi_config/pipeline.json
python tools/cabi_generator/pipeline.py bindings/cabi_config/pipeline.json --check
```

The batch now includes twelve modules: eleven distribution families and the
existing lower-triangular solver. The [chi-square manifest](chi2.json) adds eight
bulk symbols and five scalar symbols while preserving its three manual scalar
symbols. Its degrees of freedom must be finite and positive, unlike Student t's
permitted infinite-DF limit. Python exposes it as `api.chi2` with required `df`
and explicit location/scale defaults 0/1.

The [combined coverage table](interface_coverage.md) reports 88 scalar/bulk/Python
method triples plus one array adapter. Constants/excluded solver exports remain
visible in the JSON. The generated development build plan discovers the
chi-square -> gamma -> incomplete-gamma/normal dependency chain and includes only
Python-selected families, not the full manual C API object. The old 56-workload
Python timing suite retains its existing scope and reports.

## Dense matrix decompositions

The reviewed `linalg.json` now generates four explicit-shape adapters:

| C symbol | Arrays | Native operation |
| --- | --- | --- |
| `scifort_solve_lower_f64` | L(n,n), b(n), x(n) | Forward substitution |
| `scifort_cholesky_lower_f64` | A(n,n), L(n,n) | Lower Cholesky factor |
| `scifort_symmetric_eigen_jacobi_f64` | A(n,n), eigenvalues(n), eigenvectors(n,n) | Ascending symmetric eigenpairs |
| `scifort_psd_decompose_f64` | A(n,n), eigenvalues(n), eigenvectors(n,n) | PSD eigenpairs, rank and log pseudo-determinant |

All require n >= 1 and contiguous column-major storage. Every pointer, including
scalar outputs, must be valid and non-NULL; writable buffers must not overlap
other buffers. Input matrices use the lower triangle as authoritative. The PSD
adapter also returns C `int *rank`, `double *log_pdet`, and `double *cutoff`.
Its final by-value C `int allow_singular` is required: 0 rejects rank deficiency,
1 permits it. Other values return status 1 and preserve all outputs except status.
Invalid orders have the same adapter-level preservation behavior.

Once the native procedure is called, its statuses and output-on-failure behavior
apply. Status 0 means success; 1 invalid shape; 2 Cholesky non-positive-definite;
3 PSD negative eigenvalue beyond tolerance; 4 eigensolver nonfinite input or lack
of convergence; 5 singular triangular factor or rejected PSD rank deficiency.
Outputs may then be zeroed, partially computed or complete; do not treat a
nonzero status as a successful factorization. The PSD routine sets rank/log_pdet/
cutoff to zero when eigendecomposition fails. For an indefinite matrix it returns
rank=0 and log_pdet=0 with the computed cutoff; for disallowed rank deficiency it
returns the computed rank/cutoff and log_pdet=0.

On PSD success, cutoff is `1e6*epsilon(binary64)*max(abs(eigenvalues))`; rank
counts eigenvalues above that cutoff, and log_pdet sums their logarithms. A zero
matrix permitted as singular has rank=0, cutoff=0 and log_pdet=0. Existing native
semantics, including tolerance of tiny negative eigenvalues, are preserved.
These three new adapters do not change the existing solve-lower ABI and are not
yet exposed through the generated Python distribution package.

The shape-documentation audit also covers every assumed-shape dummy declaration
in [multivariate normal](../../src/stats/multivariate_normal.f90),
[multivariate t](../../src/stats/multivariate_t.f90) and
[matrix-distribution helpers](../../src/stats/matrix_distribution_helpers.f90).
It records dimensions, optional-array relationships, column-oriented samples,
empty sample/right-hand-side batches and private helper capacity requirements.
Student t's zero-draw batch skips parameter validation, whereas the normal batch
validates its dimension and parameters even with no draws; comments preserve
that distinction. The broader repository audit remains ongoing.

Validation for this extension (CMD):

```text
build\cabi-generator-venv\Scripts\python.exe -m unittest discover -s tools/cabi_generator/tests -v
build\cabi-generator-venv\Scripts\python.exe -m unittest discover -s bindings/cabi_config/tests -p "test_generated_linalg*.py" -v
build\cabi-generator-venv\Scripts\python.exe tools/cabi_generator/pipeline.py bindings/cabi_config/pipeline.json --check
```

All 48 focused tests passed: 45 reusable generator/workflow tests and three
matrix ABI tests. The new decomposition tests were also rerun after adding
one-dimensional analytic checks. Compiled fixtures and matrix adapters used GNU
Fortran/GCC/G++ 17.0.0 20260510 (experimental) on Windows, checked with `-O0 -g
-fcheck=all -fbacktrace` and optimized with `-O3`; the small independent logical
fixture used `-O0 -fcheck=all` and `-O3`. All builds used `-std=f2018
-fimplicit-none -Werror=implicit-interface`. C clients compiled/ran under C11
with `-Wall -Wextra -Werror`, and the header/client syntax passed C++17 checks.
The read-only workflow verified all 75 artifacts. The documentation audit
covered 81 assumed-shape dummy declarations across the three native modules and
confirmed unchanged executable code. Empty right-hand-side batches skip
triangular diagonal evaluation. Source encoding/line limits, local links and
whitespace checks passed. No full FPM suite was run: numerical implementations
and global build configuration are unchanged, and the changed optional ABI/tool
paths were tested directly. Other compilers/platforms remain unverified.

## Matrix-distribution helpers and structured shapes

[matrix_helpers.json](matrix_helpers.json) adds four explicit-shape C entry points
that call the original native matrix-distribution helpers:

| Symbol | Dimensions | Array arguments |
| --- | --- | --- |
| `scifort_spd_factor_f64` | n >= 1 | A(n,n), L(n,n) |
| `scifort_invert_lower_f64` | n >= 1 | L(n,n), inverse(n,n) |
| `scifort_solve_lower_matrix_f64` | n >= 1, k >= 0 | L(n,n), B(n,k), X(n,k) |
| `scifort_solve_spd_matrix_f64` | n >= 1, k >= 0 | L(n,n), B(n,k), X(n,k) |

The last routine solves `(L L^T) X = B` using a supplied lower Cholesky factor;
it does not factor the matrix itself. `spd_factor` also returns `double *logdet`.
All storage is column-major and caller-owned. Buffers must be valid, aligned,
non-NULL and nonoverlapping, including scratch pointers for empty RHS arrays.
Scalar output/status pointers must also be valid and non-NULL. There is no
allocation or copying in the adapters; native routines may allocate their own
workspaces (for example `solve_spd_matrix`).

Adapter-invalid extents return status 1 with all outputs except status unchanged.
After a native call, `matrix_status_*` values apply: 0 success, 1 invalid shape,
2 invalid parameter, 3 factorization/triangular-diagonal failure. These differ
from the linalg module's status values; use the contract of the called routine.
Native failure outputs may be zeroed or partially computed and are preserved.
`spd_factor` zeroes L/logdet before checking all matrix entries for finiteness;
it returns status 2 for any nonfinite entry, including the upper triangle.
It returns status 3 on Cholesky failure, with zero/partial L and logdet=0.
Triangular routines use the lower triangle and reject zero/nonfinite diagonals
when evaluated. Empty RHS batches (k=0) still validate array extents, then succeed
without reading triangular diagonals or modifying the zero-size output.

The native procedure headers now carry experimental `@dimension`, `@require`
and `@when` annotations in addition to their trailing shape comments. See the
[generated contract tables](../../docs/matrix_helper_shapes.md) and the
[prototype grammar and limitations](../../tools/cabi_generator/README.md).
Structured contracts are opt-in, fail generation on manifest disagreement and do
not change numerical implementations. The new functions are currently C-facing;
the generated Python distribution package does not expose them.

Validation for the shape-contract pilot (CMD):

```text
build\cabi-generator-venv\Scripts\python.exe -m unittest discover -s tools/cabi_generator/tests -v
build\cabi-generator-venv\Scripts\python.exe -m unittest discover -s bindings/cabi_config/tests -p test_generated_matrix_helpers.py -v
build\cabi-generator-venv\Scripts\python.exe tools/cabi_generator/pipeline.py bindings/cabi_config/pipeline.json --check
```

All 51 focused tests passed: 48 generator/workflow tests and three matrix-helper
ABI tests. New cases cover strict source/manifest agreement, malformed or missing
annotations, Markdown drift, non-writing failures, native/C outputs at orders 1
and 2, three RHS columns, nonfinite/singular/indefinite inputs, invalid extents,
and empty RHS batches with otherwise-invalid diagonals. The actual native helper
module and its dependencies compiled in checked (`-O0 -g -fcheck=all -fbacktrace`)
and optimized (`-O3`) configurations with GNU Fortran/GCC/G++ 17.0.0 20260510
(experimental) on Windows. The compiler flags include `-std=f2018 -fimplicit-none
-Werror=implicit-interface`; C clients compiled and ran with C11, and C++17 syntax
checks passed, both with `-Wall -Wextra -Werror`. Independent C residual and
log-determinant checks supplement direct native-oracle comparisons.

The read-only pipeline check passed for all 79 artifacts. Source encoding,
Fortran ASCII/line lengths, local links and whitespace passed; native executable
code was verified unchanged. No full FPM suite was run because the changed scope
is optional interface/documentation tooling and adapters; no global build or
numerical implementation changed. Existing distribution artifacts remain current.
Other compilers/platforms and automatic FORD/Sphinx rendering remain unverified.

## Python matrix groups and linalg shape contracts

The Python policy's `array_groups.linalg` selects both matrix manifests. All eight
adapters are available through `api.linalg`, with explicit boolean defaults,
shape/layout/alias validation and named outputs alongside the native status.
See the [Python matrix API](../python/README.md). No new C symbols or numerical
implementations are introduced by this Python expansion.

The four native linalg procedures now also opt into structured shape contracts,
with definitions/minimums matching their manifests. Their generated
[shape tables](../../docs/linalg_shapes.md) complement the existing
[matrix-helper tables](../../docs/matrix_helper_shapes.md). The broader native
documentation audit and automatic FORD/Sphinx integration remain future work.

## Multivariate normal density adapters

[multivariate_normal.json](multivariate_normal.json) selects only PDF and log-PDF
using the new `array_functions` pattern. Their generated signatures are:

```c
void scifort_multivariate_normal_pdf_f64(int n, const double *x,
    const double *mean, const double *cov, int allow_singular,
    double *value, int *status);
void scifort_multivariate_normal_logpdf_f64(int n, const double *x,
    const double *mean, const double *cov, int allow_singular,
    double *value, int *status);
```

Inputs have shapes x(n), mean(n), cov(n,n), n >= 1, with contiguous column-major
storage, valid/aligned/non-NULL/disjoint buffers and output pointers. Every input
is required; flags are 0 false / 1 true. Invalid dimension/flag returns status 1
and preserves `value`. Otherwise status is 0 and `value` is the unchanged native
numerical result, including NaN for invalid covariance/mean or rejected singular
covariance. The lower covariance triangle is authoritative; all entries must
still be finite. These call-status codes are not `multivariate_status_*` numerical
failure codes. There are no null-pointer checks or inferred buffer lengths.

The original assumed-shape functions retain their Fortran optional arguments.
Their reviewed annotations generate [native shape tables](../../docs/multivariate_normal_shapes.md).
Python adds zero-mean/identity-covariance defaults through explicit policy and
returns scalar floats; see the [Python contract](../python/README.md#multivariate-normal-densities).
Entropy and marginalization are now also generated; see the contract below.
CDF, fitting and random sampling remain outside the reviewed interfaces.

## Multivariate Student t density adapters

[multivariate_t.json](multivariate_t.json) reuses the existing array-function
pattern for these signatures:

```c
void scifort_multivariate_t_pdf_f64(int n, const double *x,
    const double *loc, const double *shape, double df, int allow_singular,
    double *value, int *status);
void scifort_multivariate_t_logpdf_f64(int n, const double *x,
    const double *loc, const double *shape, double df,
    double *value, int *status);
```

x(n), loc(n), shape(n,n), n >= 1; every native input is required in C. Storage
is column-major, aligned and contiguous. All buffers and scalar output pointers
must be valid, non-NULL and disjoint. PDF flags are 0 false / 1 true; log-PDF
always permits singular shape and exposes no flag. Invalid n/PDF flag returns
status 1, leaving `value` unchanged. Valid calls return status 0 even for native
NaN, such as invalid df or a rejected shape matrix.

Shape is not covariance: for finite df > 2 covariance is df/(df-2) times shape.
All shape entries must be finite; the lower triangle is authoritative. Finite-df
singular formulas use projected distances, while infinite df delegates to the
native normal support check. Existing numerical behavior and limitations remain
unchanged. These are call-status codes, not `multivariate_t_status_*` failures.

See the [generated native shape tables](../../docs/multivariate_t_shapes.md)
and [Python defaults/singularity contract](../python/README.md#multivariate-student-t-densities).
Entropy and marginalization are now also generated; see the contract below.
CDF and random sampling remain outside the reviewed interfaces.

## Dirichlet density, entropy and moment adapters

[dirichlet.json](dirichlet.json) selects six native functions:

```c
void scifort_dirichlet_pdf_f64(int k, int m, const double *x,
    const double *alpha, double *value, int *status);
void scifort_dirichlet_logpdf_f64(int k, int m, const double *x,
    const double *alpha, double *value, int *status);
void scifort_dirichlet_entropy_f64(int k, const double *alpha,
    double *value, int *status);
void scifort_dirichlet_mean_f64(int k, const double *alpha,
    double *value, int *status);
void scifort_dirichlet_var_f64(int k, const double *alpha,
    double *value, int *status);
void scifort_dirichlet_cov_f64(int k, const double *alpha,
    double *value, int *status);
```

alpha(k), k>=1; densities accept x(m), m>=1, with m equal to k or k-1.
For the shortened form the native function infers the final coordinate as
1-sum(x), requiring k>=2. Scalar results use one double; mean/variance need k
doubles and covariance needs k*k in column-major order. Buffers must be valid,
aligned, contiguous, non-NULL and disjoint, including all output/status pointers.
Callers provide sufficient storage; adapters do not inspect actual buffer sizes.

Status 1 rejects invalid extents or point-length choices, preserving all result
storage. Status 0 retains native numerical results, including NaN for invalid
concentrations/simplex values. Native simplex tolerance is 1e-9; zero coordinates
with concentration below 1 produce native NaN. These call statuses are distinct
from native `dirichlet_status_*` numerical/sampling statuses. Numerical formulas
and limitations are unchanged.

See [generated native shape tables](../../docs/dirichlet_shapes.md) and
[Python return/buffer contracts](../python/README.md#dirichlet-distributions).
All assumed-shape procedures in the native module, including sampling and private
validation, now document their dimensions. Sampling remains outside this C batch
because explicit RNG state needs a separately reviewed handle design.

## Multinomial and multivariate hypergeometric adapters

[multinomial.json](multinomial.json) selects PMF, log-PMF, entropy, mean and
covariance. [multivariate_hypergeom.json](multivariate_hypergeom.json) selects
PMF, log-PMF, mean, variance and covariance. Their symbols are named
`scifort_<family>_<method>_f64`; all are void subroutines with leading C-int
category count k>=1, native argument order, result pointer and status pointer.
Scalar native integers are by-value C int; integer arrays use `const int *`.
Native real arrays and scalar/array results use double. Storage is contiguous,
column-major, aligned, non-NULL, disjoint and caller-owned, with sufficient
lengths. See their generated headers for exact prototypes.

All methods limit n to [-2147483648,2147483646]. Hypergeometric integer x/m arrays
also require sum(abs(array))<=2147483646, evaluated in widened arithmetic.
These explicit call guards prevent overflow in native sums and count increments.
Invalid extents or safety bounds return status 1 with results unchanged. Otherwise
status is 0, including native NaN for invalid parameters or zero/-infinity for
impossible observations. These are not native sampling status enumerations.
The original Fortran API and numerical formulas are unchanged.

Python integer values are checked before conversion to C int; floats and booleans
are rejected for native integer arguments. Multinomial x is a native real array
and retains its integer-valued support checks. Finite probabilities use the native
final-entry adjustment before adjusted-range validation. See the detailed
[Python input/buffer contract](../python/README.md#multinomial-and-multivariate-hypergeometric-distributions),
[multinomial shape tables](../../docs/multinomial_shapes.md), and
[hypergeometric shape tables](../../docs/multivariate_hypergeom_shapes.md).
Sampling remains native-only pending explicit RNG handle design. Both modules'
private and sampling array arguments now have readable shape documentation.

## Multivariate entropy and marginal adapters

The normal and Student-t manifests also select entropy and marginal parameter
selection. Entropy C signatures take `(n, mean, cov, value, status)` and
`(n, loc, shape, df, value, status)` respectively; the arrays are always required.
Valid scalar calls return status 0, retaining native numerical NaN/infinity;
invalid extents return 1 and preserve the result.

Marginal C signatures take `(n, m, dimensions, mean, cov, marginal_mean,
marginal_cov, status)` and the corresponding `loc/shape` names for Student t.
Indices are unique **one-based** C integers in `[1,n]`. Inputs/outputs use shapes
`dimensions(m)`, location `(n)`, matrix `(n,n)`, selected location `(m)` and
selected matrix `(m,m)`, with n,m >= 1. All pointers must be valid, aligned,
nonnull, contiguous, disjoint and sufficiently sized; matrices are column-major.
The native routine validates uniqueness and ranges (status 2); invalid extents
return 1. Failure leaves both output buffers unchanged. Native lower-triangle
and finite-entry rules are retained; selection does not test definiteness.

Python separately resolves optional entropy dimensions and converts zero-based
indices. Read the [full interface contract](../python/README.md#multivariate-entropy-and-marginal-parameters)
and the generator's [dimension-source](../../tools/cabi_generator/README.md#ordered-optional-dimension-sources)
and [index policies](../../tools/cabi_generator/README.md#python-zero-based-index-inputs).
