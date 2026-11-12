## `scipy.optimize.shgo`: multiple independent nonlinear constraints (2026-11-10)

SHGO now accepts `constraints=blocks`, an array of
`shgo_nonlinear_constraint` records with **independent** vector callbacks,
per-component lower/upper bounds, and optional analytic Jacobians. The
feasibility restorer and SLSQP adapter concatenate constraint values in
input order, preserve mixed exact/numerical derivatives, and enforce the
existing objective budget and feasibility checks. The legacy `nonlinear=`
interface remains supported unchanged (but cannot be combined with the new
`constraints=` argument). A runnable example and 40 new core assertions,
5 independent SciPy numerical comparisons, and mixed-Jacobian SLSQP interface
tests are included. Full FPM/SLSQP numerical integration remains unverified.

See [SHGO guide](docs/optimize_shgo.md),
[example](example/optimize_shgo_multiconstraints_example.f90), and
[validation](validation/optimize_shgo_multiconstraints_2026-11-10.md).

## `scipy.optimize.shgo`: analytic constraints and feasibility recovery (2026-11-10)

The bounded SHGO implementation now accepts `nonlinear_jac(x,a)` with
one Jacobian row per nonlinear constraint component, forwards properly
signed/reordered analytic derivatives to SLSQP, and optionally projects
promising mesh starts onto narrow equality and inequality regions via
`restore_feasibility=.true.`. Default SHGO behavior is unchanged. New
`ncjev`, `nrestore`, and `nrestored` counters document Jacobian/projection
work. The projection is **not an optimizer or convergence guarantee**.
See the [API guide](docs/optimize_shgo.md),
[example](example/optimize_shgo_restoration_example.f90), and
[validation](validation/optimize_shgo_restoration_2026-11-10.md).

## `scipy.optimize.shgo`: linear and nonlinear constraints (2026-11-10)

SHGO now accepts linear and nonlinear lower/upper bound constraint blocks,
including equalities, using the existing `differential_evolution` convention.
Constrained polishing defaults to SciFort SLSQP, and infeasible mesh samples
provide restoration starts without being reported as local minima. New
`ncev`, `maxcv`, `feasible`, and `shgo_infeasible` diagnostics are available.
See [API](docs/optimize_shgo.md),
[example](example/optimize_shgo_constraints_example.f90), and
[validation](validation/optimize_shgo_constraints_2026-11-10.md).
This remains a restricted SciFort SHGO variant, not a full SciPy SHGO port.

## `scipy.optimize.shgo`: conforming adaptive tetrahedralization in 3D (2026-11-10)

The existing opt-in `adaptive=.true., conforming=.true.` mode now supports
**three-dimensional** bounded optimization. An octree of dyadic leaf cubes
is tetrahedralized using shared face-center vertices, exact-integer face
triangulation and coning to each sampled cell center. It checks positive
oriented volumes, face coverage, and shared triangle incidences across
coarse-fine interfaces. Results include `ntetrahedra` and `conforming_valid`;
all new face samples obey the strict objective and vertex budgets. The
original 2D construction and default nonconforming modes are unchanged.
See the [API guide](docs/optimize_shgo.md),
[3D example](example/optimize_shgo_conforming_3d_example.f90) and
[validation](validation/optimize_shgo_conforming_3d_2026-11-10.md).
This is a **SciFort extension**, not a complete port of SciPy's SHGO.

## `scipy.optimize.shgo`: conforming adaptive mesh in 2D (2026-11-10)

The bounded SHGO subset adds `conforming=.true.` alongside
`adaptive=.true.` for **two-dimensional** optimization. Each adaptive leaf
is triangulated from a sampled center to its fully subdivided boundary;
exact dyadic area and edge checks validate conformity across coarse-fine
interfaces. The default behavior and higher-dimensional support are
unchanged. Results include `ntriangles` and `conforming_valid`. See the
[SHGO guide](docs/optimize_shgo.md),
[example](example/optimize_shgo_conforming_example.f90), and
[validation](validation/optimize_shgo_conforming_2026-11-10.md).

## `scipy.optimize.shgo`: sampled homology-rank stabilization (2026-11-10)

The optional `hgr_stability_iters` stopping rule monitors the number of
sampled simplicial-graph minima over completed rounds, requiring both
consecutive unchanged counts and no material objective improvement.
`hgr_min_iters` and `hgr_f_tol` provide safeguards against premature stopping.
Connected equal-value plateaus are now deduplicated in the adaptive graph.
This is a heuristic **SciFort extension**, not SciPy's exact `minhgrd`
stopping algorithm or proof of convergence. See the
[guide](docs/optimize_shgo.md), [example](example/optimize_shgo_hgr_example.f90),
and [validation](validation/optimize_shgo_hgr_2026-11-10.md).

## `scipy.optimize.shgo`: selective adaptive refinement (2026-11-10)

The bounded SHGO subset now also accepts `adaptive=.true.` and
`adaptive_fraction=0.25_dp`. It scores leaf cells, selectively bisects the
promising ones, reuses exact dyadic vertices, and retains the existing
local-polishing and strict-budget logic. Uniform sampling is still the default;
adaptive leaf graphs are **nonconforming**, so this is not a complete port of
SciPy's SHGO or its convergence guarantees. See
[guide](docs/optimize_shgo.md),
[example](example/optimize_shgo_adaptive_example.f90) and
[validation](validation/optimize_shgo_adaptive_2026-11-10.md).

## `scipy.optimize.shgo`: deterministic simplicial global optimization (2026-11-10)

A native `shgo` variant adds nested Kuhn-Freudenthal triangulations,
sublevel-set graph minimum discovery, bounded Nelder-Mead/Powell polishing,
multiple sorted local minima, and strict objective-call budgets. It is not a
complete port of SciPy SHGO; conforming adaptive triangulation, general
constraints, and homology-rank stopping are not yet implemented. See the
[API guide](docs/optimize_shgo.md),
[standalone validation](scripts/validate_optimize_shgo.sh), and
[example](example/optimize_shgo_example.f90).

## `scipy.optimize.direct`: deterministic global optimization (2026-11-10)

`scifort_optimize_direct` implements the original DIRECT rectangle-partitioning
method and a locally biased variant, using no random-number generator or
new dependencies. Supports finite bounds, strict evaluation budgets, early
stopping by best-rectangle volume/length or known optimum, and callbacks.
See [API guide](docs/optimize_direct.md),
[standalone validation](scripts/validate_optimize_direct.sh), and
[example](example/optimize_direct_example.f90).

## Dual annealing global optimization (`scipy.optimize`, 2026-11-10)

SciFort now provides `dual_annealing` with generalized heavy-tail visiting
steps, Tsallis acceptance, cooling and restarts, explicit bounds, reproducible
local seeds, optional callbacks, and bounded local polishing through existing
L-BFGS-B or Powell methods. `dual_annealing_core` is dependency-light and can
run without any local optimizer. See [API](docs/optimize_dual_annealing.md),
[example](example/optimize_dual_annealing_example.f90), and
[validation](validation/optimize_dual_annealing_2026-11-10.md). The method
has no new runtime dependencies or copied third-party implementation code.


## Basin-hopping global optimization (`scipy.optimize`, 2026-11-10)

SciFort now provides `basinhopping` combining randomized global steps,
Metropolis acceptance and existing Nelder-Mead/Powell/BFGS local methods.
`basinhopping_core` additionally supports independently supplied local optimizers.
The engine provides reproducible RNG, adaptive step-size, custom step/acceptance
callbacks, early termination and detailed diagnostics. See the
[API guide](docs/optimize_basinhopping.md),
[example](example/optimize_basinhopping_example.f90), and
[validation](validation/optimize_basinhopping_2026-11-10.md).
No new dependencies or copied third-party code.

## Linear programming (`scipy.optimize.linprog`, 2026-11-10)

SciFort now provides a standalone, real-double-precision **two-phase simplex**
solver for dense continuous linear programs with linear inequality/equality
constraints, finite/infinite bounds, and unrestricted variables. It detects
unbounded and infeasible problems and reports SciPy-style status codes.
The method is independent and introduces no new FPM dependencies. See the
[API and limitations](docs/optimize_linprog.md),
[example](example/optimize_linprog_example.f90), and
[validation](validation/optimize_linprog_2026-11-10.md).
This initial solver is suitable for small-to-moderate dense problems, not a
replacement for HiGHS on industrial-scale, sparse, or mixed-integer LPs.

## Nonnegative and bounded linear least squares (`scipy.optimize`, 2026-11-10)

SciFort now provides `nnls` and `lsq_linear(method="bvls")` for dense real
matrices, including unbounded, mixed-bound, and rank-deficient systems. A
separate adapter exposes the BVLS routine from SciFort's existing pinned
SLSQP dependency (`nnls_slsqp`, `lsq_linear_slsqp`). See the
[API and provenance guide](docs/optimize_linear_least_squares.md),
[standalone example](example/optimize_linear_least_squares_example.f90), and
[validation](validation/optimize_linear_least_squares_2026-11-10.md).
The tested native active-set method uses the existing SciFort Jacobi SVD;
`trf`, sparse/linear-operator inputs, and arbitrary batch dimensions are not
yet included. No new dependencies were added.

## Matrix exponential (`scipy.linalg.expm`, 2026-11-10)

`scifort_linalg_expm` adds a standalone real/complex dense matrix exponential
with one leading batch dimension, using a self-contained [13/13] Pade
scaling-and-squaring algorithm and the existing pivoted LU solver. It avoids
computing a Fréchet derivative when only `exp(A)` is required. See the
[API guide](docs/linalg_expm.md), [coverage audit](docs/linalg_coverage_audit_2026-11-10.md),
[example](example/linalg_expm_example.f90), and
[validation](validation/linalg_expm_2026-11-10.md).
No existing numerical source modules were changed or dependencies added.

### Differential evolution global optimization

The `scifort_optimize` facade now exports `differential_evolution`, a bounded,
derivative-free optimizer with twelve common mutation/crossover strategies,
Latin-hypercube initialization, reproducible seeds, and progress/status
reporting. See [the differential-evolution guide](docs/optimize_differential_evolution.md)
for the Fortran API, examples, supported SciPy features, testing, and attribution.


## Matrix-free randomized interpolative SVD (`scipy.linalg.interpolative.svd`, 2026-11-10)

`scifort_linalg_interpolative_svd_operator` computes real and complex
fixed-rank or adaptive thin SVDs of implicit matrices from forward and
adjoint callbacks, without materializing the original dense matrix. It
reuses the existing matrix-free interpolative decomposition and the
small-core ID-to-SVD conversion. Adaptive error estimates use independent
random probes and a relative **Frobenius** criterion. See the
[API guide](docs/linalg_interpolative_svd_operator.md),
[implicit-operator example](example/linalg_interpolative_svd_operator_example.f90),
and [validation](docs/linalg_interpolative_svd_operator_validation_2026-11-10.md).
No existing numerical modules or dependencies were modified.

## Matrix-free randomized interpolative decomposition (`scipy.linalg.interpolative`, 2026-11-10)

`scifort_linalg_interpolative_operator` adds real and complex fixed-rank and
adaptive interpolative decompositions from forward and adjoint callbacks.
The routine selects original columns and refits interpolation coefficients
without storing the full matrix. Adaptive residual checks use probabilistic
Gaussian probes; unlike the dense ID, they are **not exact certificates**.
See the [API guide](docs/linalg_interpolative_operator.md),
[implicit-operator example](example/linalg_interpolative_operator_example.f90),
and [validation report](docs/linalg_interpolative_operator_validation_2026-11-10.md).

## Randomized interpolative decomposition (`scipy.linalg.interpolative`, 2026-11-10)

`scifort_linalg_interpolative_random` adds real and complex randomized
interpolative decompositions using Gaussian row sketches, pivoted QR, and
least-squares coefficient fitting against the original matrix. It offers
fixed-rank and adaptive relative-Frobenius-error APIs with local seeded RNG,
and returns one-based selected-column permutations and interpolation weights
without computing an SVD. See the [API guide](docs/linalg_interpolative_random.md),
[example](example/linalg_interpolative_random_example.f90), and
[validation](docs/linalg_interpolative_random_validation_2026-11-10.md).
No existing interfaces or runtime dependencies were changed.

## Randomized interpolative SVD (`scipy.linalg.interpolative.svd`, 2026-11-10)

`scifort_linalg_interpolative_svd_random` adds real and complex randomized,
ID-based thin SVD with fixed-rank and adaptive relative-Frobenius-tolerance
interfaces. It uses a small Gaussian row sketch to select columns, refits
interpolation coefficients to the original matrix, and reuses the existing
ID-to-SVD conversion. Local `int64` seeds make results reproducible; no new
dependencies or existing source changes were required. See the
[API guide](docs/linalg_interpolative_svd_random.md),
[example](example/linalg_interpolative_svd_random_example.f90), and
[validation report](docs/linalg_interpolative_svd_random_validation_2026-11-10.md).

## Interpolative SVD (`scipy.linalg.interpolative.svd`, 2026-11-09)

`scifort_linalg_interpolative_svd` adds fixed-rank and relative-tolerance
ID-based thin SVD for real and complex matrices. The deterministic algorithm
reuses SciFort's existing pivoted QR, ID, and small-core SVD without adding a
dependency or changing existing APIs. The right singular-vector factor is V,
not V^H. See [API guide](docs/linalg_interpolative_svd.md),
[example](example/linalg_interpolative_svd_example.f90), and
[validation](docs/linalg_interpolative_svd_validation_2026-11-09.md).

## Randomized matrix-rank estimation (`scipy.linalg.interpolative`, 2026-11-09)

`scifort_linalg_estimate_rank` adds real/complex numerical-rank estimates for
dense matrices and matrix-free forward-product callbacks, with optional local
seeds and oversampling controls. The adaptive randomized range finder does not
add SciPy's typical dense-rank overshoot. See the
[API guide](docs/linalg_estimate_rank.md),
[example](example/linalg_estimate_rank_example.f90), and
[reference validation](docs/linalg_estimate_rank_validation_2026-11-09.md).

## Randomized spectral-norm estimates (`scipy.linalg.interpolative`, 2026-11-09)

`scifort_linalg_spectral_norm` adds `estimate_spectral_norm` and
`estimate_spectral_norm_diff` for real and complex dense matrices, as well as
callback-based operators (`A*x`, `A^H*y`) for matrix-free computations.
The power-method estimators support deterministic seeds, caller-supplied
starting vectors, rectangular/empty matrices, and failure status reporting.
See the [API guide](docs/linalg_spectral_norm.md),
[validation](docs/linalg_spectral_norm_validation_2026-11-09.md), and
[example](example/linalg_spectral_norm_example.f90).

## Randomized Clarkson-Woodruff matrix sketches (2026-11-09)

`scifort_linalg_clarkson_woodruff` adds SciPy-style CountSketch random
projections for dense real/complex matrices, one leading batch dimension,
and CSR data stored with one-based offsets/indices. Sketches can reuse an
explicit signed bucket mapping or a local reproducible int64 seed, without
new dependencies. See the [API guide](docs/linalg_clarkson_woodruff.md),
[example](example/linalg_clarkson_woodruff.f90), and
[independent SciPy validation](docs/linalg_clarkson_woodruff_validation_2026-11-09.md).

## Interpolative decomposition (`scipy.linalg.interpolative`, 2026-11-09)

`scifort_linalg_interpolative` adds deterministic pivoted-QR interpolative
matrix decompositions, selected-column and interpolation-matrix reconstruction,
and conversion of an ID to a thin SVD, for real and complex matrices. The
fixed-rank and relative Frobenius-tolerance interfaces reuse the existing
SciFort QR and SVD code without new dependencies. See the
[API guide](docs/linalg_interpolative.md),
[example](example/linalg_interpolative_example.f90), and
[validation report](docs/linalg_interpolative_validation_2026-11-09.md).

## Hermitian pseudoinverse (`scipy.linalg.pinvh`, 2026-11-09)

`scifort_linalg_pinvh` computes Moore--Penrose pseudoinverses of real
symmetric and complex Hermitian matrices using SciFort's existing `eigh`
solver. It includes SciPy-compatible absolute and relative eigenvalue
cutoffs, numerical rank reporting, upper/lower triangle selection, and
one leading batch dimension. See the [API guide](docs/linalg_pinvh.md),
[example](example/linalg_pinvh_example.f90), and
[validation](docs/linalg_pinvh_validation_2026-11-09.md).

## Real block eigenvalue conversion (`scipy.linalg.cdf2rdf`, 2026-11-09)

`scifort_linalg_cdf2rdf` converts complex conjugate eigenpairs of a real
matrix into a real block eigenvalue matrix and a real invariant basis.
It supports individual or one-batch-dimension matrices and accepts
complex eigenvectors with independently phased conjugate partners.
See [API guide](docs/linalg_cdf2rdf.md),
[example](example/linalg_cdf2rdf_example.f90), and
[reference validation](docs/linalg_cdf2rdf_validation_2026-11-09.md).

## Pivoted LU, determinant, and inverse (`scipy.linalg.lu/det/inv`, 2026-11-09)

New `scifort_linalg_lu` offers rectangular and singular real/complex LU
factorizations, a one-based permutation-index form, and the `P*L` form.
`scifort_linalg_det_inv` adds determinant and general-matrix inverse APIs
reusing the dense LU backend. See the [API guide](docs/linalg_lu_det_inv.md),
[example](example/linalg_lu_det_inv_example.f90), and [independent SciPy
validation](docs/linalg_lu_det_inv_validation_2026-11-09.md).

## RQ factorization (`scipy.linalg.rq`, 2026-11-09)

New `scifort_linalg_rq` implements real/complex full and economic RQ
factorizations, an efficient R-only mode, and one leading batch dimension.
It reuses SciFort's existing Householder QR implementation with no new
external dependencies. Includes [API](docs/linalg_rq.md),
[example](example/linalg_rq_example.f90), and [SciPy validation](docs/linalg_rq_validation_2026-11-09.md).

## Direct Householder QR multiplication (`scipy.linalg.qr_multiply`, 2026-11-09)

`scifort_linalg_qr_multiply` adds real/complex `qr_multiply` (including mixed
real/complex multipliers, vector or matrix input, both multiplication sides,
optional pivoting and complex conjugation). It also provides reusable
`qr_apply` for previously computed QR factors, without explicitly forming Q.
[API](docs/linalg_qr_multiply.md),
[example](example/linalg_qr_multiply_example.f90), and
[SciPy validation](docs/linalg_qr_multiply_validation_2026-11-09.md).

## Cosine–sine decomposition (`scipy.linalg.cossin`, 2026-11-09)

New [`scifort_linalg_cossin`](docs/linalg_cossin.md) performs full real
orthogonal and complex unitary cosine–sine decompositions for all valid
`p`, `q` partitions. It supports both SciPy sign conventions and separate
subfactors/angles, using SciFort's existing dependency-free Jacobi SVD.
[Example](example/linalg_cossin_example.f90) and
[SciPy validation](docs/linalg_cossin_validation_2026-11-09.md).

New [`scifort_linalg_structure`](docs/linalg_structure_diagsvd.md) and
[`scifort_linalg_diagsvd`](docs/linalg_structure_diagsvd.md) implement
`scipy.linalg.bandwidth`, `issymmetric`, `ishermitian` and `diagsvd`.
These routines support 32/64-bit integers, real/complex matrices, and a
single batch dimension, with precise SciPy 1.17 edge-case behavior. See the
[native example](example/linalg_structure_diagsvd_example.f90) and
[reference validation](docs/linalg_structure_diagsvd_validation_2026-11-09.md).

## Block-diagonal and Khatri–Rao operations (2026-11-09)

New [`scifort_linalg_block_diag`](docs/linalg_block_diag_khatri_rao.md) and
[`scifort_linalg_khatri_rao`](docs/linalg_block_diag_khatri_rao.md) implement
`scipy.linalg.block_diag` and `scipy.linalg.khatri_rao`, with scalar/vector/
matrix blocks, mixed integer/real/complex types, and one broadcastable batch
dimension. The [example](example/linalg_block_diag_khatri_rao_example.f90)
and [SciPy reference validation](docs/linalg_block_diag_khatri_rao_validation_2026-11-09.md)
show how to use the caller-allocated output interfaces. No new dependencies.

## Hadamard, Leslie and DFT matrix constructors (2026-11-09)

New modules [`scifort_linalg_hadamard`](docs/linalg_hadamard_leslie_dft.md),
[`scifort_linalg_leslie`](docs/linalg_hadamard_leslie_dft.md), and
[`scifort_linalg_dft`](docs/linalg_hadamard_leslie_dft.md) implement
`scipy.linalg.hadamard`, `scipy.linalg.leslie`, and `scipy.linalg.dft`.
These provide exact Sylvester matrices, age-structured population matrices
(including one broadcastable batch dimension), and explicitly constructed
DFT matrices with SciPy normalization options. See the
[example](example/linalg_hadamard_leslie_dft_example.f90) and
[validation results](docs/linalg_hadamard_leslie_dft_validation_2026-11-09.md).
No new dependencies were introduced.

## Pascal and inverse Pascal matrices (2026-11-09)

New [`scifort_linalg_pascal`](docs/linalg_pascal.md) implements
`scipy.linalg.pascal` and `scipy.linalg.invpascal` for symmetric, lower and
upper matrices. Outputs include real64, overflow-checked exact int64, and
arbitrary-precision decimal strings without additional dependencies.
See the [example](example/linalg_pascal_example.f90) and
[474-case SciPy reference tests](docs/linalg_pascal_validation_2026-11-09.md).

## Helmert and Hilbert matrix constructors (2026-11-09)

New [`scifort_linalg_helmert`](docs/linalg_special_matrices.md) and
[`scifort_linalg_hilbert`](docs/linalg_special_matrices.md) implement
`scipy.linalg.helmert`, `hilbert`, and `invhilbert`. Inverse Hilbert
results support floating point, exact `int64` through order 14, and an
arbitrary-precision decimal-string option for larger orders. See the
[example](example/linalg_special_matrices_example.f90) and
[validation](docs/linalg_special_matrices_validation_2026-11-09.md).

## Polynomial companion and Fiedler matrices (2026-11-09)

New [`scifort_linalg_companion`](docs/linalg_companion.md) implements
`scipy.linalg.companion` and `scipy.linalg.fiedler_companion` for real,
complex and integer polynomials, including one batch dimension.
[`scifort_linalg_fiedler`](docs/linalg_companion.md) provides symmetric
Fiedler distance matrices. Both modules have no new dependencies.
The [example](example/linalg_companion_example.f90) and
[validation](docs/linalg_companion_validation_2026-11-09.md) include 3,420
independent SciPy comparisons across checked and optimized builds.

## Convolution matrices and FFT multiplication (2026-11-09)

New [`scifort_linalg_convolution_matrix`](docs/linalg_convolution_matrix.md)
implements `scipy.linalg.convolution_matrix` for all three modes (`full`,
`same`, `valid`), including one batch dimension and integer coefficients.
The [`matmul_convolution_matrix`](docs/linalg_convolution_matrix.md) SciFort
extension uses the existing Toeplitz FFT backend for structured products
without allocating the full matrix. See the [example](example/linalg_convolution_matrix_example.f90)
and [independent SciPy checks](docs/linalg_convolution_matrix_validation_2026-11-09.md).

## Hankel construction and FFT-based multiplication (2026-11-09)

New [`scifort_linalg_hankel`](docs/linalg_hankel.md) implements
`scipy.linalg.hankel`, including mixed real/complex last rows. The
additional SciFort [`matmul_hankel`](docs/linalg_hankel.md) extension
multiplies real or complex Hankel matrices without constructing them,
reusing the Toeplitz FFT implementation. The [example](example/linalg_hankel_example.f90)
and [validation report](docs/linalg_hankel_validation_2026-11-09.md) cover
3,600 independent numerical cases and 452 direct constructor comparisons.

## FFT-based Toeplitz matrix multiplication (2026-11-09)

New [`scifort_linalg_matmul_toeplitz`](docs/linalg_matmul_toeplitz.md)
implements `matmul_toeplitz` for square and rectangular Toeplitz products
with vectors or multiple columns, real/complex/mixed inputs, and SciPy's
optional conjugate first row. It uses an FFT circulant embedding without
densely constructing the matrix. See the
[example](example/linalg_matmul_toeplitz_example.f90) and
[3,600 independent SciPy/NumPy checks](docs/linalg_matmul_toeplitz_validation_2026-11-09.md).

## Circulant construction and FFT solves (2026-11-09)

The new [`scifort_linalg_circulant`](docs/linalg_circulant.md) module
adds `circulant` construction and FFT-based `solve_circulant`, including
real/complex and mixed inputs, multiple right-hand sides, and singular
least-squares solutions. See the [example](example/linalg_circulant_example.f90)
and [validation report](docs/linalg_circulant_validation_2026-11-09.md).

## Toeplitz construction and structured linear solves (2026-11-09)

New [`scifort_linalg_toeplitz`](docs/linalg_toeplitz.md) provides
`toeplitz` (including rectangular matrices) and `solve_toeplitz` for real,
complex and mixed-type vector or matrix right-hand sides. An O(n^2)
structured recurrence handles ordinary systems; difficult leading pivots
trigger a fallback to SciFort's existing pivoted dense LU. The result was
checked against SciPy and NumPy with 3,200 independent cases in checked and
optimized builds. See the [example](example/linalg_toeplitz_example.f90)
and [validation report](docs/linalg_toeplitz_validation_2026-11-09.md).

## LDL-based symmetric/Hermitian linear-system solves (2026-11-09)

The new [`scifort_linalg_indefinite_solve`](docs/linalg_indefinite_solve.md)
module implements `solve_indefinite` (the `sym`/`her` specializations of
`scipy.linalg.solve`) and `ldl_solve` for reusing SciFort's Bunch–Kaufman
factors with vector or matrix right-hand sides. It supports real symmetric,
complex Hermitian, and complex symmetric matrices, both triangle choices,
and 1x1/2x2 pivots. See the [executable example](example/linalg_indefinite_solve_example.f90)
and [SciPy reference validation](docs/linalg_indefinite_solve_validation_2026-11-09.md).

## Symmetric and Hermitian indefinite LDL factorization (`scipy.linalg.ldl`, 2026-11-09)

The new [`scifort_linalg_ldl`](docs/linalg_ldl.md) module provides pivoted
real symmetric, complex Hermitian, and complex symmetric `ldl` decompositions,
including 1x1 and 2x2 pivots, one-based permutation arrays, and upper/lower
triangle selection. It does not require positive definiteness and preserves
the SciPy reconstruction convention. See the [example](example/linalg_ldl_example.f90)
and [independent SciPy checks](docs/linalg_ldl_validation_2026-11-09.md).

## Ordered generalized complex Schur (`scipy.linalg.ordqz`, 2026-11-09)

The new [`scifort_linalg_ordqz`](docs/linalg_ordqz.md) module adds complex
QZ reordering for real, complex and mixed-kind matrix pairs. It supports
`lhp/rhp/iuc/ouc` and custom homogeneous eigenvalue selectors, including
infinite eigenvalues, and reorders previously computed Schur factors without
rerunning the QZ algorithm. The factorization still has the underlying
order-40 QZ limit; real generalized Schur output is not yet supported.

## Generalized nonsymmetric eigenproblems (`scipy.linalg.eig/eigvals`, 2026-11-08)

`scifort_linalg_generalized_eigen` implements real/complex/mixed matrix-pencil
`eig(a,b,w,status,vr=...,vl=...,alpha=...,beta=...)` and eigenvalue-only
`eigvals(a,b,w,status,...)`, including infinite roots and homogeneous eigenvalue
pairs. It reuses the two-sided QZ algorithm, with a stable triangular-pair
shortcut and an order-40 limit. Ordinary `eig`/`eigvals` remain unchanged in
`scifort_linalg_general_eigen`. See the [API](docs/linalg_generalized_eigen.md),
[example](example/linalg_generalized_eigen_example.f90), and
[validation](docs/linalg_generalized_eigen_validation_2026-11-08.md).

## Generalized complex Schur factors (`qz`, partial, 2026-11-08)

`scifort_linalg_qz.qz` provides unitary generalized Schur factors of real,
complex, or mixed coefficient pairs, including many singular-B and regular
both-singular cases. The new `method='iterative'` uses inverse-free
Hessenberg–triangular QZ iteration. The previous inverse-based construction
remains available with `method='inverse'`; the default `method='auto'`
tries iterative first, then falls back on numerical failure. Both modes
retain the order-40 limit and limitations on difficult pencils.
See the [API](docs/linalg_qz.md), [example](example/linalg_qz_example.f90),
and [SciPy comparison report](docs/linalg_qz_two_sided_validation_2026-11-08.md).

## Real-to-complex Schur conversion (`scipy.linalg.rsf2csf`, 2026-11-08)

`scifort_linalg_schur.rsf2csf` converts existing real quasi-triangular
Schur factors (`Tr`, `Zr`) into complex triangular factors (`Tc`, `Zc`),
without rerunning Schur QR. It handles complex-conjugate 2x2 blocks using
local unitary transformations. The converter can process valid precomputed
Schur factors of orders greater than 40, though the underlying SciFort
real-Schur solver is still limited to order 40. See
[API](docs/linalg_schur.md),
[example](example/linalg_rsf2csf_example.f90), and
[validation](docs/linalg_rsf2csf_validation_2026-11-08.md).

## Real and complex Schur decomposition and reordering (2026-11-08)

`scifort_linalg_schur.schur` now supports **real quasi-triangular output**
for real inputs (real `T` and `Z` arguments), in addition to its existing
complex Schur output for real or complex inputs. The real form preserves
complex-conjugate pairs as 2x2 diagonal blocks; both forms support built-in
region sorting and a custom eigenvalue predicate. The interfaces are
backward compatible. See [API and limitations](docs/linalg_schur.md) and
[real Schur validation](docs/linalg_real_schur_validation_2026-11-08.md).

## General nonsymmetric eigenproblems (`scipy.linalg`, 2026-11-08)

`scifort_linalg_general_eigen` adds `eig` and `eigvals` for ordinary
real/complex nonsymmetric matrices, including optional left and right complex
eigenvectors. The dependency-free implementation reuses SciFort's Schur QR
engine and has a 40-by-40 limit. See the [API](docs/linalg_general_eigen.md),
[example](example/linalg_general_eigen_example.f90), and
[reference validation](docs/linalg_general_eigen_validation_2026-11-08.md).

## Real banded bulge kernel and automatic full-spectrum eigenvalues (2026-11-08)

Real symmetric band-to-tridiagonal reduction now uses a separate real-arithmetic
Givens kernel rather than complex working arrays. For `eigvals_banded` and
values-only `eig_banded`, `method='auto'` selects band-preserving reduction
for narrow bands and orders at least 192, retaining dense fallback on numerical
failure. Full eigenvectors and selected eigenpairs keep their previous defaults.
See the [validation and benchmark](docs/linalg_banded_bulge_real_auto_validation_2026-11-08.md)
and [API](docs/linalg_banded_eigen.md).

## Selected symmetric-tridiagonal eigenvalues (2026-11-08)

`scifort_linalg_banded_eigen` now computes selected tridiagonal eigenvalues
using Sturm-sequence bisection in O(n) workspace. It supports one-based
`subset_by_index` and half-open `subset_by_value` modes, with an optional
`num_selected` output. Full-spectrum QL and existing interfaces are preserved;
selected eigenvectors use inverse iteration with a small-problem QL fallback. See [API](docs/linalg_banded_eigen.md),
[example](example/linalg_tridiagonal_selection_example.f90), and
[validation](docs/linalg_tridiagonal_selection_validation_2026-11-08.md).

## Dense symmetric/Hermitian eigenproblems (`scipy.linalg`, 2026-11-08)

`scifort_linalg_dense_eigen` implements real/complex `eigh` and
`eigvalsh`, including positive-definite generalized types 1–3. It supports
upper/lower triangle selection, one-based inclusive `subset_by_index`,
and half-open `subset_by_value` intervals with returned `num_selected`.
The routines reuse SciFort's Householder and tridiagonal QL algorithms,
without additional numerical dependencies. See the
[API](docs/linalg_dense_eigen.md),
[value-selection guide](docs/linalg_dense_eigen_value_selection.md),
[example](example/linalg_eigen_by_value_example.f90),
[interval validation](docs/linalg_dense_eigen_value_validation_2026-11-08.md),
and [full-spectrum validation](docs/linalg_dense_eigen_validation_2026-11-07.md).

## Symmetric tridiagonal and Hermitian band eigenvalue routines

`scifort_linalg_banded_eigen` now provides `eigh_tridiagonal`,
`eigvalsh_tridiagonal`, `eig_banded`, and `eigvals_banded` for real symmetric
tridiagonal and real/complex Hermitian band matrices in SciPy-compatible
packed storage. General-band calls now support band-preserving selected eigenvalues and
selected eigenvectors, using Sturm bisection and pivoted banded inverse
iteration; full-spectrum calls default to the dense method, with an opt-in
`method='inverse'` implementation that avoids expanding the input to dense.
See [the API and limitations](docs/linalg_banded_eigen.md).

`scifort_linalg_balance.matrix_balance` adds real/complex similarity
balancing with optional zero-coupling permutations and power-of-two scaling.
It returns either a full `T` or one-based permutation and scale vectors.
See [API](docs/linalg_balance.md),
[example](example/linalg_balance_example.f90), and
[validation](docs/linalg_balance_validation_2026-11-07.md).

## Dense upper Hessenberg reduction (`scipy.linalg`, 2026-11-07)

`scifort_linalg_hessenberg.hessenberg` now computes real and complex
upper Hessenberg reductions with optional orthogonal/unitary similarity
matrices. The dependency-free Householder implementation supports empty and
scaled inputs. See [API](docs/linalg_hessenberg.md),
[example](example/linalg_hessenberg_example.f90), and
[SciPy comparisons](docs/linalg_hessenberg_validation_2026-11-07.md).

## Schur-based Sylvester and Lyapunov solvers (2026-11-07)

`scifort_linalg_matrix_equations` now supports optional `method='schur'` for
`solve_sylvester`, `solve_continuous_lyapunov`, and `solve_discrete_lyapunov`,
using compact Schur reductions and triangular recurrences for matrices up to
order 40. The original pivoted Kronecker solver remains available via
`method='direct'`; automatic selection is the default. See the
[API](docs/linalg_matrix_equations.md), [example](example/linalg_matrix_equations_schur_example.f90)
and [independent SciPy tests](docs/linalg_matrix_equations_schur_validation_2026-11-07.md).

## Orthogonal and unitary Procrustes alignment (2026-11-07)

`scifort_linalg_procrustes` adds real/complex and mixed-type
`orthogonal_procrustes` for matrix alignment by an orthogonal or unitary
factor. It reuses the existing SVD with scaled cross-products, including
rank-deficient and empty matrices. See the [API](docs/linalg_procrustes.md),
[example](example/linalg_procrustes_example.f90), and
[SciPy comparisons](docs/linalg_procrustes_validation_2026-11-07.md).

## Spectral diagnostics and faster `svdvals` (2026-11-07)

`scifort_linalg_svd.svdvals` now uses a dedicated singular-value-only Jacobi
path rather than computing and discarding SVD vectors. New real/complex
`cond2` and `matrix_rank` interfaces provide 2-norm condition estimates and
numerical ranks (SciFort extensions corresponding to NumPy methods). See
[SVD API](docs/linalg_svd.md), [validation](docs/linalg_svd_condition_validation_2026-11-07.md),
and [example](example/linalg_svd_condition_example.f90).

## Full dense SVD, singular values and rectangular diagonal matrices (2026-11-07)

`scifort_linalg_svd.svd` now supports `full_matrices=.true.` for full square
real/complex U and Vh factors; the economy-size default remains unchanged.
The new real/complex `svdvals` and real `diagsvd` utilities follow SciPy
mathematical conventions without adding external numerical dependencies.
See the [SVD API](docs/linalg_svd.md),
[example](example/linalg_svd_full_example.f90), and
[independent SciPy checks](docs/linalg_svd_full_validation_2026-11-07.md).

## Dense column spaces, null spaces and principal angles (2026-11-07)

`scifort_linalg_subspaces` adds real/complex `orth`, `null_space`, and
`subspace_angles` (also supporting mixed real/complex inputs). It uses
SciFort's existing economy SVD, with orthonormal completion for wide null
spaces and a stable sine/cosine principal-angle formulation. See the
[API](docs/linalg_subspaces.md), [runnable example](example/linalg_subspaces_example.f90),
and [SciPy validation](docs/linalg_subspaces_validation_2026-11-07.md).

## Dense economy SVD, pseudoinverse, and rank-deficient least squares (2026-11-07)

`scifort_linalg_svd` adds real/complex `svd` (economy shapes), `pinv` and
`svd_lstsq`. The existing `scifort_linalg_qr.lstsq` now accepts optional
`method='auto'|'qr'|'svd'` (default `'auto'`): full-rank problems continue
using QR, while a QR rank-deficiency detection invokes the minimum-norm
SVD solver. No new BLAS/LAPACK dependency is introduced. See the
[SVD API](docs/linalg_svd.md), [example](example/linalg_svd_example.f90)
and [SciPy validation](docs/linalg_svd_validation_2026-11-07.md).

## QR updates, row and column edits (2026-11-09)

`scifort_linalg_qr_updates` implements real/complex `qr_update` (rank-one and
rank-k), `qr_insert`, and `qr_delete` for complete and economy QR factors.
Low-rank updates and column edits reuse QR factors through reduced-core
factorizations; row edits use a Householder refactorization. Fortran indices
are one-based and all outputs are caller-allocated. See the
[API](docs/linalg_qr_updates.md),
[example](example/linalg_qr_updates_example.f90), and
[independent SciPy checks](docs/linalg_qr_updates_validation_2026-11-09.md).

## Dense Householder QR and least-squares solves (2026-11-07)

`scifort_linalg_qr` provides real/complex reduced or complete `qr`, compact
`qr_factor`/`qr_unpack`, reusable pivoted `qr_solve`, and QR/SVD `lstsq`
(including minimum-norm underdetermined solves). Optional column pivoting uses
one-based final column permutations. See the [QR API](docs/linalg_qr.md),
[runnable example](example/linalg_qr_example.f90), and
[independent SciPy checks](docs/linalg_qr_validation_2026-11-07.md).

## Dense Cholesky factorization and reusable positive-definite solves (2026-11-07)

`scifort_linalg_dense` now supports real/complex `cholesky`, `cho_factor`,
and `cho_solve` for Hermitian/symmetric positive-definite matrices. Upper and
lower triangular factors are supported, with both vector and matrix RHSs and
reusable factors. The ignored output triangle is zeroed. See the
[API](docs/linalg_dense.md), [example](example/linalg_dense_cholesky_example.f90),
and [SciPy comparisons](docs/linalg_dense_cholesky_validation_2026-11-07.md).

## Dense triangular systems and reusable LU (`scipy.linalg`, 2026-11-07)

`scifort_linalg_dense` adds real and complex `solve_triangular`, `lu_factor`,
and `lu_solve`: normal, transpose, and conjugate-transpose solves with vectors
or multiple right-hand sides, plus reusable LU factors with one-based pivots.
See [API](docs/linalg_dense.md) and [SciPy validation](docs/linalg_dense_validation_2026-11-07.md).

## Banded linear systems (`scipy.linalg`, 2026-11-07)

`scifort_linalg_banded` adds real and complex `solve_banded`,
`solveh_banded`, `cholesky_banded`, and `cho_solve_banded`. The general solver
uses pivoted banded LU; the positive-definite solver uses Cholesky and allows
factor reuse for multiple right-hand sides. It preserves packed storage and
solves a 50,000-order tridiagonal system without creating a dense matrix.
See [API](docs/linalg_banded.md) and [independent SciPy validation](docs/linalg_banded_validation_2026-11-07.md).

## General analytic dense matrix functions (2026-11-07)

`scifort_linalg_funm` adds real/complex `funm` for general analytic scalar
callbacks on dense matrices up to order 40. It uses the existing complex Schur
kernel and a block Schur--Parlett/Sylvester recurrence. Global unitary Schur
reordering supports diagonalizable repeated eigenvalues **without** derivative
callbacks, and optional derivatives handle coupled Jordan blocks. A real-only
`funm_real` wrapper and explicit status/commutator diagnostics are included. See the
[API and limits](docs/linalg_funm.md) and [latest validation](docs/linalg_funm_block_parlett_validation_2026-11-07.md).

## Exponential matrix condition number (2026-11-07)

`scifort_linalg_expm_condition` adds real/complex `expm_cond`, an estimate of
SciPy's relative Frobenius condition number of `exp(A)`. It reuses SciFort's
Fréchet derivative kernel with matrix-free Lanczos and exact diagonal
special cases. See [API and numerical limits](docs/linalg_expm_condition.md)
and [validation](docs/linalg_expm_condition_validation_2026-11-07.md).

## Dense matrix trigonometric functions (2026-11-07)

`scifort_linalg_trigonometric` adds real/complex `sinm`, `cosm`, `tanm`,
`sinhm`, `coshm`, and `tanhm`, using paired Taylor series and double-angle
scaling, plus singularity screening for tangent functions. See
[API and limitations](docs/linalg_trigonometric.md) and
[validation](docs/linalg_trigonometric_validation_2026-11-07.md).

## Dense matrix logarithms and fractional powers (2026-11-07)

`scifort_linalg_logarithm` implements real/complex `logm` and
`fractional_matrix_power`, with real-only output wrappers, inverse
scaling-and-squaring and integer powers. It reuses the existing matrix
square-root and matrix-exponential modules. See the
[API and limitations](docs/linalg_logarithm.md) and
[SciPy comparison results](docs/linalg_logarithm_validation_2026-11-07.md).

## Dense polar decomposition and exponential derivatives (2026-11-06)

Real and complex `polar` (rectangular, both sides) and `expm_frechet`
(directional derivative of the dense matrix exponential) are available as
self-contained modern Fortran implementations.
See [API](docs/linalg_polar_frechet.md) and
[validation](docs/linalg_polar_frechet_validation_2026-11-06.md).

## Dense matrix square root and matrix sign (2026-11-06)

`scifort_linalg_matrix_functions` provides real/complex binary64 `sqrtm`,
`signm`, and real-only `sqrtm_real` for small dense matrices. It uses a
self-contained complex Schur method and scaled Newton iteration. See
[API, numerical limits, and reproducible tests](docs/linalg_matrix_functions.md),
including [comparison results](docs/linalg_matrix_functions_validation_2026-11-06.md)
and a [review of reusable Fortran matrix-function libraries](docs/linalg_external_library_review_2026-11-06.md).

## Dense algebraic Riccati solvers (2026-11-06)

`scifort_linalg_riccati` adds real/complex `solve_continuous_are` and
`solve_discrete_are`, including generalized descriptor and state-control cross
terms. The standalone matrix-sign and doubling implementations work on up to
20 state variables, validate the original Riccati residual, and avoid new
BLAS/LAPACK dependencies. See [API, assumptions, and examples](docs/linalg_riccati.md)
and [independent SciPy validation](validation/linalg_riccati_2026-11-06.md).

## Complex general matrix-function actions (2026-11-06)

`scifort_sparse_linalg_complex_funm` adds `funm_multiply_krylov_complex`
for complex CSR/CSC matrices and matrix-free operators, with complex scaling,
Arnoldi/Lanczos restarts, multiple right-hand sides, and caller-supplied
small dense matrix-function callbacks. See the
[API and limitations](docs/sparse_complex_funm.md) and
[independent SciPy comparisons](docs/sparse_complex_funm_validation_2026-11-06.md).

## Complex sparse matrix exponential actions (2026-11-06)

The `scifort_sparse_linalg_complex_expm` module computes `exp(t*A)*b` or
`exp(t*A)*B` for complex CSR/CSC matrices and matrix-free operators, with
sampling over a time interval and independently validated Taylor-tail
convergence. See [API, examples and limitations](docs/sparse_complex_expm.md).

## Matrix-free complex eigenvalues and reusable shift-invert LU (2026-11-06)

`scifort_sparse_linalg_complex_eigs` supports complex matrix-free Arnoldi
callbacks, optional matrix-free shift-invert operators, and validated reusable
sparse LU shift caches for CSR/CSC. See
[API, examples, and limitations](docs/sparse_complex_eigs_operator.md).

**Complex nonsymmetric sparse eigenpairs (2026-11-06):**
`scifort_sparse_linalg_complex_eigs` adds `eigs_complex` for complex CSR/CSC,
with six eigenvalue-selection modes, explicit restarted Arnoldi, complex
projected Schur decomposition, and true eigenvector residual checks. See the
[API and limitations](docs/sparse_complex_eigs.md) and [validation
results](docs/sparse_complex_eigs_validation_2026-11-06.md).

**Complex sparse partial SVD (2026-11-06):**
`scifort_sparse_linalg_complex_svds` adds `svds_complex` for canonical
complex CSR/CSC matrices, with largest and smallest singular triplets,
complex Golub--Kahan and one-sided Jacobi projection, and true two-sided
residual checks. See [complex sparse SVD API](docs/sparse_complex_svds.md)
and [validation](docs/sparse_complex_svds_validation_2026-11-06.md).

**Complex sparse LSMR:** `scifort_sparse_linalg_complex_lsmr` provides
`lsmr_complex` for rectangular CSR/CSC systems, regularization, and true
residual diagnostics. See [the API guide](docs/sparse_complex_lsmr.md).

**Complex sparse least squares (2026-11-06):**
`scifort_sparse_linalg_complex_lsqr` adds `lsqr_complex` for rectangular
complex binary64 CSR/CSC matrices, including underdetermined and rank-deficient
systems, damping around a nonzero starting point, and true residual diagnostics.
It uses conjugate-transpose Golub-Kahan bidiagonalization without constructing
normal equations. See [complex LSQR API](docs/sparse_complex_lsqr.md) and
[the randomized SciPy/NumPy comparisons](tools/sparse_reference/compare_sparse_complex_lsqr_scipy.py).

**Complex sparse direct solvers (2026-11-06):**
`scifort_sparse_linalg_complex_direct` adds pivoted sparse LU for complex
CSR/CSC matrices, reusable solves for normal/transpose/adjoint systems,
multiple right-hand sides, and one-shot `spsolve_complex`. The implementation
is self-contained and uses sparse fill-in rather than dense intermediates.
It does not implement SuperLU's column-ordering and supernodal optimizations.
See the [complex sparse direct API](docs/sparse_complex_direct.md) and
[reference comparisons](tools/sparse_reference/compare_sparse_complex_direct_scipy.py).

**Complex N-dimensional COO arrays (2026-11-05):**
`scifort_sparse_coo_nd_complex` and `scifort_sparse_coo_nd_complex_ops`
add complex binary64 sparse tensors, including coordinate coalescing,
conjugation, Hermitian transpose, dense/real conversion, broadcasting,
arithmetic, axis reductions, sparse tensor contraction and conjugated inner
products. One-based coordinates and finite-input checks remain unchanged.
See [complex COO API](docs/sparse_coo_nd_complex.md) and
[randomized SciPy/NumPy tests](tools/sparse_reference/compare_sparse_coo_nd_complex_scipy.py).

**Sparse Boolean broadcasting and reductions (2026-11-05):**
`scifort_sparse_coo_nd_bool_reduce` adds NumPy-compatible logical mask
broadcasting, broadcasted AND/OR/XOR/difference, scalar and axis-wise
`any`/`all`, and exact axis-wise `integer(int64)` counts with a sparse
implicit baseline. Large complement masks require no dense allocation.
See [API](docs/sparse_coo_nd_bool_reduce.md) and
[randomized NumPy tests](tools/sparse_reference/compare_sparse_coo_nd_bool_reduce_numpy.py).

**Complement-aware N-D sparse Boolean masks (2026-11-05):**
`scifort_sparse_coo_nd_bool` adds complement-aware sparse logical sets,
AND/OR/XOR/difference, scalar comparisons including equality to implicit
zeros, true-counting and COO filtering without dense masks. See the
[API](docs/sparse_coo_nd_bool.md) and
[randomized SciPy comparisons](tools/sparse_reference/compare_sparse_coo_nd_bool_scipy.py).

**Sparse N-dimensional predicates and coordinate filters (2026-11-05):**
`scifort_sparse_coo_nd_filter` adds scalar comparison masks,
nonzero counting, sparse COO mask filtering, and repeated-coordinate set
selection with optional complement. No dense mask is needed, even for huge
sparse tensor shapes. See the [predicate/filter API](docs/sparse_coo_nd_filter.md)
and [randomized SciPy comparisons](tools/sparse_reference/compare_sparse_coo_nd_filter_scipy.py).

**N-dimensional COO boolean masks and scatter updates (2026-11-05):**
`scifort_sparse_coo_nd_edit` provides full and axis boolean selection,
transactional coordinate assignment and accumulation, and scalar masked
updates, with stable handling of duplicates and stored zeros. See the
[mask/scatter API](docs/sparse_coo_nd_edit.md) and
[randomized comparisons](tools/sparse_reference/compare_sparse_coo_nd_edit_scipy.py).

**N-dimensional COO indexing and slicing (2026-11-05):**
`scifort_sparse_coo_nd_index` adds one-based scalar, inclusive strided slice,
index-array (orthogonal) and pointwise gathered-coordinate selection for
sparse tensors, preserving duplicate entries and explicitly stored zeros.
See the [indexing API](docs/sparse_coo_nd_index.md) and
[randomized comparisons](tools/sparse_reference/compare_sparse_coo_nd_index_scipy.py).

**N-dimensional sparse reductions and tensor operations (2026-11-05):**
`scifort_sparse_coo_nd_ops` adds sparse axis reductions, singleton-dimension
broadcasting, elementwise addition and multiplication, and sparse tensor
contractions with selected axes. Operations respect duplicate coordinates,
explicit zeros and one-based indices, with configurable entry caps and no
dense intermediate tensors. See [N-D COO operations](docs/sparse_coo_nd_ops.md)
and [randomized comparison tests](tools/sparse_reference/compare_sparse_coo_nd_ops_scipy.py).

**N-dimensional COO sparse arrays (2026-11-05):**
`scifort_sparse_coo_nd` introduces arbitrary-rank, real-binary64 coordinate
sparse arrays with one-based coordinates. It preserves duplicate entries and
stored zeros, supports axis permutation, C/F-order reshaping, canonicalization,
dense conversion and two-dimensional COO interoperability. See the
[N-dimensional COO API](docs/sparse_coo_nd.md) and
[SciPy comparison tests](tools/sparse_reference/compare_sparse_coo_nd_scipy.py).

**Incremental DOK and LIL sparse assembly (2026-11-05):**
`scifort_sparse_dok_lil` now provides mutable, one-based dictionary-of-keys
(hash table) and list-of-lists (sorted row lists) sparse arrays. Both support
single-coordinate assignment, addition and deletion, CSR/CSC/COO and dense
conversions, transpose and matvec. They preserve stored zeros during sparse
conversion but remove entries explicitly assigned zero. See the
[DOK/LIL API](docs/sparse_dok_lil.md) and the independent SciPy comparison
script `tools/sparse_reference/compare_sparse_dok_lil_scipy.py`.

**COO sparse arrays (2026-11-05):**
`scifort_sparse_coo` adds two-dimensional coordinate sparse matrices for
finite-element-style assembly with unsorted, duplicated entries and stored
zeros. Includes canonical CSR/CSC conversion, transpose, matrix-vector
products, `sum_duplicates`, `eliminate_zeros`, and `sparse_find` for COO,
CSR and CSC. See [COO API](docs/sparse_coo.md) and the independent SciPy
comparison script `tools/sparse_reference/compare_sparse_coo_scipy.py`.

**Matrix-function Krylov action (2026-11-05):**
`scifort_sparse_linalg.funm_multiply_krylov` now evaluates `f(tA)*b`
approximately for real canonical CSR matrices or matrix-free operators via
Arnoldi/Lanczos projection. Matrix functions are user-defined dense Fortran
callbacks on the small Krylov projection. The interface provides restart
controls, convergence diagnostics, and explicit failure statuses. See the
[API](docs/sparse_linalg.md#matrix-function-krylov-action-funm_multiply_krylov)
and [validation](docs/sparse_funm_krylov_validation_2026-11-05.md).

`scifort_sparse_laplacian_nd` now provides a matrix-free and CSR
**N-dimensional grid Laplacian** corresponding to
`scipy.sparse.linalg.LaplacianNd`, with Dirichlet/Neumann/periodic boundary
conditions and analytic eigenvalues/eigenvectors.
`scifort_sparse_linalg.spbandwidth` and `is_sptriangular` support both CSR
and CSC structural queries (including stored zeros). See
[grid Laplacian documentation](docs/sparse_laplacian_nd.md) and
[validation](docs/sparse_laplacian_nd_validation_2026-11-05.md).

**Recycled Krylov solvers (2026-11-05):**
`scifort_sparse_linalg` adds **LGMRES** and **GCROT(m,k)** for real
nonsymmetric canonical CSR systems. Both support right Jacobi/ILU(0)
preconditioning and true-residual convergence; users can retain augmentation
vectors across related solves, even if the matrix changes. See the
[API](docs/sparse_linalg.md#augmentedrecycled-krylov-solvers-lgmres-gcrotmk)
and [validation](docs/sparse_recycling_validation_2026-11-05.md).

**Further sparse Krylov solvers (2026-11-05):**
`scifort_sparse_linalg` now has **QMR** and **TFQMR** for real nonsymmetric
CSR systems. QMR evaluates both the matrix and its sparse transpose, with
optional right Jacobi preconditioning. TFQMR uses forward products only and
supports Jacobi or reusable ILU(0) right preconditioning. Both check the
original true residual, accept initial guesses, and report iteration counts.
See the [API](docs/sparse_linalg.md#quasi-minimal-residual-qmr-and-transpose-free-qmr-tfqmr)
and [validation](docs/sparse_qmr_tfqmr_validation_2026-11-05.md).

**Additional SciPy sparse algorithms (2026-11-05):**
SciFort now provides nonsymmetric **BiCG** and **Conjugate Gradient Squared**
(`scifort_sparse_linalg.bicg`, `cgs`), with true-residual stopping criteria,
initial guesses and optional Jacobi preconditioning; CGS also supports an
existing ILU(0) factor. `scifort_sparse.csr_kronsum` directly assembles
canonical CSR Kronecker sums without dense intermediate matrices. See the
[solver API](docs/sparse_linalg.md#biconjugate-gradient-bicg-and-conjugate-gradient-squared-cgs),
[Kronecker-sum API](docs/sparse.md#kronecker-sum-csr_kronsum), and
[validation](docs/sparse_bicg_kronsum_validation_2026-11-05.md).

**New sparse matrix formats and norms (2026-11-05):**
`scifort_sparse` now provides canonical **DIA** and rectangular-block **BSR**
matrices with direct CSR conversion, dense conversion and matrix-vector products.
`scifort_sparse_linalg.norm` covers Frobenius, ±1, ±infinity and approximate
spectral 2-norms on CSR/CSC; `norm_axis` supplies vector norms along a
one-based axis. See [matrix formats](docs/sparse.md#diagonal-dia-and-block-csr-bsr-formats)
and [norm API](docs/sparse_linalg.md#sparse-matrix-and-axis-norms).

**Sparse matrix powers and inverse (2026-11-05):**
`scifort_sparse_linalg` now offers `matrix_power(A,power,info)` for
nonnegative integer powers and `inv(A,info,...)` for explicit sparse
inverses. Both accept canonical one-based real CSR and CSC matrices and
preserve the input storage class. Powers use repeated squaring; inverses
reuse pivoted sparse LU on blocks of identity columns, with a diagonal
fast path. See [API](docs/sparse_linalg.md#sparse-matrix-powers-and-inverse)
and [validation](docs/sparse_matrix_functions_validation_2026-11-05.md).

**Block sparse eigenproblems (2026-11-05):**
`scifort_sparse_linalg.lobpcg` now handles ordinary and generalized real
symmetric sparse eigenproblems, with block initial guesses, constraints,
Jacobi or user-supplied CSR preconditioning, and explicit residual diagnostics.
Its Rayleigh--Ritz solve involves at most `3*k` dense dimensions; the input
matrix stays in CSR. See [API](docs/sparse_linalg.md#locally-optimal-block-preconditioned-conjugate-gradients-lobpcg)
and [validation](docs/sparse_lobpcg_validation_2026-11-05.md).

`scifort_sparse_linalg` now includes **`spilu`**, a dependency-free
natural-order incomplete LU (ILUT) preconditioner with relative threshold
dropping, per-row fill budgets, and optional diagonal shifts. CSR and CSC
inputs are supported; the factors work with existing GMRES and BiCGSTAB
preconditioning through the shared `ilu0_factor` type. See the
[ILUT API and limitations](docs/sparse_linalg.md#threshold-incomplete-lu-spilu--ilut-subset).

**Sparse direct linear solvers (2026-11-05):**
`scifort_sparse_linalg` now includes CSR/CSC `spsolve_triangular`, pivoted
sparse LU (`splu`), factor reuse (`splu_solve`), and general `spsolve` for
vector or matrix right-hand sides. Factorization keeps sparse elimination
fill-in without forming a dense matrix. Tested against SciPy and by
reconstructing `P*A=L*U`. See [API](docs/sparse_linalg.md#sparse-triangular-and-pivoted-direct-solvers-spsolve-splu).
No additional libraries or imported third-party implementations.

**Sparse matrix norm estimation and `expm` (2026-11-04):**
`scifort_sparse_linalg` now includes deterministic `onenormest` for real
canonical CSR matrices, with optional norm-witness vectors, and a full CSR
`expm` computed by applying the existing Taylor exponential action to blocks
of identity vectors. Diagonal matrices use an O(n+nnz) shortcut. These
routines require no additional dependencies; unlike SciPy's `expm`, SciFort
uses exponential action rather than Padé rational approximation, and the
result may fill in. See [API](docs/sparse_linalg.md#matrix-norm-estimation-and-sparse-exponential)
and [validation](docs/sparse_onenorm_expm_validation_2026-11-04.md).

**Sparse matrix-exponential action (2026-11-04):**
`scifort_sparse_linalg` now includes `expm_multiply` and
`expm_multiply_interval` for real CSR matrices acting on vectors or multiple
right-hand sides. A trace-shifted, scaled Taylor series computes the action
without forming `expm(A)`. Time grids are propagated from the point closest
to zero to avoid underflow when an interval crosses into a representable
region. Checked and optimized native tests and SciPy/dense-reference
comparisons are included. See [API and examples](docs/sparse_linalg.md#sparse-matrix-exponential-action-expm_multiply)
and [validation](docs/sparse_expm_validation_2026-11-04.md).

`scifort_sparse_linalg` now also provides a dependency-free **`svds`**
subset for rectangular CSR matrices. It uses Golub--Kahan bidiagonalization
and a projected Jacobi SVD, supports largest singular triplets with a bounded
Krylov subspace and smallest singular values with a full trial basis, and
checks true two-sided singular-vector residuals. See
[partial sparse SVD](docs/sparse_linalg.md#partial-singular-value-decomposition-svds)
and [numerical validation](docs/sparse_svds_validation_2026-11-04.md).

**Sparse iterative solvers (2026-11-04):**
[`scifort_sparse_linalg`](docs/sparse_linalg.md) (including optional thick-restarted `eigs_arnoldi` and symmetric `eigsh_lanczos` with shift-invert and Ritz retention) adds CSR conjugate gradients
(CG), BiCGSTAB, restarted GMRES, Jacobi preconditioning and a rectangular
LSMR adapter reusing SciFort's existing BSD-3-Clause matrix-free kernel.
LSQR for rectangular least squares and an independent ILU(0) factorization
with optional BiCGSTAB/GMRES preconditioning are also included. MINRES
handles symmetric indefinite shifted systems, while an initial real-symmetric
`eigsh` subset finds symmetric Ritz eigenpairs with optional MINRES-based
interior shift-invert (`sigma=`) with LM/LA/SA/SM selection; `eigs`
computes complex eigenpairs of real nonsymmetric sparse matrices, including
reduced-subspace SM via zero-shift inversion and complex `sigma=` via a
restarted complex GMRES inner solve. Both restarted eigenvalue methods
now accept normal-mode shift-invert without changing their defaults.
Nonsymmetric Arnoldi now also provides Schur invariant-subspace locking
via `eigs(method='arnoldi_schur')` or the direct `schur_locking=.true.` flag.
Arnoldi now accepts optional `block_size=1..k` to preserve multiple
independent Krylov search directions and capture repeated or clustered
nonsymmetric eigenvalues, including complex shift-invert. The original
single-vector behavior remains the default. See the
[block Arnoldi tests](docs/sparse_arnoldi_block_validation_2026-11-04.md).
Symmetric Lanczos also offers opt-in hard locking via `eigsh(method='lanczos_locked')`
and multi-seed deflation through `eigsh_lanczos(...,hard_locking=.true.,block_size=...)`.
Both checked and optimized regression suites compare solutions and true
residuals against SciPy. No BLAS/LAPACK dependency was added.

**Sparse random and diagonal-array constructors (2026-11-05):**
[`scifort_sparse_constructors`](docs/sparse.md#sparse-random-matrices-and-diags_array-style-construction-2026-11-05) adds
`csr_random`, `csr_random_array`, `csr_rand`, and `csr_diags_array` with reproducible explicit-RNG
sampling, canonical CSR, scalar diagonal broadcasting and shape inference.
Construction uses O(nnz + rows) memory, not a full dense matrix.

**Sparse structural editing and reshaping (2026-11-05):**
`scifort_sparse` adds `csr_delete_entries`, `csr_delete_rows`,
`csr_delete_columns`, and C-order `csr_reshape`. Sparse-to-sparse
`csr_assign` now clears implicit source zeros without expanding the selected
Cartesian product, supporting very large sparse selections. See
[API details](docs/sparse.md#structural-deletion-sparse-to-sparse-assignment-and-csr-reshape-2026-11-05).

**Sparse general block grids and assignment (2026-11-05):**
`scifort_sparse` adds `csr_bmat`/`csr_block` for two-dimensional grids with
missing blocks, and `csr_set`, `csr_assign_coo`, and overloaded `csr_assign`
for scalar, dense, and CSR assignments.  Updates are transactional on
validation errors and use deterministic last-write-wins semantics for repeated
coordinates.  See `docs/sparse.md` for block-dimension inference and examples.

**Sparse CSR assembly and slicing (2026-11-05):**
`scifort_sparse` now provides `csr_hstack`, `csr_vstack`, `csr_block_diag`,
`csr_kron`, `csr_get`, `csr_select`, and `csr_slice`, with one-based indexing,
canonical sparse storage, and no dense intermediates. New native and SciPy
randomized tests cover rectangular inputs, repeated indices, signed strides,
empty blocks and very wide sparse matrices.

**Rectangular sparse matrices (updated 2026-11-04):**
[`scifort_sparse`](docs/sparse.md) adds one-based canonical CSR storage,
COO/dense conversion, bidirectional CSC conversion, shifted diagonals,
triangular extraction, diagonal replacement, transpose, matvec,
sparse matrix multiplication, scaled addition, entrywise multiplication,
and stored-zero elimination. Standalone regression and randomized SciPy
1.17 tests cover rectangular, empty, and duplicate-coordinate inputs.
No new runtime dependencies were introduced. See also the
[sparse-library license survey](docs/sparse_reuse_repositories.md).

**Sparse graph algorithms (2026-11-03):**
[`scifort_sparse_csgraph`](docs/sparse_csgraph.md) adds one-based CSR graph
constructors; weak/strong connectivity; BFS/DFS traversal and search trees;
Dijkstra, Bellman–Ford, Floyd–Warshall, Johnson, Yen K-shortest
simple paths, and unified shortest-path dispatch; negative-cycle detection; weighted graph Laplacians; predecessor
reconstruction; Kruskal minimum spanning forests, reverse Cuthill–McKee ordering,
Dinic maximum flows, rectangular Hopcroft–Karp matching/structural rank,
minimum-weight full bipartite matching, and sparse CSR Laplacians.
Checked/optimized
native tests and randomized SciPy 1.17 comparisons are included. No new
dependencies or third-party source were added.

**Vector quantization / K-means clustering (2026-11-03):**
[`scifort_cluster_vq`](docs/cluster_vq.md) adds `whiten`, nearest-code `vq`,
fixed-iteration `kmeans2` (supplied centers, random points, K-means++,
or Gaussian initialization), and restart-based `kmeans`. All labels are one-based; stochastic entry points
use SciFort's existing explicit RNG state. Checked/optimized native and SciPy
numerical tests are provided; no third-party code or dependencies were added.

**Hierarchical agglomerative clustering (2026-11-03):**
[`scifort_cluster_hierarchy`](docs/cluster_hierarchy.md) adds a new
`scipy.cluster.hierarchy` subset with seven agglomerative linkage methods
(and seven named wrappers),
observation/condensed-distance inputs, flat distance/maximum-cluster cuts,
cophenetic distances and leaf ordering. It also supports inconsistency
statistics, subtree maxima, monotone flat-cut criteria and full/selected
cut trees, optimal leaf ordering, indexed tree representation, flat
cluster leaders, observation-to-label `fclusterdata`, unified `fcluster`,
cophenetic correlation, one-based MATLAB linkage conversion, validation,
flat-partition equivalence, and numeric dendrogram coordinates. Checked/optimized
native regressions and independent SciPy 1.17 comparisons are included;
no new dependencies.

# SciFort


SciFort is a proposed modern Fortran statistics library with an idiomatic
Fortran API, an FPM build, and a stable C ABI for future Python, R, MATLAB, and
Octave bindings. The code has been written by OpenAI and Anthropic LLMs.

This repository is an independent project. It is not affiliated with or
endorsed by SciPy, NumFOCUS, or the Python Software Foundation.

## SciPy compatibility and SciFort extensions

SciFort aims to follow SciPy's distribution names, mathematical notation,
parameter names and ordering, parameterizations, and mathematical behavior
where both libraries provide the same functionality, while keeping an idiomatic
Fortran interface. Intentional differences and compatibility limitations are
documented in the API reference.

SciFort also welcomes additional distributions, methods, and capabilities beyond
those in SciPy. Functionality absent from SciPy is labeled **SciFort extension**
in the API documentation. Extensions to a shared distribution are labeled at the
method or option level.

## Reuse and attribution

SciFort prefers reusing established Fortran implementations when their numerical
behavior, portability, maintainability, and license meet project requirements.
Adopted code is refactored to SciFort's style and API conventions, with upstream
attribution and applicable copyright and license notices preserved in the source
and documentation. Adaptations are validated using upstream tests where available
and SciFort's numerical tests.

See [AGENTS.md](AGENTS.md#reuse-of-established-fortran-implementations) for the
policy, [CODE_PROVENANCE.md](CODE_PROVENANCE.md) for sources and modifications,
and [THIRD_PARTY_LICENSES.md](THIRD_PARTY_LICENSES.md) and [LICENSES/](LICENSES/)
for license texts.

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

**Arbitrary-rank spatial transform batches (2026-11-03):**
[`scifort_spatial_transform_nd`](docs/spatial_transform_nd.md) supports
runtime-rank quaternion, rotation-vector, and rigid-transform arrays with
shape-aware broadcasting for rotation/SE(3) composition and point application.
Fortran-order flattened outputs preserve arbitrary batch shapes without
requiring language extensions or new dependencies.

**Two-dimensional SE(3) transform grids (2026-11-03):**
[`scifort_spatial_rigid_transform_grid`](docs/spatial_rigid_transform_grid.md)
adds SciPy-style 2D batch broadcasting for rigid transforms, including matrix,
dual-quaternion and twist conversions, point application, transform composition
and weighted reduction over either axis. Standalone SciPy reference tests are
included; the module has no new dependencies.

**Rigid transformations in SE(3) (2026-11-02):**
[`scifort_spatial_rigid_transform`](docs/spatial_rigid_transform.md) adds
`RigidTransform`-style translation, quaternion rotation, homogeneous matrices,
point-cloud transforms, composition, inverses, exponential coordinates, dual
quaternions, fractional screw powers and 1D stacks with weighted means.
Includes independent SciPy comparisons, with no new dependencies.

**SO(3) symmetry groups and reduction (2026-11-02):**
[`scifort_spatial_rotation_group`](docs/spatial_rotation_group.md) implements
SciPy-style `Rotation.create_group` for `Cn`, `Dn`, `T`, `O` and `I`, and
minimum-angle `Rotation.reduce` with optional left/right symmetries for
scalar and stacked rotations. The routines use original, MIT-licensed
Fortran without dependencies and are independently compared with SciPy.

**Smooth 3D attitude splines (2026-11-02):**
[`scifort_spatial_rotation_spline`](docs/spatial_rotation_spline.md) provides
`RotationSpline`-style cubic SO(3) interpolation with continuous angular
velocity and acceleration, from two or more nonuniform orientation keyframes.
Includes scalar/batched evaluation and extrapolation, a standalone 3x3
block-tridiagonal solver, independent SciPy comparison tests, and explicit
rejection of catastrophic near-antipodal nonconvergence; no new dependencies.

**Random spatial rotations and two-dimensional rotation grids (2026-11-02):**
[`scifort_spatial_rotation_grid`](docs/spatial_rotation_grid.md) and
[`scifort_spatial_rotation_stack`](docs/spatial_rotation_stack.md) provide
Haar-uniform `Rotation.random`-style sampling with reproducible caller-owned
streams, plus rank-two rotation arrays, axis-wise broadcasting of rotations
and vectors, rotation-representation conversion and weighted means. No new
runtime dependencies.

**Batched spatial rotations and stacks (2026-11-02):**
[`scifort_spatial_rotation_stack`](docs/spatial_rotation_stack.md) extends the
scalar 3D rotation implementation with 1D collections, vector and stack
broadcasting, quaternion/matrix/Euler/MRP representations, composition,
inverse, powers, indexing, concatenation and weighted means.

**Davenport-angle generalized Euler rotations (2026-11-03):**
[`scifort_spatial_transform`](docs/spatial_transform.md#davenport-angles-generalized-euler-rotations)
provides scalar, rotation-stack, and rank-two-grid generalized Euler angles
with arbitrary perpendicular adjacent axes, intrinsic/extrinsic conventions,
angle extraction, and gimbal-lock reporting.

**3D spatial rotations and orientation Slerp (2026-11-02):**
[`scifort_spatial_transform`](docs/spatial_transform.md) provides native
`scipy.spatial.transform.Rotation`-style conversions among scalar-last
quaternions, rotation vectors, Euler angles (all 24 intrinsic/extrinsic
conventions), modified Rodrigues parameters and proper orthogonal matrices,
plus inversion, composition, 3D vector application, weighted `align_vectors`,
chordal rotation means, and shortest-path pairwise/piecewise `Slerp`.
Uses scaled quaternion normalization, explicit status codes, and independent
SciPy 1.17 numerical comparisons, with no new dependencies.

**Bounded halfspace intersections (2026-11-02):**
[`scifort_spatial`](docs/spatial_halfspace.md) now provides
`halfspace_intersection` for 2D/3D bounded convex regions described by
`A*x+b<=0` and a strictly feasible interior point. The routine finds unique
intersection vertices, returns polar dual points/facets and measurements,
and handles redundant inequalities. It uses original dependency-free Fortran
geometry and is independently compared against SciPy.

**Convex hull construction (2026-11-02):**
[`scifort_spatial`](docs/spatial_convex_hull.md) now includes the native
`scipy.spatial.ConvexHull` subset for planar polygons and 3D polyhedra.
The `convex_hull` factory reports extreme vertices, oriented boundary
simplices, adjacent facets, unit supporting hyperplanes, and 2D/3D
perimeter/surface area and enclosed area/volume. Tested independently
against SciPy 1.17 and implemented without Qhull or new dependencies.

**Spherical Voronoi diagrams (2026-11-02):**
[`scifort_spatial`](docs/spatial_voronoi.md) now supports native
`scipy.spatial.SphericalVoronoi`-style vertex and region construction in 2D
and 3D, with non-origin centers, arbitrary radii, ordered generator regions,
and circular/spherical region measures. An independent incremental 3D convex
hull requires no new dependencies; geometry and areas are compared with
SciPy 1.17 in checked and optimized builds.

**Procrustes spatial alignment (2026-11-02):**
[`scifort_spatial`](docs/spatial_procrustes.md) now exports the native
`scipy.spatial.procrustes` similarity fit and the related real-valued
`scipy.linalg.orthogonal_procrustes` rotation/reflection fit. These use a
self-contained one-sided Jacobi algorithm, handle rank-deficient data,
and include randomized SciPy 1.17 comparisons in checked and optimized
GNU Fortran builds without adding an external dependency.

**Hausdorff distances and spherical interpolation (2026-11-02):**
[`scifort_spatial`](docs/spatial_geometry.md) now exports
`directed_hausdorff` with deterministic one-based witness indices and exact
Euclidean max-min distance, plus the SciFort extension `hausdorff_distance`
for the symmetric value. `geometric_slerp` supports scalar and vector time
parameters, arbitrary dimensions, and extrapolation on the unit sphere.
Antipodal ambiguities are explicitly rejected rather than choosing a random
arc. Both features are dependency-free and checked against SciPy 1.17.

**Expanded `scipy.spatial.distance` coverage (2026-11-02):**
[`scifort_spatial`](docs/spatial_distance.md) supplies `point_distance`,
`pdist`, `cdist`, and `squareform` for 19 SciPy-style distances, including
weighted Euclidean/Minkowski, cosine, correlation, Hamming, Jaccard, Canberra,
Bray--Curtis, standardized Euclidean, Mahalanobis, and Jensen--Shannon. Real-valued samples organized in rows (Fortran arrays of
shape `(observations,features)`) use SciPy's condensed pair ordering, with
explicit status codes rather than Python exceptions. It is self-contained,
without external linear algebra dependencies, and has independent SciPy 1.17
regression tests. The covariance-based metrics accept explicit variances or
inverse covariance, or estimate them automatically in `pdist` and `cdist`.
Five weighted Boolean dissimilarities -- Dice, Rogers--Tanimoto, Russell--Rao,
Sokal--Sneath, and Yule -- now work with all three distance APIs. Finite
nonzero numeric inputs are interpreted as Boolean true, and undefined
zero-denominator cases follow SciPy's bulk conventions.

**New `scipy.spatial` KD-tree pair queries (2026-11-01):**
[`scifort_spatial`](docs/spatial_kdtree.md) now also exports the native
`kd_tree` type and its exact `query_pairs`, `query_ball_tree`,
`count_neighbors`, and `sparse_distance_matrix` methods, with one-based
CSR/COO output, dual-tree pair counting, and SciPy 1.17 comparisons. This
reuses the balanced tree already used by `NearestNDInterpolator` and RBF
queries, rather than introducing another spatial-index dependency. All tree
searches now also support exact Minkowski `p=1..infinity` and optional
periodic `boxsize` (including mixed periodic/nonperiodic dimensions),
without changing existing default Euclidean APIs. `count_neighbors` also
supports multiple radii and noncumulative distance histograms; the new
`weighted_count_neighbors` supports scalar/vector radii and optionally
weighted site contributions on either side.

**New `scipy.interpolate` coverage (2026-10-31):**
[`scifort_interpolate`](docs/interpolate_piecewise_cubic.md) implements PCHIP
shape-preserving interpolation, Akima and modified Akima (`'makima'`), and
CubicHermiteSpline-style interpolation with supplied slopes. Newly added
[`CubicSpline`](docs/interpolate_cubic_spline.md) supports not-a-knot, natural,
clamped, asymmetric and periodic endpoints, including vector-valued ordinate
matrices. Reusable native objects provide evaluation, derivatives, exact
integration, extension controls and cubic coefficients. For additional SciPy
interpolation capabilities, see the [upstream Fortran package and license
review](docs/interpolation_upstream_reuse.md).

**FITPACK-style weighted least-squares B-splines (2026-10-31):**
[`lsq_univariate_spline` and `interpolated_univariate_spline`](docs/interpolate_lsq_fitpack.md)
provide SciPy-style explicit-knot weighted least-squares fits and exact
not-a-knot interpolation (degrees 1–5). A banded QR implementation supports
real multichannel observations without BLAS/LAPACK or a dense design matrix.

**Complete-knot B-spline least-squares fitting (2026-11-01):**
[`make_lsq_spline`](docs/interpolate_make_lsq_spline.md) adds a SciPy-style
weighted QR fit to caller-supplied *complete* knot vectors (degrees 0–15),
with real scalar/multichannel outputs, repeated interior knots, zero weights,
and rank-deficiency diagnostics. It reuses the existing banded QR kernel
without a new Fortran dependency.

**Natural cubic smoothing splines (2026-10-31):**
[`make_smoothing_spline`](docs/interpolate_smoothing.md) provides
weighted roughness-penalized natural cubic smoothing, either at a fixed
regularization parameter or with GCV-selected smoothing. Real scalar and
multichannel arrays are supported, returning SciFort's existing `PPoly`
object for evaluation, derivatives and exact integration. Banded Cholesky
and selected-inverse GCV require linear work and storage without new
BLAS/LAPACK dependencies.

**Scattered-data interpolation (2026-10-31):**
Native [`NearestNDInterpolator`, 2D–5D `LinearNDInterpolator`, and `griddata`
interfaces](docs/interpolate_scattered.md) add arbitrary-dimensional nearest
neighbors and 2D–5D Delaunay/barycentric linear interpolation with real/complex
multichannel data, axis rescaling, and configurable out-of-hull fills. The
local RBF neighbor selection has also been accelerated with a max-heap.
Higher-dimensional linear interpolation is available by supplying a
precomputed, 1-based simplex tessellation via `linear_nd_interpolator(...,
simplices=...)`; **automatic** Delaunay construction now supports 2D–5D
for moderate input sizes using a dimension-generic 4D/5D builder. Cached circumspheres also accelerate native 3D mesh construction. The public `griddata(..., method='cubic')` is supported
in 2D through the Clough–Tocher module.

**Rectilinear grid interpolation (2026-10-31):**
The independent [`RegularGridInterpolator` and `interpn` API](docs/interpolate_regular_grid.md)
provides 1–6D nearest, multilinear, **recursive PCHIP**, and degree-1/3/5
tensor-product B-spline interpolation (`'pchip'`, `'slinear'`, `'cubic'`,
`'quintic'`) over uneven, ascending or descending coordinate axes, including
multichannel data, fills, and extrapolation. Complex data and mixed partial
derivatives (`nu`) are supported by the B-spline methods; PCHIP is real-only
and has no derivative option, matching SciPy. The tensor B-spline fit now
uses an `O(n*k²)` banded solver rather than dense LU. SciFort requires no
new FPM dependencies. `finterp` remains an optional future back-end rather
than being vendored into SciFort.

**SciPy-style B-splines (2026-10-31):**
[`BSpline` and `make_interp_spline`](docs/interpolate_bspline.md) now support
arbitrary nonuniform knot vectors, degrees 0–15 for interpolation, scalar or
multichannel real ordinates, analytic derivatives and integrals, dense basis
matrices, knot insertion, periodic evaluation and not-a-knot/natural/clamped
boundary conditions where applicable, including periodic **construction**
with matched endpoint derivatives for degrees 2–15. See the
[periodic B-spline reference](docs/interpolate_periodic_bspline.md).
An independent SciPy reference suite and
native tests are included. This implementation requires no external runtime
libraries; `bspline-fortran` is reviewed as a possible optional high-performance
backend for tensor-product interpolation but has **not** been copied or
introduced as a dependency.

**General piecewise polynomial bases (2026-10-31):**
The independent [`PPoly`, `BPoly`, and Hermite interpolator API](docs/interpolate_piecewise_basis.md)
adds scalar and multichannel local-power and Bernstein polynomials, analytic
calculus, exact integration, periodic extension, reciprocal basis conversion,
variable-constraint `BPoly.from_derivatives`, and lossless conversion from
existing SciFort cubic spline objects. Checked and optimized native and
SciPy 1.17 reference tests are included; no external dependency is needed.
The piecewise API supports `%extend` for real/complex PPoly/BPoly,
with degree promotion and optional Hermite `orders` arguments. An
[upstream root-finding reuse review](docs/root_finding_upstream_reuse.md)
identifies permissively licensed candidate FPM backends without vendoring.
The piecewise API now also supports real-root finding and equation solving,
conversion from scalar or multichannel BSplines, and complex PPoly/BPoly
polynomial evaluation and calculus.

**Confluent polynomial and rational interpolation (2026-10-31):**
The new [`KroghInterpolator` and `FloaterHormannInterpolator`
interfaces](docs/interpolate_global.md) provide repeated-node Hermite
interpolation, value/derivative evaluation, and pole-free real-axis rational
interpolation for real or complex multichannel ordinates. These implementations
are independent of upstream Fortran libraries and introduce no new dependency.

**Barycentric polynomial interpolation and optional root backends (2026-10-31):**
The independent [`BarycentricInterpolator` and `barycentric_interpolate`
API](docs/interpolate_barycentric.md) adds real/complex multichannel global
polynomial interpolation, derivatives, exact-node evaluation, and
transactional insertion of new interpolation nodes. The `PPoly.roots` and
`PPoly.solve` methods now accept an optional `root_solver=` callback;
[the separately pinned `polyroots-fortran` adapter](docs/interpolate_polyroots_optional_backend.md)
can supply Jenkins--Traub `rpoly` when installed by the user. **No external
root-finding dependency is needed in the default SciFort build.**

**SciPy signal API additions (2026-10-31):**
[`band_stop_obj`](docs/signal_band_stop_objective.md) evaluates continuous
analog band-stop order objectives for Butterworth, Chebyshev and elliptic
families, with stable complementary-elliptic calculations.
[`get_window('exponential',...)`](docs/signal_exponential_center.md) now
supports a displaced center with `parameter2=center`, including the
`'poisson'` alias. Both additions have checked/optimized GNU Fortran tests
and independent SciPy comparisons.

The multichannel signal-analysis module also includes conditional
**group spectral coherence**: it evaluates a determinant summary of linear
association between two groups of channels, optionally conditioned on a third.
The API supports static and moving-window real/complex multitaper spectra,
explicit/automatic/adaptive-CV regularization and frequency-band integration.
This is an original SciFort extension, not a SciPy API. See
[Group spectral coherence](docs/signal_group_coherence.md). The complementary
[canonical spectral coherence modes](docs/signal_canonical_coherence.md)
return the individually sorted squared canonical correlations between channel
groups, rather than only the determinant summary. The optional
`return_vectors=.true.` also returns covariance-normalized complex canonical
coefficients, including warnings for zero and near-tied roots, across the same
static/time-frequency, real/complex and integrated-band pathways.

The optional [canonical mode tracking](docs/signal_canonical_tracking.md)
postprocessor `track_canonical_modes` uses globally optimal matching and
paired complex-phase alignment to follow modes through frequency-wise
strength crossings. Tracking is separate from the original descending-mode
results and reports matching confidence heuristics and segment boundaries.

The complementary [canonical subspace tracker](docs/signal_canonical_subspace_tracking.md) `track_canonical_subspaces` groups nearly equal canonical strengths and follows their joint coefficient subspaces using principal angles, global block matching, and common complex unitary Procrustes alignment. It reports persistent block IDs and handles splits and merges without claiming unique eigenvectors. It is a SciFort extension; rotated block bases must not be interpreted as individual canonical modes. **New:** opt-in `return_covariances=.true.` on canonical estimators and `covariance_metric=.true.` on the subspace tracker compare modes using conditioned spectral-covariance inner products, reducing dependence on channel units; see [covariance-metric tracking](docs/signal_canonical_covmetric_tracking.md).

**New:** [canonical coherence window jackknife](docs/signal_canonical_jackknife.md) provides delete-one-nonoverlapping-window standard errors, bias estimates, and exploratory pointwise Fisher-scale 95% intervals for isolated conditional canonical coherence modes and integrated bands. Both real and complex multitaper signals are supported, including adaptive taper weighting within each window. It does not assume independence of neighboring FFT frequencies, and its intervals are not calibrated hypothesis tests.

**New:** [canonical coherence block jackknife](docs/signal_canonical_block_jackknife.md) extends the window estimator with equal-length contiguous deletion blocks and real/complex `multitaper_canonical_block_jackknife` wrappers. It reports both the total spectral windows and the number of deletion groups, supporting exploratory short-memory-aware variability estimates without treating neighboring frequency bins or adaptive DPSS tapers as independent. The intervals remain heuristic and are **not calibrated**.

**New:** [statistic-specific influence-based block selection](docs/signal_canonical_influence_block_length.md) provides an opt-in `selection_method='influence'` for canonical moving-block bootstrap. It estimates dependence from delete-one-window changes in the requested conditional group coherence or individual canonical mode, with optional frequency/band targeting. The original full-spectral-matrix selector remains the default. Checked/optimized Fortran and independent NumPy regressions and a reproducible VAR(1) coverage study document the behavior and limitations; this is an exploratory heuristic, not a calibrated test.

**New:** [automatic moving-bootstrap block length](docs/signal_canonical_auto_block_length.md) estimates temporal dependence from cross-spectral windows, provides selection diagnostics, and can be invoked by omitting `block_length=` in existing canonical moving-bootstrap calls. It is an exploratory heuristic with simulated coverage checks, not a guarantee of calibrated confidence intervals. Explicit lengths retain their behavior.

**New:** [conditional multimode VAR(p) bootstrap coverage study](docs/signal_varp_conditional_coverage.md) tests five moving-block interval methods against exact finite-taper population targets in stationary five-channel VAR(2)/VAR(3) models, including two canonical modes, an optional conditioning group, and point-frequency/exact band-average spectra. A checked/optimized test runner and 180-row reproducible coverage/bias CSV are included. Near-tied population modes expose unreliable mode-specific inference despite numerical validity flags; these are exploratory diagnostics, not calibrated confidence intervals.

**New:** [exact-band fixed-taper multichannel spectrogram](docs/signal_fixed_taper_spectrogram.md) provides `fixed_taper_spectral_matrix_spectrogram` for real/complex data, arbitrary deterministic taper banks, unequal fixed weights, analytic continuous frequency-band integration, and nonoverlapping or overlapping windows. Its output plugs directly into canonical moving-block bootstrap. The [DPSS multitaper conditional VAR(p) coverage study](docs/signal_varp_multitaper_coverage.md) compares three matched population estimands across three block lengths, five interval methods and four models. Checked/optimized and SciPy references plus an exploratory CSV are included; these are SciFort extensions and do not establish interval calibration.

**New:** [fixed-weight multitaper VAR(p) population targets](docs/signal_varp_multitaper_targets.md) combine any bank of deterministic tapers, including SciFort DPSS tapers, at the complex spectral-matrix level with separate energy normalization, exact frequency-band integrals, and optional conditioning on a third channel group. The four APIs are independently checked against SciPy Lyapunov solutions, time-pair Fourier sums and canonical SVD. Data-dependent adaptive weighting is deliberately excluded. This is a SciFort extension, not a `scipy.signal` API.

**New:** [exact VAR(p) finite-window population spectra and conditional canonical targets](docs/signal_varp_population_targets.md) extend the VAR(1) spectral calculation to stable VAR(p), including analytically integrated (not FFT-grid approximated) cross-spectral bands, complex phase and two-group conditional canonical coherence. The native core has no LAPACK dependency and is independently checked using SciPy/NumPy Lyapunov, time-pair and SVD calculations plus VAR(2) Monte Carlo. This is a SciFort extension rather than a `scipy.signal` function; population plug-in coherence does not imply calibrated confidence intervals.

**New:** [exact finite-window VAR(1) spectral targets and corrected canonical-bootstrap coverage](docs/signal_var1_smoothed_coverage.md) provide `var1_expected_taper_spectrum` as a SciFort extension. Its stationary Lyapunov/taper-autocorrelation formula computes the expected finite-window spectral matrix, including complex cross-channel phases, without external LAPACK dependencies. Independent SciPy/NumPy references and a four-scenario VAR coverage study compare interval methods to the **actual smoothed population estimand**, rather than population coherence at frequency zero. The corrected target changes apparent coverage substantially; 95% calibration is *not* established.

**New:** [block-aware BCa canonical bootstrap intervals](docs/signal_canonical_bca.md)
add `interval_method='bca'` to the precomputed and real/complex multitaper
moving-block canonical bootstrap APIs. Original contiguous delete-block
refits estimate a skewness acceleration for each isolated squared mode and
aggregate group coherence, with explicit NaNs/validity flags for undefined
adjustments. Existing percentile, basic, BC and studentized behavior and
seeded draws are unchanged. Weighted incomplete deletion blocks are
heuristic; pointwise coverage remains uncalibrated under temporal dependence.

**New:** [studentized canonical moving-bootstrap intervals](docs/signal_canonical_studentized.md)
add `interval_method='studentized'` for bootstrap-t confidence limits with
block-jackknife standard errors estimated both from the original spectral
windows and independently inside each bootstrap replicate. It supports
real/complex, conditional, adaptive-multitaper and band-integrated workflows,
keeps all earlier draws and interval defaults unchanged, and flags
unidentifiable denominators or fewer than four deletion blocks rather than
reporting unsupported intervals. These remain exploratory, pointwise
uncertainty estimates, not calibrated hypothesis tests.

**New:** [bias-aware canonical bootstrap intervals](docs/signal_canonical_interval_methods.md)
add opt-in `interval_method='basic'` and `'bc'` to the existing precomputed and
real/complex multitaper moving-block bootstrap, preserving the percentile
default and bootstrap draws. The basic intervals reflect quantiles; BC adjusts
quantiles using a normal-score bias term. Neither basic nor BC is BCa or a
calibrated confidence procedure. Independent SciPy/NumPy numerical references,
public API tests and a small VAR(1) coverage diagnostic document the limitations.

[Canonical coherence moving-block bootstrap](docs/signal_canonical_moving_bootstrap.md) adds reproducible resampling of overlapping possible *contiguous windows*, including optional circular wrap, with deterministic local seeding, bootstrap standard errors/biases, and direct bounded pointwise percentile intervals. Real/complex adaptive multitaper and integrated-band estimators are supported; `canonical_moving_block_bootstrap` also operates on existing window spectral matrices. This is a SciFort extension, not a SciPy function; the intervals are exploratory and not calibrated tests. The previous jackknife APIs are unchanged.

The `scipy.differentiate` layer is exposed by `scifort_differentiate`. It provides
adaptive finite-difference `derivative`, `jacobian`, and `hessian` interfaces with
SciPy 1.17-style nested stencils, error estimates, convergence statuses, iteration
counts, and evaluation counts. Scalar and vector Jacobian callbacks use an
idiomatic Fortran API, and SciFort additionally provides complex-step derivative,
gradient, and Jacobian helpers for analytic functions. The adaptive controller is
adapted from SciPy 1.17.0 under BSD-3-Clause; the scalar complex-step formula is
adapted with attribution from ForDiff 0.1.1 (BSD-3-Clause), while NumDiff 1.7.0 was
consulted as a modern-Fortran API/architecture reference without copying its
implementation code.

MIMO state-space frequency responses also offer reusable Hessenberg and
opt-in complex-Schur plans (`prepare_mimo_frequency`, `solver='schur'`),
including one-based channel selection. See [Signal processing](docs/signal.md).

The rank-1/2/3 `filtfilt` function and rank-1..15 caller-output
`filtfilt_nd` subroutine support Gustafsson initialization with
`method='gust'` and optional `irlen`, in addition to their unchanged
padding default; see [Gustafsson filtering](docs/signal_filtfilt_gust.md).

The shared signal-processing FFT now uses radix-2 transforms for powers of
two, a direct DFT for very small lengths, mixed-radix Cooley-Tukey FFTs for
larger lengths factorizable by 2, 3, 5, and 7, and Bluestein convolution
for all other larger lengths (including prime lengths). The mixed-radix
kernel uses specialized radix-3/5/7 butterflies with pairwise sums and
differences to reduce repeated complex multiplications. This accelerates
spectral estimation and Fourier resampling without requiring FFTW.
Reusable `fft_plan` objects precompute twiddle factors or Bluestein chirps
for repeated same-size transforms
via `forward_fft_planned`, and are used internally by batch spectral estimators,
streaming STFT/inverse STFT, multitaper estimation, and FFT convolution.
For multidimensional convolution, a cost-aware planner can select more compact
mixed-radix FFT shapes when the estimated reduction in workspace offsets the
slower per-point transform; overlap-add keeps its conservative block sizing.
The original `forward_fft` is unchanged; see
[FFT implementation](docs/signal.md#fft-performance-for-arbitrary-lengths).

The `firwin_2d` interface also supports SciFort-only anisotropic
separable designs (independent per-axis cutoffs, sampling frequencies and
passband types) and rotated elliptical radial profiles. See
[Two-dimensional FIR design](docs/signal.md#recent-scipy-signal-algorithm-additions)
for the API, limitations and SciPy numerical references.

The `scipy.signal` compatibility layer is exposed by `scifort_signal`. Real/complex `resample`, `resample_poly`,
`upfirdn`, and `decimate` also support axis-selected rank-two and rank-three
arrays (one-based positive axis indices or negative indices from the end).
For rank-1 through rank-15 arrays, `resample_into`, `resample_poly_into`,
`upfirdn_into`, and `decimate_into` provide assumed-rank interfaces with
caller-allocated outputs, avoiding one wrapper per input rank.
Fourier `resample` also accepts `frequency_window=` (a real FFT-order weight
array), and complex input accepts `complex_window=`. For callable weights,
use `window_function=` (real response) or `complex_window_function=` (complex
input only), evaluated once per invocation at normalized FFT-bin frequencies.
All options are available in the caller-output rank-15 interface as well.
Real-input weights are folded across positive/negative frequencies as in SciPy.
See [frequency-bin resampling](docs/signal_resample_frequency_window.md).
`resample_poly` also supports parameterized symmetric FIR windows (Gaussian,
Kaiser, Tukey, generalized Gaussian, Taylor, general cosine) via optional
`window_parameter=`, `window_parameter2=`, and `window_coefficients=`,
including real/complex vectors and multidimensional `resample_poly_into`;
see [parameterized polyphase windows](docs/signal_resample_poly_parameter_windows.md)
and [arbitrary transfer-function decimation](docs/signal_decimate_transfer_function.md).
`resample_poly` supports complex custom FIR taps and statistical padding;
`upfirdn` additionally supports anti-reflective and anti-symmetric padding.
Discrete-time LTI system objects (`dlti_tf`, `dlti_zpk`, `dlti_ss`,
`dlti_system`) now support representation conversions, frequency analysis,
transient responses, MIMO series/parallel/feedback connections,
Kalman rank reduction, discrete Lyapunov Gramians, and square-root
balanced truncation (a SciFort extension). Continuous-time LTI objects
(`clti_tf`, `clti_zpk`, `clti_ss`, `clti_system`) provide analog simulation,
frequency response, MIMO interconnection, sample-rate conversion, continuous
Lyapunov Gramians, H2 norms, and balanced truncation. See [Discrete-time LTI system objects](docs/signal.md#discrete-time-lti-system-objects)
and [Continuous-time LTI system objects](docs/signal.md#continuous-time-lti-system-objects).
Reusable `czt_plan` and `zoom_fft_plan` derived types cache FFT chirps for repeated transforms; `czt_into`, `zoom_fft_into`, `czt_plan_into`, and `zoom_fft_plan_into` also process arbitrary signal axes of real or complex arrays through rank 15 using caller-provided complex outputs. See [Signal processing](docs/signal.md#reusable-czt-and-zoom-fft-plans).
For irregularly sampled signals, `lombscargle_into` and
`lombscargle_amplitude_into` extend weighted Lomb-Scargle analysis to
real multichannel arrays through rank 15, with a selected sample axis,
per-channel `weights_nd` or shared `weights`, optional `nan_policy='omit'`
for incomplete observations, and atomic caller-output semantics. Shared
weights reuse trigonometric calculations across channels. The reusable
`lombscargle_plan` caches the weighted trigonometric basis for repeated
calls on the same observation-time, frequency, and shared-weight grids;
`lombscargle_plan_into` and `lombscargle_plan_amplitude_into` evaluate
new signals without repeating the setup. The rank-one `lombscargle` also
accepts `nterms=` for weighted multi-harmonic periodic models, with a
`lombscargle_multiharmonic` entry point (QR-based least squares). See
[multiharmonic periodograms](docs/signal_multiharmonic.md). Coefficient recovery and reconstruction are
now available through `lombscargle_harmonic_fit` and
`lombscargle_harmonic_predict`, which return weighted Fourier amplitudes,
frequency-specific offsets, residual RMS errors, and waveform predictions;
see [harmonic fitting and prediction](docs/signal_harmonic_fit.md).
Optional `estimate_covariance=.true.` also returns coefficient covariances,
standard errors, and residual variance; `lombscargle_harmonic_predict_std`
evaluates fitted-mean uncertainty at new times. See
[harmonic uncertainty](docs/signal_harmonic_uncertainty.md).
These are SciFort extensions;
see [multichannel Lomb-Scargle](docs/signal.md#multichannel-lomb-scargle-through-rank-15-scifort-extension).
Recent signal coverage adds Whittaker-Henderson weighted smoothing,
transfer-function and state-space normalization, and FFT/direct method
selection for one-dimensional convolution and correlation. The smoother
supports explicit penalties or automatic REML tuning.
See [Signal processing](docs/signal.md) for interfaces and tests.

The spline APIs include real/complex `cspline1d`, `qspline1d` and
`cspline1d_eval`/`qspline1d_eval`, plus two-dimensional `cspline2d`,
`qspline2d` and `spline_filter`, plus `cspline2d_eval`, `qspline2d_eval`,
`cspline2d_grid` and `qspline2d_grid` for one-based scalar, paired-array and
rectangular-grid coordinates. See
[Two-dimensional splines](docs/signal_spline2d.md) for examples, supported
smoothing modes, and differences from SciPy's complex-valued implementation.

The `scipy.signal` B-spline interface now includes `cspline1d`,
`qspline1d`, `cspline1d_eval`, and `qspline1d_eval`, with real/complex
coefficients and mirrored boundary evaluation. See
[Signal B-splines](docs/signal_spline1d.md) for examples and the cubic
smoothing constraints.

Additional extrema facilities now include SciPy-style `argrelmax`, `argrelmin`,
and `argrelextrema` for one-dimensional real signals. These return one-based
indices and accept the SciPy `order` and `clip`/`wrap` boundary conventions.
The `argrelextrema` comparator is a user-supplied **pure logical function**
accepting two `real(dp)` values. The new `argrelmax_mask_into`,
`argrelmin_mask_into`, and `argrelextrema_mask_into` accept real arrays through
rank 15, operate along a selected one-based (or negative) axis, and return
caller-owned logical masks instead of SciPy's tuples of coordinate arrays.
The mask interface is a **SciFort extension**; see
[Multidimensional relative extrema](docs/signal_extrema_nd.md).
`peak_detect` is a separate SciFort extension,
inspired by the excursion-threshold detector in Jason Christopherson's
MIT-licensed [`peaks`](https://github.com/jchristopherson/peaks) (credited to
Eli Billauer):

```fortran
integer, allocatable :: imx(:), imn(:)
real(dp), allocatable :: mx(:), mn(:)
call peak_detect([0.0_dp, 2.0_dp, 0.0_dp, 3.0_dp, -1.0_dp], &
                 1.0_dp, imx, mx, imn, mn)
```

It returns confirmed alternating maxima and minima after a strictly positive
`delta` excursion, omitting an unconfirmed final swing. The inputs must be
finite; invalid `delta` or nonfinite samples return empty arrays and optional
nonzero `status`. No upstream package is required at build or run time.

The signal layer also provides `find_peaks_cwt` with default Ricker or
caller-supplied real wavelets, `cwt_ricker`, `cwt_wavelet`, and `ricker_wavelet`.
Custom wavelets are Fortran callbacks of `(points,width)` returning an
allocatable one-dimensional array; wavelet peak positions are one-based.
Complex CWT analysis is available through `cwt_complex_wavelet` and
`cwt_morlet2` for real or complex signals, with normalized Morlet wavelets,
caller-defined complex wavelets and SciPy-style conjugate-reversed convolution.
`find_peaks_cwt_complex` detects complex-wavelet ridges using `abs(CWT)`
maxima and percentile magnitude noise, while `find_peaks_cwt_magnitude`
accepts precomputed coefficients. The older real Ricker peak detector retains
its signed-maxima convention. See `docs/signal.md`.
`place_poles` now supports state-feedback pole assignment for controllable
single-input, fully actuated, and partially actuated MIMO systems, including
real targets, complex-conjugate pairs, and feasible repeated poles. General
MIMO uses a pivoted-nullspace eigenspace solver plus optional
`method='YT'` real-pair/complex-pair and `method='KNV0'` real rank-one
transfer-matrix updates. Optional `rtol` and `maxiter` control determinant
optimization; the result reports `nb_iter`, `rtol` and `det_transfer`. It
follows the Yang--Tits pair-update order but does not reproduce SciPy's
initialization or every degeneracy-handling detail. Numerically difficult or
defective assignments can still return a nonzero status. See
[Signal processing](docs/signal.md) for supported cases, validation and
attribution.

It provides one-dimensional peak detection (`find_peaks`, `peak_prominences`,
`peak_widths`), including per-sample height/threshold/prominence/width/plateau
bounds and `*_open` flags for property-only measurements, real/complex rank-1–3 `detrend` with selected axes, and direct
real/complex `convolve`, `correlate`, and `correlation_lags`.
Rank-independent `convolve_nd` and `correlate_nd` add direct, FFT and
automatic methods through rank 15, with selected axes, broadcasting and
caller-supplied outputs. Rank-two `convolve2d` and `correlate2d` support
real, complex and mixed inputs, with fill, periodic and symmetric boundaries.
The APIs use one-based sample positions; SciPy-compatible
contour/width semantics and intentional limitations are documented in
[Signal processing](docs/signal.md). Digital filter design and response
routines (`firwin`, `firwin2`, `firls`, `remez`, `minimum_phase`, `butter`, `cheby1`, `cheby2`, `ellip`, `bessel`, `iirfilter`, `iirdesign`,
`buttord`, `cheb1ord`, `cheb2ord`, `ellipord`, `iirnotch`, `iirpeak`,
`iircomb`, `gammatone`, `freqz`,
`sosfreqz`, `sos2tf`) are also available. Polynomial/roots/sections
conversions (`tf2sos`, `tf2zpk`, `zpk2tf`, and digital `zpk2sos`/`sos2zpk`) are available with
complex-capable result types and status codes. Weighted `firls` design now
uses stable narrow-band integration and a dependency-free rank-revealing
fallback for ill-conditioned Gram matrices. Butterworth, Chebyshev I/II,
elliptic, and Bessel/Thomson filter design support `analog=.true.` with rad/s cutoff frequencies,
low/high/bandpass/bandstop analog BA coefficients, and evaluation using
`freqs` or `freqs_zpk`. Minimum-order `buttord`, `cheb1ord`, `cheb2ord`,
`ellipord`, and `iirdesign` now support both analog rad/s specifications and two-edge
bandpass/bandstop designs; order results include a `cutoffs(:)` array.
Analog `iirfilter_zpk` and `iirfilter_sos` provide zero-pole-gain and
s-plane SOS output for all four classical filter families; `zpk2sos_analog`
and `sosfreqs` convert and evaluate analog sections without modifying the
existing BA-returning filter interfaces. The analog lowpass frequency
transformations `lp2lp`, `lp2hp`, `lp2bp`, `lp2bs` and their `_zpk`
variants are available for independently transforming existing analog filters.
Normalized elliptic prototypes (`ellipap`) and digital bilinear transforms
(`bilinear`, `bilinear_zpk`) are also exported from `scifort_signal`;
Bessel prototypes now support orders through 20; `bessel` and `iirfilter(ftype='bessel')` expose analog/digital low/high/band transforms. The analog `iirfilter_zpk` and `iirfilter_sos` selectors also accept Bessel filters. Analog SOS coefficients are in
descending powers of `s` and cannot be used with digital `sosfilt`.
New `deconvolve` and `envelope` APIs cover polynomial impulse-response
removal and FFT-based amplitude envelopes for real/complex signals.
`deconvolve` returns a quotient and full-length remainder;
`envelope` supports ranks 1-3, selected one-based axes, signed frequency-bin
bands, optional Fourier resampling, squared magnitudes and lowpass/all/no
residual selection. The new caller-output `envelope_into` extends the same
analysis to real or complex arrays of rank 1-15, with negative axis indexing,
noncontiguous outputs, and transactional residual/envelope writes. See [Signal processing](docs/signal.md).
The new independent `order_filter` supports masked, zero-padded rank selection on
real and integer arrays of rank 1-3, with logical, integer, and real footprint
masks. `sepfir2d` adds mirror-boundary two-dimensional separable FIR filtering
for real, complex, and integer images, including promoted complex taps; see
[Signal processing](docs/signal.md).
The independently implemented `medfilt`, `medfilt2d`, and `wiener`
functions add zero-padded median and Wiener denoising for real or integer
arrays of ranks 1-3, plus complex-valued Wiener filtering with real/complex
noise parameters. Rectangular windows, supplied or estimated noise
power, and status validation are supported; see [Signal processing](docs/signal.md).
The signal module also provides
Savitzky-Golay smoothing/derivatives and real/complex direct-form FIR/IIR
filtering (`savgol_coeffs`, `savgol_filter`, `savgol_filter_nd`, `savgol_irregular`,
`lfilter`, `lfilter_nd`, `lfilter_zi`).
The new SciFort-specific `savgol_irregular` performs weighted local-polynomial
smoothing and derivatives at arbitrary times on nonuniform sample grids, including
complex observations; see [Signal processing](docs/signal.md).
Savitzky-Golay filtering supports rank-1–4 array-returning functions and a
new caller-output, assumed-rank subroutine `savgol_filter_nd` for ranks 1–15
(real or complex, with one-based axis selection). The older functions now
reuse the same rank-independent engine. Center and interpolation-edge
weights are reused across channels for efficient multidimensional filtering. `lfilter` also promotes mixed real/complex
coefficient and signal arrays to complex results and streaming states.
The scalar/elemental `gauss_spline`, `kaiser_beta`, and `kaiser_atten`
functions and the `kaiserord` FIR design helper now cover additional SciPy
algorithms. `kaiserord` returns a typed tap-count/beta/status result, and
`firwin` accepts `window_parameter=beta` for parameterized Kaiser designs.
It also accepts SciPy-style `width=` to select a Kaiser taper from a physical
transition-band width. Scalar or vector cutoffs, logical or string
`pass_zero`, and odd-edge multiband filters are now supported. `firwin2`
now supports parameterized Gaussian, Kaiser,
Tukey and other symmetric windows, with optional explicit window arrays.
See [FIR-width and window options](docs/signal_fir_width_windows.md) and
[Gaussian splines and Kaiser design](docs/signal_kaiser_spline.md).
`get_window` now covers twenty-five named window families, including
Chebyshev, Taylor, Kaiser-Bessel-derived, coefficient-vector `general_cosine`,
and DPSS (`parameter=NW`) windows. Optional Kaiser, Gaussian, Tukey, and
other shape parameters and SciPy-compatible symmetric/periodic choices
are supported. For usage and the list of remaining windows, see
[`docs/signal_window_families.md`](docs/signal_window_families.md).
The spectral-analysis layer also provides `get_window`, `periodogram`, `welch`,
`csd`, `coherence`, `stft`, `istft`, `check_cola`, `check_nola`,
`spectrogram`, and Slepian/DPSS tapers (`dpss`), including SciPy's
`norm='subsample'` half-sample normalization. DPSS also has optional
`solver='auto'|'ql'|'partial'`: the partial tridiagonal eigensolver reduces
workspace to O(n*k) for a few tapers on long records. The optional
`ratio_method='auto'|'direct'|'fft'` controls exact concentration ratios:
for lengths of at least 256 the default uses O(k*n*log(n)) FFT-based
correlations rather than the original O(k*n*n) direct method; see
[large-window DPSS methods](docs/signal_dpss_partial.md). The new
`periodogram_into`, `welch_into`, `csd_into`, and `coherence_into` expose
caller-allocated rank-1–15 spectral results with selectable axes and real,
complex, and mixed-type input; see
[rank-independent spectral estimators](docs/signal_spectral_nd.md).
`stft_into` and `spectrogram_into` additionally handle rank-1–14 arrays
with a selectable signal axis and rank-2–15 time-frequency outputs.
`istft_into` reconstructs real/complex rank-1–14 signals from rank-2–15
STFT arrays, with independently selectable frequency and time axes; see
[rank-independent STFT and inverse STFT](docs/signal_timefreq_nd.md).
Legacy spectral estimators and the multidimensional wrappers now also accept
[parameterized or caller-supplied window arrays](docs/signal_windowed_estimators.md),
with real/complex support and SciPy-compatible validation rules.
The [streaming forward STFT and spectrogram](docs/signal_stft_stream.md)
analyze long or live input in arbitrary chunks using a bounded rolling window,
and the [streaming inverse STFT](docs/signal_istft_stream.md) reconstructs
frames incrementally with bounded overlap/FFT storage. The
[multichannel streaming API](docs/signal_stream_multi.md) adds synchronized
STFT, power spectrogram, and inverse STFT operations on `(sample,channel)`
chunks with `(frequency,channel,frame)` forward outputs, without storing a
complete recording.
The SciFort-specific `multitaper_partial_coherence` and
`partial_coherence_matrix` functions estimate conditional coherence between
channels from spectral precision matrices, with optional relative ridge
loading and explicit singular-frequency diagnostics. See
[Partial coherence](docs/signal_partial_coherence.md). Moving-window
`multitaper_tf_partial_coherence` and `partial_coherence_tf_matrix`
further provide **time-frequency partial coherence** and per-cell singularity
masks; see [Time-frequency partial coherence](docs/signal_tf_partial_coherence.md).
Both interfaces now optionally expose complex partial coherency (including
conditional phase) and diagonal-target spectral shrinkage via
`return_coherency=.true.` and `shrinkage=`. Opt-in `auto_shrinkage=.true.`
also selects frequency-local diagonal shrinkage using a Gaussian-Wishart
moment approximation, with coefficients returned in `shrinkage_used`.
Automatic shrinkage supports equal and eigenvalue-weighted DPSS tapers,
with their effective count available as `effective_tapers`; adaptive weights
are excluded **from this Gaussian-Wishart mode**. A separate opt-in
`cv_shrinkage=.true.` method refits weights leaving out each taper and selects
predictive diagonal shrinkage from held-out spectra; it also supports
**adaptive** weights. See [automatic shrinkage](docs/signal_auto_shrinkage.md)
and [cross-validated shrinkage](docs/signal_cv_shrinkage.md) for distinct
assumptions and limitations.
These are SciFort extensions.

**Multiple coherence** (`multitaper_multiple_coherence`,
`multitaper_tf_multiple_coherence`, `multiple_coherence_matrix`,
`multiple_coherence_tf_matrix`) measures how well each channel is jointly
explained by all other channels, returning conditional and explained spectral
power. It supports the existing ridge, manual, automatic and cross-validated
shrinkage modes for static and time-frequency estimates. See
[Multiple coherence](docs/signal_multiple_coherence.md). This is a SciFort
extension, not a `scipy.signal` API.
**Band-integrated spectral coherence** (`integrate_spectral_bands`,
`integrate_tf_spectral_bands`, `band_multiple_coherence_matrix`,
`band_partial_coherence_matrix`, their time-frequency variants and direct
real/complex multitaper wrappers) estimates multivariate coherence from
cross-spectral matrices **integrated over user-selected frequency intervals**.
This is different from averaging coherence values across FFT bins. It supports
nonuniform and wrapped two-sided grids, exact off-grid band endpoints and
manual post-integration regularization. See
[Band-integrated coherence](docs/signal_band_coherence.md). This is a SciFort
extension, not a `scipy.signal` function.


`multitaper_psd` combines these tapers into an equal-weight,
concentration-weighted, or adaptive Thomson power spectrum for real/complex 1-D signals or
rank-2 multichannel arrays. `multitaper_spectral_matrix`, `multitaper_csd`,
and `multitaper_coherence` add cross-spectra and interchannel analysis;
`multitaper_psd_ci` supplies pointwise log-jackknife confidence intervals;
real inputs use one-sided spectra and complex inputs use two-sided spectra. These additions use a self-contained
tridiagonal eigensolver and require no additional Fortran packages.

`multitaper_stream` incrementally accepts real or complex chunks and emits
multichannel multitaper cross-spectrogram frames without retaining the full
signal. It caches DPSS tapers, supports arbitrary chunk boundaries and an
on-demand running mean spectral matrix. See `docs/multitaper_streaming.md`
and `test/compare_signal_multitaper_stream.py` for examples and numerical checks.
The underlying FFT is now **planned once per multitaper operation** and
reused across tapers, channels, and windows. `multitaper_stream` retains the
plan through `reset()`, replacing it on `initialize()` when `nfft` changes.
The internal N-dimensional FFT convolution and tiled overlap-add kernels
similarly reuse a separate plan per selected axis (ranks 1-15), without new
API parameters. See [FFT plan reuse](docs/signal.md),
`test/test_signal_fft_plan_reuse_consumers.f90`, and the local before/after
benchmark in `benchmarks/bench_signal_multitaper_fft_plan_reuse.f90`.

Streaming multitaper estimates now support independent **linear detrending of
each completed window**; whole-record multitaper estimators support linear
detrending of the complete signal. See the SciPy comparison harness
`tools/signal_reference/compare_signal_multitaper_whole_stream_detrend.py`.

The related SciFort `multitaper_harmonic_scan` extension adds
Bonferroni/Holm/BH/BY-adjusted p-values and line significance indicators,
with configurable families of frequency bins. `multitaper_noise_floor` is
another SciFort extension, returning the local-median colored-noise background
and normalized PSD for real and complex signals. The matching
`multitaper_noise_floor_spectrogram` and streaming
`multitaper_stream%push_noise_floor`/`mean_noise_floor` APIs provide frequency-local
backgrounds for each window and channel, or for running mean PSDs, without
retaining the full record. The ratios are descriptive, not calibrated
line-significance probabilities. See
[Multitaper time-frequency analysis](docs/multitaper_time_frequency.md).

Two additional SciFort spectral extensions are `multitaper_spectrogram`
and `multitaper_harmonic_ftest`. The spectrogram uses overlapping DPSS
windows and **caches the tapers** across frames. It supports real one-sided
or complex two-sided PSDs, adaptive/equal/eigen weighting, and rank-2
multichannel input with optional one-based `axis=1`/`axis=2` selection.
Windowed multitaper time-frequency estimators additionally accept
`detrend_mode='none'`, `'constant'` or `'linear'`, with independent
least-squares linear detrending of every segment and channel. This new
optional trailing parameter overrides the legacy logical `detrend=`;
`multitaper_stream` and the whole-record multitaper PSD, spectral matrix,
CSD, coherence, jackknife CI, harmonic test and harmonic scan interfaces now
accept the same trailing `detrend_mode` option while retaining their
existing logical `detrend` argument and defaults. The
output has separate frequency and segment-center time axes, with power
indexed `(frequency, time)` for 1-D signals or `(frequency, time, channel)`
for 2-D input. Segments are unpadded (no partial end frames), and the default
overlap is half the segment length.

Multitaper cross-spectrograms, time-frequency coherence, and full
multichannel spectral-matrix spectrograms are available through
`multitaper_cross_spectrogram`, `multitaper_tf_coherence`, and
`multitaper_spectral_matrix_spectrogram`. They preserve the spectrogram's
frequency/time conventions and support equal/eigen/adaptive weights,
including independent per-channel adaptive weights for cross-spectra.
See `docs/multitaper_time_frequency.md` and
`test/compare_signal_multitaper_tf_cross.py` for the API and 120 randomized
SciPy/NumPy reference comparisons.

The harmonic test estimates a complex sinusoidal coefficient at each FFT
frequency by regression across the DPSS eigencoefficients, then returns its
F-statistic and pointwise p-value. It uses the approximate null distribution
`F(2,2*K-2)` for interior frequency bins; DC/Nyquist bins and any frequency
scan need particular care. P-values are **not** multiple-testing adjusted.
At an isolated positive interior bin of a real sinusoid, its fitted complex
amplitude is approximately half the physical sinusoidal amplitude. Both
algorithms are independent SciFort implementations, not direct SciPy APIs.

```fortran
use scifort_signal, only : dp, multitaper_spectrogram, &
    multitaper_harmonic_ftest, multitaper_spectrogram_result, &
    multitaper_harmonic_result
real(dp) :: x(256)
type(multitaper_spectrogram_result) :: tf
type(multitaper_harmonic_result) :: lines
! ... fill x with sampled data ...
tf = multitaper_spectrogram(x,3.5_dp,nperseg=64,noverlap=32, &
    kmax=5,fs=100.0_dp,nfft=128,weighting='adaptive')
lines = multitaper_harmonic_ftest(x,3.5_dp,kmax=5,fs=100.0_dp)
```

The Fortran regressions are `test/test_signal_multitaper_time_frequency*.f90`.
Run `python test/compare_signal_multitaper_time_frequency.py` for independent
SciPy/NumPy comparisons (90 spectrogram and 100 harmonic F-test cases).

The original periodogram and Welch routines follow the spectral-estimation
organization in Jason Christopherson's MIT-licensed `spectrum` project with
attribution. All routines compile without new dependencies; the shared FFT now chooses
radix-2, small direct DFT, mixed-radix, or Bluestein transforms according to
length. Real and
complex 1D spectral inputs are supported, including two-sided complex
spectra; `fftconvolve` supports rank-1/2/3 arrays, including selected
one-based convolution axes and singleton broadcasting, while the new
caller-output `fftconvolve_nd(in1,in2,y,mode,axes,status)` supports ranks
1-15, all real/complex input combinations, and noncontiguous sections.
`oaconvolve` continues to provide blocked 1D FFT convolution;
`oaconvolve_nd` extends this to ranks 1-15 with reusable kernel FFTs and
multi-axis overlap-add blocks. `lfilter` supports real and complex rank-2/3
filtering along an optional one-based axis, including initial and final filter states. For ranks 1-15,
`call lfilter_nd(b,a,x,y,zi=...,zf=...,axis=...)` accepts caller-supplied
shape-matched output and states; mixed real/complex coefficients and inputs
are supported and promote the output and states to complex.
The new `short_time_fft` derived type supplies centered, invertible
short-time Fourier transforms with canonical/tight synthesis windows,
real or complex rank-one through rank-fourteen data, four FFT layouts,
selectable signal axes, detrending, and (cross-)spectrograms.
It includes `initialize_window`, `initialize_from_dual`,
`initialize_equal_dual`, `stft`, `stft_detrend`, and `istft` methods,
and helpers for time/frequency coordinates, slice ranges and padding boundaries,
with options and SciPy comparisons documented in [Signal processing](docs/signal.md).
The public `freqz_sos` alias matches the newer SciPy naming for `sosfreqz`.
The new `czt` and `zoom_fft` routines perform complex chirp Z-transforms and
frequency-selective transforms for real or complex 1D signals. They use
Bluestein FFT convolution on unit-circle contours and Horner evaluation for
nonunit contours; no additional dependency is needed. See [Signal processing](docs/signal.md).
Cascaded second-order and forward/backward digital filtering are available
through `sosfilt`, `sosfilt_nd`, `sosfilt_zi`, `filtfilt`, and `sosfiltfilt`.
The three array-returning filtering functions accept real/complex rank-1/2/3
arrays with optional one-based axes; SOS states stream independently per
channel. `call sosfilt_nd(sos,x,y,zi=...,zf=...,axis=...)` additionally
supports assumed-rank 1-15 filtering with real/complex mixed types, using
leading-section state arrays through rank 14 or packed rank-one states for
rank 15. `sosfilt_zi` retains its rank-2 steady-state result. New caller-output
`filtfilt_nd(b,a,x,y,...)` and `sosfiltfilt_nd(sos,x,y,...)` support
forward/backward zero-phase filtering of rank-1 through rank-15 arrays,
all mixed real/complex coefficient/signal combinations, one-based axes,
and existing `'odd'`, `'even'`, `'constant'`, and `'none'` padding modes.
`filtfilt_nd` and `sosfiltfilt_nd` additionally support `method='gust'` and
`irlen` without padding, including signals shorter than the default pad
length. The rank-1/2/3 array-returning `filtfilt` and `sosfiltfilt`
interfaces also offer this method, preserving the original default and
positional calls. The new `sosfiltfilt_gust` kernel operates directly on
second-order-section states without forming high-order transfer polynomials;
it returns a numerical-failure status if initial-state estimation is ill-conditioned.
SciPy's own `sosfiltfilt` API only offers padding; the SOS Gustafsson option
is a SciFort extension. See [SOS Gustafsson filtering](docs/signal_filtfilt_gust.md).
These signal routines remain independent SciFort implementations.
Polyphase multirate filtering (`upfirdn`, `resample_poly`), FIR/IIR
`decimate`, and periodic Fourier `resample` are available for real/complex
one-dimensional signals with no additional FPM dependencies. FFT-based analytic
signals (`hilbert` for real 1D data, `hilbert2` for real 2D data) are also
implemented; each supports optional transform lengths and returns complex(dp)
results. For arrays through rank 15, the new caller-output `hilbert_into` and
`hilbert2_into` accept selected axes, optional FFT lengths, noncontiguous
sections, and preserve outputs on failure; see
[rank-independent Hilbert transforms](docs/signal_hilbert_nd.md). Frequency sweep generation (`chirp`, `chirp_complex`, `sweep_poly`)
and Gaussian-modulated pulses (`gausspulse`, `gausspulse_cutoff`) are also
available with tests against SciPy 1.17. Periodic `sawtooth`/`square`
waveforms, real 1D/2D `unit_impulse` arrays, and `max_len_seq` binary
shift-register sequences (2-32 bits, resumable state, and `int64`
forward/rewind skip-ahead) are supported. Irregularly sampled
weighted spectra are available through `lombscargle` (power, normalized,
and complex-amplitude variants). For multiple irregularly sampled signals,
`lombscargle_harmonic_prepare` and `lombscargle_harmonic_fit_batch` reuse
weighted multiharmonic QR factors and provide per-channel coefficients and
optional uncertainty estimates; see
[Prepared harmonic fitting](docs/signal_harmonic_plan.md).
`vectorstrength` computes event-phase synchrony for one or more periods.
See the signal documentation for
supported sweep methods, waveform parameters, boundary modes, filter-family
limitations, and Hilbert transform conventions.

The `scipy.optimize` layer additionally provides `linear_sum_assignment` for
rectangular minimum-/maximum-cost bipartite matching, including forbidden
edges, and `isotonic_regression` for weighted, increasing/decreasing monotone
least squares. Both are independent SciFort implementations with SciPy-compatible
numerical semantics and one-based indices. Assignment reports status rather
than raising an exception; isotonic regression returns a result with fitted
values, pooled weights, and one-based block boundaries. See `docs/api.md`.

`scipy.optimize.quadratic_assignment` now offers original, deterministic
Fortran implementations of **FAQ** (Frank-Wolfe relaxation) and **2-opt**
(local permutation swaps). Both support minimization/maximization and fixed
seeded matches; FAQ accepts a doubly stochastic `p0`, and 2-opt accepts
`partial_guess` as a warm start. Output permutations and input match pairs
use **one-based indices**. These methods are approximate and do not certify
global optimality. See `docs/optimize_quadratic_assignment.md` and
`example/optimize_quadratic_assignment_example.f90`.

`scipy.optimize.milp` now has a **dense mixed-integer branch-and-bound**
implementation reusing SciFort's existing two-phase simplex LP solver. It
handles continuous, integer, semi-continuous, and semi-integer variables,
two-sided linear constraints, explicit node limits, and a reported primal-dual
gap. This is a small-problem baseline, **not a HiGHS replacement**; it lacks
sparse matrices, cutting planes, and advanced presolve. See
`docs/optimize_milp.md` and `example/optimize_milp_example.f90`.

The initial `scipy.optimize` layer is exposed by `scifort_optimize`. Scalar root
finding provides `root_scalar`, `bisect`, `brentq`, `brenth`, `ridder`,
`toms748`, and `newton`, including secant, derivative Newton, and Halley modes.
`root_result` reports the approximate root, convergence flag/status, method,
iteration count, and function-call count. Bisection, Brent, and Ridders kernels
are adapted with attribution from Jacob Williams' roots-fortran 1.0.0
(BSD-3-Clause); SciPy 1.17.0 supplies the compatible dispatch/result semantics,
Newton/secant/Halley behavior, and TOMS748 organization under BSD-3-Clause.
SciFort returns explicit status values for invalid brackets and nonfinite
objective values instead of raising Python exceptions. Scalar minimization adds
`minimize_scalar`, `brent`, `bounded`, and `golden`, with SciPy 1.17-compatible
default tolerances, bracketing, convergence/evaluation accounting, and automatic
selection of bounded minimization when bounds are supplied. The bounded Brent
kernel is adapted with attribution from Jacob Williams' fmin 1.1.2
(BSD-3-Clause); the unbounded Brent, golden-section, bracketing, and dispatch
semantics are adapted from SciPy 1.17.0 under BSD-3-Clause. The public
`bracket` routine also exposes the downhill bracketing search already shared
by Brent and golden-section minimization. For example:

```fortran
use scifort_optimize, only : bracket, bracket_result, dp
type(bracket_result) :: br
br = bracket(objective, xa=0.0_dp, xb=1.0_dp, grow_limit=20.0_dp)
if (br%success) print *, br%xa, br%xb, br%xc, br%function_calls
```

`bracket_result` contains the three abscissae, three function values, work
count, expansion iterations, success flag and status. As in SciPy, the output
triple can be descending when the initial downhill direction is to the left;
it is not automatically sorted. An invalid or unbracketable function returns
a failure result rather than raising an exception. The default growth limit
is 110 and the default iteration limit is 1000. The tests in
`test/test_optimize_bracket.f90` are generated from
`tools/generate_optimize_bracket_reference.py` using SciPy 1.17.

The scalar `bracket_root` routine adapts
`scipy.optimize.elementwise.bracket_root` (SciPy 1.17) to a real scalar
function callback. Unlike `root_scalar`, it *finds* a sign-changing interval
by expanding from the initial endpoints. For example:

```fortran
use scifort_optimize, only : bracket_root, bracket_root_result, dp
type(bracket_root_result) :: br
br = bracket_root(f, xl0=0.0_dp, xmin=-10.0_dp, xmax=10.0_dp)
if (br%success) print *, br%xl, br%xr, br%nfev
```

If omitted, `xr0` defaults to `xl0 + 1`, reduced to
`xl0 + min((xmax-xl0)/8, 1)` when a finite `xmax` is supplied. The
`factor` default is 2 and `maxiter` default is 1000. The result contains
`xl`, `xr`, `fl`, `fr`, `nit`, `nfev`, `status`, and `success`.
A status of 0 means a bracket was found; -1 means an expansion limit was
reached, -2 means the iteration cap was reached, -3 indicates nonfinite
values, and -5 indicates invalid input. The input is checked before
calling `f`; unlike Python, invalid settings return a result instead of
raising an exception. In addition to the scalar interface, rank-one input
and rank-two input grids are supported as described below.
SciPy reference points, statuses and work counts are generated with
`tools/generate_optimize_bracket_root_reference.py`.

`find_root` refines a valid two-point, sign-changing root bracket using
Chandrupatla's interpolation/bisection algorithm from SciPy 1.17
`elementwise.find_root`. A bracket returned by `bracket_root` can be passed
directly to it:

```fortran
use scifort_optimize, only : bracket_root, bracket_root_result, &
    find_root, find_root_result, dp
type(bracket_root_result) :: br
type(find_root_result) :: sol
br = bracket_root(f, xl0=0.0_dp)
if (br%success) then
    sol = find_root(f, [br%xl, br%xr], xatol=1.0e-12_dp)
    if (sol%success) print *, sol%x, sol%f_x, sol%nfev
end if
```

The optional `xatol`, `xrtol`, `fatol`, and `frtol` control convergence;
`maxiter` limits refinement iterations (default 2046 for double precision).
Defaults match SciPy double precision: `4*tiny`, `4*epsilon`, `tiny`,
and zero, respectively. A `callback=` accepts a `find_root_result` and
an output logical `stop`; it runs initially and after each iteration.
`find_root_result` reports the root estimate, sorted bracket and function
values, `nit`, `nfev`, `status`, and `success`. Statuses are 0 for convergence,
-1 for an invalid sign-changing bracket, -2 for an iteration limit, -3 for a
nonfinite function value, -4 for callback termination, and -5 for malformed
inputs (SciFort convention). Rank-one independent problem arrays are supported;
fully general SciPy/Numpy broadcasting and keyword/extra-argument
vectorization are not yet provided. References generated by
`tools/generate_optimize_find_root_reference.py` are checked in
`test/test_optimize_find_root.f90`.

`bracket_minimum` provides a scalar adaptation of
`scipy.optimize.elementwise.bracket_minimum` for finding a three-point
bracket around a local minimum before calling `brent`. For example:

```fortran
use scifort_optimize, only : dp, bracket_minimum, bracket_minimum_result
type(bracket_minimum_result) :: br
br = bracket_minimum(objective, xm0=0.0_dp, xmin=-10.0_dp, xmax=10.0_dp)
if (br%success) print *, br%xl, br%xm, br%xr, br%nfev
```

The optional arguments are `xl0`, `xr0`, `xmin`, `xmax`, `factor`, and
`maxiter`. The defaults are half-unit initial offsets (reduced to 1/16
of the distance to a finite bound), geometric growth factor 2, and
1000 iterations. The ordered triple satisfies `fl >= fm <= fr` with at
least one strict inequality. Status values are 0 for success, -1 for
reaching the limit, -2 for the iteration cap, -3 for nonfinite values,
and -5 for invalid input. As with `bracket_root`, malformed inputs return
before evaluating the objective rather than raising Python exceptions.
Rank-one problem arrays and scalar/elementwise argument broadcasting are
supported, but not higher-rank NumPy broadcasting. SciPy-generated scalar
regressions are in `test/test_optimize_bracket_minimum.f90`.

`find_minimum` refines a valid three-point minimum bracket using the scalar
Chandrupatla minimizer from SciPy 1.17 `elementwise.find_minimum`:

```fortran
use scifort_optimize, only : dp, bracket_minimum, bracket_minimum_result, &
    find_minimum, find_minimum_result
type(bracket_minimum_result) :: br
type(find_minimum_result) :: sol
br = bracket_minimum(objective, xm0=0.0_dp)
if (br%success) then
    sol = find_minimum(objective, [br%xl, br%xm, br%xr], xatol=1.0e-10_dp)
    if (sol%success) print *, sol%x, sol%f_x, sol%nfev
end if
```

Optional scalar tolerances `xatol`, `xrtol`, `fatol`, and `frtol` control
convergence; `maxiter=100` is the default. `callback=` receives the current
`find_minimum_result` and an output logical `stop`. The result contains the
minimizer, bracket points and function values, iteration/evaluation counts,
and termination status. Independent rank-one bracket arrays are also supported.
SciPy 1.17 numerical references are generated by
`tools/generate_optimize_find_minimum_reference.py`.


### Rank-one through rank-fifteen elementwise root finding and minimization

The same four names `bracket_root`, `find_root`, `bracket_minimum`, and
`find_minimum` also accept rank-one through maximum-rank batches. They return shape-preserving
arrays of their existing result types, with **per-element** statuses, iteration
counts and objective evaluation counts. The scalar overloads are unchanged.
For refinement, provide an array with shape **`(2,n)`** to `find_root` or
**`(3,n)`** to `find_minimum`; column `i` is the bracket for problem `i`.
For bracketing, `xl0(:)` and `xm0(:)` each describe `n` independent problems.
Optional endpoints and search bounds may be scalars (broadcast to `n`) or
rank-one arrays of length 1 or `n`. The four refinement tolerances may likewise
be scalars or arrays of length 1 or `n`, while `factor` and `maxiter` are shared
scalars.

```fortran
use scifort_optimize, only : dp, bracket_root, bracket_root_result, &
    find_root, find_root_result
real(dp) :: initial(3), roots_bracket(2,3)
type(bracket_root_result) :: br(3)
type(find_root_result) :: sol(3)
integer :: i

initial = [-4.0_dp, 0.0_dp, 3.0_dp]
br = bracket_root(f, initial, xmin=-10.0_dp, xmax=[-1.0_dp, 10.0_dp, 10.0_dp])
! Use a degenerate invalid bracket if an individual search failed.
do i = 1, 3
    if (br(i)%success) then
        roots_bracket(:,i) = [br(i)%xl, br(i)%xr]
    else
        roots_bracket(:,i) = [0.0_dp, 0.0_dp]
    end if
end do
sol = find_root(f, roots_bracket, xatol=1.0e-12_dp)
```

Each independent problem invokes the **scalar** `f(x)` callback, in input
order. Refinement callbacks are also invoked one problem at a time; they do
not receive the full batch. An invalid bracket or a failed search is recorded
for that element without terminating the other problems. Mismatched optional
argument lengths or unsupported ranks return an array of invalid-input
results without evaluating `f`. The additional rank-two interface preserves a two-dimensional problem grid.
Pass `xl0(m,n)` to `bracket_root`, `xm0(m,n)` to `bracket_minimum`,
`bracket(2,m,n)` to `find_root`, or `bracket(3,m,n)` to `find_minimum`.
Results have shape `(m,n)` and retain independent status, iteration, and
function-evaluation counts. Optional bounds, initial endpoints, and
refinement tolerances accept scalars, trailing-axis vectors of length `n`
(or length 1), and matrices with shape `(m,n)`, `(m,1)`, `(1,n)`, or `(1,1)`.
Incompatible shapes return all-invalid results without evaluating `f`.

```fortran
real(dp) :: brackets(2,2,3), xatol(3)
type(find_root_result) :: roots(2,3)

! Fill brackets(:,i,j) with the endpoints for problem (i,j).
xatol = [1.0e-12_dp, 1.0e-9_dp, 1.0e-7_dp]
roots = find_root(f, brackets, xatol=xatol)
```

Rank-three overloads likewise preserve an `(m,n,p)` problem grid:
`bracket_root` and `bracket_minimum` accept `(m,n,p)` starting points, while
`find_root` and `find_minimum` accept `(2,m,n,p)` and `(3,m,n,p)` brackets,
respectively. Optional real arguments broadcast from a scalar, a rank-one
vector aligned to `p`, a rank-two array aligned to `(n,p)`, or a rank-three
array aligned to `(m,n,p)`. Any aligned dimension may be a singleton.
Incompatible shapes return invalid-input result records without evaluating
the function. A bracket's leading endpoint dimension must be exactly 2 or 3.

```fortran
real(dp) :: root_brackets(2,2,3,4)
real(dp) :: x_tolerances(3,4)
type(find_root_result) :: roots(2,3,4)

! Each root_brackets(:,i,j,k) is a separate root search.
! x_tolerances(j,k) broadcasts across i.
roots = find_root(f, root_brackets, xatol=x_tolerances)
```

All four routines now use a common flattened batch engine and thin rank-specific
result wrappers. **`bracket_root` and `bracket_minimum` support problem grids
through rank 15**, the maximum rank in Fortran 2018. **`find_root` and
`find_minimum` support problem grids through rank 14** because their starting
brackets have an additional leading dimension (2 or 3), making the input
rank 15. Ranks 1 through 3 retain their original call shapes.

For a grid with shape `(n1,...,nr)`, an optional real argument may be scalar or
have rank `k <= r`. Its dimensions must match the last `k` grid dimensions,
except that dimensions of length 1 broadcast. The result has the same rank
and shape as the grid. An incompatible argument or a bracket with the wrong
leading endpoint count produces invalid-input result records without objective
calls. Empty grids return empty result arrays. For example:

```fortran
real(dp) :: initial(2,1,1,1,1,1,1,1,1,1,1,1,1,1,2)
real(dp) :: maximums(2)
type(bracket_root_result) :: br(2,1,1,1,1,1,1,1,1,1,1,1,1,1,2)

initial = 0.0_dp
maximums = [10.0_dp, 1.0_dp]
br = bracket_root(f, initial, xr0=1.0_dp, xmax=maximums)
! The rank-one 'maximums' argument broadcasts along the last axis.
```

A separate **vector-evaluation path**, `find_root_vectorized`,
uses a vector callback to evaluate all active root candidates together. It
runs the same Chandrupatla interpolation/bisection iteration as scalar
`find_root`, with independent statuses, iteration counts, and per-point `nfev`
counts. The callback receives *packed* active abscissae and returns exactly
one function value for each input (the batch length shrinks on convergence).
There is no stable lane-index argument: the function must apply the **same
scalar function elementwise** to every supplied value. A `callback=` for
monitoring and stopping still receives scalar `find_root_result` states.

```fortran
use scifort_optimize, only : dp, find_root_vectorized, find_root_result
real(dp) :: brackets(2,3)
type(find_root_result) :: roots(3)

brackets(:,1) = [0.0_dp, 1.0_dp]
brackets(:,2) = [0.0_dp, 2.0_dp]
brackets(:,3) = [-4.0_dp, 4.0_dp]
roots = find_root_vectorized(cubic_batch, brackets, xatol=1.0e-12_dp)

contains
subroutine cubic_batch(x, fx)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: fx(:)
    fx = x**3 - 1.0_dp
end subroutine cubic_batch
```

`find_root_vectorized` supports scalar/length-one/length-`n` tolerances, an
optional common `maxiter`, callback termination, and independently invalid
brackets or nonfinite outputs. `nfev` counts **values evaluated**, not the
number of vector callback invocations. Its API is deliberately separate from
`find_root` because Fortran cannot choose a generic overload solely by the
procedure argument's scalar-versus-vector signature. Rank-two grids are also supported; grids of rank three and above remain
future work for vector callbacks. The
regression and SciPy 1.17 reference generator are in
`test/test_optimize_elementwise_vectorized_root.f90` and
`tools/generate_optimize_elementwise_vectorized_root_reference.py`.

A corresponding **vector-evaluation path**, `find_minimum_vectorized`,
accepts `(3,n)` starting triples for `n` independent Chandrupatla minimizations.
The objective callback is `subroutine f(x, fx)` with rank-one real inputs
and outputs. The first call packs all three points of every valid bracket;
later calls contain only the trial abscissae of active problems. A scalar
`callback=` receives each lane's `find_minimum_result` and may stop that
lane alone. Scalar, singleton-vector, and per-lane tolerances broadcast,
and `nfev` counts the number of **individual objective values** evaluated.

```fortran
use scifort_optimize, only : dp, find_minimum_vectorized, find_minimum_result
real(dp) :: brackets(3,2)
type(find_minimum_result) :: minima(2)

brackets(:,1) = [0.0_dp, 1.0_dp, 3.0_dp]
brackets(:,2) = [-1.0_dp, 2.0_dp, 6.0_dp]
minima = find_minimum_vectorized(objective_batch, brackets, xatol=1.0e-12_dp)

contains
subroutine objective_batch(x, fx)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: fx(:)
    fx = (x-1.5_dp)**2 + 0.02_dp*(x-1.5_dp)**4
end subroutine objective_batch
```

Both vectorized refinement routines require the *same elementwise objective*
on every lane; callback order and batch size can change as lanes finish, so
the objective should not depend on stable lane indices. The new regression
and SciPy 1.17 reference generator are
`test/test_optimize_elementwise_vectorized_minimum.f90` and
`tools/generate_optimize_elementwise_vectorized_minimum_reference.py`.


**Vectorized bracket discovery** is also available in separate APIs:
`bracket_root_vectorized` and `bracket_minimum_vectorized`. Each accepts an
array of initial coordinates and the same packed callback signature
`subroutine f(x, fx)` as the vectorized refinement functions. The first call
evaluates two root endpoints or three minimum candidates per valid lane;
subsequent calls pack only pending left/right root searches or active minimum
searches. Scalar, singleton-vector, and length-`n` values are accepted for
`xr0`, `xl0`, `xmin`, `xmax`, and `factor` as applicable. One malformed
broadcast argument rejects the entire batch; an invalid value in a single
lane leaves other lanes running. A scalar `maxiter` is shared across lanes.

```fortran
use scifort_optimize, only : dp, bracket_root_vectorized, bracket_root_result
real(dp) :: starts(3)
type(bracket_root_result) :: intervals(3)

starts = [0.0_dp, 2.0_dp, -3.0_dp]
intervals = bracket_root_vectorized(linear_batch, starts, &
    xmin=-10.0_dp, xmax=10.0_dp)

contains
subroutine linear_batch(x, fx)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: fx(:)
    fx = x - 1.0_dp
end subroutine linear_batch
```

`bracket_minimum_vectorized` uses `xm0(:)` instead of `xl0(:)` and supports
optional `xl0`, `xr0`, `xmin`, `xmax`, and `factor` arguments. For both bracketing
routines, `nfev` counts **individual function values**, not calls to the
vector callback. Packed input positions do not carry stable lane indices, so
all lanes must use the same mathematical function. The independent scalar
`bracket_root` and `bracket_minimum` overloads, including high-rank shapes,
remain unchanged. The tests and generated SciPy 1.17 reference data are in
`test/test_optimize_elementwise_vectorized_bracket_root.f90`,
`test/test_optimize_elementwise_vectorized_bracket_minimum.f90`, and
`tools/generate_optimize_elementwise_vectorized_brackets_reference.py`.

All four `*_vectorized` functions also accept **rank-two problem grids**.
The bracketing routines take `xl0(m,n)` or `xm0(m,n)` and return `(m,n)`
results. The refinement routines take `bracket(2,m,n)` or `bracket(3,m,n)`
and return `(m,n)` results. Optional real endpoints, limits, growth factors,
and tolerances broadcast from scalars, length-`n` vectors, or arrays with
shapes `(m,n)`, `(m,1)`, `(1,n)`, or `(1,1)`. The objective callback still
takes `f(x(:), fx(:))` on **packed active points**, not a rank-two array;
result states and per-lane stopping callbacks retain `(m,n)` indexing.

```fortran
real(dp) :: root_brackets(2,2,3), atol(3)
type(find_root_result) :: roots(2,3)

! root_brackets(:,i,j) holds two endpoints for problem (i,j).
atol = [1.0e-12_dp, 1.0e-9_dp, 1.0e-7_dp]
roots = find_root_vectorized(cubic_batch, root_brackets, xatol=atol)
```

The numeric controller is the same as for rank-one vectorization, and it
keeps the same input validation, active-lane packing, per-element `nfev`,
and empty-grid behavior. All four packed-callback solvers now also accept
**rank-three and rank-four problem grids**. Starting-point arrays have shapes
`(m,n,k)` or `(m,n,k,p)`, and refinement brackets have leading endpoint
axes, e.g. `(2,m,n,k)` for root finding and `(3,m,n,k,p)` for minimization.
Results retain the input problem-grid dimensions. Optional real arguments
broadcast by aligning their axes with the trailing problem-grid axes, with
singleton axes expanded as necessary. The callback still accepts **packed
rank-one active-point arrays**, not an n-dimensional grid.

Malformed optional argument shapes reject the whole batch without objective
evaluation. Invalid input values affect only their corresponding lanes.
For **rank-five through rank-fifteen** packed-callback problem grids, use
four Fortran 2018 assumed-rank subroutines: `find_root_vectorized_into`,
`bracket_root_vectorized_into`, `find_minimum_vectorized_into`, and
`bracket_minimum_vectorized_into`. Each takes an input grid and a caller-
allocated result grid of matching shape. Refinement accepts problem grids
through rank 14 (its leading bracket axis uses the fifteenth dimension),
whereas bracketing starting-point grids can have rank 15. The existing
rank-one through rank-four function interfaces remain source-compatible,
calling the new dispatchers for ranks two through four. All packed objective
callbacks still receive rank-one arrays of active points.

```fortran
! Given a packed callback root_batch(x, fx):
real(dp) :: pair(2,2,1,2,1,2), atol(2)
type(find_root_result) :: roots(2,1,2,1,2)

atol = [1.0e-12_dp, 1.0e-10_dp]
! pair(:,i,j,k,l,m) contains the two starting endpoints.
call find_root_vectorized_into(root_batch, pair, roots, xatol=atol)
```

With a scalar result, `*_into` also accepts a scalar problem (refinement
brackets then have rank one). The output is `intent(inout)` solely to avoid a
GNU Fortran 14.2 compiler error for default-initialized derived-type arrays
with `intent(out)` and assumed rank. Results are reset on entry, including
when the input shape or broadcasting is invalid; unmodified neighboring
cells of a noncontiguous output slice remain untouched.

**Broadcasted objective parameters** are supported by the four
`*_vectorized_args_into` subroutines. Their callback has signature
`f(x(:), args(:,:), fx(:))`, with `args` shaped `(nargs, size(x))` for
exactly the currently active evaluation points. The caller supplies an
assumed-rank stack whose leading axis is the number of real parameters;
its remaining axes align with the **trailing problem-grid axes** and may
have size one for broadcasting. Parameters are repacked with their
original lane identifiers whenever the active set shrinks. Existing
`f(x, fx)` functions and `*_vectorized_into` calls are unchanged.

```fortran
! Solve x = center for two-by-three independent problems.
real(dp) :: intervals(2,2,3), params(2,1,3)
type(find_root_result) :: answers(2,3)
! Each interval in intervals(:,i,j) must straddle params(1,1,j).
params(1,1,:) = [0.2_dp, 0.8_dp, 1.4_dp]  ! centers
params(2,1,:) = 1.0_dp                  ! scale factors
call find_root_vectorized_args_into(root_with_params, intervals, params, answers)
! root_with_params(x, args, fx): fx = args(2,:) * (x-args(1,:))
```

A rank-one parameter stack broadcasts constants across every problem.
A scalar (`rank(0)`) or empty parameter stack is invalid. As with the
existing tolerances, incompatible shapes return invalid-input results
without calling the objective. Higher-rank grids obey Fortran's rank-15
limit, including the parameter axis. See
`test/test_optimize_elementwise_vectorized_args.f90` and
`test/test_optimize_elementwise_vectorized_assumed_rank.f90`.

**Independently shaped real objective arguments** are also accepted by
those same four `*_vectorized_args_into` generics. Construct one
`elementwise_argument` per value, using the Fortran 2018 assumed-rank
`elementwise_argument_from(...)` helper. Unlike a stacked numeric array,
each argument may have its *own* broadcast shape (or be scalar):

```fortran
! Problem grid (2,3), objective f(x, a, fx):
real(dp) :: centers(3), offset(2,1), scale(1,3), intervals(2,2,3)
type(elementwise_argument) :: params(3)
type(find_root_result) :: answers(2,3)
params = [elementwise_argument_from(centers), &
          elementwise_argument_from(offset), &
          elementwise_argument_from(scale)]
! Set intervals(:,i,j) to a root-straddling pair.
call find_root_vectorized_args_into(f, intervals, params, answers)
! In f, a(1,:) gives center(j), a(2,:) gives offset(i,1),
! and a(3,:) gives scale(1,j) for active lane (i,j).
```

The overloaded API retains the existing packed callback signature
`f(x(:), args(:,:), fx(:))` and handles active-lane compaction without
parameter drift. All four solvers support these independently shaped
arguments, including scalar, noncontiguous and maximum-rank grids.
The numeric argument list may not be empty; malformed or incompatible
shapes yield invalid-input results without calling the objective. This
existing overload continues to use `real(dp)` arguments and packed rank-two
callback arrays. See `test/test_optimize_elementwise_vectorized_arg_list.f90`.

**Independently shaped mixed-type arguments** use the same four generic
`*_vectorized_args_into` names with `type(elementwise_typed_argument)` and
`elementwise_typed_argument_from(...)`. Each argument may contain an integer,
logical, character, real, complex, or user-defined type, and may have a
separate scalar/array shape of rank 0--15. The callback receives `x(:)` and
an array of typed arguments, plus `fx(:)`; each active lane's original
argument is available as `args(i)%values(k)%value` (argument `i`, active
point `k`) for `select type`. For example:

```fortran
! For a grid of two independent roots, solve x = [2, 3].
real(dp) :: pairs(2,2)
integer :: targets(2)
type(elementwise_typed_argument) :: mixed(2)
type(find_root_result) :: answers(2)

targets = [2, 3]
pairs(:,1) = [0.0_dp, 4.0_dp]
pairs(:,2) = [0.0_dp, 4.0_dp]
mixed = [elementwise_typed_argument_from(targets), &
         elementwise_typed_argument_from(.true.)]
call find_root_vectorized_args_into(root_with_mixed, pairs, mixed, answers)

! The callback root_with_mixed(x, args, fx) contains:
! do k = 1, size(x)
!     select type(v => args(1)%values(k)%value)
!     type is (integer)
!         fx(k) = x(k) - real(v, dp)
!     end select
! end do
```

All four root/minimum refinement/bracketing controllers preserve the
correspondence between the typed arguments and each packed active point,
including when initial endpoints repeat a lane ID or converged lanes drop
out. Fortran 2018 `select rank` handles input arrays through rank 15.
Malformed typed objects, empty lists, and incompatible broadcast shapes
return invalid-input results without calling the objective. Each typed
scalar is boxed and copied on every evaluation: prefer the existing
`real(dp)` stack/list overloads in numerical hot loops when sufficient.
`test/test_optimize_elementwise_typed_args.f90` covers mixed intrinsic and
user-defined types, strings, scalar/high-rank broadcasting, lane compaction,
noncontiguous outputs, and invalid inputs.

The four existing generic routines use independently evaluated
**scalar-kernel wrappers**, not SciPy's vectorized evaluation loop. A scalar objective and scalar refinement
callback are invoked for each element in Fortran column-major order; these scalar overloads do not accept vector-valued callbacks or additional function arguments.

`test/test_optimize_elementwise_arrays.f90` covers mixed statuses, bounds,
scalar/vector broadcasting, callbacks and invalid shapes; its generated
SciPy 1.17 reference data are in
`test/optimize_elementwise_arrays_reference.f90`.
The rank-two regression and SciPy 1.17 generator are
`test/test_optimize_elementwise_rank2.f90` and
`tools/generate_optimize_elementwise_rank2_reference.py`. The rank-three
regression and reference generator are
`test/test_optimize_elementwise_rank3.f90` and
`tools/generate_optimize_elementwise_rank3_reference.py`, including broadcasting,
callback isolation, empty dimensions and malformed argument shapes. The
`test/test_optimize_elementwise_highrank.f90` regression exercises **every
batch rank from 1 through 15**, including the rank-14 refinement ceiling,
rank-15 trailing-axis broadcasting, callbacks, failed elements, and empty
grids. It is based on scalar reference comparisons using the same existing
SciPy-validated kernels.

Multidimensional nonlinear equations add `root` with MINPACK `hybr` and `lm` methods plus
`fsolve`, using the community-maintained `fortran-lang/minpack` FPM dependency.
SciPy-style inexact-Newton methods `broyden1`, `broyden2`, `anderson`,
`diagbroyden`, `linearmixing`, `excitingmixing`, and `krylov` are also available
through `root`; the first six have direct convenience functions and Krylov is
also exposed as `newton_krylov`. Their iteration controller, Armijo and strong-
Wolfe line searches, autoscaling, and update rules are adapted from SciPy 1.17.0.
Broyden 1/2 support `max_rank=` with `restart`, `simple`, or `svd` reduction.
The Broyden rank-one formula was cross-checked against Jacob Williams'
nlesolver-fortran 1.1.0 with attribution. Krylov uses matrix-free finite-
difference Jacobian products with SciPy-compatible LGMRES (default) or GMRES;
Ivan Pribec's BSD-3-Clause `kiss` project was consulted as a modern-Fortran
Krylov API/reference without copying its implementation. `root(method='df-sane')`
adds SciPy's derivative-free spectral residual method with the default Cruz or
optional Cheng nonmonotone line search. Because Fortran is case-insensitive and
`m=` is already used for LM residual dimension, the DF-SANE history-window option
is named `sane_m=` rather than SciPy's `M`. The seven inexact-Newton/Krylov
methods now support an optional `tol_norm=` callback for convergence tests;
it returns a nonnegative real(dp) scalar from an input real(dp) vector.
The default is the maximum absolute component norm. The selected norm is used
for residual, iterate and accepted-step stopping tests, including the initial
residual norm used in relative tolerance checks. Nonfinite or negative norms
return an unsuccessful result rather than silently declaring convergence.
Krylov supports `inner_m=` preconditioners, `inner_m_update=` adaptive
preconditioner refreshes, and custom `inner_solver=` callbacks, along with
built-in GMRES, LGMRES, BiCGSTAB, CGS, MINRES, and TFQMR. MINRES assumes a
symmetric Jacobian and positive-definite preconditioner. **SciFort extension:**
`inner_m_diag_rank_one=.true.` can add a verified Sherman-Morrison
correction to exact automatic Krylov diagonal preconditioning for
`J = diag(a) + u*ones^T` mean-field coupling. A fixed
`inner_m_diag_rank_one_right=v` additionally supports known weighted
couplings `J = diag(a) + u*v^T`, including mixed-sign right factors. For
unknown rank-one right weights, `inner_m_diag_rank_one_discover=.true.`
recovers the coupling from n exact coordinate Jacobian products, reuses
them for diagonal scaling, checks every off-diagonal entry, and falls
back safely on unidentifiable or incompatible matrices. This optional
full-Jacobian diagnostic costs O(n^2) memory rather than O(n).
A separate opt-in checked rank-two Woodbury correction supports known
`J = diag(a) + U*transpose(V)` couplings with
`inner_m_diag_rank_two=.true., inner_m_diag_rank_two_right=V`, where `V`
has two columns. When both right factors are unknown, the new
`inner_m_diag_rank_two_discover=.true.` option (without
`inner_m_diag_rank_two_right`) fits a diagonal-plus-rank-two model from
n exact coordinate Jacobian products, checks every off-diagonal entry,
and falls back safely when the fit is not identifiable. Discovery stores
the complete Jacobian (O(n^2) memory), so it targets small and moderate
systems. **SciFort extension:** `inner_m_diag_rank_two_stream=.true.`
with discovery uses at most eight pilot columns and processes the rest
one at a time. It still uses n Jv products and O(n^2) arithmetic for
exhaustive off-diagonal checking, but only O(n) temporary discovery
storage. If pilot columns are uninformative or the fit is unsafe, it
retains exact diagonal scaling. See `docs/api.md` for limitations. When the uncoupled diagonal `a` is already known and `J-diag(a)` is
symmetric rank two, **SciFort extension**
`inner_m_diag_rank_two_symmetric_base=a` with the rank-two discovery
flags fits and verifies the correction with four Jv products using
O(n) workspace. Invalid or rejected fits fall back to the exact
Jacobian diagonal. This assumption-dependent, heuristically verified
path is documented in `docs/api.md` and is disabled by default.
With a known uncoupled diagonal but **nonsymmetric** rank-two coupling,
`inner_m_diag_rank_two_nonsymmetric_base=a` and an analytic
`jacobian_transpose_product=jtvec` callback enable a separate SciFort
extension. It uses four forward Jacobian-vector products and two
transpose products (including verification), O(n) discovery memory,
and exact-diagonal fallback on rejection. Neither the coupling factors
nor a dense Jacobian need be supplied. The additional opt-in
`inner_m_diag_rank_two_nonsymmetric_verify_transpose=.true.` checks two
held-out transpose products, increasing accepted setup to four forward
plus four transpose actions. This detects some inconsistent adjoints
invisible to the original fitting probes; rejection still uses the exact
diagonal. Alternatively, the opt-in
`inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost=c` selects
one or two held-out adjoint checks based on the declared relative
`J^T*v` cost, system size, and fit quality. For example, `c=16.0_dp`
uses one held-out check on well-conditioned 32-variable problems but
two on 128-variable problems. Neither option changes defaults. The
heuristic paths require an adjoint matching the forward Jacobian;
see `docs/api.md` and `benchmarks/benchmark_krylov_rank_two_discovery.f90`.
Work-aware Krylov refresh can measure the setup cost for **all** automatic
diagonal and low-rank preconditioners, including known-right rank-one,
rank-two sketch, streamed/dense discovery, symmetric and nonsymmetric
low-product fits, and ordinary randomized/exact diagonal probes. Enable
`inner_m_diag_update_measured_cost=.true.` with adaptive work-growth
scheduling, and optionally set `inner_m_diag_update_transpose_cost=` to
weight adjoint products. Counts include fitting, verification, retries,
and exact-diagonal fallback. This removes the need to guess a fixed rebuild
cost, but the comparison with inner preconditioner applications is still
heuristic; see `docs/api.md`,
`test/test_optimize_root_krylov_measured_setup_modes.f90`, and the
132-case `benchmarks/benchmark_krylov_work_refresh.f90`.
As an alternative, the opt-in
`inner_m_diag_update_timed_cost=.true.` directly measures the CPU time
of preconditioner construction and each inner linear solve. It compares
current inner solve time with the first solve after a rebuild, accounting
for the preceding rebuild time, and applies a minimum timing interval
(default 0.001 seconds) to avoid noise-induced refreshes. This mode
requires adaptive work-growth scheduling and cannot be combined with
manual or product-count rebuild costs. Process CPU time is not wall time;
see `docs/api.md` and `test/test_optimize_root_krylov_timed_refresh.f90`.
For OpenMP or other parallel callbacks, the opt-in
`inner_m_diag_update_wall_clock=.true.` switch, used with
`inner_m_diag_update_timed_cost=.true.`, instead compares wall seconds
measured by `system_clock`. The existing minimum sample duration still
applies. It rejects unavailable/nonmonotone clock samples, retains fixed
interval fallback, and does not affect default CPU-timed behavior. See
`docs/api.md` and `test/test_optimize_root_krylov_wall_clock_refresh.f90`.
A separate 32-case CPU-time comparison is available in
`benchmarks/benchmark_krylov_timed_refresh.f90` and
`benchmarks/krylov_timed_refresh_gfortran14_o2.csv`. It compares fixed
refresh intervals, manual and measured work costs, and CPU-timed refresh
on the same nonlinear rank-two systems. Two *synthetic*, single-threaded
callback workloads model expensive forward products or expensive transpose
products, at 128 variables and two nonlinearity rates. Each case checks
convergence and the known solution, while a disabled timing threshold is
compared to the initial-only schedule. Timings depend on the machine and
should not be extrapolated to application callbacks; see the accompanying
CSV and `ROADMAP.md` for remaining validation work.
An additional benchmark, `benchmark_krylov_quadrature_refresh.f90`,
replaces synthetic CPU delays with actual midpoint quadrature of local
nonlinear exponential responses. It supports serial compilation or optional
OpenMP threading of the residual and analytic Jacobian-product callbacks,
and reports both `cpu_time` and `system_clock` wall time for the same
nonlinear rank-two problems. The original CSV combines 36 GNU Fortran 14.2 `-O2`
runs at one thread and 36 at four threads. The additional
`krylov_quadrature_refresh_wall_gfortran14_o2.csv` contains 42 runs per
thread count, including the wall-timed refresh mode. Process CPU time counts thread
work and does not necessarily track wall-clock speedups; the timing
heuristic may therefore schedule different preconditioner rebuilds.
The additional opt-in
`inner_m_diag_rank_two_sketch=.true.` reconstructs the diagonal from three
probe directions and uses the same four structural Jv checks, reducing
accepted rank-two setup from `n+4` to 7 Jacobian-vector products. Unreliable
sketches retry the exact-diagonal path. Both models revert to diagonal
scaling on failed checks; see `docs/api.md` for restrictions and caveats. The four added
inner methods have SciPy 1.17.0 nonlinear work-count regressions. TFQMR offers
an opt-in `inner_tfqmr_right_preconditioned=.true.` variant that uses `J*M`
with `M`-applied iterate updates to improve accuracy with approximate inverse
preconditioners; the original SciPy-compatible `M*J` recurrence remains the
default. It can be combined with `inner_verify_residual=.true.` to check the
original linear residual. See `docs/api.md` for details and SciPy references.
**SciFort extension -- GMRES/LGMRES reorthogonalization:** pass
`inner_orthogonalization='selective'` for a second modified Gram-Schmidt pass
when the new Arnoldi vector loses at least half its norm in the first pass,
or `'twice'` for an unconditional second pass. The `'once'` default retains
existing SciPy-compatible GMRES/LGMRES behavior. Both extra passes reuse the
same `J*v` and therefore add no function evaluations. For example:

```fortran
res = root(fun, x0, method='krylov', inner_method='lgmres', &
    inner_orthogonalization='selective')
```

The option works with built-in GMRES and LGMRES (including preconditioning)
but not with custom inner solvers or BiCGSTAB/CGS/MINRES/TFQMR. It improves
basis orthogonality for ill-conditioned systems; it does not guarantee that
the nonlinear residual will be smaller. See `docs/api.md` for details.

A SciFort-specific `jacobian_diff_scheme='3-point'` option is available for
`root(method='krylov')` and `newton_krylov`: a centered directional
Jacobian-vector difference with an automatic cube-root-epsilon step rather than
the default SciPy-style two-point, square-root-epsilon difference. It uses two
function evaluations per nonzero Jacobian-vector product, improves accuracy for
smooth residuals, and detects perturbations too small to be represented.
All built-in methods and custom `inner_solver=` callbacks receive the selected
stencil; `rdiff=` still overrides the automatic step. See `docs/api.md`.
**SciFort extension -- adaptively scaled Krylov Jacobian differences:**
For strongly differently scaled variables, the default single directional
step can be inaccurate or leave a residual function's domain. Pass
`jacobian_step_mode='adaptive'` to `root(method='krylov')` or `newton_krylov`:

```fortran
res = root(fun, x0, method='krylov', inner_method='lgmres', &
    jacobian_step_mode='adaptive', jacobian_diff_scheme='3-point')
```

The default mode, `'global'`, retains SciPy-compatible directional stepping
and work counts. `'adaptive'` retains global stepping for well-scaled
Jacobian-vector products but switches to per-coordinate perturbations when
an active coordinate is much smaller than the largest state component, when
the proposed perturbation rounds back to the same point, or when it is far
smaller than the local `rdiff` scale. `'componentwise'` always uses these
coordinatewise perturbations. With either `'2-point'` or `'3-point'`, the
coordinatewise step is `max(rdiff, 8*epsilon) * max(1,abs(x(j)))`, with
`rdiff` automatically selected for the stencil if omitted. This uses one
(`'2-point'`) or two (`'3-point'`) residual evaluations per *nonzero*
component of the direction, so it can be expensive for high-dimensional
dense systems; evaluation counts intentionally differ from SciPy when the
coordinatewise path is selected. A negative-domain or nonfinite function
sample fails the inner solve without accepting a false convergence result.
The option works with built-in or custom inner solvers and preconditioners;
it is only valid for Krylov root finding.

**SciFort extension -- exact Krylov Jacobian products:** If the application
can multiply the Jacobian by a vector analytically (or through automatic
differentiation), pass `jacobian_product=` to avoid finite-difference
cancellation and its residual evaluations:

```fortran
res = newton_krylov(fun, x0, jacobian_product=my_jacobian_product)

! Callback interface:
! subroutine my_jacobian_product(x, vec, product, ok)
!     real(dp), intent(in) :: x(:), vec(:)
!     real(dp), intent(out) :: product(:)
!     logical, intent(out) :: ok
! end subroutine
```

Here `product` is `J(x)*vec`, at the *current* nonlinear iterate; `ok` must be
false if the product cannot be computed. The callback is supported by all six
built-in inner solvers and custom `inner_solver=` implementations, including
preconditioned/verified TFQMR. It is **not** a SciPy keyword: SciFort adds it
as a Fortran-specific capability. The default remains SciPy-compatible
finite differences. An analytic product does not increment the residual
function evaluation count `nfev`. Explicit `rdiff=`, `jacobian_diff_scheme=`,
and `jacobian_step_mode=` cannot be supplied simultaneously with
`jacobian_product=`; see `docs/api.md` for the full interface and a test.

The MINPACK adapters support numerical or analytic Jacobians, `hybr` band
structure, variable scaling, overdetermined LM residual systems, and SciPy-style
termination status and work-count reporting. For overdetermined `lm` problems,
pass `m=` explicitly because a Fortran callback does not expose its output size
before the call. MINPACK pivot indices remain one-based in the Fortran result
object.

Multivariate minimization provides `minimize`, `nelder_mead`, `powell`, `bfgs`,
`cg`, `newton_cg`, `dogleg`, `trust_ncg`, `trust_exact`, `trust_krylov`, `trust_constr`, `tnc`, `lbfgsb`, and `slsqp`. With neither bounds nor constraints,
`minimize` defaults to BFGS; bounds select L-BFGS-B, while an explicit combined
constraint callback selects SLSQP. With an explicit method, Nelder-Mead and Powell also
accept variable bounds. Nelder-Mead supports adaptive coefficients, a caller-supplied
initial simplex, and SciPy-style clipping/reflection at finite bounds. Powell accepts a
caller-supplied direction set and uses bounded or one-sided transformed line searches so
feasible starts remain feasible. BFGS/CG accept either an analytic gradient or SciPy-style forward differences.
Newton-CG requires an analytic gradient and accepts a dense Hessian,
Hessian-vector product, or finite-difference Hessian-vector products. TNC adds
bound-constrained truncated Newton minimization with scaled variables, an active-set
projected gradient, BFGS-preconditioned truncated CG with finite-difference
Hessian-vector products for analytic gradients, and a bounds-aware finite-difference/BFGS
fallback when no gradient callback is supplied.
It accepts `maxfun` and `max_cg_iter` controls; exact TNC work counts can differ from
SciPy's legacy C implementation. `dogleg`
implements SciPy's classical dog-leg trust-region method for a positive-definite
analytic Hessian, using the Newton point, Cauchy point, and their trust-boundary
intersection. `trust_ncg` implements SciPy's Newton-CG trust-region method with an
analytic gradient plus either a dense Hessian or Hessian-vector product. `trust_exact`
implements SciPy's nearly-exact trust-region method for an analytic gradient and
dense analytic Hessian, including indefinite-Hessian hard cases and a
`subproblem_maxiter` control for the inner damping iteration. `trust_krylov`
implements SciPy's GLTR-style truncated Krylov trust-region method with an
analytic gradient plus either a dense Hessian or Hessian-vector product;
`inexact=.true.` uses SciPy/TRLIB's default forcing sequence while
`inexact=.false.` requests tighter reduced-subproblem convergence. Its Krylov
basis is retained across rejected steps so trust-radius retries hot-start as in
SciPy. These trust-region methods expose `initial_trust_radius`,
`max_trust_radius`, `eta`, `gtol`, and `maxiter`; `trust_ncg` also includes
Steihaug negative-curvature boundary steps. BFGS and Newton-CG expose strong-Wolfe controls. L-BFGS-B is backed by Jacob Williams'
BSD-3 `lbfgsb` 1.0.0 FPM package and supports analytic or bounds-aware numerical
gradients, projected-gradient/function tolerances, correction-history control,
and SciPy-style iteration/evaluation limits. SLSQP is backed by Jacob Williams'
BSD-3-Clause `slsqp` 1.6.0 FPM package; its combined constraint callback stores
`meq=` equalities first and then `mineq=` inequalities, with inequalities feasible
when `c(x) >= 0`. Fortran cannot infer the constraint-vector length before a
callback, so `meq=`/`mineq=` are explicit. COBYLA is a native modern-Fortran
PRIMA-derived trust-region/simplex implementation. It uses the same combined
constraint callback: `meq=` leading equalities are converted to paired
inequalities and the following `mineq=` entries are feasible when `c(x) >= 0`.
COBYLA also supports finite bounds, `rhobeg`, `tol`, `catol`, `f_target`, and
`maxiter` (maximum objective evaluations), and reports maximum constraint
violation in `multivariate_minimize_result%maxcv`. `multivariate_minimize_result`
reports the final point/value, function/gradient/Hessian work counts, status,
and method-specific diagnostics. `trust_constr` implements dense
Byrd-Omojokun-style equality-constrained SQP plus a dense trust-region
interior-point layer for inequalities and finite bounds. The existing combined
constraint callback uses `meq=` leading equalities and `mineq=` trailing
inequalities feasible when `c(x) >= 0`. Alternatively, `constraint_lower=` and
`constraint_upper=` give SciPy-style per-component bounds on the raw callback
values, supporting equalities, either one-sided inequality, and two-sided intervals;
this bound-vector form is not combined with `meq`/`mineq`. Nonfixed
inequalities/bounds are converted to positive slack equalities and solved through a decreasing log-barrier
sequence, while fixed bounds become equalities. Analytic or finite-difference
objective/constraint Jacobians are supported, together with an analytic objective
Hessian, finite-difference Hessians of analytic gradients, or BFGS/SR1 Hessian updates,
pseudoinverse projections,
and SciPy-style second-order correction (SOC) for rejected equality-SQP steps,
plus `gtol`/`xtol` stopping, `initial_tr_radius`,
`initial_constr_penalty`, `barrier_tol`, `initial_barrier_parameter`,
`initial_barrier_tolerance`, `maxcv`, and final `constr_penalty` reporting.
Redundant equality rows are supported even if `meq > size(x0)` (or if
fixed bounds make the interior-point equality system overdetermined). The
dense solver uses a Gram-eigendecomposition pseudoinverse to obtain
minimum-norm projection steps and multipliers, similar in purpose to SciPy's
`factorization_method='SVDFactorization'`. This is not a sparse SVD solver;
it is appropriate for modest dense systems, and near-dependent constraints
may be less robust than a direct SVD. Incompatible equations do not count as
success: remaining feasibility error produces termination status `4`. Tests
and SciPy reference generation are in
`test/test_optimize_trust_constr_redundant.f90` and
`tools/generate_optimize_trust_constr_redundant_reference.py`.

With no constraints and no variable bounds, `trust_constr` (or
`minimize(method='trust-constr')`) uses the same equality-SQP algorithm with
zero constraint rows. For example,

```fortran
r = trust_constr(objective, x0, jac=gradient, hess=hessian)
```

The objective Hessian may alternatively be supplied as a Hessian-vector
product `hessp(x, p, hp)`, with `hp = matmul(H(x), p)`:

```fortran
r = trust_constr(objective, x0, jac=gradient, hessp=hessian_product)
```

`hessp=` is supported in both direct `trust_constr` and generic
`minimize(method='trust-constr')`, including unconstrained, equality-constrained,
nonlinear constraint-object, and inequality/bound barrier paths. Since SciFort's
SQP implementation is currently dense, it assembles the full symmetric
objective Hessian from `size(x0)` products per Hessian evaluation rather than
using a matrix-free solver. `nhev` counts these dense Hessian assemblies, not
the individual Hessian-vector products. If both `hess=` and `hessp=` are given,
`hess=` takes precedence as in SciPy. Exact multiplier-weighted constraint
Hessians can be combined with the assembled objective Hessian.

When the objective has an analytic gradient but no analytic Hessian, SciFort
also supports SciPy's finite-difference Hessian conventions. Specify
`objective_hess_fd_scheme='2-point'` or `'3-point'` for either `trust_constr`
or `minimize(method='trust-constr')`:

```fortran
r = minimize(objective, x0, method='trust-constr', jac=gradient, &
    objective_hess_fd_scheme='3-point')
```

The Hessian is reconstructed columnwise from the analytic `jac=` callback.
Two-point differences use first-order forward/backward probes, following
the sign of the coordinate or an explicit signed relative step; three-point
differences use centered or second-order one-sided probes, with feasible
steps at finite bounds and zero columns for variables fixed by bounds.
`objective_finite_diff_rel_step(:)` sets the per-variable relative steps
for this calculation, otherwise the two-point default is
`sqrt(epsilon(1.0_dp))` and the three-point default is
`epsilon(1.0_dp)**(1.0_dp/3.0_dp)`, scaled by `max(1,abs(x(i)))`.
The resulting Hessian is symmetrized. An analytic gradient is required;
`objective_hess_fd_scheme` cannot be combined with objective `hess=`,
`hessp=`, or `objective_hess_update=`. Other `minimize` methods reject the
option. As in SciPy, `njev` includes the gradient evaluations used to form
these Hessians, while `nhev` stays zero. Exact constraint Hessians are
still recorded independently in `constr_nhev`. Since SciFort uses a dense
SQP implementation, the Hessian is formed explicitly rather than stored
as a matrix-free operator. Reference results and reproducible tests are in
`tools/generate_optimize_trust_constr_hess_fd_reference.py` and
`test/test_optimize_trust_constr_hess_fd.f90`.

When no exact objective Hessian is supplied, SciFort defaults to BFGS. To
select SciPy's symmetric-rank-one (SR1) Hessian-update strategy instead, pass
`objective_hess_update='sr1'` to `trust_constr` or
`minimize(method='trust-constr')`:

```fortran
r = trust_constr(objective, x0, jac=gradient, &
    objective_hess_update='sr1')
```

`'bfgs'` and `'sr1'` are accepted case-insensitively. SR1 applies a symmetric
rank-one secant correction with SciPy's `1.0e-8` denominator safeguard and
first-step automatic scaling; unlike BFGS, SR1 can represent negative
curvature. The SR1 objective update is kept separate from a supplied exact
constraint Hessian or a BFGS approximation for missing constraint curvature.
It works with equality constraints, independent constraint blocks, inequality
barriers, variable bounds, and numerical objective gradients. An explicitly
selected update cannot be combined with an exact objective `hess=` or
`hessp=`; an invalid update name is rejected. Without a supplied Hessian,
`nhev` remains zero, since the updates do not call an objective Hessian.
Tests and reproducible references are in
`test/test_optimize_trust_constr_sr1.f90` and
`tools/generate_optimize_trust_constr_sr1_reference.py`.

For a numerical objective gradient with a second-order stencil, specify
`objective_fd_scheme='3-point'` instead of a `jac=` callback:

```fortran
r = trust_constr(objective, x0, hess=hessian, &
    objective_fd_scheme='3-point')
```

The default `objective_fd_scheme='2-point'` retains the existing forward
finite differences. Three-point differences use central probes in the
interior and second-order one-sided probes at strict `bounds_keep_feasible`
endpoints (including positive barrier slacks). The default three-point
relative step is `epsilon(1.0_dp)**(1.0_dp / 3.0_dp)`; optional `eps=`
overrides its size, scaled by `max(1, abs(x(i)))`. This objective option is
forwarded through `minimize(method='trust-constr')` and `constraint_objects=`.
Each extra objective probe contributes to `nfev`, with no analytic `jac=`
callback required. Objective finite-difference scheme selection is currently
specific to `trust-constr`; other minimizers reject that keyword rather than
silently ignoring it. Each constraint object's `%fd_scheme` separately
controls its own numerical Jacobian.

To control the objective-gradient step separately for each variable, pass
`objective_finite_diff_rel_step(:)` to either `trust_constr` or
`minimize(method='trust-constr')`. This also works with
`constraint_objects=`:

```fortran
r = trust_constr(objective, x0, hess=hessian, &
    objective_fd_scheme='3-point', &
    objective_finite_diff_rel_step=[1.0e-4_dp, 5.0e-5_dp])
```

For an explicitly supplied relative step, the trial displacement in variable
`i` is `rel_step(i) * sign(x(i)) * abs(x(i))`, matching SciPy's numerical
differentiation convention. Three-point differences use its absolute value;
two-point differences preserve the sign. If a zero coordinate or sub-ULP step
would give no displacement, the default (or `eps=`) step is substituted.
Retained variable bounds can reverse or shorten a probe. The step vector must
have `size(x0)` finite entries; negative and zero entries are allowed (the
latter triggers the usual fallback). If an analytic `jac=` is supplied, it
takes precedence. This setting does not affect the constraint-block Jacobians
and is rejected by other `minimize` methods. In the bounded `trust-constr`
interior-point solver, a negative two-point objective step is preserved even
when a positive direction is also feasible; only a conflicting bound forces
reversal or shortening. This holds for both the initial barrier gradient and
later SQP iterations. Probe directions and gradient values are verified against
SciPy 1.17 in `test/test_optimize_trust_constr_objective_signed_bounds.f90`,
with references generated by
`tools/generate_optimize_trust_constr_objective_signed_bounds_reference.py`.

The Hessian may be omitted to use the existing BFGS approximation, and the
gradient may be omitted to use numerical differences. An empty
`constraint_objects=` array is equivalent to specifying no constraint blocks;
if variable bounds are also given, the existing bounds-only solver is used.
Unconstrained results have zero-length `constr`, `v`, and per-constraint work
count arrays, and a zero-row `constr_jac` matrix, matching SciPy's reporting.
The older SciFort spelling `initial_trust_radius` remains accepted by
`trust_constr`; specifying both trust-radius spellings is invalid. As in SciPy, generic `tol=`
also supplies `barrier_tol` for `method='trust-constr'` unless that option is
given explicitly. For nonlinear equality and inequality constraints, optional
`constraint_hess(x, v, hmat)` supplies the Hessian of the multiplier-weighted
original constraint combination, following SciPy's `NonlinearConstraint.hess(x, v)`
convention. The slack-variable interior-point path maps its equality multipliers
back to the original constraint vector; simple/fixed bounds contribute zero
constraint Hessian. When the objective has no supplied Hessian or `hessp=`,
its BFGS curvature is combined with the exact constraint Hessian, rather than
approximating the constraint Hessian a second time. SOC is
enabled for equality and barrier/slack subproblems; the latter use SciPy-style
scaled slack coordinates and the `0.995` fractional slack trust box. The current
bound-vector constraint form also accepts `constraint_keep_feasible=` per raw
constraint component, and variable bounds accept `bounds_keep_feasible=` per
variable. Equality components ignore the flag, infeasible kept inequality starts
are rejected, and kept barrier slacks are tied to the corresponding constraint
gap during trial evaluation. Numerical objective-gradient and nonlinear-constraint
Jacobian probes also use one-sided or shortened finite-difference steps so retained
variable bounds are not crossed. Multiple independent linear and nonlinear
constraint blocks are now accepted by `trust_constr` and
`minimize(method='trust-constr')` using `constraint_objects=` and an array of
`nonlinear_constraint` values. Nonlinear blocks have a procedure pointer
`%fun` and optional `%jac`, `%hess`, `%hessp` callbacks. Every block has allocatable
`%lower`, `%upper`, and `%keep_feasible` arrays. At least one bound array is
required to determine each block's output size; a missing endpoint means
unbounded. For example:

```fortran
type(nonlinear_constraint) :: blocks(2)
type(multivariate_minimize_result) :: r

blocks(1)%fun => equality_values
blocks(1)%jac => equality_jacobian
blocks(1)%lower = [2.0_dp]
blocks(1)%upper = [2.0_dp]
blocks(2)%fun => inequality_values
blocks(2)%lower = [0.5_dp]
blocks(2)%keep_feasible = [.true.]
r = minimize(objective, [0.6_dp, 1.4_dp], method='trust-constr', &
    constraint_objects=blocks)
```

To supply the Hessian of the multiplier-weighted nonlinear constraint as
products rather than as a full matrix, set `%hessp` on a nonlinear block:

```fortran
blocks(1)%hessp => weighted_constraint_product
! weighted_constraint_product(x, v, p, hp) computes
! hp = matmul(Hessian(dot_product(v, constraint_values(x))), p)
r = trust_constr(objective, x0, jac=gradient, hess=hessian, &
    constraint_objects=blocks)
```

This corresponds to using a SciPy `LinearOperator` returned by a
`NonlinearConstraint.hess(x, v)` callable. It is **not** a new SciPy keyword;
`%hessp` is SciFort's explicit Fortran representation of that operator.
Constraint Hessian products receive only their own block's multipliers.
Different nonlinear blocks can use full `%hess` or product `%hessp` callbacks
in the same problem, and linear blocks continue to contribute zero Hessians.
If both callbacks are associated on a block, `%hess` takes precedence.
Since the SQP implementation is dense, each Hessian assembly invokes the
block's product `size(x0)` times, then symmetrizes the assembled matrix.
`constr_nhev` counts Hessian assemblies, *not* individual products.
Exact constraint Hessians can be combined with either an analytic objective
Hessian (possibly assembled from `hessp=`) or an objective BFGS/SR1 update.
Blocks with exact and missing constraint Hessians **can now be mixed**: a missing Hessian is
reconstructed by central differences of that block's Jacobian, with second-order
one-sided stencils near retained variable bounds. A block without `%jac` uses
its configured finite-difference Jacobian inside this approximation. This
is a SciFort dense numerical-Hessian fallback, rather than SciPy's BFGS
Hessian-update strategy; it can require substantially more evaluations. The
other blocks' callbacks are not called while assembling that Hessian.
`constr_nhev` is zero for blocks without an exact Hessian callback, while the
underlying Jacobian/function probes contribute to `constr_njev`/`constr_nfev`.

Alternatively, set a nonlinear block's `%hess_update = 'bfgs'` or `'sr1'`
to approximate its constraint Hessian from successive Jacobians without
numerically differentiating the Jacobian at every Hessian request. This is
an opt-in counterpart of SciPy's `NonlinearConstraint(hess=BFGS())` and
`NonlinearConstraint(hess=SR1())` strategies:

```fortran
blocks(1)%fun => nonlinear_values
blocks(1)%jac => nonlinear_jacobian
blocks(1)%lower = [1.0_dp]
blocks(1)%hess_update = 'sr1'
r = trust_constr(objective, x0, jac=gradient, hess=objective_hessian, &
    constraint_objects=blocks)
```

Each scalar constraint output has its own dense BFGS or SR1 approximation;
its Hessian is weighted by the **current** Lagrange multiplier when forming
the Lagrangian Hessian. This avoids false curvature when the multiplier
changes. SR1 can retain negative curvature; BFGS uses damped updates to
preserve positive definiteness. The method accepts analytic or numerical
block Jacobians, can be combined with an exact Hessian on another block,
and supports equalities, inequalities, bounds, and objective Hessian updates.
`constr_nhev` remains zero for the updated block, as for SciPy Hessian-update
strategies. Actual Jacobian/function evaluations are counted as usual.
`%hess_update` conflicts with `%hess`, `%hessp`, or `%hess_fd_scheme`, and
cannot be assigned to a linear block. With no `%hess_update`, SciFort retains
its existing numerical-derivative fallback; its dense per-output update
matrices and iteration paths need not match SciPy's weighted-block updates.
The SciPy comparison and Fortran regression are in
`tools/generate_optimize_trust_constr_constraint_updates_reference.py` and
`test/test_optimize_trust_constr_constraint_updates.f90`.

For constant-matrix SciPy `LinearConstraint` blocks, use the
`linear_constraint(matrix, lower=..., upper=..., keep_feasible=...)` factory.
It returns a `nonlinear_constraint` value so it can be mixed with nonlinear
blocks in one Fortran array. For example:

```fortran
type(nonlinear_constraint) :: blocks(2)
type(multivariate_minimize_result) :: r

blocks(1) = linear_constraint(reshape([1.0_dp, 1.0_dp], [1, 2]), &
    lower=[2.0_dp], upper=[2.0_dp])
blocks(2)%fun => inequality_values
blocks(2)%lower = [0.5_dp]
r = trust_constr(objective, [0.6_dp, 1.4_dp], constraint_objects=blocks)
```

Equivalently, assign `%matrix` directly to represent `matrix*x`. Its rows
must match the bound array length, columns must match `size(x0)`, and all
entries must be finite. Do not combine `%matrix` with `%fun`, `%jac`, `%hess`, `%hessp`,
or finite-difference settings. Matrix blocks have constant exact Jacobians,
zero constraint Hessians, and zero user callback counts (`constr_nfev`,
`constr_njev`, `constr_nhev`), matching SciPy. When the objective Hessian is
analytic, a collection of linear blocks uses the exact zero constraint
Hessian rather than a quasi-Newton approximation. Linear inequalities may
use `%keep_feasible`, just like nonlinear blocks.

Constraint values and Jacobians appear row-stacked in input block order in the
result; `constr_nfev`, `constr_njev`, and `constr_nhev` report per-block work
counts, followed by a zero-work variable-bounds block if supplied. Blocks may
mix analytic and numerical Jacobians. A block without `%jac` uses forward
(`'2-point'`) differences by default; `%fd_scheme = '3-point'` selects a
second-order central stencil and one-sided second-order probes near retained
variable bounds. The optional `%finite_diff_rel_step` array has length
`size(x0)` and supplies signed relative steps per input variable (for example,
`[-1.0e-4_dp, 2.0e-4_dp]`). As in SciPy, for two-point differences an explicit
relative step `r(i)` requests a displacement `r(i)*sign(x(i))*abs(x(i))`.
This allows negative steps to reverse the probe direction; zero or
unrepresentably small steps fall back to the automatic scale. Three-point
differences use the magnitude of the requested step and centered or one-sided
second-order stencils. At retained bounds, probes reverse direction or shorten
as necessary. Dense and `%finite_diff_jac_sparsity`-grouped Jacobians follow
the same rules. Reference cases and evaluation counts are reproduced in
`tools/generate_optimize_trust_constr_constraint_rel_step_reference.py` and
`test/test_optimize_trust_constr_constraint_rel_step.f90`.
Without this array, the two-point default is
`sqrt(epsilon(1.0_dp))` and the three-point default is
`epsilon(1.0_dp)**(1.0_dp/3.0_dp)`; `eps=` overrides those default steps.

For a nonlinear block with an analytic `%jac`, specify
`%hess_fd_scheme = '2-point'` or `'3-point'` to form its
multiplier-weighted constraint Hessian numerically by differentiating that
block's Jacobian. This is independent of `%fd_scheme`, which controls
numerical **Jacobians**. The Hessian calculation uses the same optional
`%finite_diff_rel_step` vector (two-point steps retain their signs), and
retained variable bounds select one-sided stencils as necessary. The dense
Hessian is symmetrized before use. The option requires `%jac` and cannot be
combined with `%hess` or `%hessp`; matrix (linear) blocks have zero Hessians
and cannot use it. Numerical Hessian work increases `constr_njev` through
analytic-Jacobian calls but does not increase `constr_nhev`, as in SciPy.
Without an explicit `%hess_fd_scheme`, the existing Hessian update/fallback
behavior is preserved. Compare against SciPy with
`tools/generate_optimize_trust_constr_constraint_hess_fd_reference.py` and
`test/test_optimize_trust_constr_constraint_hess_fd.f90`.

The trust-region constrained solver accepts an optional `callback=` in both
`trust_constr` and `minimize(method='trust-constr')`. This callback receives the
current physical decision variables (not interior-point slack variables) and
an iteration number starting at 1 for the initial state. A callback can request
early termination by setting `stop = .true.`; the result has `status = 3`,
`success = .false.`, and contains the final reported point and diagnostics.
The iteration number is cumulative across interior-point barrier stages,
including each stage's initial state. A callback returning `.false.` does
not change the numerical algorithm or its work counts. For example:

```fortran
r = minimize(objective, x0, method='trust-constr', &
    jac=gradient, callback=iteration_hook)

! An internal subroutine or module procedure with this interface:
subroutine iteration_hook(x, nit, stop)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: nit
    logical, intent(out) :: stop
    stop = nit >= 50
end subroutine iteration_hook
```

`callback=` is intentionally specific to trust-constr; other `minimize`
methods reject it instead of silently ignoring it. SciFort does not yet pass
SciPy's full `OptimizeResult` state into this callback. The regression is
`test/test_optimize_trust_constr_callback.f90`; the corresponding SciPy 1.17
comparison script is `tools/generate_optimize_trust_constr_callback_reference.py`.

For sparse nonlinear constraint Jacobians, allocate the logical
`%finite_diff_jac_sparsity(noutputs, size(x0))` pattern and set entries to
`.true.` wherever a constraint output can depend on an input variable. For
example, for independent output functions of three variables:

```fortran
blocks(1)%fun => three_independent_constraints
blocks(1)%lower = [1.0_dp, 4.0_dp, 9.0_dp]
blocks(1)%upper = blocks(1)%lower
allocate(blocks(1)%finite_diff_jac_sparsity(3, 3))
blocks(1)%finite_diff_jac_sparsity = .false.
blocks(1)%finite_diff_jac_sparsity(1, 1) = .true.
blocks(1)%finite_diff_jac_sparsity(2, 2) = .true.
blocks(1)%finite_diff_jac_sparsity(3, 3) = .true.
r = trust_constr(objective, x0, constraint_objects=blocks)
```

Columns with disjoint nonzero row patterns share a single forward-difference
probe or a pair of three-point probes. A zero column requires no probes. The
pattern must include every structurally nonzero derivative; omissions can
produce an incorrect Jacobian and solution. Masks are validated for shape,
ignored when an exact `%jac` callback is supplied, and rejected on linear
constraint blocks. This reduces function calls only: SciFort's `trust-constr`
still forms dense Jacobians and uses dense SQP projections, unlike SciPy's
sparse `LSMR`-based solve. It also works with exact constraint Hessian and
Hessian-product callbacks and per-block evaluation counts.
Numerical probes respect retained bounds. Barrier trials outside strict
variable bounds are rejected before user objective/constraint evaluations.
If every nonlinear block omits `%hess` and `%hessp`, the solver still uses
the existing BFGS constraint-Hessian approximation. When some blocks provide
exact Hessians, the others use the per-block finite-difference fallback just
described; an exact objective Hessian is not required because BFGS or SR1
objective curvature can be combined with exact constraint Hessians. Do not
combine `constraint_objects` with the legacy combined-constraint arguments.

Nonlinear least squares provides `least_squares`, legacy `leastsq`, and
`curve_fit`. `least_squares` defaults to dense TRF, with dense dogbox and
MINPACK-backed LM also available. TRF/dogbox support finite bounds,
bounds-aware two-point/three-point numerical Jacobians, complex-step Jacobians,
`jac_sparsity` grouped finite differences, or analytic Jacobians,
variable scales, evaluation limits, and SciPy's built-in `linear`, `soft_l1`, `huber`,
`cauchy`, and `arctan` losses with `f_scale`. Bounded `curve_fit` defaults to
TRF and can also use dogbox; its numerical model Jacobian accepts the same
`diff_scheme='2-point'`, `'3-point'`, or `'cs'` selection plus `diff_step`, with
`complex_model=` supplying the analytic complex continuation needed by complex-step
fits. Robust losses propagate through the trust-region fit and covariance uses the
robust-scaled Jacobian, matching SciPy's convention.
TRF and dogbox support `tr_solver='lsmr'`. TRF uses a matrix-free LSMR Gauss-Newton step with SciPy-style two-dimensional trust-region subspace solving; `tr_regularize=.false.` disables SciPy's default TRF regularization term. Dogbox uses LSMR for the scaled approximate Gauss-Newton step in its rectangular dogleg model. Analytic matrix-free Jacobians can be supplied with paired `jac_matvec=` and `jac_rmatvec=` callbacks for TRF; this automatically selects LSMR, avoids allocating the dense `m x n` Jacobian, and marks `least_squares_result%jac_is_operator=.true.`. As in SciPy, `x_scale='jac'` is unavailable for operator Jacobians. `jac_sparsity=` accepts a logical residual-by-parameter structural mask, batches nonconflicting numerical-difference columns, and automatically selects LSMR for TRF and dogbox; SciFort uses deterministic natural-column greedy coloring rather than SciPy's seeded permutation. Sparse stored Jacobians and parallel workers remain roadmap items. TRF and dogbox
also accept a Fortran-native `callback=` function invoked after each completed outer
iteration with the current parameters, true residuals, cost, iteration number, and
residual-evaluation count; returning `.true.` requests SciPy-style status `-2`
termination. As in SciPy, LM does not invoke the callback.
`x_scale` accepts either a positive real array or SciPy-style `x_scale='jac'`;
the latter initializes scales from Jacobian column norms and updates them
monotonically after accepted steps. Numerical Jacobians select `diff_scheme='2-point'` (default), `'3-point'`, or
`'cs'`. Three-point differences use central differences in the interior and
second-order one-sided differences near active bounds. Complex-step differences
require `complex_fun=`, an analytic complex continuation of the residual callback;
bounds do not alter the imaginary perturbation, matching SciPy's convention.

For `least_squares`, a custom robust loss can also be supplied with `loss_fun=`. The loss callback receives `z=(f/f_scale)**2` and returns vector values for `rho`, its first derivative, and its second derivative, matching SciPy's callable-loss contract.

The `scipy.integrate` layer is exposed by `scifort_integrate`. It provides
`trapezoid`, `cumulative_trapezoid`, `simpson`, `cumulative_simpson`, `romb`,
`fixed_quad`, `newton_cotes`, adaptive scalar `quad` with finite/infinite bounds,
break points, and QUADPACK weighted modes, vector `quad_vec`, nested `dblquad`,
`tplquad`, and `nquad`-style helpers, explicit-state `qmc_quad`, independent
double-exponential `tanhsinh`, adaptive multidimensional `cubature` with
Genz-Malik and product Gauss-Kronrod rules, and a `solve_ivp` layer with
adaptive RK23/RK45/DOP853, implicit Radau, variable-order BDF, and a
real-system LSODA path with automatic Adams/BDF switching; `t_eval`,
method-native dense output for RK45/DOP853/Radau/BDF/LSODA, Hermite interpolation
for RK23, optional Radau/BDF/LSODA dense Jacobian callbacks, constant dense or CSC
Jacobians for Radau/BDF, state-dependent CSC callbacks for Radau/BDF,
dependency-free sparse LU for CSC Jacobians, direct CSC numerical Jacobians and
grouped finite differences from Jacobian sparsity, LSODA banded-Jacobian and
minimum-retry-step controls, scalar or componentwise-vector `rtol`/`atol`, and scalar or vector event detection with direction filters, Nth-occurrence termination,
and per-event result views. A legacy `odeint` facade over LSODA provides
time-major output on monotone grids (including repeated times), the traditional
`func(y,t)` callback order, a type-safe `odeint_tfirst` companion, Jacobian/band
controls, scalar or componentwise-vector tolerances, hard `tcrit` stop bounds,
legacy Adams/BDF order caps (`mxordn`/`mxords`), legacy-style step controls,
optional `full_output` LSODA diagnostic histories/workspace metadata, and the
legacy `ixpr`/`mxhnil`/`printmessg` reporting controls.
RK23, RK45, DOP853, and BDF are also
overloaded for complex-valued systems, with
complex dense output, optional dense or CSC complex BDF Jacobian callbacks, and
real-valued event functions over complex states.
Adaptive scalar integration uses the BSD-3-Clause modern QUADPACK package by
Jacob Williams as a pinned FPM dependency; SciFort does not vendor that dependency.
The RK23/RK45 tableaux are adapted with attribution from Jacob Williams' BSD-3-Clause
`rklib` implementation. DOP853 is adapted from Jacob Williams' permissively licensed
modern Fortran edition of the Hairer/Wanner DOP853 solver. Radau IIA stepping,
BDF/NDF stepping, dense output, and finite-difference Jacobian logic are adapted
from SciPy 1.17.0 under its BSD-3-Clause license. The LSODA compatibility layer combines SciFort orchestration and BDF code with an
Adams Nordsieck path adapted from the public-domain ODEPACK `DSTODA` and `DCFODE`
routines. LSODA accepts `lband`/`uband` for packed analytic Jacobians or grouped
banded finite differences and `min_step` as a retry floor. The nonstiff branch now
starts at Adams method order 1 and keeps the accepted scaled-derivative Nordsieck
state persistently through orders 1--12. Trial-step changes rescale derivative
columns by powers of the step ratio, endpoint prediction uses the Pascal-triangle
shift, functional iteration uses ODEPACK's roundoff/contraction tests, and accepted
steps apply the `EL` correction directly to the predicted Nordsieck array. Successful
step/order reconsideration follows the ODEPACK `IALTH`/`TESCO` lower/current/upper
order comparison, including the 10 percent change threshold and order-dependent
`SM1` stability limits. The retained Nordsieck polynomial also supplies native dense
output on nonstiff LSODA intervals instead of the previous cubic-Hermite fallback.
SciFort retains nonuniform endpoint-derivative history as a reconstruction/fallback
path and for compatibility diagnostics rather than rebuilding the active Adams state
after every accepted step.

Method switching uses ODEPACK's 20-successful-step hysteresis, the order-5
Adams-to-BDF eligibility guard, and the common finite-error ideal-step comparison
built from ODEPACK's Adams/BDF error constants. If accepted-step error information
is dominated by roundoff, the ODEPACK-style fallback selects BDF only when the Adams
proposal is already stability-limited. Rejected Adams error-test trials do not probe
Jacobians or switch methods; after four consecutive failures SciFort conservatively
restarts low-order history and cuts the retry step by ten. Corrector-convergence
failure retries at one quarter of the failed step. Accepted BDF steps return the
endpoint derivative implied by the converged BDF equation instead of performing a
redundant RHS call. Representative nonstiff and stiff LSODA regressions now match the targeted
SciPy/ODEPACK work counts and accepted-step/controller transitions documented in
the validation notes. The separate `odeint` facade now enforces the legacy per-output-call
`mxstep` limit (default 500) without restarting LSODA history and can return
per-output `hu`/`tcur`/`tolsf`/`tsw`/`nst`/`nfe`/`nje`/`nqu`/`mused` histories
plus `imxer`, `lenrw`, and `leniw`. The facade also supports `ixpr` method-switch
messages, `mxhnil` warning limits, and `printmessg` final-status reporting; these
controls do not change LSODA numerical work.


The dependency-free special-function layer provides regularized incomplete
gamma and beta functions (`gammainc`, `gammaincc`, `betainc`, `betaincc`) with
their inverses. `scifort_special` also provides error/inverse-error
functions, real/complex `gamma`/`rgamma`/`loggamma`, real `gammasgn`/`gammaln`,
`multigammaln`, real `beta`/`betaln`, generalized `binom`, `poch`, integer-order
`polygamma`, floating `factorial`/`factorial2`/`factorialk`, `comb`/`perm`,
`bernoulli`/`euler` number sequences and `stirling2`,
`wright_bessel`/`log_wright_bessel` for nonnegative real parameters,
Jahnke-Emden `lmbda`, weighted Bessel integral `besselpoly`, complex `erf_zeros`,
real/complex `spence`, real/complex `lambertw` plus real/complex `wrightomega`,
real-degree/integer-order Ferrers `lpmv`, integer-degree real/complex
`legendre_p`/`assoc_legendre_p`, their current SciPy-style `*_all` table APIs,
Legendre-Q sequences `lqn`/`lqmn`, normalized `sph_legendre_p`/`sph_harm_y`,
and derivative helpers,
real/complex hypergeometric `hyp0f1`/`hyp1f1`/`hyp2f1`, real `hyperu`,
parabolic-cylinder `pbdv`/`pbvv`/`pbwa` plus
`pbdv_seq`/`pbvv_seq`/`pbdn_seq`, real ellipsoidal harmonics
`ellip_harm`/`ellip_harm_2`/`ellip_normal`, real `struve`/`modstruve` and the
order-zero Struve integral helpers, order-zero Bessel integral pairs
`itj0y0`/`it2j0y0` and `iti0k0`/`it2i0k0`, order-zero Kelvin
`ber`/`bei`/`ker`/`kei` functions with derivatives and zero tables, integer-order
Mathieu `mathieu_a`/`mathieu_b`, Fourier coefficients, angular
`mathieu_cem`/`mathieu_sem`, and radial `mathieu_modcem1`/`mathieu_modcem2`/
`mathieu_modsem1`/`mathieu_modsem2`, prolate/oblate spheroidal characteristic
values `pro_cv`/`obl_cv` plus `pro_cv_seq`/`obl_cv_seq`, first-kind angular
`pro_ang1`/`obl_ang1`, first-kind radial `pro_rad1`/`obl_rad1`, and second-kind
radial `pro_rad2`/`obl_rad2` functions,
Riccati-Bessel
`riccati_jn`/`riccati_yn`,
`agm`, `owens_t`, `voigt_profile`,
`huber`/`pseudo_huber`, convenience helpers such as `cbrt`, `exp10`, `exp2`,
`powm1`, `sinc`, degree trigonometric functions, `softplus`, and `diric`,
SciPy-compatible distribution helpers `bdtr`/`nbdtr`/`pdtr`,
`chdtr`/`fdtr`/`gdtr`/`stdtr` with direct complements, quantile inverses, and
legacy parameter inverses such as `btdtria`/`btdtrib`, `chdtriv`, `gdtrib`, and
`stdtridf`, plus noncentral `chndtr`/`ncfdtr`/`nctdtr` families with variate
and parameter inverses, historical `smirnov`/`smirnovi`,
`kolmogorov`/`kolmogi`, and `tklmbda` wrappers, integer-degree
Legendre/Jacobi/Gegenbauer/Laguerre/
Hermite/Chebyshev evaluators and their shifted SciPy variants, plus matching
Gaussian-quadrature root/weight generators and the legacy `assoc_laguerre`/
`*_roots` compatibility names, Riemann/Hurwitz
`zeta` and `zetac`, real/complex
`digamma`/`psi`, normal-CDF transforms including `ndtri_exp`, logistic
transforms, entropy helpers, Box-Cox transforms, and stable real/complex
`logsumexp`, `softmax`, and `log_softmax` reductions for vectors through rank-15
arrays, including tuple-style dimension lists, keepdims companions, and
sign/phase-returning log-sum-exp. Ranks 7 through 15 use a shared flattened
reduction engine rather than rank-specific numerical kernels. Weighted `logsumexp`
reductions also extend through rank 15, with scalar/trailing-rank/singleton-dimension
broadcasting that can expand singleton input dimensions and real/complex weights.
It also provides common
order-zero/one Bessel functions, real-order `jv`/`jve`/`yv`/`yve` (with SciPy-compatible
`jn` alias), integer-order `yn`/`kn`, real-order `iv`/`ive`/`kv`/`kve`, real-axis
`jvp`/`yvp`/`ivp`/`kvp` derivatives of arbitrary nonnegative integer order,
real/complex `hankel1`/`hankel2` plus scaled `hankel1e`/`hankel2e` and
`h1vp`/`h2vp` derivatives,
positive real Bessel-zero tables `jn_zeros`/`jnp_zeros`/`yn_zeros`/`ynp_zeros`
and the combined `jnyn_zeros`,
real/complex Airy `airy`/`airye` functions, `itairy` integrals, real
`ai_zeros`/`bi_zeros`, Carlson symmetric elliptic `elliprc`/`elliprd`/`elliprf`/
`elliprg`/`elliprj`, Legendre elliptic `ellipk`/`ellipkm1`/`ellipe`/`ellipkinc`/
`ellipeinc`, Jacobi elliptic `ellipj`, real/complex exponential integrals
`exp1`/`expi`, real `expn`, real/complex sine/cosine integrals `sici`, hyperbolic
integrals `shichi`, and real spherical
`spherical_jn`/`spherical_yn`/
`spherical_in`/`spherical_kn` functions with optional first derivatives. The
special-function layer exposes scaled modified-Bessel `i0e` and `i1e` kernels
shared by the Rice and von Mises
families, plus the public real Hurwitz-zeta kernel and its exponent derivative used by
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
layer adds `directional_stats` mean-direction/mean-resultant summaries for rank-2
observation matrices and rank-3 batches, unit-sphere and von Mises-Fisher
sampling/densities, Haar orthogonal/special-orthogonal/unitary generators,
prescribed-spectrum random correlation matrices, and fixed-margin random contingency
tables. `directional_stats` uses one-based observation axes and normalizes each
observation by default; `normalize=.false.` preserves input magnitudes.

The binned-statistics layer provides `binned_statistic`, `binned_statistic_2d`, and
`binned_statistic_dd` with SciPy-style `mean`, `std`, `median`, `count`, `sum`, `min`,
and `max` statistics, equal-width or explicit edges, rightmost-edge inclusion,
outlier bin numbers, and pure custom-statistic callbacks. Two-dimensional results
return both SciPy-compatible C-order flattened and expanded bin numbers. Because a
Fortran derived-type component cannot have runtime rank, the arbitrary-dimensional
result stores the statistic in C-order flattened form together with
`statistic_shape`; padded edge columns are accompanied by `edge_counts`. The
current native surface accepts one values vector per call; multiple value sets can be
processed by repeated calls while reusing the same edge specification.

The first-class distribution-infrastructure layer provides `make_distribution` and
`truncate`. `make_distribution` adapts pure scalar PDF/CDF/PPF callbacks, a support
interval, frozen parameter data, endpoint-inclusion flags, and optional stable
log/tail callbacks into a `continuous_distribution` object. The object exposes
PDF/log-PDF, CDF/log-CDF, survival/log-survival, inverse lower/upper tails, support,
median, and explicit-`rng_state` inverse-transform sampling. `truncate` wraps such an
object on the intersection of its support with user bounds, using stable normalized
CDF/SF formulas and preserving endpoint inclusivity. Unlike SciPy's Python
`make_distribution`, the native Fortran adapter does not inspect `rv_continuous`
classes or synthesize a new parameterized class: parameters are frozen when the
adapter is constructed, and generic moments/entropy/fitting are not yet part of this
small first-class interface.

The statistics layer also includes `gaussian_kde`, with weighted or
unweighted data, Scott/Silverman/constant bandwidth factors, PDF/log-PDF evaluation,
Gaussian and box integration, KDE-product integration, marginals, inverse covariance,
and explicit-state resampling. The `scipy.stats.qmc` layer provides sample scaling
and integer mapping, centered/wrap-around/mixture/L2-star discrepancies, geometric
minimum-distance/MST discrepancy, iterative discrepancy updates, Van der Corput,
Halton, and Joe-Kuo Sobol sequences, ordinary or strength-two orthogonal-array Latin
hypercubes, Poisson-disk sampling, centered-discrepancy `random-cd` and
centroidal-Voronoi `lloyd` post-processing, and QMC multinomial and
multivariate-normal samplers. Sobol supports up to 21,201
dimensions, 64 direction-number bits, reset, fast-forward, base-two blocks, and
SciFort-state LMS+digital-shift scrambling. QMC samples follow SciPy's `(n,d)`
layout, unlike the column-oriented random-variate arrays used by multivariate
distributions. `sobol_indices` adds Saltelli 2010 first- and total-order global
sensitivity indices. The native API accepts either precomputed `f_A`, `f_B`, and
`f_AB` arrays or a vectorized model callback plus independent inverse-CDF
callbacks; the latter uses the existing 64-bit scrambled Sobol engine and an
explicit `rng_state`. SciPy's callable custom estimator and `SobolResult.bootstrap`
conveniences are not yet exposed; scalar-output results retain the explicit
Fortran shape `(1,d)` rather than being squeezed to rank one.

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
correlation, and SciPy-style `rankdata` tie methods. It also provides vector
forms of SciPy's weighted geometric/harmonic/power means, standard error,
coefficient of variation, IQR, median absolute deviation, moment, skewness,
and Fisher/Pearson kurtosis helpers. All ten helpers support native Fortran
arrays through rank 15, the standard Fortran array-rank ceiling. Ranks 5 through
15 use shared assumed-rank reduction engines plus thin fixed-result-rank
wrappers, so the statistical formulas remain centralized in the vector kernels.
A one-based Fortran `axis` reduces the corresponding dimension (default 1 where
SciPy defaults to Python axis 0), while separate `*_keepdims` companions retain
the reduced dimension because a Fortran function cannot vary result rank at
runtime. Rank-2 weighted means accept rank-1 or rank-2 trailing-broadcast
weights. Rank-3 `gmean`, `hmean`, and `pmean` accept rank-1/rank-2/rank-3
trailing-broadcast weights. Rank-4 through rank-15 mean-family calls additionally
accept scalar or trailing-rank weights up to the input rank. Both forms can
expand singleton data dimensions. For `iqr`, omitting `axis` on rank-2 through
rank-15 arrays reduces all elements, matching SciPy's `axis=None` default. The
ten SciPy-style helpers accept native `nan_policy='propagate'`,
`nan_policy='omit'`, or `nan_policy='raise'` through rank 15, including weighted
mean reductions. Weighted means omit an observation when either its value or its
corresponding broadcast weight is NaN. With `'raise'`, any such NaN terminates
with Fortran 2018 `error stop`, the native analogue of SciPy's exception policy.
Infinities are not omitted. Unsupported policy names return NaN. `iqr` and
`median_abs_deviation` accept either a numeric `scale` divisor or the
case-insensitive named scale `scale='normal'` through rank 15; this uses the
standard-normal 75th-percentile conversion used by SciPy. SciFort also
provides `zscore`, `gzscore`, `zmap`, and `gstd` through rank 15. These use
one-based axes, `ddof`, and `nan_policy='propagate'` / `'omit'` / `'raise'`;
`gstd_keepdims`
retains the reduced dimension. The explicit `zscore_all`, `gzscore_all`,
`zmap_all`, and `gstd_all` companions are the Fortran analogue of SciPy
`axis=None`. `zmap` takes its location and scale from a separate comparison
array, permits a different extent along the reduced axis, and supports singleton
broadcasting on non-reduced dimensions. `gstd_keepdims(x,dims)` accepts a unique
list of one-based dimensions as the fixed-rank Fortran representation of SciPy's
tuple-axis reduction. `gstd` follows SciPy on its documented finite, strictly
positive input domain.

SciFort also provides the value-limit trimmed-statistics family `tmean`,
`tvar`, `tstd`, `tsem`, `tmin`, and `tmax` through rank 15. Optional
`lower_limit` / `upper_limit` arguments represent SciPy's possibly open
`limits` tuple, and `inclusive` controls whether exact boundary values are
retained. The variance, standard-deviation, and standard-error routines
accept `ddof`; all six accept the same three NaN policies as the other
descriptive helpers. `tmean` without `axis` reduces the whole array, matching
SciPy's `axis=None` default; its axis-specific overload takes an explicit
one-based `axis`. The other routines default to axis 1. Explicit `*_all` and
fixed-rank `*_keepdims` companions are available for all six. Empty trimmed
slices return NaN, which keeps array reductions composable in native Fortran.

The proportion-based trimming helpers `trim_mean`, `trimboth`, and `trim1` are
also available through rank 15. `trim_mean` defaults to one-based axis 1 and
supports `nan_policy='propagate'`, `'omit'`, or `'raise'`; `trim_mean_all`
represents SciPy `axis=None`, and `trim_mean_keepdims` preserves the reduced
axis as extent one. `trimboth` removes `floor(proportiontocut*n)` observations
from each value-sorted tail, while `trim1` removes that many from either the
`'left'` or `'right'` tail. Their `trimboth_all` and `trim1_all` companions
flatten before trimming, matching SciPy's `axis=None` behavior. The ordering of
returned trimmed observations is intentionally unspecified, as in SciPy.

SciFort also provides `kstat` (orders 1 through 4) and `kstatvar` (orders 1
and 2) through rank 15. With no `axis` argument, `kstat` and `kstatvar` reduce
the whole array, matching SciPy's `axis=None` default; an explicit one-based
`axis` requests slice-wise reduction. The `kstat_all` / `kstatvar_all` names
make whole-array intent explicit, while `kstat_keepdims` / `kstatvar_keepdims`
preserve an explicitly reduced axis. Both families accept
`nan_policy='propagate'`, `'omit'`, or `'raise'`. The implementation evaluates
algebraically equivalent centered-moment forms of the k-statistics to reduce
cancellation when observations share a large common offset.

Circular descriptive statistics `circmean`, `circvar`, and `circstd` are available
through rank 15. The default angular interval is `[0, 2*pi]`; `high` and `low`
can define any nonzero period, including a negative period for compatibility with
SciPy's scaling convention. Calls without `axis` reduce the whole array, while
explicit one-based axis overloads and fixed-rank `*_keepdims` companions support
slice-wise reductions. All three functions accept the standard three NaN policies,
and `circstd(..., normalize=.true.)` returns the unscaled angular standard
deviation. Explicit `circmean_all`, `circvar_all`, and `circstd_all` aliases make
whole-array intent clear.

Bayesian mean/variance/standard-deviation inference is available through
`mvsdist` and `bayes_mvs`. `mvsdist` returns three frozen `mvs_posterior` objects
for the population mean, variance, and standard deviation; each object provides
`pdf`, `logpdf`, `cdf`, `logcdf`, `sf`, `logsf`, `ppf`, `isf`, `mean`, `var`,
`std`, and `interval` methods. As in SciPy 1.17, samples of at most 1000
observations use Student-t, inverse-gamma, and generalized-gamma posteriors under
Jeffreys' prior, while larger samples use Gaussian approximations. `bayes_mvs`
returns posterior means and central probability intervals, defaults to `alpha=0.9`,
and both entry points accept assumed-rank real arrays that are flattened before
analysis. Optional integer status outputs report invalid sample size or interval
probability without relying on exceptions.

Frequency and percentile summaries are also available. `percentileofscore(a, score)`
accepts a one-dimensional reference sample and scalar or rank-1 through rank-15 score
arrays, with SciPy's `rank`, `weak`, `strict`, and `mean` tie conventions and all three
NaN policies. `scoreatpercentile` accepts rank-1 through rank-15 data, scalar or vector
percentiles, whole-array reduction by default, one-based axis reductions, optional
closed value limits, and `fraction`, `lower`, or `higher` interpolation; vector
percentiles form the leading result dimension. `cumfreq` and `relfreq` flatten input
through rank 15 and return `cumfreq_result` / `relfreq_result` records containing the
histogram values, lower limit, bin size, and extrapoint count. Optional histogram
weights and explicit real limits are supported.

Ranking helpers now match the multidimensional SciPy interface more closely. `rankdata`
supports ranks 1 through 15, whole-array Fortran-order flattening by default, explicit
one-based axis ranking, all five SciPy tie methods (`average`, `min`, `max`, `dense`,
and stable `ordinal`), and `nan_policy='propagate'` / `'omit'` / `'raise'`.
`tiecorrect` computes the standard Mann-Whitney/Kruskal-Wallis tie correction factor
from a one-dimensional rank vector.

Core multi-output summaries `mode` and `describe` are also available through rank 15. Because Python returns named result objects whose component rank depends on the runtime axis, SciFort exposes these as Fortran subroutines with separate allocatable `intent(out)` arguments. `mode` defaults to one-based axis 1, supports `mode_all` and `mode_keepdims`, and implements all three NaN policies; under `'propagate'`, NaN is counted as an ordinary modal value, matching SciPy 1.17. `describe` defaults to one-based axis 1 and reports observation count, minimum, maximum, mean, variance, skewness, and Fisher kurtosis; `describe_all` flattens the input. Arbitrary integer `ddof`, bias correction, and the three NaN policies are supported.

Information-theory summaries `entropy` and `differential_entropy` are available
through rank 15. `entropy` normalizes each probability slice before computing
Shannon entropy or, when `qk` is present, Kullback-Leibler relative entropy.
Lower-rank `qk` arrays use trailing NumPy-style broadcasting when their rank does
not exceed the rank of `pk`; same-rank singleton broadcasting is also supported.
The default reduction is one-based axis 1, with `entropy_all` representing
SciPy `axis=None` and `entropy_keepdims` retaining the reduced axis. Optional
`base` selects information units, and `nan_policy='propagate'`, `'omit'`, or
`'raise'` is supported; paired omission drops an event when either `pk` or `qk`
is NaN. `differential_entropy` implements the Vasicek, Van Es, Ebrahimi, and
Correa spacing estimators plus SciPy's `auto` selection, optional window length
and logarithm base, one-based axis reductions, `differential_entropy_all`,
`differential_entropy_keepdims`, and the same three NaN policies.

Additional one-sample summaries include `expectile` and `sigmaclip`. `expectile`
uses all elements of a rank-0 through rank-15 input and accepts an optional level
`alpha` in `[0,1]` plus NumPy-style trailing/singleton-broadcast case weights;
`alpha=0.5` is the weighted arithmetic mean. Invalid levels, empty input, invalid
weight broadcasting, NaNs, or zero total weight return NaN. `sigmaclip` flattens
its input and iteratively removes observations outside `mean-low*std` and
`mean+high*std`, preserving the order of retained values and returning the final
lower/upper thresholds alongside the clipped vector.

Sample L-moments are provided by `lmoment`. A call without `order` returns the
first four L-moments, while a scalar `order` returns the corresponding reduced
array and a rank-1 `order` vector forms the leading result dimension. The default
reduction is one-based axis 1; `lmoment_all` represents SciPy `axis=None`, and
`lmoment_keepdims` retains an explicitly reduced axis. `sorted`, `standardize`,
and `nan_policy='propagate'` / `'omit'` / `'raise'` follow SciPy 1.17 semantics.
Orders 3 and higher are L-moment ratios by default and raw L-moments when
`standardize=.false.`. Axis reductions and scalar-order keepdims support rank-15
input. Default/vector-order keepdims support rank-14 input because a rank-15
input plus the order dimension would require a rank-16 Fortran array.

Power-transform statistics are available through `boxcox`, `boxcox_llf`,
`boxcox_normmax`, `yeojohnson`, `yeojohnson_llf`, and `yeojohnson_normmax`.
The transform functions are elemental; log-likelihood reductions support one-based
axes through rank 15, explicit `*_all` companions for SciPy `axis=None`, fixed-rank
`*_keepdims`, and `nan_policy='propagate'` / `'omit'` / `'raise'`.
`boxcox_normmax` supports the SciPy `pearsonr` and `mle` criteria, with
`boxcox_normmax_all` returning both estimates. `boxcox_fit` and `yeojohnson_fit`
provide the SciPy omitted-lambda workflow as explicit Fortran subroutines returning
the transformed vector and fitted parameter.

Probability-plot diagnostics are available through `probplot`, `ppcc_plot`,
`ppcc_max`, `boxcox_normplot`, and `yeojohnson_normplot`. `probplot` returns
Filliben theoretical quantiles and ordered responses through allocatable output
vectors and optionally reports the least-squares slope, intercept, and correlation.
Built-in distribution-name dispatch covers `norm`, `tukeylambda`, `weibull_min`,
`t`, `loggamma`, `gamma`, `logistic`, `expon`, `laplace`, and `cauchy`; a
`probability_plot_ppf` callback allows arbitrary user-defined quantile functions.
`ppcc_plot` and `ppcc_max` support one-shape families, including SciPy's default
Tukey-lambda family, while the Box-Cox and Yeo-Johnson normality-plot helpers reuse
the power-transform layer and normal probability-plot correlation calculation.
Plot rendering itself is intentionally left to the caller.

Empirical-distribution distances are available through `wasserstein_distance`,
`wasserstein_distance_nd`, and `energy_distance`. The one-dimensional routines accept
independent sample/support vectors with optional nonnegative weights; weights are
normalized within each sample. The 1-D implementation integrates the absolute
(Wasserstein-1) or squared (energy) difference between the empirical CDFs.
`wasserstein_distance_nd` accepts either vectors (delegating to the 1-D routine) or
rank-2 support arrays whose rows are observations and whose columns are coordinates.
Its multidimensional path solves the balanced Euclidean optimal-transport problem with
a project-local residual-network min-cost-flow solver, so no external LP package is
required. Empty samples, weight lengths, nonnegative weights, positive finite weight
sums, dimensionality, and nonfinite observations are validated explicitly.

Empirical-distribution estimation, survival comparison, and exact quantile inference
are available through `ecdf`, `logrank`, `quantile_test`, and `order_statistic`.
`ecdf(sample)` returns CDF/SF objects with unique quantiles, probabilities, elemental
`evaluate` methods, and Greenwood or log-log (exponential-Greenwood) confidence
intervals. Right-censored data can be supplied either as `ecdf(uncensored, right)` or
through the public `censored_data` type and `ecdf(censored_data)`, both using the
Kaplan-Meier product-limit estimator. `logrank` compares two survival samples with the
signed standard-normal logrank statistic and supports `two-sided`, `less`, and
`greater` alternatives; uncensored vectors and `censored_data` objects can be mixed.
`quantile_test` returns the exact binomial T1/T2 statistic and p-value plus a type-bound
nonparametric confidence interval. `order_statistic` constructs a small callback-backed
distribution object whose PDF/CDF/SF/PPF/ISF use the beta law of continuous order
statistics; the callback parameter vector lets it wrap SciFort or user-defined
continuous families.

`scifort_hypothesis`
provides t-tests, including elemental `ttest_ind_from_stats` for independent
two-sample tests from means, corrected sample standard deviations, and sample
counts, Pearson/Spearman correlation tests, and asymptotic Mann-Whitney U.
`ttest_ind_from_stats` supports pooled or Welch standard errors and all three
SciPy alternatives; negative standard deviations have the same result as their
absolute values, and the SciFort result additionally exposes the computed degrees
of freedom. `scifort_hypothesis_extended` adds Wilcoxon rank-sum and
signed-rank tests, Kruskal-Wallis and Friedman tests, Cressie-Read power
divergence and chi-square tests, 2-by-2 Fisher exact tests, contingency-table
chi-square, classical and Welch one-way ANOVA, and Bartlett, Levene, and
Fligner-Killeen variance tests. It also provides `obrientransform`: a single
vector returns its transformed observations directly, while a `sample_group`
array returns transformed groups of the same (possibly unequal) lengths for
subsequent `f_oneway` testing. Samples shorter than three observations or with
nonfinite values return NaNs instead of triggering floating-point exceptions.
`scifort_nonparametric_extended` adds the
Ansari-Bradley and Mood scale tests, Baumgartner-Weiss-Schindler two-sample
test, Epps-Singleton two-sample test, k-sample Anderson-Darling test, and Mood
median test. `scifort_association_extended`
adds Kendall tau-b/tau-c, weighted Kendall tau, the newer `spearmanrho` test,
Somers D ordinal association, Chatterjee xi correlation, point-biserial correlation,
ordinary linear regression,
Theil-Sen and Siegel robust slopes, the Brunner-Munzel test, and Page's
ordered-alternative trend test. Somers D accepts paired rankings or
integer-valued contingency tables and returns the table used for the asymmetric
`D(Y|X)` calculation. Weighted Kendall tau uses SciPy's additive hyperbolic
weights by default; `rank_mode='lexicographic'` corresponds to SciPy `rank=None`,
`rank_mode='index'` to `rank=False`, and integer or real rank arrays can be passed
as the third argument. Custom rank-to-weight procedures are accepted with `weigher=`
and multiplicative exchange weights with `additive=.false.`. `spearmanrho` provides the
paired-vector SciPy 1.17 statistic with asymptotic inference plus exact/randomized pairing
permutations and Monte Carlo inference; randomized methods use explicit `rng_state`.
Chatterjee xi supports
the tie-aware asymptotic null variance, the continuous-response shortcut, and
exact-or-randomized pairing permutations with explicit `rng_state` for randomized
inference. `multiscale_graphcorr` adds the SciPy 1.17 MGC independence test with
Euclidean or callback-defined pairwise distances, precomputed distance matrices,
automatic/forced unpaired two-sample reduction, the full local-correlation map and
optimal scale, and plus-one permutation inference using explicit `rng_state`.
Parallel `workers` execution is intentionally not exposed; permutations run in the
calling process. `scifort_contingency_meta` adds exact binomial tests,
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
families. It also provides SciPy-style Monte Carlo `goodness_of_fit` testing
for the `norm`, `expon`, and `uniform` location-scale families with Anderson-Darling,
Kolmogorov-Smirnov, Cramer-von Mises, and Filliben statistics. The native API uses
explicit `rng_state`; `known_loc`/`known_scale` stay fixed in every Monte Carlo
refit, while `fit_loc`/`fit_scale` define the null distribution but are refit for
the statistic just like SciPy's `fit_params`. `guess_loc`/`guess_scale` are accepted
for interface parity but do not affect these closed-form fits. `scifort_resampling`
adds explicit-state percentile/basic/BCa
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

### Multidimensional spline interpolation (ranks 1-15)

`scifort_signal` now exports `spline_filter_nd`, `spline_filter1d`,
`cspline_nd`, `qspline_nd` and `map_coordinates_nd` with Fortran 2018
assumed-rank inputs. Real/complex orders 0-5 are supported; coordinates
and axes are one-based. See [ND spline guide](docs/signal_spline_nd.md). The same assumed-rank core now supports all eight
SciPy `ndimage` boundary modes, plus `affine_transform_nd`, `shift_nd`, `zoom_nd` and `rotate_nd`
for rank-1 through rank-15 arrays. See
[Geometric transforms and boundary modes](docs/signal_spline_transforms.md)
for geometric examples and padded high-order boundary handling.
The new [`scifort_ndimage` filtering and morphology guide](docs/ndimage.md)
covers Gaussian smoothing and derivatives, full-dimensional convolution
and correlation, median/rank/percentile filters, Sobel/Prewitt/Laplacian,
exact Euclidean distance transforms, binary morphology, geodesic propagation,
hole filling, connected-component labeling, complex one-axis convolution,
custom callback neighborhood/label reductions, and anisotropic or directed
weighted chamfer distance transforms.

### Scattered-data interpolation (2026-10-31)

The `scifort_interpolate` facade now exports native 2D–5D Delaunay/barycentric
`linear_nd_interpolator`, dimension-agnostic `nearest_nd_interpolator`,
and `griddata` for supported methods. Real/complex multichannel data,
rescaling, and linear fill values are supported. See the
[scattered interpolation guide](docs/interpolate_scattered.md).

### Additional interpolation algorithms (2026-10-31)

The `scifort_interpolate` facade now also exports SciPy-style `RBFInterpolator`
for scattered points in multiple dimensions and the legacy `interp1d` API
with nine interpolation kinds. Both support real/complex multichannel output
without new dependencies. See
[the RBF/interp1d guide](docs/interpolate_rbf_interp1d.md).

### FITPACK-compatible spline tuples

`scifort_interpolate` provides dependency-free `splrep` (`s=0` or
`task=-1`), `splev`, `splint`, `splder`, `splantider`, and `spalde`.
The same module now provides nonperiodic parametric `splprep` (`s=0` or
fixed-knot `task=-1`) with automatic chord-length parameterization and
multichannel `splev` evaluation.
Full adaptive-knot residual-target smoothing is available as an optional
`perazz/fitpack` integration project (not a mandatory dependency). See
[the FITPACK interoperability guide](docs/interpolate_fitpack_tck_backend.md).

Rectangular-grid **exact** tensor-product interpolation is available through
`rect_bivariate_spline(x,y,z,kx,ky,status)` (SciPy `RectBivariateSpline(s=0)`),
with degrees 1–5 and **exact rectangle integration** via
`tck%integral(xa,xb,ya,yb)` or `bispline_integral(...)`; no FITPACK dependency.
The same TCK now supports a reusable mixed-partial derivative spline through
`%partial_derivative(dx,dy)` (SciPy `BivariateSpline.partial_derivative`),
and exact x/y **knot refinement without refitting** through
`%insert_knot(axis,x,m)`. Functional aliases are `bispline_partial_derivative`
and `bispline_insert_knot`. See [bivariate spline calculus and knot refinement](docs/interpolate_bivariate_calculus.md) and
[rectangular bivariate splines](docs/interpolate_rect_bivariate.md).

Bivariate tensor-product `fitpack_bivariate_tck_type`, `bisplev`, and
`bisplev_points` are available without external dependencies. Automatic-knot
bivariate smoothing and fixed-knot surface least squares are in the optional
FITPACK 2.0.0 adapter. See
[the bivariate FITPACK guide](docs/interpolate_fitpack_bivariate.md).

## 2026-11-10 — Mixed-integer and constrained differential evolution

`scifort_optimize_differential_evolution` now accepts `integrality`, bound-constrained
linear matrices, nonlinear constraint callbacks, and an optional local polishing
callback. `scifort_optimize_differential_evolution_polish` bridges to the existing
L-BFGS-B implementation for unconstrained continuous/mixed-integer polishing.
Native and independent SciPy comparison tests are documented in
`docs/optimize_differential_evolution.md` and
`validation/optimize_differential_evolution_constraints_2026-11-10.md`.
