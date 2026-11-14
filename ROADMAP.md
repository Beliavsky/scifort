## Spherical spline interfaces (2026-11-13)

- [x] Add native fixed-knot `lsq_sphere_bivariate_spline` with weighted fitting,
  periodic and polar coefficient constraints, derivatives, and rank diagnostics.
- [x] Validate against SciPy and analytic spherical data with focused tests.
- [x] Add native rectangular-grid spherical interpolation (`s=0`) with pole options.
- [x] Add native adaptive rectangular-grid spherical smoothing (`s>0`).
- [x] Add adaptive scattered-data spherical smoothing (`SmoothSphereBivariateSpline`).
- [ ] Extend the dense solver limits and audit rank-deficient FITPACK parity.

## Complex Bessel Y zero coverage (2026-11-13)

- [x] Add complex-output overloads of `y0_zeros`, `y1_zeros`, and `y1p_zeros`
  for up to 1000 upper-half-plane roots, retaining existing real APIs.
- [x] Validate the SciPy root sequence and actual derivative/value companions
  against multiprecision references; document SciPy discrepancies explicitly.

## Sparse weighted matching compatibility (2026-11-13)

- [x] Expose the existing sparse weighted assignment solver under the
  SciPy-compatible `min_weight_full_bipartite_matching` name, retaining the
  original public API and sparse storage.
- [x] Validate both entry points against exhaustive small matchings, dense
  assignment, and randomized SciPy comparisons; correct the stale gap entry.

## Standalone FFT coverage (2026-11-12)

The user's request to work on FFT gaps approves this roadmap expansion beyond
probability distributions. Prioritize native Fortran APIs; defer C ABI work.

- [x] Expose 1D complex, real, and Hermitian FFT pairs with SciPy normalization.
- [x] Add frequency bins, FFT shifts, and FFT-backed DCT/DST types I-IV.
- [x] Add multidimensional transforms and axis selection, including real and
  Hermitian pairs and multidimensional DCT/DST.
- [x] Add fast Hankel transforms (`fht`, `ifht`, `fhtoffset`) and axis-selectable batches.
- [x] Add multidimensional shifts and two-dimensional convenience names.
- [x] Add backend-specific fast-length helpers and audit batching support.
- [ ] Add worker controls, public multidimensional plans, and backend option parity.

## Censored distribution fitting (2026-11-12)

- [x] Extend the existing censored-data type to left and interval censoring.
- [x] Add stable normal, exponential, and Weibull censored likelihoods and
  bounded fitting overloads, with independent high-precision references.
- [x] Add an analytic fixed-location right-censored exponential MLE with
  explicit status reporting when no finite identifiable estimate is available.
- [x] Add stable lognormal censored likelihoods and fixed-location bounded fits.
- [x] Add gamma censored likelihoods and fixed-location bounded fits.
- [x] Add beta censored likelihoods and fixed-support bounded fits.
- [x] Add stable logistic censored likelihoods and bounded location/scale fits.
- [x] Add stable censored likelihoods and bounded location/scale fits for both Gumbel families.
- [ ] Extend censored fitting to additional continuous families.
- [ ] Design observation-weight semantics as a separate SciFort extension.

## Mixed-integer linear programming (2026-11-10)

- [x] Add `scipy.optimize.milp`-style dense branch-and-bound with continuous,
  integer, semi-continuous, and semi-integer variable domains.
- [x] Reuse the existing two-phase simplex LP solver; add two-sided row bounds,
  feasibility verification, incumbent/dual-gap diagnostics, and node limits.
- [x] Add native tests, independent randomized SciPy comparisons, small
  exhaustive integer enumeration, a runnable example, and provenance.
- [x] Correct a large-RHS simplex Phase I ratio-tie bug discovered by MILP
  testing; retain all 1,450 independent LP comparison cases.
- [ ] Test the complete FPM integration against all real dependencies.
- [ ] Evaluate stronger MILP backends (HiGHS bindings) for sparse or large
  problems, presolve, cuts, heuristics, and proof of integer unboundedness.
- [ ] Add a warm-start feasible incumbent and variable bound propagation.

## Quadratic assignment (2026-11-10)

- [x] Implement SciPy-style FAQ and 2-opt methods with deterministic starting
  points, one-based fixed partial matches, and optional warm starts.
- [x] Reuse the existing Hungarian solver in FAQ; avoid third-party code.
- [x] Add independent SciPy comparisons, exhaustive pairwise-local tests,
  example, documentation, and complete status/error reporting.
- [ ] Add explicit-state randomized starts and multistart strategies.
- [ ] Benchmark larger/sparse graphs and consider a limited-memory FAQ path.
- [ ] Run repository-wide `fpm test` with all external dependencies.

## SHGO independent nonlinear constraint objects (2026-11-10)

- [x] Add a public type with per-block bounds, value callbacks and optional
  analytic Jacobians; forward it through the SHGO core, public adapter, and
  `scifort_optimize` umbrella exports.
- [x] Support mixed analytic/finite-difference derivatives, mixed equality
  and inequality components, deterministic feasibility restoration, and
  combined linear/nonlinear constraints without altering the old API.
- [x] Add independent invalid-input tests, native regressions, canonical
  SLSQP Jacobian-layout smoke checks and SciPy comparisons.
- [ ] Install/use FPM to link and test against the actual SLSQP implementation.
- [ ] Improve restoration for rank-deficient, highly scaled and disconnected
  feasible sets; explore alternative SHGO sampling policies.

## SHGO analytic derivatives and restoration (2026-11-10)

- [x] Pass exact nonlinear constraint Jacobians to SLSQP with the correct
  equality/inequality row order and upper-bound derivative signs.
- [x] Implement an opt-in, bounded damped Gauss-Newton feasibility projection
  for difficult equality manifolds, with finite-difference fallback and
  objective/constraint/Jacobian evaluation accounting.
- [ ] Run full `fpm test` with the real SLSQP dependency; the adapter is
  currently validated using an interface/derivative-checking stand-in only.
- [ ] Explore reliable multistart restoration for disconnected/narrow feasible
  regions without imposing false global feasibility or optimality guarantees.

## SHGO constraint support (2026-11-10)

- Implemented optional vector-valued linear and nonlinear constraint blocks
  with bound pairs, equalities, and feasibility-first graph ranking.
- Constrained SHGO polishing defaults to SciFort SLSQP. Added independently
  validated feasible results, constraint counters, and infeasible status.
- Follow-up priorities: verify full FPM/SLSQP integration and develop robust
  restoration for disconnected regions. Independent nonlinear block objects
  with optional per-block Jacobians are now implemented.

## SHGO conforming 3D adaptive mesh status (2026-11-10)

- Completed: opt-in tetrahedralization of dyadic adaptive cubes using shared
  face-center samples, exactly triangulated coarse-fine faces, positive
  integer-lattice volume checks and complete face-incidence validation.
- Added 47 native tests and five SciPy objective-quality comparisons, tested
  under checked and optimized Fortran builds.
- Remaining: conforming 4D+ complexes, rigorous homology persistence,
  improved constrained restoration and FPM-backed full integration testing.

## SHGO conforming adaptive triangulation status (2026-11-10)

- Completed: optional conforming 2D triangulation through sampled cell-center
  fans and nested edge subdivisions; exact area/edge validation; native tests.
- Remaining: conforming 4D+ complexes, topological persistence beyond graph
  minima, improved constrained restoration and
  adaptive basin selection, and full FPM integration validation.

## SHGO homology-growth stopping status (2026-11-10)

Implemented: optional conservative stabilization of the sampled graph's
local-minimum birth count, minimum-round guard, best-value improvement guard,
and connected-plateau handling on adaptive meshes. This is an approximate
H_0-birth-count proxy, **not** a SciPy-equivalent `minhgrd` implementation.
Remaining: more rigorous homology tracking, refined constrained sampling,
and full FPM-backed integration tests.

## SHGO adaptive extension status (2026-11-10)

Implemented: selective dyadic *cell* refinement with exact vertex caching,
Kuhn cell-edge graphs, deterministic ranking and strict per-cell budgets.
Remaining: globally conforming adaptive triangulation (hanging-node
resolution), rigorous homology-rank computation, richer constraint objects,
alternative sampling methods, and full FPM-backed local-polisher tests.

## 2026-11-10 — Simplicial global optimization (SHGO)

- [x] Add a deterministic, bounded simplex-graph global optimizer with
  dyadic vertex reuse and graph local-minimum (H_0 birth) detection.
- [x] Reuse SciFort's Nelder-Mead/Powell local solvers; report ordered local
  minima, strict evaluation budgets, callbacks, and status/error objects.
- [x] Add standalone native tests and independent SciPy comparison cases.
- [x] Add adaptive dyadic cell refinement and heuristic H_0-birth-count
  stabilization as optional stopping.
- [x] Add conforming 2D/3D adaptive meshes with exact geometry and
  shared-simplex incidence tests.
- [ ] Add rigorous persistent-homology analysis; the current stopping
  heuristic only observes graph-minimum births.
- [x] Add linear/nonlinear inequality and equality constraints with SLSQP
  local polishing, using differential-evolution-style constraint arguments.
- [x] Add independent nonlinear constraint objects with optional Jacobians.
- [ ] Add alternative Sobol/Halton sampling and additional constraint methods.
- [ ] Run complete `fpm test` with Git-provided numerical dependencies.

## 2026-11-09 — Hankel structured matrices

- [x] Add `scipy.linalg.hankel` equivalent, including complex and
  rectangular construction and default zero-padded last row.
- [x] Add FFT-based Hankel matrix products as a SciFort extension,
  reusing the Toeplitz implementation, with mixed real/complex operands.
- [x] Add native and SciPy comparison tests for both routines.
- [ ] Consider `scipy.linalg.convolution_matrix` as a distinct next
  matrix-construction algorithm.

## 2026-11-09 — Toeplitz linear algebra

- [x] Implement `scipy.linalg.toeplitz` with real/complex and rectangular
  outputs, and `scipy.linalg.solve_toeplitz` with O(n^2) general bordering.
- [x] Support complex Hermitian/default rows, arbitrary real or complex
  first rows, vector and multiple RHS, mixed types, and dense pivot fallback.
- [x] Validate against SciPy/NumPy in checked and optimized builds.
- [x] Add FFT-based `circulant` and `solve_circulant` (2026-11-09)
  without expanding the existing QZ/Delaunay algorithms.

## 2026-11-08 — Real-to-complex Schur conversion

- [x] Add `rsf2csf` using local complex unitary similarities on real Schur
  2x2 blocks, reusing precomputed real Schur vectors and preserving inputs.
- [x] Validate against SciPy for real, conjugate, clustered, tiny/large
  scale and large-order precomputed Schur factors.
- [ ] Investigate scalable generalized Schur/QZ algorithms and generalized
  nonsymmetric eigenproblems, without importing incompatibly licensed code.

## 2026-11-08 — Real and complex Schur decomposition

- [x] Expose `schur` for real and complex matrices with complex upper-triangular
  output, four standard eigenvalue sorting regions, optional custom selectors,
  and returned selected invariant-subspace dimension.
- [x] Add SciPy comparisons and a standalone native example.
- [x] Add real quasi-triangular Schur output for real matrices with pair-aware
  sorting and independent SciPy validation.
- [ ] Consider generalized Schur/QZ for ordinary nonsymmetric matrix pairs.
- [ ] Consider larger unblocked Schur QR limits, error estimates, and optional
  balancing without changing existing matrix-function interfaces.

## 2026-11-08 — General nonsymmetric matrix eigenproblems

- [x] Add ordinary real and complex `eig` and `eigvals` with optional
  left/right eigenvectors, Schur-based computation, and SciPy comparisons.
- [ ] Generalized nonsymmetric `(A,B)` eigenproblems and homogeneous
  eigenvalues; larger problems require a stronger QZ/Schur backend.
- [ ] Optional matrix balancing and refined ill-conditioned eigenvector solves.

## 2026-11-07 — Dense symmetric/Hermitian eigenproblems

- [x] Add real and complex dense `eigh` and `eigvalsh` with selected
  upper/lower triangle and optional one-based inclusive eigenvalue subsets.
- [x] Validate eigenvalues, eigenpair residuals, and orthogonality against
  SciPy for real/complex, repeated, clustered and scaled matrices.
- [x] Generalized symmetric/Hermitian definite eigenproblems using optional
  `b`, with type 1/2/3 Cholesky reductions and SciPy-reference tests.
- [x] Add `(low,high]` `subset_by_value` for real/complex standard and
  generalized type 1/2/3 with capacity buffers, returned selection counts,
  and native/SciPy-reference tests.
- [ ] Mixed real/complex matrix pairs and faster partial-spectrum solvers.

## 2026-11-07 — Symmetric tridiagonal and Hermitian band eigensolvers

- [x] Add `eigh_tridiagonal` / `eigvalsh_tridiagonal` (scaled
  implicit-shift QL) and `eig_banded` / `eigvals_banded` (dense-backed
  symmetric/Hermitian packed-band solver).
- [x] Add index- and value-selected tridiagonal eigenvalues with scaled
  Sturm-count bisection and O(n) memory, preserving the existing QL path.
- [x] Add selected tridiagonal eigenvectors via pivoted inverse iteration,
  reorthogonalization, O(n*k) output storage, and a full-QL fallback.
- [x] Add an opt-in full-spectrum symmetric/Hermitian band-to-tridiagonal
  Givens bulge-chasing reduction, with no dense input expansion when
  eigenvectors are not requested.
- [x] Add banded eigenpair selection via Sturm/inverse iteration.
- [ ] Optimize the real-arithmetic reduction and add blocked variants for
  larger bandwidths; consider automatic method selection after benchmarking.

## 2026-11-07 — Dense matrix balancing

- [x] Add real/complex `matrix_balance` with optional isolated-row/column
  permutation and power-of-two scaling, including separate transform factors.
- [x] Validate similarity invariants and numerical quality against SciPy.
- [ ] Consider optional balanced-Schur/Hessenberg integration and safer
  componentwise scaling for very heterogeneous input magnitudes.

## 2026-11-07 — Upper Hessenberg reduction

- [x] Add independent real/complex two-sided Householder `hessenberg`,
  with optional explicit orthogonal/unitary similarity transformation.
- [x] Validate structure, reconstructability, trace and Frobenius norms,
  nonfinite/error paths, empty shapes and extreme scales against SciPy.
- [ ] Add optional Schur input balancing and efficient blocked Hessenberg
  reduction if future eigensolvers require them; the current unblocked
  implementation targets moderate orders.

## 2026-11-07 — Scalable dense matrix equations

- [x] Add Schur/Bartels–Stewart `solve_sylvester` alongside direct Kronecker elimination.
- [x] Add single-Schur continuous and discrete Lyapunov solvers with independent tests.
- [x] Preserve original interfaces; select `auto`, `direct`, or `schur` explicitly.
- [x] Support order-40 matrices without forming order-1600 Kronecker systems.
- [ ] Consider LAPACK-backed Schur at larger orders, robust scaling and
  iterative refinement for nearly resonant operators.

## 2026-11-07 — Orthogonal/unitary Procrustes alignment

- [x] Add SciPy-compatible `orthogonal_procrustes` for real, complex, and
  mixed-input matrices, using SciFort's SVD without new dependencies.
- [x] Support empty and rank-deficient matrices, and scaled covariance
  multiplication with overflow-aware restoration of the singular-value sum.
- [x] Add native Fortran tests, SciPy randomized comparisons, a runnable
  example, and numerical API/provenance documentation.
- [ ] Consider batch-leading dimensions and optional high-performance
  LAPACK-backed SVD for larger problems.

## 2026-11-07 — Dense subspace geometry

- [x] Add real/complex `orth`, `null_space`, and `subspace_angles` based on
  economy SVD, including mixed real/complex principal angles.
- [x] Complete the right-singular basis for null spaces of wide matrices,
  with pivoted, reorthogonalized canonical vectors.
- [x] Preserve small principal angles with projection-residual sines;
  validate against SciPy using projectors, angles, and residuals.
- [x] Add optional full-size SVD bases for rectangular real/complex matrices.
- [ ] Consider optional LAPACK-backed large-scale SVD pathways.

## Dense SVD and least-squares coverage (2026-11-07)

- [x] Real and complex economy `svd` using one-sided Jacobi sweeps.
- [x] Moore-Penrose `pinv` and SVD-based rank-deficient minimum-norm `lstsq`.
- [x] Optional `method` for QR-only, SVD-only and QR-then-SVD fallback.
- [x] Full-size U/Vh output, singular-value-only `svdvals`, and `diagsvd`.
- [x] Dedicated values-only real/complex Jacobi path and SVD-based `cond2`/`matrix_rank` diagnostics.
- [ ] Add rank-deficient reusable QR solves via complete orthogonal decomposition.
- [ ] Large-scale LAPACK-backed SVD option for efficiency and difficult
  clustered singular values; preserve an intrinsic-only fallback.

## 2026-11-07 — Dense Householder QR and least squares

- [x] Add real/complex, reduced/complete Householder QR with optional column pivoting.
- [x] Add compact factor/unpack, reusable QR solves, and full-rank `lstsq` for tall and wide matrices.
- [x] Test native failure paths, exact/nonexact least squares, and independent SciPy references.
- [ ] Add complete orthogonal decomposition or SVD for rank-deficient minimum-norm least squares.
- [ ] Consider a blocked Householder implementation and optional BLAS/LAPACK backend.

## 2026-11-07 — Dense Cholesky and positive-definite systems

- [x] Add real/complex upper/lower `cholesky`, `cho_factor`, `cho_solve`.
- [x] Support factor reuse with vectors and multiple right-hand sides.
- [x] Validate factors and solutions against SciPy in checked and `-O3` builds.
- [x] Add dense QR (Householder) and reusable QR/least-squares solves.

## 2026-11-07 — Dense triangular solves and reusable LU

- [x] Add real/complex upper/lower `solve_triangular`, unit diagonal, multiple
  right-hand sides, and normal/transpose/adjoint operations.
- [x] Add dense partial-pivoted `lu_factor` and reusable `lu_solve` with the
  same `N/T/C` operations and one-based pivot indices.
- [x] Validate against SciPy in checked and optimized GNU Fortran builds.
- [ ] Consider QR dense decomposition reuse, iterative refinement,
  and optional BLAS/LAPACK backends for larger dense problems.

## 2026-11-07 — Banded dense linear systems

- [x] Add real and complex `scipy.linalg.solve_banded`-style Gaussian elimination
  with local row pivoting and compact fill-band workspace.
- [x] Add upper/lower Hermitian `solveh_banded`, `cholesky_banded`, and
  `cho_solve_banded`, for vector and multiple right-hand sides.
- [x] Compare with SciPy using independent randomized tests in checked/O3
  builds, stress row pivoting, and verify order-50,000 tridiagonal behavior.
- [ ] Consider reusable general banded LU factors, iterative refinement,
  equilibrating badly scaled systems and batch-leading dimensions.

## 2026-11-07 — General dense analytic matrix functions

- [x] Add real/complex `funm` and `funm_real` using Schur--Parlett recurrence.
- [x] Support derivative callbacks for contiguous Jordan/repeated blocks,
  explicit failure for unresolved couplings, and independent SciPy tests.
- [x] Add global Schur reordering and blocked Sylvester/Parlett recurrence
  for noncontiguous repeated eigenvalue clusters, including semisimple cases
  without derivative callbacks.
- [ ] Consider more robust diagonal-block evaluation, LAPACK Schur backend,
  certified error estimates and 1-norm condition estimation.

## 2026-11-07 — Frobenius condition number of matrix exponentials

- [x] Add real and complex `expm_cond` with Lanczos/Fréchet estimation.
- [x] Special-case diagonal matrices with divided differences, normalize
  scalar shifts, and independently check extreme shifts and large orders.
- [x] Add native tests and randomized SciPy comparisons in both build modes.
- [x] Add general dense `funm` using the existing Schur/Parlett backend.
- [ ] Consider faster block Lanczos/1-norm condition estimators for large dense matrices.

## 2026-11-07 — Dense matrix trigonometric and hyperbolic functions

- [x] Add real and complex `sinm`, `cosm`, `tanm`, `sinhm`, `coshm`, `tanhm`.
- [x] Preserve tiny arguments with paired Taylor series and double-angle
  identities; report singular tangent/cosh poles explicitly.
- [x] Add native regression and independent SciPy comparisons in checked/-O3.
- [x] Add Fréchet-derived `expm_cond` estimate (see section above).
- [ ] Consider reordered block Schur/Parlett and optional high-performance
  LAPACK Schur backend for larger matrices.

## 2026-11-06 — Polar factorization and exponential sensitivity

- [x] Add real/complex rectangular `polar` with right/left PSD factors,
  rank-deficient handling, and independent SciPy comparisons.
- [x] Add real/complex `expm_frechet` with checked/optimized comparisons.
- [ ] Consider specialized Padé derivatives and Schur-based polar backends
  for large dense matrices (optional BLAS/LAPACK); retain lightweight core.

## 2026-11-06 — Dense matrix functions

- [x] Implement real and complex `sqrtm` and `signm` for small dense matrices.
- [x] Validate principal branch-cut handling and nonlinear matrix residuals.
- [ ] Add optional high-performance LAPACK Schur backend for larger matrices.
- [x] Implement small dense `logm` and `fractional_matrix_power` with complex branches.
- [ ] Upgrade `logm`/fractional powers to Schur-Pade for difficult matrices and
  cover non-integer powers of singular PSD matrices.

## 2026-11-06 — Dense matrix-equation solvers

- [x] Add real and complex Sylvester and continuous/discrete Lyapunov solvers,
  with checked and optimized SciPy reference comparisons.
- [ ] Replace the small-system Kronecker method with scalable Schur-based
  Bartels–Stewart routines where suitable LAPACK support is available.
- [x] Add real/complex continuous/discrete algebraic Riccati solvers for
  standard and generalized equations (matrix sign / doubling methods).
- [ ] Replace the Riccati methods with robust Schur/QZ or a selectable LAPACK
  backend for larger and poorly separated problems.
- [ ] Add optional multiple-matrix batch processing.

## 2026-11-06 — General complex matrix-function actions

- [x] Add complex `funm_multiply_krylov` for CSR/CSC and matrix-free operators.
- [x] Support complex time scaling, multiple RHS, Hermitian Lanczos and
  general Arnoldi, plus independent randomized SciPy comparisons.
- [ ] Consider adaptive Krylov dimension selection and rigorous a posteriori
  error estimators for analytic matrix functions.
- [ ] Add reusable projected matrix-function callbacks and more efficient
  block Arnoldi for multiple related right-hand sides.

## Complex sparse matrix exponential actions (2026-11-06)

Completed: complex CSR/CSC and callback-based `expm_multiply`, including block
right-hand sides and intervals. Complex general matrix-function Krylov actions have now been added.
Remaining opportunities: time-varying operators and adaptive norm
estimation for opaque callbacks. See [API](docs/sparse_complex_expm.md).

## 2026-11-06 — Nonsymmetric complex sparse eigenpairs

- [x] Add `eigs_complex` for complex CSR/CSC with six eigenvalue-selection modes.
- [x] Validate restarted Arnoldi true residuals, repeated eigenvalues,
  zero cases, and explicit exhaustion status against independent references.
- [x] Add complex sparse-LU shift-invert mode (`sigma`) for interior spectra.
- [x] Add matrix-free complex eigen-operators with optional inverse callbacks,
  and a matrix-matching reusable LU shift cache for CSR/CSC calls.
- [ ] Add generalized complex eigenproblems, user-supplied linear solves for
  arbitrary preconditioned operators, and block Arnoldi/Schur locking.
- [ ] Consider block restarted Arnoldi/Schur locking for repeated or clustered spectra.

## 2026-11-06 — Complex LSMR least-squares extension

- [x] Implement `lsmr_complex` for CSR and CSC complex binary64 rectangular matrices.
- [x] Add regularized, initial-guess, rank-deficient and over/underdetermined cases.
- [x] Validate against SciPy LSMR and dense NumPy SVD references in checked/-O3 builds.
- [x] Implement complex sparse `svds` (Golub--Kahan projected SVD, largest/smallest triplets).

## 2026-11-06 — Complex sparse least squares

- [x] Add complex binary64 `lsqr_complex` for rectangular CSR/CSC matrices,
  preserving conjugate-transpose arithmetic, regularization and initial guesses.
- [x] Add independent SciPy/NumPy randomized comparisons and native tests.
- [x] Implement complex `lsmr` with external BSD-3-Clause attribution.
- [x] Add complex `svds`; remaining: restarted and shift-invert methods for challenging partial/SM searches.

## 2026-11-05 -- Complex sparse CSR/CSC and Krylov methods

- [x] Complex one-based CSR/CSC matrices, COO coalescing, dense and complex
  N-D COO interoperability, transpose and conjugate-transpose, matvec/matmat.
- [x] Complex CG, BiCGSTAB and restarted GMRES with optional Jacobi
  preconditioning, CSR/CSC overloads and true-residual reporting.
- [x] Extend to complex ILU(0) with adjoint solves and BiCG/QMR, including
  optional ILU preconditioning of BiCGSTAB/GMRES.
- [x] Add baseline complex sparse direct LU (`splu_complex`, `spsolve_complex`) with row pivoting and factor reuse.
- [ ] Extend to optimized complex sparse direct ordering/pivoting and larger-spectrum eigenproblems.
- [ ] Consider complex sparse eigenproblems, sparse matrix-function actions,
  and nonbinary64 precision where supported by the broader library.

## 2026-11-05 -- Complex N-D COO storage and tensor operations

Implemented a separate `coo_nd_complex_array` storing `complex(dp)` values,
with construction, canonicalization, reshape, conjugation, broadcasting,
reductions, products, tensor contractions and Hermitian inner products.
Next candidates: complex sparse direct solvers and additional eigensolvers;
indexed complex COO assignment; additional `einsum`-style contractions.

## 2026-11-05: COO sparse storage and extraction

- [x] Add one-based, real binary64, two-dimensional coordinate sparse storage.
- [x] Preserve unsorted duplicates and explicit zeros during construction;
  provide canonical CSR/CSC conversions, transpose, matvec, and data queries.
- [x] Add SciPy-style `coo_sum_duplicates`, `coo_eliminate_zeros`, and
  `sparse_find` (COO/CSR/CSC) with a separate SciPy comparison harness.
- [x] Add DOK and LIL sparse formats with editable single coordinates,
  structure-preserving CSR/CSC/COO conversion and independent SciPy tests.
- [x] Add a separate arbitrary-rank COO sparse-array type with duplicate and
  explicit-zero preservation, reshape, transpose, dense conversion, and lossless
  interoperability with the existing 2D COO type; verify against SciPy 1.17.
- [x] Add N-dimensional axis reductions (`sum_all`, `sum_axis`, `sum_axes`),
  broadcasted sparse COO arithmetic, and contracted-axis `tensordot` with
  sorted-key matching and independent SciPy/NumPy randomized tests.
- [x] Add N-D sparse integer/range/index-array slicing with Cartesian multi-axis
  selection and paired-coordinate gathering, validated against SciPy/NumPy.
- [x] Add N-D full/axis boolean masks and transactional coordinate scatter
  (replacement and accumulation), with exact duplicate/zero semantics.
- [x] Add scalar sparse predicates, COO-valued Boolean masks, coordinate-set
  filtering and inverted filters without full-shape allocation. Explicitly
  reject comparisons that are true for implicit zeros.
- [x] Add complement-aware sparse Boolean mask algebra (AND/OR/XOR/NOT),
  sparse scalar predicates true for implicit zeros, and sparse COO filtering.
- [x] Add right-aligned Boolean mask broadcasting, broadcast-aware logical
  operations, axis-wise `ANY`/`ALL`, and exact sparse 64-bit axis counts with
  implicit-true backgrounds and no full-volume allocation.
- [x] Implement complex-valued COO arrays and complex tensor contractions.
- [ ] Consider full NumPy advanced-index broadcasting and generalized
  Einstein summation.

## 2026-11-05: General matrix-function Krylov actions

- [x] Add `funm_multiply_krylov` with a dense real matrix-function callback,
  CSR and matrix-free operator overloads and separate Arnoldi/Lanczos paths.
- [x] Add adaptive projected Hessenberg storage, restart controls, explicit
  diagnostics, breakdown handling, and 600 randomized dense SciPy
  comparisons per build. All 32 native sparse tests passed in both builds.
- [ ] Consider complex-valued operators/functions and specialized numerical
  kernels for common `f` (`exp`, `phi`, `sqrt`, resolvents), where dedicated
  error estimators may outperform generic update-based stopping.
- [ ] Consider symmetric-matrix verification, adaptive Krylov dimensions,
  and convergence bounds for particular matrix-function families.

## 2026-11-05: Grid Laplacians and sparse structure (SciPy 1.17)

- [x] Add the matrix-free real N-dimensional `LaplacianNd` operator on
  arbitrary positive grid shapes (C-order flattened vectors), with
  Dirichlet, Neumann, and periodic boundary conditions.
- [x] Add matrix/multiple-RHS application, canonical CSR construction,
  and exact analytic eigenvalues/eigenvectors with `m` largest selection.
- [x] Add CSR/CSC structural `spbandwidth` and `is_sptriangular`, counting
  stored zeros without converting the sparse format.
- [x] Compare with SciPy 1.17 in checked and optimized GNU Fortran builds,
  including non-square grids, duplicated periodic neighbors, singleton axes,
  and rectangular sparse inputs with explicitly stored zeros.
  All 31 native sparse regression tests passed in each mode.
- [ ] Consider general grid spacing, arbitrary stencil order, and
  operator-based matrix-free integration with other SciFort solvers.

## 2026-11-05: Augmented and recycled sparse Krylov solvers

- [x] Add CSR LGMRES with reusable full-cycle correction vectors and
  configurable Krylov/outer dimensions.
- [x] Add GCROT(m,k) with orthonormal recycled images, oldest-first pruning,
  right preconditioning, and projection-based warm starts across changed A.
- [x] Validate both with native edge cases, indefinite systems, and 360
  independent SciPy comparisons per compiler mode.
- [ ] Consider complex and general `LinearOperator` variants, the GCROT
  smallest-singular-value truncation strategy, and callback hooks.

## 2026-11-05: More scipy.sparse algorithms

- [x] Implement CSR `kronsum(A,B)` by sorted sparse row-merge, retaining
  explicit zeros and rejecting nonsquare operands or unrepresentable output.
- [x] Add real sparse BiCG with transposed matrix-vector products and Jacobi.
- [x] Add real sparse CGS with Jacobi and reusable ILU(0) preconditioning.
- [x] Compare both solvers and Kronecker sums to independent SciPy 1.17
  references in bounds-checked and optimized GNU Fortran builds.
- [ ] Add BiCG general transpose preconditioners and matrix-free operators.
- [x] Add real CSR QMR (Jacobi) and TFQMR (Jacobi/ILU(0)) with true residual checks.
- [ ] Add complex and matrix-free sparse QMR/TFQMR, and other missing scipy.sparse.linalg algorithms.
- [ ] Consider algorithms using integer-valued sparse matrices and complex
  CSR/CSC storage after the numerical core's type strategy is established.

## 2026-11-03: scipy.cluster.vq core

- [x] Add `whiten` (vector and matrix), `vq` with one-based centroid indices,
  `kmeans2` (fixed Lloyd steps), and `kmeans` (mean-distance convergence).
- [x] Add explicit seeded initialization for sampling observations and
  K-means++; test provided-centroid behavior against SciPy.
- [ ] Add SciPy's covariance-based Gaussian `kmeans2(minit='random')` seeding.
- [x] Add one-dimensional overloads for K-means and K-means2.
- [ ] Consider optimized batched distance updates, and a documented optional
  Hartigan-Wong algorithm with verified
  permissive-source attribution if it is adopted.

## 2026-11-03: scipy.cluster.hierarchy foundation

- [x] Add seven independent Lance-Williams agglomerative linkage methods
  with observation/condensed-distance input and one-based node indexing.
- [x] Add `fcluster_distance`, `fcluster_maxclust`, `cophenetic_distances`,
  `leaves_list`, native regression tests, and independent SciPy comparisons.
- [x] Add `inconsistent`, `maxdists`, `maxrstat`, `maxinconsts`,
  `fcluster_inconsistent`, monotone flat-cluster criteria, and `cut_tree`
  with exact cluster-count and height selectors.
- [x] Validate new hierarchy diagnostics with randomized SciPy 1.17
  comparisons and native checked/optimized regression tests.
- [x] Add `optimal_leaf_ordering`, `to_tree`, `leaders` and `correspond`.
- [x] Add unified `fcluster`, observation-to-label `fclusterdata`, and
  scale-normalized `cophenetic_correlation`, with independent SciPy tests.
- [ ] Add faster O(n^2) single-linkage and nearest-neighbor-chain algorithms,
  and remaining tree validation/conversion helpers.

## 2026-11-02: scipy.spatial.transform rotation foundations

- [x] Implement finite SO(3) `Rotation.create_group` symmetries (`Cn`, `Dn`,
  `T`, `O`, `I`) and scalar/stack `Rotation.reduce` with one-based witness
  indices; validate closure and randomized SciPy symmetry reduction.
- [x] Add SE(3) `RigidTransform` scalar and 1D stack APIs: 4x4 matrices,
  translation/rotation composition, inversion, point application, exp/log
  coordinates, dual quaternions, screw powers, weighted means and broadcasts.
- [x] Add rank-two `RigidTransform` grids with broadcasting, constructors,
  matrix/twist/dual-quaternion conversions, axis-wise weighted means and
  independent SciPy 1.17 tests.
- [x] Add scalar, stack, and rank-two-grid Davenport-angle conversions;
  include generalized axes, gimbal-lock handling, and independent SciPy checks.
- [x] Add runtime-rank SO(3)/SE(3) transform batches with NumPy-like
  right-aligned broadcasting, component-last constructors, and shape metadata.
- [ ] Add higher-rank slicing/reductions, SVD projection of imperfect matrices,
  and symmetry reduction across arbitrary-rank batches.

- [x] Add native `rotation_type` with quaternion, matrix and rotation-vector
  construction and conversions, composition, inverse and 3D vector application.
- [x] Implement shortest-path two-keyframe and piecewise nonuniform Slerp.
- [x] Validate under native checked/optimized builds and independent SciPy 1.17
  randomized comparisons, including near-pi and near-identity orientations.
- [x] Add 24 three-axis Euler-angle conventions, 1/2-axis constructors,
  gimbal-lock reporting, and modified Rodrigues parameter conversions.
- [x] Validate Euler/MRP conversions for all conventions, singularities,
  non-principal MRPs, and extreme scales against independent SciPy 1.17 results.
- [x] Add weighted chordal quaternion rotation means, including nonnegative
  weights, antipodal quaternion invariance, and SciPy numerical comparisons.
- [x] Add `Rotation.align_vectors` with finite weighted Wahba fits,
  an infinite-primary-weight constraint, root residuals, and sensitivity.
- [x] Add SO(3) cubic `RotationSpline` with continuous angular rate and
  acceleration; use a block-tridiagonal native Fortran solver and SciPy
  reference comparisons.
- [x] Reject catastrophic spline coefficient divergence with a distinct
  `spatial_rotation_no_convergence` status and a near-antipodal regression.
- [x] Add `rotation_stack_type` with batched vector application, representation
  conversion, composition/inversion, powers, slicing, concatenation, and
  weighted means, including size-one broadcasting and SciPy reference tests.
- [x] Add Haar-uniform `Rotation.random` sampling with explicit reproducible
  caller-owned integer PRNG streams and native/SciPy statistical tests.
- [x] Implement a two-batch-axis `rotation_grid_type` with SciPy-compatible
  component-last representations, broadcasting, and weighted means by axis.
- [x] Generalize rotation batch composition and application to arbitrary
  NumPy-compatible shapes via `scifort_spatial_transform_nd`.
- [ ] Add richer rank-N slicing and reductions and more robust
  ill-conditioned attitude spline solvers.

## 2026-11-02: 2D/3D halfspace intersections

- [x] Add a native `scipy.spatial.HalfspaceIntersection` subset for bounded
  convex polygons and polyhedra with a supplied strictly feasible point.
- [x] Handle duplicate/redundant inequalities, dual geometry, 3D coplanar
  polar facets, and output unique primal vertices with one-based dual indices.
- [x] Test independent SciPy references, checked/optimized native Fortran,
  unbounded and invalid geometries, translated and anisotropic polytopes.
- [ ] Consider scalable halfspace intersection, a Chebyshev-center solver,
  adaptive precision predicates, and dimensions greater than three.

## 2026-11-02: 2D and 3D convex hulls

- [x] Implement native `scipy.spatial.ConvexHull` geometry for 2D/3D
  (one-based vertices/simplices, adjacency, outward equations).
- [x] Compute polygon perimeter and area and polyhedron surface area and
  volume, handling interior and duplicate points.
- [x] Verify using randomized SciPy oracle comparisons and native tests.
- [x] Fix the reproducible 26-point near-coplanar polar hull failure:
  validate topology/supporting planes against all input sites; retry with
  deterministic point orders; use a supporting-plane fallback if necessary.
  Add native and SciPy regressions covering permutations and cube-face clouds.
- [ ] Add robust adaptive predicates and scalable hull construction for
  large or near-degenerate sets; consider n-dimensional hulls separately.

## 2026-11-02: Spherical Voronoi construction

- [x] Add native `scipy.spatial.SphericalVoronoi`-style circle/sphere diagrams
  in 2D and 3D with region vertices and one-based CSR region indexing.
- [x] Add incremental 3D convex hull, topological cyclic region ordering,
  3D spherical polygon areas and SciPy-style 2D arc lengths.
- [x] Validate against SciPy with two randomized seeds and native edge cases
  in bounds-checked and optimized GNU Fortran builds.
- [ ] Consider robust geometric predicates for near-degenerate/coplanar facets,
  higher-dimensional spherical vertices, and scalable hull construction.

## 2026-11-02: Procrustes analysis

- [x] Implement `scipy.spatial.procrustes` for centered, Frobenius-normalized
  point-cloud alignment with rotation, reflection, and optimal scaling.
- [x] Implement the related real `scipy.linalg.orthogonal_procrustes`
  without translation or rescaling, using a standalone Jacobi polar factor.
- [x] Cover rank-deficient configurations, duplicates, nonunique rotations,
  invalid inputs, and scale/translation extremes in native regression tests.
- [x] Add independent randomized SciPy comparisons, complete spatial API
  docs, and an MIT provenance record.
- [x] Add `SphericalVoronoi` in 2D/3D without further Delaunay optimization.

## 2026-11-02: Hausdorff and spherical geometry

- [x] Implement `directed_hausdorff` with exact finite max-min Euclidean
  distances, deterministic one-based witness indices and early abandonment.
- [x] Add `hausdorff_distance` (SciFort extension) for symmetric distances
  and witness direction reporting.
- [x] Add arbitrary-dimensional `geometric_slerp` with scalar/vector path
  parameters, extrapolation, and explicit ambiguous-antipode handling.
- [x] Add native regressions and two independent SciPy comparison seeds for
  both checked and optimized builds.
- [x] Add `scipy.spatial.procrustes` and `SphericalVoronoi` geometry algorithms.
- [ ] Consider further new `scipy.spatial` algorithms.

## 2026-11-02: Weighted Boolean dissimilarities

- [x] Implement `dice`, `rogerstanimoto`, `russellrao`, `sokalsneath`, and
  `yule` in the established `point_distance`, `pdist`, and `cdist` API.
- [x] Support per-feature weights without introducing dependencies, and
  normalize the weights to avoid overflow of finite-weight sums.
- [x] Match SciPy's Boolean interpretation and bulk zero-denominator
  conventions; document intentional differences for nonbinary numeric inputs.
- [x] Add independent native regression tests and randomized SciPy `pdist`
  and `cdist` comparisons in both GNU Fortran build modes.
- [ ] Add output buffers, reuse of covariance factorizations, and other
  missing `scipy.spatial` algorithms.

## 2026-11-02: Covariance and Jensen--Shannon distance metrics

- [x] Implement `seuclidean`, `mahalanobis`, and `jensenshannon` in the
  existing `point_distance`, `pdist`, and `cdist` Fortran APIs.
- [x] Accept user-specified `variance`/`inv_cov`, estimate default sample
  statistics from `pdist` or combined `cdist` observations, and reject
  insufficient/rank-deficient covariance matrices.
- [x] Add optional Jensen--Shannon `base` and stable probability
  normalization; preserve old positional argument compatibility.
- [x] Test all three algorithms against SciPy in checked and optimized GNU
  Fortran; add scalar/degenerate-input regression tests.
- [x] Add five weighted Boolean distance metrics (2026-11-02).
- [ ] Add output buffers and faster batched Mahalanobis factor application
  for high-dimensional point clouds.

## 2026-11-01: Weighted and multi-radius KD-tree counting

- [x] Overload `count_neighbors` with a vector-radius `integer(int64)` result;
  preserve original scalar `integer(int64)` calls without source changes.
- [x] Add `weighted_count_neighbors` for scalar/array radii and optional
  source/target point weights; support signed weights and unit defaults.
- [x] Add cumulative or noncumulative radius bins, repeated radii, internally
  sorted arbitrary-order cumulative radii, and shared traversal with
  subtree-sum and interval-update pruning.
- [x] Test signed/one-sided weights, periodic axes, Minkowski norms and
  edge cases against SciPy in checked and optimized builds.
- [ ] Add approximate `eps` searches, and further metric families in
  `scifort_spatial_distance`.

## 2026-11-01: Periodic and Minkowski KD-tree queries

- [x] Add scalar and per-axis periodic `boxsize` constructors with periodic
  site validation, wrapped query points, and mixed periodic/nonperiodic axes.
- [x] Support exact Minkowski `p>=1` and `p=+infinity` in nearest,
  query-ball, pair enumeration, cross-tree, count, and sparse-matrix queries.
- [x] Generalize point/node and dual-tree bounding-box pruning to periodic
  p-norms, with scaled distance calculation to resist underflow/overflow.
- [x] Add Fortran native regression and SciPy comparisons across two seeds
  and both checked/optimized GNU Fortran builds; maintain existing APIs.
- [x] Support multiple radius bins and weighted pair counting.
- [ ] Support approximate `eps` searches.

## 2026-11-01: KD-tree pair and cross-tree queries

- [x] Extend the existing independent balanced KD-tree with exact
  `query_pairs`, `query_ball_tree`, `count_neighbors`, and
  `sparse_distance_matrix`, exporting the same tree via `scifort_spatial`.
- [x] Return deterministic one-based pair arrays and CSR/COO triplets,
  including duplicate points, zero-radius matches, and explicit zeros.
- [x] Add dual-tree bounding-box pruning and fully included subtree counting
  with `integer(int64)` counts. Use a scaled norm for extreme magnitudes.
- [x] Add deterministic native and independent SciPy 1.17 validation in
  checked and optimized GNU Fortran configurations.
- [x] Add periodic coordinate boxes and non-Euclidean Minkowski pair searches
  without changing Euclidean nearest-query defaults.
- [x] Add separate weighted counting and radius histograms.
- [ ] Add approximate `eps` searches.

## 2026-11-01: `scipy.spatial.distance` introduction

- [x] Add `scifort_spatial` with `point_distance`, `pdist`, `cdist`, and
  `squareform` (condensed/square) for real binary64 data and 11 named metrics.
- [x] Support weights, arbitrary positive Minkowski exponents and SciPy's
  pair-index convention; validate with 495 randomized cases per seed in
  both checked and optimized GNU Fortran configurations.
- [x] Add standardized Euclidean, Mahalanobis, and Jensen--Shannon distances.
- [ ] Add more Boolean distances, reusable covariance factors, output buffers,
  and callback metrics; investigate stable extreme-value summation.

## 2026-11-01: `make_lsq_spline` complete-knot weighted fitting

- [x] Add real scalar and multichannel `scipy.interpolate.make_lsq_spline`
  with user-specified **complete** knot vectors, degrees 0–15, repeated
  interior knots, nondecreasing observation positions, and zero weights.
- [x] Reuse SciFort's compact banded Givens QR kernel with shared factorization
  across output channels, preserving legacy `LSQUnivariateSpline` restrictions.
  Reject unsupported `method='norm-eq'` instead of silently changing methods.
- [x] Add independent SciPy 1.17 reference values, 480 randomized SciPy
  coefficient comparisons across two seeds, checked/optimized native tests,
  and singular-/invalid-input regressions. No new Fortran dependencies.
- [ ] Extend `make_lsq_spline` to complex data and arbitrary ordinate axes;
  assess additional numerical scaling for extremely ill-conditioned knots.

## 2026-11-01: Periodic `make_interp_spline` interpolation

- [x] Construct periodic B-spline coefficients for degrees 2–15 using cyclic
  knots and endpoint derivative matching, supporting real scalar and
  multichannel ordinates with a shared collocation factorization.
- [x] Preserve SciPy's special degree-0/1 behavior, validate periodic
  endpoint data, reject explicit knots, and expose periodic wrapping,
  derivatives, integration, and shape-preserving knot insertion.
- [x] Add two independent SciPy comparison seeds (150 + 120 cases per
  build), and all 30 checked/optimized native interpolation tests.
- [ ] Optimize large cyclic interpolation systems with a banded solver;
  positive-s smoothing on rectangular grids still needs the optional
  FITPACK backend and separate genuine-upstream validation.

## 2026-11-01: Bivariate FITPACK surface calculus

- [x] Add exact non-refitting `BivariateSpline.partial_derivative`-style TCK
  transformations for both axes, with independent x/y degrees and statuses.
- [x] Add tensor-product Boehm knot insertion in either axis, accepting TCKs
  from rectangular interpolation and the optional FITPACK backend.
- [x] Validate coefficients, knot vectors, mixed derivatives and integrals
  against SciPy/BSpline references in checked and optimized native builds.
- [ ] Add optional FITPACK rectangular-grid smoothing (`RectBivariateSpline`
  with positive `s`) and explore efficient bivariate multi-channel support.

## 2026-10-31: PPoly / BPoly / Hermite interpolation

- [x] Add real one-dimensional `PPoly` and `BPoly` with multichannel output,
  Horner and de Casteljau evaluation, analytic calculus, integration, periodic
  extension, and reciprocal conversion between the two polynomial bases.
- [x] Implement `BPoly.from_derivatives` with nonuniform endpoint derivative
  counts and exact degree elevation; bridge SciFort's cubic objects to `PPoly`.
- [x] Extend to complex coefficients in PPoly/BPoly.
- [ ] Extend to arbitrary-rank ordinate arrays.
- [x] Add `PPoly.from_spline`, real polynomial `%roots()` / `%solve()` with
  tangency and discontinuity detection, and complex PPoly/BPoly calculus.
- [x] Add `.extend()` and `BPoly.from_derivatives(orders=...)` handling.
- [x] Implement optional real-polynomial root callback and an opt-in
  version-pinned `polyroots-fortran` adapter (ABI tested with a mock).
- [ ] Build, benchmark, and validate the adapter with the **real upstream**
  `polyroots-fortran` package, and assess its system BLAS/LAPACK requirements.
- [x] Add global real/complex barycentric interpolation with multichannel
  data, derivatives and incremental node insertion.
- [x] Implement `KroghInterpolator` and Floater--Hormann rational interpolation.
- [ ] Add multiaxis ordinate support for all 1D interpolators.

## 2026-10-31: Tensor-product regular-grid B-spline interpolation

- [x] Add `RegularGridInterpolator`-style `slinear`, `cubic` and `quintic`
  methods, complex/trailing-channel support, 1–6D mixed ascending/descending
  coordinates, per-axis not-a-knot collocation, and extrapolation.
- [x] Support mixed partial derivatives via optional `nu` with native 6D and
  independent randomized SciPy tests across two seeds and two compiler builds.
- [ ] Replace dense per-axis LU by a banded collocation solver for large grids,
  and benchmark against an optional pinned `bspline-fortran` backend.
- [x] Implement 1–6D `pchip` interpolation.
- [ ] Optimize sparse support contractions for many query coordinates.

## 2026-10-31: General nonuniform B-spline interpolation

- [x] Expose native `BSpline(t,c,k)` / `make_interp_spline` for real scalar and
  multichannel arrays with degrees 0–15, analytic derivatives and integrals,
  arbitrary repeated knots, knot insertion and a dense design matrix.
- [x] Implement odd/even degree not-a-knot placement and natural/clamped cubic
  boundary constraints. Validate checked/optimized builds against SciPy across
  two random seeds, retaining all previous interpolation regressions.
- [ ] Benchmark a pinned optional FPM `bspline-fortran` backend against the
  independently implemented 2D–6D tensor-product spline support (retain
  upstream notices, verify the Windows build).
- [ ] Add complex coefficients, arbitrary ordinate axes and user-provided
  derivative boundary constraints. Periodic interpolation is now available.
- [ ] Replace the moderate-sized dense collocation solver with an efficient
  pivoted banded method and add sparse B-spline design matrices.

## 2026-10-31: First `scipy.interpolate` piecewise-cubic algorithms

- [x] Implement PCHIP monotone cubic Hermite knot derivatives, Akima and
  modified Akima knot slopes, and user-specified `CubicHermiteSpline` slopes.
- [x] Add native reusable polynomial evaluator with derivatives, extrapolation
  policy, exact piecewise integration and independent SciPy regression tests.
- [ ] Extend real 1-D interpolators to multidimensional ordinate arrays and
  arbitrary interpolation axes, retaining SciPy's coefficient conventions.
- [x] Add `CubicSpline` (`not-a-knot`, natural, clamped and periodic boundary
  conditions) and corresponding native reference and regression coverage.
- [x] Add `PPoly` root solving and periodic extrapolation, and consider inversion of
  strictly monotone PCHIP interpolants.

## 2026-10-31: Fixed-weight multitaper VAR(p) population targets

- [x] Extend exact finite-window VAR(p) point and analytic band spectra to arbitrary deterministic taper banks, including unequal fixed weights and per-taper energy normalization.
- [x] Add conditional two-group canonical population targets for fixed-weight multitaper spectral matrices and independent SciPy Lyapunov/direct-pair/Cholesky-SVD validation.
- [x] Validate checked/optimized native tests, invalid-input behavior and the full public signal dependency closure without adding runtime dependencies.
- [x] Add a native exact-band multichannel fixed-taper spectrogram and exploratory DPSS VAR(p) moving-block bootstrap coverage grid varying taper count, time-bandwidth product, weighting and block length. Larger replication studies remain necessary.
- [ ] Consider a shared-covariance, pooled-taper-autocorrelation kernel to avoid repeated Lyapunov solves for larger banks; maintain exact per-taper normalization.
- [ ] Do not conflate the deterministic weighted-target formula with the nonlinear expectation under adaptive, data-dependent taper weights.

## 2026-10-30: Multimode conditional VAR(p) bootstrap coverage

- [x] Compare five bootstrap interval methods for conditional two-mode canonical coherence in four stable VAR(2)/VAR(3) models, using stationary initialization and exact Hann-tapered window/band spectra.
- [x] Add independent Gaussian-quadrature and SVD reference checks, 64 independent datasets per model, bias/RMSE reporting, Wilson Monte Carlo uncertainty, and a reproducible 180-row results CSV.
- [ ] Increase replication counts and tune bootstrap block length, window duration, and sample size to assess calibration sensitivity; develop subspace-level inference for nearly tied canonical modes.

## 2026-10-30: Exact VAR(p) population spectra and conditional canonical bands

- [x] Extend finite-taper population spectral targets from real VAR(1) to stable real VAR(p), using a companion-state discrete Lyapunov solve with semidefinite state innovations, a positive-definite observed innovation covariance, and explicit lagged covariance propagation.
- [x] Add **analytic** integration of finite-window cross spectra over multiple possibly overlapping frequency bands; preserve complex phase, frequency periodicity, exact band additivity, energy normalization, and optional band-mean spectral density.
- [x] Expose exact-point and exact-band population conditional canonical-coherence targets through the public `scifort_signal` API, reusing the original regularization and canonical solver, with arbitrary channel groups and optional conditioning.
- [x] Validate against independent SciPy/NumPy companion Lyapunov, direct sample-pair quadrature and conditional Cholesky/SVD references; simulate VAR(2) window ensembles, retain prior public regressions, and document mathematical/statistical limits.
- [x] Use these targets in a conditional VAR(2)/VAR(3) two-mode coverage study with exact finite-window Hann targets; add explicit regularization and DPSS taper alternatives in future work.
- [ ] Develop statistically justified inference for near-degenerate canonical subspaces.

## 2026-10-30: Exact tapered VAR(1) population target and inference diagnostics

- [x] Provide an analytic finite-taper VAR(1) expected spectral matrix using stationary covariance and taper autocorrelations, with phase-aware real/complex CSD output, input validation, SciPy numerical references, and Monte Carlo validation.
- [x] Run a reproducible four-scenario study against `C(E[tapered spectrum])` rather than `C(S(0))`, comparing five interval methods and identifying invalid intervals. Preserve the old exploratory zero-frequency studies separately.
- [ ] Increase Monte Carlo replication counts and vary taper lengths, block sizes, VAR persistence, adaptive weighting and regularization; quantify point-estimator bias relative to the proper functional.
- [ ] Study multivariate group conditioning, multi-mode canonical coherence and bootstrap confidence sets for near-degenerate subspaces.

## 2026-10-30: Canonical bootstrap interval construction

- [x] Add opt-in basic/reverse and bias-corrected (BC, not BCa) percentile
  intervals to existing canonical moving-block bootstrap APIs, preserving all
  default percentile results, draws and public-ABI signatures.
- [x] Handle numerical ties in BC ranks, bound basic intervals to the valid
  coherence parameter space, and report unidentified bias adjustments.
- [x] Compare all three interval types against independent SciPy normal
  quantiles and NumPy canonical refits and run a small reproducible VAR(1)
  coverage diagnostic. Document undercoverage and scenario sensitivity.
- [x] Add opt-in studentized/bootstrap-t intervals with separate block-jackknife
  variances computed for the original and each moving-block resample, including
  weighted pseudovalues for non-divisible block lengths, validity checks,
  independent NumPy reference, a diagnostic VAR(1) coverage comparison and
  complete public API tests.
- [x] Implement and independently validate block-aware BCa acceleration,
  including original-data delete-group influences, adjusted tail quantiles,
  invalid-interval handling and equal/unequal block sizes.
- [ ] Improve bootstrap-t calibration near canonical-mode boundaries and
  repeated-root regions; compare empirical coverage with BCa on the proper
  smoothed population estimand.
- [ ] Expand block/bootstrap interval coverage experiments to more VAR
  parameters, independent random seeds, and a carefully defined spectral
  estimand, without claiming nominal coverage from small trials.

## 2026-10-30: Canonical moving-bootstrap block-length pilot

- [x] Implement exploratory autocovariance-based automatic block-length selection from full complex Hermitian window spectral matrices; preserve explicit block lengths and deterministic seeds.
- [x] Add standalone pilot diagnostics and automatic real/complex static/band multitaper bootstrap support, with independent NumPy checks and a reproducible bivariate AR(1) coverage example.
- [x] Add opt-in statistic-specific delete-one-window influence pilot, validate against independent NumPy canonical fits, and compare to full-matrix selection in a small three-channel VAR(1) coverage study.
- [ ] Expand VAR/AR processes, parameter grids, block lengths, bootstrap interval types, and Monte Carlo coverage studies; compare bias-corrected inference and dependence diagnostics.
- [ ] Consider studentized or bias-aware intervals, with calibration diagnostics for boundary and nearly tied canonical modes.

## 2026-10-30: Moving-block canonical spectral bootstrap

- [x] Add ordinary/circular moving-block bootstrap on complete nonoverlapping
  spectral windows, including local reproducible RNG, non-divisible blocks,
  direct real/complex adaptive multitaper and band-integrated estimation.
- [x] Provide standard errors, percentile bounds, bias estimates, aggregate
  group coherence and explicit validity masks for non-isolated modes.
- [x] Test checked/optimized Fortran builds, 56 randomized NumPy/SciPy
  reference cases and the full 100-module umbrella API closure.
- [ ] Study block-length selection and finite-sample empirical coverage under
  representative autoregressive and multichannel time-series processes.
- [ ] Develop uncertainty summaries for canonical subspaces and simultaneous
  frequency-band inference without claiming generic calibration.

## 2026-10-30: Contiguous-block canonical coherence resampling

- [x] Add optional equal-length contiguous block deletion to the existing
  canonical window jackknife, plus descriptive block-specific public APIs.
  Preserve the delete-one-window default and all old tests.
- [x] Validate metadata, all refit formulas, real/complex adaptive multitaper
  signals, integrated spectral bands, input rejection and edge cases against
  independent native and NumPy/SciPy references.
- [x] Add overlapping moving-block bootstrap with reproducible local seeds,
  ordinary/circular sampling, pointwise percentile bounds and integrated bands.
- [ ] Investigate unequal-length deletion blocks, automatic block-length
  sensitivity analysis and empirical coverage studies.
- [ ] Add uncertainty estimation for canonical coefficient **subspaces** that
  remains meaningful when canonical strengths are nearly degenerate.

## 2026-10-30: Window-jackknife canonical coherence uncertainty

- [x] Implement pointwise delete-one-**window** jackknife SEs and bias for
  conditional canonical squared coherence strengths and aggregate group
  coherence. Require at least four nonoverlapping, monotonically labeled
  windows; do not treat neighboring Fourier bins as independent.
- [x] Add approximate Fisher-transformed normal 95% intervals and validity
  flags rejecting zero, perfect, nearly tied or singular modes. Support
  post-window band integration, fixed ridge/shrinkage and real/complex
  equal/eigen/adaptive multitaper signal entry points.
- [x] Pass checked and optimized native tests, 100 independent NumPy/SciPy
  references (357 modes, max difference 5.88e-15), and the complete
  99-module public API closure with 14 public tests.
- [ ] Calibrate coverage under temporal dependence, consider longer block
  resampling and uncertainty for canonical subspace orientations; avoid
  interpreting pointwise jackknife limits as significance tests.

## 2026-10-30: Covariance-metric canonical subspace tracking

- [x] Expose conditioned, regularized A/B autospectral blocks via opt-in
  `return_covariances` on canonical estimators and a corresponding opt-in
  `covariance_metric` on the subspace tracker. Compare adjacent modes using
  a shared midpoint Hermitian metric. Preserve unit/coordinate invariance
  for fixed group transforms, validate invalid-metric and reset behavior,
  and cross-check complex Procrustes tracking against NumPy/SciPy.
- [ ] Investigate smoothing and inferential uncertainty for mode strengths and
  subspace directions under dependent frequency bins and adaptive weighting.
  Pointwise nonoverlapping-window jackknife is a first exploratory step,
  not a calibrated solution (see section above).

## 2026-10-29: Conditional coherence between channel groups

- [x] Add determinant-based group coherence for two disjoint channel groups,
  optionally conditioning on a third group, using Hermitian Cholesky
  log-determinants and stable small-coherence evaluation.
- [x] Support real/complex static and time-frequency DPSS spectral matrices,
  existing diagonal and ridge shrinkage, LOO adaptive coefficients and
  exact post-integration frequency-band group coherence.
- [x] Verify partial-coherence and multiple-coherence reduction identities,
  196 randomized NumPy canonical-SVD scenarios (1,418 cells) in checked and
  optimized builds, and the complete public `scifort_signal` dependency closure.
- [x] Add individually sorted squared canonical coherence modes between
  disjoint channel groups, optionally conditioned on a third, for static,
  time-frequency and integrated-band spectra. Cross-check 240 randomized
  NumPy SVD cases and native/public integration tests.
- [x] Add opt-in conditional canonical coefficient vectors for both channel
  groups, with unit conditional spectral variances and positive paired cross
  spectra. Report heuristic flags for non-identifiable zero and tied modes.
  Verify checked/optimized builds and 150 independent NumPy SVD scenarios
  across 798 cells per build, with complete public-module regression tests.
- [x] Track isolated canonical directions across ascending neighboring
  frequencies with globally optimal assignment and paired phase alignment.
  Provide separate continuity-sorted results, original rank mappings,
  match scores and safe resets at FFT wraps, invalid bins and ambiguous modes;
  validate checked/optimized builds against SciPy Hungarian assignments.
- [x] Extend isolated-direction continuity tracking to near-degenerate
  canonical **subspaces**: group almost equal roots, compare principal angles,
  match blocks globally, and align paired coefficient bases with a complex
  unitary Procrustes rotation. Preserve true canonical eigenvalue ordering;
  treat rotated columns as subspace bases, not individual eigenvectors.
  Validate split/merge interruption, 100 seeded NumPy/SciPy reference cases
  and the complete public signal dependency closure.
- [ ] Investigate calibrated statistical uncertainty for mode strengths and
  eigenvectors under dependent frequency bins and adaptive taper weighting;
  numerical direction flags are not inferential confidence estimates.

## 2026-10-29: Band-integrated multivariate spectral coherence

- [x] Integrate complex Hermitian CSDs by piecewise-linear trapezoidal
  quadrature on nonuniform and FFT-wrapped frequency grids, including
  exact off-grid endpoints, overlapping bands and time-frequency cubes.
- [x] Add band-level partial and multiple coherence, conditional power and
  conditional phase, with real/complex direct DPSS interfaces and explicit
  post-integration ridge/shrinkage. Do not transfer per-bin CV shrinkage
  or effective taper counts into integrated covariance inference.
- [x] Verify 180 randomized NumPy-oracle scenarios / 1,416 band/time cells,
  native checked/optimized builds and the public `scifort_signal` API.
- [ ] Investigate uncertainty calibration for band-integrated conditional
  coherence under correlated frequency bins and adaptive taper weights.

## 2026-10-29: Spectral multiple coherence and conditional power

- [x] Add squared multiple coherence, residual conditional spectral power,
  and explained power for every channel jointly predicted by all others.
- [x] Support real/complex static and time-frequency DPSS spectral matrices,
  existing Cholesky validity masks, trace ridge, explicit and fixed-weight
  automatic shrinkage, and leave-one-out adaptive shrinkage.
- [x] Compare 700 five-channel NumPy Schur complements in checked and
  optimized GNU Fortran builds, validate direct signal entry points, and
  compile the 93-module public API closure with six public regression tests.
- [x] Add band-integrated multivariate coherence (see section above); investigate more informative
  uncertainty diagnostics only after a justified statistical calibration.

## 2026-10-29: Data-driven diagonal spectral shrinkage

- [x] Add Gaussian-Wishart moment plug-in shrinkage estimated per spectral
  frequency and time-frequency cell from the unregularized CSD matrix.
- [x] Preserve existing explicit shrinkage and optional ridge semantics,
  expose used coefficients and validate effective taper count.
- [x] Add NumPy precision/coherency comparisons, native checked/-O3
  tests, public API smoke tests and documentation of distributional limits.
- [x] Extend automatic spectral shrinkage to deterministic eigenvalue-weighted
  DPSS tapers using effective sample counts. Test real/complex static and
  time-frequency cases and expose fixed-weight metadata.
- [x] Add optional leave-one-taper-out predictive cross-validation of
  diagonal shrinkage from individual Fourier vectors for equal, eigen and
  adaptive DPSS weights; refit adaptive weights on remaining tapers and
  test against independent SciPy/NumPy calculations.
- [ ] Develop and validate inference or oracle-risk calibrations before
  reporting any statistical confidence intervals for shrinkage.

## 2026-10-29: Phase-aware conditional spectra and diagonal shrinkage

- [x] Preserve the phase of conditional cross spectra through opt-in complex
  partial coherency, available for static and moving-window estimates.
- [x] Implement diagonal-target shrinkage prior to existing trace ridge,
  retaining per-frequency/time singularity and numerical validation.
- [x] Add independent NumPy/SciPy comparisons for real/complex signals,
  frequency and time-frequency matrices, and prior-regression coverage.
- [x] Add a separate opt-in Gaussian-Wishart automatic shrinkage mode with
  explicit assumptions; do not claim calibrated confidence intervals.

## 2026-10-29: Time-frequency partial spectral coherence

- [x] Extend conditional spectral precision calculations to moving-window
  cross-spectral matrix spectrograms, retaining all existing FFT, DPSS,
  detrending and weighting options.
- [x] Offer direct-matrix and real/complex signal interfaces with time axes,
  optional trace-scaled ridge regularization and per-cell validity flags.
- [x] Verify 120 independent NumPy/SciPy cases and retain prior multitaper
  and 252-case partial-coherence numerical regressions.
- [ ] Consider online/streaming conditional coherence and uncertainty
  estimates if supported by a clear numerical/statistical use case.

## 2026-10-29: Multichannel conditional spectral analysis

- [x] Add inverse-spectral partial coherence from Hermitian multitaper
  cross-spectral matrices; use a stable Cholesky solve per frequency.
- [x] Support real and complex multichannel input, relative diagonal loading,
  precision matrices, per-frequency singularity masks, and public export.
- [x] Compare NumPy precision inverses and SciPy DPSS/FFT spectral matrices
  in checked and targeted optimized builds, retaining previous multitaper
  regressions and the public module's API.
- [ ] Consider partial coherence for streaming/time-frequency spectral matrices
  and statistical uncertainty estimates if practical use warrants them.

## 2026-10-29: Prepared weighted multiharmonic plans

- [x] Cache QR factors and weighted Fourier design columns for arbitrary
  irregular sample times, weights, harmonic counts, and frequency grids.
- [x] Fit multiple channels efficiently from the same plan, retaining
  per-channel statistics, covariance, and fitted-mean uncertainty functions.
- [x] Test mixed valid/singular frequency bins, zero-weight samples, strided
  arrays, malformed inputs, checked and optimized builds, and public export.
- [ ] Consider streaming large channel banks, plan serialization, or further
  numerical algorithms based on measured user needs. Very large frequency
  grids should continue to use unprepared fits to control cached memory.

## 2026-10-29: Harmonic fit covariance and curve uncertainty

- [x] Estimate parameter covariance and standard errors from weighted QR,
  with floating-intercept correlations and zero-weight residual dof.
- [x] Add fitted-mean prediction standard errors, with well-defined failure
  behavior and allocation-free defaults when uncertainty is not requested.
- [x] Verify 256 independent NumPy weighted-regression references and prior
  harmonic regressions in checked and optimized builds.
- [x] Prepared weighted multiharmonic QR plans now amortize construction
  over multiple channels; see `docs/signal_harmonic_plan.md`.

## 2026-10-28: Lomb–Scargle harmonic model extraction

- [x] Add `lombscargle_harmonic_fit` for cosine/sine amplitudes, intercept,
  per-bin normalized/unnormalized power, fit RMS, and rank validity.
- [x] Add `lombscargle_harmonic_predict` for waveform reconstruction and
  evaluation at arbitrary times; restore the precentering shift in outputs.
- [x] Independently verify 256 weighted NumPy least-squares fits, SciPy
  one-harmonic results, the previous periodogram, and the public-module API.
- [x] A reusable multiharmonic plan now handles repeated-channel fits;
  coefficient covariance and prediction uncertainties are also implemented.

## 2026-10-28: Multiharmonic irregular-sampling periodogram

- [x] Add a joint weighted Fourier-series Lomb–Scargle algorithm with
  configurable harmonic count, optional weighted floating mean, and normalized
  or unnormalized power.
- [x] Preserve the existing one-harmonic API and numerical behavior; export the
  new routine through the top-level `scifort_signal` module.
- [x] Test independent NumPy least squares and SciPy's one-harmonic reference,
  invalid data, degenerate frequencies, noncontiguous arrays, and both builds.
- [x] Independent multichannel fits now share harmonic design QR factors
  through a reusable prepared plan; continue broader numerical algorithms.

## 2026-10-28: Arbitrary-input SOS simulation

- [x] Add stable SOS-native `dlsim_sos` and `dlsim(sos,u,dt,...)` for SISO
  signals without polynomial conversion, including mixed real/complex types.
- [x] Return per-sample DF-II state histories, final states for streaming,
  and support interpolation from irregular input times.
- [x] Verify 240 SciPy comparison cases in checked and targeted optimized
  builds; preserve and rerun 288 earlier SOS response comparisons.
- [ ] Investigate multichannel SOS filter banks with independently specified
  sections per channel, or other missing numerical signal algorithms.

## 2026-10-28: Stable SOS discrete-time responses

- [x] Implement `dstep_sos` / `dimpulse_sos` for real and complex cascades
  and overload the existing `dstep` / `dimpulse` public generics.
- [x] Preserve per-section DF-II `zi` initial states and avoid high-order
  coefficient expansion; validate defaults, errors, and noncontiguous arrays.
- [x] Check 288 independent SciPy cases (15,153 samples), native status tests,
  checked and optimized compilers, and public API compilation.
- [ ] Extend SOS-native time-domain simulations to multichannel models or
  further missing algorithms, without expanding maximum array ranks.

## 2026-10-28: Root-wise ZPK group-delay analysis

- [x] Add a stable group-delay algorithm for digital ZPK filters, including
  mixed real/complex zeros, poles and gains, with no polynomial expansion.
- [x] Handle unequal root counts, whole/half/custom grids, positive `fs`,
  unit-circle singularities, and public `group_delay` overloading.
- [x] Cross-check 221 SciPy/reference cases, including long root arrays,
  checked and optimized GNU Fortran builds, existing group-delay regressions,
  and public 91-source compilation.
- [ ] Continue substantive signal-processing algorithms and filtering
  improvements rather than expanding maximum array rank.

## 2026-10-28: Stable SOS group delay

- [x] Add real/complex second-order-section group delay without forming
  ill-conditioned high-order numerator/denominator polynomials.
- [x] Expose `group_delay_sos` and rank-two `group_delay` overloads through
  `scifort_signal`; preserve existing vector transfer-function overloads.
- [x] Check singular-bin NaNs, invalid SOS/FFT grids, `fs` units,
  high-order cascades and separately scaled numerator coefficients.
- [x] Verify 124 SciPy/reference cases (2,612 frequencies, up to 128
  sections) with checked and `-O3` builds; preserve old group-delay tests.
- [ ] Continue substantive numerical signal analysis and filtering gaps.

## 2026-10-28: Mixed-type rank-one convolution

- [x] Add both real/complex operand orders to 1D `convolve`, `correlate`,
  `fftconvolve`, `oaconvolve`, and `choose_conv_method`, preserving the
  established direct/FFT/overlap-add kernels and public positional APIs.
- [x] Check conjugation in mixed cross-correlation, all convolution modes,
  valid/invalid axes, zero-length inputs, strided sections, status handling,
  and measured convolution algorithm selection.
- [x] Independently verify 960 SciPy 1.17 cases / 182,304 complex values in
  checked and targeted optimized GNU Fortran builds (worst absolute error
  below `5.2e-13`); retain same-type regression tests.
- [ ] Continue missing signal-processing algorithms and substantive options
  rather than adding further rank overloads solely for completeness.

## 2026-10-28: ShortTimeFFT callable detrending

- [x] Add procedure-based real/complex `ShortTimeFFT.stft_detrend` for
  one-dimensional, rank-2 and rank-3 signals and existing signal axes.
- [x] Match SciPy's per-slice callback ordering (after padding, before
  windowing) and reject NaN/Inf callback outputs without partial results.
- [x] Add independently generated SciPy coefficients and frequency/time-grid
  comparisons for FFT modes, phase shifts, window scaling and padding;
  retain previous named-detrend behavior and high-rank source interfaces.
- [ ] Continue substantive algorithms and filtering/time-frequency options;
  do not increase supported array rank solely for completeness.

## 2026-10-28: Mixed spectral estimators

- [x] Add rank-one `csd` and `coherence` for real-complex and complex-real
  operand pairs without changing same-type generic resolution.
- [x] Add rank-one `short_time_fft%cross_spectrogram` for mixed inputs with
  both two-sided and centered FFT modes, and explicit invalid-mode rejection.
- [x] Validate the independent implementations against 384 SciPy comparisons,
  plus native callback, strided-input and legacy regression tests.
- [ ] Investigate mixed-type higher-rank `short_time_fft%cross_spectrogram`
  only if demand justifies additional generic overloads; prioritize genuinely
  missing algorithms and signal-processing options instead.

## 2026-10-28: Complex-valued Wiener denoising

- [x] Extend `wiener` for complex rank-1/2/3 arrays, real/complex explicit
  noise, automatic noise estimation, even-sized windows, and strided inputs.
- [x] Match SciPy's complex pseudo-variance convention and validate 27 new
  SciPy reference cases alongside 1,050 previous denoising comparisons in
  checked and optimized builds.
- [ ] Continue filling numerical filtering and time-frequency algorithm
  gaps rather than raising supported ranks further.

## 2026-10-28: Per-segment spectral callbacks

- [x] Provide real and complex in-place callback interfaces for all six
  legacy spectral estimator families and their multidimensional `*_into`
  adapters, without changing existing positional or string-detrend APIs.
- [x] Validate user-provided transformations, mutually exclusive arguments,
  nonfinite results, strided/rank-15 data and unchanged output on failure.
- [x] Verify 48 SciPy 1.17 reference cases and the older spectrogram and
  two-sided spectral checked tests; compile public 91-source dependency closure.
- [ ] Continue substantive filtering, spectral-analysis and time-frequency
  algorithm coverage; defer further expansion of array-rank limits.

## 2026-10-28: Shaped Fourier-domain resampling windows

- [x] Accept Gaussian, Kaiser, Tukey, generalized Gaussian, Taylor, DPSS
  and general-cosine parameters in FFT-domain `resample` for real/complex
  data, rank-2/3 arrays and `resample_into` through rank 15.
- [x] Preserve frequency-bin callback and explicit-weight alternatives,
  reject invalid/competing specifications and keep outputs unchanged on error.
- [x] Check 252 SciPy cases / 3,696 values in checked and targeted `-O3`
  builds, plus native strided/high-rank and previous multirate regressions.
- [ ] Continue substantive signal-processing algorithm/option coverage.

## 2026-10-28: Complex second-order-section decimation

- [x] Support complex `sos=` coefficients through `decimate` and
  `decimate_into`, reusing the complex SOS cascade kernels and preserving
  existing real-coefficient and transfer-function calls.
- [x] Support real-to-complex output, rank-2/3 axis selection, rank-15
  assumed-rank inputs, noncontiguous slices and unchanged output on failure.
- [x] Compare 288 seeded SciPy cases / 7,434 complex samples and rerun
  earlier decimation/multirate tests in checked and targeted `-O3` builds.
- [ ] Continue filling signal-processing algorithm gaps, especially
  higher-level filter structures and spectral-analysis options.

## 2026-10-28: Real and complex transfer-function decimation

- [x] Extend `decimate` and `decimate_into` with `b=`/`a=` while keeping
  explicit FIR/SOS and designed filters backward compatible.
- [x] Implement complex transfer-function IIR decimation via existing
  `lfilter`/`filtfilt`, scalar-denominator FIR normalization, mixed
  real-to-complex output and rank-2/3 plus assumed-rank-through-15 adapters.
- [x] Verify 576 SciPy 1.17 cases, 16,320 samples (max error `8.951e-16`),
  checked/targeted `-O3` builds and old multirate regression suites.
- [ ] Explore stable direct complex SOS decimation and more signal-analysis
  algorithm modes; avoid unnecessary array-rank expansion.

## 2026-10-28: Caller-supplied FIR and SOS decimation

- [x] Allow caller-provided real FIR `taps=` and normalized real IIR `sos=`
  in real/complex `decimate` (rank 1-3) and `decimate_into` (rank 1-15),
  while preserving existing designed-filter options and positional calls.
- [x] Match custom FIR delay compensation and causal output, and custom IIR
  causal/zero-phase filtering against SciPy `dlti` reference filters.
- [x] Run 336 independent real/complex SciPy 1.17 comparisons, 10,500 samples,
  new checked and targeted optimized native tests, and previous 53+184
  multirate regression cases. Reject invalid options without changing
  caller-allocated outputs.
- [x] Support arbitrary real and complex transfer-function coefficient arrays
  (`b=`, `a=`) including real-to-complex IIR, numeric `dlti` equivalents,
  causal/zero-phase modes, and all existing rank-aware interfaces.
- [ ] Keep extending substantive signal-processing algorithm coverage.

## 2026-10-28: Parameterized polyphase-resampling prototypes

- [x] Add Gaussian, generalized Gaussian, Kaiser, Tukey, Taylor and
  general-cosine shape arguments to real/complex `resample_poly`, including
  multidimensional array-returning and rank-independent `*_into` interfaces.
- [x] Preserve existing positional APIs and legacy Kaiser `beta`, validate
  conflicting window/tap arguments, and preserve caller output on failures.
- [x] Verify 256 SciPy reference cases, targeted checked and optimized
  Fortran builds, and earlier multirate/rank-15 independent regressions.
- [ ] Continue filling substantive missing algorithm modes and filtering
  capabilities rather than increasing supported tensor rank.
## 2026-10-28: Complete legacy spectrogram mode coverage

- [x] Add `mode='complex'`, `'magnitude'`, `'angle'`, and `'phase'` to the
  real and complex one-dimensional `spectrogram` functions while preserving
  the existing `mode='psd'` default and typed result components.
- [x] Reuse the same implementation for existing rank-independent
  `spectrogram_into` modes, eliminating duplicated STFT setup/phase logic.
- [x] Check 440 independent SciPy 1.17 configurations / 58,254 coefficients
  in checked and optimized GNU Fortran 14.2, with new native error,
  signed-zero, strided-array and multidimensional tests. An independent
  396-case (727,408 values) multidimensional SciPy comparison and the full
  91-module public build also pass in checked mode.
- [ ] Continue filling substantive signal-processing algorithm/mode gaps
  rather than increasing Fortran rank limits further.

## 2026-10-27: FIR scalar inputs and arbitrary multiband parity

- [x] Accept scalar or array cutoff frequencies and both logical and string
  `pass_zero` through the public `firwin` generic, without breaking existing
  positional interfaces.
- [x] Handle three/five/other odd-count band edges; derive Nyquist pass/stop
  selection from cutoff parity and DC mode, preserving SciPy's even-tap
  restrictions and first-passband gain scaling.
- [x] Compare 443 valid FIR designs/10,370 coefficients and 362 invalid
  designs against SciPy 1.17 in checked and optimized builds. Retest FIR
  transition-width and window-shape functionality.
- [ ] Continue extending filter design and signal analysis algorithm options.

## 2026-10-27: FIR transition-band design and parameterized tapers

- [x] Implement SciPy-compatible positive `width` for `firwin` by deriving
  Kaiser beta from the transition width expressed in physical sample-rate
  units. Preserve existing filter families, coefficient normalization and
  positional APIs; ignore the named window when `width` is specified.
- [x] Extend `firwin2` with named symmetric window parameters, two-parameter
  families, general-cosine coefficients and direct finite symmetric tapers.
  Reuse common spectral window-family validation and reject conflicts.
- [x] Verify 320 SciPy 1.17 comparisons (6,736 coefficients) in checked and
  optimized builds, plus the previous 1,620 design, 1,070 advanced-FIR and
  192 two-sided spectral comparisons; document SciPy's current
  `firwin2(general_cosine)` reference limitation.
- [ ] Continue extending remaining `scipy.signal` filter-design and signal
  analysis modes, including multichannel efficiency improvements.

## 2026-10-27: Real-input two-sided legacy spectral analysis

- [x] Implement SciPy-style `return_onesided=.false.` for real
  `periodogram`, `welch`, `csd`, `stft`, and `spectrogram`. Use the full FFT
  frequency grid (including negative Nyquist for even `nfft`) without
  one-sided PSD doubling, and keep one-sided output as the default.
- [x] Thread the optional setting through existing rank-independent
  `*_into` adapters, including all supported spectrogram output modes, while
  preserving positional compatibility and transactional output semantics.
- [x] Validate 192 deterministic SciPy 1.17 reference cases (7,248
  coefficients) and native strided/rank-14/rank-15 paths with both checked
  and optimized GNU Fortran; rerun previous spectral and inverse regressions.
- [ ] Continue filling missing `scipy.signal` options and algorithms,
  including multi-channel performance, without increasing tensor-rank limits.

## 2026-10-27: Parameterized inverse STFT and overlap checks

- [x] Extend legacy real and complex `istft` and the rank-independent
  `istft_into` interface to accept the shape parameters and general cosine
  coefficients already accepted by forward STFT. Share window-family
  validation and preserve positional call compatibility.
- [x] Extend the named-window overloads of `check_cola` and `check_nola`
  with shape parameters. Reject unused/mismatched options; keep explicit
  window arrays available as before.
- [x] Verify real/complex, rank-14, strided and invalid-option recovery
  natively and compare 192 independently generated SciPy 1.17 STFT inputs
  (5,613 reconstructed samples) with real/complex Fortran inverse output.
  Checked and optimized GNU Fortran builds pass.
- [ ] Extend remaining `scipy.signal` algorithm modes and continue improving
  efficiency on repeated/multichannel transforms.

## 2026-10-27: Parameterized and explicit windows for legacy spectra

- [x] Propagate named-window shape parameters, two-parameter families, general
  cosine coefficients, and caller-provided real window samples through the
  six legacy estimators (`periodogram`, `welch`, `csd`, `coherence`, `stft`,
  and `spectrogram`) and their assumed-rank `*_into` overloads for real and
  complex inputs. Preserve earlier positional signatures.
- [x] Infer segment length from explicit windows and reject mismatched
  lengths, conflicting window requests, invalid parameter families, and
  nonfinite window weights with predictable status codes.
- [x] Correct the omitted-window spectrogram default to SciPy's Tukey
  alpha=0.25, while leaving explicit Tukey's default at alpha=0.5.
- [x] Check 96 independent SciPy 1.17 cases (3,264 spectral coefficients),
  validate rank-14/rank-15 and strided inputs, preserve old STFT boundary
  tests, and run checked and O3 standalone builds.
- [ ] Continue adding missing `scipy.signal` algorithms and improve spectral
  performance, especially for many channels and repeated window shapes.

## 2026-10-27: Legacy STFT boundary-mode completion

- [x] Implement SciPy-style `even`, `odd`, `constant`, `zeros` and `none`
  endpoint conventions for real and complex 1-D `stft` and multidimensional
  `stft_into`, without changing the established public signatures.
- [x] Compare 1,080 deterministic cases and 54,993 STFT coefficients against
  SciPy 1.17, including odd/even FFT lengths, full and partial final segments,
  real/complex spectra, and spectrum/PSD scaling; test strided/rank-14 inputs.
- [x] Implement configurable window parameters for legacy spectral APIs
  (see the 2026-10-27 parameterized-window work above).
- [ ] Continue expanding multidimensional modes and performance.

## 2026-10-27: Multidimensional Hilbert analytic signals

- [x] Add `hilbert_into` for arbitrary-rank one-axis analytic signals and
  `hilbert2_into` for two arbitrary distinct axes, through rank 15.
  Preserve one-based/negative axes, optional FFT length changes, noncontiguous
  arrays, and transactional output updates.
- [x] Reuse explicit FFT plans and the established SciPy 1.17 frequency-mask
  conventions, including the special even-length two-dimensional Nyquist rule.
- [x] Check 500 independently generated SciPy reference cases and 191,284
  complex coefficients, and native checked/optimized public-module tests.
- [ ] Consider array-returning convenience overloads only if they can avoid
  excessive module/template duplication; continue filling `scipy.signal`
  method/algorithm gaps.

## 2026-10-27: Arbitrary-rank relative-extrema masks

- [x] Add `argrelmax_mask_into`, `argrelmin_mask_into`, and
  `argrelextrema_mask_into` for ranks 1-15 with selected one-based/negative
  axes, SciPy-equivalent `order` and `clip`/`wrap` comparison semantics, and
  immutable input / transactional logical mask outputs.
- [x] Verify 720 independent SciPy reference cases (including NaNs, ties,
  singleton and zero-size dimensions, custom non-strict comparators) and native
  public-module tests for rank 15, noncontiguous arrays and error recovery.
- [ ] Consider direct coordinate-tuple output or streamed extrema events for
  very large arrays; logical masks require O(product(shape(x))) workspace.

## 2026-10-27: Multichannel irregular-sampling spectra

- [x] Add caller-output Lomb-Scargle power and complex-amplitude routines
  for ranks 1-15, with shared sample times and weights, selected sample
  axis, and all-or-nothing output updates.
- [x] Validate against SciPy 1.17 and the original rank-one implementation;
  exercise rank 15, negative axes, noncontiguous arrays, and error recovery.
- [x] Add per-channel `weights_nd` and `nan_policy='omit'` for irregular
  missing data, with atomic outputs and SciPy comparisons after masking.
- [x] Share frequency-dependent trigonometric terms across all channels
  for common weights and benchmark against independent-lane evaluation.
- [x] Add reusable `lombscargle_plan` for repeated observations with the
  same time grid, frequency grid, and shared weights, plus separate power
  and complex-amplitude adapters. Cache weighted trigonometric columns
  and Gram denominators, validate in SciPy comparison tests, and benchmark
  repeated evaluations.
- [ ] Consider an optional missing-data mask shared by all channels, allowing
  reuse of a plan on a common subset of observations; independent channel
  masks and per-channel weights still require separate weighted fits.

## 2026-10-26: Fourier resampling with caller-supplied FFT weights

- [x] Extend `resample` (ranks 1-3) and `resample_into` (ranks 1-15) to
  accept FFT-order real spectral weights and, for complex inputs, complex
  spectral weights. Preserve all existing positional and named-window calls.
- [x] Match SciPy's Hermitian frequency-weight folding for real inputs;
  preserve unchanged outputs after errors in the assumed-rank adapter.
- [x] Verify 900 independent SciPy 1.17 comparisons including asymmetric
  windows, complex weights, original named windows, `domain='freq'`, and
  uneven down/up sampling. Check rank-15 and strided sections in Fortran.
- [x] Support SciPy-style callable real and complex frequency windows with
  explicit Fortran vector callbacks on ranks 1-3 and 1-15 (`resample_into`).
  Test even/odd FFT frequencies, one-time evaluation per multichannel call,
  rank 15, input validation, and 450 independent SciPy 1.17 cases.
- [x] Accelerate prime and other non-power-of-two FFT lengths using
  Bluestein convolution in the shared spectral FFT kernel. Check 116
  independent SciPy FFT cases and 24 direct native regressions; benchmark
  lengths 257, 1009, and 4093 with GNU Fortran `-O3`.
- [x] Introduce reusable `fft_plan` with chirp, transformed convolution
  kernel and workspace; integrate it into batch spectral estimation and
  streaming STFT/inverse STFT. Verify 198 SciPy FFT cases and native tests.
- [x] Implement mixed-radix Cooley-Tukey FFTs with factors 2, 3, 5 and 7
  for composite lengths, including reusable plan roots and scratch arrays.
  Compare 246 transforms with SciPy and check existing public FFT regressions.
- [x] Reuse plans across `fftconvolve` (ranks 1-3), `fftconvolve_nd`
  and tiled `oaconvolve_nd` (ranks 1-15), whole-record/windowed/harmonic
  multitaper estimation, and the persistent `multitaper_stream` state.
  Add prime/smooth-length tests, reset/reinitialization tests, and benchmarks.
- [x] Introduce cost-aware N-dimensional FFT convolution padding, comparing
  supported 2/3/5/7-smooth shapes with radix-2 shapes. Retain radix-2 on
  workloads where the recursive mixed-radix implementation is slower.
  Validate full/same/valid, rank-15, complex and overlap-add paths and
  benchmark unchanged and reduced-padding workloads.
- [x] Replace the generic mixed-radix output summation with radix-2, -3,
  -5 and -7 butterflies. Compare against a direct DFT, independent SciPy
  FFTs, and the previous implementation's same-length `-O3` timings.
- [ ] Consider radix-2 twiddle caching, cache-friendly radix ordering,
  vectorization, calibrated or workload-specific convolution cost models,
  and further FFT-heavy callers. Resume missing `scipy.signal` methods.

## 2026-10-26: Windowed and streaming local multitaper noise floors

- [x] Reuse the same median-frequency boundary and zero-background conventions
  in `multitaper_noise_floor_spectrogram` for rank-1 and multichannel signals.
- [x] Add `multitaper_stream%push_noise_floor` and `mean_noise_floor` without
  recording any new window history. Check invalid widths before consuming data.
- [x] Validate against independent SciPy 1.17 DPSS/detrend/FFT/median
  calculations and both debug and optimized public-module regressions.
- [ ] Examine colored-noise harmonic-line inference under locally correlated
  bins with calibrations distinct from the descriptive noise-floor ratio.
- [ ] Resume missing `scipy.signal` algorithms and modes.

## 2026-10-26: Multitaper colored-noise background and Holm correction

- [x] Add local-median `multitaper_noise_floor` with one-sided reflected and
  two-sided periodic spectral boundaries, colored-AR(1) reference data and
  explicit non-statistical interpretation of background-normalized power.
- [x] Add multiplicity-robust `method='holm'` to multitaper harmonic scans;
  retain BH, BY, Bonferroni and unadjusted alternatives. Validate against
  Statsmodels on real/complex Fourier tests and selected bins.
- [ ] Explore formal spectrum-dependent colored-noise null calibration for
  harmonic tests, including adjacent-bin correlations and finite taper count.
- [ ] Continue missing `scipy.signal` algorithms once spectral diagnostics
  are adequately tested.

# Roadmap

## Rank-independent signal filtering

- [x] Refactor existing rank-1 through rank-4 `savgol_filter` functions to
  use one shared axis-wise kernel and introduce real/complex Fortran 2018
  assumed-rank `savgol_filter_nd` with caller-supplied rank-1 through rank-15
  output, explicit one-based axis selection and rank/shape validation.
- [x] Add `lfilter_nd` using Fortran 2018 assumed-rank real/complex arrays
  (ranks 1-15), one-based axis selection, caller-supplied output and states,
  and six mixed-type promotions; preserve the rank-1/2/3 `lfilter` functions.
  Validate 300 seeded rank-5/15 cases against SciPy, including final states.
- [x] Add rank-independent `sosfilt_nd` for real/complex and mixed signals,
  preserving the leading SOS-section state dimension for ranks 1-14, with
  packed rank-one states when rank 15 requires a conceptual rank-16 state;
  preserve existing rank-1/2/3 `sosfilt` interfaces and validate 320 seeded
  SciPy outputs and streaming final states.
- [x] Add caller-output, assumed-rank `filtfilt_nd` and `sosfiltfilt_nd`
  for ranks 1-15, independent one-based axes, all real/complex input and
  coefficient combinations, and SciPy-style padding. Keep existing
  array-returning rank-1/2/3 functions unchanged. Verify 684 SciPy 1.17
  comparisons, zero-channel behavior, strided sections and failures.
- [x] Add `fftconvolve_nd` for ranks 1-15 with caller-allocated output,
  selected convolution axes, singleton broadcasting, SciPy `full`/`same`/
  `valid` modes, and all real/complex operand combinations. Preserve the
  rank-1/2/3 `fftconvolve` functions; compare 576 rank-15 cases to SciPy 1.17.
- [x] Add caller-output `convolve_nd` and `correlate_nd` for ranks 1-15,
  implementing actual direct N-dimensional summation and conjugate-reversed
  correlation kernels alongside FFT dispatch. Support `full`/`same`/`valid`,
  selected axes, singleton broadcasting, mixed real/complex inputs, output
  shape discovery, and transactional error handling. Compare 1008 rank-15
  randomized cases (67,609 elements) with SciPy 1.17 in both checked and
  optimized builds; pass 10 neighboring checked convolution regressions.
- [ ] Tune automatic direct-versus-FFT selection with compiler-specific
  timings; consider memory-bounded N-D direct correlation for large kernels.
- [x] Add caller-output `oaconvolve_nd` for ranks 1-15, sharing
  `fftconvolve_nd` axis/broadcast/shape rules and four real/complex type
  combinations. Use blocked overlap-add on the most asymmetric selected
  dimension, with separable FFTs for the remaining convolution axes and
  output cropping performed during accumulation. Validate 400 seeded
  SciPy 1.17 reference cases and deterministic Fortran regressions; keep
  existing one-dimensional `oaconvolve` overloads unchanged.
- [x] Reuse the transformed kernel across overlap-add blocks and tile
  several selected axes when one operand is much longer along each axis.
  Preserve mixed-orientation correctness by tiling one operand per call,
  and add SciPy comparisons, optimization checks and before/after timings.
- [ ] Investigate portable alternatives to C-interoperable storage remapping
  for compilers with unusual double-precision kinds, and test across GNU,
  Intel, NAG and LLVM Fortran compilers.


## Additional SciPy signal algorithms

- [x] Extend legacy inverse STFT to `istft_into`, reconstructing
  rank-1..14 real/complex signals from rank-2..15 STFTs with configurable
  frequency axis and appended time, complete-output error semantics and
  240 independently compared SciPy 1.17 cases. See `docs/signal_timefreq_nd.md`.
- [x] Allow `istft_into` to select independent frequency and time axes in
  any order, preserve output-axis compatibility, and verify 1,296 SciPy 1.17
  inverse cases, including rank-15 and strided-array regressions.
- [x] Add single-channel streaming inverse STFT with fixed-size overlap-add
  buffers and streaming forward STFT/spectrogram with fixed-size input buffers,
  flexible chunk boundaries and real/complex transforms. Verify end-to-end
  reconstruction and independent SciPy comparisons.
- [x] Add synchronized real/complex multichannel streaming STFT, power
  spectrogram and inverse STFT interfaces with sample/channel chunk input,
  frequency/channel/frame output and bounded state. Verify 256 independent
  SciPy 1.17 comparisons, noncontiguous sections and streaming round trips.
  See `docs/signal_stream_multi.md`.
- [ ] Extend streaming interfaces beyond fixed sample/channel and
  frequency/channel axes, and investigate additional boundary modes,
  cross-channel processing and parameterized window defaults.

- [x] Add `periodogram_into`, `welch_into`, `csd_into`, and `coherence_into`
  for real/complex rank-1–15 arrays, mixed cross-spectra, selected axes,
  and SciPy-compatible Welch median averaging. Confirm 396 numerical
  reference cases (446,736 values) and public/strided-array regressions.

- [x] Add complex custom FIR taps to `resample_poly` (rank 1-3),
  plus median/minimum/maximum padding and antisymmetric/antireflect
  boundary modes to polyphase filtering. Validate against SciPy.
- [x] Generalize the multirate resampling APIs to rank 15 using
  assumed-rank caller-provided outputs, following `lfilter_nd`.

- [x] Add `whittaker_henderson` with fixed lambda, weights, and a banded
  solver; add `normalize`, `abcd_normalize`, and direct/FFT convolution-method
  selection with explicit `method='fft'` and automatically selected methods.
- [x] Add SciPy-style automatic REML smoothing-penalty selection to
  Whittaker-Henderson, using a banded Cholesky profiled-likelihood criterion.
- [ ] Investigate higher-order conditioning and sparse/masked constraints,
  and adaptive log-penalty bounds when the optimum lies beyond the fixed
  REML search interval.
- [ ] Extend method selection/dispatch to rank-independent convolution
  inputs and benchmark compiler-specific switching thresholds.

- [x] Implement `closest_stft_dual_window`, `firwin_2d`, `symiirorder1`,
  `symiirorder2`, `findfreqs`, `unique_roots`, and `czt_points`; check 499
  SciPy 1.17 reference cases and 11,610 numerical values. Reuse the existing
  spline boundary recurrences and analog frequency/root methods.
- [x] Extend ShortTimeFFT forward/inverse/spectrogram operations to real and
  complex rank-two and rank-three signals with selectable axes.
- [x] Extend the ShortTimeFFT derived type to rank-four and rank-five
  real/complex inputs, including forward/inverse, detrending, and spectra.
- [x] Extend `ShortTimeFFT` to rank-six through rank-fourteen input signals
  (rank-15 output, the Fortran limit), including inverse, detrending and
  (cross-)spectrograms. Add coordinate, slice-range, and boundary helpers.
  Checked against SciPy for 216 high-rank cases and 216 coordinate cases.
- [ ] Consider a flattened or caller-output interface for rank-15 input
  signals, since conventional STFT output would require rank 16.
- [ ] Add reusable CZT plans and remaining STFT scaling conveniences.
- [x] Extend `firwin_2d` with independent second-axis cutoff bands, sampling
  rates, passband type and normalization, plus rotated elliptical radial
  interpolation as SciFort extensions. Check 239 independent SciPy/NumPy
  references, invalid-argument handling, and unchanged SciPy-style defaults.
- [ ] Explore alternative 2-D FIR design algorithms with direct frequency-
  response optimization rather than coefficient-domain radial interpolation.


- [x] Add `cspline2d`, `qspline2d` and `spline_filter` with symmetric IIR
  boundary initialization, smoothed cubic coefficients, real/complex
  Fortran interfaces, and 44 independent SciPy numerical comparisons.
- [x] Extend spline interpolation with explicit 2-D query-coordinate
  evaluation for cubic/quadratic real and complex coefficient arrays, supporting
  scalar, paired vector/matrix, and outer-grid queries with mirror boundaries.
  Verify 256 SciPy comparisons (17,920 coordinates) and public regressions.
- [x] Add rank-1 through rank-15 assumed-rank spline filtering and coordinate
  interpolation, orders 0-5, real/complex data, selected-axis filtering,
  mirror/constant boundaries, optional prefiltering, and 288 SciPy comparisons.
  See `docs/signal_spline_nd.md`.
- [x] Expand `map_coordinates_nd` to all eight SciPy boundary modes; add
  assumed-rank, real/complex `affine_transform_nd` and `shift_nd` for ranks
  1-15. Confirm 528 new SciPy interpolation/affine/shift comparisons and
  preserve the 288-case prior spline reference. Document the short-array
  reflect-prefilter numeric differences.
- [x] Complete order-2..5 nearest/grid-constant prefiltering with 12-sample
  halos, singleton-axis optimization, explicit allocation guard, and SciPy
  interpolation/affine/zoom/rotate comparisons.
- [ ] Speed up periodic spline prefiltering and profile tensor-product
  interpolation in dimensions above five.
- [x] Add assumed-rank `zoom_nd` and `rotate_nd` for real/complex arrays,
  with SciPy-compatible shape helpers, grid-mode zoom, arbitrary rotation
  planes, 1-15 and 2-15 rank coverage, and 1,156 SciPy reference cases.
- [ ] Improve short-array reflect recursive-boundary parity.
- [x] Add `scifort_ndimage` real/complex rank-1..15 Gaussian smoothing and
  derivatives, Gaussian gradient magnitude/Laplacian, rectangular mean/min/max,
  Sobel and discrete Laplacian; verify against SciPy in eight boundary modes.
- [x] Add rank-1..15 logical binary dilation/erosion and connected-component
  labeling with configurable Manhattan connectivity and SciPy comparisons.
- [x] Add custom image footprints, rank-independent median and percentile
  filtering, EDT/chamfer distance transforms and watershed segmentation.
- [x] Add real generic neighborhood/labeled callbacks, 1-D real/complex
  correlation/convolution and weighted taxicab/custom directed chamfer paths.
- [ ] Add per-axis mixed boundary modes and complex generic callbacks.

- [x] Implement `gauss_spline`, `kaiser_beta`, `kaiser_atten`, `kaiserord`,
  and parameterized Kaiser support for `firwin`. Validate 2,388 seeded
  numerical cases against SciPy and both 102-test standalone signal suites.
- [x] Extend B-spline interpolation coverage with `cspline1d`,
  `qspline1d`, `cspline1d_eval`, and `qspline1d_eval`, including real and
  complex data, mirror-symmetric boundaries, nonunit knot spacing and
  valid cubic smoothing parameters. Validate 480 SciPy reference cases.


- [x] Implement auditory `gammatone` FIR/IIR designs and periodic `iircomb`
  notch/peak filters, preserving SciFort's real transfer-function API.
  Validate 1,350 randomized cases against SciPy 1.17 and test malformed
  frequencies, comb divisibility, DC behavior, and invalid designs.
- [ ] Extend auditory filter-bank support with stable cascaded gammatone
  sections and multiband auditory feature extraction as SciFort extensions.


- [x] Expand `scipy.signal.get_window` compatibility to twenty named
  families, with SciPy-style periodic/symmetric conventions and optional
  shape parameters for Kaiser, Gaussian, Tukey, generalized Gaussian,
  generalized Hamming, and exponential windows. Verify 1,092 seeded/window
  cases against SciPy 1.17 (max error 1.43e-14), plus native tests.
- [x] Add Dolph-Chebyshev `chebwin`, Taylor and Kaiser-Bessel-derived windows,
  with direct/public interfaces, status validation and 1,718 SciPy comparisons.
- [x] Add coefficient-vector `general_cosine` (direct and `get_window`)
  and integrate the existing DPSS single taper into `get_window` with the
  approximate SciPy normalization. Verify 1,092 seeded SciPy 1.17 cases
  (maximum absolute error 2.64e-12) and 100/100 checked/optimized signal
  regressions. Maintain existing `get_window` positional arguments.

- [x] Implement `scipy.signal.windows.dpss`: symmetric/periodic Slepian
  tapers, one or more orthogonal sequences, `l2` and `approximate`
  normalization, and spectral concentration ratios using an independent
  implicit-shift tridiagonal QL solver. Compare 344 seeded cases against
  SciPy 1.17 (maximum taper error 1.46e-11).
- [x] Add SciFort-specific `multitaper_psd` for 1-D real signals,
  optional detrending, sampling rate and zero-padded FFT length, and
  equal/concentration weighting. Compare 240 seeded spectra with a
  SciPy DPSS + NumPy FFT reference (maximum PSD error 6.89e-13).
  Add checked native tests and an optimized targeted regression.
- [x] Implement SciPy's `norm='subsample'` DPSS peak correction for even
  internal lengths using an O(m) Dirichlet-kernel interpolation identity.
  Validate periodic/symmetric, single/multiple taper, and randomized cases.
- [x] Add a selected-eigenpair DPSS solver using Sturm bisection and
  pivoted inverse iteration with O(n*k) eigenvector workspace, including
  explicit QL/partial selection and automatic use for long windows with a
  few tapers. Validate up to 16,384 points against SciPy 1.17.
- [x] Accelerate DPSS concentration ratios with zero-padded FFT
  autocorrelations, reducing the Toeplitz Rayleigh-quotient work from
  O(k*n*n) to O(k*n*log(n)) with O(n) scratch. Retain the direct method
  and add explicit `ratio_method` selection; check against SciPy 1.17.
- [x] Add complex/two-sided and 2-D multichannel multitaper PSD, spectral
  matrices, cross-spectrum and coherence estimates. Support selected sample
  axes and real/complex mixed CSD; validate 160 randomized spectral matrices
  against independent SciPy DPSS + NumPy FFT estimates.
- [x] Add Thomson-style adaptive multitaper weighting to real/complex,
  channelwise and cross-spectral estimators, with pointwise log-jackknife
  confidence intervals and frequency-dependent effective degrees of freedom.
  Compare 140 randomized cases against independent SciPy DPSS / NumPy FFT /
  Student-t oracles. The original equal/eigen paths remain unchanged.
- [x] Add overlapped real/complex multichannel DPSS multitaper spectrograms
  with cached windows, adaptive weighting, FFT zero padding, and selectable
  sample axis. Add Thomson's harmonic F-test, complex line amplitude and
  pointwise F(2,2K-2) p-values. Verify 190 independent SciPy/NumPy reference
  cases. Both checked and optimized standalone signal suites pass 87/87.
- [x] Implement time-frequency cross-spectra, multitaper coherence, and
  multichannel spectral matrix spectrograms, including mixed real/complex
  signal pairs and independent per-channel adaptive weighting. Check
  120 randomized cases against SciPy DPSS/NumPy FFTs; both 89-test checked
  and optimized signal suites pass. No new dependencies.
- [x] Add incremental real/complex multitaper cross-spectral matrices with
  cached DPSS tapers, arbitrary chunks, selectable rank-2 sample axes, and
  an input-length-independent rolling window plus running spectral matrix.
  Verify chunk-boundary invariance and 80 independent SciPy/NumPy comparisons.
- [x] Extend public `detrend` to real/complex rank 1-3 arrays with one-based
  selected axes and piecewise-linear breakpoints. Verify 720 seeded SciPy 1.17
  comparisons (maximum absolute error approximately 5.35e-15).
- [x] Add `multitaper_harmonic_scan` with Bonferroni, BH and BY multiplicity
  adjustments, selected frequency families, edge-bin exclusion and boolean
  significance decisions. Verify 160 SciPy/NumPy reference cases.
- [x] Add `detrend_mode='linear'` alongside `none` and `constant` to the
  windowed multitaper power, spectral matrix, cross-spectrogram and
  coherence estimators. Validate 64 SciPy 1.17 linear-detrended reference
  cases (36,787 coefficients, maximum difference 8.92e-15).
- [x] Extend `detrend_mode='linear'` to whole-record multitaper PSD, matrix,
  CSD, coherence, jackknife CI, harmonic test/scan, and the stateful
  multichannel `multitaper_stream` estimator. Validate 96 independent
  SciPy/NumPy cases (67,127 complex coefficients) and 8 focused checked
  and optimized multitaper regression programs.
- [ ] Investigate better-calibrated harmonic-line null models under colored
  noise, and more efficient multitaper FFT batching.
- [ ] Add higher-rank or streaming spectral estimation, time-frequency
  multitaper methods and more robust intervals for highly correlated tapers.


- [x] Implement `find_peaks_cwt` with default Ricker wavelets, a reusable
  `cwt_ricker` transform, ridge tracking, percentile noise floors, and
  SciPy-style minimum SNR and ridge length thresholds. Compare 300 seeded
  nonflat cases to SciPy 1.17; note one-sample numerical ties on constant data.
- [x] Add analytic `place_poles` for controllable SISO and fully actuated
  MIMO real state-space systems, with complex-conjugate target pole support.
  Verify 280 seeded controllable systems against SciPy 1.17.
- [x] Extend `place_poles` to partially actuated and rank-deficient-input
  MIMO systems with a real feedback gain, complex-conjugate pairs, and
  feasible repeated roots. Implemented through pivoted nullspaces and
  determinant/independence improving transfer-matrix updates, not YT.
  Verify 470 seeded SciPy-feasible cases, including two published examples.
- [x] Add YT-style rank-two updates for real pairs and complex conjugates,
  KNV0 real rank-one updates, `method`/`rtol`/`maxiter` controls, plus
  iteration/relative-change/transfer-determinant diagnostics. Use the
  Yang--Tits update ordering with SciFort's own nullspace initialization;
  verify 666 randomized method-specific comparisons against SciPy.
- [ ] Improve parity with SciPy's YT/KNV0 singular-value degeneracy handling,
  transfer-matrix initialization and exact diagnostic conventions; add
  generalized Jordan-chain placement for feasible defective assignments.
- [x] Add caller-selectable real-valued wavelets to `find_peaks_cwt` and
  provide a reusable `cwt_wavelet` transform. Validate callback lengths and
  finiteness; verify custom Gaussian, shifted, oscillatory and derivative
  wavelets against SciPy 1.17 for fractional widths and boundary padding.
- [x] Add `cwt_complex_wavelet` with type-safe complex callbacks, plus
  normalized `morlet2_wavelet` and `cwt_morlet2` with configurable center
  frequency for real and complex inputs. Verify 240 SciPy 1.17 complex
  convolution comparisons (28,320 coefficients) in checked and optimized
  builds, including malformed callbacks and one-based public-module tests.
- [x] Add `find_peaks_cwt_complex` for real/complex inputs with Morlet or
  caller-supplied complex wavelets, plus `find_peaks_cwt_magnitude` for
  precomputed complex CWT matrices. Track modulus maxima, use a modulus
  percentile noise floor, and preserve the original signed real peak finder.
  Verify 280 independent magnitude-ridge cases against a SciPy 1.17
  convolution/ridge oracle in checked and optimized builds, retaining exact
  agreement on 300 original Ricker peak cases.

- [x] Add `scipy.signal.deconvolve` (real, complex, and mixed polynomial
  long division with quotient/remainder and explicit invalid-divisor status)
  and `scipy.signal.envelope` (real/complex rank 1-3, one-based time axis,
  FFT-bin bandpass, residual modes, resampling and squared magnitude).
  Verify with 1,900 seeded SciPy 1.17 comparisons and native tests.
- [x] Extend analytic-envelope filtering to real and complex assumed-rank
  `envelope_into` (ranks 1-15), one-based/negative sample axes, real/complex
  residuals, optional resampling, and atomic output after all channels succeed.
  Verify with native rank-15/strided/failing-channel tests and 360 SciPy cases.
- [ ] Consider streaming envelope estimators for long-duration signals.

- [x] Add `order_filter` for masked, zero-padded rank-1/2/3 order statistics
  and `sepfir2d` with mirror-symmetric separable FIR convolution. Provide
  real/integer masked ranks, complex and mixed separable taps, direct public
  exports, deterministic tests, and SciPy 1.17 comparison coverage. Preserve
  the documented sparse-footprint semantics when a SciPy optimized median
  shortcut disagrees with the generic reference.
- [x] Implement independent zero-padded `medfilt`, `medfilt2d`, and `wiener`
  algorithms for rank-1/2/3 arrays, with rectangular windows, real and
  default-integer inputs, and optional Wiener noise power. Check 1,050 seeded
  comparisons against SciPy 1.17, including even-sized Wiener windows.
- [ ] Extend local rank/order and denoising algorithms to higher ranks and
  consider sliding-window acceleration for very large median kernels.

- [x] Implement the inverse STFT for real one-sided and complex two-sided
  spectra, including odd/zero-padded FFT lengths, boundary removal, scaling,
  custom windows and noninvertible overlap diagnostics. Add SciPy-style
  `check_cola` and `check_nola` for named and explicit windows. Verify 480
  randomized inverse-transform cases and 560 window checks against SciPy.
- [ ] Consider newer `ShortTimeFFT` features (e.g. canonical dual windows,
  alternative analysis/synthesis windows and streaming frame subsets).
- [ ] Extend spectral processing to multichannel arrays and selected axes.

## LTI frequency responses and Bode plots

- [x] Add continuous/discrete `freqresp`, `bode`, `dfreqresp`, `dbode`
  with real/complex SISO transfer-function and state-space inputs,
  plus separate ZPK entry points.
- [x] Support SciPy default frequency-bin counts, caller-supplied frequency
  arrays, discrete full-circle grids, phase unwrapping, and digital sample
  interval scaling. Compare 720 seeded cases against SciPy 1.17.
- [x] Add direct resolvent evaluation for high-order state-space systems to
  avoid sensitivity when expanding the characteristic polynomial; compare
  560 seeded systems across 5,600 frequency points against SciPy linear solves.
- [x] Add MIMO state-space analog/digital frequency-response and Bode APIs
  with per-channel response arrays, independent phase unwrapping, and a
  common pivoted factorization for multiple inputs. Compare 480 seeded
  systems (22,648 input-output responses) against SciPy 1.17 linear solves.
- [x] Implement reusable unitary Hessenberg reduction for MIMO frequency
  sweeps, with an adjacent-pivot solver, automatic selection, explicit
  direct/Hessenberg overrides, deterministic regressions, and benchmarks.
- [x] Add one-based input/output channel selection for all four MIMO
  frequency/Bode routines and a public `mimo_frequency_plan` that caches the
  unitary Hessenberg reduction across repeated analog/digital sweeps.
  Validate 280 seeded SciPy reference systems with 18,753 selected values.
- [x] Add an optional dependency-free complex Schur-form solver for
  frequency sweeps and reusable plans, with bounded shifted-QR iteration,
  fallback to Hessenberg, and reference tests.
- [ ] Evaluate optional BLAS/LAPACK backends for large state-space systems
  and compare robust Schur routines with the internal shifted-QR algorithm.

## Discrete-time LTI simulation

- [x] Implement `dlsim`, `dstep`, and `dimpulse` with real/complex state-space,
  SISO transfer-function and ZPK inputs, multi-input response arrays, and
  optional initial states and input-time interpolation.
- [x] Add deterministic Fortran regressions and seeded cross-language
  comparison (498 cases) against SciPy 1.17 and complex state recurrences.
- [x] Add continuous-time `lsim`, `step`, and `impulse`, including exact
  matrix-exponential input interpolation and the analytical impulse response.
  Check 522 seeded SciPy/reference cases and the signal standalone suite.
- [x] Implement SciPy-like default continuous response time grids and extend
  mixed real/complex ZPK gain and input overloads. Verify 184 SciPy 1.17
  response cases and add checked Fortran regression tests.
- [ ] Evaluate state-matrix eigensolvers for stable automatic time-grid
  estimation in high-order and nearly defective systems.
- [ ] Evaluate additional optimized/BLAS-backed state-space simulation paths
  for large dense or sparse state matrices.

- [x] Add optional measured forward/transpose product counts to the
  work-aware refresh cost model for nonsymmetric rank-two discovery,
  including exact-diagonal fallback; validate against manual thresholds
  and nonlinear benchmark cases.
- [x] Extend measured setup-product counting to all automatic Krylov
  diagonal and rank-one/rank-two preconditioners, including progressive
  probes, rank-two sketches, streaming, failed fits, and exact fallback;
  validate conditional rebuilding and numerical differences in checked
  and optimized compiler builds.
- [x] Add an opt-in process-CPU-time cost model for preconditioner rebuilds
  and inner linear solves, with a minimum reliable timing interval, input
  validation, and preserved product-count/default scheduling.
- [x] Compare CPU-time-based, manual-cost, measured-product, and fixed
  refresh scheduling on 32 cases with two synthetic callback workloads
  (forward- and transpose-expensive), including no-timing guard checks.
  See `benchmarks/benchmark_krylov_timed_refresh.f90` and its CSV.
- [x] Repeat the benchmark on actual (rather than delayed) numerical kernels:
  a distributed exponential-response nonlinear model evaluated by midpoint
  quadrature, with analytic Jacobian and adjoint callbacks and OpenMP over
  variables. Compare process CPU with wall time using one and four threads,
  and validate numerical consistency and convergence in checked and optimized
  GNU Fortran builds. See `benchmark_krylov_quadrature_refresh.f90`.
- [ ] Evaluate the CPU-time refresh heuristic on independent application
  workloads (e.g. PDE/optimization solvers) and on multithreaded callbacks
  from libraries other than OpenMP. CPU time is not wall time; do not infer
  general wall-time speedups from synthetic or quadrature-only benchmarks.


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
- [x] SciPy-style vector descriptive helpers `gmean`, `hmean`, `pmean`, `sem`,
  `variation`, `iqr`, `median_abs_deviation`, `moment`, `skew`, and `kurtosis`;
- [x] rank-2 axis reductions and fixed-rank `*_keepdims` companions for those
  descriptive helpers, including trailing broadcast weights for the mean family;
- [x] vector/rank-2 `nan_policy='propagate'` and `nan_policy='omit'` overloads for
  the SciPy-style descriptive helpers, including omission of NaN weights;
- [x] unweighted rank-3 axis reductions and fixed-rank `*_keepdims` companions,
  including `nan_policy='propagate'` / `'omit'`;
- [x] rank-3 `gmean`, `hmean`, and `pmean` trailing broadcast weights, including
  joint singleton expansion and weighted NaN omission;
- [x] unweighted rank-4 axis reductions and fixed-rank `*_keepdims` companions,
  including whole-array IQR and `nan_policy='propagate'` / `'omit'`;
- [x] `gmean`, `hmean`, and `pmean` axis reductions through rank 15 using shared
  assumed-rank engines plus thin fixed-result-rank wrappers; rank-4 through
  rank-15 calls accept scalar and trailing broadcast weights with joint
  singleton expansion;
- [x] `sem`, `variation`, `iqr`, `median_abs_deviation`, `moment`, `skew`, and
  `kurtosis` axis/keepdims support through rank 15 using a shared assumed-rank
  slice engine, including `nan_policy='propagate'` / `'omit'`;
- [x] SciPy-compatible named `scale='normal'` for `iqr` and
  `median_abs_deviation` through rank 15, including keepdims and NaN policies;
- [x] SciPy-style `zscore`, `gzscore`, and `gstd` through rank 15 with
  one-based axes, `ddof`, NaN propagation/omission, and `gstd_keepdims`;
- [x] `zmap` through rank 15, including separate comparison-axis lengths and
  singleton broadcasting on non-reduced dimensions;
- [x] whole-array `axis=None` standardization via `zscore_all`, `gzscore_all`,
  `zmap_all`, and `gstd_all`, plus tuple-axis `gstd_keepdims(x,dims)`;
- [x] `nan_policy='raise'` via Fortran 2018 `error stop` on NaN input;
- [x] SciPy-style trimmed value-limit statistics `tmean`, `tvar`, `tstd`, `tsem`,
  `tmin`, and `tmax` through rank 15, including one-based axes, whole-array
  companions, fixed-rank keepdims, half-open limits, `ddof`, and all three
  NaN policies;
- [x] SciPy-style proportion trimming via `trim_mean`, `trimboth`, and `trim1`
  through rank 15, with one-based axes, explicit `*_all` companions for
  `axis=None`, `trim_mean_keepdims`, and all three `trim_mean` NaN policies;
- [x] SciPy-style `kstat` and `kstatvar` through rank 15, including orders
  1-4 / 1-2, whole-array defaults, explicit one-based axis reductions, fixed-rank
  keepdims companions, and all three NaN policies;
- [x] SciPy-style circular `circmean`, `circvar`, and `circstd` through rank 15,
  including arbitrary angular periods, whole-array defaults, explicit one-based axis
  reductions, fixed-rank keepdims, all three NaN policies, and normalized `circstd`;
- [x] SciPy-style Bayesian `mvsdist` / `bayes_mvs` inference for the mean,
  variance, and standard deviation, including frozen posterior objects, Jeffreys-prior
  Student-t/inverse-gamma/generalized-gamma posteriors through n=1000, the SciPy
  large-sample Gaussian switch, assumed-rank flattening, and central intervals;
- [x] SciPy-style frequency statistics `percentileofscore`, `scoreatpercentile`,
  `cumfreq`, and `relfreq`, including all percentile tie conventions, rank-15 score
  arrays and axis reductions, histogram limits/weights, and SciPy-compatible edge cases;
- [x] SciPy-style multidimensional `rankdata` through rank 15, including all five tie
  methods, whole-array flattening, one-based axis ranking, all three NaN policies, and
  `tiecorrect`;
- [x] SciPy-style `mode` and `describe` through rank 15 via Fortran multi-output
  subroutines, including one-based axis reductions, whole-array companions, `mode`
  keepdims, all three NaN policies, arbitrary `ddof`, and skew/kurtosis bias control;
- [x] SciPy-style `entropy` and `differential_entropy` through rank 15, including
  Shannon/KL normalization, trailing broadcast `qk`, one-based axis reductions,
  whole-array and keepdims companions, all three NaN policies, and the Vasicek,
  Van Es, Ebrahimi, Correa, and automatic differential-entropy estimators;
- [x] SciPy-style `expectile` and `sigmaclip`, including rank-15/flattened
  inputs, broadcast expectile weights, negative-weight compatibility, and iterative
  population-standard-deviation clipping with final thresholds;
- [x] SciPy-style sample `lmoment`, including scalar/vector/default orders,
  rank-15 axis and whole-array reductions, sorted-input and standardized-ratio
  modes, all three NaN policies, and fixed-rank keepdims where representable;
- [x] SciPy-style Box-Cox/Yeo-Johnson power transforms, profile log-likelihoods
  through rank 15, PPCC/MLE parameter estimation, whole-array and keepdims
  companions, fit helpers, and all three NaN policies for likelihood reductions;
- [x] SciPy-style probability-plot diagnostics `probplot`, `ppcc_plot`,
  `ppcc_max`, `boxcox_normplot`, and `yeojohnson_normplot`, including Filliben
  order-statistic medians, least-squares fit summaries, extensible PPF callbacks,
  one-shape PPCC search, and power-transform normality grids;
- [x] empirical-distribution distances `wasserstein_distance`,
  `wasserstein_distance_nd`, and `energy_distance`, including optional sample weights,
  one-dimensional CDF integration, multidimensional Euclidean optimal transport,
  SciPy-compatible input validation, and IEEE-safe nonfinite edge behavior;
- [x] empirical distribution and quantile inference via `ecdf`, `quantile_test`,
  and callback-backed `order_statistic`, including uncensored/right-censored
  Kaplan-Meier estimates, Greenwood/log-log confidence intervals, exact binomial
  quantile tests and confidence bounds, and beta-law order-statistic transforms;
- [x] SciPy-style right-censored survival comparison via public `censored_data`,
  `ecdf(censored_data)`, and `logrank`, including signed normal statistics, all
  three alternatives, mixed censored/uncensored overloads, and trap-safe zero-variance cases;
- [x] t-tests, including `ttest_ind_from_stats`, correlation tests, and
  Mann-Whitney U;
- [x] rank-sum/signed-rank, Kruskal-Wallis, and Friedman tests;
- [x] power-divergence/chi-square/Fisher contingency tests;
- [x] classical/Welch one-way ANOVA, O'Brien variance transformation, and
  Bartlett/Levene/Fligner variance tests;
- [x] one-/two-sample KS and Cramer-von Mises goodness-of-fit tests;
- [x] Shapiro-Wilk, D'Agostino/Pearson, Jarque-Bera, and Anderson-Darling normality/fit tests;
- [x] Generalized Monte Carlo `goodness_of_fit` for normal, exponential, and
  uniform location-scale families with AD, KS, CVM, and Filliben statistics,
  explicit RNG state, and SciPy-style known-versus-prefit parameter semantics;
- [x] Ansari-Bradley, Mood, Baumgartner-Weiss-Schindler, Epps-Singleton, k-sample Anderson-Darling, and median tests;
- [x] Kendall/weighted-Kendall/Spearman-rho/Somers-D/Chatterjee-xi/point-biserial
  association, linear/robust slopes, Brunner-Munzel, and Page trend tests;
- [x] Multiscale graph correlation (`multiscale_graphcorr`) for independence and
  unpaired two-sample testing, including raw/custom/precomputed distances, the full
  local-correlation map/optimal scale, and explicit-state permutation inference;
- [x] First-class continuous-distribution adapters and generic truncation via
  `make_distribution` and `truncate`, including PDF/CDF/tail/inverse-tail methods,
  endpoint-aware support clipping, stable normalization, and explicit-state sampling;
- [x] binomial, Barnard/Boschloo exact, odds-ratio/relative-risk, contingency association, p-value combination, and FDR control;
- [x] Alexander-Govern, Tukey/Tukey-Kramer/Games-Howell, Dunnett, and Poisson means E-tests;
- [x] resampling infrastructure: percentile/basic/BCa bootstrap, exact/randomized permutation tests, Monte Carlo tests, and power simulation.

Breadth should expand only after the numerical foundations and provenance
process are demonstrably reliable.

## Milestone 6: core SciPy-special compatibility

- [x] error functions and inverse error functions;
- [x] real Dawson integral, imaginary/scaled complementary error functions, and Fresnel S/C;
- [x] complex Faddeeva `wofz` foundation plus complex `erf`, `erfc`, `erfcx`,
  `erfi`, and Dawson overloads;
- [x] complex Fresnel S/C overload reusing the complex-error foundation and
  the existing real-axis kernel;
- [x] first-quadrant complex Fresnel zero tables `fresnel_zeros`, `fresnelc_zeros`, and `fresnels_zeros`;
- [x] real-argument modified Fresnel tail pairs `modfresnelp` and `modfresnelm`;
- [x] log-gamma, positive-shape log-beta, and real/complex digamma/psi;
- [x] real/complex normal-CDF `ndtr`/`log_ndtr`, real inverse-normal helpers, and real logistic transforms;
- [x] entropy, relative-entropy, and stable elementary helpers, including real/complex `xlogy`/`xlog1py`;
- [x] Box-Cox transforms and inverses;
- [x] fast real order-zero/one cylindrical Bessel J/Y/I/K functions and scaled K;
- [x] order-zero Bessel integral helpers `itj0y0`/`it2j0y0` and `iti0k0`/`it2i0k0`;
- [x] integer-order Bessel Y/K plus real-order modified-Bessel I/K and scaled I/K;
- [x] real/complex spherical Bessel J/Y and modified I/K functions with first derivatives;
- [x] real Struve/modified-Struve `struve`/`modstruve`, order-zero integral
  helpers `itstruve0`/`it2struve0`/`itmodstruve0`, and Riccati-Bessel
  sequences `riccati_jn`/`riccati_yn`;
- [x] real-order cylindrical Bessel J/Y, real-axis scaled variants, and SciPy-compatible `jn` alias;
- [x] real-axis arbitrary-order derivatives `jvp`, `yvp`, `ivp`, and `kvp`;
- [x] real/complex Hankel `hankel1`/`hankel2`, exponentially scaled
  `hankel1e`/`hankel2e`, and arbitrary-order `h1vp`/`h2vp`, including
  principal-branch continuation across the negative real axis;
- [x] real/complex cylindrical and modified Bessel `jv`/`yv`/`iv`/`kv`,
  scaled `jve`/`yve`/`ive`/`kve`, complex `jn` alias, and arbitrary-order
  `jvp`/`yvp`/`ivp`/`kvp` overloads on the principal branch;
- [x] positive real integer-order Bessel zero tables `jn_zeros`, `jnp_zeros`,
  `yn_zeros`, `ynp_zeros`, combined `jnyn_zeros`, and magnitude-ordered
  `jnjnp_zeros` with order/serial/type metadata, plus real-branch
  `y0_zeros`, `y1_zeros`, and `y1p_zeros` companion-value helpers;
- [x] order-zero Kelvin functions `ber`/`bei`/`ker`/`kei`, their first
  derivatives, combined complex `kelvin` evaluator, and positive zero tables for all
  eight functions through the documented binary64 root-search ceiling;
- [x] integer-order Mathieu characteristic values `mathieu_a`/`mathieu_b`, Fourier
  coefficients `mathieu_even_coef`/`mathieu_odd_coef`, degree-argument angular
  `mathieu_cem`/`mathieu_sem`, and radial `mathieu_modcem1`/`mathieu_modcem2`/
  `mathieu_modsem1`/`mathieu_modsem2`;
- [x] prolate/oblate spheroidal characteristic values `pro_cv`/`obl_cv` and
  sequence helpers `pro_cv_seq`/`obl_cv_seq`, using DLMF tridiagonal eigenproblems;
- [x] first-kind prolate/oblate spheroidal angular functions `pro_ang1`/`obl_ang1`
  and precomputed-characteristic variants `pro_ang1_cv`/`obl_ang1_cv`;
- [x] first-kind prolate/oblate spheroidal radial functions `pro_rad1`/`obl_rad1`
  and precomputed-characteristic variants `pro_rad1_cv`/`obl_rad1_cv`;
- [x] second-kind prolate/oblate spheroidal radial functions `pro_rad2`/`obl_rad2`
  and precomputed-characteristic variants `pro_rad2_cv`/`obl_rad2_cv`;
- [x] real/complex Airy `airy` functions and derivatives plus real/complex
  exponentially scaled `airye`;
- [x] real Airy integrals `itairy` and real Ai/Ai'/Bi/Bi' zero tables via `ai_zeros`/`bi_zeros`;
- [x] real/complex Carlson symmetric elliptic integrals `elliprc`, `elliprd`,
  `elliprf`, `elliprg`, and `elliprj`;
- [x] real Legendre complete/incomplete elliptic integrals `ellipk`, `ellipkm1`,
  `ellipe`, `ellipkinc`, and `ellipeinc`;
- [x] real Jacobi elliptic `ellipj` (`sn`, `cn`, `dn`, and amplitude);
- [x] real/complex exponential integrals `exp1` and `expi`, plus real integer-order
  `expn`, plus real/complex `sici` and `shichi`;
- [x] real/complex `gamma`, `rgamma`, and principal-branch `loggamma`, plus real
  `gammasgn`, `gammaln`, and `multigammaln`;
- [x] real gamma-combinatoric `beta`/`betaln`, generalized `binom`, `poch`,
  and integer-order `polygamma`;
- [x] floating factorial/combinatorial `factorial`, `factorial2`, `factorialk`,
  `comb`, and `perm` (`exact=False`, real zero-extension subset);
- [x] Bernoulli-number and Euler/secant-number sequences `bernoulli`/`euler`, plus
  floating Stirling numbers of the second kind `stirling2`;
- [x] Wright generalized Bessel `wright_bessel` and stable logarithm
  `log_wright_bessel` on SciPy's nonnegative real domain;
- [x] Jahnke-Emden Lambda sequences `lmbda`, weighted Bessel integral `besselpoly`,
  and first-quadrant complex error-function zeros `erf_zeros`;
- [x] real/complex dilogarithm/Spence `spence`, real/complex Riemann and real Hurwitz `zeta`, and
  `zetac` with analytic continuation and stable large-positive complement;
- [x] Add bounded-work native complex Riemann zeta with explicit binary64
  range limits, retained upstream notices, and high-precision validation;
- [x] one-dimensional log-sum-exp, softmax, and log-softmax reductions.
- [x] real rank-2 through rank-15 reduction extensions: whole-array and
  single-dimension `logsumexp`/`softmax`/`log_softmax`, sign-returning
  `logsumexp_signed`, with weighted `logsumexp` through rank 15.
- [x] complex vector through rank-15 `logsumexp`, `softmax`, and `log_softmax`,
  including real/complex weights through rank 15 and complex-unit phases from
  `logsumexp_signed`.
- [x] reduction axis/keepdims parity through rank 15, including weighted
  broadcasting through rank 15: scalar and trailing-rank weights, singleton-dimension
  joint broadcasting that can expand the input shape, complex weights with real
  inputs, tuple-style dimension lists, and fixed-rank keepdims companions. Ranks
  7 through 15 share a flattened group/index engine while retaining fixed-rank
  public wrappers required by Fortran's static result-rank rules; weighted ranks
  6 through 15 use a shared broadcast-map/reduction engine.
- [x] generalized hypergeometric `hyp0f1` for real denominator parameter and real/complex arguments.
- [x] confluent hypergeometric `hyp1f1` for real parameters and real/complex arguments.
- [x] Gauss hypergeometric `hyp2f1` for real parameters and real/complex arguments.
- [x] Tricomi confluent hypergeometric `hyperu` for real parameters and nonnegative real arguments.
- [x] Real parabolic-cylinder `pbdv`/`pbvv`/`pbwa` values and first derivatives,
  plus `pbdv_seq`/`pbvv_seq` real-order sequences and complex integer-order `pbdn_seq`.
- [x] Scalar SciPy-special helpers `agm`, `huber`, `pseudo_huber`, `ndtri_exp`,
  `owens_t`, and `voigt_profile`, with skew-normal sharing the public Owen-T kernel.
- [x] SciPy convenience helpers `cbrt`, `exp10`, `exp2`, real/complex `expm1`/`log1p`,
  `powm1`, `radian`, banker-rounding `round`, normalized `sinc`, `sindg`/`cosdg`/
  `tandg`/`cotdg`, `softplus`, and periodic `diric`.
- [x] Legacy SciPy distribution-special helpers `bdtr`/`bdtrc`/`bdtri`,
  `nbdtr`/`nbdtrc`/`nbdtri`, `pdtr`/`pdtrc`/`pdtri`, `chdtr`/`chdtrc`/`chdtri`,
  `fdtr`/`fdtrc`/`fdtri`, `gdtr`/`gdtrc`/`gdtria`/`gdtrix`, and `stdtr`/`stdtrit`.
- [x] Legacy inverse-parameter helpers `bdtrik`, `bdtrin`, `btdtria`, `btdtrib`,
  `chdtriv`, `fdtridfd`, `gdtrib`, `nbdtrik`, `nbdtrin`, `pdtrik`, `stdtridf`,
  `nrdtrimn`, and `nrdtrisd`.
- [x] Legacy noncentral-distribution special helpers `chndtr`/`chndtrix`/
  `chndtridf`/`chndtrinc`, `ncfdtr`/`ncfdtri`/`ncfdtridfd`/`ncfdtridfn`/
  `ncfdtrinc`, and `nctdtr`/`nctdtrit`/`nctdtridf`/`nctdtrinc`, reusing
  SciFort's noncentral chi-square, F, and Student-t distribution kernels.
- [x] Historical statistical special wrappers `smirnov`/`smirnovi`,
  `kolmogorov`/`kolmogi`, and `tklmbda`, reusing the `ksone`, `kstwobign`,
  and Tukey-lambda distribution kernels.
- [x] Nonnegative-integer-degree orthogonal-polynomial evaluators
  `eval_legendre`, `eval_jacobi`, `eval_gegenbauer`, `eval_laguerre`,
  `eval_genlaguerre`, `eval_hermite`, `eval_hermitenorm`, `eval_chebyt`,
  `eval_chebyu`, `eval_chebyc`, `eval_chebys`, and shifted Legendre/Jacobi/
  Chebyshev variants, with real and complex arguments where SciPy exposes them.
- [x] Classical orthogonal-polynomial coefficient constructors `legendre`, `jacobi`,
  `gegenbauer`, `laguerre`, `genlaguerre`, `hermite`, `hermitenorm`, Chebyshev
  `T/U/C/S`, and shifted Legendre/Jacobi/Chebyshev variants, with SciPy-style
  optional monic normalization and descending-power coefficient arrays.
- [x] Classical Gaussian-quadrature roots and weights `roots_legendre`,
  `roots_jacobi`, `roots_gegenbauer`, `roots_laguerre`, `roots_genlaguerre`,
  `roots_hermite`, `roots_hermitenorm`, `roots_chebyt`, `roots_chebyu`,
  `roots_chebyc`, `roots_chebys`, and shifted Legendre/Jacobi/Chebyshev variants.
- [x] Legacy orthogonal compatibility names `assoc_laguerre`, `p_roots`,
  `j_roots`, `cg_roots`, `l_roots`, `la_roots`, `h_roots`, `he_roots`,
  `t_roots`, `u_roots`, `c_roots`, `s_roots`, `ps_roots`, `js_roots`,
  `ts_roots`, and `us_roots`, delegating to the canonical evaluator/root APIs.
- [x] Complex Lambert-W `lambertw` on arbitrary integer branches and real/complex
  Wright omega `wrightomega`, with DLMF branch-point/asymptotic starts and
  Halley refinement.
- [x] Real-degree/integer-order Ferrers associated Legendre `lpmv` plus
  normalized `sph_legendre_p` and complex `sph_harm_y`, using DLMF
  hypergeometric formulas and stable integer-degree normalized recurrences.
- [x] Integer-degree real/complex `assoc_legendre_p` with branch-cut choices 2/3,
  spherical normalization, and Fortran derivative helpers for associated and
  spherical Legendre functions and spherical harmonics.
- [x] Current SciPy first-kind Legendre table APIs `legendre_p_all`,
  `assoc_legendre_p_all`, `sph_legendre_p_all`, and `sph_harm_y_all`, plus
  scalar real/complex `legendre_p`, with native signed-order Fortran bounds.
- [x] Real ellipsoidal harmonics `ellip_harm`, `ellip_harm_2`, and
  normalization `ellip_normal`, with SciPy-compatible Lame-class indexing
  and sign controls plus endpoint-regularized native quadrature.
- [x] Real/complex Legendre second-kind sequences `lqn` and `lqmn`, including
  derivatives, with Ferrers recurrences inside the unit disk and DLMF
  hypergeometric anchoring outside it.

The reduction API now covers real/complex vectors through rank-15 arrays,
including a single Fortran `dim`, fixed-rank keepdims companions, tuple-style
dimension lists, and signed/phase-returning log-sum-exp. This reaches the standard
Fortran maximum array rank; ranks 7 through 15 use a shared flattened group/index
engine. Weighted `logsumexp` broadcasting also extends through rank 15, including
scalar/trailing-rank weights and singleton-dimension joint broadcasting that can
expand the input shape. Broader complex special functions beyond the
Faddeeva/error/Fresnel/gamma/digamma/
E1/Ei/Si/Ci/Shi/Chi/Spence/Airy/Carlson/hyp0f1/hyp1f1/hyp2f1/hyperu/
normal-CDF/elementary and parabolic-cylinder/cylindrical/modified/spherical-Bessel
foundations remain future work. Python-only `scipy.special` error-state/test helpers are
intentionally outside the numerical Fortran API.

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
- [x] Add directional summary statistics via `directional_stats`, including rank-2
  observation matrices, rank-3 batches, arbitrary one-based observation axes, and
  optional per-observation normalization.
- [x] Add directional and random-matrix distributions: von Mises-Fisher, uniform
  direction, random correlation, orthogonal/special-orthogonal/unitary groups, and
  fixed-margin random tables.



## Milestone 9: density estimation and quasi-Monte Carlo

- [x] Add `binned_statistic`, `binned_statistic_2d`, and `binned_statistic_dd`
  with the seven SciPy built-in statistics, equal-width and explicit edges, custom
  pure callbacks, SciPy-compatible outlier/right-edge bin numbering, and
  arbitrary-dimensional flattened results with explicit shape metadata.
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
- [x] Add optional post-sampling `random-cd` optimization for Halton, Sobol,
  Latin-hypercube, and Poisson-disk draws, with an explicit-state standalone
  optimizer and centered-discrepancy exchange updates.
- [x] Add Lloyd post-sampling optimization for Halton, Sobol, Latin-hypercube,
  and Poisson-disk draws using a dependency-free N-dimensional Bowyer-Watson
  Delaunay/Voronoi layer; document deterministic no-Qhull-joggle behavior for
  degenerate geometries.
- [x] Add `MultinomialQMC` and `MultivariateNormalQMC` using Sobol as the default engine, including singular PSD covariance support and inverse-transform/Box-Muller normal generation.
- [x] Add `sobol_indices` with Saltelli 2010 first/total sensitivity indices, precomputed-evaluation and vectorized-model APIs, independent inverse-CDF transforms, and explicit-state scrambled Sobol sampling.


## Milestone 10: numerical integration

- [x] Add one-dimensional sampled-data `trapezoid`, `cumulative_trapezoid`,
  `simpson`, `cumulative_simpson`, and Romberg `romb`.
- [x] Add fixed-order Gauss-Legendre `fixed_quad` and Newton-Cotes rule
  construction.
- [x] Add scalar adaptive `quad` for finite, semi-infinite, and whole-line
  intervals using the pinned BSD-3-Clause modern QUADPACK 2.1.0 FPM dependency.
- [x] Add `quad` break-point and weighted-integral options using the corresponding
  QUADPACK routines where they map cleanly to a Fortran API.
- [x] Add `dblquad`, `tplquad`, and `nquad`-style nested integration helpers.
- [x] Add `quad_vec`, `qmc_quad`, and `tanhsinh`.
- [x] Add multidimensional `cubature`.
- [x] Add ODE initial-value solvers and a `solve_ivp`-style Fortran interface.
  - [x] Add adaptive real-system `RK23`, `RK45`, and `DOP853`, `t_eval`,
    native seventh-order DOP853 dense output, RK23 Hermite interpolation,
    SciPy-compatible RK45 quartic dense output, and scalar or vector event callbacks
    with per-event terminal/direction controls.
  - [x] Add real-system `BDF` with variable order 1--5, NDF enhancement,
    finite-difference or callback Jacobians, and native dense output.
  - [x] Add real-system `Radau` with fifth-order Radau IIA stepping, cubic
    collocation dense output, and finite-difference or callback Jacobians.
  - [x] Add SciPy-style positive-integer event occurrence limits and per-event
    result grouping/count helpers while preserving the existing flat event storage.
  - [x] Add complex-system `RK23`, `RK45`, and `DOP853`, including `t_eval`,
    retained dense output, and real-valued event functions over complex states.
  - [x] Add complex-system `BDF` with analytic or dense finite-difference
    Jacobians and native dense output.
  - [x] Add SciPy-style scalar-or-vector `atol` handling for all supported real
    and complex methods, including componentwise error scales and numerical
    Jacobian perturbations.
  - [x] Add scalar-or-vector `rtol` handling for all supported real and complex
    methods, including LSODA and implicit BDF/Radau error weights; use a
    conservative scalar summary only for the Newton stopping heuristic.
  - [x] Add SciPy-style constant dense Jacobian matrices for real Radau/BDF and
    complex BDF, plus logical `jac_sparsity` patterns with grouped finite
    differences; LSODA retains its separate `lband`/`uband` interface.
  - [x] Add Fortran-native vectorized RHS callbacks for real Radau/BDF and
    complex BDF finite-difference Jacobians, including sparsity-group batching
    and SciPy-compatible hidden-probe `nfev` accounting.
  - [x] Add native one-based CSC constant sparse Jacobian containers for real
    Radau/BDF and complex BDF, retaining sparse storage in solver state.
  - [x] Add dependency-free sparse LU factorization with partial row pivoting for
    constant CSC Jacobians, including real, complex, and Radau complex-shifted
    Newton systems without allocating the ordinary dense LU workspace.
  - [x] Add sparse-valued Jacobian callbacks and direct CSC numerical-Jacobian
    storage from `jac_sparsity`, with native sparse LU throughout the real
    Radau/BDF and complex BDF sparse paths.
  - [x] Add the `LSODA` method.
    - [x] Add an initial real-system `method='LSODA'` path with automatic
      Adams/BDF switching, dense or callback Jacobian stiffness checks, `t_eval`,
      backward integration, retained dense output, and event support.
    - [x] Add SciPy-style `lband`/`uband` controls for packed analytic Jacobians
      and grouped banded finite differences, plus the LSODA `min_step` retry floor.
    - [x] Extend the independent equal-step Adams path through an order-12
      Adams-Moulton corrector with generated predictor/corrector coefficients.
    - [x] Retain nonuniform Adams derivative/time history across accepted step-size
      changes, adapt Adams order independently of history depth, and use a cheap
      secant stiffness indicator so ordinary nonstiff steps avoid Jacobian probes.
    - [x] Add ODEPACK-style Adams stability limits, a 20-successful-step method-switch
      hysteresis interval, and the order-5 Adams-to-BDF switching guard, using a
      predictor/corrector directional Lipschitz estimate without extra RHS calls.
    - [x] Add the common finite-error ODEPACK ideal-step Adams/BDF switch
      comparison and avoid the redundant endpoint RHS call after accepted BDF steps.
    - [x] Add an Adams Nordsieck predictor state reconstructed from SciFort's
      nonuniform derivative history, with ODEPACK-style power rescaling when the
      trial step changes and Pascal-triangle prediction at the new endpoint.
    - [x] Replace the normal RK45 history startup with low-order Nordsieck Adams
      stepping using ODEPACK `DCFODE` correction/error constants, establishing the
      order-1 startup later extended into the persistent all-order Adams path.
    - [x] Replace the fixed low-order startup corrections with an ODEPACK-style
      functional-iteration convergence test, including roundoff convergence, contraction
      rate estimation, and a factor-of-four retry when the corrector loses contraction.
    - [x] Add the ODEPACK roundoff-pollution Adams-to-BDF switch branch and a
      conservative repeated-error-test recovery path; rejected Adams trials no longer
      probe Jacobians or switch methods.
    - [x] Carry the accepted Adams Nordsieck state persistently through orders 1--12,
      apply the ODEPACK `EL` correction directly, use native Nordsieck dense output,
      and adapt successful-step `IALTH`/`TESCO` order/step reconsideration.
    - [x] Match the ODEPACK rejected-step lower/current-order retry selection and
      repeated-error-test step cap while retaining SciFort's conservative restart.
    - [x] Transfer the current BDF/NDF multistep polynomial into Adams Nordsieck
      state on BDF-to-Adams switches and add ODEPACK-style stiff Jacobian refresh
      checks for 20-step age and substantial Newton-coefficient changes.
    - [x] Match DLSODA's automatic initial-step formula, accepted-state error
      weights, weighted maximum norm, and Adams `RMAX` startup/failure limits;
      reproduce SciPy 1.17.0 nonstiff work counts on scalar and vector regressions.
    - [x] Preserve the accepted Adams polynomial, switch order, and ideal BDF step
      ratio on Adams-to-BDF transitions; defer Jacobian work until the first stiff
      trial; and use ODEPACK-style pure-BDF coefficients, accepted-state maximum
      norms, `DPRJA` finite differences, and finite-difference `nfev` accounting
      only inside LSODA mode.
    - [x] Replace the remaining SciPy-derived stiff step/order and Newton convergence
      controller with a native BDF Nordsieck/chord path using `DCFODE`/`DSTODA`
      constants and `IALTH` selection; reproduce SciPy 1.17.0 work counts and key
      accepted-time/order transitions on a representative stiff regression.
- [x] Add a legacy `odeint`-style LSODA facade.
  - [x] Add time-grid output with repeated monotone times, legacy `func(y,t)`
    callbacks, a type-safe `odeint_tfirst` t-first companion, scalar/vector
    `rtol` and `atol`, SciPy-like default tolerances, dense or packed-banded
    `Dfun` with `col_deriv`, `h0`/`hmax`/`hmin`, and the legacy `mxstep` control.
  - [x] Add `tcrit` hard-stop semantics and LSODA order ceilings through
    `mxordn`/`mxords`, preserving solver history rather than restarting at each
    critical point.
  - [x] Add exact per-output-call `mxstep` semantics with the legacy default of 500,
    preserving LSODA history while resetting the accepted-step allowance at each requested output.
  - [x] Add `full_output`-style per-output LSODA diagnostic histories (`hu`,
    `tcur`, `tolsf`, `tsw`, `nst`, `nfe`, `nje`, `nqu`, `mused`) plus
    `imxer`, `lenrw`, and `leniw`.
  - [x] Add `ixpr`, `mxhnil`, and `printmessg` compatibility controls. `ixpr`
    reports Adams/BDF switches, `mxhnil=0` selects the DLSODA default message cap
    of 10, and `printmessg` reports the final integration status. SciFort treats a
    too-small LSODA step as terminal, so at most one `mxhnil` warning is currently
    observable per integration rather than continuing after `T+H=T` as DLSODA can.

## Milestone 11: numerical differentiation

- [x] Add `scifort_differentiate` with SciPy-style adaptive scalar `derivative`,
  including nested central/one-sided finite-difference stencils, configurable
  order/step controls, error estimates, status/success flags, and `nit`/`nfev`.
- [x] Add scalar-output gradients and vector-output `jacobian` for real vector
  arguments, including per-coordinate one-sided direction controls.
- [x] Add `hessian` by nested adaptive differentiation with SciPy-compatible
  tolerance tightening and original-function evaluation accounting.
- [x] Add SciFort extension complex-step scalar derivative, gradient, and Jacobian
  helpers for analytic functions, adapting the scalar formula from ForDiff 0.1.1.
- [x] Add SciPy 1.17 reference generation and Rosenbrock/Jacobian/one-sided
  regression coverage under runtime checks and IEEE traps.
- [ ] Consider NumDiff-style sparse Jacobian coloring/partitioning when sparse
  optimization APIs require it.

## Milestone 12: optimization

- [x] Add `scipy.optimize.linear_sum_assignment` for rectangular minimum and
  maximum bipartite matching, with forbidden edges, infeasibility statuses,
  O(n+m) workspace, and SciPy 1.17 reference objectives.
- [x] Add `scipy.optimize.isotonic_regression` via weighted O(n) pool-adjacent-
  violators, increasing/decreasing fits, pooled weights, one-based block
  boundaries, and SciPy 1.17 regression values.

- [x] Add `scifort_optimize` scalar root finding with SciPy-style `root_scalar`
  dispatch and `root_result` convergence/evaluation accounting.
- [x] Add `bisect`, `brentq`, `brenth`, and `ridder`, adapting the numerical
  kernels from roots-fortran 1.0.0 under BSD-3-Clause attribution.
- [x] Add `toms748` plus secant, derivative Newton, numerical-derivative Newton,
  and Halley open methods with SciPy 1.17-compatible defaults and method
  selection.
- [x] Add SciPy 1.17 reference generation and strict regression coverage for
  roots, iteration counts, function-call counts, invalid brackets, nonfinite
  objectives, and nonconvergence.
- [x] Add scalar minimization (`minimize_scalar`, Brent, bounded, golden), with
  SciPy 1.17 reference generation and convergence/evaluation-count regressions.
- [x] Export standalone `scipy.optimize.bracket`-style downhill bracketing,
  with configurable growth and iteration limits, explicit Fortran results,
  and SciPy 1.17 point/value/evaluation-count comparisons.
- [x] Add scalar `scipy.optimize.elementwise.bracket_root`-style sign-change
  expansion with finite/infinite limits, configurable growth and iteration cap,
  status/count reporting, and SciPy 1.17 numerical regressions.
- [x] Add scalar `scipy.optimize.elementwise.bracket_minimum`-style minimum
  bracketing, including bidirectional downhill expansion, finite limits,
  counts/status, and SciPy 1.17 reference comparisons.
- [x] Add scalar `scipy.optimize.elementwise.find_minimum` with Chandrupatla
  quadratic/golden-section bracketing, configurable tolerances, callback,
  work/status reporting and SciPy 1.17 regressions.
- [x] Add scalar `scipy.optimize.elementwise.find_root` with Chandrupatla
  interpolation/bisection, configurable tolerances, iteration callbacks,
  work/status accounting, and SciPy 1.17 reference tests.
- [x] Add rank-one independent-problem overloads for the four elementwise
  bracket/root/minimum helpers, with scalar/rank-one bound and tolerance
  broadcasting, per-element results, and SciPy 1.17 regression values.
- [x] Add rank-two grids for all four elementwise helpers, with rank-zero,
  trailing-axis rank-one, and singleton-axis rank-two parameter broadcasting.
  Match SciPy 1.17 reference values and preserve independent counts/statuses.
- [x] Extend all four elementwise helpers to rank-three independent-problem
  grids with scalar/rank-one/rank-two/rank-three trailing broadcasting,
  SciPy 1.17 reference comparisons, and trap-safe regression tests.
- [x] Extend elementwise problem grids to Fortran's maximum rank (15 for
  bracketing and 14 for root/minimum refinement because bracket inputs add
  an endpoint axis), sharing a flattened broadcasting/dispatch engine.
- [x] Add genuinely batched vector callback evaluation for rank-one
  `find_root_vectorized`, with active-lane packing and per-element convergence
  and work accounting; compare status/counts with SciPy 1.17.
- [x] Add `find_minimum_vectorized` for rank-one batched Chandrupatla
  refinement, matching SciPy 1.17 status/work counts and scalar per-lane
  results while packing only active trial points into vector callbacks.
- [x] Add rank-one `bracket_root_vectorized` and `bracket_minimum_vectorized`
  with packed initial endpoints/triples, only active expansion probes in
  later callbacks, independent status and work counts, and SciPy 1.17
  comparisons.
- [x] Extend all four packed vector-callback elementwise helpers to rank-two
  problem grids, preserving active-lane packing, shape, independent status and
  work accounting, and NumPy-style trailing-axis parameter broadcasting.
- [x] Extend packed vector-callback evaluation to rank-three and rank-four
  problem grids for all four routines, with trailing-axis broadcasting and
  parity tests against flattened rank-one callbacks.
- [x] Extend packed vector callbacks to rank 5-14 for refinement and 5-15
  for bracketing via four assumed-rank `*_vectorized_into` subroutines,
  reusing the rank-one controllers and broadcasting while retaining the
  rank-one through rank-four function overloads. High-rank/empty/invalid/
  noncontiguous-array regression added.
- [x] Add explicit trailing-axis broadcasting of a real-valued objective-
  parameter stack to all four packed-callback controllers via assumed-rank
  `*_vectorized_args_into` entry points. Repack args by active-lane ID during
  refinement and expansion, preserving invalid-input/shape isolation and work
  accounting; verify scalar/rank-five/rank-14/15 calls. The original
  `find_root`, `bracket_root`, `find_minimum`, and `bracket_minimum` generics
  still invoke scalar callbacks per element.
- [x] Support independently shaped real-valued objective argument lists in
  all four packed solvers via Fortran 2018 assumed-rank constructors and
  overloads of `*_vectorized_args_into`; retain the numeric stack API,
  original callback signature, active-lane indexing and work accounting.
  Cover rank 0/14/15, scalar/multiple-axis broadcasting, invalid lists,
  noncontiguous outputs and stack/list numerical parity.
- [x] Add independently broadcast mixed-type objective arguments with
  original Fortran 2018 typed scalar boxes, preserving integer, logical,
  real, character, complex and user-defined dynamic types through rank 15.
  Expose typed callback overloads for all four packed root/minimum solvers;
  reuse active-lane controllers, retain numeric APIs, and test compaction,
  high ranks, bad lists, type preservation and noncontiguous grids.
- [x] Offer an allocation-avoiding indexed callback for all four typed
  elementwise packed solvers. The original heterogeneous payload boxes are
  passed without per-evaluation copies; a compact integer index map selects
  source values for each active abscissa. Retain the existing copied typed
  callback API and test exact result/work-count parity, mixed-type broadcasting,
  invalid lanes and ranks 14/15.
- [x] Benchmark indexed and copied callbacks on derived-type payloads for
  all four elementwise controllers. Add `benchmarks/benchmark_elementwise_typed_args.f90`
  with command-line problem/payload/repetition sizes, warmups, counter/result
  parity checks, and alternated timing order. Reuse a preallocated integer
  gather buffer across active evaluations instead of a vector-subscripted
  temporary; experiments show payload-copy avoidance dominates the speedup,
  with index-gather differences small on GNU Fortran 14.2.

- [x] Add the core multidimensional nonlinear-equation family (`root`, `fsolve`)
  using modern MINPACK plus attributed SciPy/modern-Fortran references.
  - [x] Add MINPACK `root(method='hybr')`, `root(method='lm')`, and `fsolve`,
    including analytic Jacobians, banded finite differences, scaling, and
    overdetermined LM systems.
  - [x] Add Broyden/Anderson/linear-mixing/Krylov-style `root` methods.
    - [x] Add dense `broyden1`, `broyden2`, `anderson`, `diagbroyden`,
      `linearmixing`, and `excitingmixing`, including SciPy-default Armijo/no-line-
      search behavior and SciPy 1.17 convergence/work-count regressions.
    - [x] Add `krylov`/`newton_krylov` with matrix-free LGMRES and GMRES,
      SciPy-compatible forcing terms, finite-difference Jacobian products, and
      generated SciPy 1.17 work-count regressions.
    - [x] Add strong-Wolfe line search plus Broyden `max_rank` controls with
      restart, simple, and SVD rank reduction.
  - [x] Add `df-sane` spectral residual solving with Cruz and Cheng nonmonotone
    line searches and SciPy 1.17 work-count/termination regressions.
  - [x] Add SciPy-style `tol_norm=` convergence callbacks to all six
    inexact-Newton methods and Krylov, including direct convenience wrappers,
    residual/step/relative tolerances, and invalid-norm handling.
  - [x] Add advanced Krylov preconditioner/custom-inner-solver hooks: optional
    matrix-free inverse-Jacobian actions, accepted-step adaptive updates, and
    user-supplied inner solvers with counted Jacobian-vector evaluations.
    Retain GMRES/LGMRES defaults and SciPy 1.17 reference regressions.
  - [x] Add built-in Krylov inner methods BiCGSTAB, CGS, MINRES, and TFQMR,
    including matrix-free preconditioning, SciPy 1.17 linear/nonlinear parity
    references, and breakdown/operator-failure regressions. MINRES requires a
    symmetric Jacobian and positive-definite preconditioner.
  - [x] Add opt-in original-residual verification for TFQMR inner solves,
    retaining the SciPy-compatible estimate-based default. Reject optimistic
    convergence and incomplete inner iterates when verification is requested;
    cover failed verification callbacks and deterministic 64-variable systems.
  - [x] Add opt-in mathematically consistent right preconditioning to
    built-in TFQMR (A*M recurrence with M-applied updates), retaining SciPy
    1.17's M*A recurrence as the default. Check 16- and 64-variable
    residuals, original-residual verification, public API validation, and
    independent SciPy reference solves of the transformed system.
  - [x] Add `jacobian_diff_scheme='3-point'` to Krylov root finding, with
    centered matrix-free Jacobian-vector products, a cube-root-epsilon
    automatic step, exact residual-evaluation accounting, and guards for
    nonfinite samples and unrepresentable perturbations. Retain the default
    SciPy-compatible two-point stencil and work counts; test every public
    wrapper and all six built-in Krylov inner methods.
  - [x] Add opt-in `jacobian_step_mode='adaptive'` and `'componentwise'` to
    matrix-free Krylov Jacobian actions. Use relative steps tailored to each
    active coordinate when a global directional step is ill-scaled or rounds
    away. Preserve SciPy-compatible global behavior and nfev counts by default;
    cover different variable/residual scales and all six inner solvers.
  - [x] Add a Fortran analytic `jacobian_product=` callback for Krylov
    `root` and `newton_krylov`, bypassing finite-difference cancellation and
    extra residual evaluations when a direct Jacobian-vector product is known.
    Support all six inner methods, user solvers, true-residual verification and
    preconditioning. Test zero directions, nonlinear refresh, and failures.
  - [x] Add optional `inner_check_residual=.true.` to certify the *actual*
    normalized matrix-free linear residual before accepting any Krylov step,
    including built-in and user inner solvers. Use overflow-safe scaled
    comparisons, preserve default incomplete-iterate acceptance and SciPy
    reference counts, and cover callback failures, all six methods, and
    residual-evaluation accounting. Harden Krylov forcing-term norm ratios
    and detect tiny but nonzero steps without underflow near machine limits.
  - [x] Add opt-in GMRES/LGMRES modified Gram-Schmidt reorthogonalization
    (`inner_orthogonalization='selective'` via DGKS, or `'twice'`), with
    consistent Hessenberg correction coefficients and near-breakdown checks.
    Preserve the SciPy-compatible one-pass default. Verify 48-variable
    ill-conditioned Arnoldi orthogonality, unchanged J*v work counts,
    preconditioned LGMRES and public wrapper option validation.
  - [x] Replace the rank-deficient LGMRES/GMRES small least-squares fallback
    based on normal equations with a directly rank-revealing one-sided Jacobi
    SVD. Preserve the fast full-rank triangular backsolve, caller APIs and
    work counts; test dense rotated nearly singular matrices, large/small
    scales, minimum-norm solutions, and public nonlinear root solvers.
  - [x] Add optional `outer_prune_tol` to prune retained LGMRES augmentation
    directions that are numerically dependent on the existing solution-space
    span, including dependencies on linear combinations (not just pairwise
    duplicates). Use two-pass orthogonalization without extra J*v calls,
    preserve FIFO eviction for new directions, and keep the default retention
    behavior and SciPy evaluation counts unchanged. Cover public wrappers,
    strong scaling, nonfinite inputs, and invalid combinations.
  - [x] Add optional `outer_image_prune_tol` for LGMRES to recheck retained
    directions by their *current* Jacobian images, discard nearly dependent
    history before Arnoldi, and test the new direction on insertion. Keep
    original solution-space directions, explicit extra J*v evaluation
    accounting, SciPy-compatible defaults, and independent finite-difference
    and analytic-product regressions.
  - [x] Compare no pruning, solution-space pruning, and current-Jacobian
    image-space pruning on deterministic dense nonlinear systems with
    changing near-singular Jacobians, recording success and actual residuals,
    analytic J*v work and finite-difference `nfev`, with separate reference
    output for the optimized GNU Fortran build. Guard evaluation accounting
    and convergence on the easier cases, not nonportable performance rankings.
  - [x] Add optional `outer_image_prune_auto` for LGMRES, with a per-image
    dependence cutoff derived from observed first-pass projection loss and
    a dimension-scaled roundoff floor. Preserve explicit fixed thresholds,
    verify preflight option handling, and compare the automatic policy on the
    deterministic changing-Jacobian benchmark without asserting that it
    always wins on work or convergence.
  - [x] Share the existing retained-image buffer with the incremental
    QR basis when `outer_image_basis_reuse` is enabled without Arnoldi
    retained-image caching. Avoid an additional n-by-history allocation,
    lift the 8 MiB cap for that no-extra-storage path, and leave the
    separate-basis cached-J*v path unchanged. Test dependent image spans,
    one-shot invalidation, full FIFO history, complete nonlinear convergence
    and evaluation counts, and paired end-to-end RSS/CPU measurements.
  - [x] Reuse the current-Jacobian image-pruning basis for the *one*
    proposed Newton step following Arnoldi, rather than reconstructing the
    retained-image QR span for candidate insertion. Offer an opt-in
    `outer_image_basis_reuse` keyword for fixed/automatic image pruning;
    rebuild it at each nonlinear iterate, invalidate after insertion, and
    cap extra basis storage at 8 MiB. Verify direct decisions, FIFO
    retention, unchanged residual/J*v counts, changing-Jacobian stress cases,
    and a standalone filtering-stage benchmark.
  - [x] Add opt-in `outer_image_reuse` to avoid reevaluating filtered
    current-iterate Jacobian images inside LGMRES Arnoldi. Keep exactly
    corresponding retained directions, invalidate the cache at each new
    nonlinear iterate, and preserve the previous default and work counts.
    Cover analytic products, preconditioners, fixed and automatic pruning,
    two- and three-point residual accounting, and the five-policy stress test.
  - [x] Bound raw Arnoldi product storage for LGMRES candidate-image
    reconstruction by an explicit `outer_candidate_cache_max_mib` budget.
    If the projected raw cache exceeds the cap, recompute the candidate
    with one exact analytic Jacobian product instead of approximating from
    Hessenberg coefficients. Preserve unlimited-cache defaults, release
    prior raw storage when the history crosses its budget, and test both
    public root entry points and resource usage on large matrix-free systems.
  - [x] Add opt-in `outer_candidate_image_reuse` for analytic LGMRES J*v
    callbacks with image pruning. Reconstruct the new Newton-step image from
    pre-orthogonalization Arnoldi products, separately from the retained
    direction cache. Reject numerical-Jacobian use because finite-difference
    directional products are not exactly linear, preserve the default, and
    compare all combinations on 16/32/64/128-variable nonlinear problems.
  - [x] Reduce the memory required by image-space LGMRES pruning and
    combined retained/candidate caching: compact retained vectors in place
    and share already validated retained J*v images with Arnoldi's candidate
    reconstruction, storing raw products only for new Krylov columns. Keep
    exact pre-orthogonalization products (rather than numerically sensitive
    Hessenberg reconstructions) and verify the ill-conditioned regressions.
    Add a large, O(n)-storage process-RSS benchmark.
  - [x] Reduce basis allocation churn in GMRES and LGMRES: reuse one
    per-solve Arnoldi workspace across nonlinear iterations, use the next
    basis column as the matvec scratch vector, and avoid duplicating retained
    augmentation directions in the stored preconditioned Krylov columns.
    Preserve cached raw Jacobian products and exact default evaluation counts;
    compare peak RSS at 250,000 and 500,000 variables.
  - [x] Avoid reserving the entire `outer_k` Arnoldi basis or image matrix when
    only a few augmentation vectors are active. Geometrically grow the
    persistent basis, initially size the per-iterate image array according to
    active directions, and write preconditioner outputs into their persistent
    Krylov columns. Preserve old solver defaults and evaluation counts; test
    preconditioned GMRES and LGMRES with large unused histories and compare
    process RSS against the previous workspace snapshot.
  - [x] Remove remaining full-length candidate and spare-column allocations
    from retained-image caching: borrow the dormant preconditioned Krylov
    workspace during screening, store exactly one image per live history
    direction (including a zero-column initial history), and skip zero-filling
    columns written before use. Check caching combinations, zero/one/full FIFO
    history, strict runtime-checking builds, and matched peak-RSS comparisons.
  - [x] Reduce repeated temporary allocations in LGMRES augmentation-space
    orthogonalization: share a solver-local pruning workspace, reuse current
    Jacobian image/candidate buffers across iterations, and orthogonalize in
    retired Arnoldi basis columns instead of reserving a second large basis.
    Preserve the allocation-per-call low-level API, match filtering decisions
    and evaluation counts, and check the memory tradeoff with large matrices.

  - [x] Reduce repeated Jacobian-image span work during LGMRES history
    filtering: build one two-pass incremental orthonormal image basis per
    nonlinear iterate instead of reconstructing it for each retained vector.
    Reset it as the Jacobian changes, preserve the previous low-level
    allocation-per-call API and solver options, and reuse Arnoldi scratch or
    active-history-sized workspace storage. Validate fixed and automatic
    thresholds, near-dependence, changing-Jacobian work counts, and a
    reproducible long-history timing comparison.

  - [x] Avoid repeated solution-space span reconstruction for LGMRES while
    the augmentation history is filling: maintain a resettable two-pass
    orthonormal span, preserve exact retained solution vectors and threshold
    semantics, and fall back to rank-revealing reconstruction for FIFO
    evictions. Limit the optional additional per-solve basis to 8 MiB; verify
    near-dependent/subnormal inputs and nonlinear work-count compatibility,
    and benchmark append-only candidate insertion.

  - [x] Support opt-in full-history FIFO QR downdates for solution-space
    LGMRES pruning: maintain a small upper-triangular factor, rotate the
    orthonormal basis after evicting the oldest direction, reorthogonalize
    the proposed direction, and revert to reconstruction when the retained
    history loses numerical rank. Rebuild the factor periodically and limit
    the additional basis to the existing 8 MiB workspace threshold. Compare
    full-history decisions and stored vectors with the original routine,
    exercise public root wrappers, and benchmark repeated evictions.

- [ ] Add unconstrained and constrained multivariate minimization methods.
  - [x] Add `minimize` plus unconstrained `Nelder-Mead`, `Powell`, `BFGS`, `CG`,
    and `Newton-CG`, including analytic/forward-difference gradients,
    dense/Hessian-product Newton-CG paths, SciPy-style strong-Wolfe fallback,
    method diagnostics, and generated SciPy 1.17 regressions.
  - [x] Add L-BFGS-B with SciPy-style bounds, finite-difference gradients,
    stopping controls, and a pinned modern-Fortran FPM backend.
  - [x] Add TNC and remaining bounded-method parity.
    - [x] Add SciPy-style bounds to Nelder-Mead and Powell, including simplex
      clipping/reflection and finite/one-sided bounded Powell line searches.
    - [x] Add native TNC with scaled active bounds, projected-gradient stopping,
      BFGS-preconditioned truncated-CG Hessian products, and bounds-aware numerical gradients.
  - [x] Add SLSQP with bounds, equality/inequality callbacks, analytic or
    bounds-aware finite-difference Jacobians, and a pinned modern-Fortran backend.
  - [x] Add `dogleg` with analytic positive-definite Hessians, Cauchy/Newton
    dog-leg subproblems, SciPy-style trust-radius updates, Cholesky failure
    handling, and generated SciPy 1.17 work-count regressions.
  - [x] Add `trust-ncg` with dense-Hessian or Hessian-product Steihaug CG,
    SciPy-style trust-radius updates, negative-curvature boundary handling, and
    generated SciPy 1.17 work-count regressions.
  - [x] Add `trust-exact` nearly-exact trust-region minimization.
  - [x] Add `trust-krylov` GLTR/Lanczos trust-region minimization.
  - [x] Add COBYLA with a PRIMA-derived trust-region/simplex kernel, bounds,
    equality-to-paired-inequality conversion, `catol`/`f_target` stopping, and
    generated SciPy 1.17 regressions.
  - [ ] Add `trust-constr`.
    - [x] Add dense equality-constrained Byrd-Omojokun SQP with normal/tangential
      trust-region steps, pseudoinverse projections, analytic/numerical Jacobians,
      BFGS fallback, and SciPy 1.17 regression targets.
    - [x] Add a dense positive-slack log-barrier interior-point layer for nonlinear
      inequalities, one/two-sided bounds, fixed bounds, mixed equality/inequality
      systems, and analytic/numerical derivatives.
    - [x] Add SciPy-style `barrier_tol`, `initial_barrier_parameter`, and
      `initial_barrier_tolerance` controls, including `tol` propagation.
    - [x] Add SciPy-style `initial_constr_penalty` merit-function control and
      final `constr_penalty` reporting for equality and interior-point paths.
    - [x] Accept SciPy's `initial_tr_radius` spelling while retaining the
      earlier SciFort `initial_trust_radius` spelling for compatibility.
    - [x] Add analytic equality-constraint Hessians with SciPy's weighted
      `hess(x, v)` convention and Lagrangian-Hessian assembly.
    - [x] Extend analytic constraint Hessians through the inequality/slack path,
      mapping barrier multipliers back to the original nonlinear constraints.
    - [x] Add SciPy-style second-order correction for rejected pure equality-SQP
      steps, with a classic Maratos-effect regression.
    - [x] Extend SOC to nonlinear-inequality barrier subproblems, clipping the
      correction against SciPy's `0.995` fractional slack boundary.
    - [x] Add SciPy's scaled-slack trust-box geometry to barrier normal/tangential
      steps, carry the trust radius between barrier stages, preserve inner SQP
      multipliers for outer optimality checks, and enable bound-only SOC.
    - [x] Reuse cached objective/gradient/objective-Hessian state across barrier
      stages and count each stage's initial stopping state, giving exact SciPy 1.17
      iteration/objective/gradient/Hessian counts on the active-inequality Maratos
      regression.
    - [x] Report basic SciPy-style `constr_violation`, `optimality`,
      `barrier_parameter`, and `barrier_tolerance` result diagnostics.
    - [x] Add SciPy-style Lagrangian-gradient reporting plus cumulative
      projected-CG iteration counts and final CG termination codes.
    - [x] Match SciPy 1.17 status `4` when an apparent `gtol`/`xtol` termination
      still leaves constraint violation above `gtol`.
    - [x] Add final constraint values/Jacobians, SciPy internal method names,
      and per-block constraint function/Jacobian/Hessian evaluation counts; the
      current single combined user-constraint callback is reported first and a
      bounds block is appended second, matching SciPy's bounds ordering.
    - [x] Add SciPy-style per-component lower/upper bounds for the current combined
      nonlinear-constraint callback, including equality, lower-only, upper-only,
      and two-sided interval canonicalization with raw-value/Jacobian diagnostics.
    - [x] Add per-component `keep_feasible` masks for the current bound-vector
      nonlinear-constraint callback and variable bounds, including infeasible-start
      validation and barrier trial-point slack enforcement.
    - [x] Clip numerical objective-gradient and nonlinear-constraint Jacobian
      probes to retained variable bounds, using feasible one-sided/shortened steps.
    - [x] Permit more equality rows than variables when the system is compatible,
      using minimum-norm Gram-pseudoinverse projections for redundant linear and
      nonlinear equalities, including variable bounds and fixed bounds; detect
      inconsistent systems with status 4. Tested against SciPy 1.17's
      `factorization_method='SVDFactorization'` option.
    - [x] Add first-class multiple SciPy-style constraint objects, with linear
      matrices, nonlinear callbacks, separate bounds and derivative controls.
    - [x] Add greedy column grouping for `%finite_diff_jac_sparsity` masks on
      independent nonlinear constraint objects, using bound-aware two-point
      and three-point probes and retaining dense SQP factorization.
    - [x] Match SciPy's signed and zero-valued `%finite_diff_rel_step` semantics
      on dense and grouped constraint Jacobians, including `sign(x)*abs(x)`
      scaling, automatic fallback, and feasible one-sided bound adjustment.
    - [x] Preserve signed two-point objective-gradient relative steps at
      retained variable bounds during both initial barrier differentiation
      and subsequent equality-SQP gradient updates; verify probe directions
      and gradient values against SciPy 1.17.
    - [x] Add numerical objective Hessians from analytic gradients with
      `objective_hess_fd_scheme='2-point'` or `'3-point'`, explicit dense
      symmetrization, bounded/one-sided probes, fixed-variable handling,
      gradient work accounting, and SciPy 1.17 comparisons.
    - [x] Add opt-in `'2-point'` / `'3-point'` weighted numerical Hessians
      for nonlinear constraint blocks with analytic Jacobians, honoring
      relative steps and retained bounds, with per-block evaluation accounting
      and SciPy 1.17 regression references.
    - [x] Add opt-in per-block constraint `%hess_update='bfgs'` and `'sr1'`
      HessianUpdateStrategy counterparts, keeping independent scalar-output
      Hessian approximations so changing multipliers do not corrupt secant
      updates, mixed exact constraint Hessians, and SciPy 1.17 tests.
    - [x] Add `trust-constr` iteration callbacks with optional early stop
      (`status=3`), including the initial state, cumulative barrier-stage
      iteration numbering, physical-variable-only callbacks, and forwarding
      through `minimize`, constraint objects, and `hessp`.
- [x] Add core nonlinear least squares and curve fitting.
  - [x] Add `least_squares` with dense TRF and dogbox, finite bounds, analytic,
    bounds-aware two-point/three-point, or complex-step Jacobians, fixed or Jacobian-based
    variable scaling, evaluation limits, SciPy built-in robust losses, and custom
    callable robust losses; retain
    MINPACK-backed `method='lm'.`
  - [x] Add legacy `leastsq` and `curve_fit`, including weighted covariance,
    bounded TRF/dogbox dispatch, analytic or two-point/three-point/complex-step
    numerical Jacobians, `diff_step`, and robust losses.
  - [x] Add TRF/dogbox per-iteration callbacks with SciPy-style early-stop status.
  - [x] Add TRF/dogbox `tr_solver='lsmr'`; TRF uses SciPy-style regularized
    two-dimensional trust-region subspace solving and dogbox uses an iterative
    scaled Gauss-Newton step.
  - [x] Add analytic matrix-free/LinearOperator-style Jacobians through paired
    `jac_matvec=` and `jac_rmatvec=` callbacks without storing a dense Jacobian.
  - [x] Add `jac_sparsity` grouped finite differences for TRF/dogbox, with
    automatic LSMR selection for both trust-region methods.
  - [ ] Add sparse stored Jacobians and worker parallelism.
- [ ] Add global optimization and linear/mixed-integer programming families.

## Krylov memory management (October 2026)

- [x] Add an opt-in restarted GMRES subspace bound (`inner_restart`) to avoid
  allocating the complete `inner_maxiter` Arnoldi basis for large systems.
- [x] Preserve the SciPy-compatible non-restarted GMRES default and add
  preconditioning, finite-difference, work-count and true-residual regressions.
- [x] Add opt-in restarted LGMRES using `inner_restart`, retaining the
  augmentation history across Arnoldi cycles and updating it only after the
  complete accumulated Newton correction. Verify true residuals per cycle,
  check preconditioning and image-pruning/caching interaction, and benchmark
  full versus bounded-basis memory for large matrix-free systems.
- [x] Add optional adaptive LGMRES restart growth after weak true-residual
  reduction, with explicit maximum restart width, unchanged fixed-width
  defaults, and regressions for zero-progress and ill-scaled systems.
- [x] Benchmark fixed-width, adaptive-width and full-width LGMRES
  on large, genuinely ill-conditioned matrix-free nonlinear systems,
  reporting true and scaled residuals, solution errors, Jacobian-vector
  callback counts, CPU time and independently measured peak process RSS.
  Add scaled-stopping-norm regression and document premature absolute-residual
  convergence on difficult systems.
- [x] Add a direct `inner_m_diag=` array interface for a static diagonal
  inverse-Jacobian preconditioner in all six Krylov methods, with validation,
  compatibility checks, and a comparison against the preconditioner callback.
- [x] Add automatic diagonal inverse-Jacobian preconditioner estimation,
  refresh intervals, zero/tiny diagonal safeguards and regression coverage
  for analytic and finite-difference Jacobian-vector products.
- [x] Add low-cost deterministic Walsh-sign probing for automatic diagonal
  estimation, with an exact-mode default, finite-difference accounting,
  explicit tests for long-range aliasing, and large-system benchmarks.
- [x] Add seeded Rademacher probing for dense/irregular Jacobian estimates,
  plus optional exact-coordinate verification, split-probe consistency
  checking and exact fallback on detected estimator bias or instability;
  document that unverified random probes can have large coordinate errors.
- [x] Add progressive verified Rademacher diagonal-probe budgets with
  exact fallback at the configured limit, retaining all previous Jv products.
  Benchmark dense nonlocal couplings and report maximum as well as average
  diagonal-estimation errors; verification remains a heuristic, not a proof.
- [x] Rotate progressive exact-coordinate spot checks by bit-reversed
  within-stratum offsets, revalidate cached coordinates after each expansion,
  and reuse all checks during exact fallback. Document missed-bias tests and
  evaluation-count tradeoffs without changing fixed-probe defaults.
- [x] Prevent narrow-stratum progressive verification from reusing all
  previously checked coordinates when the bit-reversed offset rounds to
  zero. Prefer unchecked coordinates via a local successor table, and add
  dense-coupling regressions for missed bias and finite-difference work.
- [x] Add optional checked mean-field rank-one Woodbury preconditioning
  for dense row-uniform coupling, with exact diagonal estimation, two
  independent structural verification probes, safe diagonal fallback,
  ill-scaled-system and one-step Krylov regression coverage.
- [x] Generalize the checked Krylov rank-one model to a caller-supplied,
  fixed right factor `v` in `J = diag(a) + u*v^T`, with normalized
  mixed-sign weights, structural checks, and safe diagonal fallback.
- [x] Extend exact-diagonal checked Krylov preconditioning to a known
  two-column right coupling factor, with row-identifiability checks,
  independent structural probes, guarded 2-by-2 Woodbury inversion,
  and diagonal fallback on incompatible or ill-conditioned systems.
- [x] Add an opt-in three-probe Jacobian-diagonal sketch for the known
  rank-two right factor, with rowwise identifiability checks, the existing
  four-product structural verification, exact fallback, and scaled
  nonlinear benchmark coverage for diagonal/rank-one/rank-two policies.
- [x] Add an opt-in full-coordinate Jacobian diagnostic that identifies
  an unknown rank-one right factor from off-diagonal observations,
  exhaustively checks the fitted structure, and reuses the n existing
  exact diagonal-coordinate products. Document O(n^2) memory and
  safe fallback when the pattern is non-identifiable or inconsistent.
- [x] Add a first opt-in discovery path for unknown rank-two right
  factors using n coordinate Jacobian-vector products, masked low-rank
  completion of off-diagonal entries, exhaustive verification and guarded
  Woodbury inversion, with diagonal fallback and focused regressions.
- [x] Add an opt-in O(n)-memory unknown rank-two factor discovery mode
  that streams n coordinate Jacobian products, verifies every off-diagonal
  entry and retains exact diagonal fallback. Include non-pilot contamination
  and large ill-scaled regression cases.
- [x] Add an opt-in four-product, O(n)-storage rank-two fit when the
  uncoupled diagonal is known and the unknown coupling is symmetric,
  with independent held-out checks and exact-diagonal fallback.
- [x] Extend low-product unknown-factor discovery beyond symmetric
  coupling using a supplied analytic transpose-Jacobian action, a known
  uncoupled diagonal, four forward and two transpose products, and
  guarded Woodbury inversion with exact-diagonal fallback.
- [x] Benchmark dense, streamed, and nonsymmetric low-product rank-two
  discovery with forward/transpose counts, weighted transpose costs,
  rank-three fallback cases and strongly scaled problems. Add opt-in
  held-out transpose verification without changing the default API.
- [x] Add cost-aware one- or two-direction adjoint checks for nonsymmetric
  low-product rank-two discovery, with stronger verification for poorly
  conditioned fits and reproducible forward/transpose cost benchmarks.
- [x] Validate nonsymmetric rank-two preconditioning on nonlinear
  systems with varying Jacobians and add an opt-in `D(x)` callback evaluated
  at each scheduled refresh, preserving fixed-array calls, exact-diagonal
  fallback, and explicit invalid-callback errors.
- [x] Benchmark fixed preconditioner refresh intervals on genuinely
  nonlinear rank-two systems with evolving diagonal and coupling factors,
  comparing convergence, forward/transpose counts, relative adjoint costs
  and CPU time. Verify the callback schedule at intervals 2/3/4 and the
  behavior of delayed callback failures. Do not change the default interval
  based on a workload-specific benchmark.
- [x] Add an opt-in, residual-progress-triggered refresh policy with a
  maximum preconditioner age, unchanged defaults, explicit parameter
  validation, and a deterministic scheduling regression.
- [x] Benchmark work-aware refresh on nonlinear systems where stale
  preconditioners increase built-in Krylov effort. Compare fixed, stall-only,
  and growth/cost-aware policies over 108 cases, test actual early refresh
  and expensive-rebuild suppression, and document the tradeoff between fewer
  forward products and more transpose products. No universal speedup or
  change to default scheduling is justified by these results.
- [ ] Investigate a truly cost-aware nonlinear refresh policy that accounts
  for both stale-preconditioner quality and refresh expense. The initial
  progress-only policy did not improve modeled work on the 84-case nonlinear
  benchmark; explore inner Krylov work and scaling signals before automating
  refresh decisions further.
- [ ] Investigate *fewer-product* unknown-factor discovery under stronger
  structural assumptions or an adjoint Jacobian callback. Arbitrary
  unknown right factors cannot generally be identified from a fixed set
  of fewer than n coordinate probes; the O(n)-memory mode still requires
  n products and O(n^2) verification arithmetic.
- [ ] Investigate stronger adaptive LGMRES controls for weighted/ill-scaled
  residuals, including robust preconditioners for dense/nonlocal couplings,
  and sufficient nonlinear iteration budgets without weakening stopping
  criteria.
- [ ] Complete a dependency-backed `fpm test` on Windows and additional
  compiler checks when the full FPM dependency set is available.


### Krylov refresh timing follow-up (October 2026)

Wall-clock refresh is available as an opt-in alternative to CPU-timed
refresh. The serial/OpenMP quadrature benchmark now covers fixed,
product-count, CPU-timed, and wall-timed policies. Compare repeated runs
on independent hardware (including Windows and macOS), with different
OpenMP runtimes, thread counts, clock resolutions, and callback costs.
A wall-timed trigger is not guaranteed to minimize elapsed time; short
measurement windows can be noisy. Full dependency-backed `fpm test` is
still required on a system with the Fortran Package Manager installed.

## Milestone 13: signal processing

- [x] Precompute Savitzky-Golay center and `interp` boundary polynomial
  weights once per multidimensional call, and add rank-4 real/complex input
  with one-based axis selection. Validate 400 additional SciPy 1.17 cases
  and benchmark multichannel speed against the previous implementation.
- [x] Extend `savgol_filter` to real/complex rank-1/2/3 inputs,
  with selected one-based axes, boundary modes, and derivative options;
  validate 600 seeded SciPy 1.17 cases.
- [x] Promote all mixed real/complex numerator, denominator and signal
  inputs in `lfilter` for rank 1/2/3, including complex initial/final
  streaming states; validate 540 seeded SciPy 1.17 cases.
- [x] Extend `fftconvolve` to ranks 2 and 3 with mode/axis control,
  real/complex/mixed inputs, and broadcasting; add 1,000 seeded SciPy 1.17
  comparisons. Preserve the 1D overlap-add interface.
- [x] Extend `lfilter` to real/complex rank-2 and rank-3 signals, with
  one-based axis selection and per-channel initial/final states; validate
  700 randomized cases against SciPy.
- [x] Add rank-two and rank-three `sosfilt`, `filtfilt`, and `sosfiltfilt`
  overloads, with one-based axes and SciPy-compatible SOS state arrays;
  compare 480 randomized real/complex cases against SciPy 1.17.
- [ ] Add higher-rank and blocked overlap-add convolution, and consider
  performance improvements for multidimensional IIR filtering.

## Milestone 13: one-dimensional signal processing

- [x] Add real/complex SciPy-style continuous-to-discrete state-space, TF
  and ZPK conversions with `cont2discrete`/`cont2discrete_zpk` including
  zero-order hold, first-order hold, generalized bilinear (Tustin/Euler/
  backward difference), and impulse invariant methods. Validate using
  deterministic Fortran tests and 390 randomized SciPy 1.17 comparisons.
- [ ] Add higher-precision or balanced state-space realizations for very
  high-order/ill-conditioned discretization problems.

- [x] Add `tf2ss`, `ss2tf`, `zpk2ss` and `ss2zpk` with real/complex
  SISO models, common-denominator SIMO numerator matrices, selected MIMO
  input/output channels and a canonical form consistent with SciPy 1.17.
  Added 832 seeded matrix/polynomial comparisons and checked Fortran tests.
- [ ] Investigate numerically stable high-order characteristic-polynomial
  computation for `ss2tf` and additional general LTI system operations.
- [x] Add `tf2zpk`, `zpk2tf`, and digital `zpk2sos` conversion APIs
  with complex coefficients, root reconstruction, real-conjugate SOS pairing,
  and 840 seeded SciPy comparisons. Document digital degree padding and the
  distinction from analog `zpk2sos_analog`.
- [x] Extend digital `zpk2sos` to `minimal` pairing without introducing
  artificial origin roots, and add `sos2zpk` with SciPy-compatible two-root
  padding per section. Tested against SciPy 1.17 on 792 minimal designs and
  759 SOS inverses; verified direct rational responses for minimal sections.
- [x] Add `tf2sos` for real/complex BA coefficients, analog/digital
  operation and SciPy-compatible pairing selection. Resolve a `keep_odd`
  pairing failure for severely unequal zero/pole degrees; validate 992
  randomized valid designs and 64 invalid minimal designs.
- [ ] Further refine pole-zero pairing, root accuracy, and section gain
  distribution for high-order and clustered-root systems.

- [x] Add analog/digital partial fractions with `residue`, `residuez`,
  `invres`, and `invresz`, including grouping of repeated roots, complex
  coefficients, a standalone Fortran regression, and 540 SciPy comparisons.
- [ ] Improve scaling and root refinement for very high-order or tightly
  clustered poles in partial-fraction decomposition.

- [x] Add weighted real `scipy.signal.lombscargle` with optional floating
  baseline, normalized periodograms, complex amplitude extraction, and
  `scipy.signal.vectorstrength` for event-phase synchrony. Compare against
  SciPy 1.17 over 900 randomized irregular-sampling cases.
- [ ] Extend irregular spectra to batched observation sets and explore fast
  periodograms for regularly spaced trial-frequency grids.

- [x] Add `scipy.signal.sawtooth`, `square` (scalar and variable duty/width),
  `unit_impulse` (1D and 2D), and `max_len_seq` (2-32 bits with custom taps,
  reusable state and streaming), plus SciPy comparisons.
- [x] Add logarithmic-time forward/backward skip-ahead for `max_len_seq`,
  including signed `int64` distances, arbitrary nonzero seeds, custom taps,
  zero-length state-only jumps, and exact SciPy reference comparisons.
- [ ] Extend unit impulses to higher ranks, support scalar-time waveform
  generics, and add more signal-sequence generation algorithms.

- [x] Add `cheby2`, `iirfilter`, `iirnotch`, `iirpeak`, and digital
  lowpass/highpass `buttord`, `cheb1ord`, `cheb2ord`, and `iirdesign`.
  Preserve MIT-licensed independent implementations, SciPy 1.17 reference
  comparisons, SOS/BA output and explicit validation statuses.
- [x] Add analog `freqs` and `freqs_zpk`, and digital `freqz_zpk`,
  supporting automatic frequency grids, explicit frequencies, complex poles
  and zeros, and SciPy randomized comparison tests.
- [x] Add analog `butter`, `cheby1`, `cheby2`, and `iirfilter` BA designs
  (lowpass, highpass, bandpass, bandstop) using rad/s cutoff frequencies.
- [x] Extend `buttord`, `cheb1ord`, `cheb2ord`, and `iirdesign` to analog
  minimum-order selection and digital/analog bandpass and bandstop; retain
  scalar-call compatibility and add `iir_order_result%cutoffs`.
- [x] Add elliptic (Cauer) IIR design (`ellip`/`ellipord`) and
  `iirfilter`/`iirdesign` dispatch, with digital BA/SOS and analog BA,
  all four filter shapes, and randomized SciPy comparisons.
- [x] Add analog ZPK/SOS output formats for all four classical IIR design
  families using `iirfilter_zpk` and `iirfilter_sos`, and add `zpk2sos_analog`
  plus `sosfreqs`. Preserve existing analog BA behavior and distinguish
  analog s-plane sections from digital delay-polynomial sections.

- [x] Add 1D real/complex `scipy.signal.upfirdn`, `resample_poly`, and
  FIR `decimate`, including polyphase up/down sampling, Kaiser(5)
  lowpass design, optional user FIR taps, and boundary handling.
- [x] Add `cheby1`, default Chebyshev-I IIR `decimate`, and real/complex
  periodic Fourier `resample` with FFT Nyquist-bin handling.
- [x] Add axis-aware rank-2/rank-3 real/complex `resample`, `resample_poly`,
  `upfirdn`, and `decimate` overloads, including negative axis selection.
- [ ] Extend multirate resampling past rank 3, add additional boundary
  modes and complex-valued user FIR taps for `resample_poly`, and accelerate
  Fourier resampling at non-power-of-two lengths.

- [x] Add `firwin2`, weighted `firls`, and homomorphic/Hilbert
  `minimum_phase` FIR design/conversion with SciPy reference comparisons,
  fixed tests, and independent transform and normal-equation implementations.
- [x] Add `remez` equiripple FIR design for standard symmetric, Hilbert,
  and differentiator modes, with stable barycentric exchange, fixed tests,
  SciPy comparisons, and explicit nonconvergence status.
- [x] Add an independent rank-revealing Jacobi pseudoinverse fallback,
  stable narrow-band cosine integration, and weight scaling for `firls`;
  validate conditioned designs against SciPy.
- [x] Improve high-order Remez reliability using extended-precision
  intermediates when available, linear-time extrema selection, and eight
  previously failing regression cases. All 240 tested practical designs now
  converge with GNU Fortran 14.2.
- [x] Prefer an available >=18-digit intermediate real kind when full
  binary128 is unavailable for Remez exchange and reconstruction.
- [ ] Further improve convergence for extremely narrow-band or very high-order
  Remez specifications, including portable dp-only compilers.
- [x] Add `czt` and `zoom_fft` transforms, with an FFT-based unit-contour
  implementation, general complex spiral support, and SciPy comparisons.
- [x] Add `czt_into` and `zoom_fft_into` through rank 15 for real and
  complex arrays with arbitrary one-based/negative axes, fixed-shape caller
  outputs, transactional error handling and prepared-plan variants. Compare
  466 independent SciPy 1.17 cases and test rank-15/strided arrays.
- [ ] Investigate efficient nonunit-circle CZT scaling for very large
  output lengths without extreme cancellation or overflow.

- [x] Add `scipy.signal.firwin`, `freqz`, `sosfreqz`, and `butter`, including
  low/high/bandpass/bandstop designs, complex frequency responses, and
  low-order transfer-function and second-order-section outputs. Validate
  against SciPy 1.17 and preserve dependency-free standalone tests.

- [x] Add `argrelmax`, `argrelmin`, `argrelextrema` and the
  separately marked `peak_detect` excursion-based extension, inspired by
  Christopherson's MIT `peaks` and Eli Billauer; verify 1,800 randomized
  SciPy/reference comparisons, without another FPM dependency.
- [ ] Assess narrowly scoped numerical algorithms from `dynamics` (e.g.
  modal response/FRF fitting) against SciFort's existing LTI functions;
  avoid adding its full heavy dependency graph unnecessarily.
- [x] Add selected axes and mean/median Welch averaging for rank-1..15
  periodogram/Welch/CSD/coherence, reusing existing spectral kernels.
  Continue with other unimplemented time-frequency variants as needed.
- [x] Add the `scifort_signal` module with `find_peaks`, `peak_prominences`,
  and `peak_widths`, including flat-peak handling, height/threshold/distance/
  prominence/width constraints, one-based indices, and independent tests.
- [x] Add `detrend` with constant, linear, and breakpoint-segmented fits.
- [x] Add direct real/complex `convolve`, `correlate`, and `correlation_lags`,
  with `full`, `same`, and `valid` modes and input validation.
- [x] Add direct rank-two `convolve2d`/`correlate2d` for real, complex
  and mixed inputs, all three modes and SciPy fill/wrap/symm boundaries;
  test 1,600 seeded cases against SciPy 1.17.
- [x] Add real/complex 1D `fftconvolve` and `oaconvolve`, including long-signal
  blocked overlap-add and SciPy 1.17 reference comparisons.
- [x] Add N-dimensional direct and FFT convolution/correlation with
  ranks through 15, selected axes, and singleton broadcasting.
- [x] Extend `find_peaks` to position-dependent height, threshold, prominence,
  width, and plateau-size lower/upper bounds with independent SciPy comparisons.
- [x] Support SciPy open-interval `find_peaks` property requests using five
  `*_open` logical flags, preserving filters, and verify 2,400 SciPy cases.
- [x] Add `savgol_irregular` for weighted local-polynomial smoothing and
  derivatives on strictly increasing nonuniform observation times, including
  complex values, arbitrary query times, singular-window diagnostics and
  independent NumPy/SciPy numerical reference comparisons.
- [x] Add 1D `savgol_coeffs` and `savgol_filter` for smoothing,
  derivatives, even windows and five SciPy-style edge modes.
- [x] Add 1D real/complex `lfilter` and `lfilter_zi` with direct-form-II
  streaming state and step-steady-state initialization.
- [x] Add `lfiltic` for real/complex direct-form-II past-state
  initialization and `group_delay` for real/complex digital filters,
  including SciPy 1.17 reference tests and singular-frequency diagnostics.
- [x] Add the Gustafsson initial-state method to rank-1/2/3 `filtfilt`
  for real and complex signals, with optional truncated impulse response,
  rank-deficient minimum-norm least squares, and SciPy reference checks.
- [x] Extend Gustafsson initialization to assumed-rank `filtfilt_nd`
  (ranks 1-15, all real/complex promotions) with optional `irlen`, 960
  rank/axis/type reference cases and short-signal/strided-array regressions.
- [x] Improve the underdetermined Gustafsson least-squares solve with
  a pivoted Householder LQ-equivalent QR factorization, retaining SVD/QR
  fallbacks for rank-deficient systems. Four previously rejected short
  high-order SOS cases are now resolved; four severely ill-conditioned
  12th-order cases remain explicitly rejected.
- [x] Add an SOS-native Gustafsson alternative for high-order zero-phase
  filtering using cascade state responses (no `sos2tf` expansion), with
  ranks 1-15, independent SciPy SOS-state-space validation, and a clear
  numerical-failure status for ill-conditioned short signals. SciPy itself
  exposes only the padding method in `sosfiltfilt`.
- [x] Add real/complex 1D `sosfilt`, `sosfilt_zi`, `filtfilt`, and
  `sosfiltfilt`, with SOS streaming states, reflected boundary padding,
  step-steady-state initial conditions, and SciPy 1.17 reference validation.
- [x] Add one-sided 1D `periodogram`, `welch`, `csd`, `coherence`, `stft`,
  `spectrogram`, and `get_window`, with Welch averaging, periodic windows,
  density/spectrum scaling, and SciPy reference comparisons. Attribution to
  the MIT-licensed `spectrum` library is recorded in provenance and licenses.
- [x] Add complex-valued two-sided `periodogram`, `welch`, `csd`,
  `coherence`, `stft`, and `spectrogram`, preserving the original real API.
- [ ] Add two-sided spectra for real inputs, mixed-type cross spectra,
  multidimensional axes, more STFT boundary modes and flexible Tukey parameters.
- [ ] Add second-order-section (`sosfilt`), forward-backward (`filtfilt`),
  and stable IIR/FIR filter design algorithms.

## 2026-11-08 — Low-memory selected Hermitian band eigenpairs

- [x] Add selected `eigvals_banded` eigenvalues through banded Sturm/inertia.
- [x] Add selected `eig_banded` eigenvectors via pivoted band LU inverse
  iteration and two-pass reorthogonalization without a dense input matrix.
- [x] Retain dense full-spectrum behavior and selected dense fallback.
- [x] Permit opt-in full-spectrum inverse iteration without a dense input
  expansion; specialize real vectors to real work arrays, not complex n-by-n.
- [x] Add opt-in full-spectrum band-to-tridiagonal Givens bulge chasing
  without dense matrix expansion; retain the former default and fallback.
- [ ] Improve near-singular inertia robustness and enable a proven
  automatic selection policy for full-spectrum banded eigensolvers.

## 2026-11-08 — Selected symmetric band eigenvectors

- [x] Compute selected eigenvectors via pivoted banded inverse iteration,
  without allocating the full dense matrix.
- [x] Specialize real symmetric systems to real LU factors and solves.
- [x] Implement optional band-to-tridiagonal bulge chasing for full-spectrum
  eigenvectors, including complex Hermitian bands.
- [ ] Improve tightly clustered eigenpair handling and specialized real QR.

## 2026-11-08 — Generalized nonsymmetric eigenproblem groundwork

- [x] Add complex generalized Schur factors for regular small matrix pairs
  by inverse-based Schur/QR, including singular-B and some both-singular cases.
- [x] Compare homogeneous eigenpairs and reconstruction with SciPy QZ;
  test both Fortran build modes and document the numerical limitations.
- [x] Implement initial inverse-free two-sided complex QZ iteration for small
  pencils, preserving the inverse-based method as a documented fallback.
- [ ] Harden complex QZ for arbitrary regular pencils with zero-diagonal
  chasing, exceptional shifts, higher-order support, and LAPACK-grade safeguards.
- [x] Add generalized nonsymmetric `eig`/`eigvals` with finite and infinite
  eigenvalues, homogeneous output, and optional left/right eigenvectors.
- [x] Add complex-output `ordqz` with homogeneous eigenvalue sorting.
- [ ] Add proper real quasi-upper-triangular generalized Schur output.

## 2026-11-09 -- Pivoted indefinite dense factorization

- [x] Add `scipy.linalg.ldl`-style real symmetric, complex Hermitian and
  complex symmetric factorization with Bunch-Kaufman 1x1/2x2 pivots.
- [x] Support upper/lower input and triangularized output, singular zero
  columns, and one-based permutations matching SciPy's factor convention.
- [x] Independently validate 2,400 checked/release cases against SciPy.
- [ ] Consider blocked BLAS acceleration, selected-pivot condition diagnostics,
  and a separate factor-based solve API once numerical testing warrants it.

## 2026-11-12 -- Approved native boundary-value integration expansion

The user approved `solve_bvp` as the next native SciPy-coverage task after
FFT/Hankel work, while continuing to defer C ABI and wrapper expansion.

- [x] Real and complex fourth-order two-point collocation with unknown parameters.
- [x] Optional singular endpoint term, analytic/estimated Jacobians, adaptive
  meshes, dense cubic output, and explicit convergence diagnostics.
- [x] Native analytic tests and checked/optimized SciPy/analytic comparisons.
- [ ] Specialized large-system factorization and continuation interfaces, if
  warranted by numerical and performance evidence.

## 2026-11-13 -- Approved native Lebedev quadrature expansion

The user approved all 32 Lebedev unit-sphere rules as the next native
SciPy-coverage task. C ABI and wrapper expansion remain deferred.

- [x] Pure Fortran node/weight API with all algebraic orders 3-131.
- [x] Explicit license review, preserved coefficients, attribution, and notices.
- [x] Validate all nodes, weights, symmetry, and inequivalent polynomial moments
  in checked and optimized builds, including signed-weight rules.
- [ ] Reconsider a direct original Fortran import only after clear permissive
  modification and redistribution permission is established.

## 2026-11-13 -- Approved native AAA approximation expansion

The user approved AAA rational approximation as the next native SciPy-coverage
task after Lebedev quadrature. C ABI and wrapper expansion remain deferred.

- [x] Real/complex sample sets, greedy support selection, native SVD fitting,
  scalar/vector evaluation, and explicit convergence diagnostics.
- [x] Finite poles/zeros, simple-pole residues, and automatic/manual cleanup,
  with explicit size and numerical-conditioning limits.
- [x] Analytic and SciPy validation, scale/offset tests, and retained notices.
- [x] Broaden spectral calculations to 512 support points with a dedicated,
  validated Reference-LAPACK QZ backend for normalized AAA arrowhead pencils.

## 2026-11-13 -- Approved native scalar smoothing expansion

The user approved native automatic-knot `splrep(s>0)` after AAA spectra.
C ABI and wrapper expansion remain deferred.

- [x] Nonperiodic scalar automatic-knot smoothing, degrees 1-5, weighted
  residual targets, and enclosing support bounds through the existing TCK API.
- [x] Review pinned upstream source/license, retain notices, and record adaptation.
- [x] Validate checked/optimized native behavior against SciPy, the separately
  compiled upstream core, and the upstream knot-insertion crash fixture.
- [x] Extend the foundation to native nonperiodic parametric smoothing.
- [x] Add periodic scalar and parametric interpolation/smoothing.
- [ ] Add modern spline constructors.

## 2026-11-13 -- Approved native parametric smoothing expansion

The user approved nonperiodic automatic-knot `splprep(s>0)` after scalar smoothing.

- [x] Shared knots and total weighted residual across channels, degrees 1-5,
  supplied/chord-length parameters, enlarged bounds, and existing tuple layout.
- [x] Preserve higher channel counts, add scaling/range/residual diagnostics,
  and retain upstream source/license attribution.
- [x] Validate checked/optimized behavior against the pinned core, upstream
  parametric data, SciPy residual contracts, and independent normal equations.
- [x] Periodic parametric interpolation and smoothing.
- [ ] Continuation and modern constructors.

## 2026-11-13 -- Approved native periodic scalar fitting

The user approved native `splrep(per=.true.)` after nonperiodic parametric
smoothing. C ABI and wrapper expansion remain deferred.

- [x] Periodic interpolation and automatic-knot smoothing for degrees 1-5,
  weights, ignored closing ordinate/weight, and existing padded tuples.
- [x] Adapt pinned licensed periodic kernels, retain notices, and validate
  physical residuals, seam continuity, ranges, and existing dependent APIs.
- [x] Check original upstream, its analytic sine/cosine fixture, SciPy
  contracts, independent normal equations, and checked/optimized builds.
- [x] Periodic parametric fitting.
- [ ] Periodic fixed knots and remaining modern selection/capacity options.

## 2026-11-13 -- Approved native periodic parametric fitting

The user approved native `splprep(per=.true.)` after periodic scalar fitting.
C ABI and wrapper expansion remain deferred.

- [x] Closed-curve interpolation and automatic-knot smoothing, degrees 1-5,
  shared knots, weighted total residuals, and preserved higher channel counts.
- [x] Supplied/chord-length parameters, ignored closing row/weight, periodic
  endpoint derivatives, range/residual checks, and retained BSD attribution.
- [x] Checked/optimized native regressions, pinned original core and closed
  curve data, SciPy comparisons, and independent cyclic normal equations.
- [ ] Periodic fixed knots, continuation, and remaining modern selection/capacity options.

## 2026-11-13 -- Approved modern scalar constructor

The user approved `make_splrep` after native periodic parametric fitting.

- [x] Native B-spline object result, default interpolation, automatic smoothing,
  positive weights, and supplied full clamped knots with simple interiors.
- [x] Fixed-knot smoothing, polynomial limits retaining knots, physical residual
  verification, and explicit capacity/invalid-input/numerical statuses.
- [x] Checked/optimized native and independent reference/equation validation.
- [x] Native `make_splprep`; see the approved parametric constructor below.
- [x] Periodic constructor interpolation and automatic smoothing.
- [ ] Broader knot forms,
  and modern SciPy automatic-selection/soft-capacity parity.

## 2026-11-14 -- Approved modern parametric constructor

The user approved native `make_splprep` after the scalar constructor.
C ABI and wrapper expansion remain deferred.

- [x] Native array B-spline, supplied/chord-length parameters, interpolation,
  automatic/fixed-knot smoothing, shared knots, and one residual target.
- [x] Polynomial limits retaining supplied knots, explicit hard capacity,
  numerical guards, and checked/optimized numerical validation.
- [x] Periodic constructor interpolation and automatic smoothing.
- [ ] Broader knot forms, exact modern
  automatic selection, and soft capacities remain separate gaps.

## 2026-11-14 -- Approved periodic native constructors

The user requested continued native SciPy gap work after `make_splprep`.
C ABI and wrapper expansion remain deferred.

- [x] Closed-input `bc_type='periodic'` interpolation and automatic smoothing
  for both constructors, degrees 1-5, shared knots, and hard capacities.
- [x] Normalized interpolation, physical verification, checked/optimized
  regressions, SciPy interpolation and independent cyclic equation checks.
- [ ] Supplied periodic knots, broader forms, exact modern automatic-selection
  parity, and soft capacities remain gaps.
