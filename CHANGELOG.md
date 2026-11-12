## 2026-11-12 - Censored Gumbel likelihoods and fitting

- Add `gumbel_r` and `gumbel_l` censored log-likelihood and negative
  log-likelihood, with stable exponential increments for interval masses.
- Extend both Gumbel fitting generics to censored samples with bounded
  location and scale, preserving the existing array APIs.
- Share kernels through reflection and reuse logistic censoring's checked
  standardization helper. Add native and high-precision reference checks.

## 2026-11-12 - Censored logistic likelihoods and fitting

- Add logistic censored log-likelihood and negative log-likelihood with stable
  interval masses, extreme tails, and overflow-safe standardization.
- Extend `logistic_fit` to censored samples with bounded location and scale,
  preserving the existing array API. C ABI work remains deferred.
- Add native regression and independent high-precision reference checks.

## 2026-11-12 - Censored beta likelihoods and fitting

- Add beta censored likelihoods and bounded shape fitting with fixed location
  and scale, excluding exact endpoint samples from the fitting overload.
- Preserve independently formed endpoint distances, tiny normalized distances,
  and narrow interval masses; retain singular endpoint densities in likelihoods.
- Prevent intermediate overflow in negative powers in the shared gamma-ratio
  remainder, found by the extreme-shape beta checks.
- Add high-precision references, SciPy fitting comparisons, endpoint and uniform
  identities, numerical regressions, documentation, and provenance.

## 2026-11-12 - Censored gamma likelihoods and fitting

- Add gamma censored log-likelihood and negative log-likelihood, including
  support endpoints, narrow intervals, and scaled arguments below binary64 range.
- Extend `gamma_fit` to censored samples with fixed location; exact support
  endpoints are excluded from this initial fitting overload.
- Handle opposing infinite censored likelihood contributions as quiet NaN.
- Add high-precision references, fitting comparisons, tests, and documentation.

## 2026-11-12 - Censored lognormal likelihoods and fitting

- Add lognormal censored log-likelihood and negative log-likelihood with stable
  log-domain distances and narrow interval probabilities.
- Extend `lognormal_fit` to censored observations with fixed location and
  bounded shape/scale, preserving its existing array interface.
- Add native, high-precision, and SciPy fit comparisons and document limits.

## 2026-11-12 - Censored likelihoods and fitting

- Extend the existing censored-data type with left and interval observations.
  Add stable normal, exponential, and Weibull censored likelihoods and bounded
  fitting overloads, preserving existing uncensored calls.
- Add an analytic fixed-location right-censored exponential MLE, including
  explicit failure reporting for absent or unrepresentable finite estimates.
- Reject unsupported censoring groups in ECDF/logrank with NaN results.
- Add high-precision likelihood references, SciPy fit comparisons, native
  numerical regressions, documentation, and provenance.

## 2026-11-11 - External package acknowledgments

- List external packages in the README, distinguishing FPM dependencies,
  adapted implementations, optional backends, and design/testing references.
- Link to upstream projects and the detailed provenance and license records.

## 2026-11-11 - Documentation validation policy

- Exempt documentation-only and contributor-instruction edits from full FPM
  builds and tests. Require focused content, link, formatting, and applicable
  example checks instead; retain full tests for executable and build changes.

## 2026-11-11 - Project overview and historical notes

- Replace the README's accumulated announcements with a project overview,
  build and test instructions, a Fortran usage example, and documentation links.
- Preserve the previous README in `README_HISTORY.md`, with its relative links
  intact and a notice distinguishing historical reports from current API guides.

## 2026-11-11 - CPU-timed Krylov test portability

- Increase the synthetic CPU workloads and their iteration allowance so the
  timed-refresh regression works with coarse Windows CPU timers and fast builds.
- Assert that workloads reach their requested duration instead of silently
  accepting the iteration cap. Solver behavior and public APIs are unchanged.

## 2026-11-11 - Repository refresh from scifort-main

- Refresh code, tests, documentation, build configuration, validation scripts,
  and license/provenance records from the sibling scifort-main working tree.
  Include its uncommitted streaming spectrogram default-window fix.
- Retain the current Git history, local notes, batch helpers, and existing build directory.
- Normalize the Frechet/Matern spellings in two Fortran files to ASCII in
  comments and a test status message; numerical code and license notices are unchanged.

## 2026-11-11 - Streaming spectrogram default window

- Correct the omitted streaming spectrogram window from Tukey alpha 0.5 to
  alpha 0.25, matching batch spectrograms and SciPy's documented default.
- Explicit `window='tukey'` continues to select alpha 0.5, and custom window
  weights retain their supplied values. No public signatures change.
- Add 48 native batch/stream comparisons and expand independent SciPy
  comparisons to cover omitted windows and explicit Tukey names.

## 2026-11-11 - Unattended FPM test discovery

- Move 20 numerical comparison producers from `test/` into
  `validation/reference_drivers/`, preserving their source and input protocols.
  Ordinary FPM tests no longer launch programs waiting for comparison input.
- Update comparison-script, generator, documentation, and provenance paths.
- Document the test/validation split and check for explicit console reads in
  automatically discovered tests in CI.

## 2026-11-11 - Optional callback compilation

- Pass the optional nonlinear callback explicitly to the internal violation
  helper, retaining its declared interface with GNU Fortran's
  `-Werror=implicit-interface` check. Public interfaces are unchanged.
- Extend constrained differential evolution tests with a fixed-bound nonlinear
  constraint and evaluation accounting.
- Apply the same explicit callback forwarding to the real and complex
  matrix-free rank estimator and spectral norm helpers, and to dual
  annealing's best-point notification helper and SHGO's nonlinear value and
  Jacobian helpers.

## 2026-11-10 — Dense mixed-integer linear programming (`milp`)

- Add `scifort_optimize_milp` with an independent best-bound branch-and-bound
  algorithm built on SciFort's original `linprog` simplex solver.
- Implement SciPy-style variable types 0/1/2/3, two-sided linear constraint
  bounds, incumbent tracking, inherited LP dual bounds, and relative gap,
  node, CPU-time, and LP-iteration stopping options.
- Preserve correct status distinctions for infeasible, invalid, and
  relaxation-unbounded integer problems; independently re-check feasibility.
- Add umbrella exports, 78 native regression assertions, 320 randomized
  comparisons per compiler build with SciPy, 220 exhaustive integer reference
  models, an example, and API/validation documentation.
- Correct the existing `linprog` Phase I minimum-ratio tie comparison on
  large-magnitude RHS values; add two large-scale simplex regression checks.
  Standalone LP references pass 1,450 SciPy comparisons in both build modes.
- No new third-party code or FPM dependency. Full `fpm test` remains
  unverified where FPM is unavailable.

## 2026-11-10 — Quadratic assignment: FAQ and 2-opt

- Add `scifort_optimize_quadratic_assignment` with a pure approximate solver
  for minimum or maximum quadratic assignment, including seeded fixed matches.
- Implement deterministic Frank-Wolfe FAQ relaxation using the existing
  Hungarian assignment solver; support user-supplied doubly stochastic `p0`.
- Implement best-improvement 2-opt with O(n) per-pair swap deltas and
  optional initial matching guesses; return one-based permutations and status.
- Add public umbrella exports, a runnable 4x4 example, 469 native assertions,
  and 288 independent comparison cases (5,156 checks) against SciPy 1.17.0.
- No new library dependencies. `fpm test` remains unverified here because the
  Fortran Package Manager is not installed.

## 2026-11-10 — SHGO independent nonlinear constraint objects

- Add public `shgo_nonlinear_constraint` with `fun`, optional `jac`, and
  allocatable `lb`/`ub` arrays; both `shgo` and `shgo_core` accept
  `constraints(:)` as a new **optional trailing argument**.
- Concatenate each block's components while preserving existing linear
  constraints and equality/inequality semantics. Reject invalid bounds,
  unassociated procedures, missing arrays, and simultaneous legacy
  `nonlinear=` arguments before any callback invocation.
- Apply analytic Jacobians for individual blocks; use bound-aware
  finite differences **only** for blocks lacking analytic derivatives in
  both restoration and SLSQP's canonical equality/inequality layout.
- Add 40 native regression assertions, five independently feasible SciPy
  comparison cases, extended SLSQP adapter smoke tests, a public example,
  documentation and provenance details. Checked and optimized builds pass.
- No source code was copied from SciPy and no new FPM dependencies were added.
  Actual SLSQP integration still requires a real FPM-backed build.

## 2026-11-10 — SHGO analytic constraint Jacobians and feasible-point recovery

- Accept optional `nonlinear_jac(x,a)` in the core and public SHGO adapter.
  Convert linear/nonlinear equalities, lower and upper bounds into exact
  equality-first, inequality-second Jacobian rows for SciFort SLSQP.
- Add opt-in `restore_feasibility`, `restoration_starts`, and
  `restoration_maxiter` options, with deterministic damped projected
  Gauss-Newton attempts and finite-difference fallback for nonlinear Jacobians.
- Verify all accepted restored points with the actual objective, count its
  calls strictly, and expose `ncjev`, `nrestore`, and `nrestored`.
- Add 45 checked/optimized native assertions, five independent SciPy SLSQP
  feasibility comparisons, adapter derivative-layout checks, and an example.
- Preserve old defaults; restoration is heuristic, not a global convergence
  or feasibility certificate. Actual FPM/SLSQP integration remains unverified.

## 2026-11-10 — Constrained SHGO global optimization

- Add linear/nonlinear constraint blocks with lower/upper bounds and equality
  handling to `shgo` and `shgo_core` without changing unconstrained defaults.
- Use feasibility-first mesh ranking, strictly accept only feasible returned
  minima, skip objective work at infeasible mesh vertices, and reuse existing
  SciFort SLSQP for constrained local polishing.
- Add `shgo_infeasible`, `ncev`, `maxcv`, and `feasible` result diagnostics.
- Add constrained native assertions, independent SciPy benchmarks, a complete
  constrained example, and an expanded SLSQP public-adapter smoke test.
- No new dependencies, bundled external code, or license obligations.

## 2026-11-10 — Conforming three-dimensional adaptive SHGO tetrahedralization

- Extend opt-in `adaptive=.true., conforming=.true.` from 2D to 3D without
  changing the public argument list or uniform/nonconforming default paths.
- Construct conforming octree-face triangulations using cached dyadic face
  centers and exact-integer triangulation of partially subdivided faces;
  cone each face triangle to a sampled leaf center.
- Validate positive tetrahedron volumes, exact coverage of each leaf face,
  triangle incidence counts (two per interior face, one on the exterior),
  and strict atomic objective/vertex budgets for all face-center samples.
- Expose `ntetrahedra`, retain the existing `ntriangles` 2D diagnostic,
  and clear conformity diagnostics on incomplete refinement batches.
- Add 47 native regression assertions, five SciPy quality comparisons,
  a 3D example, public-adapter smoke check and standalone validation.
- Preserve all prior 2D, adaptive, HGR and DIRECT behavior. This geometric
  construction does not supply SciPy SHGO's homology guarantees or constraints.

## 2026-11-10 — Conforming two-dimensional adaptive SHGO triangulation

- Add opt-in `conforming=.true.` to `shgo_core` and the public `shgo`
  adapter, requiring `adaptive=.true.` and a two-dimensional bounds array.
- Form positive-area triangle fans from each sampled leaf center; subdivide
  coarse faces at all cached hanging vertices. Verify the resulting mesh using
  exact dyadic area sums and interior/exterior edge multiplicities.
- Expose `ntriangles` and `conforming_valid` result diagnostics, and account
  for all added cell-center samples in the strict vertex/evaluation budgets.
- Add 47 native assertions, five independent SciPy quality comparisons,
  a public-adapter forwarding smoke test, and a complete demonstration.
- Preserve previous uniform/adaptive/HGR defaults and higher-dimensional
  behavior; topology conformity does not establish global convergence.

## 2026-11-10 — Optional SHGO sampled homology-growth stopping

- Add opt-in `hgr_stability_iters`, `hgr_min_iters` and `hgr_f_tol` to the
  SHGO core and public local-minimizer adapter, preserving all earlier defaults.
- Track estimated graph local-minimum births and stable-count streaks; stop
  only after complete rounds with no material best-objective improvement.
  Return a distinct `shgo_hgr_stable` status and useful diagnostics.
- Ensure deferred local polishing runs before accepting a tentative HGR stop,
  and never mislabel incomplete adaptive rounds as stable observations.
- Correct equal-value plateau detection on nonconforming adaptive leaf graphs
  via connected-component union/find, eliminating insertion-order artifacts.
- Add dedicated tests, adapter smoke tests, example, validation, and guide.
  No new dependencies or third-party implementation source were introduced.

## 2026-11-10 — Adaptive cell refinement for SHGO

- Add opt-in `adaptive` / `adaptive_fraction` controls at the end of both the
  standalone `shgo_core` and public `shgo` argument lists, preserving uniform
  defaults and existing callers.
- Refine selected dyadic leaf cells using a deterministic best/mean-corner
  score, exact integer-lattice vertex caching, Kuhn cell-edge candidate graphs,
  and atomic per-cell evaluation and memory limits.
- Track `nleaf` and `nrefined` alongside existing SHGO counters; retain
  existing local minimizers, target stopping and callback semantics.
- Add 46 native adaptive assertions and six independent SciPy quality checks,
  each exercised in checked and `-O3` builds; retain the prior 61 native and
  33 SciPy-uniform comparisons. See the new adaptive example and documentation.
- Adaptive triangulations have hanging nodes and do *not* establish SciPy's
  theoretical SHGO convergence guarantees. Conforming triangulation,
  constraints and homology-rank stopping remain future work.

## 2026-11-10 — Deterministic SHGO subset

- Add an independent bounded simplicial global optimizer (`shgo_core`) using
  nested Kuhn-Freudenthal grids and H_0 sublevel-set component births to find
  candidate local minima. Cache vertices across refinements and enforce a
  strict counted objective budget.
- Expose `shgo` with bounded Nelder-Mead/Powell local polishing, custom
  minimizer and callback hooks, and sorted deduplicated local minima.
- Add native regression tests, standalone checked/optimized SciPy comparisons,
  an example, API guide, and provenance; no new dependencies or copied code.
- This is not a full SHGO port: adaptive meshes, nonlinear constraints,
  homology-rank stopping and Sobol/Halton strategies remain future work.

- Add independent Tsallis generalized simulated annealing core and SciFort
  `dual_annealing` public wrapper for `scipy.optimize.dual_annealing`, with
  optional bounded L-BFGS-B/Powell local minimization, 48 SciPy benchmark
  comparisons in each of checked and optimized configurations, and full
  attribution and tests (2026-11-10).

## 2026-11-10 — Deterministic DIRECT global optimizer

- Add `src/optimize/direct.f90` implementing slope-based potentially optimal
  rectangle selection, dimension-wise trisecting, and deterministic bounded
  global search.
- Provide original DIRECT and locally biased size conventions, evaluation,
  iteration, target-objective, normalized volume and length stop criteria,
  and callback termination.
- Export `direct` and statuses through `scifort_optimize`. No new dependencies.
- Add native edge-case tests, SciPy comparison harness, standalone checked/
  optimized validation script, executable example, and API documentation.
- End-to-end FPM dependency integration is still pending because FPM and its
  Git dependencies were not available in the offline validation environment.


## 2026-11-10 — Basin-hopping global optimization

- Add original `scifort_optimize_basinhopping` stochastic basin-hopping core
  and `scifort_optimize_basinhopping_minimize` adapter to existing Nelder-Mead,
  Powell and BFGS local minimizers. Export public entry points through
  `scifort_optimize` without changing existing numerical APIs.
- Include adaptive random displacements, independent seeded RNG, Metropolis
  acceptance, user-supplied local minimizers, trial and accept callbacks,
  callback-based termination, and detailed evaluation/acceptance counters.
- Add checked/optimized native tests, independent SciPy benchmark comparisons,
  example, guide, validation report and reproducible standalone test script.
- No new FPM dependencies or copied third-party implementation code.

## 2026-11-10 — Dense linear programming (linprog)

- Add `scifort_optimize_linprog` with an original, self-contained two-phase
  tableau simplex for dense continuous LPs, inequalities/equalities, variable
  lower/upper bounds, free variables and redundant constraints.
- Export `linprog` and its result/status types through `scifort_optimize`.
- Add native tests, SciPy independent comparisons in checked and optimized
  builds, an executable example, API guide and reproducible test script.
- DualSimplex (MIT) was evaluated as an alternative backend but not copied,
  because its dense dual-only formulation needs extra integration work.
- No changes to existing numerical kernels or FPM dependencies.

## 2026-11-10 — Nonnegative and bounded linear least squares

- Add native `nnls` and `lsq_linear` (BVLS active-set algorithm) for dense
  real matrices with box constraints and rank-deficient input support.
- Reuse SciFort's Jacobi SVD for constrained free-variable subproblems.
- Add optional `nnls_slsqp` and `lsq_linear_slsqp` wrappers around the existing
  SLSQP 1.6.0 `bvls_module` (BSD-style upstream license). No third-party code
  copied and no new dependency.
- Export public APIs from `scifort_optimize` without changing existing calls.
- Add native and independent SciPy tests in checked/optimized builds; the real
  SLSQP backend still needs FPM validation.

## 2026-11-10 — SciPy-style differential evolution global optimization

- Add `scifort_optimize_differential_evolution` with all 12 standard mutation
  strategy names, binomial/exponential crossover, bound handling, seeded local
  RNG, random/LHS initializations, immediate/deferred population updates,
  convergence and iteration/evaluation accounting, optional callback and
  evaluation budget. Export through `scifort_optimize` without changing its
  previous public routines.
- Record algorithmic attribution to Ian Morrall's MIT-licensed Fortran
  DE/rand/1/bin example, and preserve his original MIT notice.
- Add native checked/optimized tests, a SciPy/NumPy numerical comparison
  harness, executable example, reproducible validation runner, and API guide.
- No new dependencies.

## 2026-11-10 — Standalone matrix exponential and dense-linalg API audit

- Add `scifort_linalg_expm` implementing real/complex matrix exponentials,
  with a single leading batch dimension, via fixed-degree [13/13] Pade
  scaling and squaring and the existing SciFort LU backend.
- Audit public SciPy `linalg` symbols against the source tree; note that
  `expm` was missing, despite existing `expm_frechet` support.
- Add reproducible native and independent SciPy tests, executable example,
  API guide, validation and coverage notes; both GNU Fortran build modes pass.
- All existing numerical modules and build dependencies are preserved.

## 2026-11-10 — Matrix-free randomized interpolative SVD

- Add `scifort_linalg_interpolative_svd_operator` with fixed-rank and
  tolerance-based real/complex SVD from forward/adjoint callbacks.
- Compose the existing matrix-free randomized ID with small-core QR/SVD;
  rescale skeleton factors to improve handling of extreme input magnitudes.
- Add independent NumPy/SciPy tests, native test, matrix-free example,
  API documentation and validation report. Both checked and optimized
  GNU Fortran builds pass all new checks and relevant regressions.
- Preserve all previous Fortran modules and runtime dependencies.

## 2026-11-10 — Matrix-free randomized interpolative decomposition

- Add `scifort_linalg_interpolative_operator`, with fixed-rank and
  tolerance-based real/complex operator interfaces and optional exact
  skeleton-column outputs. No full dense input matrix is necessary.
- Use randomized adjoint sketching, pivoted QR and basis-space least-squares
  refitting; sample adaptive Frobenius residuals with independent probes.
- Add a true implicit low-rank example and native/SciPy comparisons.
- Preserve all previous source modules, interfaces and dependencies.

## 2026-11-10 — Randomized interpolative decomposition

- Add `scifort_linalg_interpolative_random` with real/complex fixed-rank and
  adaptive randomized ID, using local reproducible Gaussian row sketches.
- Select original matrix columns by pivoted QR of a sketch; fit coefficients
  against the original scaled data using the existing SVD least-squares solver.
- Verify adaptive Frobenius-norm error with sketch expansion and a full-rank
  fallback, including empty, zero, dependent, and extreme-scale matrices.
- Add native and independent NumPy/SciPy tests, executable example, and API
  documentation; preserve previous Fortran APIs and runtime dependencies.

## 2026-11-10 — Randomized interpolative SVD

- Add `scifort_linalg_interpolative_svd_random`, a new module for fixed-rank
  and adaptive real/complex randomized ID-based SVD.
- Use reproducible Gaussian sketches, original-matrix least-squares
  interpolation and SciFort's existing QR/SVD implementations.
- Verify adaptive relative Frobenius error against the original input.
- Add 50,408 native checks and 730 independent reference cases per build;
  both checked and optimized GNU Fortran builds passed, along with the
  existing interpolative, deterministic ID-SVD, and SVD tests.
- Include a reference driver, runnable example, API guide and validation.
- Preserve all existing Fortran modules and runtime dependencies.

## 2026-11-09 — Deterministic SVD via interpolative decomposition

- Add `scifort_linalg_interpolative_svd` with fixed-rank and adaptive relative
  Frobenius tolerance real/complex interfaces; produce U, s, V.
- Reuse SciFort pivoted QR, interpolative decomposition, and Jacobi SVD.
- Test checked/optimized builds, independent SciPy deterministic reference,
  native edge cases, QR/SVD regressions, and executable example.
- Document the distinction from SciPy's default randomized implementation.
- No existing source modules or dependencies were changed.

## 2026-11-09 — Randomized numerical rank estimation

- Add `scifort_linalg_estimate_rank` for real/complex dense and forward-only
  matrix-free operators, using adaptive random range probes and two-pass
  modified Gram–Schmidt.
- Support optional `int64` seeds, oversampling, rectangular and zero-sized
  matrices, and explicit shape/nonfinite/callback failure statuses.
- Add checked/optimized native tests, NumPy-SVD rank verification,
  independent SciPy comparisons, executable example, and API documentation.
- Document different stochastic output semantics from SciPy's rank estimator.
- Preserve all existing Fortran sources and build dependencies.

## 2026-11-09 — Randomized spectral-norm estimates

- Implement `scipy.linalg.interpolative.estimate_spectral_norm` and
  `estimate_spectral_norm_diff` for real and complex dense matrices.
- Add matrix-free callbacks for both real/complex operator norms and matrix
  differences, with conjugate-adjoint products and no dense allocation.
- Add seeded and explicit-start reproducibility, scaled dense arithmetic,
  explicit error status codes and zero-size handling.
- Validate against independent power iterations, SciPy randomized estimates
  and exact NumPy SVD norms in checked and optimized GNU Fortran builds.
- Preserve all existing SciFort source modules and dependency declarations.

## 2026-11-09 — Clarkson-Woodruff / CountSketch random projection

- Add `scifort_linalg_clarkson_woodruff` implementing real/complex dense,
  one-leading-batch, and one-based CSR sparse CountSketch accumulation.
- Allow reuse of a signed bucket mapping across input matrices, with
  independent int64 seeded or intrinsic Fortran random plan generation.
- Add native tests, independently seeded SciPy comparisons (990 cases in each
  checked/optimized build), an executable example, and API/validation docs.
- Keep all existing Fortran modules and dependencies unchanged.

## 2026-11-09 — Deterministic interpolative decomposition

- Add `scifort_linalg_interpolative` with real/complex fixed-rank and
  relative-Frobenius-tolerance `interp_decomp` interfaces and one-based
  permutation indices.
- Add skeleton and interpolation matrix reconstruction, direct matrix
  reconstruction, and small-core ID-to-SVD conversion.
- Reuse existing pivoted QR and Jacobi SVD, with an SVD fallback for
  numerical rank deficiency; no new runtime dependencies.
- Add 14,944 native checks and 2,016 independent SciPy reference cases
  per checked/optimized build, an example, and API/validation documentation.

## 2026-11-09 — Hermitian pseudoinverse (`pinvh`)

- Add `scifort_linalg_pinvh` for symmetric/Hermitian Moore--Penrose
  pseudoinverses, using positive and negative eigenvalues with configurable
  `atol`/`rtol`, numerical rank, upper/lower triangle selection, and a
  single batch dimension.
- Add native and independent SciPy comparisons, an executable example,
  API documentation, and reproducible validation instructions.
- No existing Fortran modules or dependencies were changed.

## 2026-11-09 — Real invariant eigenbasis (`cdf2rdf`)

- Add `scifort_linalg_cdf2rdf`, converting conjugate eigenpairs of real
  matrices to real block eigenvalues and a real invariant vector basis.
- Support real/complex eigenvalues, real/complex vectors, one batch
  dimension, nonadjacent pair members, and eigenvector phase differences.
- Add a fallback for tiny spurious imaginary components returned by
  SciFort `eig`, plus native and independent SciPy comparison tests.
- Preserve all existing source APIs and build dependencies.

## 2026-11-09 — Rectangular LU, determinant, and inverse

- Add `scifort_linalg_lu` for singular-safe rectangular, real/complex `lu`,
  `lu_indices`, and `lu_permuted` with one-based permutation conventions.
- Add `scifort_linalg_det_inv` for real/complex `det` and `inv`, reusing
  SciFort's existing dense LU factoring/solving backend for inverses.
- Add native tests, SciPy comparison driver and checker, executable example,
  API/validation documentation. No new dependencies or existing source edits.

## 2026-11-09 — Real/complex RQ factorization

- Add `scifort_linalg_rq` (`rq` and `rq_r`) for full, economic, R-only,
  empty core dimensions and one leading batch dimension.
- Compute the factorization from reverse-conjugate-transpose Householder
  QR, reusing existing internal routines; R-only skips constructing Q.
- Add native tests, an independent SciPy reference driver/checker,
  executable example and API/validation documentation. No new dependencies.

## 2026-11-09 — Householder QR multiplication

- Add `scifort_linalg_qr_multiply`, providing `qr_multiply` and reusable
  `qr_apply` for real, complex, vector, and matrix multipliers; no explicit
  Q construction or new numerical dependencies.
- Add native tests, SciPy phase-aligned reference comparisons, example and docs.

## 2026-11-09 — QR rank-k updates, insertion and deletion

- Add `scifort_linalg_qr_updates` with real/complex `qr_update`,
  `qr_insert`, and `qr_delete`, for complete and economy QR factors.
- Reuse existing QR factors for updates and column edits via reduced-core
  factorizations; use complete refactorization for row edits and dependent
  column insertions, preserving explicit status reporting.
- Add native and independent SciPy comparison tests, executable example,
  API documentation, and validation report. No new dependencies.

## 2026-11-09 — Complete cosine–sine decomposition

- Implement `scipy.linalg.cossin` in `scifort_linalg_cossin`, with real and
  complex full-factor and separate-factor APIs, all valid partitions, and
  both `swap_sign` layouts. Uses internal Jacobi SVD without LAPACK linkage.
- Add native tests (7,800 partition/sign cases in checked and optimized
  builds), independent SciPy comparisons (3,070 per build), example, and docs.
- Reject nonunitary/nonfinite input, document differences in phase gauge,
  batch and subblock input support, and optional-factor computation.

## 2026-11-09 — Matrix inspection and SVD diagonal construction

- Add `scifort_linalg_structure.bandwidth`, `issymmetric`, `ishermitian`.
  Support int32, int64, real64 and complex64, including one leading batch
  dimension and SciPy 1.17's exact integer and rectangular edge cases.
- Add `scifort_linalg_diagsvd.diagsvd`, constructing real, integer, and
  complex rectangular SVD diagonals, with one-batch overloads.
- Add native tests, reproducible seeded SciPy comparison script and driver,
  executable example, API documentation and validation report. All 7,248
  independent SciPy comparisons passed across checked and optimized builds.
  No existing source module was modified; no new dependencies were added.

## 2026-11-09 — Block-diagonal matrices and Khatri–Rao products

- Add `scifort_linalg_block_diag.as_block` and `block_diag`, supporting
  heterogeneous scalar, vector, matrix, and one-batch matrix blocks with
  integer32/integer64/real64/complex64 output kinds. Match SciPy's
  empty-vector row convention and preserve exact `int64` values.
- Add `scifort_linalg_khatri_rao.khatri_rao`, supporting all 16 pairs of
  input numeric kinds and all four 2-D/3-D input rank combinations,
  including singleton batch broadcasting and zero-dimensional core arrays.
- Add native tests, independent seeded SciPy comparisons, reference driver,
  executable example, API docs, provenance, and validation report. Existing
  Fortran source remains unchanged; no additional dependencies are needed.

## 2026-11-09 — Hadamard, Leslie and DFT matrix constructors

- Add `scifort_linalg_hadamard.hadamard` with Sylvester ordering and
  integer32/integer64/real64/complex64 output overloads.
- Add `scifort_linalg_leslie.leslie` with integer, real/complex/mixed
  vector inputs and one leading broadcastable batch dimension for real/complex.
- Add `scifort_linalg_dft.dft` with unnormalized, 1/sqrt(n) and 1/n
  Fourier matrices. Preserve existing source and runtime dependencies.
- Add native tests, independent SciPy comparison driver/harness,
  runnable example, API documentation and provenance. Both checked and
  optimized builds passed 1473 comparisons each and prior Pascal,
  Helmert/Hilbert and companion/Fiedler regression tests.

## 2026-11-09 — Pascal and inverse Pascal matrix constructors

- Add `scifort_linalg_pascal.pascal` and `invpascal`, with symmetric, lower,
  and upper forms, real64 and overflow-checked exact int64 overloads.
- Add SciFort exact-decimal-string constructors using a portable base-10^9
  integer representation for coefficients exceeding machine integer ranges.
- Validate 948 independent SciPy comparisons across checked and optimized
  builds, exact inverse identities, overflow boundaries, and prior matrix
  constructor regression tests. Add tests, runnable example and API guide.
  No existing Fortran source modules or external dependencies were changed.

## 2026-11-09 — Helmert, Hilbert and exact inverse Hilbert matrices

- Add `scifort_linalg_helmert.helmert` for full and contrast-only
  orthonormal matrices, including the single-observation edge case.
- Add `scifort_linalg_hilbert.hilbert` and `invhilbert`, with real64
  closed-form inverse coefficients and exact int64 results through order 14.
- Add `invhilbert_exact_decimal` as a SciFort arbitrary-precision option,
  using internally generated decimal strings to avoid third-party dependencies.
- Add checked and optimized native tests, 956 independent SciPy comparison
  cases, reference driver, executable example, documentation and provenance.
  No existing source files are modified, and no new dependencies are needed.

## 2026-11-09 — Companion, Fiedler companion and Fiedler matrices

- Add `scifort_linalg_companion.companion` and `fiedler_companion`,
  matching the polynomial matrix conventions of `scipy.linalg`.
- Add `scifort_linalg_fiedler.fiedler`, a symmetric absolute-distance
  matrix constructor with real, complex, and integer inputs.
- Implement one batch dimension, explicit validation statuses and
  zero-filled outputs on errors; retain existing SciFort dependencies.
- Include native checked/optimized tests, 3,420 independent SciPy matrix
  comparisons, 116 polynomial root checks, runnable example, docs and provenance.

## 2026-11-09 — Convolution matrix construction and FFT-based products

- Add `scifort_linalg_convolution_matrix.convolution_matrix` matching
  `scipy.linalg.convolution_matrix` for `full`, `same`, and `valid` modes,
  including short inputs, integer/real/complex coefficients, and one batch dimension.
- Add `scifort_linalg_matmul_convolution_matrix.matmul_convolution_matrix`
  for structured FFT products, with real/complex/mixed vector and matrix RHS;
  reuse the existing Toeplitz FFT implementation and avoid materializing dense matrices.
- Add native checked/optimized tests, independent SciPy comparison harness,
  example, API documentation and validation report. No new dependencies.

## 2026-11-09 — Hankel construction and FFT-based multiplication

- Add `scifort_linalg_hankel.hankel`, corresponding to `scipy.linalg.hankel`,
  supporting rectangular and complex matrices, mixed real/complex inputs,
  zero-filled optional last rows, and zero-size dimensions.
- Add SciFort extension `scifort_linalg_matmul_hankel.matmul_hankel`, with
  FFT-based products for real/complex vectors and multiple RHS columns.
  Reuse `matmul_toeplitz` after reversing the implied Hankel columns;
  no new FFT implementation or third-party dependencies.
- Add a native regression test, independent SciPy/NumPy reference harness,
  example, user documentation, and a validation report. All 3,600 numerical
  comparisons and 452 constructor comparisons passed in checked/optimized
  GNU Fortran builds, alongside existing structured-matrix regressions.

## 2026-11-09 — FFT-based Toeplitz multiplication

- Add `scifort_linalg_matmul_toeplitz.matmul_toeplitz`, using the existing
  FFT plan with circulant embedding to multiply square or rectangular
  Toeplitz matrices in O((m+n) log(m+n)) time per right-hand side.
- Support real, complex and mixed-type vectors and matrices, an optional
  explicitly supplied first row, conjugate default and ignored `r(1)`.
- Add checked/optimized native tests, 3,600 independent SciPy/NumPy
  reference comparisons, runnable example, API documentation and provenance.
- No new dependencies, third-party source code or license obligations.

## 2026-11-09 — FFT-based circulant construction and solves

- Add `scifort_linalg_circulant.circulant` and `solve_circulant`, supporting
  real, complex, mixed-type, vector and multiple right-hand sides.
- Reuse the existing reusable FFT plan for O(n log n) solves; implement
  SciPy-compatible singularity thresholding and spectral minimum-norm
  least squares (`singular='lstsq'`), plus one-based rank-two axis options.
- Add checked and optimized native tests, 2,400 independent SciPy 1.17.0
  reference cases, example, documentation, and provenance. No new dependencies.

## 2026-11-09 — Toeplitz construction and quadratic-time solves

- Add `scifort_linalg_toeplitz.toeplitz` for real/complex rectangular
  Toeplitz matrices and `solve_toeplitz` for real, complex and mixed-type
  vector/multiple RHS, with SciPy's optional conjugate row convention.
- Use O(n^2) general Toeplitz bordering recurrence and existing pivoted LU
  fallback for zero/small leading Schur pivots; normalize coefficients to
  handle uniformly extreme scales. No new runtime dependencies.
- Add native checked/optimized tests and 3,200 independent SciPy/NumPy
  reference comparisons, including 128 fallback cases, with executable
  example, API documentation, and provenance.

## 2026-11-09 — LDL-based symmetric/Hermitian linear solves

- Add `scifort_linalg_indefinite_solve` with `solve_indefinite` for the
  `sym`/`her` specializations of SciPy's dense `solve` and `ldl_solve` to
  reuse pivoted LDL factors, both for real/complex vector and matrix RHS.
- Handle 1x1 and 2x2 Bunch–Kaufman blocks, one-based permutations, upper
  and lower factors, and complex symmetric versus Hermitian transposes.
- Check invalid factors, singular blocks, nonfinite RHS and zero-size cases;
  add native, independent SciPy, and NumPy-reference tests (4,500 comparisons).
- Add a runnable example, API and validation guides. No new dependencies.

## 2026-11-09 -- Symmetric/Hermitian indefinite LDL factorization

- Add `scifort_linalg_ldl.ldl` for real symmetric, complex Hermitian and
  complex symmetric matrices, with Bunch-Kaufman 1x1/2x2 pivot selection.
- Follow SciPy's permuted outer-factor and block-diagonal output convention,
  accepting either selected input triangle and one-based permutations.
- Add native randomized tests, 2,400 checked/optimized SciPy reference cases,
  Fortran example, API documentation and provenance. No new dependencies.

## 2026-11-09 — Ordered complex generalized Schur decomposition (`ordqz`)

- Add `scifort_linalg_ordqz.ordqz` with four real/complex input overloads,
  homogeneous `alpha/beta` outputs, SciPy's four sorting rules, optional
  custom selector and `sdim`, and existing QZ method selection.
- Add `ordqz_factors` to reorder precomputed triangular generalized Schur
  forms using adjacent two-sided unitary swaps, without rerunning QZ.
- Normalize the separate coefficient factors and verify reconstruction;
  classify roundoff-sized denominators as numerical infinite eigenvalues.
- Add native tests, independent SciPy comparisons, example, API guide and
  provenance. No third-party code or new dependencies introduced.

## 2026-11-08 — Generalized nonsymmetric eig/eigvals

- Add `scifort_linalg_generalized_eigen.eig/eigvals` for real/complex and
  mixed-kind regular matrix pairs, with optional left/right eigenvectors.
- Return homogeneous `(alpha,beta)` eigenvalues, including infinite roots,
  and avoid overflow in finite eigenvalue ratio conversion.
- Reuse the two-sided QZ implementation with an exact triangular-pair
  shortcut; preserve ordinary eig/eigvals interfaces and dependency set.
- Add native tests, SciPy-projective-spectrum comparison harness, example,
  API documentation and provenance record.

## 2026-11-08 — Inverse-free two-sided complex QZ iteration

- Add `qz(...,method='iterative')` with Hessenberg–triangular reduction,
  complex two-sided Givens bulge chasing, and trailing 2×2 generalized shifts.
- Keep `method='inverse'` for the previous construction and add conservative
  iterative-first `method='auto'` dispatch with failure fallback.
- Add native tests of real/complex pencils, conjugate pairs, extreme scaling,
  singular denominators, and method error handling; extend the independent
  SciPy reference harness to compare all three methods.
- Document limitations relative to robust LAPACK QZ; no new dependencies
  or third-party source code.

## 2026-11-08 — Generalized complex Schur factors for regular pencils

- Add `scifort_linalg_qz.qz` with real, complex, and mixed real/complex
  inputs, complex unitary output bases, and upper-triangular generalized
  Schur factors, plus optional alpha/beta diagonal outputs.
- Build factors from ordinary Schur, LU and QR; handle invertible-B,
  singular-B/invertible-A, and many regular both-singular cases through a
  nonsingular shifted pencil. Preserve explicit failure for unsupported
  singular and poorly conditioned cases.
- Add native tests, SciPy comparison driver, executable example, API guide,
  numerical validation, and provenance. No new library dependencies.
- Note explicitly that this is not a full LAPACK-style iterative QZ solver.

## 2026-11-08 — Conversion of real Schur factors to complex (`rsf2csf`)

- Add `scifort_linalg_schur.rsf2csf` for precomputed real Schur factors and
  real orthogonal bases; construct complex triangular factors with local
  two-dimensional unitary similarities rather than rerunning Schur QR.
- Support empty inputs and real/conjugate Schur blocks, with internal input
  scaling for very small or large matrix entries and explicit error statuses.
- Add native tests, independent SciPy reference checks (including SciFort's
  own real Schur input), a runnable example, documentation and provenance.
- The conversion is not subject to the existing Schur engine's order-40 bound.
  No additional numerical library or third-party source code is required.

## 2026-11-08 — Real quasi-triangular Schur decomposition

- Extend the public `schur` generic with real `T`/`Z` output for real inputs;
  retain the established complex interfaces unchanged.
- Group complex-conjugate eigenvalues, reorder them with unitary Schur swaps,
  and construct real orthogonal invariant-subspace bases for 1x1/2x2 blocks.
- Support the four built-in sorting regions and a custom selector that moves
  whole conjugate pairs together, with `sdim` counting eigenvalues.
- Add native tests, randomized checked/-O3 SciPy comparisons, a revised
  example, API/validation notes, and provenance. No new dependencies.

## 2026-11-08 — Public complex Schur decomposition and reordering

- Add `scifort_linalg_schur.schur` for real and complex coefficient matrices,
  returning a complex upper-triangular Schur factor and a unitary basis.
- Support SciPy's four built-in eigenvalue sort regions and optional custom
  predicate callbacks, with the reordered invariant-subspace dimension.
- Add native tests, two independent checked/optimized SciPy comparison seeds,
  a standalone example, and API/validation/provenance documentation.
- Real quasi-upper-triangular Schur output remains future work. No external
  numerical dependency or imported third-party code was added.

## 2026-11-08 — General nonsymmetric dense eig/eigvals

- Add real/complex `scifort_linalg_general_eigen.eig` and `eigvals` for
  ordinary nonsymmetric matrices up to order 40.
- Recover normalized left and right eigenvectors with triangular back
  substitution after unitary Schur reduction; stabilize triangular matrices
  with a direct similarity path, avoiding unnecessary QR perturbations.
- Add native checked/optimized tests and independent SciPy comparisons,
  a runnable example, API documentation and provenance notes. No external
  numerical dependencies or third-party implementation code added.

## 2026-11-08 — Real banded Givens kernel and conservative auto dispatch

- Replace the complex conversion in `band_to_tridiagonal` for real symmetric
  packed bands with a dedicated real Givens similarity and bulge chase.
- Use `method='bulge'` automatically for full-spectrum *eigenvalues-only*
  problems when `n >= 192` and `kd <= min(8,max(1,n/16))`, falling back to
  the established dense path on nonconvergence. Eigenvectors and selected
  eigenvalue requests retain their previous automatic paths.
- Add native auto-dispatch tests, independent SciPy comparison cases and a
  deterministic microbenchmark. Preserve the no-new-dependencies policy.

## 2026-11-08 — Banded bulge-chasing tridiagonalization

- Add `scifort_linalg_banded_tridiagonal.band_to_tridiagonal`, a Givens
  bulge-chasing reduction for real symmetric and complex Hermitian bands.
- Keep the matrix in a band plus one transient bulge diagonal; full
  eigenvalues-only `method='bulge'` avoids dense input expansion.
- Add opt-in `method='bulge'` to `eig_banded` and `eigvals_banded`;
  preserve defaults and support existing subset/capacity conventions.
- Add native similarity tests, independent SciPy comparisons, example,
  validation notes, and an algorithm provenance entry.

## 2026-11-08 — Opt-in full-spectrum banded inverse iteration

- Allow `eig_banded(...,method='inverse')` without a subset, computing every
  eigenpair via banded Sturm bisection and shifted LU inverse iteration.
- Implement real-only eigenvector work arrays for real symmetric band inputs,
  avoiding the previous complex n-by-k temporary and complex LU overhead.
- Preserve dense full-spectrum mode as the default; explicitly report
  nonconvergence and insufficient caller output capacity.
- Add native tests, real/complex SciPy comparisons, a runnable example,
  API notes and provenance. A band-to-tridiagonal reduction remains future work.

## 2026-11-08 — Banded selected eigenvectors: real-LU specialization

- Extended `eig_banded` with a selected-eigenvector path based on Sturm bisection,
  pivoted band LU inverse iteration, and two-pass reorthogonalization.
- Real symmetric input uses real LU factorization and solves, avoiding a complex
  LU allocation; complex Hermitian input uses the complex solver.
- Preserved dense full-spectrum behavior and the existing positional API.
- Added native and independent SciPy tests, including large narrow bands and
  real-vs-complex consistency checks. The full FPM suite is still untested.

## 2026-11-08 — Selected symmetric/Hermitian band eigenvectors

- Add opt-in low-memory selected `eig_banded` eigenpairs for real symmetric
  and complex Hermitian input using pivoted band LU inverse iteration.
- Expose one-based inclusive index and half-open value selection, with
  explicit `num_selected`, method and tolerance controls.
- Add selected `eigvals_banded` Sturm bisection for O(n*bandwidth) workspace.
- Preserve original full-spectrum dense path, with small-size auto fallback.
- Include tests, SciPy comparisons, example, documentation, and provenance.

## 2026-11-08 — Selected tridiagonal eigenvectors by inverse iteration

- Extend selected `eigh_tridiagonal` eigenpairs with Sturm bisection followed by
  partial-pivoted tridiagonal inverse iteration and two-pass reorthogonalization.
- Reduce selected-eigenvector storage from a full n-by-n orthogonal matrix to
  the n-by-k requested output plus O(n) workspace; retain full-QL fallback
  for smaller difficult problems and explicit `eigenvector_method` controls.
- Add a 10000-by-10000 narrow-spectrum native test, independent SciPy
  eigenpair comparisons, error-case checks, documentation and a runnable example.
- Preserve existing positional call patterns; introduce no new dependencies.

## 2026-11-08 — Selected tridiagonal eigenvalues via Sturm bisection

- Extend `eigh_tridiagonal` and `eigvalsh_tridiagonal` with optional inclusive
  one-based `subset_by_index`, half-open `(low,high]` `subset_by_value`,
  `num_selected`, and absolute bisection `tol`.
- Use scaled Sturm-sequence counts and Gershgorin brackets to compute only the
  requested eigenvalues in O(n) memory without constructing eigenvectors.
- Preserve the original full-spectrum QL algorithm and all positional calls;
  selected eigenvectors currently use full-spectrum QL followed by extraction.
- Add checked and optimized native tests, an order-10000 Toeplitz stress case,
  independent SciPy comparisons, example, API documentation, and provenance.
- No third-party implementation copied or new numerical dependency introduced.

## 2026-11-08 — Half-open eigenvalue selection by value

- Extend real and complex `eigh`/`eigvalsh` with optional `subset_by_value`
  and `num_selected` for standard and generalized definite problems.
- Select eigenvalues in `(low,high]`, including repeated values, empty
  selections, and IEEE infinity bounds; reject invalid and conflicting ranges.
- Preserve all existing positional arguments and index-subset behavior.
  Selection by value uses caller-owned output capacity buffers.
- Add native tests, independent SciPy comparison harness, example, API notes,
  and provenance; no new numerical dependency or borrowed implementation.

## 2026-11-08 — Generalized symmetric/Hermitian definite eigenproblems

- Extend existing real/complex `eigh` and `eigvalsh` with optional `b` and
  `problem_type=1|2|3`, retaining the old positional interface and default.
- Reduce generalized pencils with SciFort's Cholesky/triangular solves,
  reuse the Householder/QL solver, and return correctly normalized vectors.
- Validate against SciPy including both input triangles, rank-deficient A,
  selected eigenpairs, extreme scaling, and invalid/non-positive-definite B.
- Add native tests, SciPy reference harness, example, API/validation notes,
  and provenance. No new numerical dependencies or borrowed code.

## 2026-11-07 — Dense symmetric/Hermitian eigenproblems

- Add real and complex `scifort_linalg_dense_eigen.eigh` and `eigvalsh`,
  including upper/lower triangle handling and inclusive one-based
  `subset_by_index` selection.
- Reuse existing SciFort Householder reduction and implicit-shift QL
  tridiagonal eigensolver, without new numerical dependencies.
- Include native tests, 1,980 independent SciPy comparisons,
  420 neighboring regression comparisons, example, API and provenance notes.
- Document O(n^3) dense complexity and absence of generalized problems.

## 2026-11-07 — Symmetric tridiagonal and Hermitian band eigensolvers

Added `scifort_linalg_banded_eigen` with `eigh_tridiagonal`,
`eigvalsh_tridiagonal`, `eig_banded`, and `eigvals_banded` for real
symmetric/complex Hermitian banded inputs. The tridiagonal path uses
scaled implicit-shift QL; packed bands are expanded to dense and
reduced by SciFort Householder Hessenberg. Native and SciPy comparison
tests, a runnable example, API documentation and provenance are included.
General banded eigensolving is currently O(n^3) time/O(n^2) storage.

## 2026-11-07 — Matrix balancing

- Add `scifort_linalg_balance.matrix_balance` (real and complex) with
  optional isolated-row/column permutation and binary-power diagonal scaling.
- Provide optional full transformation or separate scale/permutation outputs.
- Include native Fortran tests, an independent SciPy comparison harness,
  example, API documentation, numerical-validation notes and provenance.
- Preserve existing APIs and introduce no additional numerical dependency.

## 2026-11-07 — Hessenberg reduction

- Add `scifort_linalg_hessenberg.hessenberg` for real and complex matrices,
  returning the upper Hessenberg form and optional orthogonal/unitary basis.
- Add native tests and randomized SciPy reference comparisons, a Fortran
  example, error-handling/validation documentation and algorithm provenance.
- Preserve all previous public interfaces and introduce no new dependencies.

## 2026-11-07 — Schur matrix-equation methods

- Extend existing real/complex `solve_sylvester`, `solve_continuous_lyapunov`,
  and `solve_discrete_lyapunov` with `method='auto'|'direct'|'schur'`.
- Implement Schur-transformed triangular recurrences, allowing order-40
  inputs without constructing dense Kronecker coefficient matrices.
- Retain direct elimination as an explicit option and for small auto-selected
  problems; add singularity, shape, and nonfinite checks.
- Add native and independent SciPy tests, example, updated documentation and
  numerical provenance; no new dependencies or borrowed source code.

## 2026-11-07 — Dense orthogonal and unitary Procrustes alignment

- Add `scifort_linalg_procrustes.orthogonal_procrustes`, a generic for
  all real/complex pairings, solving `min ||A R - B||_F` with orthogonal
  or unitary `R`, and returning the sum of singular values of `A^H B`.
- Reuse the independent SciFort one-sided Jacobi SVD; normalize before
  cross-product multiplication and restore the scale with binary exponents.
- Support empty arrays and rank deficiency, with explicit invalid-input and
  numerical-failure statuses. Add native tests, streaming SciPy comparisons,
  a Fortran example, and detailed API/validation documentation.
- No third-party code or runtime dependency added.

## 2026-11-07 — Efficient singular values and spectral diagnostics

- Replace `svdvals` full-SVD allocation with independent values-only Jacobi
  sweeps; avoid singular-vector construction and use a native real arithmetic
  path for real inputs. Preserve existing signatures and rank-floor semantics.
- Add real/complex `cond2` spectral condition number and `matrix_rank` with
  optional absolute tolerance, explicitly labeled SciFort extensions.
- Add failure-path native tests, 1,260 checked/optimized SciPy/NumPy comparisons,
  reproducible benchmark, API/provenance/validation documentation and example.
- No new third-party code or numerical dependencies.

## 2026-11-07 - Full dense SVD, singular values, and rectangular diagonal matrices

- Extended `scifort_linalg_svd.svd` with optional `full_matrices=.true.`
  for full square real/complex singular vector bases, while preserving the
  economy-size default and existing positional argument order.
- Added `svdvals` (real/complex) and `diagsvd`, using existing Jacobi SVD
  infrastructure without external BLAS/LAPACK dependencies.
- Added native tests, full-SVD SciPy reference driver and comparison harness,
  runnable example, and validation documentation.

## 2026-11-07 — SVD-based range, null space and principal angles

- Add independent MIT-licensed `scifort_linalg_subspaces`, providing real and
  complex `orth` and `null_space` and all four real/complex combinations for
  `subspace_angles`.
- Complete economy right-singular bases for wide-matrix null spaces without
  forming a dense pseudoinverse. Use two-pass orthogonalization, canonical
  pivoting, and an updated residual projector diagonal.
- Use both sine and cosine formulations for stable small/large principal angles.
- Add native tests, SciPy projector/angle comparisons, a runnable example,
  API and validation documentation. No new numerical dependencies or licenses.


## 2026-11-07 — Dense economy SVD, pseudoinverse, and rank-deficient least squares

- Add original real/complex one-sided Jacobi economy `svd`, Moore-Penrose
  `pinv`, and reusable `svd_lstsq` minimum-norm solver for rectangular,
  rank-deficient matrices, with singular-value rank cutoff and status checks.
- Extend the existing `lstsq` generic with optional `method='auto'|'qr'|'svd'`;
  automatic mode preserves the previous full-rank QR path and falls back to
  SVD when pivoted QR detects rank deficiency. Keep `qr_solve` unchanged.
- Add native real/complex tests, streaming SciPy comparison driver,
  runnable example, API/validation documentation, and provenance.
- Preserve standalone Fortran builds without an additional numerical library.

## 2026-11-07 — Real/complex dense QR and full-rank least squares

- Add `scifort_linalg_qr` with reduced/complete Householder `qr`, compact
  `qr_factor`/`qr_unpack`, optional norm-based column pivoting, and one-based
  final permutation indices.
- Add reusable `qr_solve` for full-column-rank overdetermined systems and
  `lstsq` for full-rank tall and wide systems (minimum-norm solutions when wide).
- Add native shape, nonfinite, pivot and rank-deficiency tests; independent
  SciPy comparison driver and reproducible checked/optimized runner.
- Document unsupported rank-deficient minimum-norm solutions and the
  unblocked `O(m*n^2)` performance limitation.
- New source is MIT-licensed, without new numerical dependencies.

## 2026-11-07 — Dense Cholesky decomposition and reusable solves

- Add real/complex `cholesky`, `cho_factor`, `cho_solve` to `scifort_linalg_dense`,
  with upper/lower triangular storage, vector/matrix right-hand sides, and
  positive-definiteness and nonfinite-input status checks.
- Factor only the selected triangle and return zero in the unused output
  triangle; reuse `solve_triangular` for forward/back substitution.
- Remove unreachable empty loops from the older dense LU factorization.
- Add native failure tests, independent SciPy comparisons in checked and `-O3`
  builds, a runnable example, and API/validation/provenance documentation.
- No new external or third-party source dependencies.

## 2026-11-07 — Dense triangular systems and reusable LU factors

- Add `scifort_linalg_dense` with real/complex `solve_triangular`, `lu_factor`,
  and `lu_solve` generics; vector/matrix right-hand sides and `'N'/'T'/'C'` modes.
- Support upper/lower triangle selection, implicit unit diagonal, and reusable
  one-based pivoted LU factors without external BLAS/LAPACK dependencies.
- Add native Fortran tests, runnable example, independent SciPy validation,
  API documentation, and numerical provenance notes.

## 2026-11-07 — Reusable packed-band LU factors and transpose solves

- Add `lu_factor_banded` and `lu_solve_banded` with real/complex binary64,
  1-based local partial pivot indices, and packed extended LU storage.
- Enable factor reuse for different vector/matrix RHSs and for normal,
  transpose and conjugate-transpose solves (`trans='N'/'T'/'C'`).
- Refactor the previous one-shot `solve_banded` to use the same LU backend,
  eliminating duplicate matrix-elimination implementations.
- Add native tests with forced pivots and an order-50,000 system,
  independent SciPy comparisons in checked and `-O3` modes, and new
  example, API, provenance and numerical-validation documentation.
- Original MIT implementation, with no new runtime dependencies.

## 2026-11-07 — General and positive-definite packed-band solvers

- Add original MIT-licensed `scifort_linalg_banded` with `solve_banded`,
  `solveh_banded`, `cholesky_banded`, and `cho_solve_banded`; real and complex
  binary64, vector and matrix RHS, SciPy-compatible packed-band storage.
- Keep banded LU fill and banded Cholesky work linear in matrix order for
  fixed bandwidth; provide explicit status codes and preserve caller inputs.
- Add native tests, sample Fortran driver, independent checked/optimized
  SciPy comparisons, API/validation documentation and provenance notes.
- No new BLAS/LAPACK dependencies or bundled third-party source code.

## 2026-11-07 — Adaptive Taylor evaluation for dense `funm` clusters

- Add opt-in `taylor_max_order` (1--128) to the real/complex `funm` and
  `funm_real` interfaces for coupled Schur eigenvalue clusters.
- Sum analytic higher derivatives until three consecutive tail terms are
  negligible, or return nonconvergence; preserve the default degree-`n-1`
  callback contract for existing callers.
- Check the scalar spectral mapping on diagonal entries and distinguish
  mathematical nilpotence from physical-scale underflow.
- Add native failure tests, a new near-cluster example, independent SciPy
  comparisons including matrices up to order 40, and validation documentation.
- Original MIT implementation; no new external dependencies.

## 2026-11-07 — Block Schur--Parlett `funm` and semisimple spectra

- Replace scalar cross-block recurrence with a block Sylvester solve on
  globally reordered spectral clusters (including unequal block orders).
- Reorder spectral clusters even without a derivative callback; allow
  diagonalizable repeated eigenvalues with negligible within-block coupling.
- Continue to reject non-diagonalizable repeated spectra without derivatives.
- Add native mixed-block tests and SciPy tests for real/complex, upper
  triangular and dense similarity-transformed repeated spectra to order 40.
- Preserve the existing public API and dependency-free MIT implementation.

## 2026-11-07 — Global Schur clustering for dense `funm`

- Extend `scifort_linalg_funm` to cluster repeated or close Schur eigenvalues
  globally rather than only when they appear in adjacent diagonal positions.
- Stably reorder clusters through adjacent unitary transformations of both
  the triangular Schur factor and its unitary basis, then evaluate Taylor
  diagonal blocks using the existing derivative callback.
- Add native real/complex regressions and independent specialized SciPy
  comparisons for noncontiguous repeated spectra through order 40.
- Keep existing API/status values, MIT provenance, and runtime dependencies.

## 2026-11-07 — General analytic dense matrix functions (`funm`)

- Add `scifort_linalg_funm` with real/complex `funm`, `funm_real`, scalar
  callbacks, and optional k-th derivative callbacks for repeated spectra.
- Reuse the existing MIT-licensed complex Schur routine and implement the
  Parlett triangular recurrence, Taylor blocks and explicit numerical
  failure/status checks. The commutator diagnostic is not an error bound.
- Include native tests, a working Fortran example, two-mode SciPy reference
  harness, stress cases, documentation and independent provenance notes.
- No extra runtime dependencies or third-party implementation code.

## 2026-11-07 — Dense exponential condition numbers

- Add real and complex `expm_cond`, estimating SciPy's relative Frobenius
  condition number via existing SciFort Fréchet derivative kernels and a
  fully reorthogonalized matrix-free Lanczos solver.
- Add exact diagonal divided-difference evaluation, scalar-shift
  normalization for large offset inputs, explicit iteration/status controls.
- Add native tests, reproducible example, randomized and analytic SciPy
  reference checks, API/validation documentation, and provenance records.
- No additional numerical dependencies or third-party source.

## 2026-11-07 — Dense matrix trigonometric/hyperbolic functions

- Add `scifort_linalg_trigonometric` with real/complex `sinm`, `cosm`,
  `tanm`, `sinhm`, `coshm`, and `tanhm`.
- Use paired matrix Taylor series with power-of-two scaling and the
  double-angle identities, avoiding exponential subtraction near zero.
- Reuse the existing small-dense pivoted inverse for tangent solves, with
  explicit near-pole failure statuses and post-solve residual checks.
- Add native tests, runnable example, independent randomized SciPy reference
  checks and API/validation/provenance documentation. No new dependencies.

## 2026-11-07 — Dense matrix logarithm and fractional matrix powers

- Add real/complex `logm` and `fractional_matrix_power`, including real-only
  output convenience wrappers, with explicit status reporting.
- Implement independent inverse scaling-and-squaring and matrix-atanh series
  using existing SciFort square-root and exponential routines, with
  singularity screening and normalized exponential reconstruction checks.
- Add exact integer-power and identity fast paths, and normalized half-power;
  cover real/complex, triangular, branch-cut, singular, and extreme-scale cases.
- Add native regression, standalone example, randomized SciPy comparison
  harness, API/limitations documentation, and validation/provenance records.
- No new external dependencies or imported third-party source.

## 2026-11-06 — Polar decomposition and matrix exponential Fréchet derivatives

- Add `scifort_linalg_polar` with real/complex `polar` for rectangular,
  square, and rank-deficient matrices, and both left/right PSD factors.
- Add `scifort_linalg_expm_frechet` with real/complex `expm_frechet`,
  returning `exp(A)` and the directional derivative via differentiated
  scaling-and-Taylor recurrences.
- Include native tests, an executable example, randomized independent SciPy
  comparisons, documentation, validation records, and provenance notes.
- No external numerical dependencies or third-party source incorporated.

## 2026-11-06 — Dense matrix square root and matrix sign functions

- Add real/complex `sqrtm` and `signm`, with real-input complex-output
  principal square-root semantics and a `sqrtm_real` convenience routine.
- Adapt SciFort's MIT-licensed complex Schur QR for triangular square roots;
  add scaled Newton matrix-sign iteration with true residual checks.
- Add focused and randomized SciPy/NumPy comparisons, executable example,
  documentation, and review of permissively licensed external libraries.

## 2026-11-06 — Dense continuous and discrete algebraic Riccati equations

- Add `scifort_linalg_riccati` with real/complex `solve_continuous_are` and
  `solve_discrete_are`, including generalized descriptor `E` and cross-cost `S`.
- Implement a scaled Hamiltonian matrix-sign method for CARE and a
  structure-preserving doubling method for DARE; use existing Sylvester and
  Lyapunov routines for optional refinement, and check the original residual.
- Add closed-loop stability screening, explicit failure statuses, native tests,
  a standalone example, and independent SciPy reference comparisons.
- Document the 20-state limit, method differences from SciPy's Schur/QZ
  implementation, and numerical limitations near stability boundaries.

## 2026-11-06 — Dense Sylvester and Lyapunov matrix equations

- Add `scifort_linalg_matrix_equations` with real and complex binary64
  `solve_sylvester`, `solve_continuous_lyapunov`, and
  `solve_discrete_lyapunov` generics, matching SciPy equation conventions.
- Use independent partial-pivoted Kronecker linearization without external
  BLAS/LAPACK; reject singular, nonfinite, mismatched and oversized systems
  with SciFort's existing `scifort_linalg` status codes.
- Add native tests, a standalone example, and randomized SciPy comparisons
  in checked and optimized GNU Fortran builds. Document the 400-unknown
  size cap and the cubic-in-number-of-unknowns performance limitation.

## 2026-11-06 — Complex restarted Krylov matrix-function actions

- Add `scifort_sparse_linalg_complex_funm` with `funm_multiply_krylov_complex`
  (CSR/CSC vector and block overloads) and
  `funm_multiply_krylov_complex_operator` for matrix-free complex callbacks.
- Use two-pass complex Arnoldi or Hermitian Lanczos, restarted projected
  matrix-function evaluations, complex time scales, update-based stopping,
  callback failure/status propagation, and matrix-vector counts.
- Add independent SciPy/NumPy comparisons of exp/sin/polynomial actions,
  native tests, a compilable example, API and validation documentation.

## 2026-11-06 — Complex sparse exponential actions and time intervals

- Add `scifort_sparse_linalg_complex_expm` with complex CSR, CSC and
  matrix-free `expm_multiply_complex` overloads for vectors and dense blocks.
- Add `expm_multiply_complex_interval` and
  `expm_multiply_complex_operator_interval` with bidirectional propagation
  from the time sample closest to zero.
- Use a shifted, bounded Taylor expansion with conservative remainder checks,
  configurable tolerances, matrix-vector counts, and callback-error reporting.
- Add a native regression, independent randomized SciPy/NumPy comparisons,
  an executable example and API/validation documentation.

## 2026-11-06 — Matrix-free complex Arnoldi and reusable shift-invert LU

- Extend the complex `eigs_complex` generic with `eigs_complex_operator(n,k,matvec,...)`,
  accepting a complex matrix-free callback and optional shifted inverse callback.
- Share Arnoldi/Rayleigh–Ritz, selection, restart, and true-residual logic across
  CSR, CSC, and matrix-free paths, rather than duplicating eigensolver kernels.
- Add `complex_eigs_shift_cache` and `eigs_complex_prepare_shift` for CSR/CSC;
  reject mismatched shifts, modified matrices, and incomplete factors, while
  reusing the same pivoted LU without another factorization.
- Add native regression, 2,500-dimensional matrix-free stress test, independent
  randomized NumPy comparisons, documentation, and callback-error validation.

## 2026-11-06 — Complex shift-invert Arnoldi eigenvalue searches

- Add optional complex `sigma` to the CSR/CSC `eigs_complex` interface, using
  one pivoted complex sparse LU factorization of `A-sigma*I`, reused for all
  Arnoldi operator applications.
- Rank Ritz values by the transformed spectrum `mu=1/(lambda-sigma)` and
  return original eigenvalues and true residuals `||A*v-lambda*v||`.
- Distinguish singular shifted matrices (`info=-3`) from invalid options or
  numerical breakdown; retain the unshifted API and source conventions.
- Add native regression tests for complex interior spectra, singular/NaN
  shifts, CSR/CSC, limited iteration budgets, and a 70-by-70 truncated-space
  banded problem; add randomized NumPy/SciPy comparison harness.

## 2026-11-06 — Nonsymmetric complex sparse eigenpairs

- Add `scifort_sparse_linalg_complex_eigs`, providing `eigs_complex` for
  square complex binary64 CSR/CSC matrices with `LM`, `SM`, `LR`, `SR`, `LI`,
  and `SI` ordering, right eigenvectors, true residuals, and restart controls.
- Reuse SciFort's MIT-licensed complex projected Schur/QR algorithm and
  implement explicit thick restarted Arnoldi without external libraries.
- Add targeted native tests (including repeated eigenvalues, restarting,
  and nonconvergence), independent randomized NumPy/SciPy comparisons,
  a documentation example and API/validation notes.
- `SM` currently requires a full basis; no shift-invert or ARPACK-style
  implicit restarting is implemented.

## 2026-11-06 — Complex sparse SVD

- Added `scifort_sparse_linalg_complex_svds` with CSR/CSC overloads of
  `svds_complex`, returning complex U, real singular values, and complex V^H.
- Implemented complex Golub--Kahan search and one-sided Jacobi projected SVD,
  full reorthogonalization, true two-sided residual tests and rank-deficient
  zero-singular-vector completion, without dense normal matrices.
- Added 4,000-by-120 sparse native stress tests, independent randomized
  NumPy SVD comparisons and SciPy ARPACK `svds` cross-checks.
- When the full smaller-side basis is requested, complete it before accepting
  convergence, guarding against missed repeated dominant singular values.
- No additional runtime dependencies or third-party source imported.

## 2026-11-06 — Complex sparse LSMR

- Add `scifort_sparse_linalg_complex_lsmr` with CSR and CSC overloads of
  `lsmr_complex`. Supports rectangular complex binary64 matrices, damping,
  nonzero initial guesses, SciPy-style stopping codes, and true residual norms.
- Register Fortran regression tests and randomized independent SciPy/NumPy
  reference comparisons, plus a BSD-3-Clause source attribution notice.
- Retain the existing real LSQR/LSMR and complex LSQR APIs unchanged.

## 2026-11-06 — Complex sparse LSQR for rectangular least squares

- Added `scifort_sparse_linalg_complex_lsqr` with `lsqr_complex` CSR/CSC
  overloads, Golub–Kahan bidiagonalization using the conjugate-transpose
  operator, optional damping and nonzero `x0`, true residual diagnostics,
  variance-diagonal estimates, and failure status codes.
- Added native regression tests and randomized comparisons against independent
  SciPy and NumPy references; registered the new test in the standalone native
  sparse runner. See `docs/sparse_complex_lsqr.md`.

## 2026-11-06 -- Complex pivoted sparse LU and direct solves

- Add `scifort_sparse_linalg_complex_direct` with complex pivoted sparse LU
  for canonical CSR and CSC, dynamic fill-in and reusable factors.
- Add `splu_complex_solve` (normal, transpose and Hermitian-adjoint;
  multiple RHS) and one-shot `spsolve_complex` for vector and dense-block RHS.
- Add native singular/pivot/error tests and streamed randomized SciPy/NumPy
  comparisons for both checked and optimized GNU Fortran builds.
- No new runtime dependencies; document differences from SuperLU ordering.

## 2026-11-05 -- Complex ILU(0), BiCG, and QMR

- Add `scifort_sparse_complex_ilu`, with no-pivot/no-fill complex CSR LU,
  stored-diagonal validation, optional real diagonal shift, and forward and
  conjugate-adjoint triangular solves.
- Add complex `bicg_complex` and `qmr_complex` to the CSR/CSC generic solver
  module; support no preconditioner, Jacobi, or a reusable complex ILU(0).
- Extend complex `bicgstab_complex` and `gmres_complex` with optional `ilu=`
  without changing positional arguments in existing calls.
- Add focused Fortran regression and streamed randomized SciPy/NumPy tests.
  Retain BSD-3-Clause SciPy provenance for the QMR-derived recurrence.

## 2026-11-05 -- Complex CSR/CSC storage and complex Krylov solvers

- Add `scifort_sparse_complex` with canonical one-based complex128 CSR/CSC
  storage, COO/complex N-D COO interoperability, dense conversions, sparse
  matrix/vector products, transpose, conjugation and Hermitian transpose.
- Implement `scifort_sparse_linalg_complex` with complex CG (Hermitian positive
  definite), BiCGSTAB and right-preconditioned restarted GMRES; support both
  CSR and CSC, optional Jacobi preconditioning, true residual tests and
  diagnostic outputs. Keep the existing real solver interfaces unchanged.
- Add bounds-checked native tests, independently checked randomized SciPy/
  NumPy comparisons, API notes and provenance; no additional dependencies.

## 2026-11-05 -- Complex N-dimensional sparse COO arithmetic

- Add `scifort_sparse_coo_nd_complex`, a complex binary64 N-D coordinate format
  with duplicate coalescing, zero elimination, C/F reshape and dense conversion,
  conjugation, Hermitian transpose and safe conversion to/from real COO.
- Add `scifort_sparse_coo_nd_complex_ops` with complex reductions, broadcasting,
  addition, elementwise multiplication, `tensordot` without conjugation, and
  conjugated `vdot`.
- Preserve one-based coordinates, large sparse shapes, finite-input and
  output-limit validation; retain the existing real-valued API unchanged.
- Add a native Fortran test, streaming SciPy/NumPy reference comparisons,
  documentation and provenance records. No new dependencies.

## 2026-11-05 -- N-D sparse Boolean broadcasting and reductions

- Add `scifort_sparse_coo_nd_bool_reduce` with singleton and leading-axis
  broadcasting, sparse AND/OR/XOR/difference, and exact axis-wise counts.
- Represent count results by `coo_nd_bool_count_array` (int64 baseline plus
  sorted coordinate deviations); add axis-wise and scalar `ANY`/`ALL`.
- Preserve empty-axis/vacuous reduction semantics and support shapes whose
  total number of entries exceeds signed 64-bit limits.
- Include focused Fortran regression and randomized NumPy comparisons in
  checked/optimized compiler builds; original MIT code, no dependencies.

## 2026-11-05 -- Complement-aware N-D COO Boolean mask algebra

- Add `scifort_sparse_coo_nd_bool`, with canonical sorted coordinate
  exceptions to an implicit Boolean background, including large complements.
- Implement AND/OR/XOR/difference/NOT by sparse coordinate-set merging,
  one-based lookup, exact 64-bit true counts where representable, and
  conversion to/from real COO nonzero masks.
- Add scalar comparisons that are true at implicit zeros (`==0`, `/=1`),
  sparse source filtering that retains duplicates/explicit zeros, strict
  validation, and configurable output caps.
- Include native regression tests, randomized NumPy/SciPy comparisons,
  API documentation and original MIT-licensed source without new dependencies.

## 2026-11-05 -- Sparse N-dimensional predicates and coordinate selection

- Add `scifort_sparse_coo_nd_filter` for scalar comparison masks, logical
  nonzero counts, sparse-valued COO masks, coordinate-set filtering and
  coordinate-value filtering without dense allocations.
- Coalesce duplicate coordinates before logical decisions, preserve source
  duplicates and explicit zeros during filtering, and support inverted masks.
- Reject predicates true on implicit zeros, invalid shapes/indices, nonfinite
  thresholds and output-size overruns; use stable coordinate sorting and
  binary search without dense intermediates.
- Include native regression tests, randomized SciPy/NumPy comparisons in
  checked/optimized configurations, and original MIT-licensed implementation.

## 2026-11-05 -- N-dimensional COO boolean indexing and transactional scatter

- Add `scifort_sparse_coo_nd_edit` for full boolean masks, axis masks,
  last-write-wins scatter assignment, additive scatter, and masked scalar
  assignment. Use one-based COO coordinates and C-order flat full masks.
- Sort and coalesce COO entries after successful nonempty updates, preserve
  untouched stored zeros, optionally retain updated zeros, and roll back
  invalid, nonfinite, overflowing or oversized edits.
- Add standalone Fortran regression coverage and randomized independent
  NumPy/SciPy comparisons in checked and optimized builds. Original MIT code
  with no new external numerical dependency.

## 2026-11-05 -- N-dimensional COO indexing and slicing

- Add `scifort_sparse_coo_nd_index` with one-based integer, inclusive range,
  and integer-array selectors, including negative strides, duplicate selection
  indices, and optional dimension retention.
- Add sparse pointwise `coo_nd_gather` for paired coordinate tuples, distinct
  from orthogonal multi-axis selection in `coo_nd_select`.
- Preserve unsorted duplicate COO values and explicitly stored zeros; validate
  selector bounds and apply an adjustable generated-entry cap without dense
  materialization, including extremely large input dimensions.
- Add independent randomized NumPy/SciPy comparisons, native Fortran tests,
  and registration with the checked/optimized sparse regression runner.
  Original MIT implementation; no new numerical dependencies.

## 2026-11-05 -- N-dimensional COO reductions, broadcasting and contraction

- Add `scifort_sparse_coo_nd_ops` for sparse scalar, single-axis and
  multiple-axis sums, with optional `keepdims` and canonical COO results.
- Add singleton-dimension N-D broadcasting, elementwise sparse addition and
  multiplication, and `tensordot` on selected axes. The contraction uses a
  sorted join of matching contracted-coordinate groups, rather than dense
  arrays or a scan of all nonzero pairs.
- Apply configurable caps to replicated entries and matched contraction
  contributions, preserving input arrays and rejecting nonfinite results.
- Include focused native tests, streaming independent SciPy/NumPy reference
  comparisons, and integration into the checked/optimized sparse test runner.
  No new numerical dependencies or imported third-party source.

## 2026-11-05 — N-dimensional sparse COO arrays

- Add `scifort_sparse_coo_nd` with one-based `coo_nd_array` coordinates for
  arbitrary positive ranks; preserve duplicate and explicitly stored entries.
- Add C/F-order reshape, axis permutations, stable canonicalization,
  duplicate summation, zero elimination, scalar lookups, dense/flat conversion,
  and lossless conversion to and from existing 2D COO matrices.
- Add focused bounds-checked tests and 1,100 independent SciPy 1.17 randomized
  cases covering ranks 1–7 in checked and optimized gfortran builds.

## 2026-11-05 — Incremental DOK and LIL sparse formats

- Add `scifort_sparse_dok_lil`, with hash-based DOK and sorted, dynamically
  allocated row-list LIL matrices. Support one-based, real binary64 entry
  get/set/add/delete, explicit zero elimination, clear, matvec, transpose,
  and conversions to/from dense, COO, CSR, CSC and each other.
- Preserve stored zero entries imported from CSR/CSC/COO (including canceled
  COO duplicates), while explicit zero assignments delete stored keys.
  Use canonical CSR output with sorted and unique column indices.
- Add focused native Fortran tests and independent randomized SciPy 1.17
  comparisons for edit histories, stored counts and compressed structure.
  Checked and optimized GNU Fortran runs passed 500 cases each.
- No new runtime dependencies or external implementation source.

## 2026-11-05 — COO sparse matrices and sparse.find

- Add `scifort_sparse_coo` with `coo_matrix` and coordinate-format construction,
  CSR/CSC/dense conversions, transpose, matrix-vector multiplication, indexing,
  in-place duplicate summation and explicit-zero elimination. Unsorted input,
  duplicate coordinates and stored zeros are preserved until canonicalized.
- Add `sparse_find` overloads for COO, CSR and CSC, matching SciPy's behavior
  of coalescing duplicates and excluding stored zeros and cancellations.
- Add independent native Fortran regressions and SciPy 1.17 randomized
  comparisons in checked and optimized GNU Fortran builds; no added dependency.

## 2026-11-05 — Krylov actions for arbitrary matrix functions

- Add `scifort_sparse_linalg.funm_multiply_krylov` for real CSR and
  caller-provided matrix-free operators, following the restarted
  Arnoldi/Lanczos construction of Eiermann--Ernst. A dense function callback
  evaluates `f(H)` for the small projected matrix.
- Support `t`, `rtol`, `atol`, Krylov cycle size/budget, symmetric Lanczos,
  iteration/matvec counting, and last-update norm. Preserve the SciFort
  `info` convention with explicit nonconvergence/failure codes.
- Handle oversized Krylov dimensions, empty/zero-vector inputs, callback
  failures, NaN/Inf intermediate values and dynamically grown projected
  matrices without dense materialization of the input operator.
- Add native test, randomized SciPy 1.17 comparisons, and native test-runner
  registration in checked and optimized configurations. All 32 native sparse
  tests and 600 independent SciPy comparisons passed per compiler mode.
  No new dependency.

## 2026-11-05 — Grid Laplacians and sparse structural queries

- Add `scifort_sparse_laplacian_nd`, including an initialized real
  `LaplacianNd`-style operator with matrix-free matvec/matmat, CSR assembly,
  and analytic Dirichlet/Neumann/periodic eigenvalues and eigenvectors.
- Use C-order grid indexing (last axis fastest), including length-one axes
  and correctly coalesced periodic neighbors for length-two axes.
- Add CSR/CSC `spbandwidth` and `is_sptriangular` to
  `scifort_sparse_linalg`, counting explicitly stored zeros.
- Add standalone native tests, two randomized SciPy comparison drivers,
  documentation, and registration with the native sparse test runner.
  Both builds passed 180 Laplacian comparisons, 300 structural comparisons,
  and all 31 sparse native regression programs.
- Do not imitate SciPy 1.17 edge-case inconsistencies in its sparse
  singleton-axis constructor or truncated eigenspectrum optimization.
- Original MIT implementations; no new runtime dependencies or imported
  third-party code.

## 2026-11-05 — Restarted augmented/recycled sparse Krylov solvers

- Add `scifort_sparse_linalg.lgmres` (LGMRES) with optional retained
  correction directions, and `gcrotmk` (GCROT(m,k)) with orthogonalized
  `C=A*U` recycling and oldest-first truncation.
- Share a right-preconditioned flexible Arnoldi/least-squares kernel with
  two-pass MGS, Givens rotations and original-system residual checks.
- Support Jacobi and existing ILU(0) preconditioning, portable reuse
  matrices `outer_v`/`recycle_u`, optional iteration/residual diagnostics,
  starting guesses and SciFort-style status codes.
- Add standalone native regression and independent SciPy 1.17 comparisons
  of 180 randomized systems with two related right-hand sides per build.
  Both checked and optimized modes passed all 360 comparisons and all
  30 native sparse regression executables. No new dependencies.
- Original MIT-licensed implementation based on published algorithmic
  descriptions; no third-party numerical source was imported.

## 2026-11-05 — Sparse QMR and transpose-free QMR

- Add `scifort_sparse_linalg.qmr` for nonsymmetric CSR systems with explicit
  transpose products and an optional right Jacobi preconditioner.
- Add `scifort_sparse_linalg.tfqmr` using the existing SciFort matrix-free
  recurrence with right Jacobi or ILU(0) preconditioning (not simultaneously).
- Check the unpreconditioned true residual for both methods rather than
  accepting a potentially optimistic recursive or quasi-minimal estimate.
- Add native edge-case tests and an independent SciPy reference driver;
  preserve existing status codes, optional diagnostics and default tolerances.
- Retain BSD-3-Clause SciPy QMR and TFQMR attributions and full license terms.

## 2026-11-05 — BiCG, CGS and sparse Kronecker sums

- Add `scifort_sparse_linalg.bicg` with explicit transposed CSR products,
  optional Jacobi preconditioning, and verified true-residual convergence.
- Add `scifort_sparse_linalg.cgs`, with Jacobi or reusable ILU(0)
  preconditioning, finite-value guards, zero-iteration and breakdown handling.
- Add `scifort_sparse.csr_kronsum` for the SciPy-compatible Kronecker sum
  `kron(I_n,A) + kron(B,I_m)` using direct sorted CSR row merging, preserving
  explicitly stored zeros and rejecting nonfinite overlaps or overlarge output.
- Register native checked/optimized tests; compare 180 randomized Kronecker
  sums and 120 nonsymmetric linear systems in each mode against SciPy 1.17.
- Retain the dependency list and all pre-existing public APIs. Implementations
  are original MIT-licensed SciFort code based on published mathematical
  algorithms; no third-party source code has been imported.

## 2026-11-05 — Sparse random and diagonal-array constructors

- Add `scifort_sparse_constructors` with `csr_random`, `csr_random_array`, `csr_rand`,
  `csr_diags_array`, and a ragged `sparse_diagonal` helper.
- Reuse the existing explicit MRG32k3a `rng_state`; sample without
  replacement using O(nnz) auxiliary storage; support optional normal
  values or caller-supplied finite values.
- Support SciPy-like shape inference, singleton diagonal broadcasting,
  diagonal truncation, explicit zeros, and strictly unique offsets.
- Add native and randomized SciPy 1.17 comparison tests, and integrate
  native test into the sparse regression runner. No new dependencies.

## 2026-11-05 — Sparse structural deletion and scalable sparse assignment

- Add `csr_delete_entries`, `csr_delete_rows`, `csr_delete_columns` and
  C-order `csr_reshape` for general canonical CSR matrices.
- Make sparse RHS `csr_assign` avoid materializing the full Cartesian
  submatrix. Implicit source zeros now remove the old target entry, while
  explicitly stored zero values remain stored. Repeated indices retain
  Fortran column-major last-write-wins semantics.
- Preserve invalid-request transactionality and one-based CSR ordering.
  Support sparse matrices with a billion logical columns and a test of a
  50,000-by-50,000 sparse submatrix replacement.
- Add standalone native regression and randomized NumPy/SciPy comparisons
  for five new/revised operations. No new dependencies or imported code.

## 2026-11-05 — General block matrices and sparse scatter assignment

- Implement `csr_bmat` and alias `csr_block` for rank-two grids of rectangular
  CSR blocks, with unallocated placeholders for missing blocks and optional
  explicit row/column sizes where inference is impossible.
- Implement `csr_set`, `csr_assign_coo`, and generic `csr_assign` for scalar,
  dense and sparse RHS assignment to Cartesian-indexed submatrices.
- Preserve canonical CSR storage and explicitly stored zeros; use stable
  last-write-wins semantics for duplicate assignments and transactional
  validation of all updates. Coordinate scatter allocates no dense matrix.
- Add standalone checked/optimized native tests and independent SciPy/NumPy
  randomized block/assignment comparisons. No new dependencies.

## 2026-11-05 — Sparse CSR assembly, Kronecker products, indexing and slicing

- Add `csr_hstack`, `csr_vstack`, `csr_block_diag`, and `csr_kron` for
  sparse matrix assembly, including rectangular and empty matrix shapes.
- Add `csr_get` (scalar lookup), `csr_select` (Cartesian gather/reordering,
  repeated indices), and `csr_slice` (signed one-based inclusive strides).
- Preserve canonical sorted CSR storage and explicitly stored zero entries.
  Guard dimensions/nnz with 64-bit arithmetic; reject nonfinite products.
- Index selection uses memory proportional to the number of selected columns,
  not the full logical matrix width; test a one-billion-column sparse input.
- Register native Fortran regression and randomized SciPy 1.17 references
  covering positive/negative strides, duplicate coordinates and zero-sized blocks.
- No additional dependencies or third-party source implementations.

## 2026-11-05 — Sparse DIA/BSR formats and matrix norms

- Add canonical diagonal (`dia_matrix`) and block CSR (`bsr_matrix`) representations,
  rectangular blocks, validation, CSR conversions, dense conversions and direct matvec.
- Add CSR/CSC `scifort_sparse_linalg.norm` supporting Frobenius, column/row
  ±1/±infinity norms, and approximate spectral 2-norm through existing `svds`.
- Add one-based `norm_axis` for vector norms by column or by row.
- Include checked/optimized native tests, a randomized SciPy/NumPy comparison
  driver and edge-case validations for empty matrices, stored zeros, duplicate
  coordinates and large sparse bands.
- Preserve prior sparse APIs, existing provenance and build dependencies.

## 2026-11-05 — Sparse integer matrix powers and inverse

- Add `scifort_sparse_linalg.matrix_power` for nonnegative integer exponents
  via exponentiation by squaring, preserving canonical real CSR/CSC storage.
- Add `scifort_sparse_linalg.inv` with pivoted sparse LU factor reuse,
  identity-block solves, a diagonal shortcut, and optional absolute dropping.
- Reject negative exponents, nonsquare/invalid storage, and singular inverses;
  surface overflow and invalid options with existing status codes.
- Register checked/optimized native tests and NumPy/SciPy comparisons for
  rectangular rejection, stored zeros, huge powers, repeated products,
  pivoting, empty matrices, singularity, and large sparse diagonal inputs.
- No new external dependencies or third-party source implementation code.

## 2026-11-05 — LOBPCG for sparse symmetric and generalized eigenproblems

- Add `scifort_sparse_linalg.lobpcg` for multiple algebraically largest/smallest
  eigenpairs using a locally optimal Rayleigh--Ritz subspace `[X,W,P]`.
- Support SPD mass matrices `B`, constraint directions, a sparse SPD
  preconditioner action or optional positive Jacobi preconditioning.
- Return ordered values, B-orthonormal vectors, and checked full or
  constraint-projected residual norms; retain the best iterate on failure.
- Add repeated eigenvalue, large CSR, singleton, generalized and constraint
  regressions; independently compare randomized problems with dense
  generalized eigensolvers and SciPy's LOBPCG.
- Preserve all existing sparse APIs. No new runtime libraries or
  third-party source code.

## 2026-11-05 — Sparse threshold ILU with controlled fill

- Add `scifort_sparse_linalg.spilu` as an original natural-order ILUT implementation for canonical real CSR and CSC matrices, with configurable relative dropping, per-row fill budgets, and optional diagonal perturbation.
- Add `spilu_solve` vector and multiple-RHS adapters and reuse existing `ilu0_factor` storage for direct interoperability with GMRES/BiCGSTAB preconditioners.
- Cover zero-fill versus full-fill behavior, invalid arguments, pivot failures and shift recovery, sparse systems up to 1,200 rows, and randomized NumPy/SciPy reference checks.
- Document the important difference from SciPy's pivoted SuperLU implementation: no permutations or equilibration; factor sparsity patterns are not expected to coincide.
- Add new native test to the sparse regression runner. No new dependencies or third-party source code.

## 2026-11-05 — Sparse triangular and pivoted direct solvers

- Add canonical CSR/CSC `spsolve_triangular` for lower/upper systems with unit
  diagonal support, accepting vector or block right-hand sides.
- Add real-binary64 sparse LU with partial row pivoting and dynamically inserted
  fill (`splu`), reusable factors (`splu_factor`, `splu_solve`), and a general
  one-call sparse direct solve (`spsolve`).
- Add explicit negative singular-pivot status and optional absolute pivot
  threshold. Preserve previous iterative solver defaults and dependencies.
- Add native regression and randomized SciPy solver comparisons, including
  factor reconstruction, CSC/CSR formats, singular matrices, and large sparse
  tridiagonal systems. Extend the standalone sparse test runner.
- No third-party solver code copied; no extra dependencies.

## 2026-11-04 — CSR 1-norm estimation and full sparse exponential

- Added `scifort_sparse_linalg.onenormest`, a deterministic block sign/adjoint
  1-norm lower-bound estimator for real square CSR matrices with optional
  certified witness vectors. For `t >= n`, exact column sums yield the
  1-norm in O(nnz), without an n-by-n temporary matrix.
- Added `scifort_sparse_linalg.expm` returning a canonical CSR matrix,
  assembled from identity-block calls to `expm_multiply`, with optional
  absolute drop tolerance, block size, and Taylor tolerance. Diagonal inputs
  are handled directly in O(n+nnz).
- Added native and randomized SciPy/NumPy regression drivers, test runner
  registration, API documentation, and provenance record.
- No new runtime dependencies or third-party source code.

## 2026-11-04 — sparse matrix-exponential action (`expm_multiply`)

- Add real-CSR `scifort_sparse_linalg.expm_multiply` for vectors and multiple
  right-hand sides, with configurable time, trace-shift hint, local Taylor-tail
  tolerance, scaling-step and sparse-product diagnostics.
- Add `expm_multiply_interval` for evenly spaced time grids with optional
  endpoint exclusion. Anchor each grid at the time closest to zero and
  propagate in either direction to avoid unrecoverable underflow.
- Use shifted/scaled truncated Taylor actions with an explicit infinity-norm
  upper bound on omitted terms. Never materialize a dense exponential.
- Add native tests for oscillatory/triangular operators, large sparse matrices,
  numerical corner cases, invalid inputs, and time grids. Add randomized
  SciPy and independent dense exponential reference comparisons.
- No new runtime dependencies or third-party source code.

## 2026-11-04 — partial sparse singular value decomposition (`svds`)

- Add `scifort_sparse_linalg.svds` for real rectangular canonical CSR matrices,
  using two-pass-reorthogonalized Golub--Kahan search directions and a stable
  one-sided Jacobi projected singular-value decomposition. Both tall and wide
  matrices are supported without constructing a dense matrix or `A^T A`.
- Support largest-magnitude partial singular values (`which='LM'`) and
  full-basis smallest values (`which='SM'`); return left and right vectors,
  one-based CSR-compatible results, true two-sided residual norms, explicit
  nonconvergence status, iteration counts and a deterministic initial vector.
- Handle zero matrices, repeated singular values, rank deficiency, explicit
  CSR zeros and duplicates from COO construction. Add native regression and
  randomized NumPy/SciPy reference comparisons plus 1,200-by-900 sparse test.
- No added external numerical dependencies or third-party source code.

## 2026-11-04 — block-restarted Arnoldi for clustered nonsymmetric spectra

- Add optional `block_size=1..k` to `eigs_arnoldi` and to the
  `eigs(method='arnoldi'/'arnoldi_schur')` dispatcher, preserving the original
  single-vector path by default.
- Expand a block Krylov space from deterministic independent seeds, restart
  from separate Ritz residuals, and retain the existing soft/Schur locking,
  complex shift-invert, and true-residual convergence checks.
- Add native tests for exact repeated real/complex eigenspaces, clustered
  roots, Schur coupling, shift-invert, validation and 1,400-row sparse CSR.
  Add 342 SciPy dense-reference randomized comparisons across two build modes,
  plus pass the earlier 14 native and Arnoldi comparison suites.
- No added external dependencies or third-party solver source code.

## 2026-11-04 — nonsymmetric Schur-subspace locking for restarted Arnoldi

- Add optional `schur_locking=.true.` and `schur_locked_count` to
  `eigs_arnoldi`; add `eigs(...,method='arnoldi_schur')` dispatcher.
- Lock only certified approximately invariant nonsymmetric subspaces, using
  orthonormal Schur vectors and the active operator (including complex
  shift-invert); preserve the coupling to the unlocked subspace in the full
  Ritz projection. Keep legacy soft locking and default eigensolvers unchanged.
- Add dedicated native regression tests with nonnormal right-eigenvector
  coupling, complex spectra, interior shifts and two-dimensional Schur locks.
  Add randomized SciPy/dense reference checks in both build modes.
- No external solver sources, BLAS/LAPACK requirements or new runtime deps.

## 2026-11-04 — hard-locked restarted symmetric sparse eigensolver

- Add opt-in `eigsh_lanczos(...,hard_locking=.true.)` and
  `eigsh(...,method='lanczos_locked')`, retaining the existing soft-locking
  default and positional API compatibility.
- Freeze truly converged symmetric Ritz vectors, explicitly project later
  search spaces onto their orthogonal complement, and track persistent
  hard locks separately from restart cycles and operator products.
- Add optional deterministic multi-seed `block_size=1..k` on the direct
  hard-locking API for repeated eigenvalue eigenspaces. The default is one
  seed; multiple seeds can require more iterations.
- Support real MINRES-based shift-invert, four existing selection modes,
  residual validation against the original CSR matrix, and independent
  checked/optimized comparisons with NumPy and SciPy.
- Keep nonsymmetric Arnoldi soft locked: orthogonal deflation of right
  eigenvectors alone is invalid for general nonnormal operators.

## 2026-11-04 — normal-mode shift-invert for restarted Lanczos and Arnoldi

- Extend `eigsh_lanczos` with real `sigma`, shifted MINRES inner solves,
  `inner_rtol` and `inner_maxiter`. Preserve inverse-spectrum `LM`, `LA`,
  `SA` and full-subspace `SM` selection with original-matrix residuals.
- Extend `eigs_arnoldi` with complex `sigma` and restarted complex-GMRES
  inverse applications, preserving complex Ritz values, true residuals and
  convergence status. Expose inner tolerance, cap, and restart options.
- Forward shifts through `eigsh(method='lanczos')` and
  `eigs(method='arnoldi')` without changing their default algorithms.
- Add dedicated native edge tests and randomized NumPy/SciPy eigenpair
  comparisons in checked and optimized builds. No new dependency or
  copied external solver code. Hard locking and user-supplied inverse
  operators remain unsupported.

## 2026-11-04 — restarted symmetric Lanczos and eigensolver method dispatch

- Add `scifort_sparse_linalg.eigsh_lanczos`, a real-symmetric restarted
  Krylov eigensolver with thick Ritz-vector retention, two-pass full
  reorthogonalization, soft locking, true eigenpair-residual checks,
  convergence diagnostics, and one-based CSR storage.
- Add optional `method='lanczos'` on `eigsh` and `method='arnoldi'` on
  `eigs` while retaining `method='subspace'` as the backward-compatible
  default. Shift-invert support is added in the subsequent update above.
- Add native Lanczos/method-selection regression tests and independent
  checked/optimized comparisons with NumPy and SciPy, including clustered
  spectral edges and explicit nonconvergence checks.
- No new external dependencies or copied third-party solver source.

## 2026-11-04 — thick-restarted Arnoldi for nonsymmetric sparse eigenvalues

- Add the independent `scifort_sparse_linalg.eigs_arnoldi` entry point,
  leaving all existing `eigs` and `eigsh` signatures unchanged.
- Build an orthonormal Arnoldi search space, extract complex Ritz pairs,
  and retain selected Ritz vectors between restart cycles. Retain converged
  Ritz vectors first (soft locking) and require true residual convergence.
- Support `LM`, `LR`, `SR`, `LI`, and `SI` with reduced search spaces;
  accept `SM` only with `ncv=n`. Return iteration, matvec, residual and
  soft-lock diagnostics; validate invalid or nonfinite inputs.
- Add native checked/optimized regressions including a 1,200-row CSR system,
  and randomized SciPy eigenvalue comparisons plus clustered-spectrum tests.
- The new solver is **not** implicitly restarted Arnoldi, does not provide
  hard invariant-subspace locking, and has no shift-invert path yet.

## 2026-11-04 — complex shift-invert and extended sparse eigenvalue selection

- Extend `scifort_sparse_linalg.eigs` with complex `sigma` shift-invert on
  real CSR inputs. A new internal complex restarted GMRES solves the shifted
  system without a dense inverse; expose `inner_rtol`, `inner_restart`, and
  `inner_maxiter`, and verify original-matrix Ritz residuals.
- Support reduced-subspace nonsymmetric `which='SM'` for nonsingular
  matrices through automatic zero-shift inversion (previously required
  `ncv=n`). Retain the full-subspace fallback for singular matrices.
- Extend symmetric `eigsh` shift-invert beyond LM to LA and SA; accept SM
  with `ncv=n`, with transformed `1/(lambda-sigma)` selection semantics.
- Add native regression and randomized SciPy/NumPy comparison scripts for
  real and complex shifts, nonsymmetric smallest-magnitude modes, and
  four symmetric shift-invert selection modes. No additional third-party
  source or libraries. Implicitly restarted Arnoldi/Lanczos and generalized
  eigenproblems remain outstanding.

## 2026-11-04 — nonsymmetric sparse eigenvalues and symmetric shift-invert

- Implement `scifort_sparse_linalg.eigs` for real nonsymmetric CSR matrices
  with complex eigenpairs. Use complex block orthogonal iteration, real sparse
  matvecs, complex Householder Hessenberg reduction, Wilkinson-shift QR Schur
  reduction, and eigenvector reconstruction. Support LM, LR, SR, LI, SI,
  and SM with a complete ncv=n subspace; return true eigenpair residuals.
- Extend the original symmetric `eigsh` with optional `sigma` shift-invert
  for `which='LM'`, selecting the eigenvalues nearest the shift. Reuse the
  previously implemented shifted MINRES as the iterative inverse operator;
  add independent inner tolerance and iteration controls.
- Add native edge-case and reduced-subspace tests and randomized NumPy/SciPy
  eigenvalue, eigenvector, and residual comparisons. This update introduces
  no BLAS/LAPACK dependency or third-party source code.

## 2026-11-04 — MINRES and first symmetric sparse eigensolver

- Add `scifort_sparse_linalg.minres` using the three-term symmetric Lanczos
  recurrence with streamed Givens QR updates for real indefinite and singular
  systems; optional shift, initial guess, SPD Jacobi preconditioning, explicit
  symmetry checking, true-residual stopping, and iteration diagnostics.
- Add an initial `eigsh` subset for symmetric real CSR matrices, implementing
  orthogonal block subspace iteration and Rayleigh--Ritz extraction with
  a small internally diagonalized symmetric projected matrix. Support the
  LA, SA, LM and SM selection criteria, optional starting vector, subspace
  width, convergence budget and actual Ritz residuals. Generalized and
  shift-invert eigenproblems and ARPACK-style implicit restart are not
  implemented; clustered spectra may converge slowly.
- Add native regressions for indefinite/shifted/singular MINRES and symmetric
  sparse eigenvectors, plus randomized SciPy and dense numerical comparisons
  for both checked and optimized builds. The MINRES recurrences carry the
  SciPy BSD-3-Clause attribution; no new build dependency was added.

## 2026-11-04 — sparse LSQR and ILU(0) preconditioning

- Add `scifort_sparse_linalg.lsqr` using original Golub--Kahan LSQR
  recurrences, with rectangular/underdetermined/rank-deficient solutions,
  regularization around optional `x0`, numerical diagnostics, and an
  approximate covariance-diagonal estimate.
- Implement original zero-fill `ilu0` factorization, `ilu0_factor` type,
  `valid_ilu0` and `ilu0_solve` for canonical one-based CSR matrices.
  Allow optional diagonal shifting; retain the input sparsity pattern.
- Extend BiCGSTAB with right ILU(0) preconditioning and restarted GMRES with
  left ILU(0) preconditioning. Preserve existing Jacobi and default APIs.
- Add native edge-case and 1,500-row tests, plus randomized SciPy 1.17 LSQR,
  independent ILU(0) factorization and triangular solve comparisons in
  checked and optimized GNU Fortran builds. No source was imported from
  upstream libraries and no new binary dependency was added.

## 2026-11-04 — sparse CSR iterative linear algebra

- Add `scifort_sparse_linalg` with original dependency-free CG, BiCGSTAB and
  restarted GMRES solvers, one-based CSR matvecs, optional Jacobi inverse
  diagonal, initial guesses, true-residual stopping, and explicit statuses.
- Expose an `lsmr` adapter for general rectangular CSR least squares and
  damped regression, reusing the already-attributed BSD-3-Clause
  `scifort_optimize_lsmr` implementation with CSR and transpose matvecs.
- Add native tests and randomized SciPy 1.17 reference comparisons, plus
  documentation of tolerance, preconditioner, and error semantics. `kiss`
  remains a license-reviewed candidate, not imported into this update.

## 2026-11-04 — canonical CSC, packed diagonals, and triangular CSR utilities

- Extend `scifort_sparse` with canonical, one-based `csc_matrix`, validation,
  COO/dense constructors, bidirectional sparse CSR/CSC conversion, CSC matvec
  and transpose. Retain explicit zeros and sorted unique row indices.
- Add `csr_diagonal`, packed multi-diagonal constructor `csr_diags`,
  in-place whole-diagonal replacement `csr_setdiag`, and sparse upper/lower
  triangular extraction `csr_triu`/`csr_tril`. Support rectangular/empty
  shapes, off-matrix diagonals, and strict input validation.
- Add two native regression programs, randomized SciPy 1.17 cross-comparisons,
  a large sparse conversion test, and new API/attribution documentation.
  No third-party source code or additional runtime dependencies were added.

## 2026-11-03 — generic rectangular CSR sparse matrix arithmetic

- Add independent `scifort_sparse` module with canonical one-based real64
  rectangular CSR matrices and COO/dense constructors and conversions.
- Implement sparse identity/shifted diagonals, transpose, matvec,
  CSR multiplication with a two-pass row accumulator, scaled addition,
  elementwise multiplication, scalar scaling, and zero elimination.
- Add checked/optimized native regressions, including a 12,000-row sparse
  product, and randomized SciPy 1.17 comparison scripts. No additional
  runtime dependencies or third-party source code.

## 2026-11-03 — Yen K-shortest simple paths in sparse graphs

- Add `scifort_sparse_csgraph.yen` to enumerate the K shortest loopless
  directed or undirected source-to-target paths; supports `unweighted`, signed
  costs via Johnson potentials, globally detected negative cycles, zero-weight
  and parallel stored arcs, source=target and unreachable targets.
- Return ordered path costs plus a one-based predecessor row for each path.
  Distinct vertex sequences are emitted only once, independent of parallel
  arc multiplicity. Implement Yen spur generation and a candidate min-heap.
- Add standalone native regression and randomized SciPy 1.17 comparisons,
  plus exhaustive DFS enumeration of all simple routes on small graphs.
  No additional runtime dependencies or third-party source were introduced.

## 2026-11-03 — weighted rectangular matching and sparse CSR Laplacians

- Add `csr_weighted_bipartite_graph`, weighted edge/dense/square CSR
  constructors, validation, and `minimum_weight_full_bipartite_matching`
  (minimize/maximize) with one-based sorted row indices and the
  `csgraph_no_full_matching` status. Includes zero, negative and parallel edges;
  sparse shortest augmenting paths avoid a dense cost matrix.
- Add `laplacian_sparse` with in-/out-degree, normalization,
  symmetrization, parallel-arc aggregation and sparse storage.
- Add native regression tests, SciPy 1.17 randomized comparisons and
  independent small-instance exhaustive matching verification. No new
  dependencies or borrowed third-party code.

## 2026-11-03 — sparse graph reordering, maximum flow and bipartite matching

- Add one-based `reverse_cuthill_mckee` with structural symmetrization,
  duplicate/zero-edge handling, pseudo-peripheral roots, and stable degree sorting.
- Add `maximum_flow` using iterative Dinic with exact integral capacities,
  64-bit flow totals, nonnegative-capacity validation, and a coalesced
  antisymmetric CSR net-flow result.
- Add independent rectangular `csr_bipartite_graph` storage, edge/dense/square
  constructors, iterative Hopcroft–Karp `maximum_bipartite_matching`, and
  maximum-transversal `structural_rank`.
- Add native regressions in both checked and optimized modes, plus randomized
  SciPy 1.17 comparisons with exhaustive small-instance min-cut and matching
  oracles. No external dependencies or borrowed code.

## 2026-11-03 — sparse graph all-pairs paths, Laplacians and search trees

- Add independent one-based `floyd_warshall` and `johnson` shortest-path
  algorithms, with unweighted/undirected options, predecessor matrices,
  signed edges, negative-cycle detection, and selected-source Johnson paths.
- Add the `shortest_path` dispatcher (`auto` / `FW` / `D` / `BF` / `J`), dense
  degree-based and symmetrically normalized `laplacian`, iterative
  `depth_first_order`, and weighted CSR `breadth_first_tree`,
  `depth_first_tree`, and `reconstruct_path` tree construction. Provide
  `path_from_predecessors` as a distinct SciFort route-reconstruction helper.
- Add checked and optimized standalone tests and deterministic randomized
  SciPy 1.17 comparisons. Tests include zero/parallel arcs, invalid inputs,
  disconnected signed graphs, negative cycles, directedness, normalized
  Laplacian variants, selected sources, and tree edges. No new dependencies.

## 2026-11-03 — sparse CSR graph algorithms and hierarchy conveniences

- Add dependency-free `scifort_sparse_csgraph` with one-based compressed sparse
  row constructors/conversion/validation, weak and strong connected components,
  breadth-first traversal, heap-based Dijkstra, Bellman-Ford with reachable
  negative-cycle diagnostics, and Kruskal minimum spanning forests.
- Extend `scifort_cluster_hierarchy` with the seven named linkage convenience
  generics (`single`, `complete`, `average`, `weighted`, `centroid`, `median`,
  `ward`) for observation and condensed-distance inputs.
- Add checked/optimized native graph and wrapper regressions; 280 random
  SciPy 1.17 graph comparisons per mode (1,828 Dijkstra source rows), plus
  405 further graph comparisons per mode (3,026 successful Bellman-Ford
  sources, 84 negative-cycle cases, and 810 BFS traversals).
- Retain SciFort's one-based indexing, all existing hierarchy/VQ interfaces,
  and the project dependency and licensing policy.

## 2026-11-03 — numeric hierarchy dendrograms and Gaussian K-means initialization

- Add `dendrogram` to `scifort_cluster_hierarchy`, returning SciPy-compatible
  numeric branch coordinates and one-based leaf order without graphics;
  implement count and distance sorting, including SciPy's tie rules.
- Add `is_isomorphic` as an O(n log n) two-way partition equivalence checker,
  `is_valid_linkage` / `is_valid_im`, and one-based MATLAB linkage export/import
  (`to_mlab_linkage`, `from_mlab_linkage`).
- Extend `scifort_cluster_vq.kmeans2` with `init='random'`, drawing from a
  Gaussian with empirical mean and unbiased covariance via centered-data
  projection. The construction handles rank-deficient covariance matrices
  without LAPACK or any new dependency; existing defaults are unchanged.
- Add native regression suites `test_cluster_hierarchy_layout` and
  `test_cluster_vq_gaussian`, with checked/optimized build runners, and an
  additional 250-case SciPy 1.17 dendrogram and label comparison per mode.
  SciFort's partition predicate detects several SciPy 1.17 false positives
  when the same right-side label appears in separated left-side groups.

## 2026-11-03 — end-to-end hierarchical clustering and cophenetic correlation

- Implement SciPy-compatible `fclusterdata` for observation matrices, using
  the existing SciFort pairwise metric, linkage and flat-cluster kernels.
- Add `fcluster` as a unified criterion selector for inconsistent, distance,
  maxclust, monocrit and maxclust_monocrit cuts, with optional precomputed
  statistics, real/integer threshold overloads, and one-based outputs.
- Add scale-normalized `cophenetic_correlation` (equivalent to the scalar
  coefficient returned by SciPy's `cophenet(z, Y)`). Return NaN for undefined
  zero-variance correlations and reject nonfinite/mismatched distances.
- Add checked/optimized native regressions and 504 randomized SciPy 1.17
  comparisons per build, spanning seven linkage methods, six pairwise
  metrics and three flat-cut criteria. No third-party source or dependencies.

## 2026-11-03 — globally optimal hierarchy leaf ordering and tree utilities

- Implement `optimal_leaf_ordering` for condensed distances and observation
  matrices; add `optimal_ordering=.true.` to both `linkage` overloads. Minimize
  the adjacent-leaf dissimilarity sum by dynamic programming with one-based
  node IDs, retaining the exact original linkage tree.
- Add indexable `cluster_tree_node` / `to_tree`, `leaders` for validated flat
  subtree assignments with arbitrary integer labels, and size-only `correspond`.
- Add checked and optimized native regression tests, randomized SciPy
  comparisons (190 optimal-ordering and 200 cluster-leader cases per build),
  and independent exhaustive leaf-orientation verification for small trees.
  No external library code or dependencies were incorporated.

## 2026-11-03 — expanded hierarchical-clustering diagnostics and cut criteria

- Extend `scifort_cluster_hierarchy` with SciPy-style `inconsistent` (depth-limited
  sample statistics), `maxdists`, `maxrstat`, `maxinconsts`, `is_monotonic`,
  and `num_obs_linkage`.
- Add `fcluster_inconsistent`, `fcluster_monocrit`, and
  `fcluster_maxclust_monocrit`, with precomputed-statistic support and strict
  monotonicity validation of custom criteria.
- Add full/selected `cut_tree` for exact cluster counts and monotone height
  thresholds. Cluster identifiers remain one-based; height cuts use SciPy's
  strict threshold convention; tied heights and nonmonotone trees are handled
  explicitly.
- Extend native regression tests and add a randomized SciPy 1.17 comparison
  driver (180 scenarios per checked/optimized build). No new dependencies or
  incorporated third-party source code.

## 2026-11-03 — independent vector quantization and K-means clustering

- Add `scifort_cluster_vq` with population-standard-deviation `whiten`,
  one-based nearest-neighbor `vq`, fixed-iteration `kmeans2`, and restart-based
  `kmeans` with SciPy's mean unsquared Euclidean distortion convention.
- Add vector and matrix overloads, explicit SciFort RNG initialization
  (`points` and K-means++), user-supplied initial centroids, configurable
  empty-cluster behavior, finite-input
  validation and scale-normalized Euclidean distances.
- Add native regression tests and randomized comparisons against SciPy 1.17;
  retain existing hierarchy APIs and introduce no third-party algorithm code.

## 2026-11-03 — new hierarchical clustering algorithms

- Add `scifort_cluster_hierarchy` with `linkage` for seven classical methods,
  using both observation matrices and SciPy-compatible condensed distances.
- Add flat cluster thresholds (`distance` and `maxclust`), cophenetic
  distances, and deterministic leaf traversal, with one-based node indices.
- Add native checked/optimized tests and independent randomized SciPy 1.17
  comparisons; no new binary dependencies or imported algorithm source.

## 2026-11-03 — arbitrary-rank rotation and rigid-transform broadcasting

- Add `scifort_spatial_transform_nd` with assumed-rank constructors for SO(3)
  quaternions/rotation vectors/matrices and SE(3) twists/matrices/dual quats.
- Store batch-shape metadata and provide NumPy right-aligned broadcasting for
  composition and point application, including scalar operands; expose
  flattened component outputs that can be reshaped by the caller.
- Provide native regression and independent SciPy 1.17 random reference tests;
  preserve the earlier scalar, 1D, and rank-two APIs without new dependencies.
- Document remaining gaps: reductions and slicing, true dynamic-rank outputs,
  and SVD correction of imperfect rotation matrices.

## 2026-11-03 — Davenport-angle construction and extraction

- Add `Rotation.from_davenport` and `Rotation.as_davenport` equivalents to
  `scifort_spatial_transform`, supporting arbitrary normalized axes with
  perpendicular neighbors, intrinsic/extrinsic conventions, and gimbal locks.
- Add corresponding constructors and accessors for one-dimensional rotation
  stacks and rank-two rotation grids, plus scalar single-axis convenience form.
- Compare against SciPy 1.17 using independent randomized scenarios, including
  badly scaled axes and nonunique equivalent angular representations; add a
  native spatial regression and preserve all existing rotation APIs.
- Uses no additional compiled dependencies. See provenance and license notices.

## 2026-11-03 — rank-two SE(3) rigid-transform grids

- Add `scifort_spatial_rigid_transform_grid` with two-dimensional batches of
  rigid transforms, scalar/stack/grid composition, per-axis broadcasting,
  Cartesian point application, inverse, fractional power and one-based slicing.
- Support matrix, exponential-coordinate and dual-quaternion representations,
  broadcasted translation/rotation construction, weighted global and axis-wise
  chordal rotation and arithmetic translation means.
- Add a checked/optimized native regression and independent randomized SciPy
  1.17 comparisons. Existing scalar and 1D SE(3) APIs are unchanged; no new
  libraries or external numerical source are required.

## 2026-11-02 — native SE(3) rigid transformations and stacks

- Add `scifort_spatial_rigid_transform`: rotation/translation composition,
  inversion, homogeneous 4x4 matrices, point-cloud application, SE(3) exp/log
  coordinates, unit dual quaternions and fractional screw powers.
- Add one-dimensional transform stacks, singleton broadcasting, matrix/twist/
  dual quaternion conversion, weighted chordal rotation/translation means.
- Add native bounds-checked and optimized tests plus independent SciPy 1.17
  comparisons. Native Fortran MIT implementation with no new dependencies.
- Explicitly document differences from SciPy: strict rotation-matrix validation
  and lack of arbitrary-rank broadcasting.

## 2026-11-02 — finite SO(3) rotational symmetry groups and reduction

- Add `scifort_spatial_rotation_group` with `rotation_create_group` for
  SciPy-compatible cyclic, proper dihedral, tetrahedral, octahedral and
  icosahedral groups, including configurable Cn/Dn rotation axes.
- Add scalar and one-dimensional-stack `rotation_reduce`, selecting a
  minimum-angle representative under optional left/right group actions;
  return one-based witness indices and explicit invalid-input statuses.
- Add group closure/inverse tests, native checked/optimized regressions,
  an independent two-seed SciPy comparator and API documentation.
- All code is original MIT-licensed Fortran; no new external dependency.

## 2026-11-02 — reproducible random SO(3) rotations and rank-two rotation grids

- Add `rotation_rng_type`, `rotation_rng`, `rotation_random`, and
  `rotation_stack_random` using Haar-uniform Shoemake quaternions. Support
  caller-owned deterministic streams using overflow-safe 64-bit arithmetic;
  default generation delegates to Fortran `random_number`.
- Add independent `scifort_spatial_rotation_grid` implementing a two-dimensional
  rotation array with quaternion/rotvec/Euler/MRP/matrix conversions, per-axis
  broadcasting of vector application and composition, inversion, angular
  magnitudes and weighted chordal means along either axis.
- Add `test_spatial_rotation_random_grid.f90` and independent SciPy 1.17
  comparison tests, including deterministic random-sequence validation,
  isotropy checks and randomized broadcasting. No new dependencies or copied
  third-party source.

## 2026-11-02 — batched spatial rotations and rotation stacks

- Add independent MIT-licensed `scifort_spatial_rotation_stack`, exposing
  one-dimensional `rotation_stack_type` and batched quaternion, rotation-vector,
  matrix, Euler and MRP constructors and representations.
- Support scalar/stack and singleton-stack broadcasting for composition and
  vector application; add inverse, powers, magnitudes, chordal means,
  one-based indexing/slicing, comparisons, and concatenation.
- Add native regression tests and independent randomized SciPy 1.17.0
  comparisons in checked and optimized builds. Preserve scalar rotation APIs
  and introduce no new runtime dependencies.

## 2026-11-02 — smooth cubic SO(3) rotation splines

- Add `scifort_spatial_rotation_spline` implementing SciPy-style cubic
  `RotationSpline` on nonuniform orientation keyframes, with rotation,
  angular-rate and angular-acceleration evaluation, extrapolation and
  batched queries.
- Implement a small 3x3 block-tridiagonal solver and the SO(3) Jacobian
  transformations; no new runtime dependencies.
- Add native continuity/edge-case regressions and independent randomized
  SciPy 1.17 comparisons, plus API documentation and attribution. Detect
  catastrophic nonlinear-iteration divergence and return explicit status 2
  instead of interpolating from numerically invalid cubic coefficients.

## 2026-11-02 — weighted rotation alignment and chordal means

- Add `rotation_align_vectors(a,b,weights=,rssd=,sensitivity=,status=)`
  with row-major vector collections and a native quaternion Wahba solver.
  Support single-vector minimal rotation, one infinite primary weight,
  arbitrary positive finite weights, and SciPy-compatible sensitivity.
- Add `rotation_mean(rotations,weights=,status=)` for weighted chordal L2 means;
  enforce quaternion sign invariance, reject all-zero weights, and normalize
  across widely different weight magnitudes.
- Validate against SciPy `Rotation.align_vectors` and `Rotation.mean`, with
  native 3D vector/weight/invalid-input tests and independent randomized
  comparison drivers. Preserve all prior source/module exports and dependencies.

## 2026-11-02 — Euler angles and modified Rodrigues parameters

- Extend `scifort_spatial_transform` with `rotation_from_euler` (1-3 axes),
  `r%as_euler` (24 three-axis intrinsic/extrinsic sequences), and gimbal-lock
  detection with SciPy's third-angle-zero convention.
- Add `rotation_from_mrp` and `r%as_mrp` with scale-stable arithmetic and
  the principal/shadow representation semantics.
- Add dedicated native tests and the independent SciPy comparison script
  `compare_rotation_euler_mrp_scipy.py`, covering singularities, proper Euler
  angles, degrees, and extreme-magnitude inputs. No imported code or runtime
  dependencies were added.

## 2026-11-02 — spatial.transform rotations and piecewise Slerp

- Add the independent `scifort_spatial_transform` module with `rotation_type`
  and conversions to/from quaternions, rotation vectors and 3x3 matrices.
- Support 3D vector application, inverse, composition, principal rotation
  magnitude and scalar-first/canonical quaternion representations.
- Add two-rotation shortest-path `rotation_slerp` and multi-keyframe
  `rotation_slerp_path` with nonuniform knots and explicit validation.
- Add two native regression programs and randomized SciPy comparison tests
  using checked and optimized GNU Fortran. No new dependencies or imports.

## 2026-11-02 — validated 3D convex-hull reconstruction

- Fix a reproducible near-coplanar 26-point 3D hull that previously failed
  during neighbor assembly. Reject incomplete/nonconvex output by validating
  triangulation topology and all input points against each supporting plane.
- Recover difficult configurations using deterministic insertion-order retries,
  followed by a supporting-plane enumeration/face-triangulation fallback.
  Preserve public indexing, volume/area semantics, and the native MIT license.
- Add `test_spatial_convex_hull_nearcoplanar.f90`, a 26-point SciPy oracle
  fixture, and a separate 100-case nearly-coplanar reference test covering
  face perturbations from zero through `1e-8`.
- Update the native spatial test runner and 3D hull documentation. No new
  external runtime dependencies or third-party source code.

## 2026-11-02 — native 2D/3D bounded halfspace intersections

- Add `scifort_spatial_halfspace` and export `halfspace_intersection` and
  `halfspace_intersection_type` through `scifort_spatial`.
- Find distinct vertices of bounded 2D/3D intersections of affine inequalities,
  with strict-interior and recession checks, scaled coordinates, and redundant
  constraint handling. Also expose polar points, simplicial facets, supporting
  equations, area, and volume.
- Enumerate plane intersections directly to avoid missing vertices when 3D
  incremental hull triangulation is numerically incomplete. Use existing 2D
  hull geometry to triangulate non-simplicial polar faces.
- Add native status/geometry tests, independent randomized SciPy comparisons
  with two seeds, and documentation. No new dependency or imported code.

## 2026-11-02 — native 2D/3D convex hull construction

- Add `scifort_spatial_convex_hull` and export `convex_hull` and
  `convex_hull_type` through the existing `scifort_spatial` facade.
- Support 2D polygons and 3D polyhedra with outward supporting planes,
  one-based extreme point indices, simplicial facets, and facet adjacency.
- Return SciPy-compatible perimeter/surface area and enclosed area/volume;
  handle interior/duplicate points, uniform scaling, translation, and
  invalid lower-dimensional geometry without external dependencies.
- Add a native regression and independent SciPy numerical tests with
  randomized data and multiple seeds in checked and optimized builds.

## 2026-11-02 — native 2D/3D spherical Voronoi diagrams

- Add `scifort_spatial_voronoi` and export the `spherical_voronoi` factory
  and `spherical_voronoi_type` through `scifort_spatial`.
- Implement independently generated 3D triangular convex hull facets, cyclic
  region-vertex connectivity, 2D circle hull arcs, and spherical polygon areas.
- Preserve SciPy's 2D hull-edge `calculate_areas()` behavior; introduce
  `region_measures()` for per-generator measures in both dimensions.
- Add native `test_spatial_voronoi`, randomized independent SciPy comparisons,
  API documentation, and no additional build/runtime dependencies.

## 2026-11-02 -- Procrustes similarity alignment and orthogonal fits

- Add `src/spatial_procrustes.f90` with independent one-sided Jacobi polar
  decomposition, normalized/centered `procrustes`, and uncentered
  `orthogonal_procrustes` for real binary64 data.
- Support rotation/reflection/scaling, deterministic rank-deficient
  orthogonal completion, stable input centering, zero singular values,
  and explicit invalid-input/nonconvergence status codes.
- Export through `scifort_spatial`, add `test_spatial_procrustes.f90`,
  randomized SciPy comparison driver, and `docs/spatial_procrustes.md`.
- No new build dependency or imported third-party implementation code.

## 2026-11-02 — directed Hausdorff distances and great-circle interpolation

- Add MIT-licensed `scifort_spatial_hausdorff` with exact
  `directed_hausdorff` (optional one-based witness indices and early
  abandonment) and symmetric `hausdorff_distance` (SciFort extension).
- Add MIT-licensed `scifort_spatial_slerp` with rank-overloaded
  `geometric_slerp` for scalar/vector path fractions in arbitrary dimension,
  stable short-arc interpolation and extrapolation.
- Use deterministic witness ordering, explicit error statuses, scaled norms,
  and documented differences from SciPy for empty clouds/antipodal arcs.
- Export the new APIs through `scifort_spatial`, add two native regression
  tests, SciPy comparison drivers, and `docs/spatial_geometry.md`.
- No new dependency or third-party implementation source was introduced.

## 2026-11-02 — five weighted Boolean `scipy.spatial.distance` metrics

- Extend `scifort_spatial` from 14 to 19 distance metrics: add `dice`,
  `rogerstanimoto`, `russellrao`, `sokalsneath`, and `yule` to `point_distance`,
  `pdist`, and `cdist` without changing existing call signatures.
- Use independently implemented Boolean contingency counting, optional
  nonnegative per-feature weights, and weight normalization that supports
  extreme finite weight magnitudes without intermediate overflow.
- Define edge behavior for zero-weight observations and undefined ratios;
  interpret nonzero numeric samples consistently as Boolean true.
- Add `test_spatial_distance_boolean.f90`, randomized independent SciPy
  comparisons, and documentation. No new dependency or third-party source.

## 2026-11-02 — standardized Euclidean, Mahalanobis and Jensen--Shannon distances

- Extend the MIT-licensed `scifort_spatial` distance API to 14 metrics, adding
  `seuclidean`, `mahalanobis`, and `jensenshannon` to `point_distance`, `pdist`,
  and `cdist` without modifying existing positional arguments.
- Accept explicit per-coordinate `variance` and matrix `inv_cov`; automatically
  estimate sample statistics for `pdist` and `cdist` when omitted, using
  Cholesky factorization of a scale-normalized covariance matrix.
- Support normalized nonnegative probability vectors and optional base>1 for
  Jensen--Shannon. Explicitly reject zero variance, singular covariance,
  invalid distributions, and mismatched parameters with status codes.
- Add dedicated native tests, randomized independent SciPy comparisons with
  two seeds in both GNU Fortran configurations, and updated API documentation.
  No third-party code or new runtime dependencies were introduced.

## 2026-11-01 — weighted and multi-radius `scipy.spatial.cKDTree.count_neighbors`

- Overload the MIT-licensed native `kd_tree_type%count_neighbors` with
  rank-one radius arrays and `cumulative=.false.` distance histogram bins.
  Preserve existing scalar, int64-returning invocations.
- Add scalar/array `weighted_count_neighbors`, with optional finite source
  and destination weights (including zero and signed values) and default
  unit weights for omitted arrays.
- Count all radius thresholds in one dual-tree traversal, using subtree
  weight sums, bounding-box pruning, sorted-radius binary searches, and
  range updates to accumulate counts for fully included subtrees.
- Compare both APIs with SciPy 1.17.0 across 400 randomized cases in checked
  and optimized GNU Fortran builds, and extend native regression tests.
- No new dependency, third-party source, or Delaunay changes.

## 2026-11-01 — periodic and Minkowski KD-tree spatial searches

- Extend the existing balanced, native Fortran KD-tree to support scalar or
  per-coordinate `boxsize`, including nonperiodic axes (`L=0`) and wrapped
  coordinates in nearest, radius, pair, and cross-tree searches.
- Add exact Minkowski orders `1<=p<=infinity` to all existing tree searches,
  including dual-tree pair counting and one-based sparse distance triplets.
  Keep the legacy `distance2` output and expose optional unsquared `distance`.
- Implement conservative periodic minimum/maximum bounding-box distances,
  reusable scaled p-norms, and input/topology validation with error statuses.
- Add randomized native and independent SciPy reference tests using checked
  and optimized GNU Fortran. No new dependencies or copied upstream source.

## 2026-11-01 — KD-tree pair and cross-tree spatial searches

- Extend the existing balanced KD-tree with four `scipy.spatial.cKDTree`
  methods: `query_pairs`, `query_ball_tree`, `count_neighbors`, and
  `sparse_distance_matrix`. Export the shared type from `scifort_spatial`.
- Add one-based sorted pair results, CSR neighbor lists, ordered `int64` pair
  counts, and COO distance triplets, including duplicates and exact zeros.
- Count dense subtree pairs in bulk with a dual-tree bbox traversal; use
  a scaled Euclidean norm for extreme-magnitude closed-ball comparisons.
- Add independent SciPy comparisons, native regressions, API documentation,
  and spatial test runners. Existing interpolation routines and external
  dependencies remain unchanged.

## 2026-11-01 — `scipy.spatial.distance` pairwise metric coverage

- Add a new native `scifort_spatial` module exporting `point_distance`,
  `pdist`, `cdist`, and rank-overloaded `squareform` for real binary64 data.
- Support 11 metrics (Euclidean, squared Euclidean, cityblock, Chebyshev,
  Minkowski, cosine, correlation, Hamming, Jaccard, Canberra and Bray--Curtis)
  with optional nonnegative per-coordinate weights and SciPy-compatible
  condensed pair ordering.
- Use scale-normalized Euclidean/Minkowski and cosine/correlation kernels,
  explicit error statuses, degenerate-distance NaN/Inf handling, and strict
  symmetry/diagonal validation for the matrix-to-vector conversion.
- Add native edge-case tests and 495 randomized SciPy comparison cases per
  seed in both checked and optimized GNU Fortran builds. No new Fortran
  dependencies or imported upstream source code.

## 2026-11-01 — `make_lsq_spline` with complete knot vectors

- Add `make_lsq_spline` to the public `scifort_interpolate` API for real
  scalar and multichannel weighted least-squares fitting with explicit full
  knots, including nonuniform/repeated interior knots and degrees 0–15.
- Support nondecreasing sample coordinates and zero weights; expose `qr`
  method selection and report rank deficiency via the existing status codes.
- Extract the existing banded Givens QR algorithm into a shared kernel used
  by both `make_lsq_spline` and `lsq_univariate_spline`, without changing
  the FITPACK-style input restrictions of the latter.
- Add independent SciPy coefficient comparisons (480 seeded cases), fixed
  SciPy 1.17 numerical fixtures, invalid-input tests, and documentation.
  Checked and optimized GNU Fortran builds pass; no new dependencies.

## 2026-11-01 — periodic interpolating B-splines (scalar and multichannel)

- Implement `make_interp_spline(x,y,k=...,bc_type='periodic')` for real data,
  with cyclic knots, value interpolation and matching endpoint derivatives
  for polynomial degrees 2–15; preserve SciPy's degree-0/1 special cases.
- Factor a single periodic collocation system across output channels, using
  SciFort's existing dense pivoted solver; enable periodic evaluation,
  derivatives, knot insertion and exact integration across multiple periods.
- Validate against SciPy 1.17.0: 150/150 and 120/120 randomized cases per
  checked/optimized build, covering degrees 1–15, coefficients, derivatives,
  and integrals. All 30 native interpolation tests pass in both builds.
  Normalize high-degree derivative rows to avoid false singularity flags.
  No mandatory dependency.
- Record the periodic-knot construction's SciPy BSD-3-Clause provenance;
  leave Delaunay and the optional FITPACK backend unchanged.

## 2026-11-01 — bivariate spline partial derivatives and exact knot refinement

- Add `fitpack_bivariate_tck_type%partial_derivative(dx,dy,status)` to construct
  the mixed-partial derivative spline (SciPy `BivariateSpline.partial_derivative`).
  Handle asymmetric x/y degrees, repeated knots, and padded FITPACK coefficients.
- Add `%insert_knot(axis,x,m,status)` with exact tensor-product Boehm insertion
  in x or y; expose free-function aliases via `scifort_interpolate`.
- Validate against SciPy bivariate/BSpline operators across 180 randomized
  cases in checked and optimized GNU Fortran builds, including derived
  coefficients, repeated knots, values, derivatives and integrals; all pass.
- Expand the native interpolation suite to 29 programs; keep all 29 passing.
  No new dependencies, vendoring, or Delaunay changes.

## 2026-11-01 — exact rectangular-grid bivariate spline interpolation and integration

- Add dependency-free `rect_bivariate_spline(x,y,z,kx,ky,status)`, covering
  `scipy.interpolate.RectBivariateSpline(s=0)` on strictly increasing real grids
  with independent spline degrees 1–5. Reuse SciFort's existing 1D not-a-knot
  collocation solver; return the established `fitpack_bivariate_tck_type`.
- Add exact signed, domain-clipped rectangular integration for **any** valid
  bivariate TCK through `%integral()` and `bispline_integral()`, including
  optional-backend FITPACK surfaces; introduce `bispline_singular` status.
- Compare 120 randomized rectangular fits per checked and `-O3` build with
  independent SciPy 1.17.0 knots, coefficients, values/derivatives, and
  integrals; all passed with maximum normalized evaluation error 3.18e-13.
- Add native analytic, validation and all 25 degree-pair regression tests.
  All 28 interpolation test programs pass in the bounds-checked build.
- No copied upstream code, added dependencies, or Delaunay changes.

## 2026-11-01 — bivariate FITPACK smoothing and fixed-knot least squares

- Add dependency-free `fitpack_bivariate_tck_type`, `bisplev` scalar/grid
  evaluation, mixed partial derivatives, and paired-query `bisplev_points`.
- Add optional `perazz/fitpack` 2.0.0 adapters for scattered-surface automatic
  smoothing and fixed-knot least squares (including asymmetric spline degrees).
- Fix degree-dependent FITPACK surface workspace and preserve user-provided
  interior knots by requesting one actual fixed-knot fitting pass.
- Verify the genuine upstream bivariate adapter, independent SciPy TCK
  evaluation and residual reconstruction, and all 27 native interpolation tests.
- Preserve the optional dependency architecture; no FITPACK code is vendored.
  Delaunay interpolation is unchanged.

## 2026-10-31 — FITPACK-style weighted LSQ and exact univariate splines

- Implement `lsq_univariate_spline` (SciPy `LSQUnivariateSpline`): weighted
  least squares with supplied interior knots, degrees 1–5, optional bbox,
  and scalar or multichannel real observations.
- Implement `interpolated_univariate_spline` (`InterpolatedUnivariateSpline`)
  using the existing not-a-knot B-spline interpolator.
- Use a banded **Givens QR** solver rather than dense normal equations; shared
  factorization across output channels, no external BLAS/LAPACK dependency.
- Add native rank/geometry tests, large-data fitting test, and independent
  SciPy comparisons of knots, coefficients, derivatives and integration.
- Preserve Delaunay implementation unchanged.

## 2026-10-31 — natural cubic smoothing splines (new algorithm coverage)

- Stop expanding native Delaunay dimensions; existing 2D–5D triangulation
  and supplied higher-dimensional simplices remain available unchanged.
- Implement `scipy.interpolate.make_smoothing_spline`-style fixed-penalty
  and GCV-selected weighted natural cubic smoothing for one or many real
  output channels, returning an existing `PPoly` object.
- Add O(n) pentadiagonal Cholesky fits and O(n) exact GCV trace calculation
  through a selected-inverse recursion; no new BLAS/LAPACK dependencies.
- Add independent SciPy numerical comparisons (fixed lambda, GCV, values,
  derivatives, coefficients and exact integrals), native regression, and
  detailed API/algorithm documentation.

## 2026-10-31 — dimension-generic native 4D/5D Delaunay and fast simplex location

- Extend `LinearNDInterpolator` and `griddata(method='linear')` to native 4D
  and 5D Delaunay construction with one rank-two geometric implementation:
  greedy affine seed, extended-precision circumsphere calculations, cavity
  cancellation, and visible-hull preservation. Retain 2D/3D specializations.
- Preserve externally supplied simplices through ten or more dimensions and
  limit native 4D/5D builds to at most 160 points until robust adaptive
  predicates and scalable high-dimensional meshing are available.
- Add cached barycentric affine transforms, simplex-facet adjacency, and
  nearest-site-guided neighbor walking; retain a generalized BVH fallback.
- In reproducible `-O3` query benchmarks versus the predecessor ZIP, reduce
  4D and 5D evaluation CPU times by approximately 6.7× and 27.8×,
  respectively, using identical supplied meshes.
- Add native 4D/5D regression, independent SciPy/Qhull geometry and complex
  value comparisons, a predecessor query benchmark, and updated documentation.
- No extra FPM/BLAS/LAPACK dependencies or copied upstream source.

## 2026-10-31 — accelerated 3D Delaunay and dimension-generic supplied-simplices interpolation

- Cache conservative circumsphere centers/radii for 3D tetrahedra, screening
  non-cavity tetrahedra before exact extended-precision in-sphere predicates.
  On a reproducible GNU Fortran `-O3` benchmark this reduced median 3D
  construction time by 6.6–12.5× for 160–440 sites.
- Add `simplices=` to the existing real/complex `linear_nd_interpolator`
  constructors. Supplied `(ndim+1,nsimplices)` tessellations support
  dimension-generic barycentric interpolation, validated through ten dimensions.
  Native Delaunay construction remains 2D/3D only.
- Preserve real and complex multichannel values, `rescale`, fill values,
  existing 2D/3D constructors, and the tetrahedron BVH path.
- Add native and SciPy randomized comparisons for externally supplied
  2D–7D simplices, source-based predecessor benchmark, and documentation.
- No copied third-party code, added FPM dependencies, BLAS or LAPACK.

## 2026-10-31 — accelerated 3D simplex location with tetrahedron BVH

- Build a median-split axis-aligned bounding-volume hierarchy over the
  independently generated 3D Delaunay tetrahedra, with per-tetrahedron AABBs.
- Replace the former full tetrahedron scan in 3D `LinearNDInterpolator` and
  `griddata(method='linear')` with indexed candidate pruning; retain the
  existing high-precision barycentric calculation, deterministic lowest-ID
  selection, fill semantics, and public APIs.  The 2D path is unchanged.
- Add a native near-face/near-edge/complex regression, a reproducible
  predecessor-versus-indexed benchmark, and updated documentation.
- Verify 3D SciPy comparisons, 2D regressions, independent centered-Qhull
  tetrahedralization checks, and the full native interpolation test set
  in GNU Fortran checked and optimized builds.
- No additional FPM dependencies, system libraries, or third-party source.

## 2026-10-31 — 3D Delaunay tetrahedralization and LinearND/griddata support

- Extend real/complex, scalar/multichannel `LinearNDInterpolator` and
  `griddata(method='linear')` from 2D to 3D without changing public calls.
- Implement an independent hull-preserving incremental 3D Bowyer–Watson
  tetrahedralization with extended-precision in-sphere/orientation tests,
  hashed hull facets, and tetrahedral barycentric evaluation.
- Preserve the 2D triangulator and Clough–Tocher APIs and dependencies.
- Add native affine, boundary, complex, fill, and degenerate-input tests,
  plus independent SciPy and Qhull geometry comparisons.
- Explicitly document Qhull numerical-offset sensitivity, native predicate
  limitations, and the absence of new BLAS/LAPACK or FPM dependencies.

## 2026-10-31 — indexed, cached Clough–Tocher evaluation

- Cache the 19 active cubic Bernstein–Bézier controls per triangle and output
  channel during construction, avoiding repeated control reconstruction.
- Introduce a bounded-memory CSR uniform-bin index over triangle bounding
  boxes to accelerate interior queries, with exact linear-search fallback.
- Improve shared-edge containment for anisotropic, large-offset coordinates
  by accounting for normalization roundoff in the geometric tolerance.
- Add a 210-site, 403-triangle native boundary/affine/complex regression, a
  reproducible old-versus-new benchmark, and documentation updates.
- Validate checked and optimized GNU Fortran builds against SciPy; retain
  existing public API, license notices, and FPM dependencies.

## 2026-10-31 — Clough–Tocher C1 scattered-data cubic interpolation

- Add SciPy-style `CloughTocher2DInterpolator` for real/complex multichannel
  scattered data, configurable gradient-estimation tolerance/iterations,
  coordinate rescaling, and out-of-convex-hull fills.
- Implement 2-D global minimum-curvature gradient estimation and affine-
  invariant cubic triangular Bézier patches with C1 edge continuity.
  Adapt the BSD-3-Clause SciPy numerical formulas with required attribution.
- Extend the public `scifort_interpolate:griddata` generic to accept
  `method='cubic'`; preserve existing nearest/linear signatures.
- Improve Bowyer–Watson triangulation using extended-precision in-circle
  predicates (with portable double fallback), conservative convex-hull notch
  repair, local Delaunay edge flips, and a final hull-area consistency check.
  Reject unresolved incomplete meshes instead of returning silent gaps.
- Add randomized SciPy/Qhull comparisons, native C1/affine/complex/error tests,
  and updated stand-alone Fortran comparison scripts. No new required linker
  library or FPM dependency.

## 2026-10-31 — balanced KD-tree for scattered interpolation and RBF fitting

- Add independent `scifort_interpolate_kdtree` with reusable balanced KD-tree,
  deterministic duplicate/tie handling, exact k-nearest and radius queries,
  subtree bounding-box pruning, and visit-count instrumentation.
- Integrate the index with `NearestNDInterpolator` and local
  `RBFInterpolator`, replacing repeated full-site distance scans; preserve
  rescaling and one-based public Fortran indexing semantics.
- Add native brute-force/tie/error/performance regressions, randomized
  `scipy.spatial.cKDTree` comparisons, and a reproducible nearest-neighbor
  benchmark. Update standalone interpolation reference build dependencies.
- No new upstream source or FPM, BLAS/LAPACK, or Qhull dependency; retain
  the existing 2D-only limit on Delaunay-based linear interpolation.

## 2026-10-31 — scattered-data Delaunay and nearest interpolation

- Add independent `scifort_interpolate_scattered` implementing
  `NearestNDInterpolator` in any dimension, `LinearNDInterpolator` in 2D,
  and a `griddata` convenience interface for nearest and 2D linear modes.
- Implement 2D Bowyer–Watson Delaunay triangulation and barycentric
  evaluation without external Qhull, BLAS/LAPACK, or FPM dependencies.
  Support real/complex multichannel values, rescaling, configurable
  convex-hull fills, vector/batch evaluation and triangle inspection.
- Replace the `RBFInterpolator` local-site selection's repeated scans with
  an `O(n log k)` maximum heap, preserving stable nearest-first output.
- Add checked/optimized native tests and SciPy 1.17 randomized comparisons,
  including sample-site queries, convex-hull fills and invalid input tests.
  Document the 2D-only linear and non-robust-predicate limitations.

## 2026-10-31 — RBFInterpolator and legacy interp1d

- Add independent `scifort_interpolate_rbf` for real/complex multichannel
  scattered-data interpolation, eight SciPy 1.17 radial kernels, optional
  Matérn kernels, global/local neighbor fitting, polynomial augmentation,
  and scalar or per-node smoothing.
- Add independent `scifort_interpolate_legacy` for `interp1d` with nine
  interpolation kinds, sorting, asymmetric out-of-range fill, extrapolation,
  and real/complex channels; reuse native B-splines for higher-order kinds.
- Export both APIs through `scifort_interpolate`. Add two native test
  programs, two SciPy reference drivers, and complete documentation.
- Keep FPM dependencies unchanged; no third-party code copied or vendored.

## 2026-10-31 — Krogh and Floater–Hormann interpolation

- Add independent `scifort_interpolate_global`, exported by the existing
  `scifort_interpolate` facade. Implement SciPy-style confluent Newton/Hermite
  (Krogh) interpolation, including repeated-node derivatives, real/complex
  multichannel outputs, and all-derivatives-array evaluation.
- Implement Floater–Hormann pole-free rational interpolation for monotone
  real nodes and real/complex multichannel ordinates, using normalized
  log-space barycentric weights and local degree `d`.
- Add checked/optimized native regressions, two-seed SciPy comparisons,
  public API documentation, and a direct offline all-native-test runner.
- Preserve the optional, nonvendored `polyroots-fortran` integration without
  introducing any new FPM or system BLAS/LAPACK dependencies.

## 2026-10-31 — genuine upstream polyroots adapter compilation

- Compile Jacob Williams's user-supplied `polyroots-fortran` master source
  against the existing optional SciFort callback, rather than a test double.
- Verify checked/optimized builds, conventional LAPACK linking, and
  GNU linker garbage collection without BLAS/LAPACK for the `rpoly` path.
- Add 120 deterministic degree-5–25 polynomial test cases per build and
  rerun the genuine upstream-versus-native benchmark at degrees 5–17.
- Add an offline build-and-test script that accepts a separate upstream ZIP
  and does not vendor any third-party source or modify main dependencies.
- Document the upstream nonstandard `COMPLEX*16` declaration, distinct
  master-versus-tag provenance, and outstanding Windows/FPM limitations.

## 2026-10-31 — barycentric interpolation and optional polynomial-root callbacks

- Add SciPy-style `BarycentricInterpolator`/`barycentric_interpolate`
  for real and complex scalar/multichannel ordinates with stable normalized
  weights, derivatives, value updates, and incremental node insertion.
- Add `root_solver=` callback to `PPoly.solve` and `PPoly.roots` without
  changing the default dependency-free algorithm. Export the native
  per-segment real-root routine for optional solver adapters.
- Add a version-pinned, opt-in `polyroots-fortran` 1.2.0 adapter using
  its Jenkins--Traub `rpoly` routine. The main FPM manifest remains
  unchanged; upstream package linking requires system BLAS/LAPACK.
- Add checked/optimized native and independent SciPy reference drivers.
  Document the SciPy 1.17 even-batch barycentric insertion discrepancy
  and distinguish interface-stub checks from real upstream validation.

## 2026-10-31 — piecewise extension, Hermite order selection and root-finding reuse review

- Add `%extend()` to real/complex `PPoly` and `BPoly`, appending/prepending
  intervals while preserving polynomial values across different local degrees.
  Use safe, transactional validation and Bernstein degree elevation.
- Add optional per-interval `orders(:)` and uniform scalar `order` to
  `bpoly_from_derivatives`, with explicit validation and asymmetric endpoint
  derivative selection.
- Add checked/optimized native tests and randomized comparisons against SciPy
  and independently solved Hermite systems; document upstream SciPy edge cases.
- Review `polyroots-fortran`, `roots-fortran`, and other catalog root-finders,
  their prospective license obligations and optional FPM backend strategy.
  No source was copied and no dependency was introduced.

## 2026-10-31 — PPoly real roots, from_spline and complex polynomial calculus

- Add `PPoly.roots` and `PPoly.solve` for real coefficients: stationary-point
  isolation, tangential roots, exact-zero sections with NaN sentinels,
  discontinuity crossing, multiple real channels, extrapolation, and
  descending breakpoints. No companion matrix or system LAPACK is required.
- Add `ppoly_from_spline` for scalar and multichannel SciFort BSpline objects,
  including repeated knot vectors and Taylor coefficients through derivative
  evaluation.
- Add `complex_ppoly_type` and `complex_bpoly_type` through overloads of
  `ppoly()` and `bpoly()`, providing complex evaluation, derivatives,
  antiderivatives, exact integration, and power/Bernstein basis conversion.
- Add independent SciPy 1.17 randomized Fortran/Python comparison drivers and
  native regressions. Retain existing Fortran compiler and FPM dependencies.

## 2026-10-31 — PPoly, BPoly, and arbitrary-order Hermite interpolation

- Implement independent real-valued `scifort_interpolate_piecewise` with
  `PPoly` local power and `BPoly` Bernstein coefficients, scalar and
  multichannel constructors, derivative and antiderivative objects,
  exact definite integration, periodic extrapolation and basis conversions.
- Add `bpoly_from_derivatives`, supporting unequal numbers of derivatives at
  consecutive knots and Bernstein degree elevation to a shared degree.
- Add `ppoly_from_cubic` to convert existing scalar or multichannel PCHIP,
  Akima, Hermite, and CubicSpline objects into a general `PPoly`.
- Add independent SciPy 1.17 oracle drivers, native Fortran regression tests,
  documentation and examples. Update existing standalone interpolator build
  scripts for the new public-facade dependency. No new FPM dependency.

## 2026-10-31 — multidimensional PCHIP and banded tensor spline collocation

- Add `RegularGridInterpolator(method='pchip')` and `interpn(...,method='pchip')`
  for real 1–6D grids, multichannel ordinates, ascending/descending axes,
  fill/bounds/extrapolation and per-call method selection. Recursion processes
  axes last-to-first as in SciPy. Match SciPy's >=4-knot rule, reject complex
  PCHIP and disallow `nu=` for PCHIP.
- Replace dense per-axis B-spline collocation storage and pivoted LU with a
  specialized banded no-pivot LU (degree 1/3/5), reducing factor storage from
  `O(n²)` to `O(n*k)` and factorization from `O(n³)` to `O(n*k²)`.
- Add randomized SciPy reference tests, native PCHIP regression and a
  reproducible old/new fit-time benchmark; update provenance and API docs.

## 2026-10-31: Tensor-product B-spline interpolation and mixed partial derivatives

- Extend `scifort_interpolate_regular_grid` with `'slinear'`, `'cubic'`, and
  `'quintic'` methods for 1–6D rectilinear grids, in real or complex, with
  optional trailing output channels. Use a separable per-axis not-a-knot
  collocation fit and local B-spline basis evaluation, including descending
  coordinate axes and extrapolation.
- Add mixed partial derivative orders via optional `nu` to reusable real and
  complex interpolation objects; support the same method override and
  fill/bounds policies as the earlier nearest/linear methods.
- Compare 189 cases and 2,997 results per seed/compiler build against SciPy
  using its sparse direct solver, covering two seeds and checked/optimized GNU
  Fortran. Add exact-polynomial and mixed partial derivative native regression
  including six-dimensional cubic and quintic grids; rerun legacy interpolation
  suites. No upstream code copied or new FPM dependency.

## 2026-10-31: Higher-order BSpline and make_interp_spline

- Add `scifort_interpolate_bspline`, re-exported by `scifort_interpolate`,
  with explicit-knot B-spline construction, generalized de Boor evaluation,
  derivatives and derivative objects, repeated antiderivatives, exact definite
  integrals, shape-preserving knot insertion and dense basis matrices.
- Add `make_interp_spline` with SciPy-compatible default not-a-knot placement
  for degrees 0–15, and natural/clamped cubic boundary conditions. Support
  real-valued multichannel `axis=0` ordinates and optional custom knot arrays.
- Add checked/optimized native/public tests and two-seed independent SciPy
  comparisons (1,390 cases per seed and per compiler build, including
  nonclamped and repeated-knot cases), preserving the
  previous PCHIP/Akima/CubicSpline/RegularGridInterpolator test suites.
- Document numerical scope, dense-solver limitations, degree-zero endpoint
  semantics and third-party reuse review. No new FPM dependency; no third-party
  source copied.

## 2026-10-31: Rectilinear RegularGridInterpolator / interpn support

- Add `scifort_interpolate_regular_grid` re-exported by `scifort_interpolate`,
  implementing SciPy-style 1–6D nearest and multilinear interpolation over
  uneven, increasing/decreasing and singleton coordinate axes. Support real
  and complex grids, trailing output channels, reusable interpolation plans,
  `interpn` / `complex_interpn` convenience wrappers, method overrides,
  fill values, bounds errors and edge extrapolation.
- Add an independent SciPy/NumPy reference driver with 540 randomized cases
  per seed across checked and optimized GNU Fortran builds (including two
  seeds), plus a native API/invalid-input regression. Retain and rerun the
  existing PCHIP/Akima/Makima/Hermite and scalar/multichannel CubicSpline
  suites with the new module in their build dependency closure.
- Document public Fortran interfaces, shape conventions, status codes and
  lack of tensor-product cubic/quintic/PCHIP spline functionality. This is
  independently written MIT code, not a copied adaptation of `finterp`.
  No new FPM or BLAS/LAPACK dependency was introduced.

## 2026-10-31: CubicSpline, matrix-valued ordinates, upstream interpolation review

- Add the `cubic_spline` generic to `scifort_interpolate`, implementing exact
  piecewise-cubic interpolation with not-a-knot, natural, clamped, periodic and
  independent left/right first- or second-derivative end conditions. Include
  SciPy's 2- and 3-knot special cases and periodic extrapolation/integration.
- Add rank-two real ordinate arrays with first dimension indexing knots and
  second indexing output channels. `piecewise_cubic_array_interpolator`
  implements channel-wise scalar/vector evaluation, derivatives, integrals and
  coefficient retrieval without changing the existing 1D interface.
- Add native API, invalid-input, and periodic tests. Compare 847 randomized
  scalar and 480 randomized multichannel cases (412,055 numerical comparisons)
  in each of two independent SciPy reference builds (checked and optimized).
- Review the interpolation projects in Fortran-code-on-GitHub for compatible
  code reuse. Document preferred BSD-style FPM dependencies and warnings about
  ACM, noncommercial, and mixed-GPL licensing. No upstream code was copied,
  and no new runtime dependency was introduced.

## 2026-10-31: PCHIP, Akima/Makima and cubic Hermite interpolation

- Add a dependency-free `scifort_interpolate` module with `pchip_interpolator`,
  `akima1d_interpolator(method='akima'/'makima')`, `cubic_hermite_spline`,
  and the `pchip_interpolate` scalar/vector convenience generic.
- Expose reusable piecewise-cubic objects with scalar/vector evaluation,
  derivatives of arbitrary nonnegative order, exact definite integration,
  native coefficients/knots, explicit extrapolation policies, and validated
  status reporting. Preserve SciPy's default extrapolation choices.
- Validate checked/optimized GNU Fortran builds against 480 independent SciPy
  PCHIP/Akima/Makima/Hermite scenarios and 149,400 numerical values per build,
  covering knots, extrapolation, flat segments, derivatives, coefficients and
  exact integrals. Add native API/invalid-input regressions.
- Document real univariate scope and absence of full SciPy interpolation
  classes, periodic extension, ND ordinates, or root solvers. No additional
  production dependency or C ABI changes.

## 2026-10-31: SciPy band-stop objective and centered exponential windows

- Add public `band_stop_obj` for continuous-order Butterworth, Chebyshev and elliptic analog band-stop objectives. Preserve SciPy's zero-based edge selector. Use stable dB log ratios and complementary elliptic integrals for very small ripple and high attenuation; return a quiet NaN and status code for invalid inputs.
- Expand the `exponential`/`poisson` window family to arbitrary periodic window centers with `parameter2=center` alongside the existing `parameter=tau`, without changing the centered tau-only behavior.
- Check both routines with GNU Fortran `-fcheck=all` and `-O3`, including 1,102 independent stable order references, 922 direct SciPy band-stop comparisons and 400 SciPy exponential window comparisons per build. Rebuild 107 public signal dependency modules and run nine focused public/legacy regressions. No new Fortran runtime dependencies.
- Add native and independent SciPy regression drivers, documentation and provenance in `docs/signal_band_stop_objective.md` and `docs/signal_exponential_center.md`.

## 2026-10-31: Native exact-band fixed-taper spectrogram and DPSS coverage

- Add `fixed_taper_spectral_matrix_spectrogram` (real/complex overloads) for deterministic weighted taper banks, exact analytic continuous-frequency-band cross-spectral matrices, positive/zero weights, per-taper energy normalization, configurable window overlap, and explicit Fourier-phase convention. Return the existing spectral-matrix spectrogram type, so it works directly with canonical moving-block bootstrap.
- Create a native VAR(p) conditional coherence coverage driver that accepts raw Gaussian VAR(2)/VAR(3) time samples and SciPy DPSS taper banks and performs both spectral estimation and all five existing bootstrap methods in Fortran. Compare full-sample estimates with independent NumPy/SciPy finite-window spectra and exact matched stationary population targets, varying DPSS bandwidth/taper count/weighting and block length on paired datasets. Provide exploratory 1,620-row CSV (8 datasets/model, 59 draws) without claiming calibrated confidence intervals.
- Checked/optimized GNU Fortran regression tests: 60 independent random real/complex cases and 788 CSD matrices per build (maximum complex entry error 3.27e-13), DPSS Fourier-quadrature/invariance tests, checked/optimized native end-to-end bootstrap smoke tests, and 106-module public signal closure with 25 public tests. No new Fortran runtime dependencies or C ABI changes.

## 2026-10-31: Exact fixed-weight multitaper VAR(p) population targets

- Add exact deterministic fixed-weight multitaper expected spectral matrices for stable VAR(p) models, including per-taper energy normalization, unequal/zero weights, complex channel phases, and analytically integrated band averages or powers.
- Expose static and band-integrated conditional canonical-coherence wrappers that combine the complete spectral matrices **before** applying nonlinear coherence and retain arbitrary group/conditioning, fixed shrinkage/ridge, and optional canonical vectors.
- Add checked and optimized native regressions, 40 independently generated SciPy/NumPy VAR(p) and DPSS/random-taper comparisons per build, and a public umbrella regression to the full 105-module signal closure; preserve all previous public signal tests.
- Document fixed-weight and population-target limitations. These routines are an original SciFort extension, not a SciPy API, and do not estimate uncertainty or assume bootstrap confidence-interval calibration.

## 2026-10-30: Conditional multimode VAR(p) bootstrap coverage study

- Add reproducible exactly stationary Gaussian VAR(2)/VAR(3) simulation with a fifth conditioning channel, two sorted canonical roots, complex point-frequency spectra, and analytically integrated band-averaged spectra from nonoverlapping tapered windows.
- Compare all five existing moving-block-bootstrap interval methods using the *matching* finite-taper VAR(p) population canonical target, including interval validity, Wilson Monte Carlo uncertainty and point-estimator bias/RMSE. Export 180 model/frequency/statistic/method results in a checked-build CSV.
- Validate 4 simulation models and exact band matrices independently by 128-point Gaussian quadrature, and compare all Fortran full-sample canonical fits against an independent SciPy/NumPy SVD. Explore near-degenerate modes without claiming calibrated inference. No production API change or extra runtime dependency.

## 2026-10-30: Exact finite-window VAR(p) spectral and canonical population targets

- Add `varp_expected_taper_spectrum` and `varp_expected_band_taper_spectrum`, including stable real VAR(p) companion covariance, explicit complex cross-spectral phase, normalized or raw finite-taper spectra and **analytically integrated** frequency bands (with optional band averages).
- Add `varp_expected_canonical_coherence` and `varp_expected_band_canonical_coherence` for conditional canonical-mode and group-coherence population functionals, without changing the existing VAR(1) and resampling APIs.
- Add checked/optimized native regressions, independent SciPy/NumPy time-pair and Cholesky/SVD reference comparisons, 20,000 simulated VAR(2) windows, a complete public-module test, and mathematical documentation. This remains a SciFort extension, not a `scipy.signal` API or a calibrated-inference claim.

## 2026-10-30: Exact finite-window VAR(1) spectral population targets

- Add SciFort extension `var1_expected_taper_spectrum` to compute raw or energy-normalized expected tapered outer-product matrices from stable real multivariate VAR(1) coefficients and innovation covariance, solving the stationary discrete Lyapunov equation by Smith doubling without new dependencies. Expose through `scifort_signal`; report invalid and nonstationary status.
- Add native edge/analytic tests, independent SciPy Lyapunov and sample-pair NumPy checks (80 randomized scenarios, 240 matrices), and a stationary-start Monte Carlo expectation test (60,000 windows).
- Replace mismatched zero-frequency comparisons in a new four-scenario exploratory canonical-bootstrap coverage diagnostic with the exact **finite-window** population functional. Evaluate percentile, basic, BC, BCa and studentized intervals for 64 datasets per scenario and report valid fractions and Wilson Monte Carlo uncertainty. Existing coverage scripts remain for historical comparison.
- Add documentation, provenance, and public-interface regressions; do not change the existing bootstrap estimator, seeds, or interval methods.

## 2026-10-30: Block-aware BCa canonical bootstrap intervals

- Add opt-in `interval_method='bca'` for original-window spectral-matrix and
  real/complex adaptive multitaper moving-block canonical bootstraps; preserve
  all previous methods and seeded bootstrap moments exactly.
- Calculate delete-one-contiguous-block jackknife acceleration from original
  full complex cross-spectra, with a weighted pseudovalue approximation for
  incomplete final blocks. Report mode and aggregate accelerations only on
  demand; reject fewer than four blocks, degenerate influences, boundary or
  nearly tied modes, nonidentifiable bias corrections and acceleration poles.
- Verify 56 independent NumPy/SciPy scenarios (1,683 numeric quantities;
  maximum error 7.8e-13), checked and optimized GNU Fortran 14.2 native
  regressions, clean 101-module umbrella build and 21 public-interface tests.
  Add reproducible two-scenario VAR(1) coverage diagnostics and
  `docs/signal_canonical_bca.md`. BCa remains exploratory under serial
  dependence, especially with spectrally smoothed estimands.

## 2026-10-30: Studentized canonical moving-block-bootstrap intervals

- Add `interval_method='studentized'` to the precomputed and real/complex
  multitaper canonical moving-block bootstrap functions. Preserve the
  percentile, basic and bias-corrected defaults, bootstrap draws, bias
  estimates and standard deviations.
- Use a delete-one-*drawn-block* jackknife for each bootstrap replicate and a
  delete-one-contiguous-block jackknife on the original spectral windows.
  Estimate uncertainty using block pseudovalues with weighted handling of
  incomplete final blocks, and invert the empirical bootstrap-t quantiles.
- Expose separate original `studentizing_stderr` and
  `coherence_studentizing_stderr` outputs only when requested. Reject
  unidentifiable/constant pivot denominators, fewer than four blocks,
  singular fits, or nearly tied individual modes with NaN intervals and
  explicit validity flags; retain available bootstrap moment estimates.
- Add native checked/optimized regressions, independent NumPy/SciPy numerical
  reference, complete 101-module umbrella/public tests, and
  `docs/signal_canonical_studentized.md`. A two-scenario 32-dataset VAR(1)
  study demonstrates possible undercoverage despite wider studentized bounds.
  No new Fortran dependencies. Studentized pointwise bounds are exploratory,
  not calibrated confidence tests.

## 2026-10-30: Bias-aware canonical moving-block bootstrap intervals

- Extend `canonical_moving_block_bootstrap` and both real/complex
  `multitaper_canonical_moving_block_bootstrap` procedures with trailing
  optional `interval_method` (`percentile`, `basic`, `bc`), retaining original
  percentile defaults, full-data estimates, RNG streams, and other outputs.
- Implement basic/reverse percentile bounds (intersected with `[0,1]`) and
  Efron bias-corrected percentile limits with numerically robust midranks;
  return NaN/invalid flags for unidentified bias adjustment. These are not
  BCa or studentized intervals.
- Add independent SciPy/NumPy reference comparisons (120 configurations,
  777 valid scalars in checked mode; 36 configurations, 222 scalars in
  optimized mode; maximum error below 6e-15), native and umbrella-API
  tests (101 modules, 19 public tests), and an exploratory two-scenario
  VAR(1) coverage study with 32 datasets per setting.
  No new Fortran runtime dependencies. See
  `docs/signal_canonical_interval_methods.md`.

## 2026-10-30: Statistic-specific influence block selection for canonical coherence

- Add `select_canonical_influence_block_length` using delete-one-complete-window
  changes in aggregate group coherence or a chosen isolated squared canonical
  mode, optionally at a particular frequency/band, with fixed ridge/shrinkage.
- Add opt-in `selection_method='influence'`, `selection_mode` and
  `selection_frequency` to precomputed and real/complex multitaper moving-block
  bootstraps. Preserve the original matrix-Frobenius selector as the default,
  and preserve manually specified block lengths and deterministic RNG behavior.
- Expose influence scores, target statistic and selection diagnostics. Reject
  unknown method, invalid mode/bin, singular refits or repeated/boundary roots.
- Pass checked and optimized native/regression tests; 40 independent NumPy
  complex/conditioned/band scenarios (1,535 delete-one scores, max observed
  error below 1e-12); the 101-module umbrella closure and 18 public tests.
- Add an exploratory three-channel VAR(1) coverage comparison showing that
  influence selection can improve targeting but does not fix percentile
  undercoverage. See `docs/signal_canonical_influence_block_length.md`.

## 2026-10-30: Automatic block-length pilot for canonical coherence bootstrap

- Add `select_canonical_block_length` with lagged complex Hermitian Frobenius autocovariances, a capped Bartlett-weighted plug-in heuristic and pilot diagnostics.
- Make `block_length=` optional in the existing precomputed and real/complex canonical moving-bootstrap APIs. Omitted lengths select automatically, while explicit lengths retain their earlier behavior (including rejection of zero).
- Add resampling metadata, independent NumPy checks, and an exploratory bivariate Gaussian AR(1) percentile-interval coverage study. Treat automatic selection as a heuristic, not a calibrated uncertainty estimator.
- Add native, independent-reference and public-module tests plus `docs/signal_canonical_auto_block_length.md`.

## 2026-10-30: Moving-block bootstrap for canonical spectral coherence

- Add `canonical_moving_block_bootstrap` and real/complex generic
  `multitaper_canonical_moving_block_bootstrap` using overlapping possible
  contiguous-window blocks drawn with replacement; support optional circular
  wrapping, non-divisible block lengths and deterministic local RNG seeds.
- Return full spectral canonical mode strengths, aggregate group coherence,
  bootstrap standard errors, mean-minus-full-estimate biases, pointwise
  2.5%/97.5% percentile bounds, per-mode validity masks and sampling metadata.
  Preserve fixed ridge/shrinkage and existing conditional-group, multitaper
  and frequency-band APIs without changing the earlier jackknife functions.
- Pass checked and optimized GNU Fortran builds, native and previous numerical
  regressions, 56 independent SciPy/NumPy reference cases (205 valid modes,
  maximum absolute error 3.0e-15) and the full 100-module `scifort_signal`
  dependency closure (16 public-interface tests). No new Fortran dependencies.
- Pointwise bootstrap percentiles are **not calibrated significance tests**;
  document stationarity, block-length sensitivity, boundary and mode-tie
  caveats in `docs/signal_canonical_moving_bootstrap.md`.

## 2026-10-30: Contiguous equal-block canonical coherence jackknife

- Add `canonical_block_jackknife` and the real/complex
  `multitaper_canonical_block_jackknife` generic. The old window jackknife
  APIs gain an optional trailing `block_length=` argument (default 1).
- Omit one disjoint contiguous group of `block_length` complete windows per
  refit. Require at least four equal-size groups and exact divisibility;
  use the number of groups in the delete-block SE and bias formulas.
- Preserve existing full-window point estimates, integrated-band behavior,
  fixed regularization, validity masks and Fisher-scale heuristic intervals.
  Record `n_windows`, `n_blocks`, and `block_length` separately.
- Pass checked and optimized native regressions, 80 independent NumPy/SciPy
  cases (284 canonical modes per build, max absolute error 4.44e-15), and
  the full 99-module public interface with 15 passing tests. No new
  Fortran runtime dependencies. See `docs/signal_canonical_block_jackknife.md`.
- This SciFort extension allows a basic short-memory-aware resampling
  strategy; equal-block jackknife intervals are **not calibrated** and do
  not establish independence or simultaneous frequency coverage.

## 2026-10-30: Nonoverlapping-window canonical coherence jackknife

- Add `canonical_window_jackknife` on existing multichannel
  time-frequency spectral matrices, and generic real/complex
  `multitaper_canonical_window_jackknife` entry points that always
  estimate nonoverlapping windows (`noverlap=0`).
- Return full-mean canonical squared strengths, delete-one-window SEs,
  Fisher-scale SEs, jackknife bias, approximate pointwise 95% Fisher
  intervals, aggregate group coherence SEs and per-mode validity masks.
- Integrate complex Hermitian CSD matrices over frequency bands
  *within each window*, not after estimating canonical coherence.
  Exclude boundary/tied roots from ordinary rankwise inference and
  preserve previous canonical and covariance-metric tracking APIs.
- Pass checked and optimized native suites, independent NumPy/SciPy
  comparisons (100 scenarios, 357 modes, max 5.88e-15) and a clean
  umbrella public build (99 modules, 14 tests). Fortran introduces
  no new dependencies; this is a SciFort extension, not SciPy API.
- Document that these heuristic intervals are not calibrated 95%
  inference, especially under serial dependence, small window counts,
  automatic regularizer selection or mode degeneracy.

## 2026-10-30: Covariance-metric-aware canonical subspace tracking

- Add opt-in `return_covariances=.true.` to all canonical-coherence matrix,
  time-frequency, band, and real/complex multitaper APIs; return the exact
  conditioned/regularized group A/B autospectral matrices used for canonical
  coefficient normalization. Preserve the unallocated default.
- Add opt-in `covariance_metric=.true.` to `track_canonical_subspaces` for
  static and time-frequency results. Compare neighboring blocks in a shared
  midpoint conditional spectral covariance metric, using Hermitian Cholesky
  transformations before the existing principal-angle and Procrustes steps.
- Validate finite Hermitian positive-definite metrics, array shapes and source
  validity masks. Preserve legacy Euclidean tracking by default; report the
  chosen metric in `covariance_metric_used`. Fixed channel rescalings and
  nonsingular group-coordinate changes preserve the metric-aware overlaps.
- Pass new native tests, 80 independent randomized NumPy/SciPy reference cases
  (1,791 modes, 503 linked, maximum error 1.68e-14), plus 120 independent
  complex non-diagonal coordinate-invariance cases, in checked and optimized
  builds, and the 98-module public dependency closure with 13 tests. The
  original native subspace regression passes in both builds.
  See `docs/signal_canonical_covmetric_tracking.md`. No new Fortran dependency.

## 2026-10-30: Principal-angle tracking of near-degenerate canonical subspaces

- Add opt-in `track_canonical_subspaces` for static and time-frequency canonical
  results with coefficient vectors, without changing existing canonical or
  individual-mode tracking APIs. Group near-equal squared coherence roots.
- Match equal-dimensional subspaces by the smallest principal cosine, use a
  globally optimal Hungarian block assignment, and align paired coefficient
  bases by a common complex unitary Procrustes rotation. Preserve original
  canonical strength ordering separately from rotated subspace bases.
- Report persistent block IDs, dimensions, per-block principal cosines, match
  masks and new-label flags; start new labels at splits, merges, zero roots,
  invalid spectra, weak overlap and nonascending/FFT-wrapped frequencies.
- Pass checked/optimized native and prior canonical regressions, 100 independent
  NumPy/SciPy cases and 2,215 comparison modes per build (maximum absolute
  error 8.88e-16), plus the 98-module public dependency closure and 12 tests.
  See `docs/signal_canonical_subspace_tracking.md`. No new runtime dependency.

## 2026-10-30: Frequency-continuous conditional canonical modes

- Add opt-in `track_canonical_modes` for static and time-frequency canonical
  results computed with coefficient vectors. Preserve all existing descending
  canonical-mode APIs and expose separately tracked mode values and phases.
- Find globally optimal cross-frequency matches with the Hungarian assignment
  algorithm on joint normalized complex-vector overlaps. Apply the same complex
  gauge to both groups, preserving their positive canonical cross-spectrum.
- Record original descending indices, match scores, per-mode link masks,
  per-bin segment starts and source validity; reset at missing/unstable roots,
  weak matches and nonascending or FFT-wrapped frequencies. TF windows are
  independent trajectories.
- Pass native checked/optimized regressions, 80 independent NumPy/SciPy
  assignment/phase scenarios (2,307 modes per build, error <= 6.3e-16),
  and the full 97-module public dependency closure with 11 passing tests.
  See `docs/signal_canonical_tracking.md`. No new Fortran dependencies.

## 2026-10-29: Optional canonical spectral coefficient vectors

- Extend all canonical coherence matrix, band, and real/complex multitaper
  interfaces with an opt-in `return_vectors=.true.` argument; leave existing
  results and allocation-free defaults unchanged.
- Report complex, conditional spectral-variance-normalized coefficients for
  each group's ordered canonical mode, with paired positive cross covariance,
  a stable phase convention, and a numerical direction-stability mask.
- Complete orthonormal null-space directions deterministically; mark zero,
  repeated, and near-tied roots as having unreliable individual orientation.
- Verify 150 independent NumPy SVD scenarios covering 798 cells in each of
  checked and optimized builds, native edge cases and legacy regressions, and
  the full 96-module public API with 10 passing public tests.
- No new Fortran dependencies or third-party implementations. This is a
  SciFort extension; see `docs/signal_canonical_coherence.md`.

## 2026-10-29: Individual conditional canonical spectral coherence modes

- Add `canonical_coherence_matrix`, `canonical_coherence_tf_matrix` and
  frequency-band variants, plus generic real/complex multitaper interfaces
  for static and moving-window spectra.
- Return ordered squared canonical correlations between two channel groups,
  optionally controlling for a third, with determinant group coherence,
  shrinkage coefficients and per-cell validity masks.
- Solve complex Hermitian Schur complements using Cholesky and whitened
  cross-spectrum singular-value squares via a small Jacobi eigensolver.
  Preserve automatic and adaptive-CV shrinkage and explicit band regularization.
- Add checked and optimized regression suites and 240 independent NumPy SVD
  comparisons covering 1,848 time-frequency cells per build, plus full public
  module and compatibility tests. This is a SciFort extension; see
  `docs/signal_canonical_coherence.md`.

## 2026-10-29: Conditional determinant group spectral coherence

- Add `group_coherence_matrix`, `group_coherence_tf_matrix` and direct real/
  complex `multitaper_group_coherence` and `multitaper_tf_group_coherence`
  for two channel groups optionally conditioned on a third group.
- Compute `1-Wilks' Lambda` from four stable Hermitian Cholesky log
  determinants; return the nonnegative negative log ratio, valid masks and
  existing fixed/automatic/adaptive-CV spectral regularization metadata.
  Single-channel reductions agree with partial or multiple coherence.
- Add `band_group_coherence_matrix`, `band_group_coherence_tf_matrix`, and
  direct real/complex static/moving-window multitaper band versions.
- Add native checked/optimized tests, 196 independent NumPy canonical-SVD
  comparison scenarios (1,418 cells, max error below 3e-14), and a full
  95-module public build (eight public tests). See
  `docs/signal_group_coherence.md`. This is an original SciFort extension,
  not a `scipy.signal` function.

## 2026-10-29: Band-integrated spectral coherence and partial coherency

- Add `integrate_spectral_bands` and `integrate_tf_spectral_bands` to form
  frequency-band spectral matrices by exact piecewise-linear quadrature.
  Accept nonuniform, FFT-wrapped or arbitrarily ordered complex spectra,
  overlapping bands and off-grid endpoints; reject extrapolation.
- Add band-level partial/multiple coherence for static and time-frequency
  matrices, plus real/complex direct multitaper entry points. Preserve
  conditional power, partial phase, explicit regularization and validity
  masks; distinguish band-level coherence from mean squared coherence.
- Add checked/optimized native tests, independent randomized NumPy
  integration/inversion comparisons (180 scenarios / 1,416 cells) and
  full public API verification. No new runtime Fortran dependencies.
  See `docs/signal_band_coherence.md`.

## 2026-10-29: Multichannel spectral multiple coherence

- Add `multiple_coherence_matrix`, `multiple_coherence_tf_matrix`,
  `multitaper_multiple_coherence`, and `multitaper_tf_multiple_coherence`
  to `scifort_signal`. These compute each channel's squared multiple coherence,
  conditional residual spectral power, and explained spectral power.
- Reuse the Hermitian Cholesky precision solver with trace-based ridge,
  explicit/automatic/CV spectral shrinkage, and frequency/time validity
  masks. Support real/complex input and static/moving-window spectra.
- Add native checked/optimized, six public API, and 700 independent NumPy
  Schur-complement comparisons. See `docs/signal_multiple_coherence.md`.
  This is a SciFort extension, not a `scipy.signal` function.

## 2026-10-29: Predictive cross-validated shrinkage with adaptive multitaper weights

- Add opt-in `cv_shrinkage=.true.` for real/complex static and moving-window
  multichannel spectral matrices and partial coherence. For each frequency
  remove one DPSS tapered spectrum, re-estimate any adaptive weights, and
  choose bounded diagonal-target shrinkage by held-out off-diagonal loss.
- Support equal, eigen and data-dependent adaptive weighting; preserve
  Gaussian-Wishart `auto_shrinkage` as a separate fixed-weight-only method.
  Store per-frequency or per-time-frequency CV estimates on spectral matrix
  results and validate conflicts, shapes and coefficient finiteness.
- Add native and public-module regressions and separate independent
  NumPy/SciPy comparisons for 12 configurations / 1,677 bins. Checked
  and optimized GNU Fortran builds pass. The method estimates predictive
  shrinkage, **not** calibrated inference or oracle covariance risk.
  See `docs/signal_cv_shrinkage.md`.

## 2026-10-29: Eigenvalue-weighted DPSS automatic spectral shrinkage

- Extend `auto_shrinkage=.true.` to deterministic eigenvalue-weighted
  tapers for real/complex static and time-frequency partial coherence.
  Use weighted-Wishart equivalent sample count `(sum(w))**2/sum(w**2)`.
- Expose `%effective_tapers` metadata from static and moving-window
  spectral matrices and add public `multitaper_effective_tapers` with
  strict validation and overflow-safe scaling. Reuse the estimated DPSS
  eigenvalues without an additional DPSS solve.
- Preserve equal-weight, explicit shrinkage and ridge behavior; continue
  rejecting observation-dependent adaptive weights in automatic mode.
- Add native, numerical regression, and public-interface tests; checked
  and targeted optimized GNU Fortran 14.2 builds pass the existing
  196-case independent NumPy comparison and multitaper regressions.
  See `docs/signal_auto_shrinkage.md`.

## 2026-10-29: Automatic diagonal-target shrinkage for partial coherency

- Add opt-in `auto_shrinkage=.true.` to static and time-frequency
  `partial_coherence_matrix`, `partial_coherence_tf_matrix`,
  `multitaper_partial_coherence` and `multitaper_tf_partial_coherence`.
  Compute a separate Gaussian-Wishart moment plug-in shrinkage estimate
  at each frequency/time cell, retaining the existing Cholesky solver.
- Return actual coefficients in `%shrinkage_used(:)` or
  `%shrinkage_used(:,:)`; keep explicit `shrinkage=` and old defaults.
  Matrix APIs require an `effective_tapers=` count; signal-level APIs
  infer it from DPSS parameters and reject non-equal taper weighting
  when auto-selection is requested. Explicit/manual options conflict.
- Compare 196 independent NumPy cases in both checked and targeted
  `-O3` builds: 949 invertible and 8 zero/singular cells per build,
  85,545 precision/coherency numerical values, max scaled error 4.27e-15.
  Three earlier 288-, 252-, and 120-case NumPy/SciPy suites, six focused
  native regressions and four public API tests pass.
- See `docs/signal_auto_shrinkage.md` and
  `tools/run_signal_partial_shrinkage_auto.py`. Statistical independence
  is an approximation for tapered spectral data; no inferential claims.

## 2026-10-29: Complex partial coherency and diagonal-target spectral shrinkage

- Add optional `return_coherency=.true.` to static and time-frequency partial
  coherence, producing complex conditional coherency while retaining phase.
  The default returns only squared coherence and allocates no complex output.
- Add `shrinkage=` with a coefficient in `[0,1]`: blend each cross-spectral
  matrix with its diagonal before existing trace-ridge loading. This offers
  individually scale-invariant regularization when ridge is zero.
- Checked against 288 independent NumPy/SciPy cases, 105,763 complex
  coherency entries and three expected singular cells; maximum complex
  coherency error `1.073e-14`. Both checked and targeted `-O3` builds pass.
  Existing 252-case static and 120-case time-frequency comparisons pass in
  a checked build; previous native multitaper regressions remain green.
- See `docs/signal_partial_coherence.md`,
  `docs/signal_tf_partial_coherence.md` and
  `tools/run_signal_partial_coherency.py`.

## 2026-10-29: Moving-window multitaper partial coherence

- Add `multitaper_tf_partial_coherence` for real and complex multichannel
  signals, with overlapping DPSS windows, selectable weighting and detrending,
  user-selected sample axis, and optional trace-scaled diagonal loading.
- Add `partial_coherence_tf_matrix` for precomputed cross-spectral matrix
  spectrograms and `partial_coherence_tf_result`, with precision and squared
  partial coherence `(frequency,time,channel,channel)` and independent
  singularity masks `(frequency,time)`.
- Reuse the existing complex Hermitian Cholesky partial-coherence kernel.
  The old static partial-coherence API remains unchanged.
- Independent NumPy/SciPy 1.17 comparison passes 120 cases, 39,786 valid
  squared-coherence entries, and 16 singular time-frequency cells. Maximum
  coherence error `9.60e-15`, maximum CSD error `6.44e-15`.
  Checked and targeted optimized GNU Fortran builds and prior multitaper
  native tests pass; earlier 252-case partial-coherence reference passes.
- See `docs/signal_tf_partial_coherence.md` and
  `tools/run_signal_tf_partial_coherence.py`.

## 2026-10-29: Multichannel partial spectral coherence

- Add `partial_coherence_matrix` and `multitaper_partial_coherence` for
  squared pairwise partial coherence, conditioned on all other channels.
  Support real or complex multichannel time series, all existing multitaper
  estimator options, and direct use of precomputed cross-spectral matrices.
- Use Hermitian Cholesky factorization with optional scale-invariant diagonal
  loading, per-frequency singularity flags, and complex precision matrices.
  Preserve existing multitaper and frequency-analysis APIs unchanged.
- Check 252 independent NumPy/SciPy cases with 35,270 conditional-coherence
  entries and 33 singular frequency bins in checked and targeted `-O3` builds;
  retain the earlier native DPSS, multichannel, and adaptive regressions.
- See `docs/signal_partial_coherence.md` and
  `tools/run_signal_partial_coherence.py`.

## 2026-10-29: Weighted irregular-grid Savitzky-Golay smoothing and derivatives

- Add `savgol_irregular` for real and complex signals sampled at strictly
  increasing nonuniform times. Support inverse-variance observation weights,
  nearest-neighbor local windows, arbitrary query locations, derivatives,
  and explicit singular/invalid status handling.
- Reuse the existing Savitzky-Golay module and public `scifort_signal` export;
  preserve all previously supported uniform-grid interfaces. Use scaled,
  twice-reorthogonalized QR rather than fitting through normal equations.
- Checked and optimized builds each pass 220 randomized comparisons with
  independent NumPy weighted least squares (3,509 outputs, max absolute
  difference <7e-11) and 60 SciPy 1.17 uniform-grid comparisons (max
  difference <2e-11), plus the native Fortran regression.
- See `docs/signal.md` and `tools/run_signal_savgol_irregular.py`.

## 2026-10-29: Prepared multiharmonic fitting for multiple channels

- Add `harmonic_fit_plan`, `lombscargle_harmonic_prepare`,
  `lombscargle_harmonic_fit_prepared`, and `lombscargle_harmonic_fit_batch`.
  Reuse weighted QR factors for many signals on the same irregular time grid,
  with independent channel coefficients, powers, and optional covariances.
- Preserve the existing harmonic-fit results and prediction APIs. Handle
  singular frequency bins, zero weights, strided arrays, invalid channels,
  and independent precentering without changing earlier interfaces.
- Independent NumPy/SciPy comparisons pass 128 randomized multichannel cases,
  94,671 finite numerical values and 54 deliberately singular combinations
  in each of checked and optimized GNU Fortran builds (maximum absolute
  error `2.84e-13`). Previous harmonic-fit/uncertainty native regressions and
  the public `scifort_signal` smoke test pass.
- See `docs/signal_harmonic_plan.md` and
  `tools/run_signal_harmonic_plan.py`.

## 2026-10-29: Covariance and uncertainty for weighted harmonic fits

- Extend `lombscargle_harmonic_fit` with optional trailing
  `estimate_covariance=.true.`; return parameter covariance, standard errors,
  per-frequency residual variance estimates, residual degrees of freedom,
  and an indicator for whether the intercept was fitted.
- Add `lombscargle_harmonic_predict_std` for uncertainty in the fitted mean
  curve (not future observation noise); handle fixed intercepts, zero-weight
  observations, partial rank failures, and invalid prediction times.
- Compute harmonic covariance using the existing twice-reorthogonalized QR
  factor, avoiding weighted normal-equation inversion. Preserve the default
  no-covariance fitting path and earlier positional argument list.
- Checked and optimized GNU Fortran builds each pass 256 independent NumPy
  weighted-regression cases (59,796 numerical values, maximum error <1.6e-15)
  and existing native regression tests; see
  `docs/signal_harmonic_uncertainty.md` and
  `tools/run_signal_harmonic_uncertainty.py`.

## 2026-10-28: Harmonic least-squares coefficients and waveform prediction

- Add `lombscargle_harmonic_fit` and `harmonic_fit_result` for stable recovery
  of weighted multi-harmonic cosine/sine coefficients, the fitted intercept,
  two periodogram normalizations, and reconstruction RMS at every frequency.
- Add `lombscargle_harmonic_predict` for evaluation on new irregular times,
  preserving the original data scale with optional precentering, supporting
  partial rank failure without discarding valid frequency bins.
- Implement the triangular back-substitution on a twice-reorthogonalized QR
  decomposition instead of forming the weighted normal equations; retain
  the existing `lombscargle` and `lombscargle_multiharmonic` APIs unchanged.
- Both checked and optimized GNU Fortran builds pass 256 independent NumPy
  coefficient/prediction comparisons (20,061 numeric values, maximum absolute
  error `7.96e-13`), the earlier 384-case periodogram regression, and the
  native and public-interface tests. See `docs/signal_harmonic_fit.md` and
  `tools/run_signal_harmonic_fit.py`.

## 2026-10-28: Weighted multiharmonic Lomb–Scargle spectral analysis

- Add `lombscargle_multiharmonic(x,y,freqs,nterms,...)` for joint weighted
  least-squares fitting of multiple sine/cosine harmonics to irregular data.
- Extend rank-one `lombscargle` with an optional trailing `nterms=` keyword;
  its omitted or one-term path retains the original SciPy-compatible kernel.
- Use reorthogonalized QR without normal-equation inversion, support weighted
  floating means and both power normalization modes, and detect singular
  trial frequencies without losing valid bins.
- Checked and optimized GNU Fortran builds pass 384 independent NumPy
  least-squares cases (2,493 frequency bins; maximum error about `5.2e-13`),
  including single-harmonic SciPy comparisons. The older 900-case irregular
  reference suite and public-module smoke test pass as well. See
  `docs/signal_multiharmonic.md` and `tools/run_signal_multiharmonic.py`.

## 2026-10-28: Fast skip-ahead and rewind for maximum-length sequences

- Extend `max_len_seq` with optional signed `integer(int64) skip=` to seek
  forward or backward before producing output, with `length=0` returning
  only the advanced state. Existing calls and bit-generation behavior remain.
- Use packed binary transition matrices and repeated squaring, requiring
  `O(nbits**2 log(abs(skip)))` operations instead of generating intermediate
  bits. Handle custom and nonprimitive taps, full 32-bit registers, and all
  signed 64-bit distances without overflow.
- Both checked and optimized GNU Fortran 14.2 builds pass the original and
  new native tests plus 698 independent SciPy 1.17 skip/rewind comparisons
  per build, with exact output bits and final states. See
  `docs/signal_sequence_skip.md` and `tools/run_signal_sequence_skip.py`.

## 2026-10-28: Parallel SOS discrete simulation for multichannel filters

- Add `dlsim_sos_bank` and overload existing `dlsim_sos`/`dlsim` for
  time-by-channel input matrices, accepting a shared real/complex SOS cascade
  or an independent real/complex SOS cascade for each input channel.
- Return channel-separated filter outputs, pre-sample DF-II state histories,
  and final SOS states to support continuation across streaming blocks.
- Support optional irregular timestamps, initial states, and noncontiguous
  inputs, with explicit validation of every channel's SOS coefficients.
- Checked and targeted `-O3` builds pass 256 SciPy 1.17 comparisons covering
  9,808 outputs and 171,680 pre-sample states (max absolute error
  `7.22e-16`), plus the earlier 240-case SOS simulation and 288-case SOS
  response suites. The public 91-source module compiles and its new smoke
  test passes. See `docs/signal_dlsim_sos_bank.md` and
  `tools/run_signal_dlsim_sos_bank.py`.

## 2026-10-28: SOS-native arbitrary-input simulation

- Add `dlsim_sos(sos,u,dt,t,zi)` and a rank-two SOS overload of `dlsim`
  for real and complex SOS arrays and arbitrary real/complex SISO inputs.
- Return sample time grids, pre-step DF-II filter state histories, and
  final SOS states for streaming continuation. Support optional irregular
  timestamps using the same interpolation as the existing `dlsim` API.
- Validate shapes, nonfinite inputs, normalized SOS, timestamp order, and
  numerical overflow with explicit statuses, preserving previous overloads.
- Compare 240 independent SciPy 1.17 cases, 9,120 outputs, and 310,080
  pre-step SOS state values in checked and targeted `-O3` builds (maximum
  absolute difference `4.34e-15`), plus 288 existing SOS response cases.
  See `tools/run_signal_dlsim_sos.py` and `docs/signal_dlsim_sos.md`.

## 2026-10-28: SOS-native discrete time responses

- Add `dstep_sos` and `dimpulse_sos` for real and complex digital SOS cascades,
  also exposed as rank-two SOS overloads of `dstep` and `dimpulse`.
  Accept optional DF-II states `zi` without expanding SOS to TF polynomials.
- Return existing discrete response type, with sample times and clear status
  on invalid dimensions, nonfinite coefficients, or unstable trajectories.
- Compare 288 random, Butterworth, high-order, and nonzero-state reference
  cases (15,153 values) against SciPy 1.17; both checked and -O3 builds pass.
  Retain all existing transfer-function and group-delay tests and public API.
  See `tools/run_signal_sos_responses.py` and `docs/signal_sos_responses.md`.

## 2026-10-28: Stable ZPK group delay

- Add `group_delay_zpk(z,p,k,...)` and the corresponding `group_delay`
  overload, accepting mixed real/complex root vectors and real/complex
  nonzero gains. Root-wise phase derivatives avoid polynomial expansion.
- Handle unequal numbers of zeros/poles (including pure sample delays),
  complex roots, sampling-rate frequency grids, singular unit-circle bins,
  and invalid parameter values. The public API retains prior overloads.
- Both checked and fully optimized GNU Fortran builds pass 221 independent
  SciPy 1.17 ZPK comparisons (6,774 frequency values; maximum absolute
  difference `2.35e-10`) as well as existing SOS and TF regressions.
  The public 91-source module compiles and its new smoke test passes.
  See `docs/signal_group_delay_zpk.md` and
  `tools/run_signal_group_delay_zpk.py`.

## 2026-10-28: Stable SOS group-delay analysis

- Add `group_delay_sos` for real and complex second-order-section filters,
  including the overload `group_delay(sos, ...)` and exports through the
  public `scifort_signal` module. It evaluates analytic group delay by
  section to avoid high-order polynomial coefficient cancellation.
- Preserve existing transfer-function group delay and SOS frequency
  response interfaces. Undefined bins return NaN with singular status,
  matching SciFort's existing error convention.
- Native tests and independently computed SciPy 1.17 comparisons cover
  124 cases and 2,612 frequency points (up to 128 SOS sections) in checked
  and `-O3` GNU Fortran builds; maximum error to the analytic reference
  is `2.004e-12`. See `tools/run_signal_group_delay_sos.py` and
  `docs/signal_group_delay_sos.md`.

## 2026-10-28: Complex-wavelet magnitude ridge peak detection

- Add `find_peaks_cwt_complex` for real and complex input signals using
  configurable Morlet2 or user-defined complex wavelets. Add
  `find_peaks_cwt_magnitude` to detect ridge peaks in a precomputed complex
  CWT. Selection uses strict local maxima of `abs(CWT)` and a percentile
  magnitude noise floor, with existing distance/gap/length/SNR controls.
- Reuse the existing ridge tracker without changing the real-valued Ricker
  detector's signed-maxima semantics. Extend the top-level `scifort_signal`
  exports and its public-module regression.
- Validate 280 independent SciPy 1.17 magnitude-ridge cases (129,961 CWT
  positions / 939 selected peaks) and all 300 preexisting Ricker ridge cases
  exactly in both checked and fully optimized builds. The older 240 complex
  CWT coefficient comparisons also pass. See
  `tools/run_signal_cwt_complex_peaks.py` and `docs/signal.md`.

## 2026-10-28: Complex Morlet and callback continuous-wavelet transforms

- Add `cwt_complex_wavelet` for complex wavelet callbacks on real or complex
  vectors, `cwt_morlet2` with adjustable angular center frequency, and the
  normalized `morlet2_wavelet` generator. Preserve the existing real Ricker
  transform and real ridge-based peak detector without API changes.
- Apply the complex conjugate and time reversal of wavelets for centered,
  zero-padded convolution, following SciPy's historical CWT convention.
  Reject malformed callback lengths, nonfinite coefficients and signals,
  invalid wavelet scales, and nonfinite Morlet frequencies.
- Verify 240 independent SciPy 1.17 convolution comparisons / 28,320 complex
  coefficients: maximum absolute difference `1.12e-14` checked and
  `1.19e-14` optimized. Checked 91-source public-module build and existing
  CWT tests passed. See `tools/run_signal_cwt_complex.py` and `docs/signal.md`.

## 2026-10-28: Mixed real/complex one-dimensional convolution

- Completed mixed-type rank-one `convolve`, `correlate`, `fftconvolve`,
  `oaconvolve`, and `choose_conv_method` overloads for both operand orders.
  Existing same-type APIs, direct/FFT/overlap-add kernels, and positional
  arguments remain intact.
- Added native Fortran tests for method/mode combinations, noncontiguous
  signals, complex conjugation in correlation, invalid inputs, and timing
  mode; independent SciPy comparison validates 960 cases / 182,304 values.
- Checked and targeted optimized GNU Fortran comparisons both pass (maximum
  absolute error `5.18e-13`). The 91-source public module, legacy convolution
  and FFT regressions, and public mixed-type smoke test also pass. See
  `docs/signal_mixed_convolution.md` and `tools/run_signal_mixed_convolution.py`.

## 2026-10-28: ShortTimeFFT callable detrending

- Added in-place real and complex callable `detrend=` overloads for
  `short_time_fft%stft_detrend`, covering vectors and rank-2/3 arrays with
  arbitrary existing transform axes. Callbacks act on each padded slice
  before windowing, with finite-output validation and error-status cleanup.
- Preserved string-based detrending, FFT scaling, phase shift and boundary
  rules, and exported the callback interface types from `scifort_signal`.
- Verified 272 SciPy 1.17 cases (50,800 complex coefficients; checked
  maximum absolute error `1.03e-13`), with native strided/channel and
  invalid-callback regressions. See `docs/signal.md` and
  `tools/run_signal_stft_callable_detrend.py`.

## 2026-10-28: Mixed-type cross-spectral analysis

- Added rank-one `csd` and `coherence` generic overloads for both real/complex
  input orders, reusing the complex two-sided Welch kernels, including median
  averaging, named and parameterized windows and per-segment callbacks.
- Extended rank-one `short_time_fft%cross_spectrogram` to mixed inputs,
  preserving its STFT(x)*conjg(STFT(y)) orientation and rejecting onesided
  FFT plans. Existing real/complex same-type and high-rank interfaces remain
  unchanged.
- Verified 384 independent SciPy 1.17 comparisons (31,008 values, maximum
  absolute error `5.15e-13`) and native/legacy spectral regressions in checked
  and targeted optimized GNU Fortran builds. SciPy's legacy mixed-dtype
  dispatcher is inconsistent, so legacy reference comparisons explicitly
  promote the real input to complex. See `docs/signal_mixed_cross.md`.

## 2026-10-28: Complex local Wiener filtering

- Extended `scipy.signal.wiener`-compatible rank-1/2/3 Fortran overloads to
  `complex(dp)` inputs, with implicit complex noise, explicit real noise, and
  explicit complex noise. Reused zero-padded local-neighborhood geometry and
  preserved all existing median, real and integer Wiener interfaces.
- Reproduced SciPy's complex pseudo-variance `mean(x*x)-mean(x)**2` and
  NumPy's real-first lexicographic comparison, rather than substituting a
  magnitude-based complex variance. Added finite-input/noise validation,
  correct even-sized-window centering and noncontiguous-array tests.
- Verified 27 independent SciPy 1.17 cases / 1,161 output samples (max
  absolute error `2.066e-15`) plus the prior 1,050 real/integer denoising
  comparisons in checked and optimized GNU Fortran builds. See
  `docs/signal_wiener_complex.md`; full `fpm test` remains unverified.

## 2026-10-28: Callback detrending for legacy spectral estimators

- Added optional `detrend_callback=` to vector real/complex `periodogram`,
  `welch`, `csd`, `coherence`, `stft`, `spectrogram` and their existing
  rank-independent `*_into` interfaces. Each callback receives a separate
  unwindowed segment and transforms it in place before the FFT; mixed real/
  complex estimators accept a complex callback.
- Preserve all positional APIs and default string-based detrending when
  omitted. Reject simultaneous `detrend=` and callback, reject nonfinite
  callback output, and keep caller-owned multidimensional outputs unchanged
  on invalid requests. Safely omit unassociated procedure pointers in the
  multidimensional dispatch to retain checked-build compatibility.
- Added new native rank-14/15, mixed-type, strided and failure-state tests;
  independently compared 48 SciPy 1.17 cases / 2,098 spectral values with
  maximum absolute error `3.122e-13` (checked) and `7.264e-13` (`-O3`).
  Existing 440-case spectrogram-mode and 192-case two-sided spectral checked
  regression suites also passed. The top-level public module and all 91
  source dependencies passed syntax checking. Full `fpm test` unavailable.
  See `docs/signal_detrend_callbacks.md`.

## 2026-10-28: Parameterized Fourier-domain resampling windows

- Add `window_parameter`, `window_parameter2`, and `window_coefficients`
  to real and complex `resample`, axis-aware rank-2/3 overloads and
  rank-independent `resample_into` up through rank 15.
- Reuse SciFort's periodic named-window generator and existing FFT-order
  weighting, real-frequency folding and complex-domain resampling kernels.
  Reject invalid or competing window forms, preserving caller outputs.
- Validate 252 SciPy 1.17 comparisons / 3,696 complex-valued samples
  (maximum error `1.523e-14`), checked and targeted optimized builds,
  native high-rank/strided/invalid-option tests, and prior 53+184 multirate
  comparisons. See `docs/signal_resample_parameter_windows.md`.

## 2026-10-28: Complex SOS coefficients for `decimate`

- Add complex `sos=` overloads to vector and axis-selected rank-2/3
  `decimate` and to the transactional rank-1--15 `decimate_into` API.
  Real signals may now produce complex output from complex SOS cascades.
- Reuse the existing independent complex `sosfilt`/`sosfiltfilt` routines,
  keeping the original real-SOS, FIR, transfer-function and designed-filter
  interfaces intact. Reject malformed or nonfinite input with status codes.
- Verify 288 independent SciPy 1.17 comparisons / 7,434 complex samples
  in checked and targeted optimized builds (maximum observed error zero
  at printed binary64 precision), plus native rank-15, strided and
  output-preservation tests. Prior 576-case transfer-function, 336-case
  custom-filter, 53-case multidimensional and 184-case complex-resampling
  comparisons also passed. See `docs/signal_decimate_complex_sos.md`.
  Full FPM test suite remains unverified.

## 2026-10-28: Arbitrary real and complex transfer functions for decimate

- Extend the existing `decimate` vector and rank-2/3 functions and rank-1–15
  `decimate_into` with trailing `b=` and `a=` coefficient arrays. Support real
  IIR systems, complex IIR systems, and real-to-complex filtering without
  disturbing previous positional or `taps=`/`sos=` calls.
- Reuse the independent `lfilter` and `filtfilt` kernels for IIR systems and
  the polyphase/`upfirdn` paths for scalar-denominator FIRs; normalize by
  `a(1)` and reject conflicting, incomplete, nonfinite or empty arguments.
- Match 576 SciPy 1.17 custom-DLTI cases / 16,320 samples to `8.951e-16`
  maximum absolute error, alongside rank-15, strided, invalid-input native
  cases and prior multirate regression suites in checked/targeted `-O3`
  builds. See `docs/signal_decimate_transfer_function.md`.

## 2026-10-28: Caller-defined FIR/SOS filters for decimate

- Add trailing optional `taps(:)` and `sos(:,:)` keywords to `decimate` and
  `decimate_into`, reusing the existing SciFort FIR polyphase and IIR SOS
  kernels. These are available for real and complex inputs, rank-2/3
  axis-selected functions and rank-1–15 assumed-rank output subroutines.
- Preserve old positional APIs and their existing default Chebyshev-I IIR
  filters. Explicit filter coefficients are mutually exclusive and cannot
  accompany the designed-filter `ftype`/`n` options; validate real SOS
  shape and normalized a0 values. Preserve preallocated outputs on errors.
- Pass 336 SciPy 1.17 independent custom-filter cases / 10,500 samples,
  max absolute error `1.288e-14`, plus new native shape/error tests and
  existing 53-case and 184-case multirate reference suites in checked and
  targeted optimized GNU Fortran builds. Full FPM tests remain unverified.
  See `docs/signal_decimate_custom_filters.md`.

## 2026-10-28: Parameterized symmetric windows for polyphase resampling

- Added trailing `window_parameter`, `window_parameter2`, and
  `window_coefficients` keywords to real/complex `resample_poly` and
  `resample_poly_into`, including the existing axis-selected rank-2/3 and
  assumed-rank-through-15 interfaces. Supports shaped Gaussian, Kaiser,
  Tukey, generalized Gaussian, Taylor and general-cosine FIR prototypes.
- Reused the common symmetric spectral window generator and existing
  polyphase/endpoint-padding implementation. Legacy `beta` and default
  Kaiser behavior remain unchanged; conflicting shape parameters and
  explicit FIR taps now fail predictably.
- New native tests and 256 independent SciPy 1.17 real/complex comparisons
  (6,976 samples, max error `2.987e-15`) passed in checked and targeted
  `-O3` builds. Prior 53 multidimensional, 184 complex-filter and 738
  high-rank resampling SciPy cases passed. Optimized compilation covered
  both changed modules, with unchanged large dependencies compiled checked.
  Full `fpm test` remains unverified. See
  `docs/signal_resample_poly_parameter_windows.md`.
## 2026-10-28: All five legacy spectrogram output modes

- Added trailing `mode=` (`'psd'`, `'complex'`, `'magnitude'`, `'angle'`,
  `'phase'`) to the real and complex rank-one `spectrogram` functions.
  An allocatable `%values` component of the existing `spectrogram_result`
  provides complex STFT values; the original `%power` component remains
  the default real PSD result and also holds magnitude/angle/phase outputs.
- Unified the rank-independent `spectrogram_into` output modes with the
  same one-dimensional kernel, preserving its existing public interface,
  strided-array support and unchanged-output validation behavior. Implemented
  frequency-wise phase unwrapping and SciPy-compatible signed-zero principal
  angles for real DC/Nyquist FFT bins; PSD scaling remains unchanged.
- Added native real/complex and multidimensional tests, a reproducible
  checked/optimized standalone build, and a seeded independent SciPy 1.17
  comparison of 440 cases / 58,254 spectral values (maximum absolute
  difference 8.31e-14). The public-module checked build compiled all
  91 dependent modules and passed `test_signal_timefreq_nd`. Its 396-case
  independent SciPy time-frequency comparison (727,408 coefficients) also
  passed (maximum absolute difference 2.13e-12). The full FPM suite remains
  unverified. See `docs/signal_spectrogram_modes.md`.

## 2026-10-27: Scalar/logical `firwin` and odd-edge multiband design

- Extend `firwin` with scalar cutoff overloads and logical `pass_zero`
  overloads, preserving existing string-based and array-based positional APIs.
- Match SciPy's multiband transition-edge parity for any number of cutoff
  frequencies >=2, including three- and five-edge designs. Reject even-length
  filters only when the selected passbands reach Nyquist.
- Verify 443 valid SciPy 1.17 reference designs (10,370 tap coefficients;
  maximum absolute error 4.44e-16) and 362 invalid designs in both checked
  and `-O3` GNU Fortran 14.2. The previous 320-case FIR width/window
  comparison suite also passes. New standalone native tests, coefficient
  driver, randomized SciPy comparison harness, and API documentation.

## 2026-10-27: FIR transition widths and shaped `firwin2` windows

- Added trailing optional `width` to `firwin`. A finite, positive transition
  width in sample-frequency units selects SciPy's Kaiser attenuation/beta
  prescription and overrides the named window. All existing arguments retain
  their positions, and invalid widths return an empty tap array and status.
- Added trailing `window_parameter`, `window_parameter2`,
  `window_coefficients`, and `window_values` to `firwin2`. The existing shared
  named-window validator now supports symmetric (`fftbins=.false.`) tapers as
  well as periodic ones; explicit finite `window_values` is a SciFort-only
  extension. Incompatible window combinations are rejected explicitly.
- Verified 320 SciPy 1.17 FIR configurations and 6,736 coefficients, maximum
  absolute difference `4.44e-16`, in both checked and optimized GNU Fortran
  14.2 builds. The previous 1,620-case design, 1,070-case advanced-FIR,
  and 192-case two-sided spectral SciPy suites also passed. Restored the
  standalone design reference script's expanded spectral dependencies.

## 2026-10-27: Full two-sided real-input legacy spectral estimators

- Added a trailing optional `return_onesided` logical argument to real
  `periodogram`, `welch`, `csd`, `stft` and `spectrogram` and their
  rank-independent `*_into` interfaces. Existing calls continue to default
  to one-sided spectra; `return_onesided=.false.` retains all `nfft` Fourier
  bins in FFT frequency order with negative frequencies and no factor-of-two
  power correction. Even-length Nyquist bins use negative frequency, as in
  SciPy/NumPy. `coherence` has no such option in the SciPy legacy interface.
- The same keyword is supported for complex spectrograms; as with other
  complex spectral estimates, complex inputs always produce two-sided output.
  No positional arguments were reordered, and the source generator for the
  rank-independent spectral wrappers now emits the additional keyword.
- Added a native rank-14/15, strided-array, failure-transaction, frequency
  parity and PSD/CSD regression; the independent reference suite compared
  192 SciPy 1.17 cases and 7,248 spectral coefficients, with maximum absolute
  error below `3.77e-15`. Checked `-fcheck=all` and optimized `-O3` GNU Fortran
  builds passed the new tests and previous spectral-window, STFT-boundary,
  inverse-window, and ISTFT regressions. No new Fortran dependency was added.

## 2026-10-27: Parameterized inverse STFT windows and overlap checks

- Extended the real/complex legacy `istft` and rank-independent `istft_into`
  interfaces with trailing `window_parameter`, `window_parameter2`, and
  `window_coefficients` optional arguments. Named Gaussian, Kaiser, Tukey,
  DPSS, general Gaussian, Taylor, general cosine, and other supported spectral
  windows now share the same family/shape validation as the forward STFT.
  Explicit `window_values` remain supported, but cannot be combined with
  incompatible named/parameterized window requests.
- Extended named-window `check_cola`/`check_nola` to accept the same optional
  shape parameters. Retained all previous positional arguments and the
  existing explicit-vector checks.
- Added a native real/complex, rank-14, strided, and failure-transaction
  regression and an independent SciPy coefficient-driven ISTFT comparison.
  Both checked and `-O3` GNU Fortran 14.2 builds passed 192 SciPy 1.17 cases
  (5,613 reconstructed samples, max absolute error `6.57e-15`), together with
  the existing ISTFT, STFT boundary, and forward parameterized-window tests.
  No additional runtime dependency or imported source code.

## 2026-10-27: Parameterized window support for legacy `scipy.signal` estimators

- Added trailing optional `window_parameter`, `window_parameter2`,
  `window_coefficients`, and `window_values` arguments to the legacy real and
  complex `periodogram`, `welch`, `csd`, `coherence`, `stft`, and `spectrogram`
  kernels and to their assumed-rank `*_into` wrappers. Existing calls and
  positional parameter ordering remain valid.
- Caller-supplied window arrays determine `nperseg` when it is omitted;
  explicit mismatched segment sizes, conflicting window specifications,
  invalid parameter combinations, and nonfinite weights are rejected.
  In the spectrogram family, an omitted window now matches SciPy's Tukey
  alpha 0.25 rather than the alpha 0.5 default of explicitly named Tukey.
- Extended the spectral wrapper generator, added native real/complex,
  rank-14/rank-15, strided-array, and transactional-error regressions,
  and added a standalone checked/O3 build-and-test runner. All 96 SciPy 1.17
  reference configurations and 3,264 coefficients agree to absolute errors
  below `3.56e-14` in both builds; the previous STFT boundary test passes.
  No new Fortran dependencies or third-party implementation code.

## 2026-10-27: SciPy-compatible legacy STFT boundary extensions

- Extended real and complex `stft` and the `stft_into` assumed-rank adapter
  to accept `boundary='even'`, `'odd'` and `'constant'` alongside existing
  `'zeros'` and `'none'`. Reflection and endpoint padding occur before any
  optional trailing zero-padding; all extended modes use aligned zero-origin
  frame times, including odd-length windows. `spectrogram` remains unchanged.
- Kept a single real/complex extension helper in the existing spectral module
  so downstream standalone compiler source lists and FPM dependencies remain
  unchanged. The wrapper preserves transactional output semantics for invalid
  padding modes and noncontiguous array sections.
- Added native real/complex, rank-14, strided-array, and invalid-option tests,
  plus a deterministic 1,080-case independent SciPy 1.17 comparison covering
  54,993 STFT coefficients. The largest checked-build absolute difference was
  `8.21e-15`. No external Fortran dependencies or imported source code.

## 2026-10-27: Arbitrary-rank Hilbert transforms

- Added caller-output `hilbert_into` and `hilbert2_into` through rank 15,
  including selected one-based/negative transform axes, FFT-length padding or
  truncation, noncontiguous inputs/outputs, and transactional error handling.
  The older vector and matrix functions are unchanged.
- The two-dimensional mask follows SciPy 1.17's single-orthant convention,
  including omission of even-length Nyquist bins. Repeated channel FFTs use
  local reusable plans; no external Fortran dependency was added.
- Added two native tests, an independent SciPy comparison driver and harness,
  API documentation, and standalone test-runner dependency integration.
  All 500 independent comparisons (191,284 coefficients) pass in checked
  and optimized GNU Fortran builds, with max absolute error `1.47e-14`.

## 2026-10-27: Rank-independent analytic envelope analysis

- Added `envelope_into` real/complex caller-output overloads through rank 15,
  preserving all existing rank-1/2/3 `envelope` calls. The selected sample axis
  accepts one-based or negative indexing, and `n_out` adjusts that dimension.
- Reused the established FFT-domain envelope kernel for signed frequency-bin
  bands, optional squared magnitude, and lowpass/all/no residual modes. A
  supplied residual array has the input kind; `residual='none'` leaves it
  unchanged. Rank 1-15 strided input/output sections are supported.
- Staged both outputs until every channel succeeds. Invalid shapes, band
  limits, modes, nonfinite samples, or failed later channels leave caller
  outputs unchanged. Native rank-15 tests and 360 independent SciPy 1.17
  rank-four comparisons passed in checked and optimized GNU Fortran builds.
- No added external Fortran dependencies or copied SciPy implementation code.

## 2026-10-27: Multidimensional relative extrema through rank 15

- Added `argrelmax_mask_into`, `argrelmin_mask_into` and
  `argrelextrema_mask_into` to the public `scifort_signal` module for
  `real(dp)` input arrays of rank 1 through 15. The new interfaces evaluate
  neighbors on a selected one-based/negative axis, with SciPy-style `order`,
  `'clip'`/`'wrap'` boundaries, and a pure custom comparator. Logical output
  masks are a SciFort extension to SciPy's coordinate-tuple return values.
- Support zero-sized axes, singleton axes, ties, NaNs, arbitrary strided
  input/output sections, and complete preservation of caller output on invalid
  axes, output shapes, modes, or orders. Existing rank-one index-returning
  functions are unchanged.
- Added native public-module and rank-15 regression tests and 720 independent
  SciPy 1.17 comparisons per build (52,661 input positions). No new Fortran
  dependencies or imported third-party implementation code.

## 2026-10-27: Reusable weighted Lomb-Scargle plans

- Added the public `lombscargle_plan` type with `initialize`, `clear`, and
  `is_initialized` methods, and `lombscargle_plan_into` and
  `lombscargle_plan_amplitude_into` for real observations through rank 15.
  The plan caches weighted, phase-rotated sine/cosine columns and their Gram
  denominators for a fixed observation-time grid, frequency grid, shared
  weights, and floating-mean setting. Calls on new channel data avoid
  repeated trigonometric setup and retain one-based/negative sample axes.
- Plan initialization is transactional (invalid reinitialization retains a
  valid plan); evaluation stages its output so invalid shapes, nonfinite
  observations, and singular channels never partially overwrite results.
  Cached plans are read-only during evaluation and may be used concurrently.
  Per-channel weights and independently omitted samples remain available
  through the unplanned `lombscargle_into`/`lombscargle_amplitude_into` API.
- Added native rank-15/strided/error-recovery tests, a deterministic Fortran
  SciPy-reference driver, a 192-case Python comparison, and an `-O3`
  benchmark. No external Fortran libraries or runtime dependencies added.

## 2026-10-27: Per-channel Lomb-Scargle weights and cached trigonometry

- Extended the rank-1..15 `lombscargle_into` and
  `lombscargle_amplitude_into` interfaces with optional `weights_nd`, a
  per-sample/per-channel real weight array matching the input shape, and
  `nan_policy='raise'|'omit'`. Omission independently masks nonfinite
  channel observations and zero-weight samples, matching explicitly masked
  rank-one SciPy calls. These options are SciFort extensions, with no change
  to the original scalar Lomb-Scargle API.
- Shared-weight multichannel calls now cache the sine/cosine frequency basis
  and weighted Gram terms across channels rather than repeating trigonometry
  for every lane; the original scalar kernel remains the reference for
  per-channel weighting and omitted observations. All-or-nothing output
  semantics, arbitrary axes, rank 15, and strided arrays are retained.
- Added a native weighted/missing-data and error-atomicity regression,
  independent SciPy reference driver and comparison script, and an
  old-versus-new performance benchmark. No new runtime dependencies.

## 2026-10-27: Multichannel Lomb-Scargle spectral estimation through rank 15

- Added `lombscargle_into` (real power or normalized power) and
  `lombscargle_amplitude_into` (complex fitted amplitudes) to the public
  `scifort_signal` module. Real observations of rank 1..15 share irregular
  sample times, optional weights, and angular frequencies; callers select
  a one-based/negative sample axis and provide the appropriately shaped
  output array. This is a SciFort extension beyond SciPy's 1-D API.
- Reused the existing weighted scalar least-squares implementation for
  each lane, keeping its normalization, precentering, floating-mean and
  status semantics. Outputs are staged for all channels, so invalid input
  or a singular later channel leaves the full caller array unchanged.
  Noncontiguous inputs/outputs and negative axes are supported.
- Added a native rank-15/strided/error-atomicity regression, a Fortran
  reference driver, and an independent comparison script against SciPy
  1.17. No new Fortran runtime dependencies or third-party code.

## 2026-10-27: Arbitrary-axis CZT and zoom FFT through rank 15

- Added `czt_into`, `zoom_fft_into`, `czt_plan_into`, and
  `zoom_fft_plan_into` to the public `scifort_signal` module. The first two
  create a transform plan once per multichannel operation; the latter two
  reuse an initialized `czt_plan`/`zoom_fft_plan` across batches.
- Real and complex assumed-rank inputs (1-15) produce complex caller-allocated
  outputs, with a one-based or negative transform axis (default last),
  scalar or interval zoom-frequency arguments, and existing CZT spiral
  parameters. Input/output extents are validated before computation.
- The adapters stage all channels before assigning to the result. They
  preserve the caller's output on invalid shapes, nonfinite inputs and
  transform failures, including errors in later channels and strided sections.
- Added a native public-module test and a 466-case SciPy 1.17 reference
  suite covering 36,060 complex coefficients (maximum magnitude-scaled
  error 2.19e-12) in both checked and optimized GNU Fortran 14.2 builds.
  All original rank-one and rank-two CZT APIs are unchanged. No dependencies
  were added. The full FPM suite was not run.

## 2026-10-27: FFT DPSS spectral-concentration ratios

- Added `ratio_method='auto'|'direct'|'fft'` to `dpss`, appended after
  `solver` to preserve existing source calls. The automatic mode uses a
  zero-padded FFT autocorrelation for internal lengths of at least 256.
  The original O(k*n*n) Toeplitz calculation remains selectable as `direct`.
- The FFT variant uses a reusable power-of-two complex workspace across
  tapers, computes the full untrimmed window autocorrelation without circular
  aliasing, and uses compensated summation of the sinc-weighted lags.
  It reduces concentration-ratio work to O(k*n*log(n)) and temporary
  workspace to O(n). A small private radix-2 kernel avoids creating a
  circular module dependency with `scifort_signal_spectral`.
- The new reference suite passes 256 independent SciPy 1.17 ratio comparisons
  in checked and optimized builds, with a maximum absolute difference of
  1.21e-14 in the optimized build. The 150 native direct/FFT comparisons
  agree within 2.78e-15. Earlier 224 partial-eigensolver, 344 DPSS, and
  226 subsample-normalization SciPy tests still pass.
- On this GNU Fortran 14.2 `-O3` system, end-to-end DPSS computation
  with six tapers at length 8192 improved from 0.197 s (`direct`) to
  0.0399 s (`fft`), about 4.94x; smaller lengths show lesser gains.
  All 23 affected DPSS/multitaper checked public regressions passed after
  incremental rebuilding of the modified modules and dependents.
  The project-wide FPM suite remains unverified because FPM is unavailable.

## 2026-10-27: Selected-eigenpair DPSS eigensolver

- Added `solver='auto'|'ql'|'partial'` to `dpss` without changing existing
  positional arguments or default normalization. The `partial` solver uses
  Sturm eigenvalue bisection, pivoted tridiagonal inverse iteration and
  twice-reorthogonalized eigenvectors. Its workspace is O(n*k) rather than
  the dense QL solver's O(n*n). Automatic dispatch selects this method for
  internal lengths >= 256 and at most 32 requested tapers.
- Verified 224 independent SciPy 1.17 cases in both checked and optimized
  builds (1,149,984 coefficients, maximum absolute taper error 8.32e-10).
  Retained the 344-case existing DPSS and 226-case subsample-normalization
  comparison suites, and added native/public API regressions.
- An optimized 2048-point, six-taper local benchmark took 0.0162 CPU seconds
  with the partial solver versus 6.54 seconds for QL (about 403x faster),
  with peak process RSS near 1.5 MB versus 34 MB. Measurements are hardware-
  and compiler-specific. Exact concentration ratios continue to be computed
  by the O(k*n*n) Toeplitz autocorrelation calculation.

## 2026-10-27: Anisotropic and elliptical two-dimensional FIR design

- Extended `firwin_2d` with optional `fc2` (second-axis cutoff vector),
  `fs2`, `pass_zero2`, and `scale2` for independent separable filter design.
  Added `ellipse_axes` and `ellipse_angle` for rotated elliptical radial
  interpolation; all new options are labeled SciFort extensions.
- Preserved existing positional arguments, default windowing, circular
  output shape, and SciPy 1.17 compatibility. Invalid axis/radial option
  combinations or cutoff/radius values return an empty array and an
  invalid-input status, rather than silently ignoring arguments.
- Added native kernel/public-module tests and a reference harness with
  239 independent NumPy/SciPy cases covering 18,758 filter coefficients.
  The original 499-case algorithm comparison also remains unchanged.
  No external libraries or dependencies were added.

## 2026-10-27: DPSS subsample normalization

- Added `norm='subsample'` to the public `dpss` interface, completing the
  three SciPy 1.17 normalization conventions (`l2`, `approximate`, and
  `subsample`) without changing defaults or other result fields.
  Even-length first tapers use bandlimited half-sample peak interpolation;
  a closed-form Dirichlet-kernel identity avoids adding FFT dependencies.
  Periodic windows use their untrimmed internal length for normalization.
- Added `test/test_signal_dpss_subsample.f90` and
  `test/compare_signal_dpss_subsample_scipy.py`. The existing standalone
  reference driver now recognizes normalization code 2; the preceding
  normalization codes and calling convention remain unchanged.
- Evaluated mixed-radix FFT stage orderings as a potential optimization,
  but no alternative consistently improved mixed-radix workloads. The
  existing FFT ordering and public FFT APIs are therefore unchanged.

## 2026-10-26: Specialized mixed-radix FFT butterflies

- Replaced the per-output `r`-term summation in `src/signal_spectral.f90`
  with radix-sized decimation-in-time butterflies. Factoring the stage
  twiddle `W_n^(r*q)` from the small radix DFT removes an integer modulus
  and redundant twiddle lookups from each output. Radix 3, 5 and 7
  butterflies use symmetric-pair sums and differences and constant
  trigonometric coefficients; radix 2 uses the usual sum/difference pair.
- Both the unplanned `forward_fft` and reusable `fft_plan` use the updated
  kernel automatically on 2/3/5/7-smooth lengths. FFT signs, output order,
  normalization, public APIs, storage bounds, and plan threading rules are
  unchanged. Power-of-two transforms still use radix 2, and other large
  lengths still use Bluestein.
- Added `test/test_signal_fft_butterflies.f90` (66 native cases including
  independent direct DFT bins, impulses, complex frequency tones, inverse
  reconstruction, and noncontiguous sections),
  `tools/signal_reference/compare_signal_fft_butterflies.py` (138 SciPy
  reference cases and 206,016 coefficients), and
  `benchmarks/bench_signal_fft_butterflies.f90` for reproducible before/after
  comparisons. The previous 246-case SciPy FFT comparison also still passes.
  Across the two comparison suites, 384 SciPy cases and 527,718 coefficients
  passed in both checked and optimized builds; maximum magnitude-scaled
  discrepancy was `3.77e-14`.
- In a local GNU Fortran 14.2 `-O3` planned-transform benchmark against the
  previous implementation, length 343 improved ~2.12x, 2401 ~2.60x,
  and 5040 ~1.90x; across ten tested smooth lengths, improvement ranged
  from ~1.36x to ~2.60x. The benchmark reports best-of-three CPU times;
  speedup is workload/compiler/hardware dependent.
- Clean checked compilation of all 88 signal dependencies passed, followed
  by 36/36 targeted native regression programs. Those same 36 programs
  passed when the FFT kernel was compiled at `-O3` and linked with checked
  dependency objects. A full optimized 88-module build timed out;
  `fpm test` was not run in this environment.

## 2026-10-26: Cost-aware FFT padding for N-dimensional convolution

- `fftconvolve_nd`, the rank-2/3 `fftconvolve` overloads, `convolve_nd`
  and `correlate_nd` FFT dispatch, and `oaconvolve_nd` now evaluate supported
  2/3/5/7-smooth FFT lengths as possible alternatives to radix-2 padding.
  A conservative local cost heuristic accounts for padded tensor volume and
  the slower per-element mixed-radix implementation. Small and many ordinary
  convolutions retain power-of-two transforms; rank-one overlap-add keeps its
  established power-of-two block lengths to avoid regressions.
- This change is internal: modes, selected axes, singleton broadcasting,
  real/complex promotion, rank-15 adapters and public signatures are intact.
  Sizing arithmetic uses 64-bit temporaries and checks the default-integer
  bound when searching for a supported transform length.
- The new native regression and seven large-shape SciPy 1.17 comparisons
  cover complex convolution, full/same/valid cropping, rank-15 axes,
  unequal operand sizes, and blocked overlap-add. The seven cases contain
  305,675 output elements and differ from SciPy by at most 3.28e-13.
  Existing 1,008 SciPy convolution/correlation and 232 overlap-add
  comparisons pass. All 35 affected signal regression programs passed in
  both checked and optimized GNU Fortran 14.2 builds.
- A representative rank-two workload (shapes 101x251 and 29x59) now uses
  135x315 instead of 256x512 FFT work arrays, reducing each complex workspace
  from 131,072 to 42,525 elements (about 68%). In eight interleaved `-O3`
  trials, median total benchmark time was 0.03359 s for the previous version
  and 0.02486 s for this version (~1.35x faster); other shapes were similar
  or showed small timing differences. These are local measurements, not a
  portable speed guarantee. Benchmark source is included.
- The full `fpm test` was not run in this environment.

## 2026-10-26: Reuse FFT plans in multitaper and N-dimensional convolution

- Reused one `fft_plan` per selected transform axis inside `fftconvolve`
  (ranks 1-3) and the assumed-rank `fftconvolve_nd`/`oaconvolve_nd`
  kernels (ranks 1-15). Multiple channels/lines and overlap-add tiles share
  the same per-call plans. Rank-one overlap-add shares the plan across tiles.
- All whole-record, multichannel, harmonic, and windowed multitaper estimators
  use one FFT plan across their DPSS tapers, channels, and windows. A
  `multitaper_stream` instance now retains a plan across successive pushes and
  resets; reinitialization for another FFT length releases the old plan.
  All existing result shapes, options, and FFT normalization are unchanged.
- Added `test/test_signal_fft_plan_reuse_consumers.f90` (prime/smooth FFTs,
  batched versus streamed data, resets, changing FFT lengths, rank-15
  convolution, and invalid reinitialization) and
  `benchmarks/bench_signal_multitaper_fft_plan_reuse.f90` for before/after
  comparisons against the previous release. Corrected dependency ordering in
  the overlap-add SciPy reference script.
- In two local GNU Fortran 14.2 `-O3` runs of the end-to-end streaming
  benchmark, length 263 took 0.884/0.904 s before and 0.641/0.623 s after
  (~1.4x faster); length 315 was essentially unchanged. Speedups vary with
  FFT size, taper count, hardware, and compiler.
- Checked and optimized public-interface regressions passed for the modified
  convolution/multitaper paths. Existing independent SciPy 1.17 comparisons
  still pass: 1,008 ND convolution/correlation cases, 232 tiled OA cases,
  120 noise-floor time-frequency cases, 80 streamed multitaper cases,
  and 90 multitaper spectrogram plus 100 harmonic-test cases.
  Full `fpm test` was not run in this environment.

## 2026-10-26: Mixed-radix FFTs for composite lengths

- Added an allocation-free execution path for reusable `fft_plan` objects
  with lengths factorable entirely into 2, 3, 5, and 7. A recursive
  decimation-in-time Cooley-Tukey implementation combines small-radix
  transforms with precomputed roots and reusable workspace. The unplanned
  `forward_fft` interface selects the same algorithm automatically.
- Existing radix-2, small direct-DFT, and Bluestein paths remain available.
  Plan reinitialization and inverse transforms retain their public behavior.
  No external FFT library or Fortran dependency is introduced.
- Added `test/test_signal_fft_mixed_radix.f90`, independent SciPy 1.17
  verification `tools/signal_reference/compare_signal_fft_mixed_radix.py`,
  and `benchmarks/bench_signal_fft_mixed_radix.f90`. The 246 reference cases
  cover 321,702 complex coefficients, with magnitude-scaled error at most
  3.77e-14 in checked and optimized GNU Fortran 14.2 builds. New native
  checks cover 38 transform sizes and 114 repeated planned transforms.
- Local optimized timings on selected 7-smooth lengths from 225 through
  5040 show roughly 1.4-2.5x speedups for reusable plans relative to the
  preceding Bluestein-plan implementation. Benchmarks are platform-specific.
- Clean 88-module public builds passed 16 FFT-related regression programs
  in each checked and optimized configuration. An additional 16 spectral,
  streaming STFT, resampling, CZT, and multitaper regressions passed in each
  configuration. Full-project FPM testing was not available here.

## 2026-10-26: Reusable Bluestein FFT plans

- Added public `fft_plan` and `forward_fft_planned` to `scifort_signal`,
  preserving the original `forward_fft` entry point. Plans cache Bluestein
  chirps, transformed convolution kernels, and scratch space across calls.
  Methods support initialization, forward/inverse transforms, length inspection,
  and explicit release. Status codes cover invalid sizes and mismatched plans;
  failed initialization does not alter an existing plan or input array.
- Real batch Welch, CSD, periodogram, STFT, and spectrogram operations and
  the forward/inverse streaming STFT kernels reuse a plan across frames.
  No implicit global cache or additional dependency is introduced.
  A plan has mutable workspace; use separate plans for concurrent threads.
- Added native plan regression, independent SciPy reference driver, and
  `benchmarks/bench_signal_fft_plans.f90`. Checked and optimized FFTs pass
  198 SciPy 1.17 comparisons (297,189 coefficients, maximum magnitude-scaled
  difference 8.96e-13). In local GNU Fortran 14.2 `-O3` measurements, plan
  reuse is 1.46-1.83x faster for lengths 257-8191 than repeated unplanned
  Bluestein transforms. Timings are environment-specific.
- The 88-module checked public build and 29 signal regression programs pass.
  Complete FPM testing is not available in this environment.

## 2026-10-26: Fast arbitrary-length shared FFT (Bluestein)

- Replaced the quadratic-time direct-DFT fallback in `scifort_signal_spectral::forward_fft`
  with Bluestein chirp convolution for all non-power-of-two lengths above 16.
  The existing radix-2 kernel remains unchanged for power-of-two sizes, and
  the direct DFT remains for very short transforms. Forward sign convention
  and normalization (unnormalized) are unchanged.
- The algorithm uses two radix-2 transforms for the chirp sequences and a
  third for their convolution, with a transform length of the next power of
  two >= 2*n-1. Chirp phases are reduced modulo 2*n using 64-bit integer
  multiplication. The implementation is original and adds no dependencies.
- Added `test/test_signal_fft_bluestein.f90`, SciPy FFT comparison driver,
  and `benchmarks/bench_signal_fft_non_power_two.f90`. In local optimized
  GNU Fortran 14.2 measurements, prime-length transforms of 257, 1009, and
  4093 samples were about 36x, 227x, and 641x faster than the old direct DFT.
  Timings are hardware-dependent.
- The native 24-length regression and 116 comparisons with SciPy 1.17 pass
  in both checked and optimized builds (133,608 output coefficients; maximum
  magnitude-scaled difference 8.93e-13, lengths through 16,381). The
  88-module checked public build passes 27 neighboring signal regressions;
  the same tests pass with the changed FFT compiled at `-O3`.
  Full-project FPM testing remains unverified in the present environment.

## 2026-10-26: Callable Fourier-resampling windows

- Added `window_function` (real-valued vector callback) and
  `complex_window_function` (complex-valued callback for complex inputs) to
  Fourier `resample` (ranks 1-3) and `resample_into` (ranks 1-15). Existing
  named and FFT-bin array windows remain backward compatible.
- Callback inputs are frequencies from `fftfreq(nx)` in cycles per sample.
  Even FFT lengths put the Nyquist frequency at -0.5. Multichannel and
  rank-independent calls sample the callback once and reuse its weights.
- Real inputs fold positive/negative-frequency weights like SciPy; callbacks
  returning nonfinite values and conflicting window specifications return
  invalid-input status. Caller-output arrays are unchanged on failure.
- 450 independent SciPy 1.17 cases (12,722 values) pass in both checked and
  optimized GNU Fortran 14.2 builds; maximum absolute deviation 5.69e-14.
  The existing 900 array-window reference cases and four resampling
  regression programs pass unchanged. FPM is unavailable in this environment.

## 2026-10-26: Frequency-bin array windows for Fourier resampling

- Added trailing optional `frequency_window(:)` (real FFT-order weights) to
  `resample` real/complex rank-1..3 interfaces and `resample_into` for ranks
  1..15. Complex inputs also accept `complex_window(:)`, allowing complex-valued
  FFT-bin weights, including complex spectra in `domain='freq'` mode.
- For real inputs, paired positive/negative-frequency weights are averaged
  before resampling, as in SciPy's real-FFT path. The rank-independent real
  wrapper now preserves this behavior with arbitrary periodic named windows.
- Invalid window lengths, nonfinite weights and simultaneous named/array
  windows return invalid-input status; `resample_into` never modifies the
  output on failures. Positional calls and named-window defaults are unchanged.
- 900 independent SciPy 1.17 cases (32,664 samples) pass at both `-O0` with
  runtime checks and `-O3`; maximum absolute difference 5.76e-14. Five
  targeted public-module filtering/resampling tests pass with optimized
  changed modules, including rank-15 and strided-output regressions.
- Full `fpm test` and complete optimized dependency rebuild remain unverified.

## 2026-10-26: Windowed and streaming multitaper noise floors

- Added `multitaper_noise_floor_spectrogram` rank-1 and rank-2 real/complex
  overloads, returning spectral density, local-median background, and
  background-normalized power for each complete window and channel.
- Added `multitaper_stream%push_noise_floor` for irregular real/complex sample
  chunks, plus `mean_noise_floor` for the median-normalized **mean PSD**.
  Both share the existing stream state and reject invalid widths without
  consuming samples. The new mean result has no artificial time dimension.
- Reused the whole-record estimator's normalization rules; preserved its
  API, numerical behavior, and zero-background convention.
- Added an original Fortran regression and an independent SciPy/NumPy
  reference driver. All 120 randomized cases (47,098 PSD bins) passed in
  checked and optimized builds. All 16 targeted multitaper/DPSS public
  regression programs passed with both runtime checks and `-O3`.
- Full-project FPM testing remains unverified (FPM unavailable).

## 2026-10-26: Colored-noise multitaper background and Holm line tests

- Added `multitaper_noise_floor` with real/complex overloads and a new
  result type reporting multitaper PSD, local-median frequency-dependent
  background, and normalized power. One-sided endpoints reflect; full FFT
  spectra wrap. A zero spectrum stays zero, invalid median widths return
  invalid status, and zero-background division is handled explicitly.
- Added Holm step-down family-wise correction (`method='holm'`) to
  `multitaper_harmonic_scan`, without changing existing default `bh` or
  supported `none`/`bonferroni`/`by` methods. Holm's control is valid under
  arbitrary dependence provided the underlying F-test p-values are valid.
- Added native tests and SciPy/NumPy independent checks: 144 colored-noise
  PSD/median comparisons (10,889 bins) and 225 harmonic-scan comparisons.
  Checked multitaper standalone suite passes 15/15 programs. Both features
  re-export through `scifort_signal`. No new Fortran dependencies.

## 2026-10-26: Whole-record and streaming multitaper linear detrending

- Added trailing optional `detrend_mode='none'|'constant'|'linear'` to
  whole-record multitaper PSD (including multichannel/complex overloads),
  spectral matrix, CSD, coherence, jackknife confidence intervals, and
  harmonic F-test/scan. The option overrides the old logical `detrend`
  argument; original defaults and positional calls remain unchanged.
- Extended `multitaper_stream%initialize` with the same option. The stream
  detrends **each complete window and channel independently**; `reset()`
  retains the selected mode. Invalid mode strings return invalid status.
- Added `test_signal_multitaper_whole_stream_detrend.f90` and an independent
  SciPy 1.17 DPSS/linear-detrend/NumPy FFT reference driver and harness.
  96 randomized comparisons (67,127 spectral coefficients) passed with
  maximum absolute difference 2.54e-13, including complex and real inputs.
  Eight focused Fortran regression programs passed in checked and optimized
  GNU Fortran builds; FPM project-wide testing remains unverified.

## 2026-10-26: Improved Gustafsson solve and multitaper linear detrending

- Added pivoted complex Householder QR on the conjugate-transposed
  initial-state design matrix for underdetermined Gustafsson problems;
  returns a minimum-norm solution without forming normal equations. The
  prior Jacobi SVD and pivoted QR implementations remain fallbacks.
- SOS-Gustafsson independent state-space reference now passes 293/297
  cases, up from 289/297. Four near-singular 12th-order cases on only
  20 observations continue to report numerical failure rather than
  silently producing unreliable samples. No new external dependencies.
- Added trailing `detrend_mode='none'|'constant'|'linear'` to the windowed
  real/complex multitaper spectrogram, spectral-matrix spectrogram,
  cross-spectrogram and time-frequency coherence APIs. The option
  overrides the original logical `detrend` when present; old calls keep
  their previous defaults and results. Regression coverage includes
  real/complex signals, mixed pairs, per-channel slopes, axis=2,
  invalid arguments and baseline compatibility.
- Added `test_signal_multitaper_linear_detrend.f90` and a SciPy
  `detrend(type='linear')` / DPSS / FFT reference driver and script.
  The independent comparison passes 64 cases and 36,787 complex
  spectral-matrix coefficients (maximum absolute error 8.92e-15).
  All 23 focused checked Fortran regressions and five targeted optimized
  regressions pass; full-project FPM testing remains unverified.

## 2026-10-26: SOS-native Gustafsson initial conditions

- Implemented `sosfiltfilt_gust` directly on second-order-section states,
  without converting cascades to high-order transfer polynomials. The public
  `sosfiltfilt` rank-1/2/3 and `sosfiltfilt_nd` rank-1..15 interfaces accept
  trailing `method='gust'` and optional `irlen`, preserving the existing
  padding default and positional compatibility.
- Supports real, complex and mixed-type SOS/sample operands, selected axes,
  short samples and noncontiguous caller arrays. Invalid parameters produce
  NaNs and an invalid-input status; unreliable ill-conditioned state solves
  explicitly return the numerical-failure status.
- Added native regression and independent Python/SciPy SOS-state-space
  comparison driver; the SciPy API does *not* implement SOS Gustafsson, so
  comparisons use its SOS streaming kernel plus independent least squares.
- No new external dependencies; `fpm test` has not been run here.

## 2026-10-26: Gustafsson method for assumed-rank filtfilt_nd

- Extended `scifort_signal::filtfilt_nd` with trailing optional `method=` and
  `irlen=` arguments, preserving all existing positional calls and padded
  defaults. Gustafsson filtering now works on ranks 1-15 and all eight
  real/complex coefficient and signal combinations, with independent
  one-based filtering axes and support for strided sections.
- Gustafsson ignores `padlen`/`padtype` and permits signals shorter than
  the default pad length. Invalid methods, nonpositive `irlen`, invalid
  coefficients, axes or shapes return the existing invalid-input status;
  per-line errors clear any partially written lanes to quiet NaNs.
- Added `test/test_signal_filtfilt_nd_gust.f90` and
  `tools/signal_reference/{driver_signal_filtfilt_nd_gust.f90,compare_signal_filtfilt_nd_gust_scipy.py}`.
  In each of checked and optimized GNU Fortran 14.2 builds, 960 reference
  cases passed (48,384 samples; 240 direct SciPy 1.17 and 720 independent
  complex least-squares comparisons), with maximum absolute discrepancy
  `2.24e-16` or less. Legacy assumed-rank padding regressions also pass.
- Updated the legacy SciPy comparison script and standalone signal test
  runner to include the Gustafsson module dependency. No new external
  dependencies; full-project FPM testing remains unverified.

## 2026-10-26: Gustafsson initial states for zero-phase filtering

- Implemented `filtfilt(..., method='gust', irlen=...)` for real/complex
  rank-1, rank-2 and rank-3 signals, without modifying default padded
  forward/backward filtering or the SOS interface.
- Added a standalone `filtfilt_gust` rank-1 generic with mixed real/complex
  coefficient/signal support and a pure-Fortran minimum-norm solver for
  short/underdetermined signals, with optional impulse-response truncation.
- Added 858 SciPy/independent-complex-reference comparison cases spanning
  41,897 values (758 directly against SciPy 1.17, 100 against a separate
  complex least-squares reference); maximum absolute deviation `2.26e-8`
  in a short, ill-conditioned high-order case, and maximum scaled deviation
  `7.89e-11`. All five related Fortran tests passed in a clean checked build
  and with the modified filtering modules optimized at `-O3`.
- Added a public native regression, error statuses and documentation with
  no external dependencies. At this stage the higher-rank `filtfilt_nd`
  still used padding only; its Gustafsson support is documented above.
  The complete optimized dependency rebuild exceeded the execution limit;
  `fpm test` was not run because FPM is unavailable.

## 2026-10-26: caller-defined wavelets for find_peaks_cwt

- Added optional real `wavelet=` procedure to `scifort_signal::find_peaks_cwt`
  without changing the existing positional arguments; the callback accepts
  `(points,width)` and returns an allocatable one-dimensional real kernel.
- Added `cwt_wavelet` and published `cwt_wavelet_callback`; retained `cwt_ricker`
  and its original numerical behavior as the default. Reject invalid kernel
  lengths and nonfinite samples with the existing empty-result/status contract.
- Added custom-kernel native regression and SciPy 1.17 reference drivers.
  Both checked and optimized kernels matched 600 cases covering 475,752 CWT
  values and 1,561 peak indices; maximum error <= 2.49e-14.
  The original 300-case SciPy CWT suite still passes.
- Updated the original native CWT peak test to allow a one-sample
  ambiguity for nearly tied symmetric peaks under different optimization
  levels; the ridge detection algorithm itself is unchanged.
- Updated signal API documentation, README, roadmap and source provenance.
  No new external dependency. Complex wavelet callbacks are not yet supported.

## 2026-10-25: open-interval peak property requests

- Extended `scifort_signal::find_peaks` with five optional `logical` keywords:
  `height_open`, `threshold_open`, `prominence_open`, `width_open`, and
  `plateau_size_open`. True requests the corresponding properties without
  imposing a new filter, analogous to SciPy `(None, None)` intervals.
  Existing scalar/per-position constraints and positional calls still work;
  `width_open` also returns prominence details.
- Added a public-interface Fortran test covering 32 flag combinations,
  bounds, plateaus, strided inputs and empty signals. In both checked and
  `-O3` peak-kernel builds, 2,400 independent SciPy 1.17 comparisons
  passed (157,239 reported values; max error `2.84e-14`).
- Updated signal documentation and README. No new dependencies or changed
  public derived types.

## 2026-10-25: position-dependent find_peaks intervals

- Extended the existing `scifort_signal::find_peaks` function with ten optional
  per-sample `real(dp)` lower/upper constraints for height, threshold,
  prominence, width, and plateau size. Bounds are indexed at the peak center,
  and scalar plus per-sample constraints can be combined without changing the
  previous API. Plateau bounds allow fractional comparisons.
- Validates dynamic bound lengths, finiteness, and nonnegative plateau sizes;
  returns the existing invalid-input status and empty result arrays.
  Strided input and bounds are supported.
- Added native public-module regressions, a standalone reference driver and
  2,200 seeded independent SciPy 1.17 cases comparing 55,439 values with
  maximum absolute error `2.84e-14`. Existing scalar peak behavior is retained.
  Also updated the standalone signal test runner to compile previously
  missing signal convolution, time-frequency and multirate dependencies. No new library dependencies.

## 2026-10-25: N-dimensional direct convolution and correlation

- Added `convolve_nd`, `correlate_nd`, and `convolution_shape_nd` in the public
  `scifort_signal` module. Fortran 2018 caller-output routines accept real,
  complex and mixed inputs of rank 1-15, independently selected one-based
  axes, singleton broadcasting and SciPy-compatible `full`/`same`/`valid`
  cropping. Output is unchanged on invalid input. All rank-1/2/3 APIs remain.
- Implemented independent direct multidimensional accumulation and FFT-based
  convolution with conjugate reversal of selected correlation axes. A small
  workload heuristic selects a method for `method='auto'`. Both FFT and direct
  modes are explicitly selectable.
- Added 21 native interface/edge cases and 1008 independent SciPy 1.17
  comparisons (67,609 output elements, max absolute error 1.16e-14) passing
  in both checked and optimized GNU Fortran 14.2 builds; 10 related checked
  convolution/overlap-add tests passed. See `docs/signal.md`.
  No additional external dependencies.

## 2026-10-25: multichannel streaming STFT, spectrogram and inverse STFT

- Added public `stft_stream_multi` and `istft_stream_multi` state types, with
  separate multichannel init/push/finish APIs for STFT, power spectrogram,
  and inverse STFT. One user-facing state handles synchronous real or complex
  `(sample,channel)` input and `(frequency,channel,frame)` output.
- Fixed persistent memory scales with channel count and window/FFT size,
  not stream length. Preserved single-channel windows, timing, scaling,
  padding, error codes, and one-sided versus two-sided conventions.
- Added 74 deterministic Fortran cases, irregular and strided chunks,
  multichannel reconstruction, NOLA and invalid-input recovery, and a
  10,000-sample/four-channel stress test. Added 256 independent SciPy 1.17
  comparisons: 12,304 forward bins (max error `1.33e-15`) and 4,352 inverse
  samples (max error `4.78e-14`). See `docs/signal_stream_multi.md`.
  No new external dependencies.

## 2026-10-25: streaming forward STFT and spectrogram

- Added public `stft_stream` state with `stft_stream_init`, `stft_stream_push`,
  `stft_stream_finish`, `stft_stream_frequencies` and complementary
  `spectrogram_stream_init`, `spectrogram_stream_push`,
  `spectrogram_stream_finish` for incrementally analyzing arbitrary
  real/complex signal chunks with bounded persistent FFT/ring storage.
- Supports separate spectral normalizations, user windows, configurable
  detrending, one-/two-sided spectra, time/frequency coordinates, and
  delayed edge and hop padding until finalization. Rejected nonfinite chunks
  or invalid reinitialization do not corrupt the active state.
- Added 1,080 batch-versus-stream comparisons, two end-to-end streaming
  forward/inverse reconstructions, 20,000-sample stress check, and 384
  independent SciPy 1.17 reference cases (49,698 complex/real spectral bins,
  maximum absolute error `4.44e-15`). Nine related checked signal tests and
  an optimized new-module run passed. See `docs/signal_stft_stream.md`.
  No new dependencies.

## 2026-10-25: incremental inverse STFT with bounded memory

- Added `istft_stream` and `istft_stream_init` / `istft_stream_push` /
  `istft_stream_finish` to the public `scifort_signal` API. Reconstructs
  streaming one-sided real or two-sided complex STFT frames without a
  full-frame or full-output allocation; keeps only one FFT workspace and
  an overlap-add ring buffer, including correct trimming of both edges.
- Supports custom and named windows, FFT padding, PSD/spectrum scaling,
  independently emitted sample times and sticky noninvertible diagnostics.
  Rejects malformed or nonfinite frames without advancing the stream.
- Added 960 deterministic comparisons with batch `istft`, a 20,000-frame
  fixed-memory stress test, and 1,350 independent SciPy 1.17 cases
  (17,100 samples, maximum absolute error 4.41e-13). See
  `docs/signal_istft_stream.md`. No new library dependencies.

## 2026-10-25: inverse STFT with independently selectable axes

- Extended `istft_into` with optional `time_axis` and `freq_axis` arguments
  accepting arbitrary independent positions in rank-2 through rank-15 input
  arrays. Retained the original `axis` argument for the output sample
  dimension, including previous positional and keyword calls; explicit
  frequency and output axes must agree after removing the input time axis.
- Added general column-major stride mapping for frequency/time permutations
  and channels; retained atomic writes on input errors and noninvertible
  channels. Verified rank-2/rank-15 inputs and strided sections.
- Added 1,296 independent SciPy 1.17 comparisons over all six axis
  permutations, 32,832 reconstructed samples (max error 1.83e-14), and a
  standalone reference driver. Eleven related public-interface signal
  regressions also passed in checked and optimized GNU Fortran builds.
  No new FPM dependencies.

## 2026-10-25: rank-independent inverse STFT

- Added caller-output `istft_into` for complex STFT inputs of rank 2-15,
  reconstructing real or complex signals of rank 1-14. The frequency axis
  replaces the sample axis and input time is last, matching `stft_into`.
- Preserved the existing `istft` weighted overlap-add kernel and its window,
  overlap, FFT, scaling, one-/two-sided, and boundary options. Accepts
  noncontiguous array sections and commits an output only when every channel
  succeeds; invalid and noninvertible cases leave results untouched.
- Added public-module native tests through rank 15, a standalone SciPy
  comparison driver, and 240 independent SciPy 1.17 inverse-STFT cases
  covering 10,008 reconstructed samples (maximum absolute error 1.76e-14).
  No new dependencies. See `docs/signal_timefreq_nd.md`.

## 2026-10-25: rank-independent STFT and spectrogram

- Added Fortran 2018 `stft_into` and `spectrogram_into` interfaces with
  selectable one-based/negative axes for real and complex rank-1–14 input,
  yielding rank-2–15 time-frequency arrays. Frequency replaces the signal
  axis and time is last. The existing legacy one-dimensional kernels are
  reused unchanged.
- Added spectrogram `mode='psd'`, `'complex'`, `'magnitude'`, `'angle'`,
  and `'phase'` (frequency-axis phase unwrapping) with typed real/complex
  output arguments. Invalid inputs leave output arrays unchanged.
- Added native regression tests, an idempotent reference-driver generator,
  and 396 SciPy comparisons covering 727,408 values and rank-14 signals
  (excluding undefined phase at essentially zero spectral magnitude).
  See `docs/signal_timefreq_nd.md`. No additional FPM dependencies.

## 2026-10-25: rank-independent spectral estimators

- Added Fortran 2018 caller-output `periodogram_into`, `welch_into`,
  `csd_into`, and `coherence_into` supporting rank-1–15 real and complex
  inputs and selectable axes. Cross-spectra also support mixed types.
- Preserved existing 1-D kernels, periodogram and Welch averaging options,
  frequency conventions, and no-partial-write semantics.
- Added public-interface native regressions, 396 SciPy numerical comparisons
  covering 446,736 spectral coefficients, and reproducible driver generators.
  See `docs/signal_spectral_nd.md`. No new Fortran dependencies.

## 2026-10-25: rank-15 assumed-rank multirate APIs

- Added `resample_into`, `resample_poly_into`, `upfirdn_into`, and
  `decimate_into` with rank-1 through rank-15 assumed-rank input/output
  arguments and caller-provided destination arrays. These reuse the existing
  1-D kernels, preserve original allocation-return function interfaces, and
  support arbitrary axes and real/complex data (including complex taps).
- Added complete-output transaction semantics: invalid axes, bad dimensions,
  or a failed lane return an error without modifying the destination. Strided
  array sections are supported. Both real and complex high-rank paths use one
  shared Fortran-order lane traversal rather than rank-specific DSP code.
- Added an idempotent generator for portable `select rank` pack/unpack clauses,
  a public-module native regression through rank 15, and 738 randomized
  SciPy 1.17 comparisons covering 27,948 values (max absolute error ~8.7e-15).

## 2026-10-25: complex FIR polyphase resampling and padding modes

- Extended `resample_poly` with custom complex FIR `taps=` for real and
  complex 1-D input and selectable-axis rank-2/3 input. Real signals
  with complex taps return complex arrays. Existing real-tap interfaces
  are unchanged, and no additional FPM dependencies are needed.
- Added `padtype='median'`, `'minimum'`, and `'maximum'` to polyphase
  resampling (including complex lexicographic statistics), and added
  `'antisymmetric'` and `'antireflect'` boundary modes to `upfirdn`
  and `resample_poly`.
- Added source-generator support for four complex-tap rank-2/3 overloads,
  native regressions, and a 184-case SciPy 1.17 numerical reference suite.

## 2026-10-25: multidimensional multirate resampling

- Added rank-2 and rank-3 real/complex, axis-selectable overloads of `resample`,
  `resample_poly`, `upfirdn`, and `decimate` to the existing signal modules.
  `upfirdn` accepts all four real/complex tap and signal combinations.
  Positive axes are one-based; negative axes count from the end. Default axes
  correspond to SciPy's first/last-axis conventions for each function.
- The wrappers reuse SciFort's existing 1-D Fourier, polyphase, and decimation
  kernels without additional dependencies. Existing rank-1 calls are unchanged.
- Added an idempotent generator, a native regression test, and an independent
  53-case SciPy reference comparison (6,105 complex values, maximum absolute
  discrepancy about `4.2e-15`).

## 2026-10-25: continuous-time Gramians and balanced truncation

- Added `clti_system%gramians(tol,max_iter)` for complex MIMO continuous-time
  controllability and observability Gramians. A scaled Cayley transform maps
  the continuous Lyapunov equations exactly onto SciFort's discrete doubling
  solver, preserving the two Gramians. The continuous H2 norm is finite for
  strictly proper systems; a nonzero direct feedthrough produces infinity.
- Added `clti_system%balanced_truncate(order,tol)` with continuous-time
  square-root balancing and the standard twice-the-Hankel-tail bound.
  Exposed the discrete balancing projections in `dlti_balanced_result` so
  the reductions apply directly to the original analog A/B/C matrices.
  Also exposed continuous-time controllability/observability matrices.
- Added native regression tests and 180 randomized SciPy/NumPy continuous
  Lyapunov/SVD comparisons, including complex nonnormal MIMO state matrices,
  zero-state realizations, full-/zero-/reduced-order models, and direct
  feedthrough. Fixed the `signal_clti` module ordering in the standalone
  signal runner and added a `--match` option for targeted tests.
- No new dependencies. These are SciFort extensions rather than new public
  functions in `scipy.signal`.

## 2026-10-25: balanced truncation and continuous-time LTI objects

- Added `dlti_system%balanced_truncate(order,tol)`, with descending Hankel
  singular values, a twice-the-tail error estimate, zero-state handling,
  complex MIMO state-space reduction, and checks for unstable or numerically
  rank-deficient systems. This is a SciFort extension, not an existing
  `scipy.signal` function. A local complex Hermitian Jacobi eigensolver
  avoids adding a LAPACK/FPM dependency; very small Hankel singular values
  can be less accurate because the SVD uses a normal matrix.
- Added `clti_system` with real/complex analog TF, ZPK, and SS constructors,
  SISO/MIMO frequency and transient methods, conversion, state-space system
  composition and feedback, numerical minimal realization, stability checks,
  and continuous-to-discrete conversion. The implementation reuses SciFort's
  existing simulation, frequency-analysis, and algebraic interconnection code.
- Added native checked and optimized regression tests and independent
  SciPy/NumPy comparisons of 100 random complex MIMO balanced reductions and
  100 random analog systems. No additional external dependencies.

## 2026-10-25: discrete-LTI feedback, Kalman reduction and Gramians

- Extended `dlti_system` with block-diagonal `append`, one-based `select_io`,
  and MIMO positive/negative `feedback`, including nonzero direct feedthrough
  and an explicit numerical error for singular algebraic loops.
- Added `controllability` and `observability` matrices and `minreal` to remove
  numerically unreachable and unobservable modes with reorthogonalized,
  column-pivoted Gram-Schmidt projections (not balanced truncation).
- Added `gramians`: stable discrete Lyapunov solutions by Smith doubling and
  the MIMO H2 impulse-response norm. Unstable, nonconvergent and invalid inputs
  produce error status; no additional FPM dependencies.
- Added standalone Fortran regressions and SciPy 1.17 numerical comparisons
  for real/complex MIMO systems, hidden unstable states, arbitrary feedthrough,
  and zero-state systems. SciPy is used for testing only.

## 2026-10-25: discrete-time LTI system objects and composition

- Added `dlti_system` with `dlti_tf`, `dlti_zpk`, and `dlti_ss` constructors (real and complex), plus the `dlti` shorthand for TF and state-space forms. A distinct ZPK constructor avoids ambiguous Fortran overloads when the third argument is a real gain versus a sample interval.
- Added `to_tf`, `to_zpk`, `to_ss`, `freqresp`, `bode`, `freqresp_mimo`, `bode_mimo`, `output`, `step`, `impulse`, `ninputs`, `noutputs`, `nstates`, and `stable` type-bound methods. Reuses existing discrete-time simulation, frequency-response and polynomial conversion kernels. State-space MIMO simulation and frequency response require no transfer-polynomial expansion.
- Added state-space series (`system%series(upstream)`) and parallel (`system%parallel(other)`) composition, supporting nonzero feedthrough and MIMO systems. Enforces identical sample intervals and compatible channel dimensions.
- Added native regression tests, a batch reference driver, and 168 deterministic SciPy system comparisons plus 110 SciPy MIMO series/parallel comparisons (54,618 checked values in total), with no additional external dependencies.

## 2026-10-24: reusable CZT and zoom-FFT plans

- Added `czt_plan` and `zoom_fft_plan` to `scifort_signal` for repeated chirp-Z and zoom-FFT transforms. Plan initialization caches the unit-contour Bluestein convolution kernel; off-unit spirals retain numerically safer direct Horner evaluation with cached evaluation ratios.
- Added typed real/complex transforms for vectors and rank-two batch matrices, with one-based selectable axes; complex spiral points and zoom frequency coordinates; plan validation and reinitialization statuses.
- Added `test_signal_czt_plans.f90` and deterministic SciPy `CZT`/`ZoomFFT` comparison driver and script. No new third-party dependencies.

## 2026-10-24 - Rank-four/five ShortTimeFFT coverage

- Extended `short_time_fft` with real/complex rank-four and rank-five inputs,
  selectable signal axes, detrending, auto- and cross-spectrograms, and inverse
  synthesis from rank-five and rank-six coefficients.
- Added exported rank-five/six STFT and spectrogram result types, reusing the
  existing one-dimensional transform kernel for every signal in the batch.
- Added native roundtrip tests with arbitrary frequency/time axis permutations,
  invalid axis handling, and a Python/SciPy 1.17 numerical reference harness:
  108 cases checking 208,000 coefficients, with maximum error 8.15e-15.
- No new FPM dependencies; rank-six and higher input remain future work.

## 2026-10-24 - Multidimensional ShortTimeFFT and REML smoothing

- Extended `short_time_fft` to rank-two and rank-three real/complex input,
  with configurable transform axis; rank-three and rank-four STFT output
  types preserve SciPy's frequency-axis/time-last layout.
- Added full forward, detrended, inverse, spectrogram, and cross-spectrogram
  overloads for each rank, with arbitrary distinct inverse frequency/time
  axes. Input and output axis indices are one-based when positive or Python-
  compatible when negative. Kept the previous rank-one interfaces intact.
- Added optional automatic REML penalty selection to `whittaker_henderson`,
  with log-penalty optimization and a shared banded Cholesky determinant/solver
  used for evaluating the profiled likelihood. Fixed-penalty calls unchanged.
- Added native and public-interface regressions and numerical comparisons:
  92 multidimensional STFT cases (62,091 coefficients, max error 6.8e-15) and
  27 REML cases (972 fitted samples, max fitted error 5.1e-8).
- Fixed the standalone signal-suite dependency order so the spectral and
  FFT-convolution modules compile before the convolution dispatcher.
- No new FPM dependencies.

## 2026-10-24 - ShortTimeFFT numerical methods and compatibility

- Added `short_time_fft`, a configurable real/complex centered STFT plan:
  forward and inverse transforms, four FFT layouts, magnitude/PSD scaling,
  arbitrary zero-based time-slice subsets, phase offsets, four padding modes,
  and constant/linear pre-window detrending.
- Added auto-/cross-spectrograms; canonical-dual, from-dual and
  equal-dual (tight/unitary) window constructors. All use existing FFT
  and spectral-window routines; no new FPM dependencies.
- Exported the new types, statuses and methods through `scifort_signal`.
  Added `freqz_sos` as a public alias of `sosfreqz`.
- Added native, public-interface and optional SciPy-reference regressions:
  312 SciPy 1.17 cases, 44,058 complex coefficients, and approximately
  2.1e-14 maximum coefficient difference in checked and optimized builds.
  Another 96 tight-window reference cases agree to roughly 2.7e-15.

## 2026-10-24 - More signal processing: smoothing, normalization, FFT selection

- Added `whittaker_henderson` for fixed-penalty, weighted difference-order
  smoothing. A banded Cholesky solver uses O(n*order) storage and
  O(n*order**2) time, including interpolation across zero-weight samples.
- Added `normalize` for real, complex, and mixed numerator/denominator
  polynomials, with monic denominator normalization and SciPy-compatible
  leading-numerator-column removal, including multi-output transfer functions.
- Added `abcd_normalize` to infer state-space dimensions, validate compatible
  A/B/C/D matrices, and supply omitted zero matrices. It accepts Fortran
  real, complex, and integer rank-2 arrays and returns the existing state-space
  result type.
- Added `choose_conv_method` for real/complex 1-D inputs, with direct/FFT
  cost selection and optional CPU-time measurement. `convolve` and `correlate`
  now support `method='fft'` and dispatch `method='auto'` to the estimated
  faster method. Removed the old FFT/convolution-module dependency cycle.
- Added native regressions and an executable SciPy/sparse-solve reference
  harness. The 615-case cross-language suite checks 57,518 numerical values.
  All seven selected existing convolution and overlap-add test programs pass.
- No new project dependencies; SciPy is used only in optional reference tests.

## 2026-10-24 - Seven additional scipy.signal algorithms

- Added `closest_stft_dual_window`, supporting optimized scaled/unscaled
  real and complex analysis/synthesis windows with invertibility checks.
- Added `firwin_2d` with separable and circular designs, common windows,
  scalar/multiband cutoffs, and window-shape options.
- Added `symiirorder1` and `symiirorder2` for mirror-boundary IIR smoothing,
  using the existing spline recurrences, for real/complex rank-1 and
  batched rank-2 inputs.
- Exposed `findfreqs`, `unique_roots`, and `czt_points`, reusing SciFort's
  analog-root solver and chirp-Z transform conventions.
- Added full public `scifort_signal` exports, two native regressions,
  499 SciPy 1.17 comparisons (11,610 values, maximum difference 2.49e-14),
  and documented SciPy compatibility quirks. No new dependencies.

## 2026-10-24 - Bessel high-order roots and analog/digital IIR coverage

- Replaced cancellation-prone reverse-Bessel-polynomial root evaluation with
  scaled complex modified-Bessel-K evaluation and a guarded Aberth-Ehrlich
  iteration. Extended `besselap` from order 13 through order 20; beyond 20
  remains explicitly unsupported in binary64.
- Added the `bessel` filter-design generic with `phase`, `delay`, and `mag`
  normalizations, scalar/vector cutoffs, analog/digital output, low/high/
  bandpass/bandstop transforms, and digital second-order sections. Added
  `iirfilter(ftype='bessel')` dispatch and the `scifort_signal` export.
- Extended analog `iirfilter_zpk` and `iirfilter_sos` to Bessel/Thomson
  so high-order analog filters can avoid unstable polynomial evaluation.
- Added 108 full-design comparisons to SciPy and extended existing prototype
  and bilinear/elliptic reference checks to higher-order Bessel roots.
- No new FPM dependency; reuse of SciFort's complex modified-Bessel function.

## 2026-10-24 - bilinear digital transforms and elliptic prototypes

- Added `bilinear` and `bilinear_zpk` for converting analog transfer
  functions and ZPK filters to digital filters without prewarping.
- Exposed the normalized elliptic lowpass prototype as `ellipap`, reusing
  SciFort's existing elliptic design core and common prototype result type.
- Extended `besselap` double-precision support through order 13, retaining
  the order guard for less reliable higher-order reverse-polynomial roots.
- Verified 421 new reference cases (6,074 output values) against SciPy and
  independent complex-gain identities, and all 122 standalone signal tests
  with checked `-O0` and optimized `-O3` GNU Fortran builds.
- Updated `scifort_signal` exports, native and SciPy tests, documentation,
  provenance, and standalone signal test-runner dependency order; no new
  FPM dependencies.

## 2026-10-24 - analog prototype ZPK functions

- Added `buttap`, `cheb1ap`, `cheb2ap`, and `besselap` with common
  `analog_prototype_result` and status codes, exported by `scifort_signal`.
- Added Bessel `phase`/`delay`/`mag` normalization, bounded root solving,
  native tests, and SciPy prototype/frequency-response comparisons.
- Updated the standalone signal test runner, API documentation, and code
  provenance notes. No added dependencies.

## 2026-10-24 - analog lowpass frequency transformations

- Added eight analog filter-frequency transformation procedures: `lp2lp`,
  `lp2hp`, `lp2bp`, `lp2bs`, and their `_zpk` variants.
- Implemented direct rational-polynomial composition for real and complex
  transfer functions, without numerical root-finding, and algebraic ZPK maps
  including complex gains and conjugate-paired poles.
- Added native status and invalid-input tests and 328 SciPy comparison cases
  (4,256 complex output elements). Exported the procedures through
  `scifort_signal`, with updated documentation and reference scripts.
- Updated the standalone signal-suite module dependency order to include the
  ndimage modules used by its spline-related tests; no new FPM dependencies.

## 2026-10-24 - affine-parameter Hessians and displacement curvature

- Added `evaluate_spline_plan_affine_with_hessian_nd` for real and complex
  prepared spline images. Returns first derivatives and complete
  matrix-matrix, matrix-offset, and offset-offset parameter Hessians.
- Added flattened and grid-shaped
  `evaluate_spline_plan_displacement_with_hessian_nd` interfaces, including
  optional affine base transforms and rank-15 flattened outputs.
- Reused the existing analytic coordinate Hessian in bounded batches;
  introduced explicit output-shape and parameter-quadratic resource checks.
- Added new native and SciPy reference tests, public API checks, and docs.
  No new dependencies or breaking interface changes.

## 2026-10-24 - spline Hessians and shaped displacement gradients

- Added `evaluate_spline_plan_with_hessian_nd` (real/complex), computing
  interpolated values, first partial derivatives, and all second partial
  derivatives together for cached spline plans, orders 0-5 and ranks 1-15.
- Added `evaluate_spline_plan_grid_with_hessian_nd` with two leading
  coordinate axes on the Hessian tensor; shaped destinations support ranks
  1-13 and noncontiguous array sections.
- Added `evaluate_spline_plan_displacement_grid_with_gradient_nd` with
  grid-shaped per-displacement derivatives for output ranks 1-14.
- Preserved all eight boundary modes, including singleton-axis
  `grid-constant` corrections, plus existing public interfaces and dependencies.
- Added a native regression and independently differentiated SciPy tests.

## 2026-10-24 - shaped interpolation gradients and affine-parameter derivatives

- Added `evaluate_spline_plan_grid_with_gradient_nd` to return sampled values
  and source-coordinate gradients in caller-shaped destination arrays, including
  real/complex plans, noncontiguous arrays, and source rank 15.
- Added `evaluate_spline_plan_affine_with_gradient_nd`, which evaluates an
  affine warp and analytically differentiates it with respect to every matrix
  and offset parameter. The source-coordinate derivative is reused in bounded
  batches with a chain-rule transformation, supporting rank-changing affine
  maps and output ranks 1-15.
- Exposed both APIs through `scifort_signal` and `scifort_ndimage` and added
  native and independent SciPy numerical comparisons. Existing APIs and
  dependencies are unchanged.

## 2026-10-24 - analytic spline and displacement sensitivities

- Added `evaluate_spline_plan_with_gradient_nd` for first derivatives of
  cached spline interpolation with respect to source coordinates. Value and
  gradient use one support traversal, real and complex overloads, and ranks
  1-15. All eight spline modes, orders 0-5 and prefiltered/unfiltered plans
  are supported, including grid-constant singleton axes.
- Added `evaluate_spline_plan_displacement_with_gradient_nd` to evaluate a
  displacement-field warp and its gradient with respect to per-pixel offsets
  in bounded batches, with optional affine mapping and rank-changing output.
- Added two native regression programs, independent SciPy central-difference
  checks and a reproducible performance benchmark. New APIs are re-exported
  via `scifort_signal` and `scifort_ndimage` without new dependencies.
- Order-one spline derivatives at nonsmooth knots are not mathematically
  unique; the previously known short-line `reflect` prefilter difference
  from SciPy remains.

## 2026-10-24 - cached displacement-field spline warps

- Added `evaluate_spline_plan_displacement_nd` for flat real displacement
  fields, with real/complex spline plans, optional affine base mappings,
  rank-changing transforms, selected block size, and rank-15 output support.
- Added `evaluate_spline_plan_displacement_grid_nd` for shaped displacement
  tensors whose leading dimension contains source-coordinate components;
  output ranks 1-14 are supported by this convenience interface.
- Both operations reuse spline coefficients and bounded coordinate batches,
  retain the earlier eight boundary modes and spline orders 0-5, and validate
  dimensions, nonfinite fields, and affine arguments.
- Added native and independent SciPy reference tests, public signal/ndimage
  exports, documentation, and provenance. No new dependencies.

## 2026-10-24 - vectorized nonlinear coordinate mapping for spline plans

- Added `evaluate_spline_plan_geometric_batch_nd` and the
  `ndimage_coordinate_batch_map` callback interface. A callback receives a
  block of one-based output coordinates and returns one-based source positions.
- Real/complex interpolation plans, rank-changing mappings, ranks 1-15,
  noncontiguous outputs, all spline orders and boundary modes are supported.
  Existing scalar callbacks and public interfaces are unchanged.
- Added a native regression, 19,008 independent SciPy comparison values,
  a reproducible CPU-time benchmark, and API/provenance documentation.
  The callback workspace remains bounded by the selected batch size;
  no new FPM dependencies were added.

## 2026-10-23 - streaming affine and geometric spline-plan transforms

- Added `evaluate_spline_plan_affine_nd` and
  `evaluate_spline_plan_geometric_nd` to sample reusable spline plans for
  arbitrary nonseparable coordinate mappings in bounded batches.
- Real and complex overloads support output ranks 1-15 (including differing
  source/output ranks), noncontiguous output arrays, all spline modes/orders,
  and configurable batch sizes, without a full coordinate-table allocation.
- Added checked native and SciPy reference tests, a reproducible repeated
  affine-transform benchmark, public exports in signal/ndimage, and API notes.
  Existing spline and geometric-transform APIs remain unchanged; no new FPM
  dependencies.

## 2026-10-23 - separable tensor-grid spline interpolation

- Added `spline_axis_grid` and `evaluate_spline_plan_tensor_nd` for Cartesian-product
  coordinate grids. The evaluator reuses prepared coefficients and contracts one
  axis at a time, avoiding repeated full tensor-kernel evaluation and full
  coordinate meshes. Exposed through `scifort_signal` and `scifort_ndimage`.
- Supports real/complex data, spline orders 0-5, all eight boundary modes,
  noncontiguous output sections, and ranks 1-15. Added an eight-million-element
  workspace guard and explicit validation for grid dimensions and coordinates.
- Added native regression tests and 22,464 real/imaginary SciPy comparisons;
  the previous short-line `reflect` prefilter difference remains documented.
- A reproducible 180x160 to 210x180 order-5 benchmark measured a 31.9x
  improvement for four regular-grid queries using the same cached spline plan;
  the performance gain is workload-dependent. No new FPM dependencies.

## 2026-10-23 - reusable spline plans and mixed-type vectorized callbacks

- Added `spline_interpolation_plan`, `prepare_spline_plan_nd`, and
  `evaluate_spline_plan_nd`/`evaluate_spline_plan_grid_nd`. Repeated coordinate
  evaluations reuse prefiltered B-spline coefficients for the same input;
  real/complex and shaped-coordinate queries retain rank-15 source support.
  This is a SciFort extension, not a direct SciPy API name.
- Extended scalar and multi-valued `vectorized_filter_nd` callbacks to
  real-to-complex and complex-to-real outputs, preserving active-footprint
  order, bounded batches, selected axes, and `valid` output shapes.
- Added independent SciPy comparison drivers, native tests, and a CPU-time
  benchmark of cached versus repeated direct spline interpolation. No new
  FPM dependencies or ABI changes.

## 2026-10-23 - batched vectorized ndimage filtering

- `scipy.ndimage`: added `vectorized_filter_multi_nd` for multi-channel real/complex
  neighborhood callbacks, including `valid` output, selected axes, bounded
  batches, channel-last outputs and rank-15 flat output. The original scalar
  batched filter now shares an interior-window fast gather path. Added native
  and SciPy 1.17 reference comparisons; no new dependencies.

- Added `vectorized_filter_nd` for rank-1-to-rank-15 real and complex arrays,
  passing rectangular or sparse neighborhoods to user callbacks in bounded
  batches, with configurable `batch_size` and a two-million-window-element cap.
- Implemented selected-axis compact footprints, per-axis boundaries, shifted
  origins, NaN-aware callbacks, and `mode='valid'` output shape reduction.
- Re-exported `map_coordinates_grid_nd` through `scifort_ndimage`.
- Added a dedicated SciPy 1.17 reference driver and native Fortran tests for
  complex inputs, noncontiguous output sections, rank-15 arrays, invalid
  arguments, and callback batching. No new runtime dependencies.

## 2026-10-23 - linear-time periodic splines and shaped coordinate grids

- Replaced the quadratic Fourier-series periodic spline prefilter with an
  O(n) causal/anticausal B-spline pole solver for orders 2-5, real and complex.
- Added optional `mode=` to the public all-axis and axis-only spline filters,
  retaining the existing whole-sample mirror default and call signatures.
- Added `map_coordinates_grid_nd` with SciPy-shaped real/complex coordinate
  grids, source ranks up to 15, one-based coordinates, and default prefiltering.
- Verified 15,428 SciPy coefficient/interpolation values and exact original-knot
  reconstruction on short periodic and reflected inputs.
- Kept exact finite-line `reflect` reconstruction; SciPy 1.17's short-array
  recursive prefilter has a documented numerical difference.
- Added native rank-15, strided-array and invalid-argument tests. No new FPM
  dependencies.

## 2026-10-23 - geometric callback transforms and indexed label callbacks

- Added `geometric_transform_nd` for `real(dp)` and `complex(dp)` images with
  one-based caller mappings, different source/output ranks, SciPy-compatible
  interpolation defaults, and a coordinate-table resource guard.
- Reworked `labeled_comprehension_nd` to group pixels once rather than scanning
  every requested label across all pixels, preserving C-order and duplicate
  callback invocations.
- Added a SciFort-specific `axes(:)` option to `label_nd` for logical, integer,
  and real foreground masks, retaining globally distinct component numbers.
- Added native high-rank/invalid-input regressions and SciPy comparisons for
  arbitrary mappings and independent connected-component slice partitions.
- Added a reproducible indexed-versus-repeated-label benchmark (13.1x on one
  80,000-pixel/1,200-label `-O2` test run; results depend on workload).
- No new package dependencies.

## 2026-10-23 - ndimage complex derivative callbacks, selected-axis flood fill, fast region lookups

- Added complex `generic_gradient_magnitude_nd` and `generic_laplace_nd`
  overloads with `ndimage_complex_derivative_func` callbacks, complex `cval`,
  and selected-axis boundary modes. Complex gradients follow SciPy's
  `sqrt(sum(derivative**2))` convention.
- Added optional `axes(:)` to `binary_propagation_nd` and
  `binary_fill_holes_nd`, using compact structuring elements and boundary
  definitions restricted to selected axes. Prior calls are unchanged.
- Reduced requested-label search time for region sums, means, variances,
  centers, counts, histograms, extrema, and bounding boxes from O(N*K) to
  O(N*log(K)+K*log(K)), preserving unsorted and duplicate requested IDs.
- Added SciPy comparison drivers, checked native rank-15 regressions,
  `benchmarks/bench_ndimage_label_lookup.f90`, and updated `docs/ndimage.md`.
  No new FPM dependencies.

## 2026-10-23 - ndimage selected-axis grayscale and generic filtering

- Added `axes(:)` to all eight grayscale morphology procedures with compact
  footprints and non-flat structures, and added optional per-axis `modes(:)`.
  Preserved scalar modes and the existing rank-1-through-rank-15 interfaces.
- Added `axes(:)` to real and complex `generic_filter_nd`, including compact
  footprints, shifted origins, and mixed per-axis boundary modes.
- Reproduced SciPy's distinct axis-order rules for separable grayscale sizes,
  nonseparable grayscale footprints, and C-order callback neighborhoods;
  even-width dilation origin shifts are preserved for reordered axes.
- Added 3,840 SciPy/independent-reference scalar comparisons and a rank-15
  checked/optimized regression covering noncontiguous output and bad inputs.
- No new dependencies.

## 2026-10-23 - ndimage selected-axis rank filters and linear-time uniform filtering

- Added optional `axes(:)` to median, rank, and percentile filters, including
  selected-axis footprints, shifted origins, compact widths, and mixed modes.
  Preserved SciPy's ascending-axis stencil assignment even for reordered axes.
- Replaced the O(n*w) box-filter scan with an O(n+w) running sum in
  real/complex `uniform_filter_nd` and `uniform_filter1d_nd`, retaining their
  boundary modes, large-window validation, and origin conventions.
- Added a 10,740-value SciPy comparison suite and a native rank-15/error-path
  regression. All 15 native `ndimage` tests passed at `-O0` and `-O1`.
- Kept the existing FPM dependencies unchanged.

## 2026-10-23 - ndimage axis subsets and fast sliding extrema

- Accelerated rectangular `minimum_filter_nd`, `maximum_filter_nd`,
  `minimum_filter1d_nd`, and `maximum_filter1d_nd` to O(n+w) monotone-deque
  line processing, including all existing boundary modes and origins.
- Added optional `axes(:)` to real/complex Gaussian and uniform filters,
  and to rectangular minimum/maximum filters. Vector parameters and
  per-axis boundary modes can be supplied in the requested axis order.
- Added SciPy reference comparisons and checked native regression tests for
  mixed modes, 50,000-wide windows, rank-15 arrays, and invalid inputs.

## 2026-10-23: ndimage binary stencils, numeric labels, selected axes and 1D extrema

- Added `generate_binary_structure_nd` and `iterate_structure_nd` for rank-1
  through rank-15 Boolean footprints, including asymmetric and even-sized
  structuring elements and optional iterated origins.
- Extended `label_nd` to default-kind integer and real(dp) sources (any nonzero
  element foreground), preserving the logical overload and label ordering.
- Added optional `axes(:)` to generic gradient-magnitude and Laplacian
  callbacks, including subset-specific boundary modes; retained old calls.
- Added `uniform_filter1d_nd` (real/complex), `minimum_filter1d_nd`, and
  `maximum_filter1d_nd` with axis-specific origin and boundary handling.
- Added two SciPy reference suites (42,600 compared scalar elements) and two
  Fortran regression programs. No additional runtime dependencies.

## 2026-10-23: ndimage brute-force distances and generic derivatives

- Added `distance_transform_bf_nd`, with Euclidean, taxicab and chessboard
  distances, sampled Euclidean metrics, one-based nearest-background sites,
  all-foreground sentinels and a resource guard for quadratic work.
- Added `generic_gradient_magnitude_nd` and `generic_laplace_nd` with an
  explicit Fortran callback contract, error propagation and mixed modes.
- Added native rank-15/error-path tests and 1,350 SciPy reference values,
  integrated into the focused `ndimage` test runner. No dependencies added.

## 2026-10-23: ndimage Fourier ellipsoid and line-oriented filters

- Added `fourier_ellipsoid_nd` with full-FFT and packed-real-FFT support,
  including a mathematically generalized rank-4-through-rank-15 response.
- Extended median, rank, and percentile filters with `modes(:)` for independent
  per-axis boundary conditions on nonseparable footprints.
- Added `generic_filter1d_nd` with one callback per padded real-valued axis line.
- Added 9,134 SciPy/independent-reference output comparisons and a native
  regression program, without new runtime dependencies.

## 2026-10-22: ndimage callbacks, one-axis filtering, weighted chamfer paths

- Added `correlate1d_nd` and `convolve1d_nd` for real/complex data, all eight
  boundary modes, asymmetric/even-width weights and shifted kernel origins.
- Added `generic_filter_nd` with a real-valued Fortran callback, custom Boolean
  footprints and C-order neighborhood arguments. Added labeled-region
  `labeled_comprehension_nd` with caller-defined reductions and missing-label
  defaults; pixel arguments use deterministic C order.
- Added optional `axis_weights(:)` to the separable taxicab distance transform
  and a new `distance_transform_weighted_cdt_nd` with directed, positive
  3^rank integer step costs and a multi-source binary-heap Dijkstra solver.
- Added 34,296 SciPy/independent-reference elements, including 192 randomized
  custom-stencil maps, and an eighth native ndimage regression program. No
  external dependencies or copied SciPy code.

## 2026-10-22: ndimage extrema footprints, chamfer distances, labeled histograms

- Extended `minimum_filter_nd` and `maximum_filter_nd` with arbitrary Boolean
  `footprint(..)` and per-axis `origin(:)`, including the footprint-only form.
  Preserved the original separable path for rectangular size-only filters.
  Implemented filter-specific (not grayscale dilation) asymmetric orientation.
- Added independent `src/ndimage_chamfer.f90` with taxicab/city-block and
  chessboard `distance_transform_cdt_nd`, optional one-based closest-background
  coordinates, rank-15 input support, and deterministic tie choices.
- Added `count_labels_nd`, `histogram_labels_nd`, and
  `standard_deviation_labels_nd` for selected integer labels in ranks 1-15.
- Added 21,900 passing new SciPy comparisons (4,620 footprint/statistical
  scalar values and 17,280 randomized chamfer distances), plus the seventh
  focused Fortran regression program. All 15 reference scripts pass when run
  separately; all seven Fortran regression programs pass at checked `-O0`
  and optimized `-O3`. No new FPM dependencies.

## 2026-10-22: ndimage labeled extrema, masked morphology, footprint rank filters

- Added labeled extrema and one-based positions (`extrema_labels_nd`,
  `minimum_labels_nd`, `maximum_labels_nd`, `minimum_position_nd`, and
  `maximum_position_nd`), with C-order tie resolution and zeroed missing labels.
- Added `bounding_boxes_nd` and `find_objects_nd`, returning inclusive
  one-based lower/upper bounds for regions of arbitrary rank 1-15.
- Added per-pixel `mask(..)` to binary dilation, erosion, opening, and closing;
  `iterations<=0` now runs each stage until stable, using a finite-step guard
  against nonconvergent custom stencils.
- Added custom Boolean `footprint(..)` and per-axis `origin(:)` to median,
  rank, and percentile filters. The vector-width overloads allow `sizes` to be
  omitted when a footprint is provided.
- Added 8,526 passing 2D/3D SciPy 1.17 comparison values, plus focused native
  Fortran regression tests for tie resolution, error paths, and rank-15 cases.
  The new routines add no external dependencies.

## 2026-10-22: ndimage shifted stencils, region statistics, and watershed IFT

- Added per-axis `origin(:)` to binary dilation/erosion/opening/closing and
  grayscale morphological operations; binary hit-or-miss supports independent
  foreground/background origins. Grayscale morphology also supports non-flat
  real-valued `structure(..)` arrays, optionally masked by `footprint(..)`.
- Added new independent Fortran modules `ndimage_measurements.f90` for labeled
  sums, means, population variances, and one-based centers of mass, and
  `ndimage_watershed.f90` for integer marker-based image foresting transform.
  All new routines accept ranks 1-15. Watershed defaults to correct spatial
  bounds and offers `legacy_linear_bounds` for SciPy's older linear-index rule.
- Added 353 passing SciPy scalar/array comparisons for the new stencil and
  region-statistic operations; 36 bounded and 36 legacy watershed reference
  cases, of which 34 legacy outputs agree with SciPy 1.17 exactly. The two
  remaining cases differ in label tie resolution.
- Added reference scripts and a native rank-15/error-case regression test;
  updated `docs/ndimage.md` and `CODE_PROVENANCE.md`. No new FPM dependencies.

## 2026-10-22: ndimage filtering, morphology, and spline boundary prepadding

- Completed high-order `nearest` and `grid-constant` spline prefiltering in
  `map_coordinates_nd`, affine/shift, zoom, and rotate with SciPy-style 12-sample
  halos and zero-cost singleton-axis reduction. Added a padded allocation
  safety limit (status 2). Expanded independent SciPy spline comparisons to
  576 boundary/affine/shift, 1,156 zoom/rotate, and 16 singleton cases.
- Added `scifort_ndimage` Gaussian smoothing/derivatives (orders 0-3),
  Gaussian gradient magnitude, Gaussian Laplacian, uniform and minimum/maximum
  filters, Sobel and discrete Laplacian, with one-based axes and ranks 1-15.
  Validated 128 cases against SciPy 1.17 across eight boundary modes.
- Added binary dilation, erosion, and connected-component labeling with
  configurable connectivity. Validated 45 cases against SciPy 1.17 and
  added native rank-15 and error-path regressions.
- See `docs/ndimage.md` and `docs/signal_spline_transforms.md` for examples,
  restrictions and numerical verification. These routines add no dependencies.

## 2026-10-21: Rank-15 zoom and rotation

- Added real/complex `zoom_nd`, `rotate_nd`, `zoom_shape_nd` and
  `rotate_shape_nd`, reusing the existing assumed-rank affine spline core.
  Zoom accepts scalar and per-axis factors and `grid_mode`; rotation supports
  selected axis planes, degrees, and shape expansion. Both preserve the
  existing SciFort mode/prefilter defaults and orders 0-5.
- Added `test_signal_zoom_rotate_nd.f90` and a 1,076-case deterministic
  SciPy 1.17 reference comparison. 1,036 cases match within 2.47e-15;
  40 short-array reflect-prefilter cases match to 1.76e-4, reflecting the
  previously documented collocation-versus-recursive initialization gap.
- Full checked standalone signal suite passed 108/108; targeted new routines
  passed `-O3` with bounds checks and floating-point traps. Complete
  optimized-suite execution timed out. FPM was unavailable.
- See `docs/signal_spline_transforms.md` for API examples and conventions.

## 2026-10-21: Rank-independent spline boundary modes and affine transforms

- Added all eight SciPy `ndimage.map_coordinates` boundary modes to
  `map_coordinates_nd`, and new public assumed-rank real/complex
  `affine_transform_nd` and `shift_nd` for ranks 1-15 and spline orders 0-5.
- Shared the existing tensor-product cardinal B-spline kernel and
  column-major packing logic. Affine sampling computes coordinates in small
  batches; periodic filtering adds a circulant spline solver.
- Added permanent rank-15/noncontiguous/public-API tests and 528 SciPy 1.17
  interpolation/affine/shift comparisons. 480 cases agree to 2.89e-15;
  48 short-array reflect-prefilter cases differ by up to 1.40e-3 because
  the exact collocation solve differs from SciPy's recursive initialization.
- Full checked standalone signal suite: 107/107; targeted `-O3` test with
  bounds and IEEE traps: passed. Full optimized suite timed out; FPM is absent.
- For orders 2-5, `prefilter=.true.` with `nearest` or `grid-constant`
  explicitly returns invalid input until SciPy-style edge prepadding is added.
  See `docs/signal_spline_transforms.md`.

## 2026-10-21: Assumed-rank N-dimensional spline interpolation

- Added `scifort_signal_spline_nd`: rank-1..15 real/complex coefficient filters
  (`spline_filter_nd`, `spline_filter1d`, `cspline_nd`, `qspline_nd`) and
  coordinate evaluator (`map_coordinates_nd`) for B-spline orders 0..5.
- Uses banded cardinal-spline collocation on each axis and mirror-reflected
  tensor-product interpolation; supports constant fill and optional prefilter.
- Added a public regression, randomized SciPy reference comparison (288 cases,
  5,760 coordinate samples), documentation, and standalone build ordering.
- Numerical validation: checked standalone signal suite passed 106/106 tests;
  see new docs for numerical tolerances and high-dimensional cost guard.

- Added two-point/three-point/complex-step model Jacobians and `diff_step` to TRF/dogbox `curve_fit`.
# Changelog

## 2026-10-23: ndimage Fourier ellipsoid and line-oriented filters

- Added `fourier_ellipsoid_nd` with full-FFT and packed-real-FFT support,
  including a mathematically generalized rank-4-through-rank-15 response.
- Extended median, rank, and percentile filters with `modes(:)` for independent
  per-axis boundary conditions on nonseparable footprints.
- Added `generic_filter1d_nd` with one callback per padded real-valued axis line.
- Added 9,134 SciPy/independent-reference component comparisons and a native
  regression program, without new runtime dependencies.


- Added `minimize(method='trust-ncg')` and direct `trust_ncg` with SciPy-style Steihaug CG trust-region steps, dense-Hessian/Hessian-product support, and generated SciPy 1.17 regressions.
- Added dogbox `tr_solver='lsmr'` with SciPy-style scaled iterative Gauss-Newton steps; `jac_sparsity` now automatically selects LSMR for dogbox as well as TRF.
- Added `jac_sparsity` grouped finite-difference Jacobians for TRF/dogbox least squares, with automatic LSMR selection for both methods.
- Added matrix-free TRF Jacobians with paired `jac_matvec=`/`jac_rmatvec=` callbacks and automatic LSMR dispatch.
- Added TRF `tr_solver='lsmr'` with a matrix-free Fong-Saunders LSMR core and SciPy-style regularized two-dimensional trust-region subspace solving.
- Added SciPy-style complex-step (`diff_scheme='cs'`) Jacobians for `least_squares`.
- Added TRF/dogbox iteration callbacks with SciPy-style early-stop status `-2`.
- Added SciPy-style `x_scale='jac'` Jacobian-column scaling for `least_squares`.
- Added SciPy-style three-point finite-difference Jacobians for `least_squares`.
- Added callable robust-loss support (`loss_fun=`) for TRF/dogbox `least_squares`.

All notable changes will be documented here.

## Unreleased

- Added `scifort_signal.cspline2d_eval`, `qspline2d_eval`, `cspline2d_grid`,
  and `qspline2d_grid`, with real and complex tensor-product spline evaluation,
  one-based scalar and paired-array coordinates, rectangular grids, and
  mirror-symmetric extrapolation. These are independent SciFort extensions
  compatible with evaluating quadratic/cubic spline coefficient arrays using
  SciPy `ndimage.map_coordinates(..., mode="mirror", prefilter=False)`.
  Validated 256 randomized SciPy/NumPy cases (17,920 coordinate evaluations;
  largest absolute difference 9.16e-16), including one-pixel image axes,
  and added a new full public-API signal regression. See
  `docs/signal_spline2d.md`; no new external Fortran dependencies.


- Added rank-two `scifort_signal.cspline2d`, `qspline2d`, and `spline_filter`
  for real and complex images, symmetric IIR boundary handling, optional
  cubic smoothing and mirror-boundary separable FIR reconstruction.
  Implemented and attributed the symmetric filtering recurrences from SciPy
  1.17's BSD-3-Clause implementation. Verified 44 SciPy/NumPy reference
  cases (maximum absolute error 2.5e-7) and 104/104 checked standalone
  signal tests. A targeted `-O3` regression passed with IEEE traps; full
  optimized-suite and project-wide FPM tests were not run for this change.
  See `docs/signal_spline2d.md`.


- Added `scifort_signal.cspline1d`, `qspline1d`, `cspline1d_eval`, and
  `qspline1d_eval` for real and complex samples, with scalar and vector
  evaluation, mirror-symmetric boundaries and non-unit knot spacing.
  Added cubic smoothing through SciPy's recursive formulation and explicit
  invalid-input statuses for unsupported quadratic smoothing or unstable
  cubic smoothing poles. The algorithms are adapted with attribution from
  SciPy 1.17 (BSD-3-Clause); no new dependencies. Verified 480 randomized
  SciPy coefficient/evaluation cases (maximum absolute errors 1.351e-15
  and 2.798e-14) and 103/103 checked standalone signal tests. The new
  targeted spline regression passed at -O3 with IEEE traps; the full
  optimized suite timed out. See `docs/signal_spline1d.md`.


- Added `scifort_signal.gauss_spline`, `kaiser_beta`, `kaiser_atten`, and
  `kaiserord`, plus optional parameterized Kaiser-window support in `firwin`.
  The Gaussian spline and elementary Kaiser helpers use elemental Fortran
  interfaces; `kaiserord` returns a typed `(numtaps, beta, status)` result.
  Independently compared 2,388 seeded cases (including Kaiser FIR designs)
  against SciPy 1.17 with maximum absolute difference 7.77e-16. The full
  signal suite passed 102/102 checked and 102/102 optimized tests with GNU
  Fortran 14.2.0; the new routines also passed targeted IEEE-trap checks.
  No new Fortran dependencies. See docs/signal_kaiser_spline.md.


- Added `scifort_signal.gammatone` with FIR (orders 1-24) and eighth-order
  IIR auditory-filter designs, and `scifort_signal.iircomb` with notch/peak
  modes, selectable DC pass behavior, and sampling-frequency validation.
  Both return the existing `butter_result` coefficient/status type. Tested
  1,350 deterministic-seed randomized cases against SciPy 1.17.0; maximum
  absolute coefficient error 2.13e-14. GNU Fortran 14.2.0 passed
  101/101 checked signal tests and 101/101 optimized signal tests; the
  focused filter regression also passed at -O3 with IEEE traps. Added permanent public-interface
  regressions and a reproducible Python reference comparison. The formulas
  follow Slaney/Orfanidis as adapted from SciPy 1.17 (BSD-3-Clause); licensing
  and provenance are recorded below. No new Fortran dependencies.


- Added coefficient-vector `general_cosine(n, coefficients, sym, status)`
  and `get_window('general_cosine', n, coefficients=...)` with SciPy's
  origin-centered cosine coefficients and symmetric/periodic modes.
  Integrated the preexisting `dpss` solver into `get_window('dpss', n,
  parameter=NW)`, using a single approximately peak-normalized taper.
  Both are public through `scifort_signal`. Updated the standalone signal
  build order and reference test build lists for the DPSS dependency.
  Checked 1,092 deterministic and seeded cases against SciPy 1.17 (maximum
  absolute error 2.64e-12), and passed 100/100 checked and optimized signal
  tests. The new cosine code is an independent MIT-licensed implementation;
  no new Fortran dependencies or vendored third-party code were introduced.


- Added SciPy-style Dolph-Chebyshev (`chebwin`), Taylor (`taylor`), and
  Kaiser-Bessel-derived (`kaiser_bessel_derived`) windows, both as direct
  functions and through `get_window` (23 supported window families total).
  Taylor exposes `nbar`, `sll`, and `norm`, and KBD enforces even symmetric
  windows and the Princen-Bradley condition. Direct functions default to
  symmetric windows; the generic dispatcher defaults to periodic windows.
  Adapted formulas and synthesis conventions from SciPy 1.17 (BSD-3-Clause)
  with attribution in `CODE_PROVENANCE.md` and `THIRD_PARTY_LICENSES.md`.
  The Fortran implementation remains dependency-free and uses a log-scaled
  DFT for Chebyshev windows. Added 1,718 SciPy 1.17 reference comparisons,
  maximum absolute error 4.82e-13, and new targeted/public regressions.


- Expanded `scifort_signal.get_window` from seven to twenty spectral-window
  families, including Parzen, Bohman, Blackman-Harris, Nuttall, Bartlett-Hann,
  triangular, cosine, Lanczos, exponential, Kaiser, Gaussian, generalized
  Gaussian, and generalized Hamming windows. Added optional `parameter` and
  `parameter2` arguments without breaking the earlier call signature.
  Corrected the no-argument Tukey alpha from 0.25 to SciPy's 0.5 default;
  explicitly pass `parameter=0.25_dp` to obtain the old choice.
  The new stand-alone Fortran regressions and 1,092 SciPy 1.17 window
  comparisons pass (maximum absolute error 1.43e-14). No additional
  Fortran dependencies. See `docs/signal_window_families.md`.

- Expanded the existing `scifort_signal.detrend` generic to real and complex
  input arrays of ranks 1–3, including selected one-based axes and segmented
  linear fitting. Existing real one-dimensional calls are unchanged.
  Added independent Fortran tests and 720 SciPy 1.17 reference comparisons
  (maximum absolute error 5.35e-15).
- Added the SciFort spectral extension `multitaper_harmonic_scan`, with
  Bonferroni, Benjamini–Hochberg and Benjamini–Yekutieli corrected p-values,
  configurable significance level, bin selection, default DC/Nyquist exclusion,
  and significant-line flags. Compared 160 randomized cases against independent
  SciPy/NumPy references (maximum corrected-p-value error 5.03e-13).
  The corrections require valid pointwise null-model p-values, which remain
  approximate for harmonic F statistics. No additional Fortran dependencies.


- Added `multitaper_stream`, an incremental real/complex multichannel
  multitaper cross-spectral estimator with cached DPSS tapers, arbitrary
  sample chunks, one-based sample-axis selection, fixed-size window storage,
  and history-free running mean cross-spectral matrices. Existing batch
  APIs remain unchanged; no new Fortran dependencies. New native tests
  compare all weighting modes, reset, empty chunks, wrong-type inputs, odd
  and even FFT lengths and chunk boundaries with batch results. An independent
  randomized SciPy/NumPy cross-spectrum comparison covers 80 cases
  (maximum absolute spectral error 1.45e-11). GNU Fortran 14.2.0 checked
  and optimized standalone signal suites both passed 91/91 tests;
  targeted optimized regression also passed with IEEE floating-point traps.
  Project-wide FPM testing remains unavailable in this environment.


- Added time-frequency `multitaper_cross_spectrogram` for all real/complex
  input combinations, `multitaper_tf_coherence` with bounded magnitude-squared
  coherence and zero-power handling, and rank-2 multichannel
  `multitaper_spectral_matrix_spectrogram` with optional sample axis. Spectral
  matrices are Hermitian and use the same one-sided/two-sided FFT/normalization
  conventions as existing multitaper PSDs, with independent adaptive weights
  for each channel. Two Fortran regressions and 120 randomized SciPy DPSS/NumPy
  FFT comparisons passed (maximum spectral error 1.30e-11, coherence error
  2.90e-11). The full GNU Fortran 14.2.0 checked and optimized standalone
  signal suites passed 89/89 tests; the new routine passed optimized tests
  with IEEE traps. These are SciFort extensions; no new dependencies were
  added. FPM project-wide testing was unavailable.

- Added `multitaper_spectrogram` (real/complex, single-channel and
  rank-2 multichannel) with one-based sample axes, cached periodic DPSS
  tapers, noverlap, FFT padding, per-segment mean subtraction and
  equal/eigen/adaptive weighting. The output records one-sided or two-sided
  frequency bins, segment times and PSDs. Added `multitaper_harmonic_ftest`
  implementing a DPSS regression F-test with complex line amplitudes,
  F-statistics and approximate pointwise probabilities for spectral lines.
  These are SciFort extensions; no dependencies were added. Public and
  internal Fortran regressions and a reproducible NumPy/SciPy comparison
  script cover 90 spectrogram and 100 harmonic test cases (maximum absolute
  errors 8.83e-12 for PSD and 2.48e-13 for probabilities). The full
  standalone signal suite passes **87/87** with runtime checks and also
  **87/87** at `-O3`; the new regression passes with IEEE traps at `-O3`.
  Project-wide `fpm test` was not run (FPM unavailable).

- Added Thomson-style `weighting='adaptive'` to real/complex one-dimensional
  and multichannel multitaper PSD, spectral matrices, cross-spectra and
  coherence. Added `multitaper_psd_ci` returning frequency-dependent
  effective degrees of freedom and 95%-default pointwise Student-t
  log-jackknife intervals. Each leave-one-taper-out estimate refits the
  adaptive weights. The API preserves `equal`/`eigen` behavior, needs no
  new Fortran dependencies, and is a SciFort extension. New Fortran
  regressions and 140 independent SciPy DPSS/NumPy FFT/Student-t randomized
  comparisons cover real/complex signals and multichannel matrices.
  The complete runtime-checked signal suite passes 85/85 tests, including
  the new tests. The three targeted optimized multitaper regressions pass
  at `-O3` with runtime checks and floating-point traps; the full optimized
  suite and FPM project-wide tests were not run. The standalone test runner
  now compiles the DPSS and multitaper source modules explicitly.

- Expanded `scifort_signal.multitaper_psd` to complex 1-D signals and real or
  complex 2-D arrays with one-based sample-axis selection. Added
  `multitaper_spectral_matrix`, `multitaper_csd` (including mixed real/complex
  pairs), and `multitaper_coherence` with FFT-consistent frequency axes and
  Hermitian interchannel cross-spectra. The original real 1-D API is
  preserved. Added public and internal Fortran regressions plus 160
  independent SciPy DPSS/NumPy FFT comparisons (maximum absolute spectral
  matrix error 4.54e-13). The complete standalone checked signal suite passes
  84/84 tests; the new internal regressions also pass with `-O3` and IEEE
  traps. The full optimized signal-suite rebuild did not complete within the
  execution limit; project-wide `fpm test` was unavailable. No new FPM
  dependencies.

- Added `scifort_signal.dpss`, a SciPy-style implementation of the discrete
  prolate spheroidal sequences (Slepian tapers), with symmetric/periodic
  modes, multiple tapers, concentration ratios and `l2`/`approximate`
  normalization. Uses an independently written implicit-shift QL symmetric
  tridiagonal eigensolver; no new FPM dependency. The Fortran result uses
  `(sample,taper)` indexing; SciPy uses `(taper,sample)`.
- Added `multitaper_psd`, a SciFort-specific real 1-D multitaper PSD
  estimator with optional sampling rate, FFT length, mean removal and
  equal/concentration weighting. Added two Fortran test executables and
  two reproducible SciPy/NumPy comparison scripts. All 344 DPSS and 240
  multitaper comparisons passed against SciPy 1.17/NumPy; max errors were
  1.46e-11 (DPSS sample), 3.77e-15 (concentration) and 6.89e-13 (PSD).
  All 82/82 standalone signal tests passed with runtime checks. The new
  regression passes at `-O3` with bounds checks and floating-point traps;
  the full optimized signal suite was not completed. Full FPM testing
  remains unverified (FPM unavailable).

- Extended `scifort_signal.place_poles` with selectable `YT` and `KNV0`
  determinant-improving iterations, SciPy-style `rtol` and `maxiter` controls,
  and `nb_iter`, `rtol`, and `det_transfer` output diagnostics. KNV0 uses
  feasible real-eigenspace projections; YT updates real pole pairs using
  singular vectors and complex pairs using a Hermitian eigenproblem. The
  update ordering follows Tits--Yang, but initialization and tie handling
  are SciFort-specific and not claimed to produce SciPy-identical gain matrices.
  Added the iterative pole-placement regression, a dedicated reference driver
  and 666 SciPy-validated method-specific comparisons. Existing 470 MIMO
  comparisons and 280 CWT/SISO/fully actuated pole cases still pass.
  The standalone GNU Fortran 14.2.0 checked signal suite passes 80/80;
  targeted pole-placement `-O3` regressions pass with floating-point traps.
  Full optimized signal suite timed out; `fpm test` was not available.


- Expanded `scifort_signal.place_poles` to partially actuated MIMO
  (`1<m<n`) and controllable rank-deficient input matrices even when `m>=n`.
  Added pivoted nullspace eigenvector assignment, rank-one independence
  refinement for real poles, conjugate-pair candidate selection, real feedback
  reconstruction, and validation of pole multiplicity against `rank(B)`.
  This is not SciPy's full Yang-Tits robustness optimizer and cannot
  guarantee all defective or poorly conditioned assignments. The public
  procedure signature and SISO/fully actuated paths remain unchanged.
  Added a permanent MIMO regression, strengthened the old pole-placement
  test, and verified 470/470 SciPy-validated MIMO systems (maximum observed
  closed-loop eigenvalue error 3.20e-9). Signal-only GNU Fortran 14.2.0
  test suites pass 79/79 in both checked and optimized builds. Full
  project `fpm test` remains unverified (FPM unavailable).

- Added the Ricker-wavelet `find_peaks_cwt` peak detector, `cwt_ricker`
  transform, and `ricker_wavelet` for real-valued 1-D signals. The ridge-linking,
  width constraints, noise percentile and SNR filtering follow SciPy 1.17's
  wavelet peak detection behavior, with one-based output indices. Added analytic
  `place_poles` for controllable single-input real systems (including conjugate
  complex closed-loop roots) and fully actuated multi-input systems with
  full-row-rank input matrices. General partially actuated MIMO YT/KNV0 is
  explicitly not implemented (`pole_status_unsupported`). These algorithms
  add no FPM dependencies. The signal-only GNU Fortran 14.2.0 regression suite
  passes 78/78 with runtime checks and 78/78 at `-O3`; both new tests also
  pass with floating-point traps. An independent reproducible comparison script
  checks 300/300 nondegenerate peak-detection cases exactly and 280/280
  controllable pole-placement cases against SciPy 1.17. The worst closed-loop
  eigenvalue error is 1.11e-7 for an ill-conditioned sixth-order system.
  Flat-signal, near-tied CWT extrema can differ by one sample due to
  floating-point accumulation; 299/300 match in a separate suite that includes
  constant signals. Full project FPM tests remain unverified (FPM unavailable).


- Added `scipy.signal.argrelextrema`, `argrelmax`, and `argrelmin` for one-dimensional
  real signals, with pure user-supplied scalar comparators, `clip`/`wrap` edge
  modes, arbitrary positive `order`, and one-based returned indices. Added the
  `peak_detect` excursion-based peak/valley routine as an **SciFort extension**
  inspired by Eli Billauer and Jason Christopherson's MIT-licensed `peaks`
  library. No third-party implementation was copied or added as a dependency.
  New regression and public exports pass, with 1,350 randomized SciPy comparisons
  and 450 independent Billauer-style state-machine comparisons. The full
  signal suite passes 76/76 tests at both `-O1` and `-O3` with bounds/runtime
  checks.


- Added independent `scifort_signal_deconvolve` with four real/complex
  input/divisor combinations, quotient/remainder outputs, and explicit invalid
  input status; added `scifort_signal_envelope` for rank-1/2/3 real/complex
  input, signed FFT-bin bands, axis selection, resampling, squared magnitude,
  and lowpass/all/none residuals. Both are re-exported through `scifort_signal`
  without changing existing interfaces. New direct/public native regressions,
  separate reference drivers and 1,900 deterministic comparisons against
  SciPy 1.17 (600 deconvolution, 700 rank-1, 380 rank-2 and 220 rank-3
  envelope cases); maximum observed absolute difference approximately
  1.78e-13. GNU Fortran 14.2 standalone signal suite passes 74/74 checked and 74/74
  optimized (`-O3`). The new regression also passes at `-O3` with bounds
  checks and floating-point traps. Project-wide FPM tests remain unverified
  because FPM was unavailable. No external Fortran runtime dependencies or
  copied SciPy source.

- Added independently implemented `scifort_signal_spatial` with public
  `order_filter` (real/integer rank 1-3, integer/real/logical odd footprint,
  zero-based rank, zero-padded edges) and `sepfir2d` (real/complex/integer
  image, odd row/column FIR taps, mixed real/complex support, reflected edges).
  Added direct/public Fortran regressions and a seeded SciPy 1.17 comparison
  driver; 888 reference cases pass, of which 501/504 masked-order cases agree
  with `scipy.signal.order_filter` and three sparse-median cases use the
  documented behavior verified by `scipy.ndimage.generic_filter` instead.
  GNU Fortran 14.2.0 standalone signal suites pass 72/72 tests with runtime
  checks and 72/72 at `-O3`; the new tests also pass with FP traps enabled.
  FPM was unavailable. No new Fortran dependencies.

- Added `scifort_signal_denoise` and public `medfilt`, `medfilt2d`, `wiener`
  APIs for rank-1/2/3 real(dp) and default-integer arrays, zero padding,
  scalar/per-axis kernels, Wiener noise estimation, even-window alignment,
  and explicit input validation. New Fortran regressions test direct and
  public APIs; 1,050 seeded comparisons against SciPy 1.17 passed with
  maximum absolute error 6.67e-15. The standalone signal suite passes
  70/70 checked and 70/70 optimized GNU Fortran tests; the new direct
  regression and all SciPy comparisons also pass with floating-point traps
  enabled. FPM was not available. No runtime dependencies added.


- Added independently implemented `istft` inverse short-time Fourier
  transforms with real one-sided and complex two-sided outputs, odd and
  padded FFT lengths, configurable overlap, named/custom windows, boundary
  trimming, `spectrum`/`psd` scaling, and a noninvertible-overlap status.
  Added `check_cola` and `check_nola` for named or explicit window samples.
  Native regression tests include the public `scifort_signal` re-export.
  Independent SciPy 1.17 comparisons: 480 inverse transforms (max abs error
  6.84e-13) and 1,120 overlap-condition results from 560 windows.
  The standalone GNU Fortran 14.2.0 signal suite passes 68/68 tests
  with runtime checks and 68/68 at `-O3`. The two new tests pass with
  floating-point traps enabled, too. FPM testing was unavailable. No new
  Fortran dependencies or third-party source code were added.


- Optimized rank-independent `oaconvolve_nd` by transforming the unchanged
  kernel only once and reusing its spectrum for every overlap-add tile. For
  strongly asymmetric inputs it now tiles multiple selected axes of the same
  operand, uses one fixed FFT workspace, and crops contributions during
  accumulation; unlike the previous implementation it does not allocate or
  perform separate full FFT convolutions for each block. For opposing long
  directions it splits one operand per call, retaining the other as a fixed
  kernel. The original public API, rank-1 through rank-15 support, output
  shape rules and 1D overloads are unchanged. Added a deterministic
  multi-axis regression, a 232-case independent SciPy 1.17 reference driver,
  and a GNU Fortran `-O3` before/after benchmark. Timings for two 5x5-kernel
  convolutions improved from 0.0789 to 0.0246 seconds at shape 513x257 and
  from 0.3021 to 0.1004 seconds at shape 1025x513. No new runtime
  dependencies. The complete standalone signal suite passes 66/66 tests
  with runtime checks and 66/66 tests at `-O3`; the previous 400-case
  rank-1/2/3/5/15 SciPy comparison passes unchanged (maximum error 4.45e-14),
  and the new 232-case multi-axis comparison passes (maximum error 1.40e-13).
  FPM was not available in this test environment.

- Added `oaconvolve_nd`, a caller-output Fortran 2018 multidimensional
  overlap-add convolution for ranks 1-15, including mixed real/complex
  operands, selected one-based axes, singleton broadcasting, and all three
  SciPy convolution modes. For highly unequal selected-axis lengths, the
  longer input is divided into FFT blocks and the block results are
  accumulated directly into the requested output crop. For balanced lengths,
  a single multidimensional FFT avoids redundant transforms. The existing
  1D `oaconvolve` functions remain unchanged. Added two standalone Fortran
  regressions and 400 independently seeded SciPy 1.17 comparisons across
  ranks 1/2/3/5/15, with maximum absolute difference 4.55e-14.
  The standalone checked and optimized signal suites both pass 65/65 tests
  with GNU Fortran 14.2.0; the optimized new regression also passes with
  bounds checking and floating-point traps. No third-party Fortran
  dependencies or SciPy-derived source code were added. FPM is unavailable
  in the test environment.

- Added caller-output Fortran 2018 `fftconvolve_nd` for rank-1 through
  rank-15 arrays, mixed real/complex data, selected one-based convolution
  axes, singleton broadcasting on remaining axes, and SciPy-compatible
  `full`/`same`/`valid` slicing. It uses strided separable FFTs in packed
  column-major storage, returns status plus NaN-filled output on invalid
  input, and supports noncontiguous actual arguments. Existing rank-1/2/3
  `fftconvolve` functions are unchanged. Added independent native tests and
  576 randomized rank-15 SciPy 1.17 comparisons (max absolute difference
  5.02e-15). GNU Fortran 14.2.0 standalone signal tests pass 63/63 with
  runtime checks and 63/63 under -O3. No added dependencies or imported
  SciPy implementation code; project-wide FPM tests were not run.

- Added Fortran 2018 `filtfilt_nd` and `sosfiltfilt_nd` with
  caller-allocated, shape-matched rank-1 through rank-15 outputs, selected
  one-based axes, SciPy-style odd/even/constant/none endpoint padding,
  eight BA and four SOS real/complex type combinations, and scalar/shape/
  coefficient/steady-state validation. Noncontiguous array sections are
  supported. The existing rank-1/2/3 array-returning functions remain
  compatible. Added deterministic Fortran tests and 684 numerical comparisons
  against SciPy 1.17 (maximum absolute error 3.33e-16). The complete checked
  standalone signal suite passed 61/61 tests. No extra dependencies.

- Added Fortran 2018 assumed-rank `sosfilt_nd` for real and complex
  second-order-section filtering on rank-1 through rank-15 arrays. Supports
  both mixed real/complex coefficient-and-sample combinations, selected
  one-based axes, independent per-channel states, and streaming across blocks.
  The state layout retains SciPy's leading SOS-section dimension up to rank
  14; rank-one packed state buffers support rank 15 without requiring a
  forbidden rank-16 Fortran array. Noncontiguous array sections use compiler
  copy-in/out. No full-array type promotion is required. Added one standalone
  regression and 320 independently seeded SciPy 1.17 reference comparisons,
  including outputs and final states. No new dependencies were added.


- Added caller-output `lfilter_nd` for real, complex and mixed real/complex
  numerator, denominator and input types. This Fortran 2018 assumed-rank
  subroutine filters ranks 1-15 along any one-based axis and accepts
  caller-allocated initial/final state arrays whose filtered-axis extent is
  `max(size(a),size(b))-1`. Coefficients are normalized once per call and
  each data channel uses a reusable direct-form-II delay line. Existing
  `lfilter` rank-1/2/3 functions remain compatible. A new Fortran regression
  and 300 rank-5/rank-15 randomized SciPy 1.17 comparisons verify outputs,
  final states and all eight real/complex type combinations. An independent
  240-case comparison passes in checked and optimized builds. The checked
  standalone signal suite passes 58/58; the targeted optimized test passes
  with floating-point traps. The complete optimized runner timed out, and
  project-wide FPM tests were not run. No external dependencies were added.


- Refactored the eight rank-1/2/3/4 Savitzky-Golay functions to use one
  common rank-independent, allocation-free per-channel filtering kernel.
  Added `savgol_filter_nd` for real/complex assumed-rank input/output arrays
  through the standard Fortran 2018 rank limit (15). It requires a
  caller-supplied shape-matched output, accepts one-based axis selection,
  and preserves existing boundary modes, derivatives and statuses.
  Added rank-5/rank-15 tests, noncontiguous sections, and 300 randomized
  SciPy comparisons. The standalone signal suite passes 57/57 in both
  runtime-checked and optimized builds; 600 existing rank-1/2/3 and 400
  rank-4 SciPy comparisons also pass. No new external dependencies.


- Refactored Savitzky-Golay filtering to precompute the center and all
  interpolated-edge polynomial weights once per call rather than per channel.
  Added rank-4 real/complex `savgol_filter` with one-based axis selection and
  unchanged five boundary modes, derivative orders, and error statuses.
  Added 400 seeded SciPy 1.17 rank-4 comparisons and a reproducible CPU
  benchmark (about 4.8x faster for a 768-channel example with GNU Fortran
  14.2 `-O3`). The runtime-checked standalone signal suite passes 56/56;
  targeted optimized rank-2/3 and rank-4 tests pass with floating-point
  traps. The complete optimized suite timed out; project-wide FPM tests
  were not run. No new dependencies.

- Extended `savgol_filter` to rank-2 and rank-3 real/complex arrays with
  one-based axis selection; added complex rank-1 filtering as a SciFort
  extension. Preserved five edge modes, derivatives, complex constant
  padding, and 1D backward compatibility. Added 600 seeded SciPy 1.17
  comparisons, including separate real/imaginary complex references.
- Extended `lfilter` with 18 mixed-type overloads for rank-1/2/3 inputs,
  automatically promoting real/complex coefficient and signal combinations
  to complex output and initial/final states. Added 540 seeded SciPy 1.17
  comparisons. The runtime-checked standalone signal suite passes 55/55;
  six affected tests also pass separately at `-O3` with floating-point traps.
  The full optimized suite timed out; project-wide FPM tests were not run.
  No new dependencies.

- Extended `sosfilt`, `filtfilt`, and `sosfiltfilt` to real/complex rank-2
  and rank-3 arrays with one-based `axis` selection (default last axis).
  SOS initial/final states are shaped `(n_sections,...)`, replacing the
  filtered axis by two delay values. Added two fixed regression programs and
  480 seeded SciPy 1.17 comparisons (maximum absolute error `5.55e-16`).
  The checked standalone signal suite passes 53/53 tests. No new dependencies.


- Extended `fftconvolve` from 1D to real, complex and mixed-type rank-2 and
  rank-3 arrays. Optional one-based `axes=[...]` supports independent
  transforms over selected dimensions and singleton broadcasting. Extended
  `lfilter` to real/complex rank-2 and rank-3 inputs, with an optional
  one-based `axis`, initial states `zi`, and final states `zf`. Added three
  fixed regression programs and independent SciPy 1.17 comparisons covering
  1,000 FFT-convolution, 400 rank-2 filtering, and 300 rank-3 filtering
  cases. All 51 standalone signal tests passed with checked compilation.
  No additional dependencies.


- Added SciPy-style rank-two `convolve2d` and `correlate2d` with real,
  complex, and mixed inputs; full/same/valid output modes; fill/wrap/symm
  boundaries; complex fill values; and status-based validation. The
  independent direct implementation preserves SciPy's even-kernel
  `correlate2d` centering and its complex valid-mode operand-swapping
  convention. Added the 48th standalone signal regression test and a
  1,600-case SciPy 1.17 comparison (maximum scaled error 7.85e-16).
  No new FPM dependencies.

- Added opt-in `solver='schur'` for all four MIMO frequency/Bode functions
  and `prepare_mimo_frequency(...,solver='schur')`. The complex unitary
  shifted-QR reduction computes a reusable upper-triangular Schur form,
  enabling triangular multi-input solves. Nonconvergence falls back to the
  existing Hessenberg path; `mimo_frequency_plan%method` reports the actual
  reduction. The default `auto` method remains unchanged. Added regression,
  150-system SciPy comparison, and 48-state benchmark. No new dependencies.

- Added `prepare_mimo_frequency` and `mimo_frequency_plan` for reusing the
  frequency-independent Hessenberg reduction across multiple continuous-
  or discrete-time MIMO frequency sweeps and Bode plots. Extended all four
  MIMO APIs with optional one-based `inputs=`/`outputs=` channel selection,
  supporting reordering and repeated indices, without altering SISO calls.
  Kept automatic analog frequency grids identical between prepared and
  ordinary calls. Added a fixed regression, a 280-system SciPy comparison
  (18,753 selected channel-frequency values), and a repeated-sweep benchmark
  that measured a 3.6x speedup for 25 short sweeps of an 80-state model.
  No additional FPM dependencies.

- Added a unitary Householder Hessenberg reduction and adjacent-pivot
  shifted solver to MIMO `freqresp_mimo`, `dfreqresp_mimo`, `bode_mimo`, and
  `dbode_mimo`. A new optional `solver='auto'|'direct'|'hessenberg'` chooses
  automatic, original dense, or reusable Hessenberg evaluation. The
  default selects Hessenberg for at least eight states and eight frequency
  values, preserving SISO APIs and earlier calls. Checked 45 standalone
  signal tests, 480 seeded SciPy 1.17 MIMO references and a separate
  direct-vs-reduced regression. An 80-state, 600-frequency GNU Fortran
  14.2 `-O3` benchmark measured 11.6x speedup (machine-dependent).
  A complete optimized standalone signal-suite run exceeded the execution
  limit; four targeted `-O3` frequency-response tests passed.


- Added `freqresp_mimo`, `dfreqresp_mimo`, `bode_mimo`, and `dbode_mimo`
  for real/complex MIMO state-space systems. A single complex pivoted
  factorization handles all input columns at each frequency. New result
  types expose arrays in (output,input,frequency) order, with independent
  Bode phase unwrapping and partial singular-frequency reporting. Kept
  all existing SISO interfaces unchanged. Added the 44th standalone signal
  test and 480 seeded SciPy 1.17 matrix-solve comparisons spanning 3,840
  frequencies and 22,648 channel values. No new FPM dependencies.


- Changed continuous/discrete SISO state-space `freqresp`, `bode`,
  `dfreqresp`, and `dbode` to evaluate the shifted state-matrix resolvent
  directly by pivoted linear solves, rather than forming characteristic
  polynomials. Added a pole-scale automatic analog frequency grid that also
  avoids polynomial expansion, nonfinite/singular status handling,
  high-order regression tests, and 560 seeded SciPy reference cases
  covering 5,600 frequency points. TF and ZPK response paths are unchanged.


- Added continuous and discrete LTI frequency-response and Bode functions
  (`freqresp`, `bode`, `dfreqresp`, and `dbode`) through
  `src/signal_lti_frequency.f90`, plus explicit `*_zpk` variants. Real and
  complex SISO transfer-function and state-space inputs are supported.
  Discrete polynomial coefficients are converted from descending powers of
  `z` to delay coefficients before invoking `freqz`, and `dbode` returns
  rad/s after scaling by `dt`. Added 42nd standalone signal test and
  720 seeded numerical comparisons with SciPy 1.17.

- Extended continuous-time `step`, `impulse`, `step_zpk`, and `impulse_zpk`
  with optional automatic SciPy-style time grids (`n=`, default 100), using
  eigenvalue time scales with zero-root handling, while retaining explicit
  irregular grids. Added mixed real/complex ZPK gains and sampled input
  overloads for continuous-time simulation. New Fortran regressions cover
  optional grids, complex poles, integrators, and mixed-type calls; a new
  seeded SciPy 1.17 comparison covers 184 state-space, TF, and ZPK cases.

- Added continuous-time `scipy.signal.lsim`, `step`, `impulse` and the
  corresponding separate ZPK entry points. Exact FOH and ZOH sampled input
  propagation reuses SciFort's Padé matrix exponential. Supports real and
  complex state-space systems, TF/ZPK representations, initial states,
  and independent MIMO step/impulse responses. Added Fortran tests and
  522 seeded SciPy/SciPy-matrix-exponential comparisons. Continuous-time
  step and impulse now also support an automatic time array.

- Added `scipy.signal.dlsim`, `dstep`, `dimpulse`, and separate ZPK entry
  points (`dlsim_zpk`, `dstep_zpk`, `dimpulse_zpk`), supporting real/complex
  discrete state-space, SISO transfer functions/ZPK, multi-input responses,
  initial states, irregular input time grids, and explicit result/status
  types. Added a new Fortran regression and 498 seeded SciPy/independent
  recurrence comparisons. Existing SciFort state-space conversion routines
  are reused without new dependencies.

- Added `scipy.signal.cont2discrete` state-space and transfer-function
  discretization and `cont2discrete_zpk` for zeros/poles/gain. Supports ZOH,
  FOH, generalized bilinear (including Tustin, Euler, backward difference),
  and impulse invariant discretization. Includes a scaled Padé matrix
  exponential, a pivoted linear solver, exact GBT ZPK root mapping,
  result types/status codes, fixed Fortran tests, 390 SciPy 1.17 comparisons,
  and updated standalone verification scripts. No new dependencies.

- Added SciPy-style `scipy.signal.tf2ss`, `ss2tf`, `zpk2ss` and
  `ss2zpk` through `src/signal_state_space.f90` with controller-canonical
  realizations, 1D/2D numerator support, real/complex matrices, optional
  one-based MIMO input/output selection and explicit status/result types.
  Includes checked Fortran regressions and 832 seeded SciPy 1.17 numerical
  comparisons. Existing signal conversion modules and FPM dependencies are
  unchanged.

- Added `scipy.signal.tf2sos` with real and complex coefficient overloads,
  digital nearest/keep_odd/minimal pairing, and analog minimal sections.
  Fixed a `keep_odd` origin-padding failure on unequal transfer-function
  degrees. Added fixed Fortran tests and seeded analog/digital frequency-
  response comparisons (992 valid, 64 rejected improper requests) against
  SciPy 1.17 or the directly evaluated transfer-function polynomial.

- Added digital `scipy.signal.sos2zpk` for real/complex second-order-section
  matrices, including leading-zero coefficients, gain normalization and
  SciPy's origin-padding convention. Extended `zpk2sos` with
  `pairing='minimal'`, which preserves roots without introducing new ones;
  documented why minimal sections may not be suitable for `sosfilt`.
  Added fixed tests and 792 minimal-design/759 SOS-inverse comparisons against
  SciPy 1.17, without additional dependencies.

- Added `scipy.signal.tf2zpk`, `zpk2tf`, and digital `zpk2sos` to
  `scifort_signal` through an independent Fortran module. Supports real and
  complex coefficients/roots, complex gains, conjugate-aware section pairing,
  and explicit status/result types. Added fixed Fortran regression tests and
  840 seeded SciPy 1.17 comparisons. Digital SOS degree-padding conventions
  and `pairing='minimal'` limitations are documented.

- Added SciPy-style analog and digital partial-fraction routines
  (`residue`, `residuez`, `invres`, `invresz`) through new
  `src/signal_partial_fraction.f90` and public `scifort_signal` exports.
  Supports complex coefficients, repeated poles, direct polynomial parts,
  configurable root clustering, and explicit result/status types. Includes
  an independent fixed Fortran regression and 540 seeded SciPy 1.17
  comparisons with inverse transfer-function round trips. The analog
  polynomial-root utility is now a public helper consumed by the new
  module; standalone build lists have been updated.


- Added SciPy-compatible `scipy.signal.czt` and `zoom_fft` routines for
  real/complex rank-one signals. The new `signal_czt.f90` module implements
  Bluestein convolution for unit-circle spirals using SciFort's existing radix-2
  FFT, and direct Horner evaluation for nonunit spirals. Supports custom
  output lengths, complex spiral parameters, scalar or two-frequency zoom
  intervals, optional sample rates and inclusive upper frequency endpoints.
  Includes a checked Fortran regression and 480 SciPy 1.17 randomized
  comparisons (maximum normalized error 2.06e-13); the complete 32-program
  signal test suite passes in checked and optimized builds.
- Extended the portable Remez intermediate-kind choice: try binary128 first,
  then an available >=18-digit real kind (e.g. x87 80-bit), then real64.
  Forcing GNU Fortran's 80-bit path reduces failures in the 240-case Remez
  comparison from 21 (real64) to 4. The pure real64 fallback remains limited
  on demanding high-order specifications.


- Improved high-order `scipy.signal.remez` convergence with extended-precision
  barycentric exchange, extrema selection and Chebyshev reconstruction where
  the compiler supports `selected_real_kind(30,400)`, retaining double-precision
  public coefficients and a portable double-precision fallback. Replaced the
  cubic extrema-selection search with a parity-prefix dynamic program. Added
  eight demanding high-order Hilbert/differentiator regression cases, fixed
  the standalone signal test runner's omitted `signal_remez` module, and
  expanded Remez verification diagnostics. All 240 randomized designs return
  coefficients on GNU Fortran 14.2; 239 can be compared directly with SciPy
  1.17 (which fails on the other). All 31 standalone signal tests pass with
  runtime checking and `-O3`.

- Stabilized `scipy.signal.firls` on narrow frequency bands using centered,
  cancellation-resistant cosine integrals, scaled band weights, and a
  rank-revealing Jacobi pseudoinverse fallback when fast Cholesky fails.
  Valid narrow-band filters now return finite taps instead of a spurious
  singular status. Added `test_signal_firls_rank.f90` and a 348-case
  conditioned-design comparison against SciPy 1.17 using weighted in-band
  errors. All 31 standalone signal tests pass in runtime-checked and
  optimized builds. The full FPM test suite remains unverified.

- Added `scipy.signal.remez` equiripple FIR design for symmetric bandpass,
  Hilbert-transform and differentiator filters (all four linear-phase types).
  The independent Fortran barycentric Remez exchange returns explicit
  nonconvergence status for ill-conditioned specifications. Added a standalone
  fixed regression test and randomized SciPy 1.17 comparison; the existing
  signal test suite and Python reference source lists include the new module.


- Added `scipy.signal.firwin2` (all four linear-phase FIR types, windowed
  frequency interpolation and discontinuous gain constraints), `firls`
  (weighted cosine least squares), and `minimum_phase` (homomorphic and
  modified Hilbert transforms). These use the existing spectral FFT/window
  code, with no additional external dependencies. Added deterministic
  Fortran tests, SciPy 1.17 randomized comparisons, and updated standalone
  signal reference build lists. An independent rank-revealing fallback for
  `firls` has since been added (above); SciFort still does not require LAPACK.
  The Hilbert method can be sensitive to spectral cancellations.


- Added independent analog IIR ZPK (`iirfilter_zpk`) and SOS
  (`iirfilter_sos`) output for Butterworth, Chebyshev I/II, and elliptic
  filters across all four band shapes. Added `zpk2sos_analog` for
  conjugate-symmetric proper analog systems and `sosfreqs` for response
  evaluation using descending s-plane polynomial coefficients. Existing
  `analog=.true.` BA results are unchanged. Verified 1,280 randomized
  designs and 12,800 response points against SciPy 1.17; added deterministic
  tests (including high-order BA-overflow avoidance) and updated the standalone
  verification build list.

- Added `scipy.signal.ellip` and `ellipord` to the signal layer, with all
  four band shapes and analog/digital support. `iirfilter` and `iirdesign`
  now accept the elliptic family (aliases: `elliptic`, `cauer`). A new
  dependency-free elliptic-prototype module uses real AGM/Jacobi/Carlson
  formulas, s-plane transformations and bilinear mapping; digital designs
  provide BA and SOS output, analog designs BA only. The existing order
  selection uses complete elliptic integrals. Added fixed reference tests,
  randomized SciPy 1.17 comparisons, and a fix to preserve the unmatched
  real root in odd-order SOS sections. Standalone signal test scripts now
  compile the new module. Full FPM validation remains outstanding.


- Extended `buttord`, `cheb1ord`, `cheb2ord`, and `iirdesign` for analog
  rad/s minimum-order filters plus two-edge digital/analog bandpass and
  bandstop specifications. Kept scalar `cutoff` and added allocatable
  `cutoffs(:)` to `iir_order_result`; analog calls reject `fs`. Bandstop
  order minimization uses independent bounded optimization. Added fixed
  edge/attenuation/invalid-input tests and a reproducible SciPy 1.17
  comparison (1,015 order/cutoff cases, 581 low-order BA coefficient
  comparisons). All 26 standalone signal tests pass with bounds checking
  and `-O3`. No new dependencies.

- Added `analog=.true.` for `butter`, `cheby1`, `cheby2`, and `iirfilter`,
  preserving existing digital calls. The new dependency-free analog IIR
  prototype and s-plane transformation module returns BA coefficients for
  lowpass, highpass, bandpass, and bandstop filters in rad/s. Added fixed
  SciPy reference tests, a 1,200-case randomized SciPy 1.17 comparison
  (maximum scaled coefficient error 4.55e-13), and updated standalone
  signal test scripts to compile the new module. Analog SOS/ZPK output
  remains unimplemented; analog minimum-order selection is now supported
  as described above.


- Added `scipy.signal.freqs`, `freqs_zpk`, and `freqz_zpk` with real/complex
  zero-pole-gain support, analog logarithmic frequency grids, digital uniform
  grids, and explicit angular or physical frequencies. The original analog
  polynomial-root grid estimator needs no new dependencies. Added a deterministic
  test and a 1,000-case SciPy 1.17 verification harness (16,519 bins, maximum
  absolute response difference 3.8e-9). Full FPM testing remains outstanding.
  The signal test suite passes all 24 programs in both runtime-checked and
  optimized GNU Fortran 14.2 builds. Existing signal-design and irregular-
  sampling SciPy comparisons also continue to pass.

- Added weighted `scipy.signal.lombscargle` power and normalized power,
  a separate complex-amplitude API `lombscargle_amplitude`, and
  `vectorstrength` for phase-locked event times. Added status reporting,
  deterministic tests, and 900 randomized SciPy 1.17 comparisons over
  11,201 frequencies, with no new dependencies.

- Added independently implemented `scipy.signal.sawtooth`, `square`, and
  `unit_impulse` (one-dimensional/two-dimensional) plus `max_len_seq` with
  default primitive taps for 2 through 32 register bits, explicit final state,
  streaming support, and input statuses. Periodic waves support scalar or
  per-sample widths/duty cycles; impulses use one-based indices. Added two
  deterministic test programs, documentation and provenance, and a reproducible
  1,300-case SciPy 1.17 comparison. All 22 standalone signal tests pass using
  GNU Fortran 14.2 with runtime checking. No new dependencies.


- Added `scipy.signal.lfiltic` for real and complex filters, reconstructing
  direct-form-II filter states from most-recent-first output/input histories,
  with optional zero input history and correct normalization for nonunit `a0`.
  Added `scipy.signal.group_delay` for real and complex transfer functions,
  with explicit or generated frequency grids, sample-rate units, and NaN plus
  status reporting for singular frequency bins. Added independent Fortran tests
  and a reproducible 600-case SciPy 1.17 comparison (maximum absolute error
  `4.97e-14`). No extra dependencies.

- Fixed the wall-clock Krylov refresh regression workload to use `system_clock`
  rather than CPU time, with a wider timing margin for coarse clocks and scheduling
  jitter. The solver's refresh policy is unchanged.
- Fixed bounds-check failures for empty relative-frequency inputs and terminal
  logrank event groups. Circular-mean regressions now compare angular distance
  across the equivalent zero/full-period endpoints.
- Fixed omitted optional array arguments in independent elementwise solvers by
  using private subroutine batch helpers; public interfaces are unchanged.
- Removed fixed COBYLA coordinates from the search space and restored them for
  objective/constraint callbacks. All-fixed problems are evaluated once. Corrected
  SLSQP iteration-limit reporting to count completed major iterations.
- Kept bounded finite-difference TRF regression accuracy checks while allowing one
  platform-dependent final evaluation; independently checked raw callback counts.
- Fixed optional callback dispatch in multivariate minimizers through explicitly
  interfaced private module procedures, avoiding GNU Fortran implicit-interface
  compilation errors first reported in SLSQP. Public APIs and numerical formulas
  are unchanged. Added final-gradient checks for analytic and finite-difference
  SLSQP callbacks.
- Made optional Krylov `inner_m` callback dispatch go through an explicitly
  interfaced module procedure, avoiding a compiler-dependent nested-callback
  implicit-interface diagnostic. Added a direct regression for the success,
  failure, and identity-preconditioner branches.
- Added `scipy.signal.sos2tf` for real second-order-section cascades, using
  polynomial convolution to recover numerator/denominator coefficients. A
  deterministic regression checks SciPy references, empty cascades, malformed
  sections, and Butterworth/Chebyshev designs. An optional verifier compares
  545 cases against SciPy 1.17.0. No new dependencies.

- Added independent `scipy.signal.cheby2`, `iirfilter`, `iirnotch`, `iirpeak`,
  `buttord`, `cheb1ord`, `cheb2ord`, and lowpass/highpass `iirdesign` to
  `scifort_signal`. The inverse-Chebyshev analog prototype supports digital
  low/high/bandpass/bandstop transformation with SOS and BA outputs;
  order-selection routines use bilinear-prewarped pass/stop edges. No new FPM
  dependencies. Added a deterministic Fortran regression program and optional
  1,350-case SciPy 1.17 comparison script. 15 standalone signal tests pass in
  GNU Fortran checked and optimized builds. Full dependency-backed `fpm test`
  remains outstanding.


- Implemented `scipy.signal.cheby1` (analog Chebyshev-I prototype, four digital
  filter types, transfer-function and SOS coefficients), the default
  Chebyshev-I IIR branch of `decimate` (zero-phase and causal, real/complex),
  and Fourier `resample` (real/complex 1D, time/frequency-domain complex input,
  optional spectral windows and timestamps). The even-length Nyquist bin is
  split or combined explicitly when resampling. The existing FIR decimator
  is retained; omitting `ftype` now selects IIR as in SciPy. No new runtime
  dependencies. Added two Fortran regression programs and an optional
  cross-language verification script. Full FPM testing remains outstanding.

- Added three SciPy-style one-dimensional multirate signal algorithms in
  `scifort_signal_resample`: `upfirdn`, `resample_poly`, and FIR
  `decimate` (with IIR decimation added in a later increment). Real/complex signals, real/complex `upfirdn` kernels,
  rational-ratio reduction, zero-phase FIR alignment, Kaiser(5) lowpass
  design, user taps, several boundary modes, and causal FIR downsampling
  are supported. The original increment required `ftype='fir'`;
  the IIR mode and default were added subsequently. Reproducible reference
  driver and 1,120 randomized SciPy 1.17 comparisons pass (maximum absolute
  difference 8e-15); 12 standalone signal tests pass in both checked and
  optimized GNU Fortran builds. This source is original SciFort MIT code;
  no external dependencies were added. Full FPM testing is outstanding.


- Added `scipy.signal.firwin`, `freqz`, `sosfreqz`, and `butter` through
  `scifort_signal_design`. Supports windowed-sinc multiband FIR design,
  real/complex IIR frequency responses, and low/high/bandpass/bandstop
  Butterworth filters with transfer-function coefficients and stable SOS
  output. No new dependencies. Added deterministic regression tests and
  1,620 randomized SciPy 1.17 comparisons; full FPM testing is outstanding.
  Fixed conjugate/real-pole pairing for wide bandstop designs.

- Added four SciPy-style digital filtering algorithms to `scifort_signal`:
  `sosfilt`, `sosfilt_zi`, `filtfilt`, and `sosfiltfilt`, including real and
  complex 1D inputs, cascaded transposed-direct-form-II delay states,
  step steady-state initialization, and odd/even/constant/no endpoint padding.
  Two new deterministic Fortran tests and an optional 2,100-case SciPy 1.17
  comparison harness validate the supported API. The previous signal tests
  are unchanged. There are no new FPM dependencies. Full dependency-backed
  `fpm test` remains outstanding.


- Extended `scifort_signal` spectral estimators with complex(dp), two-sided
  one-dimensional `periodogram`, `welch`, `csd`, `coherence`, `stft`, and
  `spectrogram` overloads. Implemented `fftconvolve` and `oaconvolve` for
  real/complex 1D arrays and full/same/valid modes, with reusable overlap-add
  transforms for long signals. No new dependencies. Two new deterministic
  Fortran regressions and optional SciPy reference scripts validate 960
  complex spectral and 1,000 FFT-convolution cases against SciPy 1.17.0.
  Full dependency-backed `fpm test` is still outstanding.


- Added seven spectral-analysis algorithms to `scifort_signal`: `get_window`,
  `periodogram`, `welch`, `csd`, `coherence`, `stft`, and `spectrogram` for
  real 1D signals, with SciPy-style windows, detrending, Welch mean/median
  averaging, one-sided scaling, and time-frequency results. The standalone
  FFT uses radix-2 butterflies and a direct DFT fallback for other lengths.
  Jason Christopherson's MIT-licensed `spectrum` was consulted for Welch
  segmentation and spectral architecture, with explicit attribution.
  Dedicated checked/optimized Fortran regressions and 1,640 randomized
  SciPy 1.17 comparisons validate the supported subset; no new dependencies.
  The full dependency-backed FPM test suite remains outstanding.


- Expanded `scifort_signal` with four independent 1D SciPy-style algorithms:
  `savgol_coeffs`, `savgol_filter`, `lfilter`, and `lfilter_zi`.
  Savitzky-Golay uses reorthogonalized polynomial least squares, supports
  even/odd windows, derivatives and five edge modes. Real/complex filtering
  provides normalized direct-form-II FIR/IIR recurrences, streamable delay
  state, and steady-state step initialization. Two new Fortran regression
  programs pass with GNU Fortran 14.2 in checked and optimized builds;
  the optional SciPy 1.17.0 comparison checks 2,380 randomized cases.
  The full dependency-backed FPM test suite remains outstanding.


- Added seven `scipy.signal` algorithms via the new `scifort_signal` module:
  `find_peaks`, `peak_prominences`, `peak_widths`, `detrend`, `convolve`,
  `correlate`, and `correlation_lags`. Features include plateau-aware peak
  properties, prominence/window and width filtering, 1D linear/piecewise
  detrending, and direct real/complex convolution/correlation in full/same/
  valid modes. All new numerical source is independently authored under MIT.
  SciPy 1.17.0 comparisons: 700 randomized distinct-height peak signals and
  450 detrending and 900 real/complex convolution/correlation cases; a
  reproducible optional Python/SciPy reference harness is included.
  Dedicated Fortran tests include real/complex signal
  arithmetic and invalid arguments. Signal APIs are documented as a subset
  of SciPy's full functionality. Full dependency-backed `fpm test` remains
  outstanding in this environment.

- Added two new `scipy.optimize` algorithms: rectangular `linear_sum_assignment`
  (minimum/maximum matching, forbidden edges and infeasibility reporting)
  and weighted `isotonic_regression` (increasing/decreasing PAVA with pooled
  weights and block boundaries). Both are independent, MIT-licensed Fortran
  implementations with SciPy 1.17-compatible numerical behavior and one-based
  indices. Dedicated tests pass in checked and optimized GNU Fortran 14.2;
  600 randomized assignment cases and 650 randomized isotonic fits were
  separately compared against SciPy 1.17.0. Full dependency-backed `fpm test`
  remains outstanding in the current environment.

- Added optional wall-clock timing for adaptive Krylov preconditioner
  refresh, selected by `inner_m_diag_update_wall_clock=.true.` alongside
  `inner_m_diag_update_timed_cost=.true.`. Standard Fortran `system_clock`
  with int64 counters replaces process `cpu_time` for setup and inner-solve
  measurements. Counter failures/rollover and sub-resolution measurements
  suppress timing-based early refresh. Existing CPU-time and default modes
  remain unchanged. Added a wall-clock regression and an OpenMP quadrature
  benchmark mode. GNU Fortran 14.2: 68/68 Krylov-focused tests in checked
  and optimized builds; all 84 optimized and 84 checked OpenMP quadrature
  cases converged (42 at one and 42 at four threads in each build).
  Benchmark timings depend on hardware, thread scheduling, and load.


- Added `benchmark_krylov_quadrature_refresh.f90`, an independently
  checked nonlinear rank-two benchmark with per-variable exponential-response
  midpoint quadrature and analytic forward/transpose Jacobian actions.
  OpenMP parallelizes the actual quadrature loops; a serial build needs no
  OpenMP dependency. It compares fixed, product-count, and CPU-timed refresh
  policies at quadrature orders 96, 640, and 2048 and nonlinearity rates 0.9
  and 2.0, recording process CPU and elapsed wall seconds. All 72 optimized
  cases (36 each with one and four OpenMP threads) and 72 checked-debug
  cases converged, and all 67 existing Krylov-focused tests passed in both
  configurations (GNU Fortran 14.2, test-only MINPACK stub). Finite-difference
  and adjoint checks verify the quadrature Jacobian. This is an original,
  synthetic numerical workload rather than an application benchmark;
  CPU-timed refresh can choose different schedules when threaded callbacks
  change process CPU accounting. No solver APIs or defaults were changed.

- Added `benchmark_krylov_timed_refresh.f90`, a 32-case, 128-variable
  nonlinear comparison of fixed, manually costed, product-count-costed,
  and CPU-timed Krylov preconditioner refresh. Two synthetic CPU-intensive
  callback profiles distinguish expensive forward from expensive transpose
  Jacobian actions. Each case checks root accuracy and callback accounting;
  a high minimum timing threshold confirms the guard disables early timed
  updates. Both 32-case benchmark runs and all 67 existing Krylov-focused
  tests passed under GNU Fortran 14.2 with optimized and checked-debug flags
  (temporary MINPACK stubs). No solver APIs, defaults, or numerical methods
  were changed. The accompanying CSV contains single-run process-CPU times;
  those results are hardware- and load-dependent.

- Added a separate, opt-in CPU-timed preconditioner rebuild cost for
  work-aware Krylov refresh (`inner_m_diag_update_timed_cost`) with a
  configurable minimum CPU measurement duration. The policy measures setup
  and inner-solve CPU seconds, compares solve-time growth with setup cost,
  and ignores absent or low-resolution timing samples. Existing measured
  product counts, manually declared costs, defaults, and APIs remain
  unchanged when the option is not enabled. Added a CPU-load-based
  regression exercising early rebuild, noise suppression, forwarding,
  and invalid options. GNU Fortran 14.2: 67/67 Krylov tests pass in both
  checked-debug and optimized builds using temporary MINPACK stubs.

- Extended `inner_m_diag_update_measured_cost` beyond nonsymmetric
  rank-two discovery to every automatic diagonal and low-rank Krylov
  preconditioner. Setup accounting now includes exact and randomized
  diagonal probes, rank-one verification, known-right rank-two fitting,
  accepted/rejected rank-two sketches, dense or streamed discovery,
  symmetric low-product fitting, and exact fallback/retry products.
  The option is still opt-in, with the same public interface and defaults.
  Added a multi-mode regression including finite-difference accounting;
  updated the existing validation case. GNU Fortran 14.2 checked and
  optimized Krylov-focused tests pass (66/66 in each configuration), and
  all 132 existing nonlinear work benchmark cases still converge with
  identical non-timing output compared with the preceding archive.

- Added opt-in measured refresh setup costs for nonsymmetric rank-two
  Krylov discovery. The solver counts actual forward and transpose
  fitting/verification products, including n additional coordinate Jv
  products when a fit is rejected and exact-diagonal fallback runs.
  `inner_m_diag_update_measured_cost` and
  `inner_m_diag_update_transpose_cost` choose a weighted cost in the
  existing work-aware refresh rule; manual cost and all defaults are
  unchanged. Added deterministic regression and two benchmark policies
  comparing 1x and 16x transpose weights on the nonlinear work suite.
  All 65 Krylov-focused test programs passed in both checked-debug and
  optimized GNU Fortran 14.2 builds; 132/132 benchmark cases converged.


- Added a reproducible 108-case benchmark of work-aware Krylov diagonal
  preconditioner refresh on more demanding nonlinear, ill-scaled rank-two
  systems, with independent forward/transpose product counts and modeled
  setup costs. Work-based refresh triggers even when residual-stall refresh
  does not; low rebuild cost can reduce forward products, but expensive
  transpose products can outweigh the savings. Added a separate nonlinear
  regression for early refresh, moderate/high rebuild costs, and accurate
  convergence through the built-in GMRES path. All 64 Krylov-focused tests
  passed in GNU Fortran 14.2 checked-debug and optimized builds. All 108
  benchmark cases converged. No solver API or default behavior changed.

- Added opt-in, work-aware adaptive refresh for Krylov diagonal
  preconditioners. `inner_m_diag_update_work_growth=g` compares inner
  preconditioner-action counts to the last refresh, while the optional
  nonnegative `inner_m_diag_update_rebuild_cost=c` suppresses refreshes
  whose estimated setup cost outweighs extra inner work. This augments,
  rather than replaces, residual-stall and maximum-age scheduling.
  Implemented for `root` and `newton_krylov` with unchanged defaults.
  Added deterministic regression for early refresh, cost suppression,
  baseline reset, forwarding, zero-work handling, and invalid arguments.
  All 63 Krylov-focused tests pass with GNU Fortran 14.2 in checked-debug
  and optimized builds using temporary dependency stubs. A workload-wide
  performance advantage has not been established.

- Added opt-in adaptive Krylov diagonal preconditioner refresh,
  `inner_m_diag_update_adaptive=.true.`, with a maximum refresh age from
  `inner_m_diag_update_interval` and early refresh on weak *accepted*
  nonlinear residual reduction. The optional validated
  `inner_m_diag_update_stall_ratio` controls sensitivity; fixed scheduling
  and all defaults are unchanged. Added deterministic schedule, forwarding,
  age-reset, and invalid-parameter regressions. Extended the nonlinear
  rank-two benchmark to 84 configurations with adaptive thresholds 0.8/0.95.
  All 62 Krylov-focused tests pass in checked and optimized builds; all
  84 benchmark cases converge. These benchmark systems did not trigger
  earlier refreshes, so no performance improvement is claimed.

- Added a deterministic 60-case nonlinear Krylov refresh benchmark with
  changing nonsymmetric rank-two coupling and `D(x)`, dimensions 32/128,
  two diagonal condition ranges, three nonlinear rates and refresh
  intervals 0/1/2/4/8. Records convergence, independent forward/transpose
  products, base callbacks, modeled cost and CPU time. Extended the
  nonlinear base regression to verify scheduled refresh counts, distinct
  accepted iterates, and delayed callback failure for longer intervals.
  All 61 focused Krylov tests pass in checked and optimized builds;
  all 60 benchmark configurations converge. No solver API or defaults changed.

- Added an opt-in point-dependent uncoupled diagonal callback,
  `inner_m_diag_rank_two_nonsymmetric_base_fun=`, to nonsymmetric rank-two
  Krylov discovery. The callback is evaluated at each scheduled preconditioner
  refresh and forwards through `root` and `newton_krylov`; it is mutually
  exclusive with the existing fixed base array. Failed or invalid callbacks
  terminate explicitly, while verified structural mismatch retains the
  exact-diagonal fallback. Added genuinely nonlinear, changing-Jacobian
  regression cases including refresh intervals, callback failures,
  inaccurate diagonals and early option validation. Defaults are unchanged.

- Added optional cost-aware held-out transpose verification for
  nonsymmetric rank-two Krylov discovery using
  `inner_m_diag_rank_two_nonsymmetric_adaptive_transpose_cost=c`.
  The solver chooses one or two additional J^T*v checks according to the
  supplied relative callback cost, problem dimension and observed fit
  quality. The existing fixed-check Boolean and default operation are
  unchanged; conflicting, nonpositive or nonfinite settings fail early.
  Added adversarial, forward-borderline, forwarding, failure and invalid
  parameter tests. Extended the discovery benchmark to 56 cases with
  adaptive policies at transpose costs 1 and 16.


- Added opt-in `inner_m_diag_rank_two_nonsymmetric_verify_transpose=.true.`
  for low-product rank-two Krylov discovery. In addition to the existing
  four forward and two transpose Jacobian products, the method now checks
  two held-out transpose directions (four transpose products total).
  A detected inconsistency reverts to exact-diagonal preconditioning;
  an unsuccessful transpose callback remains a solver error. Default
  product counts and existing behavior are unchanged. Added an adversarial
  callback regression that is invisible along both fitting directions.
  Added a 40-case benchmark comparing exact-diagonal, dense, streamed,
  low-product, and transpose-verified discovery, with forward/adjoint counts
  and cost models for relative transpose expenses of 1, 4, and 16.


- Added opt-in low-product **nonsymmetric** rank-two preconditioning
  for Krylov `root` and `newton_krylov`. With a known uncoupled
  diagonal, `inner_m_diag_rank_two_nonsymmetric_base=a` and an analytic
  `jacobian_transpose_product=` callback, two Jv and two J^T v products
  fit the correction and two independent Jv products check it. Accepted
  systems use O(n) discovery memory; rejected or unsafe fits revert to
  the exact-coordinate diagonal. Added nonsymmetric, 128-variable
  ill-scaled, finite-difference, callback-error and option-validation
  regressions. The new options are SciFort extensions and preserve defaults.

- Added opt-in low-product symmetric rank-two preconditioning for Krylov
  `root` and `newton_krylov`. Supplying a known uncoupled diagonal as
  `inner_m_diag_rank_two_symmetric_base=a` enables a deterministic
  Nystrom fit with two Jacobian-vector products and two independent
  holdout checks (four products in accepted cases). The O(n)-storage
  algorithm rejects incompatible or unsafe fits and falls back to the
  exact Jacobian diagonal. Added analytic/finite-difference, 128-variable
  ill-scaled, rank-three contamination, callback-failure, and validation
  tests. Existing defaults and discovery modes remain unchanged.

- Added opt-in `inner_m_diag_rank_two_stream=.true.` for Krylov rank-two
  unknown-factor discovery. Retains no dense Jacobian: up to eight evenly
  spaced pilot coordinate columns identify independent off-diagonal
  anchors, and a single streamed pass through n columns learns both
  factors and checks all off-diagonal entries. The exact diagonal is
  reused on every rejected fit. Setup needs n Jacobian-vector products
  and O(n^2) arithmetic as before, but only O(n) discovery storage.
  Added tests for scaled 128-variable systems, contamination outside
  the pilot columns, finite differences, rank-deficient fits and API
  validation. Default discovery remains unchanged.

- Added opt-in `inner_m_diag_rank_two_discover=.true.` to discover both
  unknown two-column coupling factors in Krylov `root` and `newton_krylov`.
  Reuses n exact coordinate Jacobian products for full off-diagonal
  factorization, tests every off-diagonal entry, and uses a guarded
  Woodbury inverse only when the completion is identifiable and safe.
  Rejected, small, singular, or nonfinite fits retain exact diagonal
  preconditioning. The new mode is mutually exclusive with supplied
  right factors and rank-two sketching, and requires exact auto diagonal
  estimation. Defaults are unchanged. Added analytic/numerical,
  small-system, ill-scaled, adversarial, and callback-failure regressions.

- Added opt-in `inner_m_diag_rank_one_discover=.true.` for Krylov `root`
  and `newton_krylov` when the right coupling vector is not known. The
  method stores the complete exact-coordinate Jacobian already needed
  for diagonal estimation, reconstructs a rank-one off-diagonal factor
  using independent pivots, verifies all off-diagonal entries, and
  applies the guarded Sherman-Morrison inverse only to safe fits.
  Degenerate/sparse, nonsingular-base failures, singular Woodbury factors,
  incompatible coupling, and unsafe fits retain the original diagonal
  preconditioner. No extra Jv products beyond n exact coordinate probes;
  O(n^2) storage and comparison work; option disabled by default. The
  existing known-right and rank-two methods are unchanged. Added tests
  for analytic and numerical Jv, mixed-sign/zero weights, ill scaling,
  fallbacks, validation, and callback failure.

- Added opt-in `inner_m_diag_rank_two_sketch=.true.` for checked
  rank-two Krylov preconditioning. Three deterministic bounded sign
  probes recover the Jacobian diagonal in the known-right-factor
  `D + U*V^T` model, subject to worst-row 3-by-3 identifiability checks.
  The existing two fit and two verification probes then validate the
  Woodbury correction, so accepted setups use seven Jv calls instead of
  `n+4`. Degenerate fits or failed checks retry the exact-diagonal path;
  Jacobian callback failures remain fatal. Defaults are unchanged.
  Added analytic, finite-difference, error/fallback tests and a coupled
  benchmark comparing five policies across four size/scaling regimes.

- Added an opt-in checked rank-two Woodbury preconditioner to Krylov
  `root` and `newton_krylov`. With exact automatic diagonal estimation,
  `inner_m_diag_rank_two=.true.` and a known `(n,2)`
  `inner_m_diag_rank_two_right=V`, it fits `J = diag(a) + U*V^T` from two
  bounded Jacobian-vector products, checks two independent directions,
  and safely applies a two-by-two Woodbury inverse. Four additional Jv
  evaluations per refresh; rejects non-identifiable or numerically unsafe
  fits and retains the existing diagonal fallback. Opt-in, mutually
  exclusive with rank one, and unsupported with built-in MINRES. Added
  analytic, finite-difference, extreme-scale, singular, invalid-input,
  small-system, and adversarial-coupling regression cases. This is an
  original SciFort extension; omitted options preserve existing behavior.

- Extended the opt-in checked Krylov rank-one preconditioner with
  `inner_m_diag_rank_one_right=v` for a known, fixed right coupling vector:
  `J = diag(a) + u*v^T`. The estimator uses a bounded sign-weighted Jv
  fit, normalizes user weights safely, and applies a weighted
  Sherman-Morrison inverse after two structural checks. For systems of at
  least three variables, checks avoid sign-duplicate probe directions. It retains the
  existing diagonal fallback, exact-diagonal and non-MINRES constraints,
  O(n) memory, and three additional Jv products per refresh. Omitting the
  new vector preserves the original mean-field behavior. Added analytic,
  finite-difference, mixed-sign, scaled-weight, shape/error, and fallback
  regression coverage; this is a SciFort extension rather than SciPy API.

- Added opt-in `inner_m_diag_rank_one=.true.` for Krylov `root` and
  `newton_krylov`: an exact-diagonal, checked mean-field
  `J = diag(a) + u*ones^T` preconditioner with a Sherman-Morrison inverse.
  Three additional Jv products fit and verify the model; unsafe or
  inconsistent models fall back to the ordinary safeguarded diagonal.
  The default is unchanged and built-in MINRES rejects the nonsymmetric
  option. Added tests for one-step GMRES, ill scaling, callback failure,
  invalid options, and analytic/finite-difference evaluation counts.

- Progressive Krylov Rademacher diagonal verification now prefers **unchecked
  coordinates** when a rotated sample overlaps a cached site. A local
  path-compressed successor table with cyclic wraparound prevents repeated
  selections and avoids costly searches for large verification sets. Cached
  exact coordinates are still revalidated against expanded estimates; exact
  fallback reuses every product. Added a 64-variable coupled-Jacobian
  counterexample where high-coverage, narrow-stratum verification previously
  missed a 30% biased diagonal at all probe stages. New analytic, two-point,
  and three-point tests cover the recovered exact diagonal and Jv work counts.
  Public APIs, defaults, and fixed-probe behavior are unchanged.

- Progressive verified Krylov diagonal probing now rotates exact-coordinate
  spot checks at each probe-budget expansion using deterministic base-2
  bit-reversed offsets within equally spaced strata. It rechecks all earlier
  exact samples from a per-estimate cache without repeated Jv calls. This
  detects localized probe bias that could pass repeated checks at the original
  coordinates; exact fallback still reuses every sampled coordinate. No public
  API or default/fixed-budget path changes. Added two- and three-stage,
  callback-failure, full-coverage and finite-difference count regressions.

- Progressive Krylov diagonal probing now caches exact coordinate spot checks
  across randomized-probe expansions and reuses them on exact fallback. This
  avoids duplicate Jacobian-vector products at a fixed nonlinear iterate,
  while preserving estimate values, default/fixed-probe paths and public APIs.
  Added deterministic biased-probe and finite-difference count regressions.

- Added opt-in `inner_m_diag_probe_max=` for **progressively verified seeded
  Rademacher diagonal estimation** in Krylov `root` and `newton_krylov`. When
  a low-budget estimate fails split-batch or coordinate verification, the
  estimator doubles the probe budget up to the supplied cap, reusing all
  earlier Jacobian-vector products. Persistent disagreement still triggers
  exact coordinate estimation. Omitting the argument preserves all previous
  behavior, including exact-fallback evaluation counts. Added tests for
  intermediate acceptance, budget exhaustion, callback failures, fixed-path
  compatibility, validation, and forward/centered finite-difference counts.
  A dense rank-one coupling benchmark and CSV document both work reductions
  and remaining coordinatewise approximation errors.

- Added optional `inner_m_diag_probe_scheme='rademacher'` with reproducible
  `inner_m_diag_probe_seed=` for approximate diagonal estimation on dense or
  irregularly coupled Krylov systems. The existing exact estimator and Walsh
  probe order remain the defaults. New `inner_m_diag_verify_samples=` and
  `inner_m_diag_verify_rtol=` controls perform coordinate spot checks and,
  for seeded probes, a split-batch consistency check without extra Jv calls;
  suspect estimates fall back to exact coordinate probing. A regression
  documents both reduced mean aliasing and potentially *larger* worst-case
  error for random probes, with validation and evaluation-count checks.
  Added `benchmark_krylov_random_diagonal.f90` and reference CSV. The
  extension is original SciFort code; complete FPM validation remains pending.

- Added opt-in `inner_m_diag_probes=` for deterministic Walsh-sign estimates
  of the inverse-Jacobian diagonal in Krylov `root` and `newton_krylov`.
  The default remains exact coordinate probing; powers of two can recover
  narrow-banded diagonals exactly using far fewer Jv calls, but nonlocal
  couplings may alias. Added regression tests for exactness, known aliasing,
  finite-difference evaluation counts, refresh intervals, and input failures.
  On a 2,048-variable tridiagonal nonlinear example, analytic Jv calls fell
  from 2,073 to 33 with 8 probes, with matching convergence and solution
  accuracy. GNU Fortran 14.2 debug and optimized focused tests passed;
  dependency-backed FPM testing remains outstanding.

- Added opt-in `inner_m_diag_auto=` for Krylov `root` and `newton_krylov`
  with configurable `inner_m_diag_update_interval=` and relative
  `inner_m_diag_floor=`. Diagonals are obtained from analytic Jv actions or
  finite differences, refreshed at accepted nonlinear iterates, and inverted
  with zero/tiny-entry safeguards. The implementation preserves solver
  defaults and rejects conflicting preconditioners before residual calls.
  New coupled-system regression covers both entry points, update cadence,
  MINRES sign safety, finite-difference counts and callback failures. All
  44 focused Krylov/root programs pass GNU Fortran 14.2 checked debug and
  optimized builds. A new benchmark and reference CSV demonstrate both
  convergence gains and the O(n)-probe estimation overhead. Full FPM
  validation remains pending.

- Added optional static diagonal inverse-Jacobian preconditioning through
  `inner_m_diag=` for Krylov `root` and `newton_krylov`. All six built-in
  Krylov methods are supported, including restarted GMRES/LGMRES; custom
  inner solvers receive the same array action through their preconditioner
  callback. Preflight validation rejects wrong lengths, nonfinite/zero
  entries, conflicting user preconditioners, and nonpositive MINRES factors;
  overflow of a preconditioner product returns a callback failure. The
  previous absent-option algorithms are unchanged. Added a regression and
  benchmark comparing this convenience option with the existing callback,
  including a changing nonlinear Jacobian and systems with diagonal condition
  number `1e8`. All 43 focused programs passed GNU Fortran 14.2 checked
  debug and optimized builds; dependency-backed FPM validation remains open.


- Added a matrix-free large-system comparison for full, fixed-restart, and
  adaptive restarted LGMRES. New benchmark, 78-case reference CSV, and
  regression cover 64-65,536 variables, condition numbers `1e4`-`1e8`,
  true residuals and solution errors, analytic Jacobian products, and
  optional scaled stopping norms. At 65,536 variables and condition `1e6`,
  adaptive restarting used 24,004 KiB peak process RSS versus 59,832 KiB
  for full-width LGMRES, with 276 versus 164 analytic Jv calls. The tests
  also document why an absolute residual can be misleading for a severely
  ill-conditioned nonlinear system. No solver defaults or implementation
  changed. All 42 focused root/Krylov tests passed in checked debug and
  optimized GNU Fortran 14.2 builds; complete FPM validation remains pending.

- Added an opt-in adaptive restarted-LGMRES extension,
  `inner_restart_adaptive=.true.`, with `inner_restart_max=` as the maximum
  Arnoldi subspace width. If a cycle reduces the true linear residual by less
  than 20%, the next cycle can double its width without exceeding the total
  `inner_maxiter` budget. By default growth is capped at four times the
  initial `inner_restart`; omitted or disabled options retain fixed-restart
  and SciPy-compatible behavior. This also rescues a genuinely stagnant
  LGMRES(1) case by growing to width two. Added a Fortran regression and a
  reproducible six-case benchmark. All 41 focused root/Krylov tests pass with
  GNU Fortran 14.2 in checked debug and optimized builds. The complete
  dependency-backed `fpm test` remains outstanding.

- Extended opt-in `inner_restart=` to built-in **LGMRES**. Restart cycles use
  no more than the specified number of new Krylov vectors and preserve the
  retained `outer_k` augmentation history until the complete Newton step is
  formed. Every cycle checks the original linear residual, and the final
  check's exact Jacobian image is reused for candidate pruning. Restart mode
  avoids holding raw Arnoldi products solely for candidate-image reuse. The
  absent and full-width option paths keep the previous algorithm and counts.
  Added a 40th focused root/Krylov regression and a large matrix-free memory
  benchmark; all 40 focused tests passed in GNU Fortran 14.2 checked debug
  and optimized builds. On the 250,000-variable benchmark, two-run peak RSS
  was roughly 288,460 KiB without restarts versus 73,500 KiB with
  `inner_restart=8`, in exchange for 105 versus 90 analytic Jv calls.

- Added optional `inner_restart=` to built-in GMRES in `root(method='krylov')`
  and `newton_krylov`. It caps the per-cycle Arnoldi subspace while
  treating `inner_maxiter` as the total step budget, and recomputes the true
  linear residual after each restart. The original single-cycle default,
  including its work counts, is unchanged. Added a checked Fortran
  regression and full-solve peak-memory benchmark. All 39 focused Krylov/
  root programs passed in GNU Fortran 14.2 debug and optimized builds.
  For 100,000 unknowns, the reproducible `inner_maxiter=64` versus
  `inner_restart=8` comparison reduced peak RSS from 111,472 to 25,612 KiB,
  with increased analytic Jacobian-vector calls (26 versus 101).

- Added an opt-in raw Arnoldi-product memory budget for built-in LGMRES:
  `outer_candidate_cache_max_mib=` is forwarded by `root(method='krylov')`
  and `newton_krylov()`. With `outer_candidate_image_reuse=.true.`, a
  projected cache exceeding the positive MiB budget triggers an exact
  analytic Jacobian-vector callback for the candidate rather than storing
  its raw Arnoldi products. Zero or omission retains unlimited caching and
  existing evaluation counts; invalid combinations are rejected before
  residual evaluation. A new regression exercises both caching modes, the
  budget boundary across nonlinear iterations, and public forwarding.
  All 39 focused root/Krylov programs pass in GNU Fortran 14.2 checked
  debug and optimized builds. On the included 250,000-variable benchmark,
  a 1 MiB raw-cache cap reduced peak RSS from approximately 174.9 MiB to
  136.6 MiB in exchange for five additional analytic Jv calls.

- Reduced cached LGMRES image-pruning memory without changing the public
  solver interface: the pre-Arnoldi direction scratch now borrows the first
  persistent Krylov `z` column, and the Jacobian-image matrix contains only
  the current history (not a spare candidate column). Removed a redundant
  full-matrix zero-fill: every live image is produced before use. The
  zero-history case uses a zero-column Fortran array. Added a short-history
  and FIFO regression. All 38 focused root/Krylov programs pass under GNU
  Fortran 14.2 debug and optimized checked builds. Matched full-solve peak
  RSS measurements with `outer_k=8` improved by approximately 3% at
  100,000-500,000 unknowns, with identical convergence and evaluation counts.

- Reused the already allocated LGMRES Jacobian-image buffer for the
  current-iterate orthogonal pruning basis when `outer_image_basis_reuse`
  is enabled without `outer_image_reuse`. This removes a separate
  `n*outer_count` allocation. The in-place path no longer needs the 8 MiB
  extra-basis cap, while raw-image caching retains that cap and the
  original algorithm. Cached candidate screening accepts an optional
  caller-owned basis and has a guarded legacy fallback if the stream
  becomes invalid. A new direct regression covers dependence, FIFO
  replacement, and candidate rejection. All 37 focused root/Krylov tests
  passed in checked debug and optimized builds. Two paired large-system
  runs at 250,000 variables reduced peak RSS from about 251.2 MiB to
  241.7 MiB (3.8%), with identical results and evaluation counts;
  full-solver CPU times were similar, with no consistent speedup.


- Added opt-in `outer_image_basis_reuse=.true.` to built-in Krylov LGMRES
  with current-Jacobian image pruning. The candidate Newton direction is
  tested against the same two-pass image basis built during history filtering,
  instead of rebuilding the basis after Arnoldi. One-shot validity and
  per-iterate resets prevent reuse across changing Jacobians; extra basis
  storage is limited to 8 MiB, with transparent legacy fallback above that
  limit. The default algorithm and Jacobian evaluation counts remain
  unchanged. New direct and public regressions, plus 16-variable
  changing-Jacobian stress comparisons, pass GNU Fortran 14.2 checked debug
  and optimized builds. A 192-vector, 48-history pruning benchmark showed
  approximately 1.8x less CPU time for the filtering/insertion stage,
  excluding Jacobian products and Arnoldi iteration.


- Added opt-in `outer_fifo_prune_fast=.true.` to the built-in LGMRES
  `root(method='krylov')` and `newton_krylov()` interfaces when
  `outer_prune_tol > 0` and `outer_k >= 2`. A full augmentation history is
  maintained by a small QR factor and Givens-rotation downdates after FIFO
  eviction, rather than repeatedly rebuilding its entire orthonormal span.
  Rank-deficient histories automatically fall back to the previous
  reconstruction method; the cached factor is refreshed periodically.
  The default and all previous root options remain unchanged. Two new tests
  cover FIFO equivalence and public forwarding/invalid combinations; all 34
  focused Krylov/root programs pass in GNU Fortran 14.2 debug and `-O2`
  checked builds. A 192-element, 48-history, 200-proposal benchmark ran
  7.8x faster for the isolated pruning loop (not overall nonlinear solving).

- Added an incremental two-pass orthogonal basis for append-only
  solution-space LGMRES augmentation pruning. The built-in root solver now
  reuses this per-solve basis while `outer_count < outer_k`, switching back
  to the original reconstructed-span/FIFO routine once the history is full.
  Extra basis storage is capped at 8 MiB (1,048,576 `real(dp)` elements);
  larger systems retain the previous path. Default SciPy-compatible root
  behavior and all public signatures remain unchanged. Added the low-level
  `retain_lgmres_direction_stream` and `reset_lgmres_direction_pruning_stream`
  entry points, a deterministic equivalence regression, and a standalone
  benchmark. All 32 focused root/Krylov programs pass in checked GNU Fortran
  14.2 debug and optimized builds; the previous changing-Jacobian comparison
  remains byte-for-byte identical.

- Made the repeated current-Jacobian LGMRES image-pruning pass incremental.
  As retained directions are screened, a two-pass orthonormal image basis is
  extended rather than reconstructed for each new image. Reset the basis at
  each new nonlinear iterate; keep the existing candidate-insertion/FIFO path
  and public `root`/`newton_krylov` options unchanged. Reuse retired Arnoldi
  storage where available, and reserve only the active history length when
  separate scratch is needed. Added direct fixed/automatic filtering parity
  tests and a repeatable long-history benchmark. All 31 focused root/Krylov
  programs passed in checked debug and optimized GNU Fortran 14.2 builds;
  the optimized changing-Jacobian stress output matched the previous archive
  byte for byte.

- Reused allocation-heavy LGMRES augmentation-pruning work arrays within
  nonlinear solves. Both solution-direction and Jacobian-image filters accept
  an optional `lgmres_pruning_workspace`; built-in LGMRES reuses the unused
  Arnoldi basis as orthogonalization scratch, avoiding a second `n*outer_k`
  workspace. Jacobian images and candidate vectors now persist across Newton
  iterations in solver-local storage. No public `root` options, default
  pruning decisions, or evaluation counts changed. Added a direct equivalence
  regression, an optional pruning microbenchmark, and a matched large-system
  process-memory comparison. All 30 focused root/Krylov tests passed in
  checked debug and optimized GNU Fortran 14.2 builds.


- Reduced large-`outer_k` memory overhead in Krylov LGMRES without changing
  its public interface. The Arnoldi basis now grows geometrically with the
  number of active retained vectors; Jacobian-image pruning allocates only
  enough images for active history plus one candidate. Preconditioners write
  directly into the stored Krylov direction instead of making two full-vector
  copies. Added a capacity-independence regression, extended the matrix-free
  memory benchmark with optional `outer_k` and preconditioning, and compared
  GNU Fortran 14.2 peak RSS with the previous snapshot. All 29 focused
  root/Krylov tests passed in both debug and optimized checked builds.

- Reduced GMRES/LGMRES Arnoldi allocation churn without changing user-visible
  solver options: reuse the basis/Hessenberg/least-squares workspace across
  Newton iterations, use the next Arnoldi basis column as J*v scratch, and
  retain only ordinary Krylov directions in the separate `z` workspace
  (augmentation directions already live in `outer_vec`). GMRES also reuses a
  finished basis column for the final linear-residual check. Added multi-size
  regression coverage and a 250k/500k-variable peak-RSS comparison against
  the previous snapshot. SciPy-compatible evaluation counts are unchanged.

- Reduced transient memory in the opt-in LGMRES Jacobian-image cache paths,
  without changing their public arguments. Compact retained augmentation
  directions in place using a single candidate vector, and when both cache
  options are active, share the previously evaluated retained `J*v` images
  instead of duplicating them in the raw Arnoldi-product workspace. Keep the
  actual (pre-Gram-Schmidt) products for new Krylov columns to preserve
  stability on nearly singular systems; no extra Jacobian callbacks are
  introduced. Added a 64-variable preconditioned regression and a matrix-free
  memory benchmark, with baseline/updated peak-RSS reference measurements.
  Existing Krylov defaults are unchanged.


- Added opt-in `outer_candidate_image_reuse=.true.` to built-in Krylov
  LGMRES with analytic `jacobian_product` and Jacobian-image pruning. Form
  the candidate step image from the raw Arnoldi operator products rather than
  repeat J*v, with no additional nonlinear residual calls. The change can
  be combined with retained-image caching; numerical finite-difference
  products remain on the previous path. Added API validation, nonlinear
  regression, expanded seven-policy 16-variable stress comparison, and
  a 32/64/128-variable timing and callback-count benchmark. Defaults are
  unchanged.

- Added opt-in `outer_image_reuse=.true.` for built-in LGMRES with active
  automatic or fixed Jacobian-image pruning. Reuse prevalidated retained
  `J(x)*v` from the pruning pass in the immediately following Arnoldi solve,
  instead of recomputing the same products. Preserve identical stored
  directions for the cache, invalidate all images after the current nonlinear
  iterate, retain the existing candidate-image check, and preserve legacy
  evaluation counts when disabled (the default). Added public API tests,
  two-/three-point function-evaluation accounting, callback/preconditioner
  coverage, and a fifth policy to the changing-Jacobian stress comparison.

- Added opt-in `outer_image_prune_auto=.true.` to `root(method='krylov')`
  and `newton_krylov` for built-in LGMRES. A new automatic Jacobian-image
  cutoff depends on how much of each normalized image is removed by the
  first modified Gram-Schmidt pass, with an epsilon- and dimension-based
  roundoff floor. Two-pass span checks and original solution-space direction
  retention are preserved. This avoids requiring a user-chosen fixed
  dependence threshold and is explicitly a heuristic, not a rigorous
  Jacobian-error estimate. The default and fixed-tolerance paths remain
  unchanged. Added direct and public regression tests, option/callback
  validation, and a fourth policy to the deterministic pruning benchmark.


- Added `test/test_optimize_root_krylov_pruning_comparison.f90`, a
  deterministic 16-variable nonlinear stress comparison of unpruned LGMRES,
  solution-space pruning, and current-Jacobian-image pruning. The test uses
  dense orthogonal factors and singular values spanning three to seven orders
  of magnitude; it records success, nonlinear iterations, residual evaluations,
  analytic J*v calls, and independent final residuals. With GNU Fortran 14.2,
  the two pruning policies converged on a condition-`1e7` instance on which
  unpruned LGMRES exhausted 80 iterations. It also checks finite-difference
  evaluation accounting and does not assert that one policy always wins.
  Added an optimized-build CSV reference and reproduction notes.

- Added opt-in `outer_image_prune_tol=` to Krylov `root` and `newton_krylov`
  for built-in LGMRES. It checks the *current* Jacobian images of retained
  solution directions before each Arnoldi solve and when inserting a new
  direction, removing directions nearly dependent in image space. Existing
  solution-space vectors remain unchanged; old images are recomputed at every
  nonlinear iterate. The default `0` preserves SciPy-compatible behavior;
  positive image and solution-space pruning thresholds cannot be combined.
  Extra numerical Jacobian products are counted in `nfev` for finite
  differences. Added direct span/FIFO/scale tests, wrapper and failure tests,
  exact extra finite-difference sample counts, and optimized/debug regressions.

- Added opt-in `outer_prune_tol=` to `root(method='krylov')` and
  `newton_krylov` for built-in LGMRES with `outer_k > 0`. When positive, each
  new normalized Newton direction is projected with two modified Gram-Schmidt
  passes onto a temporary orthonormal basis of retained augmentation
  directions. Directions whose relative orthogonal component is at or below
  the tolerance are not stored. Independent new directions retain the
  established FIFO history policy; the original augmentation vectors are not
  modified. The default `0` leaves SciPy-compatible retention and evaluation
  counts unchanged. No extra Jacobian-vector or residual evaluations are
  needed. Added direct span-dependence, extreme-scale and FIFO regressions,
  public wrapper and option validation, and nonlinear LGMRES convergence tests.


- Replaced LGMRES/GMRES rank-deficient triangular least-squares fallback
  (normal equations plus eigen-decomposition of `H**T*H`) with a directly
  rank-revealing one-sided Jacobi SVD of the small rotated Hessenberg matrix.
  The new internal `scifort_optimize_krylov_rank_ls` module scales the matrix,
  rotates columns directly, and truncates singular values at a relative
  `max(m,n)*epsilon` cutoff to compute a minimum-norm solution. It avoids
  squaring the condition number and requires no BLAS/LAPACK dependency.
  Ordinary nonsingular back-substitution and all public Krylov options are
  unchanged; no additional Jacobian evaluations are needed. New standalone
  rank-deficient, near-dependent, scaled and rotated-matrix tests, and public
  LGMRES nonlinear regressions exercise this fallback.

- Added opt-in `inner_orthogonalization='selective'` (DGKS norm-loss
  criterion) and `'twice'` (two-pass modified Gram-Schmidt) to built-in GMRES
  and LGMRES through `root(method='krylov')` and `newton_krylov`.
  Accumulate the second projection into the Arnoldi Hessenberg column,
  re-evaluate its remaining norm and detect near-breakdown without extra
  Jacobian products. Preserve one-pass `'once'` as the SciPy-compatible default;
  reject unsupported Krylov/custom solvers and invalid options before any
  residual evaluation. Added a focused ill-conditioned 48-variable Arnoldi
  orthogonality regression, preconditioned LGMRES, multi-step nonlinear
  augmentation reuse, and public API tests.

- Added SciFort-specific `inner_check_residual=.true.` to Krylov `root`
  and `newton_krylov`. This opt-in check evaluates the final unpreconditioned
  inner linear residual against the requested forcing tolerance, rejecting
  inaccurate steps from all six built-in or custom inner solvers before a
  nonlinear update. It uses overflow-safe scaled arithmetic and a floating-
  point roundoff floor; analytic products do not add residual evaluations,
  while finite-difference checks count their probes in `nfev`. Separately
  fixed an underflowing zero-step check and squared-norm over/underflow in
  the inexact-Newton forcing update for extreme finite residuals. Added a
  focused regression for default compatibility, failure propagation, and
  1e-180/1e180/1e308 test cases; unchanged default paths retain SciPy parity.

- Added a SciFort-specific optional `jacobian_product=` analytic Jacobian-vector
  callback to `root(method='krylov')` and `newton_krylov`. This supplies
  `J(x)*v` directly, without finite-difference residual probes or erroneous
  cancellation near a large residual offset, and works with all six inner
  methods, custom inner solvers, preconditioning, and verified TFQMR.
  User products that fail or return nonfinite values reject the inner solve;
  zero directions bypass the callback. Explicit finite-difference controls
  cannot be combined with this option. An analytic product does not increment
  `nfev` (which counts residual evaluations). Added dedicated cancellation,
  nonlinear-refresh, dispatch, and failure-mode regressions. The default
  finite-difference API and SciPy reference counts remain unchanged.

- Added SciFort-specific `jacobian_step_mode='adaptive'` and
  `'componentwise'` to Krylov `root` and `newton_krylov`. The adaptive mode
  retains the SciPy-compatible global stencil on well-scaled products but
  uses individually scaled finite differences for disparate active coordinates
  or insufficient directional perturbations. The forced componentwise option
  supports two-point and centered three-point stencils, custom inner solvers,
  and all six built-in Krylov methods; exact extra `nfev` counts and guards
  against nonfinite function samples are covered by a new focused regression.
  The default remains `'global'`, with unchanged SciPy reference counts.

- Added opt-in `jacobian_diff_scheme='3-point'` to Krylov `root` and
  `newton_krylov`. It uses two-sided directional Jacobian differences with
  the cube-root-epsilon automatic step, detects unrepresentable perturbations,
  and avoids intermediate overflow for large finite centered samples. The
  SciPy-compatible default `'2-point'` and its evaluation counts are unchanged.
  New regression tests compare against analytic Jacobian products, exercise
  all solver dispatch paths, and cover invalid and nonfinite stencil inputs.

- Added opt-in `inner_tfqmr_right_preconditioned=.true.` for built-in TFQMR
  through `root(method='krylov')`, `newton_krylov` and a low-level
  `right_preconditioned=` switch. The matrix-free recurrence now uses `J*M`
  with `M`-applied iterate updates, consistently solving the right-preconditioned
  system when selected. The original SciPy 1.17 `M*J` default is unchanged.
  Added diagonal/nondiagonal approximate-inverse tests, original-residual
  verification, public-wrapper/invalid-option tests, and an independent SciPy
  transformed-system reference generator with work-count parity.

- Added SciFort-specific opt-in `inner_verify_residual=.true.` to
  `root(method='krylov')` and `newton_krylov` with built-in TFQMR. It checks
  the actual (unpreconditioned) `norm2(rhs - J*step)` at candidate stopping
  points and the iteration limit, continues if the estimator is optimistic,
  and fails instead of accepting an inaccurate finite iterate. Unsupported
  method/custom-solver combinations are rejected. The default remains the
  SciPy 1.17-compatible estimate-based stopping with unchanged work counts.
  Added false-convergence, delayed operator failure, and 64-variable
  robustness tests and a reproducible SciPy reference script.

- Added built-in `inner_method='bicgstab'`, `'cgs'`, `'minres'`, and `'tfqmr'`
  for matrix-free `root(method='krylov')` and `newton_krylov`. All four support
  existing inverse-preconditioner callbacks, bounded Krylov iterations,
  finite/operator-breakdown guards, and SciPy-compatible acceptance of finite
  incomplete inner iterates. MINRES assumes a symmetric Jacobian and
  positive-definite preconditioner. Added standalone matrix-free tests,
  public-root linear/nonlinear regressions, and a SciPy 1.17.0 reference
  generator with exact iteration/evaluation-count checks.

- Added `inner_m=`, `inner_m_update=`, and `inner_solver=` callbacks to
  `root(method='krylov')` and direct `newton_krylov`, including a public
  matrix-free linear-operator interface. GMRES and LGMRES support native
  right preconditioning; custom inner solvers receive the current Jacobian
  product, inverse-Jacobian action, linear tolerance and iteration budget.
  Invalid operators and unsuccessful callbacks produce a failed result.
  Added adaptive-update/invalid-output tests, SciPy 1.17.0 work-count parity
  cases and a standalone reference generator. Existing defaults are unchanged.

- Added `tol_norm=` custom vector stopping norms to `root(method='broyden1',
  'broyden2', 'anderson', 'linearmixing', 'diagbroyden', 'excitingmixing',
  'krylov')` and their direct convenience wrappers. The selected norm applies
  consistently to residual, initial relative-residual scale, iterates and
  accepted steps; the original max norm remains the default. Corrected
  premature convergence at the initial point when the residual is nonzero,
  guarded invalid custom norms, and added SciPy 1.17.0 parity regressions.


- Benchmarked all four elementwise copied versus indexed mixed-type callbacks
  with large allocatable derived-type payloads and repeatable CPU-time runs.
  Indexed callbacks now gather active argument offsets into a reusable integer
  work array rather than forming vector-subscripted temporaries on every call.
  All public signatures and numerical controllers are unchanged. Added tests
  for reordered/repeated/empty gather maps and a standalone benchmark runner.

- Added allocation-avoiding `*_vectorized_indexed_args_into` variants for
  all four packed `scipy.optimize.elementwise` helpers. Mixed-type callbacks
  receive the original `elementwise_typed_argument` objects plus per-active-
  point indices, avoiding per-evaluation allocation/copying of polymorphic
  scalar payloads. Existing copied typed/numeric entry points remain intact.
  Added parity tests for statuses, evaluations and results, independent
  trailing-axis broadcasting, noncontiguous outputs, invalid inputs, and
  ranks 14/15.

- Added mixed-type argument-list overloads to all four
  `*_vectorized_args_into` solvers. `elementwise_typed_argument_from`
  preserves integers, logicals, character strings, complex values and
  user-defined types alongside real parameters, with independent
  trailing-axis broadcasts through rank 15. Each callback receives
  boxed polymorphic scalars for precisely its active lanes. Added tests
  for numerical correctness, type preservation, compaction, invalid inputs,
  noncontiguous grids, and ranks 14/15. The existing numeric interfaces
  and numerical controllers are unchanged.


- Added independently shaped real(dp) parameter lists to all four packed
  `*_vectorized_args_into` routines. New `elementwise_argument` and
  `elementwise_argument_from(scalar_or_array)` preserve the shape of each
  argument; each broadcasts separately against the result grid. Existing
  stacked numeric args remain source-compatible via generic overloads.
  Added active-lane, stack/list parity, invalid-shape, noncontiguous and
  rank-14/15 regression coverage without duplicating numerical kernels.

- Added all four `*_vectorized_args_into` assumed-rank entry points for
  broadcasted objective parameter stacks (`args(nargs, grid...)`). Parameter-
  aware callbacks receive `(x(:), args(:,:), fx(:))` and retain the correct
  per-lane arguments even as active evaluations are compacted. The existing
  packed `(x, fx)` callbacks and rank-one through rank-four function APIs
  remain source-compatible. Added rank-zero/rank-five/rank-14/15,
  invalid-shape, noncontiguous, status/work-count and bracketing regressions.

- Added assumed-rank `find_root_vectorized_into`, `find_minimum_vectorized_into`,
  `bracket_root_vectorized_into`, and `bracket_minimum_vectorized_into` for
  Fortran 2018 grids through the highest usable ranks (14 for refinement,
  15 for bracketing), including scalar problems. Shared the flattening helper
  with trailing-axis broadcasting and refactored existing rank-two through
  rank-four function overloads into thin adapters. Results are explicitly
  initialized for malformed input or shape mismatches, with noncontiguous
  section support. Added a rank-0/rank-5/rank-14/rank-15 regression checking
  packed rank-one parity, statuses, counters, shape errors, and empty grids.

- Extended the packed vector-callback `find_root_vectorized`,
  `bracket_root_vectorized`, `find_minimum_vectorized`, and
  `bracket_minimum_vectorized` overloads to rank-three and rank-four problem
  grids. The new adapters retain Fortran element order, expand trailing-axis
  singleton broadcasts, forward optional callbacks and bounds, and reuse the
  rank-one active-lane solver. Added checked rank-three/rank-four equivalence
  tests for all four solvers, explicit bracketing endpoints, status/count
  parity, callbacks, invalid shapes, nonfinite lane inputs, and empty axes.

- Added rank-two packed vector-callback overloads for
  `find_root_vectorized`, `bracket_root_vectorized`,
  `find_minimum_vectorized`, and `bracket_minimum_vectorized`.
  Starting-point grids `(m,n)` and refinement brackets `(2,m,n)` or
  `(3,m,n)` yield `(m,n)` result arrays. Scalars, trailing vectors,
  and singleton-axis matrices broadcast optional parameters. Active-point
  packing, invalid-input isolation, callback stops, and per-element work
  counts reuse the rank-one controllers. Added a four-solver rank-two
  regression covering shape, broadcasts, callbacks, and IEEE-trap safety.

- Added `bracket_root_vectorized` and `bracket_minimum_vectorized` with
  packed `f(x(:), fx(:))` callbacks, broadcast bounds/growth factors,
  independent endpoint/expansion states, SciPy-style statuses and per-value
  `nfev` counts. Completed the vectorized rank-one bracketing/refinement
  quartet. New regressions compare with SciPy 1.17 and scalar SciFort,
  including active-lane packing, bad inputs, nonfinite outputs, and limits.


- Added `find_minimum_vectorized` using a packed `f(x(:), fx(:))`
  callback for valid starting triples and active Chandrupatla trial points.
  It retains independent statuses, per-point `nfev`, scalar stopping callbacks,
  broadcast convergence tolerances, and isolation of invalid or nonfinite
  lanes. Added a SciPy 1.17 reference generator and regression against scalar
  results under strict `-O0`, `-O2`, and `-O3` GNU Fortran builds.

- Added `find_root_vectorized` with an explicit `f(x(:), fx(:))` callback,
  batched endpoint evaluation, and dynamic active-lane packing at each
  Chandrupatla iteration. The previous `find_root` generic remains unchanged.
  Scalar state callbacks, per-element tolerance broadcasting, independent
  convergence statuses/point counts, invalid lane isolation, and zero-sized
  batches are retained. Added SciPy 1.17 status/iteration/evaluation references
  and a strict vectorized-root regression at `-O0`, `-O2`, and `-O3`.

- Extended all four `scipy.optimize.elementwise` solvers using a shared
  Fortran 2018 assumed-rank broadcasting engine. `bracket_root` and
  `bracket_minimum` now support problem grids through rank 15; `find_root`
  and `find_minimum` support rank 14 because bracket inputs include one
  endpoint dimension. Generated thin rank-specific result wrappers preserve
  grid shapes, scalar callbacks, independent iteration/function counts,
  singleton/trailing-axis broadcasting, and invalid-input handling. Added a
  rank-1-through-15 regression with high-rank broadcasting, callback, and
  empty-grid tests. Older rank-one invalid-shape tests now use truly
  incompatible arrays, since singleton rank-one inputs broadcast.


- Added rank-three `(m,n,p)` independent-problem overloads for `bracket_root`,
  `find_root`, `bracket_minimum`, and `find_minimum`. Optional arguments support
  NumPy-style trailing-axis scalar/rank-one/rank-two/rank-three broadcasting,
  preserving independent statuses and work counts. Added a SciPy 1.17 reference
  generator and rank-three regression covering shapes, singleton expansion,
  callbacks and invalid/empty grids. Existing scalar and lower-rank APIs remain
  unchanged; function evaluations still use scalar callbacks.
- Fixed a pre-existing GNU Fortran `-O3` IEEE-invalid trap in scalar `find_root`
  with a NaN tolerance: validate each tolerance for finiteness before performing
  negative-tolerance comparisons. The scalar and array solver regressions now
  pass with `-O3` and floating-point exception traps enabled.

- Added rank-two independent-problem overloads for `bracket_root`, `find_root`,
  `bracket_minimum`, and `find_minimum`. Bracketing accepts `(m,n)` initial
  points, refinement accepts `(2,m,n)` / `(3,m,n)` initial brackets, and
  results preserve `(m,n)` shape. Scalars, trailing-axis vectors, and
  singleton-axis matrices broadcast optional bounds and tolerances. SciPy
  1.17 generated reference cases and a dedicated regression check mixed
  successes, work counts, callbacks, malformed shapes, and empty grids.
  Scalar callbacks remain independent rather than vectorized.

- Added rank-one independent-problem overloads for all four
  `scipy.optimize.elementwise` bracketing/refinement helpers:
  `bracket_root`, `find_root`, `bracket_minimum`, and `find_minimum`.
  Scalar bounds, endpoints and convergence tolerances broadcast over the
  batch, or rank-one values can be supplied per problem. Each result retains
  its own status, iteration and evaluation counts; callbacks run separately
  per problem. Added independent scalar comparisons, shape/validation cases,
  and generated SciPy 1.17 rank-one numerical/status/work-count references.
  Higher-rank broadcasting and true vectorized objective callbacks remain open.

- Added scalar `find_root` (SciPy 1.17 `elementwise.find_root`), using the
  Chandrupatla hybrid interpolation/bisection iteration. Supports four
  tolerances, reversed brackets, iteration/evaluation counts, per-iteration
  stop callback (status -4), and handoff from `bracket_root`. Added 18
  generated SciPy numerical references and trap-safe Fortran regressions.
  Higher-rank broadcasting is not yet supported.

- Added scalar `find_minimum` based on SciPy 1.17
  `elementwise.find_minimum` and Chandrupatla's three-point minimization
  algorithm. Supports four optional tolerances, a configurable iteration cap,
  ordered bracket/function-value results, evaluation counts, and callbacks
  that may stop the solver with status -4. Generated SciPy references cover
  smooth/nonsmooth objectives, initial bracket ordering, convergence tolerances,
  invalid/nonfinite inputs, and iteration limits. This is not the full
  general array-broadcasting SciPy interface.

- Added scalar `bracket_minimum` with SciPy 1.17
  `elementwise.bracket_minimum`-style geometric search, optional initial
  endpoints and finite search limits, work counts, and explicit termination
  statuses. Reuses `scalar_objective_function` and can supply a three-point
  bracket to `brent`. Generated SciPy cases cover both downhill directions,
  bounded searches, plateaus, iteration limits, and exhaustion of a bound.
  A separate rank-one overload now supports independent elementwise problems.

- Added scalar `bracket_root` in `scifort_optimize_root_scalar` and the umbrella
  `scifort_optimize` module, following SciPy 1.17 `elementwise.bracket_root`.
  Supports independent left/right geometric expansion, optional interval
  limits, finite-bound initial-gap scaling, configurable `factor`/`maxiter`,
  status, function counts, and bounded/no-root failure cases. SciPy-generated
  regression cases cover both expansion directions, endpoints, one-sided
  searches, and limited growth. A rank-one overload now also exists.

- Exposed SciPy 1.17-style `bracket(f, xa, xb, grow_limit, maxiter)` as a public
  scalar-minimization helper, reusing the existing downhill expansion algorithm.
  Added `bracket_result` with the bracket points, corresponding function values,
  function-call count, iteration count, and status. Added SciPy-generated
  reference comparisons for reversed direction and limited-growth searches.
  The shared bracketing kernel now checks for NaNs before ordered comparisons
  so invalid-objective cases work with floating-point exception traps enabled.

- Fixed signed two-point **objective-gradient** differences at retained
  variable bounds in `trust-constr`: negative relative steps retain their
  directions unless the requested probe would leave the bound interval.
  The rule now holds both for initial interior-point gradient construction
  and later equality-SQP iterations. Added SciPy 1.17-generated gradient
  references and regression probes for interior, active-bound, linear-block,
  generic `minimize`, and accepted-step cases. Default forward steps and
  three-point behavior remain unchanged.

- Added opt-in `%hess_update='bfgs'` and `'sr1'` strategies for individual
  `trust-constr` nonlinear constraint blocks. Independent dense per-output
  secant matrices are weighted by the current constraint multipliers;
  the implementation supports mixed exact/updated Hessians and analytic or
  numerical block Jacobians without changing the numerical-Hessian default.
  Duplicate Jacobian work is avoided by reusing per-block iterate caches.
  These strategies count no exact constraint-Hessian evaluations, following
  SciPy's `constr_nhev` convention. Added validation, SciPy 1.17 comparison
  values, and regression coverage for equality, barrier and bounded solves.

- Added `%hess_fd_scheme='2-point'` / `'3-point'` for individual nonlinear
  `trust-constr` constraint blocks with analytic `%jac`. The dense SQP solver
  differentiates only the selected block's multiplier-weighted Jacobian,
  symmetrizes the result, obeys signed relative steps and retained bounds,
  and counts numerical Hessian probes as Jacobian calls. Missing options
  preserve the existing BFGS/mixed-Hessian behavior; explicit finite-difference
  Hessians are rejected if the Jacobian is also numerical or if `%hess` / `%hessp`
  was already supplied. Added SciPy 1.17 reference solutions and regressions.

- Aligned independent `trust-constr` nonlinear constraint objects' signed
  `%finite_diff_rel_step` vectors with SciPy 1.17's `_numdiff` behavior.
  Explicit two-point steps now use `r*sign(x)*abs(x)`, retaining negative
  directions unless variable bounds force reversal; three-point steps use
  absolute magnitudes. Zero and sub-ULP relative steps fall back to an
  automatic difference. Dense and grouped sparse Jacobians share the behavior.
  Added SciPy-generated probe/Jacobian/evaluation-count regressions and updated
  the former zero-step invalid-input test.

- Added opt-in `objective_hess_fd_scheme='2-point'` / `'3-point'` for
  `trust_constr` and `minimize(method='trust-constr')`. Numerically reconstruct
  a dense symmetric objective Hessian from the analytic objective gradient,
  using shortened/one-sided bound-safe differences and optional per-coordinate
  `objective_finite_diff_rel_step`; two-point probes preserve signed step
  directions whenever feasible. Reject conflicting exact Hessians/Hessian
  products/SR1/BFGS overrides and missing `jac`. Count gradient probes in
  `njev` but keep `nhev=0`, matching SciPy's numerical-Hessian reporting.
  SciPy 1.17 references and a new regression cover Rosenbrock, nonlinear
  equality, box barriers, fixed bounds, callbacks, and invalid input.

- Added an optional `callback=` to `trust_constr` and
  `minimize(method='trust-constr')`. The callback receives the physical
  variable vector and 1-based cumulative SQP iteration number and may request
  early termination with status 3; the full `OptimizeResult` callback argument
  is not yet implemented. Supports direct, objective `hessp`, constraint
  objects, nonlinear equalities, interior-point inequalities, and bounds.
  Added SciPy 1.17 reference verification and targeted tests for stop at the
  initial point and later iterations, convergence without stopping, status,
  shape, and invalid-method rejection.

- Added `%finite_diff_jac_sparsity(:, :)` logical masks to independent
  `trust-constr` nonlinear constraint blocks. A greedy grouping of columns
  with disjoint structural nonzero rows shares two-point or three-point
  function probes; one-sided and shortened steps honor retained variable
  bounds, all-zero columns are skipped, and exact Jacobian callbacks take
  precedence. No sparse linear algebra is implied: SQP projections are still
  dense. Added SciPy 1.17 comparison data and regression tests for colored
  and uncolored evaluation counts, overlapping row patterns, Hessian callbacks,
  fixed sparsity zeros, strict bounds, and invalid mask shapes.

- Accepted overdetermined, rank-deficient equality systems in dense
  `trust_constr` and `minimize(method='trust-constr')`, including redundant
  nonlinear equalities, constraint-object blocks, interior-point variable
  bounds, and fixed bounds. The existing Gram-eigendecomposition pseudoinverse
  is used to compute minimum-norm projected steps and multipliers rather than
  rejecting `meq > size(x0)`; inconsistent constraints still terminate as
  infeasible (status 4), not successful. Added five SciPy 1.17
  `factorization_method='SVDFactorization'` reference cases and regressions for
  exact constraint Hessians, BFGS, and SR1 objective updates.


- Added `objective_hess_update='sr1'` as an alternative to the default BFGS
  objective Hessian approximation for `trust_constr` and
  `minimize(method='trust-constr')`. The symmetric-rank-one update uses SciPy's
  first-secant scaling and small-denominator skip rule; negative curvature is
  retained. It propagates through unconstrained, nonlinear equality,
  independent constraint blocks, inequality/slack-barrier, and bounds-only
  paths. Objective SR1 updates remain separate from exact/approximated
  constraint Hessians, while explicit `hess=`/`hessp=` conflicts and invalid
  strategy names are rejected. Added a SciPy 1.17 reference generator and
  tests for all paths, analytic/numerical gradients, BFGS compatibility,
  Hessian work accounting, and option validation.

- Enabled mixed exact and missing Hessians for independent `trust-constr`
  nonlinear constraint blocks. Analytic `%hess` and `%hessp` callbacks keep
  precedence; blocks lacking both use bound-aware, symmetrized finite
  differences of their own (analytic or numerically computed) Jacobians.
  Other blocks' evaluation counts remain independent, and approximated
  Hessians do not increment `constr_nhev`. Added SciPy 1.17 reference regression
  for analytic/missing and numerical/missing combinations, generic `minimize`
  dispatch, objective `hessp=`, mixed-stencil retained bounds, and
  exact-Hessian precedence.

- Added objective `hessp(x, p, hp)` Hessian-vector products to direct
  `trust_constr` and `minimize(method='trust-constr')`. Dense reconstruction
  from coordinate-basis products reuses the existing analytic-Hessian path,
  including independent constraint blocks, multiplier-weighted constraint
  Hessians, barrier slacks, and Hessian evaluation accounting. Explicit `hess=`
  takes precedence. Added SciPy 1.17-generated references and tests for
  unconstrained, linear/nonlinear, bounds, and numerical-gradient cases.

- Added `objective_finite_diff_rel_step(:)` to `trust_constr` and generic
  `minimize(method='trust-constr')`, accepting different signed relative steps
  for each original variable. This follows SciPy's explicit relative-step
  convention, including the near-zero representable-step fallback, in both
  two-point and three-point gradients. Constraint-object and slack-barrier
  paths preserve the vector and existing gradient/function work accounting;
  strict retained bounds still clip probes. Added SciPy 1.17 comparison tests
  for mixed steps, signed two-point probes, one-sided bounded probes, zero
  coordinates, exact-gradient precedence, validation, and a complete solve.

- Added opt-in `objective_fd_scheme='3-point'` to `trust_constr` and
  `minimize(method='trust-constr')`, including independent constraint objects.
  Objective gradients now use centered second-order differences in the interior
  and second-order one-sided/shortened differences at retained variable bounds.
  The explicit `eps=` step is honored, while the default three-point step is
  `epsilon(1.0_dp)**(1.0_dp/3.0_dp)` independently of two-point constraint
  Jacobian steps. The existing two-point default remains unchanged. The
  numerical derivative code counts all objective probes, including initial
  barrier stages; positive-slack domains are protected from negative probes.
  A finite slack floor also prevents an overflow in the analytic barrier
  Hessian when a retained bound is active at the starting point. Generated
  SciPy 1.17 references and new tests cover work-count parity for unconstrained,
  linear, and nonlinear equality cases, strict bounds, one-sided stencils,
  explicit steps, analytic-gradient precedence, and invalid options.

- Enabled `trust_constr` and `minimize(method='trust-constr')` without any
  constraints or finite variable bounds, including an empty
  `constraint_objects=` array. The zero-row equality-SQP path now safely
  handles empty projections, multipliers, and constraint diagnostics without
  evaluating nonexistent callbacks. Added generated SciPy 1.17 reference tests
  for analytic/numerical gradients, an analytic Hessian, Rosenbrock, BFGS,
  stationary starts, and iteration limits.

- Added SciPy-style three-point numerical Jacobians for independently defined
  `trust-constr` nonlinear constraints. Each block now accepts
  `%fd_scheme = '2-point'` (default) or `'3-point'`, plus optional positive
  `%finite_diff_rel_step(:)` values per optimization variable. Central and
  second-order one-sided stencils use the actual representable displacements,
  including shortened probes at retained variable bounds. SciPy 1.17 reference
  tests cover the exact evaluation counts with mixed analytic/numerical blocks,
  nonlinear constraint derivatives, per-coordinate steps, and invalid inputs.
  Barrier trials outside `bounds_keep_feasible` bounds now return an infinite
  objective and do not invoke user constraints, including for an initial point
  exactly on a strict bound.

- Added first-class multiple nonlinear constraint blocks to `trust_constr` and
  `minimize` through `constraint_objects=` and the public `nonlinear_constraint`
  type. Each block has independently bound-valued components, feasibility
  retention, optional analytic Jacobians and weighted Hessians, and per-block
  evaluation counts; mixed analytic/numerical Jacobians are supported. The
  existing combined callback interface remains unchanged, and the generic
  minimizer now selects `trust-constr` when constraint objects are supplied.
  The analytic/Jacobian accounting matches SciPy 1.17 for separate equality and
  inequality objects, including finite-difference function probes and zero
  explicit Jacobian calls for numerical blocks. New tests cover nonlinear
  Hessian multiplier routing, bounds-block ordering, strict finite-difference
  clipping, invalid shapes, infeasible retained starts, and legacy conflicts.
  Partial provision of analytic constraint Hessians remains unsupported and is
  explicitly rejected rather than silently dropping curvature.

- Added strict finite-difference clipping for `trust-constr` variable bounds marked
  `bounds_keep_feasible=.true.`. Numerical objective-gradient and nonlinear-
  constraint Jacobian probes now reverse or shorten their one-sided step when the
  default forward probe would cross a retained lower/upper bound, following the
  `strict_bounds`/bounded numerical-differentiation behavior used by SciPy 1.17.
  The clipping is threaded through barrier-stage cached derivatives and subsequent
  equality-SQP numerical gradient evaluations. A new SciPy-backed regression starts
  immediately below a retained upper bound with both objective and constraint
  derivatives numerical and verifies that neither callback is evaluated above the
  bound. Multiple first-class constraint objects remain the outstanding
  `trust-constr` roadmap item.

- Added SciPy-style per-component feasibility retention to the current
  `trust-constr` APIs. `constraint_keep_feasible=` applies to raw components used
  with `constraint_lower=`/`constraint_upper=`, while `bounds_keep_feasible=`
  applies per optimization variable. Equality components ignore the flag; an
  infeasible initial point for any retained inequality is rejected. Interior-point
  subproblems carry the canonical keep mask into the slack system and tie each
  retained slack to its current constraint gap before trial evaluation, matching
  SciPy's barrier feasibility mechanism. Generated SciPy 1.17 regressions match
  the retained nonlinear-lower-bound case at 9/6/6/6/5
  iteration/objective/gradient/Hessian/CG counts and a retained variable lower
  bound at 8/5/5/5/4 (with `barrier_tol=1e-3`), including generic `minimize`
  forwarding and infeasible-start validation. Strict finite-difference clipping
  is covered by the entry above; first-class multiple constraint objects remain
  roadmap work.

- Added SciPy-style per-component lower/upper bounds to the existing combined
  `trust-constr` nonlinear-constraint callback through `constraint_lower=` and
  `constraint_upper=`. Components are canonicalized as equalities, lower-only,
  upper-only, or two-sided intervals using SciPy 1.17 ordering/sign conventions,
  while final `constr`/`constr_jac` diagnostics remain in the callback's original
  value/Jacobian representation. Weighted analytic constraint Hessians map the
  canonical SQP multipliers back to the raw callback multiplier vector before
  invoking `constraint_hess(x, v)`. The legacy `meq`/`mineq` interface remains
  unchanged and is rejected only when mixed ambiguously with the new bound-vector
  form. Generated SciPy regressions cover all four component types in one callback,
  exact 12/8/7 iteration/objective/CG accounting, direct and `minimize` dispatch,
  equality-only bound-vector dispatch, and upper-bounded analytic Hessian mapping.

- Added rich SciPy-style `trust-constr` result diagnostics: final constraint
  values, row-stacked constraint Jacobians, per-block constraint
  function/Jacobian/Hessian evaluation counts, and the internal method names
  `equality_constrained_sqp` / `tr_interior_point`. SciFort's current combined
  user-constraint callback is one reporting block and any variable bounds form
  a second block after it, preserving SciPy's bounds-after-user-constraints
  order. Barrier stages now carry cached constraint values/Jacobians across
  unchanged restart points and outer stopping checks consume that cached state;
  analytic constraint Hessians also perform SciPy's initial zero-multiplier
  setup call. The equality Maratos constraint counts now match SciPy 1.17 at
  9/7/8 and the barrier Maratos counts at 15/12/19 for function/Jacobian/Hessian
  calls, with a mixed user-constraint-plus-bounds regression covering block
  ordering and zero-cost bounds diagnostics.

- Extended `trust-constr` result parity with SciPy 1.17. The result now reports
  the full Lagrangian gradient, cumulative projected-CG iteration count, and final
  projected-CG stop condition. The dense tangential subproblem now uses the same
  SciPy-style projected-CG loop for pure equality and slack-barrier problems,
  including the `n-m` iteration cap and stop codes; equality and barrier Maratos
  regressions match SciPy at 7 and 17 cumulative CG iterations respectively.
  Apparent `gtol`/`xtol` success is also converted to status `4` with
  `success=.false.` when constraint violation remains above `gtol`, matching
  SciPy's post-solve success semantics.

- Closed the representative `trust-constr` interior-point work-count gap. Barrier
  stages now reuse the already-computed objective, gradient, and analytic objective
  Hessian at an unchanged restart point, matching SciPy's cached barrier-state
  handoff instead of reevaluating user callbacks. Equality SQP iteration accounting
  also counts each stage's initial stopping state, including the `maxiter` boundary.
  The nonlinear-inequality Maratos regression now matches SciPy 1.17 exactly at
  19 iterations, 15 objective evaluations, 12 gradient evaluations, and 12 Hessian
  evaluations. `multivariate_minimize_result` also reports SciPy-style
  `constr_violation`, `optimality`, `barrier_parameter`, and `barrier_tolerance`
  diagnostics; generated references cover the new counts and barrier diagnostics.

- Completed the next `trust-constr` interior-point fidelity slice. Barrier SQP now
  uses SciPy-style scaled slack coordinates, the `-0.995` slack trust-box boundary
  in normal/tangential subproblems, and a fivefold trust-radius carry between
  barrier stages. Inner SQP multipliers are retained for outer original-problem
  optimality checks, and missing nonlinear-constraint Hessians use a separate
  SciPy-style BFGS approximation instead of being treated as zero when the
  objective Hessian is analytic. Pure-bound barriers avoid premature `gtol`
  termination while the barrier parameter is still coarse, and bound-only SOC is
  now enabled. The nonlinear-inequality Maratos regression remains at the SciPy
  1.17 solution while dropping from 95 to 22 objective evaluations (SciPy: 15).

- Extended second-order correction (SOC) into nonlinear-inequality `trust-constr`
  barrier subproblems. Barrier SQP now clips each SOC correction against SciPy's
  `0.995` fractional slack boundary before evaluating the corrected point, so a
  correction cannot consume the remaining positive slack. A new active-inequality
  Maratos regression matches the SciPy 1.17 solution while reducing SciFort's
  objective evaluations from 172 in the earlier checkpoint to 95 at that stage.
  The later scaled-slack work above reduces this further to 22 versus SciPy's 15.

- Added SciPy-style second-order correction (SOC) to the dense equality-constrained
  `trust-constr` SQP kernel. Rejected steps with a small normal component now use
  the current Jacobian pseudoinverse to correct nonlinear constraint curvature
  before merit-function acceptance. A classic Maratos-effect regression reaches
  the solution in 9 objective evaluations, matching SciPy 1.17 and improving on
  the previous SciFort checkpoint's 14 evaluations. Nonlinear-inequality barrier
  SOC is covered by the entries above; bound-only barrier SOC is now enabled by
  the scaled-slack trust-box implementation.

- Extended analytic constraint-Hessian support through the `trust-constr`
  inequality/slack interior-point path. `constraint_hess(x, v, hmat)` continues
  to follow SciPy's `NonlinearConstraint.hess(x, v)` convention; barrier-problem
  multipliers are mapped back to the original equality/inequality constraint
  vector before evaluating the weighted Hessian, while bound rows contribute
  zero curvature. Added SciPy 1.17 nonlinear-inequality reference coverage plus
  direct and generic `minimize(method='trust-constr')` callback checks. An
  analytic objective Hessian is still required when `constraint_hess` is used.

- Added equality-only analytic constraint-Hessian support to `trust-constr`.
  The new `constraint_hess(x, v, hmat)` callback follows SciPy's
  `NonlinearConstraint.hess(x, v)` convention and is added to the analytic
  objective Hessian to form the Lagrangian Hessian. Direct and generic
  `minimize(method='trust-constr')` paths are covered by SciPy 1.17 reference
  regressions.

- Added SciPy's `trust-constr` spelling `initial_tr_radius` as an alias for the
  existing SciFort `initial_trust_radius` control. Direct `trust_constr` and
  generic `minimize(method='trust-constr')` both accept the SciPy spelling;
  supplying both spellings is rejected as ambiguous.

- Added SciPy-style `trust-constr` `initial_constr_penalty` handling. The dense
  SQP merit function now starts from the requested penalty, interior-point
  barrier subproblems receive the same initial penalty, and
  `multivariate_minimize_result%constr_penalty` reports the final value. Added
  SciPy 1.17 reference and direct/generic/interior-point regression coverage.

- Added SciPy-style `trust-constr` barrier controls: `barrier_tol`,
  `initial_barrier_parameter`, and `initial_barrier_tolerance`. Generic `tol=`
  now also sets the trust-constr barrier tolerance, matching SciPy's minimize
  dispatch. Added a generated SciPy 1.17 regression with nondefault barrier
  settings and direct/generic-dispatch coverage.

- Extended `trust-constr` from dense equality SQP to a dense trust-region
  interior-point path for inequalities and finite bounds. Nonfixed constraints use
  positive slack variables with SciPy-style initial barrier parameter/tolerance
  `0.1` and `0.2` decay, fixed bounds become equalities, and each barrier problem
  reuses the Byrd-Omojokun normal/tangential SQP kernel. The path supports mixed
  equalities/inequalities, one/two-sided and fixed bounds, analytic or numerical
  objective/constraint derivatives, analytic objective Hessians or BFGS, global
  KKT-style `gtol` stopping, `xtol`, and `maxcv`. Generated SciPy 1.17 references
  now cover active inequalities, bounds, mixed constraints, and fixed bounds in
  addition to the earlier equality cases.

- Added `minimize(method='cobyla')` and direct `cobyla` with a native
  PRIMA-derived trust-region/simplex kernel, SciFort equality/inequality callback
  conversion, finite and fixed bounds, SciPy-style bound-only start projection,
  `rhobeg`, `tol`, `catol`, `f_target`, evaluation limits, and `maxcv` result
  reporting. Generated SciPy 1.17 regressions cover nonlinear inequalities,
  equality conversion, mixed/fixed and one-sided bounds, incompatible constraints,
  objective-target stopping, generic dispatch, and evaluation-limit termination.

- Added `minimize(method='tnc')` and direct `tnc` with scaled finite bounds, active-set projected-gradient stopping, BFGS-preconditioned truncated-CG finite-difference Hessian-vector products for analytic gradients, and a bounds-aware numerical-gradient/BFGS fallback. Generated SciPy 1.17 regressions cover analytic and numerical gradients, Rosenbrock active-bound cases, generic dispatch, and evaluation-limit termination.

- Added SciPy 1.17-style bounds to Nelder-Mead and Powell minimization. Nelder-Mead clips trial vertices and reflects an upper-bound-crossing initial simplex back into the feasible box; Powell now uses finite bounded scalar searches or SciPy's tangent transform for one-sided line bounds, including bounded extrapolation and out-of-bounds status handling. Generated regressions cover box bounds, clipped starts, one-sided bounds, and exact iteration/function-evaluation counts.

- Added `minimize(method='trust-krylov')` and direct `trust_krylov` with a native
  GLTR/Lanczos trust-region subproblem solver, dense-Hessian or Hessian-product
  support, SciPy/TRLIB inexact and exact forcing sequences, hot-started Krylov
  state across rejected steps, and generated SciPy 1.17 solution/work-count
  regressions including a five-variable Rosenbrock case.

- Added `minimize(method='trust-exact')` and direct `trust_exact` with SciPy 1.17-style
  nearly-exact trust-region subproblems, analytic gradient/Hessian support,
  indefinite-Hessian handling, `subproblem_maxiter`, and generated reference/work-count
  regressions.

- Added `minimize(method='dogleg')` and direct `dogleg` with SciPy 1.17-style Cauchy/Newton dog-leg trust-region steps, positive-definite Hessian checks, trust-radius controls, and generated reference/work-count regressions.

### Added

- Added SciPy 1.17-style `minimize(method='trust-ncg')` and direct `trust_ncg`.
  The solver accepts an analytic gradient with either a dense Hessian or
  Hessian-vector product, solves trust-region subproblems with Steihaug CG,
  handles nonpositive curvature by selecting the best trust-boundary
  intersection, and exposes SciPy's `initial_trust_radius`, `max_trust_radius`,
  `eta`, `gtol`, and `maxiter` controls. Generated regressions cover Rosenbrock
  dense/Hessian-product paths, custom trust controls, iteration limits, convex
  quadratic convergence, indefinite-Hessian boundary steps, and SciPy-compatible
  objective/gradient/Hessian work counts.

- Added the core nonlinear least-squares family to `scifort_optimize`:
  `least_squares` now supports dense SciPy-style TRF (default), dogbox, and
  MINPACK LM; legacy `leastsq` and `curve_fit` are also public. TRF/dogbox add
  finite bounds, strict-feasible TRF steps, Coleman-Li scaling, dense exact
  trust-region/dogleg solves, bounds-aware two-point/three-point and complex-step
  finite differences, analytic Jacobians, fixed or Jacobian-column variable scales,
  evaluation-limit semantics, and the built-in
  `soft_l1`, `huber`, `cauchy`, and `arctan` robust losses with `f_scale`.
  `curve_fit` now defaults to TRF for bounded problems and propagates robust
  trust-region losses/covariance scaling. Added deterministic SciPy 1.17
  references and regressions for bounded/unbounded TRF/dogbox, robust outlier
  fits, curve-fit covariance, work counts, and termination limits.
  TRF also supports `tr_solver='lsmr'` with optional `tr_regularize=.false.`;
  analytic Jacobians may now be supplied matrix-free through paired
  `jac_matvec=`/`jac_rmatvec=` callbacks, avoiding dense Jacobian storage and
  automatically selecting LSMR. Numerical Jacobians also accept a logical `jac_sparsity` structural mask for grouped two-point, three-point, and complex-step differences; TRF and dogbox automatically select LSMR for this path. Dogbox also accepts explicit `tr_solver='lsmr'` for dense analytic or numerical Jacobians. Dense analytic and numerical Jacobians remain supported.

- Added bounded/constrained multivariate minimization to `scifort_optimize`.
  `minimize` now dispatches to L-BFGS-B when bounds are supplied and to SLSQP
  when constraints are supplied; direct `lbfgsb` and `slsqp` entry points are
  also public. L-BFGS-B uses the pinned Jacob Williams `lbfgsb` 1.0.0 FPM
  backend with SciPy-style clipping, projected-gradient/function tolerances,
  correction history, bounds-aware finite differences, and work-limit
  accounting. SLSQP uses Jacob Williams `slsqp` 1.6.0 with equality-first
  combined constraints, SciPy-style nonnegative inequality convention, variable
  bounds, analytic or bounds-aware numerical objective/constraint Jacobians,
  and status/work-count mapping. Added generated SciPy 1.17 references and a
  dedicated dependency-backed regression program for bounded, equality,
  inequality, mixed-Jacobian, default-dispatch, and iteration-limit cases.

- Added the first unconstrained multivariate-minimization slice to
  `scifort_optimize`: SciPy 1.17-style `minimize`, `nelder_mead`, `powell`,
  `bfgs`, `cg`, and `newton_cg`. Nelder-Mead includes adaptive coefficients and
  custom initial simplices; Powell supports custom direction sets; BFGS and
  Polak-Ribiere+ conjugate gradient accept analytic or forward-difference
  gradients; Newton-CG accepts dense Hessians, Hessian-vector products, or
  finite-difference Hessian-vector products. Gradient-based methods use
  SciPy-compatible strong-Wolfe search with the second Wolfe-search fallback.
  `multivariate_minimize_result` exposes SciPy-style solution/status/work counts,
  including `nhev`, plus method-specific gradient, inverse-Hessian, final-simplex,
  and direction diagnostics. Added generated SciPy 1.17 Rosenbrock/quadratic
  references and strict regression coverage for convergence, evaluation
  accounting, generic dispatch, adaptive Nelder-Mead, Hessian/Hessian-product
  Newton-CG, and iteration/evaluation limits.

- Completed the core SciPy 1.17 multidimensional `root` method family with
  `root(method='krylov')` and direct `newton_krylov`. The matrix-free solver
  implements SciPy-compatible finite-difference Jacobian products, nonlinear
  forcing terms, LGMRES recycling (`outer_k`) and GMRES inner solves, with
  `rdiff`, inner-dimension, tolerance, line-search, and iteration-limit controls.
  Generated regressions cover default LGMRES, GMRES, disabled recycling, custom
  `rdiff`, stricter `fatol`, and failure at the iteration limit.

- Added SciPy-compatible strong-Wolfe line search for the inexact-Newton and
  Krylov root methods. `src/optimize/wolfe.f90` adapts SciPy 1.17's
  DCSRCH/DCSTEP translation of the MINPACK-2 safeguarded strong-Wolfe search;
  Broyden and Krylov regressions match SciPy's iterates and residual-evaluation
  counts. Broyden 1/2 now also accept `max_rank` and `reduction_method` with
  `restart`, `simple`, and SVD rank reduction, including convergence and
  nonconvergence reference cases.

- Added `root(method='df-sane')`, adapting SciPy 1.17's derivative-free
  spectral residual iteration and both Cruz and Cheng nonmonotone line searches.
  SciFort exposes the relative/absolute residual tolerances, evaluation cap,
  spectral safeguards, initial spectral coefficient, and nonmonotone history
  window (`sane_m=` in Fortran because `m=` is already the LM residual count).
  Generated SciPy 1.17 regressions cover both line searches, generic `tol`, and
  function-evaluation-limit termination.

- Expanded multidimensional `scifort_optimize.root` with SciPy 1.17-style
  `broyden1`, `broyden2`, `anderson`, `diagbroyden`, `linearmixing`, and
  `excitingmixing` methods plus matching direct convenience functions. The dense
  inexact-Newton controller reproduces SciPy's default Jacobian autoscaling,
  Armijo line search, convergence tests, iteration counts, and residual-evaluation
  counts; `line_search='none'`, generic `tol`, and iteration-limit behavior are
  regression-tested against generated SciPy 1.17 references. The Broyden update
  is also cross-checked against Jacob Williams' nlesolver-fortran 1.1.0 under
  BSD-3-Clause attribution.

- Added the first multidimensional nonlinear-equation slice to `scifort_optimize`: SciPy-style `root` methods `hybr` and `lm`, plus `fsolve`, backed by the community-maintained `fortran-lang/minpack` FPM dependency. The adapters support finite-difference or analytic Jacobians, MINPACK variable scaling, banded finite differences for `hybr`, overdetermined Levenberg-Marquardt systems, SciPy-compatible default tolerances and evaluation accounting, and diagnostic QR work arrays. Added SciPy 1.17 reference generation and regressions for square, banded, analytic-Jacobian, and overdetermined systems.

- Added SciPy 1.17-style scalar minimization to `scifort_optimize`:
  `minimize_scalar`, `brent`, `bounded`, and `golden`, with two/three-point
  bracketing, automatic bounded dispatch when bounds are supplied, explicit
  status/result reporting, and SciPy-compatible iteration/function-call counts.
  The bounded Brent kernel is adapted from Jacob Williams' fmin 1.1.2 under
  BSD-3-Clause; unbounded Brent, golden-section, bracketing, and dispatch
  semantics are adapted from SciPy 1.17.0. Added generated SciPy references and
  runtime-check regressions for convergence, nonconvergence, invalid input, and
  nonfinite objectives.

- Added the first `scifort_optimize` slice with SciPy 1.17-style scalar root
  finding: `root_scalar`, `bisect`, `brentq`, `brenth`, `ridder`, `toms748`,
  and `newton` (secant/Newton/Halley). `root_result` records convergence,
  status/flag, method, iterations, and function-call counts. Bisection, Brent,
  and Ridders kernels are adapted from roots-fortran 1.0.0 by Jacob Williams
  under BSD-3-Clause; root-solver dispatch, open-method behavior, and TOMS748
  organization are adapted from SciPy 1.17.0 under BSD-3-Clause. Added SciPy
  reference generation and strict IEEE-trap/runtime-check regressions.

- Added the `scifort_differentiate` subpackage with SciPy 1.17-style adaptive
  `derivative`, `jacobian`, and `hessian` APIs. The nested finite-difference
  stencil controller reports derivative/error estimates, convergence status,
  iteration counts, and function-evaluation counts; Jacobian/Hessian interfaces
  support scalar/vector Fortran callbacks. Added SciFort extension complex-step
  derivative/gradient/Jacobian helpers for analytic functions, with the scalar
  formula adapted from ForDiff 0.1.1 under BSD-3-Clause. NumDiff 1.7.0 was
  consulted as an architectural reference without copying implementation code.

- Added first-class continuous-distribution infrastructure with native Fortran
  `make_distribution` and `truncate`. Pure scalar PDF/CDF/PPF callbacks can now be
  frozen into `continuous_distribution` objects with optional stable log/tail
  callbacks, support metadata, inverse tails, medians, and explicit-state sampling;
  `truncate` provides normalized support clipping with SciPy 1.17-compatible
  density/CDF/survival/inverse-transform semantics. The Fortran adapter deliberately
  freezes parameters at construction rather than dynamically wrapping Python-style
  `rv_continuous` classes.

### Fixed

- Fixed `odeint` and `odeint_tfirst` Jacobian adapters for compilers that enforce
  `-Werror=implicit-interface`: optional Jacobian callbacks are now captured in
  explicitly typed procedure pointers before the internal adapter invokes them.
  This preserves the public callback API while keeping full compile-time interface
  checking.

- Made the inverse-gamma and generalized-gamma positive-affine quantile guards
  safe under IEEE overflow traps when `scale < 1`. The previous guard could overflow
  while computing its own threshold even when the requested quantile was finite.

- Corrected complex `log_softmax` max-shift branch handling to match SciPy when
  subtracting the complex maximum changes the logarithm branch by an exact multiple
  of `2*pi*i`.

- Fixed the Radau sparse-Jacobian rejection path so the second embedded-error
  solve uses the sparse real factorization instead of the intentionally
  unallocated dense LU workspace.

- Use FPM monolithic library mode for SciFort so dependency objects are bundled
  into the root archive. This avoids unresolved `quadpack_double` symbols seen
  with per-package static archives on Windows/MinGW while preserving the
  external Quadpack dependency and public SciFort API.

- Restored FPM dependency resolution for Quadpack. The released Quadpack 2.1.x
  manifests use spaced array-table headers (`[[ test ]]`) that some current FPM
  parsers reject before compilation. SciFort now follows Quadpack `master`, whose
  manifest uses the parser-compatible `[[test]]` spelling while retaining the
  required public double-precision QUADPACK interfaces.

- Guarded the extreme negative tail of `normal_logcdf` before forming `z*z`.
  Values whose logarithmic probability is below binary64 range now return
  negative infinity without raising an overflow exception in trapping builds.

- Hardened large-integer-order real `J_v` evaluation when `x` is small relative
  to the order by using the existing convergent power-series kernel instead of
  the platform `BESSEL_JN` intrinsic. This prevents a Windows gfortran/libm
  internal overflow from turning `jvp(0,1,1000)` into NaN.
- Reworked nonnegative real `factorial` evaluation as a guarded Gamma recurrence
  from an anchor in `[1,2)`, avoiding near-overflow `exp(loggamma)` error
  amplification and restoring portable accuracy at `factorial(170.5)`.

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

- Added SciPy-style `multiscale_graphcorr` (MGC) independence/two-sample testing.
  The native Fortran API accepts vector or observation-by-feature matrix inputs,
  Euclidean or callback-defined distances, and precomputed square distance matrices;
  it returns the full local-correlation map, one-based optimal scale, statistic,
  permutation null distribution, and plus-one p-value. Randomized permutations use
  explicit `rng_state`; SciPy's `workers` parallelism is intentionally omitted.
  Added SciPy 1.17.0 full-map/statistic/scale references for raw, custom-distance,
  tied-precomputed, and automatic two-sample cases plus deterministic replay tests.

- Added SciPy-style generalized Monte Carlo `goodness_of_fit` testing for
  normal, exponential, and uniform location-scale families. The implementation
  supports Anderson-Darling, Kolmogorov-Smirnov, Cramer-von Mises, and Filliben
  statistics, uses explicit `rng_state`, returns the fitted null parameters and
  Monte Carlo null distribution, applies the plus-one p-value correction, and
  preserves SciPy's distinction between parameters known a priori and parameters
  supplied only to define the null distribution. Added SciPy 1.17.0 deterministic
  statistic/fitted-parameter references and deterministic replay tests.

- Added SciPy-style `sobol_indices` global sensitivity analysis with the Saltelli
  2010 first- and total-order estimators. The native Fortran API supports both
  precomputed `f_A`/`f_B`/`f_AB` evaluations and a vectorized model callback with
  independent inverse-CDF distribution callbacks. Model sampling reuses SciFort's
  64-bit scrambled Sobol engine with explicit `rng_state`, preserves SciPy's
  mean-normalization and constant-output-to-zero behavior, and enforces power-of-two
  sample counts. Added SciPy 1.17.0 scalar/multi-output reference regressions plus
  deterministic additive-model replay tests.

- Added SciPy-style `mvsdist` and `bayes_mvs` Bayesian inference for a normal
  sample with unknown mean and variance. `mvsdist` returns frozen posterior objects
  for the mean, variance, and standard deviation with density/CDF/quantile/moment/
  interval methods; `bayes_mvs` returns posterior-mean centers and central intervals.
  The implementation preserves SciPy 1.17's Student-t, inverse-gamma, and
  generalized-gamma posteriors for n <= 1000 and its Gaussian approximations for
  larger samples, with assumed-rank flattening and explicit Fortran status codes.
  Added SciPy 1.17.0 regressions for the documented example, small-sample infinite
  moments, the n > 1000 switch, matrix flattening, constant/NaN samples, and invalid
  interval probabilities.

- Added SciPy-style `obrientransform` for O'Brien's homogeneity-of-variance
  transformation. A real vector overload returns one transformed sample and a
  `sample_group` overload preserves unequal group lengths for direct use with
  `f_oneway`. Added SciPy 1.17.0 hand-check, unequal-group, constant-sample,
  short-sample, and transformed-ANOVA regressions.

- Added SciPy-style `ttest_ind_from_stats` as an elemental native-Fortran
  summary-statistics t-test. The interface accepts sample means, corrected sample
  standard deviations, and integer observation counts; supports pooled and Welch
  tests plus `two-sided`/`less`/`greater` alternatives; treats negative standard
  deviations by magnitude; preserves SciPy's one-observation pooled/Welch edge
  behavior; and returns SciFort's existing `ttest_result` with the computed degrees
  of freedom. Added SciPy 1.17.0 regression values for scalar, one-sided,
  negative-standard-deviation, one-observation, and elemental vector cases.

- Added the SciPy 1.17 `binned_statistic` family for one-, two-, and
  arbitrary-dimensional samples. Built-in mean, standard deviation, median, count,
  sum, minimum, and maximum statistics support equal-width or explicit edges,
  SciPy-style outlier/rightmost-edge bin numbering, custom pure callbacks, 2-D
  expanded bin numbers, and flattened C-order arbitrary-dimensional results with an
  explicit shape vector. Added regenerated SciPy 1.17.0 references for built-in
  statistics, explicit/outlier bins, callbacks, 2-D numbering, and 3-D
  equal/unequal-edge cases.

- Added SciPy 1.17-style `directional_stats` for vector observations. Rank-2
  inputs return a unit mean direction and scalar mean resultant length, while rank-3
  inputs provide batched results for any one-based observation axis. The default
  normalizes observations before averaging and `normalize=.false.` preserves their
  magnitudes. Added SciPy 1.17.0 references for both matrix axes, all rank-3 axes,
  raw-vector averaging, zero-resultant behavior, and invalid-axis handling.

- Added SciPy 1.17-style `spearmanrho` for paired vectors, including all three
  alternatives, asymptotic rank-correlation inference, exact permutation p-values
  when the full pairing null fits within `n_resamples`, randomized pairing tests, and
  continuous-null Monte Carlo inference with explicit `rng_state`. Added SciPy 1.17.0
  references for the documented tied-rank example and all exact alternatives, plus
  deterministic randomized replay and invalid-input regression coverage.

- Added SciPy-style Somers D ordinal association via `somersd`, with paired-ranking
  and direct contingency-table overloads, `two-sided`/`less`/`greater` alternatives,
  SciPy-compatible asymmetric `D(Y|X)` normalization, returned contingency tables,
  and trap-safe degenerate/perfect-association handling. Added SciPy 1.17.0 references
  for the canonical documentation table, transposed asymmetry, and tied sample ranks.

- Added SciPy-style Chatterjee xi correlation and independence testing via `chatterjeexi`.
  The scalar/vector interface implements the tie-aware asymptotic variance, the
  `y_continuous` shortcut, exact pairing-permutation inference when the full null
  distribution fits within `n_resamples`, and randomized Fisher-Yates pairings with
  explicit `rng_state` otherwise. Added SciPy 1.17.0 references for tied responses,
  strongly negative xi, continuous-response normalization, exact permutation p-values,
  deterministic randomized replay, and invalid-input/trap-safe edge cases.

- Added SciPy-style Baumgartner-Weiss-Schindler two-sample testing via `bws_test`,
  with `two-sided`, `less`, and `greater` alternatives, exact permutation
  enumeration for small sample splits, randomized permutation inference for larger
  problems, tied-rank support, and explicit `rng_state` for randomized draws. Added
  SciPy 1.17.0 exact references for canonical and tied cases plus trap-safe regression
  coverage.

- Added SciPy-style right-censored survival comparison with public `censored_data`,
  `ecdf(censored_data)`, and `logrank`. The logrank result preserves the signed
  standard-normal statistic, supports `two-sided`, `less`, and `greater` alternatives,
  and accepts uncensored vectors and right-censored objects in any pairing. Added the
  canonical malignant-glioma example, randomized fixed-reference cases, IEEE-trap
  regressions, and child-process validation/error checks.

- Added SciPy-style empirical-distribution and quantile-inference helpers: `ecdf`
  for uncensored and right-censored/Kaplan-Meier samples, Greenwood and log-log
  confidence intervals, exact `quantile_test` with type-bound confidence intervals,
  and callback-backed `order_statistic` PDF/CDF/SF/PPF/ISF evaluation via the beta
  order-statistic law. Added SciPy 1.17.0 examples, IEEE-trap regressions, and
  child-process validation/error checks.

- Added SciPy-style `wasserstein_distance_nd` for vector and multidimensional
  empirical distributions. Rank-2 inputs use rows as observations, support optional
  independently normalized weights, and solve the balanced Euclidean optimal-transport
  problem with a project-local successive-shortest-augmenting-path min-cost-flow
  solver. Added fixed SciPy 1.17.0 examples, randomized weighted/unweighted validation,
  IEEE-trap regressions, and child-process validation/error checks.

- Added SciPy-style one-dimensional `wasserstein_distance` and `energy_distance`
  for empirical distributions, including optional independently normalized sample
  weights, repeated support values, SciPy-compatible validation, and trap-safe
  NaN/infinity behavior.

- Added SciPy-style probability-plot diagnostics `probplot`, `ppcc_plot`,
  `ppcc_max`, `boxcox_normplot`, and `yeojohnson_normplot`. The implementation uses
  Filliben order-statistic medians, returns native Fortran output vectors and fit
  summaries, supports the common SciFort named PPF families plus custom PPF
  callbacks, and reproduces SciPy's two-point PPCC bracket-expansion semantics.
  Added fixed SciPy 1.17.0 references, randomized optimizer validation, IEEE-trap
  regressions, and child-process argument/error checks.

- Added SciPy-style Box-Cox and Yeo-Johnson power-transform statistics. `boxcox` and
  `yeojohnson` are elemental transforms; `boxcox_llf` / `yeojohnson_llf` support
  one-based axis reductions through rank 15, `*_all`, `*_keepdims`, and all three
  NaN policies. Added Box-Cox PPCC/MLE and Yeo-Johnson MLE parameter searches plus
  explicit vector fit helpers. Added SciPy 1.17.0 reference, rank-15, optimizer,
  error-path, and floating-point-trap regressions.

- Added SciPy-style sample `lmoment`. Scalar and vector orders are supported along
  one-based axes through rank 15, calls without `order` return the first four
  L-moments, `lmoment_all` represents `axis=None`, and `lmoment_keepdims` preserves
  reduced axes while matching SciPy 1.17's default-order versus explicit-vector
  layout. Added `sorted`, `standardize`, and all three NaN policies, rank-15 and
  floating-point-trap regressions, plus child-process validation of order and NaN
  error paths. Vector/default keepdims is available through rank 14 because rank-15
  input would require a rank-16 Fortran result.

- Added SciPy-style `expectile` and `sigmaclip`. `expectile` accepts scalar through
  rank-15 inputs, levels in `[0,1]`, and NumPy-style trailing/singleton-broadcast
  case weights; finite nonnegative weights use an exact piecewise-linear solution,
  with a secant fallback for SciPy-compatible negative-weight cases. `sigmaclip`
  flattens rank-15 input, iteratively removes observations outside population-
  standard-deviation thresholds, preserves retained observation order, and returns
  the final lower/upper limits. Added SciPy 1.17.0 weighted, broadcast, rank-15,
  edge-case, and floating-point-trap regressions.

- Added SciPy-style `mode` and `describe` through rank 15 using Fortran multi-output
  subroutines. `mode` supports one-based axis reductions, whole-array and keepdims
  companions, NaN-as-a-category propagation, omission, and raising. `describe` reports
  nobs/min/max/mean/variance/skewness/Fisher-kurtosis with arbitrary `ddof`, optional
  bias correction, whole-array reduction, and all three NaN policies. Added rank-15,
  floating-point-trap, empty-mode, masked-invalid, and child-process error regressions.

- Added SciPy-style `entropy` and `differential_entropy` through rank 15. `entropy`
  normalizes probability slices and computes Shannon or Kullback-Leibler entropy,
  supports lower-rank trailing and same-rank singleton broadcasting of `qk`, optional
  logarithm bases, one-based axis reductions, whole-array/keepdims companions, and
  paired `nan_policy='propagate'` / `'omit'` / `'raise'` semantics. `differential_entropy`
  implements the Vasicek, Van Es, Ebrahimi, Correa, and SciPy-compatible automatic
  estimators with optional window length/base, axis/whole-array/keepdims reductions,
  and all three NaN policies. Added SciPy 1.17.0 numerical, rank-15, broadcasting,
  floating-point-trap, and child-process error regressions.

- Expanded `rankdata` from its original vector-only implementation to rank-1 through
  rank-15 arrays, with whole-array flattening, explicit one-based axis ranking, all five
  SciPy tie methods, and `nan_policy='propagate'` / `'omit'` / `'raise'`. Added
  `tiecorrect` for Mann-Whitney/Kruskal-Wallis tie factors, plus SciPy 1.17.0 numerical,
  NaN, rank-15, floating-point-trap, and child-process error regressions.

- Added SciPy-style frequency statistics `percentileofscore`, `scoreatpercentile`,
  `cumfreq`, and `relfreq`. `percentileofscore` supports `rank`, `weak`, `strict`,
  and `mean` tie conventions, NaN policies, and score arrays through rank 15.
  `scoreatpercentile` supports whole-array and one-based axis reductions, scalar or
  vector percentiles, `fraction` / `lower` / `higher` interpolation, and closed
  value limits through rank 15. `cumfreq` and `relfreq` support flattened rank-15
  input, explicit/default real limits, weights, and SciPy-compatible extrapoint and
  normalization behavior. Added SciPy 1.17.0 numerical, NaN, weighted-histogram,
  rank-15, floating-point-trap, and child-process error regressions.

- Added SciPy-style circular descriptive statistics `circmean`, `circvar`, and
  `circstd` through rank 15. The new shared assumed-rank engine supports arbitrary
  nonzero `high-low` angular periods, whole-array defaults, explicit one-based
  axis reductions, fixed-rank keepdims, `nan_policy='propagate'` / `'omit'` /
  `'raise'`, and `circstd(normalize=.true.)`. Added SciPy 1.17.0 numerical,
  custom-range, NaN, empty-slice, rank-15, and child-process error regressions.

- Added SciPy-style `kstat` and `kstatvar` through rank 15. `kstat` supports
  orders 1-4 and `kstatvar` orders 1-2, with whole-array `axis=None` behavior,
  explicit one-based axis reductions, fixed-rank keepdims companions, and
  `nan_policy='propagate'` / `'omit'` / `'raise'`. The shared vector kernel
  uses centered-moment identities algebraically equivalent to the standard
  k-statistic power-sum formulas. Added SciPy 1.17.0 whole-array, axis, NaN,
  insufficient-sample, rank-15, and error-path regressions.

- Added SciPy-style proportion trimming `trim_mean`, `trimboth`, and `trim1`
  through rank 15. Shared assumed-rank kernels sort each reduction slice by
  value, round each requested cut down with `floor(proportiontocut*n)`, and
  provide explicit `*_all` companions for SciPy `axis=None`; `trim_mean` also
  has fixed-rank keepdims and `nan_policy='propagate'` / `'omit'` / `'raise'`.
  Added SciPy 1.17.0 numerical, NaN, axis, whole-array, and rank-15 regressions.

- Added SciPy-style trimmed value-limit statistics `tmean`, `tvar`, `tstd`,
  `tsem`, `tmin`, and `tmax` through rank 15. Shared slice kernels implement
  lower/upper value bounds, inclusive/exclusive endpoints, `ddof`, and
  `nan_policy='propagate'` / `'omit'` / `'raise'`; rank-split wrappers add
  one-based axis reduction, fixed-rank keepdims, and explicit `*_all`
  whole-array companions. `tmean` itself follows SciPy's `axis=None` default,
  while the other five default to Fortran axis 1. Added SciPy 1.17.0
  numerical, rank-15, facade, NaN, and floating-point-trap regressions.

- Added SciPy-style `nan_policy='raise'` across the descriptive-statistics
  compatibility layer. The ten descriptive helpers (`gmean` through `kurtosis`)
  now terminate with Fortran 2018 `error stop` when NaNs are encountered,
  including NaN broadcast weights through rank 15. `zscore`, `gzscore`, `zmap`,
  `gstd`, whole-array companions, and tuple-axis `gstd_keepdims` implement the
  same policy; `zmap` raises only for NaNs in the comparison data, matching
  SciPy. Added a child-process regression that verifies expected nonzero
  termination without causing the FPM parent test to fail.

- Added whole-array standardization companions `zscore_all`, `gzscore_all`,
  `zmap_all`, and `gstd_all` through rank 15 as the fixed-signature Fortran
  analogue of SciPy `axis=None`. Added tuple-axis geometric-standard-deviation
  support through `gstd_keepdims(x,dims)`, where `dims` is a unique list of
  one-based dimensions retained with extent one. SciPy 1.17.0 whole-array,
  tuple-axis, NaN-policy, invalid-dimension, and rank-15 regressions were added.
  Also explicitly re-exported `zmap` from the private `scifort_stats` facade.

- Added SciPy-style `zmap` through rank 15. The shared standardization engine
  computes location and scale from a separate comparison array, accepts a
  different comparison extent along the reduced axis, broadcasts singleton
  non-reduced dimensions, and preserves `ddof` plus NaN propagation/omission.
  Added SciPy 1.17.0 vector, broadcast, NaN, constant-reference, invalid-shape,
  rank-15, and floating-point-trap regressions.

- Added SciPy-style `zscore`, `gzscore`, and `gstd` through rank 15, plus
  fixed-rank `gstd_keepdims`. A shared assumed-rank engine handles one-based
  axis slicing, scaled standardization, `ddof`, NaN propagation/omission, and
  zero extents; thin rank-specific wrappers preserve static Fortran result
  ranks. Added SciPy 1.17.0 references, all-rank smoke tests, and IEEE
  invalid/zero/overflow trapping coverage through the public `scifort_stats`
  facade.

- Added SciPy-compatible character `scale='normal'` overloads for `iqr` and
  `median_abs_deviation` through rank 15, including whole-array IQR, one-based
  axes, fixed-rank keepdims, and `nan_policy='propagate'` / `'omit'`. The
  overloads delegate to the existing numeric-scale kernels using the standard
  normal 75th-percentile conversion. Added SciPy 1.17.0 reference coverage for
  every rank, custom IQR ranges, invalid names, and case-insensitive `normal`.

- Extended `sem`, `variation`, `iqr`, `median_abs_deviation`, `moment`, `skew`,
  and `kurtosis` from rank 4 through the standard Fortran rank-15 ceiling.
  Ranks 5-15 share one assumed-rank slice/index engine and thin fixed-result-rank
  wrappers, with one-based axes, fixed-rank keepdims, whole-array IQR, and
  `nan_policy='propagate'` / `'omit'`. Added all-rank smoke coverage and
  SciPy 1.17.0 rank-15 references for non-default axes, nonlinear moments/skew,
  SEM, NaN policies, invalid axes, and floating-point-trap builds.

- Extended weighted `gmean`, `hmean`, and `pmean` from rank 3 through the
  standard Fortran rank-15 ceiling. Rank-4 through rank-15 wrappers accept an
  assumed-rank `weights(..)` dummy and use a shared trailing joint-broadcast
  engine, so scalar through full-rank weights can expand singleton data
  dimensions without materializing broadcast arrays. Added weighted keepdims,
  `nan_policy='propagate'` / `'omit'`, invalid-broadcast/rank checks, SciPy
  1.17.0 rank-4/rank-15 references, and smoke coverage across every rank.

- Extended unweighted `gmean`, `hmean`, and `pmean` from rank 4 through the
  standard Fortran rank-15 ceiling. Ranks 5-15 use one assumed-rank flatten/index
  engine and thin fixed-result-rank wrappers, with one-based axes, fixed-rank
  keepdims companions, and `nan_policy='propagate'` / `'omit'`. Split the
  descriptive implementation into low-rank and high-rank modules behind the
  unchanged `scifort_descriptive` facade to keep optimized compile times practical.
  Added SciPy 1.17.0 high-rank references, zero-extent coverage, and public-facade
  tests.

- Extended all ten SciPy-style descriptive helpers to unweighted rank-4 arrays
  with one-based axes 1 through 4, fixed-rank `*_keepdims` companions,
  whole-array `iqr(axis=None)` behavior, and `nan_policy='propagate'` /
  `'omit'`. Added SciPy 1.17.0 rank-4 axis, keepdims, invalid-axis, NaN,
  and floating-point-trap regressions through the public `scifort_stats` facade.

- Extended rank-3 `gmean`, `hmean`, and `pmean` with rank-1/rank-2/rank-3
  NumPy-style trailing broadcast weights, including joint expansion of singleton
  data dimensions, weighted `nan_policy='propagate'` / `'omit'`, and fixed-rank
  `*_keepdims` companions. Added SciPy 1.17.0 broadcast, expansion, invalid-shape,
  weighted-NaN, and floating-point-trap regressions through `scifort_stats`.

- Extended the SciPy-style descriptive helpers to unweighted rank-3 arrays with
  one-based axis selection, fixed-rank `*_keepdims` companions, whole-array
  `iqr(axis=None)` behavior, and `nan_policy='propagate'` / `'omit'`. Added
  SciPy 1.17.0 rank-3 axis, keepdims, invalid-axis, propagation, omission, and
  floating-point-trap regressions through the public `scifort_stats` facade.

- Added `nan_policy='propagate'` and `nan_policy='omit'` overloads for vector and
  rank-2 `gmean`, `hmean`, `pmean`, `sem`, `variation`, `iqr`,
  `median_abs_deviation`, `moment`, `skew`, and `kurtosis`, including their
  fixed-rank `*_keepdims` companions. Weighted means omit an observation when
  either its value or its broadcast weight is NaN. Added SciPy 1.17.0 vector,
  axis, keepdims, weighted-NaN, all-NaN, invalid-policy, and FPE-trap regressions.
  Also sequenced `pmean`, `skew`, and `kurtosis` NaN/finite checks explicitly so
  propagation does not depend on compiler short-circuit evaluation.

- Extended the SciPy-style descriptive helpers to rank-2 arrays with native
  one-based `axis` selection and fixed-rank `*_keepdims` companions. Weighted
  `gmean`, `hmean`, and `pmean` accept rank-1/rank-2 trailing broadcast weights;
  matrix `iqr` without an axis reduces all elements like SciPy `axis=None`. The
  mean kernels now isolate NaN tests before ordered comparisons so propagation is
  safe under floating-point exception traps. Added SciPy 1.17.0 axis, keepdims,
  broadcasting, invalid-axis, and slice-NaN regressions.

- Added vector SciPy-style descriptive helpers `gmean`, `hmean`, `pmean`,
  `sem`, `variation`, `iqr`, `median_abs_deviation`, `moment`, `skew`, and
  `kurtosis`, re-exported through `scifort_stats`. Weighted mean helpers accept
  finite same-length weights; skew/kurtosis include SciPy's bias-correction
  options, and power-mean overflow/underflow is handled without floating-point
  exceptions in trapping builds.

- Extended weighted real/complex `logsumexp`, signed/phase-returning companions,
  and fixed-rank keepdims from rank 5 through the standard Fortran rank-15 limit.
  Scalar through rank-15 real/complex weights use NumPy-style trailing-dimension
  joint broadcasting and may expand singleton input dimensions. Weighted ranks
  6 through 15 share a broadcast-map/reduction engine split into rank-6--10 and
  rank-11--15 wrapper modules behind the unchanged public reduction facade.

- Extended unweighted real/complex `logsumexp`, signed/phase-returning companions,
  fixed-rank keepdims, `softmax`, and `log_softmax` from rank 6 through the
  standard Fortran rank-15 limit. Ranks 7 through 15 use a shared flattened
  reduction layout/grouping engine, preserving arbitrary unique tuple-axis order,
  empty-axis semantics, zero extents, and invalid/duplicate-axis rejection without
  duplicating rank-specific numerical kernels. SciPy 1.17.0 references cover the
  rank-15 real and complex paths, with smoke coverage for every intermediate rank.
  The implementation is split into rank-1--6 and rank-7--15 modules behind the
  existing `scifort_special_reductions` facade so optimized builds do not force one
  very large compilation unit; the public module and generic procedure names are
  unchanged.

- Added unweighted real/complex rank-6 `logsumexp`, signed/phase-returning
  companions, fixed-rank keepdims, `softmax`, and `log_softmax`. Rank-6 tuple
  dimension lists support any unique subset/order of the six Fortran dimensions,
  including empty-axis, zero-extent, and invalid-axis behavior.

- Added weighted real/complex rank-5 `logsumexp`, signed/phase-returning
  companions, and keepdims variants. Scalar through rank-5 weights use
  NumPy-style trailing-dimension joint broadcasting, including singleton input
  expansion, zero extents, complex weights for real inputs, and real weights for
  complex inputs.

- Added unweighted real/complex rank-5 `logsumexp`, signed/phase-returning
  companions, fixed-rank keepdims, `softmax`, and `log_softmax`. Rank-5 tuple
  dimension lists support any unique subset/order of the five Fortran dimensions.

- Added weighted real/complex rank-4 `logsumexp`, signed/phase-returning
  companions, and keepdims variants. Scalar through rank-4 weights use
  NumPy-style trailing-dimension joint broadcasting and can expand singleton
  dimensions of the rank-4 input; real inputs accept complex weights and complex
  inputs accept either real or complex weights.

- Added unweighted real/complex rank-4 `logsumexp`, signed/phase-returning
  companions, fixed-rank keepdims, `softmax`, and `log_softmax`. Rank-4 tuple
  dimension lists support any subset/order of the four Fortran dimensions.

- Extended weighted `logsumexp` broadcasting so scalar/rank-1/rank-2/rank-3
  weights can expand singleton input dimensions for vectors and rank-2/rank-3
  arrays, matching NumPy/SciPy joint broadcasting. This applies to real and
  complex inputs/weights and to signed/phase-returning and keepdims variants.

- Added rank-3 real/complex `logsumexp`, signed/keepdims companions, `softmax`,
  and `log_softmax`. Rank-3 weights support NumPy-style scalar and trailing
  rank-1/rank-2/rank-3 broadcasting, and tuple dimension lists can reduce any
  subset of the three Fortran dimensions.

- Extended real/complex `logsumexp` to scalar, trailing-vector, and singleton-
  dimension weight broadcasting for vectors/rank-2 arrays, including complex
  weights with real inputs. Added `logsumexp_keepdims` and
  `logsumexp_signed_keepdims` fixed-rank companions plus tuple-style rank-2
  dimension lists for `softmax`/`log_softmax`, including empty-axis semantics.

- Added complex vector/rank-2 `logsumexp`, `softmax`, and `log_softmax`.
  Complex `logsumexp` accepts real or complex same-shape weights, and
  `logsumexp_signed` returns a real log-magnitude plus a unit complex phase,
  matching SciPy's complex `return_sign=True` convention. The complex reduction
  kernel uses magnitude shifting and compensated real/imaginary summation, and
  zero weights are skipped before exponentiation.

- Expanded `logsumexp`, `softmax`, and `log_softmax` beyond their original
  one-dimensional forms. Real rank-2 arrays now support whole-array or single-
  dimension reductions/normalization; `logsumexp` accepts same-shape real
  weights, and `logsumexp_signed` exposes SciPy-style log-absolute-value and
  sign outputs for subtraction-safe weighted sums. Zero weights are removed
  before exponentiation so masked/disabled terms cannot overflow.

- Added SciPy-style QMC `lloyd` post-sampling optimization. A public
  `qmc_lloyd` routine performs centroidal Voronoi relaxation with SciPy's
  decay/stopping convention, while a new dependency-free N-dimensional
  Bowyer-Watson Delaunay layer supplies finite Voronoi vertices. Halton,
  Sobol, LatinHypercube, and unit-cube PoissonDisk now accept
  `optimization='lloyd'`. Generic nondegenerate geometry is validated against
  Qhull without `QJ`; degenerate/cospherical cases can differ from SciPy because
  SciFort deliberately avoids Qhull's randomized joggle.

- Added SciPy-style QMC `random-cd` post-sampling optimization. A public
  `qmc_random_cd` routine performs explicit-state pairwise coordinate exchanges
  with centered-discrepancy acceptance and optional work counters. Halton, Sobol,
  LatinHypercube, and unit-cube PoissonDisk now accept
  `optimization='random-cd'`; optimized Halton/Sobol draws take an explicit
  draw-time RNG state because SciFort does not keep hidden mutable RNG state.

- Completed the listed legacy `odeint` facade controls with `ixpr`, `mxhnil`, and
  `printmessg`. Method-switch reporting is driven by actual LSODA Adams/BDF state
  changes, `mxhnil=0` selects the DLSODA default cap of 10 warning messages, and
  `printmessg=.true.` emits the final integration status. Printing controls leave
  the solution and work counters unchanged. Because SciFort treats an LSODA
  too-small-step condition as terminal, only one hnil warning can currently be
  observed in a single integration.

- Added `odeint(..., full_output=.true.)` diagnostic histories compatible with
  SciPy/ODEPACK names: `hu`, `tcur`, `tolsf`, `tsw`, `nst`, `nfe`, `nje`, `nqu`,
  and `mused`, plus scalar `imxer`, `lenrw`, and `leniw`. LSODA now records the
  step size, order, method actually used, and latest method-switch time at each
  requested output boundary. Regressions cover exact SciPy 1.17.0 nonstiff work
  counts/orders, stiff BDF histories, repeated initial outputs, and signed
  backward-step diagnostics.

- Replaced the initial pooled `odeint` `mxstep` approximation with exact legacy
  per-output-call accounting. The LSODA driver now resets a successful-step counter
  whenever the next requested output is satisfied while preserving its internal
  Adams/BDF history; `mxstep=0` selects the DLSODA/SciPy default of 500. A SciPy
  1.17.0 parity regression matches the `mxstep=6` excess-work boundary at 10
  accepted steps/21 RHS calls and the `mxstep=7` success case at 11/23.

- Extended the legacy `odeint` facade with `tcrit`, `mxordn`, and `mxords`.
  LSODA now accepts a hard critical-time bound without resetting its accumulated
  Adams/BDF history, and its nonstiff/stiff order controllers honor configurable
  ceilings up to the ODEPACK maxima 12 and 5. Added regressions that verify no
  RHS evaluation crosses an active critical time, first-order caps increase work
  on representative nonstiff/stiff systems, and invalid controls are rejected.

- Added scalar-or-componentwise-vector `rtol` across the real and complex
  `solve_ivp` stack. RK23/RK45/DOP853, Radau, BDF, complex BDF, and both LSODA
  branches now use per-component relative tolerances in adaptive error weights
  and startup selection; implicit solvers use `maxval(rtol)` only for their
  scalar Newton stopping heuristic. The `odeint` facade inherits the same vector
  `rtol` support. Added shape/value validation and multiscale real, complex, and
  `odeint` regressions.

- Added a legacy `odeint`-style LSODA facade. `odeint` preserves the historical
  `func(y,t)` callback order and returns time-major `(nt,nstate)` output on
  monotone increasing or decreasing grids with repeated times; `odeint_tfirst`
  provides the type-safe `func(t,y)` form. The facade supports SciPy-like default
  tolerances, scalar or vector `atol`, dense or packed-banded `Dfun`,
  `col_deriv`, `ml`/`mu`, `h0`, `hmax`, `hmin`, and an initial pooled `mxstep`
  budget that was later replaced by exact per-output-call accounting. Added SciPy
  1.17.0 numerical-reference regressions for the
  documented pendulum example, backward/repeated-time output, banded Jacobians,
  t-first callbacks, controls, and invalid arguments.

- Added state-dependent CSC Jacobian callbacks and direct sparse numerical
  Jacobian storage for implicit `solve_ivp`. Real `BDF`/`Radau` and complex
  `BDF` accept `jac_sparse=` callbacks returning the existing one-based CSC
  containers; every callback result is validated before factorization and each
  invocation contributes to `njev`. When `jac_sparsity=` is supplied without an
  analytic or constant Jacobian, finite-difference values are now written
  directly into CSC storage and solved with the native sparse LU path rather
  than materializing an `n` by `n` Jacobian. Added regressions verify callback
  accounting, malformed callback rejection, and absence of dense Jacobian
  storage in sparsity-driven BDF/Radau initialization.

- Added native one-based CSC constant Jacobian containers and sparse direct
  factorization for `solve_ivp`. Real `BDF`/`Radau` accept `ivp_sparse_matrix`
  and complex `BDF` accepts `ivp_complex_sparse_matrix` through
  `jac_sparse_matrix=`; `ivp_csc_matrix(...)` constructs CSC values that
  `solve_ivp` validates. Constant sparse Jacobians remain sparse in solver state,
  report `njev=0`, and now factor their shifted Newton matrices with an original
  SciFort sparse-row LU implementation using partial row pivoting and sparse
  fill-in. The CSC path no longer allocates or assembles the ordinary dense LU
  workspace.

- Added optional batched RHS callbacks for implicit `solve_ivp` numerical
  Jacobians. Real `BDF`/`Radau` and complex `BDF` accept
  `fun_vectorized(t,y(:, :),dydt(:, :))`; dense differences perturb all columns
  in one batch, while `jac_sparsity=` batches one column per independent group.
  The scalar RHS remains the stepping callback, Jacobian probes remain excluded
  from `nfev`, and unsupported explicit/LSODA use is rejected.

- Added constant dense Jacobian matrices and finite-difference Jacobian sparsity
  patterns to `solve_ivp`. Real `BDF`/`Radau` and complex `BDF` accept
  `jac_matrix=` as SciFort's Fortran spelling of SciPy's constant array-like
  Jacobian form; the matrix is stored once and does not increment `njev`.
  `jac_sparsity=` accepts an `n` by `n` logical structural pattern and greedily
  groups columns with disjoint affected rows so one RHS evaluation can perturb
  several state components. LSODA continues to use `jac=` with optional
  `lband`/`uband`, matching SciPy's method-specific restriction. Invalid shapes,
  unsupported methods, and simultaneous `jac=`/`jac_matrix=` are rejected.

- Added SciPy-style componentwise absolute tolerances to real and complex
  `solve_ivp`. The `atol=` argument now accepts either a scalar or a rank-one
  array of length `size(y0)`; explicit Runge-Kutta methods, Radau, real/complex
  BDF, and both LSODA branches use the componentwise scales consistently,
  including finite-difference Jacobian perturbations and initial-step selection.
  Invalid ranks, lengths, nonfinite entries, and negative entries are rejected.
  Existing scalar calls retain their public API and regression behavior.

- Made the LSODA ODEPACK-style BDF next-step cap overflow-safe when
  `max_step` is effectively unlimited. The previous mathematically harmless
  `max_step/h` expression could overflow under `-ffpe-trap=overflow`; the new
  comparison applies the same cap without forming an unbounded quotient.

- Completed the remaining LSODA stiff-controller parity work with an LSODA-only
  ODEPACK-style BDF Nordsieck recurrence. The stiff branch now applies `DCFODE`
  BDF correction and `TESCO` constants directly, advances the Nordsieck history
  with Pascal prediction, uses `DSTODA` chord convergence with up to three
  corrections and first-correction acceptance, and follows the ODEPACK
  `IALTH`/lower-current-upper order and step-size selection rules. Convergence
  failures use the ODEPACK factor-of-four retry. The representative stiff cosine
  regression now reproduces SciPy 1.17.0 LSODA's 116 RHS evaluations, 10 Jacobian
  evaluations, 10 LU factorizations, and 59 stored points; the numerical-Jacobian
  version reports 126 RHS evaluations, accounting exactly for ten scalar
  perturbations. Standalone `method='BDF'` remains on the existing SciPy/NDF
  implementation and retains its reference work counts.

- Continued the LSODA stiff-method parity work at the Adams-to-BDF handoff.
  A successful switch now carries the accepted Adams Nordsieck polynomial into
  the SciFort BDF backward-difference state at the ODEPACK-selected order and
  ideal BDF step ratio instead of restarting at order 1. Jacobian formation is
  deferred until the first actual BDF trial, matching `DSTODA` switch timing.
  While the LSODA stiff branch is active it now uses pure BDF coefficients
  rather than SciPy's NDF enhancement, freezes tolerance weights from the
  accepted state, and uses ODEPACK's weighted maximum norm. Internally generated
  dense and banded Jacobians use the `DPRJA` perturbation rule, reuse the
  predictor RHS in the first Newton iteration, and include every perturbation
  callback in `nfev`; standalone `method='BDF'` retains its existing SciPy/NDF
  behavior and work counts. Added `test/test_integrate_lsoda_stiff_switch.f90`
  and strengthened finite-difference accounting checks in
  `test/test_integrate_ivp.f90`. The later native stiff-controller batch extends
  this handoff work to the complete representative ODEPACK/SciPy trajectory.

- Matched the nonstiff LSODA startup and error-control path more closely to
  public-domain ODEPACK `DLSODA`/`DSTODA`. Automatic LSODA startup now uses
  DLSODA's initial-step formula from the already evaluated initial RHS instead
  of the generic `solve_ivp` startup probe, so it no longer spends an extra RHS
  evaluation. Adams error weights are frozen from the accepted state at the
  start of each trial, tolerance-scaled controller norms use DLSODA's weighted
  maximum norm, and `RMAX` follows the ODEPACK startup/normal/error-failure
  limits of 10000/10/2. These changes reproduce SciPy 1.17.0 LSODA work counts
  on scalar exponential regressions at both `rtol=1e-6` and `rtol=1e-8`, and
  on an eight-component nonstiff growth regression. Added
  `test/test_integrate_lsoda_work_counts.f90`. A later stiff-controller batch
  adds the corresponding representative stiff work-count parity.

- Closed the LSODA BDF-to-Adams state-transfer and stiff-Jacobian refresh gaps.
  A successful BDF-to-Adams switch now converts the current SciPy-style BDF/NDF
  backward-difference polynomial into an Adams Nordsieck power-series state at
  the current endpoint, preserves the selected method order through order 5,
  and enforces the Adams `h*f` first-derivative column before continuing. This
  avoids the previous order-1 history restart after returning to a nonstiff
  region. LSODA's BDF branch now also enables ODEPACK-style Jacobian freshness
  checks: reevaluate after 20 accepted stiff steps or when the Newton coefficient
  changes by more than 30 percent. Standalone `BDF` retains its existing SciPy
  refresh policy. Added `test/test_integrate_lsoda_transfer.f90` for multistep
  state transfer, the first post-switch Adams step, periodic refresh, and the
  coefficient-change refresh trigger. Exact LSODA work-counter parity remains.

- Brought LSODA Adams error-test retries closer to public-domain ODEPACK `DSTODA`.
  A failed step now compares the ODEPACK current-order and one-order-lower retry
  factors using the retained Nordsieck column and `TESCO(1,q)`, applies the `SM1`
  stability limits, forbids an order increase, and caps the retry factor at 0.2
  after the second consecutive error-test failure. A converged-but-rejected
  functional correction also retains its Lipschitz estimate for the stability
  test. The existing conservative four-failure restart remains unchanged. Added
  `test/test_integrate_lsoda_rejections.f90` to exercise lower-order retry choice
  and the repeated-failure step cap.

- Completed the successful-step persistent Nordsieck path for LSODA's Adams branch.
  The solver now starts at method order 1, carries scaled derivatives persistently
  through orders 1--12, power-rescales them on step changes, applies the Pascal
  prediction and ODEPACK `EL` correction, and uses `IALTH`/`TESCO` to delay and
  compare lower/current/upper order choices. The retained Nordsieck polynomial now
  provides native nonstiff LSODA dense output instead of cubic Hermite interpolation.
  The functional-corrector contraction tests apply at all Adams orders and retain the
  converged functional-iteration derivative without an extra accepted-endpoint RHS
  refresh. Nonuniform derivative/time history remains available for reconstruction and
  conservative recovery but is no longer rebuilt into active state after every
  successful step. Regressions cover delayed order changes, a direct order-12
  Nordsieck step/dense state, variable-step history, signed-step rescaling, backward
  integration, and a deliberately divergent corrector.

- Moved LSODA method switching closer to public-domain ODEPACK `DSTODA`: Adams
  proposals now respect the order-dependent `SM1` stability envelope, method-switch
  tests are delayed for 20 successful steps after initialization or a switch, and
  Adams-to-BDF switching is skipped above order 5. The common finite-error switch
  path now compares the ideal Adams and BDF step factors using ODEPACK's `CM1`/`CM2`
  error constants and ratio test. A directional local Lipschitz estimate from the
  existing Adams predictor/corrector evaluations supplies the stability estimate
  without extra RHS or Jacobian calls. Accepted BDF steps also return the endpoint
  derivative implied by the converged BDF equation, removing a redundant LSODA RHS
  evaluation. Direct regressions check the hysteresis, stability cap, BDF activation,
  and endpoint-derivative reuse. The accepted-step switch path now also retains the
  latest Adams correction-based Lipschitz estimate for the ODEPACK-style
  roundoff-pollution branch: when local error information is dominated by roundoff,
  Adams switches to BDF only if its proposed step is already stability-limited.
  Rejected Adams trials no longer probe a Jacobian or switch methods; after four
  consecutive error-test failures SciFort conservatively restarts low-order history
  and cuts the retry step by a factor of ten. Error-test retries before that restart
  now use the ODEPACK lower/current-order comparison and repeated-failure step cap.
  Remaining parity work includes BDF-to-Adams state transfer/Jacobian refresh and
  exact work counters.

- Extended the SciFort LSODA nonstiff branch to variable-step Adams history.
  Accepted derivatives are retained with their time offsets, Adams-Bashforth and
  Adams-Moulton weights are regenerated for the actual nonuniform nodes, and order
  adapts independently of history depth instead of rebuilding through RK45 whenever
  the step changes. Nonstiff switching uses directional secant and correction-based
  stiffness indicators without routine Jacobian sampling, eliminating Jacobian work
  on smooth nonstiff problems. The later persistent-Nordsieck work above supersedes
  the original rebuild-after-accept design; the retained nonuniform history is now a
  reconstruction and recovery aid.

- Added an initial real-system `method='LSODA'` compatibility path to `solve_ivp`.
  The predictor/corrector implementation is SciFort code: it establishes an
  equal-step derivative history with adaptive RK45 startup and advances nonstiff regions with
  Adams-Bashforth predictors that ramp from order 4 through 11 and Adams-Moulton
  correctors through order 12, and automatically switches to SciFort's existing
  variable-order BDF/NDF solver when
  a dense Jacobian stiffness indicator crosses the switching threshold. Analytic or
  finite-difference dense Jacobians are supported for switching/BDF work. The initial
  implementation used endpoint Hermite interpolation on Adams steps; the persistent
  Nordsieck work above supersedes that startup/dense-output design while reusing the
  existing `t_eval`, backward-integration, dense-output, and event machinery. LSODA
  now also accepts `lband`/`uband`: analytic callbacks use SciPy's
  packed-band layout, while finite differences group nonoverlapping band columns to
  reduce RHS work. `min_step` supplies a floor for rejected Adams/BDF retry steps.
  Predictor/corrector weights are generated by integrating Lagrange basis
  polynomials, avoiding imported coefficient tables. ODEPACK parity work remains
  for Nordsieck step/order adaptation, remaining ideal-step switching comparisons,
  and exact work counters.

- Added complex-system `method='BDF'` to `solve_ivp`, completing the complex
  method set supported by SciPy among SciFort's current IVP methods. The new
  stepper adapts SciPy 1.17.0's variable-order 1--5 BDF/NDF algorithm, Newton
  iteration, native dense output, complex dense LU solves, and real-axis
  finite-difference Jacobian convention. Complex BDF accepts an optional
  trailing `jac(t,y,J)` callback and supports `t_eval`, retained dense output,
  backward integration, and existing scalar/vector event machinery.

- Added complex-system overloads of `solve_ivp` for `RK23`, `RK45`, and
  `DOP853`. The generic `scifort_integrate%solve_ivp` now resolves from the
  type of `y0`; complex runs support `t_eval`, backward integration, retained
  dense output, scalar/vector real-valued events, direction filters, and
  Nth-occurrence termination. New `solve_ivp_complex_result`,
  `ivp_complex_dense_solution`, and `ivp_complex_event_set` types preserve
  complex state/event storage while keeping the existing real API unchanged.

- Extended `solve_ivp` event compatibility with SciPy 1.17.0 positive-integer
  terminal semantics. Scalar `event_max_occurrences` and vector
  `ivp_event_set%max_occurrences(:)` can terminate on the Nth occurrence while
  preserving existing logical terminal controls. `solve_ivp_result` now reports
  per-event counts and provides `event_times_for(i)` / `event_states_for(i)`
  grouping helpers over the existing flat event arrays. Same-step roots are
  ordered and truncated at the first event that reaches its occurrence limit,
  and initial endpoint direction filtering follows SciPy's inclusive crossing
  tests.

- Added real-system `method='Radau'` to `solve_ivp`, adapting SciPy 1.17.0's
  fifth-order three-stage Radau IIA collocation method, embedded third-order
  error estimator, transformed Newton iteration, cubic dense output, and dense
  finite-difference Jacobian strategy under BSD-3-Clause. Analytic Jacobian
  callbacks, `njev`/`nlu` accounting, `t_eval`, retained dense solutions,
  backward integration, and event-root interpolation are supported.

- Added real-system `method='BDF'` to `solve_ivp`, adapting SciPy 1.17.0's
  variable-order 1--5 quasi-constant-step BDF/NDF algorithm, Newton iteration,
  dense output, and dense finite-difference Jacobian strategy under BSD-3-Clause.
  An optional Fortran Jacobian callback is supported; `solve_ivp_result` now also
  reports `njev` and `nlu`. BDF dense output is used by `t_eval`, retained
  continuous solutions, and event root location.

- Replaced RK45's cubic-Hermite interpolation with SciPy 1.17.0-compatible
  quartic dense output for `t_eval`, retained dense solutions, and event-root
  states. The continuous extension reuses the existing seven Dormand-Prince
  stage derivatives, so dense interpolation adds no right-hand-side evaluations.

- Added native seventh-order DOP853 dense output for `t_eval`, retained dense
  solutions, and event root location. Added `ivp_event_set` for multiple event
  functions through one vector callback, with per-event terminal/direction controls,
  flattened `t_events`/`y_events`, and `event_indices` identifying each root.
- Added `method='DOP853'` to `solve_ivp`, adapting the Hairer/Wanner DOP853
  eighth-order tableau and embedded 5th/3rd-order error estimator from Jacob
  Williams' permissively licensed modern Fortran implementation. The method
  reuses SciFort's existing adaptive driver and result/event API; a subsequent
  batch adds the method's native seventh-order dense interpolant.
- Started the ODE portion of `scipy.integrate` with a Fortran-native `solve_ivp`
  result/API and adaptive real-system `RK23`/`RK45`. The Bogacki-Shampine 3(2)
  and Dormand-Prince 5(4) tableaux are adapted from Jacob Williams' BSD-3-Clause
  `rklib`; SciFort supplies the controller, storage, `t_eval`, cubic-Hermite
  dense interpolation, backward integration, and scalar event location.
- Added adaptive multidimensional `cubature` with scalar and vector-valued
  callbacks, Genz-Malik degree-7/degree-5 and product `gk15`/`gk21` embedded
  rules, finite/infinite/reversed bounds, user break points, global error-driven
  subdivision, and an explicit max-subdivision status. Genz-Malik formulas and
  Gauss-Kronrod tables are adapted from SciPy 1.17.0 under its BSD-3-Clause
  license; the Fortran region manager and infinite-bound transformations are
  SciFort code.
- Added `quad_vec`, `qmc_quad`, and `tanhsinh` to `scifort_integrate`.
  `quad_vec` is a correctness-first vector facade over the scalar adaptive
  integrator, `qmc_quad` uses independently scrambled Halton replicates with an
  explicit `rng_state`, and `tanhsinh` is an independent Takahasi-Mori
  double-exponential implementation for finite and infinite real intervals.
  The tanh-sinh convergence check includes a transformed-tail guard so slowly
  decaying infinite integrals do not report false convergence at the fixed
  binary64 truncation boundary.
- Started `scipy.integrate` compatibility with one-dimensional `trapezoid`,
  `cumulative_trapezoid`, `simpson`, `cumulative_simpson`, `romb`, `fixed_quad`,
  `newton_cotes`, and scalar adaptive `quad`. The adaptive facade uses the
  BSD-3-Clause modern QUADPACK 2.1.0 package as a pinned FPM dependency and
  supports finite, semi-infinite, whole-line, and reversed integration limits.
- Extended `scipy.integrate` adaptive quadrature with finite-interval break
  points; cosine, sine, algebraic-endpoint, logarithmic-endpoint, and Cauchy
  weights; and nested `dblquad`, `tplquad`, and `nquad`-style integration. The
  specialized modes dispatch through Quadpack2 `dqagp`, `dqawo`, `dqawf`,
  `dqaws`, and `dqawc` without vendoring QUADPACK source.
- Added SciPy-style classical orthogonal-polynomial coefficient constructors
  `legendre`, `jacobi`, `gegenbauer`, `laguerre`, `genlaguerre`, `hermite`,
  `hermitenorm`, `chebyt`, `chebyu`, `chebyc`, `chebys`, `sh_legendre`,
  `sh_jacobi`, `sh_chebyt`, and `sh_chebyu`. They return descending-power
  coefficient arrays matching SciPy `orthopoly1d.coeffs`, support optional
  monic normalization, and retain evaluation/quadrature in the existing
  `eval_*`/`roots_*` interfaces.
- Added complex overloads for `ndtr`, `log_ndtr`, `expm1`, and `log1p`,
  plus mixed real/complex overloads for `xlogy` and `xlog1py`. Complex normal-CDF
  evaluation reuses SciFort's Faddeeva `erfc`/`erfcx` foundation, with scaled
  left-tail and reflected `log1p` branches for a stable principal `log_ndtr`;
  the elementary helpers preserve cancellation-sensitive and zero-multiplier
  semantics while leaving the established real APIs unchanged.
- Added complex-argument overloads for `spherical_jn`, `spherical_yn`,
  `spherical_in`, and `spherical_kn`, including first derivatives. Nonreal values
  use the DLMF half-integer definitions over the complex cylindrical Bessel layer,
  while nonzero real-axis complex inputs delegate to the established real kernels.
  High-order accuracy is hardened with direct complex-I series evaluation, a
  backward Miller recurrence for spherical J near its turning region, and tighter
  scaled-Hankel continuation.
- Added complex-argument overloads for SciPy-style cylindrical and modified
  Bessel `jv`/`yv`/`iv`/`kv`, scaled `jve`/`yve`/`ive`/`kve`, the `jn`
  alias, and arbitrary-order `jvp`/`yvp`/`ivp`/`kvp`. The implementation
  reuses scaled Hankel continuation, adds a direct J power-series branch to
  avoid large-order/small-argument cancellation, and preserves the existing
  real-axis APIs.
- Added real/complex SciPy-style Hankel functions `hankel1`, `hankel2`, scaled
  `hankel1e`/`hankel2e`, and arbitrary-order derivatives `h1vp`/`h2vp`.
  Real arguments reuse the existing `jv`/`yv` kernels. Complex arguments use
  scaled-K connection formulas, scaled Bessel-ODE continuation, and DLMF
  principal-branch continuation across the negative real axis.
- Added real ellipsoidal harmonics `ellip_harm`, `ellip_harm_2`, and
  `ellip_normal`. First-kind Lame-polynomial coefficient construction and
  scaling are adapted from SciPy 1.17.0's BSD-2-Clause kernel, replacing
  its LAPACK eigenproblem with SciFort's existing pure symmetric eigensolver.
  Second-kind and normalization integrals use independently written
  endpoint-removing substitutions with adaptive Gauss-Legendre quadrature.
- Added the current SciPy first-kind Legendre API layer: real/complex
  `legendre_p`, `legendre_p_all`, `assoc_legendre_p_all`,
  `sph_legendre_p_all`, and `sph_harm_y_all`. Table routines allocate native
  zero-based degree/derivative axes and a signed `-m:m` order axis; harmonic
  derivatives use one uniform mixed-derivative tensor. The implementation
  reuses SciFort's independent DLMF recurrences and Taylor-jet derivatives.
- Added the real parabolic-cylinder function `pbwa(a,x,w,wprime)` for DLMF
  `W(a,x)`. The independent implementation uses exact DLMF origin values from
  complex log-gamma ratios and Taylor continuation of
  `y'' + (x*x/4-a)*y = 0`. The public domain follows SciPy's legacy
  `abs(a) <= 5`, `abs(x) <= 5` interface; outside-domain inputs return NaN.
- Added second-kind prolate/oblate spheroidal radial functions `pro_rad2`/
  `obl_rad2` and precomputed-characteristic variants `pro_rad2_cv`/`obl_rad2_cv`.
  The independent implementation uses the DLMF spherical-Bessel-Y expansion at
  a stable large-argument anchor and adaptive Dormand-Prince continuation of the
  spheroidal differential equation toward the prolate boundary or oblate origin.
- Added first-kind prolate/oblate spheroidal radial functions `pro_rad1`/
  `obl_rad1` and precomputed-characteristic variants `pro_rad1_cv`/`obl_rad1_cv`.
  The independent implementation uses the DLMF spherical-Bessel expansion with
  logarithmically scaled negative-order coefficients and an analytic oblate-origin
  limit.
- Added first-kind prolate/oblate spheroidal angular functions `pro_ang1`/
  `obl_ang1` and precomputed-characteristic variants `pro_ang1_cv`/`obl_ang1_cv`.
  The original implementation reconstructs the DLMF Ferrers expansion from the
  tridiagonal eigenvector by split continued-fraction recurrences and enforces the
  Flammer/SciPy center normalization.
- Added prolate/oblate spheroidal characteristic values `pro_cv`/`obl_cv` and
  sequence APIs `pro_cv_seq`/`obl_cv_seq`. The independent implementation uses
  the DLMF finite tridiagonal eigenproblems, a symmetric similarity transform,
  and targeted Sturm bisection rather than importing SciPy's legacy backend.
- Added the four modified/radial Mathieu functions `mathieu_modcem1`,
  `mathieu_modcem2`, `mathieu_modsem1`, and `mathieu_modsem2`, using DLMF
  Bessel cross-product expansions with adaptive anchoring for numerical stability.
- Added the integer-order angular Mathieu family: `mathieu_a`/`mathieu_b`,
  `mathieu_even_coef`/`mathieu_odd_coef`, and degree-argument
  `mathieu_cem`/`mathieu_sem`. The independent implementation solves the four
  DLMF Fourier parity eigenproblems with SciFort's existing symmetric eigensolver
  and applies exact negative-`q` coefficient symmetries.
- Added the order-zero Kelvin family `ber`, `bei`, `ker`, `kei`, their first
  derivatives, combined complex `kelvin` evaluator, and all eight positive zero-table
  routines. The independent implementation uses the DLMF Bessel definitions,
  convergent complex series, large-argument asymptotics, and safeguarded sign bracketing.
- Added `lmbda`, `besselpoly`, and `erf_zeros` to `scifort_special`. The Lambda
  sequence uses its exact `0F1` representation; the weighted Bessel integral uses
  endpoint-aware composite Gauss-Legendre quadrature; and error-function zeros use
  DLMF asymptotic seeds followed by Newton refinement against SciFort's complex `erf`.
- Added `wright_bessel` and `log_wright_bessel` for SciPy-compatible nonnegative real
  parameters. The implementation evaluates the DLMF defining series around its largest
  term in the log domain, with an asymptotic saddle branch for exceptionally large modes.
- Added `bernoulli`, `euler`, and `stirling2` to `scifort_special`. Bernoulli and Euler
  APIs return zero-based number sequences with explicit invalid-argument status; Stirling
  numbers use the direct positive recurrence in binary64 rather than SciPy's documented
  large-`n` approximation.
- Added order-zero Bessel integral pairs `itj0y0`, `it2j0y0`, `iti0k0`,
  and `it2i0k0`. The independent implementation uses DLMF defining series,
  composite Gauss-Legendre continuation, and large-argument Bessel-tail
  asymptotics, reusing SciFort's existing J/Y/K kernels. Divergent endpoint
  values use IEEE infinities rather than legacy finite sentinels.
- Added real Struve `struve(v,x)` and modified-Struve `modstruve(v,x)`,
  together with the historical order-zero integral helpers `itstruve0`,
  `it2struve0`, and `itmodstruve0`. The independent implementation combines
  DLMF power/Bessel series, large-argument asymptotics, exact negative
  half-integer identities, and composite Gauss-Legendre continuation for the
  oscillatory H0 integrals. The large-x `itstruve0` path follows the defining
  integral rather than reproducing observed SciPy 1.17 legacy failures.
- Added SciPy-style Riccati-Bessel sequences `riccati_jn` and `riccati_yn`,
  returning zero-based values and first derivatives through order `n`. They
  reuse SciFort's spherical-Bessel kernels and use the mathematical x=0 limits;
  singular higher-order `riccati_yn` values are IEEE infinities rather than
  legacy finite sentinels.

- Added real/complex Legendre second-kind sequence APIs `lqn` and `lqmn`,
  including first derivatives. Ferrers values inside the unit disk use exact
  Legendre recurrences, while exterior order-zero values use the DLMF 14.3.7
  hypergeometric representation before recurrence in order. Nonreal points on
  the unit circle use the finite canonical complex-Q continuation rather than
  reproducing the infinities/NaNs from SciPy 1.17.0's legacy backend.

- Added real-degree/integer-order Ferrers `lpmv`, normalized
  `sph_legendre_p`, and complex `sph_harm_y`. Noninteger `lpmv` uses
  DLMF hypergeometric representations, integer degrees use recurrence to
  avoid gamma-pole cancellation, and spherical functions use a directly
  normalized recurrence with the Condon-Shortley phase.
- Added integer-degree real/complex `assoc_legendre_p` with SciPy-style
  branch-cut choices 2/3 and optional spherical normalization, plus
  `assoc_legendre_p_derivatives`, `sph_legendre_p_derivatives`, and
  `sph_harm_y_derivatives`. The derivative routines propagate truncated Taylor
  jets through the recurrence instead of using finite differences.

- Added SciPy-style complex `lambertw(z[,k,tol])` for arbitrary integer branches
  and real/complex `wrightomega(z)`. The independent implementation uses DLMF
  branch-point and large-argument starts followed by Halley refinement, with a
  direct Wright-omega path outside the safe exponential range. Exact `-1/e` on
  the real Lambert branch point returns the mathematical value `-1`.

- Added legacy SciPy orthogonal compatibility names: `assoc_laguerre` plus
  `p_roots`, `j_roots`, `cg_roots`, `l_roots`, `la_roots`, `h_roots`,
  `he_roots`, `t_roots`, `u_roots`, `c_roots`, `s_roots`, `ps_roots`,
  `js_roots`, `ts_roots`, and `us_roots`. They delegate directly to the
  existing generalized-Laguerre evaluator and canonical quadrature routines.

- Added historical statistical `scipy.special` wrappers `smirnov`, `smirnovi`,
  `kolmogorov`, `kolmogi`, and `tklmbda`. The wrappers delegate to SciFort's
  existing `ksone`, `kstwobign`, and Tukey-lambda distribution kernels rather
  than duplicating probability algorithms.

- Added Gaussian-quadrature root/weight generators `roots_legendre`,
  `roots_jacobi`, `roots_gegenbauer`, `roots_laguerre`, `roots_genlaguerre`,
  `roots_hermite`, `roots_hermitenorm`, `roots_chebyt`, `roots_chebyu`,
  `roots_chebyc`, `roots_chebys`, `roots_sh_legendre`, `roots_sh_jacobi`,
  `roots_sh_chebyt`, and `roots_sh_chebyu`. General classical families use an
  independent Golub-Welsch construction with SciFort's symmetric eigensolver;
  Chebyshev families use their closed-form nodes and weights.

- Added integer-degree classical orthogonal-polynomial evaluators `eval_legendre`,
  `eval_jacobi`, `eval_gegenbauer`, `eval_laguerre`, `eval_genlaguerre`,
  `eval_hermite`, `eval_hermitenorm`, `eval_chebyt`, `eval_chebyu`,
  `eval_chebyc`, `eval_chebys`, `eval_sh_legendre`, `eval_sh_jacobi`,
  `eval_sh_chebyt`, and `eval_sh_chebyu`. Real/complex argument overloads use
  independent DLMF three-term recurrences; Hermite follows SciPy's real-argument
  interface.

- Added the legacy noncentral-distribution `scipy.special` families `chndtr`,
  `chndtridf`, `chndtrinc`, `chndtrix`, `ncfdtr`, `ncfdtri`, `ncfdtridfd`,
  `ncfdtridfn`, `ncfdtrinc`, `nctdtr`, `nctdtridf`, `nctdtrinc`, and `nctdtrit`.
  Forward CDFs and variate inverses reuse SciFort's existing `ncx2`, `ncf`, and
  `nct` kernels; parameter inverses use safeguarded scalar searches. Noncentral-F
  degree inverses return the smallest positive solution when the CDF equation has
  two roots.

- Added SciPy-style convenience special helpers `cbrt`, `exp10`, `exp2`,
  `expm1`, `log1p`, `powm1`, `radian`, `round`, normalized `sinc`,
  `sindg`/`cosdg`/`tandg`/`cotdg`, `softplus`, and `diric`. The degree
  trigonometric functions use explicit modulo-360 reduction for finite inputs,
  and `cotdg` guards tiny angles before reciprocal overflow.

- Added the remaining legacy inverse-parameter distribution helpers `bdtrik`,
  `bdtrin`, `btdtria`, `btdtrib`, `chdtriv`, `fdtridfd`, `gdtrib`, `nbdtrik`,
  `nbdtrin`, `pdtrik`, `stdtridf`, `nrdtrimn`, and `nrdtrisd`. The parameter
  inverses use safeguarded geometric bracketing around SciFort's existing
  incomplete-beta/gamma and normal kernels. `fdtridfd` returns the smaller
  positive root when the denominator-df equation has two solutions.

- Added SciPy-compatible legacy distribution-special helpers for binomial, negative-binomial, Poisson, chi-square, F, gamma, and Student-t probabilities/inverses: `bdtr`/`bdtrc`/`bdtri`, `nbdtr`/`nbdtrc`/`nbdtri`, `pdtr`/`pdtrc`/`pdtri`, `chdtr`/`chdtrc`/`chdtri`, `fdtr`/`fdtrc`/`fdtri`, `gdtr`/`gdtrc`/`gdtria`/`gdtrix`, and `stdtr`/`stdtrit`. The wrappers reuse SciFort's regularized incomplete beta/gamma and normal kernels, with stable F-coordinate and extreme Student-t handling.

- Added SciPy-style scalar special-function helpers `agm`, `huber`,
  `pseudo_huber`, `ndtri_exp`, `owens_t`, and `voigt_profile`. Owen's T uses
  a 32-point angular quadrature after reciprocal-slope reduction; `ndtri_exp`
  directly inverts the log-normal CDF in the extreme lower tail; and the Voigt
  profile reuses SciFort's complex Faddeeva kernel. Skew-normal CDF evaluation
  now reuses the public `owens_t` implementation instead of a private duplicate.

- Extended the parabolic-cylinder family with SciPy-style real `pbvv(v,x)` values and derivatives, real `pbdv_seq`/`pbvv_seq` sequence wrappers, and analytic complex integer-order `pbdn_seq`. `pbvv` combines exact reciprocal-gamma origin data, ODE continuation, a large-positive-argument DLMF expansion, and an exact negative-axis connection formula. The sequence APIs preserve SciPy's minimum-two-output shape convention; complex `pbdn_seq` follows the analytic DLMF `D_n(z)` definition. Added high-precision scalar, sequence, complex, compatibility-shape, invalid-input, asymptotic, and negative-argument regressions.

- Added SciPy-style real `pbdv(v,x)` functionality through the elemental Fortran subroutine `pbdv(v,x,d,dprime)`, returning `D_v(x)` and its derivative. Positive arguments reuse the new `hyperu` kernel through the DLMF relation; negative arguments use Taylor continuation of the parabolic-cylinder ODE from exact origin values, with a Hermite recurrence for nonnegative integer orders. Added ordinary, positive/negative, high-order, large-argument, elemental-array, origin, and IEEE regressions.

- Added SciPy-style real `hyperu(a,b,x)` for Tricomi's confluent hypergeometric function U. The implementation combines DLMF large-positive-argument asymptotics, inward Taylor continuation of Kummer's differential equation, the small-argument 1F1 connection formula with near-integer-`b` Richardson regularization, and exact negative-integer-`a` terminating polynomials. Added ordinary, small/large argument, difficult continuation, Kummer-transformation, elemental-array, branch-cut, endpoint, and IEEE regressions.

- Added SciPy-style `hyp2f1(a,b,c,z)` for real parameters and real or complex
  binary64 arguments. The implementation evaluates the defining Gauss series
  near the origin and analytically continues the regular solution with local
  Taylor series whose recurrence is derived from the hypergeometric
  differential equation. Complex continuation detours around `z=1` and uses
  the sign of zero on the real-axis branch cut to select the upper or lower
  principal value. Added Gauss-at-one, terminating-polynomial, denominator-pole
  cancellation, real branch-cut, elemental-array, IEEE, and difficult complex
  continuation regressions.

- Added SciPy-style `hyp1f1(a,b,z)` for real parameters and real or complex binary64 arguments. The implementation combines the defining Kummer series, Taylor continuation of the Kummer differential equation, Kummer left-half-plane transformation, terminating-polynomial special cases, and DLMF 13.7 large-argument asymptotics. Added regressions for large positive and negative real arguments, complex continuation and asymptotics, removable real denominator-pole cases, SciPy's complex denominator-pole convention, IEEE values, and elemental arrays.

- Added SciPy-style `hyp0f1(v,z)` for real denominator parameter `v` and real or complex binary64 `z`. Near the origin it uses the defining 0F1 series; larger real arguments use the existing real Bessel kernels through the DLMF connection formula, while larger complex arguments use outward Taylor continuation of the generalized-hypergeometric differential equation. Added SciPy 1.17.0 regressions for oscillatory/large real values, complex continuation, negative noninteger parameters, IEEE values, poles, and elemental arrays.

- Extended the Carlson symmetric elliptic interfaces `elliprc`, `elliprd`,
  `elliprf`, `elliprg`, and `elliprj` to complex binary64 arguments while
  preserving their existing real specifics and real Cauchy-principal-value
  behavior. Complex RC/RF/RD use principal-branch duplication formulas, RG
  uses the Carlson RF/RD identity with branch-preserving square-root factors,
  and RJ uses a scaled adaptive integral to avoid unverified complex recurrence
  branches. Added SciPy 1.17.0 regression coverage for cuts, conjugate cases,
  unbalanced arguments, IEEE values, homogeneity, and elemental arrays.

- Extended `airye` to complex binary64 arguments while preserving the existing
  real specific and its negative-real NaN convention. Moderate recessive Ai/Ai'
  values use Airy-ODE Taylor continuation from a scaled asymptotic anchor to
  avoid power-series cancellation; outer-plane Ai and Bi scaling is applied
  before rotated-Ai connection formulas so exponentially large intermediates
  are not formed. Added SciPy 1.17.0 regression coverage through the Stokes
  sectors, negative complex real axis, elemental arrays, IEEE values, and a
  `-200+300i` scaled stress case.

- Extended `airy` to complex binary64 arguments while preserving the existing
  real specific. The complex kernel uses the Airy-ODE power series for moderate
  arguments and NIST DLMF asymptotic plus rotated-Ai connection formulas in the
  outer complex plane. The paired outputs remain Ai, Ai-prime, Bi, and Bi-prime;
  exact real-axis complex inputs reuse the existing real implementation.

- Extended `spence` to complex binary64 arguments while preserving the existing
  real specific. Complex evaluation follows the principal `Li_2(1-z)` branch
  from DLMF dilogarithm series, inversion, and logarithmic expansion formulas;
  exact negative-real complex inputs preserve signed-zero side selection on the
  branch cut. No complex C ABI or new dependency was added.

- Extended the top-level `sici` and `shichi` interfaces to complex binary64 arguments while
  preserving the existing real specifics. Small complex arguments use the defining DLMF
  power series; larger `sici` values reuse the validated principal E1 kernel and larger
  `shichi` values reuse principal Ei. Complex `Ci`/`Chi` follow the principal-logarithm
  negative-real branch cut, while `Si`/`Shi` remain entire. No complex C ABI or new
  dependency was added.

- Extended the top-level `exp1` and `expi` interfaces to complex binary64 arguments while
  preserving the existing real specifics. Complex `exp1` follows the principal E1 branch
  with SciPy-compatible signed-zero behavior on the negative-real cut; complex `expi`
  follows SciPy's exact-cut convention and analytic continuation off the cut. The original
  SciFort implementation uses the DLMF Ein series and E1 continued fraction with guarded
  complex exponential scaling. No complex C ABI or new dependency was added.

- Extended the top-level `digamma` and `psi` aliases to complex binary64 arguments while
  preserving the existing real specifics. The complex kernel uses DLMF recurrence,
  reflection, and asymptotic expansion with a stabilized complex `cot(pi*z)` evaluation
  near the negative-real poles, including finite behavior for imaginary offsets as small
  as `1e-300`. No complex C ABI or new dependency was added.

- Extended `gamma`, `rgamma`, and `loggamma` to complex binary64 arguments while
  preserving the existing real specifics. Complex `loggamma` follows the
  principal branch with its single negative-real-axis cut, using DLMF recurrence,
  reflection, and Stirling expansion with explicit sheet correction; `gamma` and
  `rgamma` exponentiate that branch with guarded complex overflow handling. No
  complex C ABI or new dependency was added.

- Added real-argument `modfresnelp` and `modfresnelm` paired-output
  interfaces for the modified Fresnel tail integrals and their scaled K forms.
  The implementation reuses the validated Faddeeva kernel for the smooth scaled
  positive branch and an error-free split square for oscillatory phase
  reconstruction. It follows the documented integral for negative lower limits,
  intentionally correcting the first returned component produced by SciPy 1.17.0
  on that branch.

- Added `fresnel_zeros`, `fresnelc_zeros`, and `fresnels_zeros` for the first-quadrant complex zeros of the Fresnel integrals. The original SciFort implementation uses NIST DLMF 7.13.5-7.13.7 asymptotic seeds followed by Newton refinement against the existing complex Fresnel kernel, and intentionally remains finite beyond the range where SciPy 1.17.0 begins returning NaNs.

- Extended the top-level `fresnel` subroutine to complex binary64 arguments while
  preserving the existing real kernel. Complex S/C evaluation follows DLMF 7.5.8
  through the validated complex-error foundation, uses the defining entire series
  in the well-conditioned diagonal sector, and maps exact real/imaginary axes to
  the existing real Fresnel implementation. No complex C ABI or dependency was added.

- Extended the top-level `erf`, `erfc`, `erfi`, and `dawsn` APIs to complex
  binary64 arguments while preserving the existing real kernels. The new
  specifics reuse the independently implemented Faddeeva `wofz` foundation,
  local entire series near the origin, conjugation/reflection symmetries, and
  scaled exponential products to avoid avoidable overflow. No complex C ABI or
  new dependency was added.

- Added complex `wofz` (Faddeeva) support and overloaded top-level `erfcx` so the existing real API now also accepts complex arguments. The independent implementation follows NIST DLMF Chapter 7 using the entire power series, real-axis Taylor continuation, upper-half-plane integral representation, large-argument continued fraction, and reflection symmetry. No third-party implementation code or dependency was added.

- Added real `dawsn`, `erfi`, `erfcx`, and paired-output `fresnel` special
  functions. The independent implementation uses DLMF definitions/asymptotics,
  standard Gauss-Legendre quadrature, stable overflow handling, and compensated
  quadratic-phase reduction for large real Fresnel arguments.

- Added real-branch `y0_zeros`, `y1_zeros`, and `y1p_zeros` interfaces. They
  return the positive real roots together with the documented derivative/value
  arrays while reusing SciFort's existing Y/Y-prime root kernels. SciFort
  returns the mathematical `Y_0'=-Y_1` companion for `y0_zeros`; this corrects
  the opposite-sign numerical companion produced by SciPy 1.17.0 despite its
  documentation. The real-only SciFort routines also remain finite beyond the
  140-entry limit observed in SciPy 1.17.0's legacy specialized routines.
- Added positive real integer-order Bessel-zero table interfaces `jn_zeros`,
  `jnp_zeros`, `yn_zeros`, `ynp_zeros`, combined `jnyn_zeros`, and
  `jnjnp_zeros`. The original SciFort implementation brackets successive simple
  roots of the existing `jv`/`yv`/`jvp`/`yvp` kernels and refines them with
  safeguarded Newton steps; negative integer orders reuse the standard
  parity-equivalent positive-order zeros. `jnjnp_zeros` merges J/J-prime roots
  across integer orders and returns SciPy-style order, serial, and family
  metadata. Regression coverage includes orders through 1000, 1200-root fixed
  order tables, and the full 1200-entry merged J/J-prime compatibility range.
- Added real `spence`, Riemann/Hurwitz `zeta`, and `zetac` interfaces to
  `scifort_special`. `spence(x)` evaluates SciPy's real `Li_2(1-x)` branch
  from DLMF dilogarithm series/connection formulas; `zeta` extends the
  existing Euler-Maclaurin Hurwitz kernel with real Riemann analytic
  continuation and integer-exponent negative-shift continuation. `zetac`
  evaluates `zeta(s)-1` directly for positive `s`, preserving tiny
  large-exponent values beyond SciPy 1.17.0's early underflow.
- Added SciPy-style floating factorial/combinatorial interfaces `factorial`,
  `factorial2`, `factorialk`, `comb`, and `perm` to `scifort_special`.
  Integer factorial/multifactorial/permutation arguments use guarded direct
  products for exact small binary64 results, while real factorial and
  generalized combination/permutation paths reuse the existing Gamma,
  `binom`, and `poch` kernels. The public Fortran API corresponds to SciPy's
  floating `exact=False`, `extend='zero'` behavior; arbitrary-precision
  `exact=True` results are intentionally not emulated.
- Added real gamma-combinatoric interfaces `beta`, generalized `binom`, `poch`,
  and integer-order `polygamma` to `scifort_special`, and expanded `betaln` from
  the positive statistical domain to SciPy-compatible real meromorphic
  continuation. The original SciFort implementation uses gamma-ratio identities,
  exact negative-integer pole cancellations, the Hurwitz-zeta polygamma identity,
  recurrence/reflection on the negative axis, and exact reduced trigonometric
  values at integer/half-integer reflection points.
- Added real gamma-family interfaces `gamma`, `rgamma`, `gammasgn`, real-axis
  `loggamma`, and `multigammaln` to `scifort_special`. The original SciFort
  implementation reuses the existing real `gammaln` reflection kernel, evaluates
  gamma/reciprocal-gamma magnitudes in the log domain, preserves signed zeros and
  subnormal values when representable, and matches SciPy 1.17.0 pole conventions.
  The batch also fixes the existing `digamma(-0.0)` sign convention.
- Added real exponential/trigonometric integral interfaces `exp1`, `expi`, `expn`,
  `sici`, and `shichi` to `scifort_special`. `exp1`/`expi`/`expn` use independent
  SciFort series, continued-fraction, and asymptotic kernels; `shichi` uses the
  standard Ei/E1 identities; `sici` adapts the Caltech/JPL MATH77 DSI/DCI
  Chebyshev kernels under the already-retained BSD-style license. Regression
  coverage includes high orders, breakpoints, large arguments, elemental arrays,
  IEEE endpoints, and a finite Ei range beyond SciPy 1.17.0 premature overflow.
- Added real Jacobi elliptic `ellipj(u,m,sn,cn,dn,ph)` to `scifort_special`,
  returning `sn`, `cn`, `dn`, and the Jacobi amplitude for `0 <= m <= 1`.
  The original SciFort implementation uses the DLMF arithmetic-geometric-mean
  inversion with period reduction, stable exact endpoint formulas, and a
  cancellation-resistant `dn` identity near `m = 1`. Regression coverage
  includes SciPy 1.17.0 values, quarter-period/incomplete-integral identities,
  elemental arrays, large arguments, invalid domains, and IEEE endpoints.
- Added real Legendre elliptic integrals `ellipk`, `ellipkm1`, `ellipe`,
  `ellipkinc`, and `ellipeinc` to `scifort_special`. The implementation is
  original SciFort code using DLMF Carlson connection formulas over the
  existing RC/RD/RF/RG kernels, with a convergent near-singularity series
  for complete E, stable amplitude reduction for arbitrary real phi, and
  explicit IEEE/domain behavior. Regression coverage includes SciPy 1.17.0
  values, tiny complementary parameters, extreme negative parameters,
  periodic identities, huge amplitudes, elemental arrays, and IEEE limits.
- Added real Carlson symmetric elliptic integrals `elliprc`, `elliprd`, `elliprf`,
  `elliprg`, and `elliprj` to `scifort_special`. The RC/RD/RF/RJ duplication
  kernels adapt Jacob Williams' BSD-3-Clause `carlson-elliptic-integrals`
  release 1.0.0 / public-domain SLATEC algorithms; SciFort adds pure elemental
  wrappers, homogeneity scaling, IEEE/domain handling, real principal-value
  continuations for negative RC/RJ pole arguments, and `elliprg` from Carlson's
  standard identity. Regression coverage includes SciPy values, identities,
  singular cases, elemental arrays, and `1e+/-200` scaling.
- Completed the real Airy compatibility slice with `itairy`, `ai_zeros`, and
  `bi_zeros`. Airy integrals use an independently derived integrated ODE series,
  moderate-range Gauss-Legendre bridge, and DLMF tail expansions; zero tables use
  DLMF asymptotic seeds plus Newton refinement through SciFort's `airy` kernel.
  Regression coverage includes signed-limit symmetry, large positive/oscillatory
  arguments, elemental invocation, invalid counts, and the first five
  zeros/cross-values of Ai/Ai'/Bi/Bi'.
- Added SciPy-style real Airy `airy` and exponentially scaled `airye` interfaces. The implementation is independently derived from DLMF Airy/Bessel identities and the Airy ODE, reuses SciFort's existing real-order Bessel kernels away from zero, and uses a local power series near zero. Regression coverage includes the series/Bessel boundary, oscillatory and exponential regions, scaled large arguments, IEEE values, elemental arrays, and a binary64 extreme-argument asymptotic path.
- Added SciPy-style real-axis Bessel derivatives `jvp`, `yvp`, `ivp`, and `kvp` with optional nonnegative derivative order `n` (default 1). The implementation uses DLMF shifted-order identities over SciFort's existing real-order Bessel kernels, normalized binomial weights, and scaled accumulation to avoid intermediate coefficient overflow. Regression coverage includes higher derivatives, negative orders/arguments, zero and IEEE endpoints, elemental arrays, and a 1000th-derivative stress case.
- Added SciPy-style real-order cylindrical Bessel `jv`, `jve`, `yv`, and `yve` to `scifort_special`, plus real-order `jn` as the SciPy-compatible alias of `jv`. The positive-order core adapts the permissively licensed MATH77 `DBESJN`, `DBESYN`, and `DBESPQ` algorithms to modern elemental Fortran, while negative orders use stable reflection identities. Regression coverage includes small/middle/asymptotic regions, near-integer orders, negative arguments, overflow/underflow edges, IEEE values, aliases, and elemental arrays.
- Added SciPy-style real spherical Bessel functions `spherical_jn`, `spherical_yn`,
  `spherical_in`, and `spherical_kn` to `scifort_special`, including optional
  first derivatives. The implementation uses DLMF series, recurrences, and
  explicit half-integer identities with signed-log evaluation for difficult
  modified-K branches; regression coverage includes negative arguments,
  tiny derivatives, high-order finite values, IEEE endpoints, and elemental arrays.
- Added SciPy-style ordered Bessel interfaces `yn`, `kn`, `iv`, `ive`, `kv`, and `kve` to `scifort_special`. Integer-order Y uses the Fortran standard intrinsic; ordered I/K reuse SciFort's existing scaled-I/log-K kernels, including the standard connection formula for negative noninteger I order. Added domain, IEEE, sign-change, scaling, and elemental-array regressions against SciPy 1.17.0.
- Added SciPy-style fast real Bessel functions `j0`, `j1`, `y0`, `y1`, `i0`,
  `i1`, `k0`, `k0e`, `k1`, and `k1e` to `scifort_special`. J/Y use the
  Fortran standard intrinsics; I/K reuse SciFort's existing scaled-I and log-K
  kernels, with regression coverage for domains, IEEE values, arrays, large
  arguments, and scaled/unscaled behavior.
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
- Added `scifort_association_extended` with Kendall tau-b/tau-c, weighted Kendall
  tau, point-biserial correlation, SciPy-style least-squares `linregress`, Theil-Sen and Siegel
  robust slope estimators, Brunner-Munzel, and Page's ordered-alternative trend
  test. Kendall includes exact untied inversion-count probabilities; weighted tau
  supports SciPy-style default/directional/index/custom importance ranks, additive
  or multiplicative exchange weights, and custom rank weighers. Page exact p-values
  use convolution of the one-row permutation distribution. Added
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

- Added `wright_bessel` and `log_wright_bessel` for SciPy-compatible nonnegative real
  parameters. The implementation evaluates the DLMF defining series around its largest
  term in the log domain, with an asymptotic saddle branch for exceptionally large modes.
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

## 2026-10-22: ndimage Fourier filters and per-axis extension modes

- Added `fourier_gaussian_nd`, `fourier_uniform_nd`, and `fourier_shift_nd`
  with full FFT and packed real FFT conventions, scalar/vector parameters,
  fractional Fourier shifts, and one-based FFT axis arguments.
- Added `modes(:)` to Gaussian, uniform, min/max, full N-D correlation,
  convolution, and generic neighborhood filters. Existing `mode=` calls remain
  valid; passing both returns invalid-input status.
- `generic_filter_nd` now supports complex-valued Fortran callbacks with
  complex input/output arrays.
- Added two numerical reference scripts (19,584 mixed-mode scalar checks,
  402 Fourier-domain complex coefficients) and ninth native ndimage test.

## 2026-11-01 — FITPACK TCK compatibility and optional adaptive-knot backend

- Add dependency-free `fitpack_tck_type`, `splrep` for s=0 and task=-1,
  `splev`, `splint`, `splder`, `splantider`, and `spalde`.
- Add real SciPy FITPACK numerical reference cases and native tests.
- Add optional pinned `perazz/fitpack` adaptive-knot `UnivariateSpline`
  adapter under `optional/fitpack_backend`, retaining BSD-3-Clause notices.
  Upstream build/link remains to be tested; no new mandatory dependency.
- Continue prioritizing new SciPy algorithms; Delaunay is unchanged.

## 2026-11-01 — Parametric FITPACK-compatible curve interpolation

- Add dependency-free, nonperiodic `splprep` with SciPy-style zero-residual
  interpolation and fixed-knot multichannel weighted least squares.
- Add automatic chord-length parameterization and `splev` evaluation/derivatives
  for parametric TCK curves, including the FITPACK four extrapolation modes.
- Include native validation and independent SciPy coefficient/value comparisons.
- Preserve positive-s and periodic smoothing as unsupported in the core rather
  than incorrectly substituting separate univariate smoothing splines.
- Make no changes to Delaunay algorithms or FPM dependencies.

## 2026-11-01 — Verified user-supplied FITPACK 2.0.0 offline

- Compiled all 14 core and 26 C-API modules from `fitpack-main.zip` with GNU
  Fortran 14.2.0; upstream Fortran/C++ tests pass 66/66 in checked mode.
- Compiled and tested the real optional SciFort adapter in checked and -O3
  modes, with 140 independent TCK and smoothing-residual reference cases each.
- Added `check_optional_fitpack_local.py`, which builds and tests directly
  from a ZIP or extracted source tree without FPM or network access.
- Corrected the optional SciPy comparison procedure to distinguish different
  adaptive knot selections from errors in SciFort's TCK conversion.
- Kept FITPACK external: no upstream numerical source was copied into SciFort;
  the existing optional FPM tag remains unchanged pending separate verification.

## 2026-11-01 — FITPACK 2.0.0 periodic and parametric smoothing

- Extend the verified optional FITPACK adapter with periodic scalar smoothing,
  adaptive-knot open parametric curves and closed parametric curves.
- Handle the final-ordinate conventions of SciPy's periodic FITPACK interfaces.
- Add checked/optimized native tests and independent SciPy-BSpline reconstruction,
  derivative, residual and direct-fitting comparisons across two random seeds.
- Add an offline FPM path-dependency setup/restoration helper for supplied
  FITPACK 2.0.0 sources without changing the main SciFort dependencies.
- Document positive-s differences from SciPy, including large extrapolation
  differences in certain higher-degree smoothing cases.
- Leave the existing Delaunay implementation unchanged.

## 2026-11-10 — Differential evolution integer constraints and polishing hook

- Add integral coordinates, linear and nonlinear constraints, feasibility-first
  trial selection, constraint evaluation counters, and infeasible-status reporting.
- Add optional polishing callback plus L-BFGS-B bridge using the existing SciFort
  optimizer; polishing is restricted to unconstrained cases.
- Add native tests, randomized comparisons against SciPy, example, and API notes.
