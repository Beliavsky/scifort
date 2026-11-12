# Candidate upstream libraries for SciFort `scipy.sparse`

Surveyed 2026-11-04 against the [Fortran-code-on-GitHub catalog](https://beliavsky.github.io/Fortran-code-on-GitHub/#linear-algebra),
repository README and license files. **These are reuse candidates, not
libraries currently imported into `scifort_sparse`.** Retain original notices
and investigate per-file provenance before copying any implementation.

| Repository | Observed project license | Potential use | Considerations |
|---|---|---|---|
| [FSPARSE](https://github.com/jalvesz/FSPARSE) | MIT ([license](https://github.com/jalvesz/FSPARSE/blob/main/LICENSE.md)) | COO/CSR/CSC/ELL/SELL-C, sparse matrix-vector kernels | FYPP preprocessing in source distribution; good reference for extra formats/precision, less need to replace existing CSR core |
| [kiss](https://github.com/ivan-pi/kiss) | BSD-3-Clause | scipy.sparse.linalg iterative methods (CG, BiCG, BiCGSTAB, GMRES, LGMRES, MINRES, etc.) | Candidate for future methods; current FPM manifest links BLAS/LAPACK, so not imported in the dependency-free CSR layer |
| [LSQR](https://github.com/jacobwilliams/LSQR) | BSD-style ([license](https://github.com/jacobwilliams/LSQR/blob/master/LICENSE.txt)) | Sparse rectangular least squares | Has FPM manifest; adapt callback/solver API to SciFort CSR |
| [LSMR](https://github.com/jacobwilliams/LSMR) | BSD-2-Clause | Sparse regularized least squares | Has FPM manifest; check Stanford upstream notices |
| [LUSOL](https://github.com/jacobwilliams/lusol) | MIT or BSD ([license](https://github.com/jacobwilliams/lusol/blob/master/LICENSE.md)) | Sparse LU and linear-system solving | Solver is more complex than a sparse storage layer; review upstream provenance |
| [SPOOF](https://github.com/samharrison7/spoof) | MIT | CSR and diagonal matrix types | Project says it draws on SPARSKIT; review each source file before transplanting |
| [SPRAL](https://github.com/ralna/spral) | BSD (repository) | LSMR, scaling, indefinite factorizations, matrix-free eigensolvers | Heavy native dependencies (BLAS/LAPACK, METIS/hwloc); optional backend rather than core dependency |
| [Sparse iterative parallel solver](https://github.com/tueboesen/Sparse-iterative-parallel-linear-solver) | MIT | Block-splitting/OpenMP iterative solver | Specialized algorithm and system assumptions; not necessary for baseline SciPy feature parity |
| [sparse_fortran](https://github.com/Euler-37/sparse_fortran) | No license identified in current repository listing | COO/CSR examples | **Do not copy** until an explicit compatible license is confirmed |

## Suggested integration order

1. Add a thin `scifort_sparse_linalg` CSR/matrix-free interface that takes
   a matvec callback or `csr_matrix`; validate with native numerical tests.
2. Integrate or independently implement iterative solvers (`kiss`) with
   license text preserved and algorithm-level provenance. Keep optional
   preconditioners and stopping criteria testable independently.
3. Add LSQR/LSMR to support rectangular least-squares problems, then evaluate
   sparse LU via LUSOL as an optional backend.
4. Evaluate alternative data formats (ELL, BSR, SELL-C) only when new
   functionality or application performance justifies the extra API surface.

All new source code in the 2026-11-04 CSC/diagonal update is independent.
A project's top-level MIT/BSD license does not automatically cover imported
third-party routines; inspect original copyrights, notices, redistribution
conditions, tests, and dependencies before including code. Legal attribution
is distinct from scholarly citation: preserve the license notice in
`LICENSES/` or `THIRD_PARTY_LICENSES.md`, record exact repository/commit and
modifications in `CODE_PROVENANCE.md`, and acknowledge upstream authors.

## 2026-11-04 solver module update

SciFort now includes `scifort_sparse_linalg` with original CG, BiCGSTAB,
GMRES and a thin adapter reusing the existing in-tree BSD-licensed LSMR.
`kiss` has **not** been imported; its upstream FPM manifest specifies
`link = ["blas", "lapack"]`. Attribution for the reused in-tree LSMR
implementation is retained at `src/optimize/lsmr.f90` and in
`CODE_PROVENANCE.md`. See [`sparse_linalg.md`](sparse_linalg.md).

## 2026-11-04 — LSQR / ILU(0) follow-up

The updated `src/sparse_linalg.f90` adds native, dependency-free LSQR and
ILU(0), with optional ILU-preconditioned BiCGSTAB and GMRES. These are
independent MIT-licensed SciFort implementations of published numerical
methods, **not** imported code from `jacobwilliams/LSQR`, `ivan-pi/kiss`,
`LUSOL`, or SuperLU. Reference tests use SciPy's public LSQR outputs and a
separate Python implementation of zero-fill ILU to verify factors and solves.
The upstream candidates remain useful for future solvers and more general
sparse factorization, but are not required by this addition.
