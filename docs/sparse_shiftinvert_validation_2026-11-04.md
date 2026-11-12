# Sparse shift-invert extension: validation (2026-11-04)

This update adds **complex-shifted nonsymmetric eigenpairs** to `eigs`,
reduced-subspace zero-shift `SM` selection, and the transformed `LA`, `SA`,
`SM` criteria to the existing symmetric `eigsh(...,sigma=...)` path.
It does **not** implement implicitly restarted Arnoldi or Lanczos.

## Algorithms and semantics

- `eigs` still uses orthogonal block subspace iteration and the existing
  projected complex Schur QR eigensolver. For a supplied complex `sigma`,
  its block update applies `(A-sigma*I)^(-1)` using a new private complex
  restarted GMRES with twice modified Gram-Schmidt, complex Givens rotations,
  and true-residual checks. `which='LM'` orders by distance to the shift;
  values and residuals refer to the *original* sparse matrix. Complex
  shifts use a **full complex inverse**, unlike SciPy's default real-matrix
  `OPpart='r'` transform; the comparison casts the SciPy matrix to complex.
- Reduced-subspace nonsymmetric `which='SM'` now implicitly applies the
  inverse at zero for nonsingular matrices; when `ncv=n`, no inverse is
  needed. Singular operators may cause inner breakdown at zero shift.
- Symmetric `eigsh` with a shift selects by the transformed eigenvalues
  `1/(lambda-sigma)`: `LM` closest, `LA` just above, `SA` just below, and
  `SM` farthest. Shifted `SM` requires `ncv=n`.
- No new build dependencies or third-party source are introduced. The
  earlier MINRES BSD-3-Clause notice is retained.

## Environment and reproduced checks

- GNU Fortran 14.2.0, `-std=f2018`; checked (`-O0 -g -fcheck=all`)
  and optimized (`-O3`).
- SciPy 1.17.0 and NumPy for independent eigenvalue/reference comparisons.
- `python tools/sparse_reference/run_native_sparse.py --mode both`: all
  eight native sparse-matrix and sparse-solver programs passed in both modes.
- `python tools/sparse_linalg_reference/compare_eigs_shift_scipy.py --mode both
  --cases 90 --seed 198765`: 90 cases per build, of which 60 have complex
  shifts. All selected eigenvalues agree with SciPy's full-complex inverse
  on these cases and with NumPy dense eigenvalues. True residuals verified.
- The same comparison with `--smallest --cases 65`: 65 reduced-subspace
  nonsymmetric `SM` cases per build, all passed.
- `python tools/sparse_linalg_reference/compare_shiftinvert_modes_scipy.py
  --mode both --cases 84 --seed 198765`: 84 symmetric shifted-mode cases
  per build (`LA`, `SA`, `SM`), all passed against dense `eigvalsh` and SciPy
  (a SciPy reference call initially exhausted its default iteration budget;
  increasing that budget to 2,000 allowed the reference to converge).
- Existing symmetric shift-invert 35 cases per build, original `eigs` 30
  cases per build, original `eigsh` 35 cases per build, original iterative
  sparse-solver/LSQR/ILU/MINRES comparisons passed in both modes.
- Existing sparse-graph, hierarchical-clustering and K-means native
  regression runners passed in both modes.

## Known limitations

- The outer eigenalgorithm remains **block subspace iteration, not implicit
  restarting (IRAM/IRLM)**. Large, highly nonnormal or tightly clustered
  problems may converge slowly or give nonconvergence.
- Restarted inner GMRES without a general complex preconditioner may be
  expensive or fail near singular shifts; failure yields a negative status.
- Current CSR input is real; `eigs` returns complex eigenpairs but cannot
  directly accept complex sparse matrices or generalized mass matrices.
- Symmetric shifted `SM` intentionally requires a full subspace and therefore
  may use quadratic memory. The reduced nonsymmetric `SM` path assumes
  invertible `A`.
- The full FPM suite was not run (FPM not installed). A pre-existing
  executable-stack linker warning from internal LSMR callback trampolines
  remains; it is not caused by these shift-invert changes.
