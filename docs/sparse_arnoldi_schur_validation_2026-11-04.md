# Nonsymmetric Arnoldi Schur invariant-subspace locking — 2026-11-04

## Algorithm and compatibility

`eigs_arnoldi(...,schur_locking=.true.)` and the `eigs` method selector
`'arnoldi_schur'` opt into a frozen orthonormal Schur subspace. A selected
Ritz vector must first pass the original-matrix true residual check and then
pass a separate active-operator invariant-subspace leakage check. Nearly
linearly dependent vectors are not frozen. A small complex Schur QR rotation
triangularizes the projected operator on the frozen span. Subsequent Arnoldi
search directions are projected away from the frozen span, but the locked
block remains inside the *full* Rayleigh--Ritz extraction with its
nonsymmetric coupling to the orthogonal complement preserved. Inverse-mode
locking uses the shifted inverse as its active operator, and returned
residuals are always checked against the original matrix.

The legacy `eigs` method, original soft-locked `eigs_arnoldi`, symmetric
`eigsh_lanczos`, and all existing one-based CSR storage rules are unchanged.
No external library or vendored source is used.

## Reproducible checks

```sh
python tools/sparse_reference/run_native_sparse.py --mode checked
python tools/sparse_reference/run_native_sparse.py --mode optimized
python tools/sparse_linalg_reference/compare_schur_locked_scipy.py --mode both --cases 91 --seed 20261104
python tools/sparse_linalg_reference/compare_schur_locked_scipy.py --mode both --cases 91 --seed 20261105
```

Each randomized reference run contains 91 two-eigenpair cases covering
normal and nonnormal real and complex spectra, plus 12 additional
three-eigenpair cases with at least two independently frozen Schur vectors.
Both returned eigenvalues and the *full complex right eigenvectors* are
checked against SciPy's dense eigenvalues and independently computed true
residuals. The direct native test also checks method dispatch,
complex shift-invert, persistent lock count, nonconvergence and bad inputs.

## Numerical boundaries

The method remains a thick-restarted Rayleigh--Ritz eigensolver, **not**
ARPACK/IRAM or a guarantee of convergence for defective or highly nonnormal
matrices. Tied or nearly tied eigenvalues can be selected differently by
independent iterative eigensolvers even if each returned eigenpair is valid;
small eigenpair residuals alone do not certify an accurate eigenvalue in a
nonnormal problem. Approximate invariant subspaces may pass the lock
threshold yet require further refinement at tighter tolerances. The opt-in
mode can spend more operator products than soft-locking on easy problems,
because it verifies and rotates each new locked Schur block.

The project-wide `fpm test` may still be unavailable without FPM. GNU
Fortran may emit the pre-existing executable-stack linker warning from
nested LSMR callbacks; no such issue is introduced by the new Arnoldi
code path.

## Results in this snapshot

- Randomized seed `20261104`: 91 two-eigenpair cases **per build** passed,
  including 49 cases with frozen nonempty Schur subspaces and 85 requiring
  more than one restart; an additional 12 three-eigenpair cases with two
  Schur locks passed per build.
- Independent seed `20261105`: 91 further cases per build passed, including
  50 with nonempty locks and 87 requiring restarts; the 12 additional
  three-eigenpair cases again passed.
- Total: **206 distinct eigenproblems per build** (182 two-eigenpair and
  24 three-eigenpair problems), in each of checked and optimized builds.
- The standalone native test additionally exercises a 900-row CSR matrix,
  complex-shift dispatcher, malformed inputs and explicit nonconvergence.
- Existing soft-locked Arnoldi randomized reference tests also pass in both
  builds (20 randomized matrices and four clustered-spectrum cases per mode).
