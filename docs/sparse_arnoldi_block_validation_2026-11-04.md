# Block restarted Arnoldi: numerical validation (2026-11-04)

## Change

`src/sparse_linalg.f90` adds an optional `block_size` argument to both
`eigs_arnoldi` and the `eigs` dispatcher (only for `method='arnoldi'` and
`method='arnoldi_schur'`). The default `block_size=1` preserves the existing
single-vector Arnoldi basis construction and residual restart. With
`block_size=b>1`, the initial subspace contains `b` deterministic independent
complex starting directions. Each new Arnoldi sweep applies the operator to
an earlier column `j-b+1`, providing a block Krylov expansion. On restart,
`b` independent residual/guard directions are retained rather than combining
all residuals into one vector. Two-pass orthogonalization remains enabled.

The Schur-locking path can be combined with block Arnoldi; its invariant-
subspace leakage checks, coupling to the unlocked Schur complement, and
original-matrix eigenpair residual checks are unchanged. The complex
shift-invert operator is supported using SciFort's existing restarted GMRES.
The numerical method remains **thick restarted**, *not* an implicitly restarted
Arnoldi method (IRAM).

`block_size` must be between 1 and `k`; for `b>1` a reduced subspace needs
`ncv >= k+b+1` so that space remains for wanted Ritz vectors, a restart
block and an expansion direction. A full-size `ncv=n` search needs no
restart and may be smaller than that bound. The `subspace` method of `eigs`
rejects `block_size` rather than silently ignoring it. A nonconvergent full
search reports a nonzero `info`, not success.

## Tests

All tests used GNU Fortran checked (`-O0 -g -fcheck=all`) and optimized (`-O3`)
builds, NumPy 2.3.5 and SciPy 1.17.0.

- `test/test_sparse_arnoldi_block.f90`: exact semisimple double real and
  complex eigenvalues, clustered triplet, real shift-invert, nonnormal
  Schur-locked restarts, direct and dispatched modes, independent right
  eigenvectors, invalid block widths, exhausted budgets, a full-size
  subspace, and a **1,400-row** CSR problem with a repeated dominant
  eigenvalue. It independently computes original-matrix residuals.
- `tools/sparse_linalg_reference/compare_arnoldi_block_scipy.py`:
  **108 + 63 = 171** independently generated matrices per build, using
  seeds 20261104 and 20261105, for **342 comparisons across both builds**.
  Nine scenario families cover repeated positive/negative real eigenvalues,
  tight real triplets, repeated imaginary eigenvalues, real/complex shift-
  invert, and nonnormal Schur locking. With these seeds, **149** cases per
  build required more than one Arnoldi cycle. All passed. The tests match
  returned eigenvalues against `scipy.linalg.eigvals`, and validate true
  eigenpair residuals, one-based-friendly eigenvector shapes, unit norms,
  and independent eigenvectors for semisimple repeated eigenvalues.
- The full standalone native sparse runner passed **all 14 tests** in both
  checked and optimized builds, including the new block regression.
- Pre-existing reference comparison scripts still passed in both modes:
  28 random Schur-locked cases and 12 multi-lock cases, plus 18 randomized
  ordinary Arnoldi cases and 4 clustered-spectrum cases per build.

Reproduce the new comparisons with:

```bash
python tools/sparse_reference/run_native_sparse.py --mode both
python tools/sparse_linalg_reference/compare_arnoldi_block_scipy.py --mode both --cases 108 --seed 20261104
python tools/sparse_linalg_reference/compare_arnoldi_block_scipy.py --mode both --cases 63 --seed 20261105
```

## Known limitations

- Independent seeds improve representation of repeated *semisimple*
  eigenspaces but do not guarantee completeness for defective/Jordan
  eigenvalues, or reliable eigenvalue accuracy for highly nonnormal systems.
  Small eigenpair residuals may coexist with poorly conditioned eigenvalues.
- Block size greater than one incurs more expansion/restart work and may
  be slower than the default for easy, well-separated eigenvalues.
- The solver still has no robust, general Schur reordering or IRAM implicit
  shift mechanism, and uses SciFort's inner complex GMRES for shift-invert.
- This update does **not** copy third-party numerical source; SciPy is used
  as a test-only reference. The pre-existing MIT/BSD license notices remain.
- The complete FPM project suite was **not run** because `fpm` was not
  available in the test environment. The pre-existing executable-stack
  warning caused by the LSMR callback design remains.
