# Hard-locking validation (`scipy.sparse.linalg.eigsh` subset)

Added `eigsh_lanczos(...,hard_locking=.true.)` and
`eigsh(...,method='lanczos_locked')` with original-operator residual checks,
intermediate hard locks and explicit symmetric-subspace deflation. All input
and output indices remain one-based. `block_size=2` on the direct hard-locked
interface also permits a block Krylov space to resolve a double eigenvalue
that the single-start soft-locking baseline misses.

## Reproduce

```sh
python tools/sparse_reference/run_native_sparse.py --mode checked
python tools/sparse_reference/run_native_sparse.py --mode optimized
python tools/sparse_linalg_reference/compare_lanczos_locked_scipy.py --mode both --cases 80
python tools/sparse_linalg_reference/compare_lanczos_locked_shift_scipy.py --mode both --cases 36
```

The native `test_sparse_lanczos_locked` program checks full-subspace
selection, an 80-row matrix with an early hard-locked eigenpair and a
second slowly converging Ritz pair, real shift-invert selection, zero
iteration budget, singular inner solves, invalid symmetry, invalid input,
Gram-matrix orthogonality, and exact double-eigenvalue recovery with
`block_size=2`.

The independent Python oracles compare returned eigenvalues with dense
`numpy.linalg.eigvalsh` and sparse `scipy.sparse.linalg.eigsh`, and check
*actual* `||A v-lambda v||` independently of solver-reported residuals.
The shift tests compare selection according to `mu=1/(lambda-sigma)`
with `which=LM`, `LA`, `SA`, `SM`. They also verify orthonormality,
convergence statuses, and cycle counts.

## Results verified on 2026-11-04

- All **12** standalone native sparse regressions passed with `-O0 -g
  -fcheck=all -fbacktrace` and `-O3` (GNU Fortran 14.2).
- The new unshifted reference suite passed **180 randomized matrices per
  build**, using two independent seeds and including separate difficult
  and exact-double-eigenvalue cases (80+100 randomized per build).
- The new shift-invert suite passed **36 random shifted matrices per build**
  with NumPy/SciPy comparisons for LM, LA, SA and SM. Additional shorter
  shifted suites passed as well.
- The original Lanczos reference suite passed **24 randomized + 2
  difficult cases per build**, and the earlier restarted shift-invert
  reference suite passed **12 cases per build**.
- One constructed 80-row system returned `info>0`, `locked_count=1`
  after its first cycle and later converged both eigenpairs with hard
  deflation. An exact double-eigenvalue case converged with two seeds
  while the earlier single-start baseline returned only one member of
  that eigenspace.
- The existing LSMR internal-callback executable-stack linker warning
  remains. The full `fpm test` suite was not run (FPM unavailable here).

## Mathematical and API limitations

- Hard deflation is **only** implemented for real symmetric `eigsh`.
  For nonsymmetric and nonnormal matrices, simply removing the span of
  converged right eigenvectors can alter the remaining eigenproblem.
- This is a thick-restarted Rayleigh--Ritz procedure with frozen vectors,
  **not** an implicitly restarted Lanczos (IRLM) or ARPACK clone.
- A single Krylov seed may miss exact eigenvalue multiplicities even when
  every returned pair has a tiny residual. Increasing the deterministic
  `block_size` can help, but no Krylov method can guarantee that an
  arbitrary initial block overlaps every wanted invariant subspace.
- Hard locking does not prove that the requested eigenvalues are the
  global top-k under the chosen ordering; it certifies residuals of the
  selected Ritz pairs. Inexact locking can affect closely spaced values.
- Real shifted MINRES inner solves near a singular shift can fail or
  become costly. These failures are not reported as converged eigenpairs.
- The default soft solver and default `eigsh` dispatcher were left
  unchanged; `block_size` is only valid with the direct hard-locking API.

The method draws on standard thick restarting and locking principles from
Wu and Simon (2000), DOI:10.1137/S0895479898334605. The implementation
contains no imported third-party numerical source.
