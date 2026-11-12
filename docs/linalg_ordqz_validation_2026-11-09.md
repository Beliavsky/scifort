# Ordered generalized Schur validation (2026-11-09)

**Numerical implementation:** `src/linalg_ordqz.f90`. Independent reference:
`scipy.linalg.ordqz(A,B,output='complex')` (SciPy 1.17.0). The new code and
reference driver use only the SciFort project dependencies already present.
The Fortran compiler is GNU Fortran 14.2.0, using `-std=f2018`.

## Completed comparisons

- `--mode both --cases 330 --max-order 25 --seed 80515`:
  **330/330 cases in each build**. Maximum normalized reconstruction error
  `1.16e-14`, maximum projective generalized-eigenvalue difference `3.33e-14`.
- `--mode both --cases 500 --max-order 40 --seed 47611`:
  **500/500 cases in each build**. Maximum normalized reconstruction error
  `1.38e-14`, maximum projective difference `7.79e-14`.
- `--mode checked --cases 160 --max-order 128 --seed 69831 --factors-only`:
  **160/160 cases** using QZ factors independently provided by SciPy. Maximum
  normalized reconstruction error `7.30e-15`. The combined run timed out
  before the optimized batch finished.
- `--mode optimized --cases 45 --max-order 128 --seed 69831 --factors-only`:
  **45/45 cases** passed in a separate optimized run. Maximum normalized
  reconstruction error `7.30e-15`.
- `--mode both --cases 250 --max-order 28 --seed 15109`:
  the checked batch **250/250 passed**, with maximum normalized reconstruction
  error `2.03e-14`. The combined command timed out before optimized results;
  it is **not** counted as a two-build pass.
- `--mode both --cases 80 --max-order 20 --seed 20261109`:
  **80/80 cases in each build**. Maximum reconstruction error `1.71e-14`,
  maximum projective eigenvalue difference `1.91e-14`.
- Native `test_linalg_ordqz.f90`, including built-in and callback selection,
  invalid sort and shape handling, infinite eigenvalues and zero-sized
  pencils, passed in both checked and optimized builds. The previous
  `test_linalg_qz.f90`, `test_linalg_qz_two_sided.f90` and
  `test_linalg_generalized_eigen.f90` native tests also passed in both modes.

**Total completed independent comparisons: 2,275**, counting completed
checked and optimized batches separately. The test cases include real,
complex, and mixed-kind matrix pairs; singular `B` matrices with regular
pencils; 1x1 and larger matrices; sparse/triangular inputs; distinct complex
and real eigenvalues; and extreme coefficient scale ratios.

The reference checks both reconstructing equations, orthogonality of the
left/right Schur bases, upper triangularity, the position and count of
selected eigenvalues, and generalized eigenvalue agreement under optimal
matching of homogeneous/projective coordinates. It intentionally does not
compare Q/Z entries or rely on a particular ordering within either selected
or unselected group; these decompositions are nonunique.

## Boundary and failure behavior

A pencil with exact zero diagonal entries in the original `B` may acquire
roundoff-sized nonzero entries after QZ, so the built-in selectors use a
scale-aware numerical-infinity check. This changes the classification of
very small but finite `beta` values relative to naive literal division,
which is a documented difference from SciPy's exact-zero handling.

Adjacent swaps of nearly coincident eigenvalues can be ill-conditioned.
The implementation returns nonconvergence on a failed swap or excessive
reconstruction residual rather than claiming a valid decomposition.
No real quasi-triangular generalized Schur output or condition estimate is
provided. The underlying `qz` factorization remains limited to order 40;
`ordqz_factors` operates on supplied triangular factors beyond this order.

The complete project-wide `fpm test` suite was not run because `fpm` was
unavailable in the execution environment. Focused GNU Fortran builds and
Python/SciPy checks do not replace full FPM integration testing.

## Reproduce

```bash
OPENBLAS_NUM_THREADS=1 python tools/linalg_reference/run_ordqz_checks.py \
    --mode both --cases 330 --max-order 25 --seed 80515

OPENBLAS_NUM_THREADS=1 python tools/linalg_reference/run_ordqz_checks.py \
    --mode optimized --cases 45 --max-order 128 --seed 69831 --factors-only
```
