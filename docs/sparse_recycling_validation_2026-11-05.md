# 2026-11-05 — LGMRES and GCROT(m,k) sparse Krylov validation

## Implementation

`scifort_sparse_linalg.lgmres` augments restarted, right-preconditioned flexible
GMRES with normalized correction directions retained from previous cycles or
passed across related solves. `scifort_sparse_linalg.gcrotmk` implements
GCROT(m,k) with a projected Arnoldi iteration and `C=A*U`, `C^T*C=I` recycling.
The retained `U` vectors are reimaged and reorthogonalized at the start of a
new solve, so `recycle_u` is useful even when the operator has changed.

Both solvers share the original SciFort CSR validation, right Jacobi/ILU(0)
application and explicit true-residual stopping rule. Their Arnoldi routine
uses two modified Gram-Schmidt passes and streamed Givens rotations for
least squares. Fortran integer column indices remain one-based. Only
real binary64 CSR is supported; neither complex matrices nor matrix-free
preconditioner callbacks are in the new public API. GCROT currently supports
FIFO (`oldest`) truncation, not SciPy's SVD-based `smallest` truncation.

These are original SciFort MIT implementations based on published Krylov
subspace algorithms. No SciPy source code was translated or imported.
SciPy 1.17.0 was used solely to compare independent numerical results.

## Commands and environment

- GNU Fortran 14.2.0, with `-std=f2018`.
- Checked: `-O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow`.
- Optimized: `-O3`.
- NumPy and SciPy 1.17.0 for independent reference comparisons.

From the repository root:

```bash
python tools/sparse_reference/compare_sparse_recycling_scipy.py --mode both --systems 180
python tools/sparse_reference/run_native_sparse.py --mode both
```

The comparison harness creates 180 random sparse nonsymmetric systems of
order 2–45 **per build** (some with extra upper-banded structure and uneven
row scaling). It solves **two related right-hand sides per system** to
exercise reuse of `outer_v` and `recycle_u`, alternating solvers, with zero
and nonzero initial guesses, Jacobi, ILU(0) and no preconditioner. Each
solution is independently checked against `numpy.linalg.solve` and the
matching SciPy `lgmres`/`gcrotmk` solve, as well as the unpreconditioned
residual and SciFort's returned residual norm. SciPy preconditioning is
omitted from the independent comparison: the Fortran runs use its own
right preconditioners and must still solve the same original problem.

| Check | Checked | Optimized |
|---|---:|---:|
| Randomized systems | 180 / 180 | 180 / 180 |
| Fortran solutions checked against NumPy and SciPy | 360 / 360 | 360 / 360 |
| Standalone sparse native regression executables | 30 / 30 | 30 / 30 |

The focused native test `test/test_sparse_recycling.f90` includes an
indefinite, nonsymmetric, nonnormal system, narrow-restart convergence,
Jacobi and ILU(0), zero-recycling operation, exact starting solutions, empty
matrices, invalid inputs, bad preconditioner diagonals, NaN RHS/history,
zero matrices, a capped nonconvergent solve, correction-vector augmentation,
and recycled GCROT projections after changing the operator. The latter two
exercises demonstrate that the recycled vectors really affect the solver;
they do not merely warm-start an otherwise plain restarted GMRES.

`fpm test` and integration tests requiring external FPM dependencies were
not executed because `fpm` is unavailable in this environment. The
standalone native runner covers the SciFort sparse dependency subset.

## Limitations and interpretation

`maxiter` counts outer/restart cycles, whereas `iterations` counts Arnoldi
columns generated (including augmented LGMRES columns). The positive
nonconvergence status is not a guarantee of a useful approximate solution;
verify `residual_norm`. LGMRES may stagnate with an excessively small
`inner_m`; GCROT can require additional storage and products for the
recycling basis. The preconditioner application order and iteration counts
can differ from SciPy, which uses a more general linear-operator API. No
claims of speed advantage or general convergence on arbitrary nonsymmetric
matrices are made by the randomized tests.
