# Sparse 1-norm estimation and matrix exponential validation

Date: 2026-11-04. Reference: Python 3.13, NumPy 2.3.5, SciPy 1.17.0;
GNU Fortran in both `-O0 -g -fcheck=all` and `-O3` configurations.

## Implemented functionality

- `onenormest(a,t,itmax,v,w,info)` implements a deterministic sign/adjoint
  block iteration, returning a certified induced-1-norm lower bound with
  normalized input `v` and `w=A*v`. For `t >= n`, a one-pass sparse column
  sum computes the exact norm instead of assembling a dense probe block.
- `expm(a,info,drop_tol,block_size,tol)` forms a canonical CSR exponential
  from block-identity applications of SciFort's existing Taylor exponential
  action. Purely diagonal inputs use an O(n+nnz) direct exponential.
  The default cutoff is zero; positive cutoffs deliberately prune entries.

The Fortran code uses one-based canonical CSR matrices and no new external
runtime libraries. It does not implement SciPy's Padé-based sparse `expm`
or a fully general matrix-free `LinearOperator`-based `onenormest`.

## Automated numerical checks

Reproduce from the project root:

```bash
python tools/sparse_reference/run_native_sparse.py --mode both
python tools/sparse_linalg_reference/compare_onenorm_expm_scipy.py --mode both --cases 240 --seed 20261104
python tools/sparse_linalg_reference/compare_onenorm_expm_scipy.py --mode both --cases 240 --seed 20261234
```

Results:

| Check | Checked (`-fcheck=all`) | Optimized (`-O3`) |
|---|---:|---:|
| Complete native sparse regression programs | 17 passed | 17 passed |
| First-seed random exponential/norm comparisons | 240 passed | 240 passed |
| Second-seed random exponential/norm comparisons | 240 passed | 240 passed |
| Exact-norm checks (`t >= n`), both seeds | 106 passed | 106 passed |
| First seed minimum norm-estimate/exact ratio | 0.741 | 0.741 |
| Second seed minimum norm-estimate/exact ratio | 0.377 | 0.377 |
| First seed maximum scaled `expm` deviation | 3.05e-13 | 3.05e-13 |
| Second seed maximum scaled `expm` deviation | 3.17e-12 | 3.17e-12 |

The worst norm estimate is still a valid mathematical lower bound. SciPy's
own randomized estimator also reports valid lower bounds, but may produce
different estimates because it samples probe directions independently.
Every norm estimate is compared against an exact 1-norm and against SciPy's
`onenormest` as an independent bound, with `sum(abs(v)) = 1`, `w = A*v`, and
`sum(abs(w)) = estimate` verified independently.

Every assembled matrix exponential is checked against both
`scipy.linalg.expm` and `scipy.sparse.linalg.expm` (only as independent
numerical references, not as copied implementations). Tolerances use
matrix-scaled maximum absolute error, rather than excessive elementwise
relative accuracy near cancellation zeros. Matrices include diagonal,
non-normal, triangular/nilpotent, skew-symmetric, symmetric, rank-deficient,
shifted, sparse and empty-pattern examples, with duplicate input coordinates
and explicitly stored zeros.

An additional 50-case optimized regression of the pre-existing
`expm_multiply` interface passed after the changes.

The native test additionally covers a 7,000-row diagonal CSR matrix,
where exponential assembly and the `t=n` exact-norm request avoid any dense
n-by-n arrays, plus a 160-row upper-bidiagonal matrix with independently
known exponential coefficients. Invalid inputs, empty matrices, overflow
reporting, positive drop tolerances, and block-size errors are tested.

## Known limitations

- For large matrices whose exponentials are dense, the result naturally
  requires O(n^2) stored entries and up to O(n) block exponential actions;
  use `expm_multiply` instead when only the action is needed.
- The `expm` algorithm uses trace-shifted scaled Taylor series, rather than
  SciPy's rational Padé algorithm. Floating-point roundoff and performance
  may differ, especially for strongly nonnormal or very large matrices.
- `onenormest` is a lower-bound heuristic, not guaranteed to equal the
  induced 1-norm unless all coordinate directions are probed (`t >= n`).
- The full FPM test suite was not run because `fpm` is unavailable in the
  execution environment. GNU ld reports the pre-existing executable-stack
  warning associated with the LSMR callback implementation.
