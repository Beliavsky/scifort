# Sparse direct solver validation (2026-11-05)

## Scope

New `scifort_sparse_linalg` entry points:

- `spsolve_triangular` for CSR/CSC, lower/upper, unit/nonunit diagonals,
  vector RHS and multiple RHS blocks.
- `splu` (partial row pivoting, dynamically growing sparse rows),
  `splu_solve`, `spsolve`, and `valid_splu`.

Both sparse formats are canonical, real binary64 and one-based. LU factors
contain `row_perm` and combined sparse `lu` such that `P*A=L*U`.

## Reproduction

```bash
python tools/sparse_reference/run_native_sparse.py --mode both
python tools/sparse_linalg_reference/compare_sparse_direct_scipy.py --mode checked --cases 240 --seed 20261105
python tools/sparse_linalg_reference/compare_sparse_direct_scipy.py --mode optimized --cases 300 --seed 20261106
```

## Results

- All **18** standalone native sparse regression programs passed with GNU
  Fortran 14.2 under checked (`-O0 -fcheck=all`) and optimized (`-O3`) builds.
- The first 240 randomized SciPy 1.17 comparisons passed in checked mode:
  **80** sparse LU cases and **160** triangular cases; **42** general
  factorizations required nontrivial row permutations.
- Another 300 randomized cases passed in optimized mode: **100** sparse LU
  cases and **200** triangular cases; **51** general factorizations used a
  nontrivial row permutation.
- Independent NumPy reconstruction checks verify `P*A = L*U` (including
  fill-in and row swaps), rather than merely trusting the returned solution.
- All vector RHS cases compare directly with SciPy SuperLU `spsolve` or
  SciPy `spsolve_triangular`. Native tests additionally cover malformed
  dimensions, singular pivots, missing diagonal entries, nonfinite RHS,
  duplicate COO entries, explicitly stored zeros, empty matrices, reusable
  factors, and a 1,500-row tridiagonal sparse system.
- A GNU linker warning that the pre-existing internal LSMR callbacks require
  an executable stack remains unchanged.

## Limitations

This is a basic serial sparse LU factorization with **row pivoting only**.
SuperLU offers further ordering, scaling, supernodal factorization, and
numerical controls. Large fill-in can make this implementation slower or
more memory-intensive than those specialized external solvers. There is no
complex matrix support, transposed factor solve, column permutation, or
iterative refinement yet. A pivot that is nonzero but numerically tiny may be
accepted at the default `pivot_tol=0`; use a positive absolute threshold to
reject it. The test suite does not establish performance parity with SciPy.

No FPM build was executed: FPM is not installed in this environment. The
standalone native sparse suite uses the project's source modules directly.
