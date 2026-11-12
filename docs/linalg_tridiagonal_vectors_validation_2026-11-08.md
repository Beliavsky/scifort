# Selected tridiagonal eigenvector validation — 2026-11-08

## Implementation

`eigh_tridiagonal` now uses Sturm bisection and tridiagonal inverse
iteration for selected eigenvectors. Two-pass reorthogonalization keeps the
selected columns orthonormal, including repeated eigenvalues. `auto`
falls back to full-QL on nonconvergence for order <= 256; `inverse`
never does so.

## Reproduction

```sh
python tools/linalg_reference/run_tridiagonal_vectors_checks.py --mode both \
    --cases 600 --max-order 512 --seed 15803
python tools/linalg_reference/run_tridiagonal_selection_checks.py --mode both \
    --cases 720 --max-order 128 --seed 20261108
```

Each `both` run builds with GNU Fortran 14.2.0 in runtime-checked and
`-O3` configurations, and compares results against SciPy 1.17.0.

## Results

| Suite | Checked | Optimized | Largest normalized residual |
|---|---:|---:|---:|
| Existing selection harness, n <= 128 | 720/720 | 720/720 | 2.293e-12 |
| New inverse harness, n <= 128, seed 1441 | 280/280 | 280/280 | 1.665e-12 |
| New inverse harness, n <= 512, seed 15803 | 600/600 | 600/600 | 2.946e-12 |

The new harness checks eigenvalue differences, eigen-equation residuals,
mutual orthogonality, and sign-invariant agreement for well-separated simple
eigenvectors. The order-10000 native Toeplitz test selects and verifies five
eigenvectors without allocating an order-10000 dense matrix. Other native
checks cover exact repeated diagonal eigenvalues, explicit QL fallback,
invalid `eigenvector_method`, and undersized output buffers.

## Limitations

- Bisection is limited to the existing relative-scale tolerance, and may
  not resolve eigenvalue differences substantially below machine precision.
- Strong eigenvalue clusters can require substantial reorthogonalization.
- For a large, especially difficult problem, inverse iteration may fail
  explicitly rather than return an inaccurate eigenvector.
- The general Hermitian banded eigensolver still constructs a dense matrix;
  this improvement applies only to symmetric tridiagonal selected eigenpairs.
- The entire FPM suite was not executed because FPM is unavailable here.
