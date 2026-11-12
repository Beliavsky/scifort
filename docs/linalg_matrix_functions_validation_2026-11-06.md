# Dense sqrtm/signm validation — 2026-11-06

Compiler: GNU Fortran (checked: `-O0 -fcheck=all`; optimized: `-O3`).
Reference: NumPy 2.3.5, SciPy 1.17.0.

The randomized tests exercise real and complex square matrices of orders 1–8:
positive-definite real/complex matrices, indefinite real symmetric matrices,
triangular matrices with positive and negative diagonals, and real/complex
nonsymmetric matrices. `sqrtm` is checked by `X@X≈A`, and `signm` by
`S@S≈I`, independent of the external reference value. Inputs and outputs are
compared in double precision.

| Run | Cases | Result |
|---|---:|---|
| Checked initial run, seed 20261106 | 300 | 300 passed |
| Checked, seed 20261108 | 700 | 700 passed |
| Checked, seed 20261110 | 700 | 700 passed |
| Optimized, seed 20261109 | 1,600 | 1,600 passed |
| Optimized, seed 20261110 | 700 | 700 passed |
| **Total** | **4,000** | **4,000 passed** |

The largest relative error in the 1,600-case optimized batch was
`7.626e-14` and largest relative square-root/sign-involution residual was
`1.432e-14`. The 700-case checked batch had largest error `7.403e-14`.

One `sqrtm` branch-cut discrepancy was identified on a real symmetric matrix
with exact negative eigenvalues. SciPy's Schur QR selected negative imaginary
roots based on insignificant `O(1e-16)` imaginary eigenvalue roundoff. The
Fortran result selected positive imaginary roots and had a true residual of
`5.13e-15`. For real symmetric matrices, comparison therefore uses an
independent `numpy.linalg.eigh` principal-branch reference.

The focused Fortran test `test/test_linalg_matrix_functions.f90` passed in
both compiler modes. Previous dense Sylvester/Lyapunov and CARE/DARE native
tests also passed in both configurations. The new source is usable without
BLAS/LAPACK. The complete project-wide `fpm test` suite was not executed.
